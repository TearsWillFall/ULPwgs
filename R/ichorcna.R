

#' IchorCNA implementation for Capture Panel Data
#' 
#'
#' @param cnr [REQUIRED] Path to tumour CNR file.
#' @param normal [OPTIONAL] Path to normal samples. Default c(0.5,0.6,0.7,0.8,0.9)
#' @param ploidy [OPTIONAL] Initial tumour ploidy; can be more than one value if additional ploidy initializations are desired. Default: 2
#' @param lambda [OPTIONAL] Initial Student's t precision; must contain 4 values (e.g. c(1500,1500,1500,1500)); if not provided then will automatically use based on variance of data
#' @param scStates [OPTIONAL]  Subclonal states to consider. Default NULL
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
#' @param lambdaScaleHyperParam [OPTIONAL] Hyperparameter (scale) for Gamma prior on Student's-t precision. Default 3
#' @param maxCN [OPTIONAL] Total clonal states. Default 7.
#' @param estimateNormal [OPTIONAL] Estimate normal. Default TRUE.
#' @param estimateScPrevalence [OPTIONAL] Estimate subclonal prevalence. Default TRUE.
#' @param estimatePloidy [OPTIONAL] Estimate tumour ploidy. Default TRUE.
#' @param maxFracGenomeSubclone [OPTIONAL] Exclude solutions with subclonal genome fraction greater than this value. Default 0.5
#' @param maxFracCNASubclone  [OPTIONAL] Exclude solutions with fraction of subclonal events greater than this value. Default 0.7
#' @param minSegmentBins [OPTIONAL]  Minimum number of bins for largest segment threshold required to estimate tumor fraction; if below this threshold, then will be assigned zero tumor fraction
#' @param altFracThreshold [OPTIONAL] Minimum proportion of bins altered required to estimate tumor fraction; if below this threshold, then will be assigned zero tumor fraction. Default: [0.05]
#' @param includeHOMD [OPTIONAL] If FALSE, then exclude HOMD state. Useful when using large bins (e.g. 1Mb). Default FALSE.
#' @param txnE [OPTIONAL] Self-transition probability. Increase to decrease number of segments. Default: [0.9999999]
#' @param txnStrength [OPTIONAL] Transition pseudo-counts. Exponent should be the same as the number of decimal places of --txnE. Default: [1e+07]
#' @param plotFileType [OPTIONAL] File format for output plots. Default pdf
#' @param plotYLim [OPTIONAL] ylim to use for chromosome plots. Default: [c(-2,2)]
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param threads [OPTIONAL] Number of threads to split the work. Default 4
#' @param ram [OPTIONAL] RAM memory to asing to each thread. Default 4
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Task EXECUTOR ID. Default "recalCovariates"
#' @param task_name Task name. Default "recalCovariates"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export



ichor_capture=function(
    cnr="R/TRAILS_TR067_I260038_BULK_CAPTURE_PCF_V2_HG19_LB1_411_HW57JDMXX_1.cnr",
    normal=0.5,
    ploidy=2,
    maxCN=7,
    includeHOMD=FALSE,
    scStates=c(1,3),
    txnE=0.9999999,
    txnStrength=1e7,
    output_name="",
    output_dir=".",
    min_cov=-15,
    lambda=NULL,
    coverage=NULL,
    minSegmentBins=50,
    minTumFracToCorrect=0.01,
    chrs=c(1:22,"X"),
    chrTrain=seq(1:22),
    gender="male",
    maxFracCNASubclone=0.7,
    maxFracGenomeSubclone=0.5,
    lambdaScaleHyperParam=3,
    altFracThreshold=0.05,
    estimateNormal=TRUE,
    estimatePloidy=TRUE,
    estimateScPrevalence=TRUE,
    plotYLim=c(-2,2),
    plotFileType="pdf",
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=4,ram=4,mode="local",
    executor_id=make_unique_id("hybridIchorCNA"),
    task_name="hybridIchorCNA",time="48:0:0",
    update_time=60,
    wait=FALSE,hold=""
){  

    options(scipen=0, stringsAsFactors=F)
    options(stringsAsFactors=FALSE)
    options(bitmapType='cairo')
    argg <- as.list(environment())
    
    id=""
    if(output_name!=""){
      id=output_name
    }else{
      id=get_file_name(cnr)
    }

    out_file_dir=set_dir(dir=output_dir,name=paste0(id,"/ichor_capture"))
    tmp_dir=set_dir(dir=out_file_dir,name="tmp")


    
    outImage <- paste0(out_file_dir,"/",id,".RData")
    names(cnr)=Vectorize(get_file_name)(cnr)
    numSamples <- length(cnr)
    read_tumour_copy=function(file="R/TRAILS_TR067_I260038_BULK_CAPTURE_PCF_V2_HG19_LB1_411_HW57JDMXX_1.cnr",
        header=TRUE,sep="\t"){
        tumour_copy=read.table(file=file,header=header,sep=sep)
        
        output <- GenomicRanges::GRanges(ranges = IRanges::IRanges(start = tumour_copy$start,
        width = tumour_copy$end-tumour_copy$start),
                        seqnames = tumour_copy$chromosome, 
                        log2 = tumour_copy$log2,
                        gene = tumour_copy$gene,
                        depth = tumour_copy$depth,
                        weight = tumour_copy$weight)
        message(paste0("Dropped ",length(output[!as.character(GenomeInfoDb::seqnames(output)) %in% chrs,])," bins with omitted chromosomes."))
        return(output[as.character(GenomeInfoDb::seqnames(output)) %in% chrs,])
    }


    drop_low_coverage=function(tumour_copy,min_cov=-15){
       message(paste0(sum( tumour_copy[[1]]$log2<min_cov | tumour_copy[[1]]$depth==0)," bins excluded from analysis due to low coverage."))
        return(!tumour_copy[[1]]$log2<min_cov | tumour_copy[[1]]$depth==0)
    }

    drop_weight=function(tumour_copy,min_weight=0){
       message(paste0(sum(tumour_copy[[1]]$weight<min_weight),
       " bins excluded from analysis due to low weight."))
        return(!tumour_copy[[1]]$weight<min_weight)
    }

    drop_outliers=function(tumour_copy,width=50,factor=10,partial=1){
        out=as.data.frame(tumour_copy[[1]]) %>% dplyr::group_by(seqnames) %>% 
        dplyr::mutate(smooth_log2=abs(log2-pracma::savgol(log2,
        ifelse(width%%2==0,width+1,width)))) %>% 
        dplyr::mutate(quant=as.vector(zoo::rollapply(smooth_log2, width = width,
         FUN = "quantile", p = .95,fill=FALSE,partial=partial))) %>% 
         dplyr::mutate(filter=abs(smooth_log2)>(quant*factor))
       
        message(paste0("Dropped ",sum(out$filter)," outlier bins."))
        return(!out$filter)
    }

    tumour_copy<-lapply(cnr,FUN=read_tumour_copy)
    low_coverage=drop_low_coverage(tumour_copy)
    weight=drop_weight(tumour_copy)
    outliers=drop_outliers(tumour_copy)


    valid <- as.character(GenomeInfoDb::seqnames(tumour_copy[[1]])) %in% chrs
    valid <- valid & low_coverage & weight & outliers
    chrInd <- as.character(GenomeInfoDb::seqnames(tumour_copy[[1]])) %in% chrTrain
    
  


    
    ### RUN HMM ###
    ## store the results for different normal and ploidy solutions ##

    ptmTotalSolutions <- proc.time() # start total timer
    results <- list()
    loglik <- as.data.frame(matrix(NA, nrow = length(normal) * length(ploidy), ncol = 7, 
                    dimnames = list(c(), c("init", "n_est", "phi_est", "BIC", 
                                                                    "Frac_genome_subclonal", "Frac_CNA_subclonal", "loglik"))))
    counter <- 1
    compNames <- rep(NA, nrow(loglik))
    mainName <- rep(NA, length(normal) * length(ploidy))
    #### restart for purity and ploidy values ####

for (n in normal){
    for (p in ploidy){
        if (n == 0.95 & p != 2) {
            next
        }

        logR <- as.data.frame(lapply(tumour_copy, function(x) { x$log2 })) # NEED TO EXCLUDE CHR X #
        param <- getDefaultParameters(logR[valid & chrInd, , drop=F], maxCN = maxCN, includeHOMD = includeHOMD, 
                    ct.sc=scStates, ploidy = floor(p), e=txnE, e.same = 50, strength=txnStrength)
        param$phi_0 <- rep(p, numSamples)
        param$n_0 <- rep(n, numSamples)
        
        ############################################
        ######## CUSTOM PARAMETER SETTINGS #########
        ############################################
        # 0.1x cfDNA #
        if (is.null(lambda)){
                logR.var <- 1 / ((apply(logR, 2, sd, na.rm = TRUE) / sqrt(length(param$ct))) ^ 2)
                param$lambda <- rep(logR.var, length(param$ct))
                param$lambda[param$ct %in% c(2)] <- logR.var 
                param$lambda[param$ct %in% c(1,3)] <- logR.var 
                param$lambda[param$ct >= 4] <- logR.var / 5
                param$lambda[param$ct == max(param$ct)] <- logR.var / 15
                param$lambda[param$ct.sc.status] <- logR.var / 10
        }else{
                param$lambda[param$ct %in% c(2)] <- lambda[2]
                param$lambda[param$ct %in% c(1)] <- lambda[1]
                param$lambda[param$ct %in% c(3)] <- lambda[3]
                param$lambda[param$ct >= 4] <- lambda[4]
                param$lambda[param$ct == max(param$ct)] <- lambda[2] / 15
                param$lambda[param$ct.sc.status] <- lambda[2] / 10
            }
            param$alphaLambda <- rep(lambdaScaleHyperParam, length(param$ct))  
        # 1x bulk tumors #
        #param$lambda[param$ct %in% c(2)] <- 2000
        #param$lambda[param$ct %in% c(1)] <- 1750
        #param$lambda[param$ct %in% c(3)] <- 1750
        #param$lambda[param$ct >= 4] <- 1500
        #param$lambda[param$ct == max(param$ct)] <- 1000 / 25
            #param$lambda[param$ct.sc.status] <- 1000 / 75
            #param$alphaLambda[param$ct.sc.status] <- 4
            #param$alphaLambda[param$ct %in% c(1,3)] <- 5
            #param$alphaLambda[param$ct %in% c(2)] <- 5
            #param$alphaLambda[param$ct == max(param$ct)] <- 4
                    
            #############################################
            ################ RUN HMM ####################
            #############################################
        hmmResults.cor <- HMMsegment(x=tumour_copy, validInd=valid, dataType = "log2", 
                                    param = param, chrTrain = chrTrain, maxiter = 50,
                                    estimateNormal = estimateNormal, estimatePloidy = estimatePloidy,
                                    estimateSubclone = estimateScPrevalence, verbose = TRUE)
      
        for (s in 1:numSamples){
            iter <- hmmResults.cor$results$iter
            id <- names(hmmResults.cor$cna)[s]

            ## convert full diploid solution (of chrs to train) to have 1.0 normal or 0.0 purity
            ## check if there is an altered segment that has at least a minimum # of bins
            segsS <- hmmResults.cor$results$segs[[s]]
            segsS <- segsS[segsS$chr %in% chrTrain, ]
            segAltInd <- which(segsS$event != "NEUT")
            maxBinLength = -Inf
            if (sum(segAltInd) > 0){
                maxInd <- which.max(segsS$end[segAltInd] - segsS$start[segAltInd] + 1)
                maxSegRD <- GenomicRanges::GRanges(seqnames=segsS$chr[segAltInd[maxInd]], 
                                    ranges=IRanges::IRanges(start=segsS$start[segAltInd[maxInd]], end=segsS$end[segAltInd[maxInd]]))
                hits <- IRanges::findOverlaps(query=maxSegRD, subject=tumour_copy[[s]][valid, ])
                maxBinLength <- length(S4Vectors::subjectHits(hits))
            }
            ## check if there are proportion of total bins altered 
            # if segment size smaller than minSegmentBins, but altFrac > altFracThreshold, then still estimate TF
            cnaS <- hmmResults.cor$cna[[s]]
            altInd <- cnaS[cnaS$chr %in% chrTrain, "event"] == "NEUT"
            altFrac <- sum(!altInd, na.rm=TRUE) / length(altInd)
            if ((maxBinLength <= minSegmentBins) & (altFrac <= altFracThreshold)){
                hmmResults.cor$results$n[s, iter] <- 1.0
            }

        # correct integer copy number based on estimated purity and ploidy
        correctedResults <- correctIntegerCN(cn = hmmResults.cor$cna[[s]],
                segs = hmmResults.cor$results$segs[[s]], 
                purity = 1 - hmmResults.cor$results$n[s, iter], ploidy = hmmResults.cor$results$phi[s, iter],
                cellPrev = 1 - hmmResults.cor$results$sp[s, iter], 
                maxCNtoCorrect.autosomes = maxCN, maxCNtoCorrect.X = maxCN, minPurityToCorrect = minTumFracToCorrect, 
                gender = gender[[s]], chrs = chrs, correctHOMD = includeHOMD)
        hmmResults.cor$results$segs[[s]] <- correctedResults$segs
        hmmResults.cor$cna[[s]] <- correctedResults$cn

            ## plot solution ##
            outPlotFile <- paste0(out_file_dir,  "/", id, "_genomeWide_", "n", n, "-p", p)
            mainName[counter] <- paste0(id, ", n: ", n, ", p: ", p, ", log likelihood: ", signif(hmmResults.cor$results$loglik[hmmResults.cor$results$iter], digits = 4))
            plotGWSolution(hmmResults.cor=hmmResults.cor, s=s, outPlotFile=outPlotFile, plotFileType=plotFileType, 
                logR.column = "logR", call.column = "Corrected_Call",
                        plotYLim=plotYLim, estimateScPrevalence=estimateScPrevalence,
                         seqinfo=NULL, main=mainName[counter])
        }
        iter <- hmmResults.cor$results$iter
        results[[counter]] <- hmmResults.cor
        loglik[counter, "loglik"] <- signif(hmmResults.cor$results$loglik[iter], digits = 4)
        subClonalBinCount <- unlist(lapply(hmmResults.cor$cna, function(x){ sum(x$subclone.status) }))
        fracGenomeSub <- subClonalBinCount / unlist(lapply(hmmResults.cor$cna, function(x){ nrow(x) }))
        fracAltSub <- subClonalBinCount / unlist(lapply(hmmResults.cor$cna, function(x){ sum(x$copy.number != 2) }))
        fracAltSub <- lapply(fracAltSub, function(x){if (is.na(x)){0}else{x}})
        loglik[counter, "Frac_genome_subclonal"] <- paste0(signif(fracGenomeSub, digits=2), collapse=",")
        loglik[counter, "Frac_CNA_subclonal"] <- paste0(signif(as.numeric(fracAltSub), digits=2), collapse=",")
        loglik[counter, "init"] <- paste0("n", n, "-p", p)
        loglik[counter, "n_est"] <- paste(signif(hmmResults.cor$results$n[, iter], digits = 2), collapse = ",")
        loglik[counter, "phi_est"] <- paste(signif(hmmResults.cor$results$phi[, iter], digits = 4), collapse = ",")

        counter <- counter + 1
    }
    }
    ## get total time for all solutions ##
    elapsedTimeSolutions <- proc.time() - ptmTotalSolutions
    message("Total ULP-WGS HMM Runtime: ", format(elapsedTimeSolutions[3] / 60, digits = 2), " min.")

    ### SAVE R IMAGE ###
    save.image(outImage)
    #save(tumour_copy, results, loglik, file=paste0(outDir,"/",id,".RData"))

    ### SELECT SOLUTION WITH LARGEST LIKELIHOOD ###
    loglik <- loglik[!is.na(loglik$init), ]
    if (estimateScPrevalence){ ## sort but excluding solutions with too large % subclonal 
        fracInd <- which(loglik[, "Frac_CNA_subclonal"] <= maxFracCNASubclone & 
                                    loglik[, "Frac_genome_subclonal"] <= maxFracGenomeSubclone)
        if (length(fracInd) > 0){ ## if there is a solution satisfying % subclonal
            ind <- fracInd[order(loglik[fracInd, "loglik"], decreasing=TRUE)]
        }else{ # otherwise just take largest likelihood
            ind <- order(as.numeric(loglik[, "loglik"]), decreasing=TRUE) 
        }
    }else{#sort by likelihood only
    ind <- order(as.numeric(loglik[, "loglik"]), decreasing=TRUE) 
    }

    #new loop by order of solutions (ind)
    outPlotFile <- paste0(out_file_dir,  "/", id, "_genomeWide_all_sols")
    for(i in 1:length(ind)) {
        hmmResults.cor <- results[[ind[i]]]
        turnDevOff <- FALSE
        turnDevOn <- FALSE
        if (i == 1){
            turnDevOn <- TRUE
        }
        if (i == length(ind)){
            turnDevOff <- TRUE
    }
    plotGWSolution(hmmResults.cor, s=s, outPlotFile=outPlotFile, plotFileType="pdf", 
                        logR.column = "logR", call.column = "Corrected_Call",
                        plotYLim=plotYLim, estimateScPrevalence=estimateScPrevalence, 
                        seqinfo = NULL,
                        turnDevOn = turnDevOn, turnDevOff = turnDevOff, main=mainName[ind[i]])
    }

    hmmResults.cor <- results[[ind[1]]]
    hmmResults.cor$results$loglik <- as.data.frame(loglik)
    hmmResults.cor$results$gender <- gender
    hmmResults.cor$results$coverage <- coverage

    outputHMM(cna = hmmResults.cor$cna, segs = hmmResults.cor$results$segs, 
                        results = hmmResults.cor$results, patientID = id, outDir=out_file_dir)
    outFile <- paste0(out_file_dir, "/",id, ".params.txt")
    outputParametersToFile(hmmResults=hmmResults.cor, file = outFile)


    job_report=build_job_report(
        job_id=job,
        executor_id=executor_id,
        exec_code=exec_code,
        task_id=task_id,
        input_args = argg,
        out_file_dir=out_file_dir,
        out_files=list(
            param=outFile,
            plot_solutions=outPlotFile,

        )
    )

    return(job_report)


}



#' IchorCNA implementation for Capture Panel Data
#' 
#'
#' @param cnr [REQUIRED] Path to tumour CNR file.
#' @param normal [OPTIONAL] Path to normal samples. Default c(0.5,0.6,0.7,0.8,0.9)
#' @param ploidy [OPTIONAL] Initial tumour ploidy; can be more than one value if additional ploidy initializations are desired. Default: 2
#' @param lambda [OPTIONAL] Initial Student's t precision; must contain 4 values (e.g. c(1500,1500,1500,1500)); if not provided then will automatically use based on variance of data
#' @param scStates [OPTIONAL]  Subclonal states to consider. Default NULL
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
#' @param lambdaScaleHyperParam [OPTIONAL] Hyperparameter (scale) for Gamma prior on Student's-t precision. Default 3
#' @param maxCN [OPTIONAL] Total clonal states. Default 7.
#' @param estimateNormal [OPTIONAL] Estimate normal. Default TRUE.
#' @param estimateScPrevalence [OPTIONAL] Estimate subclonal prevalence. Default TRUE.
#' @param estimatePloidy [OPTIONAL] Estimate tumour ploidy. Default TRUE.
#' @param maxFracGenomeSubclone [OPTIONAL] Exclude solutions with subclonal genome fraction greater than this value. Default 0.5
#' @param maxFracCNASubclone  [OPTIONAL] Exclude solutions with fraction of subclonal events greater than this value. Default 0.7
#' @param minSegmentBins [OPTIONAL]  Minimum number of bins for largest segment threshold required to estimate tumor fraction; if below this threshold, then will be assigned zero tumor fraction
#' @param altFracThreshold [OPTIONAL] Minimum proportion of bins altered required to estimate tumor fraction; if below this threshold, then will be assigned zero tumor fraction. Default: [0.05]
#' @param includeHOMD [OPTIONAL] If FALSE, then exclude HOMD state. Useful when using large bins (e.g. 1Mb). Default FALSE.
#' @param txnE [OPTIONAL] Self-transition probability. Increase to decrease number of segments. Default: [0.9999999]
#' @param txnStrength [OPTIONAL] Transition pseudo-counts. Exponent should be the same as the number of decimal places of --txnE. Default: [1e+07]
#' @param plotFileType [OPTIONAL] File format for output plots. Default pdf
#' @param plotYLim [OPTIONAL] ylim to use for chromosome plots. Default: [c(-2,2)]
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param threads [OPTIONAL] Number of threads to split the work. Default 4
#' @param ram [OPTIONAL] RAM memory to asing to each thread. Default 4
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Task EXECUTOR ID. Default "recalCovariates"
#' @param task_name Task name. Default "recalCovariates"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export



parallel_sample_ichor_capture=function(
    cnrs="",
    normal=0.5,
    ploidy=2,
    maxCN=7,
    includeHOMD=FALSE,
    scStates=c(1,3),
    txnE=0.9999999,
    txnStrength=1e7,
    output_dir=".",
    min_cov=-15,
    lambda=NULL,
    coverage=NULL,
    minSegmentBins=50,
    minTumFracToCorrect=0.01,
    chrs=c(1:22,"X"),
    chrTrain=seq(1:22),
    gender="male",
    maxFracCNASubclone=0.7,
    maxFracGenomeSubclone=0.5,
    lambdaScaleHyperParam=3,
    altFracThreshold=0.05,
    estimateNormal=TRUE,
    estimatePloidy=TRUE,
    estimateScPrevalence=TRUE,
    plotYLim=c(-2,2),
    plotFileType="pdf",
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=4,ram=4,mode="local",
    executor_id=make_unique_id("hybridIchorCNA"),
    task_name="hybridIchorCNA",time="48:0:0",
    update_time=60,
    wait=FALSE,hold=""
){  

    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir,name="ichor_capture")

    jobs_report=build_job_report(
          job_id=job,
          executor_id=executor_id,
          exec_code=list(), 
          task_id=task_id,
          input_args=argg,
          out_file_dir=out_file_dir,
          out_files=list(
            )
      )
    
    cnr_list=cnr
    names(cnr_list)=Vectorize(get_file_name)(cnrs)

    if(mode=="local"){
        jobs_report[["steps"]][["par_sample_ichor_capture"]]<-
        parallel::mclapply(cnr_list,FUN=function(cnr){
        job_report <-ichor_capture(
            cnr=cnr,
            normal=normal,
            ploidy=ploidy,
            maxCN=maxCN,
            includeHOMD=includeHOMD,
            scStates=scStates,
            txnE=txnE,
            txnStrength=txnStrength,
            output_name=get_file_name(cnr),
            output_dir=out_file_dir,
            min_cov=min_cov,
            lambda=lambda,
            coverage=coverage,
            minSegmentBins=minSegmentBins,
            minTumFracToCorrect=minTumFracToCorrect,
            chrs=chrs,
            chrTrain=chrTrain,
            gender=gender,
            maxFracCNASubclone=maxFracCNASubclone,
            maxFracGenomeSubclone= maxFracGenomeSubclone,
            lambdaScaleHyperParam=lambdaScaleHyperParam,
            altFracThreshold=altFracThreshold,
            estimateNormal=estimateNormal,
            estimatePloidy=estimatePloidy,
            estimateScPrevalence=estimateScPrevalence,
            plotYLim=plotYLim,
            plotFileType=plotFileType,
            verbose=verbose,
            batch_config=batch_config,
            threads=threads,
            ram=ram,mode=mode,
            executor_id=task_id,
            time=time,
            hold=hold
    )
    },mc.cores=threads)
    
    }else if(mode=="batch"){
          rdata_file=paste0(tmp_dir,"/",job,".samples.RData")
          output_dir=out_file_dir
          save(cnr_list,normal,
            ploidy,
            maxCN,
            includeHOMD,
            scStates,
            txnE,
            txnStrength,
            min_cov,
            lambda,
            coverage,
            minSegmentBins,
            minTumFracToCorrect,
            chrs,
            chrTrain,
            gender,
            maxFracCNASubclone,
            maxFracGenomeSubclone,
            lambdaScaleHyperParam,
            altFracThreshold,
            estimateNormal,
            estimatePloidy,
            estimateScPrevalence,
            plotYLim,
            plotFileType,verbose,file = rdata_file)
          exec_code=paste0("Rscript -e \"ULPwgs::ichor_capture(rdata=\\\"",
          rdata_file,"\\\",selected=$SGE_TASK_ID)\"")
          out_file_dir2=set_dir(dir=out_file_dir,name="batch")
          batch_code=build_job_exec(job=job,time=time,ram=ram,
          threads=1,output_dir=out_file_dir2,
          hold=hold,array=length(cnrs_list))
          exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)

          if(verbose){
              print_verbose(job=job,arg=argg,exec_code=exec_code)
          }
          error=execute_job(exec_code=exec_code)
          if(error!=0){
              stop("cnvkit failed to run due to unknown error.
              Check std error for more information.")
          }
         
         jobs_report[["steps"]][["par_sample_ichor_capture"]]<- build_job_report(
              job_id=job,
              executor_id=executor_id,
              exec_code=exec_code, 
              task_id=task_id,
              input_args=argg,
              out_file_dir=out_file_dir,
              out_files=list(
                  param=paste0(out_file_dir,"/",names(cnrs_list),"/ichor_capture/",names(cnrs_list),".params.txt"),
                  plot=paste0(out_file_dir,"/",names(cnrs_list),"/ichor_capture/",names(cnrs_list),"_genomeWide_all_sols.",plotFileType)
              )
        )
    }

      return(jobs_report)



}
