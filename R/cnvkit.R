
#' Variant Calling using Mutect2
#'
#' This function functions calls Mutect2 for variant calling.
#' If a vector of tumour samples are provided these will be processed in multi-sample mode.
#' To run in tumour-normal mode suppply a single tumour and normal sample.
#' If no normal is supplied this will run in tumour only.
#' TO DO// Implement mitochondrial mode feature
#' 
#' For more information read:
#' https://cnvkit.readthedocs.io/en/stable/pipeline.html
#'
#'
#' @param sif_cnvkit [REQUIRED] Path to cnvkit sif file.
#' @param tumour [REQUIRED] Path to tumour BAM file.
#' @param normal [OPTIONAL] Path to normal BAM file.
#' @param ref_genome [REQUIRED] Path to reference genome fasta file.
#' @param rdata [OPTIONAL] Import R data information with list of BAM.
#' @param selected [OPTIONAL] Select BAM from list.
#' @param pon [OPTIONAL] Path to Panel of Normals. Default none
#' @param access [OPTIONAL] Path to reference genome accessibility. Default none
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
#' @param pon [OPTIONAL] Path to panel of normal.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param seg_method [OPTIONAL] Method to use to generate segmentation calls.Default  cbs. Options ["cbs","flasso","haar","none","hmm","hmm-tumor","hmm-germline"]
#' @param seq_method [OPTIONAL] Sequenced methods used. Default hybrid. Options ["hybrid","amplicon","wgs"]
#' @param trend_scatter [OPTIONAL] Show trendline. Default TRUE
#' @param range_scatter [OPTIONAL] Range to show in chr:start-end format. Default none
#' @param range_list_scatter [OPTIONAL] Path to BED file with range list. Default none
#' @param genes_scatter [OPTIONAL] List of genes to focus on. Default none.
#' @param margin_width_scatter [OPTIONAL] Default gene margin width. Default 1000000
#' @param plot_by_bin_scatter [OPTIONAL] Assume all bins are the same size. Default FALSE
#' @param trend_scatter [OPTIONAL] Show trendline. Default TRUE.
#' @param antitarget_symbol_scatter [OPTIONAL] Show antitarget probes with this symbol. Default "@"
#' @param segment_colour_scatter [OPTIONAL] Show antitarget probes with this symbol. Default "@"
#' @param y_max_scatter [OPTIONAL] Y max range. Default NULL.
#' @param y_min_scatter [OPTIONAL] Y min range. Default NULL.
#' @param cn_thr_diagram [OPTIONAL] Copy number threshold. Default 0.5
#' @param min_probes_diagram [OPTIONAL] Minimum number of probes to show a gene. Default 3
#' @param male_reference_diagram [OPTIONAL]  Assume inputs were normalized to a male reference. Default FALSE
#' @param gender_diagram [OPTIONAL] Assume sample gender. Default male
#' @param shift_diagram [OPTIONAL] Adjust the X and Y chromosomes according to sample sex. Default TRUE
#' @param normal_ichor [OPTIONAL] Path to normal samples. Default c(0.5,0.6,0.7,0.8,0.9)
#' @param ploidy_ichor [OPTIONAL] Initial tumour ploidy; can be more than one value if additional ploidy initializations are desired. Default: 2
#' @param lambda_ichor [OPTIONAL] Initial Student's t precision; must contain 4 values (e.g. c(1500,1500,1500,1500)); if not provided then will automatically use based on variance of data
#' @param scStates_ichor [OPTIONAL]  Subclonal states to consider. Default NULL
#' @param output_name_ichor [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
#' @param lambdaScaleHyperParam_ichor [OPTIONAL] Hyperparameter (scale) for Gamma prior on Student's-t precision. Default 3
#' @param maxCN_ichor [OPTIONAL] Total clonal states. Default 7.
#' @param estimateNormal_ichor [OPTIONAL] Estimate normal. Default TRUE.
#' @param estimateScPrevalence_ichor [OPTIONAL] Estimate subclonal prevalence. Default TRUE.
#' @param estimatePloidy_ichor [OPTIONAL] Estimate tumour ploidy. Default TRUE.
#' @param maxFracGenomeSubclone_ichor [OPTIONAL] Exclude solutions with subclonal genome fraction greater than this value. Default 0.5
#' @param maxFracCNASubclone_ichor  [OPTIONAL] Exclude solutions with fraction of subclonal events greater than this value. Default 0.7
#' @param minSegmentBins_ichor [OPTIONAL]  Minimum number of bins for largest segment threshold required to estimate tumor fraction; if below this threshold, then will be assigned zero tumor fraction
#' @param altFracThreshold_ichor [OPTIONAL] Minimum proportion of bins altered required to estimate tumor fraction; if below this threshold, then will be assigned zero tumor fraction. Default: [0.05]
#' @param includeHOMD_ichor [OPTIONAL] If FALSE, then exclude HOMD state. Useful when using large bins (e.g. 1Mb). Default FALSE.
#' @param txnE_ichor [OPTIONAL] Self-transition probability. Increase to decrease number of segments. Default: [0.9999999]
#' @param txnStrength_ichor [OPTIONAL] Transition pseudo-counts. Exponent should be the same as the number of decimal places of --txnE. Default: [1e+07]
#' @param plotFileType_ichor [OPTIONAL] File format for output plots. Default pdf
#' @param plotYLim_ichor [OPTIONAL] ylim to use for chromosome plots. Default: [c(-2,2)]
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


process_cnvkit=function(
    rdata=NULL,selected=FALSE,
    target=build_default_reference_list()$HG19$panel$PCF_V3$binned$target,
    antitarget=build_default_reference_list()$HG19$panel$PCF_V3$binned$antitarget,
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    access=build_default_reference_list()$HG19$reference$access_5k,
    sif_cnvkit=build_default_sif_list()$sif_cnvkit,
    pon=build_default_reference_list()$HG19$panel$PCF_V3$variant$pon_cn_male,
    tumour="",
    seg_method="cbs",
    seq_method="hybrid",
    output_name="",
    output_dir=".",
    diagram=TRUE,
    scatter=TRUE,
    verbose=FALSE,
    read_count=FALSE,
    min_mapq=0,
    gc=TRUE,
    edge=TRUE,
    rmask=TRUE,
    smooth=FALSE,
    drop_low_coverage=TRUE,
    drop_outliers=10,
    range_scatter="",
    range_list_scatter="",
    genes_scatter="",
    margin_width_scatter=1000000,
    plot_by_bin_scatter=FALSE,
    trend_scatter=FALSE,
    antitarget_symbol_scatter="@",
    segment_colour_scatter="red",
    y_max_scatter=NULL,
    y_min_scatter=NULL,
    cn_thr_diagram=0.5,
    min_probes_diagram=3,
    male_reference_diagram=FALSE,
    gender_diagram="male",
    shift_diagram=TRUE,
    batch_config=build_default_preprocess_config(),
    threads=4,ram=4,mode="local",
    executor_id=make_unique_id("processCNVkit"),
    task_name="processCNVkit",time="48:0:0",
    update_time=60,
    wait=FALSE,hold=NULL
  ){

    if(!is.null(rdata)){
        load(rdata)
        if(!is.null(selected)){
            tumour=region_list[selected]
        }
    }
      
  id=""
    if(output_name!=""){
      id=output_name
    }else{
      id=get_file_name(tumour)
    }



    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir,name=id)
    job=build_job(executor_id=executor_id,task=task_id)



    jobs_report=build_job_report(
      job_id=job,
      executor_id=executor_id,
      exec_code=list(), 
      task_id=task_id,
      input_args = argg,
      out_file_dir=out_file_dir,
      out_files=list()
    )


 
  jobs_report[["steps"]][["targetCoverageCNVkit"]]<-coverage_cnvkit(
        sif_cnvkit=sif_cnvkit,
        ref_genome=ref_genome,
        bed=target,
        bam=tumour,
        output_name=paste0(id,".targetcoverage"),
        output_dir=out_file_dir,
        read_count=read_count,
        min_mapq=min_mapq,
        verbose=verbose,
        batch_config=batch_config,
        threads=threads,ram=ram,mode=mode,
        executor_id=task_id,
        time=time,
        hold=hold
  )

  jobs_report[["steps"]][["antitargetCoverageCNVkit"]]<-coverage_cnvkit(
        sif_cnvkit=sif_cnvkit,
        ref_genome=ref_genome,
        bed=antitarget,
        bam=tumour,
        output_name=paste0(id,".antitargetcoverage"),
        output_dir=out_file_dir,
        read_count=read_count,
        min_mapq=min_mapq,
        verbose=verbose,
        batch_config=batch_config,
        threads=threads,ram=ram,mode=mode,
        executor_id=task_id,
        time=time,
        hold=hold
  )

  jobs_report[["steps"]][["fixCNVkit"]]<-fix_cnvkit(
        sif_cnvkit=sif_cnvkit,
        pon=pon,
        target=jobs_report[["steps"]][["targetCoverageCNVkit"]]$out_files$cnn,
        antitarget=jobs_report[["steps"]][["antitargetCoverageCNVkit"]]$out_files$cnn,
        output_name=id,
        output_dir=out_file_dir,
        gc=gc,
        edge=edge,
        rmask=rmask,
        verbose=verbose,
        batch_config=batch_config,
        threads=1,ram=1,mode=mode,
        executor_id=task_id,
        time=time,
        hold=c(jobs_report[["steps"]][["targetCoverageCNVkit"]]$job_id,
        jobs_report[["steps"]][["antitargetCoverageCNVkit"]]$job_id)
  )

  jobs_report[["steps"]][["segmentCNVkit"]]<-segment_cnvkit(
        sif_cnvkit=sif_cnvkit,
        cnr=jobs_report[["steps"]][["fixCNVkit"]]$out_files$cnr,
        seg_method=seg_method,
        output_name=id,
        output_dir=out_file_dir,
        smooth=smooth,
        drop_low_coverage=drop_low_coverage,
        drop_outliers=drop_outliers,
        verbose=verbose,
        batch_config=batch_config,
        threads=threads,ram=ram,mode=mode,
        executor_id=task_id,
        time=time,
        hold=jobs_report[["steps"]][["fixCNVkit"]]$job_id
  )


  if(scatter){
    jobs_report[["steps"]][["scatterCNVkit"]]<-scatter_cnvkit(
      sif_cnvkit=sif_cnvkit,
      cnr=jobs_report[["steps"]][["fixCNVkit"]]$out_files$cnr,
      cns=jobs_report[["steps"]][["segmentCNVkit"]]$out_files$cns,
      output_name=id,
      output_dir=out_file_dir,
      range=range_scatter,
      range_list=range_list_scatter,
      genes=genes_scatter,
      margin_width=margin_width_scatter,
      plot_by_bin=plot_by_bin_scatter,
      trend=trend_scatter,
      antitarget_symbol=antitarget_symbol_scatter,
      segment_colour=segment_colour_scatter,
      title=id,
      y_max=y_max_scatter,
      y_min=y_max_scatter,
      verbose=verbose,
      batch_config=batch_config,
      threads=1,ram=1,mode=mode,
      executor_id=task_id,
      time=time,
      hold=jobs_report[["steps"]][["segmentCNVkit"]]$job_id
    )
  }

  if(diagram){
      jobs_report[["steps"]][["diagramCNVkit"]]<-diagram_cnvkit(
        sif_cnvkit=sif_cnvkit,
        cnr=jobs_report[["steps"]][["fixCNVkit"]]$out_files$cnr,
        cns=jobs_report[["steps"]][["segmentCNVkit"]]$out_files$cns,
        output_name=id,
        output_dir=out_file_dir,
        cn_thr=cn_thr_diagram,
        min_probes=min_probes_diagram,
        male_reference=male_reference_diagram,
        gender=gender_diagram,
        shift=shift_diagram,
        title=id,
        verbose=verbose,
        batch_config=batch_config,
        threads=1,ram=1,mode=mode,
        executor_id=task_id,
        time=time,
        hold=jobs_report[["steps"]][["segmentCNVkit"]]$job_id
      )
  }


  return(jobs_report)

 
}

#' Multiregion parallelization of Haplotypecaller Gatk Variant Calling for multiple-samples
#'
#' This function functions calls Haplotypecaller across multiple regions in parallel across multiple sanmples
#' If a vector of normal samples are provided these will be processed in co-joint calling  mode.
#' To run in normal mode suppply a normal sample.
#' 
#' //TODO validate co-joint calling mode
#' 
#' For more information read:
#' https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller
#'
#' @param rdata [OPTIONAL] Import R data information with list of BAM.
#' @param selected [OPTIONAL] Select BAM from list.
#' @param sif_gatk [REQUIRED] Path to gatk sif file.
#' @param bin_bcftools [REQUIRED] Path to bcftools binary file.
#' @param bin_bgzip [REQUIRED] Path to bgzip binary file.
#' @param bin_tabix [REQUIRED] Path to tabix binary file.
#' @param bin_samtools [REQUIRED] Path to samtools binary file.
#' @param normal [OPTIONAL] Path to normal BAM file.
#' @param patient_id [OPTIONAL] Patient ID.
#' @param ref_genome [REQUIRED] Path to reference genome fasta file.
#' @param indel_db [OPTIONAL] Path to database with indel info.
#' @param haplotype_db [OPTIONAL] Path to database with haplotype info.
#' @param info_key [OPTIONAL] Info Key for CNN model training. Default CNN_1D. Options ["CNN_1D","CNN_2D"]
#' @param snp_tranche [OPTIONAL] SNP tranche cut-off.
#' @param indel_tranche [OPTIONAL] INDEL tranche cut-off.
#' @param keep_previous_tranche [OPTIONAL] Remove previous filter information.
#' @param regions [OPTIONAL] Regions to parallelize through. If not given will be infered from BAM file.
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param threads [OPTIONAL] Number of threads to split the work. Default 4
#' @param filter [OPTIONAL] Filter variants. Default TRUE
#' @param ram [OPTIONAL] RAM memory to asing to each thread. Default 4
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param batch_config [REQUIRED] Additional batch configuration if batch mode selected.
#' @param executor_id Task EXECUTOR ID. Default "recalCovariates"
#' @param task_name Task name. Default "recalCovariates"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export


parallel_samples_process_cnvkit=function(
    target=build_default_reference_list()$HG19$panel$PCF_V3$binned$target,
    antitarget=build_default_reference_list()$HG19$panel$PCF_V3$binned$antitarget,
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    access=build_default_reference_list()$HG19$reference$access_5k,
    sif_cnvkit=build_default_sif_list()$sif_cnvkit,
    pon=build_default_reference_list()$HG19$panel$PCF_V3$variant$pon_cn_male,
    tumour="",
    patient_id="",
    seg_method="cbs",
    seq_method="hybrid",
    output_name="",
    output_dir=".",
    diagram=TRUE,
    scatter=TRUE,
    verbose=FALSE,
    read_count=FALSE,
    min_mapq=0,
    gc=TRUE,
    edge=TRUE,
    rmask=TRUE,
    smooth=FALSE,
    drop_low_coverage=TRUE,
    drop_outliers=10,
    range_scatter="",
    range_list_scatter="",
    genes_scatter="",
    margin_width_scatter=1000000,
    plot_by_bin_scatter=FALSE,
    trend_scatter=FALSE,
    antitarget_symbol_scatter="@",
    segment_colour_scatter="red",
    y_max_scatter=NULL,
    y_min_scatter=NULL,
    cn_thr_diagram=0.5,
    min_probes_diagram=3,
    male_reference_diagram=FALSE,
    gender_diagram="male",
    shift_diagram=TRUE,
    batch_config=build_default_preprocess_config(),
    threads=4,ram=4,mode="local",
    executor_id=make_unique_id("parSampleProcessCNVkit"),
    task_name="parSampleProcessCNVkit",time="48:0:0",
    update_time=60,
    wait=FALSE,hold=NULL
){

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir,name=patient_id)
  tmp_dir=set_dir(dir=out_file_dir,name="tmp")
 

  job=build_job(executor_id=executor_id,task_id=task_id)


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


    tumour_list=tumour
    names(tumour_list)=Vectorize(get_file_name)(tumour)

    if(mode=="local"){
      jobs_report[["steps"]][["par_sample_processs_cnvkit"]]<-
      parallel::mclapply(tumour_list,FUN=function(tumour){
      job_report <- process_cnvkit(
          target=target,
          antitarget=antitarget,
          ref_genome=ref_genome,
          access=access,
          sif_cnvkit=sif_cnvkit,
          pon=pon,
          tumour=tumour,
          seg_method=seg_method,
          seq_method=seq_method,
          output_name=get_file_name(tumour),
          output_dir=out_file_dir,
          diagram=diagram,
          scatter=scatter,
          verbose=verbose,
          read_count=read_count,
          min_mapq=min_mapq,
          gc=gc,
          edge=edge,
          rmask=rmask,
          smooth=smooth,
          drop_low_coverage=drop_low_coverage,
          drop_outliers=drop_outliers,
          range_scatter=range_scatter,
          range_list_scatter=range_list_scatter,
          genes_scatter=genes_scatter,
          margin_width_scatter=margin_width_scatter,
          plot_by_bin_scatter=plot_by_bin_scatter,
          trend_scatter=trend_scatter,
          antitarget_symbol_scatter=antitarget_symbol_scatter,
          segment_colour_scatter=segment_colour_scatter,
          y_max_scatter=y_max_scatter,
          y_min_scatter=y_min_scatter,
          cn_thr_diagram=cn_thr_diagram,
          min_probes_diagram=min_probes_diagram,
          male_reference_diagram=male_reference_diagram,
          gender_diagram=gender_diagram,
          shift_diagram=shift_diagram,
          batch_config=batch_config,
          threads=1,ram=4,
          executor_id=task_id,
          time=time,
          hold=hold)

        },mc.cores=threads)
    }else if(mode=="batch"){

            rdata_file=paste0(tmp_dir,"/",job,".samples.RData")
            output_dir=out_file_dir
            executor_id=task_id
            save(
              tumour_list,
              target,
              antitarget,
              ref_genome,
              access,
              sif_cnvkit,
              pon,
              seg_method,
              seq_method,
              output_dir,
              diagram,
              scatter,
              verbose,
              read_count,
              min_mapq,
              gc,
              edge,
              rmask,
              smooth,
              drop_low_coverage,
              drop_outliers,
              range_scatter,
              range_list_scatter,
              genes_scatter,
              margin_width_scatter,
              plot_by_bin_scatter,
              trend_scatter,
              antitarget_symbol_scatter,
              segment_colour_scatter,
              y_max_scatter,
              y_min_scatter,
              cn_thr_diagram,
              min_probes_diagram,
              male_reference_diagram,
              gender_diagram,
              shift_diagram
            )
            exec_code=paste0("Rscript -e \"ULPwgs::process_cnvkit(rdata=\\\"",
            rdata_file,"\\\",selected=$SGE_TASK_ID)\"")
            out_file_dir2=set_dir(dir=out_file_dir,name="batch")
            batch_code=build_job_exec(job=job,time=time,ram=ram,
            threads=1,output_dir=out_file_dir2,
            hold=hold,array=length(tumour_list))
            exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)

            if(verbose){
                print_verbose(job=job,arg=argg,exec_code=exec_code)
            }
            error=execute_job(exec_code=exec_code)
            if(error!=0){
                stop("gatk failed to run due to unknown error.
                Check std error for more information.")
            }
    
          jobs_report[["steps"]][["par_sample_processs_cnvkit"]]<- build_job_report(
                job_id=job,
                executor_id=executor_id,
                exec_code=exec_code, 
                task_id=task_id,
                input_args=argg,
                out_file_dir=out_file_dir,
                out_files=list(
                  diagram=ifelse(diagram,paste0(out_file_dir,"/",
                  names(tumour_list),"/cnvkit_reports/",
                  names(tumour_list),".diagram.pdf"),""),
                  scatter=ifelse(diagram,paste0(out_file_dir,"/",
                  names(tumour_list),"/cnvkit_reports/",
                  names(tumour_list),".scatter.pdf"),""),
                  target_cnn=paste0(out_file_dir,"/",
                  names(tumour_list),"/cnvkit_reports/",
                  names(tumour_list),".target.cnn"),
                  antitarget_cnn=paste0(out_file_dir,"/",
                  names(tumour_list),"/cnvkit_reports/",
                  names(tumour_list),".antitarget.cnn"),
                  corrected_cnr=paste0(out_file_dir,"/",
                  names(tumour_list),"/cnvkit_reports/",
                  names(tumour_list),".cnr"),
                  segments_cns=paste0(out_file_dir,"/",
                  names(tumour_list),"/cnvkit_reports/",
                  names(tumour_list),".cns"),
                )
            )
    }

  if(wait&&mode=="batch"){
    job_validator(job=unlist_lvl(jobs_report[["steps"]],var="job_id"),time=update_time,
    verbose=verbose,threads=threads)
  }

  return(jobs_report)

}



#' Multiregion parallelization of Haplotypecaller Gatk Variant Calling for multiple-samples for single and multiple samples
#'
#' This function functions calls Haplotypecaller across multiple regions in parallel across multiple sanmples
#' If a vector of normal samples are provided these will be processed in co-joint calling  mode.
#' To run in normal mode suppply a normal sample.
#' 
#' //TODO validate co-joint calling mode
#' 
#' For more information read:
#' https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller
#' 
#' 
#' @param rdata [OPTIONAL] Import R data information with list of BAM.
#' @param selected [OPTIONAL] Select BAM from list.
#' @param sif_gatk [REQUIRED] Path to gatk sif file.
#' @param bin_bcftools [REQUIRED] Path to bcftools binary file.
#' @param bin_bgzip [REQUIRED] Path to bgzip binary file.
#' @param bin_tabix [REQUIRED] Path to tabix binary file.
#' @param bin_samtools [REQUIRED] Path to samtools binary file.
#' @param sample_sheet [OPTIONAL] Path to sheet with sample information.
#' @param normal [OPTIONAL] Path to normal BAM file.
#' @param patient_id [OPTIONAL] Patient ID.
#' @param bam_dir [OPTIONAL] Path to directory with BAM files.
#' @param pattern [OPTIONAL] Pattern to use to search for BAM files in BAM directory.
#' @param ref_genome [REQUIRED] Path to reference genome fasta file.
#' @param indel_db [OPTIONAL] Path to database with indel info.
#' @param haplotype_db [OPTIONAL] Path to database with haplotype info.
#' @param info_key [OPTIONAL] Info Key for CNN model training. Default CNN_1D. Options ["CNN_1D","CNN_2D"]
#' @param snp_tranche [OPTIONAL] SNP tranche cut-off.
#' @param indel_tranche [OPTIONAL] INDEL tranche cut-off.
#' @param keep_previous_tranche [OPTIONAL] Remove previous filter information.
#' @param regions [OPTIONAL] Regions to parallelize through. If not given will be infered from BAM file.
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param threads [OPTIONAL] Number of threads to split the work. Default 4
#' @param filter [OPTIONAL] Filter variants. Default TRUE
#' @param ram [OPTIONAL] RAM memory to asing to each thread. Default 4
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param batch_config [REQUIRED] Additional batch configuration if batch mode selected.
#' @param executor_id Task EXECUTOR ID. Default "recalCovariates"
#' @param task_name Task name. Default "recalCovariates"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export


multisample_process_cnvkit=function(
    target=build_default_reference_list()$HG19$panel$PCF_V3$binned$target,
    antitarget=build_default_reference_list()$HG19$panel$PCF_V3$binned$antitarget,
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    access=build_default_reference_list()$HG19$reference$access_5k,
    sif_cnvkit=build_default_sif_list()$sif_cnvkit,
    pon=build_default_reference_list()$HG19$panel$PCF_V3$variant$pon_cn_male,
    sample_sheet=NULL,
    bam_dir="",
    patient_id="",
    pattern="bam$",
    header=TRUE,
    sep="\t",
    seg_method="cbs",
    seq_method="hybrid",
    output_name="",
    output_dir=".",
    diagram=TRUE,
    scatter=TRUE,
    verbose=FALSE,
    read_count=FALSE,
    min_mapq=0,
    gc=TRUE,
    edge=TRUE,
    rmask=TRUE,
    smooth=FALSE,
    drop_low_coverage=TRUE,
    drop_outliers=10,
    range_scatter="",
    range_list_scatter="",
    genes_scatter="",
    margin_width_scatter=1000000,
    plot_by_bin_scatter=FALSE,
    trend_scatter=FALSE,
    antitarget_symbol_scatter="@",
    segment_colour_scatter="red",
    y_max_scatter=NULL,
    y_min_scatter=NULL,
    cn_thr_diagram=0.5,
    min_probes_diagram=3,
    male_reference_diagram=FALSE,
    gender_diagram="male",
    shift_diagram=TRUE,
    batch_config=build_default_preprocess_config(),
    threads=4,ram=4,mode="local",
    executor_id=make_unique_id("parSampleProcessCNVkit"),
    task_name="parSampleProcessCNVkit",time="48:0:0",
    update_time=60,
    wait=FALSE,hold=NULL
){

  argg <- as.list(environment())

        task_id=make_unique_id(task_name)
        out_file_dir=set_dir(dir=output_dir)

        job=build_job(executor_id=executor_id,task_id=task_id)

        job_report=build_job_report(
                job_id=job,
                executor_id=executor_id,
                exec_code=list(), 
                task_id=task_id,
                input_args=argg,
                out_file_dir=out_file_dir,
                out_files=list(
                )
        )



    columns=c(
        "seg_method",
        "seq_method",
        "output_dir",
        "diagram",
        "scatter",
        "verbose",
        "read_count",
        "min_mapq",
        "gc",
        "edge",
        "rmask",
        "smooth",
        "drop_low_coverage",
        "drop_outliers",
        "range_scatter",
        "range_list_scatter",
        "genes_scatter",
        "margin_width_scatter",
        "plot_by_bin_scatter",
        "trend_scatter",
        "antitarget_symbol_scatter",
        "segment_colour_scatter",
        "y_max_scatter",
        "y_min_scatter",
        "cn_thr_diagram",
        "min_probes_diagram",
        "male_reference_diagram",
        "gender_diagram",
        "shift_diagram",
        "batch_config",
        "threads",
        "ram",
        "mode",
        "time",
        "hold"
    )


    if(!is.null(sample_sheet)){
      
        if(!is.data.frame(sample_sheet)){
                file_info=read.csv(sample_sheet,header=header,sep=sep,stringsAsFactors=FALSE)
                if(!header){
                    names(file_info)=columns
                }
        }else{
                file_info=sample_sheet
        }
        
        file_info=file_info %>% dplyr::group_by(dplyr::across(-tumour)) %>% dplyr::summarise(normal=list(tumour))

        job_report[["steps"]][["multisample_process_cnvkit"]]=parallel::mclapply(seq(1,nrow(file_info)),FUN=function(x){
            
            lapply(columns,FUN=function(col){
                if(is.null(file_info[[col]])){
                    file_info[[col]]<<-get(col)
                }

                if(is.null(file_info[[x,col]])){
                    file_info[[x,col]]<<-get(col)
                }
            
            })
        
            job_report<-parallel_samples_process_cnvkit(     
                target=file_info[x,]$target,
                antitarget=file_info[x,]$antitarget,
                ref_genome=file_info[x,]$ref_genome,
                access=file_info[x,]$access,
                sif_cnvkit=file_info[x,]$sif_cnvkit,
                pon=file_info[x,]$pon,
                tumour=file_info[x,]$tumour,
                seg_method=file_info[x,]$seg_method,
                seq_method=file_info[x,]$seq_method,
                output_dir=file_info[x,]$output_dir,
                diagram=file_info[x,]$diagram,
                scatter=file_info[x,]$scatter,
                verbose=file_info[x,]$verbose,
                read_count=file_info[x,]$read_count,
                min_mapq=file_info[x,]$min_mapq,
                gc=file_info[x,]$gc,
                edge=file_info[x,]$edge,
                rmask=file_info[x,]$rmask,
                smooth=file_info[x,]$smooth,
                drop_low_coverage=file_info[x,]$drop_low_coverage,
                drop_outliers=file_info[x,]$drop_outliers,
                range_scatter=file_info[x,]$range_scatter,
                range_list_scatter=file_info[x,]$range_list_scatter,
                genes_scatter=file_info[x,]$genes_scatter,
                margin_width_scatter=file_info[x,]$margin_width_scatter,
                plot_by_bin_scatter=file_info[x,]$plot_by_bin_scatter,
                trend_scatter=file_info[x,]$trend_scatter,
                antitarget_symbol_scatter=file_info[x,]$antitarget_symbol_scatter,
                segment_colour_scatter=file_info[x,]$segment_colour_scatter,
                y_max_scatter=file_info[x,]$x_min_scatter,
                y_min_scatter=file_info[x,]$y_min_scatter,
                cn_thr_diagram=file_info[x,]$cn_thr_diagram,
                min_probes_diagram=file_info[x,]$min_probes_diagram,
                male_reference_diagram=file_info[x,]$male_reference_diagram,
                gender_diagram=file_info[x,]$gender_diagram,
                shift_diagram=file_info[x,]$shift_diagram,
                batch_config=file_info[x,]$batch_config,
                threads=file_info[x,]$threads,ram=file_info[x,]$ram,
                mode=file_info[x,]$mode,
                executor_id=task_id,
                time=file_info[x,]$time,
                hold=file_info[x,]$hold
              )
            },mc.cores=ifelse(mode=="local",1,3))

    }else{
        bam_dir_path=system(paste("realpath",bam_dir),intern=TRUE)
        tumour=system(paste0("find ",bam_dir_path,"| grep ",pattern),intern=TRUE)


            job_report<-parallel_samples_process_cnvkit(     
                target=target,
                antitarget=antitarget,
                ref_genome=ref_genome,
                access=access,
                sif_cnvkit=sif_cnvkit,
                pon=pon,
                tumour=tumour,
                seg_method=seg_method,
                seq_method=seq_method,
                output_dir=output_dir,
                diagram=diagram,
                scatter=scatter,
                verbose=verbose,
                read_count=read_count,
                min_mapq=min_mapq,
                gc=gc,
                edge=edge,
                rmask=rmask,
                smooth=mooth,
                drop_low_coverage=drop_low_coverage,
                drop_outliers=drop_outliers,
                range_scatter=range_scatter,
                range_list_scatter=range_list_scatter,
                genes_scatter=genes_scatter,
                margin_width_scatter=margin_width_scatter,
                plot_by_bin_scatter=plot_by_bin_scatter,
                trend_scatter=trend_scatter,
                antitarget_symbol_scatter=antitarget_symbol_scatter,
                segment_colour_scatter=segment_colour_scatter,
                y_max_scatter=x_min_scatter,
                y_min_scatter=y_min_scatter,
                cn_thr_diagram=cn_thr_diagram,
                min_probes_diagram=min_probes_diagram,
                male_reference_diagram=male_reference_diagram,
                gender_diagram=gender_diagram,
                shift_diagram=shift_diagram,
                batch_config=file_info[x,]$batch_config,
                threads=threads,ram=ram,
                mode=mode,
                executor_id=task_id,
                time=time,
                hold=hold
              )
    }


    if(wait&&mode=="batch"){
        job_validator(job=unlist_level(named_list=job_report[["steps"]][["multisample_haplotypecaller"]],var="job_id"),
        time=update_time,verbose=verbose,threads=threads)
    }

    return(job_report)

}











#' Wrapper around access command in CNVkit
#'
#' Create a BED file with non-mappable regions to exclude from the genome.
#' Additional regions can be supplied in a BED format using the exclude_regions argument
#' 
#' For more information read:
#' https://cnvkit.readthedocs.io/en/stable/pipeline.html
#'
#' @param sif_cnvkit [REQUIRED] Path to cnvkit sif file.
#' @param ref_genome [REQUIRED] Path to reference genome fasta file.
#' @param exclude_regions [REQUIRED] Additional regions to exclude in BED format. Default none.
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
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


access_cnvkit=function(
    sif_cnvkit=build_default_sif_list()$sif_cnvkit,
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    output_name="access",
    output_dir=".",
    exclude_regions="",
    gap_size=5000,
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=1,mode="local",
    executor_id=make_unique_id("accessCNVkit"),
    task_name="accessCNVkit",time="48:0:0",
    update_time=60,
    wait=FALSE,hold=NULL
  ){


    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir)
    job=build_job(executor_id=executor_id,task_id=task_id)


    out_file=paste0(out_file_dir,"/",output_name,".bed")

    if(exclude_regions!=""){
      exclude_regions=paste0(" -x ",exclude_regions)
    }
      
    exec_code=paste("singularity exec -H ",paste0(getwd(),":/home "),sif_cnvkit,
    " cnvkit.py access -o ",out_file,exclude_regions," -s ",gap_size, ref_genome)


    if(mode=="batch"){
        out_file_dir2=set_dir(dir=out_file_dir,name="batch")
        batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
        threads=threads,output_dir=out_file_dir2)
        exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
    }

    if(verbose){
        print_verbose(job=job,arg=argg,exec_code=exec_code)
    }

    error=execute_job(exec_code=exec_code)
    
    if(error!=0){
      stop("cnvkit failed to run due to unknown error.
      Check std error for more information.")
    }

    job_report=build_job_report(
      job_id=job,
      executor_id=executor_id,
      exec_code=exec_code, 
      task_id=task_id,
      input_args = argg,
      out_file_dir=out_file_dir,
      out_files=list(
        access=out_file)
    )


    if(wait&&mode=="batch"){
      job_validator(job=job_report$job_id,time=update_time,
      verbose=verbose,threads=threads)
    }

    return(job_report)

  }







#' Wrapper around target function from CNVkit
#'
#' This function wraps around target function for CNVkit
#' This function generates an target BED file from an input target file. 
#' Additional parameters can be used to exclude regions and modify the average bin size
#' 
#' 
#' For more information read:
#' https://cnvkit.readthedocs.io/en/stable/pipeline.html
#' 
#' 
#' @param sif_cnvkit [REQUIRED] Path to cnvkit sif file.
#' @param ref_genome [REQUIRED] Path to reference genome.
#' @param bin_size_target [OPTIONAL] Size of bins for targets. Default 75
#' @param bin_size_antitarget [OPTIONAL] Size of bins for antitargets. Default 500000
#' @param min_bin_size_antitarget [OPTIONAL] Size of bins for antitargets. Default NULL
#' @param access [OPTIONAL] Path to reference genome accessibility. Default none
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
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


  create_pon_cnvkit=function(
    sif_cnvkit=build_default_sif_list()$sif_cnvkit,
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    access=build_default_reference_list()$HG19$reference$access_5k,
    output_name="reference",
    normals="",
    output_dir=".",
    target="",
    bin_size_target=75,
    bin_size_antitarget=500000,
    min_bin_size_antitarget=NULL,
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=1,mode="local",
    executor_id=make_unique_id("createPonCNVkit"),
    task_name="createPonCNVkit",time="48:0:0",
    update_time=60,
    wait=FALSE,hold=NULL
  ){

    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir)
    job=build_job(executor_id=executor_id,task_id=task_id)


    out_file=paste0(out_file_dir,"/",output_name,".cnn")


    if(target!=""){
      target=paste0(" -t ",target)
    }


    if(access!=""){
      paste0(" -g ",access)
    }
    if(!is.null(min_bin_size_antitarget)){
      min_bin_size_antitarget=paste0(" --antitarget-min-size ",min_bin_size_antitarget)
    }

    exec_code=paste("singularity exec -H ",paste0(getwd(),":/home "),sif_cnvkit,
    " cnvkit.py batch -p ",threads, " -n ",paste0(normals,collapse=" "),
     " --output-reference ",out_file," -f ", ref_genome,access,target,
     " --target-avg-size ",bin_size_target,
     " --antitarget-avg-size ",bin_size_antitarget,min_bin_size_antitarget
     )


    if(mode=="batch"){
        out_file_dir2=set_dir(dir=out_file_dir,name="batch")
        batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
        threads=threads,output_dir=out_file_dir2)
        exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
    }

    if(verbose){
        print_verbose(job=job,arg=argg,exec_code=exec_code)
    }

    error=execute_job(exec_code=exec_code)
    
    if(error!=0){
      stop("cnvkit failed to run due to unknown error.
      Check std error for more information.")
    }

    job_report=build_job_report(
      job_id=job,
      executor_id=executor_id,
      exec_code=exec_code, 
      task_id=task_id,
      input_args = argg,
      out_file_dir=out_file_dir,
      out_files=list(
        access=out_file)
    )


    if(wait&&mode=="batch"){
      job_validator(job=job_report$job_id,time=update_time,
      verbose=verbose,threads=threads)
    }

    return(job_report)

  }







#' Wrapper around target function from CNVkit
#'
#' This function wraps around target function for CNVkit
#' This function generates an target BED file from an input target file. 
#' Additional parameters can be used to exclude regions and modify the average bin size
#' 
#' 
#' For more information read:
#' https://cnvkit.readthedocs.io/en/stable/pipeline.html
#'
#' @param sif_cnvkit [REQUIRED] Path to cnvkit sif file.
#' @param bed [REQUIRED] Path to input BED file with target regions. Default none
#' @param annotation [OPTIONAL] BED file with region annotation information. Default none
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param short_names [OPTIONAL] Use short annotation names. Default FALSE
#' @param bin_size [OPTIONAL] Average bin size. Default 100.
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


  create_target_cnvkit=function(
    sif_cnvkit=build_default_sif_list()$sif_cnvkit,
    bed="",
    annotation="",
    output_name="",
    output_dir=".",
    split=FALSE,
    short_names=FALSE,
    bin_size=100,
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=1,mode="local",
    executor_id=make_unique_id("targetCNVkit"),
    task_name="targetCNVkit",time="48:0:0",
    update_time=60,
    wait=FALSE,hold=NULL
  ){

    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir)
    job=build_job(executor_id=executor_id,task_id=task_id)

      
    id=""
    if(output_name!=""){
      id=output_name
    }else{
      id=get_file_name(bed)
    }

    if(annotation!=""){
      annotation=paste0(" --annotate ",annotation)
    }

    add=""
    if(short_names){
      add=paste(add," --short-names ")
    }

    out_file=paste0(out_file_dir,"/",id,".binned.targets.bed")

    exec_code=paste("singularity exec -H ",paste0(getwd(),":/home "),sif_cnvkit,
    " cnvkit.py target --split -a ",bin_size," -o ",out_file, add, bed)


    if(mode=="batch"){
        out_file_dir2=set_dir(dir=out_file_dir,name="batch")
        batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
        threads=threads,output_dir=out_file_dir2)
        exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
    }

    if(verbose){
        print_verbose(job=job,arg=argg,exec_code=exec_code)
    }

    error=execute_job(exec_code=exec_code)
    
    if(error!=0){
      stop("cnvkit failed to run due to unknown error.
      Check std error for more information.")
    }

    job_report=build_job_report(
      job_id=job,
      executor_id=executor_id,
      exec_code=exec_code, 
      task_id=task_id,
      input_args = argg,
      out_file_dir=out_file_dir,
      out_files=list(
        target=out_file)
    )


    if(wait&&mode=="batch"){
      job_validator(job=job_report$job_id,time=update_time,
      verbose=verbose,threads=threads)
    }

    return(job_report)

  }



#' Wrapper around antitarget function from CNVkit
#'
#' This function wraps around antitarget function for CNVkit
#' This function generates an antitarget BED file from an input target file. 
#' Additional parameters can be used to exclude regions and modify the average bin size
#' 
#' 
#' For more information read:
#' https://cnvkit.readthedocs.io/en/stable/pipeline.html
#'
#' @param sif_cnvkit [REQUIRED] Path to cnvkit sif file.
#' @param bed [REQUIRED] Path to input BED file with target regions. Default none
#' @param access [OPTIONAL] Path to non-accessible regions to exclude. Default none
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param bin_size [OPTIONAL] Average bin size. Default 100.
#' @param min_bin_size [OPTIONAL] Minimum average bin size. Bins smaller that this will be dropped. Default NULL.
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

  create_antitarget_cnvkit=function(
    sif_cnvkit=build_default_sif_list()$sif_cnvkit,
    access=build_default_reference_list()$HG19$reference$access_5k,
    bed="",
    output_name="",
    output_dir=".",
    bin_size=500000,
    min_bin_size=NULL,
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=1,mode="local",
    executor_id=make_unique_id("antitargetCNVkit"),
    task_name="antitargetCNVkit",time="48:0:0",
    update_time=60,
    wait=FALSE,hold=NULL
  ){

    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir)
    job=build_job(executor_id=executor_id,task_id=task_id)


    id=""
    if(output_name!=""){
      id=output_name
    }else{
      id=get_file_name(bed)
    }

    add=""
    if(!is.null(min_bin_size)){
      add=paste0(" -m ",min_bin_size)
    }

    out_file=paste0(out_file_dir,"/",id,".binned.antitargets.bed")


    if(access!=""){
      access=paste0(" -g ",access)
    }


    exec_code=paste("singularity exec -H ",paste0(getwd(),":/home "),sif_cnvkit,
    " cnvkit.py antitarget -a ",bin_size,access," -o ",out_file, add, bed)

    if(mode=="batch"){
        out_file_dir2=set_dir(dir=out_file_dir,name="batch")
        batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
        threads=threads,output_dir=out_file_dir2)
        exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
    }

    if(verbose){
        print_verbose(job=job,arg=argg,exec_code=exec_code)
    }

    error=execute_job(exec_code=exec_code)
    
    if(error!=0){
      stop("cnvkit failed to run due to unknown error.
      Check std error for more information.")
    }

    job_report=build_job_report(
      job_id=job,
      executor_id=executor_id,
      exec_code=exec_code, 
      task_id=task_id,
      input_args = argg,
      out_file_dir=out_file_dir,
      out_files=list(
        antitarget=out_file)
    )


    if(wait&&mode=="batch"){
      job_validator(job=job_report$job_id,time=update_time,
      verbose=verbose,threads=threads)
    }

    return(job_report)

  }


#' Create binned target and anitarget beds from target BED file
#'
#' This function wraps around target and antitarget functiosn for CNVkit
#' This function generates an target BED file from an input target file. 
#' Additional parameters can be used to exclude regions and modify the average bin size
#' 
#' 
#' For more information read:
#' https://cnvkit.readthedocs.io/en/stable/pipeline.html
#'
#' @param sif_cnvkit [REQUIRED] Path to cnvkit sif file.
#' @param access [OPTIONAL] Path to non-accessible regions to exclude. Default none
#' @param bed [REQUIRED] Path to input BED file with target regions. Default none
#' @param annotation [OPTIONAL] BED file with region annotation information. Default none
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param short_names [OPTIONAL] Use short annotation names. Default FALSE
#' @param bin_size_target [OPTIONAL] Average bin size for target. Default 120.
#' @param bin_size_antitarget [OPTIONAL] Average bin size for target. Default 500000.
#' @param min_bin_size_antitarget [OPTIONAL] Average bin size for target. Default null.
#' @param split [OPTIONAL] Average bin size for target. Default null.
#' @param output_dir [OPTIONAL] Split annotation names.
#' @param output_dir [OPTIONAL] Split annotation names.
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



bin_targets_cnvkit=function(
    sif_cnvkit=build_default_sif_list()$sif_cnvkit,
    access=build_default_reference_list()$HG19$reference$access_5k,
    bed="",
    output_name="",
    annotation="",
    output_dir=".",
    bin_size_antitarget=100000,
    bin_size_target=100,
    min_bin_size_antitarget=NULL,
    split=FALSE,
    short_names=FALSE,
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=1,mode="local",
    executor_id=make_unique_id("bintargetCNVkit"),
    task_name="bintargetCNVkit",time="48:0:0",
    update_time=60,
    wait=FALSE,hold=NULL
){

    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir)
    job=build_job(executor_id=executor_id,task_id=task_id)

    
    jobs_report=build_job_report(
      job_id=job,
      executor_id=executor_id,
      exec_code=list(), 
      task_id=task_id,
      input_args = argg,
      out_file_dir=out_file_dir,
      out_files=list()
    )


    jobs_report[["steps"]][["TargetsCNVkit"]]<-create_target_cnvkit(
      sif_cnvkit=sif_cnvkit,
      bed=bed,
      annotation=annotation,
      output_name=output_name,
      output_dir=out_file_dir,
      split=split,
      short_names=short_names,
      bin_size=bin_size_target,
      verbose=verbose,
      batch_config=batch_config,
      threads=1,ram=1,mode=mode,
      executor_id=task_id,
      time=time,
      hold=hold
  )


   jobs_report[["steps"]][["AntitargetsCNVkit"]]<-create_antitarget_cnvkit(
      sif_cnvkit=sif_cnvkit,
      access=access,
      bed=jobs_report[["steps"]][["TargetsCNVkit"]]$out_files$target,
      output_name=output_name,
      output_dir=out_file_dir,
      bin_size=bin_size_antitarget,
      min_bin_size=min_bin_size_antitarget,
      verbose=verbose,
      batch_config=batch_config,
      threads=1,ram=1,mode=mode,
      executor_id=task_id,
      time=time,
      hold=jobs_report[["steps"]][["TargetsCNVkit"]]$job_id
  )

  return(jobs_report)

}


#' Wrapper around autobin function from CNVkit
#'
#' This function wraps around target function for CNVkit
#' This function generates an target BED file from an input target file. 
#' Additional parameters can be used to exclude regions and modify the average bin size
#' 
#' 
#' For more information read:
#' https://cnvkit.readthedocs.io/en/stable/pipeline.html
#'
#' @param sif_cnvkit [REQUIRED] Path to cnvkit sif file.
#' @param ref_genoma [REQUIRED] Path to reference genoms.
#' @param bed [REQUIRED] Path to input BED file with target regions. Default none
#' @param bams [REQUIRED] Path to BAM files. Default none
#' @param annotation [OPTIONAL] BED file with region annotation information. Default none
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param short_names [OPTIONAL] Use short annotation names. Default FALSE
#' @param seq_method [OPTIONAL] Sequenced methods used. Default hybrid. Options ["hybrid","amplicon","wgs"]
#' @param seq_type [OPTIONAL] Use short annotation names. Default FALSE
#' @param min_bin_size_target [OPTIONAL] Mininimum target bin size. Default 20.
#' @param max_bin_size_target [OPTIONAL] Mininimum target bin size. Default 20000.
#' @param max_bin_size_antitarget [OPTIONAL] Mininimum target bin size. Default 500000.
#' @param min_bin_size_antitarget [OPTIONAL] Mininimum target bin size. Default 500.
#' @param bp_per_bin [OPTIONAL] Bases per bin. Default 100000.
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


  autobin_cnvkit=function(
    sif_cnvkit=build_default_sif_list()$sif_cnvkit,
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    access=build_default_reference_list()$HG19$reference$access_5k,
    bed="",
    bams="",
    annotation="",
    output_name="",
    output_dir=".",
    seq_method="hybrid",
    split=FALSE,
    short_names=FALSE,
    bp_per_bin=100000,
    min_bin_size_target=20,
    max_bin_size_target=20000,
    min_bin_size_antitarget=500,
    max_bin_size_antitarget=500000,
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=1,mode="local",
    executor_id=make_unique_id("autobinCNVkit"),
    task_name="autobinCNVkit",time="48:0:0",
    update_time=60,
    wait=FALSE,hold=NULL
  ){

    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir)
    job=build_job(executor_id=executor_id,task_id=task_id)


    if(annotation!=""){
      annotation=paste0(" --annotate ",annotation)
    }

    if(access!=""){
      access=paste0(" -g ",access)
    }

    add=""
    if(short_names){
      add=paste(add," --short-names ")
    }

    exec_code=paste("singularity exec -H ",paste0(getwd(),":/home "),sif_cnvkit,
    " cnvkit.py autobin -t ",bed," -f ",ref_genome,
    " -b ",bp_per_bin, " -m ",seq_method,access,
    "  --target-max-size ",max_bin_size_target,
    "  --target-min-size ",min_bin_size_target,
    "  --antitarget-min-size ",min_bin_size_antitarget,
    "  --antitarget-max-size ",max_bin_size_antitarget,add,annotation,
     paste0(bams,collapse=" "))
de
    if(mode=="batch"){
        out_file_dir2=set_dir(dir=out_file_dir,name="batch")
        batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
        threads=threads,output_dir=out_file_dir2)
        exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
    }

    if(verbose){
        print_verbose(job=job,arg=argg,exec_code=exec_code)
    }

    error=execute_job(exec_code=exec_code)
    
    if(error!=0){
      stop("cnvkit failed to run due to unknown error.
      Check std error for more information.")
    }

    job_report=build_job_report(
      job_id=job,
      executor_id=executor_id,
      exec_code=exec_code, 
      task_id=task_id,
      input_args = argg,
      out_file_dir=out_file_dir,
      out_files=list(
        target=out_file)
    )


    if(wait&&mode=="batch"){
      job_validator(job=job_report$job_id,time=update_time,
      verbose=verbose,threads=threads)
    }

    return(job_report)

  }



#' Wrapper around autobin function from CNVkit
#'
#' This function wraps around target function for CNVkit
#' This function generates an target BED file from an input target file. 
#' Additional parameters can be used to exclude regions and modify the average bin size
#' 
#' 
#' For more information read:
#' https://cnvkit.readthedocs.io/en/stable/pipeline.html
#'
#' @param sif_cnvkit [REQUIRED] Path to cnvkit sif file.
#' @param ref_genoma [REQUIRED] Path to reference genoms.
#' @param bed [REQUIRED] Path to input BED file with target regions. Default none
#' @param bam [REQUIRED] Path to BAM files. Default none
#' @param read_count [OPTIONAL] Alternative method for coverage. Default FALSE
#' @param min_mapq [OPTIONAL] Minimum mapping quality to count a read for coverage. Default 0.
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
#' @param output_dir [OPTIONAL] Path to the output directory.
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


  coverage_cnvkit=function(
    rdata=NULL,
    selected=NULL,
    sif_cnvkit=build_default_sif_list()$sif_cnvkit,
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    bed="",
    bam="",
    output_name="",
    output_dir=".",
    read_count=FALSE,
    min_mapq=0,
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=1,mode="local",
    executor_id=make_unique_id("coverageCNVkit"),
    task_name="coverageCNVkit",time="48:0:0",
    update_time=60,
    wait=FALSE,hold=NULL
  ){

    if(!is.null(rdata)){
      load(rdata)
      if(!is.null(selected)){
        bam=bam_list[selected]
      }
    }

    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir)
    job=build_job(executor_id=executor_id,task_id=task_id)

         
    id=""
    if(output_name!=""){
      id=output_name
    }else{
      id=get_file_name(bam)
    }
  
    add=""
    if(read_count){
      add=paste(add," -c ")
    }

    out_file=paste0(out_file_dir,"/",id,".cnn")

    exec_code=paste("singularity exec -H ",paste0(getwd(),":/home "),sif_cnvkit,
    " cnvkit.py coverage -p ",threads, "-q ",min_mapq,
    " -f ",ref_genome," -o ",out_file, add, bam, bed)

    if(mode=="batch"){
        out_file_dir2=set_dir(dir=out_file_dir,name="batch")
        batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
        threads=threads,output_dir=out_file_dir2)
        exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
    }

    if(verbose){
        print_verbose(job=job,arg=argg,exec_code=exec_code)
    }

    error=execute_job(exec_code=exec_code)
    
    if(error!=0){
      stop("cnvkit failed to run due to unknown error.
      Check std error for more information.")
    }

    job_report=build_job_report(
      job_id=job,
      executor_id=executor_id,
      exec_code=exec_code, 
      task_id=task_id,
      input_args = argg,
      out_file_dir=out_file_dir,
      out_files=list(
        cnn=out_file)
    )


    if(wait&&mode=="batch"){
      job_validator(job=job_report$job_id,time=update_time,
      verbose=verbose,threads=threads)
    }

    return(job_report)

  }


#' Wrapper around autobin function from CNVkit
#'
#' This function wraps around target function for CNVkit
#' This function generates an target BED file from an input target file. 
#' Additional parameters can be used to exclude regions and modify the average bin size
#' 
#' 
#' For more information read:
#' https://cnvkit.readthedocs.io/en/stable/pipeline.html
#'
#' @param sif_cnvkit [REQUIRED] Path to cnvkit sif file.
#' @param ref_genoma [REQUIRED] Path to reference genoms.
#' @param bed [REQUIRED] Path to input BED file with target regions. Default none
#' @param bam [REQUIRED] Path to BAM files. Default none
#' @param read_count [OPTIONAL] Alternative method for coverage. Default FALSE
#' @param min_mapq [OPTIONAL] Minimum mapping quality to count a read for coverage. Default 0.
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
#' @param output_dir [OPTIONAL] Path to the output directory.
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


  parallel_sample_coverage_cnvkit=function(
    sif_cnvkit=build_default_sif_list()$sif_cnvkit,
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    bed="",
    bams="",
    output_dir=".",
    read_count=FALSE,
    min_mapq=0,
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=1,mode="local",
    executor_id=make_unique_id("parSamplecoverageCNVkit"),
    task_name="parSamplecoverageCNVkit",time="48:0:0",
    update_time=60,
    wait=FALSE,hold=NULL
  ){

    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir,name="cnn")
    tmp_dir=set_dir(dir=output_dir,name="tmp")
    job=build_job(executor_id=executor_id,task_id=task_id)

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

    bam_list=bams
    names(bam_list)=Vectorize(get_file_name)(bams)
    

    if(mode=="local"){
    jobs_report[["steps"]][["par_sample_coverage_cnvkit"]]<-
    parallel::mclapply(bam_list,FUN=function(bam){
      job_report <- coverage_cnvkit(
                sif_cnvkit=sif_cnvkit,
                ref_genome=ref_genome,
                bed=bed,
                bam=bam,
                output_name=get_file_name(bam),
                output_dir=out_file_dir,
                read_count=read_count,
                min_mapq=min_mapq,
                verbose=verbose,
                batch_config=batch_config,
                executor_id=task_id
              )
    },mc.cores=threads)
    
    }else if(mode=="batch"){

          rdata_file=paste0(tmp_dir,"/",job,".samples.RData")
          output_dir=out_file_dir
          save(bam_list,bed,sif_cnvkit,ref_genome,
          output_dir,read_count,min_mapq,output_dir,verbose,file = rdata_file)
          exec_code=paste0("Rscript -e \"ULPwgs::coverage_cnvkit(rdata=\\\"",
          rdata_file,"\\\",selected=$SGE_TASK_ID)\"")
          out_file_dir2=set_dir(dir=out_file_dir,name="batch")
          batch_code=build_job_exec(job=job,time=time,ram=ram,
          threads=1,output_dir=out_file_dir2,
          hold=hold,array=length(bam_list))
          exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)

          if(verbose){
              print_verbose(job=job,arg=argg,exec_code=exec_code)
          }
          error=execute_job(exec_code=exec_code)
          if(error!=0){
              stop("cnvkit failed to run due to unknown error.
              Check std error for more information.")
          }
         
         jobs_report[["steps"]][["par_sample_coverage_cnvkit"]]<- build_job_report(
              job_id=job,
              executor_id=executor_id,
              exec_code=exec_code, 
              task_id=task_id,
              input_args=argg,
              out_file_dir=out_file_dir,
              out_files=list(
                  cnn=paste0(out_file_dir,"/",names(bam_list),".cnn")
              )
        )
    }

      return(jobs_report)


  }














#' Wrapper around reference function from CNVkit
#'
#' This function wraps around target function for CNVkit
#' This function generates an target BED file from an input target file. 
#' Additional parameters can be used to exclude regions and modify the average bin size
#' 
#' 
#' For more information read:
#' https://cnvkit.readthedocs.io/en/stable/pipeline.html
#'
#' @param sif_cnvkit [REQUIRED] Path to cnvkit sif file.
#' @param ref_genoma [REQUIRED] Path to reference genoms.
#' @param cnn [REQUIRED] Path to normal coverage profiles. Default none
#' @param gender [OPTIONAL] Sample gender. Default male.
#' @param target [OPTIONAL] Path to BED file with target regions. Default none.
#' @param antitarget [OPTIONAL] Path to BED file with target regions. Default none.
#' @param cluster [OPTIONAL] Calculate and store summary stats for clustered subsets of the normal samples with similar coverage profiles. Default FALSE
#' @param min_cluster_size [OPTIONAL] Minimum size to keep in reference profiles.
#' @param gc [OPTIONAL] Disble GC bias correction. Default FALSE.
#' @param edge [OPTIONAL] Disble edge correction. Default FALSE.
#' @param rmask [OPTIONAL] Disble repeat mask correction. Default FALSE.
#' @param male_reference [OPTIONAL] Adjust X chromosome to log2 0. Default FALSE.
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
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


  reference_cnvkit=function(
    sif_cnvkit=build_default_sif_list()$sif_cnvkit,
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    cnn="",
    output_name="reference",
    output_dir=".",
    gender="male",
    target="",
    antitarget="",
    gc=TRUE,
    edge=TRUE,
    rmask=TRUE,
    verbose=FALSE,
    cluster=FALSE,
    min_cluster_size=10,
    male_reference=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=1,mode="local",
    executor_id=make_unique_id("referenceCNVkit"),
    task_name="referenceCNVkit",time="48:0:0",
    update_time=60,
    wait=FALSE,hold=NULL
  ){

    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir)
    job=build_job(executor_id=executor_id,task_id=task_id)

    if(!is.null(gender)){
      gender=paste0(" -x ",gender)
    }

    add=""

    if(!gc){
      add=paste(add," --no-gc ")
    }

    if(!edge){
      add=paste(add," --no-edge ")
    }

    if(!rmask){
      add=paste(add," --no-rmask ")
    }


    if(male_reference){
      add=paste(add," -y ")
    }

    if(cluster){
      add=paste(add," -c ")
      min_cluster_size=paste0(" --min-cluster-size ",min_cluster_size) 
    }else{
      min_cluster_size=""
    }
  


    out_file=paste0(out_file_dir,"/",output_name,".cnn")

    exec_code=paste("singularity exec -H ",paste0(getwd(),":/home "),sif_cnvkit,
    " cnvkit.py reference ", gender, min_cluster_size,
    " -f ",ref_genome," -o ",out_file, add, paste0(cnn,collapse=" "))

    if(mode=="batch"){
        out_file_dir2=set_dir(dir=out_file_dir,name="batch")
        batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
        threads=threads,output_dir=out_file_dir2)
        exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
    }

    if(verbose){
        print_verbose(job=job,arg=argg,exec_code=exec_code)
    }

    error=execute_job(exec_code=exec_code)
    
    if(error!=0){
      stop("cnvkit failed to run due to unknown error.
      Check std error for more information.")
    }

    job_report=build_job_report(
      job_id=job,
      executor_id=executor_id,
      exec_code=exec_code, 
      task_id=task_id,
      input_args = argg,
      out_file_dir=out_file_dir,
      out_files=list(
        cnn=out_file)
    )


    if(wait&&mode=="batch"){
      job_validator(job=job_report$job_id,time=update_time,
      verbose=verbose,threads=threads)
    }

    return(job_report)

  }





#' Wrapper around reference function from CNVkit
#'
#' This function wraps around target function for CNVkit
#' This function generates an target BED file from an input target file. 
#' Additional parameters can be used to exclude regions and modify the average bin size
#' 
#' 
#' For more information read:
#' https://cnvkit.readthedocs.io/en/stable/pipeline.html
#'
#' @param sif_cnvkit [REQUIRED] Path to cnvkit sif file.oms.
#' @param target [REQUIRED] Path to CNN file for target regions. Default none.
#' @param antitarget [REQUIRED] Path to BED file with target regions. Default none.
#' @param gc [OPTIONAL] Disble GC bias correction. Default FALSE.
#' @param edge [OPTIONAL] Disble edge correction. Default FALSE.
#' @param rmask [OPTIONAL] Disble repeat mask correction. Default FALSE.
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
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


 fix_cnvkit=function(
    rdata=NULL,
    selected=NULL,
    sif_cnvkit=build_default_sif_list()$sif_cnvkit,
    pon=build_default_reference_list()$HG19$panel$PCF_V3$variant$pon_cn_male,
    target="",
    antitarget="",
    output_name="",
    output_dir=".",
    gc=TRUE,
    edge=TRUE,
    rmask=TRUE,
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=1,mode="local",
    executor_id=make_unique_id("fixCNVkit"),
    task_name="fixCNVkit",time="48:0:0",
    update_time=60,
    wait=FALSE,hold=NULL
  ){

    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir)
    job=build_job(executor_id=executor_id,task_id=task_id)

    if(!is.null(rdata)){
      load(rdata)
      if(!is.null(selected)){
        target=sample_list[selected,]$targets
        antitarget=sample_list[selected,]$antitargets
      }
    }

    id=""
    if(output_name!=""){
      id=output_name
    }else{
      id=get_file_name(target)
    }

    add=""

    if(!gc){
      add=paste(add," --no-gc ")
    }

    if(!edge){
      add=paste(add," --no-edge ")
    }

    if(!rmask){
      add=paste(add," --no-rmask ")
    }


    out_file=paste0(out_file_dir,"/",id,".cnr")

    exec_code=paste("singularity exec -H ",paste0(getwd(),":/home "),sif_cnvkit,
    " cnvkit.py fix -o ",out_file, add, target, antitarget, pon)

    if(mode=="batch"){
        out_file_dir2=set_dir(dir=out_file_dir,name="batch")
        batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
        threads=threads,output_dir=out_file_dir2)
        exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
    }

    if(verbose){
        print_verbose(job=job,arg=argg,exec_code=exec_code)
    }

    error=execute_job(exec_code=exec_code)
    
    if(error!=0){
      stop("cnvkit failed to run due to unknown error.
      Check std error for more information.")
    }

    job_report=build_job_report(
      job_id=job,
      executor_id=executor_id,
      exec_code=exec_code, 
      task_id=task_id,
      input_args = argg,
      out_file_dir=out_file_dir,
      out_files=list(
        cnr=out_file)
    )


    if(wait&&mode=="batch"){
      job_validator(job=job_report$job_id,time=update_time,
      verbose=verbose,threads=threads)
    }

    return(job_report)

  }




#' Wrapper around parallel intergration for fix function from CNVkit
#'
#' This function wraps around fix function for CNVkit
#' This function correct sample coverage using a reference panel of normals files.
#' Additional parameters can be used to exclude regions and modify the average bin size.
#' 
#' 
#' For more information read:
#' https://cnvkit.readthedocs.io/en/stable/pipeline.html
#'
#' @param sif_cnvkit [REQUIRED] Path to cnvkit sif file.
#' @param ref_genoma [REQUIRED] Path to reference genome.
#' @param pon [REQUIRED] Path to panel of normals location.
#' @param targets [REQUIRED] Path to target CNN files.
#' @param antitargets [REQUIRED] Path to antitarget CNN files.
#' @param gc [OPTIONAL] Disble GC bias correction. Default FALSE.
#' @param edge [OPTIONAL] Disble edge correction. Default FALSE.
#' @param rmask [OPTIONAL] Disble repeat mask correction. Default FALSE.
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
#' @param output_dir [OPTIONAL] Path to the output directory.
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

  parallel_sample_fix_cnvkit=function(
    sif_cnvkit=build_default_sif_list()$sif_cnvkit,
    pon=build_default_reference_list()$HG19$panel$PCF_V3$variant$pon_cn_male,
    targets="",
    antitargets="",
    output_name="",
    output_dir=".",
    gc=TRUE,
    edge=TRUE,
    rmask=TRUE,
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=1,mode="local",
    executor_id=make_unique_id("parFixCNVkit"),
    task_name="parFixCNVkit",time="48:0:0",
    update_time=60,
    wait=FALSE,hold=NULL
  ){

    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir,name="cnr")
    tmp_dir=set_dir(dir=output_dir,name="tmp")
    job=build_job(executor_id=executor_id,task_id=task_id)

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

    sample_list=data.frame(targets=targets,
    antitargets=antitargets,
    names=Vectorize(get_file_name)(targets))
    row.names(sample_list)= sample_list$names
    
    
    if(mode=="local"){
    jobs_report[["steps"]][["par_sample_fix_cnvkit"]]<-
    parallel::mclapply(seq(1,nrow(sample_list)),FUN=function(row){
      job_report <- fix_cnvkit(
                sif_cnvkit=sif_cnvkit,
                pon=pon,
                target=sample_list[row,]$targets,
                antitarget=sample_list[row,]$antitargets,
                output_name=sample_list[row,]$names,
                output_dir=out_file_dir,
                gc=gc,
                edge=edge,
                rmask=rmask,
                verbose=verbose,
                batch_config=batch_config,
                executor_id=task_id
              )
    },mc.cores=threads)
    
    }else if(mode=="batch"){

          rdata_file=paste0(tmp_dir,"/",job,".samples.RData")
          output_dir=out_file_dir
          save(sample_list,bed,sif_cnvkit,
          gc,edge,rmask,output_dir,verbose,file = rdata_file)
          exec_code=paste0("Rscript -e \"ULPwgs::fix_cnvkit(rdata=\\\"",
          rdata_file,"\\\",selected=$SGE_TASK_ID)\"")
          out_file_dir2=set_dir(dir=out_file_dir,name="batch")
          batch_code=build_job_exec(job=job,time=time,ram=ram,
          threads=1,output_dir=out_file_dir2,
          hold=hold,array=nrow(sample_list))
          exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)

          if(verbose){
              print_verbose(job=job,arg=argg,exec_code=exec_code)
          }
          error=execute_job(exec_code=exec_code)
          if(error!=0){
              stop("gatk failed to run due to unknown error.
              Check std error for more information.")
          }
         
         jobs_report[["steps"]][["par_sample_fix_cnvkit"]]<- build_job_report(
              job_id=job,
              executor_id=executor_id,
              exec_code=exec_code, 
              task_id=task_id,
              input_args=argg,
              out_file_dir=out_file_dir,
              out_files=list(
                  cnr=paste0(out_file_dir,"/",sample_list$names,".cnr")
              )
        )
    }

      return(jobs_report)


  }






#' Wrapper around reference function from CNVkit
#'
#' This function wraps around target function for CNVkit
#' This function generates an target BED file from an input target file. 
#' Additional parameters can be used to exclude regions and modify the average bin size
#' 
#' 
#' For more information read:
#' https://cnvkit.readthedocs.io/en/stable/pipeline.html
#'
#' @param sif_cnvkit [REQUIRED] Path to cnvkit sif file.oms.
#' @param seg_method [OPTIONAL] Method to use for segmentation. Default cbs. Options ["cbs","flasso","haar","none","hmm","hmm-tumor","hmm-germline"]
#' @param smooth [OPTIONAL] Smooth before CBS. Default TRUE
#' @param drop_low_coverage [OPTIONAL] Drop low coverage bins before segmentation. Default TRUE
#' @param drop_outliers [OPTIONAL] Drop outlier bins before segmentation. Default TRUE
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
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

 segment_cnvkit=function(
    rdata=NULL,
    selected=NULL,
    sif_cnvkit=build_default_sif_list()$sif_cnvkit,
    cnr="",
    seg_method="cbs",
    output_name="",
    output_dir=".",
    smooth=FALSE,
    drop_low_coverage=FALSE,
    drop_outliers=10,
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=1,mode="local",
    executor_id=make_unique_id("segmentCNVkit"),
    task_name="segmentCNVkit",time="48:0:0",
    update_time=60,
    wait=FALSE,hold=NULL
  ){

    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir)
    job=build_job(executor_id=executor_id,task_id=task_id)

    if(!is.null(rdata)){
      load(rdata)
      if(!is.null(selected)){
        cnr=sample_list[selected]
      }
    }

    id=""
    if(output_name!=""){
      id=output_name
    }else{
      id=get_file_name(cnr)
    }

    add=""

    if(smooth){
      add=paste(add," --smooth-cbs ")
    }


     if(drop_low_coverage){
      add=paste(add," --drop-low-coverage ")
    }


    if(drop_outliers){
        drop_outliers=paste0("  --drop-outliers ",drop_outliers)
    }
  
    out_file=paste0(out_file_dir,"/",id,".cns")

    exec_code=paste("singularity exec -H ",paste0(getwd(),":/home "),sif_cnvkit,
    " cnvkit.py segment -p ",threads,drop_outliers,paste0(" -o ",out_file),add, cnr)

    if(mode=="batch"){
        out_file_dir2=set_dir(dir=out_file_dir,name="batch")
        batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
        threads=threads,output_dir=out_file_dir2)
        exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
    }

    if(verbose){
        print_verbose(job=job,arg=argg,exec_code=exec_code)
    }

    error=execute_job(exec_code=exec_code)
    
    if(error!=0){
      stop("cnvkit failed to run due to unknown error.
      Check std error for more information.")
    }

    job_report=build_job_report(
      job_id=job,
      executor_id=executor_id,
      exec_code=exec_code, 
      task_id=task_id,
      input_args = argg,
      out_file_dir=out_file_dir,
      out_files=list(
        cns=out_file)
    )


    if(wait&&mode=="batch"){
      job_validator(job=job_report$job_id,time=update_time,
      verbose=verbose,threads=threads)
    }

    return(job_report)

  }


#' Wrapper around reference function from CNVkit
#'
#' This function wraps around target function for CNVkit
#' This function generates an target BED file from an input target file. 
#' Additional parameters can be used to exclude regions and modify the average bin size
#' 
#' 
#' For more information read:
#' https://cnvkit.readthedocs.io/en/stable/pipeline.html
#'
#' @param sif_cnvkit [REQUIRED] Path to cnvkit sif file.oms.
#' @param ref_genome [REQUIRED] Path to CNN file for target regions. Default none.
#' @param seg_method [OPTIONAL] Method to use for segmentation. Default cbs. Options ["cbs","flasso","haar","none","hmm","hmm-tumor","hmm-germline"]
#' @param smooth [OPTIONAL] Smooth before CBS. Default TRUE
#' @param drop_low_coverage [OPTIONAL] Drop low coverage bins before segmentation. Default TRUE
#' @param drop_low_coverage [OPTIONAL] Drop low coverage bins before segmentation. Default TRUE
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
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

 parallel_samples_segment_cnvkit=function(
    sif_cnvkit=build_default_sif_list()$sif_cnvkit,
    cnrs="",
    seg_method="cbs",
    output_name="",
    output_dir=".",
    smooth=FALSE,
    drop_low_coverage=FALSE,
    drop_outliers=10,
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=1,mode="local",
    executor_id=make_unique_id("segmentCNVkit"),
    task_name="segmentCNVkit",time="48:0:0",
    update_time=60,
    wait=FALSE,hold=NULL
  ){

    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir,name="cns")
    tmp_dir=set_dir(dir=output_dir,name="tmp")
    job=build_job(executor_id=executor_id,task_id=task_id)

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

    sample_list=cnrs
    names(sample_list)=Vectorize(get_file_name)(sample_list)
    
    
    if(mode=="local"){
    jobs_report[["steps"]][["par_sample_segment_cnvkit"]]<-
    parallel::mclapply(sample_list,FUN=function(sample){
      job_report <-segment_cnvkit(
                sif_cnvkit=sif_cnvkit,
                cnr=sample,
                seg_method=seg_method,
                smooth=smooth,
                drop_low_coverage=drop_low_coverage,
                drop_outliers=drop_outliers,
                output_name=get_file_name(sample),
                output_dir=out_file_dir,
                verbose=verbose,
                batch_config=batch_config,
                executor_id=task_id
              )
    },mc.cores=threads)
    
    }else if(mode=="batch"){
          rdata_file=paste0(tmp_dir,"/",job,".samples.RData")
          output_dir=out_file_dir
          save(sample_list,seg_method,sif_cnvkit,smooth,
          drop_low_coverage,drop_outliers,output_dir,verbose,file = rdata_file)
          exec_code=paste0("Rscript -e \"ULPwgs::segment_cnvkit(rdata=\\\"",
          rdata_file,"\\\",selected=$SGE_TASK_ID)\"")
          out_file_dir2=set_dir(dir=out_file_dir,name="batch")
          batch_code=build_job_exec(job=job,time=time,ram=ram,
          threads=1,output_dir=out_file_dir2,
          hold=hold,array=length(sample_list))
          exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)

          if(verbose){
              print_verbose(job=job,arg=argg,exec_code=exec_code)
          }
          error=execute_job(exec_code=exec_code)
          if(error!=0){
              stop("gatk failed to run due to unknown error.
              Check std error for more information.")
          }
         
         jobs_report[["steps"]][["par_sample_segment_cnvkit"]]<- build_job_report(
              job_id=job,
              executor_id=executor_id,
              exec_code=exec_code, 
              task_id=task_id,
              input_args=argg,
              out_file_dir=out_file_dir,
              out_files=list(
                  cnn=paste0(out_file_dir,"/",sample_list$names,".cnr")
              )
        )
    }
      return(jobs_report)
  }







#' Wrapper around reference scatter from CNVkit
#'
#' This function wraps around target function for CNVkit
#' This function generates an target BED file from an input target file. 
#' Additional parameters can be used to exclude regions and modify the average bin size
#' 
#' 
#' For more information read:
#' https://cnvkit.readthedocs.io/en/stable/pipeline.html
#'
#' @param sif_cnvkit [REQUIRED] Path to cnvkit sif file.oms.
#' @param trend [OPTIONAL] Smooth before CBS. Default TRUE
#' @param range [OPTIONAL] Range to show in chr:start-end format. Default none
#' @param range_list [OPTIONAL] Path to BED file with range list. Default none
#' @param genes [OPTIONAL] List of genes to focus on. Default none.
#' @param margin_width [OPTIONAL] Default gene margin width. Default 1000000
#' @param plot_by_bin [OPTIONAL] Assume all bins are the same size. Default FALSE
#' @param trend [OPTIONAL] Show trendline. Default TRUE.
#' @param antitarget_symbol [OPTIONAL] Show antitarget probes with this symbol. Default "@"
#' @param segment_colour [OPTIONAL] Show antitarget probes with this symbol. Default "@"
#' @param title [OPTIONAL] Sample title. Default NULL.
#' @param y_max [OPTIONAL] Y max range. Default NULL.
#' @param y_min[OPTIONAL] Y min range. Default NULL.
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
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

  scatter_cnvkit=function(
    rdata=NULL,
    selected=NULL,
    sif_cnvkit=build_default_sif_list()$sif_cnvkit,
    chrs=build_default_chr_list()$canonical,
    cnr="",
    cns="",
    output_name="",
    output_dir=".",
    range="",
    range_list="",
    genes="",
    margin_width=1000000,
    plot_by_bin=FALSE,
    trend=TRUE,
    antitarget_symbol="@",
    segment_colour="red",
    title=NULL,
    y_max=NULL,
    y_min=NULL,
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=1,mode="local",
    executor_id=make_unique_id("scatterCNVkit"),
    task_name="scatterCNVkit",time="48:0:0",
    update_time=60,
    wait=FALSE,hold=NULL
  ){

    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir)
    job=build_job(executor_id=executor_id,task_id=task_id)

    if(!is.null(rdata)){
      load(rdata)
      if(!is.null(selected)){
        cnr=sample_list[selected]
      }
    }

    id=""
    if(output_name!=""){
      id=output_name
    }else{
      id=get_file_name(cnr)
    }


    if(range!=""){
      range=paste0(" -c ",range)
    }

     if(segment_colour!=""){
     segment_colour=paste0(" --segment-color ",segment_colour)
    }

    
    if(genes!=""){
      genes=paste0(" --gene ",genes)
    }

    if(range_list!=""){
      range_list=paste0(" -l ",range_list)
    }

    add=""

     if(!is.null(title)){
      title=paste0(" --title ",title)
    }
     if(!is.null(y_max)){
      y_max=paste0(" --y-max ",y_max)
    }


    if(!is.null(y_min)){
      y_min=paste0(" --y-min ",y_min)
    }
    

    if(trend){
      add=paste(add," --trend ")
    }


     if(plot_by_bin){
      add=paste(add," --by-bin ")
    }


    if(margin_width!=""){
      margin_width=paste0(" --width ",margin_width)
    }


    out_file=paste0(out_file_dir,"/",id,".scatter.pdf")

    cnr_tmp=paste0(cnr,".tmp")
    cns_tmp=""
    exec_code=paste0("cat ",cnr," | head -n 1 > ",cnr_tmp," && cat ",cnr," | grep \"",
    paste0(paste0("^",chrs),collapse="\\|"),"\" >> ",cnr_tmp)
    if(cns!=""){
      cns_tmp=paste0(cns,".tmp")
      cns_code=paste(" -s ",cns_tmp)
      exec_code=paste0(exec_code," && cat ",cns," |  head -n 1 > ",cns_tmp," && ",paste0("cat ",cns," | grep \" ",
      paste0(paste0("^",unlist(chrs),collapse="\\|")," \" >> ",cns_tmp)))
    }

    exec_code=paste(exec_code,";singularity exec -H ",paste0(getwd(),":/home "),sif_cnvkit,
    " cnvkit.py scatter -o ",out_file,cns_code,add,title,segment_colour,
    y_max,y_min,range,margin_width,range_list,cnr_tmp,
    "&& rm ",cnr_tmp,cns_tmp)

    if(mode=="batch"){
        out_file_dir2=set_dir(dir=out_file_dir,name="batch")
        batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
        threads=threads,output_dir=out_file_dir2)
        exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
    }

    if(verbose){
        print_verbose(job=job,arg=argg,exec_code=exec_code)
    }

    error=execute_job(exec_code=exec_code)
    
    if(error!=0){
      stop("cnvkit failed to run due to unknown error.
      Check std error for more information.")
    }

    job_report=build_job_report(
      job_id=job,
      executor_id=executor_id,
      exec_code=exec_code, 
      task_id=task_id,
      input_args = argg,
      out_file_dir=out_file_dir,
      out_files=list(
        scatter=out_file)
    )


    if(wait&&mode=="batch"){
      job_validator(job=job_report$job_id,time=update_time,
      verbose=verbose,threads=threads)
    }

    return(job_report)

  }






#' Wrapper around reference scatter from CNVkit
#'
#' This function wraps around target function for CNVkit
#' This function generates an target BED file from an input target file. 
#' Additional parameters can be used to exclude regions and modify the average bin size
#' 
#' 
#' For more information read:
#' https://cnvkit.readthedocs.io/en/stable/pipeline.html
#'
#' @param sif_cnvkit [REQUIRED] Path to cnvkit sif file.oms.
#' @param cnr [REQUIRED] Path to CNR as produced by fix_cnvkit command.
#' @param cns [OPTIONAL] Path to CNS as produced by segment_cnvkit command.
#' @param cn_thr [OPTIONAL] Copy number threshold. Default 0.5
#' @param min_probes [OPTIONAL] Minimum number of probes to show a gene. Default 3
#' @param male_reference [OPTIONAL]  Assume inputs were normalized to a male reference. Default FALSE
#' @param gender [OPTIONAL] Assume sample gender. Default male
#' @param shift [OPTIONAL] Adjust the X and Y chromosomes according to sample sex. Default TRUE
#' @param title [OPTIONAL] Sample ID. Default NULL.
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
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

  diagram_cnvkit=function(
    rdata=NULL,
    selected=NULL,
    chrs=build_default_chr_list()$canonical,
    sif_cnvkit=build_default_sif_list()$sif_cnvkit,
    cnr="",
    cns="",
    output_name="",
    output_dir=".",
    cn_thr=0.5,
    min_probes=3,
    male_reference=FALSE,
    gender="male",
    shift=TRUE,
    title=NULL,
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=1,mode="local",
    executor_id=make_unique_id("diagramCNVkit"),
    task_name="diagramCNVkit",time="48:0:0",
    update_time=60,
    wait=FALSE,hold=NULL
  ){

    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir)
    job=build_job(executor_id=executor_id,task_id=task_id)

    if(!is.null(rdata)){
      load(rdata)
      if(!is.null(selected)){
        cnr=sample_list[selected]
      }
    }

    id=""
    if(output_name!=""){
      id=output_name
    }else{
      id=get_file_name(cnr)
    }

    if(!is.null(gender)){
      gender=paste0(" -x ",gender)
    }





    if(cn_thr!=""){

      cn_thr=paste0(" -t ",cn_thr)

    }

    if(min_probes!=""){

      min_probes=paste0(" -m ",min_probes)

    }

    if(!is.null(title)){
      title=paste0(" --title ",title)
    }


    add=""

    if(!shift){
      add=paste(add," --no-shift-xy ")
    }

    if(male_reference){
      add=paste(add," -y ")
    }
    tmp_cns=paste0(cns,".tmp")
    tmp_cnr=paste0(cnr,".tmp")
 

    out_file=paste0(out_file_dir,"/",id,".diagram.pdf")
    cnr_tmp=paste0(cnr,".tmp")
    cns_tmp=""
    exec_code=paste0("cat ",cnr," | head -n 1 > ",cnr_tmp," && cat ",cnr," | grep \"",
    paste0(paste0("^",chrs),collapse="\\|"),"\" >> ",cnr_tmp)
    if(cns!=""){
      cns_tmp=paste0(cns,".tmp")
      cns_code=paste(" -s ",cns_tmp)
      exec_code=paste0(exec_code," && cat ",cns," |  head -n 1 > ",cns_tmp," && ",paste0("cat ",cns," | grep \" ",
      paste0(paste0("^",unlist(chrs),collapse="\\|")," \" >> ",cns_tmp)))
    }
    exec_code=paste(exec_code,"; singularity exec -H ",paste0(getwd(),":/home "),sif_cnvkit,
    " cnvkit.py diagram -o ",out_file,cns_code,gender,add,title,min_probes,cn_thr,
    cnr_tmp,"&& rm ",cnr_tmp,cns_tmp)

    if(mode=="batch"){
        out_file_dir2=set_dir(dir=out_file_dir,name="batch")
        batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
        threads=threads,output_dir=out_file_dir2)
        exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
    }

    if(verbose){
        print_verbose(job=job,arg=argg,exec_code=exec_code)
    }

    error=execute_job(exec_code=exec_code)
    
    if(error!=0){
      stop("cnvkit failed to run due to unknown error.
      Check std error for more information.")
    }

    job_report=build_job_report(
      job_id=job,
      executor_id=executor_id,
      exec_code=exec_code, 
      task_id=task_id,
      input_args = argg,
      out_file_dir=out_file_dir,
      out_files=list(
        diagram=out_file)
    )


    if(wait&&mode=="batch"){
      job_validator(job=job_report$job_id,time=update_time,
      verbose=verbose,threads=threads)
    }

    return(job_report)

  }

