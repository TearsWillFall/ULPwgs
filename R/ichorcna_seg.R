# file:   segmentation.R
# author: Gavin Ha, Ph.D.
#         Justin Rhoades
#               Dana-Farber Cancer Institute
#               Broad Institute
# contact: <gavinha@broadinstitute.org>
# ULP-WGS website: http://www.broadinstitute.org/~gavinha/ULP-WGS/
# HMMcopy website: http://compbio.bccrc.ca/software/hmmcopy/ and https://www.bioconductor.org/packages/release/bioc/html/HMMcopy.html
# date:   Oct 26, 2016
# description: Hidden Markov model (HMM) to analyze Ultra-low pass whole genome sequencing (ULP-WGS) data.
# This script is the main script to run the HMM.

runHMMsegment <- function(x,
    validInd = NULL, dataType = "log2", param = NULL, 
    chrTrain = c(1:22), maxiter = 50, estimateNormal = TRUE, estimatePloidy = TRUE, 
    estimatePrecision = TRUE, estimateSubclone = FALSE, estimateTransition = TRUE,
    estimateInitDist = TRUE, logTransform = FALSE, verbose = TRUE) {
    chr <- as.factor(GenomeInfoDb::seqnames(x[[1]]))
	# setup columns for multiple samples #
	dataMat <- as.matrix(as.data.frame(lapply(x, function(y) { GenomicRanges::mcols(x[[1]])[, dataType] })))
	
	# normalize by median and log data #
	if (logTransform){
    dataMat <- apply(dataMat, 2, function(x){ log(x / median(x, na.rm = TRUE)) })
	}else{
	  dataMat <- log(2^dataMat)
	}
	## update variable x with loge instead of log2
  for (i in 1:length(x)){
    GenomicRanges::mcols(x[[i]])[, dataType] <- dataMat[, i]
  }
  if (!is.null(chrTrain)) {
		chrInd <- chr %in% chrTrain
  }else{
  	chrInd <- !logical(length(chr))
  }

  if (!is.null(validInd)){
    chrInd <- chrInd & validInd
  }  

	if (is.null(param)){
		param <- getDefaultParameters(dataMat[chrInd])
	}
	#if (param$n_0 == 0){
	#	param$n_0 <- .Machine$double.eps
	#}
	####### RUN EM ##########
  convergedParams <- runEMs(dataMat, chr, chrInd, param, maxiter, 
      verbose, estimateNormal = estimateNormal, estimatePloidy = estimatePloidy, 
      estimateSubclone = estimateSubclone, estimatePrecision = estimatePrecision, 
      estimateTransition = estimateTransition, estimateInitDist = estimateInitDist)
  # Calculate likelihood using converged params
 # S <- param$numberSamples
 # K <- length(param$ct)
 # KS <- K ^ S
 # py <- matrix(0, KS, nrow(dataMat))
 # iter <- convergedParams$iter
  # lambdasKS <- as.matrix(expand.grid(as.data.frame(convergedParams$lambda[, , iter])))
  # for (ks in 1:KS) {
  #   probs <- tdistPDF(dataMat, convergedParams$mus[ks, , iter], lambdasKS[ks, ], param$nu)
  #   py[ks, ] <- apply(probs, 1, prod) # multiply across samples for each data point to get joint likelihood.
  # }
  # 
  viterbiResults <- runViterbis(convergedParams, chr)

  # setup columns for multiple samples #
  segs <- segmentData(dataGR=x, validInd=validInd, 
  states=viterbiResults$states, convergedParams=convergedParams)
  #output$segs <- processSegments(output$segs, chr, start(x), end(x), x$DataToUse)
  names <- c("HOMD","HETD","NEUT","GAIN","AMP","HLAMP",paste0(rep("HLAMP", 8), 2:25))
  #if (c(0) %in% param$ct){ #if state 0 HOMD is IN params#
  	#names <- c("HOMD", names)
  	# shift states to start at 2 (HETD)
    #tmp <- lapply(segs, function(x){ x$state <- x$state + 1; x})
    #viterbiResults$states <- as.numeric(viterbiResults$states) + 1
	#}
	### PUTTING TOGETHER THE COLUMNS IN THE OUTPUT ###
  cnaList <- list()
  S <- length(x)
  for (s in 1:S){
    id <- names(x)[s]
    copyNumber <- param$jointCNstates[viterbiResults$state, s]
    subclone.status <- param$jointSCstatus[viterbiResults$state, s]
  	cnaList[[id]] <- data.frame(cbind(sample = as.character(id), 
                  chr = as.character(GenomeInfoDb::seqnames(x[[s]])),	
                  start = start(x[[s]]), end = end(x[[s]]), 
                  copy.number = copyNumber,
                  event = names[copyNumber + 1], 
                  logR = round(log2(exp(dataMat[,s])), digits = 4),
                  subclone.status = as.numeric(subclone.status)
  	))

    cnaList[[id]] <- transform(cnaList[[id]], 
                              start = as.integer(as.character(start)),
                              end = as.integer(as.character(end)), 
                              copy.number = as.numeric(copy.number),
                              logR = as.numeric(as.character(logR)),
                              subclone.status = as.numeric(subclone.status))
  
  	## order by chromosome ##
  	chrOrder <- unique(chr) #c(1:22,"X","Y")
  	cnaList[[id]] <- cnaList[[id]][order(match(cnaList[[id]][, "chr"],chrOrder)),]
  	## remove MT chr ##
    cnaList[[id]] <- cnaList[[id]][cnaList[[id]][,"chr"] %in% chrOrder, ]
    
    ## segment mean loge -> log2
    #segs[[s]]$median.logR <- log2(exp(segs[[s]]$median.logR))
    segs[[s]]$median <- log2(exp(segs[[s]]$median))
    ## add subclone status
    segs[[s]]$subclone.status <-  param$jointSCstatus[segs[[s]]$state, s]
  }	
  convergedParams$segs <- segs
  return(list(cna = cnaList, results = convergedParams, viterbiResults = viterbiResults))
}

getTransitionMatrix <- function(K, e, strength){
  A <- matrix(0, K, K)
  for (j in 1:K) {
    A[j, ] <- (1 - e[1]) / (K - 1)
    A[j, j] <- e[1]
  }
  A <- normalize(A)
  A_prior <- A
  dirPrior <- A * strength[1]
  return(list(A=A, dirPrior=dirPrior))
}

getDefaultParameters <- function(x, maxCN = 5, ct.sc = NULL, ploidy = 2, e = 0.9999999, e.sameState = 10, strength = 10000000, includeHOMD = FALSE){
  if (includeHOMD){
    ct <- 0:maxCN
  }else{
    ct <- 1:maxCN
  }
	param <- list(
		strength = strength, e = e,
		ct = c(ct, ct.sc),
		ct.sc.status = c(rep(FALSE, length(ct)), rep(TRUE, length(ct.sc))),
		phi_0 = 2, alphaPhi = 4, betaPhi = 1.5,
		n_0 = 0.5, alphaN = 2, betaN = 2,
		sp_0 = 0.5, alphaSp = 2, betaSp = 2,
		lambda = as.matrix(rep(100, length(ct)+length(ct.sc)), ncol=1),
		nu = 2.1,
		kappa = rep(75, length(ct)), 
		alphaLambda = 5
	)
	K <- length(param$ct)
  ## initialize hyperparameters for precision using observed data ##
	if (!is.null(dim(x))){ # multiple samples (columns)
    param$numberSamples <- ncol(x)
    #betaLambdaVal <- ((apply(x, 2, function(x){ sd(diff(x), na.rm=TRUE) }) / sqrt(length(param$ct))) ^ 2)
    betaLambdaVal <- ((apply(x, 2, sd, na.rm = TRUE) / sqrt(length(param$ct))) ^ 2)   
	}else{ # only 1 sample
	  param$numberSamples <- 1
	  betaLambdaVal <- ((sd(x, na.rm = TRUE) / sqrt(length(param$ct))) ^ 2)
	}
	param$betaLambda <- matrix(betaLambdaVal, ncol = param$numberSamples, nrow = length(param$ct), byrow = TRUE)
  param$alphaLambda <- rep(param$alphaLambda, K)
  
	# increase prior precision for -1, 0, 1 copies at ploidy
	#param$lambda[param$ct %in% c(1,2,3)] <- 1000 # HETD, NEUT, GAIN
	#param$lambda[param$ct == 4] <- 100 
	#param$lambda[which.max(param$ct)] <- 50 #highest CN
	#param$lambda[param$ct == 0] <- 1 #HOMD
	S <- param$numberSamples
	logR.var <- 1 / ((apply(x, 2, sd, na.rm = TRUE) / sqrt(length(param$ct))) ^ 2)
	if (!is.null(dim(x))){ # multiple samples (columns)
		param$lambda <- matrix(logR.var, nrow=K, ncol=S, byrow=T, dimnames=list(c(),colnames(x)))
	}else{ # only 1 sample    
		#logR.var <- 1 / ((sd(x, na.rm = TRUE) / sqrt(length(param$ct))) ^ 2)
    param$lambda <- matrix(logR.var, length(param$ct))
    param$lambda[param$ct %in% c(2)] <- logR.var 
    param$lambda[param$ct %in% c(1,3)] <- logR.var 
    param$lambda[param$ct >= 4] <- logR.var / 5
    param$lambda[param$ct == max(param$ct)] <- logR.var / 15
    param$lambda[param$ct.sc.status] <- logR.var / 10
  }
  # define joint copy number states #
  param$jointCNstates <- expand.grid(rep(list(param$ct), S))
  param$jointSCstatus <- expand.grid(rep(list(param$ct.sc.status), S))
  colnames(param$jointCNstates) <- paste0("Sample.", 1:param$numberSamples)
  colnames(param$jointSCstatus) <- paste0("Sample.", 1:param$numberSamples)
  
	# Initialize transition matrix to the prior
	txn <- getTransitionMatrix(K ^ S, e, strength)
  ## set higher transition probs for same CN states across samples ##
  # joint states where at least "tol" fraction of samples with the same CN state
	#apply(param$jointCNstates, 1, function(x){ sum(duplicated(as.numeric(x))) > 0 })
  cnStateDiff <- apply(param$jointCNstates, 1, function(x){ (abs(max(x) - min(x)))})
  if (e.sameState > 0 & S > 1){
		txn$A[, cnStateDiff == 0] <- txn$A[, cnStateDiff == 0] * e.sameState * K 
		txn$A[, cnStateDiff >= 3] <- txn$A[, cnStateDiff >=3]  / e.sameState / K
	}
  for (i in 1:nrow(txn$A)){
    for (j in 1:ncol(txn$A)){
      if (i == j){
        txn$A[i, j] <- e
      }
    }
  }
  txn$A <- normalize(txn$A)
	param$A <- txn$A
	param$dirPrior <- txn$A * strength[1] 
  param$A[, param$ct.sc.status] <- param$A[, param$ct.sc.status] / 10
  param$A <- normalize(param$A)
  param$dirPrior[, param$ct.sc.status] <- param$dirPrior[, param$ct.sc.status] / 10
  
  if (includeHOMD){
    K <- length(param$ct)
    param$A[1, 2:K] <- param$A[1, 2:K] * 1e-5; param$A[2:K, 1] <- param$A[2:K, 1] * 1e-5;
    param$A[1, 1] <- param$A[1, 1] * 1e-5
    param$A <- normalize(param$A); param$dirPrior <- param$A * param$strength
  }

  param$kappa <- rep(75, K ^ S)
  param$kappa[cnStateDiff == 0] <- param$kappa[cnStateDiff == 0] + 125
	param$kappa[cnStateDiff >=3] <- param$kappa[cnStateDiff >=3] - 50
	param$kappa[which(rowSums(param$jointCNstates==2) == S)] <- 800
  
  return(param)
}


segmentData <- function(dataGR, validInd, states, convergedParams, dataType="log2"){
  if (sum(convergedParams$param$ct == 0) ==0){
  	includeHOMD <- FALSE
  }else{
  	includeHOMD <- TRUE
  }
  if (!includeHOMD){
    names <- c("HETD","NEUT","GAIN","AMP","HLAMP",paste0(rep("HLAMP", 8), 2:25))
  }else{
    names <- c("HOMD","HETD","NEUT","GAIN","AMP","HLAMP",paste0(rep("HLAMP", 8), 2:25))
  }
  states <- states[validInd]
  S <- length(dataGR)
  jointStates <- convergedParams$param$jointCNstates
  jointSCstatus <- convergedParams$param$jointSCstatus
  colNames <- c("seqnames", "start", "end", dataType)
  segList <- list()
 
  for (i in 1:S){
  	id <- names(dataGR)[i]
  	dataIn <- dataGR[[i]][validInd, ]
    rleResults <- t(sapply(S4Vectors::runValue(GenomeInfoDb::seqnames(dataIn)), function(x){
      ind <- as.character(GenomeInfoDb::seqnames(dataIn)) == x
      r <- rle(states[ind])
    }))
    rleLengths <- unlist(rleResults[, "lengths"])
    rleValues <- unlist(rleResults[, "values"])
    sampleDF <- as.data.frame(dataIn)
    numSegs <- length(rleLengths)
    segs <- as.data.frame(matrix(NA, ncol = 7, nrow = numSegs, 
                   dimnames = list(c(), c("chr", "start", "end", "state", "event", "median", "copy.number"))))
    prevInd <- 0
    
   
    for (j in 1:numSegs){
      start <- prevInd + 1
      end <- prevInd + rleLengths[j]
      segDF <- sampleDF[start:end, colNames]
      prevInd <- end
      numR <- nrow(segDF)
      segs[j, "chr"] <- as.character(segDF[1, "seqnames"])
      segs[j, "start"] <- segDF[1, "start"]
      segs[j, "state"] <- rleValues[j]
      segs[j, "copy.number"] <- jointStates[rleValues[j], i]
      if (segDF[1, "seqnames"] == segDF[numR, "seqnames"]){
        segs[j, "end"] <- segDF[numR, "end"]
        segs[j, "median"] <- round(median(segDF[,dataType], na.rm = TRUE), digits = 6)
        if (includeHOMD){
        	segs[j, "event"] <- names[segs[j, "copy.number"] + 1]
        }else{
        	segs[j, "event"] <- names[segs[j, "copy.number"]]
        }
      }else{ # segDF contains 2 different chromosomes
        print(j)
      }                                      
    }
    segList[[id]] <- segs
  }
  return(segList)
}
    

runViterbis <- function(convergedParams, chr){
  message("runViterbi: Segmenting and classifying")
  chrs <- levels(chr)
  chrsI <- vector('list', length(chrs))
  # initialise the chromosome index and the init state distributions
  for(i in 1:length(chrs)) {
    chrsI[i] <- list(which(chr == chrs[i]))
  }
  segs <- vector('list', length(chrs))
  py <- convergedParams$py
  N <- ncol(py)
  Z <- rep(0, N)
  convergeIter <- convergedParams$iter
  piG <- convergedParams$pi[, convergeIter]
  A <- convergedParams$A


  for(c in 1:length(chrsI)) {
    I <- chrsI[[c]]
    output <- .Call("viterbi", log(piG), log(A), log(py[, I]), PACKAGE = "HMMcopy")
    Z[I] <- output$path
    segs[[c]] <- output$seg
  }
  return(list(segs=segs, states=Z))
}


# Normalize a given array to sum to 1
normalize <- function(A) {
	vectorNormalize <- function(x){ x / (sum(x) + (sum(x) == 0)) }
	if (length(dim(A)) < 2){
  	M <- vectorNormalize(A)
  }else{
  	M <- t(apply(A, 1, vectorNormalize))
  }
  return(M);
}


# processSegments <- function(seg, chr, start, end, copy) {
#   segment <- data.frame()
#   chromosomes <- levels(chr)
#   for (i in 1:length(chromosomes)) {
#     seg_length = dim(seg[[i]])[1]
#     chr_name <- rep(chromosomes[i], seg_length)
#     chr_index <- which(chr == chromosomes[i])
#     chr_start <- start[chr_index][seg[[i]][, 1]]
#     chr_stop <- end[chr_index][seg[[i]][, 2]]
#     chr_state <- seg[[i]][, 3]
#     chr_median <- rep(0, seg_length)
#     for(j in 1:seg_length) {
#       chr_median[j] <-
#         median(na.rm = TRUE, log2(exp(copy[chr_index][seg[[i]][j, 1]:seg[[i]][j, 2]])))
#     }
#     segment <- rbind(segment, cbind(chr = chr_name,
#       start = as.numeric(chr_start), end = chr_stop, state = chr_state,
#       median = chr_median))
#   }
#   segment <- transform(segment, start = as.numeric(as.character(start)),
#     end = as.numeric(as.character(end)), as.numeric(as.character(state)),
#     median = as.numeric(as.character(median)))
#   return(segment)
# }


# file:   EM.R
# author: Gavin Ha, Ph.D.
#         Justin Rhoades
#               Dana-Farber Cancer Institute
#               Broad Institute
# contact: <gavinha@broadinstitute.org>
# ULP-WGS website: http://www.broadinstitute.org/~gavinha/ULP-WGS/
# HMMcopy website: http://compbio.bccrc.ca/software/hmmcopy/ and https://www.bioconductor.org/packages/release/bioc/html/HMMcopy.html
# date:   Oct 26, 2016
# description: Hidden Markov model (HMM) to analyze Ultra-low pass whole genome sequencing (ULP-WGS) data.
# This script is the main script to run the HMM.

runEMs <- function(copy, chr, chrTrain, param, maxiter, verbose = TRUE, 
									estimateNormal = TRUE, estimatePloidy = TRUE, estimatePrecision = TRUE,
                  estimateTransition = TRUE, estimateInitDist = TRUE, estimateSubclone = TRUE,
                  likChangeConvergence = 1e-3) {
  
  library(HMMcopy)
  if (nrow(copy) != length(chr) || nrow(copy) != length(chrTrain)) {
    stop("runEM: Length of inputs do not match for one of: copy, chr, chrTrain")
  }
  
  if (is.null(param$ct) || is.null(param$lambda) || is.null(param$nu) ||
      is.null(param$kappa)) {
    stop("runEM: Parameter missing, ensure all parameters exist as columns in",
         "data frame: ct, lambda, nu, kappa")
  }
  
  S <- param$numberSamples
  K <- length(param$ct)
  Z <- sum(param$ct.sc) #number of subclonal states (# repeated copy number states)
  KS <- K ^ S
  N <- nrow(copy)
  rho <- matrix(0, KS, N)
  py <- matrix(0, KS, N)            # Local evidence
  mus <- array(0, dim = c(KS, S, maxiter))     # State means
  lambdas <- array(0, dim = c(K, S, maxiter))  # State Variances
  phi <- matrix(NA, length(param$phi_0), maxiter)					 # Tumour ploidy
  n <- matrix(NA, S, maxiter)						 # Normal contamination
  sp <- matrix(NA, S, maxiter)     # cellular prevalence (tumor does not contain event)
  piG <- matrix(0, KS, maxiter)     # Initial state distribution
  converged <- FALSE               # Flag for convergence
  Zcounts <- matrix(0, KS, KS)
  loglik <- rep(0, maxiter)
  
  ptmTotal <- proc.time() # start total timer
  
  # SET UP
  # Set up the chromosome indices and make cell array of chromosome indicies
  chrs <- levels(chr)              # WARNING, gets resorted here...
  chrsI <- vector('list', length(chrs))
  # initialise the chromosome index and the init state distributions
  for(i in 1:length(chrs)) {
    chrsI[i] <- list(which(chr == chrs[i]))
  }
  
  # INITIALIZATION
  if (verbose) { message("runEM: Initialization") }
  i <- 1
  piG[, i] <- normalize(param$kappa)
  n[, i] <- param$n_0
  sp[, i] <- param$sp_0
  phi[, i] <- param$phi_0
  lambdas[, , i] <- param$lambda #matrix(param$lambda, nrow = K, ncol = S, byrow = TRUE)
  lambdasKS <- as.matrix(expand.grid(as.data.frame(lambdas[, , i]))) #as.matrix(expand.grid(rep(list(param$lambda), S)))
  mus[, , i] <- as.matrix(get2and3ComponentMixture(param$jointCNstates, param$jointSCstatus, n[, i], sp[, i], phi[, i]))
  
  # Likelihood #
  for (ks in 1:KS) {
    probs <- tdistPDF(copy, mus[ks, , i], lambdasKS[ks, ], param$nu)
    py[ks, ] <- apply(probs, 1, prod) # multiply across samples for each data point to get joint likelihood.
  }
  
  # initialize transition prior #
  A <- normalize(param$A)
  A_prior <- A
  dirPrior <- param$dirPrior
  
  loglik[i] <- -Inf
  
  while(!converged && (i < maxiter)) {
    ptm <- proc.time()
    #if (verbose) { message(paste("runEM: EM iteration:", i, "Log likelihood:", loglik[i])) }
    i <- i + 1
    
    ################ E-step ####################
    if (verbose) { message("runEM iter", i-1 , ": Expectation") }
    Zcounts <- matrix(0, KS, KS)
    for (j in 1:length(chrsI)) {
      I <- intersect(chrsI[[j]], which(chrTrain))
      if (length(I) > 0){   
        output <- .Call("forward_backward", piG[, i - 1], A, py[, I],PACKAGE = "HMMcopy")
        rho[, I] <- output$rho
        loglik[i] <- loglik[i] + output$loglik
        Zcounts <- Zcounts + t(colSums(aperm(output$xi, c(3, 2, 1))))
      }
    }
    
    ## marginalize rho by sample ##
    #rhoS <- list()
    #for (s in 1:S){
    #  rhoS[[s]] <- matrix(0, nrow = K, ncol = N)
    #  for (k in 1:K){
    #    rhoS[[s]][k, ] <- colSums(rho[param$jointCNstates[, s] == k, chrTrain])
    #  }
    #}
    
    ################ M-step ####################
    if (verbose) { message("runEM iter", i-1 , ": Maximization") }
    output <- estimateParamsMap(copy[chrTrain, ], n[, i - 1], sp[, i - 1], phi[, i - 1], 
                                lambdas[, , i - 1], piG[, i - 1], A, param,
                                rho[, chrTrain], Zcounts, 
                                estimateNormal = estimateNormal, estimatePloidy = estimatePloidy,
                                estimatePrecision = estimatePrecision, estimateTransition = estimateTransition,
                                estimateInitDist = estimateInitDist, estimateSubclone = estimateSubclone)
    if (verbose == TRUE) {
      for (s in 1:S){
        message("Sample", s, " n=", signif(output$n[s], digits = 4), ", sp=", signif(output$sp[s], digits = 4), 
                ", phi=", signif(output$phi[s], digits = 4), 
                ", lambda=", paste0(signif(output$lambda[, s], digits=5), collapse = ","),
                ", F=", signif(output$F, digits=4))
      }     
    }
    n[, i] <- output$n
    sp[, i] <- output$sp
    phi[, i] <- output$phi
    lambdas[, , i] <- output$lambda
    piG[, i] <- output$piG
    A <- output$A
    estF <- output$F
    
    # Recalculate the likelihood
    lambdasKS <- as.matrix(expand.grid(as.data.frame(lambdas[, , i])))
    mus[, , i] <- as.matrix(get2and3ComponentMixture(param$jointCNstates, param$jointSCstatus, n[, i], sp[, i], phi[, i]))
    for (ks in 1:KS) {
      probs <- tdistPDF(copy, mus[ks, , i], lambdasKS[ks, ], param$nu)
      py[ks, ] <- apply(probs, 1, prod) # multiply across samples for each data point to get joint likelihood.
    }
    
    prior <- priorProbs(n[, i], sp[, i], phi[, i], lambdas[, , i], piG[, i], A, param, 
                        estimateNormal = estimateNormal, estimatePloidy = estimatePloidy,
                        estimatePrecision = estimatePrecision, estimateTransition = estimateTransition,
                        estimateInitDist = estimateInitDist, estimateSubclone = estimateSubclone)
    
    
    # check converence 
    loglik[i] <- loglik[i] + prior$prior
    elapsedTime <- proc.time() - ptm
    if (verbose) { 
      message(paste("runEM iter", i-1, " Log likelihood:", loglik[i])) 
      message("runEM iter", i-1, " Time: ", format(elapsedTime[3] / 60, digits = 2), " min.")
    }
    if ((abs(loglik[i] - loglik[i - 1]) / abs(loglik[i])) < likChangeConvergence){#} && loglik[i] > loglik[i - 1]) {
      converged = 1
    }
    if (loglik[i] < loglik[i - 1]){
      #message("LIKELIHOOD DECREASED!!!")
      #message("Using previous iteration ", i-2)
      i <- i - 1
      converged = 1
    }
  }

 
  if (converged) {
    # Perform one last round of E-step to get latest responsibilities
    #E-step
    if (verbose) { message("runEM iter", i-1 ,": Re-calculating responsibilties from converged parameters.") }   
    for (j in 1:length(chrsI)) {
        if (length(I) > 0){   
            I <- intersect(chrsI[[j]], which(chrTrain))
            output <- .Call("forward_backward", piG[, i], A, py[, I], PACKAGE = "HMMcopy")
            rho[, I] <- output$rho
        }
    }
  }
  
  if (verbose) {
    totalTime <- proc.time() - ptmTotal
    message("runEM: Using optimal parameters from iter", i-1)
    message("runEM: Total elapsed time: ", format(totalTime[3] / 60, digits = 2), "min.")
  }
  
  #### Return parameters ####
  n <- n[, 1:i, drop=FALSE]
  sp <- sp[, 1:i, drop=FALSE]
  phi <- phi[, 1:i, drop=FALSE]
  mus <- mus[, , 1:i, drop=FALSE]
  lambdas <- lambdas[, , 1:i, drop=FALSE]
  piG <- piG[, 1:i, drop=FALSE]
  loglik = loglik[1:i]
  
  output <- vector('list', 0);
  output$n <- n
  output$sp <- sp
  output$phi <- phi
  output$mus <- mus
  output$lambdas <- lambdas
  output$pi <- piG
  output$A <- A
  output$loglik <- loglik
  output$rho <- rho
  output$param <- param
  output$py <- py
  output$iter <- i
  
  return(output)
}

get2and3ComponentMixture <- function(ct, ct.sc.status, n, sp, phi){
  S <- length(n)
  cn <- 2
  mu <- matrix(NA, nrow = nrow(ct), ncol = ncol(ct))
  for (s in 1:S){
    #subclonal 3 component
    mu[ct.sc.status[, s] == TRUE, s] <- (((1 - n[s]) * (1 - sp[s]) * ct[ct.sc.status[, s] == TRUE, s]) + 
                                           ((1 - n[s]) * sp[s] * cn) + (n[s] * cn)) / ((1 - n[s]) * phi[s] + n[s] * cn)
    #clonal 2 component
    mu[ct.sc.status[ ,s] == FALSE, s] <- ((1 - n[s]) * ct[ct.sc.status[, s] == FALSE, s] + n[s] * cn) / ((1 - n[s]) * phi[s] + n[s] * cn)
  }
  return(log(mu))
}

get2ComponentMixture <- function(ct, n, phi){
  if (length(phi) > 1 && length(phi) != length(n)){
    stop("get2ComponentMixture: length(n) not equal to length(phi)")
  }
  S <- length(n)
  #if (length(phi) == 1){
  #  phi <- rep(phi, length(n))
  #}
  cn <- 2
  #mu <-  ((1 - n) * ct + n * cn) / ((1 - n) * phi + n * cn)
  mu <- NULL
  for (s in 1:S){
    mu <- cbind(mu, ((1 - n[s]) * ct[, s] + n[s] * cn) / ((1 - n[s]) * phi[s] + n[s] * cn))
  }
  return(log(mu))
}


# Dirichlet probability density function, returns the probability of vector
# x under the Dirichlet distribution with parameter vector alpha
# Author: David Ross http://www.cs.toronto.edu/~dross/code/dirichletpdf.m
dirichletpdf <- function(x, alpha) {
  if (any(x < 0)) {
    return(0);
  }
  if (abs(sum(x) - 1) > 1e-3) {
    stop("Dirichlet PDF: probabilities do not sum to 1")
    return(0);
  }
  p <- exp(lgamma(sum(alpha)) - sum(lgamma(alpha))) * prod(x ^ (alpha - 1))
  return(p)
}

dirichletpdflog <- function(x, k) {
  c <- lgamma(sum(k, na.rm = TRUE)) - sum(lgamma(k), na.rm = TRUE)  #normalizing constant
  l <- sum((k - 1) * log(x), na.rm = TRUE)  #likelihood
  y <- c + l
  return(y)
}

gammapdflog <- function(x, a, b) { #rate and scale parameterization
  c <- a * log(b) - lgamma(a)  # normalizing constant  
  l <- (a - 1) * log(x) + (-b * x)  #likelihood  
  y <- c + l
  return(y)
}

betapdflog <- function(x, a, b) {
  y = -lbeta(a, b) + (a - 1) * log(x) + (b - 1) * log(1 - x)
  return(y)
}

tdistPDF <- function(x, mu, lambda, nu) {
  S <- ncol(x)
  if (!is.null(S)){
    p <- NULL
    for (s in 1:S){
      tpdf <- (gamma(nu / 2 + 0.5)/gamma(nu / 2)) * ((lambda[s] / (pi * nu)) ^ (0.5)) *
        (1 + (lambda[s] * (x[, s] - mu[s]) ^ 2) / nu) ^ (-0.5 * nu - 0.5)
      p <- cbind(p, tpdf)
    }
  }else{
    nu <- rep(nu, length(mu))
    p <- (gamma(nu / 2 + 0.5)/gamma(nu / 2)) * ((lambda / (pi * nu)) ^ (0.5)) *
      (1 + (lambda * (x - mu) ^ 2) / nu) ^ (-0.5 * nu - 0.5)
  }
  p[is.na(p)] <- 1
  return(p)
}


estimateParamsMap <- function(D, n_prev, sp_prev, phi_prev, lambda_prev, pi_prev, A_prev, 
                              params, rho, Zcounts, 
                              estimateNormal = TRUE, estimatePloidy = TRUE,
                              estimatePrecision = TRUE, estimateInitDist = TRUE, 
                              estimateTransition = TRUE, estimateSubclone = TRUE,
                              verbose = TRUE) {
  KS <- nrow(rho)
  K <- length(params$ct)
  T <- ncol(rho)
  S <- length(n_prev)
  #dimen <- ncol(D)
  intervalNormal <- c(1e-6, 1 - 1e-6)
  intervalSubclone <- c(1e-6, 1 - 1e-6)
  intervalPhi <- c(.Machine$double.eps, 10)
  intervalLambda <- c(1e-5, 1e4)
  
  # initialize params to be estimated
  n_hat <- n_prev
  sp_hat <- sp_prev
  phi_hat <- phi_prev
  lambda_hat <- lambda_prev
  pi_hat <- pi_prev
  A_hat <- A_prev
  
  # Update transition matrix A
  if (estimateTransition){
    for (k in 1:KS) {
      A_hat[k, ] <- Zcounts[k, ] + params$dirPrior[k, ]
      A_hat[k, ] <- normalize(A_hat[k, ])
    }
  }
  # map estimate for pi (initial state dist)
  if (estimateInitDist){
    pi_hat <- estimateMixWeightsParamMap(rho, params$kappa)
  }
  
  # map estimate for normal 
  if (estimateNormal){
    suppressWarnings(
      estNorm <- optim(n_prev, fn = completeLikelihoodFun, pType = rep("n", S), n = n_prev, sp = sp_prev, phi = phi_prev, 
                       lambda = lambda_prev, piG = pi_hat, A = A_hat,
                       params = params, D = D, rho = rho, Zcounts = Zcounts, 
                       estimateNormal = estimateNormal, estimatePloidy = estimatePloidy,
                       estimatePrecision = estimatePrecision, estimateInitDist = estimateInitDist,
                       estimateTransition = estimateTransition, estimateSubclone = estimateSubclone,
                       method = "L-BFGS-B",
                       lower = intervalNormal[1], upper = intervalNormal[2],
                       control = list(trace = 0, fnscale = -1))
    )
    n_hat <- estNorm$par
  }
  if (estimateSubclone){
    suppressWarnings(
      estSubclone <- optim(sp_prev, fn = completeLikelihoodFun, pType = rep("sp", S), n = n_hat, sp = sp_prev, 
                           phi = phi_prev, lambda = lambda_prev, piG = pi_hat, A = A_hat,
                           params = params, D = D, rho = rho, Zcounts = Zcounts, 
                           estimateNormal = estimateNormal, estimatePloidy = estimatePloidy,
                           estimatePrecision = estimatePrecision, estimateInitDist = estimateInitDist,
                           estimateTransition = estimateTransition, estimateSubclone = estimateSubclone,
                           method = "L-BFGS-B",
                           lower = intervalNormal[1], upper = intervalNormal[2],
                           control = list(trace = 0, fnscale = -1))
    )
    sp_hat <- estSubclone$par
  }
  if (estimatePloidy){
    suppressWarnings(
      estPhi <- optim(phi_prev, fn = completeLikelihoodFun, pType = rep("phi", length(phi_prev)), n = n_hat, sp = sp_hat, 
                      phi = phi_prev, lambda = lambda_prev, piG = pi_hat, A = A_hat,
                      params = params, D = D, rho = rho, Zcounts = Zcounts, 
                      estimateNormal = estimateNormal, estimatePloidy = estimatePloidy,
                      estimatePrecision = estimatePrecision, estimateInitDist = estimateInitDist,
                      estimateTransition = estimateTransition, estimateSubclone = estimateSubclone,
                      method = "L-BFGS-B",
                      lower = intervalPhi[1], upper = intervalPhi[2],
                      control = list(trace = 0, fnscale = -1))
    )
    phi_hat <- estPhi$par
  }
  if (estimatePrecision){
    suppressWarnings(
      estLambda <- optim(c(lambda_prev), fn = completeLikelihoodFun, pType = rep("lambda", K*S), n = n_hat, sp = sp_hat,
                         phi = phi_hat, lambda = lambda_prev, piG = pi_hat, A = A_hat,
                         params = params, D = D, rho = rho, Zcounts = Zcounts, 
                         estimateNormal = estimateNormal, estimatePloidy = estimatePloidy,
                         estimatePrecision = estimatePrecision, estimateInitDist = estimateInitDist,
                         estimateTransition = estimateTransition, estimateSubclone = estimateSubclone,
                         method = "L-BFGS-B", lower = intervalLambda[1],
                         control = list(trace = 0, fnscale = -1))
    )
    lambda_hat <- matrix(estLambda$par, ncol = S, byrow = FALSE)
  }
  
  #}
  #}
  
  #map estimate for phi (ploidy) 
  # if (estimatePloidy) {
  #    suppressWarnings(
  #   estPhi <- optim(phi_prev, fn = phiLikelihoodFun, n = n_prev, lambda = lambda_prev, 
  #  								params = params, D = D, rho = rhoS, method = "L-BFGS-B",
  #                 control = list(fnscale = -1, ndeps = 1e-5), lower = intervalPhi[1])#, upper = intervalPhi[2])
  #  )
  #  phi_hat <- estPhi$par
  #}
  
  # map estimate for lambda (Student's t precision)
  # if (estimatePrecision){
  #for (s in 1:S){
  #  paramsS <- params
  #  paramsS$betaLambda <- paramsS$betaLambda[s]
  #  for (k in 1:K){
  #		  suppressWarnings(
  #	  estLambda <- optim(lambda_prev[k, s], fn = lambdaLikelihoodFun, n = n_prev[s], phi = phi_prev, 
  #										params = paramsS, D = D[, s], rho = rhoS[[s]][k, , drop = FALSE], k = k,
  #										control = list(fnscale = -1), method = "L-BFGS-B",
  #										lower = intervalLambda[1])#, upper = intervalLambda[2])
  
  
  #lambda_hat <- matrix(estLambda$par, ncol = S, byrow = TRUE)
  #  }
  #}
  #}
  #lambda_hat <- estimatePrecisionParamMap(lambda_prev, n_prev, phi_prev, params, D, rho)
  
  
  #rm(a, b, c, d, e)
  gc(verbose = FALSE, reset = TRUE)
  #if (verbose == TRUE) {
  #  message("n=", format(n_hat, digits = 2), ", phi=", format(phi_hat, digits = 4), ", lambda=", paste0(format(lambda_hat, digits=2), collapse = " ", sep = ","))
  #}
  return(list(n = n_hat, sp = sp_hat, phi = phi_hat, lambda = lambda_hat, piG = pi_hat, A = A_hat, F = estLambda$value))
}

# Student's t likelihood function #
stLikelihood <- function(n, sp, phi, lambda, params, D, rho){
  KS <- nrow(rho)
  lik <- 0
  # Recalculate the likelihood
  lambdaKS <- as.matrix(expand.grid(as.data.frame(lambda)))
  mus <- as.matrix(get2and3ComponentMixture(params$jointCNstates, params$jointSCstatus, n, sp, phi))
  for (ks in 1:KS) {
    probs <- log(tdistPDF(D, mus[ks, ], lambdaKS[ks, ], params$nu))
    # multiply across samples for each data point to get joint likelihood.
    l <- rho[ks, ] %*% rowSums(as.data.frame(probs)) 
    lik <- lik + as.numeric(l)
  }
  return(lik)
}

## length of x will depend on which pararmeters are being estimated
## x values should be same as n, phi, lambda, pi, but may not necessarily
## include all of those variables
completeLikelihoodFun <- function(x, pType, n, sp, phi, lambda, piG, A, params, D, rho, Zcounts,
                                  estimateNormal = TRUE, estimatePloidy = TRUE,
                                  estimatePrecision = TRUE, estimateTransition = FALSE,
                                  estimateInitDist = TRUE, estimateSubclone = TRUE){
  KS <- nrow(rho)
  N <- ncol(rho)
  S <- params$numberSamples
  K <- length(params$ct)
  
  lambda <- matrix(lambda, ncol = S, byrow = FALSE)
  
  if (estimatePrecision && sum(pType == "lambda") > 0){
    lambda <- matrix(x[pType == "lambda"], ncol = S, byrow = FALSE)
  } 
  if (estimateNormal && sum(pType == "n") > 0){
    n <- x[pType == "n"]
  }  
  if (estimateSubclone && sum(pType == "sp") > 0){
    sp <- x[pType == "sp"]
  }  
  if (estimatePloidy && sum(pType == "phi") > 0){
    phi <- x[pType == "phi"]
  }
  
  ## prior probabilities ##
  prior <- priorProbs(n, sp, phi, lambda, piG, A, params, 
                      estimateNormal = estimateNormal, estimatePloidy = estimatePloidy,
                      estimatePrecision = estimatePrecision, estimateTransition = estimateTransition,
                      estimateInitDist = estimateInitDist, estimateSubclone = estimateSubclone)
  ## likelihood terms ##
  likObs <- stLikelihood(n, sp, phi, lambda, params, D, rho)
  likInit <- rho[, 1] %*% log(piG)
  likTxn <- 0
  for (ks in 1:KS){
    likTxn <- likTxn + Zcounts[ks, ] %*% log(A[ks, ])
  }
  
  ## sum together ##
  lik <- likObs + likInit + likTxn
  prior <- prior$prior
  f <- lik + prior
  return(f)
}

priorProbs <- function(n, sp, phi, lambda, piG, A, params, 
                       estimateNormal = TRUE, estimatePloidy = TRUE,
                       estimatePrecision = TRUE, estimateTransition = TRUE,
                       estimateInitDist = TRUE, estimateSubclone = TRUE){
  S <- params$numberSamples
  K <- length(params$ct)
  KS <- K ^ S
  ## prior terms ##
  priorA <- 0
  if (estimateTransition){
    for (ks in 1:KS){
      priorA <- priorA + dirichletpdflog(A[ks, ], params$dirPrior[ks, ])
    }
  }
  priorLambda <- 0
  if (estimatePrecision){
    for (s in 1:S){
      for (k in 1:K){
        priorLambda <- priorLambda + 
          gammapdflog(as.data.frame(lambda)[k, s], params$alphaLambda[k], params$betaLambda[k, s])
      }
    } 
  }
  priorLambda <- as.numeric(priorLambda)
  priorN <- 0
  if (estimateNormal){
    priorN <- sum(betapdflog(n, params$alphaN, params$betaN))
  }
  priorSP <- 0
  if (estimateNormal){
    priorSP <- sum(betapdflog(sp, params$alphaSp, params$betaSp))
  }
  priorPhi <- 0
  if (estimatePloidy){
    priorPhi <- sum(gammapdflog(phi, params$alphaPhi, params$betaPhi))
  }
  priorPi <- 0
  if (estimateInitDist){
    priorPi <- dirichletpdflog(piG, params$kappa)
  }
  prior <- priorA + priorLambda + priorN + priorSP + priorPhi + priorPi  
  return(list(prior = prior, priorA = priorA, priorLambda = priorLambda, 
              priorN = priorN, priorSP = priorSP, priorPhi = priorPhi, priorPi = priorPi))
}

estimatePrecisionParamMap <- function(lambda, n, phi, params, D, rho){
  mu <- get2ComponentMixture(params$ct, n, phi)
  mu_0 <- get2ComponentMixture(params$ct, params$n_0, params$phi_0)
  yr <- t(matrix(D, length(D), length(mu)))        # Vectorize parameter
  u <- (1 + params$nu) / (((yr - mu) ^ 2) * lambda + params$nu); # scale parameter
  
  # Calculate the precision
  lambda_hat <- (rowSums(rho, na.rm = TRUE) + params$alphaLambda + 1) /
    (rowSums(rho * u * ((yr - mu) ^ 2), na.rm = TRUE) +
       (params$eta * (mu - mu_0) ^ 2 + params$betaLambda))
  return(lambda_hat)
}

estimateMixWeightsParamMap <- function(rho, kappa) {
  K <- nrow(rho)
  #pi <- (rho[, 1] + kappa - 1) / (sum(rho[, 1]) + sum(kappa) - K)
  pi <- (rowSums(rho, na.rm = TRUE) + kappa - 1) / 
    #	(ncol(rho) + sum(kappa) - K)
    (sum(rowSums(rho, na.rm = TRUE)) + sum(kappa) - K)
  
  return(pi)
} 