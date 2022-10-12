# file:   utils.R
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

####################################
##### FUNCTION TO FILTER CHRS ######
####################################
# updated for GRanges #


keepChr <- function(tumour_reads, chrs = c(1:22,"X","Y")){	
	tumour_reads <- keepSeqlevels(tumour_reads, chrs, pruning.mode="coarse")
	sortSeqlevels(tumour_reads)
	return(sort(tumour_reads))
}

filterEmptyChr <- function(gr){
	require(plyr)
	ind <- daply(as.data.frame(gr), .variables = "seqnames", .fun = function(x){
	  rowInd <- apply(x[, 6:ncol(x), drop = FALSE], 1, function(y){
	    sum(is.na(y)) == length(y)
	  })
	  sum(rowInd) == nrow(x)
	})	
	return(keepSeqlevels(gr, value = names(which(!ind))))
}

####################################
##### FUNCTION GET SEQINFO ######
####################################
getSeqInfo <- function(genomeBuild = "hg19", genomeStyle = "NCBI"){
	bsg <- paste0("BSgenome.Hsapiens.UCSC.", genomeBuild)
	if (!require(bsg, character.only=TRUE, quietly=TRUE, warn.conflicts=FALSE)) {
		seqinfo <- Seqinfo(genome=genomeBuild)
	} else {
		seqinfo <- seqinfo(get(bsg))
	}
	seqlevelsStyle(seqinfo) <- genomeStyle
	seqinfo <- keepSeqlevels(seqinfo, value = chrs)
	#seqinfo <- cbind(seqnames = seqnames(seqinfo), as.data.frame(seqinfo))
	return(seqinfo)	
}

##################################################
##### FUNCTION TO FILTER CENTROMERE REGIONS ######
##################################################
excludeCentromere <- function(x, centromere, flankLength = 0, genomeStyle = "NCBI"){
	require(GenomeInfoDb)
	colnames(centromere)[1:3] <- c("seqnames","start","end")
	centromere$start <- centromere$start - flankLength
	centromere$end <- centromere$end + flankLength
	centromere <- as(centromere, "GRanges")
	seqlevelsStyle(centromere) <- genomeStyle
	centromere <- sort(centromere)	
	hits <- findOverlaps(query = x, subject = centromere)
	ind <- queryHits(hits)
	message("Removed ", length(ind), " bins near centromeres.")
	if (length(ind) > 0){
		x <- x[-ind, ]
	}
	return(x)
}

##################################################
##### FUNCTION TO USE NCBI CHROMOSOME NAMES ######
##################################################
## deprecated ##
setGenomeStyle <- function(x, genomeStyle = "NCBI", species = "Homo_sapiens"){
        require(GenomeInfoDb)
        #chrs <- genomeStyles(species)[c("NCBI","UCSC")]
        if (!genomeStyle %in% seqlevelsStyle(as.character(x))){
        x <- suppressWarnings(mapSeqlevels(as.character(x),
                                        genomeStyle, drop = FALSE)[1,])
    }

    autoSexMChr <- extractSeqlevelsByGroup(species = species,
                                style = genomeStyle, group = "all")
    x <- x[x %in% autoSexMChr]
    return(x)
}

wigToGRanges <- function(wigfile, verbose = TRUE){
  if (verbose) { message(paste("Slurping:", wigfile)) }
  input <- readLines(wigfile, warn = FALSE)
  breaks <- c(grep("fixedStep", input), length(input) + 1)
  temp <- NA
  span <- NA
  for (i in 1:(length(breaks) - 1)) {
    data_range <- (breaks[i] + 1):(breaks[i + 1] - 1)
    track_info <- input[breaks[i]]
    if (verbose) { message(paste("Parsing:", track_info)) }
    tokens <- strsplit(
      sub("fixedStep chrom=(\\S+) start=(\\d+) step=(\\d+) span=(\\d+)",
          "\\1 \\2 \\3 \\4", track_info, perl = TRUE), " ")[[1]]
    span <- as.integer(tokens[4])
    chr <- rep.int(tokens[1], length(data_range))
    pos <- seq(from = as.integer(tokens[2]), by = as.integer(tokens[3]),
               length.out = length(data_range))
    val <- as.numeric(input[data_range])
    temp <- c(temp, list(data.frame(chr, pos, val)))
  }
  if (verbose) { message("Sorting by decreasing chromosome size") }
  lengths <- as.integer(lapply(temp, nrow))
  temp <- temp[order(lengths, decreasing = TRUE)]
  temp = do.call("rbind", temp)
  output <- GenomicRanges::GRanges(ranges = IRanges(start = temp$pos, width = span),
                       seqnames= temp$chr, value = temp$val)
  return(output)
}


loadReadCountsFromWig <- function(counts, chrs = c(1:22, "X", "Y"), gc = NA, map = NA, centromere = NA, flankLength = 100000, targetedSequences = NA, genomeStyle = "NCBI", applyCorrection = TRUE, mapScoreThres = 0.9, chrNormalize = c(1:22, "X", "Y"), fracReadsInChrYForMale = 0.002, chrXMedianForMale = -0.5, useChrY = TRUE){
	require(HMMcopy)
	require(GenomeInfoDb)
	seqlevelsStyle(counts) <- genomeStyle
	counts.raw <- counts	
	counts <- keepChr(counts, chrs)
	
	if (!is.na(gc)){ 
		seqlevelsStyle(gc) <- genomeStyle
		counts$gc <- keepChr(gc, chrs)$value
	}
	if (!is.na(map)){ 
		seqlevelsStyle(map) <- genomeStyle
		counts$map <- keepChr(map, chrs)$value
	}
	colnames(values(counts))[1] <- c("reads")
	
	# remove centromeres
	if (!is.na(centromere)){ 
		counts <- excludeCentromere(counts, centromere, flankLength = flankLength, genomeStyle=genomeStyle)
	}
	# keep targeted sequences
	if (!is.na(targetedSequences)){
		colnames(targetedSequences)[1:3] <- c("chr", "start", "end")
		targetedSequences.GR <- as(targetedSequences, "GRanges")
		seqlevelsStyle(targetedSequences.GR) <- genomeStyle
		countsExons <- filterByTargetedSequences(counts, targetedSequences.GR)
		counts <- counts[countsExons$ix,]
	}
	gender <- NA
	if (applyCorrection){
		## correct read counts ##
		counts <- correctReadCounts(counts, chrNormalize = chrNormalize)
		if (!is.na(map)) {
		  ## filter bins by mappability
		  counts <- filterByMappabilityScore(counts, map=map, mapScoreThres = mapScoreThres)
		}
		## get gender ##
		gender <- getGender(counts.raw, counts, gc, map, fracReadsInChrYForMale = fracReadsInChrYForMale, 
							chrXMedianForMale = chrXMedianForMale, useChrY = useChrY,
							centromere=centromere, flankLength=flankLength, targetedSequences = targetedSequences,
							genomeStyle = genomeStyle)
    }
  return(list(counts = counts, gender = gender))
}

filterByMappabilityScore <- function(counts, map, mapScoreThres = 0.9){
	message("Filtering low uniqueness regions with mappability score < ", mapScoreThres)
	counts <- counts[counts$map >= mapScoreThres, ]
	return(counts)
}

filterByTargetedSequences <- function(counts, targetedSequences){
 ### for targeted sequencing (e.g.  exome capture),
    ### ignore bins with 0 for both tumour and normal
    ### targetedSequence = GRanges object
    ### containing list of targeted regions to consider;
    ### 3 columns: chr, start, end
					
	hits <- findOverlaps(query = counts, subject = targetedSequences)
	keepInd <- unique(queryHits(hits))    

	return(list(counts=counts, ix=keepInd))
}

selectFemaleChrXSolution <- function(){
	
}

##################################################
### FUNCTION TO DETERMINE GENDER #################
##################################################
getGender <- function(rawReads, normReads, gc, map, fracReadsInChrYForMale = 0.002, chrXMedianForMale = -0.5, useChrY = TRUE,
					  centromere=NA, flankLength=1e5, targetedSequences=NA, genomeStyle="NCBI"){
	chrXStr <- grep("X", runValue(GenomeInfoDb::seqnames(normReads)), value = TRUE)
	chrYStr <- grep("Y", runValue(GenomeInfoDb::seqnames(rawReads)), value = TRUE)
	chrXInd <- as.character(GenomeInfoDb::seqnames(normReads)) == chrXStr
	if (sum(chrXInd) > 1){ ## if no X 
		chrXMedian <- median(normReads[chrXInd, ]$copy, na.rm = TRUE)
		# proportion of reads in chrY #
		tumY <- loadReadCountsFromWig(rawReads, chrs=chrYStr, genomeStyle=genomeStyle,
				gc=gc, map=map, applyCorrection = FALSE, centromere=centromere, flankLength=flankLength, 
				targetedSequences=targetedSequences)$counts
		chrYCov <- sum(tumY$reads) / sum(rawReads$value)
		if (chrXMedian < chrXMedianForMale){
			if (useChrY && (chrYCov < fracReadsInChrYForMale)){ #trumps chrX if using chrY
					gender <- "female"  
			}else{
				gender <- "male" # satisfies decreased chrX log ratio and/or increased chrY coverage
			}
		}else{
			gender <- "female" # chrX is provided but does not satisfies male critera
		}
	}else{
		gender <- "unknown" # chrX is not provided
		chrYCov <- NA
		chrXMedian <- NA
	}
	return(list(gender=gender, chrYCovRatio=chrYCov, chrXMedian=chrXMedian))
}
	
	
normalizeByPanelOrMatchedNormal <- function(tumour_copy, chrs = c(1:22, "X", "Y"), 
      normal_panel = NA, normal_copy = NA, gender = "female", normalizeMaleX = FALSE){
    genomeStyle <- seqlevelsStyle(tumour_copy)[1]
    seqlevelsStyle(chrs) <- genomeStyle
 	### COMPUTE LOG RATIO FROM MATCHED NORMAL OR PANEL AND HANDLE CHRX ###
	## NO PANEL
	# matched normal but NO panel, then just normalize by matched normal (WES)
	## WHY DO WE NOT NORMALIZE BY NORMAL WITH PANEL? ##
	chrXInd <- grep("X", as.character(GenomeInfoDb::seqnames(tumour_copy)))
	chrXMedian <- median(tumour_copy[chrXInd, ]$copy, na.rm = TRUE)
	if (!is.na(normal_copy) && is.na(normal_panel)){
			message("Normalizing Tumour by Normal")
			tumour_copy$copy <- tumour_copy$copy - normal_copy$copy
			rm(normal_copy)
	}
	# matched normal and panel and male, then compute normalized chrX median (WES)
	if (!is.na(normal_copy) && !is.na(normal_panel) && gender=="male"){
			message("Normalizing by matched normal for ChrX")
			chrX.MNnorm <- tumour_copy$copy[chrXInd] - normal_copy$copy[chrXInd]
			chrXMedian.MNnorm <- median(chrX.MNnorm, na.rm = TRUE)
	}
	# if male, then just normalize chrX to median (ULP and WES)
	if (is.na(normal_copy) && gender=="male" && !gender.mismatch && normalizeMaleX){
			tumour_copy$copy[chrXInd] <- tumour_copy$copy[chrXInd] - chrXMedian
	}
	# PANEL, then normalize by panel instead of matched normal (ULP and WES)
	if (!is.na(normal_panel)){
		## load in IRanges object, then convert to GRanges
		panel <- readRDS(normal_panel)
		seqlevelsStyle(panel) <- genomeStyle
		panel <- keepChr(panel, chr = chrs)
        # intersect bins in sample and panel
        hits <- findOverlaps(tumour_copy, panel, type="equal")
        tumour_copy <- tumour_copy[queryHits(hits),]
        panel <- panel[subjectHits(hits),]
        # subtract out panel median
		tumour_copy$copy <- tumour_copy$copy - panel$Median
		# if male, then shift chrX by +chrXMedian.MNnorm
		if (gender == "male" && exists("chrXMedian.MNnorm")){
			tumour_copy$copy[chrXInd] <- tumour_copy$copy[chrXInd] + chrXMedian.MNnorm
		}
	}
	return(tumour_copy)
}

##################################################
###### FUNCTION TO CORRECT GC/MAP BIASES ########
##################################################
correctReadCounts <- function(x, chrNormalize = c(1:22), mappability = 0.9, samplesize = 50000,
    verbose = TRUE) {
  if (length(x$reads) == 0 | length(x$gc) == 0) {
    stop("Missing one of required columns: reads, gc")
  }
  chrInd <- as.character(GenomeInfoDb::seqnames(x)) %in% chrNormalize
  if(verbose) { message("Applying filter on data...") }
  x$valid <- TRUE
  x$valid[x$reads <= 0 | x$gc < 0] <- FALSE
  x$ideal <- TRUE
  routlier <- 0.01
  range <- quantile(x$reads[x$valid & chrInd], prob = c(0, 1 - routlier), na.rm = TRUE)
  doutlier <- 0.001
  domain <- quantile(x$gc[x$valid & chrInd], prob = c(doutlier, 1 - doutlier), na.rm = TRUE)
  if (length(x$map) != 0) {
    x$ideal[!x$valid | x$map < mappability | x$reads <= range[1] |
      x$reads > range[2] | x$gc < domain[1] | x$gc > domain[2]] <- FALSE
  } else {
    x$ideal[!x$valid | x$reads <= range[1] |
      x$reads > range[2] | x$gc < domain[1] | x$gc > domain[2]] <- FALSE
  }

  if (verbose) { message("Correcting for GC bias...") }
  set <- which(x$ideal & chrInd)
  select <- sample(set, min(length(set), samplesize))
  rough = loess(x$reads[select] ~ x$gc[select], span = 0.03)
  i <- seq(0, 1, by = 0.001)
  final = loess(predict(rough, i) ~ i, span = 0.3)
  x$cor.gc <- x$reads / predict(final, x$gc)

  if (length(x$map) != 0) {
    if (verbose) { message("Correcting for mappability bias...") }
    coutlier <- 0.01
    range <- quantile(x$cor.gc[which(x$valid & chrInd)], prob = c(0, 1 - coutlier), na.rm = TRUE)
    set <- which(x$cor.gc < range[2] & chrInd)
    select <- sample(set, min(length(set), samplesize))
    final = approxfun(lowess(x$map[select], x$cor.gc[select]))
    x$cor.map <- x$cor.gc / final(x$map)
  } else {
    x$cor.map <- x$cor.gc
  }
  x$copy <- x$cor.map
  x$copy[x$copy <= 0] = NA
  x$copy <- log(x$copy, 2)
  return(x)
}

## Recompute integer CN for high-level amplifications ##
## compute logR-corrected copy number ##
correctIntegerCN <- function(cn, segs, callColName = "event", 
		purity, ploidy, cellPrev, maxCNtoCorrect.autosomes =NA, 
		maxCNtoCorrect.X = NA, correctHOMD = TRUE,
		 minPurityToCorrect = 0.2, gender = "male", chrs = c(1:22, "X")){
	names <- c("HOMD","HETD","NEUT","GAIN","AMP","HLAMP", rep("HLAMP", 1000))

	## set up chromosome style
	autosomeStr <- grep("X|Y", chrs, value=TRUE, invert=TRUE)
	chrXStr <- grep("X", chrs, value=TRUE)
	
	if (is.na(maxCNtoCorrect.autosomes)){
		maxCNtoCorrect.autosomes <- max(segs[segs$chr %in% autosomeStr, "copy.number"], na.rm = TRUE)
	}
	if (is.na(maxCNtoCorrect.X) & gender == "female" & length(chrXStr) > 0){
		maxCNtoCorrect.X <- max(segs[segs$chr == chrXStr, "copy.number"], na.rm=TRUE)
	}
	## correct log ratio and compute corrected CN
	cellPrev.seg <- rep(1, nrow(segs))
	cellPrev.seg[as.logical(segs$subclone.status)] <- cellPrev
	segs$logR_Copy_Number <- logRbasedCN(segs[["median"]], purity, ploidy, cellPrev.seg, cn=2)
	if (gender == "male" & length(chrXStr) > 0){ ## analyze chrX separately
		ind.cnChrX <- which(segs$chr == chrXStr)
		segs$logR_Copy_Number[ind.cnChrX] <- logRbasedCN(segs[["median"]][ind.cnChrX], purity, ploidy, cellPrev.seg[ind.cnChrX], cn=1)
	}

	## assign copy number to use - Corrected_Copy_Number
	# 1) ame ichorCNA calls for autosomes - no change in copy number
	segs$Corrected_Copy_Number <- as.integer(segs$copy.number)
	segs$Corrected_Call <- segs[[callColName]]

	ind.change <- c()
	if (purity >= minPurityToCorrect){
		# 2) ichorCNA calls adjusted for >= copies - HLAMP
		# perform on all chromosomes
		ind.cn <- which(segs$copy.number >= maxCNtoCorrect.autosomes | 
						(segs$logR_Copy_Number >= maxCNtoCorrect.autosomes * 1.2 & !is.infinite(segs$logR_Copy_Number)))
		segs$Corrected_Copy_Number[ind.cn] <- as.integer(round(segs$logR_Copy_Number[ind.cn]))
		segs$Corrected_Call[ind.cn] <- names[segs$Corrected_Copy_Number[ind.cn] + 1]
		ind.change <- c(ind.change, ind.cn)
		
		# 3) ichorCNA calls adjust for HOMD
		if (correctHOMD){
			ind.cn <- which(segs$chr %in% chrs & 
				(segs$copy.number == 0 | segs$logR_Copy_Number == 1/2^6))
			segs$Corrected_Copy_Number[ind.cn] <- as.integer(round(segs$logR_Copy_Number[ind.cn]))
			segs$Corrected_Call[ind.cn] <- names[segs$Corrected_Copy_Number[ind.cn] + 1]
			ind.change <- c(ind.change, ind.cn)
		}
		# 4) Re-adjust chrX copy number for males (females already handled above)
		if (gender == "male" & length(chrXStr) > 0){
			ind.cn <- which(segs$chr == chrXStr & 
				(segs$copy.number >= maxCNtoCorrect.X | segs$logR_Copy_Number >= maxCNtoCorrect.X * 1.2))
			segs$Corrected_Copy_Number[ind.cn] <- as.integer(round(segs$logR_Copy_Number[ind.cn]))
			segs$Corrected_Call[ind.cn] <- names[segs$Corrected_Copy_Number[ind.cn] + 2]
			ind.change <- c(ind.change, ind.cn)
		}
	}
	## adjust the bin level data ##
	# 1) assign the original calls
	cn$Corrected_Copy_Number <- as.integer(cn$copy.number)
	cn$Corrected_Call <- cn[[callColName]]
	cellPrev.cn <- rep(1, nrow(cn))
	cellPrev.cn[as.logical(cn$subclone.status)] <- cellPrev
	cn$logR_Copy_Number <- logRbasedCN(cn[["logR"]], purity, ploidy, cellPrev.cn, cn=2)
	if (gender == "male" & length(chrXStr) > 0){ ## analyze chrX separately
		ind.cnChrX <- which(cn$chr == chrXStr)
		cn$logR_Copy_Number[ind.cnChrX] <- logRbasedCN(cn[["logR"]][ind.cnChrX], purity, ploidy, cellPrev.cn[ind.cnChrX], cn=1)
	}
	if (purity >= minPurityToCorrect){
		# 2) correct bins overlapping adjusted segs
		ind.change <- unique(ind.change)
		ind.overlapSegs <- c()
		if (length(ind.change) > 0){		
			cn.gr <- as(cn, "GRanges")
			segs.gr <- as(segs, "GRanges")
			hits <- IRanges::findOverlaps(query = cn.gr, subject = segs.gr[ind.change])
			cn$Corrected_Copy_Number[S4Vectors::queryHits(hits)] <- segs$Corrected_Copy_Number[ind.change][S4Vectors::subjectHits(hits)]
			cn$Corrected_Call[S4Vectors::queryHits(hits)] <- segs$Corrected_Call[ind.change][S4Vectors::subjectHits(hits)]
			ind.overlapSegs <- S4Vectors::queryHits(hits)
		}
		# 3) correct bins that are missed as high level amplifications
		ind.hlamp <- which(cn$copy.number >= maxCNtoCorrect.autosomes | 
	 					(cn$logR_Copy_Number >= maxCNtoCorrect.autosomes * 1.2 & !is.infinite(cn$logR_Copy_Number)))
		ind.cn <- unique(ind.hlamp, ind.overlapSegs)
	 	cn$Corrected_Copy_Number[ind.cn] <- as.integer(round(cn$logR_Copy_Number[ind.cn]))
	 	cn$Corrected_Call[ind.cn] <- names[cn$Corrected_Copy_Number[ind.cn] + 1]
	 }
	 
	return(list(cn = cn, segs = segs))
}


## compute copy number using corrected log ratio ##
logRbasedCN <- function(x, purity, ploidyT, cellPrev=NA, cn = 2){
	if (length(cellPrev) == 1 && is.na(cellPrev)){
		cellPrev <- 1
	}else{ #if cellPrev is a vector
		cellPrev[is.na(cellPrev)] <- 1
	}
	ct <- (2^x 
		* (cn * (1 - purity) + purity * ploidyT * (cn / 2)) 
		- (cn * (1 - purity)) 
		- (cn * purity * (1 - cellPrev))) 
	ct <- ct / (purity * cellPrev)
	ct <- sapply(ct, max, 1/2^6)
	return(ct)
}


computeBIC <- function(params){
  iter <- params$iter
  KS <- nrow(params$rho) # num states
  N <- ncol(params$rho) # num data points
  NP <- nrow(params$n) + nrow(params$phi) # normal + ploidy
  L <- 1 # precision (lambda)
  numParams <- KS * (KS - 1) + KS * (L + NP) - 1
  b <- -2 * params$loglik[iter] + numParams * log(N)
  return(b)
} 


# file:   output.R
# author: Gavin Ha, Ph.D.
#		Fred Hutchinson Cancer Research Center
# contact: <gha@fredhutch.org>
# author: Justin Rhoades
#       Dana-Farber Cancer Institute
#       Broad Institute
# https://github.com/broadinstitute/ichorCNA
# date:  July 12, 2019
# description: Hidden Markov model (HMM) to analyze Ultra-low pass whole genome sequencing (ULP-WGS) data.

##################################################
###### FUNCTION TO GET OUTPUT HMM RESULTS ########
##################################################






outputHMM <- function(cna, segs, results, patientID = NA, outDir = "."){
  names <- c("HOMD","HETD","NEUT","GAIN","AMP","HLAMP",paste0(rep("HLAMP", 8), 2:25))

  S <- results$param$numberSamples
  
  segout <- NA
  shuffle <- NA
  if (is.na(patientID)){
    patientID <- names(cna)[1]
  }
  ### PRINT OUT SEGMENTS FOR EACH SAMPLE TO SEPARATE FILES ###
  for (s in 1:S){
    id <- names(cna)[s]
    bin_size <- as.numeric(cna[[s]][1,"end"]) - as.numeric(cna[[s]][1,"start"]) + 1
    markers <- (segs[[s]]$end - segs[[s]]$start + 1) / bin_size  
    segTmp <- cbind(sample = as.character(id), segs[[s]][, 1:3],
                    event = names[segs[[s]]$copy.number + 1],
                    copy.number = segs[[s]]$copy.number, 
                    bins = markers, median = segs[[s]]$median,
                    subclone.status=segs[[s]]$subclone.status,
                    segs[[s]][, 9:ncol(segs[[s]])])
    segout <- rbind(segout, segTmp[, 1:8])
    ### Re-ordering the columns for output ###
    shuffleTmp <- segTmp[, c(1, 2, 3, 4, 7, 8, 6, 5, 9:ncol(segTmp))]
    colnames(shuffleTmp)[1:9] <- c("ID", "chrom", "start", "end", 
                              "num.mark", "seg.median.logR", "copy.number", "call", "subclone.status")
    shuffle <- rbind(shuffle, shuffleTmp) 
  }
  
  shu_out = paste(outDir,"/", patientID,".seg.txt",sep="")
  seg_out = paste(outDir,"/", patientID,".seg",sep="")
  message("Writing segments to ", seg_out)
  write.table(shuffle, file = shu_out, quote = FALSE, sep = "\t", row.names = FALSE)
  write.table(segout, file = seg_out, quote = FALSE, sep = "\t", row.names = FALSE)
  ### PRINT OUT BIN-LEVEL DATA FOR ALL SAMPLES IN SAME FILE ##
  cnaout <- cbind(cna[[1]][, -c(1)])
  colnames(cnaout)[4:ncol(cnaout)] <- paste0(names(cna)[1], ".", colnames(cnaout)[4:ncol(cnaout)])
  if (S >= 2) {
    for (s in 2:S){
      id <- names(cna)[s]
      cnaTmp <- cna[[s]][, -c(1)]
      colnames(cnaTmp)[4:ncol(cnaTmp)] <- paste0(id, ".", colnames(cnaTmp)[4:ncol(cnaTmp)])
      cnaout <- merge(cnaout, cnaTmp, by = c("chr", "start", "end"))
    }
  }
  cnaout$chr <- factor(cnaout$chr, levels = unique(cna[[1]]$chr))
  cnaout <- cnaout[order(cnaout[, 1], cnaout[, 2], cnaout[, 3]), ]
  cna_out = paste(outDir,"/", patientID,".cna.seg",sep="")
  message("Outputting to bin-level results to ", cna_out)
  write.table(cnaout, file = cna_out, quote = FALSE, sep = "\t", row.names = FALSE)
}

outputParametersToFile <- function(hmmResults, file){
  S <- hmmResults$results$param$numberSamples
  x <- hmmResults$results
  i <- x$iter
  fc <- file(file, "w+")
  outMat <- as.data.frame(cbind(`Tumor Fraction` = signif(1 - x$n[, i], digits = 4), Ploidy = signif(x$phi[, i], digits = 4)))
  outMat <- cbind(Sample = names(hmmResults$cna), outMat)
  rownames(outMat) <- names(hmmResults$cna)
  write.table(outMat, file = fc, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  write.table("\n", file = fc, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
  for (s in 1:S){
    id <- names(hmmResults$cna)[s]
    write.table(id, file = fc, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
    write.table(paste0("Gender:\t", x$gender), file = fc, col.names = FALSE, 
                row.names = FALSE, quote = FALSE, sep = "\t")
    write.table(paste0("Tumor Fraction:\t", signif(1 - x$n[s, i], digits = 4)), file = fc, col.names = FALSE, 
                row.names = FALSE, quote = FALSE, sep = "\t")
    ploidy <- (1 - x$n[i]) * x$phi[i] + x$n[i] * 2
    write.table(paste0("Ploidy:\t", signif(x$phi[s, i], digits = 4)), file = fc, col.names = FALSE, 
                row.names = FALSE, quote = FALSE, sep = "\t")
   	subcloneGenomeFrac <- sum(hmmResults$cna[[s]]$subclone.status) / nrow(hmmResults$cna[[s]])
   	subcloneCNAFrac <- sum(hmmResults$cna[[s]]$subclone.status) / sum(hmmResults$cna[[s]]$copy.num != 2)
   	scFrac <- 1 - x$sp[s, i]
   	if (subcloneGenomeFrac == 0){
   		scFrac <- NA
   	}else{
   		scFrac <- signif(scFrac, digits = 4)
   	}
    #if (sum(x$param$ct.sc.status) != 0){    	
      write.table(paste0("Subclone Fraction:\t", scFrac), file = fc, col.names = FALSE, 
                  row.names = FALSE, quote = FALSE, sep = "\t")
      write.table(paste0("Fraction Genome Subclonal:\t", signif(subcloneGenomeFrac, digits = 4)), file = fc, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
      write.table(paste0("Fraction CNA Subclonal:\t", signif(subcloneCNAFrac, digits = 4)), file = fc, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
    #}
    if (is.na(x$coverage)){
    	coverage <- NA
    }else{
    	coverage <- signif(coverage, digits = 4)
    }
		write.table(paste0("Coverage:\t", coverage), file = fc, col.names = FALSE, 
								row.names = FALSE, quote = FALSE, sep = "\t")
	
    if (!is.na(x$chrYCov)){
      write.table(paste0("ChrY coverage fraction:\t", signif(x$chrYCov[s], digits = 4)), file = fc, col.names = FALSE, 
                  row.names = FALSE, quote = FALSE, sep = "\t")
    }
    if (!is.na(x$chrXMedian)){
      write.table(paste0("ChrX median log ratio:\t", signif(x$chrXMedian[s], digits = 4)), file = fc, col.names = FALSE, 
                  row.names = FALSE, quote = FALSE, sep = "\t")
    }
    write.table(paste0("Student's t mean: ", paste0(signif(x$mus[,s,i], digits = 2), collapse = ", ")), 
                file = fc, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
    write.table(paste0("Student's t precision: ", paste0(signif(x$lambdas[,s,i], digits = 2), collapse = ", ")), file = fc, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
    write.table(paste0("Gamma Rate Init:\t", signif(hmmResults$results$param$betaLambda[1], digits=2)), file=fc, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
    write.table(paste0("GC-Map correction MAD:\t", format(mad(diff(2^as.numeric(hmmResults$cna[[s]][,"logR"])), na.rm=T), digits=4)), file = fc, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
    write.table("\n", file = fc, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
  }
  write.table(x$loglik, file = fc, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  close(fc)
  invisible()
}



