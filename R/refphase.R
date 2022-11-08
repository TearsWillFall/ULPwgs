
### Modified version of refphase_load_ascat to fix issue with ASCAT naming

modified_refphase_load_ascat <- function(samples, ascat_input, ascat_output, het_only = FALSE, homozygous_cutoff = 0.7) {

  validate_ascat_data(samples, ascat_input, ascat_output)
  data <- list("segs" = list(), "snps" = list(), "purity" = list(), "ploidy" = list())

    #' Format chromosomes to fit numerical
    #'
    #' Formatting such that chromosom 1 is "1". Chromosome X and Y are labelled as "X" and "Y", respectively
    format_chromosomes <- function(chromosomes) {
    chromosomes <- gsub("(23)|x", "X", chromosomes)
    chromosomes <- gsub("(24)|y", "Y", chromosomes)
    chromosomes <- gsub("(c|C)hr(om)?([0-9XYxy]{1,2})", "\\3", chromosomes)

    stopifnot(gsub("[0-9]{1,2}|[XYxy]", "", chromosomes) == "")

    chromosomes
    }

    format_snps <- function(snps, samples, het_only = FALSE) {
    for (sample in samples) {
        # Mask non-hz SNPs with NA
        GenomicRanges::mcols(snps[[sample]])$baf[GenomicRanges::mcols(snps[[sample]])$germline_zygosity != "het"] <- NA
        # Only keep heterozygous positions
        if (isTRUE(het_only)) {
        snps[[sample]] <- snps[[sample]][GenomicRanges::mcols(snps[[sample]])$germline_zygosity == "het"]
        }
    }

    snps <- GenomicRanges::GRangesList(snps, compress = FALSE)
    snps <- intersect_all_snps(snps)

    names(snps) <- samples

    snps
    }

    format_segs <- function(segs, samples) {
    if (any(segs$start > segs$end)) {
        stop("Error: Start position in segmentation must be <= end position")
    }

    segs <- GenomicRanges::GRangesList(segs, compress = FALSE)
    names(segs) <- samples

    segs
    }

    intersect_all_snps <- function(snps) {

    intersected_snps <- IRanges::subsetByOverlaps(snps[[1]], snps[[2]])
    if (length(snps) >= 3) {
        for (sample in names(snps)[3:length(snps)]) {
        intersected_snps <- IRanges::subsetByOverlaps(intersected_snps, snps[[sample]])
        }
    }
    for (sample in names(snps)) {
        snps[[sample]] <- IRanges::subsetByOverlaps(snps[[sample]], intersected_snps)
    }

    snps
    }
    sample_data_from_ploidy_purity <- function(samples, ploidy, purity) {
        sample_data <- data.frame("sample_id" = samples,
                                    "ploidy" = unname(unlist(ploidy[samples])),
                                    "purity" = unname(unlist(purity[samples])))
        sample_data$segmentation <- ""
        sample_data$snps <- ""

        sample_data
    }



  for (sample in samples) {
    cur_ascat_input <- ascat_input[[sample]]
    cur_ascat_output <- ascat_output[[sample]]

    # 1) SNPs
    cur_chroms <- format_chromosomes(cur_ascat_input$SNPpos$Chromosome)
    cur_snps <- GPos(
      seqnames = Rle(cur_chroms),
      pos = cur_ascat_input$SNPpos$Position,
      baf = cur_ascat_input$Tumor_BAF[[1]],
      logr = cur_ascat_input$Tumor_LogR[[1]],
      germline_zygosity = ifelse(is.na(cur_ascat_input$Germline_BAF[[1]]) |
        abs(0.5 - cur_ascat_input$Germline_BAF[[1]]) >= (homozygous_cutoff / 2), "hom", "het")
    )

    data$snps[[sample]] <- cur_snps

    # 2) Segments
    cur_segs <- cur_ascat_output$segments
    names(cur_segs) <- c("sample_id", "chrom", "start", "end", "cn_major", "cn_minor")
    cur_segs$chrom <- format_chromosomes(cur_segs$chrom)
    cur_segs$sample_id <- sample

    cur_segs <- GenomicRanges::makeGRangesFromDataFrame(cur_segs,
      keep.extra.columns = TRUE,
      ignore.strand = TRUE,
      seqnames.field = "chrom",
      start.field = "start",
      end.field = "end",
      starts.in.df.are.0based = FALSE
    )

    data$segs[[sample]] <- cur_segs

    # 3) Purity and Ploidy
    data$purity[[sample]] <- cur_ascat_output$aberrantcellfraction[[1]]
    data$ploidy[[sample]] <- cur_ascat_output$ploidy[[1]]
  }

  
  data$segs <- format_segs(data$segs, samples)
  data$snps <- format_snps(data$snps, samples, het_only = het_only)

  data$sample_data <- sample_data_from_ploidy_purity(samples, data$ploidy, data$purity)
  class(data) <- "RefphaseInput"

  data
}


