% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gatk.R
\name{mutect_filter_gatk}
\alias{mutect_filter_gatk}
\title{Filter Mutect2 Calls wrapper for R}
\usage{
mutect_filter_gatk(
  sif_gatk = build_default_sif_list()$sif_gatk,
  bin_bcftools = build_default_tool_binary_list()$bin_bcftools,
  bin_bgzip = build_default_tool_binary_list()$bin_bgzip,
  bin_tabix = build_default_tool_binary_list()$bin_tabix,
  vcf = "",
  stats = NULL,
  contamination_table = NULL,
  extract_pass = TRUE,
  segmentation_table = NULL,
  orientation_model = NULL,
  output_name = "",
  ref_genome = build_default_reference_list()$HG19$reference$genome,
  output_dir = ".",
  verbose = FALSE,
  clean = FALSE,
  batch_config = build_default_preprocess_config(),
  threads = 4,
  ram = 4,
  mode = "local",
  executor_id = make_unique_id("filterMutect2Gatk"),
  task_name = "filterMutect2Gatk",
  time = "48:0:0",
  update_time = 60,
  wait = FALSE,
  hold = NULL
)
}
\arguments{
\item{sif_gatk}{\link{REQUIRED} Path to gatk sif file.}

\item{bin_bcftools}{\link{REQUIRED} Path to bcftools binary file.}

\item{bin_bgzip}{\link{REQUIRED} Path to bgzip binary file.}

\item{bin_tabix}{\link{REQUIRED} Path to tabix binary file.}

\item{vcf}{\link{REQUIRED} Path to VCF file.}

\item{stats}{\link{OPTIONAL} Path to stats information generated by Mutect2 for supplied VCF.}

\item{contamination_table}{\link{OPTIONAL} Path to contamination tables for each tumour samples in VCF.}

\item{extract_pass}{\link{OPTIONAL} Extract PASSing variants. Default TRUE.}

\item{segmentation_table}{\link{OPTIONAL} Path to segmentation tables for each tumour samples in VCF.}

\item{orientation_model}{\link{OPTIONAL} Path to orientation model generated F1R2 read information.}

\item{output_name}{\link{OPTIONAL} Name for the output. If not given the name of one of the samples will be used.}

\item{output_dir}{\link{OPTIONAL} Path to the output directory.}

\item{verbose}{\link{OPTIONAL} Enables progress messages. Default False.}

\item{clean}{\link{OPTIONAL} Remove unfiltered VCF after completion. Default FALSE.}

\item{threads}{\link{OPTIONAL} Number of threads to split the work. Default 4}

\item{ram}{\link{OPTIONAL} RAM memory to asing to each thread. Default 4}

\item{mode}{\link{OPTIONAL}  Where to parallelize. Default local. Options \link{"local","batch"}}

\item{executor_id}{\link{OPTIONAL} Task EXECUTOR ID. Default "recalCovariates"}

\item{task_name}{\link{OPTIONAL} Task name. Default "recalCovariates"}

\item{time}{\link{OPTIONAL} If batch mode. Max run time per job. Default "48:0:0"}

\item{update_time}{\link{OPTIONAL} If batch mode. Job update time in seconds. Default 60.}

\item{wait}{\link{OPTIONAL} If batch mode wait for batch to finish. Default FALSE}

\item{hold}{\link{OPTIONAL} HOld job until job is finished. Job ID.}
}
\description{
This function filters Mutect2 callsets under multiple conditions
}
\details{
For more information:

https://gatk.broadinstitute.org/hc/en-us/articles/360036856831-FilterMutectCalls
}
