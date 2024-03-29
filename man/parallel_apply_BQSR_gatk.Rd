% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gatk.R
\name{parallel_apply_BQSR_gatk}
\alias{parallel_apply_BQSR_gatk}
\title{Multiregion parallelization of apply_BQSR function}
\usage{
parallel_apply_BQSR_gatk(
  bin_samtools = build_default_tool_binary_list()$bin_samtools,
  sif_gatk = build_default_sif_list()$sif_gatk,
  bin_picard = build_default_tool_binary_list()$bin_picard,
  bam = "",
  regions = NULL,
  ref_genome = "",
  rec_table = "",
  clean = TRUE,
  output_dir = ".",
  verbose = FALSE,
  tmp_dir = ".",
  batch_config = build_default_preprocess_config(),
  mode = "local",
  executor_id = make_unique("par_applyBQSR"),
  task_name = "par_applyBQSR",
  time = "48:0:0",
  threads = 4,
  ram = 4,
  update_time = 60,
  wait = FALSE,
  hold = NULL
)
}
\arguments{
\item{bin_samtools}{\link{REQUIRED} Path to samtools executable. Default tools/samtools/samtools.}

\item{sif_gatk}{\link{REQUIRED} Path to gatk executable. Default tools/gatk/gatk.}

\item{bin_picard}{\link{REQUIRED} Path to picard executable. Default tools/picard/build/libs/picard.jar}

\item{bam}{\link{REQUIRED} Path to the BAM file.}

\item{regions}{\link{OPTIONAL} Regions to parallelize through.}

\item{ref_genome}{\link{REQUIRED} Path to reference genome}

\item{rec_table}{\link{REQUIRED} Path to the recalibratio table.}

\item{clean}{Clean intermediary files Default TRUE}

\item{output_dir}{\link{OPTIONAL} Path to the output directory.}

\item{verbose}{\link{OPTIONAL} Enables progress messages. Default False.}

\item{tmp_dir}{\link{OPTIONAL} Path to the temporary directory.}

\item{mode}{\link{REQUIRED} Where to parallelize. Default local. Options \link{"local","batch"}}

\item{executor_id}{\link{OPTIONAL} Task executor name. Default "par_applyBQSR"}

\item{task_name}{\link{OPTIONAL} Task name. Default "par_applyBQSR"}

\item{time}{\link{OPTIONAL} If batch mode. Max run time per job. Default "48:0:0"}

\item{threads}{\link{OPTIONAL} Number of threads for the main job. Default 4}

\item{ram}{\link{OPTIONAL} If batch mode. RAM memory in GB per job. Default 1}

\item{update_time}{\link{OPTIONAL} If batch mode. Show job updates every update time. Default 60}

\item{wait}{\link{OPTIONAL} If batch mode wait for batch to finish. Default FALSE}

\item{hold}{\link{OPTIONAL} HOld job until job is finished. Job ID.}
}
\description{
Recalibrates
Applies numerical corrections to each individual basecall based on the covariates analyzed before.
For more information about this function: https://gatk.broadinstitute.org/hc/en-us/articles/360050814312-ApplyBQSR
}
