% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gatk.R
\name{recal_gatk}
\alias{recal_gatk}
\title{Function for base quality recalibration}
\usage{
recal_gatk(
  bin_samtools = build_default_tool_binary_list()$bin_samtools,
  sif_gatk = build_default_sif_list()$sif_gatk,
  bin_picard = build_default_tool_binary_list()$bin_picard,
  bam = "",
  ref_genome = "",
  dbsnp = "",
  ram = 4,
  threads = 4,
  output_dir = ".",
  tmp_dir = ".",
  verbose = FALSE,
  batch_config = build_default_preprocess_config(),
  executor_id = make_unique_id("recalGATK"),
  task_name = "recalGATK",
  clean = TRUE,
  mode = "local",
  time = "48:0:0",
  update_time = 60,
  wait = FALSE,
  hold = NULL
)
}
\arguments{
\item{bin_samtools}{\link{REQUIRED} Path to picard executable. Default path tools/samtools/samtools.}

\item{sif_gatk}{\link{REQUIRED} Path to picard executable. Default path tools/gatk/gatk.}

\item{bin_picard}{\link{REQUIRED} Path to picard executable. Default path tools/picard/build/libs/picard.jar.}

\item{bam}{\link{REQUIRED}  Path to BAM file.}

\item{ref_genome}{\link{REQUIRED}  Path to reference genome.}

\item{dbsnp}{\link{REQUIRED} Known variant database.Requires atleast 1.}

\item{ram}{\link{OPTIONAL} RAM memory per thread.}

\item{threads}{\link{OPTIONAL} Number of threads to split the work.}

\item{output_dir}{\link{OPTIONAL} Path to the output directory.}

\item{tmp_dir}{\link{OPTIONAL} Path to the temporary directory.}

\item{verbose}{\link{OPTIONAL} Enables progress messages. Default False.#}

\item{executor_id}{Task EXECUTOR ID. Default "recalGATK"}

\item{clean}{Clean input files. Default TRUE.}

\item{mode}{\link{REQUIRED} Where to parallelize. Default local. Options \link{"local","batch"}}

\item{time}{\link{OPTIONAL} If batch mode. Max run time per job. Default "48:0:0"}

\item{update_time}{\link{OPTIONAL} If batch mode. Job update time in seconds. Default 60.}

\item{wait}{\link{OPTIONAL} If batch mode wait for batch to finish. Default FALSE}

\item{hold}{\link{OPTIONAL} HOld job until job is finished. Job ID.}

\item{task_id}{Task nam. Default "recalGATK"}
}
\description{
This function recalibrates the base quality of the reads in two steps process based on GATK best practices guides.
}
