% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gatk.R
\name{estimate_contamination_gatk}
\alias{estimate_contamination_gatk}
\title{Estimate sample contamination Gatk}
\usage{
estimate_contamination_gatk(
  sif_gatk = build_default_sif_list()$sif_gatk,
  rdata = NULL,
  selected = NULL,
  tumour = NA,
  normal = NA,
  tumour_pileup = "",
  normal_pileup = "",
  output_name = "",
  output_dir = ".",
  tmp_dir = ".",
  biallelic_db = build_default_reference_list()$HG19$variant$biallelic_reference,
  db_interval = build_default_reference_list()$HG19$variant$biallelic_reference,
  verbose = FALSE,
  batch_config = build_default_preprocess_config(),
  threads = 1,
  ram = 4,
  mode = "local",
  executor_id = make_unique_id("estimateContaminationGatk"),
  task_name = "estimateContaminationGatk",
  time = "48:0:0",
  update_time = 60,
  wait = FALSE,
  hold = NULL
)
}
\arguments{
\item{sif_gatk}{Path to gatk sif file.}

\item{output_name}{\link{OPTIONAL} Name for the output. If not given the name of one of the samples will be used.}

\item{output_dir}{Path to the output directory.}

\item{verbose}{\link{OPTIONAL} Enables progress messages. Default False.}

\item{threads}{\link{OPTIONAL} Number of threads to split the work. Default 4}

\item{ram}{\link{OPTIONAL} RAM memory to asing to each thread. Default 4}

\item{mode}{\link{REQUIRED} Where to parallelize. Default local. Options \link{"local","batch"}}

\item{executor_id}{Task EXECUTOR ID. Default "recalCovariates"}

\item{task_name}{Task name. Default "recalCovariates"}

\item{time}{\link{OPTIONAL} If batch mode. Max run time per job. Default "48:0:0"}

\item{update_time}{\link{OPTIONAL} If batch mode. Job update time in seconds. Default 60.}

\item{wait}{\link{OPTIONAL} If batch mode wait for batch to finish. Default FALSE}

\item{hold}{\link{OPTIONAL} HOld job until job is finished. Job ID.}

\item{f1r2}{Path to f1r2 files.}

\item{ref_genome}{Path to reference genome fasta file.}
}
\description{
This function estimates the sample cross-contamination
}
