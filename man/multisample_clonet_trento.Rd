% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trento.R
\name{multisample_clonet_trento}
\alias{multisample_clonet_trento}
\title{Process multiple tumour samples using CLONET}
\usage{
multisample_clonet_trento(
  sample_sheet = NA,
  bam_dir = "",
  normal_id = "",
  patient_id = "",
  version = "V3",
  tc = NULL,
  ploidy = NULL,
  tmp_dir = ".",
  header = TRUE,
  sep = "",
  threads = 3,
  ram = 4,
  output_dir = ".",
  verbose = FALSE,
  batch_config = build_default_preprocess_config(),
  executor_id = make_unique_id("multi_clonet"),
  task_name = "multi_clonet",
  mode = "local",
  time = "48:0:0",
  update_time = 60,
  wait = FALSE,
  hold = NULL
)
}
\arguments{
\item{sample_sheet}{Path to sample sheet or data.frame.}

\item{bam_dir}{Path to bam directory. Only if sample sheet is not provided.}

\item{normal_id}{Normal sample identifier. Only if sample sheet is not provided.}

\item{patient_id}{Patient id. Only if sample sheet is not provided.}

\item{tc}{Pre-computed tumour content. Default NULL.}

\item{ploidy}{Pre-computed ploidy. Default NULL.}

\item{tmp_dir}{Path to temporary file directory.}

\item{threads}{Number of CPU cores to use. Default 3.}

\item{ram}{RAM memory for batched job. Default 4}

\item{output_dir}{Path to the output directory.}

\item{verbose}{Enables progress messages. Default False.}

\item{executor_id}{Executor ID. Default "preprocess_trento"}

\item{task_name}{Name of the task. Default "preprocess_trento"}

\item{mode}{\link{REQUIRED} Where to parallelize. Default local. Options \link{"local","batch"}}

\item{time}{\link{OPTIONAL} If batch mode. Max run time per job. Default "48:0:0"}

\item{update_time}{\link{OPTIONAL} If batch mode. Job update time in seconds. Default 60.}

\item{wait}{\link{OPTIONAL} If batch mode wait for batch to finish. Default FALSE}

\item{hold}{Job to hold on in batched mode.}
}
\description{
This function identifies a set of BAM files as tumour and normal
and processes them using the CLONET pipeline in parallel
If sample sheet is provided data has to be supplied in the following format:
patient_id   tumour  normal  version
}
\details{
Header can be ommitted if data is given in the order above and header argument is set to FALSE.
Sample sheet separator can be set using the sep argument.
}
