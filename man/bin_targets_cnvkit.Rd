% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cnvkit.R
\name{bin_targets_cnvkit}
\alias{bin_targets_cnvkit}
\title{Create binned target and anitarget beds from target BED file}
\usage{
bin_targets_cnvkit(
  sif_cnvkit = build_default_sif_list()$sif_cnvkit,
  access = build_default_reference_list()$HG19$reference$access_5k,
  bed = NULL,
  annotation = NULL,
  bin_size_antitarget = 1e+05,
  bin_size_target = 100,
  min_bin_size_antitarget = NULL,
  split = FALSE,
  short_names = FALSE,
  ...
)
}
\arguments{
\item{sif_cnvkit}{\link{REQUIRED} Path to cnvkit sif file.}

\item{access}{\link{OPTIONAL} Path to non-accessible regions to exclude. Default none}

\item{bed}{\link{REQUIRED} Path to input BED file with target regions. Default none}

\item{annotation}{\link{OPTIONAL} BED file with region annotation information. Default none}

\item{bin_size_antitarget}{\link{OPTIONAL} Average bin size for target. Default 500000.}

\item{bin_size_target}{\link{OPTIONAL} Average bin size for target. Default 120.}

\item{min_bin_size_antitarget}{\link{OPTIONAL} Average bin size for target. Default null.}

\item{split}{\link{OPTIONAL} Average bin size for target. Default null.}

\item{short_names}{\link{OPTIONAL} Use short annotation names. Default FALSE}

\item{output_name}{\link{OPTIONAL} Name for the output. If not given the name of the first tumour sample of the samples will be used.}

\item{output_dir}{\link{OPTIONAL} Split annotation names.}

\item{threads}{\link{OPTIONAL} Number of threads to split the work. Default 4}

\item{ram}{\link{OPTIONAL} RAM memory to asing to each thread. Default 4}

\item{verbose}{\link{OPTIONAL} Enables progress messages. Default False.}

\item{mode}{\link{REQUIRED} Where to parallelize. Default local. Options \link{"local","batch"}}

\item{executor_id}{Task EXECUTOR ID. Default "recalCovariates"}

\item{task_name}{Task name. Default "recalCovariates"}

\item{time}{\link{OPTIONAL} If batch mode. Max run time per job. Default "48:0:0"}

\item{update_time}{\link{OPTIONAL} If batch mode. Job update time in seconds. Default 60.}

\item{wait}{\link{OPTIONAL} If batch mode wait for batch to finish. Default FALSE}

\item{hold}{\link{OPTIONAL} HOld job until job is finished. Job ID.}
}
\description{
This function wraps around target and antitarget functiosn for CNVkit
This function generates an target BED file from an input target file.
Additional parameters can be used to exclude regions and modify the average bin size
}
\details{
For more information read:
https://cnvkit.readthedocs.io/en/stable/pipeline.html
}
