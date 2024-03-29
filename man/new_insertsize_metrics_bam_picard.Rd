% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/picard.R
\name{new_insertsize_metrics_bam_picard}
\alias{new_insertsize_metrics_bam_picard}
\title{Generate BAM Insert Size Metrics}
\usage{
new_insertsize_metrics_bam_picard(
  bin_picard = build_default_tool_binary_list()$bin_picard,
  bam = NULL,
  deviations = NULL,
  min_width = NULL,
  width = NULL,
  ...
)
}
\arguments{
\item{bin_picard}{Path to bwa executable. Default path tools/samtools/samtools.}

\item{bam}{Path to the input file with the sequence.}

\item{output_dir}{Path to the output directory.}

\item{verbose}{Enables progress messages. Default False.}

\item{ram}{RAM memory to use in GB. Default 4.}

\item{tmp_dir}{Path to TMP directory. Default .}

\item{executor_id}{\link{OPTIONAL} Task executor name. Default "recalCovariates"}

\item{task_name}{\link{OPTIONAL} Task name. Default "recalCovariates"}

\item{time}{\link{OPTIONAL} If batch mode. Max run time per job. Default "48:0:0"}

\item{update_time}{\link{OPTIONAL} If batch mode. Job update time in seconds. Default 60.}

\item{wait}{\link{OPTIONAL} If batch mode wait for batch to finish. Default FALSE}

\item{hold}{\link{OPTIONAL} HOld job until job is finished. Job ID.}
}
\description{
Generate BAM Insert Size Metrics
}
