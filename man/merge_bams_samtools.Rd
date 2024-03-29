% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{merge_bams_samtools}
\alias{merge_bams_samtools}
\title{Merge BAM files in directory}
\usage{
merge_bams_samtools(
  bin_samtools = build_default_tool_binary_list()$bin_samtools,
  bams = "",
  output_name = "",
  verbose = TRUE,
  threads = 3,
  ram = 4,
  executor = make_unique_id("mergeBAMs"),
  task = "mergeBAMs",
  mode = "local",
  time = "48:0:0",
  update_time = 60,
  wait = FALSE,
  hold = NULL
)
}
\arguments{
\item{bin_samtools}{Path to samtools binary.}

\item{output_name}{Output file name}

\item{verbose}{Enables progress messages. Default False.}

\item{threads}{Number of threads to use.Default 3.}

\item{bam}{Path to the input bam file with the sequence.}

\item{bam_dir}{Path to directory with BAM files to merge.}
}
\description{
This function takes a BAM file and merges it with others found within a
directory. This function is still Work In Progress (WIP).
}
