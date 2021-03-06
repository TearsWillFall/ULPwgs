% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{fastqc}
\alias{fastqc}
\title{Generate a quality control (QC) report from a fastaqc file}
\usage{
fastqc(
  bin_path = "tools/FastQC/bin/fastqc",
  file_R1 = "",
  file_R2 = "",
  n_cores = 3,
  output_dir = "",
  verbose = FALSE
)
}
\arguments{
\item{bin_path}{Path to fastQC executable. Default path tools/FastQC/bin/fastqc.}

\item{file_R1}{Path to the input file with the sequence.}

\item{file_R2}{\link{Optional} Path to the input with the reverse read sequence.}

\item{n_cores}{Number of CPU cores to use. Default 3.}

\item{output_dir}{Path to the output directory.}

\item{verbose}{Enables progress messages. Default False.}
}
\description{
This function takes a set of sequence files (fastq,SAM,BAM...) and
returns a report in HTML format.
}
