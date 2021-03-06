% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{read_counter}
\alias{read_counter}
\title{Generate a WIG file}
\usage{
read_counter(
  bin_path = "tools/samtools/samtools",
  bin_path2 = "tools/hmmcopy_utils/bin/readCounter",
  chrs = c(1:22, "X", "Y"),
  win = 5e+05,
  bam = "",
  output_dir = "",
  verbose = FALSE,
  threads = 3
)
}
\arguments{
\item{bin_path}{Path to readCounter executable. Default path tools/samtools/samtools.}

\item{bin_path2}{Path to readCounter executable. Default path tools/hmmcopy_utils/bin/readCounter.}

\item{chrs}{String of chromosomes to include. c()}

\item{win}{Size of non overlaping windows. Default 500000.}

\item{bam}{Path to the BAM file .}

\item{output_dir}{Path to the output directory.}

\item{verbose}{Enables progress messages. Default False.}

\item{threads}{Number of threads to use. Default 3}
}
\description{
This function generates a WIG file.
}
