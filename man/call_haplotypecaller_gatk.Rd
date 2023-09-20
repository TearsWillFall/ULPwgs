% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gatk.R
\name{call_haplotypecaller_gatk}
\alias{call_haplotypecaller_gatk}
\title{Variant Calling using HaplotypeCaller}
\usage{
call_haplotypecaller_gatk(
  sif_gatk = build_default_sif_list()$sif_gatk,
  bin_samtools = build_default_tool_binary_list()$bin_samtools,
  ref_genome = build_default_reference_list()$HG19$reference$genome,
  bam = NULL,
  chromosomes = c(1:22, "X", "Y"),
  region = NULL,
  score = TRUE,
  score_CNN = "CNN_1D",
  filter = TRUE,
  ...
)
}
\arguments{
\item{sif_gatk}{\link{REQUIRED} Path to gatk sif file.}

\item{ref_genome}{\link{REQUIRED} Path to reference genome fasta file.}

\item{bam}{\link{OPTIONAL} Path to normal BAM file.}

\item{region}{\link{REQUIRED} Genomic position in samtools format chr:start-end.}
}
\description{
This function functions calls HaplotypeCaller for variant calling.
If a vector of normal samples are provided these will be processed in multi-sample mode.
To run in normal mode suppply a single normal sample.
TO DO// Implement mitochondrial mode feature
}
\details{
For more information read:
https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller
}