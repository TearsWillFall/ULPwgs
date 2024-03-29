% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/batch_constructor.R
\name{build_job_exec}
\alias{build_job_exec}
\title{Build job executor in SGE}
\usage{
build_job_exec(
  job = "",
  time = "48:0:0",
  ram = 3,
  threads = 1,
  output_dir = ".",
  hold = NULL,
  bypass = FALSE,
  wd = getwd(),
  array = ""
)
}
\arguments{
\item{job}{Name of job or jobs.}

\item{time}{Time in seconds between checks. Default 10.}

\item{ram}{\link{OPTIONAL} RAM memory per thread requested. Default 4.}

\item{threads}{\link{OPTIONAL} Number of threads to use to. Default 3.}

\item{output_dir}{\link{OPTIONAL} PATH to output directory.}

\item{hold}{Job IDs to hold job.}

\item{verbose}{\link{OPTIONAL} Enables progress messages. Default False.}
}
\description{
Build job executor in SGE
}
