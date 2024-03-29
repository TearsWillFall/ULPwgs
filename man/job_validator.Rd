% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/batch_constructor.R
\name{job_validator}
\alias{job_validator}
\title{Validate job submited in batch on sungrid based cluster}
\usage{
job_validator(job = "", time = 10, verbose = FALSE, threads = 3)
}
\arguments{
\item{job}{Name of job or jobs.}

\item{time}{Time in seconds between checks. Default 10.}

\item{verbose}{\link{OPTIONAL} Enables progress messages. Default False.}

\item{threads}{\link{OPTIONAL} Number of threads to use to. Default 3.}
}
\description{
Validate job submited in batch on sungrid based cluster
}
