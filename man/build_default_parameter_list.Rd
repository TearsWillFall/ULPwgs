% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/config.R
\name{build_default_parameter_list}
\alias{build_default_parameter_list}
\title{Build default tool configuration}
\usage{
build_default_parameter_list(
  parameters = list(parameter = c("name", "step", "mode", "time", "threads", "ram",
    "batch_config", "args"), text = c("Step Name: ", "Run: ", "Run Mode: ", "Run Time: ",
    "Threads: ", "Ram: ", "Batch Config: ", "Tool Parameters : "))
)
}
\arguments{
\item{config}{List with name,threads,ram,step,mode,verbose,time,args config}
}
\description{
Build default tool configuration
}
