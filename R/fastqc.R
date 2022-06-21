#' Generate a quality control (QC) report from a fastaqc file
#'
#' This function takes a set of sequence files (fastq,SAM,BAM...) and
#' returns a report in HTML format.
#'
#'
#' @param file_R1 Path to the input file with the sequence.
#' @param file_R2 [Optional] Path to the input with the reverse read sequence.
#' @param bin_path Path to fastQC executable. Default path tools/FastQC/bin/fastqc.
#' @param threads Number of CPU cores to use. Default 3.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor Name of the executor. Default "fastQC"
#' @param task Name of the task. Default "fastQC"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @export


qc_fastqc=function (bin_path="tools/FastQC/bin/fastqc",file_R1="",file_R2="",threads=3,
output_dir="",verbose=FALSE,executor=make_unique_id("fastQC"),task="fastQC",mode="local",
time="48:0:0",update_time=60,wait=FALSE,hold=""){

  out_file_dir=set_dir(dir=output_dir,name="fastqc_reports")


  if (!file_R2==""){

    exec_code=paste(bin_path,"-o ", out_file_dir,"-t ",threads,"--noextract",file_R1,file_R2)


  }else{
    exec_code=paste(bin_path,"-o ", out_file_dir,"-t ",threads,"--noextract",file_R1)
  }

  job=build_job(task=make_unique_id(task),executor=executor)
  if(mode=="batch"){
    out_file_dir2=set_dir(dir=out_file_dir,name="batch")
    batch_code=build_job_exec(job_name=job_name,time=time,ram=ram,threads=threads,output_dir=out_file_dir2)
    exec_code=paste0("echo 'source ~/.bashrc;",exec_code,"'|",batch_code)
  }

  if(verbose){
      print(exec_code)
  }

  error=system(exec_code)
  if(error!=0){
    stop("fastqc failed to run due to unknown error.
    Check std error for more information.")
  }

  if(wait&&mode=="batch"){
    job_validator(job=job,
    time=update_time,verbose=verbose,threads=threads)
  }

  return(job)


}
