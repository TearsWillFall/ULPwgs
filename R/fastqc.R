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
#' @param executor_id Executor ID. Default "fastQC"
#' @param task_name Name of the task. Default "fastQC"
#' @param ram RAM memory for batched job. Default 4
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param hold Job to hold on in batched mode.
#' @export


qc_fastqc=function (bin_path="tools/FastQC/bin/fastqc",file_R1="",file_R2="",
output_dir="",verbose=FALSE,executor_id=make_unique_id("fastQC"),
task_name="fastQC",mode="local",threads=3,ram=4,
time="48:0:0",update_time=60,wait=FALSE,hold=""){

  task_id=make_unique_id(task_name)

  out_file_dir=set_dir(dir=output_dir,name="fastqc_reports")


  if (!check_missing(file_R2)){

    exec_code=paste(bin_path,"-o ", out_file_dir,"-t ",threads,"--noextract",file_R1,file_R2)


  }else{
    exec_code=paste(bin_path,"-o ", out_file_dir,"-t ",threads,"--noextract",file_R1)
  }

  job=build_job(executor_id=executor_id,task_id=task_id)

  if(mode=="batch"){
    out_file_dir2=set_dir(dir=out_file_dir,name="batch")
    batch_code=build_job_exec(job=job,time=time,ram=ram,threads=threads,output_dir=out_file_dir2,hold=hold)
    exec_code=paste0("echo 'source ~/.bashrc;",exec_code,"'|",batch_code)
  }

  if(verbose){
    print_verbose(exec_code=exec_code)
  }

  error=system(exec_code)
  if(error!=0){
    stop("fastqc failed to run due to unknown error.
    Check std error for more information.")
  }


  job_report=build_job_report(job_id=job,executor_id=executor_id,
  task_id=task_id,job_order=1,out_files=list(
    R1=list( 
      ZIP=paste0(out_file_dir,"/",get_file_name(file_R1),"_fastqc.zip"),
      HTML=paste0(out_file_dir,"/",get_file_name(file_R2),"_fastqc.html"),
    R2=list(
      ZIP=paste0(out_file_dir,"/",get_file_name(file_R2),"_fastqc.zip"),
      HTML=paste0(out_file_dir,"/",get_file_name(file_R2),"_fastqc.html")
        ) 
      )
    )
  )
 
  if(wait&&mode=="batch"){
    job_validator(job=job_report$job_id,
    time=update_time,verbose=verbose,threads=threads)
  }

  return(job_report)
  
  }



