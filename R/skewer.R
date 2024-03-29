#' Adapter sequence trimmer using skewer
#'
#' Detects and removes adapter sequences from single and paired read data
#'
#' @param bin_skewer Path to skewer executable. Default path tools/skewer/skewer.
#' @param file_R1 Path to the input file with the sequence.
#' @param file_R2 [Optional] Path to the input with the reverse read sequence.
#' @param xadapt [Optional] Adapter sequence/file.
#' @param yadapt [Optional] Adapter sequence/file.
#' @param threads Number of CPU cores to use. Default 3.
#' @param mean_quality The lowest mean quality value allowed before trimming. Default 0
#' @param min_length Minimum length of reads allowed after trimming. Default 18
#' @param max_length Maximum length of reads allowed after trimming. Default No limit.
#' @param output_dir Path to the output directory.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Executor ID. Default "trimmingSkewer"
#' @param task_name Name of the task. Default "trimmingSkewer"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param verbose Enables progress messages. Default False.
#' @export


trimming_skewer=function(
  bin_skewer=build_default_tool_binary_list()$bin_skewer,
  file_R1="",file_R2="",xadapt="",
  yadapt="",mean_quality=0,min_length=18,max_length="",
  threads=3,ram=4,output_dir=".",verbose=FALSE,
  batch_config=build_default_preprocess_config(),
  output_name="",mode="local",
  executor_id=make_unique_id("trimmingSkewer"),
  task_name="trimmingSkewer",time="48:0:0",
  update_time=60,wait=FALSE,hold=NULL
){


  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir,name="skewer_reports")

  func=paste(bin_skewer,"-m tail -t",threads,"-Q",mean_quality,"-l",min_length)
  
  if(!check_missing(var=xadapt) &&!check_missing(var=yadapt)){
    func=paste(func,"-x", xadapt,"-y", yadapt)
  }

  if(!check_missing(max_length)){
    func=paste(func,"-L",max_length)
  }

  if (!check_missing(file_R2)){
    out_file=paste0(out_file_dir,"/",ifelse(output_name=="",
    intersect_file_name(file_R1,file_R2),output_name))
    exec_code=paste(func,"-z -f sanger --quiet -o",out_file,file_R1,file_R2)
  }else{
    out_file=paste0(out_file_dir,"/",ifelse(output_name=="",
    get_file_name(file_R1),output_name))
    exec_code=paste(func,"-z -f sanger --quiet -o",out_file,file_R1)
  }


  job=build_job(executor_id=executor_id,task_id=task_id)
  
  if(mode=="batch"){
    out_file_dir2=set_dir(dir=out_file_dir,name="batch")
    batch_code=build_job_exec(job=job,time=time,ram=ram,threads=threads,output_dir=out_file_dir2)
    exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
  }
  

  if(verbose){
      print_verbose(job=job,arg=argg,exec_code=exec_code)
    }
  error=execute_job(exec_code=exec_code)
  if(error!=0){
    stop("skewer failed to run due to unknown error.
    Check std error for more information.")
  }

  job_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    exec_code=exec_code,
    task_id=task_id,
    input_args=argg,
    out_file_dir=out_file_dir,
    out_files=list(
      r1=paste0(out_file,"-trimmed-pair1.fastq.gz"),
      r2=paste0(out_file,"-trimmed-pair2.fastq.gz"),
      log=paste0(out_file,"-trimmed.log")
    )
  )

  if(wait&&mode=="batch"){
    batch_validator(job=job_report$job_id,
    time=update_time,verbose=verbose,threads=threads)
  }

  return(job_report)

}