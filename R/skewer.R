#' Adapter sequence trimmer using skewer
#'
#' Detects and removes adapter sequences from single and paired read data
#'
#' @param file_R1 Path to the input file with the sequence.
#' @param file_R2 [Optional] Path to the input with the reverse read sequence.
#' @param xadapt [Optional] Adapter sequence/file.
#' @param yadapt [Optional] Adapter sequence/file.
#' @param bin_path Path to skewer executable. Default path tools/skewer/skewer.
#' @param threads Number of CPU cores to use. Default 3.
#' @param mean_quality Minimum mean quality of reads to be kept.Dedault 0
#' @param min_length Minimum length of reads to keep.Default 35
#' @param max_length Maximum length of reads to keep.Default NA.
#' @param output_dir Path to the output directory.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param verbose Enables progress messages. Default False.
#' @export


trimming_skewer=function(bin_path="tools/skewer/skewer",file_R1="",file_R2="",xadapt=NA,
yadapt=NA,threads=3,output_dir="",verbose=FALSE,mean_quality=0,min_length=35,max_length=NA,
mode="local",time="48:0:0",update_time=60,wait=FALSE){

  out_file_dir=set_dir(dir=output_dir,name="skewer_reports")

  if ((!is.na(xadapt)) & (!is.na(yadapt))){
    func=paste(bin_path,"-m tail -t",threads,"-x", xadapt,"-y", yadapt,"-Q",mean_quality,"-l",min_length)
  }

  else{
    func=paste(bin_path,"-m tail -t",threads,"-Q",mean_quality,"-l",min_length)
  }
  if(!is.na(max_length)){
    func=paste(func,"-L",max_length)
  }

  if (!file_R2==""){
    exec_code=paste(func,"-z -f sanger --quiet -o",paste0(out_file_dir,"/",intersect_file_name(file_R1,file_R2)),file_R1,file_R2)
  }else{
    exec_code=paste(func,"-z -f sanger --quiet -o",paste0(out_file_dir,"/",get_file_name(file_R1)),file_R1)
  }

  job_name="skewer"

  if(mode=="batch"){
    out_file_dir2=set_dir(dir=out_file_dir,name="batch")
    job_name=build_job(task="skewer",input_name=get_file_name(bam))
    batch_code=build_job_exec(job_name=job_name,time=time,ram=ram,threads=1,output_dir=out_file_dir2)
    exec_code=paste0("echo 'source ~/.bashrc;",exec_code,"'|",batch_code)
  }

  if(verbose){
      print(exec_code)
    }
  error=system(exec_code)
  if(error!=0){
    stop("skewer failed to run due to unknown error.
    Check std error for more information.")
  }

   if(wait&&mode=="batch"){
    batch_validator(job=paste(job_name,"_batch"),
    time=update_time,verbose=verbose,threads=threads)
  }

  return(job_name)
  
}