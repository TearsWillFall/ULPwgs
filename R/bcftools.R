
#' Concatenate a list of vcf files 
#'
#' This function functions concatenates a list of vcf files
#'
#' @param vcf [Optional] Path to vcf files to merge.
#' @param output_name [OPTIONAL] Name for the output. If not given the name of one of the samples will be used.
#' @param pon [Optional] Path to panel of normal.
#' @param output_dir Path to the output directory.
#' @param format Output file format. vcf/bcf
#' @param compress Output compressed file. Default file.
#' @param threads [OPTIONAL] Number of threads to split the work. Default 1
#' @param ram [OPTIONAL] RAM memory to asing to each thread. Default 4
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Task EXECUTOR ID. Default "recalCovariates"
#' @param task_name Task name. Default "recalCovariates"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export


concat_vcf=function(
    bin_bcftools=build_default_tool_binary_list()$bin_bcftools,
    bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
    bin_tabix=build_default_tool_binary_list()$bin_tabix,
    vcf="",sort=TRUE,compress=FALSE,clean=TRUE,format="vcf",
    output_name="",output_dir=".",tmp_dir=".",
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=4,mode="local",
    executor_id=make_unique_id("concatVCF"),
    task_name="concatVCF",time="48:0:0",
    update_time=60,wait=FALSE,hold=""
){

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(name=output_dir)
  job=build_job(executor_id=executor_id,task_id=task_id)
  
  
  id=""
  if(output_name!=""){
    id=output_name
  }else{
    id=get_file_name(vcf[1])
  }

  vcf=paste(vcf, collapse=' ' )



if(format=="vcf"){
    out_type=" v "
    out_file=paste0(out_file_dir,"/",id,".merged.vcf.gz")
}else if(format=="bcf"){
    out_type=" t "
    out_file=paste0(out_file_dir,"/",id,".merged.bcf.gz")
}else{
    stop("Not supported format supplied")
}

if(compressed){
    if(format=="vcf"){
        out_type=" z "
    }else if(format=="bcf"){
        out_type=" b "
    }
    out_file=paste0(out_file,".gz")
}



exec_code=paste(bin_bcftools,"concat -o",out_file, " - O ", out_type, vcf)

if(clean){
    exec_code=paste(exec_code," && rm",vcf)
}
  
  if(mode=="batch"){
       out_file_dir2=set_dir(dir=out_file_dir,name="batch")
       batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
       threads=threads,output_dir=out_file_dir2)
       exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
  }

  error=execute_job(exec_code=exec_code)
  
  if(error!=0){
    stop("gatk failed to run due to unknown error.
    Check std error for more information.")
  }

  job_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    exec_code=exec_code, 
    task_id=task_id,
    input_args = argg,
    out_file_dir=out_file_dir,
    out_files=list(
      concat_vcf=out_file)
  )
   if(sort){
        job_report[["step"]][["sortVcf"]]<-sort_vcf(
            bin_bcftools=bin_bcftools,
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            vcf=job_report$out_files$concat_vcf,
            compress=compress,format=format,
            clean=clean,
            output_name=output_name,
            output_dir=out_file_dir,tmp_dir=tmp_dir,
            batch_config=batch_config,
            threads=threads,ram=ram,mode=mode,
            executor_id=task_id,time=time,
            hold=job)
    }

  if(wait&&mode=="batch"){
    job_validator(job=job_report$job_id,time=update_time,
    verbose=verbose,threads=threads)
  }

  return(job_report)

}



#' Concatenate a list of vcf files 
#'
#' This function functions concatenates a list of vcf files
#'
#' @param vcf [Optional] Path to vcf files to merge.
#' @param output_name [OPTIONAL] Name for the output. If not given the name of one of the samples will be used.
#' @param pon [Optional] Path to panel of normal.
#' @param output_dir Path to the output directory.
#' @param format Output file format. vcf/bcf
#' @param compress Output compressed file. Default file.
#' @param threads [OPTIONAL] Number of threads to split the work. Default 1
#' @param ram [OPTIONAL] RAM memory to asing to each thread. Default 4
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Task EXECUTOR ID. Default "recalCovariates"
#' @param task_name Task name. Default "recalCovariates"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export


sort_vcf=function(
    bin_bcftools=build_default_tool_binary_list()$bin_bcftools,
    bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
    bin_tabix=build_default_tool_binary_list()$bin_tabix,
    vcf="",compress=FALSE,clean=TRUE,format="vcf",
    output_name="",output_dir=".",tmp_dir=".",
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=4,mode="local",
    executor_id=make_unique_id("concatVCF"),
    task_name="concatVCF",time="48:0:0",
    update_time=60,wait=FALSE,hold=""
){

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(name=output_dir)
  job=build_job(executor_id=executor_id,task_id=task_id)
  
  
  id=""
  if(output_name!=""){
    id=output_name
  }else{
    id=get_file_name(vcf[1])
  }


 if (!tmp_dir==""){
    tmp=paste0("-T ",tmp_dir)
  }


  

if(format=="vcf"){
    out_type=" v "
    out_file=paste0(out_file_dir,"/",id,".sorted.vcf.gz")
}else if(format=="bcf"){
    out_type=" t "
    out_file=paste0(out_file_dir,"/",id,".sorted.bcf.gz")
}else{
    stop("Not supported format supplied")
}

if(compressed){
    if(format=="vcf"){
        out_type=" z "
    }else if(format=="bcf"){
        out_type=" b "
    }
    out_file=paste0(out_file,".gz")
}



exec_code=paste0(bin_bcftools," sort ",vcf , paste0(" -m ",ram,"G "), tmp_dir, " -o ",out_file, " -O ",out_type)

if(clean){
    exec_code=paste(exec_code," && rm",vcf)
}
  
  if(mode=="batch"){
       out_file_dir2=set_dir(dir=out_file_dir,name="batch")
       batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
       threads=threads,output_dir=out_file_dir2)
       exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
  }

  error=execute_job(exec_code=exec_code)
  
  if(error!=0){
    stop("gatk failed to run due to unknown error.
    Check std error for more information.")
  }

  job_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    exec_code=exec_code, 
    task_id=task_id,
    input_args = argg,
    out_file_dir=out_file_dir,
    out_files=list(
      sorted_vcf=out_file)
  )


  if(wait&&mode=="batch"){
    job_validator(job=job_report$job_id,time=update_time,
    verbose=verbose,threads=threads)
  }

  return(job_report)

}


