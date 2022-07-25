
#' Mark duplicated reads
#'
#' This function marks duplicated reads (artifacts) found in aligned sequences.
#'
#' @param bin_path Path to picard executable. Default path tools/picard/build/libs/picard.jar.
#' @param bam Path to the input file with the aligned sequence.
#' @param output_dir Path to the output directory.
#' @param tmp_dir Path to tmp directory.
#' @param verbose Enables progress messages. Default False.
#' @param remove_duplicates Do not write duplicates to the output file. Default FALSE
#' @param hnd Maximum number of file handles. Default 1000.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor [OPTIONAL] Task executor name. Default "recalCovariates"
#' @param task [OPTIONAL] Task name. Default "recalCovariates"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export

markdups_picard=function(bin_path="tools/picard/build/libs/picard.jar",bam="",
output_dir="",verbose=FALSE,hnd=1000,threads,ram=4,tmp_dir="",remove_duplicates=TRUE,
mode="local",executor=make_unique_id("markDups"),task="markDups",
time="48:0:0",update_time=60,wait=FALSE,hold=""){

    argg <- as.list(environment())
  task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir,name="markdups_reports")

    tmp=""
    if (!tmp_dir==""){
      tmp=paste0("TMP_DIR=",tmp_dir)
    }

    if(remove_duplicates){
      remove_duplicates=" REMOVE_DUPLICATES=true "
    }else{
      remove_duplicates=" REMOVE_DUPLICATES=false "
    }


    job=build_job(executor_id=executor_id,task_id=task_id)

    exec_code=paste0("java -Xmx",ram,"g", " -Djava.io.tmpdir=",tmp_dir," -jar ",
      bin_path," MarkDuplicates I=",bam, " O=",paste0(out_file_dir,"/",
      get_file_name(bam),".rmdup.",get_file_ext(bam)),
      " M=",paste0(out_file_dir,"/",
      get_file_name(bam),".picard_rmdup.txt"),
      remove_duplicates, " AS=true VALIDATION_STRINGENCY=LENIENT ",
      paste0("MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=",hnd)," ",tmp)
    
    
  if(mode=="batch"){
       out_file_dir2=set_dir(dir=out_file_dir,name="batch")
       exec_batch=build_job_exec(job=job,hold=hold,time=time,ram=ram,
       threads=threads,output_dir=out_file_dir2)
       exec_code=paste("echo 'source ~/.bashrc;",exec_code,"'|",exec_batch)
  }
    
  if(verbose){
       print_verbose(job=job,arg=argg,exec_code=exec_code)
    }

    error=system(exec_code)
    
    if(error!=0){
      stop("markdups failed to run due to unknown error.
      Check std error for more information.")
    }

  job_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    task_id=task_id, 
    input_args = argg,
    out_file_dir=out_file_dir,
      out_file=list(
        )
  )


   if(wait&&mode=="batch"){
    job_validator(job=job_report$job_id,
    time=update_time,verbose=verbose,threads=threads)
   }

    return(job_report)

  }

#' Generate BAM General Summary Metrics
#'
#'
#' @param bam Path to the input file with the sequence.
#' @param bin_path Path to bwa executable. Default path tools/samtools/samtools.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param ram RAM memory to use in GB. Default 4.
#' @param tmp_dir Path to TMP directory. Default .
#' @param executor [OPTIONAL] Task executor name. Default "recalCovariates"
#' @param task [OPTIONAL] Task name. Default "recalCovariates"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export

summary_metrics_bam_picard=function(bin_path="tools/picard/build/libs/picard.jar",bam="",output_dir="",
verbose=FALSE,tmp_dir=".",threads=3,ram=4,mode="local",executor=make_unique_id("summaryMetrics"),task="summaryMetrics",
time="48:0:0",update_time=60,wait=FALSE,hold=""){


  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir,name="summary")

  tmp=""
  if (!tmp_dir==""){
    tmp=paste0(" TMP_DIR=",tmp_dir)
  }

  exec_code=paste0("java -Xmx",ram,"g", " -Djava.io.tmpdir=",tmp_dir,
        " -jar ",bin_path," CollectAlignmentSummaryMetrics ",
        "VALIDATION_STRINGENCY=SILENT I=",bam," O=",paste0(out_file_dir,"/",get_file_name(bam),".picard_summary.txt "),tmp)
  
  job=build_job(executor_id=executor_id,task_id=task_id)

  if(mode=="batch"){
       out_file_dir2=set_dir(dir=out_file_dir,name="batch")
       exec_batch=build_job_exec(job=job,hold=hold,time=time,ram=ram,
       threads=threads,output_dir=out_file_dir2)
       exec_code=paste("echo 'source ~/.bashrc;",exec_code,"'|",exec_batch)
  }


 if(verbose){
       print_verbose(job=job,arg=argg,exec_code=exec_code)
    }
  error=system(exec_code)
  if(error!=0){
    stop("picard failed to run due to unknown error.
    Check std error for more information.")
  }

  job_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    task_id=task_id, 
    input_args = argg,
    out_file_dir=out_file_dir,
      out_file=list(
        )
)


 if(wait&&mode=="batch"){
    job_validator(job=job_report$job_id,
    time=update_time,verbose=verbose,threads=threads)
  }
  return(job_report)
}



#' Generate BAM Insert Size Metrics
#'
#'
#' @param bam Path to the input file with the sequence.
#' @param bin_path Path to bwa executable. Default path tools/samtools/samtools.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param ram RAM memory to use in GB. Default 4.
#' @param tmp_dir Path to TMP directory. Default .
#' @param executor [OPTIONAL] Task executor name. Default "recalCovariates"
#' @param task [OPTIONAL] Task name. Default "recalCovariates"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export

insertsize_metrics_bam_picard=function(bin_path="tools/picard/build/libs/picard.jar",
bam="",output_dir="",verbose=FALSE,tmp_dir=".",threads=1,ram=4,
mode="local",executor=make_unique_id("insertsizeMetrics"),task="insertsizeMetrics",
time="48:0:0",update_time=60,wait=FALSE,hold=""){



  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir,name="insertsize")


  tmp=""
  if (!tmp_dir==""){
    tmp=paste0(" TMP_DIR=",tmp_dir)
  }

  exec_code=paste0("java -Xmx",ram,"g", " -Djava.io.tmpdir=",tmp_dir," -jar ",
      bin_path," CollectInsertSizeMetrics ","VALIDATION_STRINGENCY=SILENT I=",
      bam," O=",paste0(out_file_dir,"/",get_file_name(bam),".picard_insert_size.txt")," H=",
      paste0(out_file_dir,"/",get_file_name(bam),".picard_insert_size.pdf "),tmp)

  job=build_job(executor_id=executor_id,task_id=task_id)
  if(mode=="batch"){
       out_file_dir2=set_dir(dir=out_file_dir,name="batch")
       exec_batch=build_job_exec(job=job,hold=hold,time=time,ram=ram,
       threads=threads,output_dir=out_file_dir2)
       exec_code=paste("echo 'source ~/.bashrc;",exec_code,"'|",exec_batch)
  }

  if(verbose){
       print_verbose(job=job,arg=argg,exec_code=exec_code)
    }
  error=system(exec_code)
  if(error!=0){
    stop("picard failed to run due to unknown error.
    Check std error for more information.")
  }

  job_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    task_id=task_id, 
    input_args = argg,
    out_file_dir=out_file_dir,
      out_file=list(
        )
)


 if(wait&&mode=="batch"){
    job_validator(job=job_report$job_id,
    time=update_time,verbose=verbose,threads=threads)
  }
  return(job_report)
}




#' Generate BAM Summary for Targeted data
#'
#'
#' @param bam Path to the input file with the sequence.
#' @param bin_path Path to bwa executable. Default path tools/samtools/samtools.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param ram RAM memory to use in GB. Default 4.
#' @param tmp_dir Path to TMP directory. Default .
#' @param bi Bait capture target interval for panel data. Interval format.
#' @param ti Primary target intervals for panel data. Interval format.
#' @param executor [OPTIONAL] Task executor name. Default "recalCovariates"
#' @param task [OPTIONAL] Task name. Default "recalCovariates"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export

tg_summary_metrics_bam_picard=function(bin_path="tools/picard/build/libs/picard.jar",bam="",output_dir="",
verbose=FALSE,tmp_dir=".",threads=1,ram=4,bi="",ti="",mode="local",executor_id=make_unique_id("TGsummaryMetrics"),
task_name="TGsummaryMetrics",time="48:0:0",update_time=60,wait=FALSE,hold=""){

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir,name="summary")


  tmp=""
  if (!tmp_dir==""){
    tmp=paste0(" TMP_DIR=",tmp_dir)
  }

  out_file_ts=paste0(out_file_dir,"/",get_file_name(bam),".picard_TS.txt")
  out_file=paste0(out_file_dir,"/",get_file_name(bam),".picard_CollectHSmetrics.txt ")
  exec_code=paste0("java -Xmx",ram,"g", " -Djava.io.tmpdir=",tmp_dir,
        " -jar ",bin_path," CollectHsMetrics VALIDATION_STRINGENCY=SILENT BI=",
        bi," TI=",ti," I=",bam," THEORETICAL_SENSITIVITY_OUTPUT=",
        out_file_ts," O=",out_file,tmp)
        
 
 
 
job=build_job(executor_id=executor_id,task_id=task_id)
 if(mode=="batch"){
       out_file_dir2=set_dir(dir=out_file_dir,name="batch")
       exec_batch=build_job_exec(job=job,hold=hold,time=time,ram=ram,
       threads=threads,output_dir=out_file_dir2)
       exec_code=paste("echo 'source ~/.bashrc;",exec_code,"'|",exec_batch)
  }


   if(verbose){
       print_verbose(job=job,arg=argg,exec_code=exec_code)
   }

  error=system(exec_code)
  if(error!=0){
    stop("picard failed to run due to unknown error.
    Check std error for more information.")
  }

  job_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    task_id=task_id, 
    input_args = argg,
    out_file_dir=out_file_dir,
      out_file=list(
        ts=out_file_ts,
        summary=out_file
        )
  )



  if(wait&&mode=="batch"){
    job_validator(job=job_report$job_id,
    time=update_time,verbose=verbose,threads=threads)
  }
  return(job_report)


}


#' Generate BAM Summary for RNAseq data
#'
#'
#' @param bam Path to the input file with the sequence.
#' @param bin_path Path to bwa executable. Default path tools/samtools/samtools.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param ram RAM memory to use in GB. Default 4.
#' @param tmp_dir Path to TMP directory. Default .
#' @param ri Ribosomal interval for RNAseq data. Interval format.
#' @param ref_flat Path to flat refrence. Interval format.
#' @param executor [OPTIONAL] Task executor name. Default "recalCovariates"
#' @param task [OPTIONAL] Task name. Default "recalCovariates"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export

rnaseq_summary_metrics_bam_picard=function(bin_path="tools/picard/build/libs/picard.jar",
bam="",output_dir="",verbose=FALSE,tmp_dir=".",threads=1,ram=4,ri="",ref_flat="",mode="local",
executor=make_unique_id("RNAsummaryMetrics"),task="RNAsummaryMetrics",time="48:0:0",
update_time=60,wait=FALSE,hold=""){


  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir,name="summary")

  tmp=""
  if (!tmp_dir==""){
    tmp=paste0(" TMP_DIR=",tmp_dir)
  }

  exec_code=paste0("java -Xmx",ram,"g", " -Djava.io.tmpdir=",tmp_dir," -jar ",bin_path,
        " CollectRnaSeqMetrics VALIDATION_STRINGENCY=SILENT STRAND_SPECIFICITY='NONE' REF_FLAT=",
         ref_flat, " RIBOSOMAL_INTERVALS=",ri,
         " I=",bam," O=",paste0(out_file_dir,"/",get_file_name(bam),".CollectRNAseqMetrics.txt "),tmp)
  
  
  
job=build_job(executor_id=executor_id,task_id=task_id)
  if(mode=="batch"){
       out_file_dir2=set_dir(dir=out_file_dir,name="batch")
       exec_batch=build_job_exec(job=job,hold=hold,time=time,ram=ram,
       threads=threads,output_dir=out_file_dir2)
       exec_code=paste("echo 'source ~/.bashrc;",exec_code,"'|",exec_batch)
  }


  if(verbose){
       print_verbose(job=job,arg=argg,exec_code=exec_code)
  }

  error=system(exec_code)
  if(error!=0){
    stop("picard failed to run due to unknown error.
    Check std error for more information.")
  }

  job_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    task_id=task_id, 
    input_args = argg,
    out_file_dir=out_file_dir,
      out_file=list(
        )
)


  if(wait&&mode=="batch"){
    job_validator(job=job_report$job_id,
    time=update_time,verbose=verbose,threads=threads)
  }
  return(job_report)
}



#' Generate BAM Summary for WGS data
#'
#'
#' @param bam Path to the input file with the sequence.
#' @param bin_path Path to bwa executable. Default path tools/samtools/samtools.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param ram RAM memory to use in GB. Default 4.
#' @param tmp_dir Path to TMP directory. Default .
#' @param executor [OPTIONAL] Task executor name. Default "recalCovariates"
#' @param task [OPTIONAL] Task name. Default "recalCovariates"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export

wgs_summary_metrics_bam_picard=function(bin_path="tools/picard/build/libs/picard.jar",
bam="",output_dir="",verbose=FALSE,tmp_dir=".",threads=1,ram=4,mode="local",
executor=make_unique_id("WGSsummaryMetrics"),task="WGSsummaryMetrics",time="48:0:0",
update_time=60,wait=FALSE,hold=""){


  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir,name="summary")


  tmp=""
  if (!tmp_dir==""){
    tmp=paste0(" TMP_DIR=",tmp_dir)
  }

  exec_code=paste0("java -Xmx",ram,"g", " -Djava.io.tmpdir=",tmp_dir," -jar ",
            bin_path," CollectWgsMetrics VALIDATION_STRINGENCY=SILENT MINIMUM_MAPPING_QUALITY=",
            mapq," I=",bam," O=",paste0(out_file_dir,"/",get_file_name(bam),".picard_wgs_q00.txt "),tmp)
  
job=build_job(executor_id=executor_id,task_id=task_id)
  if(mode=="batch"){
       out_file_dir2=set_dir(dir=out_file_dir,name="batch")
       exec_batch=build_job_exec(job=job,hold=hold,time=time,ram=ram,
       threads=threads,output_dir=out_file_dir2)
       exec_code=paste("echo 'source ~/.bashrc;",exec_code,"'|",exec_batch)
  }


 if(verbose){
       print_verbose(job=job,arg=argg,exec_code=exec_code)
    }
  error=system(exec_code)
  if(error!=0){
    stop("picard failed to run due to unknown error.
    Check std error for more information.")
  }


  job_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    task_id=task_id, 
    input_args = argg,
    out_file_dir=out_file_dir,
      out_file=list(
        )
)

 if(wait&&mode=="batch"){
    job_validator(job=job_report$job_id,
    time=update_time,verbose=verbose,threads=threads)
  }
  return(job_report)
}