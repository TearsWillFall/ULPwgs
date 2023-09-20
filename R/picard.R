

#' Wrapper around gatk GatherBamFiles function
#'
#' This functions collects the Recalibration reports generated from scattered parallel_apply_BQSR output
#' This function wraps around gatk GatherBamFiles function.
#' For more information about this function: https://gatk.broadinstitute.org/hc/en-us/articles/360037055512-GatherBamFiles-Picard-
#'
#' @param bin_picard [REQUIRED] Path to picard executable. Default tools/picard/build/libs/picard.jar.
#' @param bam [REQUIRED] Path to BAM file/s.
#' @param output_name [OPTIONAL] Name for the output file name.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param clean Clean input files. Default TRUE.
#' @param tmp_dir [OPTIONAL] Path to the temporary directory.
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param threads [OPTIONAL] Number of threads to split the work. Default 4
#' @param ram [OPTIONAL] RAM memory to asing to each thread. Default 4
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Task EXECUTOR ID. Default "gatherBAM"
#' @param task_name Task name. Default "gatherBAM"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Show job updates every update time. Default 60
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export

gather_bam_files_picard=function(
  bin_picard=build_default_tool_binary_list()$bin_picard,
  bam="",output_name="File",output_dir=".",
  verbose=FALSE,batch_config=build_default_preprocess_config(),
  threads=4, ram=4, mode="local",clean=FALSE,
  executor_id=make_unique_id("gatherBAM"),task_name="gatherBAM",
  time="48:0:0",update_time=60,wait=FALSE,hold=NULL
){

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir)

  
  out_file=paste0(out_file_dir,"/",output_name,".bam")
  bam=bam[order(as.numeric(lapply(lapply(lapply(lapply(lapply(lapply(basename(bam),
    FUN=strsplit,split="\\."),FUN="[[",index=1),FUN="[",index=2),
    FUN=strsplit,split="__"),FUN="[[",index=1),FUN="[",index=1)))]


  exec_code=paste0("java -jar ",bin_picard," GatherBamFiles ",
    paste0(" I=",bam,collapse=" ")," O=",out_file)

  if(clean){
    exec_code=paste(exec_code," && rm",paste(paste0(bam,"*"),collapse=" "))
  }

 
  job=build_job(executor_id=executor_id,task_id=task_id)

  if(mode=="batch"){
       out_file_dir2=set_dir(dir=out_file_dir,name="batch")
       batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
       threads=threads,output_dir=out_file_dir2)
       exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
  }

  if(verbose){
   print_verbose(job=job,arg=argg,exec_code=exec_code)
  }

  error=execute_job(exec_code=exec_code)

  if(error!=0){
    stop("picard failed to run due to unknown error.
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
      bam=out_file)
    )

  if(wait&&mode=="batch"){
    job_validator(job=job_report$job_id,time=update_time,
    verbose=verbose,threads=threads)
  }
  return(job_report)
}



#' Mark duplicated reads
#'
#' This function marks duplicated reads (artifacts) found in aligned sequences.
#'
#' @param bin_picard Path to picard executable. Default path tools/picard/build/libs/picard.jar.
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

markdups_picard=function(
  bin_picard=build_default_tool_binary_list()$bin_picard,bam="",
  output_dir=".",verbose=FALSE,
  batch_config=build_default_preprocess_config(),
  hnd=1000,threads,ram=4,tmp_dir=".",
  remove_duplicates=TRUE,mode="local",executor_id=make_unique_id("markDups"),
  task_name="markDups",time="48:0:0",
  update_time=60,wait=FALSE,hold=NULL
){

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

    out_file=paste0(out_file_dir,"/",get_file_name(bam),".rmdup.",get_file_ext(bam))
    out_file_md=paste0(out_file_dir,"/",get_file_name(bam),".picard_rmdup.txt")

    exec_code=paste0("java -Xmx",ram,"g", " -Djava.io.tmpdir=",tmp_dir," -jar ",
      bin_picard," MarkDuplicates I=",bam, " O=",out_file,
      " M=",out_file_md, remove_duplicates, " AS=true VALIDATION_STRINGENCY=LENIENT ",
      paste0("MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=",hnd)," ",tmp)
    
  if(mode=="batch"){
       out_file_dir2=set_dir(dir=out_file_dir,name="batch")
       batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
       threads=threads,output_dir=out_file_dir2)
       exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
  }
    
  if(verbose){
       print_verbose(job=job,arg=argg,exec_code=exec_code)
    }

    error=execute_job(exec_code=exec_code)
    
    if(error!=0){
      stop("markdups failed to run due to unknown error.
      Check std error for more information.")
    }

  job_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    exec_code=exec_code, 
    task_id=task_id, 
    input_args = argg,
    out_file_dir=out_file_dir,
      out_file=list(
        bam=out_file,
        log=out_file_md
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
#' @param bin_picard Path to picard executable. Default tools/picard/build/libs/picard.jar.
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

summary_metrics_bam_picard=function(
  bin_picard=build_default_tool_binary_list()$bin_picard,
  bam="",output_dir=".",verbose=FALSE,
  batch_config=build_default_preprocess_config(),tmp_dir=".",
  threads=3,ram=4,mode="local",executor_id=make_unique_id("summaryMetrics"),
  task_name="summaryMetrics",
  time="48:0:0",update_time=60,wait=FALSE,hold=NULL
){
  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir,name="summary")

  tmp=""
  if (!tmp_dir==""){
    tmp=paste0(" TMP_DIR=",tmp_dir)
  }

  out_file=paste0(out_file_dir,"/",get_file_name(bam),".picard_summary.txt")
  exec_code=paste0("java -Xmx",ram,"g", " -Djava.io.tmpdir=",tmp_dir,
        " -jar ",bin_picard," CollectAlignmentSummaryMetrics ",
        "VALIDATION_STRINGENCY=SILENT I=",bam," O=",out_file,tmp)
  
  job=build_job(executor_id=executor_id,task_id=task_id)

  if(mode=="batch"){
       out_file_dir2=set_dir(dir=out_file_dir,name="batch")
       batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
       threads=threads,output_dir=out_file_dir2)
       exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
  }

 if(verbose){
       print_verbose(job=job,arg=argg,exec_code=exec_code)
  }

  error=execute_job(exec_code=exec_code)
  if(error!=0){
    stop("picard failed to run due to unknown error.
    Check std error for more information.")
  }

  job_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    exec_code=exec_code, 
    task_id=task_id, 
    input_args = argg,
    out_file_dir=out_file_dir,
    out_file=list(
      summary=out_file
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
#' @param bin_picard Path to bwa executable. Default path tools/samtools/samtools.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param ram RAM memory to use in GB. Default 4.
#' @param tmp_dir Path to TMP directory. Default .
#' @param executor_id [OPTIONAL] Task executor name. Default "recalCovariates"
#' @param task_name [OPTIONAL] Task name. Default "recalCovariates"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export

insertsize_metrics_bam_picard=function(
  bin_picard=build_default_tool_binary_list()$bin_picard,
  bam="",output_dir=".",verbose=FALSE,
  batch_config=build_default_preprocess_config(),tmp_dir=".",threads=1,ram=4,
  mode="local",executor_id=make_unique_id("insertsizeMetrics"),
  task_name="insertsizeMetrics",
  time="48:0:0",update_time=60,wait=FALSE,hold=NULL){

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir,name="insertsize")


  tmp=""
  if (!tmp_dir==""){
    tmp=paste0(" TMP_DIR=",tmp_dir)
  }

  out_file=paste0(out_file_dir,"/",get_file_name(bam),".picard_insert_size.txt")
  out_file_pdf=paste0(out_file_dir,"/",get_file_name(bam),".picard_insert_size.pdf")
  exec_code=paste0("java -Xmx",ram,"g", " -Djava.io.tmpdir=",tmp_dir," -jar ",
      bin_picard," CollectInsertSizeMetrics ","VALIDATION_STRINGENCY=SILENT I=",
      bam," O=",out_file," H=",out_file_pdf,tmp)

  job=build_job(executor_id=executor_id,task_id=task_id)
  if(mode=="batch"){
       out_file_dir2=set_dir(dir=out_file_dir,name="batch")
       batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
       threads=threads,output_dir=out_file_dir2)
       exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
  }

  if(verbose){
       print_verbose(job=job,arg=argg,exec_code=exec_code)
    }
  error=execute_job(exec_code=exec_code)
  if(error!=0){
    stop("picard failed to run due to unknown error.
    Check std error for more information.")
  }

  job_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    exec_code=exec_code, 
    task_id=task_id, 
    input_args = argg,
    out_file_dir=out_file_dir,
      out_file=list(
        insert=out_file,
        pdf=out_file_pdf
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
#' @param bin_picard Path to bwa executable. Default path tools/samtools/samtools.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param ram RAM memory to use in GB. Default 4.
#' @param tmp_dir Path to TMP directory. Default .
#' @param bi Bait capture target interval for panel data. Interval format.
#' @param ti Primary target intervals for panel data. Interval format.
#' @param executor_id [OPTIONAL] Task executor name. Default "recalCovariates"
#' @param task_name [OPTIONAL] Task name. Default "recalCovariates"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export

tg_summary_metrics_bam_picard=function(
  bin_picard=build_default_tool_binary_list()$bin_picard,bam="",output_dir=".",
  verbose=FALSE,batch_config=build_default_preprocess_config(),
  tmp_dir=".",ref_genome="",threads=1,ram=4,bi="",ti="",
  mode="local",executor_id=make_unique_id("TGsummaryMetrics"),
  task_name="TGsummaryMetrics",time="48:0:0",update_time=60,
  wait=FALSE,hold=NULL){

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir,name="summary")


  tmp=""
  if (!tmp_dir==""){
    tmp=paste0(" TMP_DIR=",tmp_dir)
  }

  ref=""
  if(ref_genome!=""){
    ref=paste0(" R=",ref_genome)
  }
 


  out_file_ts=paste0(out_file_dir,"/",get_file_name(bam),".picard_TS.txt")
  out_file=paste0(out_file_dir,"/",get_file_name(bam),".picard_CollectHSmetrics.txt")
  exec_code=paste0("java -Xmx",ram,"g", " -Djava.io.tmpdir=",tmp_dir,
        " -jar ",bin_picard," CollectHsMetrics VALIDATION_STRINGENCY=SILENT BI=",
        bi," TI=",ti," I=",bam," THEORETICAL_SENSITIVITY_OUTPUT=",
        out_file_ts," O=",out_file,ref,tmp)


  job=build_job(executor_id=executor_id,task_id=task_id)

  if(mode=="batch"){
        out_file_dir2=set_dir(dir=out_file_dir,name="batch")
        batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
        threads=threads,output_dir=out_file_dir2)
        exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
  }


   if(verbose){
       print_verbose(job=job,arg=argg,exec_code=exec_code)
   }

  error=execute_job(exec_code=exec_code)
  if(error!=0){
    stop("picard failed to run due to unknown error.
    Check std error for more information.")
  }

  job_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    exec_code=exec_code, 
    task_id=task_id, 
    input_args = argg,
    out_file_dir=out_file_dir,
      out_file=list(
        ts=out_file_ts,
        summary=out_file
        )
  )

  
  job_report[["steps"]][["parse_summary"]]<-parse_picard_metrics(
    summary=out_file,
    output_dir=out_file_dir,
    mode=mode,executor_id=task_id,
    output_name=get_file_name(bam),
    time=time,update_time=60,
    wait=FALSE,hold=job_report$job_id
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
#' @param bin_samtools Path to picard executable. Default path tools/picard/build/libs/picard.jar.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param ram RAM memory to use in GB. Default 4.
#' @param tmp_dir Path to TMP directory. Default .
#' @param ri Ribosomal interval for RNAseq data. Interval format.
#' @param ref_flat Path to flat refrence. Interval format.
#' @param executor_id [OPTIONAL] Task executor name. Default "recalCovariates"
#' @param task_name [OPTIONAL] Task name. Default "recalCovariates"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export

rnaseq_summary_metrics_bam_picard=function(
  bin_picard=build_default_tool_binary_list()$bin_picard,
  bam="",output_dir=".",verbose=FALSE,
  batch_config=build_default_preprocess_config(),tmp_dir=".",threads=1,ram=4,
  ri="",ref_flat="",mode="local",executor_id=make_unique_id("RNAsummaryMetrics"),
  task_name="RNAsummaryMetrics",time="48:0:0",
  update_time=60,wait=FALSE,hold=NULL){


  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir,name="summary")

  tmp=""
  if (!tmp_dir==""){
    tmp=paste0(" TMP_DIR=",tmp_dir)
  }

  out_file=paste0(out_file_dir,"/",get_file_name(bam),".CollectRNAseqMetrics.txt")
  exec_code=paste0("java -Xmx",ram,"g", " -Djava.io.tmpdir=",tmp_dir," -jar ",bin_picard,
        " CollectRnaSeqMetrics VALIDATION_STRINGENCY=SILENT STRAND_SPECIFICITY='NONE' REF_FLAT=",
         ref_flat, " RIBOSOMAL_INTERVALS=",ri,
         " I=",bam," O=",out_file,tmp)
  
  
  
job=build_job(executor_id=executor_id,task_id=task_id)
  if(mode=="batch"){
       out_file_dir2=set_dir(dir=out_file_dir,name="batch")
       batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
       threads=threads,output_dir=out_file_dir2)
       exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
  }


  if(verbose){
       print_verbose(job=job,arg=argg,exec_code=exec_code)
  }

  error=execute_job(exec_code=exec_code)
  if(error!=0){
    stop("picard failed to run due to unknown error.
    Check std error for more information.")
  }

  job_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    exec_code=exec_code, 
    task_id=task_id, 
    input_args = argg,
    out_file_dir=out_file_dir,
      out_file=list(
        summary=out_file
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
#' @param bin_picard Path to picard executable. Default path tools/samtools/samtools.
#' @param ref_genome Path to reference genome.
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

wgs_summary_metrics_bam_picard=function(
  bin_picard=build_default_tool_binary_list()$bin_picard,
  bam="",mapq=0,ref_genome=build_default_reference_list()$HG19$reference$genome,
  output_dir=".",verbose=FALSE,
  batch_config=build_default_preprocess_config(),
  tmp_dir=".",threads=1,ram=4,mode="local",
  executor_id=make_unique_id("WGSsummaryMetrics"),
  task_name="WGSsummaryMetrics",time="48:0:0",
  update_time=60,wait=FALSE,hold=NULL){


  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir,name="summary")


  tmp=""
  if (!tmp_dir==""){
    tmp=paste0(" TMP_DIR=",tmp_dir)
  }

  out_file=paste0(out_file_dir,"/",get_file_name(bam),".picard_wgs_q00.txt")
  exec_code=paste0("java -Xmx",ram,"g", " -Djava.io.tmpdir=",tmp_dir," -jar ",
            bin_picard," CollectWgsMetrics VALIDATION_STRINGENCY=SILENT MINIMUM_MAPPING_QUALITY=",
            mapq," I=",bam," O=",out_file, " R=",ref_genome,tmp)
  
job=build_job(executor_id=executor_id,task_id=task_id)
  if(mode=="batch"){
       out_file_dir2=set_dir(dir=out_file_dir,name="batch")
       batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
       threads=threads,output_dir=out_file_dir2)
       exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
  }


 if(verbose){
       print_verbose(job=job,arg=argg,exec_code=exec_code)
    }
  error=execute_job(exec_code=exec_code)
  if(error!=0){
    stop("picard failed to run due to unknown error.
    Check std error for more information.")
  }


  job_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    exec_code=exec_code, 
    task_id=task_id, 
    input_args = argg,
    out_file_dir=out_file_dir,
      out_file=list(
        summary=out_file
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
#' @param bin_picard Path to bwa executable. Default path tools/samtools/samtools.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param ram RAM memory to use in GB. Default 4.
#' @param tmp_dir Path to TMP directory. Default .
#' @param executor_id [OPTIONAL] Task executor name. Default "recalCovariates"
#' @param task_name [OPTIONAL] Task name. Default "recalCovariates"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export

new_insertsize_metrics_bam_picard=function(
  bin_picard=build_default_tool_binary_list()$bin_picard,
  bam=NULL,
  deviations=NULL,
  min_width=NULL,
  width=NULL,
  ...
){



 run_main=function(
    .env
  ){
    options(scipen=999)
    .this.env=environment()
    append_env(to=.this.env,from=.env)

    set_main(.env=.this.env)

    if(!is.null(width)){
        width=paste0(" W=",width)
    }

     if(!is.null(deviations)){
        deviations=paste0(" W=",deviations)
     }

    if(!is.null(min_width)){
        min_width=paste0(" MW=",min_width)
    }
  
    .main$out_files$insert_size=paste0(out_file_dir,"/",input_id,".picard_insert_size.txt")
    .main$out_files$pdf=paste0(out_file_dir,"/",input_id,".picard_insert_size.pdf")
    .main$exec_code=paste0("java -Xmx",ram,"g", " -Djava.io.tmpdir=",tmp_dir," -jar ",
        bin_picard," CollectInsertSizeMetrics ","VALIDATION_STRINGENCY=SILENT I=",
        input," O=",.main$out_files$insert_size," H=",.main$out_files$pdf,
      " TMP_DIR=",tmp_dir, " DEVIATIONS=",deviations,width,min_width)

 

    run_job(.env=.this.env)

    .main.step=.main$steps[[fn_id]]

    .env$.main <- .main
  }

  .base.env=environment()
  list2env(list(...),envir=.base.env)
  set_env_vars(
    .env= .base.env,
    vars="bam"
  )
 
   launch(.env=.base.env)


}


#' Extract Histogram and Metric Information from Picard Summary File
#'
#'
#' @param summary Path to Picard summary file
#' @export

parse_summary_picard=function(
  summary=NULL,
  ...
  ){

  run_main=function(
    .env
  ){
    options(scipen=999)
    .this.env=environment()
    append_env(to=.this.env,from=.env)

    set_main(.env=.this.env)

    .main$out_files$histogram=paste0(out_file_dir,"/",input_id,".histogram.txt")
    .main$out_files$metrics=paste0(out_file_dir,"/",input_id,".metrics.txt")
    .main$exec_code=paste0(
      "awk -v RS= \"{if(NR==2) {out=\\\"",.main$out_files$metrics,
      "\\\"} else if(NR==3){ out=\\\"",
    .main$out_files$histogram,
    "\\\"}; if (NR>1)  print > (out)}\" ",input
    )

    run_job(.env=.this.env)

    .main.step=.main$steps[[fn_id]]

    .env$.main <- .main
  }
  
  .base.env=environment()
  list2env(list(...),envir=.base.env)
  set_env_vars(
    .env= .base.env,
    vars="summary"
  )

  launch(.env=.base.env)

}

