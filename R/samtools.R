#' Sort and index a sequence file
#' 
#'
#' Wrapper around index_bam_samtools and sort_bam_samtools functions
#' 
#' 
#' @param bam Path to the input file with the sequence.
#' @param bin_samtools Path to samtools executable. Default path tools/samtools/samtools.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param threads Number of threads. Default 3
#' @param ram Ram memory to use per thread in GB. Default 1GB
#' @param sort Sort BAM file. Default TRUE.
#' @param coord_sort Generate a coord sorted file. Otherwise queryname sorted. Default TRUE
#' @param index Generate an index file for sorted BAM. Default TRUE
#' @param clean Remove input files. Default FALSE
#' @param stats Generate BAM stats. Default all. Options ["all","flag","index",""]
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Task EXECUTOR id. Default "mardupsGATK"
#' @param task_name Name of the task. Default "mardupsGATK"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] Hold job until job is finished. Job ID. 
#' @export

sort_and_index_bam_samtools=function(
  bin_samtools=build_default_tool_binary_list()$bin_samtools,
  bam="",output_dir=".",verbose=FALSE,
  batch_config=build_default_preprocess_config(),threads=3,ram=1,sort=TRUE,
  coord_sort=TRUE,index=TRUE,stats="all", clean=FALSE,
  mode="local",executor_id=make_unique_id("sortANDindex"),
  task_name="sortANDindex",time="48:0:0",
  update_time=60,wait=FALSE,hold=NULL
){

  argg <- as.list(environment())

  task_id=make_unique_id(task_name)

  out_file_dir=set_dir(dir=output_dir)

  job=build_job(executor_id=executor_id,task=task_id)

  job_report=build_job_report(
    job_id=job,
    executor_id=executor_id, 
    exec_code=list(), 
    task_id=task_id,
    input_args=argg,
    out_file_dir=out_file_dir,
      out_files=list()
  )

  if(sort){
      job_report[["steps"]][["sort"]]=sort_bam_samtools(
        bin_samtools=bin_samtools,bam=bam,output_dir=out_file_dir,
        ram=ram,verbose=verbose,threads=threads,coord_sort=coord_sort,clean=clean,
        executor_id=task_id,mode=mode,time=time,batch_config = batch_config,
        update_time=update_time,wait=FALSE,hold=hold)

      out_file_dir=set_dir(dir=output_dir,name="sorted")
      bam=job_report[["steps"]][["sort"]]$out_files$bam
      hold=job_report[["steps"]][["sort"]]$job_id

      if (coord_sort){
        if(index){
            job_report[["steps"]][["index"]] <-index_bam_samtools(
              bin_samtools=bin_samtools,
              bam=bam,verbose=verbose,threads=threads,ram=ram,
              executor_id=task_id,mode=mode,time=time,batch_config = batch_config,
              update_time=update_time,wait=FALSE,hold=hold,
              output_dir = out_file_dir)

          if(stats=="index"|stats=="all"){
              job_report[["steps"]][["index_stats"]]<- stats_bam_samtools(
                bin_samtools=bin_samtools,bam=bam,
                output_dir=out_file_dir,batch_config = batch_config,
                verbose=verbose,threads=threads,stats="index",executor_id=task_id,
                mode=mode,time=time,update_time=update_time,wait=FALSE,hold=hold)
          }
        }
      }
  }else{
     job_report[["steps"]][["index"]] <- index_bam_samtools(
      bin_samtools=bin_samtools,bam=bam,verbose=verbose,threads=threads,
      executor_id=task_id,mode=mode,time=time,batch_config = batch_config,
      update_time=update_time,wait=FALSE,hold=hold)
  }

  if(stats=="flag"|stats=="all"){
    job_report[["steps"]][["flag_stats"]] <- stats_bam_samtools(
      bin_samtools=bin_samtools,bam=bam,output_dir=out_file_dir,
      verbose=verbose,threads=threads,stats="flag",executor_id=task_id,
      mode=mode,time=time,update_time=update_time,batch_config = batch_config,
      wait=FALSE,hold=hold)
  }

  if(wait&&mode=="batch"){
      job_validator(job=unlist_lv(job_report[["steps"]],var="job_id"),
      time=update_time,verbose=verbose,threads=threads)
  }

  return(job_report)
}


#' Sort a BAM file
#'
#' This function sorts a genome sequence file (BAM/SAM)
#'
#' @param bam Path to the input file with the sequence.
#' @param bin_samtols Path to samtools executable. Default path tools/samtools/samtools.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param threads Number of threads. Default 3
#' @param ram Ram memory to use per thread in GB. Default 1GB
#' @param output_dir Path to the output directory.
#' @param coord_sort Generate a coord sorted file. Otherwise queryname sorted. Default TRUE
#' @param clean Remove input files. Default FALSE
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Task EXECUTOR ID. Default "mardupsGATK"
#' @param task_name Name of the task. Default "mardupsGATK"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] Hold job until job is finished. Job ID. 
#' @export

sort_bam_samtools=function(
  bin_samtools=build_default_tool_binary_list()$bin_samtools,
  bam="",output_dir=".",verbose=FALSE,
  batch_config=build_default_preprocess_config(),
  threads=3,ram=1,coord_sort=TRUE,mode="local",
  executor_id=make_unique_id("sortBAM"),clean=FALSE,task_name="sortBAM",
  time="48:0:0",update_time=60,wait=FALSE,hold=NULL){

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir,name="sorted")

  sort_type=""
  
  out_file=paste0(out_file_dir,"/",get_file_name(bam),".sorted.",get_file_ext(bam))
  if(!coord_sort){
    sort_type=" -n "
  }
  exec_code=paste0(bin_samtools," sort ",sort_type, bam," -@ ",threads," -m ",ram,"G"," -o ",
  out_file)

  if(clean){
    exec_code=paste(exec_code," && rm",paste(bam,collapse=" "))
  }

  job=build_job(executor_id=executor_id,task=task_id)

  if(mode=="batch"){
       out_file_dir2=set_dir(dir=out_file_dir,name="batch")
       batch_code=build_job_exec(job=job,time=time,ram=ram,threads=threads,
       output_dir=out_file_dir2,hold=hold)
       exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
  }

    
  if (verbose){
    print_verbose(job=job,arg=argg,exec_code=exec_code)
  }
  error=execute_job(exec_code=exec_code)
  if(error!=0){
    stop("samtools failed to run due to unknown error.
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
      bam=out_file)
    )

  if(wait&&mode=="batch"){
      job_validator(job=job_report$job_id,
      time=update_time,verbose=verbose,threads=threads)
  }

  return(job_report)
}


index_bam_samtools=function(
  bin_samtools=build_default_tool_binary_list()$bin_samtools,
  bam="",verbose=FALSE,batch_config=build_default_preprocess_config(),
  threads=3,ram=4,mode="local",executor_id=make_unique_id("indexBAM"),
  task_name="indexBAM",time="48:0:0",update_time=60, output_dir=".",
  wait=FALSE,hold=NULL
){

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir)
  exec_code=paste(bin_samtools," index",bam," -@ ",threads)
  job=build_job(executor_id=executor_id,task=task_id)

  if(mode=="batch"){
       out_file_dir2=set_dir(dir=out_file_dir,name="batch")
       batch_code=build_job_exec(job=job,time=time,ram=ram,threads=threads,
       output_dir=out_file_dir2,hold=hold)
       exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
  }

    
  if (verbose){
    print_verbose(job=job,arg=argg,exec_code=exec_code)
  }

  error=execute_job(exec_code=exec_code)
  if(error!=0){
    stop("samtools failed to run due to unknown error.
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
      bai=paste0(bam,".bai"))
      )


  if(wait&&mode=="batch"){
      job_validator(job=job_report,
      time=update_time,verbose=verbose,threads=threads)
  }

  return(job_report)
}






#' Generate BAM file flag and index stats
#'
#'
#' @param bam Path to the input file with the sequence.
#' @param bin_samtools Path to bwa executable. Default path tools/samtools/samtools.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param threads Number of threads. Default 3
#' @param ram RAM per thread to use. Default 4.
#' @param stats Generate BAM stats. Default all. Options ["all","flag","index",""]
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Task EXECUTOR ID. Default "mardupsGATK"
#' @param task_name Name of the task. Default "mardupsGATK"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] Hold job until job is finished. Job ID.
#' @export


stats_bam_samtools=function(
  bin_samtools=build_default_tool_binary_list()$bin_samtools,
  bam="",output_dir=".",verbose=FALSE,
  batch_config=build_default_preprocess_config(),
  threads=3,ram=4,stats="all",
  mode="local",executor_id=make_unique_id("statsBAM"),
  task_name="statsBAM",time="48:0:0",
  update_time=60,wait=FALSE,hold=NULL
){

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir,name="stats")

  job=build_job(executor_id=executor_id,task=task_id)

  job_report=build_job_report(
    job_id=job, 
    executor_id=executor_id,
    exec_code=list(),  
    task_id=task_id,
    input_args=argg,
    out_file_dir=out_file_dir,
    out_files=list()
  )


  if(stats=="all"|stats=="flag"){
    job_report[["steps"]][["flag"]]=stats_flag_samtools(
      bin_samtools=bin_samtools,bam=bam,output_dir=out_file_dir,
      verbose=verbose,threads=threads,mode=mode,time=time,
      update_time=update_time,wait=FALSE,hold=hold)

  }

  if(stats=="all"|stats=="index"){
    job_report[["steps"]][["index"]]=stats_index_samtools(
      bin_samtools=bin_samtools,
      bam=bam,output_dir=out_file_dir,
      verbose=verbose,threads=threads,mode=mode,time=time,
      update_time=update_time,wait=FALSE,hold=hold)
  }


  if(wait&&mode=="batch"){
      job_validator(job=c(
      job_report[["steps"]][["index"]]$job_id,
      job_report[["steps"]][["flag"]]$job_id),
      time=update_time,verbose=verbose,threads=threads)
  }

  return(job_report)
}

#' Generate BAM file flagstats
#'
#'
#' @param bam Path to the input file with the sequence.
#' @param bin_samtools Path to samtools executable. Default path tools/samtools/samtools.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param threads Number of threads. Default 3
#' @param ram RAM per thread to use. Default 4.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor Job EXECUTOR ID. Default "mardupsGATK"
#' @param task_name Name of the task. Default "mardupsGATK"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] Hold job until job is finished. Job ID. 
#' @export



stats_flag_samtools=function(
  bin_samtools=build_default_tool_binary_list()$bin_samtools,
  bam="",output_dir=".",verbose=FALSE,
  batch_config=build_default_preprocess_config(),
  threads=3,ram=4,mode="local",executor_id=make_unique_id("statsFlag"),
  task_name="statsFlag",time="48:0:0",update_time=60,
  wait=FALSE,hold=NULL
){

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir,name="flag")

  out_file=paste0(out_file_dir,"/",get_file_name(bam),".flagstat.txt")
  exec_code=paste0(bin_samtools," flagstat ",bam," -@ ",threads," > ",out_file)

  job=build_job(executor_id=executor_id,task_id=task_id)
  if(mode=="batch"){
       out_file_dir2=set_dir(dir=out_file_dir,name="batch")
       batch_code=build_job_exec(job=job,time=time,ram=ram,threads=threads,
       output_dir=out_file_dir2,hold=hold)
       exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
  }

    
  
  if (verbose){
    print_verbose(job=job,arg=argg,exec_code=exec_code)
  }


  error=execute_job(exec_code=exec_code)
  if(error!=0){
    stop("samtools failed to run due to unknown error.
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
      flag_stat=out_file)
    )

  if(wait&&mode=="batch"){
      job_validator(job=job_report$job_id,
      time=update_time,verbose=verbose,threads=threads)
  }

  return(job_report)
}


#' Generate BAM file indexstats
#'
#'
#' @param bam Path to the input file with the sequence.
#' @param bin_samtools Path to samtools executable. Default path tools/samtools/samtools.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param threads Number of threads. Default 3
#' @param ram RAM per thread to use. Default 4.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Job EXECUTOR ID. Default "mardupsGATK"
#' @param task_name Name of the task. Default "mardupsGATK"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] Hold job until job is finished. Job ID. 
#' @export


stats_index_samtools=function(
  bin_samtools=build_default_tool_binary_list()$bin_samtools,
  bam="",output_dir=".",verbose=FALSE,
  batch_config=build_default_preprocess_config(),threads=3,ram=4,mode="local",
  executor_id=make_unique_id("statsINDEX"),task_name="statsINDEX",
  time="48:0:0",update_time=60,wait=FALSE,hold=NULL
){

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir,name="index")

  out_file=paste0(out_file_dir,"/",get_file_name(bam),".idxstats.txt")
  exec_code=paste0(bin_samtools," idxstats ",bam," > ",out_file)


  job=build_job(executor_id=executor_id,task_id=task_id)
  if(mode=="batch"){
       out_file_dir2=set_dir(dir=out_file_dir,name="batch")
       batch_code=build_job_exec(job=job,time=time,ram=ram,threads=threads,
       output_dir=out_file_dir2,hold=hold)
       exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
  }

    
  if (verbose){
    print_verbose(job=job,arg=argg,exec_code=exec_code)
  }
     error=execute_job(exec_code=exec_code)
  if(error!=0){
    stop("samtools failed to run due to unknown error.
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
      idx_stat=out_file)
    )

  if(wait&&mode=="batch"){
      job_validator(job=job_report$job_id,
      time=update_time,verbose=verbose,threads=threads)
  }

  return(job_report)

}


#' Generate BAM MapQ metrics
#'
#'
#' @param bam Path to the input file with the sequence.
#' @param bin_samtools Path to samtools executable. Default path tools/samtools/samtools.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param threads Number of threads. Default 3
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Job EXECUTOR ID. Default "mardupsGATK"
#' @param task_name Name of the task. Default "mardupsGATK"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] Hold job until job is finished. Job ID. 
#' @export

mapq_metrics_bam_samtools=function(
  bin_samtools=build_default_tool_binary_list()$bin_samtools,
  bam="",
  output_dir=".",
  verbose=FALSE,batch_config=build_default_preprocess_config(),
  threads=3,ram=4,mode="local",
  executor_id=make_unique_id("metricsMAPQ"),
  task_name="metricsMAPQ",time="48:0:0",update_time=60,wait=FALSE,hold=NULL
){

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir,name="mapq")


  out_file=paste0(out_file_dir,"/",get_file_name(bam),".mapq_dist.txt")
  exec_code=paste(bin_samtools,"view",bam," -@ ",threads, 
  " | awk \"{c[\\$5]++} END { for (i in c) printf(\\\"%s\\t%s\\n\\\",i,c[i])}\"",
  " | sort -k 1 -g >", out_file)

  job=build_job(executor_id=executor_id,task_id=task_id)

  if(mode=="batch"){
       out_file_dir2=set_dir(dir=out_file_dir,name="batch")
       batch_code=build_job_exec(job=job,time=time,ram=ram,threads=threads,
       output_dir=out_file_dir2,hold=hold)
       exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
  }

  
  if (verbose){
    print_verbose(job=job,arg=argg,exec_code=exec_code)
  }

  error=execute_job(exec_code=exec_code)
  if(error!=0){
    stop("samtools failed to run due to unknown error.
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
      mapq_metric=out_file)
    )

  
  if(wait&&mode=="batch"){
      job_validator(job=job_report$job_id,
      time=update_time,verbose=verbose,threads=threads)
  }

  return(job_report)
}




#' Filter BAM file by size using samtools
#'
#'
#' @param bam Path to the input file with the sequence.
#' @param bin_samtools Path to samtools executable. Default path tools/samtools/samtools.
#' @param region Genomic region to search
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param threads Number of threads. Default 3
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Job EXECUTOR ID. Default "mardupsGATK"
#' @param task_name Name of the task. Default "mardupsGATK"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] Hold job until job is finished. Job ID. 
#' @export


filter_bam_by_size_samtools=function(
  rdata=NULL,
  selected=NULL,
  bin_samtools=build_default_tool_binary_list()$bin_samtools,
  bam="",
  region=NULL,
  output_name="",
  min_frag_size=0,
  max_frag_size=180,
  output_name="",
  output_dir=".",
  include=TRUE, 
  verbose=FALSE,
  batch_config=build_default_preprocess_config(),
  threads=1,ram=4,mode="local",
  executor_id=make_unique_id("filterBAMbyInsertSize"),
  task_name="filterBAMbyInsertSize",time="48:0:0",
  update_time=60,wait=FALSE,hold=NULL
){

  if(!is.null(rdata)){
    load(rdata)
    if(!is.null(selected)){
      region=region_list[selected]
    }
  }

      
    id=""
    if(output_name!=""){
      id=output_name
    }else{
      id=get_file_name(bam)
    }


  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir)
  job=build_job(executor_id=executor_id,task_id=task_id)

  position=""
  if(!is.null(region)){
    position=strsplit(region,split="__")[[1]][2]
    out_file=paste0(out_file_dir,"/",id,".",min_frag_size,"_",max_frag_size,".",region,".bam")
  }else{
    out_file=paste0(out_file_dir,"/",id,".",min_frag_size,"_",max_frag_size,".bam")
  }


  if(include){
    exec_code=paste(bin_path,"view -h ",bam,," | \ awk 'substr($0,1,1)==\"@\""," || ($9>=",
    min_frag_size,"&& $9<=",max_frag_size,") ||", "($9<=-",min_frag_size,"&& $9>=-",max_frag_size,
    ")'|",bin_path, "view -b >",out_file)
  }else{
    exec_code=paste(bin_path,"view -h ",bam,," | \ awk 'substr($0,1,1)==\"@\""," || ($9=<",
    min_frag_size,"&& $9>=",max_frag_size,") ||", "($9>=-",min_frag_size,"&& $9<=-",max_frag_size,
    ")'|",bin_path, "view -b >",out_file)

  }




  job=build_job(executor_id=executor_id,task_id=task_id)
  if(mode=="batch"){
      out_file_dir2=set_dir(dir=out_file_dir,name="batch")
      batch_code=build_job_exec(job=job,time=time,ram=ram,threads=threads,
      output_dir=out_file_dir2,hold=hold)
      exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
  } 

  if(verbose){
    print_verbose(job=job,arg=argg,exec_code=exec_code)
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
        frag_bam=out_file
    )
  )


  if(wait&&mode=="batch"){
      job_validator(job=job_report$job_id,
      time=update_time,verbose=verbose,threads=threads)
  }

  return(job_report)

}


#' Filter BAM file by size using samtools in parallel
#'
#'
#' @param bam Path to the input file with the sequence.
#' @param bin_samtools Path to samtools executable. Default path tools/samtools/samtools.
#' @param region Genomic region to search
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param threads Number of threads. Default 3
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Job EXECUTOR ID. Default "mardupsGATK"
#' @param task_name Name of the task. Default "mardupsGATK"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] Hold job until job is finished. Job ID. 
#' @export




parallel_region_filter_bam_by_size_samtools=function(
  bin_samtools=build_default_tool_binary_list()$bin_samtools,
  bin_picard=build_default_tool_binary_list()$bin_picard,
  bam="",bed=NULL,
  min_frag_size=0,
  max_frag_size=180,
  sep="\t",
  header=FALSE,
  verbose=FALSE,
  output_dir=".",
  include=TRUE,
  index=TRUE,
  clean=TRUE,
  batch_config=build_default_preprocess_config(),
  threads=3,ram=4,mode="local",
  executor_id=make_unique_id("parRegionfilterBAMbyInsertSize"),
  task_name="parRegionfilterBAMbyInsertSize",time="48:0:0",
  update_time=60,wait=FALSE,hold=NULL

){

  options(scipen=999)

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir)
  out_file_dir_tmp=set_dir(dir=out_file_dir,name="tmp")
  job=build_job(executor_id=executor_id,task_id=task_id)

  
  jobs_report=build_job_report(
    job_id=job,
    executor_id=list(),
    exec_code=list(), 
    task_id=task_id,
    input_args=argg,
    out_file_dir=out_file_dir,
    out_files=list(
      )
    )

  jobs_report[["steps"]][["getChr"]] <- get_bam_reference_chr(
    bin_samtools=bin_samtools,
    bam=bam,verbose=verbose,output_dir=out_file_dir_tmp,
    executor_id=task_id,mode="local",threads=threads,ram=ram,
    time=time,update_time=update_time,wait=FALSE,hold=hold)

  bam_chr=read.table(jobs_report[["steps"]][["getChr"]]$out_files$ref,
  sep="\t",header=TRUE)
  bam_chr$start=bam_chr$start+1
  bam_chr$order=as.numeric(as.factor(bam_chr$chr))

  if(!is.null(bed)){
      regions=read.table(bed,sep=sep,header=header)
      names(regions)[1:3]=c("chr","start","end")
      miss_chr=setdiff(unique(regions[1,]), bam_chr)
    
      ### Validate chromosome format in BED file

      if(length(miss_chr)>0){
        stop(paste0(paste(miss_chr,collapse=",")," are not present in the BAM file."))
      }

      regions=dplyr::left_join(regions,bam_chr[,c("chr","order")])
      regions=regions %>% arrange(order,start,end)
  }else{
      regions=bam_chr
  }

  regions$pos=1:nrow(regions)
  regions=regions %>% dplyr::mutate(region=paste0(pos,"__",chr,":",start,"-",end))
    
  region_list=regions$region
  names(region_list)=regions$region
  if(mode=="local"){
        jobs_report[["steps"]][["par_sample_fragment_length"]]<-
        parallel::mclapply(region_list,FUN=function(region){
        job_report <-filter_bam_by_size_samtools(
          bin_samtools=bin_samtools,
          bam=bam,
          region=region,
          min_frag_size=min_frag_size,
          max_frag_size=max_frag_size,
          include=include, 
          verbose=verbose,
          output_name=get_file_name(bam),
          output_dir=out_file_dir_tmp,
          batch_config=batch_config,
          threads=1,ram=ram,mode="local",
          executor_id=task_id,
          time=time,
          hold=hold)
    },mc.cores=threads)
    
  }else if(mode=="batch"){
          rdata_file=paste0(out_file_dir,"/",job,".regions.RData")
          output_dir=out_file_dir_tmp
          save(
            region_list,
            bam,
            min_frag_size,
            max_frag_size,
            output_dir,
            verbose,
            exclude,
            file = rdata_file)
          exec_code=paste0("Rscript -e \"ULPwgs::filter_bam_by_size_samtools(rdata=\\\"",
          rdata_file,"\\\",selected=$SGE_TASK_ID)\"")
          out_file_dir2=set_dir(dir=out_file_dir,name="batch")
          batch_code=build_job_exec(job=job,time=time,ram=ram,
          threads=1,output_dir=out_file_dir2,
          hold=hold,array=length(region_list))
          exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)

          if(verbose){
              print_verbose(job=job,arg=argg,exec_code=exec_code)
          }
          error=execute_job(exec_code=exec_code)
          if(error!=0){
              stop("cnvkit failed to run due to unknown error.
              Check std error for more information.")
          }
         
         jobs_report[["steps"]][["par_region_ichor_capture"]]<- build_job_report(
              job_id=job,
              executor_id=executor_id,
              exec_code=exec_code, 
              task_id=task_id,
              input_args=argg,
              out_file_dir=out_file_dir,
              out_files=list(
                  frag_bam=paste0(out_file_dir,"/",get_file_name(bam),".",
                  min_frag_size,"_",max_frag_size,".",names(region_list),".bam")
              )
        )
    }

    jobs_report[["steps"]][["gather_bam"]]<-gather_bam_files(
      bin_picard=bin_picard,
      bam=unlist_lvl(jobs_report[["steps"]][["par_region_fragment_length"]],var="frag_bam"),
      output_dir=out_file_dir,
      output_name=paste0(get_file_name(bam),".",min_frag_size,"_",max_frag_size,".bam"),
      executor_id=task_id,mode=mode,time=time,threads=threads,ram=ram,
      update_time=update_time,wait=FALSE,
      clean=clean,
      hold=unlist_lvl(jobs_report[["steps"]][["par_region_fragment_length"]],var="job_id",recursive=TRUE)
    )

    if(index){
        jobs_report[["steps"]][["sort_and_index"]] <- sort_and_index_bam_samtools(
        bin_samtools=bin_samtools,
        bam=jobs_report[["steps"]][["gather_bam"]]$out_files$bam,
        output_dir=out_file_dir,batch_config=batch_config,
        ram=ram,verbose=verbose,threads=threads,sort=FALSE,
        stats="",index=TRUE,clean=clean,
        mode=mode,executor_id=task_id,time=time,
        update_time=update_time,wait=FALSE,
        hold=unlist_lvl(jobs_report[["steps"]][["gather_bam"]],var="job_id")
      )
  
    }
    
    return(jobs_report)

}



