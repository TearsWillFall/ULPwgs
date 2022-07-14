#' Sort and index a sequence file
#' 
#'
#' Wrapper around index_bam_samtools and sort_bam_samtools functions
#' 
#' 
#' @param bam Path to the input file with the sequence.
#' @param bin_path Path to bwa executable. Default path tools/samtools/samtools.
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

sort_and_index_bam_samtools=function(bin_path="tools/samtools/samtools",bam="",output_dir="",
verbose=FALSE,threads=3,ram=1,sort=TRUE,coord_sort=TRUE,index=TRUE,stats="all", clean=FALSE,
mode="local",executor_id=make_unique_id("sortANDindex"),task_name="sortANDindex",time="48:0:0",
update_time=60,wait=FALSE,hold=""){

  argg <- as.list(environment())

  task_id=make_unique_id(task_name)

  out_file_dir=set_dir(dir=output_dir)

  job=build_job(executor_id=executor_id,task=task_id)

  job_report=build_job_report(
    job_id=job,
    executor_id=executor_id, 
    task_id=task_id,
    input_args=argg,
    out_file_dir=out_file_dir,
      out_files=list()
  )

  if(sort){
      job_report[["steps"]][["sort"]]=sort_bam_samtools(bin_path=bin_path,bam=bam,output_dir=out_file_dir,
      ram=ram,verbose=verbose,threads=threads,coord_sort=coord_sort,clean=clean,
      executor_id=task_id,mode=mode,time=time,
      update_time=update_time,wait=FALSE,hold=hold)

      out_file_dir=set_dir(dir=output_dir,name="sorted")
  
      if (coord_sort){
      
        bam=paste0(out_file_dir,"/",get_file_name(bam),".sorted.",get_file_ext(bam))

        if(index){
            job_report[["steps"]][["index"]]=index_bam_samtools(bin_path=bin_path,
            bam=bam,verbose=verbose,threads=threads,ram=ram,
            executor_id=task_id,mode=mode,time=time,
            update_time=update_time,wait=FALSE,hold=job,output_dir = out_file_dir)
          if(stats=="index"|stats=="all"){
              job_report[["stats"]]=stats_bam_samtools(bin_path=bin_path,bam=bam,output_dir=out_file_dir,
              verbose=verbose,threads=threads,stats="index",executor_id=task_id,
              mode=mode,time=time,update_time=update_time,wait=FALSE,hold=job)
          }
        }
      }
  }else{
     job_report[["steps"]][["index"]]=index_bam_samtools(bin_path=bin_path,bam=bam,verbose=verbose,threads=threads,
     executor_id=task_id,mode=mode,time=time,update_time=update_time,wait=FALSE,hold=hold)
  }


  if(stats=="flag"|stats=="all"){
    job_report[["steps"]][["stats"]]=append(job_report[["steps"]][["stats"]],
    stats_bam_samtools(bin_path=bin_path,bam=bam,output_dir=out_file_dir,
      verbose=verbose,threads=threads,stats="flag",executor_id=task_id,
      mode=mode,time=time,update_time=update_time,
      wait=FALSE,hold=hold))
  }

  if(wait&&mode=="batch"){
      job_validator(job=c(
      job_report[["steps"]][["sort"]]$job_id,
      job_report[["steps"]][["index"]]$job_id,
      job_report[["steps"]][["stats"]][["steps"]][["stats"]]$job_id,
      job_report[["steps"]][["stats"]][["steps"]][["index"]]$job_id),
      time=update_time,verbose=verbose,threads=threads)
  }

  return(job_report)
}


#' Sort a BAM file
#'
#' This function sorts a genome sequence file (BAM/SAM)
#'
#' @param bam Path to the input file with the sequence.
#' @param bin_path Path to bwa executable. Default path tools/samtools/samtools.
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

sort_bam_samtools=function(bin_path="tools/samtools/samtools",bam="",output_dir="",
verbose=FALSE,threads=3,ram=1,coord_sort=TRUE,mode="local",
executor_id=make_unique_id("sortBAM"),clean=FALSE,task_name="sortBAM",
time="48:0:0",update_time=60,wait=FALSE,hold=""){

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir,name="sorted")

  sort_type=""
  
  out_file=paste0(out_file_dir,"/",get_file_name(bam),".sorted.",get_file_ext(bam))
  if(!coord_sort){
    sort_type=" -n "
  }
  exec_code=paste0(bin_path," sort ",sort_type, bam," -@ ",threads," -m ",ram,"G"," -o ",
  out_file)

  if(clean){
    exec_code=paste(exec_code," && rm",paste(bam,collapse=" "))
  }

  job=build_job(executor_id=executor_id,task=task_id)

  if(mode=="batch"){
       out_file_dir2=set_dir(dir=out_file_dir,name="batch")
       exec_batch=build_job_exec(job=job,time=time,ram=ram,threads=threads,
       output_dir=out_file_dir2,hold=hold)
       exec_code=paste("echo 'source ~/.bashrc;",exec_code,"'|",exec_batch)
  }

    
  if (verbose){
    print_verbose(job=job,arg=argg,exec_code=exec_code)
  }
  error=system(exec_code)
  if(error!=0){
    stop("samtools failed to run due to unknown error.
    Check std error for more information.")
  }

  job_report=build_job_report(
    job_id=job, 
    executor_id=executor_id, 
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


index_bam_samtools=function(bin_path="tools/samtools/samtools",
bam="",verbose=FALSE,threads=3,ram=4,
mode="local",executor_id=make_unique_id("indexBAM"),
task_name="indexBAM",time="48:0:0",update_time=60, output_dir="",
wait=FALSE,hold=""){

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir)
  exec_code=paste(bin_path," index",bam," -@ ",threads)
  job=build_job(executor_id=executor_id,task=task_id)

  if(mode=="batch"){
       out_file_dir2=set_dir(dir=out_file_dir,name="batch")
       exec_batch=build_job_exec(job=job,time=time,ram=ram,threads=threads,
       output_dir=out_file_dir2,hold=hold)
       exec_code=paste("echo 'source ~/.bashrc;",exec_code,"'|",exec_batch)
  }

    
  if (verbose){
    print_verbose(job=job,arg=argg,exec_code=exec_code)
  }

  error=system(exec_code)
  if(error!=0){
    stop("samtools failed to run due to unknown error.
    Check std error for more information.")
  }

  
  job_report=build_job_report(
    job_id=job,
    executor_id=executor_id, 
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
#' @param bin_path Path to bwa executable. Default path tools/samtools/samtools.
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


stats_bam_samtools=function(bin_path="tools/samtools/samtools",bam="",output_dir="",
verbose=FALSE,threads=3,ram=4,stats="all",mode="local",executor_id=make_unique_id("statsBAM"),
task_name="statsBAM",time="48:0:0",update_time=60,wait=FALSE,hold=""){

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir,name="stats")

  job=build_job(executor_id=executor_id,task=task_id)

  job_report=build_job_report(
    job_id=job, 
    executor_id=executor_id, 
    task_id=task_id,
    input_args=argg,
    out_file_dir=out_file_dir,
    out_files=list()
  )


  if(stats=="all"|stats=="flag"){
    job_report[["steps"]][["flag"]]=stats_flag_samtools(bin_path=bin_path,bam=bam,output_dir=out_file_dir,
    verbose=verbose,threads=threads,mode=mode,time=time,update_time=update_time,wait=FALSE,hold=hold)

  }

  if(stats=="all"|stats=="index"){
    job_report[["steps"]][["index"]]=stats_index_samtools(bin_path=bin_path,bam=bam,output_dir=out_file_dir,
    verbose=verbose,threads=threads,mode=mode,time=time,update_time=update_time,wait=FALSE,hold=hold)
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
#' @param bin_path Path to bwa executable. Default path tools/samtools/samtools.
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



stats_flag_samtools=function(bin_path="tools/samtools/samtools",bam="",output_dir="",
verbose=FALSE,threads=3,ram=4,mode="local",executor_id=make_unique_id("statsFlag"),
task_name="statsFlag",time="48:0:0",update_time=60,wait=FALSE,hold=""){

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir,name="flag")

  out_file=paste0(out_file_dir,"/",get_file_name(bam),".flagstat.txt")
  exec_code=paste0(bin_path," flagstat ",bam," -@ ",threads," > ",out_file)

  job=build_job(executor_id=executor_id,task_id=task_id)
  if(mode=="batch"){
       out_file_dir2=set_dir(dir=out_file_dir,name="batch")
       exec_batch=build_job_exec(job=job,time=time,ram=ram,threads=threads,
       output_dir=out_file_dir2,hold=hold)
       exec_code=paste("echo 'source ~/.bashrc;",exec_code,"'|",exec_batch)
  }

    
  
  if (verbose){
    print_verbose(job=job,arg=argg,exec_code=exec_code)
  }


  error=system(exec_code)
  if(error!=0){
    stop("samtools failed to run due to unknown error.
    Check std error for more information.")
  }


  job_report=build_job_report(
    job_id=job,
    executor_id=executor_id, 
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
#' @param bin_path Path to bwa executable. Default path tools/samtools/samtools.
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


stats_index_samtools=function(bin_path="tools/samtools/samtools",bam="",output_dir="",
verbose=FALSE,threads=3,ram=4,mode="local",executor_id=make_unique_id("statsINDEX"),
task_name="statsINDEX",time="48:0:0",update_time=60,wait=FALSE,hold=""){

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir,name="index")

  out_file=paste0(out_file_dir,"/",get_file_name(bam),".idxstats.txt")
  exec_code=paste0(bin_path," idxstats ",bam," > ",out_file)


  job=build_job(executor_id=executor_id,task_id=task_id)
  if(mode=="batch"){
       out_file_dir2=set_dir(dir=out_file_dir,name="batch")
       exec_batch=build_job_exec(job=job,time=time,ram=ram,threads=threads,
       output_dir=out_file_dir2,hold=hold)
       exec_code=paste("echo 'source ~/.bashrc;",exec_code,"'|",exec_batch)
  }

    
  if (verbose){
    print_verbose(job=job,arg=argg,exec_code=exec_code)
  }
     error=system(exec_code)
  if(error!=0){
    stop("samtools failed to run due to unknown error.
    Check std error for more information.")
  }

  
  job_report=build_job_report(
    job_id=job,
    executor_id=executor_id, 
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
#' @param bin_path Path to bwa executable. Default path tools/samtools/samtools.
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

mapq_metrics_bam_samtools=function(bin_path="tools/samtools/samtools",bam="",output_dir="",
verbose=FALSE,threads=3,ram=4,mode="local",executor_id=make_unique_id("metricsMAPQ"),
task_name="metricsMAPQ",time="48:0:0",update_time=60,wait=FALSE,hold=""){

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir,name="mapq")


  out_file=paste0(out_file_dir,"/",get_file_name(bam),".mapq_dist.txt")
  exec_code=paste(bin_path,"view",bam," -@ ",threads, " | awk -F", "'\\t'",
    "'{c[$5]++} END { for (i in c) catf(\"%s\\t%s\\n\",i,c[i]) }'",
    " | sort -t$'\\t' -k 1 -g >", out_file)

  job=build_job(executor_id=executor_id,task_id=task_id)

  if(mode=="batch"){
       out_file_dir2=set_dir(dir=out_file_dir,name="batch")
       exec_batch=build_job_exec(job=job,time=time,ram=ram,threads=threads,
       output_dir=out_file_dir2,hold=hold)
       exec_code=paste("echo 'source ~/.bashrc;",exec_code,"'|",exec_batch)
  }

  
  if (verbose){
    print_verbose(exec_code=exec_code)
  }

    error=system(exec_code)
  if(error!=0){
    stop("samtools failed to run due to unknown error.
    Check std error for more information.")
  }

  job_report=build_job_report(
    job_id=job,
    executor_id=executor_id, 
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


