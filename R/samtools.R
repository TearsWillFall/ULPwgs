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
#' @param executor Name of the executor. Default "mardupsGATK"
#' @param task Name of the task. Default "mardupsGATK"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] Hold job until job is finished. Job ID. 
#' @export

sort_and_index_bam_samtools=function(bin_path="tools/samtools/samtools",bam="",output_dir="",
verbose=FALSE,threads=3,ram=1,sort=TRUE,coord_sort=TRUE,index=TRUE,stats="all", clean=FALSE,
mode="local",executor=make_unique_id("sortANDindex"),task="sortANDindex",time="48:0:0",
update_time=60,wait=FALSE,hold=""){

  out_file_dir=set_dir(dir=output_dir)

  if(sort){
      job=sort_bam_samtools(bin_path=bin_path,bam=bam,output_dir=out_file_dir,
      ram=ram,verbose=verbose,threads=threads,coord_sort=coord_sort,clean=clean,
      executor=executor,mode=mode,time=time,
      update_time=update_time,wait=FALSE,hold=hold)

      out_file_dir=set_dir(dir=output_dir,name="sorted")
  
      if (coord_sort){
      
        bam=paste0(out_file_dir,"/",get_file_name(bam),".sorted.",get_file_ext(bam))

        if(index){
            job=index_bam_samtools(bin_path=bin_path,bam=bam,verbose=verbose,threads=threads,ram=ram,
            executor=executor,mode=mode,time=time,
            update_time=update_time,wait=FALSE,hold=job,output_dir = out_file_dir)
          if(stats=="index"|stats=="all"){
              job=stats_bam_samtools(bin_path=bin_path,bam=bam,output_dir=out_file_dir,
              verbose=verbose,threads=threads,stats="index",executor=executor,
              mode=mode,time=time,update_time=update_time,wait=FALSE,hold=job)
          }
        }
      }
  }else{
     job=index_bam_samtools(bin_path=bin_path,bam=bam,verbose=verbose,threads=threads,
     executor=executor,mode=mode,time=time,update_time=update_time,wait=FALSE,hold=hold)
  }


  if(stats=="flag"|stats=="all"){
    job=stats_bam_samtools(bin_path=bin_path,bam=bam,output_dir=out_file_dir,
      verbose=verbose,threads=threads,stats="flag",executor=executor,
      mode=mode,time=time,update_time=update_time,
      wait=FALSE,hold=hold)
  }

  if(wait&&mode=="batch"){
      job_validator(job=job,
      time=update_time,verbose=verbose,threads=threads)
  }

  return(job)
  
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
#' @param executor Name of the executor. Default "mardupsGATK"
#' @param task Name of the task. Default "mardupsGATK"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] Hold job until job is finished. Job ID. 
#' @export

sort_bam_samtools=function(bin_path="tools/samtools/samtools",bam="",output_dir="",
verbose=FALSE,threads=3,ram=1,coord_sort=TRUE,mode="local",
executor=make_unique_id("sortBAM"),clean=FALSE,task="sortBAM",
time="48:0:0",update_time=60,wait=FALSE,hold=""){

  out_file_dir=set_dir(dir=output_dir,name="sorted")

  sort_type=""
  
  if(!coord_sort){
    sort_type=" -n "
  }
  exec_code=paste0(bin_path," sort ",sort_type, bam," -@ ",threads," -m ",ram,"G"," -o ",
  paste0(out_file_dir,"/",get_file_name(bam),".sorted.",get_file_ext(bam)))

  if(clean){
    exec_code=paste(exec_code," && rm",paste(bam,collapse=" "))
  }

  job=build_job(executor=executor,task=make_unique_id(task))
  if(mode=="batch"){
    out_file_dir2=set_dir(dir=out_file_dir,name="batch")
    batch_code=build_job_exec(job=job,time=time,ram=ram,threads=threads,output_dir=out_file_dir2)
    exec_code=paste0("echo 'source ~/.bashrc;",exec_code,"'|",batch_code)
  }
    
    
  if (verbose){
    print(exec_code)
  }
  error=system(exec_code)
  if(error!=0){
    stop("samtools failed to run due to unknown error.
    Check std error for more information.")
  }

  if(wait&&mode=="batch"){
      job_validator(job=job,
      time=update_time,verbose=verbose,threads=threads)
  }

  return(job)
}






index_bam_samtools=function(bin_path="tools/samtools/samtools",bam="",verbose=FALSE,threads=3,ram=4,
mode="local",executor=make_unique_id("indexBAM"),task="indexBAM",time="48:0:0",update_time=60, output_dir="",
wait=FALSE,hold=""){
  out_file_dir=set_dir(dir=output_dir)
  exec_code=paste(bin_path," index",bam," -@ ",threads)
  job=build_job(executor=executor,task=make_unique_id(task))
  if(mode=="batch"){
    out_file_dir2=set_dir(dir=out_file_dir,name="batch")
    batch_code=build_job_exec(job=job,time=time,ram=ram,threads=threads,output_dir=out_file_dir2)
    exec_code=paste0("echo 'source ~/.bashrc;",exec_code,"'|",batch_code)
  }
    
  if (verbose){
    print(exec_code)
  }
     error=system(exec_code)
  if(error!=0){
    stop("samtools failed to run due to unknown error.
    Check std error for more information.")
  }

  if(wait&&mode=="batch"){
      job_validator(job=job,
      time=update_time,verbose=verbose,threads=threads)
  }

  return(job)
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
#' @param executor Name of the executor. Default "mardupsGATK"
#' @param task Name of the task. Default "mardupsGATK"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] Hold job until job is finished. Job ID.
#' @export


stats_bam_samtools=function(bin_path="tools/samtools/samtools",bam="",output_dir="",
verbose=FALSE,threads=3,ram=4,stats="all",mode="local",executor=make_unique_id("statsBAM"),
task="statsBAM",time="48:0:0",update_time=60,wait=FALSE,hold=""){


  out_file_dir=set_dir(dir=output_dir,name="stats")

  if(stats=="all"|stats=="flag"){
    job1=stats_flag_samtools(bin_path=bin_path,bam=bam,output_dir=out_file_dir,
    verbose=verbose,threads=threads,mode=mode,time=time,update_time=update_time,wait=FALSE,hold=hold)
  }


  if(stats=="all"|stats=="index"){
    job2=stats_index_samtools(bin_path=bin_path,bam=bam,output_dir=out_file_dir,
    verbose=verbose,threads=threads,mode=mode,time=time,update_time=update_time,wait=FALSE,hold=hold)
  }

  
  if(wait&&mode=="batch"){
      job_validator(job=c(job1,job2),
      time=update_time,verbose=verbose,threads=threads)
  }

  return(c(job1,job2))
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
#' @param executor Name of the executor. Default "mardupsGATK"
#' @param task Name of the task. Default "mardupsGATK"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] Hold job until job is finished. Job ID. 
#' @export



stats_flag_samtools=function(bin_path="tools/samtools/samtools",bam="",output_dir="",
verbose=FALSE,threads=3,ram=4,mode="local",executor=make_unique_id("statsFlag"),
task="statsFlag",time="48:0:0",update_time=60,wait=FALSE,hold=""){

  out_file_dir=set_dir(dir=output_dir,name="flag")

  exec_code=paste0(bin_path," flagstat ",bam," -@ ",threads," > ",
    paste0(out_file_dir,"/",get_file_name(bam),".flagstat.txt"))

  job=build_job(executor=executor,task=make_unique_id(task))
  if(mode=="batch"){
    out_file_dir2=set_dir(dir=out_file_dir,name="batch")
    batch_code=build_job_exec(job=job,time=time,ram=ram,threads=threads,output_dir=out_file_dir2)
    exec_code=paste0("echo 'source ~/.bashrc;",exec_code,"'|",batch_code)
  }
    
  
  if (verbose){
    print(exec_code)
  }


  error=system(exec_code)
  if(error!=0){
    stop("samtools failed to run due to unknown error.
    Check std error for more information.")
  }

  
  if(wait&&mode=="batch"){
      job_validator(job=job,
      time=update_time,verbose=verbose,threads=threads)
  }

  return(job)
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
#' @param executor Name of the executor. Default "mardupsGATK"
#' @param task Name of the task. Default "mardupsGATK"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] Hold job until job is finished. Job ID. 
#' @export


stats_index_samtools=function(bin_path="tools/samtools/samtools",bam="",output_dir="",
verbose=FALSE,threads=3,ram=4,mode="local",executor=make_unique_id("statsINDEX"),
task="statsINDEX",time="48:0:0",update_time=60,wait=FALSE,hold=""){

  out_file_dir=set_dir(dir=output_dir,name="index")

  exec_code=paste0(bin_path," idxstats ",bam," > ",paste0(out_file_dir,"/",get_file_name(bam),".idxstats.txt"))


  job=build_job(executor=executor,task=make_unique_id(task))
  if(mode=="batch"){
    out_file_dir2=set_dir(dir=out_file_dir,name="batch")
    batch_code=build_job_exec(job=job,time=time,ram=ram,threads=threads,output_dir=out_file_dir2)
    exec_code=paste0("echo 'source ~/.bashrc;",exec_code,"'|",batch_code)
  }
    
  if (verbose){
    print(exec_code)
  }
     error=system(exec_code)
  if(error!=0){
    stop("samtools failed to run due to unknown error.
    Check std error for more information.")
  }

  
  if(wait&&mode=="batch"){
      job_validator(job=job,
      time=update_time,verbose=verbose,threads=threads)
  }

  return(job)

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
#' @param executor Name of the executor. Default "mardupsGATK"
#' @param task Name of the task. Default "mardupsGATK"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] Hold job until job is finished. Job ID. 
#' @export

metrics_mapq_bam_samtools=function(bin_path="tools/samtools/samtools",bam="",output_dir="",
verbose=FALSE,threads=3,ram=4,mode="local",executor=make_unique_id("metricsMAPQ"),
task="metricsMAPQ",time="48:0:0",update_time=60,wait=FALSE,hold=""){


  out_file_dir=set_dir(dir=output_dir,name="mapq")


  exec_code=paste(bin_path,"view",bam," -@ ",threads, " | awk -F", "'\\t'",
    "'{c[$5]++} END { for (i in c) printf(\"%s\\t%s\\n\",i,c[i]) }'",
    " | sort -t$'\\t' -k 1 -g >>", paste0(out_file_dir,"/",get_file_name(bam),".mapq_dist.txt"))

  job=build_job(executor=executor,task=make_unique_id(task))
  if(mode=="batch"){
    out_file_dir2=set_dir(dir=out_file_dir,name="batch")
    batch_code=build_job_exec(job=job,time=time,ram=ram,threads=threads,output_dir=out_file_dir2)
    exec_code=paste0("echo 'source ~/.bashrc;",exec_code,"'|",batch_code)
  }
    
  
  if (verbose){
    print(exec_code)
  }

    error=system(exec_code)
  if(error!=0){
    stop("samtools failed to run due to unknown error.
    Check std error for more information.")
  }
  
  if(wait&&mode=="batch"){
      job_validator(job=job,
      time=update_time,verbose=verbose,threads=threads)
  }

  return(job)
}


