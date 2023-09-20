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

new_sort_and_index_bam_samtools=function(
  bin_samtools=build_default_tool_binary_list()$bin_samtools,
  bam=NULL,
  sort=TRUE,
  index=TRUE,
  coord_sort=TRUE,
  stats=TRUE,
  ...
){


  run_main=function(
    .env
  ){
    .this.env=environment()
    append_env(to=.this.env,from=.env)

    set_main(.env=.this.env)


    .main$steps[[fn_id]]<-.this.env
    .main.step=.main$steps[[fn_id]]
  
    if(sort){
        .main.step$steps <-append(
          .main.step$steps ,
            new_sort_bam_samtools(
              bin_samtools=bin_samtools,
              bam=input,
              output_dir=out_file_dir,
              tmp_dir=tmp_dir,
              env_dir=env_dir,
              batch_dir=batch_dir,
              ram=ram,
              index=index,
              stats=stats,
              verbose=verbose,
              threads=threads,
              coord_sort=coord_sort,
              err_msg=err_msg,
              clean=clean,
              executor_id=task_id
        )
      )
          .this.step=.main.step$steps$new_sort_bam_samtools
          .main.step$out_files=append(.main.step$out_files,.this.step$out_files)

    }else{
        .main.step$steps  <-append(
          .main.step$steps ,
            new_index_bam_samtools(
              bin_samtools=bin_samtools,
              bam=input,
              stats=stats,
              tmp_dir=tmp_dir,
              env_dir=env_dir,
              batch_dir=batch_dir,
              verbose=verbose,
              err_msg=err_msg,
              threads=threads,
              executor_id=task_id
          )
      )
        .this.step=.main.step$steps$new_index_bam_samtools
        .main.step$out_files=append(.main.step$out_files,.this.step$out_files)
    }

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

new_sort_bam_samtools=function(
  bin_samtools=build_default_tool_binary_list()$bin_samtools,
  bam=NULL,
  index=TRUE,
  stats=TRUE,
  coord_sort=TRUE,
  ...
){



    run_main=function(
      .env
    ){
      

      .this.env=environment()
      append_env(to=.this.env,from=.env)

      set_main(.env=.this.env)
      
     
      .main$out_file=paste0(out_file_dir,"/",input_id,".sorted.",input_ext)

    
      .main$exec_code=paste0(
        bin_samtools," sort ",
        ifelse(coord_sort,""," -n "), 
        input,
        " -@ ",threads,
        " -m ",ram,
        "G"," -o ",
        .main$out_file
      )
    
      run_job(.env=.this.env)

      .main.step=.main$steps[[fn_id]]
      .main.step$out_files$srt_bam=.main.step$out_file
      

      if(index & coord_sort){

          .main.step$steps <-append(
            .main.step$steps ,
              new_index_bam_samtools(
                  bin_samtools=bin_samtools,
                  bam=.main$out_file,
                  stats=stats,
                  tmp_dir=tmp_dir,
                  env_dir=env_dir,
                  batch_dir=batch_dir,
                  verbose=verbose,
                  threads=threads,
                  err_msg=err_msg,
                  ram=ram,
                  executor_id=task_id
                )
            )

          .this.step=.main.step$steps$new_index_bam_samtools
          .main.step$out_files=append(.main.step$out_files,.this.step$out_files)
      }

      if(stats){
          .main.step$steps <-append(
           .main.step$steps ,
              new_stats_bam_samtools(
                  bin_samtools=bin_samtools,
                  bam=.main$out_file,
                  stats="flag",
                  verbose=verbose,
                  threads=threads,
                  err_msg=err_msg,
                  ram=ram,
                  executor_id=task_id,
                  output_dir=paste0(out_file_dir,"/stats"),
                  tmp_dir=tmp_dir,
                  env_dir=env_dir,
                  batch_dir=batch_dir
                )
            )

          .this.step=.main.step$steps$new_stats_bam_samtools
          .main.step$out_files=append(.main.step$out_files,.this.step$out_files) 
      }

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






#' Index a BAM file
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




#' Index a BAM file
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



new_index_bam_samtools=function(
  bin_samtools=build_default_tool_binary_list()$bin_samtools,
  bam=NULL,
  stats=TRUE,
  ...
){


  run_main=function(
    .env
  ){


    .this.env=environment()
    append_env(to=.this.env,from=.env)

    set_main(.env=.this.env)

    .main$out_file=paste0(input,".bai")
    .main$exec_code=paste(
      bin_samtools," index",
      input," -@ ",threads
    )
   
    run_job(.env=.this.env)

    .main.step=.main$steps[[fn_id]]
    .main.step$out_files$index_bam=.main.step$out_file

    if(stats){
      .main.step$steps<-append(
       .main.step$steps,
        new_stats_bam_samtools(
                  bin_samtools=bin_samtools,
                  bam=input,
                  output_dir=paste0(dirname(input),"/stats"),
                  tmp_dir=tmp_dir,
                  env_dir=env_dir,
                  batch_dir=batch_dir,
                  verbose=verbose,
                  threads=threads,
                  err_msg=err_msg,
                  stats="index",
                  executor_id=task_id
        )
      )
      
      .this.step=.main.step$steps$new_stats_bam_samtools
      .main.step$out_files=append(.main.step$out_files,.this.step$out_files) 
    } 

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


new_stats_bam_samtools=function(
  bin_samtools=build_default_tool_binary_list()$bin_samtools,
  bam=NULL,
  stats="all",
  ...
){

  run_main=function(
    .env
  ){


    .this.env=environment()
    append_env(to=.this.env,from=.env)
    set_main(.env=.this.env)

    .main$steps[[fn_id]]<-.this.env
    .main.step=.main$steps[[fn_id]]
    

    if(stats=="all"|stats=="flag"){
         .main.step$steps<-append(
           .main.step$steps,
            new_flag_stats_samtools(
              bin_samtools=bin_samtools,
              bam=input,
              output_dir=out_file_dir,
              tmp_dir=tmp_dir,
              env_dir=env_dir,
              batch_dir=batch_dir,
              verbose=verbose,
              threads=threads,
              err_msg=err_msg,
              ram=ram,
              executor_id=task_id
          )
        )
        .this.step=.main.step$steps$new_flag_stats_samtools
        .main.step$out_files=append(.main.step$out_files,.this.step$out_files) 
    } 

    if(stats=="all"|stats=="index"){
       .main.step$steps<-append(
         .main.step$steps,
          new_index_stats_samtools(
            bin_samtools=bin_samtools,
            bam=input,
            output_dir=out_file_dir,
            tmp_dir=tmp_dir,
            env_dir=env_dir,
            batch_dir=batch_dir,
            verbose=verbose,
            threads=threads,
            err_msg=err_msg,
            ram=ram,
            executor_id=task_id
          )
        )
        .this.step=.main.step$steps$new_index_stats_samtools
        .main.step$out_files=append(.main.step$out_files,.this.step$out_files) 
      }

    .env$.main <-.main

  }


  .base.env=environment()
  list2env(list(...),envir=.base.env)
  set_env_vars(
    .env= .base.env,
    vars="bam"
  )
 
  launch(.env=.base.env)

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



new_flag_stats_samtools=function(
  bin_samtools=build_default_tool_binary_list()$bin_samtools,
  bam=NULL,
  ...
){

  
  run_main=function(
    .env
  ){

    .this.env=environment()
    append_env(to=.this.env,from=.env)


    set_main(.env=.this.env)
   

    .main$out_file=paste0(
      out_file_dir,"/",
      input_id,".flagstat.txt"
    )
    .main$exec_code=paste0(
      bin_samtools," flagstat ",
      input," -@ ",
      threads," > ",
      .main$out_file
    )

    run_job(.env=.this.env)

    .main.step=.main$steps[[fn_id]]

    .main.step$out_files$flag_stats <- .main$out_file

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


new_index_stats_samtools=function(
  bin_samtools=build_default_tool_binary_list()$bin_samtools,
  bam=NULL,
  ...
){  


   run_main=function(
    .env
  ){

    .this.env=environment()
    append_env(to=.this.env,from=.env)
    
    
    set_main(.env=.this.env)
  

    .main$out_file=paste0(
      out_file_dir,"/",input_id,".idxstats.txt"
    )

    .main$exec_code=paste0(
        bin_samtools," idxstats ",
        input," > ",.main$out_file
    )

    run_job(.env=.this.env)

    .this.step=.main$steps[[fn_id]]
    
    .this.step$out_files$index_stats <- .main$out_file

    .env$.main <- .main

  }

  .base.env=environment()
  list2env(list(...),envir=.base.env)
  set_env_vars(
    .env=.base.env,
    vars="bam"
  )

  launch(.env=.base.env)



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
  region=NA,
  min_frag_size=1,
  max_frag_size=167,
  flags=c(99, 147, 83, 163),
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

  fg=""
  if(!is.na(flags)){
    fg=paste0(" -f ",paste0(flags,collapse=","))
  }
  

  position=""
  if(!is.na(region)){
    position=strsplit(region,split="__")[[1]][2]
    out_file=paste0(out_file_dir,"/",id,".",ifelse(include,"include_","exclude_"),
    min_frag_size,"_",max_frag_size,".",region,".bam")
  }else{
    out_file=paste0(out_file_dir,"/",id,".",ifelse(include,"include_","exclude_"),
    min_frag_size,"_",max_frag_size,".bam")
  }


  if(include){
    exec_code=paste(bin_samtools,"view ",fg,bam,position," | \ awk 'substr($0,1,1)==\"@\""," || ($9>=",
    min_frag_size,"&& $9<=",max_frag_size,") ||", "($9<=-",min_frag_size,"&& $9>=-",max_frag_size,
    ")'|",bin_samtools, "view -b >",out_file)
  }else{
    exec_code=paste(bin_samtools,"view ",fg,bam,position," | \ awk 'substr($0,1,1)==\"@\""," || ($9=<",
    min_frag_size,"&& $9>=",max_frag_size,") ||", "($9>=-",min_frag_size,"&& $9<=-",max_frag_size,
    ")'|",bin_samtools, "view -b >",out_file)

  }

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
#' @param flags Read flags to filter by. Default c(99, 147, 83, 163)
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
  bam="",bed=NA,
  min_frag_size=1,
  max_frag_size=167,
  flags=c(99, 147, 83, 163),
  sep="\t",
  header=FALSE,
  verbose=FALSE,
  output_dir=".",
  output_name="",
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


  id=""
  if(output_name!=""){
    id=output_name
  }else{
    id=get_file_name(bam)
  }

  
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
  sep="\t",header=TRUE,stringsAsFactors = FALSE)
  bam_chr$order=as.numeric(as.factor(bam_chr$chr))

  if(!is.na(bed)){
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
        jobs_report[["steps"]][["par_region_fragment_length"]]<-
        parallel::mclapply(region_list,FUN=function(region){
        job_report <- filter_bam_by_size_samtools(
          bin_samtools=bin_samtools,
          bam=bam,
          region=region,
          min_frag_size=min_frag_size,
          max_frag_size=max_frag_size,
          include=include, 
          verbose=verbose,
          flags=flags,
          output_name=id,
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
            flags,
            include,
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
         
         jobs_report[["steps"]][["par_sample_fragment_length"]]<- build_job_report(
              job_id=job,
              executor_id=executor_id,
              exec_code=exec_code, 
              task_id=task_id,
              input_args=argg,
              out_file_dir=out_file_dir,
              out_files=list(
                  frag_bam=paste0(out_file_dir,"/",get_file_name(bam),".",
                  ifelse(include,"include_","exclude_"),
                  min_frag_size,"_",max_frag_size,".",names(region_list),".",
                  ".bam")
              )
        )
    }

    jobs_report[["steps"]][["gather_bam"]]<-gather_bam_files_picard(
      bin_picard=bin_picard,
      bam=unlist_lvl(jobs_report[["steps"]][["par_region_fragment_length"]],var="frag_bam"),
      output_dir=out_file_dir,
      output_name=paste0(id,".",ifelse(include,"include_","exclude_"),
      min_frag_size,"_",max_frag_size),
      executor_id=task_id,mode=mode,time=time,threads=threads,ram=ram,
      update_time=update_time,wait=FALSE,
      clean=clean,
      hold=unlist_lvl(jobs_report[["steps"]][["par_region_fragment_length"]],
      var="job_id",recursive=TRUE)
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





#' Filter BAM file by size using samtools
#'
#'
#' @param bam Path to the input file with the sequence.
#' @param bin_samtools Path to samtools executable. Default path tools/samtools/samtools.
#' @param mapq Mapping quality of the read to analyze. Default 60.
#' @param mapq Flags of the reads to read. Default c(99, 147, 83, 163)
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


get_insert_size_samtools=function(
  rdata=NULL,
  selected=NULL,
  bin_samtools=build_default_tool_binary_list()$bin_samtools,
  bam="",
  region=NA,
  mq=0,
  flags=c(99, 147, 83, 163),
  output_name="",
  output_dir=".",
  verbose=FALSE,
  batch_config=build_default_preprocess_config(),
  threads=1,ram=4,mode="local",
  executor_id=make_unique_id("getInsertSizeBAM"),
  task_name="getInsertSizeBAM",time="48:0:0",
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

 fg=""
  if(!is.null(flags)){
    fg=paste0(" -f ",paste0(flags,collapse=","))
  }


  mapq=""
  if(!is.null(mq)){
    mapq=paste0(" -q ",mq)
  }
  
  


  position=""
  if(!is.na(region)){
    if(grepl("__",region)){
      position=strsplit(region,split="__")[[1]][2]
    }else{
      position=region
  }
    
    out_file=paste0(out_file_dir,"/",id,".",region,".fragments.txt")
  }else{
    out_file=paste0(out_file_dir,"/",id,".fragments.txt")
  }


  exec_code=paste(bin_samtools,"view ",fg,mapq,bam,position,
   " | awk '{
      mot = substr($10, 1, 4);
      fl=($9^2)^(1/2);
      fl_count[NR] = fl;
      fl_dist[fl] = fl_dist[fl]+1;
      motif_dist[mot] = motif_dist[mot]+1;

    }END{
        fl_str_dist=\"\";
        fl_median = 0;
        fl_average = 0;
        fl_sd = 0;
        fl_mode= 0;
        fl_max= 0;
        for( fl in fl_dist ) {
            if(fl_str_dist!=\"\"){
              fl_str_dist = fl_str_dist\"|\"fl\":\"fl_dist[fl];
            } else{
              fl_str_dist = fl\":\"fl_dist[fl];
            }
            
            if(fl_max<=fl_dist[fl]){
              fl_max=fl_dist[fl];
              fl_mode=fl;
            }

        }

        motif_str_dist=\"\";
        motif_max= 0;
        motif_mode= 0;

        for( mot in motif_dist ) {
            if(motif_str_dist!=\"\"){
              motif_str_dist  = motif_str_dist\"|\"mot\":\"motif_dist[mot];
            }else{
               motif_str_dist =mot\":\"motif_dist[mot];
            }

            if(motif_max<=motif_dist[mot]){
              motif_max=motif_dist[mot];
              motif_mode=mot;
            }
        }

        if (NR > 1) {
            if ((NR % 2) == 1) {
                fl_median = fl_count[(NR + 1) / 2];
                motif_median=motif_count[(NR +1)/2];
            } else {
                fl_median = (fl_count[NR / 2] + fl_count[(NR / 2) + 1]) / 2.0;
        
            }
            fl_sum = 0;
            for( i = 1; i <= length( fl_count ); i++ ) {
                fl_sum += fl_count[i];
            }
            fl_average = fl_sum / NR
            fl_sumsd = 0;
            for( i = 1; i <= length( fl_count ); i++ ) {
                fl_sumsd += (fl_count[i] - fl_average)^2;
            }
            fl_sd = (fl_sumsd /(NR - 1))^(1/2);
        } else {
            if (NR == 1) {
                fl_median = fl_count[1];
                fl_average = fl_count[1];
                fl_sd = 0;
            } else {
                fl_median = 0;
                fl_average = 0;
                fl_sd = 0;
            }
        };
      printf(\"ID\\tFLAGS\\tMAPQ\\tREGION\\tTOTAL\\tfl_median\\tfl_mode\\tfl_max\\tfl_average\\tfl_sd\\tmotif_mode\\tmotif_max\\tfl_str_dist\\tmotif_str_dist\\n\");
      printf(\"",id,"\\t",paste0(flags,collapse=","),"\\t",mq,"\\t",ifelse(position=="","GENOME",position),
      "\\t%d\\t%d\\t%d\\t%d\\t%d\\t%d\\t%s\\t%s\\t%s\\t%s\\n\", NR , fl_median, fl_mode, fl_max, fl_average , fl_sd , motif_mode , motif_max, fl_str_dist , motif_str_dist);}'> ",out_file
  )


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
    stop("samtools failed to run due to unknown error.
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
#' @param flags Read flags to filter by. Default c(99, 147, 83, 163)
#' @param threads Number of threads. Default 3
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Job EXECUTOR ID. Default "mardupsGATK"
#' @param task_name Name of the task. Default "mardupsGATK"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] Hold job until job is finished. Job ID. 
#' @export




parallel_region_get_insert_size_samtools=function(
  bin_samtools=build_default_tool_binary_list()$bin_samtools,
  bin_picard=build_default_tool_binary_list()$bin_picard,
  bam="",
  bed=NULL,
  mq=0,
  flags=c(99, 147, 83, 163),
  sep="\t",
  header=FALSE,
  verbose=FALSE,
  output_dir=".",
  output_name="",
  include=TRUE,
  index=TRUE,
  clean=TRUE,
  batch_config=build_default_preprocess_config(),
  threads=3,ram=4,mode="local",
  executor_id=make_unique_id("parRegionGetInsertSizeBAM"),
  task_name="parRegionGetInsertSizeBAM",time="48:0:0",
  update_time=60,wait=FALSE,hold=NULL
){

  options(scipen=999)

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir)
  out_file_dir_tmp=set_dir(dir=out_file_dir,name="tmp")
  job=build_job(executor_id=executor_id,task_id=task_id)


  id=""
  if(output_name!=""){
    id=output_name
  }else{
    id=get_file_name(bam)
  }

  
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
  sep="\t",header=TRUE,stringsAsFactors = FALSE)
  bam_chr$order=as.numeric(as.factor(bam_chr$chr))

  if(!is.null(bed)){
      regions=read.table(bed,sep=sep,header=header)
      names(regions)[1:3]=c("chr","start","end")
      miss_chr=setdiff(unique(regions[,1]), bam_chr$chr)
    
      ### Validate chromosome format in BED file

      if(length(miss_chr)>0){
        stop(paste0(paste(miss_chr,collapse=",")," are not present in the BAM file."))
      }

      regions=dplyr::left_join(regions,bam_chr[,c("chr","order")])
      regions=regions %>% dplyr::arrange(order,start,end)
  }else{
      regions=bam_chr
  }

  regions$pos=1:nrow(regions)
  regions=regions %>% dplyr::mutate(region=paste0(pos,"__",chr,":",start,"-",end))
    
  region_list=regions$region
  names(region_list)=regions$region

  if(mode=="local"){
        jobs_report[["steps"]][["par_region_insert_size"]]<-
        parallel::mclapply(region_list,FUN=function(region){
        job_report <- get_insert_size_samtools(
          bin_samtools=bin_samtools,
          bam=bam,
          region=region,
          mq=mq,
          verbose=verbose,
          flags=flags,
          output_name=id,
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
            output_dir,
            verbose,
            mq,
            flags,
            file = rdata_file)
          exec_code=paste0("Rscript -e \"ULPwgs::get_insert_size_samtools(rdata=\\\"",
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
              stop("samtools failed to run due to unknown error.
              Check std error for more information.")
          }
         
         jobs_report[["steps"]][["par_region_insert_size"]]<- build_job_report(
              job_id=job,
              executor_id=executor_id,
              exec_code=exec_code, 
              task_id=task_id,
              input_args=argg,
              out_file_dir=out_file_dir,
              out_files=list(
                  frag_bam=paste0(out_file_dir,"/",get_file_name(bam),".",names(region_list),".txt")
              )
        )
    }

    return(jobs_report)

}



#' Generate mpileup file from BAM file
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


mpileup_samtools=function(
  bam=NULL,
  gpos=NULL
){
   run_main=function(
    .env
  ){
    .this.env=environment()
    append_env(to=.this.env,from=.env)
   
    set_main(.env=.this.env)
    
    gpos=read_gpos(gpos)

    .main$out_file=paste0(
      out_file_dir,"/",input_id,".pileup.txt"
    )

    .main$exec_code=paste0(
        bin_samtools," mpileup ",
        input," > ",.main$out_file
    )

    run_job(.env=.this.env)

    .this.step=.main$steps[[fn_id]]
    
    .this.step$out_files$index_stats <- .main$out_file

    .env$.main <- .main

  }

  .base.env=environment()
  list2env(list(...),envir=.base.env)
  set_env_vars(
    .env=.base.env,
    vars="bam"
  )

  launch(.env=.base.env)


}






#' Filter BAM file by size using samtools
#'
#'
#' @param bam Path to the input file with the sequence.
#' @param bin_samtools Path to samtools executable. Default path tools/samtools/samtools.
#' @param mapq Mapping quality of the read to analyze. Default 60.
#' @param mapq Flags of the reads to read. Default c(99, 147, 83, 163)
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


new_get_insert_size_samtools=function(
  bin_samtools=build_default_tool_binary_list()$bin_samtools,
  bam=NULL,
  region=NULL,
  mapq=NULL,
  flags=c(99, 147, 83, 163),
  min_fl=1,
  max_fl=1000,
  ignore_N_bases=TRUE,
  ...
){


   run_main=function(
        .env
    ){
        .this.env=environment()
        append_env(to=.this.env,from=.env)
   
        set_main(.env=.this.env)


        add=""
          if(!is.null(flags)){
            add=paste0(" -f ",paste0(flags,collapse=","))
          }

        if(!is.null(mapq)){
          add=paste0(" -q ",mapq)
        }
       
        if(!is.null(input)){
          .main$out_files$frags=paste0(out_file_dir,"/",input_id,".",input,".fragments.txt")
        }else{
          .main$out_files$frags=paste0(out_file_dir,"/",input_id,".fragments.txt")
        }

        position="GENOME"
        if(!is.null(input)){
          position=input
        }
        reg=""
        if(ignore_N_bases){
          reg=" && index(mot,\"N\")==0"
        }

      
      .main$exec_code=paste(bin_samtools,"view ",add,bam,input," -@ ",threads,
      " | awk '{
          mot = substr($10, 1, 4);
          fl=($9^2)^(1/2);",
          paste0("if(fl >= ",min_fl,"&& fl <= ",max_fl,reg,")"),
          "{
            fl_count[NR] = fl;
            fl_dist[fl] = fl_dist[fl]+1;
            motif_dist[mot] = motif_dist[mot]+1;
          };

        }END{
            fl_str_dist=\"\";
            fl_median = 0;
            fl_average = 0;
            fl_sd = 0;
            fl_mode= 0;
            fl_max= 0;
            for( fl in fl_dist ) {
                if(fl_str_dist!=\"\"){
                  fl_str_dist = fl_str_dist\"|\"fl\":\"fl_dist[fl];
                } else{
                  fl_str_dist = fl\":\"fl_dist[fl];
                }
                
                if(fl_max<=fl_dist[fl]){
                  fl_max=fl_dist[fl];
                  fl_mode=fl;
                }

            }

            motif_str_dist=\"\";
            motif_max= 0;
            motif_mode= 0;

            for( mot in motif_dist ) {
                if(motif_str_dist!=\"\"){
                  motif_str_dist  = motif_str_dist\"|\"mot\":\"motif_dist[mot];
                }else{
                  motif_str_dist =mot\":\"motif_dist[mot];
                }

                if(motif_max<=motif_dist[mot]){
                  motif_max=motif_dist[mot];
                  motif_mode=mot;
                }
            }

            if (NR > 1) {
                if ((NR % 2) == 1) {
                    fl_median = fl_count[(NR + 1) / 2];
                    motif_median=motif_count[(NR +1)/2];
                } else {
                    fl_median = (fl_count[NR / 2] + fl_count[(NR / 2) + 1]) / 2.0;
            
                }
                fl_sum = 0;
                for( i = 1; i <= length( fl_count ); i++ ) {
                    fl_sum += fl_count[i];
                }
                fl_average = fl_sum / NR
                fl_sumsd = 0;
                for( i = 1; i <= length( fl_count ); i++ ) {
                    fl_sumsd += (fl_count[i] - fl_average)^2;
                }
                fl_sd = (fl_sumsd /(NR - 1))^(1/2);
            } else {
                if (NR == 1) {
                    fl_median = fl_count[1];
                    fl_average = fl_count[1];
                    fl_sd = 0;
                } else {
                    fl_median = 0;
                    fl_average = 0;
                    fl_sd = 0;
                }
            };
          printf(\"ID\\tFLAGS\\tMAPQ\\tREGION\\tTOTAL\\tfl_median\\tfl_mode\\tfl_max\\tfl_average\\tfl_sd\\tmotif_mode\\tmotif_max\\tfl_str_dist\\tmotif_str_dist\\n\");
          printf(\"",input_id,"\\t",paste0(flags,collapse=","),"\\t",mapq,"\\t",position,
          "\\t%d\\t%d\\t%d\\t%d\\t%d\\t%d\\t%s\\t%s\\t%s\\t%s\\n\", NR , fl_median, fl_mode, fl_max, fl_average , fl_sd , motif_mode , motif_max, fl_str_dist , motif_str_dist);}'> ",.main$out_files$frags
      )
      run_job(.env=.this.env)
      .env$.main <- .main
    }
    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env=.base.env,
      output_name=get_file_name(bam),
      vars="region"
    )

  launch(.env=.base.env)


}
 




#' Sort and index a sequence file
#' 
#'
#' Wrapper around index_bam_samtools and sort_bam_samtools functions
#' 
#' @param bin_samtools Path to samtools executable. Default path tools/samtools/samtools.
#' @param bam Path to BAM file.
#' @param read_fraction Fraction of reads to subsample.
#' @export

sample_bam_samtools=function(
  bin_samtools=build_default_tool_binary_list()$bin_samtools,
  bam=NULL,
  read_fraction=NULL,
  ...
){


  run_main=function(
    .env
  ){
    .this.env=environment()
    append_env(to=.this.env,from=.env)

    set_main(.env=.this.env)

    .main$out_files$sampled_bam<- paste0(out_file_dir,"/",input_id,".subsampled_",read_fraction,".bam")
     .main$exec_code=paste0(
      bin_samtools," view -b -s ",
      read_fraction,
      input," -@ ",
      threads," > ",
      .main$out_file
    )

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



#' Sort and index a sequence file
#' 
#'
#' Wrapper around index_bam_samtools and sort_bam_samtools functions
#' 
#' @param bin_samtools Path to samtools executable. Default path tools/samtools/samtools.
#' @param bam Path to BAM file.
#' @param read_fraction Fraction of reads to subsample.
#' @export

sample_bam_samtools=function(
  bin_samtools=build_default_tool_binary_list()$bin_samtools,
  bam=NULL,
  read_fraction=NULL,
  ...
){


  run_main=function(
    .env
  ){
    .this.env=environment()
    append_env(to=.this.env,from=.env)

    set_main(.env=.this.env)

    .main$out_files$sampled_bam<- paste0(out_file_dir,"/",input_id,".subsampled_",read_fraction,".bam")
     .main$exec_code=paste0(
      bin_samtools," view -b -s ",
      read_fraction,
      input," -@ ",
      threads," > ",
      .main$out_file
    )

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