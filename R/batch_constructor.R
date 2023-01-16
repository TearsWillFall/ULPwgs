#' Build job for SGE executor
#' Job constructor only uses the top level Executor ID 
#' and the bottom level task ID.
#' 
#' Middle level task are omited in Job ID for simplicity
#'
#' @param executor Job EXECUTOR ID
#' @param task Task ID
#' @export




build_job=function(executor_id="executor",task_id="task"){
options(scipen = 999)
  executor=paste0("EXECUTOR_",executor_id)
  task=paste0("TASK_",task_id)
  job=paste0(c(executor,task),collapse=".")
  return(job)
}



#' Build a job report data structure
#' 
#' Data structure for job_report:
#' job_id job_order out_files
#'
#'
#' @param job_id Job ID
#' @param executor_id Task executor ID
#' @param task_id Task ID
#' @param exec_code Execution code
#' @param input Input arguments for job
#' @param job_order Job execution order
#' @param out_file_fir Output dir for 
#' @param out_files Output files for job
#' @export

build_job_report=function(job_id="job_1",executor_id="",task_id="",exec_code="",
job_order=1, input_args="",out_file_dir="", out_files=list(file="file")){
  options(scipen = 999)
  job_report=list(job_id=job_id,executor_id=executor_id,task_id=task_id,
  job_order=job_order,input_args=input_args,exec_code=exec_code,
  out_file_dir=out_file_dir,out_files=out_files)
  return(job_report)
}


#' Build job executor in SGE
#'
#' @param job Name of job or jobs.
#' @param time Time in seconds between checks. Default 10.
#' @param output_dir [OPTIONAL] PATH to output directory. 
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param threads [OPTIONAL] Number of threads to use to. Default 3.
#' @param ram  [OPTIONAL] RAM memory per thread requested. Default 4.
#' @param hold Job IDs to hold job.
#' 
#' @export

build_job_exec=function(job="",time="48:0:0",ram=3,threads=1,
output_dir=".",hold=NULL,wd=getwd(),array=""){

  if(is.null(job)){ 
    stop("A job ID has to be supplied to uniquely identify the job.")
  }
  
  if(!is.null(hold)){
    hold=paste0(" -hold_jid ",paste0(hold,collapse=","))
  }
  if(array!=""){
      array=paste0(" -t 1-",array)
  }
  exec_code=paste("qsub -V -N ",job,array,paste0(" -l h_rt=",time),
  paste0(" -l mem=",ram,"G"), paste0(" -pe smp ",threads), paste0(" -wd ",wd), 
  paste0(" -o ",output_dir,"/",job,".std_out"),
  paste0(" -e ",output_dir,"/",job,".std_error"),hold)
}


#' Validate job submited in batch on sungrid based cluster
#'
#' @param job Name of job or jobs.
#' @param time Time in seconds between checks. Default 10.
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param threads [OPTIONAL] Number of threads to use to. Default 3.
#' @export


job_validator=function(job="",time=10,verbose=FALSE,threads=3){
  error=FALSE
  col_names=c("job_id","job_priority","job_name","user","status","start_time","nodes")
  exec_code="qstat -xml | tr '\n' ' ' | sed 's#<job_list[^>]*>#\\n#g'   | sed 's#<[^>]*>##g' | grep \" \" | column -t"
  
  while(!error){
   
    dat_info=try({
      dat_info=read.table(text=system(exec_code,intern=TRUE),fill=TRUE)
      names(dat_info)=col_names
      dat_info=dat_info[dat_info$job_name %in% job,]
      jobs=nrow(dat_info)
      },silent=TRUE)

    if (inherits(dat_info, "try-error")|jobs==0) {
      break
    }else{
      if(verbose){
          cat("----------------------------------")
          cat(dat_info)
          cat("----------------------------------")
      }
      if(any(grepl("E",dat_info$status))){
        error=TRUE
      }
      Sys.sleep(time)
    }
  }

  if(error){
      cat("----------------------------------")
      cat(dat_info)
      cat("----------------------------------")
      parallel::mclapply(1:nrow(dat_info),FUN=function(x){
        system(paste0("qdel ",dat_info[x,]$job_id))
      },mc.cores=threads)
    
      stop("One or more jobs failed")
  }

return()
    

}

#' Create an unique ID
#'
#' @param name Name of the ID
#' @param id Unique identfier. If not given selected randomly
#' @param sep Separator between name and ID
#' @export

make_unique_id=function(name,id=sample(1:100000000000000,1),sep="_"){
  options(scipen = 999)
  unique_name=paste0(c(name,id),collapse=sep)
  return(unique_name)
}


#' Check for job limit in Myriad SGE
#'
#' @param exec_code Execution code
#' @param delay Time to wait for job limit checks in seconds
#' @export


check_job_limit=function(job_limit=1000000){
  n_jobs=system("qstat -u $USER | wc -l",intern=TRUE)
  return(as.numeric(n_jobs)>=job_limit)
}


#' Execute job with delay if job limit in Myriad SGE reached
#' 
#' @param exec_code Execution code
#' @param delay Time to wait for job limit checks in seconds
#' @export

execute_job=function(exec_code){
  error=system(exec_code)
  return(error)
}




#' Build RData object to inherit in main_function
#' 
#' @param objects List of objects to save in RData object
#' @param job_id Job Identifier
#' @param output_dir Path to output directory 
#' @export


build_self=function(
  .env=environment()
){        
      .this.env=environment()
      append_env(to=.this.env,from=.env)
      self_file=paste0(out_file_dir_tmp,"/",job_id,".self.RData")
      saveRDS(.this.env,file = self_file)
      append_env(to=.env,from=.this.env)
}




      



#' Build execution innit for function to execute
#' 
#' @param env List of objects to save in RData object
#' @export

build_main=function(
  .env=environment()
  ){
  .this.env=environment()
  append_env(to=.this.env,from=.env)
  main_file=paste0(out_file_dir_tmp,"/",job_id,".main.RData")
  saveRDS(.env$steps,file=main_file)
  append_env(to=.env,from=.this.env)
}


#' Build execution innit for function to execute
#' 
#' @param .env List of objects to save in RData object
#' @export

wait_scheduler=function(.env){

    .this.env=environment()
    append_env(to=.this.env,from=.env)
    check=TRUE

    while(check){
      Sys.sleep(60)
      jobs_in_queue=suppressWarnings(system("qstat -r | grep  \"jobname\"",intern=TRUE))
      if(length(jobs_in_queue)>0){
         jobs_in_queue=gsub(".* ","",jobs_in_queue)
         if(!any(jobs_in_queue %in% job_id)){
            check=FALSE
         }
      }else{
        check=FALSE
      }
     
    }
    return()

}




#' Build execution innit for function to execute
#' 
#' @param .env List of objects to save in RData object
#' @export

read_main=function(
  .env=environment()
  ){
    
    .this.env=environment()
    append_env(to=.this.env,from=.env)

 
    ### Reads mains and updates values for enviroments with data

    if(n_inputs>1){
      main.envs=lapply(1:n_inputs,function(n){
          main.env=readRDS(main.envs[[n]]$main_file)
        
        })
    }else{
      main.envs=readRDS(main.envs[[1]]$main_file)
    }
    
    .env$main.envs<-main.envs
}





#' Build execution innit for function to execute
#' 
#' @param .env List of objects to save in RData object
#' @export



build_exec_innit=function(
        .env=environment()
){
        .this.env=environment()
        append_env(to=.this.env,from=.env)

        ### Use SGE TASK ID if mode is set to batch otherwise use value
        
        if(mode=="local"){
             
              exec_code=paste0("Rscript -e \" invisible(lapply(1:",n_inputs,
              ",FUN=function(select){",ns,"::",fn,"(inherit=\\\"",
              self_file,"\\\",select=select)}))\"")
        
        }else if(mode=="batch"){
              exec_code=paste0("Rscript -e \" invisible(",
              ns,"::",fn,"(inherit=\\\"",self_file,
              "\\\",select=$SGE_TASK_ID))\"")
        }else{
          stop("Unkown mode type")
        }

        append_env(to=.env,from=.this.env)
}





#' Run job from batch constructor
#' 
#' @param .env Inherit current enviroment.
#' 
#' 
#' @export



build_clean_exec=function(
  .env
){
  .env$exec_code<- paste0(.env$exec_code," && rm",paste(input))
}



#' Build execution innit for batch
#' 
#' @param .env Inherit envoment variables
#' @export






build_batch_exec_innit=function(
  .env=environment()
){  
    n_inputs<-1
    .this.env=environment()
    append_env(to=.this.env,from=.env)

  
    out_file_dir_batch=set_dir(dir=out_file_dir_tmp,name="batch")

    batch_code=build_job_exec(
      job=job_id,time=time,ram=ram,
      threads=threads,wd=getwd(),
      output_dir=out_file_dir_batch,
      hold=hold,array=n_inputs
    )

    exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,
    ";",exec_code,"'|",batch_code)

    append_env(to=.env,from=.this.env)
}


#' Set vars for self run
#' 
#' @param env Inherit current enviroment.
#' 
#' 
#' @export


set_self=function(
  .env
){
    .this.env=environment()
    append_env(to=.this.env,from=.env)

    error=0
    self_file=""
    main_code=""
    exec_code=""
    batch_code=""
    .env$.self=environment()
    return()

}



#' Run job from batch constructor
#' 
#' @param env Inherit current enviroment.
#' 
#' 
#' @export

run_self=function(
  .env=environment()
){  


    .this.env=environment()
    append_env(to=.this.env,from=.env)
    set_self(.env=.env)
    
    .self=.env$.self

    ### Create RData object to inherit vars

    build_self(.env=.self)

    build_exec_innit(
            .env=.self
    )


    if(mode=="batch"){
        build_batch_exec_innit(
              .env=.self
        )
    }


    if(verbose){
          print_verbose(job=.self$job_id,
            arg=as.list(.self),
            exec_code=.self$exec_code
          )
    }

    .self$error=execute_job(exec_code=.self$exec_code)

    if(.self$error!=0){
        stop(.self$err_msg)
    }


    if(mode=="batch"){
      if(wait){
         wait_scheduler(.env=.self)
         read_main(.env=.self)
      }
    }else{
      read_main(.env=.self)
    }

    build_self(.env=.self)
    
    if(is.null(.self$main.envs)){
      return(.self)
    }else{
      return(.self$main.envs)
    }
  
}

#' Run job from batch constructor
#' 
#' @param env Inherit current enviroment.
#' 
#' 
#' @export


 run_job=function(
  .env=environment()
 ){
      .this.env=environment()
      append_env(to=.this.env,from=.env$.main)

      if(clean){
        build_clean_exec(.env=.this.env)
      }

      if(verbose){
        print_verbose(
          job=job_id,
          arg=as.list(.this.env),
          exec_code=exec_code
        )
      }

      error=execute_job(
          exec_code=exec_code
        )
      
      if(error!=0){
          stop(err_msg)
      }
      
      .env$.main$steps[[fn]] <- .this.env
 }
   







#' Wrapper around qstat call for SGE
#' 
#' 
#' @export


qstat=function(){
  system("qstat")
}

#' Wrapper around qdel call for SGE
#' 
#' @param jobs Inherit current enviroment.
#' 
#' 
#' @export

qdel=function(jobs){
  system("qdel ",paste(jobs,collapse=" "))
}









#' Set default enviromental variables based on input for all
#' functions in package.
#' 
#' @param envir Inherit current enviroment.
#' @param id Input identifier. If not given name of input will be used.
#' @param name Name of output file directory
#' @export

set_env_vars=function(
  .env=environment(),
  vars=NULL,
  fn=NULL,
  output_dir=".",
  output_name=NULL,
  verbose=FALSE,
  bgzip_index=FALSE,
  index=TRUE,
  index_format="tbi",
  clean=FALSE,
  batch_config=build_default_preprocess_config(),
  threads=1,
  ram=4,
  ns="ULPwgs",
  mode="local",
  time="48:0:0",
  sheet=NULL,
  inherit=NULL,
  select=NULL,
  err_msg=NULL,
  executor_id=NULL,
  wait=TRUE,
  hold=NULL
){  

    .this.env=environment()
    append_env(to=.this.env,from=.env)
    .env$n_jobs <- 1

       
    if(!is.null(inherit)){
        if(!is.environment(inherit)){
          inherit <-readRDS(file=inherit)
        }
        append_env(to=.this.env,from=inherit)
        main.envs<-inherit$main.envs
        .renv=main.envs[[select]]
        .env$self.envs <- .this.env
        return()
    }
    
    if(is.null(fn)){

      ## GET CALLER FUNCTION NAME IF NOT GIVEN

    fn <- sub(".*::","",sub("\\(.*","",
        paste0(deparse(sys.calls()[[sys.nframe()-1]]),collapse=","))
      )
    }


    if(is.null(executor_id)){
      executor_id <- make_unique_id(fn)
    }

    if (!is.null(sheet)){
        read_sheet(.env=.this.env)
        set_ss_env(.env=.this.env)
        return()
    }


    if(!is.null(vars)){
      inputs <- get(vars)
      n_inputs <- length(inputs)
      inputs_id <- set_input_id(
        inputs=inputs,
        ids=output_name
      )
      inputs_ext <- unname(Vectorize(get_file_ext)(inputs))
    }


    out_file_dir <- set_dir(
          dir=output_dir
    )

    out_file_dir_tmp <- set_dir(
      dir=out_file_dir,
      name="tmp"
    )
    task_id <- make_unique_id(fn)

    job_id <- build_job(
      executor_id=executor_id,
      task_id=task_id
    )

    if(is.null(err_msg)){
        err_msg <- paste0("CRITICAL ERROR: ",fn," (",job_id,") "," -> ")
    }
    
    set_main_env(.env=.this.env)

    .env$self.envs<-.this.env

  
}


#' Set steps enviroment for use
#' 
#' @param .env Environment
#' @export


run=function(.env){

  .this.env=environment()
  append_env(to=.this.env,from=.env)

  if(is.null(select)){
    run_self(.env=.env)
  }else{
    .renv=.env$.renv
    run_main(.env=.renv)
    build_main(.env=.renv$.main)
  }
}


#' Set steps enviroment for use
#' 
#' @param .env Environment
#' @export

launch=function(.env){
      .this.env=environment()
      append_env(to=.this.env,from=.env)

      if(n_jobs>1){
        reports=lapply(1:n_jobs,function(n){
            run(.env=self.envs[[n]])
        })
      }else{
          reports=run(.env=self.envs)
      }
      return(reports)
  }


#' Set steps enviroment for use
#' 
#' @param .env Environment
#' @export


read_sheet=function(.env){

    .this.env=environment()
    append_env(to=.this.env,from=.env)
    .base.env=.env$.env

    sheet=read.delim(sheet,header=TRUE)
    n_total<-nrow(sheet)
    sheet=sheet %>% dplyr::distinct()
    n_jobs<-nrow(sheet)
    n_vars<-ncol(sheet)
    n_dup=n_total-n_jobs
    
    if(n_dup>0){
      warning(paste0(n_dup, " were duplicated in sheet"))
    }

    sheet=sheet %>% 
      dplyr::group_by(dplyr::across(-c(vars))) %>%
      dplyr::summarise(!! vars := list(!! rlang::sym(vars)))
    
    print(sheet)
    .env$n_jobs<- .base.env <- n_jobs
    .env$n_vars<- .base.env <- n_vars
    .env$sheet<- .base.env <- sheet
}






#' Set steps enviroment for use
#' 
#' @param .env Environment
#' @export


set_ss_env=function(.env){

    .this.env=environment()
    append_env(to=.this.env,from=.env)
    .base.env <- .env$.env
    this.sheet=sheet
    vars_sheet=colnames(this.sheet)

    lapply(seq(1,n_jobs),
      FUN=function(row,.env){
        .this.env=environment()
        append_env(to=.this.env,from=.env)
        sheet <- NULL

        ### ASSIGN VARS IN SHEET TO ENVIROMENT
      
        invisible(
            lapply(seq(1,n_vars),
              FUN=function(col)
            {
              .this.env[[vars_sheet[col]]] <- this.sheet[[row,col]]
            }
          )
        )
        
        set_env_vars(
            .env=.this.env
        )

        .env$self.envs[[row]] <- self.envs
      },.env=.this.env
    )
    .base.env$self.envs<-self.envs      
}




#' Set steps enviroment for use
#' 
#' @param .env Environment
#' @export


set_main_env=function(.env){
   
    .this.env=environment()
    append_env(to=.this.env,from=.env)
    lapply(
        1:n_inputs,
        function(n,.env){
          .this.env=environment()
          append_env(to=.this.env,from=.env)

          .env$self.envs <- .this.env
          input<-inputs[n]
          n_inputs<-1
          input_id<-inputs_id[n]
          input_ext<-inputs_ext[n]

          task_id <- make_unique_id(fn)

          job_id <- build_job(
            executor_id=executor_id,
            task_id=task_id
          )

          err_msg <- paste0(err_msg ,fn," (",job_id,") "," -> ")
          
          
          build_main(.env=.this.env)
        
          .env$main.envs[[n]]<-.this.env
        },
        .env=.this.env
    )
  .env$main.envs<-main.envs
}



#' Set steps enviroment for use
#' 
#' @param envir Environment
#' @export




set_main=function(
  .env=environment()
){     
      .this.env=environment()
      append_env(to=.this.env,from=.env)
      exec_code=""
      error=0
      out_file=""
      steps=list()
      out_files=list()
      .env$.main=environment()
      return()
  }



#' Append two enviroments
#' 
#' @param to Enviroment to append to
#' @param from Enviroment to append from
#' @export

append_env = function(to=environment(), from=NULL) {
      
      from_list = ls(from)
      for(var in from_list) {
        if(!grepl("\\.",var)){
          if(is.null(to[[var]])){
             to[[var]] <- NULL
          }
        }

        if(!is.null(from[[var]])){
           to[[var]] <- from[[var]]
        }
      }
}



#' Set default environment ID
#' 
#' @param inputs Input file/s
#' @param ids Input file identifier
#' @export

set_input_id=function(inputs,ids=NULL){
      if(!is.null(ids)){
        my_id=ids
      }else{
        my_id=unname(Vectorize(get_file_name)(inputs))
      }
      return(my_id)
}




#' @export

print_verbose=function(exec_code,arg=NULL,job,ws=1){
      rep(cat("    \n"),ws)
      cat(crayon::blue("Job:"))
      rep(cat("    \n"),ws)
      cat(paste0(crayon::red(job,"\n")))
      rep(cat("    \n"),ws)
      if(!is.null(arg)){
         cat(crayon::blue("Arguments:"))
         rep(cat("    \n"),ws)
         lapply(names(arg),FUN=function(ag){
          tryCatch({
            rep(cat("    \n"),ws);
            cat(paste0(crayon::silver(ag),": ",arg[[ag]],"\n"))},error=function(e){
              cat(paste0(crayon::silver(ag),": ","\n"))
            })
      })
      }
      rep(cat("    \n"),ws)
      cat(crayon::blue("Command:"))
      rep(cat("    \n"),ws)
      cat(paste0(crayon::green(exec_code,"\n")))
      rep(cat("    \n"),ws)
}






