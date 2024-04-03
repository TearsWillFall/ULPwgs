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
  job=lapply(task_id,FUN=function(id){
    task=paste0("TASK_",id)
    job=paste0(c(executor,task),collapse=".")
  })
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

build_job_report=function(job_id="job_1",
  executor_id="",
  task_id="",
  exec_code="",
  job_order=1,
  input_args="",
  out_file_dir="",
  out_files=list(file="file")
){
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

build_job_exec=function(
  job="",
  time="48:0:0",
  ram=3,
  threads=1,
  output_dir=".",
  hold=NULL,
  bypass=FALSE,
  wd=getwd(),
  array=""
){

  if(is.null(job)){ 
    stop("A job ID has to be supplied to uniquely identify the job.")
  }
  
  if(!is.null(hold)){
    hold=paste0(" -hold_jid ",paste0(hold,collapse=","))
  }
  if(array!=""){
      array=paste0(" -t 1-",array)
  }

  bps=""
  if(bypass){
    bps=" -P crag7day "
  }
  exec_code=paste("qsub -V -N ",job,array,bps,paste0(" -l h_rt=",time),
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


job_validator=function(
  job="",
  time=10,
  verbose=FALSE,
  threads=3
){

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

make_unique_id=function(
  name,
  id=sample(1:100000000000000,1,replace=FALSE),
  sep="_"
){
  options(scipen = 999)
  unique_name=paste0(name,"_",id)
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
  #Read .bashrc to import all envriomental variables
  error=system(paste0(". $HOME/.bashrc;",exec_code))
  return(error)
}




#' Build RData object to inherit in main_function
#' 
#' @param objects List of objects to save in RData object
#' @param job_id Job Identifier
#' @param output_dir Path to output directory 
#' @export


build_self=function(){   
      .base.env=parent.frame()     
      .this.env=environment()
      append_env(to=.this.env,from=.base.env)
      self_file=paste0(env_dir,"/",job_id,".self.RData")
      saveRDS(.this.env,file = self_file)
      append_env(to=.base.env,from=.this.env)
}




      



#' Build execution innit for function to execute
#' 
#' @param env List of objects to save in RData object
#' @export

build_main=function(){
  .base.env=parent.frame()
  .this.env=environment()
  append_env(to=.this.env,from=.base.env)
  main_file=paste0(env_dir,"/",job_id,".main.RData")
  saveRDS(.base.env$steps,file=main_file)
  append_env(to=.base.env,from=.this.env)
}


#' Build execution innit for function to execute
#' 
#' @param .env List of objects to save in RData object
#' @export

wait_scheduler=function(){
    .base.env=parent.frame()
    .this.env=environment()
    append_env(to=.this.env,from=.base.env)
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

read_main=function(){
    .base.env=parent.frame()
    .this.env=environment()
    append_env(to=.this.env,from=.base.env)

 
    ### Reads mains and updates values for enviroments with data
      main.envs=lapply(1:n_inputs,function(n){
          main.env=readRDS(main.envs[[n]]$main_file)}
    )
  
    append_env(to=.base.env,from=.this.env)

}





#' Build execution innit for function to execute
#' 
#' @param .env List of objects to save in RData object
#' @export



build_exec_innit=function(
){
        .base.env=parent.frame()
        .this.env=environment()
        append_env(to=.this.env,from=.base.env)

        ### Use SGE TASK ID if mode is set to batch otherwise use value
        
        if(mode=="local"){
              exec_code=paste0("Rscript -e \" invisible(lapply(1:",n_inputs,
              ",FUN=function(select){",ns,"::",fn,"(inherit=\\\"",
              self_file,"\\\",select=select)}))\"")
        }else if(mode=="local_parallel"){
              exec_code=paste0("Rscript -e \" invisible(parallel::mclapply(1:",n_inputs,
              ",FUN=function(select){",ns,"::",fn,"(inherit=\\\"",
              self_file,"\\\",select=select)},mc.cores=",threads-1,"))\"")
        }else if(mode=="batch"){
              exec_code=paste0("Rscript -e \" invisible(",
              ns,"::",fn,"(inherit=\\\"",self_file,
              "\\\",select=$SGE_TASK_ID))\"")
        }else{
          stop("Unkown mode type")
        }

        append_env(to=.base.env,from=.this.env)
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
  .this.env=environment()
  append_env(to=.this.env,from=.env)
  .env$exec_code<- paste0(.env$exec_code," && rm ",paste(input))
}







#' Build execution innit for batch
#' 
#' @param .env Inherit envoment variables
#' @export



build_batch_exec_innit=function(){ 
    .base.env=parent.frame()
    .this.env=environment()
    append_env(to=.this.env,from=.base.env)

    batch_code=build_job_exec(
      job=job_id,time=time,ram=ram,
      threads=threads,wd=getwd(),
      output_dir=batch_dir,
      bypass=bypass,
      hold=hold,array=n_inputs
    )

    exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,
    ";",exec_code,"'|",batch_code)

    append_env(to=.batch.env,from=.this.env)
}


#' Set vars for self run
#' 
#' @param env Inherit current enviroment.
#' 
#' 
#' @export

set_self=function(){
    
    .base.env=parent.frame()
    .this.env=environment()
    append_env(to=.this.env,from=.base.env)
    error=0
    self_file=""
    main_code=""
    exec_code=""
    batch_code=""
    append_env(to=.base.env,from=.this.env)
}




#' Run job from batch constructor
#' 
#' @param env Inherit current enviroment.
#' 
#' 
#' @export

run_self=function(){  
    .base.env=parent.frame()
    .this.env=environment()
    append_env(to=.this.env,from=.base.env$self_env)
    

    ### Initiate 
    set_self()
    ### Create RData object to inherit vars
    build_self()
    ### Create caller function
    build_exec_innit()

    ### UPDATE CALLER FUNCTION IF SCHEDULER IS USED
    if(mode=="batch"){
        build_batch_exec_innit()
    }

    ### PRODUCE VERBOSE
    if(verbose){
          print_verbose(job=job_id,
            arg=as.list(.this.env),
            exec_code=exec_code
          )
    }

    ### EXECUTE JOB
    error=execute_job(exec_code=exec_code)

    ### RETURN ERROR MESSAGE
    if(error!=0){
        stop(err_msg)
    }

    ### WAIT FOR SCHEDULER TO FINISH
    if(mode=="batch"){
      if(wait){
         wait_scheduler()
         read_main()
      }
    }else{
      read_main()
    }

    return(main.envs)
 
  
}


#' Run job from batch constructor
#' 
#' @param env Inherit current enviroment.
#' 
#' 
#' @export


 run_job=function(){
      .base.env=parent.frame()
      .this.env=environment()
      append_env(to=.this.env,from=.base.env)
    
      if(clean){
        build_clean_exec()
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

      .base.env$steps[[fn_id]] <- .this.env
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

#' Function to consolidate the types within main function
#' 
#' @param .env Environment
#' 
#' 
#' @export

consolidate_type<-function(){
      .base.env=parent.frame()
      .this.env=environment()
      append_env(to=.this.env,from=.base.env)

      ## WE LOOP THROUGH ALL VARIABLES FOR MAIN FUNCTION
      for(var in fn_vars){
        var_value=get(var)
        ### CHECK VARIABLE TYPE
        var_type=typeof(var_value)
        ### IF VARIABLE TYPE IS CHARACTER
        if(var_type=="character"){
          ## WE CHECK IF VARIABLE CONTAINS A PATH
          if(tryCatch({file.exists(var_value)},
            error=function(e){
              FALSE
          })){
            var_value=normalizePath(var_value)
            ## IF FILE EXISTS LOCALLY WE CREATE A SYMLINK IN THE 
            ## TEMP DIRECTORY FOR EACH VARIABLE
            var_dir=set_dir(dir=ln_dir,name=var)
            system(paste("ln -fs ",var_value, var_dir)," > /dev/null 2>&1 ")
            ### UPDATE THE VARIABLE TO THE SYMLINK
            .this.env[[var]]=paste0(var_dir,"/",basename(var_value))

          }else{
            ## CHECK MISSING CASES
            ## CHECK IF REMOTE NODE IS GIVEN
            if(remote){

                ### CHECK IF VARIABLE PATH EXIST IN REMOTE

                check=suppressWarnings(system(paste(
                "sshpass -f ",password,
                " ssh ",paste0(user,
                  ifelse(!is.null(user),"@",""),node),
                  "\" realpath -e ",var_value,"\""),intern=TRUE
                ))

                if(length(check)!=0){
                  var_dir=set_dir(dir=rmt_dir,name=var)
                  ### COPY REMOTE FILE TO LOCAL TMP DIR IF REMOTE FILE EXISTS
                  system(paste(
                  "sshpass -f ",password,
                  " ssh ",paste0(user,
                    ifelse(!is.null(user),"@",""),node),
                    "\" cp -rn ",check," -t ",
                    var_dir, "\"")
                  )                
                  ### UPDATE THE VARIABLE TO THE REMOTE FILE
                  .this.env[[var]]=paste0(var_dir,"/",basename(check))
                }
              }
            }
          }
        }
  append_env(to=.base.env,from=.this.env)  
}









#' Set default enviromental variables based on input for all
#' functions in package.
#' 
#' @param envir Inherit current enviroment.
#' @param id Input identifier. If not given name of input will be used.
#' @param name Name of output file directory
#' @export

set_env_vars=function(
  vars=NULL,
  fn=NULL,
  fn_id=NULL,
  fn_vars=NULL,
  output_dir=".",
  tmp_dir=NULL,
  ln_dir=NULL,
  rmt_dir=NULL,
  env_dir=NULL,
  batch_dir=NULL,
  output_name=NULL,
  verbose=FALSE,
  bgzip_index=FALSE,
  index_format="tbi",
  license_dir=build_default_license_list()$dir,
  clean=FALSE,
  batch_config=build_default_preprocess_config(),
  threads=1,
  ram=4,
  ns="ULPwgs",
  mode="local",
  time="48:0:0",
  get=FALSE,
  bypass=FALSE,
  node=NULL,
  user=NULL,
  password=NULL,
  sheet=NULL,
  inherit=NULL,
  select=NULL,
  err_msg=NULL,
  remote=FALSE,
  executor_id=NULL,
  wait=FALSE,
  hold=NULL
){  

    ## START WITH BASE.ENV VARIABLES
    .base.env=parent.frame()
    .this.env=environment()
    append_env(to=.this.env,from= .base.env)
    
    ### ONE JOB FOR TOP FUNCTION
    .base.env$n_jobs <- 1

    ### IF WE INHERIT A RDS FILE READ AND APPEND
    ### WE SAVE RDS FILES WHEN SUBMITTING JOBS TO SCHEDULER OR JOBS RUN LOCALLY IN PARALLEL
    ### SELECT VARIABLE DEFINES WHICH ENV WE ARE RUNNING
    ### WE APPEND THIS ENV AND STOP
    if(!is.null(inherit)){
        if(!is.environment(inherit)){
          inherit <-readRDS(file=inherit)
        }
        append_env(to=.this.env,inherit$main.envs[[select]])
        .base.env$self.envs <- .this.env
        return()
    }
    
    #### IF FUNCTION ID IS NOT GIVE WE USE FN NAME AS ID
    #### OTHERWISE WE APPEND FUNCTION ID
    if(is.null(fn_id)){
      fn_id<-fn
    }else{
      fn_id<-paste0(fn,".",fn_id)
    }

    ### IF EXECTUCTOR ID IS NOT GIVEN WE CREATE AN UNIQUE NAME USING THE FUNCTION ID
    if(is.null(executor_id)){
      executor_id <- make_unique_id(fn)
    }

    ### CHECK IF FILES PATH CAN BE LOCATED REMOTELY
    if(!is.null(node)){
      remote=TRUE
    }

    ### CREATE WORK DIRECTORIES
    set_work_dir()
  
    ### CHECK IF VARIABLE SHEET IS GIVEN
    ### IF SHEET IS GIVEN WE ASSIGN THE VARIABLES AND QUIT
    if (!is.null(sheet)){
        ### READ VARIABLE SHEET
        read_sheet()
    }else{
    ### SET CREATE A SHEET FROM THE FUNCTION VARIABLES
        sheet=as.data.frame(as.list(.base.env)[fn_vars])
    }

    ascertain_sheet()

    ### WE ASSIGN A TASK ID TO THE CALLER FUNCTION

    task_id <- make_unique_id(fn)

    job_id <- build_job(
      executor_id=executor_id,
      task_id=task_id
    )
    
    ### WE CREATE AN ERROR MESSAGE TO TRACK WHERE JOB FAILS
    if(is.null(err_msg)){
        err_msg <- paste0("CRITICAL ERROR: ",fn," (",job_id,") "," -> ")
    }

    ### FROM THE PARENT ENVIRONMENT WE CREATE A NEW ENVIRONMENT FOR EACH INPUT
    set_main_env()
    
    ### 
    .base.env$self.envs<-.this.env

  
}

#' Set working directory to use
#' 
#' @param output_dir Environment
#' @export


set_work_dir=function(){
      .base.env=parent.frame()
      .this.env=environment()
      append_env(to=.this.env,from=.base.env)
  
      ### CREATE MAIN WORKING DIRECTORY
      out_file_dir <- set_dir(
          dir=output_dir
      )


      ### CREATE TMP DIRECTORY
      ### WE WILL STORE TMP FILES FOR ALL FUNCTIONS HERE 
      if(is.null(tmp_dir)){
          tmp_dir <- set_dir(
            dir=out_file_dir,
            name="tmp"
        )
      }
      
      ### WITHIN TMP DIRECTORY CREATE DIRECTORY TO STORE SYMLINK OF FILES
      ### WE WILL STORE SYMLINK FILES FOR LOCAL FILES
      if(is.null(ln_dir)){
          ln_dir <- set_dir(
            dir=tmp_dir,
            name="ln"
        )
      }

      ### WITHIN TMP DIRECTORY CREATE DIRECTORY TO STORE REMOTE FILES
      ### WE WILL STORE REMOTE DOWNLOAD FILES HERE
      if(is.null(rmt_dir)){
          rmt_dir <- set_dir(
            dir=tmp_dir,
            name="rmt"
        )
      }

      ### CREATE DIRECTORY TO STORE ENVIROMENT DATA
      ### WE WILL STORE R ENVIRONMENT DATA HERE
      if(is.null(env_dir)){
            env_dir<- set_dir(
              dir=out_file_dir,
              name="env"
          )
        }

      ### CREATE DIRECTORY TO BATCH DATA
      ### WE WILL STORE SCHEDULER DATA HERE
      if(is.null(batch_dir)){
          batch_dir<- set_dir(
            dir=out_file_dir,
            name="batch"
        )
      }

    ### WE APPEND NEW VARIABLES BACK TO THE MAIN FUNCTION
    append_env(to=.base.env,from=.this.env)
}

#' Set main enviroment inputs
#'
#' @export

set_inputs<-function(){
       .base.env=parent.frame()
       .this.env=environment()
        append_env(to=.this.env,from=.base.env)

        n_inputs <- 1
        inputs<-NULL
        inputs_id <- output_name
        inputs_ext <- NULL

      ### IF VARS ARGUMENT IS ASSIGNED WE PERMIT MUTLIPLE INPUTS FOR SPECIFIC VARIABLE
        if(!is.null(vars)){
          inputs <- get(vars)
          n_inputs <- length(inputs)
         
        }

        append_env(to=.base.env,from=.this.env)
}


#' Ascertain values in variable sheet
#'
#' @export

 ascertain_sheet<-function(){
      .base.env=parent.frame()
      .this.env=environment()
      append_env(to=.this.env,from=.base.env)
      n_total<-nrow(sheet)
      sheet=sheet %>% dplyr::distinct()
      n_inputs<-nrow(sheet)
      n_vars<-ncol(sheet)
      n_dup=n_total-n_inputs
      if(n_dup>0){
        warning(paste0(n_dup, " were duplicated in sheet"))
      }
      .base.env=parent.frame()
      .this.env=environment()
      append_env(to=.base.env,from=.this.env)
}



#' Set steps enviroment for use
#' 
#' @param .env Environment
#' @export


run=function(){
  append_env()
  if(is.null(select)){
    run_self()
  }else{
    consolidate_type()
    FUN()
    build_main()
  }
}



#' Set steps enviroment for use
#' 
#' @param .env Environment
#' @export

launch=function(){
  
      append_env()
      reports=run()
      return(reports)
  }


#' Set steps enviroment for use
#' 
#' @param .env Environment
#' @export


read_sheet=function(){
    .base.env=parent.frame()
    .this.env=environment()
    append_env(to=.this.env,from=.base.env)
    sheet=read.delim(sheet,header=TRUE)
    append_env(to=.base.env,from=.this.env)
}





#' Set steps enviroment for use
#' 
#' @export


set_main_env=function(){
    .base.env=parent.frame()
    .this.env=environment()
    append_env(to=.this.env,from=.base.env)
    ### WE CREATE AN ENVIRONMENT FOR EACH INPUT VALUE
    main.envs=parallel::mclapply(
        1:n_inputs,
        function(row){
        .this.env=environment()
        append_env(to=.this.env,from=.base.env)
        
        ### ASSIGN VARS IN SHEET TO ENVIROMENT
        for (col in n_vars){
          .this.env[[names(sheet)[col]]]<-sheet[row,col]
        }

        ### WE CREATE A TASK ID FOR EACH JOB
        executor_id=task_id
        task_id <- make_unique_id(fn)
        job_id <- build_job(
          executor_id=executor_id,
          task_id=task_id
        )

        ## WE TRACE ERROR MESSAGE
        err_msg <- paste0(err_msg ,fn," (",job_id,") "," -> ")
        
        ### TO TRACK EACH JOB WE SAVE EACH ENVIROMENT IN RDS FORMAT
        build_main()

        return(.this.env)
        },
        mc.cores=parallel::detectCores()
    )
    append_env(to=.base.env,from=.this.env)
}



#' Set steps enviroment for use
#' 
#' @param envir Environment
#' @export

set_main=function(){     
      .base.env=parent.frame()
      .this.env=environment()
      append_env(to=.this.env,from=.base.env)
      exec_code=""
      error=0
      out_file=""
      steps=list()
      out_files=list()
      append_env(to=.base.env,from=.this.env)
}



#' Append two enviroments
#' 
#' @param to Enviroment to append to
#' @param from Enviroment to append from
#' @export

append_env = function(to=environment(), from=parent.frame()) {
      
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
        my_id=rep(ids,length(inputs))
      }else if(!all(sapply(inputs,typeof)=="character")){
        my_id <- ifelse(!is.null(names(inputs)),names(inputs),seq(1,length(inputs)))
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

#' Caller function within routine
#' 
#' @param FUN Function to call
#' @export

 call_function=function(
    FUN=NULL,
    ...
 ){
       if(is.null(FUN)){
          stop("Please define a function to run")
       }

       .base.env=parent.frame()
       ### ADD OTHER VARIABLES TO BASE ENV
       list2env(x=list(...),envir=.base.env)
       ## CREATE FUNCTION VARIABLE NAMES
       .base.env$fn_vars=names(.base.env)[!grepl("\\.",names(.base.env))]
       ## IF NOT SET UP BY THE USER WE WILL GET THE MAIN FUNCTION NAME
       if(is.null(.base.env$fn)){
          ## GET CALLER FUNCTION NAME
          fn <- sub(".*::","",sub("\\(.*","",
            paste0(deparse(sys.calls()[[sys.nframe()-1]]),collapse=","))
          )
       }

      ## WE WILL DEFINE THE ENVIROMENTAL VARIABLES
      set_env_vars()
      print(as.list(.base.env))
      ## WE WILL LAUNCH THE MAIN FUNCTION
      launch()
    }



