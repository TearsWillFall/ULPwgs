#' Build job for SGE executor
#' Job constructor only uses the top level Executor ID 
#' and the bottom level task ID.
#' 
#' Middle level task are omited in Job ID for simplicity
#'
#' @param executor Job EXECUTOR ID
#' @param task Task ID
#' @export

build_job=function(
  parent_id=NULL,
  child_id=NULL
){

  options(scipen = 999)
  if(is.null(parent_id)){
    stop("parent_id argument is required to allocate an id for job")
  }
  parent_id=paste0("parent.",parent_id)
  if(!is.null(child_id)){
    job=lapply(child_id,FUN=function(id){
      child_id=paste0("child.",id)
      job=paste0(c(parent_id,child_id),collapse=".")
    })
  }else{
    job=parent_id
  }

  return(job)
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



#' Build a job report data structure
#' 
#' Data structure for job_report:
#' job_id job_order out_files
#'
#'
#' @param job_id Job ID
#' @param parent_id Task executor ID
#' @param child_id Task ID
#' @param exec_code Execution code
#' @param input Input arguments for job
#' @param job_order Job execution order
#' @param out_file_fir Output dir for 
#' @param out_files Output files for job
#' @export

build_job_report=function(job_id="job_1",
  parent_id="",
  child_id="",
  exec_code="",
  job_order=1,
  input_args="",
  out_file_dir="",
  out_files=list(file="file")
){
  options(scipen = 999)
  job_report=list(job_id=job_id,parent_id=parent_id,child_id=child_id,
  job_order=job_order,input_args=input_args,exec_code=exec_code,
  out_file_dir=out_file_dir,out_files=out_files)
  return(job_report)
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
  sep="."
){
  options(scipen = 999)
  unique_name=paste0(name,sep,id)
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




#' Build execution innit for function to execute
#' 
#' @param .env List of objects to save in RData object
#' @export

wait_scheduler=function(){
    append_env(to=environment(),from=parent.frame())
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


#' Set main enviroment inputs
#'
#' @export

set_inputs<-function(){
        append_env(to=environment(),from=parent.frame())
        n_inputs <- 1
        inputs<-NULL
        inputs_id <- output_name
        inputs_ext <- NULL

      ### IF VARS ARGUMENT IS ASSIGNED WE PERMIT MUTLIPLE INPUTS FOR SPECIFIC VARIABLE
        if(!is.null(vars)){
          inputs <- get(vars)
          n_inputs <- length(inputs)
         
        }

        append_env(from=environment(),to=parent.frame())
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


callFUN<-function(){
  UseMethod("callFUN")
}

callFUN.call<-function(
  ...,
  args
  ){
    .this.env=environment()
    .base.env=parent.frame()
    ### ADD OTHER VARIABLES TO BASE ENV
    list2env(x=list(...),envir=.base.env)

    ## GET VARIABLE NAMES
    fn_vars=names(.base.env)[!grepl("\\.|def_args|FUN",names(.base.env))]

    append_env(to=.this.env,from=.base.env)

    ## WE SET THE ENVIROMENT TO CALLER FUNCTION
    callFUN.setEnv()

    ## WE WILL RUN SELF OR A PROCESS
    if(is.null(select)){
    
      callFUN.runSelf()
  
    }else{
      callFUN.runProcess()
    }

    return(.this.env)
}


callFUN.setEnv<-function(){

    append_env(to=environment(),from=parent.frame())
    
    ## WE CHECK IF WE ARE INHERITING A PARENT ENVIROMENT
    ## ELSE WE CREATE A PARENT ENVIROMENT
    ## WHEN INHERITING WE DON'T CHECK ARG TYPES AGAIN

    if(exists("env")){
      if(!is.environment(env)){
        env <-readRDS(file=env)
        append_env(to=environment(),from=env$child.envs[[select]])
      }else{
        append_env(to=environment(),from=env)
      }
      env<-NULL

    }else{
      self<-TRUE

      ## WE VALIDATE USER DEFINED VARIABLE FOR PARENT FUNCTION
      callFUN.checkArgs()

      ### WE BUILD THE PARENT ENVIRONMENT
      callFUN.buildParent()

      ### WE SET THE DUMPSTER WHERE TO PUT CHILDREN INFO
      callFUN.dumpInfo()

      ### FROM THE PARENT ENVIRONMENT WE CREATE A NEW ENVIRONMENT FOR EACH UNIQUE INPUT
      child.envs=parallel::mclapply(
          1:n_inputs,
          function(row,.env){

            append_env(to=environment(),from=.env)

            ### WE BUILD THE CHILD ENV
            callFUN.buildChild()

            ### WE DUMP CHILDREN INFO
            callFUN.dumpInfo()

            ### WE WRITE CHILDREN ENV
            callFUN.writeEnv()

            return(environment())
            },
            .env=environment(),
            mc.cores=ifelse(
              (parallel::detectCores()-1)<1,
              1,
              (parallel::detectCores()-1)
        )
      )
    }

    append_env(from=environment(),to=parent.frame())

}



callFUN.checkArgs<-function(){
   .this.env=environment()
    append_env(to=.this.env,from=parent.frame())

    if(exists("def_args")){
      types=def_args$types
      required=def_args$required

      for (arg in names(types)){
        if(!exists(arg)){
          stop(paste0("Variable : ",arg,
          " ( type : ",types[[arg]]," ) -> Value: Not defined. Define a value"))
        }

        if(is.null(.this.env[[arg]])&required[[arg]]){
          stop(paste0("Variable : ",arg,
          " ( type : ",types[[arg]]," ) -> Value: NULL (type: NULL). Define a non-NULL value"))
        }

        if(typeof(.this.env[[arg]])!=types[arg]&required[[arg]]){
          stop(paste0("Variable : ",arg,
          " ( type : ",types[[arg]]," ) -> Value: ",.this.env[[arg]],
          " ( type : ",typeof(.this.env[[arg]])," ). Invalid type"))
        }
      } 
  }
    
}



callFUN.runSelf=function(){
  append_env(to=environment(),from=parent.frame())
 

  ### WE WRITE PARENT ENVIROMENT
  callFUN.writeEnv()

  ### CREATE CALLER
  callFUN.buildCall()

  ### RUN CALL
  callFUN.runCall()
  
  ### WAIT FOR SCHEDULER TO FINISH
  if(mode=="batch"){
    if(wait){
        wait_scheduler()
    }
  }

  ### READ CHILD ENVIRONMENTS
  callFUN.readEnv()

  ### REMOVE WORK DIRECTORIES IF PROCESS FINISHED
  callFUN.rmDir()
    
  append_env(from=environment(),to=parent.frame())
}


callFUN.rmDir=function(){
    append_env(to=environment(),from=parent.frame())
    if(!self){
    if(preserve=="partial"){
        if(error!=0){
            system(paste0("rm -rf ",parent_dir))
          }
    }else if (preserve=="none"){
        system(paste0("rm -rf ",parent_dir))
      }
    }

  append_env(from=environment(),to=parent.frame())
}


callFUN.runProcess=function(){

  append_env(to=environment(),from=parent.frame())
  
  ### GET ALL DATA BEFORE RUNNING ENV
  callFUN.collectData()

  ### CREATE CALLER
  callFUN.buildCall()

  ### WE RUN THE ENVIROMENT
  callFUN.runCall()
  
  ### WE WRITE RESULTS FOR CHILD TO RDS
  callFUN.writeEnv()

  append_env(from=environment(),to=parent.frame())

}

callFUN.runCall=function(FUN=NULL){
  append_env(to=environment(),from=parent.frame())
  ### PRODUCE VERBOSE
  if(verbose){
        print_verbose(job=job_id,
          arg=as.list(environment())[names(environment()) %in% dump_names],
          exec_code=exec_code
        )
  }

  #Read .bashrc to import all envriomental variables
  error=system(paste0(". $HOME/.bashrc;",exec_code),wait=wait)
  ### RETURN ERROR MESSAGE
  if(error!=0){
      stop(err_msg)
  }

  append_env(from=environment(),to=parent.frame())
}



callFUN.collectData.createLink=function(){
    append_env(to=environment(),from=parent.frame())
    var_value=normalizePath(var_value)
    ## IF FILE EXISTS LOCALLY WE CREATE A SYMLINK IN THE 
    ## TEMP DIRECTORY FOR EACH VARIABLE
    var_dir=set_dir(dir=ln_dir,name=var)
    system(paste("ln -fs ",var_value, var_dir," > /dev/null 2>&1 "))
    ### UPDATE THE VARIABLE TO THE SYMLINK
    append_env(from=environment(),to=parent.frame())
}


callFUN.collectData.remoteCheck=function(){
    append_env(to=environment(),from=parent.frame())
    check=suppressWarnings(system(paste(
    "sshpass -f ",password,
    " ssh ",paste0(user,
      ifelse(!is.null(user),"@",""),node),
      "\" realpath -e ",var_value,"\" "),intern=TRUE
    ))
    append_env(from=environment(),to=parent.frame())
}

callFUN.collectData.remoteGet=function(){
  append_env(to=environment(),from=parent.frame())
  var_dir=set_dir(dir=rmt_dir,name=var)
  ### COPY REMOTE FILE TO LOCAL TMP DIR IF REMOTE FILE EXISTS
  system(paste(
  "sshpass -f ",password,
  " ssh ",paste0(user,
    ifelse(!is.null(user),"@",""),node),
    "\" cp -rn ",check," -t ",
    var_dir, "\"")
  )
  append_env(from=environment(),to=parent.frame())
}

callFUN.collectData=function(){
    ## WE LOOP THROUGH ALL VARIABLES FOR MAIN FUNCTION
    .base.env=parent.frame()
    append_env(to=environment(),from=.base.env)
    parallel::mclapply(
      fn_vars,
      FUN=function(
        var,
        env
      ){
        append_env(to=environment(),from=env)
        var_value=get(var)
        ### CHECK VARIABLE TYPE
        var_type=typeof(var_value)
        ### IF VARIABLE TYPE IS CHARACTER
        if(var_type=="character"){
          ## WE CHECK IF VARIABLE CONTAINS A PATH LOCALLY
          if(file.exists(var_value)){
            callFUN.collectData.createLink()
          }else{
            ## CHECK IF REMOTE NODE IS GIVEN
            if(is.null(node)){
              return()
            }
          

            ### CHECK IF VARIABLE PATH EXIST IN REMOTE
            callFUN.collectData.remoteCheck()

            ### WE STOP FUNCTION IF FILE IS NOT FOUND REMOTELY
            ### WE STOP FUNCTION IF FILE IS NOT FOUND REMOTELY
            if(length(check)==0){
              return()
            }

            ### GET REMOTE PATH
            callFUN.collectData.remoteGet()
            var_value=paste0(var_dir,"/",basename(check))     
      
            ### We CREATE A SYMLINK TO THE DATA
            callFUN.collectData.createLink()
          }
          env[[var]]=normalizePath(paste0(var_dir,"/",basename(var_value)))

          }
      },
      env=.base.env,
      mc.cores=ifelse(
                (parallel::detectCores()-1)<1,
                1,
                (parallel::detectCores()-1)
      )
    )
}




callFUN.buildCall=function(){
    append_env(to=environment(),from=parent.frame())
    ### Use SGE TASK ID if mode is set to batch otherwise use value
   
    callFUN.setCall()
  
    ### CHECK IF WE ARE IN A CHILD ENVIROMENT
    if(!exists("child_id")){

      if(self){
          exec_code=paste0("Rscript -e \"",
            ns,"::",fn,"(env=\\\"",
            env_file,"\\\")})\""
          )
      }
    
      ## WE ASSUME WE HAVE INFINITE CORES AND CAN RUN INFINITE JOBS 
      if(mode=="local"){
            ### LOCALLY WE ARE LIMITED IN NUMBER OF CORES
            cores=parallel::detectCores()-1
            ### WE ASSIGN A REASONABLE NUMBER OF JOBS FoR THE REQUESTED NUMBER OF THREADS 
            rjobs=floor(cores/threads)
            exec_code=paste0("Rscript -e \" invisible(parallel::mclapply(1:",n_inputs,
            ",FUN=function(select){",ns,"::",fn,"(env=\\\"",
            env_file,"\\\",select=select)},mc.cores=",rjobs,"))\"")
      
      }else if(mode=="batch"){
            cores=Inf
            rjobs=Inf
            ### IN BATCH WE ASSUME INFINITE RESOURCES 
            exec_code=paste0("Rscript -e \" invisible(",
            ns,"::",fn,"(env=\\\"",env_file,
            "\\\",select=$SGE_TASK_ID))\"")
          

            ### WE CREATE SCHEDULER SPECIFIC VARIABLES

            batch_code=paste(
                  "qsub -V ", 
            paste0("-N ",job_id),
            paste0(" -t 1-",n_inputs),
            paste0(" -l h_rt=",time),
            paste0(" -l mem=",ram,"G"),
            paste0(" -pe smp ",threads), 
            paste0(" -wd ",getwd()), 
            paste0(" -o ",batch_dir,"/",job_id,".std_out"),
            paste0(" -e ",batch_dir,"/",job_id,".std_error"))

            if(!is.null(hold)){
              batch_code=paste0(batch_code,paste0(" -hold_jid ",paste0(hold,collapse=",")))
            }

            if(bypass){
              batch_code=paste0(batch_code," -P crag7day ")
            }

            ### WE APPEND CONFIG EXEC CODE AND BATCH CODE DATA

            exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,
            ";",exec_code,"'|",batch_code)
      }else{
        stop("Unkown mode type")
      }

       
    }else{
      ### OF THIS IS A CHILD ENVIROMENT WE DIRECTLY RUN USER DEFINED FUNCTION
      FUN()
    }

    append_env(from=environment(),to=parent.frame())
}


callFUN.setCall=function(){
  append_env(to=environment(),from=parent.frame())
  cores=1
  rjobs=1
  error=0
  exec_code=NULL
  batch_code=NULL
  steps=list()
  out_files=list()
  append_env(from=environment(),to=parent.frame())
}


callFUN.writeEnv=function(){
  append_env(to=environment(),from=parent.frame())
  if(exists("child_id")){
      env_file=paste0(env_dir,"/",child_id,".child.RData")
      saveRDS(environment(),file= env_file)
  }else{
      self<-FALSE
      env_file=paste0(env_dir,"/",parent_id,".parent.RData")
      saveRDS(environment(),file = env_file)
  }
  append_env(from=environment(),to=parent.frame())
}


callFUN.readEnv=function(){
  append_env(to=environment(),from=parent.frame())
  child.envs=lapply(1:n_inputs,function(n){
        child.env=readRDS(child.envs[[n]]$env_file)}
  )
  append_env(from=environment(),to=parent.frame())
}



callFUN.buildId=function(){
  append_env(to=environment(),from=parent.frame())

  if(!exists("job_id")){
    ### IF parent ID IS NOT GIVEN WE CREATE AN UNIQUE NAME USING THE FUNCTION ID
    if(is.null(parent_id)){
      parent_id <- make_unique_id(fn)
    }

    ### WE ASSIGN A JOB ID TO THE CALLER FUNCTION

    job_id <- build_job(
      parent_id=parent_id
    )

  }else{
    ### WE CREATE A JOB ID FOR EACH CHILD
    ### CHILDREN SHALL WORK!!

    child_id <- make_unique_id(fn)
    
    job_id <- build_job(
      parent_id=parent_id,
      child_id=child_id
    )
  }
  
  append_env(from=environment(),to=parent.frame())
}




callFUN.buildError<-function(){
  append_env(to=environment(),from=parent.frame())
  if(!is.null(err_msg)){
     err_msg <- paste0(err_msg ,fn," (",job_id,") "," -> ")
  }else{
     err_msg <- paste0("CRITICAL ERROR: ",fn," (",job_id,") "," -> ")
  }
  append_env(from=environment(),to=parent.frame())
}




callFUN.reassignArgs<-function(){
  append_env(to=environment(),from=parent.frame())
  ### REASSIGN VAR VALUES IN SHEET TO ENVIRONMENT
  for (col in 1:n_vars){
    assign(names(sheet)[col],sheet[row,col])
  }
  append_env(from=environment(),to=parent.frame())

}




callFUN.buildSheet=function(){
  append_env(to=environment(),from=parent.frame())
  ### CHECK IF VARIABLE SHEET IS GIVEN
  ### IF SHEET IS GIVEN WE ASSIGN THE VARIABLES AND QUIT
  if (!is.null(sheet)){
      ### READ VARIABLE SHEET
       sheet=read.delim(sheet,header=TRUE)
  }else{
  ### SET CREATE A SHEET FROM THE FUNCTION VARIABLES
      sheet=as.data.frame(as.list(parent.frame())[fn_vars])
  }

  callFUN.setSheet()

  append_env(from=environment(),to=parent.frame())

}

callFUN.setSheet=function(){
  append_env(to=environment(),from=parent.frame())
  n_total<-nrow(sheet)
  sheet=sheet %>% dplyr::distinct()
  n_inputs<-nrow(sheet)
  n_vars<-ncol(sheet)
  n_dup=n_total-n_inputs
  if(n_dup>0){
    warning(paste0(n_dup, " were duplicated in sheet"))
  }
  append_env(from=environment(),to=parent.frame())
}







callFUN.buildChild=function(){
  append_env(to=environment(),from=parent.frame())

  ### USING THE VAR SHEET WE REASSIGN THE VALUES OF EACH VARIABLE
  callFUN.reassignArgs()

  ### CREATE CHILD JOB ID
  callFUN.buildId()

  ## CREATE ERROR MESSAGE FOR EACH CHILD
  callFUN.buildError()

  append_env(from=environment(),to=parent.frame())

}



callFUN.buildParent=function(
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
  bypass=FALSE,
  preserve="partial",
  node=NULL,
  user=NULL,
  password=NULL,
  sheet=NULL,
  env=NULL,
  select=NULL,
  err_msg=NULL,
  complimentary=FALSE,
  parent_id=NULL,
  wait=FALSE,
  hold=NULL
){
    append_env(to=environment(),from=parent.frame())
    ### IF WE INHERIT A RDS FILE READ AND APPEND
    ### WE SAVE RDS FILES WHEN SUBMITTING JOBS TO SCHEDULER OR JOBS RUN LOCALLY IN PARALLEL
    ### SELECT VARIABLE DEFINES WHICH ENV WE ARE RUNNING
    ### WE APPEND THIS ENV AND STOP

    ### CREATE VARIABLES FOR THE ENVIRONMENT
    callFUN.buildSelf()
      
    ### CREATE JOB ID FOR THE PARENT FUNCTION
    callFUN.buildId()

    ### CREATE WORK DIRECTORIES
    callFUN.buildDirs()
    
    ### CREATE SHEET WITH VARIABLES
    callFUN.buildSheet()
    
    ### CREATE AN ERROR MESSAGE TO TRACK WHERE JOB FAILS
    callFUN.buildError()

    append_env(from=environment(),to=parent.frame())
}


callFUN.buildSelf=function(){
  append_env(to=environment(),from=parent.frame())
   ## IF NOT SET UP BY THE USER WE WILL GET THE MAIN FUNCTION NAME
    
  if(is.null(fn)){
    ## GET CALLER FUNCTION NAME
    fn <- sub(".*::","",sub("\\(.*","",
      paste0(deparse(sys.calls()[[sys.nframe()-4]]),collapse=","))
    )
  }

  #### IF FUNCTION ID IS NOT GIVE WE USE FN NAME AS ID
  #### OTHERWISE WE APPEND FUNCTION ID

  if(is.null(fn_id)){
    fn_id<-fn
  }else{
    fn_id<-paste0(fn,".",fn_id)
  }

  append_env(from=environment(),to=parent.frame())
}




callFUN.buildDirs=function(){
      append_env(
        to=environment(),
        from=parent.frame()
      )

      ### CREATE MAIN WORKING DIRECTORY
      out_file_dir <- set_dir(
          dir=output_dir
      )


      ### CREATE MAIN WORKING DIRECTORY
      parent_dir <- set_dir(
          dir=out_file_dir,
          name=parent_id
      )

      ### CREATE TMP DIRECTORY
      ### WE WILL STORE TMP FILES FOR ALL FUNCTIONS HERE 
      if(is.null(tmp_dir)){
          tmp_dir <- set_dir(
            dir=parent_dir,
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
      }else{
        ln_dir <- set_dir(
          dir=ln_dir,
          name=child_id
        )
      }

      ### WITHIN TMP DIRECTORY CREATE DIRECTORY TO STORE REMOTE FILES
      ### WE WILL STORE REMOTE DOWNLOAD FILES HERE
      if(is.null(rmt_dir)){
          rmt_dir <- set_dir(
            dir=tmp_dir,
            name="rmt"
        )
      }else{
          rmt_dir <- set_dir(
            dir=rmt_dir,
            name=child_id
        )
      }

      ### CREATE DIRECTORY TO STORE ENVIROMENT DATA
      ### WE WILL STORE R ENVIRONMENT DATA HERE
      if(is.null(env_dir)){
            env_dir<- set_dir(
              dir=parent_dir,
              name="env"
          )
      }else{
            env_dir<- set_dir(
              dir=env_dir,
              name=child_id
          )
      }

      ### CREATE DIRECTORY TO BATCH DATA
      ### WE WILL STORE SCHEDULER DATA HERE
      if(is.null(batch_dir)){
          batch_dir<- set_dir(
            dir=parent_dir,
            name="batch"
        )
      }else{
          batch_dir<- set_dir(
            dir=batch_dir,
            name=child_id
          )

      }

    ### WE APPEND NEW VARIABLES BACK TO THE MAIN FUNCTION
    append_env(from=environment(),to=parent.frame())
}



callFUN.dumpInfo<-function(){
  append_env(to=environment(),from=parent.frame())



  if(!exists("dump_file")){
    ### TRACING CHILDREN CAN BE DIFFICULT
    ### WE DEFINE THE VARIABLES THAT WILL HELP IDENTIFY CHILDREN IN JOB HIERARCHY

    dump_names=c("parent_id","child_id","child_order",fn_vars)
    
    ## WE DEFINE A FILE WERE WE WILL DUMP THIS INFO
    
    dump_file=paste0(parent_dir,"/",parent_id,".dump")
    ### WE CREATE AN ENVIRONMENT FOR EACH INPUT VALUE

    write.table(
          file=dump_file,
          x=as.data.frame(matrix(dump_names,nrow=1)),
          sep="\t",
          col.names=FALSE,
          row.names=FALSE,
          quote=FALSE
    )
  }else{

    ### WE APPEND THE CHILDREN INFO TO THE DUMP FILE

    write.table(
      file=dump_file,
      x=data.frame(
      parent_id=parent_id,
      child_id=child_id,
      child_order=row,
      sheet[row,]),
      sep="\t",
      append=TRUE,
      col.names=FALSE,
      row.names=FALSE,
      quote=FALSE
    )          
  }
  append_env(from=environment(),to=parent.frame())
}


dumpInfo.append=function(){
  append_env(to=environment(),from=parent.frame())
  
  
  
}










