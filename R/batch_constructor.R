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
    fn_vars=names(.base.env)[!grepl("\\.|args|FUN",names(.base.env))]

    append_env(to=.this.env,from=.base.env)

    ## WE SET THE ENVIROMENT TO CALLER FUNCTION
    callFUN.setEnv()

    ## WE WILL RUN SELF OR A PROCESS
    if(name_env!="child"){
      callFUN.runSelf()
    }else{
      callFUN.runProcess()
    }
    return(.this.env)
}


callFUN.setEnv<-function(
   
){

    append_env(to=environment(),from=parent.frame())
    
    ## WE CHECK IF WE ARE INHERITING A PARENT ENVIROMENT
    ## ELSE WE CREATE A PARENT ENVIROMENT
    ## WHEN INHERITING WE DON'T CHECK ARG TYPES AGAIN

    if(exists("env")){
      if(!is.environment(env)){
        env <-readRDS(file=env)
        if(exists("select")){
          env <- env$child.env[[select]]
        }
      }
      append_env(to=environment(),from=env)
      env<-NULL
    }else{
     
      args=appendList(args,build_default_variable_list())

      ### CREATE VARIABLES FOR THE ENVIRONMENT
      callFUN.buildSelf()
     
      ### WE WRITE PARENT ENVIROMENT
      callFUN.writeEnv()

    }

    append_env(from=environment(),to=parent.frame())

}




callFUN.createLink=function(){
    append_env(to=environment(),from=parent.frame())
    arg_value=normalizePath(arg_value)
    ## IF FILE EXISTS LOCALLY WE CREATE A SYMLINK IN THE 
    ## TEMP DIRECTORY FOR EACH VARIABLE
    arg_dir=set_dir(dir=ln_dir,name=arg)
    system(paste("ln -fs ",arg_value, arg_dir," > /dev/null 2>&1 "))
    ### UPDATE THE VARIABLE TO THE SYMLINK
    append_env(from=environment(),to=parent.frame())
}



callFUN.remoteCheck=function(){
  append_env(to=environment(),from=parent.frame())
  arg_dir=set_dir(dir=rmt_dir,name=arg)
  ### COPY REMOTE FILE TO LOCAL TMP DIR IF REMOTE FILE EXISTS
  check=system(paste(
  "sshpass -e ssh -Y ",ip_address,
    "\" realpath -e ",arg_value,"\""),intern=TRUE
  )  
  append_env(from=environment(),to=parent.frame())
}


callFUN.remoteGet=function(){
  append_env(to=environment(),from=parent.frame())
  arg_dir=set_dir(dir=rmt_dir,name=arg)
  ### COPY REMOTE FILE TO LOCAL TMP DIR IF REMOTE FILE EXISTS
  system(paste(
  "sshpass -e ssh -Y ",ip_address,
    "\" cp -r ",arg_value," -t ",
    arg_dir, "\"")
  )
  arg_value=paste0(arg_dir,"/",basename(arg_value))   
  append_env(from=environment(),to=parent.frame())
}





callFUN.checkTypes<-function(){
   .this.env=environment()
   .base.env=parent.frame()
    append_env(to=.this.env,from=.base.env)

    if(exists("args")){
      types=args$types
      subtypes=args$subtypes
      required=args$required

      for (arg in names(types)){
        arg_value=.this.env[[arg]]
        arg_type=types[[arg]]
        arg_subtype=subtypes[[arg]]
        arg_required=required[[arg]]

        if(!exists(arg)){
          mssg=paste0("Variable : ",arg,
          " ( type : ",arg_type,
          " ) -> Value: Not defined . Define a value")
          callFUN.callError()
        }

        if(typeof(arg_value)!=arg_type){
          mssg=paste0("Variable : ",arg,
          " ( type : ",arg_type," ) -> Value: ",arg_value,
          " ( type : ",typeof(arg_value),
          " ) . Invalid type")
          callFUN.callError()
        }
      }
    }
}
     
  


callFUN.checkSubtypes=function(){
    .this.env=environment()
    append_env(to=environment(),from=parent.frame())
 

    if(exists("args")){
      types=args$types
      subtypes=args$subtypes
      required=args$required

      for (arg in names(subtypes[!is.null(subtypes)])){
        arg_value=.this.env[[arg]]
        arg_type=types[[arg]]
        arg_subtype=subtypes[[arg]]
        arg_required=required[[arg]]

        if(!is.null(arg_subtype)){
           if(arg_subtype=="path"){
              if(!is.null(arg_value)){
                if(exists("remote")){
                  if(!any(remote %in% arg)){
                    print(arg_value)
                    if(!file.exists(arg_value)){
                       mssg=paste0("Variable : ",arg,
                              " ( type : ",arg_type," ) [ subtype : path ] -> Value: ",arg_value,
                              " ( type : ",typeof(arg_value),
                              " ) [ subtype : NULL ] . Path doesn't exist locally")
                           callFUN.callError()
                           
                    }
                  }else{
                    callFUN.remoteCheck()
                    if(length(check)==0){
                      mssg=paste0("Variable : ",arg,
                              " ( type : ",arg_type," ) [ subtype : path ] -> Value: ",arg_value,
                              " ( type : ",typeof(arg_value),
                              " ) [ subtype : NULL ] . Path doesn't exist remotely")
                      callFUN.callError()
                    }
                    
                  }
                }
            }
          }
        }
      }
    }
  }


callFUN.verbose=function(){
  append_env(to=environment(),from=parent.frame())
  if(name_env=="self"){
    if(verbose){
          print_verbose(job=env_id,
            arg=as.list(environment())[names(environment()) %in% fn_vars],
            exec_code=exec_code
          )
    }
  }else{
    if(verbose){
          print_verbose(job=env_id,
            arg=as.list(environment())[names(environment()) %in% dump_names],
            exec_code=exec_code
          )
    }
  }
}

    
callFUN.runSelf=function(){
  append_env(to=environment(),from=parent.frame())

  ### CREATE CALLER
  callFUN.buildCall()


  ### RUN CALL
  callFUN.runCall()
  
  ### WAIT FOR SCHEDULER TO FINISH
  if(rmode=="batch"){
    if(await){
        await_scheduler()
    }
  }

  ### READ CHILD ENVIRONMENTS
  callFUN.readEnv()
    
  append_env(from=environment(),to=parent.frame())
}




callFUN.rmDir=function(){
    append_env(to=environment(),from=parent.frame())
    if(!self){
      if(preserve=="partial"){
          if(error!=0){
              system(paste0("rm -rf ",work_dir))
            }
      }else if (preserve=="none"){
          system(paste0("rm -rf ",work_dir))
        }
    }

  append_env(from=environment(),to=parent.frame())
}


callFUN.runProcess=function(){

  append_env(to=environment(),from=parent.frame())
  
  ### CREATE CALLER
  callFUN.buildCall()

  ### WE RUN THE ENVIROMENT
  callFUN.runCall()
  
  ### WE WRITE RESULTS FOR CHILD TO RDS
  callFUN.writeEnv()

  ### WE MOVE DATA TO OUT_FILE_DIR DIRECTORY
  callFUN.moveData()

  append_env(from=environment(),to=parent.frame())

}



callFUN.moveData=function(){
  append_env(to=environment(),from=parent.frame())

  if(!exists("remote")){
    if(length(out_files)!=0){
      check=system(paste("mv ",ifelse(overwrite," -n ",""),paste0(out_dir,"/*"),out_file_dir, " 2> /dev/null ; echo $?"),intern=TRUE)
    }
  }else{
     check=system(
      paste("sshpass -e ssh -Y ",ip_address, " \" mv ",
      ifelse(overwrite," -n ",""),
      paste0(out_dir,"/*"),
      out_file_dir,";echo 0\"" ),intern=TRUE
    )
  }
  if(check!=0){
    mssg=paste0("Files exist. To overwrite the files in `output_dir` set variable `overwrite` to : TRUE ")
    callFUN.callError()
  }

  mssg=paste0("Files moved from : ",work_dir," to ", ifelse(exists("remote"),"remote", "local")," directory  `output_dir` : `",out_file_dir,"'")
  append_env(from=environment(),to=parent.frame())
}



callFUN.runCall=function(FUN=NULL){
  append_env(to=environment(),from=parent.frame())

  callFUN.verbose()

  #Read .bashrc to import all envriomental variables
  error=system(paste0(". $HOME/.bashrc;",exec_code),wait=await)
  ### RETURN ERROR MESSAGE
  if(error!=0){
    mssg=paste0("Execution of environment halted due to internal error ")
    callFUN.callError()
  }

  append_env(from=environment(),to=parent.frame())
}





callFUN.buildCall=function(){
    append_env(to=environment(),from=parent.frame())
    ### Use SGE TASK ID if mode is set to batch otherwise use value
   
    callFUN.setCall()

    if(name_env=="self"){
       exec_code=paste0("Rscript -e \" invisible(",
            ns,"::",fn,"(env=\\\"",
            parent.env$env_file,"\\\"))\""
      )
    }else if (name_env=="parent"){
      if(rmode=="local"){
              ### LOCALLY WE ARE LIMITED IN NUMBER OF CORES
              cores=parallel::detectCores()-1
              ### WE ASSIGN A REASONABLE NUMBER OF JOBS FoR THE REQUESTED NUMBER OF THREADS 
              rjobs=floor(cores/threads)
              exec_code=paste0("Rscript -e \" invisible(parallel::mclapply(1:",n_inputs,
              ",FUN=function(select){",ns,"::",fn,"(env=\\\"",
              env_file,"\\\",select=select)},mc.cores=",rjobs,"))\"")
        }else if(rmode=="batch"){
              cores=Inf
              rjobs=Inf

              ### IN BATCH WE ASSUME INFINITE RESOURCES 
              exec_code=paste0("Rscript -e \" invisible(",
              ns,"::",fn,"(env=\\\"",env_file,
              "\\\",select=$SGE_TASK_ID))\"")
  
              ### WE CREATE SCHEDULER SPECIFIC VARIABLES

              batch_code=paste(
                    "qsub -V ", 
              paste0("-N ",env_id),
              paste0(" -t 1-",n_inputs),
              paste0(" -l h_rt=",time),
              paste0(" -l mem=",ram,"G"),
              paste0(" -pe smp ",threads), 
              paste0(" -wd ",getwd()), 
              paste0(" -o ",work_dir,"/",env_id,".std_out"),
              paste0(" -e ",work_dir,"/",env_id,".std_error"))

              if(exists("hold")){
                batch_code=paste0(batch_code,paste0(" -hold_jid ",paste0(hold,collapse=",")))
              }

              if(bypass){
                batch_code=paste0(batch_code," -P crag7day ")
              }

              ### WE APPEND CONFIG EXEC CODE AND BATCH CODE DATA

              exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,
              ";",exec_code,"'|",batch_code)
        }else{
           mssg="Incorrect run mode selected. Available run modes include: `local` and `batch` "
           callFUN.callError()
        }
    }else if(name_env=="child"){
      ### OF THIS IS A CHILD ENVIROMENT WE DIRECTLY RUN USER DEFINED FUNCTION
      FUN()
    }else{
      mssg=paste0("Unknown environment : " , name_env)
      callFUN.callError()
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
  if(name_env=="self"){
      env_file=paste0(work_dir,"/",env_id,".self.RData")
      saveRDS(environment(),file = env_file)
  }else if (name_env=="parent"){
      env_file=paste0(work_dir,"/",env_id,".parent.RData")
      saveRDS(environment(),file = env_file)
  }else if(name_env=="child"){
      env_file=paste0(work_dir,"/",env_id,".child.RData")
      saveRDS(environment(),file= env_file)
  }else{
    mssg=paste0("Unknown environment : " , name_env)
    callFUN.callError()
  }
  append_env(from=environment(),to=parent.frame())
}


callFUN.readEnv=function(){
  append_env(to=environment(),from=parent.frame())
  if(name_env=="self"){
    parent.env=readRDS(parent.env$env_file)
  }else if(name_env=="parent"){
    child.env=lapply(1:n_inputs,function(n){
        env=readRDS(child.env[[n]]$env_file)}
    )
  }

  append_env(from=environment(),to=parent.frame())
}


make_hash_id=function(name,fn,vars=NULL){
  hash=rlang::hash(c(name,fn,vars))
  hash=paste0(fn,".",hash)
  return(hash)
}






callFUN.buildId=function(){
  append_env(to=environment(),from=parent.frame())
  env_id<-make_hash_id(name=name_env,fn=fn,vars=sheet)
  append_env(from=environment(),to=parent.frame())
}




callFUN.buildError<-function(){
  append_env(to=environment(),from=parent.frame())
  if(exists("err_msg")){
     err_msg <- paste0(err_msg ,name_env," ( ",env_id," ) "," -> ")
  }else{
     err_msg <- paste0("CRITICAL ERROR: ",name_env," ( ", env_id," ) "," -> ")
  }
  append_env(from=environment(),to=parent.frame())
}


callFUN.callError<-function(){
  append_env(to=environment(),from=parent.frame())
  light_red <- crayon::make_style("tomato")
  cat(crayon::red(paste0("[",Sys.time(),"]", "[",name_env,"]","[",env_id,"]\n")))
  cat(light_red(paste0(mssg, "\n\n")))
   opt <- options(show.error.messages=FALSE)
  on.exit(options(opt))
  stop()

}

callFUN.callWarning<-function(){
  append_env(to=environment(),from=parent.frame())
  orange <- crayon::make_style("orange")
  cat(orange(paste0("[",Sys.time(),"]", "[",name_env,"]","[",env_id,"]\n")))
  cat(crayon::yellow(paste0(mssg,"\n\n")))
}



callFUN.buildSheet=function(){
  append_env(to=environment(),from=parent.frame())
  ### CHECK IF VARIABLE SHEET IS GIVEN
  ### IF SHEET IS GIVEN WE ASSIGN THE VARIABLES AND QUIT
  if (exists("sheet")){
    if(file.exists(sheet)){
      ### READ VARIABLE SHEET
      sheet=read.delim(sheet,header=TRUE)
    }else{
      stop(
        paste0("sheet: ", sheet, "doesn't exist")
      )
    }
  }else{
  ### SET CREATE A SHEET FROM THE FUNCTION VARIABLES
      sheet=as.data.frame(as.list(parent.frame())[fn_vars])
  }

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

callFUN.assignSheetChild=function(){
  .this.env=environment()
  .base.env=parent.frame()
  append_env(to=environment(),from=)
  sheet=sheet[row,]
  assign("sheet",sheet)
  for (col in 1:n_inputs){
    assign(names(sheet)[col],sheet[,col],envir=.base.env)
  }
}


callFUN.setProcess=function(){
  append_env(to=environment(),from=parent.frame())

  ## CREATE ERROR MESSAGE FOR EACH CHILD
  callFUN.buildError()

  ## CREATE DIRECTORIES FOR HANDLE PROCESS
  callFUN.buildDir()

  append_env(from=environment(),to=parent.frame())
}



callFUN.buildChild=function(){
  .base.env=parent.frame()
  append_env(to=environment(),from=.base.env)
  
  name_env<-"child"
  await<-TRUE
  ### CREATE CHILD JOB ID
  callFUN.buildId()

  ### CREATE SHEET
  callFUN.assignSheetChild()

  print(sheet)

  ### SET PROCESSS
  callFUN.setProcess()
  
  
  ### CHECK VARIABLE SUBTYPES
  callFUN.checkSubtypes()

  ### WE DUMP CHILDREN INFO
  callFUN.dumpInfo()

  ### WE WRITE CHILDREN ENV
  callFUN.writeEnv()

  return(environment())

}



callFUN.buildParent=function(){
  
    append_env(to=environment(),from=parent.frame())
    ### IF WE INHERIT A RDS FILE READ AND APPEND
    ### WE SAVE RDS FILES WHEN SUBMITTING JOBS TO SCHEDULER OR JOBS RUN LOCALLY IN PARALLEL
    ### SELECT VARIABLE DEFINES WHICH ENV WE ARE RUNNING
    ### WE APPEND THIS ENV AND STOP

    name_env<-"parent"

    await<-TRUE
    
    ### CREATE CHILD JOB ID
    callFUN.buildId()

    ### WE BUILD THE PROCESS AND THE ERROR MESSAGES
    callFUN.setProcess()
    
    ### WE CHECK VARIABLE TYPES
    callFUN.checkTypes()

  
    ### CHECK IF WE REQUIRE REMOTE DATA
    callFUN.remoteValidate()


    ### WE SET THE DUMPSTER WHERE TO PUT CHILDREN INFO
    callFUN.dumpInfo()

    print(sheet)

    ### FROM THE PARENT ENVIRONMENT WE CREATE A NEW ENVIRONMENT FOR EACH UNIQUE INPUT
    child.env=list()
    for(row in 1:n_inputs){
        child.env[[row]]<-callFUN.buildChild()
    }

    ### WE WRITE PARENT ENV
    callFUN.writeEnv()
    
    return(environment())
}


callFUN.setSelf<-function(){
  append_env(to=environment(),from=parent.frame())


  if(!exists("verbose")){
    verbose<-FALSE
  }


  ### CHECK IF RUN MODE VARIABLE EXISTS
  if(!exists("rmode")){
    mssg="Variable: `rmode` has not been provided. Setting default run mode to : `local`"
    callFUN.callWarning()
    ## IF NOT SET TO LOCAL
    rmode<-"local"
  }


  if(!exists("threads")){
     mssg="Variable: `threads` has not been provided. Setting default threads to : `1` "
     callFUN.callWarning()
    threads<-1
  }

  if(rmode=="batch"){
    if(!exists("time")){
       mssg="Variable: `time` has not been provided. Setting default run time to : `48:00:00` "
       callFUN.callWarning()
      time<-"48:0:0"
    }
    
    if(!exists("ram")){
       mssg="Variable: `ram` has not been provided. Setting default ram (Gb) to : `1` "
       callFUN.callWarning()
      ram<-1
    }

    if(!exists("bypass")){
       mssg="Variable: `bypass` has not been provided. Setting default wallclock bypass to : `FALSE` "
       callFUN.callWarning()
    }
  }

   
  if(!exists("await")){
     mssg="Variable: `await` has not been provided. Setting default strategy to await for parent : `TRUE`"
     callFUN.callWarning()
    await<-TRUE
  }


  if(!exists("lic_dir")){
    mssg=text=paste0("Variable: `lic_dir` has not been provided. Setting default license directory to : `",build_default_license_list()$dir)
     callFUN.callWarning()
    lic_dir<-build_default_license_list()$dir
  }


  if(!exists("batch_cfg")){
    mssg=paste0("Variable: `batch_cfg` has not been provided. Setting default config for batch mode to : `",build_default_preprocess_config())
     callFUN.callWarning()
    batch_cfg<-build_default_preprocess_config()
  }

  if(!exists("preserve")){
     mssg="Variable: `preserve` has not been provided. Setting default strategy to deal with working directory  to : `partial` "
     callFUN.callWarning()
    preserve<-"partial"
  }

  if(!exists("compl")){
     mssg="Variable: `compl` has not been provided. Setting default strategy to deal with complementary files : [ `bai` , `tbi` ] "
     callFUN.callWarning()
    compl<-c(".bai",".tbi")
  }

  if(!exists("work_dir")){
        mssg=paste0("Variable: `work_dir` has not beed provided. Setting default working directory to `",getwd(),"`")
        callFUN.callWarning()
        work_dir<-getwd()
      }

  if(!exists("output_dir")){
        mssg=paste0("Variable: `output_dir` has not beed provided. Setting default output directory to `",getwd(),"`")
        callFUN.callWarning()
        output_dir<-getwd()
  }

  if(!exists("overwrite")){
     mssg="Variable: `overwrite` has not been provided. Setting default to overwrite the output_dir: FALSE "
     callFUN.callWarning()
    overwrite=FALSE
  }

  if(!exists("rds")){
      rds=list()
      rds$user=Sys.getenv("RDS_USER")
      rds$password=Sys.getenv("RDS_PASS")
      rds$node=Sys.getenv("RDS_NODE")
      mssg=paste0("Variable: `rds` has not been provided. Setting default RDS config using enviromental variables :\n", 
     "\t`$RDS_USER` : ", rds$user,"\n",
     "\t`$RDS_PASS` : ", rds$password,"\n",
     "\t`$RDS_NODE` : ", rds$node)
      callFUN.callWarning()
  }else{
      if(is.null(rds$user)){
        rds$user<-Sys.getenv("RDS_USER")
        mssg=paste0("Variable: `rds` has been provided but `user` not given. Setting default RDS user from enviromental variable :\n",
        "\t`$RDS_USER` : ",  rds$user)
        callFUN.callWarning()
      }
      if(is.null(rds$password)){
          rds$password<-Sys.getenv("RDS_PASS")
        mssg=paste0("Variable: `rds` has been provided but `password` not given. Setting default RDS user from enviromental variable :\n", 
        "\t`$RDS_PASS` : ", rds$password)
        callFUN.callWarning()
      }

      if(is.null(rds$node)){
        rds$node<-Sys.getenv("RDS_NODE")
         mssg=paste0("Variable: `rds` has been provided but `node` not given. Setting default RDS user from enviromental variable :\n", 
        "\t`$RDS_NODE` : ",rds$node)
        callFUN.callWarning()
      }
  }

 append_env(from=environment(),to=parent.frame())
}


callFUN.buildSelf=function(){
  append_env(to=environment(),from=parent.frame())

  name_env<-"self"

  ## SET NAMESPACE
  ns <- "ULPwgs"
  
  ## WE DEFINE THE MAIN NAME FOR THE RUNNING FUNCTION
  fn <- sub(".*::","",sub("\\(.*","",
    paste0(deparse(sys.calls()[[1]]),collapse=","))
  )

  ### CREATE SHEET WITH VARIABLES
  callFUN.buildSheet()

  ### CREATE CHILD JOB ID
  callFUN.buildId()

  ### WE SET THE DEFAULT VARIABLES
  callFUN.setSelf()

  ### WE BUILD THE PROCESS AND THE ERROR MESSAGES
  callFUN.setProcess()

  ### WRITE SELF ENV
  callFUN.writeEnv()

  ### WE BUILD THE PARENT ENVIRONMENT
  parent.env=callFUN.buildParent()

  append_env(from=environment(),to=parent.frame())
}

callFUN.buildDir=function(){
    append_env(
        to=environment(),
        from=parent.frame()
      )

        work_dir=set_dir(
                dir=work_dir,
                name=env_id
          )

      if (name_env=="child"){
          if(!exists("remote")){
            ### CREATE OUTPUT DIR
          
            out_file_dir <- set_dir(
                    dir=output_dir
            )
          }else if(any(remote %in% "output_dir")){
            callFUN.remoteCreateDir()
          }
              
          tmp_dir <- set_dir(
              dir=work_dir,
              name="tmp"
          )


          ### CREATE TMP DIRECTORY
          ### WE WILL STORE TMP FILES FOR ALL FUNCTIONS HERE 
        
          out_dir <- set_dir(
              dir=tmp_dir,
              name="out"
          )
    
          ### WITHIN TMP DIRECTORY CREATE DIRECTORY TO STORE SYMLINK OF FILES
          ### WE WILL STORE SYMLINK FILES FOR LOCAL FILES
      
          ln_dir <- set_dir(
            dir=tmp_dir,
            name="ln"
          )
    
        

          ### WITHIN TMP DIRECTORY CREATE DIRECTORY TO STORE REMOTE FILES
          ### WE WILL STORE REMOTE DOWNLOAD FILES HERE
        
          rmt_dir <- set_dir(
            dir=tmp_dir,
            name="rmt"
          )
      }


      ### WE APPEND NEW VARIABLES BACK TO THE MAIN FUNCTION
        append_env(from=environment(),to=parent.frame())
}





callFUN.remoteCreateDir=function(){
    append_env(
      to=environment(),
      from=parent.frame()
    )

    check=system(paste0("sshpass -e ssh -Y ",ip_address,
              " \" mkdir -p ",output_dir," 2>/dev/null ;echo 0\""
      ),intern=TRUE
    )

    if(check!=0){
       mssg=paste0("Failed to create remote output directory : `",output_dir, "`" )
       callFUN.callError()
    }

    mssg=paste0("Succesfully created remote output directory : `",output_dir,"`")
    callFUN.callWarning()
    out_file_dir<-output_dir
    
    append_env(from=environment(),to=parent.frame())
}


callFUN.remoteScpDir=function(){
    append_env(
      to=environment(),
      from=parent.frame()
    )
    check=system(paste0("sshpass -e ",password,
            " scp -r ",paste0(work_dir,"/*"),ip_address,":",out_file_dir, "; echo $?" 
      ),intern=TRUE
    )

    if(check!=0){
        mssg=paste0("Failed to copy data to remote output directory : [ ",out_file_dir, " ] " )
        callFUN.callError()
    }
   

    append_env(from=environment(),to=parent.frame())
}


callFUN.remoteValidate=function(){
    append_env(
      to=environment(),
      from=parent.frame()
    )
    if(exists("remote")){

      node=ifelse(rds$node=="",NULL,rds$node)
      user=ifelse(rds$user=="",NULL,rds$user)
      password=ifelse(rds$password=="",NULL,rds$password)
      ip_address=paste0(user,ifelse(!is.null(user),"@",""),node)
    
      Sys.setenv("SSHPASS"=password)
      ### WE PING THE REMOTE SERVER TO SEE IF AVAILABLE

      check=system(paste0("sshpass -e ",
              " ssh -Y ", ip_address," \"test;echo 0\""
        ),intern=TRUE
      )

      callFUN.checkRemote()
    }
    
    append_env(from=environment(),to=parent.frame())
}




callFUN.dumpInfo<-function(){
  append_env(to=environment(),from=parent.frame())

  if(!exists("dump_file")){
    ### TRACING CHILDREN CAN BE DIFFICULT
    ### WE DEFINE THE VARIABLES THAT WILL HELP IDENTIFY CHILDREN IN JOB HIERARCHY

    dump_names=c("id","order",fn_vars)
    
    ## WE DEFINE A FILE WERE WE WILL DUMP THIS INFO
    
    dump_file=paste0(work_dir,"/",env_id,".dump")
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
      id=env_id,
      order=row,
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


callFUN.setOutput=function(...){
      ### WE CREATE ALL OUTPUT VARIABLES WITH THE OUT_FILE_DIR
      append_env(to=environment(),from=parent.frame())
      for(var in names(list(...))){
        out_files[[var]]=paste0(out_file_dir,"/",list(...)[[var]])
      }
      append_env(from=environment(),to=parent.frame())
  }


callFUN.checkRemote=function(){
      append_env(
        to=environment(),
        from=parent.frame()
      )
    ### IF SERVER DOESN'T RESPOND WE RETURN ERROR 
    if(check==255){
      mssg=paste0("Remote server not accessible : `",ip_address,"`")
      callFUN.callError()
    }else if(check==5){
      mssg=paste0("Failed to log in remote server :`",ip_address,"`")
      callFUN.callError()
    }else if(check!=0){
      mssg=paste0("Unknown error with remote server :`",ip_address,"`")
      callFUN.callError()
    }
    mssg=paste0("Remote server available : `", ip_address,"`")
    callFUN.callWarning()

    append_env(from=environment(),to=parent.frame())
  }