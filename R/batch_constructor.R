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
    if(!exists("select")){
      callFUN.runSelf()
  
    }else{
      callFUN.runProcess()
    }

    return(.this.env)
}


callFUN.setEnv<-function(
    fn_id=NULL,
    fn_vars=NULL,
    output_dir=".",
    parent_dir=NULL,
    output_name=NULL,
    verbose=FALSE,
    bgzip_idx=FALSE,
    vcf_idx_fmt="tbi",
    lic_dir=build_default_license_list()$dir,
    batch_cfg=build_default_preprocess_config(),
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
    remotes=NULL,
    compl=FALSE,
    parent_id=NULL,
    wait=FALSE,
    hold=NULL
){

    append_env(to=environment(),from=parent.frame())
    
    ## WE CHECK IF WE ARE INHERITING A PARENT ENVIROMENT
    ## ELSE WE CREATE A PARENT ENVIROMENT
    ## WHEN INHERITING WE DON'T CHECK ARG TYPES AGAIN

    if(exists("env")){
      if(!is.environment(env)){
        env <-readRDS(file=env)
        if(exists("select")){
          env <- env$child.envs[[select]]
        }
      }
      append_env(to=environment(),from=env)
      env<-NULL
    }else{
      self<-FALSE

      args=appendList(args,build_default_variable_list())


      ### CREATE VARIABLES FOR THE ENVIRONMENT
      callFUN.buildSelf()

      ### WE BUILD THE PARENT ENVIRONMENT
      callFUN.buildParent()




     
      ### WE WRITE PARENT ENVIROMENT
      callFUN.writeEnv()

      self<-TRUE
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
    check=suppressWarnings(system(paste(
    "sshpass -f ",password,
    " ssh ",ip_remote,
      "\" realpath -e ",arg_value,"\" "),intern=TRUE
    ))
    append_env(from=environment(),to=parent.frame())
}


callFUN.remoteGet=function(){
  append_env(to=environment(),from=parent.frame())
  arg_dir=set_dir(dir=rmt_dir,name=arg)
  ### COPY REMOTE FILE TO LOCAL TMP DIR IF REMOTE FILE EXISTS
  system(paste(
  "sshpass -f ",password,
  " ssh ",ip_remote,
    "\" cp -rn ",check," -t ",
    arg_dir, "\"")
  )
  arg_value=paste0(arg_dir,"/",basename(check))   
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
          stop(paste0("Variable : ",arg,
          " ( type : ",arg_type,
          " ) -> Value: Not defined. Define a value"))
        }

        if(is.null(arg_value)&arg_required){
          stop(paste0("Variable : ",arg,
          " ( type : ",arg_type,
          " ) -> Value: NULL (type: NULL). Define a non-NULL value"))
        }

        if(typeof(arg_value)!=arg_type&arg_required){
          stop(paste0("Variable : ",arg,
          " ( type : ",arg_type," ) -> Value: ",arg_value,
          " ( type : ",typeof(arg_value),
          " ). Invalid type"))
        }
      }
    }
}
     
  
          #       #   ###CHECK IF REMOTE LOCATION HAS BEEN DEFINED
          #       #   if(is.null(node)){
           
          #       #     )
          #       #     }else{
          #       #         ## WE CHECK IF PATH IS IN REMOTE
          #       #         callFUN.remoteCheck()

          #       #         if(length(check)==0){
          #       #           stop(paste0("Variable : ",arg,
          #       #             " ( type : ",arg_type," ) [ subtype : path ] -> Value: ",arg_value,
          #       #             " ( type : ",typeof(arg_value),
          #       #             " ) [ subtype : NULL ] . Path doesn't locally and remotely")
          #       #           )
                        
          #       #         }
          #       #         ### WE GET THE REMOTE PATH
          #       #         callFUN.remoteGet()
          #       #     }
          #       # }
          
          #   callFUN.createLink()

          #   ## WE REASSIGN THE VALUE OF ARGUMENT TO THE LOCAL/REMOTE
          #   .base.env[[arg]]=normalizePath(paste0(arg_dir,"/",basename(arg_value)))
          #   }
          # }


callFUN.checkSubtypes=function(){
   .this.env=environment()
   .base.env=parent.frame()
    append_env(to=.this.env,from=.base.env)

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
                if(!is.null(remotes)){
                  if(!any(remotes %in% arg)){
                    if(!file.exists(arg_value)){
                              stop(paste0(err_msg," -> Variable : ",arg,
                            " ( type : ",arg_type," ) [ subtype : path ] -> Value: ",arg_value,
                            " ( type : ",typeof(arg_value),
                            " ) [ subtype : NULL ] . Path doesn't exist locally"))
                    }
                  }else{
                    callFUN.remoteCheck()
                    if(length(check)==0){
                        stop(paste0(err_msg," -> Variable : ",arg,
                              " ( type : ",arg_type," ) [ subtype : path ] -> Value: ",arg_value,
                              " ( type : ",typeof(arg_value),
                              " ) [ subtype : NULL ] . Path doesn't exist remotely"))
                    }
                  }
                }
            }
          }
        }
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
  
  ### CREATE CALLER
  callFUN.buildCall()

  ### WE RUN THE ENVIROMENT
  callFUN.runCall()
  
  ### WE WRITE RESULTS FOR CHILD TO RDS
  callFUN.writeEnv()

  ### WE MOVE DATA TO RESULTS DIRECTORY
  callFUN.moveData()


  append_env(from=environment(),to=parent.frame())

}


callFUN.moveData=function(){
  append_env(to=environment(),from=parent.frame())
  system(paste0("mv ",child_dir,out_file_dir))
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
      stop(paste0(err_msg," -> Execution halted due to internal error"))
  }

  append_env(from=environment(),to=parent.frame())
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
            env_file,"\\\")\""
          )
      }else{
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
          stop(err_msg, " -> Unknown process mode ")
        }

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

  if(!exists("self_id")){
    self_id=make_unique_id(fn)
  }else if(!exists("job_id")&exists("self_id")){
    
    ### IF parent ID IS NOT GIVEN WE CREATE AN UNIQUE NAME USING THE FUNCTION ID
    parent_id <- make_unique_id(fn)
  
    ### WE ASSIGN A JOB ID TO THE CALLER FUNCTION

    job_id <- build_job(
      parent_id=self_id,
      child_id=parent_id
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
  if(exists("err_msg")){
     err_msg <- paste0(err_msg ,fn," (",job_id,") "," -> ")
  }else{
     err_msg <- paste0("CRITICAL ERROR: ",fn," (",self_id,") "," -> ")
  }
  append_env(from=environment(),to=parent.frame())
}






callFUN.buildSheet=function(){
  append_env(to=environment(),from=parent.frame())
  ### CHECK IF VARIABLE SHEET IS GIVEN
  ### IF SHEET IS GIVEN WE ASSIGN THE VARIABLES AND QUIT
  if (exists("sheet")){
    if(file.exists(sheet)){
      ### READ VARIABLE SHEET
      sheet=read.delim(sheet,header=TRUE)
    }
  }else{
  ### SET CREATE A SHEET FROM THE FUNCTION VARIABLES
      sheet=as.data.frame(as.list(parent.frame())[fn_vars])
  }
  append_env(from=environment(),to=parent.frame())

}

callFUN.assignSheetParent=function(){
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



callFUN.assignSheetChild=function(){
  append_env(to=environment(),from=parent.frame())
  for (col in 1:n_inputs){
    assign(names(sheet)[col],sheet[row,col])
  }  
  append_env(from=environment(),to=parent.frame())
}


callFUN.setProcess=function(){
  append_env(to=environment(),from=parent.frame())
  ### CREATE CHILD JOB ID
  callFUN.buildId()

  ## CREATE ERROR MESSAGE FOR EACH CHILD
  callFUN.buildError()

  ## CREATE DIRECTORIES FOR HANDLE PROCESS
  callFUN.buildDir()

  append_env(from=environment(),to=parent.frame())
}


callFUN.buildChilds<-function(){

      append_env(to=environment(),from=parent.frame())

      ### WE SET THE DUMPSTER WHERE TO PUT CHILDREN INFO
      callFUN.dumpInfo()

      ### FROM THE PARENT ENVIRONMENT WE CREATE A NEW ENVIRONMENT FOR EACH UNIQUE INPUT

      child.envs=list()
      for(row in 1:n_inputs){
        get_child_env=function(row){
           append_env(to=environment(),from=parent.frame())
           ### USING THE VAR SHEET WE REASSIGN THE VALUES OF EACH VARIABLE
           callFUN.assignSheetChild()
           callFUN.buildChild()
           return(environment())
        }
        child.envs[[row]]<-get_child_env(row)
       }

      append_env(from=environment(),to=parent.frame())

}

callFUN.buildChild=function(){
  append_env(to=environment(),from=parent.frame())

  ### SET PROCESSS
  callFUN.setProcess()
  
  ### CHECK VARIABLE SUBTYPES
  callFUN.checkSubtypes()

  ### WE DUMP CHILDREN INFO
  callFUN.dumpInfo()

  ### WE WRITE CHILDREN ENV
  callFUN.writeEnv()

  append_env(from=environment(),to=parent.frame())

}



callFUN.buildParent=function(){
  
    append_env(to=environment(),from=parent.frame())
    ### IF WE INHERIT A RDS FILE READ AND APPEND
    ### WE SAVE RDS FILES WHEN SUBMITTING JOBS TO SCHEDULER OR JOBS RUN LOCALLY IN PARALLEL
    ### SELECT VARIABLE DEFINES WHICH ENV WE ARE RUNNING
    ### WE APPEND THIS ENV AND STOP

    ### WE BUILD THE PROCESS AND THE ERROR MESSAGES
    callFUN.buildProcess()
    

    ### WE IMPORT AND/OR CREATE SHEET TO APOINT VARIABLES
    callFUN.assignSheetParent()

    ### WE CHECK VARIABLE TYPES
    callFUN.checkTypes()

    ### CHECK IF WE REQUIRE REMOTE DATA
    callFUN.remoteValidate()

 


    ### CREATE WORK DIRECTORIES
    callFUN.buildDir()

    append_env(from=environment(),to=parent.frame())
}


callFUN.buildSelf=function(){
  append_env(to=environment(),from=parent.frame())
   ## IF NOT SET UP BY THE USER WE WILL GET THE MAIN FUNCTION NAME
    
  if(!exists("fn")){
    ## GET CALLER FUNCTION NAME
    fn <- sub(".*::","",sub("\\(.*","",
      paste0(deparse(sys.calls()[[sys.nframe()-3]]),collapse=","))
    )
  }

  #### IF FUNCTION ID IS NOT GIVE WE USE FN NAME AS ID
  #### OTHERWISE WE APPEND FUNCTION ID

  if(is.null(fn_id)){
    fn_id<-fn
  }else{
    fn_id<-paste0(fn,".",fn_id)
  }

  ### WE BUILD THE PROCESS AND THE ERROR MESSAGES
  callFUN.buildProcess()


  ### CREATE SHEET WITH VARIABLES
  callFUN.buildSheet()
  

  append_env(from=environment(),to=parent.frame())
}

callFUN.buildDir=function(){
    append_env(
        to=environment(),
        from=parent.frame()
      )

      if(!exists("child_id")&!exists("parent_id")){
        
        self_dir <- set_dir(
                dir=work_dir,
                name=self_id
        )
        
      }else if(!exists("child_id")&exists("parent_id")){

        ### CREATE MAIN WORKING DIRECTORY
        
        parent_dir <- set_dir(
                dir=self_dir,
                name=parent_id
        )

      
        ### CREATE TMP DIRECTORY
        ### WE WILL STORE TMP FILES FOR ALL FUNCTIONS HERE 
      
        tmp_dir <- set_dir(
            dir=parent_dir,
            name="tmp"
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

        ### CREATE DIRECTORY TO STORE ENVIROMENT DATA
        ### WE WILL STORE R ENVIRONMENT DATA HERE
        
        env_dir<- set_dir(
          dir=parent_dir,
          name="env"
        )


        ### CREATE DIRECTORY TO BATCH DATA
        ### WE WILL STORE SCHEDULER DATA HERE
      
        batch_dir<- set_dir(
          dir=parent_dir,
          name="batch"
        )


      }else if (exists("child_id")){

        if(!remote_output){
           ### CREATE OUTPUT DIR
        
          out_file_dir <- set_dir(
                  dir=output_dir
          )
        }else{
  
          callFUN.remoteCreateDir()
        }
       

        child_dir <- set_dir(
            dir=tmp_dir,
            name=child_id
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

    check=system(paste0("sshpass -f ",password,
            " ssh ",paste0(user,
              ifelse(!is.null(user),"@",""),node),
              "\" mkdir ",out_file_dir,"; echo $? \""
      )
    )

    if(check!=0){
       stop(paste0(err_msg," -> Failed to create remote output directory : [ ",out_file_dir, " ] " ))
    }
   
    append_env(from=environment(),to=parent.frame())
}


callFUN.remoteScpDir=function(){
    append_env(
      to=environment(),
      from=parent.frame()
    )
    check=system(paste0("sshpass -f ",password,
            " scp -r ",paste0(child_dir,"/*"),ip_remote,":",out_file_dir, "; echo $?" 
      ),intern=TRUE
    )

    if(check!=0){
       stop(paste0(err_msg," -> Failed to move data to remote output directory : [ ",out_file_dir, " ] " ))
    }
   

    append_env(from=environment(),to=parent.frame())
}



callFUN.remoteValidate=function(){

    append_env(
      to=environment(),
      from=parent.frame()
    )

    ip_remote=paste0(user,paste0(ifelse(!is.null(user),"@",""),node))

    ### WE PING THE REMOTE SERVER TO SEE IF AVAILABLE

    check=system(paste0("sshpass -f ",password,
            " ssh -q ", ip_remote," exit"
      ),intern=TRUE
    )

    ### IF SERVER DOESN'T RESPOND WE RETURN ERROR 
    if(check==255){
      stop(paste0(err_msg," -> "," [ ",node," ] " ,"  Server  not accessible. Wrong IP or server currently unavailable"))
    }else if(
      check==5
    ){
      stop(paste0(err_msg," -> "," [ ",node," ] " ," -> Failed to log into remote server. Incorrect login details"))
    }

    ## WE CHECK IF THE REMOTE VARIABLES PROVIDED ARE INPUTS OR/AND OUTPUT DIR
    remotes_input=remotes[!remotes %in% "output_dir"]
    remotes_output=remotes[remotes %in% "output_dir"]
    
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












