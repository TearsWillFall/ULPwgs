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


storeEnv=function(){
  UseMethod("storeEnv")
}

storeEnv.parent.write=function(){
  append_env(to=environment(),from=parent.frame())
  parent_file=paste0(env_dir,"/",job_id,".parent.RData")
  saveRDS(environment(),file = parent_file)
  append_env(from=environment(),to=parent.frame())
}


storeEnv.child.write=function(){
  append_env(to=environment(),from=parent.frame())
  child_file=paste0(env_dir,"/",job_id,".child.RData")
  saveRDS(environment(),file=child_file)
  append_env(from=environment(),to=parent.frame())
}

storeEnv.child.read=function(){
  append_env(to=environment(),from=parent.frame())
  ### Reads mains and updates values for enviroments with data
    child.envs=lapply(1:n_inputs,function(n){
        child.env=readRDS(child.envs[[n]]$child_file)}
  )
  append_env(from=environment(),to=parent.frame())
}


storeEnv.parent.read=function(){
  append_env(from=readRDS(parent_file),to=parent.frame())
}


 append_env(to=environment(),from=parent.frame())

  

  append_env(from=environment(),to=parent.frame())



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





buildCall=function(){
  UseMethod("buildCaller")
}


buildCall.init=function(){
    append_env(to=environment(),from=parent.frame())
    ### Use SGE TASK ID if mode is set to batch otherwise use value
    
    if(mode=="local"){
          exec_code=paste0("Rscript -e \" invisible(lapply(1:",n_inputs,
          ",FUN=function(select){",ns,"::",fn,"(.env=\\\"",
          parent_file,"\\\",select=select)}))\"")
    }else if(mode=="local_parallel"){
          exec_code=paste0("Rscript -e \" invisible(parallel::mclapply(1:",n_inputs,
          ",FUN=function(select){",ns,"::",fn,"(.env=\\\"",
          parent_file,"\\\",select=select)},mc.cores=",threads-1,"))\"")
    }else if(mode=="batch"){
          exec_code=paste0("Rscript -e \" invisible(",
          ns,"::",fn,"(.env=\\\"",parent_file,
          "\\\",select=$SGE_TASK_ID))\"")
          buildCaller.batch()

    }else{
      stop("Unkown mode type")
    }

    append_env(from=environment(),to=parent.frame())
}

buildCall.batch=function(){
    append_env(to=environment(),from=parent.frame())
    exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,
    ";",exec_code,"'|",batch_code)
    append_env(from=environment(),to=parent.frame())
}

buildCall.batch.ini=function(){
    append_env(to=environment(),from=parent.frame())
    buildCall.batch.init()
    exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,
    ";",exec_code,"'|",batch_code)
    append_env(from=environment(),to=parent.frame())
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

buildCall.batch.init=function(){
  append_env(to=environment(),from=parent.frame())
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
  
}


setVars=function(){
  UseMethod("setVars")
}

setVars.parent=function(){
    append_env(to=environment(),from=parent.frame())
    error=0
    parent_file=""
    main_code=""
    exec_code=""
    batch_code=""
    append_env(from=environment(),to=parent.frame())
}

setVars.child=function(){
    append_env(to=environment(),from=parent.frame())
    exec_code=""
    error=0
    out_file=""
    steps=list()
    out_files=list()
    append_env(from=environment(),to=parent.frame())
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




runEnv=function(){
  UseMethod("runEnv")
}


runEnv.parent=function(){
  append_env(to=environment(),from=parent.frame())
 
  ### WE SET VARS TO COLLECT INFO FROM PARENT
  setVars.parent()

  ### CREATE RDS OBJECT TO STORE ENVIRONMENT
  storeEnv.parent.write()

  ### CREATE CALLER
  buildCall.init()

  runEnv.run()
 
  ### WAIT FOR SCHEDULER TO FINISH
  if(mode=="batch"){
    if(wait){
        wait_scheduler()
    }
  }
  storeEnv.child.read()
  append_env(from=environment(),to=parent.frame())
}


runEnv.child=function(){

  append_env(to=environment(),from=parent.frame())
  
  ### WE SET VARS TO COLLECT INFO FROM CHILD
  setVars.child()

  runEnv.consolidate()
  FUN()
  storeEnv.child.read()

  if(clean){
    build_clean_exec()
  }

  runEnv.run()

  append_env(from=environment(),to=parent.frame())

}

runEnv.run=function(){
  append_env(to=environment(),from=parent.frame())
 
   ### PRODUCE VERBOSE
  if(verbose){
        print_verbose(job=job_id,
          arg=as.list(environment())[dump_names],
          exec_code=exec_code
        )
  }

  #Read .bashrc to import all envriomental variables
  error=system(paste0(". $HOME/.bashrc;",exec_code))
  ### RETURN ERROR MESSAGEchild_order
  if(error!=0){
      stop(err_msg)
  }

  append_env(from=environment(),to=parent.frame())
}


runEnv.consolidate=function(){
    .base.env=environment()
    append_env(to=environment(),from=parent.frame())
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
          system(paste("ln -fs ",var_value, var_dir," > /dev/null 2>&1 "))
          ### UPDATE THE VARIABLE TO THE SYMLINK
          .base.env[[var]]=paste0(var_dir,"/",basename(var_value))

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
                var_value=paste0(var_dir,"/",basename(check))     
                var_dir=set_dir(dir=ln_dir,name=var)

                ### UPDATE THE VARIABLE TO THE SYMLINK
                system(paste("ln -fs ",var_value, var_dir," > /dev/null 2>&1 "))
                .base.env[[var]]=paste0(var_dir,"/",basename(var_value))
                
              }
            }
          }
        }
      }
}


















buildDirs=function(){
    UseMethod("buildDirs")
}

buildDirs.parent=function(){
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
              dir=parent_dir,
              name="env"
          )
        }

      ### CREATE DIRECTORY TO BATCH DATA
      ### WE WILL STORE SCHEDULER DATA HERE
      if(is.null(batch_dir)){
          batch_dir<- set_dir(
            dir=parent_dir,
            name="batch"
        )
      }

    ### WE APPEND NEW VARIABLES BACK TO THE MAIN FUNCTION
    append_env(from=environment(),to=parent.frame())
}


buildDirs.child=function(){
      append_env(
        to=environment(),
        from=parent.frame()
      )

      ### WITHIN TMP DIRECTORY CREATE DIRECTORY TO STORE SYMLINK OF FILES
      ### WE WILL STORE SYMLINK FILES FOR LOCAL FILES
    
      ln_dir <- set_dir(
        dir=ln_dir,
        name=child_id
      )


      ### WITHIN TMP DIRECTORY CREATE DIRECTORY TO STORE REMOTE FILES
      ### WE WILL STORE REMOTE DOWNLOAD FILES HERE
      
      rmt_dir <- set_dir(
          dir=rmt_dir,
          name=child_id
      )
      

      ### CREATE DIRECTORY TO STORE ENVIROMENT DATA
      ### WE WILL STORE R ENVIRONMENT DATA HERE
    
      env_dir<- set_dir(
          dir=env_dir,
          name=child_id
      )
        

      ### CREATE DIRECTORY TO BATCH DATA
      ### WE WILL STORE SCHEDULER DATA HERE
   
      batch_dir<- set_dir(
        dir=batch_dir,
        name=child_id
      )


    ### WE APPEND NEW VARIABLES BACK TO THE MAIN FUNCTION
    append_env(from=environment(),to=parent.frame())
}




buildId=function(){
      UseMethods("buildId")
}

buildId.child=function(){
  append_env(to=environment(),from=parent.frame())
  
  ### WE CREATE A JOB ID FOR EACH CHILD
  ### CHILDREN SHALL WORK!!

  child_id <- make_unique_id(fn)
  
  job_id <- build_job(
    parent_id=parent_id,
    child_id=child_id
  )

   append_env(from=environment(),to=parent.frame())
}

buildId.parent=function(){
  append_env(to=environment(),from=parent.frame())
  
  ### IF parent ID IS NOT GIVEN WE CREATE AN UNIQUE NAME USING THE FUNCTION ID
    if(is.null(parent_id)){
      parent_id <- make_unique_id(fn)
    }

    ### WE ASSIGN A JOB ID TO THE CALLER FUNCTION

    job_id <- build_job(
      parent_id=parent_id
    )

   append_env(from=environment(),to=parent.frame())
}



buildEnv<-function(){
          UseMethods("buildEnv")
    }

buildEnv.reassign<-function(){
  append_env(to=environment(),from=parent.frame())

  ### REASSIGN VAR VALUES IN SHEET TO ENVIRONMENT
  for (col in 1:n_vars){
    assign(names(sheet)[col],sheet[row,col])
  }
  append_env(from=environment(),to=parent.frame())

}


buildEnv.child=function(){
  append_env(to=environment(),from=parent.frame())

  ### USING THE VAR SHEET WE REASSIGN THE VALUES OF EACH VARIABLE
  buildEnv.reassign()

  ### CREATE CHILD JOB ID
  buildId.child()

  ### CREATE WORKING DIRECTORIES FOR EACH CHILD

  buildDirs.child()

  ## CREATE ERROR MESSAGE FOR EACH CHILD

  buildErrorMessage.child()

  ### CREATE RDS FORMAT
  storeEnv.child.write()

  append_env(from=environment(),to=parent.frame())

}

buildEnv.parent=function(
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
  node=NULL,
  user=NULL,
  password=NULL,
  sheet=NULL,
  env=NULL,
  select=NULL,
  err_msg=NULL,
  remote=FALSE,
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
    buildEnv.parent.set()
      
    ### CREATE JOB ID FOR THE PARENT FUNCTION
    buildId.parent()

    ### CREATE WORK DIRECTORIES
    buildDirs.parent()
    
    ### CREATE SHEET WITH VARIABLES
    buildSheet.create()
    
    ### CREATE AN ERROR MESSAGE TO TRACK WHERE JOB FAILS
    buildErrorMessage.parent()

    append_env(from=environment(),to=parent.frame())
}



buildEnv.parent.set=function(){
  append_env(to=environment(),from=parent.frame())
   ## IF NOT SET UP BY THE USER WE WILL GET THE MAIN FUNCTION NAME
    
  if(is.null(fn)){
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


  ### CHECK IF FILES PATH CAN BE LOCATED REMOTELY
  if(!is.null(node)){
    remote=TRUE
  }

  append_env(from=environment(),to=parent.frame())
}




dumpInfo<-function(){
  UseMethod("dumpInfo")
}


dumpInfo.set<-function(){
  append_env(to=environment(),from=parent.frame())

  ### TRACING CHILDREN CAN BE DIFFICULT
  ### WE DEFINE THE VARIABLES THAT WILL HELP IDENTIFY CHILDREN IN JOB HIERARCHY

  dump_names=c("parent_id","child_order","child_id",fn_vars)
  
  ## WE DEFINE A FILE WERE WE WILL DUMP THIS INFO
  
  dump_file=paste0(parent_dir,"/",parent_id,".dump")
  ### WE CREATE AN ENVIRONMENT FOR EACH INPUT VALUE

  write.table(
        file=dump_file,
        x=data.frame(dump_names),
        sep="\t",
        col.names=FALSE,
        row.names=FALSE,
        quote=FALSE
  )

  append_env(from=environment(),to=parent.frame())
}


dumpInfo.append=function(){
  append_env(to=environment(),from=parent.frame())
  
  
  ### WE APPEND THE CHILDREN INFO TO THE DUMP FILE

  info=data.frame(
      parent_id=parent_id,
      child_order=row,
      child_id=child_id,
      sheet[row,]
  )

    write.table(
      file=dump_file,
      x=info,
      sep="\t",
      append=TRUE,
      col.names=FALSE,
      row.names=FALSE,
      quote=FALSE
  )          
}


buildSheet=function(){
  UseMethod("buildSheet")
}

buildSheet.create=function(){
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

  buildSheet.ascertain()

  append_env(from=environment(),to=parent.frame())

}

buildSheet.ascertain=function(){
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










buildErrorMessage<-function(){
    UseMethod("buildErrorMessage")
}

buildErrorMessage.child<-function(){
  append_env(to=environment(),from=parent.frame())
  err_msg <- paste0(err_msg ,fn," (",job_id,") "," -> ")
  append_env(from=environment(),to=parent.frame())
}

buildErrorMessage.parent<-function(){
  append_env(to=environment(),from=parent.frame())
  if(is.null(err_msg)){
    err_msg <- paste0("CRITICAL ERROR: ",fn," (",job_id,") "," -> ")
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
      .this.env=environment()
      .base.env=parent.frame()
      ### ADD OTHER VARIABLES TO BASE ENV
      list2env(x=list(...),envir=.base.env)
      
      ## GET VARIABLE NAMES
      fn_vars=names(.base.env)[!grepl("\\.|FUN",names(.base.env))]

      append_env(to=environment(),from=parent.frame())

      ## WE CHECK IF WE ARE INHERITING A PARENT ENVIROMENT
      ## ELSE WE CREATE A PARENT ENVIROMENT

      if(exists(".env")){
        if(!is.environment(.env)){
          .env <-readRDS(file=.env)
        }else{
          append_env(to=environment(),from=.env)
        }
        append_env(to=environment(),from=.env$child.envs[[select]])
      }else{
        ### WE BUILD THE PARENT ENVIRONMENT
        
        buildEnv.parent()

        ### WE SET THE DUMPSTER WHERE TO PUT CHILDREN INFO
        dumpInfo.set()

        ### FROM THE PARENT ENVIRONMENT WE CREATE A NEW ENVIRONMENT FOR EACH UNIQUE INPUT
        child.envs=parallel::mclapply(
            1:n_inputs,
            function(row,.env){

              append_env(to=environment(),from=.env)


              ### WE BUILD THE CHILD ENV
              buildEnv.child()

              ### WE DUMP CHILDREN INFO
              dumpInfo.append()

              return(environment())
              },
              .env=environment(),
              mc.cores=parallel::detectCores()
        )
      }
      ## WE WILL LAUNCH THE PARENT OR THE CHILD FUNCTION
      if(is.null(select)){
        runEnv.parent()
      }else{
        runEnv.child()
      }

      append_env(from=environment(),to=parent.frame())
      

    }



