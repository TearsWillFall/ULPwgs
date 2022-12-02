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


build_rdata_object=function(
  envir=NULL,
  job_id=NULL,
  output_dir="."
){
          out_file_dir=set_dir(dir=output_dir,name="RData")
          rdata_file=paste0(out_file_dir,"/",job_id,".RData")
          save(list = ls(envir,all.names = TRUE),envir=envir,file = rdata_file)
          return(rdata_file)
  }
      




#' Build execution innit for function to execute
#' 
#' @param objects List of objects to save in RData object
#' @param job_id Job Identifier
#' @param nspace Function namespace. Default ULPwgs
#' @param fun Function name. Default NULL
#' @param inherit_scheduler Inherit submission order from scheduler. Default FALSE
#' @param output_dir Path to output directory 
#' @export


build_exec_innit=function(
        envir=NULL,
        job_id=NULL,
        output_dir=".",
        nspace="ULPwgs",
        fun=NULL,
        inherit_scheduler=FALSE,
        selected=NULL
      ){

        if(is.null(fun)){
            stop("fun argument can't be type NULL.")
        }

        if(is.null(selected)){
          stop("selected argument can't be of type NULL.")
        }
        ### Create RData object to inherit vars

        rdata=build_rdata_object(
          envir=envir,
          job_id=job_id,
          output_dir=output_dir
        )
      
        ### Use SGE TASK ID if mode is set to batch otherwise use value

        selected=ifelse(inherit_scheduler,"$SGE_TASK_ID",selected)
         
        exec_code=paste0("Rscript -e \" ",nspace,"::",fun,"(rdata=\\\"",
        rdata,"\\\",selected=$SGE_TASK_ID)\"")
        
        return(exec_code)

}




#' Build execution innit for batch
#' 
#' @param objects List of objects to save in RData object
#' @param job_id Job Identifier
#' @param nspace Function namespace. Default ULPwgs
#' @param fun Function name. Default NULL
#' @param inherit_scheduler Inherit submission order from scheduler. Default FALSE
#' @param output_dir Path to output directory 
#' @export


build_batch_exec_innit=function(
  exec_code=NULL,
  job_id=NULL,
  time="48:0:0",
  ram=3,
  threads=1,
  wd=getwd(),
  hold=NULL,
  array=NULL,
  output_dir=".",
  batch_config=build_default_preprocess_config()
){  

    if(is.null(exec_code)){
      stop("exec_code argument can't be of type NULL")
    }


    batch_code=build_job_exec(
      job=job_id,time=time,ram=ram,
      threads=threads,wd=wd,
      output_dir=output_dir,
      hold=hold,array=array
    )

    out_file_dir=set_dir(dir=output_dir,name="batch")
    exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,
    ";",exec_code,"'|",batch_code)

    return(exec_code)
}


#' Run job from batch constructor
#' 
#' @param envir Current enviroment.
#' @param job_id Job Identifier
#' @param time Current job time
#' @param time Job time
#' @param output_dir Path to output directory
#' @param nspace Function namespace. Default ULPwgs
#' @param fun Function name. Default NULL
#' @param hold Job identifier to hold on. Default NULL
#' @param output_dir Path to output directory 
#' 
#' 
#' @export


run_job=function(
  envir=environment(),
  slist=NULL,
  job_id=NULL,
  output_dir=".",
  verbose=FALSE,
  nspace="ULPwgs",
  fun=NULL,
  hold=NULL,
  mode="local",
  time="48:00:00",
  batch_config=build_default_preprocess_config(),
  threads=1,
  ram=1,
  error_mssg="Job failed"
){

 
    if(mode=="local"){
        lapply(seq(1,length(slist)),
            FUN=function(selected){
                exec_code=build_exec_innit(
                      envir=envir,
                      job_id=job_id,
                      output_dir=output_dir,
                      nspace=nspace,
                      fun=fun,
                      inherit_scheduler=FALSE,
                      selected=selected
                )

                if(verbose){
                    print_verbose(
                      job=job_id,
                      arg=as.list(envir),
                      exec_code=exec_code)
                }


                error=execute_job(exec_code=exec_code)
                
                if(error!=0){
                    stop(error_mssg)
                }
            }
        )
        
    }else if(mode=="batch"){
        
        
      n_jobs=length(slist)


      exec_code=build_exec_innit(
            envir=envir,
            job_id=job_id,
            output_dir=output_dir,
            nspace=nspace,
            fun=fun,
            inherit_scheduler=TRUE
      )

      exec_code=build_batch_exec_innit(
            exec_code=exec_code,
            job_id=job_id,
            time=time,
            threads=threads,
            ram=ram,
            hold=hold,
            array=n_jobs,
            output_dir=output_dir,
            batch_config=batch_config
          )
      }

      if(verbose){
          print_verbose(job=job_id,
            arg=as.list(envir),
            exec_code=exec_code
          )
      }


      error=execute_job(exec_code=exec_code)

      if(error!=0){
          stop(error_message)
      }

}



