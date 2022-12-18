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
  envir=NULL
){
          envir$rdata_file=paste0(envir$out_file_dir_tmp,"/",envir$job_id,".RData")
          saveRDS(envir,file = envir$rdata_file)
  }
      




#' Build execution innit for function to execute
#' 
#' @param envir List of objects to save in RData object
#' @export


build_exec_innit=function(
        envir=NULL
      ){

        ### Create RData object to inherit vars

        build_rdata_object(envir=envir)
        
        ### Use SGE TASK ID if mode is set to batch otherwise use value
        
        if(envir$mode=="local"){
             
              envir$exec_code=paste0("Rscript -e \" lapply(1:",envir$n_inputs,
              ",FUN=function(select){",envir$ns,"::",envir$fn,"(inherit=\\\"",
              envir$rdata_file,"\\\",select=select)})\"")
        
        }else if(envir$mode=="batch"){
    
              envir$exec_code=paste0("Rscript -e \" ",
              envir$ns,"::",envir$fn,"(inherit=\\\"",envir$rdata_file,
              "\\\",select=$SGE_TASK_ID)\"")

        }else{

          stop("Unkown mode type")
        }


}




#' Build execution innit for batch
#' 
#' @param envir Inherit enviroment variables
#' @param exec_code Execution code
#' @export


build_batch_exec_innit=function(
  envir=environment()
){  

    if(is.null(envir$exec_code)){
      stop("exec_code: Invalid argument type NULL")
    }

    envir$out_file_dir_batch=set_dir(dir=envir$out_file_dir_tmp,name="batch")

    envir$batch_code=build_job_exec(
      job=envir$job_id,time=envir$time,ram=envir$ram,
      threads=envir$threads,wd=getwd(),
      output_dir=envir$out_file_dir_batch,
      hold=envir$hold,array=length(envir$input)
    )

    envir$exec_code=paste0("echo '. $HOME/.bashrc;",envir$batch_config,
    ";",envir$exec_code,"'|",envir$batch_code)

}


#' Run job from batch constructor
#' 
#' @param envir Inherit current enviroment.
#' 
#' 
#' @export


run_self=function(
  envir=environment()
){  

    build_exec_innit(
            envir=envir
    )
  
    if(envir$mode=="batch"){
        build_batch_exec_innit(
              envir=envir
        )
    }

    if(envir$verbose){
          print_verbose(job=envir$job_id,
            arg=as.list(envir),
            exec_code=envir$exec_code
          )
    }

    envir$error=execute_job(exec_code=envir$exec_code)

    if(envir$error!=0){
        stop(envir$err_msg)
    }

}

#' Run job from batch constructor
#' 
#' @param envir Inherit current enviroment.
#' 
#' 
#' @export


 run_job=function(
  envir=environment()
 ){

      this.envir=environment()
      append_envir(this.envir,envir)


      
      if(verbose){
        print_verbose(
          job=job_id,arg=as.list(this.envir),
          exec_code=steps[[fn]]$exec_code
        )
      }


      steps[[fn]]$error=execute_job(
          exec_code=steps[[fn]]$exec_code
        )
      
      if(steps[[fn]]$error!=0){
          stop(err_mssg)
      }
      envir$steps <- steps
 }
   




#' Set default enviromental variables based on input for all
#' functions in package.
#' 
#' @param envir Inherit current enviroment.
#' @param id Input identifier. If not given name of input will be used.
#' @param name Name of output file directory
#' @export

set_envir_vars=function(
  envir=environment(),
  vars=NULL,
  fn=NULL,
  sheet=NULL,
  executor_id=NULL,
  err_mssg=NULL,
  ns="ULPwgs",
  dir_name=""
){

    this.envir=environment()
    append_envir(this.envir,envir)
    
    if(is.null(fn)){
        ## GET CALLER FUNCTION NAME IF NOT GIVEN
        
   
      this.envir$fn <- sub(".*::","",sub("\\(.*","",
          paste0(deparse(sys.calls()[[sys.nframe()-1]]),collapse=","))
        )
      }else{
      this.envir$fn <- fn
    }
      
    
    this.envir$ns <- ns

    if(is.null(executor_id)){
      this.envir$executor_id <- make_unique_id(this.envir$fn)
    }else{
      this.envir$executor_id <- executor_id
    }

    this.envir$task_id <- make_unique_id(this.envir$fn)

    this.envir$job_id <- build_job(executor_id=this.envir$executor_id,task_id=this.envir$task_id)

    if(is.null(err_mssg)){
        this.envir$err_msg <- paste0("CRITICAL ERROR: ",fn," -> ")
    }else{
      this.envir$err_msg <- paste0(this.envir$err_msg ,fn," -> ")
    }

    if (!is.null(sheet)){
        envir$envirs=set_ss_envir(envir=this.envir)
        return()
    }

    
    if(!is.null(inherit)){
        if(!is.environment(inherit)){
          inherit <-readRDS(file=inherit)
        }else{
          append_envir(envir,inherit)
        }
    }else{
      if(!is.null(inputs)){
        this.envir$inputs <- this.envir[[vars]]
        this.envir$n_inputs <- length(this.envir$inputs)
        this.envir$inputs_id <- set_input_id(inputs=this.envir$inputs,ids=this.envir$output_name)
        this.envir$inputs_ext <- Vectorize(get_file_ext)(this.envir$inputs)
      }

      this.envir$out_file_dir <- set_dir(
            dir=this.envir$output_dir,
            name=ifelse(this.envir$n_inputs>1,paste0(dir_name,"/",this.envir$n_inputs)
          )
      )

      this.envir$out_file_dir_tmp <- set_dir(dir=envir$out_file_dir,name="tmp")

      envir$envirs<-this.envir
    }

 
    return()
 
}



#' Set steps enviroment for use
#' 
#' @param envir Environment
#' @export


set_ss_envir=function(envir){
          dat=read.delim(envir$sheet,header=TRUE)
          dat_filt=dat %>% dplyr::distinct()
          nrows_dup=nrow(dat)-nrow(dat_filt)
          if(nrows_dup>0){
            warning(paste0(nrows_dup, " were duplicated in sheet"))
          }
          


          dat_filt=dat_filt %>% 
          dplyr::group_by(dplyr::across(-c(this.envir$vars))) %>%
          summarise(!! this.envir$vars := list(!! rlang::sym(this.envir$vars)))


          envirs=lapply(seq(1,nrow(dat_filt)),
            FUN=function(row){
              this.envir=environment()
              append_envir(this.envir,envir)
              this.envir$sheet <- NULL

              ### ASSIGN VARS IN SHEET TO ENVIROMENT
              invisible(lapply(seq(1,nrows(dat_row)),FUN=function(col){
                    this.envir[[col]] <- dat_row[row,col]
                    set_envir_vars(
                        envir=this.envir,
                        vars=this.envir$vars,
                        fn=this.envir$fn,
                        executor_id=this.envir$executor_id,
                        err_mssg=this.envir$err_msg
                    )
              }))

              return(this.envir)
            }
          )
          return(envirs)
}



#' Run enviroment
#' 
#' @param envir Environment
#' @export


run_envir=function(envir){

  this.envir=environment()
  append_envir(this.envir,envir)

  if(is.null(select)){
      runs=lapply(
          1:length(this.envir$envirs),FUN=function(n){
                  run_self(
                      envir=this.envir$envirs[n]
                  )
              }   
          )
          envir$runs <- runs
      }else{
          set_envir_inputs(envir=this.envir)
          run_main(
                  envir=this.envir
          )
          envir$steps <- steps
      }
}









#' Set steps enviroment for use
#' 
#' @param envir Environment
#' @export




set_steps_vars=function(
  envir=environment()
){      
      this.envir=environment()
      append_envir(this.envir,envir)
      steps=list()
      steps[[fn]]$job_id=job_id
      steps[[fn]]$executor_id=executor_id
      steps[[fn]]$task_id=task_id
      steps[[fn]]$select=select
      steps[[fn]]$input=input
      steps[[fn]]$input_id=input_id
      steps[[fn]]$exec_code=""
      steps[[fn]]$error=0
      steps[[fn]]$output_dir=out_file_dir
      steps[[fn]]$out_file=""
      return(steps)
  }



#' Append two enviroments
#' 
#' @param to Enviroment to append to
#' @param from Enviroment to append from
#' @export

append_envir = function(to=environment(), from=NULL) {
      
      from_list = ls(from)
      for(var in from_list) {
           if(!is.null(from[[var]])&!is.environment(from[[var]])){
              to[[var]] = from[[var]]
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
        my_id=Vectorize(get_file_name)(inputs)
      }
      return(my_id)
}

#' Set inputs withib enviroment from input list 
#' 
#' @param envir Enviroment to work within
#' @export


set_envir_inputs=function(envir){
      envir$input<-envir$inputs[envir$select]
      envir$input_id<-envir$inputs_id[envir$select]
      envir$input_ext<-envir$inputs_ext[envir$select]
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






