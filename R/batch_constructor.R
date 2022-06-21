build_job=function(executor="executor",task="task"){
options(scipen = 999)

  executor=paste0("E:",executor)
  job_name=paste0(c(executor,task),collapse="|")
  return(job_name)
}



build_job_exec=function(job="",time="48:0:0",ram=3,threads=1,output_dir="",hold=""){
  
  if(hold!=""){
    hold=paste0(" -hold_jid",paste0(hold,collapse=","))
  }
  exec_code=paste("qsub -N ",job, paste0(" -l h_rt=",time),
  paste0(" -l mem=",ram,"G"), paste0(" -pe smp ",threads), paste0(" -wd ",getwd()), 
  paste0(" -o ",output_dir,"/",full_name,".std_out"),
  paste0(" -e ",output_dir,"/",full_name,".std_error"),hold)
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
          print("----------------------------------")
          print(dat_info)
          print("----------------------------------")
      }
      if(any(grepl("E",dat_info$status))){
        error=TRUE
      }
      Sys.sleep(time)
    }
  }

  if(error){
      print("----------------------------------")
      print(dat_info)
      print("----------------------------------")
      parallel::mclapply(1:nrow(dat_info),FUN=function(x){
        system(paste0("qdel ",dat_info[x,]$job_id))
      },mc.cores=threads)
    
      stop("One or more jobs failed")
  }

return()
    

}

make_unique_id=function(name,id=sample(1:100000000000000,1),sep="_"){
  options(scipen = 999)
  unique_name=paste0(c(name,id),sep=sep)
  return(unique_name)
}