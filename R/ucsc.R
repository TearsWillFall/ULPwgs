
#' Function to convert bigBed to bed file
#'
#' This function takes a bigBed file and converts it to a bed file
#'
#' @param bin_ucsc Path to bigBedtoBed script. Default path tools/OtherToolsr/bigBedToBed.
#' @param bigBed Path to directory with BAM files to merge.
#' @param header Create column header. Default TRUE
#' @param verbose Enables progress messages. Default False.
#' @param threads Number of threads . Default 4
#' @param ram RAM memory. Default 4
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Task EXECUTOR ID. Default "mardupsGATK"
#' @param task_name Task name. Default "mardupsGATK"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param verbose [OPTIONAL] Enables progress messages. Default False.#
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] Hold job until job is finished. Job ID. 
#' @export

bigBedToBed_ucsc=function(
    rdata=NULL,
    selected=NULL,
    bin_ucsc=build_default_tool_binary_list()$bin_ucsc$bigBedToBed,
    bigBed="",output_name="",
    output_dir=".",verbose=FALSE,batch_config=build_default_preprocess_config(),
    executor_id=make_unique_id("bigBedToBed"),task_name="bigBedToBed",
    mode="local",time="48:0:0",
    threads=4,ram=4,update_time=60,wait=FALSE,hold=NA
  ){
 
    options(scipen = 999)
    

    if(!is.null(rdata)){
        load(rdata)
        if(!is.null(selected)){
            bigBed=bigBed_list[selected]
        }
    }


    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir)

    job=build_job(executor_id=executor_id,task_id=task_id)


    id=""
    if(output_name!=""){
      id=output_name
    }else{
      id=get_file_name(bigBed)
    }


    out_file=paste0(out_file_dir,id,".bed")
    exec_code=paste(bin_ucsc, bigBed,out_file)
    
    
    
    if(mode=="batch"){
        out_file_dir2=set_dir(dir=out_file_dir,name="batch")
        batch_code=build_job_exec(job=job,
        time=time,ram=ram,threads=threads,
        output_dir=out_file_dir2,hold=hold)
        exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
    }


    if(verbose){
      print_verbose(job=job,arg=argg,exec_code=exec_code)
    }


  error=execute_job(exec_code=exec_code)
  if(error!=0){
    stop("ucsc tools failed to run due to unknown error.
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
        bed=out_file)
  )




if(wait&&mode=="batch"){
    batch_validator(job=job_report$job_id,
    time=update_time,verbose=verbose,threads=threads)
  }


  return(job_report)
}





#' Parallel implementation of bigBedToBed function
#'
#' This function takes a bigBed file and converts it to a bed file
#'
#' @param bin_ucsc Path to bigBedtoBed script. Default path tools/OtherToolsr/bigBedToBed.
#' @param bigBeds Path to directory with BAM files to merge.
#' @param header Create column header. Default TRUE
#' @param verbose Enables progress messages. Default False.
#' @param threads Number of threads . Default 4
#' @param ram RAM memory. Default 4
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Task EXECUTOR ID. Default "mardupsGATK"
#' @param task_name Task name. Default "mardupsGATK"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param verbose [OPTIONAL] Enables progress messages. Default False.#
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] Hold job until job is finished. Job ID. 
#' @export

parallel_sample_bigBedToBed_ucsc=function(
    bin_ucsc=build_default_tool_binary_list()$bin_ucsc$bigBedToBed,
    bigBeds="",
    output_dir=".",verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    executor_id=make_unique_id("parSampleBigBedToBed"),
    task_name="parSampleBigBedToBed",
    mode="local",time="48:0:0",
    threads=4,ram=4,update_time=60,wait=FALSE,hold=NA
  ){
 
    options(scipen = 999)
    
    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir)

    job=build_job(executor_id=executor_id,task_id=task_id)

    job_report=build_job_report(
        job_id=job,
        executor_id=executor_id,
        exec_code=list(), 
        task_id=task_id,
        input_args=argg,
        out_file_dir=out_file_dir,
        out_files=list(
      )
    )


    
    bigBed_list=bigBeds
    names(bigBed_list)=Vectorize(get_file_name)(bigBeds)

    if(mode=="local"){
            job_report[["steps"]][["parSampleBigBedToBed"]]<-parallel::mclapply(
                bigBed_list,FUN=function(bed){
                job_report <- bigBedToBed_ucsc(
                    bigBed=bed,batch_config=batch_config,
                    bin_ucsc=bin_ucsc,
                    output_name=get_file_name(bed),
                    output_dir=out_file_dir,
                    verbose=verbose,
                    executor_id=task_id
            )
        },mc.cores=threads)

    }else if(mode=="batch"){
            rdata_file=paste0(tmp_dir,"/",job,".samples.RData")
            output_dir=out_file_dir
            save(bigBed_list,output_dir,verbose,file = rdata_file)
            exec_code=paste0("Rscript -e \"ULPwgs::bigBedToBed_ucsc(rdata=\\\"",
            rdata_file,"\\\",selected=$SGE_TASK_ID)\"")
            out_file_dir2=set_dir(dir=out_file_dir,name="batch")
            batch_code=build_job_exec(job=job,time=time,ram=ram,
            threads=1,output_dir=out_file_dir2,
            hold=hold,array=length(bigBed_list))
            exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)


            if(verbose){
                print_verbose(job=job,arg=argg,exec_code=exec_code)
            }
            error=execute_job(exec_code=exec_code)
            if(error!=0){
                stop("ucsc tools failed to run due to unknown error.
                Check std error for more information.")
            }

            job_report[["steps"]][["parSampleBigBedToBed"]]<- build_job_report(
                job_id=job,
                executor_id=executor_id,
                exec_code=exec_code, 
                task_id=task_id,
                input_args=argg,
                out_file_dir=out_file_dir,
                out_files=list(
                    bed=paste0(out_file_dir,names(bigBed_list),".bed")
                    )
            )
    }

  return(job_report)
}




#' Function to liftOver files between versions
#'
#' This function takes a bigBed file and converts it to a bed file
#'
#' @param bin_ucsc Path to bigBedtoBed script. Default path tools/OtherToolsr/bigBedToBed.
#' @param bigBed Path to directory with BAM files to merge.
#' @param header Create column header. Default TRUE
#' @param verbose Enables progress messages. Default False.
#' @param threads Number of threads . Default 4
#' @param ram RAM memory. Default 4
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Task EXECUTOR ID. Default "mardupsGATK"
#' @param task_name Task name. Default "mardupsGATK"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param verbose [OPTIONAL] Enables progress messages. Default False.#
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] Hold job until job is finished. Job ID. 
#' @export

liftOver_ucsc=function(
    rdata=NULL,
    selected=NULL,
    bin_ucsc=build_default_tool_binary_list()$bin_ucsc$bigBedToBed,
    bigBed="",output_name="",
    output_dir=".",verbose=FALSE,batch_config=build_default_preprocess_config(),
    executor_id=make_unique_id("bigBedToBed"),task_name="bigBedToBed",
    mode="local",time="48:0:0",
    threads=4,ram=4,update_time=60,wait=FALSE,hold=NA
  ){
 
    options(scipen = 999)
    

    if(!is.null(rdata)){
        load(rdata)
        if(!is.null(selected)){
            bigBed=bigBed_list[selected]
        }
    }


    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir)

    job=build_job(executor_id=executor_id,task_id=task_id)


    id=""
    if(output_name!=""){
      id=output_name
    }else{
      id=get_file_name(bigBed)
    }


    out_file=paste0(out_file_dir,id,".bed")
    exec_code=paste(bin_ucsc, bigBed,out_file)
    
    
    
    if(mode=="batch"){
        out_file_dir2=set_dir(dir=out_file_dir,name="batch")
        batch_code=build_job_exec(job=job,
        time=time,ram=ram,threads=threads,
        output_dir=out_file_dir2,hold=hold)
        exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
    }


    if(verbose){
      print_verbose(job=job,arg=argg,exec_code=exec_code)
    }


  error=execute_job(exec_code=exec_code)
  if(error!=0){
    stop("ucsc tools failed to run due to unknown error.
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
        bed=out_file)
  )




if(wait&&mode=="batch"){
    batch_validator(job=job_report$job_id,
    time=update_time,verbose=verbose,threads=threads)
  }


  return(job_report)
}

