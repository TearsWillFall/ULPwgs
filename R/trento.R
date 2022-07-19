#' Generate a quality control (QC) report from a fastaqc file
#'
#' This function takes a set of sequence files (fastq,SAM,BAM...) and
#' returns a report in HTML format.
#'
#'
#' @param sif_path Path to singularity image file. 
#' @param fastq_dir Path to fastq directory. 
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Executor ID. Default "preprocess_trento"
#' @param task_name Name of the task. Default "preprocess_trento"
#' @param threads Number of CPU cores to use. Default 3.
#' @param ram RAM memory for batched job. Default 4
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param hold Job to hold on in batched mode.
#' @export



preprocess_seq_trento=function(
    sif_path="Singularity_Images/preProcess_latest.sif",
    fastq_dir="", threads=3,ram=4,ref_genome="",output_dir="",verbose=FALSE,
    executor_id=make_unique_id("preprocess_trento"),tmp_dir="",
    task_name="reprocess_trento",mode="local",time="48:0:0",
    update_time=60,wait=FALSE,hold=""){

    argg <- as.list(environment())

    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir,name="preprocess")
    
 
    exec_code=paste(" singularity run --app preProcess ",
    sif_path, " -i ", fastq_dir ," -o ", out_file_dir," -t ",tmp_dir,
    " -n " , threads, " -m ", ram)
    
    
    out_file_dir2=set_dir(dir=out_file_dir,name="batch")
    job=build_job(executor=executor,task=make_unique_id(task))

    if(mode=="batch"){
        exec_batch=build_job_exec(job=job,hold=hold,time=time,ram=ram,threads=threads,
        output_dir=out_file_dir2)
        exec_code=paste("echo 'source ~/.bashrc;",exec_code,"'|",exec_batch)
    }

    if(verbose){
        print_verbose(job=job,arg=argg,exec_code=exec_code)
    }


    job_report=build_job_report(
        job_id=job,
        executor_id=executor_id,
        task_id=task_id,
        input_args=argg,
        out_file_dir=out_file_dir,
        out_files=list(
            NA
        )
    )

    error=system(exec_code)
    if(error!=0){
        stop("Preprocess failed to run due to unknown error.
        Check std error for more information.")
    }
    
    if(wait&&mode=="batch"){
        batch_validator(job=job_report$job_id,
        time=update_time,verbose=verbose,threads=threads)
    }

    return(job_report)


}
