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
    " -n " , threads, " -m ", ram, 
    "-f .+_R1_.+[.]f.+[.]gz,.+_1[.].+[.]gz,R1[.].+[.]gz -r .+_R2_.+[.]f.+[.]gz,.+_2[.].+[.]gz,R2[.].+[.]gz")
    
    
    out_file_dir2=set_dir(dir=out_file_dir,name="batch")
    job=build_job(executor_id=executor_id,task_id=make_unique_id(task_id))

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



#' Process multiple tumour samples using CLONET
#'
#' This function identifies a set of BAM files as tumour and normal
#' and processes them using the CLONET pipeline
#'
#' @param sif_path Path to singularity image file. 
#' @param bam_dir Path to bam directory. 
#' @param normal_id Normal sample identifier. 
#' @param patient_id Patient id. 
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Executor ID. Default "preprocess_trento"
#' @param task_name Name of the task. Default "preprocess_trento"
#' @param threads Number of CPU cores to use. Default 3.
#' @param ram RAM memory for batched job. Default 4
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param output_dir Path to the output directory.
#' @param tmp_dir Path to temporary file directory.
#' @param verbose Enables progress messages. Default False.
#' @param hold Job to hold on in batched mode.
#' @export



multisample_clonet_trento=function(
    sif_path="Singularity_Images/preProcess_latest.sif",
    bam_dir="",normal_id="",patient_id="",tmp_dir="",threads=3,
    ram=4,output_dir="",verbose=FALSE,
    executor_id=make_unique_id("multi_clonet"),
    task_name="multi_clonet",mode="local",time="48:0:0",
    update_time=60,wait=FALSE,hold=""){

        argg <- as.list(environment())

        task_id=make_unique_id(task_name)
        out_file_dir=set_dir(dir=output_dir)

        job=build_job(executor_id=executor_id,task_id=task_id)

        job_report=build_job_report(
                job_id=job,
                executor_id=executor_id,
                task_id=task_id,
                input_args=argg,
                out_file_dir=out_file_dir,
                out_files=list(
                )
        )
        bam_dir_path=system(paste("realpath",bam_dir),intern=TRUE)
        bam_files=system(paste0("find ",bam_dir_path,"| grep bam$"),intern=TRUE)
        t_files=bam_files[!grepl(normal_id,bam_files)]
        normal=bam_files[grepl(normal_id,bam_files)]

        parallel::mclapply(t_files,FUN=function(tumour){
            job_report[["steps"]][["clonet"]][[ULPwgs::get_file_name(tumour)]]<<-clonet_trento(
                sif_path=sif_path,
                tumour=tumour,normal=normal,
                patient_id=patient_id,
                threads=threads,
                ram=ram,output_dir=paste0(out_file_dir,patient_id),verbose=verbose,
                executor_id=task_id,mode=local,time=time,
                update_time=60,wait=FALSE,hold="")
        },mc.cores=ifelse(mode=="local",1,3))


    if(wait&&mode=="batch"){
        job_validator(job=unlist_level(named_list=job_report[["steps"]][["clonet"]],var="job_id"),
        time=update_time,verbose=verbose,threads=threads)
    }

    return(job_report)
    }





#' Process a pair of tumour-normal samples using CLONET
#'
#' This function takes a pair of tumour of tumour and 
#' normal BAMS and applies the CLONET pipeline
#'
#'
#' @param sif_path Path to singularity image file. 
#' @param tumour Path to tumour BAM file 
#' @param normal Path to normal BAM file
#' @param patient_id Patient id. 
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Executor ID. Default "preprocess_trento"
#' @param task_name Name of the task. Default "preprocess_trento"
#' @param threads Number of CPU cores to use. Default 3.
#' @param ram RAM memory for batched job. Default 4
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param output_dir Path to the output directory.
#' @param tmp_dir Path to temporary file directory.
#' @param verbose Enables progress messages. Default False.
#' @param hold Job to hold on in batched mode.
#' @export




clonet_trento=function(
    sif_path="Singularity_Images/preProcess_latest.sif",
    tumour="",normal="",patient_id="",tmp_dir="",threads=3,
    ram=4,output_dir="",verbose=FALSE,
    executor_id=make_unique_id("clonet"),
    task_name="clonet",mode="local",time="48:0:0",
    update_time=60,wait=FALSE,hold=""){

    argg <- as.list(environment())

    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir,name="clonet_reports")
    out_file_dir_tmp=set_dir(dir=output_dir,name="clone_tmp")
    out_file_dir=set_dir(dir=out_file_dir,name=get_file_name(tumour))


    file_info=data.frame(Patient=patient_id,Tumour=tumour,Normal=normal)

    sample_sheet=paste0(getwd(),"/",out_file_dir_tmp,"/tmp.txt")
    write.table(file_info,file=sample_sheet,quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")

    exec_code=paste(paste0(" export SINGULARITY_BINDPATH=",getwd()),"; singularity run --app PCFS ",
    sif_path, " -s ", sample_sheet ," -o ", paste0(getwd(),"/",out_file_dir)," -t ",out_file_dir_tmp,
    " -n " , threads)
    
    out_file_dir2=set_dir(dir=out_file_dir,name="batch")
    job=build_job(executor_id=executor_id,task_id=make_unique_id(task_id))

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
        stop("Clonet failed to run due to unknown error.
        Check std error for more information.")
    }
    
    if(wait&&mode=="batch"){
        batch_validator(job=job_report$job_id,
        time=update_time,verbose=verbose,threads=threads)
    }

    return(job_report)
}
