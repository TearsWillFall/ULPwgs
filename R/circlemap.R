
#' Extract Circular DNA Candidates using Circlemap Realign tool
#'
#' This function generates a BED file with circular DNA candidates
#' 
#' 
#' //TODO validate co-joint calling mode
#' 
#' For more information read:
#' https://github.com/iprada/Circle-Map
#' 
#' 
#' 
#' @param env_circlemap [REQUIRED] Conda enviroment for circlemap tool.
#' @param bam [OPTIONAL] Path to BAM file
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param threads [OPTIONAL] Number of threads to split the work. Default 4
#' @param ram [OPTIONAL] RAM memory to asing to each thread. Default 4
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param batch_config [REQUIRED] Additional batch configuration if batch mode selected.
#' @param executor_id Task EXECUTOR ID. Default "recalCovariates"
#' @param task_name Task name. Default "recalCovariates"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export


realign_circlemap=function(
        env_circlemap=build_default_python_enviroment_list()$env_circlemap,
        bin_samtools=build_default_tool_binary_list()$bin_samtools,
        bam=NULL,output_dir=".",verbose=FALSE,
        batch_config=build_default_preprocess_config(),
        threads=3,ram=1, mode="local",
        tmp_dir=NULL,
        executor_id=make_unique_id("realignCircleMap"),
        task_name="realignCircleMap",time="48:0:0",
        update_time=60,wait=FALSE,hold=NULL
){

    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir,name="realign_call")

    if(!is.null(tmp_dir)){
        out_file_dir_tmp=tmp_dir
    }else {
        out_file_dir=set_dir(dir=out_file_dir,name="tmp")
    }  

    job=build_job(executor_id=executor_id,task_id=task_id)

    if(is.null(bam)){
        stop("bam argument is required")
    }

    id=""

    if(output_name!=""){
        id=output_name
    }else{
        id=get_file_name(bam[1])
    }



    jobs_report=build_job_report(
        job_id=job,
        executor_id=executor_id,
        exec_code=list(), 
        task_id=task_id,
        input_args = argg,
        out_file_dir=out_file_dir,
        out_files=list(
            circ_bed=out_file
        )
    )



    jobs_report[["steps"]][["sort_bam_name"]]<-sort_and_index_bam_samtools(
            bin_samtools=bin_samtools,
            bam=bam,output_dir=out_file_dir_tmp,verbose=verbose,
            batch_config=batch_config,threads=threads,ram=ram,sort=TRUE,
            coord_sort=FALSE,index=FALSE,stats="all", clean=TRUE,
            mode=mode,executor_id=task_id,
            time=time,
            hold=hold
    )

    jobs_report[["steps"]][["extract_circular_reads"]]<- read_extractor_circlemap(
                env_circlemap=env_circlemap,
                bin_samtools=bin_samtools,
                bam=jobs_report[["steps"]][["sort_and_index"]][["steps"]][["sort"]]$out_files$bam,
                output_dir=out_file_dir_tmp,verbose=verbose,
                batch_config=batch_config,
                threads=threads,ram=ram,
                sort=TRUE,
                coord_sort=TRUE,
                index=TRUE,stats="all", clean=TRUE,
                mode=mode,executor_id=task_id,
                time=time,
                hold=hold
    )



    hold=unlist_lvl(job_report[["steps"]][["extract_circular_reads"]],var="job_id")


    out_file=paste0(out_file_dir,"/",id,".circular_candidates.bed")
    exec_code=paste(paste("conda ",env_circlemap),";Circle-Map Realign -sbam ",bam, " -qbam ", 
    jobs_report[["steps"]][["sort_and_index"]][["steps"]][["sort"]]$out_files$bam," -i ",
    jobs_report[["steps"]][["extract_circular_reads"]][["steps"]][["sort_and_index"]][["steps"]][["sort"]]$out_files$bam,
    " -o ",out_file," -t ",threads," -dir ",out_file_dir)
    
     
    if(mode=="batch"){
        out_file_dir2=set_dir(dir=out_file_dir,name="batch")
        batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
        threads=threads,output_dir=out_file_dir2)
        exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,
        ";",exec_code,"'|",batch_code)
    }


    jobs_report$exec_code=exec_code

    if(verbose){
        print_verbose(job=job,arg=argg,exec_code=exec_code)
    }

    error=execute_job(exec_code=exec_code)
    
    
    if(error!=0){
        stop("circlemap failed to run due to unknown error.
        Check std error for more information.")
    }



    if(wait&&mode=="batch"){
        job_validator(
            job=job_report$job_id,
            time=update_time,
            verbose=verbose,
            threads=threads
        )
    }

    return(jobs_report)
}




#' Extract Circular Reads using Circlemap tool
#'
#' This function generates a BAM file with circular DNA candidate reads
#' 
#' 
#' //TODO validate co-joint calling mode
#' 
#' For more information read:
#' https://github.com/iprada/Circle-Map
#' 
#' 
#' 
#' @param env_circlemap [REQUIRED] Conda enviroment for circlemap tool.
#' @param bin_samtools [REQUIRED] Path to samtools binary file.
#' @param bam [OPTIONAL] Path to BAM file
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param threads [OPTIONAL] Number of threads to split the work. Default 4
#' @param ram [OPTIONAL] RAM memory to asing to each thread. Default 4
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param sort Sort BAM file. Default TRUE.
#' @param coord_sort Generate a coord sorted file. Otherwise queryname sorted. Default TRUE
#' @param index Generate an index file for sorted BAM. Default TRUE
#' @param clean Remove input files. Default FALSE
#' @param stats Generate BAM stats. Default all. Options ["all","flag","index",""]
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param batch_config [REQUIRED] Additional batch configuration if batch mode selected.
#' @param executor_id Task EXECUTOR ID. Default "recalCovariates"
#' @param task_name Task name. Default "recalCovariates"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export


read_extractor_circlemap=function(
    env_circlemap=build_default_python_enviroment_list()$env_circlemap,
    bin_samtools=build_default_tool_binary_list()$bin_samtools,
    bam=NULL,output_dir=".",output_name="",verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=3,ram=1,sort=TRUE,coord_sort=FALSE,
    index=FALSE,stats="all", clean=FALSE,
    mode="local",executor_id=make_unique_id("readExtractCircleMap"),
    task_name="readExtractCircleMap",time="48:0:0",
    update_time=60,wait=FALSE,hold=NULL

  ){

    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir)
    job=build_job(executor_id=executor_id,task_id=task_id)

    if(is.null(bam)){
        stop("bam argument is required")
    }

    id=""
    if(output_name!=""){
        id=output_name
    }else{
        id=get_file_name(bam[1])
    }

    out_file=paste0(out_file_dir,"/",id,".circular_read_candidates.bam")
    exec_code=paste(paste("conda ",env_circlemap),";Circle-Map ReadExtractor -i ",
    bam, " -o ", out_file," -dir ",out_file_dir)
    
    if(mode=="batch"){
        out_file_dir2=set_dir(dir=out_file_dir,name="batch")
        batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
        threads=threads,output_dir=out_file_dir2)
        exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,
        ";",exec_code,"'|",batch_code)
    }

    if(verbose){
        print_verbose(job=job,arg=argg,exec_code=exec_code)
    }

    error=execute_job(exec_code=exec_code)
    
    
    if(error!=0){
        stop("circlemap failed to run due to unknown error.
        Check std error for more information.")
    }

    jobs_report=build_job_report(
        job_id=job,
        executor_id=executor_id,
        exec_code=exec_code, 
        task_id=task_id,
        input_args = argg,
        out_file_dir=out_file_dir,
        out_files=list(
            circ_bam=out_file
        )
    )

    jobs_report[["steps"]][["sort_bam_name"]]<-sort_and_index_bam_samtools(
                    bin_samtools=bin_samtools,
                    bam=out_file,
                    output_dir=out_file_dir,
                    verbose=verbose,
                    batch_config=batch_config,
                    threads=threads,ram=ram,sort=sort,
                    coord_sort=coord_sort,index=index,
                    stats=stats,clean=clean,
                    mode=mode,executor_id=task_id,
                    time=time,
                    hold=job
    )


    if(wait&&mode=="batch"){
        job_validator(
            job=job_report$job_id,
            time=update_time,
            verbose=verbose,
            threads=threads
        )
    }

    return(jobs_report)

}


#' Extract Circular DNA Repeat Candidates using Circlemap tool
#'
#' This function generates a BED file with circular DNA repeat candidates
#' 
#' 
#' //TODO validate co-joint calling mode
#' 
#' For more information read:
#' https://github.com/iprada/Circle-Map
#' 
#' 
#' 
#' @param env_circlemap [REQUIRED] Conda enviroment for circlemap tool.
#' @param bam [OPTIONAL] Path to BAM file
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param threads [OPTIONAL] Number of threads to split the work. Default 4
#' @param ram [OPTIONAL] RAM memory to asing to each thread. Default 4
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param batch_config [REQUIRED] Additional batch configuration if batch mode selected.
#' @param executor_id Task EXECUTOR ID. Default "recalCovariates"
#' @param task_name Task name. Default "recalCovariates"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export




repeat_caller_circlemap=function(
    env_circlemap=build_default_python_enviroment_list()$env_circlemap,
    bam=NULL,output_dir=".",verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=1,
    mode="local",
    executor_id=make_unique_id("RepeatsCircleMap"),
    task_name="RepeatsCircleMap",time="48:0:0",
    update_time=60,wait=FALSE,hold=NULL
){

    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir,name="repeat_call")
    job=build_job(executor_id=executor_id,task_id=task_id)

    if(is.null(bam)){
        stop("bam argument is required")
    }

    id=""
    if(output_name!=""){
        id=output_name
    }else{
        id=get_file_name(bam[1])
    }


    out_file=paste0(out_file_dir,"/",id,".circular_repeat_candidates.bed")
    exec_code=paste(paste("conda ",env_circlemap),";Circle-Map Repeats -i ", bam, " -o ",
     out_file, " -dir ",out_file_dir)
    
    if(mode=="batch"){
        out_file_dir2=set_dir(dir=out_file_dir,name="batch")
        batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
        threads=threads,output_dir=out_file_dir2)
        exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,
        ";",exec_code,"'|",batch_code)
    }

    if(verbose){
        print_verbose(job=job,arg=argg,exec_code=exec_code)
    }

    error=execute_job(exec_code=exec_code)
    
    
    if(error!=0){
        stop("circlemap failed to run due to unknown error.
        Check std error for more information.")
    }

    jobs_report=build_job_report(
        job_id=job,
        executor_id=executor_id,
        exec_code=exec_code, 
        task_id=task_id,
        input_args = argg,
        out_file_dir=out_file_dir,
        out_files=list(
            circ_repeat_bed=out_file
        )
    )

    if(wait&&mode=="batch"){
        job_validator(
            job=job_report$job_id,
            time=update_time,
            verbose=verbose,
            threads=threads
        )
    }

    return(jobs_report)

}



#' Extract Circular DNA Candidates using Circlemap tool
#'
#' This function generates a BED file with circular DNA candidates
#' 
#' 
#' //TODO validate co-joint calling mode
#' 
#' For more information read:
#' https://github.com/iprada/Circle-Map
#' 
#' 
#' 
#' @param env_circlemap [REQUIRED] Conda enviroment for circlemap tool.
#' @param bam [OPTIONAL] Path to BAM file
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param threads [OPTIONAL] Number of threads to split the work. Default 4
#' @param ram [OPTIONAL] RAM memory to asing to each thread. Default 4
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param batch_config [REQUIRED] Additional batch configuration if batch mode selected.
#' @param executor_id Task EXECUTOR ID. Default "recalCovariates"
#' @param task_name Task name. Default "recalCovariates"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export


circdna_circlemap=function(
        env_circlemap=build_default_python_enviroment_list()$env_circlemap,
        bin_samtools=build_default_tool_binary_list()$bin_samtools,
        bam=NULL,output_dir=".",output_name="",verbose=FALSE,
        batch_config=build_default_preprocess_config(),
        threads=3,ram=1, mode="local",
        executor_id=make_unique_id("circdnaCircleMap"),
        task_name="circdnaCircleMap",time="48:0:0",
        update_time=60,wait=FALSE,hold=NULL
){

    id=""

    if(output_name!=""){
        id=output_name
    }else{
        id=get_file_name(bam[1])
    }



    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir,name=paste0(id,"/circlemap"))
    out_file_dir_tmp=set_dir(dir=out_file_dir,name="tmp")
    job=build_job(executor_id=executor_id,task_id=task_id)

    if(is.null(bam)){
        stop("bam argument is required")
    }



    jobs_report=build_job_report(
        job_id=job,
        executor_id=executor_id,
        exec_code=list(), 
        task_id=task_id,
        input_args = argg,
        out_file_dir=out_file_dir,
        out_files=list(
        )
    )

    jobs_report[["steps"]][["realign_circlemap"]]<-realign_circlemap(
        env_circlemap=env_circlemap,
        bin_samtools=bin_samtools,
        bam=bam,output_dir=out_file_dir,verbose=verbose,
        tmp_dir=out_file_dir_tmp,
        batch_config=batch_config,
        threads=threads,ram=ram, mode=mode,
        executor_id=task_id,
        time=time,
        hold=hold
    )



    jobs_report[["steps"]][["repeats_circlemap"]]<-repeat_caller_circlemap(
        env_circlemap=env_circlemap,
        bam=bam,output_dir=out_file_dir,verbose=verbose,
        batch_config=batch_config,
        threads=threads,ram=ram,
        mode=mode,
        executor_id=task_id,
        time=time,
        hold=hold
    )



 
    if(wait&&mode=="batch"){
        job_validator(
            job=job_report$job_id,
            time=update_time,
            verbose=verbose,
            threads=threads
        )
    }




    return(jobs_report)


}





