
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
        bam=NULL,output_dir=".",output_name="",verbose=FALSE,
        ref_genome=build_default_reference_list()$HG19$reference$genome,
        batch_config=build_default_preprocess_config(),
        threads=3,ram=1, mode="local",
        tmp_dir=NULL,
        executor_id=make_unique_id("realignCircleMap"),
        task_name="realignCircleMap",time="48:0:0",
        update_time=60,wait=FALSE,hold=NULL
){

    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir,name="realign_reports")

    if(!is.null(tmp_dir)){
        out_file_dir_tmp=tmp_dir
    }else {
        out_file_dir_tmp=set_dir(dir=out_file_dir,name="tmp")
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

    
    out_file=paste0(out_file_dir,"/",id,".circular_candidates.bed")

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



    jobs_report[["steps"]][["sort_and_index"]]<-sort_and_index_bam_samtools(
            bin_samtools=bin_samtools,
            bam=bam,output_dir=out_file_dir_tmp,verbose=verbose,
            batch_config=batch_config,threads=threads,ram=ram,sort=TRUE,
            coord_sort=FALSE,index=FALSE,stats="all", clean=FALSE,
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



    hold=unlist_lvl(jobs_report[["steps"]][["extract_circular_reads"]],var="job_id")


    exec_code=paste(set_conda_enviroment(env_circlemap),
    "Circle-Map Realign -sbam ",normalizePath(bam), " -qbam ", 
    jobs_report[["steps"]][["sort_and_index"]][["steps"]][["sort"]]$out_files$bam," -i ",
    jobs_report[["steps"]][["extract_circular_reads"]][["steps"]][["sort_and_index"]][["steps"]][["sort"]]$out_files$bam,
    " -o ",out_file," -t ",threads," -dir /", " -fasta ",normalizePath(ref_genome), " -tdir ",out_file_dir_tmp)
    
     
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
    exec_code=paste(set_conda_enviroment(env_circlemap),"Circle-Map ReadExtractor -i ",
    normalizePath(bam), " -o ", out_file," -dir /")
    
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

    jobs_report[["steps"]][["sort_and_index"]]<-sort_and_index_bam_samtools(
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
    bam=NULL,output_dir=".",output_name="",verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=1,
    mode="local",
    executor_id=make_unique_id("RepeatsCircleMap"),
    task_name="RepeatsCircleMap",time="48:0:0",
    update_time=60,wait=FALSE,hold=NULL
){

    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir,name="repeat_reports")
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
    exec_code=paste(set_conda_enviroment(env_circlemap),"Circle-Map Repeats -i ",
     normalizePath(bam), " -o ",
    out_file, " -dir /")
    




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
#' @param bin_samtools [REQUIRED] Path to samtools binary file.
#' @param bam [OPTIONAL] Path to BAM file
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param tmp_dir [OPTIONAL] Path to the temporary directory.
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
        ref_genome=build_default_reference_list()$HG19$reference$genome,
        bam=NULL,output_dir=".",tmp_dir=NULL,output_name="",verbose=FALSE,
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
    out_file_dir=set_dir(dir=output_dir,name=paste0(id,"/circlemap_reports"))
    out_file_dir_tmp=set_dir(dir=out_file_dir,name="tmp")
    job=build_job(executor_id=executor_id,task_id=task_id)


    if(!is.null(tmp_dir)){
        out_file_dir_tmp=tmp_dir
    }else {
        out_file_dir_tmp=set_dir(dir=out_file_dir,name="tmp")
    }  

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
        bam=normalizePath(bam),
        ref_genome=normalizePath(ref_genome),
        output_dir=out_file_dir,verbose=verbose,
        tmp_dir=out_file_dir_tmp,
        batch_config=batch_config,
        threads=threads,ram=ram, mode=mode,
        executor_id=task_id,
        time=time,
        hold=hold
    )


    
    jobs_report[["steps"]][["repeats_circlemap"]]<-repeat_caller_circlemap(
        env_circlemap=env_circlemap,
        bam=normalizePath(bam),output_dir=out_file_dir,verbose=verbose,
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
#' @param tmp_dir [OPTIONAL] Path to the temporary directory.
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




parallel_samples_circdna_circlemap=function(
        env_circlemap=build_default_python_enviroment_list()$env_circlemap,
        bin_samtools=build_default_tool_binary_list()$bin_samtools,
        ref_genome=build_default_reference_list()$HG19$reference$genome,
        bam=NULL,output_dir=".",output_name="",verbose=FALSE,
        patient_id=NULL,
        batch_config=build_default_preprocess_config(),
        threads=3,ram=1, mode="local",
        executor_id=make_unique_id("parcircdnaCircleMap"),
        task_name="parcircdnaCircleMap",time="48:0:0",
        update_time=60,wait=FALSE,hold=NULL
    ){


        argg <- as.list(environment())
        task_id=make_unique_id(task_name)
        out_file_dir=set_dir(dir=output_dir,name=patient_id)
        out_file_dir_tmp=set_dir(dir=out_file_dir,name="tmp")
        job=build_job(executor_id=executor_id,task_id=task_id)

        if(is.null(bam)){
            stop("bam argument is required")
        }

        if(is.null(patient_id)){
            stop("patient_id argument is required")
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


    bam_list=bam
    names(bam_list)=Vectorize(get_file_name)(bam)

    if(mode=="local"){
      jobs_report[["steps"]][["par_sample_call_circdna"]]<-
      lapply(bam_list,FUN=function(b){
            job_report <-circdna_circlemap(
            env_circlemap=env_circlemap,
            bin_samtools=bin_samtools,
            ref_genome=normalizePath(ref_genome),
            bam=normalizePath(bam),output_dir=out_file_dir,
            output_name=get_file_name(bam),
            tmp_dir=out_file_dir_tmp,
            verbose=verbose,
            batch_config=batch_config,
            threads=threads,ram=ram,
            mode=mode,
            executor_id=task_id,
            time=time,
            hold=hold
    )
      })
    }else if(mode=="batch"){
            rdata_file=paste0(tmp_dir,"/",job,".samples.RData")
            output_dir=out_file_dir
            executor_id=task_id
            tmp_dir=out_file_dir_tmp
            save(
              bam_list,
              bin_samtools,
              env_circlemap,
              normalizePath(ref_genome),
              output_dir,
              tmp_dir,
              executor_id,
              clean,
              verbose
            )
            exec_code=paste0("Rscript -e \"ULPwgs::circdna_circlemap(rdata=\\\"",
            rdata_file,"\\\",selected=$SGE_TASK_ID)\"")
            out_file_dir2=set_dir(dir=out_file_dir,name="batch")
            batch_code=build_job_exec(job=job,time=time,ram=ram,
            threads=2,output_dir=out_file_dir2,
            hold=hold,array=length(bam_list))
            exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)

            if(verbose){
                print_verbose(job=job,arg=argg,exec_code=exec_code)
            }
            error=execute_job(exec_code=exec_code)
            if(error!=0){
                stop("circlemap failed to run due to unknown error.
                Check std error for more information.")
            }
    
          jobs_report[["steps"]][["par_sample_call_circdna"]]<- build_job_report(
                job_id=job,
                executor_id=executor_id,
                exec_code=exec_code, 
                task_id=task_id,
                input_args=argg,
                out_file_dir=out_file_dir,
                out_files=list(
                  circ_bed=paste0(out_file_dir,"/",names(bam_list),"/circlemap_reports/",
                  "repeat_reports/",names(bam_list),".circular_candidates.bed"),
                  circ_repeat_bed=paste0(out_file_dir,"/",names(bam_list),"/circlemap_reports/",
                  "realign_reports/",names(bam_list),".circular_repeat_candidates.bed"))
                )
        
    }

    return(jobs_report)

}





#' Read Circle-Map circDNA BED file output
#'
#' This function reads the circle-map circDNA output file
#' 
#' 
#' //TODO validate co-joint calling mode
#' 
#' For more information read:
#' https://github.com/iprada/Circle-Map
#' 
#' 
#' 
#' @param bed [REQUIRED] Path to bed file.
#' @param type [REQUIRED] circleMap output type. Default ["repeat"].
#' @param id [OPTIONAL] Sample identifier. If not given the name of the first tumour sample of the samples will be used.
#' @export


read_bed_circlemap=function(
    bed=NULL,type="repeat",
    id=NULL,
    sep="\t"
){
   
    if(is.null(id)){
        id=get_file_name(bed)
    }

    if(type=="repeat"){

        dat=read.table(file=bed,sep=sep)
        names(dat)=c("chr","start","end","n_repeats",
        "circ_score","mean_cov","sd_cov","cov_incr_start",
        "cov_incr_end","cov_continuity")
        ### There is an empty column in the BED file generated by circlemap
        dat$`NA`=NULL
        ### Redundant column for repeat data
        dat$circ_score=NULL
    }else if (type=="align"){

        dat=read.table(file=bed,sep=sep)
        names(dat)=c("chr","start","end","disc_reads","split_reads",
        "circ_score","mean_cov","sd_cov","cov_incr_start",
        "cov_incr_end","cov_continuity")
        ### There is an empty column in the BED file generated by circlemap
        dat$`NA`=NULL
       
    }

    dat$id=id

    return(dat)

}


#' Annotate Circle-Map circDNA BED file output
#'
#' This function reads the circle-map circDNA output file
#' 
#' 
#' //TODO validate co-joint calling mode
#' 
#' For more information read:
#' https://github.com/iprada/Circle-Map
#' 
#' 
#' 
#' @param bed [REQUIRED] Path to bed file.
#' @param annotation_ref [REQUIRED] Path to annotation file.
#' @param type [REQUIRED] circleMap output type. Default ["repeat"].
#' @param id [OPTIONAL] Sample identifier. If not given the name of the first tumour sample of the samples will be used.
#' @export

annotate_bed_circlemap=function(
    bed=NULL,
    id=NULL,
    type="repeat",
    annotation_ref=build_default_reference_list()$HG19$panel$PCF_V3$annotation$genes,
    sep="\t",
    threads=8
){
    dat=read_bed_circlemap(bed=bed,id=id,type=type,sep=sep)
    annotation=read.table(annotation_ref,sep="\t",header=TRUE) %>% 
    dplyr::select(chr,start,end,gene_id)


    summarised_dat=parallel::mclapply(unique(dat$chr),function(chrom){
        tmp_dat=dat %>% dplyr::filter(chr==chrom)
        tmp_annotation=annotation %>% dplyr::filter(chr==chrom)

        full_gene=fuzzyjoin::fuzzy_left_join(tmp_dat,tmp_annotation,
            by=c("chr"="chr","start"="start","end"="end"),
            match_fun=c(`==`,`<=`,`>=`)
            )

        full_gene$annot_type="COMPLETE"

        partial_left=fuzzyjoin::fuzzy_left_join(tmp_dat,tmp_annotation,
            by=c("chr"="chr","start"="start","end"="start"),
            match_fun=c(`==`,`<=`,`>=`)
        )

        partial_left$annot_type="PARTIAL"
        partial_left=dplyr::anti_join(partial_left,full_gene,
            by=c("chr.x"="chr.x",
                "start.x"="start.x",
                "end.x"="end.x",
                "gene_id"="gene_id"
            ))

        partial_right=fuzzyjoin::fuzzy_left_join(tmp_dat,tmp_annotation,
            by=c("chr"="chr","start"="end","end"="end"),
            match_fun=c(`==`,`<=`,`>=`)
        )

        partial_right$annot_type="PARTIAL"
        partial_right=dplyr::anti_join(partial_right,full_gene,
            by=c("chr.x"="chr.x",
                "start.x"="start.x",
                "end.x"="end.x",
                "gene_id"="gene_id"
            ))


        complete_dat=rbind(full_gene,partial_left,partial_right)
        summarised_dat=complete_dat %>% 
            dplyr::rename(chr=chr.x,start=start.x,end=end.x) %>%
            dplyr::select(chr:gene_id,annot_type)  %>% 
            dplyr::distinct(chr:annot_type) %>%
            dplyr::mutate(name=paste0(gene_id,":",annot_type)) %>%
            dplyr::group_by(dplyr::across(chr:id)) %>% 
            dplyr::summarise(genes=paste0(name,collapse=";")
        )
        return(summarised_dat)
    },mc.cores=threads)

    summarised_dat=dplyr::bind_rows(summarised_dat) %>% 
    dplyr::arrange(gtools::mixedsort(chr),start)

    return(summarised_dat)
}

