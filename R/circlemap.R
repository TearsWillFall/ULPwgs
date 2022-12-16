
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
        inherit=NULL,
        select=NULL,
        env_circlemap=build_default_python_enviroment_list()$env_circlemap,
        bin_samtools=build_default_tool_binary_list()$bin_samtools,
        bam=NULL,
        output_dir=".",
        output_name=NULL,
        verbose=FALSE,
        ref_genome=build_default_reference_list()$HG19$reference$genome,
        batch_config=build_default_preprocess_config(),
        threads=3,
        ram=1,
        ns="ULPwgs",
        mode="local",
        tmp_dir=NULL,
        executor_id=make_unique_id("realignCircleMap"),
        task_name="realignCircleMap",
        time="48:0:0",
        update_time=60,
        wait=FALSE,
        hold=NULL,
        err_mssg="realign failed to run"
){

    this.envir=environment()
    set_envir_vars(envir=this.envir,inputs=bam,ids=output_name,dir_name="realign_reports")


    out_file=paste0(out_file_dir,"/",id,".circular_candidates.bed")

    main_realign_circlemap=function(envir){



        this.envir=environment()
        append_envir(this.envir,envir)

        steps=list()
        steps$realign$job_id=job_id
         
        steps$realign <- append(steps$realign,new_sort_and_index_bam_samtools(
            bin_samtools=bin_samtools,
            bam=input,
            output_dir=out_file_dir_tmp,
            verbose=verbose,
            batch_config=batch_config,
            threads=threads,
            ram=ram,
            coord_sort=FALSE,
            index=FALSE,
            clean=FALSE,
            executor_id=task_id
        ))

        steps$realign <- append(steps$realign,read_extractor_circlemap(
                env_circlemap=env_circlemap,
                bin_samtools=bin_samtools,
                bam=steps$sort_and_index$sort$out_file,
                output_dir=out_file_dir_tmp,
                verbose=verbose,
                batch_config=batch_config,
                threads=threads,ram=ram,
                sort=TRUE,
                coord_sort=TRUE,
                index=TRUE,stats="all", 
                clean=TRUE,
                executor_id=task_id
            )
        )

        steps$realign$exec_code <- paste(
            set_conda_enviroment(env_circlemap),
            " Circle-Map Realign -sbam ",normalizePath(input),
            " -qbam ", steps$sort_and_index$out_file,
            " -i ",steps$read_extractor$sort_and_index$sort$out_file,
            " -o ",steps$realign$out_file,
            " -t ",threads," -dir /", 
            " -fasta ",normalizePath(ref_genome), 
            " -tdir ",out_file_dir_tmp
        )

        envir$steps <- steps

    }

    if(is.null(select)){
        run_job(
            envir=this.envir
        )
    }else{

       set_envir_inputs(envir=this.envir)
       main_repeat_caller_circlemap(
            envir=this.envir
       )
       return(steps)
    }

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
    inherit=NULL,
    select=NULL,
    env_circlemap=build_default_python_enviroment_list()$env_circlemap,
    bin_samtools=build_default_tool_binary_list()$bin_samtools,
    bam=NULL,
    output_dir=".",
    output_name=NULL,
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=3,
    ram=1,
    sort=TRUE,
    coord_sort=FALSE,
    index=FALSE,
    stats="all",
    clean=FALSE,
    ns="ULPwgs",
    mode="local",
    executor_id=make_unique_id("readExtractCircleMap"),
    task_name="readExtractCircleMap",
    time="48:0:0",
    update_time=60,
    wait=FALSE,
    hold=NULL,
    err_mssg="read_extractor"

  ){



    this.envir=environment()
    set_envir_vars(envir=this.envir,inputs=bam,ids=output_name,dir_name="read_extractor")



    main_read_extractor_circlemap=function(
        envir
    ){
    
        this.envir=environment()
        append_envir(this.envir,envir)

        steps=list()
        steps$read_extractor$job_id=job_id
    
        steps$read_extractor$out_file=paste0(
            out_file_dir,"/",input_id,".circular_read_candidates.bam"
        )

        steps$read_extractor$exec_code=paste(
            set_conda_enviroment(env_circlemap),
            "Circle-Map ReadExtractor -i ",normalizePath(input), 
            " -o ", steps$read_extractor$out_file," -dir /"
        )

        if(verbose){
            print_verbose(
                job=job_id,arg=as.list(this.envir),
                exec_code=steps$read_extractor$exec_code)
        }
        
        steps$read_extractor$error=execute_job(
            exec_code=steps$read_extractor$exec_code
        )
    
        if(steps$read_extractor$error!=0){
            stop(err_mssg)
        }


        steps$read_extractor<-append(
            steps$read_extractor,
                sort_and_index_bam_samtools(
                bin_samtools=bin_samtools,
                bam=steps$read_extractor$out_file,
                output_dir=out_file_dir,
                verbose=verbose,
                batch_config=batch_config,
                threads=threads,
                ram=ram,
                sort=TRUE,
                coord_sort=FALSE,
                index=TRUE,
                clean=FALSE,
                executor_id=task_id
            )
        )

        envir$steps <- steps
        
    }

    if(is.null(select)){
        run_job(
            envir=this.envir
        )
    }else{
       set_envir_inputs(envir=this.envir)
       main_read_extractor_circlemap(
            envir=this.envir
       )
       return(steps)
    }


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
    inherit=NULL,
    select=NULL,
    env_circlemap=build_default_python_enviroment_list()$env_circlemap,
    bam=NULL,
    output_dir=".",
    output_name=NULL,
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,
    ram=1,
    ns="ULPwgs",
    mode="local",
    executor_id=make_unique_id("RepeatsCircleMap"),
    task_name="RepeatsCircleMap",
    time="48:0:0",
    update_time=60,
    wait=FALSE,
    hold=NULL,
    err_mssg="repeat_caller failed to run"
){

  
    this.envir=environment()
    set_envir_vars(envir=this.envir,inputs=bam,ids=output_name,dir_name="repeat_reports")



    main_repeat_caller_circlemap=function(
        envir
    ){
            this.envir=environment()
            append_envir(this.envir,envir)

            steps=list()
            steps$repeat_caller$job_id<-job_id

            steps$repeat_caller$out_file=paste0(
                out_file_dir,"/",input_id,".circular_repeat_candidates.bed"
            )
            steps$repeat_caller$exec_code=paste(
                set_conda_enviroment(env_circlemap),
                " Circle-Map Repeats -i ",normalizePath(input), " -o ",
                steps$repeat_caller$out_file, " -dir /"
            )

            if(verbose){
                print_verbose(job=job_id,arg=as.list(this.envir),
                exec_code=steps$repeat_caller$exec_code
                )
            }
    
            steps$repeat_caller$error=execute_job(exec_code=steps$repeat_caller$exec_code)
        
            if(steps$repeat_caller$error!=0){
                stop(err_mssg)
            }

            envir$steps <- steps
    }

    if(is.null(select)){
        run_job(
            envir=this.envir
        )
    }else{

       set_envir_inputs(envir=this.envir)
       main_repeat_caller_circlemap(
            envir=this.envir
        )
       return(steps)
    }



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
        inherit=NULL,
        select=NULL,
        env_circlemap=build_default_python_enviroment_list()$env_circlemap,
        bin_samtools=build_default_tool_binary_list()$bin_samtools,
        ref_genome=build_default_reference_list()$HG19$reference$genome,
        bam=NULL,
        output_dir=".",
        output_name=NULL,
        tmp_dir=NULL,
        verbose=FALSE,
        batch_config=build_default_preprocess_config(),
        threads=3,
        ram=1,
        mode="local",
        ns="ULPwgs",
        executor_id=make_unique_id("circdnaCircleMap"),
        task_name="circdnaCircleMap",
        time="48:0:0",
        update_time=60,
        wait=FALSE,
        hold=NULL,
        err_mssg="circdna failed to run"
){

 
    this.envir=environment()
    set_envir_vars(
        envir=this.envir,
        inputs=bam,
        ids=output_name,
        dir_name="repeat_reports"
    )


    main_circdna_circlemap=function(
        envir
    ){

        this.envir=environment()
        append_envir(this.envir,envir)

        steps=list()
        steps$circdna$job_id=job_id
        steps$circdna <- append(steps$circdna,realign_circlemap(
            env_circlemap=env_circlemap,
            bin_samtools=bin_samtools,
            bam=normalizePath(input),
            ref_genome=normalizePath(ref_genome),
            output_dir=out_file_dir,
            verbose=verbose,
            tmp_dir=out_file_dir_tmp,
            batch_config=batch_config,
            threads=threads,
            ram=ram,
            executor_id=task_id,
        ))


        steps$circdna <- append(steps$circdna,repeat_caller_circlemap(
            env_circlemap=env_circlemap,
            bam=normalizePath(input),
            output_dir=out_file_dir,
            verbose=verbose,
            batch_config=batch_config,
            threads=threads,ram=ram,
            executor_id=task_id
        ))

        envir$steps<-steps
        
    }

   
    if(is.null(select)){
        run_job(
            envir=this.envir
        )
    }else{

       set_envir_inputs(envir=this.envir)
       main_circdna_circlemap(
            envir=this.envir
        )
       return(steps)
    }

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
#' @param output_dir Path to the output directory.
#' @param threads [OPTIONAL] Number of threads to split the work. Default 4
#' @param ram [OPTIONAL] RAM memory to asing to each thread. Default 4
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Task EXECUTOR ID. Default "recalCovariates"
#' @param task_name Task name. Default "recalCovariates"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export

annotate_bed_circlemap=function(
    inherit=NULL,
    select=NULL,
    bed=NULL,
    id=NULL,
    output_name=NULL,
    type="repeat",
    write=TRUE,
    annotation_ref=build_default_reference_list()$HG19$panel$PCF_V3$annotation$genes,
    ns="ULPwgs",
    output_dir=".",
    mode="local",
    time="48:0:0",
    threads=threads,
    ram=1,
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    executor_id=make_unique_id("annotateBEDcirclemap"),
    task_name="annotateBEDcirclemap",
    hold=NULL,
    err_mssg="annotate_bed failed to run"
){
    
    this.envir=environment()
    set_envir_vars(envir=this.envir,inputs=bed,ids=output_name)

    main_annotate_bed_circlemap=function(
        envir
    ){
        this.envir=environment()
        append_envir(this.envir,envir)


        steps=list()
        steps$annotate_bed$job_id <- job_id

        dat=read_bed_circlemap(bed=input,id=id,type=type,sep="\t")
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
                dplyr::mutate(annot_name=paste0(gene_id,":",annot_type)) %>%
                dplyr::select(chr:id,annot_name) %>%
                dplyr::distinct()%>%
                dplyr::group_by(dplyr::across(c(-annot_name))) %>%
                dplyr::summarise(genes=paste0(annot_name,collapse=";")
            )
            return(summarised_dat)
        },mc.cores=threads)

        summarised_dat=dplyr::bind_rows(summarised_dat) %>% 
        dplyr::arrange(gtools::mixedsort(chr))

        if(write){
            steps$annotate_bed$out_file=paste0(out_file_dir,"/",id,".circular_",type,".annotated.bed")
            write.table(file=steps$annotate_bed$out_file,summarised_dat,sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
        }

        envir$steps<-steps
       
    }


    if(is.null(select)){
        run_job(
            envir=this.envir
        )
    }else{
        
        set_envir_inputs(envir=this.envir)
        main_annotate_bed_circlemap(
           envir=this.envir
        )
        return(steps)
    }


  

}

