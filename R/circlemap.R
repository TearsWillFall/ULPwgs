
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
        ref_genome=build_default_reference_list()$HG19$reference$genome,
        bam=NULL,
        output_dir=".",
        output_name=NULL,
        verbose=FALSE,
        batch_config=build_default_preprocess_config(),
        threads=3,
        ram=1,
        ns="ULPwgs",
        mode="local",
        tmp_dir=NULL,
        time="48:0:0",
        err_mssg="realign failed to run",
        inherit=NULL,
        select=NULL,
        executor_id=NULL,
        hold=NULL

){

    this.envir=environment()
    set_envir_vars(
        envir=this.envir,
        vars="bam",
        executor_id = executor_id,
        dir_name="realign_reports"
    )


    run_main=function(envir){



        this.envir=environment()
        append_envir(this.envir,envir)
        set_steps_vars(envir=this.envir)
    
         
        steps[[fn]] <- append(steps[[fn]],new_sort_and_index_bam_samtools(
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

        steps[[fn]] <- append(steps[[fn]],read_extractor_circlemap(
                env_circlemap=env_circlemap,
                bin_samtools=bin_samtools,
                bam=steps$new_sort_and_index_bam_samtools$out_file,
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


        steps[[fn]]$out_file=paste0(out_file_dir,"/",input_id,".circular_candidates.bed")

        steps[[fn]]$exec_code <- paste(
            set_conda_envir(env_circlemap),
            " Circle-Map Realign -sbam ",normalizePath(input),
            " -qbam ", steps[[fn]]$new_sort_and_index_bam_samtools$sort$out_file,
            " -i ",steps[[fn]]$new_sort_and_index_bam_samtools$sort$out_file,
            " -o ",steps[[fn]]$out_file,
            " -t ",threads," -dir /", 
            " -fasta ",normalizePath(ref_genome), 
            " -tdir ",out_file_dir_tmp
        )
        

        run_job(envir=this.envir)

        envir$steps <- steps

    }

    run_envir(envir=this.envir)
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
    time="48:0:0",
    err_mssg="read_extractor",
    inherit=NULL,
    select=NULL,
    executor_id=NULL,
    hold=NULL

  ){



    this.envir=environment()
    set_envir_vars(
        envir=this.envir,
        vars="bam",
        executor_id = executor_id,
        dir_name="read_extractor"
    )



    run_main=function(
        envir
    ){
    
        this.envir=environment()
        append_envir(this.envir,envir)
        set_steps_vars(envir=this.envir)
    
        steps[[fn]]$out_file=paste0(
            out_file_dir,"/",input_id,".circular_read_candidates.bam"
        )

        steps[[fn]]$exec_code=paste(
            set_conda_envir(env_circlemap),
            "Circle-Map ReadExtractor -i ",normalizePath(input), 
            " -o ", steps[[fn]]$out_file," -dir /"
        )

        run_job(envir=this.envir)

        steps[[fn]]<-append(
            steps[[fn]],
                new_sort_and_index_bam_samtools(
                    bin_samtools=bin_samtools,
                    bam=steps[[fn]]$out_file,
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



    run_envir(envir=this.envir)

 
    


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
    bam=NULL,
    output_dir=".",
    output_name=NULL,
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,
    ram=1,
    ns="ULPwgs",
    mode="local",
    time="48:0:0",
    err_mssg="repeat_caller failed to run",
    inherit=NULL,
    select=NULL,
    executor_id=NULL,
    hold=NULL
){

  
    this.envir=environment()
    set_envir_vars(
        envir=this.envir,
        vars="bam",
        executor_id = executor_id,
        dir_name="repeat_reports"
    )



    run_main=function(
        envir
    ){
            this.envir=environment()
            append_envir(this.envir,envir)
            set_steps_vars(envir=this.envir)
          
            steps[[fn]]$out_file=paste0(
                out_file_dir,"/",input_id,".circular_repeat_candidates.bed"
            )
            steps[[fn]]$exec_code=paste(
                set_conda_envir(env_circlemap),
                " Circle-Map Repeats -i ",normalizePath(input), " -o ",
                steps[[fn]]$out_file, " -dir /"
            )

            run_job(envir=this.envir)

            envir$steps <- steps
    }

    run_envir(envir=this.envir)



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
        bam=NULL,
        output_dir=".",
        output_name=NULL,
        verbose=FALSE,
        batch_config=build_default_preprocess_config(),
        threads=3,
        ram=1,
        mode="local",
        ns="ULPwgs",
        time="48:0:0",
        err_mssg="circdna failed to run",
        inherit=NULL,
        select=NULL,
        executor_id=NULL,
        hold=NULL
){

 
    this.envir=environment()
    set_envir_vars(
        envir=this.envir,
        vars="bam",
        executor_id = executor_id,
        dir_name="repeat_reports"
    )


    run_main=function(
        envir
    ){

        this.envir=environment()
        append_envir(this.envir,envir)
        set_steps_vars(envir=this.envir)

  
        steps[[fn]] <- append(steps[[fn]],realign_circlemap(
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


        steps[[fn]] <- append(steps[[fn]],repeat_caller_circlemap(
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

    run_envir(envir=this.envir)

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
    annotation_ref=build_default_reference_list()$HG19$panel$PCF_V3$annotation$genes,
    bed=NULL,
    output_dir=".",
    output_name=NULL,
    type="repeat",
    write=TRUE,
    ns="ULPwgs",
    mode="local",
    time="48:0:0",
    threads=8,
    ram=1,
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    err_mssg="annotate_bed failed to run",
    inherit=NULL,
    select=NULL,
    executor_id=NULL,
    hold=NULL
){
    
    this.envir=environment()
    set_envir_vars(
        envir=this.envir,
        vars="bed",
        executor_id = executor_id
    )

    run_main=function(
        envir
    ){
        this.envir=environment()
        append_envir(this.envir,envir)
        set_steps_vars(envir=this.envir)


        dat=read_bed_circlemap(bed=input,id=output_name,type=type,sep="\t")
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
            steps[[fn]]$out_file=paste0(out_file_dir,"/",input_id,".circular_",type,".annotated.bed")
            write.table(
                file=steps[[fn]]$out_file,
                summarised_dat,
                sep="\t",
                col.names=TRUE,
                row.names=FALSE,
                quote=FALSE
            )
        }

        envir$steps<-steps
       
    }


    run_envir(envir=this.envir)
  

}

