
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
        ...

){

    run_main=function(
        .env
    ){

        .this.env=environment()
        append_env(to=.this.env,from=.env)

        set_main(.env=.this.env)


        .main$out_file=paste0(out_file_dir,"/",input_id,".realign.circ_candidates.bed")
        .main$out_files$realign_bed <- .main$out_file

            
        .main$steps <- append(.main$steps,
            new_sort_and_index_bam_samtools(
                bin_samtools=bin_samtools,
                bam=input,
                output_dir=paste0(out_file_dir,"/sorted"),
                tmp_dir=tmp_dir,
                env_dir=env_dir,
                batch_dir=batch_dir,
                verbose=verbose,
                threads=threads,
                ram=ram,
                err_msg=err_msg,
                coord_sort=FALSE,
                stats=FALSE,
                index=FALSE,
                clean=FALSE,
                executor_id=task_id
            )
        )

        .this.step=.main$steps$new_sort_and_index_bam_samtools
        .main$out_files$srt_qbam=.this.step$out_files$srt_bam
       
        .main$steps <- append(.main$steps,
            read_extractor_circlemap(
                env_circlemap=env_circlemap,
                bin_samtools=bin_samtools,
                bam=.main$out_files$srt_qbam,
                output_dir=paste0(out_file_dir,"/read_extractor"),
                tmp_dir=tmp_dir,
                env_dir=env_dir,
                batch_dir=batch_dir,
                verbose=verbose,
                threads=threads,
                ram=ram,
                err_msg=err_msg,
                sort=TRUE,
                coord_sort=TRUE,
                index=TRUE,
                stats="all",
                executor_id=task_id
            )
        )
        
        .this.step=.main$steps$read_extractor_circlemap
        .main$out_files=append(.main$out_files,.this.step$out_files)
        

        .main$exec_code <- paste(
            set_conda_envir(env_circlemap),
            " Circle-Map Realign -sbam ",normalizePath(input),
            " -qbam ", .main$out_files$srt_qbam,
            " -i ",.main$out_files$srt_cbam,
            " -o ",.main$out_files$realign_bed,
            " -t ",threads," -dir /", 
            " -fasta ",normalizePath(ref_genome), 
            " -tdir ",out_file_dir
        )

        run_job(.env=.this.env)

        .env$.main <- .main

    }

    .base.env=environment()
    list2env(list(...),env=.base.env)
    set_env_vars(
        .env= .base.env,
        vars="bam"
    )

    launch(.env=.base.env)

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
    ...
  ){


    

    run_main=function(
       .env
    ){
    
        .this.env=environment()
        append_env(to=.this.env,from=.env)

        set_main(.env=.this.env)
    
        .main$out_file=paste0(
            out_file_dir,"/",input_id,".circular_read_candidates.bam"
        )

       .main$exec_code=paste(
            set_conda_envir(env_circlemap),
            "Circle-Map ReadExtractor -i ",normalizePath(input), 
            " -o ", .main$out_file," -dir /"
        )

        run_job(.env=.this.env)

        .main.step=.main$steps[[fn]]
        .main.step$out_files$unsrt_cbam=.main.step$out_file

        .main$steps[[fn]]$steps<-append(
           .main$steps[[fn]]$steps,
                new_sort_and_index_bam_samtools(
                    bin_samtools=bin_samtools,
                    bam=.main$out_file,
                    output_dir=out_file_dir,
                    tmp_dir=tmp_dir,
                    env_dir=env_dir,
                    batch_dir=batch_dir,
                    verbose=verbose,
                    threads=threads,
                    ram=ram,
                    err_msg=err_msg,
                    stats=FALSE,
                    sort=TRUE,
                    index=TRUE,
                    executor_id=task_id
            )
        )

        .this.step=.main.step$steps$new_sort_and_index_bam_samtools
        .main.step$out_files=append(.main.step$out_files,.this.step$out_files)
        names(.main.step$out_files)[-1]<-c("srt_cbam","index_cbam")
         
        .env$.main <- .main   
    }



       
    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="bam"
    )

    launch(.env=.base.env)

 
    


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
    ...
){

  
    run_main=function(
        .env
    ){
        
            .this.env=environment()
            append_env(to=.this.env,from=.env)

            set_main(.env=.this.env)
    
          
            .main$out_files$realign_bed=paste0(
                out_file_dir,"/",input_id,".repeat.circ_candidates.bed"
            )
            .main$exec_code=paste(
                set_conda_envir(env_circlemap),
                " Circle-Map Repeats -i ",normalizePath(input), " -o ",
                .main$out_file, " -dir /"
            )

            run_job(.env=.this.env)

            .env$.main <- .main   
    }

         
    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="bam"
    )

    launch(.env=.base.env)




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
        annotation_ref=build_default_reference_list()$HG19$panel$PCF_V3$annotation$genes,
        env_circlemap=build_default_python_enviroment_list()$env_circlemap,
        bin_samtools=build_default_tool_binary_list()$bin_samtools,
        ref_genome=build_default_reference_list()$HG19$reference$genome,
        bam=NULL,
        ...
){

    run_main=function(
       .env
    ){

        .this.env=environment()
        append_env(to=.this.env,from=.env)

        set_main(.env=.this.env)

        .main$steps[[fn]]<-.this.env
        
        .main.step=.main$steps[[fn]]
    
  
        .main$steps[[fn]]$steps <- append(
            .main$steps[[fn]]$steps,realign_circlemap(
                env_circlemap=env_circlemap,
                bin_samtools=bin_samtools,
                bam=normalizePath(input),
                ref_genome=normalizePath(ref_genome),
                output_dir=paste0(out_file_dir,"/realign_reports"),
                tmp_dir=tmp_dir,
                env_dir=env_dir,
                batch_dir=batch_dir,
                verbose=verbose,
                threads=threads,
                err_msg=err_msg,
                ram=ram,
                executor_id=task_id
            )
        )
        

        .this.step=.main.step$steps$realign_circlemap
        .main.step$out_files=append(.main.step$out_files,.this.step$out_files)



        .main$steps[[fn]]$steps  <- append(
            .main$steps[[fn]]$steps,
            repeat_caller_circlemap(
                env_circlemap=env_circlemap,
                bam=normalizePath(input),
                output_dir=paste0(out_file_dir,"/repeat_reports"),
                tmp_dir=tmp_dir,
                env_dir=env_dir,
                batch_dir=batch_dir,
                verbose=verbose,
                batch_config=batch_config,
                threads=threads,
                ram=ram,
                err_msg=err_msg,
                executor_id=task_id
            )
        )

        .this.step=.main.step$steps$repeat_caller_circlemap
        .main.step$out_files=append(.main.step$out_files,.this.step$out_files)


        .env$.main <- .main
        
    }

    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="bam"
    )

    launch(.env=.base.env)


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
    type="repeat",
    ...
){
  
  
    run_main=function(
        .env
    ){
   
        .this.env=environment()
        append_env(to=.this.env,from=.env)
        set_main(.env=.this.env)

        .main$steps[[fn]]<-.this.env
  
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

        
        .main$out_files$annot_circ_bed=paste0(
            out_file_dir,"/",input_id,".",type,
            ".circ_candidates.annotated.bed")


        
        write.table(
            file=.main$out_files$annot_circ_bed,
            summarised_dat,
            sep="\t",
            col.names=TRUE,
            row.names=FALSE,
            quote=FALSE
        )
   
        .env$.main <-.main

    }


    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
        .env=.base.env,
        vars="bed"
    )

    launch(.env=.base.env)


}

