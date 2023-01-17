

#' Strelka wrapper for SNV variant calling
#'
#' This function wraps the STRELKA functions for variant calling
#' 
#' @param bin_strelka_somatic Path to strelka somatic workflow binary
#' @param bin_strelka_germline Path to strelka germline workflow binary
#' @param tumour [OPTIONAL] Path to tumour BAM file. If not given will assume germline variant calling.
#' @param normal [REQUIRED] Path to tumour BAM file.
#' @param ref_genome [REQUIRED] Path to reference genome FASTA
#' @param variants [REQUIRED] Variants types to call. Default all. Options ["snv","sv","all"]
#' @param indel_candidates [OPTIONAL] Path to indel candidates file produced by MANTA.
#' @param targeted [REQUIRED] Remove coverage filtering for exome/targeted data. Default TRUE
#' @param output_dir [OPTIONAL] Path to the output directory. Default current directory
#' @param threads [OPTIONAL] Number of threads to split the work. Default 4
#' @param batch_config [OPTIONAL] Default configuration for job submission in batch.
#' @param ram [OPTIONAL] RAM memory to asing to each thread. Default 4
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Task EXECUTOR ID. Default "recalCovariates"
#' @param task_name Task name. Default "recalCovariates"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export


call_sv_manta=function( 
    bin_bcftools=build_default_tool_binary_list()$bin_bcftools,
    bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
    bin_tabix=build_default_tool_binary_list()$bin_tabix,
    bin_vep=build_default_tool_binary_list()$bin_vep,
    bin_manta=build_default_tool_binary_list()$bin_manta,
    cache_vep=build_default_cache_list()$cache_vep,
    tumour=NULL,
    normal=NULL,
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    targeted=TRUE,
    annotate=TRUE,
    tabulate=TRUE,
    ...
){


    run_main=function(
        .env
    ){
      .this.env=environment()
      append_env(to=.this.env,from=.env)
      set_main(.env=.this.env)
      

      if(!is.null(normal)){

        if(!is.null(tumour)){

            .main.step$steps<-append(
              .main.step$steps,
              call_somatic_sv_manta(
                bin_manta=bin_manta,
                bin_bcftools=bin_bcftools,
                bin_bgzip=bin_bgzip,
                bin_tabix=bin_tabix,
                bin_vep=bin_vep,
                cache_vep=cache_vep,
                tumour=input,
                normal=normal,
                annotate=annotate,
                tabulate=tabulate,
                ref_genome=ref_genome,
                output_dir=paste0(out_file_dir,"/somatic/sv"),
                tmp_dir=tmp_dir,
                env_dir=env_dir,
                batch_dir=batch_dir,
                err_msg=err_msg,
                targeted=targeted,
                verbose=verbose,
                threads=threads,
                ram=ram,
                executor_id=task_id

              )
            )

          .this.step=.main.step$steps$call_somatic_sv_manta
          .main.step$out_files$somatic=.this.step$out_files
        }

      }else{
           .main.step$steps<-append(
              .main.step$steps,
              call_germline_sv_manta(
                bin_manta=bin_manta,
                bin_bcftools=bin_bcftools,
                bin_bgzip=bin_bgzip,
                bin_tabix=bin_tabix,
                bin_vep=bin_vep,
                cache_vep=cache_vep,
                normal=normal,
                annotate=annotate,
                tabulate=tabulate,
                ref_genome=ref_genome,
                output_dir=paste0(out_file_dir,"/germline/sv"),
                tmp_dir=tmp_dir,
                env_dir=env_dir,
                batch_dir=batch_dir,
                err_msg=err_msg,
                targeted=targeted,
                verbose=verbose,
                threads=threads,
                ram=ram,
                executor_id=task_id
            )
           )
          .this.step=.main.step$steps$call_somatic_sv_manta
          .main.step$out_files$germline=.this.step$out_files
      }

      .env$.main <- .main

    }
        
    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="tumour"
    )

    launch(.env=.base.env)

    
}




#' Manta  wrapper for structural variant calling
#'
#' This function wraps the Manta workflow functions for structural variant calling 
#' 
#' @param bin_manta Path to manta pipeline binary
#' @param tumour [OPTIONAL] Path to tumour BAM file. If not given will assume germline variant calling.
#' @param normal [REQUIRED] Path to tumour BAM file.
#' @param ref_genome [REQUIRED] Path to reference genome FASTA
#' @param targeted [REQUIRED] Remove coverage filtering for exome/targeted data. Default TRUE
#' @param output_dir [OPTIONAL] Path to the output directory. Default current directory
#' @param threads [OPTIONAL] Number of threads to split the work. Default 4
#' @param ram [OPTIONAL] RAM memory to asing to each thread. Default 4
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param batch_config [OPTIONAL] Default configuration for job submission in batch.
#' @param executor_id Task EXECUTOR ID. Default "recalCovariates"
#' @param task_name Task name. Default "recalCovariates"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export



call_germline_sv_manta=function(
    bin_bcftools=build_default_tool_binary_list()$bin_bcftools,
    bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
    bin_tabix=build_default_tool_binary_list()$bin_tabix,
    bin_vep=build_default_tool_binary_list()$bin_vep,
    cache_vep=build_default_cache_list()$cache_vep,
    bin_manta=build_default_tool_binary_list()$bin_manta,
    normal=NULL,
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    targeted=TRUE,
    annotate=TRUE,
    tabulate=TRUE,
    ...
){

    run_main=function(
        .env
    ){
          .this.env=environment()
          append_env(to=.this.env,from=.env)
      
          set_main(.env=.this.env)


         .main$exec_code=paste0(
          bin_manta,
          " --bam ",normalizePath(input),
          " --referenceFasta ", normalizePath(ref_genome) ,
          " --runDir ", out_file_dir, 
          ifelse(targeted," --exome ",""), "; ",
          paste0(
            out_file_dir,
            "/runWorkflow.py -m local -j ",
            threads)
          )

        .main$out_files$strelka$workflow=paste0(out_file_dir,"/runWorkflow.py")
        .main$out_files$strelka$stats=list(
              tsv=paste0(out_file_dir,"/results/stats/runStats.tsv"),
              xml=paste0(out_file_dir,"/results/stats/runStats.xml")
        )
        .main$out_files$strelka$variants=list(
              indel_candidate=paste0(out_file_dir,"/results/variants/candidateSmallIndels.vcf.gz"),
              indel_candidate_idx=paste0(out_file_dir,"/results/variants/candidateSmallIndels.vcf.gz.tbi"),
              diploid_sv=paste0(out_file_dir,"/results/variants/diploidSV.vcf.gz"),
              diploid_sv_idx=paste0(out_file_dir,"/results/variants/diploidSV.vcf.gz.tbi")
        )

        run_job(
            .env=.this.env
        )

        .main.step=.main.step

    
        .main.step$steps <- append(
          .main.step$steps, 
          add_af_strelka_vcf(
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            vcf= .main.step$out_files$variants$diploid_sv,
            output_dir=paste0(out_file_dir,"/annotated"),
            tmp_dir=tmp_dir,
            env_dir=env_dir,
            batch_dir=batch_dir,
            err_msg=err_msg,
            verbose=verbose, 
            executor_id=task_id,
            threads=threads,
            ram=ram
          )
        )
        .this.step=.main.step$steps$add_af_strelka_vcf
        .main.step$out_files$annotated$af=.this.step$out_files
      
        .main.step$steps <-append(
          .main.step$steps, 
          extract_pass_variants_strelka_vcf(
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            vcf=.main.step$out_files$annotated$af$sv$bgzip_vcf,
            type="sv",
            output_dir=paste0(out_file_dir,"/annotated"),
            tmp_dir=tmp_dir,
            env_dir=env_dir,
            batch_dir=batch_dir,
            err_msg=err_msg,
            verbose=verbose,
            executor_id=task_id,
            threads=threads,
            ram=ram
    
          )
        )

        .this.step=.main.step$steps$extract_pass_variants_strelka_vcf.sv
        .main.step$out_files$annotated$filter=.this.step$out_files
      
        if(annotate){
            .main.step$steps<-append(
              .main.step$steps,
              annotate_strelka_vep(
                bin_vep=bin_vep,
                bin_bgzip=bin_bgzip,
                bin_tabix=bin_tabix,
                cache_vep=cache_vep,
                vcf=.main.step$out_files$annotated$filter$sv$bgzip_vcf,
                type="sv",
                output_dir=paste0(out_file_dir,"/annotated"),
                tmp_dir=tmp_dir,
                env_dir=env_dir,
                batch_dir=batch_dir,
                err_msg=err_msg,
                verbose=verbose,
                executor_id=task_id,
                threads=threads,
                ram=ram
            )
          )

        .this.step=.main.step$steps$annotate_strelka_vep.sv
        .main.step$out_files$annotated$vep=.this.step$out_files
        
          if(tabulate){
              .main.step$steps<-append(
                .main.step$steps,
                tabulate_strelka_vcf(
                  vcf=.main.step$out_files$annotated$vep$sv$bgzip_vcf,
                  type="sv",
                  output_dir=paste0(out_file_dir,"/annotated"),
                  tmp_dir=tmp_dir,
                  env_dir=env_dir,
                  batch_dir=batch_dir,
                  err_msg=err_msg,
                  threads=threads,
                  ram=ram,
                  verbose=verbose,
                  executor_id=task_id
              )
            )

            .this.step=.main.step$steps$tabulate_strelka_vcf.sv
            .main.step$out_files$annotated$tabulate=.this.step$out_files
          }
        }

        .env$.main <- .main
    }


  .base.env=environment()
  list2env(list(...),envir=.base.env)
  set_env_vars(
    .env= .base.env,
    vars="normal"
  )

  launch(.env=.base.env)


}




#' Manta  wrapper for structural variant calling
#'
#' This function wraps the Manta workflow functions for structural variant calling 
#' 
#' @param bin_manta Path to manta pipeline binary
#' @param tumour [OPTIONAL] Path to tumour BAM file. If not given will assume germline variant calling.
#' @param normal [REQUIRED] Path to tumour BAM file.
#' @param ref_genome [REQUIRED] Path to reference genome FASTA
#' @param targeted [REQUIRED] Remove coverage filtering for exome/targeted data. Default TRUE
#' @export



call_somatic_sv_manta=function(
    bin_bcftools=build_default_tool_binary_list()$bin_bcftools,
    bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
    bin_tabix=build_default_tool_binary_list()$bin_tabix,
    bin_vep=build_default_tool_binary_list()$bin_vep,
    cache_vep=build_default_cache_list()$cache_vep,
    bin_manta=build_default_tool_binary_list()$bin_manta,
    tumour=NULL,
    normal=NULL,
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    targeted=TRUE,
    annotate=TRUE,
    tabulate=TRUE,
    ...
){
   
    run_main=function(
        .env
    ){  
      
        .this.env=environment()
        append_env(to=.this.env,from=.env)
        set_main(.env=.this.env)


         .main$exec_code=paste0(
          bin_manta,
          " --tumorBam ",input,
          " --normalBam ",normal,
          " --referenceFasta ", ref_genome ,
          " --runDir ", out_file_dir, 
          ifelse(targeted," --exome ",""),"; ",
          paste0(
            out_file_dir,
            "/runWorkflow.py -m local -j ",
            threads)
          )

        .main$out_files$strelka$workflow=paste0(out_file_dir,"/runWorkflow.py")
        .main$out_files$strelka$stats=list(
              tsv=paste0(out_file_dir,"/results/stats/runStats.tsv"),
              xml=paste0(out_file_dir,"/results/stats/runStats.xml")
        )
       .main$out_files$strelka$variants=list(
              sv=paste0(out_file_dir,"/results/variants/somaticSV.vcf.gz"),
              sv_idx=paste0(out_file_dir,"/results/variants/somaticSV.vcf.gz"),
              indel_candidates=paste0(out_file_dir,"/results/variants/candidateSmallIndels.vcf.gz"),
              indel_candidates_idx=paste0(out_file_dir,"/results/variants/candidateSmallIndels.vcf.gz.tbi"),
              sv_candidate=paste0(out_file_dir,"/results/variants/candidateSV.vcf.gz"),
              sv_candidate_index=paste0(out_file_dir,"/results/variants/candidateSV.vcf.gz.tbi"),
              diploid_sv=paste0(out_file_dir,"/results/variants/diploidSV.vcf.gz"),
              diploid_sv_idx=paste0(out_file_dir,"/results/variants/diploidSV.vcf.gz.tbi")
        )

        run_job(
          .env=.this.env
        )
  
        .main.step=.main$steps[[fn]]

    
        .main.step$steps <- append(
          .main.step$steps, 
          add_af_strelka_vcf(
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            vcf=.main.step$out_files$strelka$variants$sv,
            type="sv",
            output_dir=paste0(out_file_dir,"/annotated"),
            tmp_dir=tmp_dir,
            env_dir=env_dir,
            batch_dir=batch_dir,
            err_msg=err_msg,
            verbose=verbose, 
            executor_id=task_id
          )
        )
        .this.step=.main.step$steps$add_af_strelka_vcf.sv
        .main.step$out_files$annotated$af=.this.step$out_files
      
       .main.step$steps <-append(
          .main.step$steps, 
          extract_pass_variants_strelka_vcf(
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            vcf=.main.step$out_files$annotated$af$sv$bgzip_vcf,
            type="sv",
            output_dir=paste0(out_file_dir,"/annotated"),
            tmp_dir=tmp_dir,
            env_dir=env_dir,
            batch_dir=batch_dir,
            err_msg=err_msg,
            verbose=verbose,
            executor_id=task_id
          )
        )

        .this.step=.main.step$steps$extract_pass_variants_strelka_vcf.sv
        .main.step$out_files$annotated$filter=.this.step$out_files

        if(annotate){
            .main.step$steps<-append(
              .main.step$steps,
              annotate_strelka_vep(
                bin_vep=bin_vep,
                bin_bgzip=bin_bgzip,
                bin_tabix=bin_tabix,
                cache_vep=cache_vep,
                vcf=.main.step$out_files$annotated$filter$sv$bgzip_vcf,
                type="sv",
                output_dir=paste0(out_file_dir,"/annotated"),
                tmp_dir=tmp_dir,
                env_dir=env_dir,
                batch_dir=batch_dir,
                err_msg=err_msg,
                verbose=verbose,
                executor_id=task_id
            )
          )
           .this.step=.main.step$steps$annotate_strelka_vep.sv
          .main.step$out_files$annotated$vep=.this.step$out_files

        
          if(tabulate){
              .main.step$steps<-append(
                .main.step$steps,
                tabulate_strelka_vcf(
                  vcf=.main.step$out_files$annotated$vep$sv$bgzip_vcf,
                  type="sv",
                  output_dir=paste0(out_file_dir,"/annotated"),
                  tmp_dir=tmp_dir,
                  env_dir=env_dir,
                  batch_dir=batch_dir,
                  err_msg=err_msg,
                  threads=threads,
                  ram=ram,
                  verbose=verbose,
                  executor_id=task_id
              )
            )
            .this.step=.main.step$steps$tabulate_strelka_vep.sv
            .main.step$out_files$annotated$tabulate=.this.step$out_files
          }
        }

        .env$.main <- .main
    }

    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="tumour"
    )
  
    launch(.env=.base.env)


}

