
#' Manta and Strelka wrapper for variant calling
#'
#' This function wraps the functions for variant calling with Strelka and Manta
#' 
#' @param bin_strelka_somatic Path to strelka somatic workflow binary
#' @param bin_strelka_germline Path to strelka germline workflow binary
#' @param bin_manta Path to manta pipeline binary
#' @param tumour [OPTIONAL] Path to tumour BAM file. If not given will assume germline variant calling.
#' @param normal [REQUIRED] Path to tumour BAM file.
#' @param ref_genome [REQUIRED] Path to reference genome FASTA
#' @param variants [REQUIRED] Variants types to call. Default all. Options ["snv","sv","all"]
#' @param indel_cnds [REQUIRED] Use indel candidates to correct SNV calls. Default TRUE
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


call_variants_strelka=function(
    bin_bcftools=build_default_tool_binary_list()$bin_bcftools,
    bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
    bin_tabix=build_default_tool_binary_list()$bin_tabix,
    bin_vep=build_default_tool_binary_list()$bin_vep,
    cache_vep=build_default_cache_list()$cache_vep,
    bin_strelka_somatic=build_default_tool_binary_list()$bin_strelka$somatic,
    bin_strelka_germline=build_default_tool_binary_list()$bin_strelka$germline,
    bin_manta=build_default_tool_binary_list()$bin_manta,
    tumour=NULL,
    normal=NULL,
    patient_id=NULL,
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    annotate=TRUE,
    tabulate=TRUE,
    targeted=TRUE,
    ...
){




    run_main=function(
      .env
    ){

      .this.env=environment()
      append_env(to=.this.env,from=.env)
  
      set_main(.env=.this.env)

      .main$steps[[fn]]<-.this.env
    
      .main$steps[[fn]]$steps <-append(
            .main$steps[[fn]]$steps,
            call_sv_manta(
                bin_bcftools=bin_bcftools,
                bin_bgzip=bin_bgzip,
                bin_tabix=bin_tabix,
                bin_manta=bin_manta,
                bin_vep=bin_vep,
                cache_vep=cache_vep,
                tumour=input,
                normal=normal,
                annotate=annotate,
                tabulate=tabulate,
                ref_genome=ref_genome,
                output_dir=paste0(
                  out_file_dir,"/",
                  input_id,
                  "/strelka_reports"
                ),
                targeted=targeted,
                verbose=verbose,
                threads=threads,
                ram=ram,
                executor_id=task_id
        )
    )
      .this.step=.main.step$steps$call_sv_manta
      .main.step$out_files$sv=.this.step$out_files

      
      .main$steps[[fn]]$steps <-append(
        .main$steps[[fn]]$steps,
            call_snvs_strelka(
              bin_bcftools=bin_bcftools,
              bin_bgzip=bin_bgzip,
              bin_tabix=bin_tabix,
              bin_vep=bin_vep,
              bin_strelka_somatic=bin_strelka_somatic,
              bin_strelka_germline=bin_strelka_germline,
              cache_vep=cache_vep,
              annotate=annotate,
              tabulate=tabulate,
              indel_candidates=.main.step$out_files$indel_candidates,
              tumour=tumour,
              normal=normal,
              ref_genome=ref_genome,
              output_dir=paste0(
                  out_file_dir,"/",
                  input_id,
                  "/strelka_reports"
              ),
              targeted=targeted,
              verbose=verbose,
              threads=threads,
              ram=ram,
              executor_id=task_id
        )
      )

      .this.step=.main.step$steps$call_snvs_manta
      .main.step$out_files$snvs=.this.step$out_files

      .env$.main <-.main
    }

    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars=ifelse(is.null(tumour),"normal","tumour")
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
          bin_strelka,
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
              diploid_sv_idx=paste0(out_file_dir,"/results/variants/diploidSV.vcf.gz.tbi"),
        )

        run_job(
          .env=.this.env
        )
  
        .main.step=.main$steps[[fn]]

    
        .main$steps[[fn]]$steps <- append(
          .main$steps[[fn]]$steps, 
          add_af_strelka_vcf(
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            vcf_sv=.main.step$out_files$strelka$variants$sv,
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
      
       .main$steps[[fn]]$steps <-append(
          .main$steps[[fn]]$steps, 
          extract_pass_variants_strelka_vcf(
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            vcf=.main.step$out_files$annotated$af$bgzip_vcf,
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
            .main$steps[[fn]]$steps<-append(
              .main$steps[[fn]]$steps,
              annotate_strelka_vep(
                bin_vep=bin_vep,
                bin_bgzip=bin_bgzip,
                bin_tabix=bin_tabix,
                cache_vep=cache_vep,
                vcf=.main.step$out_files$annotated$filter$bgzip_vcf,
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
              steps[[fn]]$steps<-append(
                steps[[fn]]$steps,
                tabulate_strelka_vcf(
                  vcf=.main.step$out_files$annotated$vep,
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
            .this.step=.main.step$steps$annotate_strelka_vep.sv
            .main.step$out_files$annotated$tabulate=.this.step$out_files
          }
        }

        .env$.main <- .main
    }

    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars=ifelse(is.null(tumour),"normal","tumour")
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
          bin_strelka,
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
              diploid_sv_idx=paste0(out_file_dir,"/results/variants/diploidSV.vcf.gz.tbi"),
        )

        run_job(
            .env=.this.env
        )

        .main.step=.main$steps[[fn]]

    
        .main$steps[[fn]]$steps <- append(
          .main$steps[[fn]]$steps, 
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
      
        .main$steps[[fn]]$steps <-append(
          .main$steps[[fn]]$steps, 
          extract_pass_variants_strelka_vcf(
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            vcf=.main.step$out_files$annotated$af$bgzip_vcf,
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
            .main$steps[[fn]]$steps<-append(
              .main$steps[[fn]]$steps,
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
              .main$steps[[fn]]$steps<-append(
                .main$steps[[fn]]$steps,
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


call_snvs_strelka=function( 
    bin_bcftools=build_default_tool_binary_list()$bin_bcftools,
    bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
    bin_tabix=build_default_tool_binary_list()$bin_tabix,
    bin_vep=build_default_tool_binary_list()$bin_vep,
    bin_strelka_somatic=build_default_tool_binary_list()$bin_strelka$somatic,
    bin_strelka_germline=build_default_tool_binary_list()$bin_strelka$germline,
    cache_vep=build_default_cache_list()$cache_vep,
    tumour=NULL,
    normal=NULL,
    indel_candidates=NULL,
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    targeted=TRUE,
    verbose=TRUE,
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

      .main$steps[[fn]]<-.this.env
      tumour<-ifelse(is.null(tumour),input,NULL)
      normal<-ifelse(is.null(tumour),normal,input)
    

      if(!is.null(normal)){

        if(!is.null(tumour)){

            .main$steps[[fn]]$steps<-append(
              .main$steps[[fn]]$steps,
              call_somatic_snvs_strelka(
                bin_strelka=bin_strelka_somatic,
                bin_bcftools=bin_bcftools,
                bin_bgzip=bin_bgzip,
                bin_tabix=bin_tabix,
                bin_vep=bin_vep,
                cache_vep=cache_vep,
                tumour=tumour,
                normal=normal,
                annotate=annotate,
                tabulate=tabulate,
                ref_genome=ref_genome,
                output_dir=paste0(out_file_dir,"/somatic/snvs_indels"),
                tmp_dir=tmp_dir,
                env_dir=env_dir,
                batch_dir=batch_dir,
                err_msg=err_msg,
                indel_candidates=indel_candidates,
                targeted=targeted,
                verbose=verbose,
                threads=threads,
                ram=ram,
                executor_id=task_id

              )
            )

          .this.step=.main.step$steps$call_somatic_snvs_strelka
          .main.step$out_files$snvs$somatic=.this.step$out_files

        }
       
      }else{
            .main$steps[[fn]]$steps<-append(
              .main$steps[[fn]]$steps,
              call_germline_snvs_strelka(
                bin_strelka=bin_strelka_germline,
                bin_bcftools=bin_bcftools,
                bin_bgzip=bin_bgzip,
                bin_tabix=bin_tabix,
                bin_vep=bin_vep,
                cache_vep=cache_vep,
                normal=normal,
                annotate=annotate,
                tabulate=tabulate,
                ref_genome=ref_genome,
                output_dir=paste0(out_file_dir,"/germline/snvs_indels"),
                tmp_dir=tmp_dir,
                env_dir=env_dir,
                batch_dir=batch_dir,
                err_msg=err_msg,
                indel_candidates=indel_candidates,
                targeted=targeted,
                verbose=verbose,
                threads=threads,
                ram=ram,
                executor_id=task_id
            )
           )
          .this.step=.main.step$steps$call_germline_snvs_strelka
          .main.step$out_files$snvs$germline=.this.step$out_files
      }
      .env$.main<-.main

    }

    
    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars=ifelse(is.null(tumour),"normal","tumour")
    )

    launch(.env=.base.env)
  
  

}




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


call_sv_strelka=function( 
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
  


    build_main=function(
        .env
    ){
      .this.env=environment()
      append_env(to=.this.env,from=.env)
      set_main(.env=.this.env)
      tumour<-ifelse(is.null(tumour),input,NULL)
      normal<-ifelse(is.null(tumour),normal,input)

      if(!is.null(normal)){

        if(!is.null(tumour)){

            .main$steps[[fn]]$steps<-append(
              .main$steps[[fn]]$steps,
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
           .main$steps[[fn]]$steps<-append(
              .main$steps[[fn]]$steps,
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
      vars=ifelse(is.null(tumour),"normal","tumour")
    )

    launch(.env=.base.env)

    
}




#' Strelka wrapper for somatic SNV variant calling
#'
#' This function wraps the STRELKA functions for somatic variant calling
#' 
#' @param bin_strelka_somatic Path to strelka somatic workflow binary
#' @param tumour [OPTIONAL] Path to tumour BAM file.
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




call_somatic_snvs_strelka=function(
    bin_bcftools=build_default_tool_binary_list()$bin_bcftools,
    bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
    bin_tabix=build_default_tool_binary_list()$bin_tabix,
    bin_vep=build_default_tool_binary_list()$bin_vep,
    bin_strelka=build_default_tool_binary_list()$bin_strelka$somatic,
    cache_vep=build_default_cache_list()$cache_vep,
    tumour=NULL,
    normal=NULL,
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    indel_candidates=NULL,
    annotate=TRUE,
    tabulate=TRUE,
    targeted=TRUE,
    ...
){

    

    run_main=function(
        .env
    ){
         .this.env=environment()
        append_env(to=.this.env,from=.env)
        set_main(.env=.this.env)


          
        .main$exec_code=paste0(
          bin_strelka,
          " --tumorBam ",input,
          " --normalBam ",normal,
          " --referenceFasta ", ref_genome ,
          " --runDir ", out_file_dir, 
          ifelse(targeted," --exome ",""),
          ifelse(indel_candidates,
            paste0("--indelCandidates ",indel_candidates),""), "; ",
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
              indel=paste0(out_file_dir,"/results/variants/somatic.indels.vcf.gz"),
              indel_idx=paste0(out_file_dir,"/results/variants/somatic.indels.vcf.gz.tbi"),
              snv=paste0(out_file_dir,"/results/variants/somatic.snvs.vcf.gz"),
              snv_idx=paste0(out_file_dir,"/results/variants/somatic.snvs.vcf.gz.tbi")
        )

        run_job(
            .env=.this.env
        )

        .main.step=.main$steps[[fn]]


        
        ### ADD AF for SNVS

         .main$steps[[fn]]$steps <- append(
          .main$steps[[fn]]$steps,
           add_af_strelka_vcf(
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            vcf=.main$steps$out_files$strelka$snv,
            type="snv",
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

        .this.step=.main.step$steps$add_snv_af_strelka_vcf.snv
        .main.step$out_files$annotated$af=.this.step$out_files


        ### ADD AF for INDELS
        .main$steps[[fn]]$steps <- append(.main$steps[[fn]]$steps,
            add_af_strelka_vcf(
              bin_bgzip=bin_bgzip,
              bin_tabix=bin_tabix,
              vcf=.main$steps$out_files$strelka$indel,
              type="indel",
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
  
        .this.step=.main.step$steps$add_indel_af_strelka_vcf.indel
        .main.step$out_files$annotated$af=append(
          .main.step$out_files$annotated$af,
          .this.step$out_files
        )


        ### FILTER SNV BY FILTERS
         
        .main$steps[[fn]]$steps <-append(
          .main$steps[[fn]]$steps, 
          extract_pass_variants_strelka_vcf(
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            vcf=.main$steps$out_files$annotated$af$snv$bgzip_vcf,
            type="snv",
            output_dir=paste0(out_file_dir,"/annotated"),
            tmp_dir=tmp_dir,
            env_dir=env_dir,
            batch_dir=batch_dir,
            err_msg=err_msg,
            verbose=verbose,
            threads=threads,
            ram=ram,
            executor_id=task_id,
          )
        )

        .this.step=.main.step$steps$extract_pass_variants_strelka_vcf.snv
        .main.step$out_files$annotated$filter=.this.step$out_files



        ### FILTER INDELS BY FILTERS

      
        .main$steps[[fn]]$steps<-append(
          .main$steps[[fn]]$steps, 
          extract_pass_variants_strelka_vcf(
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            vcf=.main$steps$out_files$annotated$af$indel$bgzip_vcf,
            type="indel",
            output_dir=paste0(out_file_dir,"/annotated"),
            tmp_dir=tmp_dir,
            env_dir=env_dir,
            batch_dir=batch_dir,
            err_msg=err_msg,
            verbose=verbose,
            threads=threads,
            ram=ram,
            executor_id=task_id
          )
        )

        .this.step=.main.step$steps$extract_pass_variants_strelka_vcf.indel
        .main.step$out_files$annotated$filter=append(
          .main.step$out_files$annotated$filter,
          .this.step$out_files
        )



        if(annotate){

            ### ANNOTATE SNV USING VEP
            .main$steps[[fn]]$steps <-append(
              .main$steps[[fn]]$steps ,
              annotate_strelka_vep(
                bin_vep=bin_vep,
                bin_bgzip=bin_bgzip,
                bin_tabix=bin_tabix,
                cache_vep=cache_vep,
                vcf=.main.step$out_files$annotated$filter$snv$bgzip_vcf,
                type="snv",
                output_dir=paste0(out_file_dir,"/annotated"),
                tmp_dir=tmp_dir,
                env_dir=env_dir,
                batch_dir=batch_dir,
                err_msg=err_msg,
                verbose=verbose,
                threads=threads,
                ram=ram,
                executor_id=task_id
            )
          )

          .this.step=.main.step$steps$annotate_snv_strelka_vep.snv
          .main.step$out_files$annotated$vep=.this.step$out_files


          ### ANNOTATE INDEL USING VEP

          .main$steps[[fn]]$steps <-append(
              .main$steps[[fn]]$steps ,
              annotate_strelka_vep(
                bin_vep=bin_vep,
                bin_bgzip=bin_bgzip,
                bin_tabix=bin_tabix,
                cache_vep=cache_vep,
                vcf=.main.step$out_files$annotated$filter$indel$bgzip_vcf,
                type="indel",
                output_dir=paste0(out_file_dir,"/annotated"),
                tmp_dir=tmp_dir,
                env_dir=env_dir,
                batch_dir=batch_dir,
                err_msg=err_msg,
                verbose=verbose,
                threads=threads,
                ram=ram,
                executor_id=task_id
            )
          )

          .this.step=.main.step$steps$annotate_indel_strelka_vep.indel
          .main.step$out_files$annotated$vep=append(
            .main.step$out_files$annotated$vep,
            .this.step$out_files
          )

        
          if(tabulate){
              .main$steps[[fn]]$steps <-append(
                .main$steps[[fn]]$steps ,
                tabulate_strelka_vcf(
                  vcf=.main.step$out_files$annotated$vep$snv$bgzip_vcf,
                  type="snv",
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

          .this.step=.main.step$steps$tabulate_strelka_vcf.snv
          .main.step$out_files$annotated$tabulate=.this.step$out_files


              .main$steps[[fn]]$steps <-append(
                .main$steps[[fn]]$steps ,
                tabulate_strelka_vcf(
                  vcf=.main.step$out_files$annotated$vep$indel$bgzip_vcf,
                  type="indel",
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


          .this.step=.main.step$steps$tabulate_strelka_vcf.indel
            .main.step$out_files$annotated$tabulate=append(
            .main.step$out_files$annotated$tabulate,
              .this.step$out_files
            )

          }
        }


          .env$.main<-.main
      }


    
    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars=ifelse(is.null(tumour),"normal","tumour")
    )

    launch(.env=.base.env)

}



#' Strelka wrapper for germline SNV variant calling
#'
#' This function wraps the STRELKA functions for germline variant calling
#' 
#' @param bin_strelka_somatic Path to strelka somatic workflow binary
#' @param tumour [OPTIONAL] Path to tumour BAM file.
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


call_germline_snvs_strelka=function(
    bin_bcftools=build_default_tool_binary_list()$bin_bcftools,
    bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
    bin_tabix=build_default_tool_binary_list()$bin_tabix,
    bin_vep=build_default_tool_binary_list()$bin_vep,
    cache_vep=build_default_cache_list()$cache_vep,
    bin_strelka=build_default_tool_binary_list()$bin_strelka$germline,
    normal=NULL,
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    indel_candidates=NULL,
    extract_pass=TRUE,
    annotate=TRUE,
    tabulate=TRUE,
    targeted=TRUE,
    ...
){
  
 
    run_main=function(
        .env
    ){
        .this.env=environment()
        append_env(to=.this.env,from=.env)
        set_main(.env=.this.env)


          

        .main$exec_code=paste0(
          bin_strelka,
          " --bam ",normalizePath(input),
          " --referenceFasta ", normalizePath(ref_genome),
          " --runDir ", out_file_dir, 
          ifelse(targeted," --exome ",""),
          ifelse(indel_candidates,
            paste0("--indelCandidates ",indel_candidates),""), "; ",
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
              indel=paste0(out_file_dir,"/results/variants/somatic.indels.vcf.gz"),
              indel_idx=paste0(out_file_dir,"/results/variants/somatic.indels.vcf.gz.tbi"),
              snv=paste0(out_file_dir,"/results/variants/somatic.snvs.vcf.gz"),
              snv_idx=paste0(out_file_dir,"/results/variants/somatic.snvs.vcf.gz.tbi")
        )
      
        run_job(
              .env=.this.env
        )


        .main.step=.main$steps[[fn]]



         
        ### ADD AF for SNVS

         .main$steps[[fn]]$steps <- append(
          .main$steps[[fn]]$steps,
           add_af_strelka_vcf(
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            vcf=.main$steps$out_files$strelka$snv,
            type="snv",
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

        .this.step=.main.step$steps$add_snv_af_strelka_vcf.snv
        .main.step$out_files$annotated$af=.this.step$out_files


        ### ADD AF for INDELS
        .main$steps[[fn]]$steps <- append(.main$steps[[fn]]$steps,
            add_af_strelka_vcf(
              bin_bgzip=bin_bgzip,
              bin_tabix=bin_tabix,
              vcf=.main$steps$out_files$strelka$indel,
              type="indel",
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
  
        .this.step=.main.step$steps$add_indel_af_strelka_vcf.indel
        .main.step$out_files$annotated$af=append(
          .main.step$out_files$annotated$af,
          .this.step$out_files
        )


        ### FILTER SNV BY FILTERS
         
        .main$steps[[fn]]$steps <-append(
          .main$steps[[fn]]$steps, 
          extract_pass_variants_strelka_vcf(
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            vcf=.main$steps$out_files$annotated$af$snv$bgzip_vcf,
            type="snv",
            output_dir=paste0(out_file_dir,"/annotated"),
            tmp_dir=tmp_dir,
            env_dir=env_dir,
            batch_dir=batch_dir,
            err_msg=err_msg,
            verbose=verbose,
            threads=threads,
            ram=ram,
            executor_id=task_id
          )
        )

        .this.step=.main.step$steps$extract_pass_variants_strelka_vcf.snv
        .main.step$out_files$annotated$filter=.this.step$out_files



        ### FILTER INDELS BY FILTERS

      
        .main$steps[[fn]]$steps<-append(
          .main$steps[[fn]]$steps, 
          extract_pass_variants_strelka_vcf(
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            vcf=.main$steps$out_files$annotated$af$indel$bgzip_vcf,
            type="indel",
            output_dir=paste0(out_file_dir,"/annotated"),
            tmp_dir=tmp_dir,
            env_dir=env_dir,
            batch_dir=batch_dir,
            err_msg=err_msg,
            verbose=verbose,
            threads=threads,
            ram=ram,
            executor_id=task_id
          )
        )

        .this.step=.main.step$steps$extract_pass_variants_strelka_vcf.indel
        .main.step$out_files$annotated$filter=append(
          .main.step$out_files$annotated$filter,
          .this.step$out_files
        )



        if(annotate){

            ### ANNOTATE SNV USING VEP
            .main$steps[[fn]]$steps <-append(
              .main$steps[[fn]]$steps ,
              annotate_strelka_vep(
                bin_vep=bin_vep,
                bin_bgzip=bin_bgzip,
                bin_tabix=bin_tabix,
                cache_vep=cache_vep,
                vcf=.main.step$out_files$annotated$filter$snv$bgzip_vcf,
                type="snv",
                output_dir=paste0(out_file_dir,"/annotated"),
                tmp_dir=tmp_dir,
                env_dir=env_dir,
                batch_dir=batch_dir,
                err_msg=err_msg,
                verbose=verbose,
                threads=threads,
                ram=ram,
                executor_id=task_id
            )
          )

          .this.step=.main.step$steps$annotate_snv_strelka_vep.snv
          .main.step$out_files$annotated$vep=.this.step$out_files


          ### ANNOTATE INDEL USING VEP

          .main$steps[[fn]]$steps <-append(
              .main$steps[[fn]]$steps ,
              annotate_strelka_vep(
                bin_vep=bin_vep,
                bin_bgzip=bin_bgzip,
                bin_tabix=bin_tabix,
                cache_vep=cache_vep,
                vcf=.main.step$out_files$annotated$filter$indel$bgzip_vcf,
                type="indel",
                output_dir=paste0(out_file_dir,"/annotated"),
                tmp_dir=tmp_dir,
                env_dir=env_dir,
                batch_dir=batch_dir,
                err_msg=err_msg,
                verbose=verbose,
                threads=threads,
                ram=ram,
                executor_id=task_id
            )
          )

          .this.step=.main.step$steps$annotate_indel_strelka_vep.indel
          .main.step$out_files$annotated$vep=append(
            .main.step$out_files$annotated$vep,
            .this.step$out_files
          )

        
          if(tabulate){
              .main$steps[[fn]]$steps <-append(
                .main$steps[[fn]]$steps ,
                tabulate_strelka_vcf(
                  vcf=.main.step$out_files$annotated$vep$snv$bgzip_vcf,
                  type="snv",
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

          .this.step=.main.step$steps$tabulate_strelka_vcf.snv
          .main.step$out_files$annotated$tabulate=.this.step$out_files


              .main$steps[[fn]]$steps <-append(
                .main$steps[[fn]]$steps ,
                tabulate_strelka_vcf(
                  vcf=.main.step$out_files$annotated$vep$indel$bgzip_vcf,
                  type="indel",
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


          .this.step=.main.step$steps$tabulate_strelka_vcf.indel
            .main.step$out_files$annotated$tabulate=append(
            .main.step$out_files$annotated$tabulate,
              .this.step$out_files
            )

          }
        }


          .env$.main<-.main

      }


    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars=ifelse(is.null(tumour),"normal","tumour")
    )

    launch(.env=.base.env)
      

}