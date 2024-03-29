
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
    bin_samtools=build_default_tool_binary_list()$bin_samtools,
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
    tumour_id=NULL,
    normal_id=NULL,
    chromosomes=c(1:22,"X","Y","MT"),
    regions=NULL,
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

      .main$steps[[fn_id]]<-.this.env
      .main.step<-.main$steps[[fn_id]]

      
      if(is.list(tumour)){
        tumour=input
      }else if(is.list(normal)){
        normal=input
      }

      ### IF NO REGION IS GIVEN SPLIT BAM PER CHROMOSOME
      if(is.null(regions)){
        .main.step$steps <-append(
          .main.step$steps,
          get_sq_bam(
            bin_samtools=bin_samtools,
            bin_bcftools=bin_bcftools,
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            bam=input,
            output_name=input_id,
            output_dir=tmp_dir,
            header=FALSE,
            index=TRUE,
            compress=TRUE,
            chromosomes=chromosomes,
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
        .this.step=.main.step$steps$get_sq_bam
        .main.step$out_files$regions=.this.step$out_files$bgzip_bed
        regions=.main.step$out_files$regions
      }

      .main.step$steps <-append(
            .main.step$steps,
            call_sv_manta(
                bin_samtools=bin_samtools,
                bin_bcftools=bin_bcftools,
                bin_bgzip=bin_bgzip,
                bin_tabix=bin_tabix,
                bin_manta=bin_manta,
                bin_vep=bin_vep,
                cache_vep=cache_vep,
                tumour=tumour,
                normal=normal,
                chromosomes=chromosomes,
                patient_id=patient_id,
                tumour_id=tumour_id,
                normal_id=normal_id,
                annotate=annotate,
                tabulate=tabulate,
                regions=regions,
                ref_genome=ref_genome,
                output_dir=paste0(
                  out_file_dir,"/",patient_id,
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

      
      .main.step$steps <-append(
        .main.step$steps,
            call_snvs_strelka(
              bin_samtools=bin_samtools,
              bin_bcftools=bin_bcftools,
              bin_bgzip=bin_bgzip,
              bin_tabix=bin_tabix,
              bin_vep=bin_vep,
              bin_strelka_somatic=bin_strelka_somatic,
              bin_strelka_germline=bin_strelka_germline,
              cache_vep=cache_vep,
              annotate=annotate,
              tabulate=tabulate,
              regions=regions,
              chromosomes=chromosomes,
              indel_candidates=.main.step$out_files$sv$indel_candidates,
              tumour=tumour,
              normal=normal,
              patient_id=patient_id,
              tumour_id=tumour_id,
              normal_id=normal_id,
              ref_genome=ref_genome,
              output_dir=paste0(
                  out_file_dir,"/",patient_id,
                  "/strelka_reports"
              ),
              targeted=targeted,
              verbose=verbose,
              threads=threads,
              ram=ram,
              executor_id=task_id
        )
      )

      .this.step=.main.step$steps$call_snvs_strelka
      .main.step$out_files$snvs=.this.step$out_files

      .env$.main <-.main
    }

    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars=ifelse(!is.null(tumour),"tumour","normal")
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
    bin_samtools=build_default_tool_binary_list()$bin_samtools,
    bin_bcftools=build_default_tool_binary_list()$bin_bcftools,
    bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
    bin_tabix=build_default_tool_binary_list()$bin_tabix,
    bin_vep=build_default_tool_binary_list()$bin_vep,
    bin_strelka_somatic=build_default_tool_binary_list()$bin_strelka$somatic,
    bin_strelka_germline=build_default_tool_binary_list()$bin_strelka$germline,
    cache_vep=build_default_cache_list()$cache_vep,
    tumour=NULL,
    normal=NULL,
    patient_id=NULL,
    tumour_id=NULL,
    normal_id=NULL,
    indel_candidates=NULL,
    chromosomes=c(1:22,"X","Y","MT"),
    regions=NULL,
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

      .main$steps[[fn_id]]<-.this.env
      .main.step<-.main$steps[[fn_id]]
   
  
      if(!is.null(normal)){
        if(is.null(normal_id)){
          normal_id=get_file_name(normal)
        }

        if(!is.null(tumour)){

          if(is.null(tumour_id)){
            tumour_id=get_file_name(tumour)
          }


            .main.step$steps<-append(
              .main.step$steps,
              call_somatic_snvs_strelka(
                bin_samtools=bin_samtools,
                bin_strelka=bin_strelka_somatic,
                bin_bcftools=bin_bcftools,
                bin_bgzip=bin_bgzip,
                bin_tabix=bin_tabix,
                bin_vep=bin_vep,
                cache_vep=cache_vep,
                tumour=input,
                normal=normal,
                chromosomes=chromosomes,
                patient_id=patient_id,
                tumour_id=tumour_id,
                normal_id=normal_id,
                annotate=annotate,
                tabulate=tabulate,
                regions=regions,
                ref_genome=ref_genome,
                output_dir=paste0(out_file_dir,"/",input_id,"/somatic/snvs_indels"),
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

        }else{
            .main.step$steps<-append(
              .main.step$steps,
              call_germline_snvs_strelka(
                bin_samtools=bin_samtools,
                bin_strelka=bin_strelka_germline,
                bin_bcftools=bin_bcftools,
                bin_bgzip=bin_bgzip,
                bin_tabix=bin_tabix,
                bin_vep=bin_vep,
                cache_vep=cache_vep,
                normal=input,
                patient_id=patient_id,
                normal_id=normal_id,
                chromosomes=chromosomes,
                annotate=annotate,
                tabulate=tabulate,
                regions=regions,
                ref_genome=ref_genome,
                output_dir=paste0(out_file_dir,"/",input_id,"/germline/snvs_indels"),
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

    }
      .env$.main<-.main

  }
    
    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars=ifelse(!is.null(tumour),"tumour","normal")
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
    bin_samtools=build_default_tool_binary_list()$bin_samtools,
    bin_bcftools=build_default_tool_binary_list()$bin_bcftools,
    bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
    bin_tabix=build_default_tool_binary_list()$bin_tabix,
    bin_vep=build_default_tool_binary_list()$bin_vep,
    bin_strelka=build_default_tool_binary_list()$bin_strelka$somatic,
    cache_vep=build_default_cache_list()$cache_vep,
    tumour=NULL,
    normal=NULL,
    patient_id=NULL,
    tumour_id=NULL,
    normal_id=NULL,
    chromosomes=c(1:22,"X","Y","MT"),
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    indel_candidates=NULL,
    regions=NULL,
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
          ifelse(!is.null(regions),paste0(" --callRegions ",regions),""),
          ifelse(targeted," --exome ",""),
          ifelse(!is.null(indel_candidates),
            paste0("--indelCandidates ",indel_candidates),""), "; ",
          paste0(
            out_file_dir,
            "/runWorkflow.py -m local -j ",
            threads
            )
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

        .main.step=.main$steps[[fn_id]]


        .main.step$steps=append(.main.step$steps,
          annotate_somatic_output_strelka(
            chromosomes=ifelse(targeted,NULL,chromosomes),
            vcf_list=list(list(
              snv=.main.step$out_files$strelka$variants$snv,
              indel=.main.step$out_files$strelka$variants$indel)),
              bin_samtools=bin_samtools,
              bin_bcftools=bin_bcftools,
              bin_bgzip=bin_bgzip,
              bin_tabix=bin_tabix,
              bin_vep=bin_vep,
              bin_strelka=bin_strelka,
              cache_vep=cache_vep,
              annotate=annotate,
              patient_id=patient_id,
              tumour_id=tumour_id,
              normal_id=normal_id,
              output_dir=out_file_dir,
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


          .env$.main<-.main
      }


    
    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="tumour"
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
    bin_samtools=build_default_tool_binary_list()$bin_samtools,
    bin_bcftools=build_default_tool_binary_list()$bin_bcftools,
    bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
    bin_tabix=build_default_tool_binary_list()$bin_tabix,
    bin_vep=build_default_tool_binary_list()$bin_vep,
    cache_vep=build_default_cache_list()$cache_vep,
    bin_strelka=build_default_tool_binary_list()$bin_strelka$germline,
    normal=NULL,
    patient_id=NULL,
    normal_id=NULL,
    chromosomes=c(1:22,"X","Y","MT"),
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
          ifelse(!is.null(regions),paste0(" --callRegions ",regions),""),
          ifelse(targeted," --exome ",""),
          ifelse(!is.null(indel_candidates),
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
              snv_and_indel=paste0(out_file_dir,"/results/variants/variants.vcf.gz"),
              snv_and_indel_index=paste0(out_file_dir,"/results/variants/variants.vcf.gz.tbi"),
              genome=paste0(out_file_dir,"/results/variants/genome.S1.vcf.gz"),
              genome_idx=paste0(out_file_dir,"/results/variants/genome.S1.vcf.gz.tbi")
        )
      
        run_job(
              .env=.this.env
        )
        run_mode="local"
        if(!targeted){
          run_mode="local_parallel"
        }else{
          chromosomes=NULL
        }

        .main.step<-.main$steps[[fn_id]]

        .main.step$steps=append(.main.step$steps,
        annotate_germline_output_strelka(
          chromosomes=chromosomes,
            vcf=.main.step$out_files$strelka$variants$snv_and_indel,
            bin_samtools=bin_samtools,
            bin_bcftools=bin_bcftools,
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            bin_vep=bin_vep,
            bin_strelka=bin_strelka,
            cache_vep=cache_vep,
            annotate=annotate,
            patient_id=patient_id,
            normal_id=normal_id,
            output_dir=out_file_dir,
            tmp_dir=tmp_dir,
            env_dir=env_dir,
            batch_dir=batch_dir,
            err_msg=err_msg,
            verbose=verbose, 
            threads=threads,
            mode=run_mode,
            ram=ram,
            executor_id=task_id
          )
        )

        .env$.main<-.main

      }


    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="normal"
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

annotate_germline_output_strelka<-function(
    bin_samtools=build_default_tool_binary_list()$bin_samtools,
    bin_bcftools=build_default_tool_binary_list()$bin_bcftools,
    bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
    bin_tabix=build_default_tool_binary_list()$bin_tabix,
    bin_vep=build_default_tool_binary_list()$bin_vep,
    bin_strelka=build_default_tool_binary_list()$bin_strelka$somatic,
    cache_vep=build_default_cache_list()$cache_vep,
    chromosomes=NULL,
    vcf=NULL,
    annotate=TRUE,
    patient_id=NULL,
    normal_id=NULL,
    ...
  ){

     run_main=function(
        .env
    ){
        .this.env=environment()
        append_env(to=.this.env,from=.env)
        set_main(.env=.this.env)
        .main.step=.main$steps[[fn_id]]

        if(!is.null(chromosomes)){
          chromosomes=input
        }else{
          vcf=input
        }
        
     
        if(mode=="local_parallel"){
          threads=1
        }

       
            ### ADD AF for SNVS and INDELS 

        .main.step$steps<- append(
         .main.step$steps,
           add_af_strelka_vcf(
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            vcf=vcf,
            type="germline",
            fn_id="germline",
            output_dir=paste0(out_file_dir,"/annotated/af"),
            chromosomes=chromosomes,
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

        .this.step=.main.step$steps$add_af_strelka_vcf.germline
        .main.step$out_files$annotated$af$germline=.this.step$out_files


        ### FILTER VARIANTS BY FILTERS

      
        .main.step$steps<-append(
          .main.step$steps, 
          extract_pass_variants_strelka_vcf(
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            vcf=.main.step$out_files$annotated$af$germline$bgzip_vcf,
            type="germline",
            fn_id="germline",
            output_dir=paste0(out_file_dir,"/annotated/filter"),
            chromosomes=chromosomes,
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

        .this.step=.main.step$steps$extract_pass_variants_strelka_vcf.germline
        .main.step$out_files$annotated$filter$germline=.this.step$out_files
        

        if(annotate){

            ### ANNOTATE SNV USING VEP
            .main.step$steps <-append(
              .main.step$steps ,
              annotate_strelka_vep(
                bin_vep=bin_vep,
                bin_bgzip=bin_bgzip,
                bin_tabix=bin_tabix,
                cache_vep=cache_vep,
                patient_id=patient_id,
                normal_id=normal_id,
                vcf=.main.step$out_files$annotated$filter$germline$bgzip_vcf,
                chromosomes=chromosomes,
                type="germline",
                fn_id="germline",
                output_dir=paste0(out_file_dir,"/annotated/vep"),
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

          .this.step=.main.step$steps$annotate_strelka_vep.germline
          .main.step$out_files$annotated$vep$germline=.this.step$out_files
        }
        .env$.main<-.main
    }

    
    .base.env=environment()
      list2env(list(...),envir=.base.env)
      set_env_vars(
        .env= .base.env,
        vars=ifelse(!is.null(chromosomes),"chromosomes","vcf")
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




annotate_somatic_output_strelka<-function(
    bin_samtools=build_default_tool_binary_list()$bin_samtools,
    bin_bcftools=build_default_tool_binary_list()$bin_bcftools,
    bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
    bin_tabix=build_default_tool_binary_list()$bin_tabix,
    bin_vep=build_default_tool_binary_list()$bin_vep,
    bin_strelka=build_default_tool_binary_list()$bin_strelka$somatic,
    cache_vep=build_default_cache_list()$cache_vep,
    chromosomes=NULL,
    vcf_list=NULL,
    annotate=TRUE,
    patient_id=NULL,
    tumour_id=NULL,
    normal_id=NULL,
    ...
){

        run_main=function(
        .env
    ){
        .this.env=environment()
        append_env(to=.this.env,from=.env)
        set_main(.env=.this.env)
        .main.step=.main$steps[[fn_id]]

        if(!is.null(chromosomes)){
          chromosomes=input
        }else{
          vcf_list=unlist(input)
        }

        if(mode=="local_parallel"){
          threads=1
        }

        ### ADD AF for SNVS

         .main.step$steps <- append(
          .main.step$steps,
           add_af_strelka_vcf(
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            vcf=vcf_list[1],
            type="snv",
            fn_id="snv",
            output_dir=paste0(out_file_dir,"/annotated/af/snv"),
            chromosomes=chromosomes,
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

        .this.step=.main.step$steps$add_af_strelka_vcf.snv
        .main.step$out_files$annotated$af$snv=.this.step$out_files


        ### ADD AF for INDELS
        .main.step$steps <- append(
          .main.step$steps,
            add_af_strelka_vcf(
              bin_bgzip=bin_bgzip,
              bin_tabix=bin_tabix,
              vcf=vcf_list[2],
              chromosomes=chromosomes,
              type="indel",
              fn_id="indel",
              output_dir=paste0(out_file_dir,"/annotated/af/indel"),
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
  
        .this.step=.main.step$steps$add_af_strelka_vcf.indel
        .main.step$out_files$annotated$af$indel=.this.step$out_files



        ### FILTER SNV BY FILTERS
         
        .main.step$steps<-append(
          .main.step$steps, 
          extract_pass_variants_strelka_vcf(
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            vcf=.main.step$out_files$annotated$af$snv$bgzip_vcf,
            chromosomes=chromosomes,
            type="snv",
            fn_id="snv",
            output_dir=paste0(out_file_dir,"/annotated/filter/snv"),
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
        .main.step$out_files$annotated$filter$snv=.this.step$out_files



        ### FILTER INDELS BY FILTERS

      
        .main.step$steps<-append(
          .main.step$steps, 
          extract_pass_variants_strelka_vcf(
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            vcf=.main.step$out_files$annotated$af$indel$bgzip_vcf,
            chromosomes=chromosomes,
            type="indel",
            fn_id="indel",
            output_dir=paste0(out_file_dir,"/annotated/filter/indel"),
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
        .main.step$out_files$annotated$filter$indel=.this.step$out_files




        if(annotate){

            ### ANNOTATE SNV USING VEP
            .main.step$steps <-append(
              .main.step$steps,
              annotate_strelka_vep(
                bin_vep=bin_vep,
                bin_bgzip=bin_bgzip,
                bin_tabix=bin_tabix,
                cache_vep=cache_vep,
                patient_id=patient_id,
                tumour_id=tumour_id,
                normal_id=normal_id,
                vcf=.main.step$out_files$annotated$filter$snv$bgzip_vcf,
                chromosomes=chromosomes,
                type="snv",
                fn_id="snv",
                output_dir=paste0(out_file_dir,"/annotated/vep/snv"),
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

          .this.step=.main.step$steps$annotate_strelka_vep.snv
          .main.step$out_files$annotated$vep$snv=.this.step$out_files


          ### ANNOTATE INDEL USING VEP

          .main.step$steps <-append(
              .main.step$steps ,
              annotate_strelka_vep(
                bin_vep=bin_vep,
                bin_bgzip=bin_bgzip,
                bin_tabix=bin_tabix,
                cache_vep=cache_vep,
                patient_id=patient_id,
                tumour_id=tumour_id,
                normal_id=normal_id,
                vcf=.main.step$out_files$annotated$filter$indel$bgzip_vcf,
                chromosomes=chromosomes,
                type="indel",
                fn_id="indel",
                output_dir=paste0(out_file_dir,"/annotated/vep/indel"),
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

          .this.step=.main.step$steps$annotate_strelka_vep.indel
          .main.step$out_files$annotated$vep$indel=.this.step$out_files
        
        }
        .env$.main<-.main
    }

    .base.env=environment()
      list2env(list(...),envir=.base.env)
      set_env_vars(
        .env= .base.env,
        vars=ifelse(!is.null(chromosomes),"chromosomes","vcf_list")
      )

    launch(.env=.base.env)


}
      