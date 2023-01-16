
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
    verbose=FALSE,
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
    
      .main$steps[[fn]]$steps <-append(
            .main$steps[[fn]]$steps,
            call_sv_manta(
                bin_bcftools=bin_bcftools,
                bin_bgzip=bin_bgzip,
                bin_tabix=bin_tabix,
                bin_manta=bin_manta,
                bin_vep=bin_vep,
                cache_vep=cache_vep,
                tumour=tumour,
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
      .main.step$out_files=append(.main.step$out_files,.this.step$out_files) 

      
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
        envir
    ){
        this.envir=environment()
        append_envir(this.envir,envir)
        set_steps_vars(envir=this.envir)



         steps[[fn]]$exec_code=paste0(
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

        steps[[fn]]$out_file$workflow=paste0(out_file_dir,"/runWorkflow.py")
        steps[[fn]]$out_file$stats=list(
              tsv=paste0(out_file_dir,"/results/stats/runStats.tsv"),
              xml=paste0(out_file_dir,"/results/stats/runStats.xml")
        )
        steps[[fn]]$out_file$variants=list(
              sv=paste0(out_file_dir,"/results/variants/somaticSV.vcf.gz"),
              sv_idx=paste0(out_file_dir,"/results/variants/somaticSV.vcf.gz"),
              indel_candidate=paste0(out_file_dir,"/results/variants/candidateSmallIndels.vcf.gz"),
              indel_candidate_idx=paste0(out_file_dir,"/results/variants/candidateSmallIndels.vcf.gz.tbi"),
              sv_candidate=paste0(out_file_dir,"/results/variants/candidateSV.vcf.gz"),
              sv_candidate_index=paste0(out_file_dir,"/results/variants/candidateSV.vcf.gz.tbi"),
              diploid_sv=paste0(out_file_dir,"/results/variants/diploidSV.vcf.gz"),
              diploid_sv_idx=paste0(out_file_dir,"/results/variants/diploidSV.vcf.gz.tbi"),
        )

        run_job(
            envir=this.envir
        )


    
        steps[[fn]] <- append(
          steps[[fn]], 
          add_af_strelka_vcf(
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            vcf_sv=steps[[fn]]$out_file$variants$sv,
            output_dir=paste0(out_file_dir,"/results"),
            verbose=verbose, 
            executor_id=task_id,
            threads=threads,ram=ram
          )
        )
      
        steps[[fn]] <-append(
          steps[[fn]], 
          extract_pass_variants_strelka_vcf(
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            vcf_sv=steps[[fn]]$add_af_strelka_vcf$out_file$variants$sv,
            output_dir=paste0(out_file_dir,"/results"),
            verbose=verbose,
            threads=threads,
            ram=ram,
            executor_id=task_id,
          )
        )

        if(annotate){
            steps[[fn]]<-append(
              steps[[fn]],
              annotate_strelka_vep(
                bin_vep=bin_vep,
                bin_bgzip=bin_bgzip,
                bin_tabix=bin_tabix,
                cache_vep=cache_vep,
                vcf_sv=steps[[fn]]$extract_pass_variants_strelka_vcf$out_file$variants$sv,
                output_dir=paste0(out_file_dir,"/results"),
                verbose=verbose,
                threads=threads,
                ram=ram,
                executor_id=task_id
            )
          )

        
          if(tabulate){
              steps[[fn]]<-append(
                steps[[fn]],
                tabulate_strelka_vcf(
                  vcf_sv=steps[[fn]]$annotate_strelka_vep$out_file$variants$sv,
                  output_dir=paste0(out_file_dir,"/results"),
                  threads=threads,
                  ram=ram,
                  verbose=verbose,
                  executor_id=task_id
              )
            )
          }
        }

        envir$steps <- steps
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
    output_dir="./germline/sv",
    output_name=NULL,
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,
    ram=4,
    mode="local",
    time="48:0:0",
    sheet=NULL,
    executor_id=NULL,
    inherit=NULL,
    select=NULL,
    hold=NULL
){
   this.envir=environment()
    set_envir_vars(
        envir=this.envir,
        vars="normal"
    )

    run_main=function(
        envir
    ){
        this.envir=environment()
        append_envir(this.envir,envir)
        set_steps_vars(envir=this.envir)



         steps[[fn]]$exec_code=paste0(
          bin_strelka,
          " --bam ",input,
          " --referenceFasta ", ref_genome ,
          " --runDir ", out_file_dir, 
          ifelse(targeted," --exome ",""), "; ",
          paste0(
            out_file_dir,
            "/runWorkflow.py -m local -j ",
            threads)
          )

        steps[[fn]]$out_file$workflow=paste0(out_file_dir,"/runWorkflow.py")
        steps[[fn]]$out_file$stats=list(
              tsv=paste0(out_file_dir,"/results/stats/runStats.tsv"),
              xml=paste0(out_file_dir,"/results/stats/runStats.xml")
        )
        steps[[fn]]$out_file$variants=list(
              indel_candidate=paste0(out_file_dir,"/results/variants/candidateSmallIndels.vcf.gz"),
              indel_candidate_idx=paste0(out_file_dir,"/results/variants/candidateSmallIndels.vcf.gz.tbi"),
              diploid_sv=paste0(out_file_dir,"/results/variants/diploidSV.vcf.gz"),
              diploid_sv_idx=paste0(out_file_dir,"/results/variants/diploidSV.vcf.gz.tbi"),
        )

        run_job(
            envir=this.envir
        )


    
        steps[[fn]] <- append(
          steps[[fn]], 
          add_af_strelka_vcf(
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            vcf_sv=steps[[fn]]$out_file$variants$diploid_sv,
            output_dir=paste0(out_file_dir,"/results"),
            verbose=verbose, 
            executor_id=task_id,
            threads=threads,ram=ram
          )
        )
      
        steps[[fn]] <-append(
          steps[[fn]], 
          extract_pass_variants_strelka_vcf(
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            vcf_sv=steps[[fn]]$add_af_strelka_vcf$out_file$variants$sv,
            output_dir=paste0(out_file_dir,"/results"),
            verbose=verbose,
            threads=threads,
            ram=ram,
            executor_id=task_id,
          )
        )

        if(annotate){
            steps[[fn]]<-append(
              steps[[fn]],
              annotate_strelka_vep(
                bin_vep=bin_vep,
                bin_bgzip=bin_bgzip,
                bin_tabix=bin_tabix,
                cache_vep=cache_vep,
                vcf_sv=steps[[fn]]$extract_pass_variants_strelka_vcf$out_file$variants$sv,
                output_dir=paste0(out_file_dir,"/results"),
                verbose=verbose,
                threads=threads,
                ram=ram,
                executor_id=task_id
            )
          )

        
          if(tabulate){
              steps[[fn]]<-append(
                steps[[fn]],
                tabulate_strelka_vcf(
                  vcf_sv=steps[[fn]]$annotate_strelka_vep$out_file$variants$sv,
                  output_dir=paste0(out_file_dir,"/results"),
                  threads=threads,
                  ram=ram,
                  verbose=verbose,
                  executor_id=task_id
              )
            )
          }
        }

        envir$steps <- steps
    }


     envirs=run_envir(envirs=this.envir$envirs)


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
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    output_dir="./strelka_reports",
    output_name=NULL,
    indel_candidates=NULL,
    targeted=TRUE,
    verbose=TRUE,
    annotate=TRUE,
    tabulate=TRUE,
    batch_config=build_default_preprocess_config(),
    threads=1,
    ram=4,
    mode="local",
    time="48:0:0",
    sheet=NULL,
    inherit=NULL,
    select=NULL,
    executor_id=NULL,
    hold=NULL
){
  
 this.envir=environment()
    set_envir_vars(
        envir=this.envir,
        vars="tumour"
    )

    run_main=function(
        envir
    ){
        this.envir=environment()
        append_envir(this.envir,envir)
        set_steps_vars(envir=this.envir)

    


      if(!is.null(normal)){

        if(!is.null(tumour)){

            steps[[fn]]<-append(
              steps[[fn]],
              call_somatic_snvs_strelka(
                bin_strelka=bin_strelka_somatic,
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
                output_dir=paste0(out_file_dir,"/somatic/snvs_indels"),
                indel_candidates=indel_candidates,
                targeted=targeted,
                verbose=verbose,
                threads=threads,
                ram=ram,
                executor_id=task_id

              )
            )

        }
       
      }else{
           steps[[fn]]<-append(
              steps[[fn]],
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
                indel_candidates=indel_candidates,
                targeted=targeted,
                verbose=verbose,
                threads=threads,
                ram=ram,
                executor_id=task_id
            )
           )
      }
      envir$steps <-steps

    }

    
  
     envirs=run_envir(envirs=this.envir$envirs)

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
    output_dir="./strelka_reports",
    output_name=NULL,
    targeted=TRUE,
    verbose=TRUE,
    annotate=TRUE,
    tabulate=TRUE,
    batch_config=build_default_preprocess_config(),
    threads=1,
    ram=4,
    mode="local",
    time="48:0:0",
    sheet=NULL,
    inherit=NULL,
    select=NULL,
    executor_id=NULL,
    hold=NULL
){
  
 this.envir=environment()
    set_envir_vars(
        envir=this.envir,
        vars="tumour"
    )

    run_main=function(
        envir
    ){
        this.envir=environment()
        append_envir(this.envir,envir)
        set_steps_vars(envir=this.envir)

    


      if(!is.null(normal)){

        if(!is.null(tumour)){

            steps[[fn]]<-append(
              steps[[fn]],
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
                targeted=targeted,
                verbose=verbose,
                threads=threads,
                ram=ram,
                executor_id=task_id

              )
            )

        }
       
      }else{
           steps[[fn]]<-append(
              steps[[fn]],
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
                targeted=targeted,
                verbose=verbose,
                threads=threads,
                ram=ram,
                executor_id=task_id
            )
           )
      }
      envir$steps <-steps

    }

    
  
     envirs=run_envir(envirs=this.envir$envirs)
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
    output_dir="./somatic/snv",
    output_name=NULL,
    indel_candidates=NULL,
    annotate=TRUE,
    tabulate=TRUE,
    targeted=TRUE,
    verbose=TRUE,
    batch_config=build_default_preprocess_config(),
    threads=1,
    ram=4,
    mode="local",
    time="48:0:0",
    sheet=NULL,
    inherit=NULL,
    select=NULL,
    executor_id=NULL,
    hold=NULL
){

    this.envir=environment()
    set_envir_vars(
        envir=this.envir,
        vars="tumour"
    )

    run_main=function(
        envir
    ){
        this.envir=environment()
        append_envir(this.envir,envir)
        set_steps_vars(envir=this.envir)


          
        steps[[fn]]$exec_code=paste0(
          bin_strelka,
          " --tumorBam ",input,
          " --normalBam ",normal,
          " --referenceFasta ", ref_genome ,
          " --runDir ", out_file_dir, 
          ifelse(targeted," --exome ",""),
          ifelse(indel_candidates,paste0("--indelCandidates ",indel_candidates),""), "; ",
          paste0(
            out_file_dir,
            "/runWorkflow.py -m local -j ",
            threads)
          )

        steps[[fn]]$out_file$workflow=paste0(out_file_dir,"/runWorkflow.py")
        steps[[fn]]$out_file$stats=list(
              tsv=paste0(out_file_dir,"/results/stats/runStats.tsv"),
              xml=paste0(out_file_dir,"/results/stats/runStats.xml")
        )
        steps[[fn]]$out_file$variants=list(
              indel=paste0(out_file_dir,"/results/variants/somatic.indels.vcf.gz"),
              indel_idx=paste0(out_file_dir,"/results/variants/somatic.indels.vcf.gz.tbi"),
              snv=paste0(out_file_dir,"/results/variants/somatic.snvs.vcf.gz"),
              snv_idx=paste0(out_file_dir,"/results/variants/somatic.snvs.vcf.gz.tbi")
        )

        run_job(
            envir=this.envir
        )


        steps[[fn]] <- append(steps[[fn]], 
          add_af_strelka_vcf(
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            vcf_snv=steps[[fn]]$out_file$variants$snv,
            vcf_indel=steps[[fn]]$out_file$variants$indel,
            output_dir=out_file_dir,
            verbose=verbose, 
            executor_id=task_id,
            threads=threads,
            ram=ram
            )
        )
      
        steps[[fn]] <-append(
          steps[[fn]], 
          extract_pass_variants_strelka_vcf(
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            vcf_snv=steps[[fn]]$add_af_strelka_vcf$out_file$variants$snv,
            vcf_indel=steps[[fn]]$add_af_strelka_vcf$out_file$variants$indel,
            output_dir=out_file_dir,
            verbose=verbose,
            threads=threads,
            ram=ram,
            executor_id=task_id,
          )
        )

        if(annotate){
            steps[[fn]]<-append(
              steps[[fn]],
              annotate_strelka_vep(
                bin_vep=bin_vep,
                bin_bgzip=bin_bgzip,
                bin_tabix=bin_tabix,
                cache_vep=cache_vep,
                vcf_snv=steps[[fn]]$extract_pass_variants_strelka_vcf$out_file$variants$snv,
                vcf_indel=steps[[fn]]$extract_pass_variants_strelka_vcf$out_file$variants$indel,
                output_dir=out_file_dir,
                verbose=verbose,
                threads=threads,
                ram=ram,
                executor_id=task_id
            )
          )

        
          if(tabulate){
              steps[[fn]]<-append(
                steps[[fn]],
                tabulate_strelka_vcf(
                  vcf_snv=steps[[fn]]$annotate_strelka_vep$out_file$variants$snv,
                  vcf_indel=steps[[fn]]$annotate_strelka_vep$out_file$variants$ndel,
                  output_dir=out_file_dir,
                  threads=threads,
                  ram=ram,
                  verbose=verbose,
                  executor_id=task_id
              )
            )
          }
        }


          envir$steps <-steps
      }


    
       envirs=run_envir(envirs=this.envir$envirs)

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
    output_dir="./germline/snv",
    output_name=NULL,
    indel_candidates=NULL,
    extract_pass=TRUE,
    annotate=TRUE,
    tabulate=TRUE,
    targeted=TRUE,
    verbose=TRUE,
    batch_config=build_default_preprocess_config(),
    mode="local",
    threads=1,
    ram=4,
    time="48:0:0",
    sheet=NULL,
    inherit=NULL,
    select=NULL,
    executor_id=NULL,
    hold=NULL
){
  
    this.envir=environment()
    set_envir_vars(
        envir=this.envir,
        vars="normal"
    )

    run_main=function(
        envir
    ){
        this.envir=environment()
        append_envir(this.envir,envir)
        set_steps_vars(envir=this.envir)

        steps[[fn]]$exec_code=paste0(
          bin_strelka,
          " --bam ",input,
          " --referenceFasta ", ref_genome ,
          " --runDir ", out_file_dir, 
          ifelse(targeted," --exome ",""),
          ifelse(indel_candidates,paste0("--indelCandidates ",indel_candidates),""), "; ",
          paste0(
            out_file_dir,
            "/runWorkflow.py -m local -j ",
            threads)
          )

        steps[[fn]]$out_file$workflow=paste0(out_file_dir,"/runWorkflow.py")
        steps[[fn]]$out_file$stats=list(
              tsv=paste0(out_file_dir,"/results/stats/runStats.tsv"),
              xml=paste0(out_file_dir,"/results/stats/runStats.xml")
        )
        steps[[fn]]$out_file$variants=list(
              indel=paste0(out_file_dir,"/results/variants/somatic.indels.vcf.gz"),
              indel_idx=paste0(out_file_dir,"/results/variants/somatic.indels.vcf.gz.tbi"),
              snv=paste0(out_file_dir,"/results/variants/somatic.snvs.vcf.gz"),
              snv_idx=paste0(out_file_dir,"/results/variants/somatic.snvs.vcf.gz.tbi")
        )
      
        run_job(
              envir=this.envir
        )

        envir$steps <-steps

      }
       envirs=run_envir(envirs=this.envir$envirs)

}