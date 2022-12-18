
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
    output_dir=".",
    output_name=NULL,
    annotate=TRUE,
    tabulate=TRUE,
    targeted=TRUE,
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,
    ram=4,
    mode="local",
    time="48:0:0",
    ss=NULL,
    inherit=NULL,
    select=NULL,
    executor_id=NULL,
    hold=NULL
){



    this.envir=environment()

    

    set_envir_vars(
        envir=this.envir,
        vars=if(tumour,"tumour","normal"),
        ids=output_name,
        executor_id = executor_id,
        dir_name = patient_id
    )

    main=function(
      envir
    ){

      this.envir=environment()
      append_envir(this.envir,envir)
      #### Overide envir vars

      set_steps_vars(envir=this.envir)


      steps[[fn]] <-append(
            steps[[fn]],
            call_sv_manta(
                bin_bcftools=bin_bcftools,
                bin_bgzip=bin_bgzip,
                bin_tabix=bin_tabix,
                bin_manta=bin_manta,
                bin_vep=bin_vep,
                cache_vep=cache_vep,
                tumour=ifelse(tumour,input,NULL),
                normal=ifelse(tumour,normal,input),
                annotate=annotate,
                tabulate=tabulate,
                ref_genome=ref_genome,
                output_dir=out_file_dir,
                targeted=targeted,
                verbose=verbose,
                threads=threads,
                ram=ram,
                executor_id=task_id
        )
    )
        


      steps[[fn]] <-append(
          steps[[fn]],
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
              indel_candidates=steps[[fn]]$call_sv_manta$variants$out_file$indel_candidates,
              tumour=ifelse(tumour,input,NULL),
              normal=ifelse(tumour,normal,input),
              ref_genome=ref_genome,
              output_dir=out_file_dir,
              targeted=targeted,
              verbose=verbose,
              threads=threads,
              ram=ram,
              executor_id=task_id
        )
      )
}

    if(is.null(select)){
        run_self(
            envir=this.envir
        )
    }else{
        set_envir_vars(envir=this.envir)
        run_main(
          envir=this.envir
        )
        return(steps)
    }
    
    
}




#' Parallel call variants per sample Strelka
#'
#' This function wraps the functions for variant calling with Strelka and Manta
#' 
#' @param bin_strelka_somatic Path to strelka somatic workflow binary
#' @param bin_strelka_germline Path to strelka germline workflow binary
#' @param bin_manta Path to manta pipeline binary
#' @param sample_sheet [OPTIONAL] Path to sheet with sample information.
#' @param bam_dir [OPTIONAL] Path to directory with BAM files.
#' @param normal_id [OPTIONAL] Path to directory with BAM files.
#' @param pattern [OPTIONAL] Pattern to search BAM files within directory.
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




multisample_call_variants_strelka=function(
    bin_bcftools=build_default_tool_binary_list()$bin_bcftools,
    bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
    bin_tabix=build_default_tool_binary_list()$bin_tabix,
    bin_vep=build_default_tool_binary_list()$bin_vep,
    cache_vep=build_default_cache_list()$cache_vep,
    bin_strelka_somatic=build_default_tool_binary_list()$bin_strelka$somatic,
    bin_strelka_germline=build_default_tool_binary_list()$bin_strelka$germline,
    bin_manta=build_default_tool_binary_list()$bin_manta,
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    sample_sheet=NULL,
    normal_id=NULL,
    bam_dir="",
    pattern="bam$",
    patient_id="",
    variants="all",
    extract_pass=TRUE,
    annotate=TRUE,
    tabulate=TRUE,
    indel_cnds=TRUE,
    output_dir=".",
    header=TRUE,
    sep="\t",
    targeted=TRUE,verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=4,mode="local",
    executor_id=make_unique_id("multisampleStrelka"),
    task_name="multisampleStrelka",time="48:0:0",
    update_time=60,wait=FALSE,hold=NULL
){

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir)
  job=build_job(executor_id=executor_id,task_id=task_id)

        job_report=build_job_report(
                job_id=job,
                executor_id=executor_id,
                exec_code=list(), 
                task_id=task_id,
                input_args=argg,
                out_file_dir=out_file_dir,
                out_files=list(
                )
        )



    columns=c(
        "bin_bcftools",
        "bin_bgzip",
        "bin_tabix",
        "bin_vep",
        "cache_vep",
        "bin_strelka_somatic",
        "bin_strelka_germline",
        "bin_manta",
        "patient_id",
        "tumour",
        "normal",
        "extract_pass",
        "annotate",
        "tabulate",
        "variants",
        "indel_cnds",
        "ref_genome",
        "output_dir",
        "verbose",
        "targeted",
        "batch_config",
        "threads",
        "ram",
        "mode",
        "time",
        "hold"
    )


    if(!is.null(sample_sheet)){
      
        if(!is.data.frame(sample_sheet)){
                file_info=read.csv(sample_sheet,header=header,sep=sep,stringsAsFactors=FALSE)
                if(!header){
                    names(file_info)=columns
                }
        }else{
                file_info=sample_sheet
        }
        
        file_info=file_info %>% dplyr::group_by(dplyr::across(-tumour)) %>% dplyr::summarise(tumour=list(tumour))

        job_report[["steps"]][["multisample_strelka"]]=parallel::mclapply(seq(1,nrow(file_info)),FUN=function(x){
            
            lapply(columns,FUN=function(col){
                if(is.null(file_info[[col]])){
                    file_info[[col]]<<-get(col)
                }

         
            })
            
        
            job_report<-parallel_samples_call_variants_strelka(   
                bin_bcftools=file_info[x,]$bin_bcftools,
                bin_bgzip=file_info[x,]$bin_bgzip,
                bin_tabix=file_info[x,]$bin_tabix,
                bin_vep=file_info[x,]$bin_vep,
                cache_vep=file_info[x,]$cache_vep,
                bin_strelka_somatic=file_info[x,]$bin_strelka_somatic,
                bin_strelka_germline=file_info[x,]$bin_strelka_germline,
                bin_manta=file_info[x,]$bin_manta,
                patient_id=file_info[x,]$patient_id,
                tumour=file_info[x,]$tumour,
                normal=file_info[x,]$normal,
                extract_pass=file_info[x,]$extract_pass,
                annotate=file_info[x,]$annotate,
                tabulate=file_info[x,]$tabulate,
                variants=file_info[x,]$variants,
                indel_cnds=file_info[x,]$indel_cnds,
                ref_genome=file_info[x,]$ref_genome,
                output_dir=file_info[x,]$output_dir,
                targeted=file_info[x,]$targeted,
                verbose=file_info[x,]$verbose,
                batch_config=file_info[x,]$batch_config,
                threads=file_info[x,]$threads,
                ram=file_info[x,]$ram,mode=file_info[x,]$mode,
                executor_id=task_id,
                time=file_info[x,]$time,
                hold=file_info[[x,"hold"]]
              )
            },mc.cores=ifelse(mode=="local",1,3))

    }else{
        bam_dir_path=system(paste("realpath",bam_dir),intern=TRUE)
        normal=system(paste0("find ",bam_dir_path,"| grep ",pattern),intern=TRUE)
        tumour=bam_files[!grepl(normal_id,bam_files)]
        normal=bam_files[grepl(normal_id,bam_files)]


              job_report<-parallel_samples_call_variants_strelka(
                bin_bcftools=bin_bcftools,
                bin_bgzip=bin_bgzip,
                bin_tabix=bin_tabix,
                bin_vep=bin_vep,
                cache_vep=cache_vep,
                bin_strelka_somatic=bin_strelka_somatic,
                bin_strelka_germline=bin_strelka_germline,
                bin_manta=bin_manta,
                patient_id=patient_id,
                tumour=tumour,
                normal=normal,
                variants=variants,
                extract_pass=extract_pass,
                annotate=annotate,
                tabulate=tabulate,
                indel_cnds=indel_cnds,
                ref_genome=ref_genome,
                output_dir=output_dir,
                targeted=targeted,verbose=verbose,
                batch_config=batch_config,
                threads=threads,ram=ram,mode=mode,
                executor_id=task_id,
                time=time,
                hold=hold
              )
    }


    if(wait&&mode=="batch"){
        job_validator(job=unlist_level(named_list=job_report[["steps"]][["multisample_strelka"]],var="job_id"),
        time=update_time,verbose=verbose,threads=threads)
    }

    return(job_report)

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
    output_dir=".",
    output_name=NULL,
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,
    ram=4,
    mode="local",
    time="48:0:0",
    executor_id=NULL,
    inherit=NULL,
    hold=NULL
){
   this.envir=environment()
    set_envir_vars(
        envir=this.envir,
        vars="tumour",
        ids=output_name,
        executor_id = executor_id,
        dir_name="somatic"
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
            vcf_snv=steps[[fn]]$out_file$variants$snv,
            vcf_indel=steps[[fn]]$out_file$variants$indel,
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
                vcf_sv=steps[[fn]]$extract_pass_variants_strelka_vcf$out_file$variants$sv,
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


    if(is.null(select)){
        run_self(
            envir=this.envir
        )
    }else{
        set_envir_vars(envir=this.envir)
        run_main(
          envir=this.envir
        )
        return(steps)
    }


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
    output_dir=".",
    output_name=NULL,
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,
    ram=4,
    mode="local",
    time="48:0:0",
    executor_id=NULL,
    inherit=NULL,
    hold=NULL
){
   this.envir=environment()
    set_envir_vars(
        envir=this.envir,
        vars="normal",
        ids=output_name,
        executor_id = executor_id,
        dir_name="germline"
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
              sv=paste0(out_file_dir,"/results/variants/somaticSV.vcf.gz"),
              sv_idx=paste0(out_file_dir,"/results/variants/somaticSV.vcf.gz"),
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
            vcf_snv=steps[[fn]]$out_file$variants$snv,
            vcf_indel=steps[[fn]]$out_file$variants$indel,
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
                vcf_sv=steps[[fn]]$extract_pass_variants_strelka_vcf$out_file$variants$sv,
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


    if(is.null(select)){
        run_self(
            envir=this.envir
        )
    }else{
        set_envir_vars(envir=this.envir)
        run_main(
          envir=this.envir
        )
        return(steps)
    }


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
    output_dir=".",
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
    inherit=NULL,
    select=NULL,
    executor_id=NULL,
    hold=NULL
){
  
 this.envir=environment()
    set_envir_vars(
        envir=this.envir,
        vars="tumour",
        ids=output_name,
        executor_id = executor_id,
        dir_name="snv_indel"
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
                output_dir=out_file_dir,
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
                output_dir=out_file_dir,
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

    
  
    if(is.null(select)){
        run_self(
            envir=this.envir
        )
    }else{
        
        set_envir_inputs(envir=this.envir)
        run_main(
          envir=this.envir
        )
        return(steps)
    }

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
    output_dir=".",
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
    inherit=NULL,
    select=NULL,
    executor_id=NULL,
    hold=NULL
){
  
 this.envir=environment()
    set_envir_vars(
        envir=this.envir,
        vars="tumour",
        ids=output_name,
        executor_id = executor_id,
        dir_name="sv"
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
                output_dir=out_file_dir,
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
              call_germline_sv_manta(
                bin_manta=bin_manta,
                bin_strelka=bin_manta_germline,
                bin_bcftools=bin_bcftools,
                bin_bgzip=bin_bgzip,
                bin_tabix=bin_tabix,
                bin_vep=bin_vep,
                cache_vep=cache_vep,
                normal=normal,
                annotate=annotate,
                tabulate=tabulate,
                ref_genome=ref_genome,
                output_dir=out_file_dir,
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

    
  
    if(is.null(select)){
        run_self(
            envir=this.envir
        )
    }else{
        
        set_envir_inputs(envir=this.envir)
        run_main(
          envir=this.envir
        )
        return(steps)
    }

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
    output_dir=".",
    output_name=NULL
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
    inherit=NULL,
    select=NULL,
    executor_id=NULL,
    hold=NULL
){

    this.envir=environment()
    set_envir_vars(
        envir=this.envir,
        vars="tumour",
        ids=output_name,
        executor_id = executor_id,
        dir_name="somatic"
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


        steps[[fn]] <- append(steps[[fn]], add_af_strelka_vcf(
          bin_bgzip=bin_bgzip,
          bin_tabix=bin_tabix,
          vcf_snv=steps[[fn]]$out_file$variants$snv,
          vcf_indel=steps[[fn]]$out_file$variants$indel,
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


    
      if(is.null(select)){
          run_self(
              envir=this.envir
          )
      }else{
          
          set_envir_inputs(envir=this.envir)
          run_main(
            envir=this.envir
          )
          return(steps)
      }


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
    output_dir=".",
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
    inherit=NULL,
    select=NULL,
    executor_id=NULL,
    hold=NULL
){
  
    this.envir=environment()
    set_envir_vars(
        envir=this.envir,
        vars=NULL,
        ids=output_name,
        executor_id = executor_id,
        dir_name="germline"
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
      if(is.null(select)){
          run_self(
              envir=this.envir
          )
      }else{
          
          set_envir_inputs(envir=this.envir)
          run_main(
            envir=this.envir
          )
          return(steps)
      }

}