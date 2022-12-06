

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
    rdata=NULL,
    selected=NULL,
    bin_bcftools=build_default_tool_binary_list()$bin_bcftools,
    bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
    bin_tabix=build_default_tool_binary_list()$bin_tabix,
    bin_vep=build_default_tool_binary_list()$bin_vep,
    cache_vep=build_default_cache_list()$cache_vep,
    bin_strelka_somatic=build_default_tool_binary_list()$bin_strelka$somatic,
    bin_strelka_germline=build_default_tool_binary_list()$bin_strelka$germline,
    bin_manta=build_default_tool_binary_list()$bin_manta,
    tumour="",normal="",variants="all",indel_cnds=TRUE,
    ref_genome="",output_dir=".",
    output_name="",
    extract_pass=TRUE,
    annotate=TRUE,
    tabulate=TRUE,
    targeted=TRUE,verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=4,mode="local",
    executor_id=make_unique_id("callSVManta"),
    task_name="callSVManta",time="48:0:0",
    update_time=60,wait=FALSE,hold=NULL
){

  if(!is.null(rdata)){
    load(rdata)
    if(!is.null(selected)){
      tumour=tumour_list[selected]
    }
  }


    id=""
    if(output_name!=""){
      id=output_name
    }else{
      id=get_file_name(tumour)
    }



    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir,name=paste0("strelka_reports/",id))
    job=build_job(executor_id=executor_id,task_id=task_id)


    jobs_report=build_job_report(
      job_id=job,
      executor_id=executor_id,
      task_id=task_id,
      input_args = argg,
      out_file_dir=out_file_dir,
      out_files=list()
    )

    

    if(variants=="sv"|variants=="all"){
         jobs_report[["steps"]][["callSVManta"]]=call_sv_manta(
            bin_bcftools=bin_bcftools,
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            bin_vep=bin_vep,
            cache_vep=cache_vep,
            bin_manta=bin_manta,
            tumour=tumour,normal=normal,
            extract_pass=extract_pass,
            annotate=annotate,
            tabulate=tabulate,
            ref_genome=ref_genome,output_dir=out_file_dir,
            targeted=targeted,verbose=verbose,
            batch_config=batch_config,threads=threads,
            ram=ram,mode=mode,executor_id=task_id,
            time=time,hold=hold
        )
        ### Call INDELs before calling SNVs if indel cnds are given
        if(indel_cnds){
          hold=jobs_report[["steps"]][["callSVManta"]]$job_id
          indel_candidates=jobs_report[["steps"]][["callSVManta"]]$out_files$vcf_small_indel_candidate
          
        }
   }


    if(variants=="snv"|variants=="all"){
        jobs_report[["steps"]][["callSNVStrelka"]]=call_snvs_strelka(
          bin_bcftools=bin_bcftools,
          bin_bgzip=bin_bgzip,
          bin_tabix=bin_tabix,
          bin_vep=bin_vep,
          cache_vep=cache_vep,
          extract_pass=extract_pass,
          annotate=annotate,
          tabulate=tabulate,
          bin_strelka_somatic=bin_strelka_somatic,
          bin_strelka_germline=bin_strelka_germline,
          indel_candidates=indel_candidates,
          tumour=tumour,normal=normal,
          ref_genome=ref_genome,output_dir=out_file_dir,
          targeted=targeted,verbose=verbose,
          batch_config=batch_config,threads=threads,
          ram=ram,mode=mode,executor_id=task_id,
          time=time,hold=hold
      )
    }
    return(jobs_report)
}





#' Parallel call variants per sample Strelka
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



parallel_samples_call_variants_strelka=function(
    bin_bcftools=build_default_tool_binary_list()$bin_bcftools,
    bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
    bin_tabix=build_default_tool_binary_list()$bin_tabix,
    bin_vep=build_default_tool_binary_list()$bin_vep,
    cache_vep=build_default_cache_list()$cache_vep,
    bin_strelka_somatic=build_default_tool_binary_list()$bin_strelka$somatic,
    bin_strelka_germline=build_default_tool_binary_list()$bin_strelka$germline,
    bin_manta=build_default_tool_binary_list()$bin_manta,
    patient_id=NA,tumour=NA,normal=NA,variants="all",indel_cnds=TRUE,
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    output_dir=".",
    extract_pass=TRUE,
    annotate=TRUE,
    tabulate=TRUE,
    targeted=TRUE,verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=4,mode="local",
    executor_id=make_unique_id("parSampleCallSVManta"),
    task_name="parSampleCallSVManta",time="48:0:0",
    update_time=60,wait=FALSE,hold=NULL
){

    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir,name=patient_id)
    job=build_job(executor_id=executor_id,task_id=task_id)


    
    jobs_report=build_job_report(
      job_id=job,
      executor_id=list(),
      task_id=task_id,
      input_args = argg,
      out_file_dir=out_file_dir,
      out_files=list()
    )


    tumour_list=tumour
    names(tumour_list)=Vectorize(get_file_name)(tumour_list)
 
    
    if(mode=="local"){
      jobs_report[["steps"]][["par_sample_call_variants"]]<-
      lapply(tumour_list,FUN=function(tumour){
        job_report <- call_variants_strelka(
            bin_bcftools=bin_bcftools,
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            bin_vep=bin_vep,
            cache_vep=cache_vep,
            bin_strelka_somatic=bin_strelka_somatic,
            bin_strelka_germline=bin_strelka_germline,
            bin_manta=bin_manta,
            extract_pass=extract_pass,
            annotate=annotate,
            tabulate=tabulate,
            tumour=tumour,normal=normal,
            variants=variants,indel_cnds=indel_cnds,
            ref_genome=ref_genome,output_dir=out_file_dir,
            targeted=targeted,verbose=verbose,
            batch_config=batch_config,
            threads=threads,ram=ram,mode=mode,
            executor_id=task_id,
            time=time,
            hold=hold
        )
      })
    }else if(mode=="batch"){
            rdata_file=paste0(tmp_dir,"/",job,".samples.RData")
            output_dir=out_file_dir
            executor_id=task_id
            save(
              tumour_list,
              bin_bcftools,
              bin_bgzip,
              bin_tabix,
              bin_vep,
              cache_vep,
              bin_strelka_somatic,
              bin_strelka_germline,
              bin_manta,
              variants,
              extract_pass,
              annotate,
              tabulate,
              indel_cnds,
              ref_genome,
              targeted,
              verbose
            )
            exec_code=paste0("Rscript -e \"ULPwgs::parallel_samples_call_variants_strelka(rdata=\\\"",
            rdata_file,"\\\",selected=$SGE_TASK_ID)\"")
            out_file_dir2=set_dir(dir=out_file_dir,name="batch")
            batch_code=build_job_exec(job=job,time=time,ram=ram,
            threads=threads,output_dir=out_file_dir2,
            hold=hold,array=length(tumour_list))
            exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)

            if(verbose){
                print_verbose(job=job,arg=argg,exec_code=exec_code)
            }
            error=execute_job(exec_code=exec_code)
            if(error!=0){
                stop("strelka failed to run due to unknown error.
                Check std error for more information.")
            }

    
          jobs_report[["steps"]][["par_sample_call_variants"]]<- build_job_report(
                job_id=job,
                executor_id=executor_id,
                exec_code=exec_code, 
                task_id=task_id,
                input_args=argg,
                out_file_dir=out_file_dir,
                out_files=list(
                    somatic=list(
                      snvs=list(
                             workflow=paste0(out_file_dir,"/",patient_id,
                          "/strelka_reports/",names(tumour_list),"/somatic_snvs/runWorkflow.py"),
                              stats=list(
                                tsv=paste0(out_file_dir,"/",patient_id,
                          "/strelka_reports/",names(tumour_list),"/somatic_snvs/results/stats/runStats.tsv"),
                                xml=paste0(out_file_dir,"/",patient_id,
                          "/strelka_reports/",names(tumour_list),"/somatic_snvs/results/stats/runStats.tsv")
                            ),
                            variants=list(
                                indel=paste0(out_file_dir,"/",patient_id,
                          "/strelka_reports/",names(tumour_list),"/somatic_snvs/results/variants/somatic.indels.vcf.gz"),
                                indel_index=paste0(out_file_dir,"/",patient_id,
                          "/strelka_reports/",names(tumour_list),"/somatic_snvs/results/variants/somatic.indels.vcf.gz.tbi"),
                                snvs=paste0(out_file_dir,"/",patient_id,
                          "/strelka_reports/",names(tumour_list),"/somatic_snvs/results/variants/somatic.snvs.vcf.gz"),
                                snvs_index=paste0(out_file_dir,"/",patient_id,
                          "/strelka_reports/",,names(tumour_list),,"/somatic_snvs/results/variants/somatic.snvs.vcf.gz.tbi")
                            )
                    ),
                    svs=list(
                        workflow=paste0(out_file_dir,"/",patient_id,
                          "/strelka_reports/",names(tumour_list),
                          "/somatic_sv/runWorkflow.py"),
                        stats=list(
                          aligment=paste0(out_file_dir,"/",patient_id,
                          "/strelka_reports/",names(tumour_list),
                          "/somatic_sv/results/stats/alignmentStatsSummary.txt"),
                          tsv=paste0(out_file_dir,"/",patient_id,
                          "/strelka_reports/",names(tumour_list),
                          "/somatic_sv/results/stats/svCandidateGenerationStats.tsv"),
                          xml=paste0(out_file_dir,"/",patient_id,
                          "/strelka_reports/",names(tumour_list),
                          "/somatic_sv/results/stats/svCandidateGenerationStats.xml"),
                          graph=paste0(out_file_dir,"/",patient_id,
                          "/strelka_reports/",names(tumour_list),
                          "/somatic_sv/results/stats/svLocusGraphStats.tsv")
                        ),
                        variants=list(
                          small_indel_candidate=paste0(out_file_dir,"/",patient_id,
                          "/strelka_reports/",names(tumour_list),
                          "/somatic_sv/results/variants/candidateSmallIndels.vcf.gz"),
                          small_indel_candidate_index=paste0(out_file_dir,"/",patient_id,
                          "/strelka_reports/",names(tumour_list),
                          "/somatic_sv/results/variants/candidateSmallIndels.vcf.gz.tbi"),
                          sv_candidate=paste0(out_file_dir,"/",patient_id,
                          "/strelka_reports/",names(tumour_list),
                          "/somatic_sv/results/variants/candidateSV.vcf.gz"),
                          sv_candidate_index=paste0(out_file_dir,"/",
                          patient_id,"/",names(tumour_list),
                          "/strelka_reports/",names(tumour_list),
                          "/somatic_sv/results/variants/candidateSV.vcf.gz.tbi"),
                          diploid_sv=paste0(out_file_dir,"/",patient_id,
                          "/strelka_reports/",names(tumour_list),
                          "/somatic_sv/results/variants/diploidSV.vcf.gz"),
                          diploid_sv_index=paste0(out_file_dir,"/",patient_id,
                          "/strelka_reports/",names(tumour_list),
                          "/somatic_sv/results/variants/diploidSV.vcf.gz.tbi"),
                          somatic_sv=paste0(out_file_dir,"/",patient_id,
                          "/strelka_reports/",names(tumour_list),
                          "/somatic_svs/results/variants/somaticSV.vcf.gz"),
                          somatic_sv_index=paste0(out_file_dir,"/",patient_id,
                          "/strelka_reports/",names(tumour_list),
                          "/somatic_svs/results/variants/somaticSV.vcf.gz.tbi")
                    )
                 )
              )
            )
          )
        
      }

    return(jobs_report)
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



call_sv_manta=function(
    bin_bcftools=build_default_tool_binary_list()$bin_bcftools,
    bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
    bin_tabix=build_default_tool_binary_list()$bin_tabix,
    bin_vep=build_default_tool_binary_list()$bin_vep,
    cache_vep=build_default_cache_list()$cache_vep,
    bin_manta=build_default_tool_binary_list()$bin_manta,
    tumour="",normal="",
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    extract_pass=TRUE,
    annotate=TRUE,
    tabulate=TRUE,
    output_dir=".",targeted=TRUE,verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=4,mode="local",
    executor_id=make_unique_id("callSVManta"),
    task_name="callSVManta",time="48:0:0",
    update_time=60,wait=FALSE,hold=NULL
){
  
    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    job=build_job(executor_id=executor_id,task_id=task_id)

    somatic_sv=NULL
    soamtic_sv_index=NULL
    if(tumour!=""){
      out_file_dir=set_dir(dir=output_dir,name="somatic_sv")
      tumour_input=paste(" --tumorBam ",tumour)
      normal_input=paste(" --normalBam ",normal)
      somatic_sv=paste0(out_file_dir,"/results/variants/somaticSV.vcf.gz")
      somatic_sv_index=paste0(out_file_dir,"/results/variants/somaticSV.vcf.gz.tbi")
    }else if(normal!=""){
      out_file_dir=set_dir(dir=output_dir,name="normal_sv")
      normal_input=paste0(" --normalBam ", paste0(normal,collapse=" --normalBam "))
    
    }

    exome=""
    if (targeted){
        exome=" --exome "
    } 

  
    jobs_report=build_job_report(
      job_id=job,
      executor_id=executor_id,
      task_id=task_id,
      input_args = argg,
      out_file_dir=out_file_dir,
      out_files=list(
        workflow=paste0(out_file_dir,"/runWorkflow.py"),
        stats=list(
          aligment=paste0(out_file_dir,"/results/stats/alignmentStatsSummary.txt"),
          tsv=paste0(out_file_dir,"/results/stats/svCandidateGenerationStats.tsv"),
          xml=paste0(out_file_dir,"/results/stats/svCandidateGenerationStats.xml"),
          graph=paste0(out_file_dir,"/results/stats/svLocusGraphStats.tsv")
        ),
        variants=list(
          vcf_small_indel_candidate=paste0(out_file_dir,"/results/variants/candidateSmallIndels.vcf.gz"),
          small_indel_candidate_index=paste0(out_file_dir,"/results/variants/candidateSmallIndels.vcf.gz.tbi"),
          vcf_sv_candidate=paste0(out_file_dir,"/results/variants/candidateSV.vcf.gz"),
          sv_candidate_index=paste0(out_file_dir,"/results/variants/candidateSV.vcf.gz.tbi"),
          vcf_diploid_sv=paste0(out_file_dir,"/results/variants/diploidSV.vcf.gz"),
          diploid_sv_index=paste0(out_file_dir,"/results/variants/diploidSV.vcf.gz.tbi"),
          vcf_somatic_sv=somatic_sv,
          somatic_sv_index=somatic_sv_index
          )
        )
      )

  
    exec_code=paste0(bin_manta,tumour_input,normal_input," --referenceFasta ", 
    ref_genome ," --runDir ", out_file_dir, exome, "; ",
    paste0(out_file_dir,"/runWorkflow.py -m local -j ",threads))

    if(mode=="batch"){
        hold=unlist_lvl(job_report[["steps"]],var="job_id")
        out_file_dir2=set_dir(dir=out_file_dir,name="batch")
        batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
        threads=threads,output_dir=out_file_dir2)
        exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
    } 

    if(verbose){
        print_verbose(job=job,arg=argg,exec_code=exec_code)
    }

    error=execute_job(exec_code=exec_code)
    
    if(error!=0){
        stop("manta failed to run due to unknown error.
        Check std error for more information.")
    }


     jobs_report[["steps"]][["addAFStrelka"]] <- add_af_strelka_vcf(
      bin_bgzip=bin_bgzip,
      bin_tabix=bin_tabix,
      vcf_sv=jobs_report$out_files$variants$vcf_somatic_sv,
      verbose=verbose, 
      executor_id=task_id,
      batch_config=batch_config,
      mode=mode,time=time,
      threads=1,ram=1,
      hold=job
    )
    
    vcf_sv=jobs_report[["steps"]][["addAFStrelka"]]$out_files$vcf_sv

    if(extract_pass){
        jobs_report[["steps"]][["extractPASS"]]<-extract_pass_variants_strelks_vcf(
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            vcf_sv=vcf_sv,
            output_dir=out_file_dir,
            verbose=verbose,
            batch_config=batch_config,
            threads=1,ram=1,mode=mode,
            executor_id=task_id,
            time=time,
            hold=unlist_lvl(jobs_report[["steps"]][["addAFStrelka"]],var="job_id")
        )
        vcf_sv=jobs_report[["steps"]][["extractPASS"]]$out_files$vcf_sv
    
    }

  if(annotate){

      jobs_report[["steps"]][["annotateVEP"]]<-annotate_strelka_vep(
          bin_vep=bin_vep,
          bin_bgzip=bin_bgzip,
          bin_tabix=bin_tabix,
          cache_vep=cache_vep,
          vcf_sv=vcf_sv,
          extract_pass=extract_pass,
          output_dir=out_file_dir,
          verbose=verbose,
          batch_config=batch_config,
          threads=threads,ram=ram,mode=mode,
          executor_id=task_id,
          time=time,
          hold=hold
      )
      
  }


    if(wait&&mode=="batch"){
        job_validator(job=job_report$job_id,time=update_time,
        verbose=verbose,threads=threads)
    }

    return(jobs_report)



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
    tumour="",normal="",
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    output_dir=".",indel_candidates="",
    targeted=TRUE,verbose=TRUE,
    extract_pass=TRUE,annotate=TRUE,
    tabulate=TRUE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=4,mode="local",
    executor_id=make_unique_id("callSNVStrelka"),
    task_name="callSNVStrelka",time="48:0:0",
    update_time=60,wait=FALSE,hold=NULL
){
  
    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir)
    job=build_job(executor_id=executor_id,task_id=task_id)

    jobs_report=build_job_report(
      job_id=job,
      executor_id=executor_id,
      task_id=task_id,
      input_args = argg,
      out_file_dir=out_file_dir,
      out_files=list()
    )

    

    if(tumour!=""&normal!=""){
       jobs_report[["steps"]][["callSomaticSNVStrelka"]]<-call_somatic_snvs_strelka(
            bin_bcftools=bin_bcftools,
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            bin_vep=bin_vep,
            cache_vep=cache_vep,
            extract_pass=extract_pass,
            annotate=annotate,
            tabulate=tabulate,
            bin_strelka=bin_strelka_somatic,
            tumour=tumour,normal=normal,
            ref_genome=ref_genome,
            output_dir=out_file_dir,indel_candidates=indel_candidates,
            targeted=targeted,verbose=verbose,
            batch_config=batch_config,
            threads=threads,ram=ram,mode=mode,
            executor_id=task_id,
            hold=hold

      )
    }else if(normal!=""){
         jobs_report[["steps"]][["callGermlineSNVStrelka"]]<-call_germline_snvs_strelka(
             bin_bcftools=bin_bcftools,
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            bin_vep=bin_vep,
            cache_vep=cache_vep,
            bin_strelka=bin_strelka_germline,
            normal=normal,
            extract_pass=extract_pass,
            annotate=annotate,
            tabulate=tabulate,
            ref_genome=ref_genome,
            output_dir=out_file_dir,
            indel_candidates=indel_candidates,
            targeted=targeted,verbose=verbose,
            batch_config=batch_config,
            threads=threads,ram=ram,mode=mode,
            executor_id=task_id,
            hold=hold
      )
    }else {
      stop("Missing tumour/normal BAM file/s.")
    }
  return(jobs_report)

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
    tumour="",normal="",
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    output_dir=".",indel_candidates="",
    extract_pass=TRUE,
    annotate=TRUE,
    tabulate=TRUE,
    targeted=TRUE,verbose=TRUE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=4,mode="local",
    executor_id=make_unique_id("callSomaticSNVStrelka"),
    task_name="callSomaticSNVStrelka",time="48:0:0",
    update_time=60,wait=FALSE,hold=NULL
){

    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir,name="somatic_snvs")
    job=build_job(executor_id=executor_id,task_id=task_id)


    if(tumour!=""&normal!=""){
      tumour_input=paste(" --tumorBam ",tumour)
      normal_input=paste(" --normalBam ",normal)
    }else {
      stop("Missing tumour/normal BAM file/s.")
    }
    

    jobs_report=build_job_report(
      job_id=job,
      executor_id=executor_id,
      task_id=task_id,
      input_args = argg,
      out_file_dir=out_file_dir,
      out_files=list(
        workflow=paste0(out_file_dir,"/runWorkflow.py"),
        stats=list(
          tsv=paste0(out_file_dir,"/results/stats/runStats.tsv"),
          xml=paste0(out_file_dir,"/results/stats/runStats.tsv")
      ),
      variants=list(
          vcf_indel=paste0(out_file_dir,"/results/variants/somatic.indels.vcf.gz"),
          indel_index=paste0(out_file_dir,"/results/variants/somatic.indels.vcf.gz.tbi"),
          vcf_snv=paste0(out_file_dir,"/results/variants/somatic.snvs.vcf.gz"),
          snv_index=paste0(out_file_dir,"/results/variants/somatic.snvs.vcf.gz.tbi")
      )
      ) 
    )

    exome=""
    if (targeted){
        exome=" --exome "
    }

    exec_code=paste0(bin_strelka,tumour_input,normal_input," --referenceFasta ", 
    ref_genome ," --runDir ", out_file_dir, exome, indel_candidates, "; ",
    paste0(out_file_dir,"/runWorkflow.py -m local -j ",threads))

    if(mode=="batch"){
        hold=unlist_lvl(job_report[["steps"]],var="job_id")
        out_file_dir2=set_dir(dir=out_file_dir,name="batch")
        batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
        threads=threads,output_dir=out_file_dir2)
        exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
    } 

    jobs_report$exec_code=exec_code

    if(verbose){
        print_verbose(job=job,arg=argg,exec_code=exec_code)
    }

    error=execute_job(exec_code=exec_code)



    if(error!=0){
        stop("gatk failed to run due to unknown error.
        Check std error for more information.")
    }

      
    jobs_report[["steps"]][["addAFStrelka"]] <- add_af_strelka_vcf(
      bin_bgzip=bin_bgzip,
      bin_tabix=bin_tabix,
      vcf_snv=jobs_report$out_files$variants$vcf_snv,
      vcf_indel=jobs_report$out_files$variants$vcf_indel,
      verbose=verbose, 
      executor_id=task_id,
      batch_config=batch_config,
      mode=mode,time=time,
      threads=1,ram=1,
      hold=job
    )
    
    vcf_snv=jobs_report[["steps"]][["addAFStrelka"]]$out_files$vcf_snv
    vcf_indel=jobs_report[["steps"]][["addAFStrelka"]]$out_files$vcf_indel

    if(extract_pass){
        jobs_report[["steps"]][["extractPASS"]]<-extract_pass_variants_strelks_vcf(
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            vcf_snv=vcf_snv,
            vcf_indel=vcf_indel,
            output_dir=out_file_dir,
            verbose=verbose,
            batch_config=batch_config,
            threads=1,ram=1,mode=mode,
            executor_id=task_id,
            time=time,
            hold=unlist_lvl(jobs_report[["steps"]][["addAFStrelka"]],var="job_id")
        )
        vcf_snv=jobs_report[["steps"]][["extractPASS"]]$out_files$vcf_snv
        vcf_indel=jobs_report[["steps"]][["extractPASS"]]$out_files$vcf_indel
       

    }

  if(annotate){

      jobs_report[["steps"]][["annotateVEP"]]<-annotate_strelka_vep(
          bin_vep=bin_vep,
          bin_bgzip=bin_bgzip,
          bin_tabix=bin_tabix,
          cache_vep=cache_vep,
          vcf_snv=vcf_snv,
          vcf_indel=vcf_indel,
          extract_pass=extract_pass,
          output_dir=out_file_dir,
          verbose=verbose,
          batch_config=batch_config,
          threads=threads,ram=ram,mode=mode,
          executor_id=task_id,
          time=time,
          hold=hold
      )
      
  }


    if(wait&&mode=="batch"){
        job_validator(job=job_report$job_id,time=update_time,
        verbose=verbose,threads=threads)
    }

    return(job_report)



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
    normal="",ref_genome=build_default_reference_list()$HG19$reference$genome,
    output_dir=".",indel_candidates="",
    extract_pass=TRUE,annotate=TRUE,
    tabulate=TRUE,targeted=TRUE,verbose=TRUE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=4,mode="local",
    executor_id=make_unique_id("callGermlineSNVStrelka"),
    task_name="callGermlineSNVStrelka",time="48:0:0",
    update_time=60,wait=FALSE,hold=NULL
){
    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir,name="germline_snvs")
    job=build_job(executor_id=executor_id,task_id=task_id)

    if(normal!=""){
      normal_input=paste(" --bam ",normal)
    }else {
      stop("Missing normal BAM file.")
    }
    
    job_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    task_id=task_id,
    input_args = argg,
    out_file_dir=out_file_dir,
    out_files=list(
      workflow=paste0(out_file_dir,"/runWorkflow.py"),
      stats=list(
        tsv=paste0(out_file_dir,"/results/stats/runStats.tsv"),
        xml=paste0(out_file_dir,"/results/stats/runStats.tsv")
    ),
    variants=list(
        indel=paste0(out_file_dir,"/results/variants/somatic.indels.vcf.gz"),
        indel_index=paste0(out_file_dir,"/results/variants/somatic.indels.vcf.gz.tbi"),
        snvs=paste0(out_file_dir,"/results/variants/somatic.snvs.vcf.gz"),
        snvs_index=paste0(out_file_dir,"/results/variants/somatic.snvs.vcf.gz.tbi")
    )
    )
  )

    exome=""
    if (targeted){
        exome=" --exome "
    }

    exec_code=paste0(bin_strelka,normal_input," --referenceFasta ", 
    ref_genome ," --runDir ", out_file_dir, exome, indel_candidates, "; ",
    paste0(out_file_dir,"/runWorkflow.py -m local -j ",threads))

    if(mode=="batch"){
        hold=unlist_lvl(job_report[["steps"]],var="job_id")
        out_file_dir2=set_dir(dir=out_file_dir,name="batch")
        batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
        threads=threads,output_dir=out_file_dir2)
        exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
    } 

    job_report$exec_code=exec_code

    if(verbose){
        print_verbose(job=job,arg=argg,exec_code=exec_code)
    }

    error=execute_job(exec_code=exec_code)
    
    if(error!=0){
        stop("gatk failed to run due to unknown error.
        Check std error for more information.")
    }

    if(wait&&mode=="batch"){
        job_validator(job=job_report$job_id,time=update_time,
        verbose=verbose,threads=threads)
    }

    return(job_report)

}