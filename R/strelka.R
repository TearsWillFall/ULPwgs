

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
    bin_strelka_somatic=build_default_tool_binary_list()$bin_strelka$somatic,
    bin_strelka_germline=build_default_tool_binary_list()$bin_strelka$germline,
    bin_manta=build_default_tool_binary_list()$bin_manta,
    tumour="",normal="",variants="all",indel_cnds=TRUE,
    ref_genome="",output_dir=".",
    targeted=TRUE,verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=4,mode="local",
    executor_id=make_unique_id("callSVManta"),
    task_name="callSVManta",time="48:0:0",
    update_time=60,wait=FALSE,hold=""
){

    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir,name="strelka_reports")
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
            bin_manta=bin_manta,
            tumour=tumour,normal=normal,
            ref_genome=ref_genome,output_dir=out_file_dir,
            targeted=targeted,verbose=verbose,
            batch_config=batch_config,threads=threads,
            ram=ram,mode=mode,executor_id=task_id,
            time=time,hold=hold
        )
        ### Call INDELs before calling SNVs if indel cnds are given
        if(indel_cnds){
          hold=jobs_report[["steps"]][["callSVManta"]]$job_id
          indel_candidates=jobs_report[["steps"]][["callSVManta"]]$out_files$indel_candidates
          indel_candidates=paste0(" --indelCandidates ",indel_candidates)
        }
   }


    if(variants=="snv"|variants=="all"){
        jobs_report[["steps"]][["callSNVStrelka"]]=call_snvs_strelka(
          bin_strelka_somatic=bin_strelka_somatic,
          bin_strelka_germline=bin_strelka_germline,
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
    bin_manta=build_default_tool_binary_list()$bin_manta,
    tumour="",normal="",ref_genome=build_default_reference_list()$HG19$reference$genome,
    output_dir=".",targeted=TRUE,verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=4,mode="local",
    executor_id=make_unique_id("callSVManta"),
    task_name="callSVManta",time="48:0:0",
    update_time=60,wait=FALSE,hold=""
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
      somatic_sv=paste0(out_file_dir,"MantaWorkflow/results/variants/somaticSV.vcf.gz")
      somatic_sv_index=paste0(out_file_dir,"MantaWorkflow/results/variants/somaticSV.vcf.gz.tbi")
    }else if(normal!=""){
      out_file_dir=set_dir(dir=output_dir,name="normal_sv")
      normal_input=paste0(" --normalBam ", paste0(normal,collapse=" --normalBam "))
    
    }

    exome=""
    if (targeted){
        exome=" --exome "
    } 

  
    job_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    task_id=task_id,
    input_args = argg,
    out_file_dir=out_file_dir,
    out_files=list(
      workflow=paste0(out_file_dir,"/MantaWorkflow/runWorkflow.py"),
      stats=list(
        aligment=paste0(out_file_dir,"/MantaWorkflow/results/stats/alignmentStatsSummary.txt"),
        tsv=paste0(out_file_dir,"/MantaWorkflow/results/stats/svCandidateGenerationStats.tsv"),
        xml=paste0(out_file_dir,"/MantaWorkflow/results/stats/svCandidateGenerationStats.xml"),
        graph=paste0(out_file_dir,"/MantaWorkflow/results/stats/svLocusGraphStats.tsv")
      ),
      variants=list(
        small_indel_candidate=paste0(out_file_dir,"/MantaWorkflow/results/variants/candidateSmallIndels.vcf.gz"),
        small_indel_candidate_index=paste0(out_file_dir,"/MantaWorkflow/results/variants/candidateSmallIndels.vcf.gz.tbi"),
        sv_candidate=paste0(out_file_dir,"/MantaWorkflow/results/variants/candidateSV.vcf.gz"),
        sv_candidate_index=paste0(out_file_dir,"/MantaWorkflow/results/variants/candidateSV.vcf.gz.tbi"),
        diploid_sv=paste0(out_file_dir,"/MantaWorkflow/results/variants/diploidSV.vcf.gz"),
        diploid_sv_index=paste0(out_file_dir,"/MantaWorkflow/results/variants/diploidSV.vcf.gz.tbi"),
        somatic_sv=somatic_sv,
        somatic_sv_index=somatic_sv_index
        )
      )
    )

  
    exec_code=paste0(bin_manta,tumour_input,normal_input," --referenceFasta ", 
    ref_genome ," --runDir ", out_file_dir, exome, "; ",
    paste0(out_file_dir,"/MantaWorkflow/runWorkflow.py -m local -j ",threads))

    if(mode=="batch"){
        hold=unlist_lvl(jobs_report[["steps"]],var="job_id")
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
    bin_strelka_somatic=build_default_tool_binary_list()$bin_strelka$somatic,
    bin_strelka_germline=build_default_tool_binary_list()$bin_strelka$germline,
    tumour="",normal="",ref_genome=build_default_reference_list()$HG19$reference$genome,
    output_dir=".",indel_candidates="",
    targeted=TRUE,verbose=TRUE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=4,mode="local",
    executor_id=make_unique_id("callSNVStrelka"),
    task_name="callSNVStrelka",time="48:0:0",
    update_time=60,wait=FALSE,hold=""
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
            bin_strelka=bin_strelka_germline,
            normal=normal,
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
    bin_strelka=build_default_tool_binary_list()$bin_strelka$somatic,
    tumour="",normal="",ref_genome=build_default_reference_list()$HG19$reference$genome,
    output_dir=".",indel_candidates="",
    targeted=TRUE,verbose=TRUE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=4,mode="local",
    executor_id=make_unique_id("callSomaticSNVStrelka"),
    task_name="callSomaticSNVStrelka",time="48:0:0",
    update_time=60,wait=FALSE,hold=""
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
      workflow=paste0(out_file_dir,"/StrelkaSomaticWorkflow/runWorkflow.py"),
      stats=list(
        tsv=paste0(out_file_dir,"/StrelkaSomaticWorkflow/results/stats/runStats.tsv"),
        xml=paste0(out_file_dir,"/StrelkaSomaticWorkflow/results/stats/runStats.tsv")
    ),
    variants=list(
        indel=paste0(out_file_dir,"/StrelkaSomaticWorkflow/results/variants/somatic.indels.vcf.gz"),
        indel_index=paste0(out_file_dir,"/StrelkaSomaticWorkflow/results/variants/somatic.indels.vcf.gz.tbi"),
        snvs=paste0(out_file_dir,"/StrelkaSomaticWorkflow/results/variants/somatic.snvs.vcf.gz"),
        snvs_index=paste0(out_file_dir,"/StrelkaSomaticWorkflow/results/variants/somatic.snvs.vcf.gz.tbi")
    )
    )
  )

    exome=""
    if (targeted){
        exome=" --exome "
    }

    exec_code=paste0(bin_strelka,tumour_input,normal_input," --referenceFasta ", 
    ref_genome ," --runDir ", out_file_dir, exome, indel_candidates, "; ",
    paste0(out_file_dir,"/StrelkaSomaticWorkflow/runWorkflow.py -m local -j ",threads))

    if(mode=="batch"){
        hold=unlist_lvl(jobs_report[["steps"]],var="job_id")
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

    if(wait&&mode=="batch"){
        job_validator(job=jobs_report$job_id,time=update_time,
        verbose=verbose,threads=threads)
    }

    return(jobs_report)



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
    bin_strelka=build_default_tool_binary_list()$bin_strelka$germline,
    normal="",ref_genome=build_default_reference_list()$HG19$reference$genome,
    output_dir=".",indel_candidates="",
    targeted=TRUE,verbose=TRUE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=4,mode="local",
    executor_id=make_unique_id("callGermlineSNVStrelka"),
    task_name="callGermlineSNVStrelka",time="48:0:0",
    update_time=60,wait=FALSE,hold=""
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
    
    jobs_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    task_id=task_id,
    input_args = argg,
    out_file_dir=out_file_dir,
    out_files=list(
      workflow=paste0(out_file_dir,"/StrelkaGermlineWorkflow/runWorkflow.py"),
      stats=list(
        tsv=paste0(out_file_dir,"/StrelkaGermlineWorkflow/results/stats/runStats.tsv"),
        xml=paste0(out_file_dir,"/StrelkaGermlineWorkflow/results/stats/runStats.tsv")
    ),
    variants=list(
        indel=paste0(out_file_dir,"/StrelkaGermlineWorkflow/results/variants/somatic.indels.vcf.gz"),
        indel_index=paste0(out_file_dir,"/StrelkaGermlineWorkflow/results/variants/somatic.indels.vcf.gz.tbi"),
        snvs=paste0(out_file_dir,"/StrelkaGermlineWorkflow/results/variants/somatic.snvs.vcf.gz"),
        snvs_index=paste0(out_file_dir,"/StrelkaGermlineWorkflow/results/variants/somatic.snvs.vcf.gz.tbi")
    )
    )
  )

    exome=""
    if (targeted){
        exome=" --exome "
    }

    exec_code=paste0(bin_strelka,normal_input," --referenceFasta ", 
    ref_genome ," --runDir ", out_file_dir, exome, indel_candidates, "; ",
    paste0(out_file_dir,"/StrelkaGermlineWorkflow/runWorkflow.py -m local -j ",threads))

    if(mode=="batch"){
        hold=unlist_lvl(jobs_report[["steps"]],var="job_id")
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

    if(wait&&mode=="batch"){
        job_validator(job=jobs_report$job_id,time=update_time,
        verbose=verbose,threads=threads)
    }

    return(jobs_report)

}