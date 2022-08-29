

call_somatic_sv_strelka=function(
    bin_strelka=build_default_tool_binary_list()$bin_strelka,
    bin_manta=build_default_tool_binary_list()$bin_manta,
    bin_bcftools=build_default_tool_binary_list()$bin_bcftools,
    bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
    bin_tabix=build_default_tool_binary_list()$bin_tabix,
    tumor="",normal="",
    ref_genome="",
    output_dir=".",output_name="",
    verbose=FALSE,targeted=FALSE,
    threads=3,exec_options="local"
){

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir,name="somatic_sv")
  job=build_job(executor_id=executor_id,task_id=task_id)


  id=""
  if(output_name!=""){
    id=output_name
  }else if(tumour!=""){
    tumor_input=paste(" --tumorBam ",tumor)
    normal_input=paste(" --normalBam ",normal)
    id=get_file_name(tumour)
  }else if(normal!=""){
    id=get_file_name(normal)
    normal_input=paste0(" --bam ", paste0(normal,collapse=" --bam "))
  }

    jobs_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    task_id=task_id,
    input_args = argg,
    out_file_dir=out_file_dir,
    out_files=list(
      workflow=paste0(out_file_dir,"StrelkaSomaticWorkflow/runWorkflow.py"),
      stats=list(
        tsv=paste0(out_file_dir,"StrelkaSomaticWorkflow/results/stats/runStats.tsv"),
        xml=paste0(out_file_dir,"StrelkaSomaticWorkflow/results/stats/runStats.tsv")
    ),
    variants=list(
        indel=paste0(out_file_dir,"StrelkaSomaticWorkflow/results/variants/somatic.indels.vcf.gz"),
        indel_index=paste0(out_file_dir,"StrelkaSomaticWorkflow/results/variants/somatic.indels.vcf.gz.tbi"),
        snvs=paste0(out_file_dir,"StrelkaSomaticWorkflow/results/variants/somatic.snvs.vcf.gz"),
        snvs_index=paste0(out_file_dir,"StrelkaSomaticWorkflow/results/variants/somatic.snvs.vcf.gz.tbi")
    )
  )
  )
  
    ### Call INDELs before calling SNVs

    if (tumor!=""){
        jobs_report[["steps"]][["callSVManta"]]=call_sv_manta(
            bin_manta=bin_manta,
            tumor_bam=tumor,normal=normal,
            ref_genome=ref_genome,output_dir=output_dir,
            verbose=verbose,targeted=targeted,
            threads=threads,
            patient_id=patient_id
        )
        
        indel_candidates=jobs_report[["steps"]][["callSVManta"]]$out_files$indel_candidates
        indel_candidates=paste0(" --indelCandidates ",indel_candidates)
    }

    exome=""
    if (targeted){
        exome=" --exome "
    }

    

    exec_code=paste0(bin_strelka,tumour_input,normal_input," --referenceFasta ", 
    ref_genome , out_file2," --runDir ", out_file_dir, exome, "; ",
    paste0(out_file_dir,"StrelkaSomaticWorkflow/runWorkflow.py -m local -j ",threads))

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




call_somatic_variants_strelka=function(
    bin_strelka=build_default_tool_binary_list()$bin_strelka,
    bin_manta=build_default_tool_binary_list()$bin_manta,
    bin_bcftools=build_default_tool_binary_list()$bin_bcftools,
    bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
    bin_tabix=build_default_tool_binary_list()$bin_tabix,
    tumor="",normal="",
    ref_genome="",
    output_dir=".",output_name="",
    verbose=FALSE,targeted=FALSE,
    threads=3,exec_options="local"
){

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir,name="strelka_reports")
  job=build_job(executor_id=executor_id,task_id=task_id)


  id=""
  if(output_name!=""){
    id=output_name
  }else if(tumour!=""){
    tumor_input=paste(" --tumorBam ",tumor)
    normal_input=paste(" --normalBam ",normal)
    id=get_file_name(tumour)
  }else if(normal!=""){
    id=get_file_name(normal)
    normal_input=paste0(" --bam ", paste0(normal,collapse=" --bam "))
  }

    jobs_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    task_id=task_id,
    input_args = argg,
    out_file_dir=out_file_dir,
    out_files=list(
      workflow=paste0(out_file_dir,"MantaWorkflow/runWorkflow.py"),
      stats=list(
        aligment=paste0(out_file_dir,"MantaWorkflow/results/stats/alignmentStatsSummary.txt"),
        tsv=paste0(out_file_dir,"MantaWorkflow/results/stats/svCandidateGenerationStats.tsv"),
        xml=paste0(out_file_dir,"MantaWorkflow/results/stats/svCandidateGenerationStats.xml"),
        graph=paste0(out_file_dir,"MantaWorkflow/results/stats/svLocusGraphStats.tsv")
      ),
    variants=list(
        small_indel_candidate=paste0(out_file_dir,"MantaWorkflow/results/variants/candidateSmallIndels.vcf.gz"),
        sv_candidate=paste0(out_file_dir,"MantaWorkflow/results/variants/candidateSV.vcf.gz"),
        diploid_sv=paste0(out_file_dir,"MantaWorkflow/results/variants/diploidSV.vcf.gz"),
        somatic_sv=paste0(out_file_dir,"MantaWorkflow/results/variants/somaticSV.vcf.gz"),
        snvs=paste0(out_file_dir,"MantaWorkflow/results/variants/somatic.snvs.vcf.gz"),
        snvs_index=paste0(out_file_dir,"MantaWorkflow/results/variants/somatic.snvs.vcf.gz.tbi")
        )
    )
  )
  
    ### Call INDELs before calling SNVs

    if (tumor!=""){
        jobs_report[["steps"]][["callSVManta"]]=call_sv_manta(
            bin_manta=bin_manta,
            tumor_bam=tumor,normal=normal,
            ref_genome=ref_genome,output_dir=output_dir,
            verbose=verbose,targeted=targeted,
            threads=threads,
            patient_id=patient_id
        )
        
        indel_candidates=jobs_report[["steps"]][["callSVManta"]]$out_files$indel_candidates
        indel_candidates=paste0(" --indelCandidates ",indel_candidates)
    }

    exome=""
    if (targeted){
        exome=" --exome "
    }

    

    exec_code=paste0(bin_strelka,tumour_input,normal_input," --referenceFasta ", 
    ref_genome , out_file2," --runDir ", out_file_dir, exome, "; ",
    paste0(out_file_dir,"StrelkaSomaticWorkflow/runWorkflow.py -m local -j ",threads))

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