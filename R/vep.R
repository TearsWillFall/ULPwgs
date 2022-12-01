#' Create Read Orientation Model for Mutect2 Filter
#'
#' This function creates read orientation model for mutect2 filtering
#'
#' @param bin_vep [REQUIRED] Path to VEP binary.
#' @param bin_bgzip [REQUIRED] Path to bgzip binary.
#' @param bin_tabix [REQUIRED] Path to tabix binary.
#' @param vep_cache [REQUIRED] Path to vep cache location.
#' @param vcf [REQUIRED] Path to VCF file.
#' @param compress [OPTIONAL] Generate a compressed VCF
#' @param output_name [OPTIONAL] Name of output file
#' @param clean [OPTIONAL]Remove extra files.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param threads [OPTIONAL] Number of threads to split the work. Default 4
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


annotate_vep=function(
    rdata=NULL,
    selected=NULL,
    bin_vep=build_default_tool_binary_list()$bin_vep,
    bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
    bin_tabix=build_default_tool_binary_list()$bin_tabix,
    cache_vep=build_default_cache_list()$cache_vep,
    vcf="",
    output_name="",
    compress=TRUE,
    clean=FALSE,
    output_dir=".",
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=1,mode="local",
    executor_id=make_unique_id("callVEP"),
    task_name="callVEP",time="48:0:0",
    update_time=60,
    wait=FALSE,hold=NULL
){

    if(!is.null(rdata)){
      load(rdata)
      if(!is.null(selected)){
        cnr=sample_list[selected]
      }
    }

    id=""
    if(output_name!=""){
      id=output_name
    }else{
      id=get_file_name(vcf)
    }


    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir)
    job=build_job(executor_id=executor_id,task_id=task_id)


  out_file=paste0(out_file_dir,"/",id,".annotated.vcf")

  exec_code=paste(bin_vep,"-i",vcf,"-o",out_file,
  "--cache --port 3337 --everything --force_overwrite --vcf --fork ",threads," --dir ",cache_vep)

  if(mode=="batch"){
       out_file_dir2=set_dir(dir=out_file_dir,name="batch")
       batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
       threads=threads,output_dir=out_file_dir2)
       exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
  }


  if(verbose){
       print_verbose(job=job,arg=argg,exec_code=exec_code)
  }

    
    jobs_report=build_job_report(
            job_id=job,
            executor_id=executor_id,
            exec_code="", 
            task_id=task_id,
            input_args = argg,
            out_file_dir=out_file_dir,
            out_files=list(
            annotated_vcf=out_file
        )
    )

    
    error=execute_job(exec_code=exec_code)
    
    if(error!=0){
        stop("vep failed to run due to unknown error.
        Check std error for more information.")
    }

    if(compress){
        jobs_report[["steps"]][["compressAndIndex"]]<-compress_and_index_vcf_htslib(
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            vcf=out_file,
            output_dir=out_file_dir,output_name=id,
            clean=clean,verbose=verbose,
            batch_config=batch_config,
            threads=1,ram=ram,mode=mode,
            executor_id=task_id,
            time=time,
            hold=job
    )
    }

return(jobs_report)



 
}



#' Annotate Strelka VCF files with VEP annotation tool
#'
#' This function annotates Strelka VCF file using VEP
#'
#' @param bin_vep [REQUIRED] Path to VEP binary.
#' @param bin_bgzip [REQUIRED] Path to bgzip binary.
#' @param bin_tabix [REQUIRED] Path to tabix binary.
#' @param vep_cache [REQUIRED] Path to vep cache location.
#' @param vcf_snv [REQUIRED] Path to VCF file.
#' @param vcf_in [REQUIRED] Path to VCF file.
#' @param compress [OPTIONAL] Generate a compressed VCF
#' @param output_name [OPTIONAL] Name of output file
#' @param clean [OPTIONAL]Remove extra files.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param threads [OPTIONAL] Number of threads to split the work. Default 4
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


anotate_strelka_vep=function(
    bin_vep=build_default_tool_binary_list()$bin_vep,
    bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
    bin_tabix=build_default_tool_binary_list()$bin_tabix,
    cache_vep=build_default_cache_list()$cache_vep,
    vcf_snv="",
    vcf_indel="",
    extract_pass=TRUE,
    output_dir=".",
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=1,mode="local",
    executor_id=make_unique_id("callVEP"),
    task_name="callVEP",time="48:0:0",
    update_time=60,
    wait=FALSE,hold=NULL
){


    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir)
    job=build_job(executor_id=executor_id,task_id=task_id)


    if(is.null(vcf_snv)|is.null(vcf_indel)){
        stop("vcf_snv/vcf_indel arguments are required")
    }

     jobs_report=build_job_report(
      job_id=job,
      executor_id=executor_id,
      exec_code=list(),
      task_id=task_id,
      input_args = argg,
      out_file_dir=out_file_dir,
      out_files=list()
    )


    jobs_report[["steps"]][["annotateSnvStrelka"]]<- annotate_vep(
          bin_vep=bin_vep,
          bin_bgzip=bin_bgzip,
          bin_tabix=bin_tabix,
          cache_vep=cache_vep,
          vcf=vcf_snv,
          output_name=paste0("somatic.snv",ifelse(extract_pass,".PASS","")),
          output_dir=out_file_dir,
          verbose=verbose,
          batch_config=batch_config,
          threads=threads,
          ram=ram,mode=mode,
          executor_id=task_id,
          time=time,
          hold=hold
      )

    jobs_report[["steps"]][["annotateIndelStrelka"]]<-annotate_vep(
          bin_vep=bin_vep,
          bin_bgzip=bin_bgzip,
          bin_tabix=bin_tabix,
          cache_vep=cache_vep,
          vcf=vcf_indel,
          output_name=paste0("somatic.indel",ifelse(extract_pass,".PASS","")),
          output_dir=out_file_dir,
          verbose=verbose,
          batch_config=batch_config,
          threads=threads,
          ram=ram,mode=mode,
          executor_id=task_id,
          time=time,
          hold=hold
      )

    jobs_report$out_files=list(
        vcf_snv=unlist_lvl(jobs_report[["steps"]][["annotateeSnvStrelka"]],var="compressed_vcf"),
        vcf_indel=unlist_lvl(jobs_report[["steps"]][["annotateIndelStrelka"]],var="compressed_vcf")
    )


    return(jobs_report)

}