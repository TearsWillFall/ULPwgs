
#' Wrapper for MarkDuplicatesSpark from gatk
#'
#' This function removes duplicated reads (artifacts) found in aligned sequences and sorts the output bam.
#'
#' @param bam Path to the input file with the aligned sequence.
#' @param sif_gatk Path to gatk executable. Default path tools/picard/build/libs/picard.jar.
#' @param output_dir Path to the output directory.
#' @param tmp_dir Path to tmp directory.
#' @param remove_duplicates Remove all sequencing duplicates from BAM file. Default TRUE.
#' @param threads Number of threads . Default 4
#' @param ram RAM memory. Default 4
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Task EXECUTOR ID. Default "mardupsGATK"
#' @param task_name Task name. Default "mardupsGATK"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param verbose [OPTIONAL] Enables progress messages. Default False.#
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] Hold job until job is finished. Job ID. 
#' @import tidyverse
#' @export


markdups_gatk=function(
  sif_gatk=build_default_sif_list()$sif_gatk,bam="",output_dir=".",
  verbose=FALSE,batch_config=build_default_preprocess_config(),
  tmp_dir=".",threads=3,ram=4,remove_duplicates=TRUE,
  executor_id=make_unique_id("markdupsGATK"),task_name="markdupsGATK",
  mode="local",time="48:0:0",update_time=60,wait=FALSE,hold=NULL){

    argg <- as.list(environment())
    task_id=make_unique_id(task_name)

    out_file_dir=set_dir(dir=output_dir,name="markdups_reports")

    tmp=""
    if (!tmp_dir==""){
    tmp=paste0(" --tmp-dir ",tmp_dir)
    }

    dups=""
    if(remove_duplicates){
        dups=" --remove-all-duplicates"
    }

    out_file=paste0(out_file_dir,"/",get_file_name(bam),".sorted.rmdup.",
    get_file_ext(bam))
    
    out_file_md=paste0(out_file_dir,"/",get_file_name(bam),".gatk_rmdup.txt")
    exec_code=paste0(" singularity run ",sif_gatk," /gatk/gatk MarkDuplicatesSpark -I ",bam, " -O ",  out_file,
    " -M ",out_file_md," ",tmp," --conf \'spark.executor.cores=",threads,"\'", dups)


    job=build_job(executor_id=executor_id,task_id=task_id)

    if(mode=="batch"){
        out_file_dir2=set_dir(dir=out_file_dir,name="batch")
        batch_code=build_job_exec(job=job,time=time,ram=ram,threads=threads,
        output_dir=out_file_dir2,hold=hold)
        exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
    }
  
    if(verbose){
         print_verbose(job=job,arg=argg,exec_code=exec_code)
    }

    error=execute_job(exec_code=exec_code)
    if(error!=0){
        stop("markdups failed to run due to unknown error.
        Check std error for more information.")
    }

  job_report=build_job_report(
    job_id=job, 
    executor_id=executor_id, 
    exec_code=exec_code, 
    task_id=task_id,
    input_args=argg,
    out_file_dir=out_file_dir,
    out_files=list(
      bam=out_file,
      log=out_file_md
    )
  )

  if(wait&&mode=="batch"){
        job_validator(job=job_report$job_id,
        time=update_time,verbose=verbose,threads=threads)
  } 


  return(job_report)
}



#' Function for base quality recalibration
#'
#' This function recalibrates the base quality of the reads in two steps process based on GATK best practices guides.
#'
#' @param bin_samtools [REQUIRED] Path to picard executable. Default path tools/samtools/samtools.
#' @param sif_gatk [REQUIRED] Path to picard executable. Default path tools/gatk/gatk.
#' @param bin_picard [REQUIRED] Path to picard executable. Default path tools/picard/build/libs/picard.jar.
#' @param bam [REQUIRED]  Path to BAM file.
#' @param ref_genome [REQUIRED]  Path to reference genome.
#' @param dbsnp [REQUIRED] Known variant database.Requires atleast 1.
#' @param threads [OPTIONAL] Number of threads to split the work.
#' @param ram [OPTIONAL] RAM memory per thread.
#' @param tmp_dir [OPTIONAL] Path to the temporary directory.
#' @param clean Clean input files. Default TRUE.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Task EXECUTOR ID. Default "recalGATK"
#' @param task_id Task nam. Default "recalGATK"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param verbose [OPTIONAL] Enables progress messages. Default False.#
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export


recal_gatk=function(
  bin_samtools=build_default_tool_binary_list()$bin_samtools,
  sif_gatk=build_default_sif_list()$sif_gatk,
  bin_picard=build_default_tool_binary_list()$bin_picard,
  bam="",ref_genome="",dbsnp="",ram=4,threads=4,output_dir=".",
  tmp_dir=".",verbose=FALSE,batch_config=build_default_preprocess_config(),
  executor_id=make_unique_id("recalGATK"),
  task_name="recalGATK",clean=TRUE,mode="local",time="48:0:0",
  update_time=60,wait=FALSE,hold=NULL
  ){



  argg <- as.list(environment())
  task_id=make_unique_id(task_name)

  out_file_dir_before=set_dir(dir=output_dir,name="recal_reports/recal_before")
  out_file_dir_after=set_dir(dir=output_dir,name="recal_reports/recal_after")
  out_file_dir_main=set_dir(dir=output_dir,name="recal_reports")



  job=build_job(executor_id=executor_id,task_id=task_id)
  
  job_report=build_job_report(
    job_id=job, 
    executor_id=executor_id, 
    exec_code=list(), 
    task_id=task_id,
    input_args=argg,
    out_file_dir=list(
      main=out_file_dir_main,
      after=out_file_dir_after,
      before=out_file_dir_before
      ),
    out_files=list(
    )
  )
  
  job_report[["steps"]][["getChr"]] <- get_fai_reference_chr(
    fasta=ref_genome,verbose=verbose,output_dir=tmp_dir,
    batch_config=batch_config,executor_id=task_id,mode="local",
    threads=threads,ram=ram,
    time=time,update_time=update_time,wait=FALSE,hold=hold)


  regions=read.table(job_report[["steps"]][["getChr"]]$out_files$ref,
  sep="\t",header=TRUE)

  
  job_report[["steps"]][["par_bqsr_before"]] <- parallel_generate_BQSR_gatk(
    bin_samtools=bin_samtools,sif_gatk=sif_gatk,bam=bam,batch_config=batch_config,
    ref_genome=ref_genome,dbsnp=dbsnp,regions=regions,
    output_dir=out_file_dir_before,clean=clean,tmp_dir=tmp_dir,
    verbose=verbose,executor_id=task_id,mode=mode,threads=threads,ram=ram,
    time=time,update_time=update_time,wait=FALSE,
    hold=hold)

  job_report[["steps"]][["par_apply_bqsr"]] <- parallel_apply_BQSR_gatk(
    bin_samtools=bin_samtools,sif_gatk=sif_gatk,bin_picard=bin_picard,
    bam=bam,ref_genome=ref_genome,batch_config=batch_config,
    rec_table=job_report[["steps"]][["par_bqsr_before"]][["steps"]][["gather_bqsr_report"]]$out_files$table,
    output_dir=tmp_dir,regions=regions,tmp_dir=tmp_dir,
    clean=clean,verbose=verbose,executor_id=task_id,
    threads=threads,mode=mode,ram=ram,time=time,
    update_time=update_time,wait=FALSE,
    hold=unlist_lvl(job_report[["steps"]][["par_bqsr_before"]],var="job_id"))

  
  job_report[["steps"]][["sort_and_index"]] <- sort_and_index_bam_samtools(
    bin_samtools=bin_samtools,
    bam=job_report[["steps"]][["par_apply_bqsr"]][["steps"]][["gather_bam"]]$out_files$bam,
    output_dir=out_file_dir_main,batch_config=batch_config,
    ram=ram,verbose=verbose,threads=threads,sort=TRUE,
    stats="",index=TRUE,clean=clean,
    mode=mode,executor_id=task_id,time=time,
    update_time=update_time,wait=FALSE,
    hold=unlist_lvl(job_report[["steps"]][["par_apply_bqsr"]],var="job_id"))
  
 

  job_report[["steps"]][["par_bqsr_after"]] <- parallel_generate_BQSR_gatk(
    bin_samtools=bin_samtools,sif_gatk=sif_gatk,tmp_dir=tmp_dir,batch_config=batch_config,
    bam=job_report[["steps"]][["sort_and_index"]][["steps"]][["sort"]]$out_files$bam,
    ref_genome=ref_genome,dbsnp=dbsnp,threads=threads,regions=regions,
    output_dir=out_file_dir_after,verbose=verbose,executor_id=task_id,clean=clean,
    mode=mode,ram=ram,time=time,update_time=update_time,
    wait=FALSE,hold=unlist_lvl(job_report[["steps"]][["sort_and_index"]],var="job_id"))
  

 job_report[["steps"]][["analyse_covariates"]] <- analyze_covariates_gatk(
  sif_gatk=sif_gatk,tmp_dir=tmp_dir,batch_config=batch_config,
  before=job_report[["steps"]][["par_bqsr_before"]][["steps"]][["gather_bqsr_report"]]$out_files$table,
  after=job_report[["steps"]][["par_bqsr_after"]][["steps"]][["gather_bqsr_report"]]$out_files$table,
  output_dir=out_file_dir_main,executor_id=task_id,mode=mode,threads=threads,
  ram=ram,time=time,update_time=update_time,wait=FALSE,
  hold=unlist_lvl(job_report[["steps"]][["par_bqsr_after"]],var="job_id"))
  

  if(wait&&mode=="batch"){
    job_validator(job=job_report$job_id,
    time=update_time,verbose=verbose,threads=threads)
  }
  return(job_report)
}





#' Wrapper of BaseRecalibrator function from gatk
#'
#' Generates a recalibration table based on various covariates.
#' This function wraps around gatk BaseRecalibrator function.
#' For more information about this function: https://gatk.broadinstitute.org/hc/en-us/articles/360036898312-BaseRecalibrator
#'
#' @param bam [REQUIRED] Path to the BAM file.
#' @param sif_gatk [REQUIRED] Path to gatk executable. Default tools/gatk/gatk.
#' @param ref_genome [REQUIRED] Path to reference genome
#' @param dbsnp [REQUIRED] Path to known snp positions in VCF format. Multiple vcf can be supplied as a vector.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param tmp_dir [OPTIONAL] Path to the temporary directory.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Task EXECUTOR ID. Default "generateBQSR"
#' @param task_id Task NAME ID. Default "generateBQSR"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param threads [OPTIONAL] Number of threads to split the work.
#' @param ram [OPTIONAL] If batch mode. RAM memory in GB per job. Default 1
#' @param update_time [OPTIONAL] If batch mode. Show job updates every update time. Default 60
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID.
#' @export

generate_BQSR_gatk=function(
  region="",rdata=NULL,selected=NULL,
  sif_gatk=build_default_sif_list()$sif_gatk,bam="",
  ref_genome="",dbsnp="",output_dir=".",tmp_dir=".",verbose=FALSE,
  batch_config=build_default_preprocess_config(),
  threads=4,ram=4,mode="local",
  executor_id=make_unique_id("generateBQSR"),task_name="generateBQSR",
  time="48:0:0",update_time=60,wait=FALSE,hold=NULL
){  

  if(!is.null(rdata)){
    load(rdata)
    if(!is.null(selected)){
      region=region_list[selected]
    }
  }

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir)

  if(tmp_dir!=""){
    tmp_dir=paste0(" --tmp-dir ",tmp_dir)
  }

  reg=""
  if (region==""){
      out_file=paste0(out_file_dir,get_file_name(bam),".recal.table")
  }else{
      reg=paste0(" -L ",region)
      out_file=paste0(out_file_dir,get_file_name(bam),".",region,".recal.table")
  }

  ## Multiple vcf with snps can be given

  if (dbsnp!=""){
    dbsnp=paste(" --known-sites ",dbsnp,collapse=" ")
  }

  exec_code=paste0("singularity run ", sif_gatk," /gatk/gatk BaseRecalibrator -I ",bam, " -R ", ref_genome, dbsnp,
  reg," -O ",out_file,tmp_dir)


  job=build_job(executor_id=executor_id,task_id=task_id)
  if(mode=="batch"){
       out_file_dir2=set_dir(dir=out_file_dir,name="batch")
       batch_code=build_job_exec(job=job,time=time,ram=ram,threads=threads,
       output_dir=out_file_dir2,hold=hold)
       exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
  }

  if(verbose){
     print_verbose(job=job,arg=argg,exec_code=exec_code)
  }
  error=execute_job(exec_code=exec_code)
  if(error!=0){
    stop("gatk failed to run due to unknown error.
    Check std error for more information.")
  }

   job_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    exec_code=exec_code, 
    task_id=task_id,
    input_args = argg,
    out_file_dir=out_file_dir,
    out_files=list(
      recal_table=out_file
    )
  )


  if(wait&&mode=="batch"){
      job_validator(job=job_report$job_id,
      time=update_time,verbose=verbose,threads=threads)
  }

  return(job_report)
}



#' Multiregion parallelization of generate_BQSR function
#'
#' Generates a recalibration table based on various covariates.
#' This function wraps around gatk BaseRecalibrator function.
#' For more information about this function: https://gatk.broadinstitute.org/hc/en-us/articles/360036898312-BaseRecalibrator
#'
#' @param bam [REQUIRED] Path to the BAM file.
#' @param bin_samtools [REQUIRED] Path to gatk executable. Default tools/samtools/samtools.
#' @param sif_gatk [REQUIRED] Path to gatk executable. Default tools/gatk/gatk.
#' @param ref_genome [REQUIRED] Path to reference genome
#' @param dbsnp [REQUIRED] Path to known snp positions in VCF format. Multiple vcf can be supplied as a vector.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param tmp_dir [OPTIONAL] Path to the temporary directory.
#' @param regions [OPTIONAL] Regions to parallelize through.
#' @param clean Remove intermediary files. Default TRUE
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Task EXECUTOR ID. Default "par_generateBQSR"
#' @param task_name Task name. Default "par_generateBQSR"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param threads Number of threads to split the work. Default 3
#' @param ram [OPTIONAL] If batch mode. RAM memory in GB per job. Default 1
#' @param update_time [OPTIONAL] If batch mode. Show job updates every update time. Default 60
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export
#' @import pbapply



parallel_generate_BQSR_gatk=function(
  bin_samtools=build_default_tool_binary_list()$bin_samtools,
  sif_gatk=build_default_sif_list()$sif_gatk,bam="",
  regions=NULL,ref_genome="", clean=TRUE, dbsnp="",
  tmp_dir=".",threads=3,ram=4,
  executor_id=make_unique_id("par_generateBQSR"),
  task_name="par_generateBQSR",output_dir=".",
  verbose=FALSE,batch_config=build_default_preprocess_config(),
  mode="local",time="48:0:0",
  update_time=60,wait=FALSE,hold=NULL){

  options(scipen = 999)
 
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

  if(is.null(regions)){

    job_report[["steps"]][["getChr"]] <- get_bam_reference_chr(
      bin_samtools=bin_samtools,
      bam=bam,verbose=verbose,output_dir=tmp_dir,
      executor_id=task_id,mode="local",threads=threads,ram=ram,
      time=time,update_time=update_time,wait=FALSE,hold=hold)


    regions=read.table(job_report[["steps"]][["getChr"]]$out_files$ref,
    sep="\t",header=TRUE)
  }

  regions$start=regions$start+1
  regions=regions %>% dplyr::mutate(region=paste0(chr,":",start,"-",end))


    region_list=regions$region
    names(region_list)=regions$region

  if(mode=="local"){
      job_report[["steps"]][["generate_bqsr_report"]]<-parallel::mclapply(
        region_list,FUN=function(region){
        job_report <- generate_BQSR_gatk(
        region=region,tmp_dir=tmp_dir,batch_config=batch_config,
        sif_gatk=sif_gatk,bam=bam,ref_genome=ref_genome,
        dbsnp=dbsnp,output_dir=out_file_dir,verbose=verbose,
        executor_id=task_id
    )
  },mc.cores=threads)

  }else if(mode=="batch"){
        rdata_file=paste0(tmp_dir,"/",job,".regions.RData")
        executor_id=task_id
        save(
          region_list,
          sif_gatk,bam,
          ref_genome,dbsnp,
          output_dir,
          executor_id,
          verbose,
          tmp_dir,
          file = rdata_file
        )
        exec_code=paste0("Rscript -e \"ULPwgs::generate_BQSR_gatk(rdata=\\\"",
        rdata_file,"\\\",selected=$SGE_TASK_ID)\"")
        out_file_dir2=set_dir(dir=out_file_dir,name="batch")
        batch_code=build_job_exec(job=job,time=time,ram=ram,
        threads=2,output_dir=out_file_dir2,
        hold=hold,array=length(region_list))
        exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)


        if(verbose){
            print_verbose(job=job,arg=argg,exec_code=exec_code)
        }
        error=execute_job(exec_code=exec_code)
        if(error!=0){
            stop("gatk failed to run due to unknown error.
            Check std error for more information.")
        }

        job_report[["steps"]][["generate_bqsr_report"]]<- build_job_report(
              job_id=job,
              executor_id=executor_id,
              exec_code=exec_code, 
              task_id=task_id,
              input_args=argg,
              out_file_dir=out_file_dir,
              out_files=list(
                  recal_table=paste0(out_file_dir,get_file_name(bam),".",region_list,".recal.table")
                )
        )

  }
  

  job_report[["steps"]][["gather_bqsr_report"]]<- gather_BQSR_reports_gatk(
    sif_gatk=sif_gatk,
    report=unlist_lvl(named_list=job_report[["steps"]][["generate_bqsr_report"]],var="recal_table"),
    executor_id=task_id,tmp_dir=tmp_dir,output_dir=out_file_dir,clean=clean,
    output_name=get_file_name(bam),verbose=verbose,mode=mode,time=time,
    threads=threads,ram=ram,update_time=update_time,wait=FALSE,
    hold=unlist_lvl(named_list=job_report[["steps"]][["generate_bqsr_report"]],var="job_id",recursive=TRUE))
  
  if(wait&&mode=="batch"){
    job_validator(job=unlist_lvl(named_list=job_report[["steps"]],var="job_id"),
    time=update_time,verbose=verbose,threads=threads)
  }


  return(job_report)
}



#' Wrapper around gatk GatherBQSRReports function
#'
#' This functions collects the Recalibration reports generated from scattered parallel_generate_BQSR output
#' This function wraps around gatk GatherBQSRReports function.
#' For more information about this function: https://gatk.broadinstitute.org/hc/en-us/articles/360036829851-GatherBQSRReports
#'
#' @param sif_gatk [REQUIRED] Path to gatk executable. Default tools/gatk/gatk.
#' @param report [REQUIRED] Path to indivual reports.
#' @param output_name [OPTIONAL] Name for the output report file.
#' @param clean Clean input files. Default TRUE.
#' @param executor_id Task EXECUTOR ID. Default "gatherBQSR"
#' @param tmp_dir [OPTIONAL] Path to the temporary directory.
#' @param task_name Task name. Default "gatherBQSR"
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param tmp_dir [OPTIONAL] Path to temporary directory.
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param threads Number of threads to split the work. Default 3
#' @param ram [OPTIONAL] If batch mode. RAM memory in GB per job. Default 1
#' @param update_time [OPTIONAL] If batch mode. Show job updates every update time. Default 60
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export

gather_BQSR_reports_gatk=function(
  sif_gatk=build_default_sif_list()$sif_gatk,
  report="",output_name="Report",clean=FALSE,
  output_dir=".",tmp_dir=".",verbose=FALSE,
  executor_id=make_unique_id("gatherBQSR"),
  task_name="gatherBQSR",
  batch_config=build_default_preprocess_config(),
  mode="local",time="48:0:0",
  threads=4,ram=4,update_time=60,wait=FALSE,hold=NULL){

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir)

  if(tmp_dir!=""){
    tmp_dir=paste0(" --tmp-dir ",tmp_dir)
  }
  out_file=paste0(out_file_dir,output_name,".recal.table")
  exec_code=paste0("singularity run ",sif_gatk,
  " /gatk/gatk GatherBQSRReports ",paste(" -I ",report,collapse=" "),
    " -O ",out_file,tmp_dir)


  if(clean){
    exec_code=paste(exec_code," && rm",paste(report,collapse=" "))
  }

  job=build_job(executor_id = executor_id,task_id=task_id)
  

  if(mode=="batch"){
  
      out_file_dir2=set_dir(dir=out_file_dir,name="batch")
       batch_code=build_job_exec(job=job,
       time=time,ram=ram,threads=threads,
       output_dir=out_file_dir2,hold=hold)
       exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
  }

  if(verbose){
     print_verbose(job=job,arg=argg,exec_code=exec_code)
  }

  error=execute_job(exec_code=exec_code)
  if(error!=0){
    stop("gatk failed to run due to unknown error.
    Check std error for more information.")
  }

  job_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    exec_code=exec_code, 
    task_id=task_id, 
    input_args = argg,
    out_file_dir=out_file_dir,
      out_files=list(
        table=out_file)
  )


  if(wait&&mode=="batch"){
    job_validator(job=job_report$job_id,time=update_time,verbose=verbose,
    threads=threads)
  }
  return(job_report)
}

#' Gather multiple BQSR reports with GATK
#'
#' Gathers multiple BQSR recalibration reports generated from different genomic regions into a single report.
#' This function wraps around the GATK GatherBQSRReports function.
#' For more information about this function: https://gatk.broadinstitute.org/hc/en-us/articles/360036901871-GatherBQSRReports
#'
#' @param sif_gatk \link{REQUIRED} Path to GATK singularity image. Default tools/gatk/gatk.
#' @param report \link{REQUIRED} Path(s) to BQSR recalibration report files to gather.
#' @param output_name \link{OPTIONAL} Name for the output recalibration report file. Default "Report"
#' @param clean_reports \link{OPTIONAL} If TRUE, remove the input report files after gathering. Default FALSE
#' @param verbose \link{OPTIONAL} Enables progress messages. Default False.
#' @param threads \link{OPTIONAL} Number of threads. Default 3
#' @param ram \link{OPTIONAL} RAM memory to assign to each thread in GB. Default 1GB
#' @param output_dir \link{OPTIONAL} Path to the output directory.
#' @param tmp_dir \link{OPTIONAL} Path to the temporary directory.
#' @param mode \link{OPTIONAL} Where to parallelize. Default local. Options \link{"local","batch"}
#' @param executor_id \link{OPTIONAL} Task EXECUTOR ID. Default "gatherBQSR"
#' @param task_name \link{OPTIONAL} Task name. Default "gatherBQSR"
#' @param time \link{OPTIONAL} If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time \link{OPTIONAL} If batch mode. Job update time in seconds. Default 60.
#' @param wait \link{OPTIONAL} If batch mode wait for batch to finish. Default FALSE
#' @param hold \link{OPTIONAL} Hold job until job is finished. Job ID.
#'
#' @export


new_gather_BQSR_reports_gatk=function(
  sif_gatk=build_default_sif_list()$sif_gatk,
  report=NULL,
  output_name="Report",
  clean_reports=FALSE,
  ...
  ){

      run_main=function(
        .env
      ){
      
        .this.env=environment()
        append_env(to=.this.env,from=.env)
        set_main(.env=.this.env)

        if(!is.null(tmp_dir)){
          tmp_dir=paste0(" --tmp-dir ",tmp_dir)
        }

        .main$out_files$rec_table=paste0(out_file_dir,"/",input_id,".recal.table")

        .main$exec_code=paste0(
          "singularity run ",sif_gatk, " /gatk/gatk GatherBQSRReports ",
          paste(" -I ",report,collapse=" "),
          " -O ", .main$out_files$rec_table,tmp_dir
        )

        if(clean_reports){
          .main$exec_code=paste(.main$exec_code, " && ",paste(paste0(report,"*"),collapse=" "))
        }
          
        run_job(.env=.this.env)
        .env$.main <- .main

      }

      .base.env=environment()
      list2env(list(...),envir=.base.env)

      set_env_vars(
        .env= .base.env,
        vars="output_name"
      )
      
      launch(.env=.base.env)
    
}





#' Wrapper of applyBQSR function gatk
#'
#' Applies numerical corrections to each individual basecall based on the covariates analyzed before.
#' This function wraps around gatk applyBQSR  function.
#' For more information about this function: https://gatk.broadinstitute.org/hc/en-us/articles/360050814312-ApplyBQSR
#'
#' @param bam [REQUIRED] Path to the BAM file.
#' @param sif_gatk [REQUIRED] Path to gatk executable. Default tools/gatk/gatk.
#' @param ref_genome [REQUIRED] Path to reference genome
#' @param rec_table [REQUIRED] Path to covariates table generated by generate_BSQR.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param tmp_dir [OPTIONAL] Path to the temporary directory.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param threads Number of threads to split the work. Default 3
#' @param ram [OPTIONAL] If batch mode. RAM memory in GB per job. Default 1
#' @param executor_id [OPTIONAL] Task executor name. Default "applyBQSR"
#' @param task_name [OPTIONAL] Task name. Default "applyBQSR"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param ram [OPTIONAL] If batch mode. RAM memory in GB per job. Default 1
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export


apply_BQSR_gatk=function(
  region="",rdata=NULL,selected=NULL,
  sif_gatk=build_default_sif_list()$sif_gatk,bam="",ref_genome="",
  rec_table="",output_dir=".",verbose=FALSE,tmp_dir=".",
  batch_config=build_default_preprocess_config(),mode="local", threads=4,ram=4,
  executor_id=make_unique_id("applyBQSR"),task_name="applyBQSR",time="48:0:0",
  update_time=60,wait=TRUE,hold=NULL){


  if(!is.null(rdata)){
    load(rdata)
    if(!is.null(selected)){
      region=region_list[selected]
    }
  }
  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir)
 
  if(tmp_dir!=""){
    tmp_dir=paste0(" --tmp-dir ",tmp_dir)
  }

  reg=""
  if (region==""){
      out_file=paste0(" ", out_file_dir,"/", get_file_name(bam),".recal.",get_file_ext(bam))
  }else{
      reg=paste0(" -L ",strsplit(region,"__")[[1]][2], " ")
      out_file=paste0(out_file_dir,"/", get_file_name(bam),".",region,".recal.",get_file_ext(bam))
  }
  exec_code=paste("singularity run ", sif_gatk,
  " /gatk/gatk ApplyBQSR -I ",bam, " -R ", ref_genome,
   " --bqsr-recal-file ",rec_table, " -O ",out_file,reg)
   
  job=build_job(executor_id=executor_id,task_id=task_id)

  if(mode=="batch"){
       out_file_dir2=set_dir(dir=out_file_dir,name="batch")
       batch_code=build_job_exec(job=job,time=time,ram=ram,threads=threads,
       output_dir=out_file_dir2,hold=hold)
       exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
  }

  if(verbose){
     print_verbose(job=job,arg=argg,exec_code=exec_code)
  }


  error=execute_job(exec_code=exec_code)
  if(error!=0){
    stop("gatk failed to run due to unknown error.
    Check std error for more information.")
  }

  job_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    exec_code=exec_code,
    task_id=task_id,
    input_args = argg,
    out_file_dir=out_file_dir,
    out_files=list(
      recal_bam=out_file,
      recal_bai=sub(".bam",".bai",out_file)
    )
  )

  if(wait&&mode=="batch"){
      job_validator(job=job,
      time=update_time,verbose=verbose,threads=threads)
  }

  return(job_report)
}

#' Multiregion parallelization of apply_BQSR function
#'
#' Recalibrates
#' Applies numerical corrections to each individual basecall based on the covariates analyzed before.
#' For more information about this function: https://gatk.broadinstitute.org/hc/en-us/articles/360050814312-ApplyBQSR
#'
#' @param bam [REQUIRED] Path to the BAM file.
#' @param bin_samtools [REQUIRED] Path to samtools executable. Default tools/samtools/samtools.
#' @param sif_gatk [REQUIRED] Path to gatk executable. Default tools/gatk/gatk.
#' @param bin_picard [REQUIRED] Path to picard executable. Default tools/picard/build/libs/picard.jar
#' @param ref_genome [REQUIRED] Path to reference genome
#' @param rec_table [REQUIRED] Path to the recalibratio table.
#' @param clean Clean intermediary files Default TRUE
#' @param tmp_dir [OPTIONAL] Path to the temporary directory.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param regions [OPTIONAL] Regions to parallelize through.
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id [OPTIONAL] Task executor name. Default "par_applyBQSR"
#' @param task_name [OPTIONAL] Task name. Default "par_applyBQSR"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param threads [OPTIONAL] Number of threads for the main job. Default 4
#' @param ram [OPTIONAL] If batch mode. RAM memory in GB per job. Default 1
#' @param update_time [OPTIONAL] If batch mode. Show job updates every update time. Default 60
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export
#' @import pbapply

parallel_apply_BQSR_gatk=function(
  bin_samtools=build_default_tool_binary_list()$bin_samtools,
  sif_gatk=build_default_sif_list()$sif_gatk,
  bin_picard=build_default_tool_binary_list()$bin_picard,
  bam="",regions=NULL,ref_genome="",rec_table="",clean=TRUE,
  output_dir=".",verbose=FALSE,tmp_dir=".",
  batch_config=build_default_preprocess_config(),mode="local",
  executor_id=make_unique("par_applyBQSR"),
  task_name="par_applyBQSR",
  time="48:0:0",threads=4,ram=4,
  update_time=60,wait=FALSE, hold=NULL){

  options(scipen = 999)
 

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir)

   if(is.null(regions)){
      job_report[["steps"]][["getChr"]] <- get_bam_reference_chr(
      bin_samtools=bin_samtools,
      bam=bam,verbose=verbose,output_dir=tmp_dir,
      executor_id=task_id,mode="local",threads=threads,ram=ram,
      time=time,update_time=update_time,wait=FALSE,hold=hold)


      regions=read.table(job_report[["steps"]][["getChr"]]$out_files$ref,
      sep="\t",header=TRUE)
  }
  regions$start=regions$start+1
  regions$pos=1:nrow(regions)
  regions=regions %>% dplyr::mutate(region=paste0(pos,"__",chr,":",start,"-",end))

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

    region_list=regions$region
    names(region_list)=regions$region


    if(mode=="local"){
      job_report[["steps"]][["apply_bqsr"]]=parallel::mclapply(
      region_list,FUN=function(region){

        job_report<-apply_BQSR_gatk(
        region=region,batch_config=batch_config,
        sif_gatk=sif_gatk,bam=bam,ref_genome=ref_genome,
        executor_id=executor_id,rec_table=rec_table,tmp_dir=tmp_dir,
        output_dir=out_file_dir,verbose=verbose)},
        mc.cores=threads)
  }else if(mode=="batch"){
        rdata_file=paste0(tmp_dir,"/",job,".regions.RData")
        executor_id=task_id
        save(
          region_list,
          sif_gatk,
          bam,
          ref_genome,
          rec_table,
          output_dir,
          executor_id,
          verbose,
          tmp_dir,
          file = rdata_file
        )
        exec_code=paste0("Rscript -e \"ULPwgs::apply_BQSR_gatk(rdata=\\\"",
        rdata_file,"\\\",selected=$SGE_TASK_ID)\"")
        out_file_dir2=set_dir(dir=out_file_dir,name="batch")
        batch_code=build_job_exec(job=job,time=time,ram=ram,
        threads=2,output_dir=out_file_dir2,
        hold=hold,array=length(region_list))
        exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)

        if(verbose){
            print_verbose(job=job,arg=argg,exec_code=exec_code)
        }

        error=execute_job(exec_code=exec_code)
        
        if(error!=0){
            stop("gatk failed to run due to unknown error.
            Check std error for more information.")
        }

        job_report[["steps"]][["apply_bqsr"]] <- build_job_report(
              job_id=job,
              executor_id=executor_id,
              exec_code=exec_code, 
              task_id=task_id,
              input_args=argg,
              out_file_dir=out_file_dir,
              out_files=list(
                  recal_bam=paste0(out_file_dir,"/", get_file_name(bam),".",region_list,
                  ".recal.",get_file_ext(bam))
                )
        )

  }
  
  job_report[["steps"]][["gather_bam"]]=gather_bam_files_picard(
    bin_picard=bin_picard,
    bam=unlist_lvl(job_report[["steps"]][["apply_bqsr"]],var="recal_bam"),
    output_dir=out_file_dir,
    output_name=paste0(get_file_name(bam),".recal.sorted.rmdup.sorted"),
    executor_id=task_id,mode=mode,time=time,threads=threads,ram=ram,
    update_time=update_time,wait=FALSE,
    clean=clean,
    hold=unlist_lvl(job_report[["steps"]][["apply_bqsr"]],var="job_id",recursive=TRUE))


    if(wait&&mode=="batch"){
      job_validator(job=unlist_lvl(named_list=job_report[["steps"]][["apply_bqsr"]],
      var="job_id"),time=update_time,verbose=verbose,threads=threads)
    }

  return(job_report)
}






#' Wrapper of AnalyzeCovariates function in gatk
#'
#' Generates a report of the recalibrated values.
#' This function wraps around gatk AnalyzeCovariates function.
#' For more information about this function: https://gatk.broadinstitute.org/hc/en-us/articles/360037066912-AnalyzeCovariates
#'
#' @param sif_gatk [REQUIRED] Path to gatk executable. Default tools/gatk/gatk.
#' @param before [REQUIRED] Recalibration table produced by generate_BQSR function.
#' @param after [OPTIONAL] Recalibration table produced by generate_BQSR function.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param tmp_dir [OPTIONAL] Path to the temporary directory.
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


analyze_covariates_gatk=function(
  sif_gatk=build_default_sif_list()$sif_gatk,before="",after="",
  output_dir=".",verbose=FALSE,tmp_dir=".",
  batch_config=build_default_preprocess_config(),
  threads=4,ram=4,mode="local",
  executor_id=make_unique_id("recalCovariates"),
  task_name="recalCovariates",time="48:0:0",
  update_time=60,wait=FALSE,hold=NULL
){

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir)



  if(tmp_dir!=""){
    tmp_dir=paste0(" --tmp-dir ",tmp_dir)
  }
 
  if (before!="" & after==""){
    out_file=paste0(out_file_dir,"/",get_file_name(before),
    "_covariates_analysis_before.pdf")

    out_file_csv=paste0(out_file_dir,"/",get_file_name(before),
    "_covariates_analysis_before.csv")

    exec_code=paste0("singularity run ",sif_gatk,
    " /gatk/gatk  AnalyzeCovariates -bqsr ",before, " -plots ",out_file,tmp_dir," -csv ",out_file_csv)
  }else if(before=="" & after!=""){

    out_file=paste0(out_file_dir,"/",get_file_name(after),
       "_covariates_analysis_after.pdf")

    out_file_csv=paste0(out_file_dir,"/",get_file_name(after),
    "_covariates_analysis_after.csv")

    exec_code=paste0("singularity run ",sif_gatk,
    " /gatk/gatk AnalyzeCovariates -bqsr ",after, " -plots ",out_file,tmp_dir," -csv ",out_file_csv)
  }else{
    out_file=paste0(out_file_dir,"/",get_file_name(before),"_covariates_analysis.pdf")
    out_file_csv=paste0(out_file_dir,"/",get_file_name(before),"_covariates_analysis.csv")
    exec_code=paste0("singularity run ",sif_gatk,
    " /gatk/gatk  AnalyzeCovariates -before ",before," -after ",after,
      " -plots ",out_file,tmp_dir," -csv ",out_file_csv)
  }

  job=build_job(executor_id=executor_id,task=task_id)

  if(mode=="batch"){
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
    stop("gatk failed to run due to unknown error.
    Check std error for more information.")
  }

  job_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    exec_code=exec_code, 
    task_id=task_id,
    input_args = argg,
    out_file_dir=out_file_dir,
    out_files=list(
      plot=out_file)
  )


  if(wait&&mode=="batch"){
    job_validator(job=job_report$job_id,time=update_time,
    verbose=verbose,threads=threads)
  }

  return(job_report)
}





#' Run GATK AnalyzeCovariates to compare BQSR covariates
#'
#' Wrapper for GATK's `AnalyzeCovariates` that runs inside the project's
#' Singularity image. Use `before` and/or `after` to provide the
#' pre- and after-recalibration tables (from `BaseRecalibrator`); the
#' function will create PDF plots and CSV summaries in the task output
#' directory. If only `before` or `after` is provided, the function will
#' produce the corresponding single plot; if both are provided, a
#' before/after comparison is produced.
#'
#' @param sif_gatk Path to the Singularity image for GATK. Defaults to
#'   `build_default_sif_list()$sif_gatk`.
#' @param before Path to the pre-recalibration table (output of
#'   `BaseRecalibrator`) for the "before" comparison. Optional.
#' @param after Path to the after-recalibration table for the "after"
#'   comparison. Optional.
#' @param output_name Base name used for output files. Defaults to
#'   "sample".
#' @param ... Additional arguments forwarded into the internal job
#'   environment (for example `output_dir`, `tmp_dir`, `mode`,
#'   `batch_config`, `threads`, `ram`, `verbose`).
#'
#' @return Invisibly returns the job report created by the internal job
#'   runner and writes PDF and CSV outputs to the task output directory.
#' @export


new_analyze_covariates_gatk=function(
  sif_gatk=build_default_sif_list()$sif_gatk,
  before=NULL,
  after=NULL,
  output_name="sample",
  ...
){

  run_main=function(
    .env
  ){
    .this.env=environment()
    append_env(to=.this.env,from=.env)
    set_main(.env=.this.env)

      if(!is.null(tmp_dir)){
          tmp_dir=paste0(" --tmp-dir ",tmp_dir)
        }

    if (!is.null(before) & is.null(after)){
      .main$out_files$before_pdf=paste0(out_file_dir,"/",input_id,
      "_covariates_analysis_before.pdf")

      .main$out_files$before_csv=paste0(out_file_dir,"/",input_id,
      "_covariates_analysis_before.csv")

      .main$exec_code=paste0("singularity run ",sif_gatk,
      " /gatk/gatk  AnalyzeCovariates -bqsr ",before, 
      " -plots ",.main$out_files$before_pdf,tmp_dir,
      " -csv ",.main$out_files$before_csv)
  }else if(is.null(before) & !is.null(after)){

    .main$out_files$after_pdf=paste0(out_file_dir,"/",input_id,
       "_covariates_analysis_after.pdf")

    .main$out_files$after_csv=paste0(out_file_dir,"/",input_id,
    "_covariates_analysis_after.csv")

    .main$exec_code=paste0("singularity run ",sif_gatk,
    " /gatk/gatk AnalyzeCovariates -bqsr ",after, 
    " -plots ",.main$out_files$after_pdf,tmp_dir,
    " -csv ",.main$out_files$after_csv)
  }else{
    .main$out_files$pdf=paste0(out_file_dir,"/",input_id,"_covariates_analysis.pdf")
    .main$out_files$csv=paste0(out_file_dir,"/",input_id,"_covariates_analysis.csv")
    .main$exec_code=paste0("singularity run ",sif_gatk,
    " /gatk/gatk  AnalyzeCovariates -before ",before," -after ",after,
      " -plots ",.main$out_files$pdf,tmp_dir," -csv ",.main$out_files$csv)
  }

    run_job(.env=.this.env)
    .env$.main <- .main
  }

   .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="output_name"
    )

    launch(.env=.base.env)

}





#' Variant Calling using Mutect2
#'
#' This function functions calls Mutect2 for variant calling.
#' If a vector of tumour samples are provided these will be processed in multi-sample mode.
#' To run in tumour-normal mode suppply a single tumour and normal sample.
#' If no normal is supplied this will run in tumour only.
#' TO DO// Implement mitochondrial mode feature
#' 
#' For more information read:
#' https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2
#'
#' @param sif_gatk [REQUIRED] Path to gatk sif file.
#' @param tumour [REQUIRED] Path to tumour BAM file.
#' @param normal [OPTIONAL] Path to normal BAM file.
#' @param ref_genome [REQUIRED] Path to reference genome fasta file.
#' @param rdata [OPTIONAL] Import R data information with list of BAM.
#' @param selected [OPTIONAL] Select BAM from list.
#' @param germ_resource [REQUIRED]Path to germline resources vcf file.
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
#' @param pon [OPTIONAL] Path to panel of normal.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param mnps [OPTIONAL] Report MNPs in vcf file.
#' @param contamination [OPTIONAL] Produce sample cross-contamination reports. Default TRUE.
#' @param orientation [OPTIONAL] Produce read orientation inforamtion. Default FALSE
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



mutect2_gatk=function(region="",
  rdata=NULL,selected=NULL,
  sif_gatk=build_default_sif_list()$sif_gatk,
  tumour="",normal=NA,output_name="",
  ref_genome=build_default_reference_list()$HG19$reference$genome,
  germ_resource=build_default_reference_list()$HG19$variant$germ_reference,
  pon="",output_dir=".",tmp_dir=".",
  verbose=FALSE,orientation=FALSE,mnps=FALSE,
  batch_config=build_default_preprocess_config(),
  threads=4,ram=4,mode="local",
  executor_id=make_unique_id("Mutect2"),
  task_name="Mutect2",time="48:0:0",
  update_time=60,wait=FALSE,hold=NULL
){

  if(!is.null(rdata)){
    load(rdata)
    if(!is.null(selected)){
      region=region_list[selected]
    }
  }

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir,name="vcf")
  job=build_job(executor_id=executor_id,task=task_id)
  
  
  if(tmp_dir!=""){
    tmp_dir=paste0(" --tmp-dir ",tmp_dir)
  }
  
  id=""
  if(output_name!=""){
    id=output_name
  }else{  
    id=get_file_name(tumour[1])
  }
  
  reg=""
  if (region==""){
      out_file=paste0(out_file_dir,"/",id,".unfilt.vcf")
  }else{
      reg=paste0(" -L ",region)
      id=paste0(id,".",region)
      out_file=paste0(out_file_dir,"/",id,".unfilt.vcf")
  }

  if (is.vector(tumour)){
    tumour=paste0(" -I ",paste(normalizePath(tumour),collapse=" -I "))
  }else{
    tumour=paste0(" -I ",normalizePath(tumour))
  }
  norm=" "
  if (!is.null(normal)){
    if (is.vector(normal)){
      norm=paste0(" -I ",paste(normalizePath(normal),collapse=" -I ")," -normal ",
      paste(Vectorize(get_file_name)(normal),collapse=" -normal "))
  }else{
      norm=paste0(" -I ",normalizePath(normal)," -normal ",get_file_name(normal))
      }
  }

  if (pon!=""){
    pon=paste0(" --panel-of-normals ",pon)
  }

  f1r2=""
  if (orientation){
      out_file_dir_ort=set_dir(dir=output_dir,name="orientation")
      out_file2=paste0(out_file_dir_ort,"/",id,".f1r2.tar.gz")
      f1r2=paste0(" --f1r2-tar-gz ",out_file2)
  }



  filter_mnps=""
  if (mnps){
    filter_mnps=" -max-mnp-distance 0 "
  }

  exec_code=paste("singularity run ",sif_gatk,
  " /gatk/gatk   Mutect2 -R ",ref_genome,tumour, norm,
   " --germline-resource ",germ_resource, pon, " -O ",out_file, reg,f1r2,filter_mnps)

  if(mode=="batch"){
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
    stop("gatk failed to run due to unknown error.
    Check std error for more information.")
  }

  job_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    exec_code=exec_code, 
    task_id=task_id,
    input_args = argg,
    out_file_dir=out_file_dir,
    out_files=list(
      vcf=out_file,
      stats=paste0(out_file,".stats"),
      idx=paste0(out_file,".idx"),
      f1r2=out_file2)
  )




  if(wait&&mode=="batch"){
    job_validator(job=job_report$job_id,time=update_time,
    verbose=verbose,threads=threads)
  }

  return(job_report)

}



#' Multiregion parallelization across Mutect2 Gatk Variant Calling
#'
#' This function functions calls Mutect2 across multiple regions in parallel.
#' If a vector of tumour samples are provided these will be processed in multi-sample mode.
#' To run in tumour-normal mode suppply a single tumour and normal sample.
#' If no normal is supplied this will run in tumour only.
#' TO DO// Implement mitochondrial mode feature
#' 
#' For more information read:
#' https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2
#'
#' @param sif_gatk [REQUIRED] Path to gatk sif file.
#' @param bin_bcftools [REQUIRED] Path to bcftools binary file.
#' @param bin_bgzip [REQUIRED] Path to bgzip binary file.
#' @param bin_vep [REQUIRED] Path to VEP binary.
#' @param bin_tabix [REQUIRED] Path to tabix binary file.
#' @param bin_samtools [REQUIRED] Path to samtools binary file.
#' @param tumour [REQUIRED] Path to tumour BAM file.
#' @param normal [OPTIONAL] Path to normal BAM file.
#' @param ref_genome [REQUIRED] Path to reference genome fasta file.
#' @param germ_resource [REQUIRED]Path to germline resources vcf file.
#' @param regions [OPTIONAL] Regions to analyze. If regions for parallelization are not provided then these will be infered from BAM file.
#' @param output_name [OPTIONAL] Name for the output. If not given the name of one of the samples will be used.
#' @param pon [OPTIONAL] Path to panel of normal.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param mnps [OPTIONAL] Report MNPs in vcf file.
#' @param extract_pass [OPTIONAL] Extract variants with a PASS filter.Default TRUE
#' @param annotate [OPTIONAL] Annotate variants. Default TRUE
#' @param contamination [OPTIONAL] Produce sample cross-contamination reports. Default TRUE.
#' @param orientation [OPTIONAL] Produce read orientation inforamtion. Default FALSE
#' @param filter [OPTIONAL] Filter Mutect2. Default TRUE.
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





parallel_regions_mutect2_gatk=function(
  rdata=NULL,selected=NULL,
  sif_gatk=build_default_sif_list()$sif_gatk,
  bin_bcftools=build_default_tool_binary_list()$bin_bcftools,
  bin_vep=build_default_tool_binary_list()$bin_vep,
  bin_samtools=build_default_tool_binary_list()$bin_samtools,
  bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
  bin_tabix=build_default_tool_binary_list()$bin_tabix,
  tumour="",normal="",output_name="",
  chr=c(1:22,"X","Y","MT"),
  ref_genome=build_default_reference_list()$HG19$reference$genome,
  germ_resource=build_default_reference_list()$HG19$variant$germ_reference,
  biallelic_db=build_default_reference_list()$HG19$variant$biallelic_reference,
  db_interval=build_default_reference_list()$HG19$variant$biallelic_reference,
  regions=NULL,pon="",output_dir=".",
  verbose=FALSE,filter=TRUE,
  extract_pass=TRUE,annotate=TRUE,
  orientation=TRUE,mnps=FALSE,
  contamination=TRUE,clean=FALSE,
  batch_config=build_default_preprocess_config(),
  threads=4,ram=4,mode="local",
  executor_id=make_unique_id("parRegionMutect2"),
  task_name="parRegionMutect2",time="48:0:0",
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
  out_file_dir=set_dir(dir=output_dir,name=paste0("mutect2_reports/",id))
  tmp_dir=set_dir(dir=out_file_dir,name="tmp")

  job=build_job(executor_id=executor_id,task_id=task_id)


  jobs_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    exec_code=list(), 
    task_id=task_id,
    input_args=argg,
    out_file_dir=out_file_dir,
    out_files=list(
      )
    )

  if(is.null(regions)){

    jobs_report[["steps"]][["getChr"]] <- get_bam_reference_chr(
      bin_samtools=bin_samtools,
      bam=normalizePath(tumour[1]),verbose=verbose,output_dir=tmp_dir,
      output_name=paste0(get_file_name(tumour[1]),"_Ref"),
      executor_id=task_id,mode="local",threads=threads,ram=ram,
      time=time,update_time=update_time,wait=FALSE,hold=hold)


    regions=read.table(jobs_report[["steps"]][["getChr"]]$out_files$ref,
    sep="\t",header=TRUE)
    hold=jobs_report[["steps"]][["getChr"]]$job_id
  }

  regions$start=regions$start+1
  regions=regions %>% dplyr::mutate(region=paste0(chr,":",start,"-",end))
  
  
  ##### Ignore regions that are found in chromosome we are not interested
  
  if(!is.null(chr)){
    regions=regions[regions$chr %in% chr,]
  }
  

  region_list=regions$region
  names(region_list)=regions$region

  if(mode=="local"){
    jobs_report[["steps"]][["par_region_call_variants"]]<-
    parallel::mclapply(region_list,FUN=function(region){
      job_report <- mutect2_gatk(
            sif_gatk=sif_gatk,
            region=region,
            tumour=normalizePath(tumour),
            normal=normalizePath(normal),
            output_name=id,
            ref_genome=normalizePath(ref_genome),
            germ_resource = normalizePath(germ_resource),
            pon=normalizePath(pon),output_dir=tmp_dir,tmp_dir=tmp_dir,
            verbose=verbose,orientation=orientation,mnps=mnps,
            executor_id=task_id)
    },mc.cores=threads)
    
    }else if(mode=="batch"){
          rdata_file=paste0(tmp_dir,"/",job,".regions.RData")
          output_dir=tmp_dir
          executor_id=task_id
          output_name=id
          save(
            region_list,
            normalizePath(tumour),
            normalizePath(normal),
            sif_gatk,
            normalizePath(ref_genome),
            output_name,
            normalizePath(germ_resource),
            normalizePath(pon),
            orientation,
            mnps,
            executor_id,
            output_dir,
            verbose,tmp_dir,
            file = rdata_file)
          exec_code=paste0("Rscript -e \"ULPwgs::mutect2_gatk(rdata=\\\"",
          rdata_file,"\\\",selected=$SGE_TASK_ID)\"")
          out_file_dir2=set_dir(dir=out_file_dir,name="batch")
          batch_code=build_job_exec(job=job,time=time,ram=ram,
          threads=2,output_dir=out_file_dir2,
          hold=hold,array=length(region_list))
          exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)

          if(verbose){
              print_verbose(job=job,arg=argg,exec_code=exec_code)
          }
          error=execute_job(exec_code=exec_code)
          if(error!=0){
              stop("gatk failed to run due to unknown error.
              Check std error for more information.")
          }
        
         jobs_report[["steps"]][["par_region_call_variants"]]<- build_job_report(
              job_id=job,
              executor_id=executor_id,
              exec_code=exec_code, 
              task_id=task_id,
              input_args=argg,
              out_file_dir=tmp_dir,
              out_files=list(
                  vcf=paste0(tmp_dir,"/vcf/",id,".",region_list,".unfilt.vcf"),
                  stats=paste0(tmp_dir,"/vcf/",id,".",region_list,".unfilt.vcf.stats"),
                  idx=paste0(tmp_dir,"/vcf/",id,".",region_list,".unfilt.vcf.idx"),
                  f1r2=paste0(tmp_dir,"/orientation/",id,".",region_list,".f1r2.tar.gz")
                )
        )

  }

  vcfs=unlist_lvl(jobs_report[["steps"]][["par_region_call_variants"]],var="vcf")
  stats=unlist_lvl(jobs_report[["steps"]][["par_region_call_variants"]],var="stats")
  f1r2=unlist_lvl(jobs_report[["steps"]][["par_region_call_variants"]],var="f1r2")
  hold=unlist_lvl(jobs_report,var="job_id")

  jobs_report[["steps"]][["gatherFilesGatk"]]<-gather_mutect2_gatk(
      sif_gatk=sif_gatk,
      bin_samtools=bin_samtools,
      bin_bcftools=bin_bcftools,
      bin_bgzip=bin_bgzip,
      bin_tabix=bin_tabix,
      vcfs=normalizePath(vcfs),
      stats=normalizePath(stats),
      f1r2=normalizePath(f1r2),
      clean=clean,
      output_name=output_name,
      output_dir=out_file_dir,tmp_dir=tmp_dir,
      verbose=verbose,orientation=orientation,
      batch_config=batch_config,
      threads=threads,ram=ram,mode=mode,
      executor_id=task_id,
      time=time,
      hold=hold
  )



  vcf=unlist_lvl(jobs_report[["steps"]][["gatherFilesGatk"]],var="sorted_vcf")
  stats=unlist_lvl(jobs_report[["steps"]][["gatherFilesGatk"]],var="merged_stats")
  if(orientation){
    orientation_model=unlist_lvl(jobs_report[["steps"]][["gatherFilesGatk"]],var="orientation_model")
  }else{
    orientation_model=""
  }
 

  if(contamination){
          jobs_report[["steps"]][["estimateContaminationGatk"]]<- parallel_estimate_contamination_gatk(
                sif_gatk=sif_gatk,
                tumours=normalizePath(tumour),
                normal=normalizePath(normal),
                output_dir=out_file_dir,
                verbose=verbose,tmp_dir=tmp_dir,
                batch_config=batch_config,
                threads=threads,ram=ram,mode=mode,
                executor_id=task_id,
                time=time,
                hold=hold
          )

        contamination_table=unlist_lvl(jobs_report[["steps"]][["estimateContaminationGatk"]],var="contamination_table")
        segmentation_table=unlist_lvl(jobs_report[["steps"]][["estimateContaminationGatk"]],var="segmentation_table")
  
  }else{
    contamination_table=""
    segmentation_table=""
  }

  if(filter){
    jobs_report[["steps"]][["filterMutectGatk"]]<- mutect_filter_gatk(
            sif_gatk=sif_gatk,
            bin_bcftools=bin_bcftools,
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            vcf=normalizePath(vcf),stats=normalizePath(stats),contamination_table=normalizePath(contamination_table),
            segmentation_table=normalizePath(segmentation_table),
            orientation_model=normalizePath(orientation_model),output_name=output_name,
            ref_genome=normalizePath(ref_genome),
            extract_pass=extract_pass,
            output_dir=out_file_dir,
            verbose=verbose,clean=clean,
            batch_config=batch_config,
            threads=1,ram=ram,mode=mode,
            executor_id=task_id,
            time=time,
            hold=unlist_lvl(jobs_report[["steps"]],var="job_id")
        )
    if(extract_pass){
        vcf=unlist_lvl(jobs_report[["steps"]][["filterMutectGatk"]],var="extract_vcf")
    }else{
        vcf=unlist_lvl(jobs_report[["steps"]][["filterMutectGatk"]],var="filtered_vcf")
    }
  }

  if(annotate){
    jobs_report[["steps"]][["filterMutectGatk"]]<- annotate_vep(
        bin_vep=bin_vep,
        bin_bgzip=bin_bgzip,
        bin_tabix=bin_tabix,
        vcf=normalizePath(vcf),
        output_name=paste0(get_file_name(vcf),ifelse(extract_pass,".PASS","")),
        output_dir=out_file_dir,
        verbose=verbose,
        batch_config=batch_config,
        threads=threads,
        ram=ram,mode=mode,
        executor_id=task_id,
        time=time,
        hold=unlist_lvl(jobs_report[["steps"]],var="job_id")
    )
  }




  if(wait&&mode=="batch"){
    job_validator(job=unlist_lvl(jobs_report[["steps"]],var="job_id"),time=update_time,
    verbose=verbose,threads=threads)
  }

  return(jobs_report)

}


#' Multiregion parallelization across Mutect2 Gatk Variant Calling
#'
#' This function functions calls Mutect2 across multiple regions in parallel.
#' If a vector of tumour samples are provided these will be processed in multi-sample mode.
#' To run in tumour-normal mode suppply a single tumour and normal sample.
#' If no normal is supplied this will run in tumour only.
#' TO DO// Implement mitochondrial mode feature
#' 
#' For more information read:
#' https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2
#'
#' @param sif_gatk [REQUIRED] Path to gatk sif file.
#' @param bin_bcftools [REQUIRED] Path to bcftools binary file.
#' @param bin_vep [REQUIRED] Path to VEP binary.
#' @param bin_bgzip [REQUIRED] Path to bgzip binary file.
#' @param bin_tabix [REQUIRED] Path to tabix binary file.
#' @param bin_samtools [REQUIRED] Path to samtools binary file.
#' @param tumour [REQUIRED] Path to tumour BAM file.
#' @param normal [OPTIONAL] Path to normal BAM file.
#' @param ref_genome [REQUIRED] Path to reference genome fasta file.
#' @param chr [OPTIONAL] Chromosomes to analyze. Default c(1:22,"X","Y","MT")
#' @param germ_resource [REQUIRED]Path to germline resources vcf file.
#' @param regions [OPTIONAL] Regions to analyze. If regions for parallelization are not provided then these will be infered from BAM file.
#' @param output_name [OPTIONAL] Name for the output. If not given the name of one of the samples will be used.
#' @param pon [OPTIONAL] Path to panel of normal.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param mnps [OPTIONAL] Report MNPs in vcf file.
#' @param extract_pass [OPTIONAL] Extract variants with a PASS filter.Default TRUE
#' @param annotate [OPTIONAL] Annotate variants. Default TRUE
#' @param method [OPTIONAL] Default variant calling method. Default single. Options ["single","multi"]
#' @param contamination [OPTIONAL] Produce sample cross-contamination reports. Default TRUE.
#' @param orientation [OPTIONAL] Produce read orientation inforamtion. Default FALSE
#' @param filter [OPTIONAL] Filter Mutect2. Default TRUE.
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



parallel_samples_mutect2_gatk=function(
  sif_gatk=build_default_sif_list()$sif_gatk,
  bin_vep=build_default_tool_binary_list()$bin_vep,
  bin_bcftools=build_default_tool_binary_list()$bin_bcftools,
  bin_samtools=build_default_tool_binary_list()$bin_samtools,
  bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
  bin_tabix=build_default_tool_binary_list()$bin_tabix,
  tumour="",normal="",patient_id="",
  chr=c(1:22,"X","Y","MT"),
  ref_genome=build_default_reference_list()$HG19$reference$genome,
  germ_resource=build_default_reference_list()$HG19$variant$germ_reference,
  biallelic_db=build_default_reference_list()$HG19$variant$biallelic_reference,
  db_interval=build_default_reference_list()$HG19$variant$biallelic_reference,
  regions=NULL,pon="",output_dir=".",
  method="single",
  verbose=FALSE,filter=TRUE,
  orientation=TRUE,mnps=FALSE,
  annotate=TRUE,extract_pass=TRUE,
  contamination=TRUE,clean=FALSE,
  batch_config=build_default_preprocess_config(),
  threads=4,ram=4,mode="local",
  executor_id=make_unique_id("parSamplesMutect2"),
  task_name="parSamplesMutect2",time="48:0:0",
  update_time=60,wait=FALSE,hold=NULL
){

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir,name=patient_id)
  tmp_dir=set_dir(dir=out_file_dir,name="tmp")
 

  job=build_job(executor_id=executor_id,task_id=task_id)


  jobs_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    exec_code=list(), 
    task_id=task_id,
    input_args=argg,
    out_file_dir=out_file_dir,
    out_files=list(
      )
    )

  if(method=="single"){

    tumour_list=tumour
    names(tumour_list)=Vectorize(get_file_name)(tumour)

    if(mode=="local"){
      jobs_report[["steps"]][["par_sample_call_variants"]]<-
      parallel::mclapply(tumour_list,FUN=function(tumour){
        job_report <- parallel_regions_mutect2_gatk(
              sif_gatk=sif_gatk,
              bin_bcftools=bin_bcftools,
              bin_samtools=bin_samtools,
              bin_bgzip=bin_bgzip,
              bin_tabix=bin_tabix,
              regions=regions,
              output_name=get_file_name(tumour),
              chr=chr,
              ref_genome=normalizePath(ref_genome),
              germ_resource=normalizePath(germ_resource),
              biallelic_db=normalizePath(biallelic_db),
              db_interval=normalizePath(db_interval),
              pon=normalizePath(pon),
              filter=filter,
              orientation=orientation,
              mnps=mnps,
              contamination=contamination,
              clean=clean,
              tumour=tumour,
              normal=normal,
              output_dir=out_file_dir,
              verbose=verbose,
              threads=threads,
              executor_id=task_id)
      },mc.cores=threads)
    }else if(mode=="batch"){
            rdata_file=paste0(tmp_dir,"/",job,".samples.RData")
            output_dir=out_file_dir
            executor_id=task_id
            save(tumour_list,normal,
              sif_gatk,
              bin_bcftools,
              bin_samtools,
              bin_bgzip,
              bin_tabix,
              chr,
              germ_resource,
              biallelic_db,db_interval,pon,ref_genome,
              filter,orientation,mnps,
              mode,ram,
              executor_id,
              output_dir,
              contamination,clean,
              verbose,file = rdata_file)
            exec_code=paste0("Rscript -e \"ULPwgs::parallel_regions_mutect2_gatk(rdata=\\\"",
            rdata_file,"\\\",selected=$SGE_TASK_ID)\"")
            out_file_dir2=set_dir(dir=out_file_dir,name="batch")
            batch_code=build_job_exec(job=job,time=time,ram=ram,
            threads=4,output_dir=out_file_dir2,
            hold=hold,array=length(tumour_list))
            exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)

            if(verbose){
                print_verbose(job=job,arg=argg,exec_code=exec_code)
            }
            error=execute_job(exec_code=exec_code)
            if(error!=0){
                stop("gatk failed to run due to unknown error.
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
                  filtered_vcf=ifelse(filter,paste0(out_file_dir,"/mutect2_reports/",names(tumour_list),"/",
                  names(tumour_list),".filtered.vcf"),""),
                  sorted_vcf=paste0(out_file_dir,"/mutect2_reports/",names(tumour_list),"/",
                  names(tumour_list),".sorted.vcf"),
                  compressed_vcf=paste0(out_file_dir,"/mutect2_reports/",names(tumour_list),"/",
                  names(tumour_list),".sorted.vcf.gz")
                )
                  )
             }
    }else if (method=="multi"){
        jobs_report[["steps"]][["par_sample_call_variants"]]<-
              parallel_regions_mutect2_gatk(
                sif_gatk=sif_gatk,
                bin_bcftools=bin_bcftools,
                bin_samtools=bin_samtools,
                bin_bgzip=bin_bgzip,
                bin_tabix=bin_tabix,
                regions=regions,
                output_name=patient_id,
                chr=chr,
                ref_genome=normalizePath(ref_genome),
                germ_resource=normalizePath(germ_resource),
                biallelic_db=normalizePath(biallelic_db),
                db_interval=normalizePath(db_interval),
                pon=normalizePath(pon),
                filter=filter,
                orientation=orientation,
                mnps=mnps,
                contamination=contamination,
                clean=clean,
                tumour=tumour,
                normal=normal,
                output_dir=out_file_dir,
                verbose=verbose,
                threads=threads,
                executor_id=task_id,
                mode=mode
              )

    }else{
      stop("Wrong method supplied. Only single or multi methods available.")
    }     


  if(wait&&mode=="batch"){
    job_validator(job=unlist_lvl(jobs_report[["steps"]],var="job_id"),time=update_time,
    verbose=verbose,threads=threads)
  }

  return(jobs_report)

}




#' Multiregion parallelization across Mutect2 Gatk Variant Calling for multiple samples
#'
#' This function functions calls Mutect2 across multiple regions in parallel.
#' If a vector of tumour samples are provided these will be processed in multi-sample mode.
#' To run in tumour-normal mode suppply a single tumour and normal sample.
#' If no normal is supplied this will run in tumour only.
#' TO DO// Implement mitochondrial mode feature
#' 
#' For more information read:
#' https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2
#'
#' @param sif_gatk [REQUIRED] Path to gatk sif file.
#' @param bin_bcftools [REQUIRED] Path to bcftools binary file.
#' @param bin_vep [REQUIRED] Path to VEP binary.
#' @param bin_bgzip [REQUIRED] Path to bgzip binary file.
#' @param bin_tabix [REQUIRED] Path to tabix binary file.
#' @param bin_samtools [REQUIRED] Path to samtools binary file.
#' @param sample_sheet [OPTIONAL] Path to sheet with sample information.
#' @param bam_dir [OPTIONAL] Path to directory with BAM files.
#' @param normal_id [OPTIONAL] Path to directory with BAM files.
#' @param chr [OPTIONAL] Chromosomes to analyze. Default c(1:22,"X","Y","MT")
#' @param ref_genome [REQUIRED] Path to reference genome fasta file.
#' @param germ_resource [REQUIRED]Path to germline resources vcf file.
#' @param regions [OPTIONAL] Regions to analyze. If regions for parallelization are not provided then these will be infered from BAM file.
#' @param output_name [OPTIONAL] Name for the output. If not given the name of one of the samples will be used.
#' @param pon [OPTIONAL] Path to panel of normal.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param mnps [OPTIONAL] Report MNPs in vcf file.
#' @param method [OPTIONAL] Default variant calling method. Default single. Options ["single","multi"]
#' @param contamination [OPTIONAL] Produce sample cross-contamination reports. Default TRUE.
#' @param orientation [OPTIONAL] Produce read orientation inforamtion. Default FALSE
#' @param filter [OPTIONAL] Filter Mutect2. Default TRUE.
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



multisample_mutect2_gatk=function(
  sif_gatk=build_default_sif_list()$sif_gatk,
  bin_vep=build_default_tool_binary_list()$bin_vep,
  bin_bcftools=build_default_tool_binary_list()$bin_bcftools,
  bin_samtools=build_default_tool_binary_list()$bin_samtools,
  bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
  bin_tabix=build_default_tool_binary_list()$bin_tabix,
  sample_sheet=NULL,
  bam_dir="",
  normal_id="",
  patient_id="",
  pattern="bam$",
  ref_genome=build_default_reference_list()$HG19$reference$genome,
  germ_resource=build_default_reference_list()$HG19$variant$germ_reference,
  biallelic_db=build_default_reference_list()$HG19$variant$biallelic_reference,
  db_interval=build_default_reference_list()$HG19$variant$biallelic_reference,
  pon=build_default_reference_list()$HG19$panel$PCF_V3$variant$pon_muts,
  chr=c(1:22,"X","Y","MT"),
  method="single",
  regions=NULL,output_dir=".",
  verbose=FALSE,filter=TRUE,
  annotate=TRUE,extract_pass=TRUE,
  orientation=TRUE,mnps=FALSE,
  contamination=TRUE,clean=FALSE,
  header=TRUE,
  sep="\t",
  batch_config=build_default_preprocess_config(),
  threads=4,ram=8,mode="local",
  executor_id=make_unique_id("multiSampleMutect2"),
  task_name="multiSampleMutect2",time="48:0:0",
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
      "sif_gatk",
      "bin_bcftools",
      "bin_samtools",
      "bin_bgzip",
      "bin_vep",
      "bin_tabix",
      "patient_id",
      "tumour",
      "normal",
      "chr",
      "ref_genome",
      "germ_resource",
      "biallelic_db",
      "db_interval",
      "pon",
      "regions",
      "filter",
      "orientation",
      "annotate",
      "extract_pass",
      "mnps",
      "method",
      "contamination",
      "output_dir",
      "clean",
      "verbose",
      "batch_config",
      "threads",
      "ram",
      "time",
      "mode",
      "hold")


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

        job_report[["steps"]][["multisample_gatk"]]=parallel::mclapply(seq(1,nrow(file_info)),FUN=function(x){
            
            lapply(columns,FUN=function(col){
                if(is.null(file_info[[col]])){
                    file_info[[col]]<<-get(col)
                }

            
            })

       
            job_report<- parallel_samples_mutect2_gatk(
              sif_gatk=sif_gatk,
              bin_vep=bin_vep,
              bin_bcftools=bin_bcftools,
              bin_samtools=bin_samtools,
              bin_bgzip=bin_bgzip,
              bin_tabix=bin_tabix,
              tumour=normalizePath(unlist(file_info[x,]$tumour)),
              normal=normalizePath(file_info[x,]$normal),
              ref_genome=normalizePath(file_info[x,]$ref_genome),
              germ_resource=normalizePath(file_info[x,]$germ_resource),
              biallelic_db=normalizePath(file_info[x,]$biallelic_db),
              db_interval=normalizePath(file_info[x,]$biallelic_db),
              patient_id=file_info[x,]$patient_id,
              pon=normalizePath(file_info[x,]$pon),
              chr=eval(parse(text=file_info[x,]$chr)),
              regions=file_info[[x,"regions"]],
              method=file_info[x,]$method,
              output_dir=file_info[x,]$output_dir,
              extract_pass=file_info[x,]$extract_pass,
              annotate=file_info[x,]$annotate,
              verbose=file_info[x,]$verbose,
              filter=file_info[x,]$filter,
              orientation=file_info[x,]$orientation,mnps=file_info[x,]$mnps,
              contamination=file_info[x,]$contamination,
              clean=file_info[x,]$clean,
              batch_config=file_info[x,]$batch_config,
              threads=file_info[x,]$threads,
              ram=file_info[x,]$ram,
              mode=file_info[x,]$mode,
              executor_id=task_id,
              time=file_info[x,]$time,
              hold=file_info[[x,"hold"]])
            },mc.cores=ifelse(mode=="local",1,3))

    }else{
        bam_dir_path=system(paste("realpath",bam_dir),intern=TRUE)
        bam_files=system(paste0("find ",bam_dir_path,"| grep ",pattern),intern=TRUE)
        tumour=bam_files[!grepl(normal_id,bam_files)]
        normal=bam_files[grepl(normal_id,bam_files)]

      
          job_report<-parallel_samples_mutect2_gatk(
                  sif_gatk=sif_gatk,
                  bin_vep=bin_vep,
                  bin_bcftools=bin_bcftools,
                  bin_samtools=bin_samtools,
                  bin_bgzip=bin_bgzip,
                  bin_tabix=bin_tabix,
                  tumour=normalizePath(tumour),
                  normal=normalizePath(normal),
                  patient_id=patient_id,
                  ref_genome=normalizePath(ref_genome),
                  germ_resource=normalizePath(germ_resource),
                  biallelic_db=normalizePath(biallelic_db),
                  db_interval=normalizePath(biallelic_db),
                  extract_pass=extract_pass,
                  chr=chr,
                  annotate=annotate,
                  regions=regions,
                  method=method,
                  pon=pon,
                  output_dir=output_dir,
                  verbose=verbose,
                  filter=filter,
                  orientation=orientation,mnps=mnps,
                  contamination=contamination,
                  clean=clean,
                  batch_config=batch_config,
                  threads=threads,ram=ram,
                  mode=mode,
                  executor_id=task_id,
                  time=time,
                  hold=hold)
    }


    if(wait&&mode=="batch"){
        job_validator(job=unlist_level(named_list=job_report[["steps"]][["multisample_gatk"]],var="job_id"),
        time=update_time,verbose=verbose,threads=threads)
    }

    return(job_report)

}



#' Filter Mutect2 Calls wrapper for R
#'
#' This function filters Mutect2 callsets under multiple conditions
#'
#' For more information: 
#' 
#' https://gatk.broadinstitute.org/hc/en-us/articles/360036856831-FilterMutectCalls
#' 
#' 
#' @param sif_gatk [REQUIRED] Path to gatk sif file.
#' @param bin_bcftools [REQUIRED] Path to bcftools binary file.
#' @param bin_bgzip [REQUIRED] Path to bgzip binary file.
#' @param bin_tabix [REQUIRED] Path to tabix binary file.
#' @param vcf [REQUIRED] Path to VCF file.
#' @param stats [OPTIONAL] Path to stats information generated by Mutect2 for supplied VCF.
#' @param contamination_table [OPTIONAL] Path to contamination tables for each tumour samples in VCF.
#' @param segmentation_table [OPTIONAL] Path to segmentation tables for each tumour samples in VCF.
#' @param orientation_model [OPTIONAL] Path to orientation model generated F1R2 read information.
#' @param clean [OPTIONAL] Remove unfiltered VCF after completion. Default FALSE.
#' @param extract_pass [OPTIONAL] Extract PASSing variants. Default TRUE.
#' @param output_name [OPTIONAL] Name for the output. If not given the name of one of the samples will be used.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param threads [OPTIONAL] Number of threads to split the work. Default 4
#' @param ram [OPTIONAL] RAM memory to asing to each thread. Default 4
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param mode [OPTIONAL]  Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id [OPTIONAL] Task EXECUTOR ID. Default "recalCovariates"
#' @param task_name [OPTIONAL] Task name. Default "recalCovariates"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export


mutect_filter_gatk=function(
  sif_gatk=build_default_sif_list()$sif_gatk,
  bin_bcftools=build_default_tool_binary_list()$bin_bcftools,
  bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
  bin_tabix=build_default_tool_binary_list()$bin_tabix,
  vcf="",stats=NULL,contamination_table=NULL,
  extract_pass=TRUE,
  segmentation_table=NULL,orientation_model=NULL,output_name="",
  ref_genome=build_default_reference_list()$HG19$reference$genome,
  output_dir=".",verbose=FALSE,clean=FALSE,
  batch_config=build_default_preprocess_config(),
  threads=4,ram=4,mode="local",
  executor_id=make_unique_id("filterMutect2Gatk"),
  task_name="filterMutect2Gatk",time="48:0:0",
  update_time=60,wait=FALSE,hold=NULL
){

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir)
  job=build_job(executor_id=executor_id,task_id=task_id)


  id=""
  if(output_name!=""){
    id=output_name
  }else{
     id=get_file_name(vcf)
  }

  if(!is.null(stats)){
      stats=paste0(" -stats ",normalizePath(stats))
  }

  if (!is.null(contamination_table)){
     contamination_table=paste0(" --contamination-table ", paste0(normalizePath(contamination_table),
     collapse=" --contamination-table "))
  }
  
  if(!is.null(segmentation_table)){
      segmentation_table=paste0(" --tumor-segmentation ",paste0(normalizePath(segmentation_table),
      collapse=" --tumor-segmentation "))
    }

  if(!is.null(orientation_model)){
      orientation_model=paste0(" --ob-priors ",normalizePath(orientation_model))
  }
  
  out_file=paste0(out_file_dir,"/",id,".filtered.vcf")


  exec_code=paste0("singularity exec -H ",getwd(),":/home ",sif_gatk,
  " /gatk/gatk   FilterMutectCalls -R ",ref_genome," -O ",out_file," -V ",normalizePath(vcf), stats,
  contamination_table,segmentation_table,orientation_model
  )

  if(clean){
    exec_code=paste(exec_code," && rm",vcf)
  }


 if(mode=="batch"){
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
    stop("gatk failed to run due to unknown error.
    Check std error for more information.")
  }




  job_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    exec_code=exec_code, 
    task_id=task_id,
    input_args = argg,
    out_file_dir=out_file_dir,
    out_files=list(
      filtered_vcf=out_file)
  )


    if(extract_pass){
      job_report[["steps"]][["extractPASSvcf"]]<-
        parallel_vcfs_variants_by_filters_vcf(
          bin_bgzip=bin_bgzip,
          bin_tabix=bin_tabix,
          vcf=out_file,filters="PASS",
          exclusive=TRUE,
          compress=FALSE,
          output_dir=out_file_dir,
          verbose=verbose,
          batch_config=batch_config,
          threads=threads,ram=ram,mode=mode,
          executor_id=task_id,
          time=time,
          hold=job
        )
    }

  if(wait&&mode=="batch"){
    job_validator(job=job_report$job_id,time=update_time,
    verbose=verbose,threads=threads)
  }

  return(job_report)

}




gather_mutect2_gatk=function(
  sif_gatk=build_default_sif_list()$sif_gatk,
  bin_samtools=build_default_tool_binary_list()$bin_samtools,
  bin_bcftools=build_default_tool_binary_list()$bin_bcftools,
  bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
  bin_tabix=build_default_tool_binary_list()$bin_tabix,
  vcfs="",stats=NULL,f1r2=NULL,output_dir=".",tmp_dir=".",
  output_name="",verbose=FALSE,orientation=FALSE,clean=TRUE,
  batch_config=build_default_preprocess_config(),
  threads=4,ram=4,mode="local",
  executor_id=make_unique_id("gatherMutect"),
  task_name="gatherMutect",time="48:0:0",
  update_time=60,wait=FALSE,hold=NULL
){

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir)
  job=build_job(executor_id=executor_id,task_id=task_id)

  jobs_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    exec_code=list(), 
    task_id=task_id,
    input_args=argg,
    out_file_dir=out_file_dir,
    out_files=list(
      )
    )

  jobs_report[["steps"]][["concatVCF"]] <- concat_vcf(
    bin_bcftools=bin_bcftools,
    bin_bgzip=bin_bgzip,
    bin_tabix=bin_tabix,
    vcfs=normalizePath(vcfs),clean=clean,sort=TRUE,
    compress=FALSE,format="vcf",
    verbose=verbose,output_name=output_name,
    output_dir=out_file_dir,
    tmp_dir=tmp_dir,batch_config=batch_config,
    threads=1,ram=ram,mode=mode,
    executor_id=task_id,
    time=time,hold=hold
  )

  jobs_report[["steps"]][["compressAndIndex"]]<-compress_and_index_vcf_htslib(
    bin_bgzip=bin_bgzip,
    bin_tabix=bin_tabix,
    vcf=unlist_lvl(jobs_report[["steps"]][["concatVCF"]],var="sorted_vcf"),
    output_dir=out_file_dir,output_name=output_name,
    clean=clean,verbose=verbose,
    batch_config=batch_config,
    threads=1,ram=ram,mode=mode,
    executor_id=task_id,
    time=time,
    hold=unlist_lvl(jobs_report[["steps"]][["concatVCF"]],var="job_id")
  )



  jobs_report[["steps"]][["mergeStatMutect"]] <- merge_mutect_stats_gatk(
    sif_gatk=sif_gatk,
    stats=stats,output_name=output_name,
    output_dir=out_file_dir,clean=clean,
    verbose=verbose,batch_config=batch_config,
    threads=1,ram=ram,mode=mode,
    executor_id=task_id,
    time=time,hold=hold
  )

  if(orientation){
      jobs_report[["steps"]][["orientationModel"]] <- learn_orientation_gatk(
        sif_gatk=sif_gatk,
        f1r2=f1r2,output_name=output_name,
        output_dir=out_file_dir,clean=clean,
        verbose=verbose,batch_config=batch_config,
        threads=1,ram=4,mode=mode,
        executor_id=task_id,
        time=time,hold=hold
      )
  }

  if(wait&&mode=="batch"){
    job_validator(job=jobs_report$job_id,time=update_time,
    verbose=verbose,threads=threads)
  }

  return(jobs_report)

}   

#' Create Read Orientation Model for Mutect2 Filter
#'
#' This function creates read orientation model for mutect2 filtering
#'
#' @param f1r2 Path to f1r2 files.
#' @param sif_gatk Path to gatk sif file.
#' @param ref_genome Path to reference genome fasta file.
#' @param output_name [OPTIONAL] Name for the output. If not given the name of one of the samples will be used.
#' @param output_dir Path to the output directory.
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


learn_orientation_gatk=function(
  sif_gatk=build_default_sif_list()$sif_gatk,
  f1r2="",output_name="",output_dir=".",clean=TRUE,
  verbose=FALSE,batch_config=build_default_preprocess_config(),
  threads=4,ram=4,mode="local",
  executor_id=make_unique_id("learnOrientationMutect2"),
  task_name="learnOrientationMutect2",time="48:0:0",
  update_time=60,wait=FALSE,hold=NULL
){

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir)
  job=build_job(executor_id=executor_id,task_id=task_id)


  id=""
  if(output_name!=""){
    id=output_name
  }else{
     id=get_file_name(f1r2[1])
  }
 
  if (!is.null(f1r2)){
    f1r2_list=paste0(" -I ",paste0(normalizePath(f1r2),collapse=" -I "))
  }
  
  out_file=paste0(out_file_dir,"/",id,".ROM.tar.gz")

  exec_code=paste0("singularity exec -H ",getwd(),":/home ",sif_gatk,
  " /gatk/gatk   LearnReadOrientationModel  ",f1r2_list," -O ",out_file)


  if(clean){
    exec_code=paste(exec_code," && rm",paste(f1r2,collapse=" "))
  }

  if(mode=="batch"){
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
    stop("gatk failed to run due to unknown error.
    Check std error for more information.")
  }

  job_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    exec_code=exec_code, 
    task_id=task_id,
    input_args = argg,
    out_file_dir=out_file_dir,
    out_files=list(
      orientation_model=out_file)
  )


  if(wait&&mode=="batch"){
    job_validator(job=job_report$job_id,time=update_time,
    verbose=verbose,threads=threads)
  }

  return(job_report)



}




#' Estimate sample contamination Gatk
#'
#' This function estimates the sample cross-contamination
#'
#' @param f1r2 Path to f1r2 files.
#' @param sif_gatk Path to gatk sif file.
#' @param ref_genome Path to reference genome fasta file.
#' @param output_name [OPTIONAL] Name for the output. If not given the name of one of the samples will be used.
#' @param output_dir Path to the output directory.
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


estimate_contamination_gatk=function(
  sif_gatk=build_default_sif_list()$sif_gatk,
  rdata=NULL,selected=NULL,
  tumour=NA,normal=NA,tumour_pileup="",
  normal_pileup="",output_name="",output_dir=".",tmp_dir=".",
  biallelic_db=build_default_reference_list()$HG19$variant$biallelic_reference,
  db_interval=build_default_reference_list()$HG19$variant$biallelic_reference,
  verbose=FALSE,batch_config=build_default_preprocess_config(),
  threads=1,ram=4,mode="local",
  executor_id=make_unique_id("estimateContaminationGatk"),
  task_name="estimateContaminationGatk",time="48:0:0",
  update_time=60,wait=FALSE,hold=NULL){

  if(!is.null(rdata)){
    load(rdata)
    if(!is.null(selected)){
      tumour=tumour_list[selected]
    }
  }

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir,name="contamination_reports")
  out_file_dir_tab=set_dir(dir=out_file_dir,name="contamination")
  out_file_dir_seg=set_dir(dir=out_file_dir,name="segmentation")
  job=build_job(executor_id=executor_id,task_id=task_id)


  if(output_name!=""){
    id=output_name
  }else if(tumour!=""){
    id=get_file_name(tumour)
  }else if(tumour_pileup!=""){
    id=get_file_name(tumour_pileup)
  }

  out_file=paste0(out_file_dir_tab,"/",id,".contamination.table")
  out_file2=paste0(out_file_dir_seg,"/",id,".segmentation.table")

  jobs_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    task_id=task_id,
    input_args = argg,
    out_file_dir=list(
      contamination=out_file_dir_tab,
      segmentation=out_file_dir_seg
      ),
    out_files=list(
      contamination_table=out_file,
      segmentation_table=out_file2
    )
  )

  if(!is.na(normal)){
      jobs_report[["steps"]][["nPileupGatk"]]<-pileup_summary_gatk(
        sif_gatk=sif_gatk,
        bam=normalizePath(normal),
        output_name=get_file_name(normal),
        output_dir=paste0(out_file_dir,"/pileup_reports/normal"),
        verbose=verbose,batch_config=batch_config,
        biallelic_db=normalizePath(biallelic_db),
        db_interval=normalizePath(db_interval),
        threads=1,ram=ram,mode=mode,
        executor_id=task_id,hold=hold
      )
      normal_pileup<-jobs_report[["steps"]][["nPileupGatk"]]$out_files$pileup_table

  }

  if(!is.na(tumour)){
        jobs_report[["steps"]][["tPileupGatk"]]<-pileup_summary_gatk(
          sif_gatk=sif_gatk,
          bam=normalizePath(tumour),
          output_name=get_file_name(tumour),
          output_dir=paste0(out_file_dir,"/pileup_reports/tumour"),
          verbose=verbose,batch_config=batch_config,
          biallelic_db=normalizePath(biallelic_db),
          db_interval=normalizePath(db_interval),
          threads=1,ram=ram,mode=mode,
          executor_id=task_id,hold=hold
        )
        tumour_pileup<-jobs_report[["steps"]][["tPileupGatk"]]$out_files$pileup_table
  }


  
  if (normal_pileup!=""){
        normal_pileup=paste0(" -matched ",normalizePath(normal_pileup))
  }




  exec_code=paste0("singularity exec -H ",getwd(),":/home ",sif_gatk,
  " /gatk/gatk   CalculateContamination  -O ",out_file, " -tumor-segmentation ", out_file2,
  " -I ",normalizePath(tumour_pileup),normal_pileup)

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





#' Estimate multis-sample contamination GATK
#'
#' This function estimates the sample cross-contamination
#'
#' @param f1r2 Path to f1r2 files.
#' @param sif_gatk Path to gatk sif file.
#' @param ref_genome Path to reference genome fasta file.
#' @param output_name [OPTIONAL] Name for the output. If not given the name of one of the samples will be used.
#' @param output_dir Path to the output directory.
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


parallel_estimate_contamination_gatk=function(
  sif_gatk=build_default_sif_list()$sif_gatk,
  tumours="",normal="",output_dir=".",tmp_dir=".",
  biallelic_db=build_default_reference_list()$HG19$variant$biallelic_reference,
  db_interval=build_default_reference_list()$HG19$variant$biallelic_reference,
  verbose=FALSE,batch_config=build_default_preprocess_config(),
  threads=1,ram=4,mode="local",
  executor_id=make_unique_id("parEstimateContaminationGatk"),
  task_name="parEstimateContaminationGatk",time="48:0:0",
  update_time=60,wait=FALSE,hold=NULL
){

    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir)
    job=build_job(executor_id=executor_id,task_id=task_id)

      
    jobs_report=build_job_report(
      job_id=job,
      executor_id=executor_id,
      exec_code=list(), 
      task_id=task_id,
      input_args=argg,
      out_file_dir=out_file_dir,
      out_files=list(
        )
      )


    tumour_list=tumours
    names(tumour_list)=Vectorize(get_file_name)(tumours)
      
    if(normal!=""){
        jobs_report[["steps"]][["nPileupGatk"]]<-pileup_summary_gatk(
          sif_gatk=sif_gatk,
          bam=normal,output_name=get_file_name(normal),
          output_dir=paste0(out_file_dir,"/contamination_reports/pileup_reports/normal"),
          verbose=verbose,batch_config=batch_config,
          biallelic_db=biallelic_db,
          db_interval=db_interval,
          threads=1,ram=ram,mode=mode,
          executor_id=task_id,hold=hold
        )
        normal_pileup<-jobs_report[["steps"]][["nPileupGatk"]]$out_files$pileup_table
        hold=jobs_report[["steps"]][["nPileupGatk"]]$job_id
    }


    if(mode=="local"){
      jobs_report[["steps"]][["parEstimateContaminationGatk"]]<-
      parallel::mclapply(tumour_list,FUN=function(tumour){
        job_report <-  estimate_contamination_gatk(
            sif_gatk= sif_gatk,
            tumour = normalizePath(tumour),
            normal="",
            normal_pileup= normal_pileup,
            output_name=get_file_name(tumour),
            output_dir=out_file_dir,
            tmp_dir=tmp_dir,
            verbose=verbose,
            executor_id=task_id)
      },mc.cores=threads)
      }else if(mode=="batch"){
            rdata_file=paste0(tmp_dir,"/",job,".pileups.RData")
            save(tumour_list,normal_pileup,sif_gatk,biallelic_db,db_interval,
            output_dir,verbose,tmp_dir,file = rdata_file)
            exec_code=paste0("Rscript -e \"ULPwgs::estimate_contamination_gatk(rdata=\\\"",
            rdata_file,"\\\",selected=$SGE_TASK_ID)\"")
            out_file_dir2=set_dir(dir=out_file_dir,name="batch")
            batch_code=build_job_exec(job=job,time=time,ram=ram,
            threads=1,output_dir=out_file_dir2,
            hold=hold,array=length(tumour_list))
            exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)

            if(verbose){
                print_verbose(job=job,arg=argg,exec_code=exec_code)
            }
            error=execute_job(exec_code=exec_code)
            if(error!=0){
                stop("gatk failed to run due to unknown error.
                Check std error for more information.")
            }
          
            jobs_report[["steps"]][["parEstimateContaminationGatk"]]<- build_job_report(
                  job_id=job,
                  executor_id=executor_id,
                  exec_code=exec_code, 
                  task_id=task_id,
                  input_args=argg,
                  out_file_dir=out_file_dir,
                  out_files=list(
                      contamination_table=paste0(out_file_dir,"/contamination_reports/contamination/",names(tumour_list),".contamination.table"),
                      segmentation_table=paste0(out_file_dir,"/contamination_reports/segmentation/",names(tumour_list),".segmentation.table")
                    )
            )
    }



    if(wait&&mode=="batch"){
      jobs_validator(job=jobs_report$job_id,time=update_time,
      verbose=verbose,threads=threads)
    }

    return(jobs_report)

}





#' Generate a pileup summary for BAM file
#'
#' This function generates a pilelup summary for a bam file
#'
#' @param bam Path to a BAM file
#' @param sif_gatk Path to gatk sif file.
#' @param output_name [OPTIONAL] Name for the output. If not given the name of one of the samples will be used.
#' @param output_dir Path to the output directory.
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


pileup_summary_gatk=function(
  sif_gatk=build_default_sif_list()$sif_gatk,
  bam="",output_name="",output_dir=".",
  rdata=NULL,selected=NULL,
  verbose=FALSE,batch_config=build_default_preprocess_config(),
  biallelic_db=build_default_reference_list()$HG19$variant$biallelic_reference,
  db_interval=build_default_reference_list()$HG19$variant$biallelic_reference,
  threads=1,ram=4,mode="local",
  executor_id=make_unique_id("pileupSummaryGatk"),
  task_name="pileupSummaryGatk",time="48:0:0",
  update_time=60,wait=FALSE,hold=NULL
){
  
  if(!is.null(rdata)){
    load(rdata)
  if(!is.null(selected)){
    bam=bam_list[selected]
  }
  }

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir)
  job=build_job(executor_id=executor_id,task_id=task_id)


  id=""
  if(output_name!=""){
    id=output_name
  }else{
     id=get_file_name(bam)
  }

  out_file=paste0(out_file_dir,"/",id,".pileup.table")

  exec_code=paste0("singularity exec -H ",getwd(),":/home ",sif_gatk,
  " /gatk/gatk  GetPileupSummaries  -O ",out_file," -V ",normalizePath(biallelic_db)," -L ",
   normalizePath(db_interval), " -I ",normalizePath(bam))


  if(mode=="batch"){
       out_file_dir2=set_dir(dir=out_file_dir,name="batch")
       batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
       threads=threads,output_dir=out_file_dir2)
       exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",
       exec_code,"'|",batch_code)
  }

  if(verbose){
       print_verbose(job=job,arg=argg,exec_code=exec_code)
  }

  error=execute_job(exec_code=exec_code)
  
  if(error!=0){
    stop("gatk failed to run due to unknown error.
    Check std error for more information.")
  }

  job_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    exec_code=exec_code, 
    task_id=task_id,
    input_args = argg,
    out_file_dir=out_file_dir,
    out_files=list(
      pileup_table=out_file)
  )


  if(wait&&mode=="batch"){
    job_validator(job=job_report$job_id,time=update_time,
    verbose=verbose,threads=threads)
  }

  return(job_report)

}



#' Generate pileup summary for BAM files 
#'
#' This function generates a pilelup summary for a bam file
#'
#' @param bam Path to a BAM file
#' @param sif_gatk Path to gatk sif file.
#' @param output_name [OPTIONAL] Name for the output. If not given the name of one of the samples will be used.
#' @param output_dir Path to the output directory.
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


parallel_pileup_summary_gatk=function(
  sif_gatk=build_default_sif_list()$sif_gatk,
  bams="",output_dir=".",tmp_dir=".",verbose=FALSE,
  batch_config=build_default_preprocess_config(),
  biallelic_db=build_default_reference_list()$HG19$variant$biallelic_reference,
  db_interval=build_default_reference_list()$HG19$variant$biallelic_reference,
  threads=1,ram=4,mode="local",
  executor_id=make_unique_id("parPileupSummaryGatk"),
  task_name="parPileupSummaryGatk",time="48:0:0",
  update_time=60,wait=FALSE,hold=NULL){


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

  bam_list=bams
  names(bam_list)=Vectorize(get_file_name)(bams)



  if(mode=="local"){
    job_report[["steps"]][["parPileupSummaryGatk"]]<-
    parallel::mclapply(bam_list,FUN=function(bam){
      job_report <-  pileup_summary_gatk(
          sif_gatk=sif_gatk,
          bam=normalizePath(bam),output_name=get_file_name(bam),output_dir=out_file_dir,
          verbose=verbose,batch_config=batch_config,
          biallelic_db=normalizePath(biallelic_db),
          db_interval=normalizePath(db_interval),
          executor_id=task_id,
          time=time,hold=hold)
    },mc.cores=threads)
    }else if(mode=="batch"){
          rdata_file=paste0(tmp_dir,"/",job,".tumours.RData")
          output_dir=out_file_dir
          save(bam_list,sif_gatk,output_dir,normalizePath(biallelic_db),normalizePath(db_interval),verbose,tmp_dir,file = rdata_file)
          exec_code=paste0("Rscript -e \"ULPwgs::pileup_summary_gatk(rdata=\\\"",rdata_file,"\\\",selected=$SGE_TASK_ID)\"")
          out_file_dir2=set_dir(dir=out_file_dir,name="batch")
          batch_code=build_job_exec(job=job,time=time,ram=ram,
          threads=1,output_dir=out_file_dir2,
          hold=hold,array=length(bam_list))
          exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)

          if(verbose){
              print_verbose(job=job,arg=argg,exec_code=exec_code)
          }
          error=execute_job(exec_code=exec_code)
          if(error!=0){
              stop("gatk failed to run due to unknown error.
              Check std error for more information.")
          }
        
          job_report[["steps"]][["parPileupSummaryGatk"]]<- build_job_report(
                job_id=job,
                executor_id=executor_id,
                exec_code=exec_code, 
                task_id=task_id,
                input_args=argg,
                out_file_dir=out_file_dir,
                out_files=list(
                    pileup_table=paste0(out_file_dir,"/pileup/",names(bam_list),".pileup.table")
                  )
          )
  }


  if(wait&&mode=="batch"){
    job_validator(job=job_report$job_id,time=update_time,
    verbose=verbose,threads=threads)
  }

  return(job_report)

}





#' Multiregion parallelization across Mutect2 Gatk Variant Calling
#'
#' This function functions calls Mutect2 across multiple regions in parallel
#' 
#' @param sif_gatk Path to gatk sif file.
#' @param stats Path to Mutect2 stats files.
#' @param output_name [OPTIONAL] Name for the output. If not given the name of one of the samples will be used.
#' @param output_dir Path to the output directory.
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


merge_mutect_stats_gatk=function(
  sif_gatk=build_default_sif_list()$sif_gatk,
  stats=NA,output_name="",output_dir=".",clean=FALSE,
  verbose=TRUE, batch_config=build_default_preprocess_config(),
  threads=1,ram=4,mode="local",
  executor_id=make_unique_id("mergeMutectStats"),
  task_name="mergeMutectStats",time="48:0:0",
  update_time=60,wait=FALSE,hold=NULL){

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir)
  job=build_job(executor_id=executor_id,
  task_id=task_id)

  id=""
  if(output_name!=""){
    id=output_name
  }else{
     id=get_file_name(stats[1])
  }
 
  if (!is.null(stats)){
    stat=paste0(" -stats ",paste0(normalizePath(stats),collapse=" -stats "))
  }
  
  out_file=paste0(out_file_dir,"/",id,".merged.vcf.stats")

  exec_code=paste0("singularity exec -H ",getwd(),":/home ",sif_gatk,
  " /gatk/gatk  MergeMutectStats -O ", out_file ,stat)

  if(clean){
    exec_code=paste(exec_code," && rm",paste(normalizePath(stats),collapse=" "))
}

  if(mode=="batch"){
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
    stop("gatk failed to run due to unknown error.
    Check std error for more information.")
  }

  job_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    exec_code=exec_code, 
    task_id=task_id,
    input_args = argg,
    out_file_dir=out_file_dir,
    out_files=list(
      merged_stats=out_file)
  )


  if(wait&&mode=="batch"){
    job_validator(job=job_report$job_id,time=update_time,
    verbose=verbose,threads=threads)
  }

  return(job_report)

}


#' Wrapper around CreateSomaticPanelOfNormals from GATK package 
#'
#' This function creates a somatic panel of normals from normal samples bam files.
#' Recommended minimum number of normal samples is around 40.
#' 
#' For more information about this function: https://gatk.broadinstitute.org/hc/en-us/articles/360037058172-CreateSomaticPanelOfNormals-BETA-
#'
#' @param sif_gatk [REQUIRED] Path to gatk sif file.
#' @param bin_bcftools [REQUIRED] Path to bcftools binary file.
#' @param bin_bgzip [REQUIRED] Path to bgzip binary file.
#' @param bin_tabix [REQUIRED] Path to tabix binary file.
#' @param bin_samtools [REQUIRED] Path to samtools binary file.
#' @param normals [OPTIONAL] Path to normal BAM file.
#' @param ref_genome [REQUIRED] Path to reference genome fasta file.
#' @param germ_resource [REQUIRED]Path to germline resources vcf file.
#' @param regions [OPTIONAL] Regions to analyze. If regions for parallelization are not provided then these will be infered from BAM file.
#' @param output_name [OPTIONAL] Name for the output. If not given the name of one of the samples will be used.
#' @param pon [OPTIONAL] Path to panel of normal.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param mnps [OPTIONAL] Report MNPs in vcf file.
#' @param contamination [OPTIONAL] Produce sample cross-contamination reports. Default TRUE.
#' @param orientation [OPTIONAL] Produce read orientation inforamtion. Default FALSE
#' @param filter [OPTIONAL] Filter Mutect2. Default TRUE.
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




create_pon_gatk=function(
  sif_gatk=build_default_sif_list()$sif_gatk,
  bin_bcftools=build_default_tool_binary_list()$bin_bcftools,
  bin_samtools=build_default_tool_binary_list()$bin_samtools,
  bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
  bin_tabix=build_default_tool_binary_list()$bin_tabix,
  output_name="PoN",
  vcfs="",
  regions=NULL,output_dir=".",
  ref_genome=build_default_reference_list()$HG19$reference$genome,
  germ_resource=build_default_reference_list()$HG19$variant$germ_reference,
  verbose=FALSE,
  batch_config=build_default_preprocess_config(),
  threads=4,ram=4,mode="local",
  executor_id=make_unique_id("createPoNGatk"),
  task_name="createPoNGATK",time="48:0:0",
  update_time=60,wait=FALSE,hold=NULL
){

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir,name="pon_gatk")
  job=build_job(executor_id=executor_id,task_id=task_id)


  out_file=paste0(out_file_dir,output_name,".vcf.gz")
  jobs_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    exec_code="", 
    task_id=task_id,
    input_args = argg,
    out_file_dir=out_file_dir,
    out_files=list(
      pon=out_file
    )
  )


  jobs_report[["steps"]][["createGenomicDbGatk"]]<-create_genomic_db_gatk(
      sif_gatk=sif_gatk,
      vcfs=normalizePath(vcfs),output_dir=out_file_dir,
      ref_genome=build_default_reference_list()$HG19$reference$genome,
      verbose=verbose,
      batch_config=batch_config,
      threads=threads,ram=ram,mode=mode,
      executor_id=task_id,
      time=time,
      hold=hold
  )

  hold=unlist_lvl(jobs_report[["steps"]][["createGenomicDbGatk"]],var="job_id")



  exec_code=paste("singularity exec -H ",paste0(getwd(),":/home "),sif_gatk,
  " /gatk/gatk  CreateSomaticPanelOfNormals  -R ",ref_genome,paste0(" -V gendb://",
  normalizePath(jobs_report[["steps"]][["createGenomicDbGatk"]]$out_file_dir))," -O ",out_file)

  if(mode=="batch"){
       out_file_dir2=set_dir(dir=out_file_dir,name="batch")
       batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
       threads=threads,output_dir=out_file_dir2)
       exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
  }


   if(verbose){
       print_verbose(job=job,arg=argg,exec_code=exec_code)
  }

  jobs_report$exec_code=exec_code

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





#' Wrapper around CreateDBImport from GATK package 
#'
#' This function creates a genomic database
#' 
#' For more information about this function: https://gatk.broadinstitute.org/hc/en-us/articles/360037058172-CreateSomaticPanelOfNormals-BETA-
#'
#' @param sif_gatk [REQUIRED] Path to gatk sif file.
#' @param vcfs [OPTIONAL] Path to VCFs files. Only required if BAM files are not given.
#' @param ref_genome [REQUIRED] Path to reference genome fasta file.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param threads [OPTIONAL] Number of threads to split the work. Default 4
#' @param ram [OPTIONAL] RAM memory to asing to each thread. Default 4
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param size [OPTIONAL] Number of VCF per thread. Default 10.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Task EXECUTOR ID. Default "recalCovariates"
#' @param task_name Task name. Default "recalCovariates"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export




create_genomic_db_gatk=function(
  sif_gatk=build_default_sif_list()$sif_gatk,
  bin_samtools=build_default_tool_binary_list()$bin_samtools,
  bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
  bin_tabix=build_default_tool_binary_list()$bin_tabix,
  vcfs="",regions=NULL,output_dir=".",
  ref_genome=build_default_reference_list()$HG19$reference$genome,
  verbose=FALSE,size=10,
  batch_config=build_default_preprocess_config(),
  threads=4,ram=4,mode="local",
  executor_id=make_unique_id("createGenomicDB"),
  task_name="createGenomicDB",time="48:0:0",
  update_time=60,wait=FALSE,hold=NULL
){

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  tmp_dir=set_dir(dir=output_dir,name="tmp")
  job=build_job(executor_id=executor_id,task_id=task_id)


  map_file=paste0(tmp_dir,"/vcfs.map")
  map=data.frame(id=Vectorize(get_file_name)(vcfs),file=vcfs)
  write.table(x=map,file=map_file,quote=FALSE,sep="\t",
  row.names=FALSE,col.names=FALSE)




  jobs_report=build_job_report(
      job_id=job,
      executor_id=executor_id,
      exec_code=list(), 
      task_id=task_id,
      input_args = argg,
      out_file_dir=paste0(output_dir,"/db"),
      out_files=list(
        map_file=map_file
      )
  )


  if(is.null(regions)){

    jobs_report[["steps"]][["getChr"]] <- get_fai_reference_chr(
      fasta=ref_genome,verbose=verbose,output_dir=tmp_dir,
      output_name="refChr",header=FALSE,
      executor_id=task_id,mode="local",threads=threads,ram=ram,
      time=time,update_time=update_time,wait=FALSE,hold=hold)

    regions=jobs_report[["steps"]][["getChr"]]$out_files$ref
  
  }

  

  exec_code=paste("singularity exec -H ",paste0(getwd(),":/home "),sif_gatk,
  " /gatk/gatk GenomicsDBImport -R ",normalizePath(ref_genome),
  " --genomicsdb-workspace-path ",paste0(output_dir,"/db"),
  " --batch-size ",size,
  " --max-num-intervals-to-import-in-parallel ",threads,
  " --tmp-dir ",normalizePath(tmp_dir)," --reader-threads ",threads,
  " --overwrite-existing-genomicsdb-workspace TRUE",
  " --sample-name-map ",normalizePath(map_file), " -L ",regions)

  if(mode=="batch"){
       out_file_dir2=set_dir(dir=output_dir,name="batch")
       batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
       threads=threads,output_dir=out_file_dir2)
       exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
  }


   if(verbose){
       print_verbose(job=job,arg=argg,exec_code=exec_code)
  }

   jobs_report$exec_code=exec_code



  error=execute_job(exec_code=exec_code)
  
  if(error!=0){
    stop("gatk failed to run due to unknown error.
    Check std error for more information.")
  }

  if(wait&&mode=="batch"){
    job_validator(job=jobs_report$job_id,
    time=update_time,verbose=verbose,threads=threads)
  }

  return(jobs_report)

}




#' Variant Calling using HaplotypeCaller
#'
#' This function functions calls HaplotypeCaller for variant calling.
#' If a vector of normal samples are provided these will be processed in multi-sample mode.
#' To run in normal mode suppply a single normal sample.
#' TO DO// Implement mitochondrial mode feature
#' 
#' For more information read:
#' https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller
#'
#' @param sif_gatk [REQUIRED] Path to gatk sif file.
#' @param normal [OPTIONAL] Path to normal BAM file.
#' @param ref_genome [REQUIRED] Path to reference genome fasta file.
#' @param rdata [OPTIONAL] Import R data information with list of BAM.
#' @param selected [OPTIONAL] Select BAM from list.
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


haplotypecaller_gatk=function(
    rdata=NULL,selected=NULL,
    sif_gatk=build_default_sif_list()$sif_gatk,
    normal="",
    region="",
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    output_dir=".",
    tmp_dir=".",
    output_name="",
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=4,ram=4,mode="local",
    executor_id=make_unique_id("parallelCallHaplotypecallerGatk"),
    task_name="parallelCallHaplotypecallerGatk",time="48:0:0",
    update_time=60,wait=FALSE,hold=NULL
){

  if(!is.null(rdata)){
    load(rdata)
    if(!is.null(selected)){
      region=region_list[selected]
    }
  }

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir,name="vcf")
  job=build_job(executor_id=executor_id,task_id=task_id)



  if(tmp_dir!=""){
    tmp_dir=paste0(" --tmp-dir ",normalizePath(tmp_dir))
  }
  
  id=""
  if(output_name!=""){
    id=output_name
  }else{  
    id=get_file_name(normal[1])
  }
  
  reg=""
  if (region==""){
      out_file=paste0(out_file_dir,"/",id,".unfilt.vcf.gz")
  }else{
      reg=paste0(" -L ",region)
      id=paste0(id,".",region)
      out_file=paste0(out_file_dir,"/",id,".unfilt.vcf.gz")
  }

  if (is.vector(normal)){
    normal=paste0(" -I ",paste(normalizePath(normal),collapse=" -I "))
  }else{
    normla=paste0(" -I ",normalizePath(normal))
  }


  exec_code=paste("singularity exec -H ",paste0(getwd(),":/home "),sif_gatk,
  " /gatk/gatk HaplotypeCaller -R ",normalizePath(ref_genome), normal," -O ",out_file,reg)

  if(mode=="batch"){
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
    stop("gatk failed to run due to unknown error.
    Check std error for more information.")
  }

  job_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    exec_code=exec_code, 
    task_id=task_id,
    input_args = argg,
    out_file_dir=out_file_dir,
    out_files=list(
      vcf=out_file,
      idx=paste0(out_file,".idx")
    )
  )


  if(wait&&mode=="batch"){
    job_validator(job=job_report$job_id,time=update_time,
    verbose=verbose,threads=threads)
  }

  return(job_report)
}



#' Variant Calling using HaplotypeCaller
#'
#' This function functions calls HaplotypeCaller for variant calling.
#' If a vector of normal samples are provided these will be processed in multi-sample mode.
#' To run in normal mode suppply a single normal sample.
#' TO DO// Implement mitochondrial mode feature
#' 
#' For more information read:
#' https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller
#'
#' @param sif_gatk [REQUIRED] Path to gatk sif file.
#' @param bam [OPTIONAL] Path to normal BAM file.
#' @param ref_genome [REQUIRED] Path to reference genome fasta file.
#' @param region [REQUIRED] Genomic position in samtools format chr:start-end.
#' @export


new_haplotypecaller_gatk=function(
    sif_gatk=build_default_sif_list()$sif_gatk,
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    indel_db=build_default_reference_list()$HG19$variant$mills_reference,
    haplotype_db=build_default_reference_list()$HG19$variant$hapmap_reference,
    bam=NULL,
    region=NULL,
    score=TRUE,
    score_CNN="CNN_1D",
    filter=TRUE,
    snp_tranche=99.95,
    indel_tranche=99.4,
    keep_previous_filters=FALSE,
    ...
){

   run_main=function(
    .env
  ){

    .this.env=environment()
    append_env(to=.this.env,from=.env)
    set_main(.env=.this.env)


  
    tmp_dir=paste0(" --tmp-dir ",normalizePath(tmp_dir))
    if(mode=="local_parallel"){
      threads=1
    }
  
    if (is.null(region)){
        .main$out_files$unfiltered_vcf=paste0(
          out_file_dir,"/",input_id,".unfilt.vcf.gz"
        )
    }else{
        .main$out_files$unfiltered_vcf=paste0(
          out_file_dir,"/",paste0(input_id,".",input),".unfilt.vcf.gz"
        )
        .main$out_files$unfiltered_vcf_index=paste0(
          out_file_dir,"/",paste0(input_id,".",input),".unfilt.vcf.gz.tbi"
        )
        region=paste0(" -L ",input)
    }

    .main$exec_code=paste(
      "singularity exec -H ",paste0(getwd(),":/home "),sif_gatk,
      " /gatk/gatk HaplotypeCaller --native-pair-hmm-threads ",threads," -R ",
      normalizePath(ref_genome), 
      paste0(" -I ",normalizePath(bam)),
      " -O ",.main$out_files$unfiltered_vcf,
      region
    )

    run_job(.env=.this.env)

    .main.step=.main$steps[[fn_id]]
    
    if(score){
        .main.step$steps<-append(
        .main.step$steps,
        cnn_score_variants_gatk(
          sif_gatk=sif_gatk,
          ref_genome=normalizePath(ref_genome),
          haplotype_db=normalizePath(haplotype_db),
          indel_db=normalizePath(indel_db),
          bam=bam,
          vcf=.main.step$out_files$unfiltered_vcf,
          filter=filter,
          info_key=score_CNN,
          snp_tranche=snp_tranche,
          indel_tranche=indel_tranche,
          keep_previous_filters=keep_previous_filters,
          output_name=paste0(input_id,".",input),
          output_dir=paste0(out_file_dir,"/scored"),
          tmp_dir=tmp_dir,
          batch_dir=batch_dir,
          env_dir=env_dir,
          err_msg=err_msg,
          verbose=verbose,
          threads=threads,
          ram=ram,
          executor_id=task_id
        )
      )
      .this.step=.main.step$steps$cnn_score_variants_gatk
      .main.step$out_files=append(.main.step$out_files,.this.step$out_files)
    }
    .env$.main<-.main
   }

  
    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      output_name=get_file_name(bam),
      vars="region"
    )

    launch(.env=.base.env)
}




#' Variant Calling using HaplotypeCaller
#'
#' This function functions calls HaplotypeCaller for variant calling.
#' If a vector of normal samples are provided these will be processed in multi-sample mode.
#' To run in normal mode suppply a single normal sample.
#' TO DO// Implement mitochondrial mode feature
#' 
#' For more information read:
#' https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller
#'
#' @param sif_gatk [REQUIRED] Path to gatk sif file.
#' @param bam [OPTIONAL] Path to normal BAM file.
#' @param ref_genome [REQUIRED] Path to reference genome fasta file.
#' @param region [REQUIRED] Genomic position in samtools format chr:start-end.
#' @export


call_haplotypecaller_gatk=function(
    sif_gatk=build_default_sif_list()$sif_gatk,
    bin_samtools=build_default_tool_binary_list()$bin_samtools,
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    indel_db=build_default_reference_list()$HG19$variant$mills_reference,
    haplotype_db=build_default_reference_list()$HG19$variant$hapmap_reference,
    bam=NULL,
    patient_id=NULL,
    chromosomes=c(1:22,"X","Y"),
    region=NULL,
    score=TRUE,
    filter=TRUE,
    score_CNN="CNN_1D",
    snp_tranche=99.95,
    indel_tranche=99.4,
    keep_previous_filters=FALSE,
    ...
){

   run_main=function(
    .env
  ){

    .this.env=environment()
    append_env(to=.this.env,from=.env)
    out_file_dir=set_dir(
        out_file_dir,
        name=paste0(patient_id,"/haplotypecaller_reports/",input_id)
    )

    set_main(.env=.this.env)


    .main$steps[[fn_id]]<-.this.env
    .main.step=.main$steps[[fn_id]]

    ### IF NO REGION IS GIVEN SPLIT BAM PER CHROMOSOME
    if(is.null(region)){
      .main.step$steps <-append(
        .main.step$steps,
        get_sq_bam(
          bin_samtools=bin_samtools,
          bam=input,
          output_name=input_id,
          output_dir=tmp_dir,
          header=TRUE,
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
      .main.step$out_files$region=.this.step$out_files$index_bed
      region=.main.step$out_files$region
     
    }

    ### ASCERTAIN REGION INPUT IF GIVEN
    
   if(length(region)==1){
      ### IF PATH READ AS BED
      if (file.exists(region)){
        region=read_bed(
          bed=region
        )$body
      }

      ### IF DATA.FRAME GENERATE GID
      if (is.data.frame(region)){
        region=region[region$chrom %in% chromosomes,]
      ### CHANGE TO BASE 1 SYSTEM
        region$gid=paste0(region$chrom,":",as.numeric(region$chromStart)+1,"-",as.numeric(region$chromEnd))
        region=unlist(region$gid)
      }
   }

    ###  ONLY RUN LOCALLY PARALLEL JOBS AT REGION LEVEL IF SAMPLE LEVEL LOCAL/BATCH MODE
    ###  LIMITATION OF MCLAPPLY

    run_mode="local_parallel"
    if(mode=="local_parallel"){
      run_mode="local"
    }
  
    .main.step$steps <-append(
      .main.step$steps,
      new_haplotypecaller_gatk(
        sif_gatk=sif_gatk,
        ref_genome=normalizePath(ref_genome),
        haplotype_db=normalizePath(haplotype_db),
        indel_db=normalizePath(indel_db),
        bam=input,
        output_name=input_id,
        region=region,
        score=score,
        filter=filter,
        score_CNN=score_CNN,
        snp_tranche=snp_tranche,
        indel_tranche=indel_tranche,
        keep_previous_filters=keep_previous_filters,
        mode=run_mode,
        output_dir=paste0(out_file_dir,"/unfiltered"),
        tmp_dir=tmp_dir,
        env_dir=env_dir,
        batch_dir=batch_dir,
        err_msg=err_msg,
        verbose=verbose,
        threads=threads,
        ram=ram,
        executor_id=task_id
    ))

    .this.step=.main.step$steps
    .main.step$out_files$region_vcfs=get_variable_env(env=.this.step)
    .env$.main<-.main

  }

  

  .base.env=environment()
  list2env(list(...),envir=.base.env)
  set_env_vars(
    .env= .base.env,
    vars="bam"
  )
  launch(.env=.base.env)

}



#' Germline Variant Scoring Using CNN
#'
#' This function functions calls CNNScoreVariants.
#' If a vector of normal samples are provided these will be processed in multi-sample mode.
#' To run in normal mode suppply a single normal sample.
#' TO DO// Implement mitochondrial mode feature
#' 
#' For more information read:
#' https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller
#'
#' @param sif_gatk [REQUIRED] Path to gatk sif file.
#' @param vcf [REQUIRE] Path to VCF file.
#' @param ref_genome [REQUIRED] Path to reference genome fasta file.
#' @param bam [OPTIONAL] Path to BAM file. If given then 2D CNN will be applied
#' @export
#' 
cnn_score_variants_gatk=function(
  sif_gatk=build_default_sif_list()$sif_gatk,
  ref_genome=build_default_reference_list()$HG19$reference$genome,
  indel_db=build_default_reference_list()$HG19$variant$mills_reference,
  haplotype_db=build_default_reference_list()$HG19$variant$hapmap_reference,
  vcf=NULL,
  bam=NULL,
  filter=TRUE,
  info_key="CNN_1D",
  snp_tranche=99.95,
  indel_tranche=99.4,
  keep_previous_filters=FALSE,
  ...
){

   run_main=function(
    .env
  ){

    .this.env=environment()
    append_env(to=.this.env,from=.env)
    set_main(.env=.this.env)

    opt=""
    .main$out_files$scored_vcf=paste0(out_file_dir,"/",input_id,".CNNscored.1D.vcf.gz")
    .main$out_files$scored_vcf_index=paste0(out_file_dir,"/",input_id,".CNNscored.1D.vcf.gz.tbi")
    if (info_key=="CNN_2D"){
       bam=paste0(" -I ",bam)
       opt=" -tensor-type read_tensor "
      .main$out_files$scored_vcf=paste0(out_file_dir,"/",input_id,".CNNscored.2D.vcf.gz")
      .main$out_files$scored_vcf_index=paste0(out_file_dir,"/",input_id,".CNNscored.2D.vcf.gz.tbi")
    }else{
      bam=NULL
    }

    .main$exec_code=paste(
      "singularity exec -H ",paste0(getwd(),":/home "),sif_gatk,
      " /gatk/gatk CNNScoreVariants -R ",
      normalizePath(ref_genome), 
      paste0(" -V ",normalizePath(vcf)),
      " -O ",.main$out_files$scored_vcf,
      bam,opt
    )

     run_job(.env=.this.env)
    
    .main.step=.main$steps[[fn_id]]

    if(filter){
      .main.step$steps <-append(
      .main.step$steps,
      filter_variant_tranches_gatk(
        sif_gatk=build_default_sif_list()$sif_gatk,
        ref_genome=normalizePath(ref_genome),
        indel_db=normalizePath(indel_db),
        haplotype_db=normalizePath(haplotype_db),
        vcf=.main$out_files$scored_vcf,
        info_key=info_key,
        snp_tranche=snp_tranche,
        indel_tranche=indel_tranche,
        keep_previous_filters=keep_previous_filters,
        output_dir=paste0(out_file_dir,"/filtered"),
        output_name=input_id,
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
      .this.step=.main.step$steps$filter_variant_tranches_gatk
      .main.step$out_files=append(.main.step$out_files,.this.step$out_files)
    }
    .env$.main<-.main
  } 
    
   .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="vcf"
    )

    launch(.env=.base.env)

}








#' Variant Trench Filtering using GATK
#'
#' This function filters variant tranches calls generated using CNNScoreVariants method
#' 
#' For more information read:
#' https://gatk.broadinstitute.org/hc/en-us/articles/360051308071-FilterVariantTranches
#'
#' @param sif_gatk [REQUIRED] Path to gatk sif file.
#' @param vcf [REQUIRED] Path to VCF file.
#' @param ref_genome [REQUIRED] Path to reference genome fasta file.
#' @param indel_db [OPTIONAL] Path to database with indel info.
#' @param haplotype_db [OPTIONAL] Path to database with haplotype info.
#' @param info_key [OPTIONAL] Info Key for CNN model training. Default CNN_1D. Options ["CNN_1D","CNN_2D"]
#' @param snp_tranche [OPTIONAL] SNP tranche cut-off.
#' @param indel_tranche [OPTIONAL] INDEL tranche cut-off.
#' @param keep_previous_tranche [OPTIONAL] Remove previous filter information.
#' @export



filter_variant_tranches_gatk=function(
  sif_gatk=build_default_sif_list()$sif_gatk,
  ref_genome=build_default_reference_list()$HG19$reference$genome,
  indel_db=build_default_reference_list()$HG19$variant$mills_reference,
  haplotype_db=build_default_reference_list()$HG19$variant$hapmap_reference,
  vcf=NULL,
  info_key="CNN_1D",
  snp_tranche=99.95,
  indel_tranche=99.4,
  keep_previous_filters=FALSE,
  ...
){



  run_main=function(
    .env
  ){
    .this.env=environment()
    append_env(to=.this.env,from=.env)
    set_main(.env=.this.env)


    prev_filters=" --invalidate-previous-filters "
    if (keep_previous_filters){
      prev_filters=" "
    }

    if(!is.null(indel_db)){
      indel_db=paste0(" --resource ",normalizePath(indel_db))
    }

    if(!is.null(haplotype_db)){
      haplotype_db=paste0(" --resource ",normalizePath(haplotype_db))
    }

    .main$out_files$filtered_vcf=paste0(out_file_dir,"/",input_id,".filteredTranches.vcf.gz")
    .main$out_files$filtered_vcf_index=paste0(out_file_dir,"/",input_id,".filteredTranches.vcf.gz.tbi")
    .main$exec_code=paste(
      "singularity exec -H ",paste0(getwd(),":/home "),sif_gatk,
      " /gatk/gatk  FilterVariantTranches -R ",normalizePath(ref_genome), 
      " -O ", .main$out_files$filtered_vcf, 
      " -V ",normalizePath(input), " --info-key ",info_key,
      " --snp-tranche ",snp_tranche,
      " --indel-tranche ",indel_tranche,
      indel_db,
      haplotype_db,prev_filters
    )

    run_job(.env=.this.env)
    .env$.main<-.main
  }

  .base.env=environment()
  list2env(list(...),envir=.base.env)
  set_env_vars(
      .env=.base.env,
      vars="vcf"
  )

  launch(.env=.base.env)

}






#' Convert FASTQ files to unmapped SAM/BAM format using GATK
#'
#' This function converts paired-end FASTQ files to unmapped SAM/BAM format.
#' The output is an unmapped BAM file that can be used as input for alignment and variant calling pipelines.
#' Sample name is automatically extracted from the intersection of R1 and R2 file names.
#' 
#' For more information read:
#' https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-FastqToSam
#'
#' @param sif_gatk [REQUIRED] Path to gatk sif file.
#' @param fastq [REQUIRED] Path to fastq files in list format with R1 and R2 paired reads.
#' @export
#' 
fastq_to_sam_gatk=function(
  sif_gatk=build_default_sif_list()$sif_gatk,
  fastq=NULL,
  ...
){

   run_main=function(
    .env
  ){

    .this.env=environment()
    append_env(to=.this.env,from=.env)
    set_main(.env=.this.env)

    .main$out_files$bam=paste0(out_file_dir,"/",input_id,".bam")


    .main$exec_code=paste(
      "singularity exec ",sif_gatk,
      " /gatk/gatk  FastqToSam",
      " -F1", input$fastq_r1,
      " -F2 ",input$fastq_r2,
      " -O ",.main$out_files$bam,
      " -SM ",input_id
    )

     run_job(.env=.this.env)

    .env$.main<-.main
  } 
    
   .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="fastq"
    )

    launch(.env=.base.env)

}





#' Convert SAM/BAM files to paired-end FASTQ format using GATK
#'
#' This function converts SAM/BAM files to paired-end FASTQ format.
#' The output are paired-end FASTQ files (R1 and R2) that can be used for downstream analysis.
#' Clipping attributes can be applied to filter reads based on quality scores.
#' 
#' For more information read:
#' https://gatk.broadinstitute.org/hc/en-us/articles/360036883491-SamToFastq
#'
#' @param sif_gatk [REQUIRED] Path to gatk sif file.
#' @param bam [REQUIRED] Path to input BAM/SAM file to convert to FASTQ format.
#' @param clipping_attribute [OPTIONAL] Apply clipping attributes during conversion. Default TRUE.
#' @export
#' 
sam_to_fastq_gatk=function(
  sif_gatk=build_default_sif_list()$sif_gatk,
  bam=NULL,
  clipping_attribute=TRUE,
  ...
){

   run_main=function(
    .env
  ){

    .this.env=environment()
    append_env(to=.this.env,from=.env)
    set_main(.env=.this.env)

    .main$out_files$fastq_r1=paste0(out_file_dir,"/",input_id,"_R1.fastq")
    .main$out_files$fastq_r2=paste0(out_file_dir,"/",input_id,"_R2.fastq")

    .main$exec_code=paste(
      "singularity exec ",sif_gatk,
      " /gatk/gatk  SamToFastq -F",.main$out_files$fastq_r1,
      " -F2 ",.main$out_files$fastq_r2,
      " -I ", input,
      ifelse(clipping_attribute," --CLIPPING_ATTRIBUTE XT --CLIPPING_ACTION 2 "," ")

    )

     run_job(.env=.this.env)

    .env$.main<-.main
  } 
    
   .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="bam"
    )

    launch(.env=.base.env)
}





#' Merge Aligned and Unmapped BAM Files using GATK
#'
#' Merge Aligned and Unmapped BAM Files using GATK
#'
#' This function merges aligned and unmapped BAM files using GATK's MergeBamAlignment tool.
#' It is typically used after UMI extraction and alignment workflows to combine the aligned
#' reads with the original unmapped reads, preserving UMI and other metadata.
#'
#' The function provides configurable options for:
#' \itemize{
#'   \item Attribute retention/removal (default: retain X0, remove NM and MD)
#'   \item BAM sort order (default: queryname)
#'   \item Including/excluding aligned reads only
#'   \item Adding mate CIGAR information for paired reads
#'   \item Primary alignment strategy: MostDistant
#'   \item Proper pair flags alignment enabled
#'   \item Overlapping reads not clipped
#'   \item No limit on insertions/deletions (-1)
#' }
#'
#' For more information read:
#' https://gatk.broadinstitute.org/hc/en-us/articles/360036452392-MergeBamAlignment
#'
#' @param sif_gatk [REQUIRED] Path to GATK Singularity container image. 
#'   Default: from build_default_sif_list()$sif_gatk.
#' @param ref_genome [REQUIRED] Path to reference genome FASTA file. 
#'   Default: HG19 reference from build_default_reference_list().
#' @param bam [REQUIRED] BAM file(s) to merge. Should contain mapped and unmapped BAM information 
#'   accessible via input$mapped and input$unmapped.
#' @param output_name [OPTIONAL] Name for the output files. Default: "sample".
#' @param attributes [OPTIONAL] Character vector of BAM tags to retain. 
#'   Default: c("XO","NM","MD"). Tags listed here will be retained from input BAM files.
#' @param sort_order [OPTIONAL] Output sort order for merged BAM. 
#'   Options: "queryname", "coordinate". Default: "queryname".
#' @param aligned_reads_only [OPTIONAL] Include only aligned reads in output. Default: TRUE.
#' @param add_mate_cigar [OPTIONAL] Add CIGAR information for mate in MC tag. Default: TRUE.
#' @param ... Additional arguments passed to internal functions for environment setup and job execution.
#'
#' @return Merged BAM file with original reads and alignments combined. Output path tracked in job report.
#'   Output BAM is accessible via .main$out_files$merged_bam.
#'
#' @seealso
#'   \link{group_by_umi_fgbio} for grouping reads by UMI before consensus calling
#'
#' @export
#' 
merge_bam_umi_gatk=function(
  sif_gatk=build_default_sif_list()$sif_gatk,
  ref_genome=build_default_reference_list()$HG19$reference$genome,
  bam=NULL,
  output_name="sample",
  attributes=c("X0","NM","MD"),
  sort_order="queryname",
  aligned_reads_only=TRUE,
  add_mate_cigar=FALSE,
  ...
){

   run_main=function(
    .env
  ){

    .this.env=environment()
    append_env(to=.this.env,from=.env)
    set_main(.env=.this.env)

    .main$out_files$bam=paste0(out_file_dir,"/",input_id,".merged.bam")

    .main$exec_code=paste(
      "singularity exec ",sif_gatk,
      " /gatk/gatk MergeBamAlignment",
      paste(" --ATTRIBUTES_TO_RETAIN ",attributes,collapse=" "),
      " --ALIGNED_BAM ", input$mapped,
      " --UNMAPPED_BAM ", input$unmapped, 
      " --OUTPUT ", .main$out_files$bam,
      " --REFERENCE_SEQUENCE ", ref_genome,
      " --SORT_ORDER ",sort_order,
      ifelse(aligned_reads_only," --ALIGNED_READS_ONLY true  ",""),
      ifelse(add_mate_cigar," --ADD_MATE_CIGAR true ",""),
      " --MAX_INSERTIONS_OR_DELETIONS -1",
      " --PRIMARY_ALIGNMENT_STRATEGY MostDistant ",
      " --ALIGNER_PROPER_PAIR_FLAGS true ",
      " --CLIP_OVERLAPPING_READS false "
      )

     run_job(.env=.this.env)

    .env$.main<-.main
  } 
    
   .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="bam"
    )

    launch(.env=.base.env)
}








#' GATK Base Quality Score Recalibration (BQSR) workflow
#'
#' High-level wrapper that runs GATK's BaseRecalibrator and ApplyBQSR steps
#' following GATK best-practices. The function constructs containerized
#' commands (via Singularity) and submits them through the package's job
#' runner. It integrates with the pipeline environment and expects common
#' runner-managed variables (for example `out_file_dir`, `input`,
#' `input_id`, and `fn_id`) when called via `launch()`.
#'
#' @param bin_samtools Path to the `samtools` binary. Defaults to
#'   `build_default_tool_binary_list()$bin_samtools`.
#' @param sif_gatk Path to the Singularity image containing GATK. Defaults
#'   to `build_default_sif_list()$sif_gatk`.
#' @param bin_picard Path to the Picard jar. Defaults to
#'   `build_default_tool_binary_list()$bin_picard`.
#' @param ref_genome Path to the reference genome FASTA (required for
#'   recalibration and realignment steps).
#' @param dbsnp Path to known-sites VCF(s) used by `BaseRecalibrator`.
#' @param bam Path to the input BAM file to recalibrate.
#' @param ... Additional named arguments forwarded into the internal job
#'   environment (commonly: `output_dir`, `tmp_dir`, `threads`, `ram`,
#'   `mode`, `batch_config`, `verbose`, `clean`, etc.).
#'
#' @details
#' The function implements a two-step BQSR process: it first calls
#' `BaseRecalibrator` to compute recalibration tables and then calls
#' `ApplyBQSR` to produce the recalibrated BAM. When invoked with a
#' `region` argument (via the runner) the function will restrict
#' recalibration to that interval and name outputs accordingly.
#'
#' @return Invisibly returns the job report produced by the internal job
#'   runner and writes recalibration tables and recalibrated BAM(s) into
#'   the configured output directory.
#' @export


new_recal_gatk=function(
  bin_samtools=build_default_tool_binary_list()$bin_samtools,
  sif_gatk=build_default_sif_list()$sif_gatk,
  bin_picard=build_default_tool_binary_list()$bin_picard,
  ref_genome=build_default_reference_list()$HG19$reference$genome,
  dbsnp=build_default_reference_list()$HG19$database$all_common,
  chromosomes=c(1:22,"X","Y"),
  bam=NULL,
  ...
  ){

      run_main=function(
            .env
      ){

              # Copy environment variables from parent scope to current environment
              .this.env=environment()
              append_env(to=.this.env,from=.env)
              set_main(.env=.this.env)

              # Store this function's environment in main steps registry
              .main$steps[[fn_id]]<-.this.env
              # Reference the current processing step for appending results
              .main.step=.main$steps[[fn_id]]

              # Record pipeline start time for elapsed time tracking
              start_time <- Sys.time()
                  
        
        # Define BQSR processing pipeline steps in logical execution order.
        # Each step name directly corresponds to conditional processing blocks below.
        # Steps are designed to handle: reference extraction → pre-recal tables → post-recal tables → recalibrated BAM → gather/sort → covariate analysis
        
        
        steps=c(
            "get_chrom",        # Step 1: Extract chromosome sizes and create genome BED from BAM header
            "before_bqsr",      # Step 2: Run BaseRecalibrator to compute pre-recalibration tables
            "before_merge_bqsr",
            "apply_bqsr",       # Step 3: Run ApplyBQSR to recalibrate BAM using pre-recal tables
            "gather_bam",       # Step 4: Gather scattered recalibrated BAMs into single file
            "sort_bam",         # Step 5: Sort and index the final recalibrated BAM
            "after_bqsr", 
            "after_merge_bqsr",      # Step 6: Run BaseRecalibrator again on recalibrated BAM for comparison
            "analyze_covariates"  # Step 7: Generate covariate analysis plots comparing before/after BQSR
        )
        
        # Total number of BQSR pipeline steps for progress reporting and loop control
        total_steps=length(steps)

        # Execute each BQSR step sequentially with progress logging
        for(step in 1:total_steps){
            
            # Display pipeline progress: current step number, total steps, and step name
            logger(paste("Running step", step, "of", total_steps, ":", steps[step]),start_time)


            tryCatch({

              # --- STEP 1: Extract chromosome sizes and genome structure ---

              if(steps[step]=="get_chrom"){
                  .main.step$steps <-append(
                            .main.step$steps,get_ref_from_bam(
                            bin_samtools=bin_samtools,
                            bam=input,
                            output_dir=tmp_dir,
                            output_name=input_id,
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
                  .this.step=.main.step$steps$get_ref_from_bam
                  .main.step$out_files$genome_bed=.this.step$out_files$genome_bed
                    
                  regions=read.table(.main.step$out_files$genome_bed,
                  sep="\t",header=TRUE) %>% dplyr::filter(chr %in% chromosomes) %>% 
                  dplyr::mutate(regions=paste0(chr,":",start+1,"-",end-1))

              }
                

              # --- STEP 2: Compute base recalibration tables (BaseRecalibrator) ---
              if(steps[step]=="before_bqsr"){

                  .main.step$steps$new_generate_BQSR_gatk.before<-
                              new_generate_BQSR_gatk(
                              sif_gatk=sif_gatk,
                              ref_genome=ref_genome,
                              dbsnp=dbsnp,
                              bam=input,
                              region=regions$regions,
                              output_dir=paste0(out_file_dir,"/before_recal/tables"),
                              output_name=paste0(input_id),
                              tmp_dir=tmp_dir,
                              env_dir=env_dir,
                              batch_dir=batch_dir,
                              err_msg=err_msg,
                              verbose=verbose,
                              threads=threads,
                              mode="local_parallel",
                              ram=ram,
                              fn_id="before",
                              executor_id=task_id
                      )

                  .this.step=.main.step$steps$new_generate_BQSR_gatk.before
                  .main.step$out_files$before_bqsr$table$scattered=get_variable_env(env=.this.step)
              }

             # --- STEP 2: Compute base recalibration tables (BaseRecalibrator) ---
              if(steps[step]=="before_merge_bqsr"){
                  .main.step$steps<-append(
                    .main.step$steps,
                              new_gather_BQSR_reports_gatk(
                                sif_gatk=sif_gatk,
                                report=.main.step$out_files$before_bqsr$table$scattered,
                                clean_reports=TRUE,
                                output_dir=paste0(out_file_dir,"/before_recal/tables"),
                                output_name=paste0(input_id),
                                tmp_dir=tmp_dir,
                                env_dir=env_dir,
                                batch_dir=batch_dir,
                                err_msg=err_msg,
                                verbose=verbose,
                                threads=threads,
                                ram=ram,
                                fn_id="before",
                                executor_id=task_id
                      )
                  )

                  .this.step=.main.step$steps$new_gather_BQSR_reports_gatk.before
                  .main.step$out_files$before_bqsr$table$merged=.this.step$out_files
                
              }



              # --- STEP 3: Apply recalibration to BAM (ApplyBQSR) ---
              if(steps[step]=="apply_bqsr"){
                .main.step$steps$new_apply_BQSR_gatk<- 
                            new_apply_BQSR_gatk(
                            sif_gatk=sif_gatk,
                            ref_genome=ref_genome,
                            dbsnp=dbsnp,
                            bam=input,
                            region=regions$regions,
                            rec_table=.main.step$out_files$before_bqsr$table$scatter,
                            output_dir=paste0(out_file_dir,"/before_recal/bam"),
                            output_name=paste0(input_id),
                            tmp_dir=tmp_dir,
                            env_dir=env_dir,
                            batch_dir=batch_dir,
                            err_msg=err_msg,
                            verbose=verbose,
                            threads=threads,
                            mode="local_parallel",
                            ram=ram,
                            executor_id=task_id
                    )

                .this.step=.main.step$steps$new_apply_BQSR_gatk
                .main.step$out_files$before_bqsr$bam=get_variable_env(env=.this.step)
              }

              # --- STEP 4: Gather scattered BAMs into single file (Picard) ---
                if(steps[step]=="gather_bam"){
                  
                  .main.step$steps <-append(
                            .main.step$steps,
                              new_gather_bam_files_picard(
                              bin_picard = bin_picard,
                              bam=.main.step$out_files$before_bqsr$bam,
                              region=regions$regions,
                              rec_table=.main.step$out_files$before_bqsr$table,
                              output_dir=paste0(out_file_dir,"/before_recal/bam"),
                              output_name=paste0(input_id,".recal"),
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

                  .this.step=.main.step$steps$new_gather_bam_files_picard
                  .main.step$out_files$before_bqsr$bam$unsorted=.this.step$out_files$bam
                
              }


              ### STEP 5

              if(steps[step]=="sort_bam"){

                  .main.step$steps <-append(
                      .main.step$steps,
                          new_sort_and_index_bam_samtools(
                                  bin_samtools=bin_samtools,
                                  bam=.main.step$out_files$before_bqsr$bam$unsorted,
                                  sort=TRUE,
                                  index=TRUE,
                                  coord_sort=TRUE,
                                  stats=TRUE,
                                  output_dir=paste0(out_file_dir,"/before_recal/bam/sorted"),
                                  output_name=paste0(input_id,".recal.sorted"),
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

                  .this.step=.main.step$steps$new_sort_and_index_bam_samtools
                  .main.step$out_files$before_bqsr$bam$sorted=.this.step$out_files
              }




              # --- STEP 6: Re-analyze recalibration on recalibrated BAM (BaseRecalibrator) ---

              if(steps[step]=="after_bqsr"){
                  .main.step$steps$new_generate_BQSR_gatk.after <-
                            new_generate_BQSR_gatk(
                            sif_gatk=sif_gatk,
                            ref_genome=ref_genome,
                            dbsnp=dbsnp,
                            bam=.main.step$out_files$before_bqsr$bam$sorted$srt_bam,
                            region=regions$regions,
                            output_dir=paste0(out_file_dir,"/after_recal/tables"),
                            output_name=paste0(input_id),
                            tmp_dir=tmp_dir,
                            env_dir=env_dir,
                            batch_dir=batch_dir,
                            err_msg=err_msg,
                            verbose=verbose,
                            threads=threads,
                            mode="local_parallel",
                            ram=ram,
                            fn_id="after",
                            executor_id=task_id
                    )
              

                .this.step=.main.step$steps$new_generate_BQSR_gatk.after
                .main.step$out_files$after_bqsr$table$scattered=get_variable_env(env=.this.step)
              }



              # --- STEP 2: Compute base recalibration tables (BaseRecalibrator) ---
              if(steps[step]=="after_merge_bqsr"){
                  .main.step$steps<-append(
                    .main.step$steps,
                              new_gather_BQSR_reports_gatk(
                                sif_gatk=sif_gatk,
                                report=.main.step$out_files$after_bqsr$table$scattered,
                                clean_reports=TRUE,
                                output_dir=paste0(out_file_dir,"/after_recal/tables"),
                                output_name=paste0(input_id),
                                tmp_dir=tmp_dir,
                                env_dir=env_dir,
                                batch_dir=batch_dir,
                                err_msg=err_msg,
                                verbose=verbose,
                                threads=threads,
                                ram=ram,
                                fn_id="after",
                                executor_id=task_id
                      )
                  )

                  .this.step=.main.step$steps$new_gather_BQSR_reports_gatk.after
                  .main.step$out_files$after_bqsr$table$merged=.this.step$out_files
                
              }


              ###STEP 7

              if(steps[step]=="analyze_covariates"){
                  .main.step$steps <-append(
                          .main.step$steps,
                            new_analyze_covariates_gatk(
                            sif_gatk=sif_gatk,
                            before=.main.step$out_files$before_bqsr$table$merged,
                            after=.main.step$out_files$after_bqsr$table$merged,
                            output_dir=paste0(out_file_dir),
                            output_name=paste0(input_id),
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

                .this.step=.main.step$steps$new_analyze_covariates_gatk
                .main.step$out_files$covariates=.this.step$out_files
              }

            # Log successful step completion
            logger(paste("Completed step", step, "of", total_steps, ":", steps[step]),start_time)
          

          }, error=function(e){
                          # Handle step execution errors with informative message
                          logger(paste("ERROR in step", step, ":", steps[step]),start_time)
                          stop(paste("Step '" , steps[step], "' failed. Error:", e$message,
                                    "\nReview input files and parameters before retrying."))
        
        })

        # Log pipeline completion with total runtime
        total_elapsed <- as.numeric(difftime(Sys.time(), start_time, units="secs"))
        total_elapsed_str <- sprintf("%.1f", total_elapsed)
        logger(paste("GATK Recalibration pipeline completed successfully."),start_time)
        logger(paste("Total steps executed:", total_steps, "| Total runtime:", total_elapsed_str, "seconds"),start_time)
          
        # Return main object to parent environment for job tracking and output reporting
        .env$.main <- .main

      }
  }



       # Setup execution environment by copying parent scope variables
    .base.env=environment()
    # Merge additional parameters passed via ... into execution environment
    list2env(list(...),envir=.base.env)
    # Configure environment variables required for batch processing and temp directories
    set_env_vars(
        .env= .base.env,
        vars="bam"
    )

    # Launch the UMI processing pipeline with fully prepared environment
    launch(.env=.base.env)
  

}




 
 
 
#' Generate GATK BaseRecalibrator (BQSR) table
#'
#' Runs GATK BaseRecalibrator inside a Singularity image to produce a BQSR
#' recalibration table from an input BAM and known-sites VCF(s). Optionally the
#' operation can be restricted to a genomic region. Additional execution
#' parameters (tmp dirs, batch settings, threads, ram, etc.) are passed via
#' the variadic arguments and handled by the internal launch environment.
#'
#' @param sif_gatk Path to the GATK Singularity image. Default build_default_sif_list()$sif_gatk.
#' @param ref_genome Path to the reference genome FASTA. Default build_default_reference_list()$HG19$reference$genome.
#' @param dbsnp Known-sites VCF(s) used by BaseRecalibrator. May be a vector. Default build_default_reference_list()$HG19$database$all_common.
#' @param bam Input BAM file to recalibrate. If missing, the function expects a value passed into the internal environment as `input`.
#' @param region Optional region string to restrict processing (e.g. "chr1:1000-2000"). If provided, the job will be limited to this region.
#' @param ... Additional arguments forwarded to the internal launch environment (e.g. `tmp_dir`, `out_file_dir`, `executor_id`, `batch_config`, `threads`, `ram`, `mode`, `time`, `update_time`, `wait`, `hold`, `verbose`).
#' @return A job report list. The generated recalibration table path is available in `out_files$recal_table`.
#' @details The function builds and runs a command equivalent to:
#' \code{singularity run <sif_gatk> /gatk/gatk BaseRecalibrator -I <bam> -R <ref_genome> --known-sites <dbsnp> -O <out.recal.table>}.
#' Multiple known-sites VCFs are passed as repeated \code{--known-sites} arguments.
#' @export




new_generate_BQSR_gatk=function(
  sif_gatk=build_default_sif_list()$sif_gatk,
  ref_genome=build_default_reference_list()$HG19$reference$genome,
  dbsnp=build_default_reference_list()$HG19$database$all_common,
  bam=NULL,
  region=NULL,
  ...
){  

   run_main=function(
      .env
    ){
      
      .this.env=environment()
      append_env(to=.this.env,from=.env)

      set_main(.env=.this.env)

      if(!is.null(tmp_dir)){
        tmp_dir=paste0(" --tmp-dir ",tmp_dir)
      }

  
      reg=paste0(" -L ",input)
      .main$out_files$recal_table=paste0(out_file_dir,"/",get_file_name(bam),".",input,".recal.table")


      ## Multiple vcf with snps can be given

      if (!is.null(dbsnp)){
        dbsnp=paste(" --known-sites ",dbsnp,collapse=" ")
      }

      .main$exec_code=paste0(
        "singularity run ",sif_gatk, " /gatk/gatk BaseRecalibrator -I ", bam,
        " -R ", ref_genome, 
        dbsnp, reg,
        " -O ", .main$out_files$recal_table,tmp_dir

      )



        
      run_job(.env=.this.env)
      .env$.main <- .main

    }



    .base.env=environment()
    list2env(list(...),envir=.base.env)

    set_env_vars(
      .env= .base.env,
      vars="region"
    )
    

    launch(.env=.base.env)
 
}





#' Apply GATK Base Quality Score Recalibration (BQSR)
#'
#' Run GATK's `ApplyBQSR` using a Singularity image. The function builds
#' the container command to apply a recalibration table to a BAM file and
#' writes the recalibrated BAM to the task output directory. When `region`
#' is provided the function will apply BQSR to the specified region and
#' name the output accordingly.
#'
#' @param sif_gatk Path to the Singularity image for GATK. Defaults to
#'   `build_default_sif_list()$sif_gatk`.
#' @param ref_genome Reference genome fasta path. Defaults to
#'   `build_default_reference_list()$HG19$reference$genome`.
#' @param dbsnp Path to dbSNP VCF (used by other BQSR steps). Defaults to
#'   `build_default_reference_list()$HG19$database$all_common`.
#' @param bam Path to the input BAM to recalibrate.
#' @param region Optional region string (e.g. "chr1:10000-20000"). If
#'   provided the command will limit processing to this interval.
#' @param rec_table Path or vector of recalibration table(s). When `region`
#'   is given a matching entry will be looked up from `rec_table`.
#'
#' @return Invisibly returns the job report created by the internal job
#'   runner and writes the recalibrated BAM (and/or recalibration table)
#'   to the task output directory.
#' @export

new_apply_BQSR_gatk=function(
  sif_gatk=build_default_sif_list()$sif_gatk,
  ref_genome=build_default_reference_list()$HG19$reference$genome,
  dbsnp=build_default_reference_list()$HG19$database$all_common,
  bam=NULL,
  region=NULL,
  rec_table=NULL,
  ...
  ){

    run_main=function(
      .env
    ){
      
      .this.env=environment()
      append_env(to=.this.env,from=.env)

      set_main(.env=.this.env)

      if(!is.null(tmp_dir)){
        tmp_dir=paste0(" --tmp-dir ",tmp_dir)
      }


      reg=paste0(" -L ",input)
      recal=rec_table[grepl(input,rec_table)]
      .main$out_files$recal_bam=paste0(out_file_dir,"/",get_file_name(bam),".",input,".recal.",get_file_ext(bam))
  

      .main$exec_code=paste0(
        "singularity run ",sif_gatk, " /gatk/gatk ApplyBQSR -I ", bam,
        " -R ", ref_genome, 
        " --bqsr-recal-file ", recal ,reg,
        " -O ", .main$out_files$recal_bam,tmp_dir
      )

        
      run_job(.env=.this.env)
      .env$.main <- .main

    }

    .base.env=environment()
    list2env(list(...),envir=.base.env)

    set_env_vars(
      .env= .base.env,
      vars="region"
    )
    

    launch(.env=.base.env)

}





