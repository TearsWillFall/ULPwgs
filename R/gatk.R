
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
  mode="local",time="48:0:0",update_time=60,wait=FALSE,hold=""){

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
    exec_code=paste0(" singularity exec -H ",getwd(),":/home ",sif_gatk," /gatk/gatk MarkDuplicatesSpark -I ",bam, " -O ",  out_file,
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
  update_time=60,wait=FALSE,hold=""
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
  time="48:0:0",update_time=60,wait=FALSE,hold=""
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

  exec_code=paste0("singularity exec -H ",getwd(),":/home ",sif_gatk," /gatk/gatk BaseRecalibrator -I ",bam, " -R ", ref_genome, dbsnp,
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
  regions="",ref_genome="", clean=TRUE, dbsnp="",
  tmp_dir=".",threads=3,ram=4,
  executor_id=make_unique_id("par_generateBQSR"),
  task_name="par_generateBQSR",output_dir=".",
  verbose=FALSE,batch_config=build_default_preprocess_config(),
  mode="local",time="48:0:0",
  update_time=60,wait=FALSE,hold=""){

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

  if(regions==""){

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
        save(region_list,sif_gatk,bam,ref_genome,dbsnp,output_dir,verbose,tmp_dir,file = rdata_file)
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
  threads=4,ram=4,update_time=60,wait=FALSE,hold=""){

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir)

  if(tmp_dir!=""){
    tmp_dir=paste0(" --tmp-dir ",tmp_dir)
  }
  out_file=paste0(out_file_dir,output_name,".recal.table")
  exec_code=paste0("singularity exec -H ",getwd(),":/home " ,sif_gatk,
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
  update_time=60,wait=TRUE,hold=""){


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
  exec_code=paste("singularity exec -H ",paste0(getwd(),":/home "),sif_gatk,
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
  bam="",regions="",ref_genome="",rec_table="",clean=TRUE,
  output_dir=".",verbose=FALSE,tmp_dir=".",
  batch_config=build_default_preprocess_config(),mode="local",
  executor_id=make_unique("par_applyBQSR"),
  task_name="par_applyBQSR",
  time="48:0:0",threads=4,ram=4,
  update_time=60,wait=FALSE, hold=""){

  options(scipen = 999)
 

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir)

   if(regions==""){
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
        save(region_list,sif_gatk,bam,ref_genome,rec_table,output_dir,verbose,tmp_dir,file = rdata_file)
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
  
  job_report[["steps"]][["gather_bam"]]=gather_bam_files(
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





#' Wrapper around gatk GatherBamFiles function
#'
#' This functions collects the Recalibration reports generated from scattered parallel_apply_BQSR output
#' This function wraps around gatk GatherBamFiles function.
#' For more information about this function: https://gatk.broadinstitute.org/hc/en-us/articles/360037055512-GatherBamFiles-Picard-
#'
#' @param bin_picard [REQUIRED] Path to picard executable. Default tools/picard/build/libs/picard.jar.
#' @param bam [REQUIRED] Path to BAM file/s.
#' @param output_name [OPTIONAL] Name for the output file name.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param clean Clean input files. Default TRUE.
#' @param tmp_dir [OPTIONAL] Path to the temporary directory.
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param threads [OPTIONAL] Number of threads to split the work. Default 4
#' @param ram [OPTIONAL] RAM memory to asing to each thread. Default 4
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Task EXECUTOR ID. Default "gatherBAM"
#' @param task_name Task name. Default "gatherBAM"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Show job updates every update time. Default 60
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export

gather_bam_files=function(
  bin_picard=build_default_tool_binary_list()$bin_picard,
  bam="",output_name="File",output_dir=".",
  verbose=FALSE,batch_config=build_default_preprocess_config(),
  threads=4, ram=4, mode="local",clean=FALSE,
  executor_id=make_unique_id("gatherBAM"),task_name="gatherBAM",
  time="48:0:0",update_time=60,wait=FALSE,hold=""
){

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir)

  
  out_file=paste0(out_file_dir,"/",output_name,".bam")
  bam=bam[order(as.numeric(lapply(lapply(lapply(lapply(lapply(lapply(basename(bam),
    FUN=strsplit,split="\\."),FUN="[[",index=1),FUN="[",index=2),
    FUN=strsplit,split="__"),FUN="[[",index=1),FUN="[",index=1)))]


  exec_code=paste0("java -jar ",bin_picard," GatherBamFiles ",
    paste0(" I=",bam,collapse=" ")," O=",out_file)

  if(clean){
    exec_code=paste(exec_code," && rm",paste(bam,collapse=" "),
    " && rm ",paste0(gsub(".bam",".bai",bam),collapse=" "))
  }

 
  job=build_job(executor_id=executor_id,task_id=task_id)

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
    stop("picard failed to run due to unknown error.
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
      bam=out_file)
    )

  if(wait&&mode=="batch"){
    job_validator(job=job_report$job_id,time=update_time,
    verbose=verbose,threads=threads)
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
  update_time=60,wait=FALSE,hold=""
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

    exec_code=paste0("singularity exec -H ",getwd(),":/home ",sif_gatk,
    " /gatk/gatk  AnalyzeCovariates -bqsr ",before, " -plots ",out_file,tmp_dir," -csv ",out_file_csv)
  }else if(before=="" & after!=""){

    out_file=paste0(out_file_dir,"/",get_file_name(after),
       "_covariates_analysis_after.pdf")

    out_file_csv=paste0(out_file_dir,"/",get_file_name(after),
    "_covariates_analysis_after.csv")

    exec_code=paste0("singularity exec -H ",getwd(),":/home ",sif_gatk,
    " /gatk/gatk AnalyzeCovariates -bqsr ",after, " -plots ",out_file,tmp_dir," -csv ",out_file_csv)
  }else{
    out_file=paste0(out_file_dir,"/",get_file_name(before),"_covariates_analysis.pdf")
    out_file_csv=paste0(out_file_dir,"/",get_file_name(before),"_covariates_analysis.csv")
    exec_code=paste0("singularity exec -H ",getwd(),":/home ",sif_gatk,
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
  tumour="",normal="",output_name="",
  ref_genome=build_default_reference_list()$HG19$reference,
  germ_resource=build_default_reference_list()$HG19$variant$germ_reference,
  pon="",output_dir=".",tmp_dir=".",
  verbose=FALSE,orientation=FALSE,mnps=FALSE,
  batch_config=build_default_preprocess_config(),
  threads=4,ram=4,mode="local",
  executor_id=make_unique_id("Mutect2"),
  task_name="Mutect2",time="48:0:0",
  update_time=60,wait=FALSE,hold=""
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
  

  if(output_name!=""&length(tumour)>1){
    fname=output_name
  }else{  
    fname=get_file_name(tumour[1])
  }
  
  reg=""
  if (region==""){
      out_file=paste0(out_file_dir,"/",fname,".unfilt.vcf")
  }else{
      reg=paste0(" -L ",region)
      fname=paste0(fname,".",region)
      out_file=paste0(out_file_dir,"/",fname,".unfilt.vcf")
  }

  if (is.vector(tumour)){
    tumour=paste0(" -I ",paste(tumour,collapse=" -I "))
  }else{
    tumour=paste0(" -I ",tumour)
  }
  norm=" "
  if (normal!=""){
    if (is.vector(normal)){
      norm=paste0(" -I ",paste(normal,collapse=" -I ")," -normal ",
      paste(Vectorize(get_file_name)(normal),collapse=" -normal "))
  }else{
      norm=paste0(" -I ",normal," -normal ",get_file_name(normal))
      }
  }

  if (pon!=""){
    pon=paste0(" --panel-of-normals ",pon)
  }

  f1r2=""
  if (orientation){
      out_file_dir_ort=set_dir(dir=output_dir,name="orientation")
      out_file2=paste0(out_file_dir_ort,"/",fname,".f1r2.tar.gz")
      f1r2=paste0(" --f1r2-tar-gz ",out_file2)
  }



  filter_mnps=""
  if (mnps){
    filter_mnps=" -max-mnp-distance 0 "
  }

  exec_code=paste("singularity exec -H ",paste0(getwd(),":/home "),sif_gatk,
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
  bin_samtools=build_default_tool_binary_list()$bin_samtools,
  bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
  bin_tabix=build_default_tool_binary_list()$bin_tabix,
  tumour="",normal="",output_name="",
  ref_genome=build_default_reference_list()$HG19$reference,
  germ_resource=build_default_reference_list()$HG19$variant$germ_reference,
  biallelic_db=build_default_reference_list()$HG19$variant$biallelic_reference,
  db_interval=build_default_reference_list()$HG19$variant$biallelic_reference,
  regions="",pon="",output_dir=".",
  verbose=FALSE,filter=TRUE,
  orientation=TRUE,mnps=FALSE,
  contamination=TRUE,clean=FALSE,
  batch_config=build_default_preprocess_config(),
  threads=4,ram=4,mode="local",
  executor_id=make_unique_id("parRegionMutect2"),
  task_name="parRegionMutect2",time="48:0:0",
  update_time=60,wait=FALSE,hold=""
){


  if(!is.null(rdata)){
    load(rdata)
    if(!is.null(selected)){
      tumour=tumours_list[selected]
      output_name=get_file_name(tumour)
    }
  }

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir,name=paste0(output_name,"/mutect2_reports"))
  tmp_dir=set_dir(dir=out_file_dir,name="mutect2_tmp")

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

  if(regions==""){

    jobs_report[["steps"]][["getChr"]] <- get_bam_reference_chr(
      bin_samtools=bin_samtools,
      bam=tumour[1],verbose=verbose,output_dir=tmp_dir,
      output_name=paste0(get_file_name(tumour[1]),"_Ref"),
      executor_id=task_id,mode="local",threads=threads,ram=ram,
      time=time,update_time=update_time,wait=FALSE,hold=hold)


    regions=read.table(jobs_report[["steps"]][["getChr"]]$out_files$ref,
    sep="\t",header=TRUE)
    hold=jobs_report[["steps"]][["getChr"]]$job_id
  }

  regions$start=regions$start+1
  regions=regions %>% dplyr::mutate(region=paste0(chr,":",start,"-",end))


  region_list=regions$region
  names(region_list)=regions$region

  if(mode=="local"){
    jobs_report[["steps"]][["par_region_call_variants"]]<-
    parallel::mclapply(region_list,FUN=function(region){
      job_report <- mutect2_gatk(
            sif_gatk=sif_gatk,
            region=region,
            tumour=tumour,
            normal=normal,
            output_name=output_name,
            ref_genome=ref_genome,
            germ_resource = germ_resource,
            pon=pon,output_dir=tmp_dir,tmp_dir=tmp_dir,
            verbose=verbose,orientation=orientation,mnps=mnps,
            executor_id=task_id)
    },mc.cores=threads)
    
    }else if(mode=="batch"){
          rdata_file=paste0(tmp_dir,"/",job,".regions.RData")
          output_dir=tmp_dir
          save(region_list,tumour,normal,sif_gatk,ref_genome,output_name,
          germ_resource,pon,orientation,mnps,output_dir,verbose,tmp_dir,file = rdata_file)
          exec_code=paste0("Rscript -e \"ULPwgs::mutect2_gatk(rdata=\\\"",
          rdata_file,"\\\",selected=$SGE_TASK_ID)\"")
          out_file_dir2=set_dir(dir=out_file_dir,name="batch")
          batch_code=build_job_exec(job=job,time=time,ram=ram,
          threads=1,output_dir=out_file_dir2,
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
           if(output_name!=""&length(tumour)>1){
              fname=output_name
           }else{  
              fname=get_file_name(tumour[1])
           }
  
         jobs_report[["steps"]][["par_region_call_variants"]]<- build_job_report(
              job_id=job,
              executor_id=executor_id,
              exec_code=exec_code, 
              task_id=task_id,
              input_args=argg,
              out_file_dir=tmp_dir,
              out_files=list(
                  vcf=paste0(tmp_dir,"/vcf/",fname,".",region_list,".unfilt.vcf"),
                  stats=paste0(tmp_dir,"/vcf/",fname,".",region_list,".unfilt.vcf.stats"),
                  idx=paste0(tmp_dir,"/vcf/",fname,".",region_list,".unfilt.vcf.idx"),
                  f1r2=paste0(tmp_dir,"/orientation/",fname,".",region_list,".f1r2.tar.gz")
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
      vcfs=vcfs,stats=stats,
      f1r2=f1r2,
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

  vcf=unlist_lvl(jobs_report[["steps"]][["gatherFilesGatk"]],var="concat_vcf")
  stats=unlist_lvl(jobs_report[["steps"]][["gatherFilesGatk"]],var="merged_stats")
  if(orientation){
    orientation_model=unlist_lvl(jobs_report[["steps"]][["gatherFilesGatk"]],var="orientation_model")
  }else{
    orientation_model=""
  }
 

  if(contamination){
          jobs_report[["steps"]][["estimateContaminationGatk"]]<- parallel_estimate_contamination_gatk(
                sif_gatk=sif_gatk,
                tumours=tumour,
                normal=normal,
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
            vcf=vcf,stats=stats,contamination_table=contamination_table,
            segmentation_table=segmentation_table,
            orientation_model=orientation_model,output_name=output_name,
            ref_genome=ref_genome,
            output_dir=out_file_dir,
            verbose=verbose,clean=clean,
            batch_config=batch_config,
            threads=1,ram=ram,mode=mode,
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
#' @param bin_bgzip [REQUIRED] Path to bgzip binary file.
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
  bin_bcftools=build_default_tool_binary_list()$bin_bcftools,
  bin_samtools=build_default_tool_binary_list()$bin_samtools,
  bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
  bin_tabix=build_default_tool_binary_list()$bin_tabix,
  tumours="",normal="",
  ref_genome=build_default_reference_list()$HG19$reference,
  germ_resource=build_default_reference_list()$HG19$variant$germ_reference,
  biallelic_db=build_default_reference_list()$HG19$variant$biallelic_reference,
  db_interval=build_default_reference_list()$HG19$variant$biallelic_reference,
  regions="",pon="",output_dir=".",
  verbose=FALSE,filter=TRUE,
  orientation=TRUE,mnps=FALSE,
  contamination=TRUE,clean=FALSE,
  batch_config=build_default_preprocess_config(),
  threads=4,ram=4,mode="local",
  executor_id=make_unique_id("parRegionMutect2"),
  task_name="parRegionMutect2",time="48:0:0",
  update_time=60,wait=FALSE,hold=""
){

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir,name="mutect2_reports")
  tmp_dir=set_dir(dir=out_file_dir,name="mutect2_tmp")
 

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

  tumours_list=tumours
  names(tumours_list)=Vectorize(get_file_name)(tumours)

  if(mode=="local"){
    jobs_report[["steps"]][["par_sample_call_variants"]]<-
    parallel::mclapply(tumours_list,FUN=function(tumour){
      job_report <- parallel_regions_mutect2_gatk(
            sif_gatk=sif_gatk,
            bin_bcftools=bin_bcftools,
            bin_samtools=bin_samtools,
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            regions=regions,
            output_name=get_file_name(tumour),
            ref_genome=ref_genome,
            germ_resource=germ_resource,
            biallelic_db=biallelic_db,
            db_interval=db_interval,
            pon=pon,
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
          save(tumours_list,normal,
          sif_gatk,
          bin_bcftools,
          bin_samtools,
          bin_bgzip,
          bin_tabix,
          germ_resource,
          biallelic_db,db_interval,pon,ref_genome,
          filter,orientation,mnps,
          mode,ram,threads, time,
          output_dir,contamination,clean,
          verbose,file = rdata_file)
          exec_code=paste0("Rscript -e \"ULPwgs::parallel_regions_mutect2_gatk(rdata=\\\"",
          rdata_file,"\\\",selected=$SGE_TASK_ID)\"")
          out_file_dir2=set_dir(dir=out_file_dir,name="batch")
          batch_code=build_job_exec(job=job,time=time,ram=ram,
          threads=1,output_dir=out_file_dir2,
          hold=hold,array=length(tumours_list))
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
                )
              )

  }


  if(wait&&mode=="batch"){
    job_validator(job=unlist_lvl(jobs_report[["steps"]],var="job_id"),time=update_time,
    verbose=verbose,threads=threads)
  }

  return(jobs_report)

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
  vcf="",stats="",contamination_table="",
  segmentation_table="",orientation_model="",output_name="",
  ref_genome=build_default_reference_list()$HG19$reference,
  output_dir=".",verbose=FALSE,clean=FALSE,
  batch_config=build_default_preprocess_config(),
  threads=4,ram=4,mode="local",
  executor_id=make_unique_id("filterMutect2Gatk"),
  task_name="filterMutect2Gatk",time="48:0:0",
  update_time=60,wait=FALSE,hold=""
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

  if(stats!=""){
      stats=paste0(" -stats ",stats)
  }

  if (contamination_table!=""){
     contamination_table=paste0(" --contamination-table ", paste0(contamination_table,
     collapse=" --contamination-table "))
  }
  
  if(segmentation_table!=""){
      segmentation_table=paste0(" --tumor-segmentation ",paste0(segmentation_table,
      collapse=" --tumor-segmentation "))
    }

  if(orientation_model!=""){
      orientation_model=paste0(" --ob-priors ",orientation_model)
  }
  
  out_file=paste0(out_file_dir,"/",id,".filtered.vcf")


  exec_code=paste0("singularity exec -H ",getwd(),":/home ",sif_gatk,
  " /gatk/gatk   FilterMutectCalls -R ",ref_genome," -O ",out_file," -V ",vcf, stats,
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
  vcfs="",stats="",f1r2="",output_dir=".",tmp_dir=".",
  output_name="",verbose=FALSE,orientation=FALSE,clean=TRUE,
  batch_config=build_default_preprocess_config(),
  threads=4,ram=4,mode="local",
  executor_id=make_unique_id("gatherMutect"),
  task_name="gatherMutect",time="48:0:0",
  update_time=60,wait=FALSE,hold=""
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
    vcfs=vcfs,clean=clean,sort=TRUE,
    compress=FALSE,format="vcf",
    verbose=verbose,output_name=output_name,
    output_dir=out_file_dir,
    tmp_dir=tmp_dir,batch_config=batch_config,
    threads=1,ram=ram,mode=mode,
    executor_id=task_id,
    time=time,hold=hold
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
  update_time=60,wait=FALSE,hold=""
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
 
  if (f1r2!=""){
    f1r2_list=paste0(" -I ",paste0(f1r2,collapse=" -I "))
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
  tumour="",normal="",tumour_pileup="",
  normal_pileup="",output_name="",output_dir=".",tmp_dir=".",
  biallelic_db=build_default_reference_list()$HG19$variant$biallelic_reference,
  db_interval=build_default_reference_list()$HG19$variant$biallelic_reference,
  verbose=FALSE,batch_config=build_default_preprocess_config(),
  threads=1,ram=4,mode="local",
  executor_id=make_unique_id("estimateContaminationGatk"),
  task_name="estimateContaminationGatk",time="48:0:0",
  update_time=60,wait=FALSE,hold=""){

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

  if(normal!=""){
      jobs_report[["steps"]][["nPileupGatk"]]<-pileup_summary_gatk(
        sif_gatk=sif_gatk,
        bam=normal,output_name=get_file_name(normal),
        output_dir=paste0(out_file_dir,"/pileup_reports/normal"),
        verbose=verbose,batch_config=batch_config,
        biallelic_db=biallelic_db,
        db_interval=db_interval,
        threads=1,ram=ram,mode=mode,
        executor_id=task_id,hold=hold
      )
      normal_pileup<-jobs_report[["steps"]][["nPileupGatk"]]$out_files$pileup_table

  }

  if(tumour!=""){
        jobs_report[["steps"]][["tPileupGatk"]]<-pileup_summary_gatk(
          sif_gatk=sif_gatk,
          bam=tumour,output_name=get_file_name(tumour),
          output_dir=paste0(out_file_dir,"/pileup_reports/tumour"),
          verbose=verbose,batch_config=batch_config,
          biallelic_db=biallelic_db,
          db_interval=db_interval,
          threads=1,ram=ram,mode=mode,
          executor_id=task_id,hold=hold
        )
        tumour_pileup<-jobs_report[["steps"]][["tPileupGatk"]]$out_files$pileup_table
  }


  
  if (normal_pileup!=""){
        normal_pileup=paste0(" -matched ",normal_pileup)
  }




  exec_code=paste0("singularity exec -H ",getwd(),":/home ",sif_gatk,
  " /gatk/gatk   CalculateContamination  -O ",out_file, " -tumor-segmentation ", out_file2,
  " -I ",tumour_pileup,normal_pileup)

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
  update_time=60,wait=FALSE,hold=""
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
            tumour = tumour,
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
  update_time=60,wait=FALSE,hold=""
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
  " /gatk/gatk  GetPileupSummaries  -O ",out_file," -V ",biallelic_db," -L ",
   db_interval, " -I ",bam)


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
  update_time=60,wait=FALSE,hold=""){


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
          bam=bam,output_name=get_file_name(bam),output_dir=out_file_dir,
          verbose=verbose,batch_config=batch_config,
          biallelic_db=biallelic_db,
          db_interval=db_interval,
          executor_id=task_id,
          time=time,hold=hold)
    },mc.cores=threads)
    }else if(mode=="batch"){
          rdata_file=paste0(tmp_dir,"/",job,".tumours.RData")
          output_dir=out_file_dir
          save(bam_list,sif_gatk,output_dir,biallelic_db,db_interval,verbose,tmp_dir,file = rdata_file)
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
  stats="",output_name="",output_dir=".",clean=FALSE,
  verbose=TRUE, batch_config=build_default_preprocess_config(),
  threads=1,ram=4,mode="local",
  executor_id=make_unique_id("mergeMutectStats"),
  task_name="mergeMutectStats",time="48:0:0",
  update_time=60,wait=FALSE,hold=""){

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
 
  if (stats!=""){
    stat=paste0(" -stats ",paste0(stats,collapse=" -stats "))
  }
  
  out_file=paste0(out_file_dir,"/",id,".merged.vcf.stats")

  exec_code=paste0("singularity exec -H ",getwd(),":/home ",sif_gatk,
  " /gatk/gatk  MergeMutectStats -O ", out_file ,stat)

  if(clean){
    exec_code=paste(exec_code," && rm",paste(stats,collapse=" "))
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
#' @param vcfs [OPTIONAL] Path to VCFs files. Only required if BAM files are not given.
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
  normals="",output_name="PoN",
  vcfs="",
  regions="",output_dir=".",
  ref_genome=build_default_reference_list()$HG19$reference,
  germ_resource=build_default_reference_list()$HG19$variant$germ_reference,
  verbose=FALSE,
  batch_config=build_default_preprocess_config(),
  threads=4,ram=4,mode="local",
  executor_id=make_unique_id("createPoNGatk"),
  task_name="createPoNGATK",time="48:0:0",
  update_time=60,wait=FALSE,hold=""
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

  ### Process NORMAL BAMS if vcfs are not given otherwise use vcfs
  if(normals!=""){

    jobs_report[["steps"]][["par_samples_mutect2"]]<-parallel_samples_mutect2_gatk(
      sif_gatk=sif_gatk,
      bin_bcftools=bin_bcftools,
      bin_samtools=bin_samtools,
      bin_bgzip=bin_bgzip,
      bin_tabix=bin_tabix,
      tumours=normals,
      ref_genome=ref_genome,
      germ_resource=germ_resource,
      biallelic_db="",
      db_interval="",
      output_dir=out_file_dir,
      verbose=verbose,filter=FALSE,
      orientation=FALSE,mnps=FALSE,
      contamination=FALSE,clean=FALSE,
      batch_config=batch_config,
      threads=threads,ram=ram,mode=mode,
      executor_id=task_id,
      time=time,
      hold=hold
    )

    vcfs=unlist_lvl(jobs_report[["steps"]][["par_samples_mutect2"]],var="concat_vcf")
    hold=unlist_lvl(jobs_report[["steps"]][["par_samples_mutect2"]],var="job_id")
  }

  vcfs=paste0(" --vcfs ",paste0(vcfs,collapse=" --vcfs "))

  exec_code=paste("singularity exec -H ",paste0(getwd(),":/home "),sif_gatk,
  " /gatk/gatk  CreateSomaticPanelOfNormals  ",vcfs," -O ",out_file)

  if(mode=="batch"){
       out_file_dir2=set_dir(dir=out_file_dir,name="batch")
       batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
       threads=threads,output_dir=out_file_dir2)
       exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
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
    ref_genome=build_default_reference_list()$HG19$reference,
    output_dir=".",
    tmp_dir=".",
    output_name="",
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=4,ram=4,mode="local",
    executor_id=make_unique_id("parallelCallHaplotypecallerGatk"),
    task_name="parallelCallHaplotypecallerGatk",time="48:0:0",
    update_time=60,wait=FALSE,hold=""
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
    tmp_dir=paste0(" --tmp-dir ",tmp_dir)
  }
  

  if(output_name!=""&length(normal)>1){
    fname=output_name
  }else{  
    fname=get_file_name(normal[1])
  }
  
  reg=""
  if (region==""){
      out_file=paste0(out_file_dir,"/",fname,".unfilt.vcf.gz")
  }else{
      reg=paste0(" -L ",region)
      fname=paste0(fname,".",region)
      out_file=paste0(out_file_dir,"/",fname,".unfilt.vcf.gz")
  }

  if (is.vector(normal)){
    normal=paste0(" -I ",paste(normal,collapse=" -I "))
  }else{
    normla=paste0(" -I ",normal)
  }


  exec_code=paste("singularity exec -H ",paste0(getwd(),":/home "),sif_gatk,
  " /gatk/gatk HaplotypeCaller -R ",ref_genome, normal," -O ",out_file,reg)

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


#' Parallel implementation per region for Variant Calling using HaplotypeCaller 
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
#' @param bin_bcftools [REQUIRED] Path to bcftools binary file.
#' @param bin_bgzip [REQUIRED] Path to bgzip binary file.
#' @param bin_tabix [REQUIRED] Path to tabix binary file.
#' @param bin_samtools [REQUIRED] Path to samtools binary file.
#' @param normal [OPTIONAL] Path to normal BAM file.
#' @param ref_genome [REQUIRED] Path to reference genome fasta file.
#' @param indel_db [OPTIONAL] Path to database with indel info.
#' @param haplotype_db [OPTIONAL] Path to database with haplotype info.
#' @param info_key [OPTIONAL] Info Key for CNN model training. Default CNN_1D. Options ["CNN_1D","CNN_2D"]
#' @param snp_tranche [OPTIONAL] SNP tranche cut-off.
#' @param indel_tranche [OPTIONAL] INDEL tranche cut-off.
#' @param keep_previous_tranche [OPTIONAL] Remove previous filter information.
#' @param rdata [OPTIONAL] Import R data information with list of BAM.
#' @param selected [OPTIONAL] Select BAM from list.
#' @param regions [OPTIONAL] Regions to parallelize through. If not given will be infered from BAM file.
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param threads [OPTIONAL] Number of threads to split the work. Default 4
#' @param filter [OPTIONAL] Filter variants. Default TRUE
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



parallel_regions_haplotypecaller_gatk=function(
  rdata=NULL,selected=NULL,
  sif_gatk=build_default_sif_list()$sif_gatk,
  bin_bcftools=build_default_tool_binary_list()$bin_bcftools,
  bin_samtools=build_default_tool_binary_list()$bin_samtools,
  bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
  bin_tabix=build_default_tool_binary_list()$bin_tabix,
  normal="",
  ref_genome=build_default_reference_list()$HG19$reference,
  regions="",
  output_dir=".",
  indel_db=build_default_reference_list()$HG19$variant$mills_reference,
  haplotype_db=build_default_reference_list()$HG19$variant$hapmap_reference,
  filter=TRUE,
  output_name="",
  info_key="CNN_1D",
  snp_tranche=99.95,
  indel_tranche=99.4,
  keep_previous_filters=FALSE,
  clean=FALSE,
  verbose=FALSE,
  batch_config=build_default_preprocess_config(),
  threads=4,ram=4,mode="local",
  executor_id=make_unique_id("parallelCallHaplotypecallerGatk"),
  task_name="parallelCallHaplotypecallerGatk",time="48:0:0",
  update_time=60,wait=FALSE,hold=""
){
  
  if(!is.null(rdata)){
    load(rdata)
    if(!is.null(selected)){
      tumour=tumours_list[selected]
      output_name=get_file_name(tumour)
    }
  }

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir,name=paste0(output_name,"/haplotypecaller_reports"))
  tmp_dir=set_dir(dir=out_file_dir,name="haplotypecaller_tmp")

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

  if(regions==""){

    jobs_report[["steps"]][["getChr"]] <- get_bam_reference_chr(
      bin_samtools=bin_samtools,
      bam=normal[1],verbose=verbose,output_dir=tmp_dir,
      output_name=paste0(get_file_name(normal[1]),"_Ref"),
      executor_id=task_id,mode="local",threads=threads,ram=ram,
      time=time,update_time=update_time,wait=FALSE,hold=hold)

    regions=read.table(jobs_report[["steps"]][["getChr"]]$out_files$ref,
    sep="\t",header=TRUE)
    hold=jobs_report[["steps"]][["getChr"]]$job_id
  }

  regions$start=regions$start+1
  regions=regions %>% dplyr::mutate(region=paste0(chr,":",start,"-",end))


  region_list=regions$region
  names(region_list)=regions$region

  if(mode=="local"){
    jobs_report[["steps"]][["par_region_call_variants"]]<-
    parallel::mclapply(region_list,FUN=function(region){
      job_report <- haplotypecaller_gatk(
            sif_gatk=sif_gatk,
            region=region,
            normal=normal,
            output_name=output_name,
            ref_genome=ref_genome,
            output_dir=tmp_dir,
            tmp_dir=tmp_dir,
            verbose=verbose,
            executor_id=task_id
      )
    },mc.cores=threads)
    
    }else if(mode=="batch"){
          rdata_file=paste0(tmp_dir,"/",job,".regions.RData")
          output_dir=tmp_dir
          save(region_list,normal,sif_gatk,ref_genome,output_name,
          output_dir,verbose,tmp_dir,file = rdata_file)
          exec_code=paste0("Rscript -e \"ULPwgs::haplotypecaller_gatk(rdata=\\\"",
          rdata_file,"\\\",selected=$SGE_TASK_ID)\"")
          out_file_dir2=set_dir(dir=out_file_dir,name="batch")
          batch_code=build_job_exec(job=job,time=time,ram=ram,
          threads=1,output_dir=out_file_dir2,
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
           if(output_name!=""&length(normal)>1){
              fname=output_name
           }else{  
              fname=get_file_name(normal[1])
           }
  
         jobs_report[["steps"]][["par_region_call_variants"]]<- build_job_report(
              job_id=job,
              executor_id=executor_id,
              exec_code=exec_code, 
              task_id=task_id,
              input_args=argg,
              out_file_dir=tmp_dir,
              out_files=list(
                  vcf=paste0(tmp_dir,"/vcf/",fname,".",region_list,".unfilt.vcf.gz"),
                  idx=paste0(tmp_dir,"/vcf/",fname,".",region_list,".unfilt.vcf.idx")
                )
        )

  }


  vcfs<-unlist_lvl(jobs_report[["steps"]][["par_region_call_variants"]],var="vcf")
  hold<-unlist_lvl(jobs_report[["steps"]][["par_region_call_variants"]],var="job_id")
  jobs_report[["steps"]][["concatVCF"]] <- concat_vcf(
    bin_bcftools=bin_bcftools,
    bin_bgzip=bin_bgzip,
    bin_tabix=bin_tabix,
    vcfs=vcfs,clean=clean,sort=TRUE,
    compress=FALSE,format="vcf",
    verbose=verbose,output_name=output_name,
    output_dir=out_file_dir,
    tmp_dir=tmp_dir,batch_config=batch_config,
    threads=1,ram=ram,mode=mode,
    executor_id=task_id,
    time=time,hold=hold
  )

  if(filter){
        vcf<-jobs_report[["steps"]][["concatVCF"]]$out_files$concat_vcf
        hold<-jobs_report[["steps"]][["concatVCF"]]$job_id
        jobs_report[["steps"]][["cnnScoreVariantsGatk"]]<-cnn_score_variants_gatk(
        sif_gatk=sif_gatk,
        vcf=vcf,bam=ifelse(info_key=="CNN_1D","",bam),
        ref_genome=ref_genome,
        output_dir=out_file_dir,output_name=output_name,
        verbose=verbose,
        batch_config=batch_config,
        threads=threads,ram=ram,mode=mode,
        executor_id=task_id,
        time=time,
        hold=hold
    )

      vcf<-jobs_report[["steps"]][["cnnScoreVariantsGatk"]]$out_files$scored_vcf
      hold<-jobs_report[["steps"]][["cnnScoreVariantsGatk"]]$job_id

      jobs_report[["steps"]][["filterVariantTranchesGatk"]]<-filter_variant_tranches_gatk(
        sif_gatk=sif_gatk,
        vcf=vcf,
        ref_genome=ref_genome,
        indel_db=indel_db,
        haplotype_db=haplotype_db,
        output_dir=out_file_dir,output_name=output_name,
        info_key=info_key,
        snp_tranche=snp_tranche,
        indel_tranche=indel_tranche,
        keep_previous_filters=keep_previous_filters,
        verbose=verbose,
        batch_config=batch_config,
        threads=threads,ram=ram,mode=mode,
        executor_id=task_id,
        time=time,
        hold=hold
    )

  }

  return(jobs_report)

}


#' CNN to score germline variants in VCF
#'
#' This function functions calls CNNScoreVariants for a VCF files and scores
#' the variants using a predictive CNN model
#' 
#' For more information read:
#' https://gatk.broadinstitute.org/hc/en-us/articles/360037226672-CNNScoreVariants
#'
#' @param sif_gatk [REQUIRED] Path to gatk sif file.
#' @param vcf [REQUIRED] Path to VCF file.
#' @param bam [OPTIONAL] Path to BAM file.
#' @param ref_genome [REQUIRED] Path to reference genome fasta file.
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
#' @param output_dir [OPTIONAL] Path to the output directory.
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




cnn_score_variants_gatk=function(

  sif_gatk=build_default_sif_list()$sif_gatk,
  vcf="",bam="",
  ref_genome=build_default_reference_list()$HG19$reference,
  output_dir=".",output_name="",
  verbose=FALSE,
  batch_config=build_default_preprocess_config(),
  threads=4,ram=4,mode="local",
  executor_id=make_unique_id("CNNScoreVariantsGatk"),
  task_name="CNNScoreVariantsGatk",time="48:0:0",
  update_time=60,wait=FALSE,hold=""
){

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir,name="cnnscores")
  job=build_job(executor_id=executor_id,task_id=task_id)
  

  id=""
  if(output_name!=""){
    id=output_name
  }else{
    id=get_file_name(vcf[1])
  }



  opt=""
  if (bam!=""){
    bam=paste0(" -I ",bam)
    opt=" -tensor-type read_tensor "
    out_file=paste0(out_file_dir,"/",id,".CNNscored.2D.vcf")
  }else{
    out_file=paste0(out_file_dir,"/",id,".CNNscored.1D.vcf")
  }

  exec_code=paste("singularity exec -H ",paste0(getwd(),":/home "),sif_gatk,
  " /gatk/gatk CNNScoreVariants -R ",ref_genome, bam," -O ",out_file, " -V ",vcf, opt, bam)

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
      scored_vcf=out_file,
      idx=paste0(out_file,".idx")
    )
  )


  if(wait&&mode=="batch"){
    job_validator(job=job_report$job_id,time=update_time,
    verbose=verbose,threads=threads)
  }

  return(job_report)

}




#' Variant Trench Filtering using GATK
#'
#' This function filters variant tranches calls generated using CNNScoreVariants method
#' 
#' For more information read:
#' https://gatk.broadinstitute.org/hc/en-us/articles/360051308071-FilterVariantTranches
#'
#' @param sif_gatk [REQUIRED] Path to gatk sif file.
#' @param vcf [OPTIONAL] Path to VCF file.
#' @param ref_genome [REQUIRED] Path to reference genome fasta file.
#' @param indel_db [OPTIONAL] Path to database with indel info.
#' @param haplotype_db [OPTIONAL] Path to database with haplotype info.
#' @param info_key [OPTIONAL] Info Key for CNN model training. Default CNN_1D. Options ["CNN_1D","CNN_2D"]
#' @param snp_tranche [OPTIONAL] SNP tranche cut-off.
#' @param indel_tranche [OPTIONAL] INDEL tranche cut-off.
#' @param keep_previous_tranche [OPTIONAL] Remove previous filter information.
#' @param regions [OPTIONAL] Regions to parallelize through. If not given will be infered from BAM file.
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
#' @param output_dir [OPTIONAL] Path to the output directory.
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



filter_variant_tranches_gatk=function(
  sif_gatk=build_default_sif_list()$sif_gatk,
  vcf="",
  ref_genome=build_default_reference_list()$HG19$reference,
  indel_db=build_default_reference_list()$HG19$variant$mills_reference,
  haplotype_db=build_default_reference_list()$HG19$variant$hapmap_reference,
  output_dir=".",output_name="",
  info_key="CNN_1D",
  snp_tranche=99.95,
  indel_tranche=99.4,
  keep_previous_filters=FALSE,
  verbose=FALSE,
  batch_config=build_default_preprocess_config(),
  threads=4,ram=4,mode="local",
  executor_id=make_unique_id("FilterVariantTranchesGatk"),
  task_name="FilterVariantTranchesGatk",time="48:0:0",
  update_time=60,wait=FALSE,hold=""
){


  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir)
  job=build_job(executor_id=executor_id,task_id=task_id)


  id=""
  if(output_name!=""){
    id=output_name
  }else{
    id=get_file_name(vcf[1])
  }


  prev_filters=" --invalidate-previous-filters "
  if (keep_previous_filters){
    prev_filters=" "
  }

  if(indel_db!=""){
    indel_db=paste0(" --resource ",indel_db)
  }

  if(haplotype_db!=""){
    haplotype_db=paste0(" --resource ",haplotype_db)
  }



  out_file=paste0(out_file_dir,"/",id,".filteredTranches.vcf")

  
  exec_code=paste(
    "singularity exec -H ",paste0(getwd(),":/home "),sif_gatk,
    " /gatk/gatk  FilterVariantTranches -R ",ref_genome, 
    " -O ", out_file, 
    " -V ",vcf, " --info-key ",info_key,
    " --snp-tranche ",snp_tranche,
    " --indel-tranche ",indel_tranche,indel_db,
    haplotype_db,prev_filters
  )

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
      scored_vcf=out_file,
      idx=paste0(out_file,".idx")
    )
  )


  if(wait&&mode=="batch"){
    job_validator(job=job_report$job_id,time=update_time,
    verbose=verbose,threads=threads)
  }

  return(job_report)

}



