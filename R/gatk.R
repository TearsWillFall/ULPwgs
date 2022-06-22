
#' Wrapper for MarkDuplicatesSpark from gatk
#'
#' This function removes duplicated reads (artifacts) found in aligned sequences and sorts the output bam.
#'
#' @param bam Path to the input file with the aligned sequence.
#' @param bin_path Path to picard executable. Default path tools/picard/build/libs/picard.jar.
#' @param output_dir Path to the output directory.
#' @param tmp_dir Path to tmp directory.
#' @param remove_duplicates Remove all sequencing duplicates from BAM file. Default TRUE.
#' @param threads Number of threads . Default 4
#' @param ram RAM memory. Default 4
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor Name of the executor. Default "mardupsGATK"
#' @param task Name of the task. Default "mardupsGATK"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param verbose [OPTIONAL] Enables progress messages. Default False.#
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] Hold job until job is finished. Job ID. 
#' @export


markdups_gatk=function(bin_path="tools/gatk/gatk",bam="",output_dir="",
verbose=FALSE,tmp_dir="",threads=3,ram=4,remove_duplicates=TRUE,
executor=make_unique_id("markdupsGATK"),task="markdupsGATK",mode="local",time="48:0:0",
update_time=60,wait=FALSE,hold=""){

    out_file_dir=set_dir(dir=output_dir,name="markdups_reports")

    tmp=""
    if (!tmp_dir==""){
    tmp=paste0("--tmp-dir ",tmp_dir)
    }

    dups=""
    if(remove_duplicates){
        dups="--remove-all-duplicates"
    }

    exec_code=paste0(bin_path," MarkDuplicatesSpark -I ",bam, " -O ",
    paste0(out_file_dir,"/",get_file_name(bam),".sorted.rmdup.",
    get_file_ext(bam))," -M ",paste0(out_file_dir,"/",get_file_name(bam),".gatk_rmdup.txt"),
    " ",tmp," --conf \'spark.executor.cores=",threads,"\'", dups)

   

    
    job=build_job(executor=executor,task=make_unique_id(task))
    if(mode=="batch"){
        out_file_dir2=set_dir(dir=out_file_dir,name="batch")
        exec_batch=build_job_exec(job=job,time=time,ram=ram,threads=threads,
        output_dir=out_file_dir2,hold=hold)
        exec_code=paste("echo 'source ~/.bashrc;",exec_code,"'|",exec_batch)
    }
    if(verbose){
        print(exec_code)
    }
    error=system(exec_code)
    if(error!=0){
        stop("markdups failed to run due to unknown error.
        Check std error for more information.")
    }

if(wait&&mode=="batch"){
    job_validator(job=job,
    time=update_time,verbose=verbose,threads=threads)
}



return(job)

}



#' Function for base quality recalibration
#'
#' This function recalibrates the base quality of the reads in two steps process based on GATK best practices guides.
#'
#' @param bin_path [REQUIRED] Path to picard executable. Default path tools/samtools/samtools.
#' @param bin_path2 [REQUIRED] Path to picard executable. Default path tools/gatk/gatk.
#' @param bin_path3 [REQUIRED] Path to picard executable. Default path tools/picard/build/libs/picard.jar.
#' @param bam [REQUIRED]  Path to BAM file.
#' @param ref_genome [REQUIRED]  Path to reference genome.
#' @param dbsnp [REQUIRED] Known variant database.Requires atleast 1.
#' @param threads [OPTIONAL] Number of threads to split the work.
#' @param ram [OPTIONAL] RAM memory per thread.
#' @param clean Clean input files. Default TRUE.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor Name of the executor. Default "recalGATK"
#' @param task Name of the task. Default "recalGATK"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param verbose [OPTIONAL] Enables progress messages. Default False.#
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export


recal_gatk=function(bin_path="tools/samtools/samtools",bin_path2="tools/gatk/gatk",
bin_path3="tools/picard/build/libs/picard.jar",bam="",ref_genome="",dbsnp="",ram=4,
threads=4,output_dir="",verbose=FALSE,executor=make_unique_id("recalGATK"),clean=TRUE,
task="recalGATK",mode="local",time="48:0:0",update_time=60,wait=FALSE,hold=""){

  out_file_dir=set_dir(dir=output_dir,name="recal_reports/recal_before")
  out_file_dir2=set_dir(dir=output_dir,name="recal_reports/recal_after")
  out_file_dir3=set_dir(dir=output_dir,name="recal_reports/recal_tmp")
  out_file_dir4=set_dir(dir=output_dir,name="recal_reports")


  job=parallel_generate_BQSR_gatk(bin_path=bin_path,bin_path2=bin_path2,bam=bam,
    ref_genome=ref_genome,dbsnp=dbsnp,
    output_dir=out_file_dir,
    verbose=verbose,executor=executor,mode=mode,threads=threads,ram=ram,clean=clean,
    time=time,update_time=update_time,wait=FALSE,hold=hold)

  job=parallel_apply_BQSR_gatk(bin_path=bin_path,bin_path2=bin_path2,bin_path3=bin_path3,
    bam=bam,ref_genome=ref_genome,rec_table=paste0(out_file_dir,"/",get_file_name(bam),".recal.table"),
    output_dir=out_file_dir4,clean=clean,verbose=verbose,executor=executor,threads=threads,mode=mode,ram=ram,time=time,
    update_time=update_time,wait=FALSE,hold=job)

  job=sort_and_index_bam_samtools(bin_path=bin_path,bam=paste0(out_file_dir4,"/",
    get_file_name(bam),".recal.",get_file_ext(bam)),output_dir=out_file_dir4,
    ram=ram,verbose=verbose,threads=threads,sort=FALSE,stats="",index=TRUE,clean=clean,
    mode=mode,executor=executor,ram=ram,time=time,update_time=update_time,wait=FALSE,hold=job)

  job=parallel_generate_BQSR_gatk(bin_path=bin_path,bin_path2=bin_path2,
    bam=paste0(out_file_dir4,"/",get_file_name(bam),".recal.",get_file_ext(bam)),
    ref_genome=ref_genome,dbsnp=dbsnp,threads=threads,clean=clean,output_dir=out_file_dir2,verbose=verbose,executor=executor,
    mode=mode,ram=ram,time=time,update_time=update_time,wait=FALSE,hold=job)

  job=recal_covariates_gatk(bin_path=bin_path2,before=paste0(out_file_dir,"/",get_file_name(bam),".recal.table"),
    after=paste0(out_file_dir2,"/",get_file_name(bam),".recal.table"),output_dir=out_file_dir4,executor=executor,
    mode=mode,threads=threads,ram=ram,time=time,update_time=update_time,wait=FALSE,hold=job)

  if(wait&&mode=="batch"){
    job_validator(job=job,time=update_time,verbose=verbose,threads=threads)
  }

  return(job)
}





#' Wrapper of BaseRecalibrator function from gatk
#'
#' Generates a recalibration table based on various covariates.
#' This function wraps around gatk BaseRecalibrator function.
#' For more information about this function: https://gatk.broadinstitute.org/hc/en-us/articles/360036898312-BaseRecalibrator
#'
#' @param bam [REQUIRED] Path to the BAM file.
#' @param bin_path [REQUIRED] Path to gatk executable. Default tools/gatk/gatk.
#' @param ref_genome [REQUIRED] Path to reference genome
#' @param dbsnp [REQUIRED] Path to known snp positions in VCF format. Multiple vcf can be supplied as a vector.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor Name of the executor. Default "generateBQSR"
#' @param task Name of the task. Default "generateBQSR"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param threads [OPTIONAL] Number of threads to split the work.
#' @param ram [OPTIONAL] If batch mode. RAM memory in GB per job. Default 1
#' @param update_time [OPTIONAL] If batch mode. Show job updates every update time. Default 60
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export

generate_BQSR_gatk=function(region="",bin_path="tools/gatk/gatk",bam="",ref_genome="",
dbsnp="",output_dir="",verbose=FALSE, threads=4,ram=4,mode="local",
executor=make_unique_id("generateBQSR"),task="generateBQSR",
time="48:0:0",update_time=60,wait=FALSE,hold=""){

  out_file_dir=set_dir(dir=output_dir)

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

  exec_code=paste0(bin_path," BaseRecalibrator -I ",bam, " -R ", ref_genome,dbsnp,
  reg," -O ",out_file)


  job=build_job(executor=executor,task=make_unique_id(task))
  if(mode=="batch"){
       out_file_dir2=set_dir(dir=out_file_dir,name="batch")
       exec_batch=build_job_exec(job=job,time=time,ram=ram,threads=threads,
       output_dir=out_file_dir2,hold=hold)
       exec_code=paste("echo 'source ~/.bashrc;",exec_code,"'|",exec_batch)
  }

  if(verbose){
    print(exec_code)
  }
  error=system(exec_code)
  if(error!=0){
    stop("gatk failed to run due to unknown error.
    Check std error for more information.")
  }

  if(wait&&mode=="batch"){
      job_validator(job=job,
      time=update_time,verbose=verbose,threads=threads)
  }

  return(job)
}



#' Multiregion parallelization of generate_BQSR function
#'
#' Generates a recalibration table based on various covariates.
#' This function wraps around gatk BaseRecalibrator function.
#' For more information about this function: https://gatk.broadinstitute.org/hc/en-us/articles/360036898312-BaseRecalibrator
#'
#' @param bam [REQUIRED] Path to the BAM file.
#' @param bin_path [REQUIRED] Path to gatk executable. Default tools/samtools/samtools.
#' @param bin_path2 [REQUIRED] Path to gatk executable. Default tools/gatk/gatk.
#' @param ref_genome [REQUIRED] Path to reference genome
#' @param dbsnp [REQUIRED] Path to known snp positions in VCF format. Multiple vcf can be supplied as a vector.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor Name of the executor. Default "generateBQSR"
#' @param task Name of the task. Default "generateBQSR"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param threads Number of threads to split the work. Default 3
#' @param ram [OPTIONAL] If batch mode. RAM memory in GB per job. Default 1
#' @param update_time [OPTIONAL] If batch mode. Show job updates every update time. Default 60
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export
#' @import pbapply



parallel_generate_BQSR_gatk=function(bin_path="tools/samtools/samtools",
bin_path2="tools/gatk/gatk",bam="",ref_genome="",dbsnp="",threads=3,ram=4,
executor=make_unique_id("par_generateBQSR"),task="par_generateBQSR",output_dir="",
verbose=FALSE,mode="local",time="48:0:0",update_time=60,wait=FALSE,hold=""){

  options(scipen = 999)
  options(warn = -1)

  out_file_dir=set_dir(dir=output_dir)

  dat=get_bam_reference_chr(bin_path=bin_path,bam=bam,verbose=verbose)
  dat$start=dat$start+1
  dat=dat %>% dplyr::mutate(region=paste0(chr,":",start,"-",end),
  report_loc=paste0(out_file_dir,get_file_name(bam),".",paste0(chr,":",start,"-",end),".recal.table"))


  jobs=parallel::mclapply(1:nrow(dat),FUN=function(x){
    tmp=dat[x,]
    id=generate_BQSR_gatk(region=tmp$region,
    bin_path=bin_path2,bam=bam,ref_genome=ref_genome,
    dbsnp=dbsnp,output_dir=out_file_dir,verbose=verbose,
    executor=executor,task="generateBQSR",mode=mode,time=time,
    threads=threads,ram=ram,update_time=update_time,wait=FALSE,hold=hold)   
  },mc.cores=ifelse(mode=="local",threads,3))
  

  job=gather_BQSR_reports_gatk(bin_path=bin_path2,report=dat$report_loc,
  executor=executor,task="gatherBQSR",output_dir=out_file_dir,
  output_name=get_file_name(bam),verbose=verbose,mode=mode,time=time,
  threads=threads,ram=ram,update_time=update_time,wait=FALSE,hold=jobs)
  
  if(wait&&mode=="batch"){
    job_validator(job=job,time=update_time,verbose=verbose,threads=threads)
  }

  return(job)

}



#' Wrapper around gatk GatherBQSRReports function
#'
#' This functions collects the Recalibration reports generated from scattered parallel_generate_BQSR output
#' This function wraps around gatk GatherBQSRReports function.
#' For more information about this function: https://gatk.broadinstitute.org/hc/en-us/articles/360036829851-GatherBQSRReports
#'
#' @param bin_path [REQUIRED] Path to gatk executable. Default tools/gatk/gatk.
#' @param reports_dir [REQUIRED] Path to the directory where reports are stored.
#' @param report [REQUIRED] Path to indivual reports.
#' @param output_name [OPTIONAL] Name for the output report file.
#' @param clean Clean input files. Default TRUE.
#' @param executor [OPTIONAL] Task executor name. Default "gatherBQSR"
#' @param task [OPTIONAL] Task name. Default "gatherBQSR"
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param threads Number of threads to split the work. Default 3
#' @param ram [OPTIONAL] If batch mode. RAM memory in GB per job. Default 1
#' @param update_time [OPTIONAL] If batch mode. Show job updates every update time. Default 60
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export

gather_BQSR_reports_gatk=function(bin_path="tools/gatk/gatk",report="",reports_dir="",
executor=make_unique_id("gatherBQSR"),task="gatherBQSR",output_name="Report", clean=TRUE,
output_dir="",verbose=FALSE,mode="local",time="48:0:0",threads=4,ram=4,update_time=60,wait=FALSE,hold=""){

  out_file_dir=set_dir(dir=output_dir)
  if(report==""){
      files=list.files(reports_dir,full.names=TRUE,pattern=":")
      files=files[grepl(".recal.table",files)]
  }else{
    files=report
  }



  exec_code=paste0(bin_path," GatherBQSRReports ",paste(" -I ",files,collapse=" "),
    " -O ",paste0(out_file_dir,output_name,".recal.table"))


  if(clean){
    exec_code=paste(exec_code," && rm",paste(files,collapse=" "))
  }

  job=build_job(executor = executor,task=make_unique_id(task))
  if(mode=="batch"){
    
      out_file_dir2=set_dir(dir=out_file_dir,name="batch")
       exec_batch=build_job_exec(job=job,
       time=time,ram=ram,threads=threads,
       output_dir=out_file_dir2,hold=hold)
       exec_code=paste("echo 'source ~/.bashrc;",exec_code,"'|",exec_batch)
  }

  if(verbose){
    print(exec_code)
  }

  error=system(exec_code)
  if(error!=0){
    stop("gatk failed to run due to unknown error.
    Check std error for more information.")
  }

  if(wait&&mode=="batch"){
    job_validator(job=job,time=update_time,verbose=verbose,
    threads=threads)
  }
  return(job)
}



#' Wrapper of applyBQSR function gatk
#'
#' Applies numerical corrections to each individual basecall based on the covariates analyzed before.
#' This function wraps around gatk applyBQSR  function.
#' For more information about this function: https://gatk.broadinstitute.org/hc/en-us/articles/360050814312-ApplyBQSR
#'
#' @param bam [REQUIRED] Path to the BAM file.
#' @param bin_path [REQUIRED] Path to gatk executable. Default tools/gatk/gatk.
#' @param ref_genome [REQUIRED] Path to reference genome
#' @param rec_table [REQUIRED] Path to covariates table generated by generate_BSQR.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param threads Number of threads to split the work. Default 3
#' @param ram [OPTIONAL] If batch mode. RAM memory in GB per job. Default 1
#' @param executor [OPTIONAL] Task executor name. Default "applyBQSR"
#' @param task [OPTIONAL] Task nam. Default "applyBQSR"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param ram [OPTIONAL] If batch mode. RAM memory in GB per job. Default 1
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export


apply_BQSR_gatk=function(region="",bin_path="tools/gatk/gatk",bam="",ref_genome="",
rec_table="",output_dir="",verbose=FALSE,mode="local", threads=4,ram=4,
executor=make_unique_id("applyBQSR"),task="applyBQSR",time="48:0:0",
update_time=60,wait=TRUE,hold=""){

  out_file_dir=set_dir(dir=output_dir)
 
  reg=""
  if (region==""){
      out_file=paste0(" ", out_file_dir,"/", get_file_name(bam),".recal.",get_file_ext(bam))
  }else{
      reg=paste0(" -L ",strsplit(region,"_")[[1]][2], " ")
      out_file=paste0(out_file_dir,"/", get_file_name(bam),".",region,".recal.",get_file_ext(bam))
  }
  exec_code=paste(bin_path," ApplyBQSR -I ",bam, " -R ", ref_genome,
   " --bqsr-recal-file ",rec_table, " -O ",out_file,reg)
   
  job=build_job(executor=executor,task=make_unique_id(task))
  if(mode=="batch"){
       out_file_dir2=set_dir(dir=out_file_dir,name="batch")
       exec_batch=build_job_exec(job=job,time=time,ram=ram,threads=threads,
       executor=executor,task=task,
       output_dir=out_file_dir2,hold=hold)
       exec_code=paste("echo 'source ~/.bashrc;",exec_code,"'|",exec_batch)
  }

  if(verbose){
    print(exec_code)
  }


  error=system(exec_code)
  if(error!=0){
    stop("gatk failed to run due to unknown error.
    Check std error for more information.")
  }

  if(wait&&mode=="batch"){
      job_validator(job=job,
      time=update_time,verbose=verbose,threads=threads)
  }

  return(job)
}

#' Multiregion parallelization of apply_BQSR function
#'
#' Recalibrates
#' Applies numerical corrections to each individual basecall based on the covariates analyzed before.
#' For more information about this function: https://gatk.broadinstitute.org/hc/en-us/articles/360050814312-ApplyBQSR
#'
#' @param bam [REQUIRED] Path to the BAM file.
#' @param bin_path [REQUIRED] Path to samtools executable. Default tools/samtools/samtools.
#' @param bin_path2 [REQUIRED] Path to gatk executable. Default tools/gatk/gatk.
#' @param bin_path3 [REQUIRED] Path to picard executable. Default tools/picard/build/libs/picard.jar
#' @param ref_genome [REQUIRED] Path to reference genome
#' @param rec_table [REQUIRED] Path to the recalibratio table.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor [OPTIONAL] Task executor name. Default "par_applyBQSR"
#' @param task [OPTIONAL] Task name. Default "par_applyBQSR"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param threads [OPTIONAL] Number of threads for the main job. Default 4
#' @param ram [OPTIONAL] If batch mode. RAM memory in GB per job. Default 1
#' @param update_time [OPTIONAL] If batch mode. Show job updates every update time. Default 60
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export
#' @import pbapply

parallel_apply_BQSR_gatk=function(bin_path="tools/samtools/samtools",bin_path2="tools/gatk/gatk",
bin_path3="tools/picard/build/libs/picard.jar",bam="",ref_genome="",rec_table="",
output_dir="",verbose=FALSE,mode="local",executor=make_unique("par_applyBQSR"),task="par_applyBQSR",
time="48:0:0",threads=4,ram=4,update_time=60,wait=FALSE, hold=""){

  options(scipen = 999)
  options(warn = -1)

  out_file_dir=set_dir(dir=output_dir)
  dat=get_bam_reference_chr(bin_path=bin_path,bam=bam,verbose=verbose)
  dat$start=dat$start+1
  dat$pos=1:nrow(dat)
  dat=dat %>% dplyr::mutate(region=paste0(pos,"_",chr,":",start,"-",end),
    bam_loc=paste0(out_file_dir,"/", get_file_name(bam),".",paste0(pos,"_",chr,":",start,"-",end),".recal.",get_file_ext(bam)))


  jobs=parallel::mclapply(1:nrow(dat),FUN=function(x){
    tmp=dat[x,]
    apply_BQSR_gatk(region=tmp$region,bin_path=bin_path2,bam=bam,ref_genome=ref_genome,
    executor=executor,task="applyBQSR",rec_table=rec_table,
    output_dir=out_file_dir,verbose=verbose,mode=mode,time=time,threads=threads,
    ram=ram,update_time=update_time,hold=hold,wait=FALSE)},mc.cores=ifelse(mode=="local",threads,3))
   

  job=gather_bam_files(bin_path=bin_path3,bam=dat$bam_loc,output_dir=out_file_dir,
  output_name=paste0(get_file_name(bam),".recal.sorted.rmdup.sorted"),executor=executor,
  task="GatherBAM",hold=jobs,mode=time,time=time,threads=threads,ram=ram,
  update_time=update_time,wait=FALSE)

  if(wait&&mode=="batch"){
    job_validator(job=job_name,time=update_time,verbose=verbose,threads=threads)
  } 
  return(job)
}





#' Wrapper around gatk GatherBamFiles function
#'
#' This functions collects the Recalibration reports generated from scattered parallel_apply_BQSR output
#' This function wraps around gatk GatherBamFiles function.
#' For more information about this function: https://gatk.broadinstitute.org/hc/en-us/articles/360037055512-GatherBamFiles-Picard-
#'
#' @param bin_path [REQUIRED] Path to gatk executable. Default tools/picard/build/libs/picard.jar.
#' @param bamr [REQUIRED] Path to BAM file/s.
#' @param bams_dir [REQUIRED] Path to the directory where BAM files are stored.
#' @param output_name [OPTIONAL] Name for the output file name.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param clean Clean input files. Default TRUE.
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param threads [OPTIONAL] Number of threads to split the work. Default 4
#' @param ram [OPTIONAL] RAM memory to asing to each thread. Default 4
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor [OPTIONAL] Task executor name. Default "gatherBAM"
#' @param task [OPTIONAL] Task name. Default "gatherBAM"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Show job updates every update time. Default 60
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export

gather_bam_files=function(bin_path="tools/picard/build/libs/picard.jar",bam="",bams_dir="",
output_name="File",output_dir="",verbose=FALSE,threads=4, ram=4, mode="local",clean=FALSE,
executor=make_unique_id("gatherBAM"),task="gatherBAM",time="48:0:0",update_time=60,
wait=FALSE,hold=""){

  out_file_dir=set_dir(dir=output_dir)
  if(bam==""){
    files=list.files(bams_dir,full.names=TRUE,pattern=":")
    files=files[grepl("bam$",files)]
  }else{
    files=bam
  }

  files=files[order(as.numeric(lapply(lapply(lapply(lapply(lapply(lapply(basename(files),
    FUN=strsplit,split="\\."),FUN="[[",index=1),FUN="[",index=2),
    FUN=strsplit,split="_"),FUN="[[",index=1),FUN="[",index=1)))]

  exec_code=paste0("java -jar ",bin_path," GatherBamFiles ",
    paste0(" I=",files,collapse=" ")," O=",paste0(out_file_dir,"/",output_name,".bam"))

  if(clean){
    exec_code=paste(exec_code," && rm",paste(files,collapse=" "))
  }

 
  job=build_job(executor=executor,task=make_unique_id(task))

  if(mode=="batch"){
       out_file_dir2=set_dir(dir=out_file_dir,name="batch")
       exec_batch=build_job_exec(job=job,hold=hold,time=time,ram=ram,
       threads=threads,output_dir=out_file_dir2)
       exec_code=paste("echo 'source ~/.bashrc;",exec_code,"'|",exec_batch)
  }

  if(verbose){
    print(exec_code)
  }

  error=system(exec_code)

  if(error!=0){
    stop("picard failed to run due to unknown error.
    Check std error for more information.")
  }

  if(wait&&mode=="batch"){
    job_validator(job=job,time=update_time,
    verbose=verbose,threads=threads)
  }
  return(job)
}

#' Wrapper of AnalyzeCovariates function in gatk
#'
#' Generates a report of the recalibrated values.
#' This function wraps around gatk AnalyzeCovariates function.
#' For more information about this function: https://gatk.broadinstitute.org/hc/en-us/articles/360037066912-AnalyzeCovariates
#'
#' @param bin_path [REQUIRED] Path to gatk executable. Default tools/gatk/gatk.
#' @param before [REQUIRED] Recalibration table produced by generate_BQSR function.
#' @param after [OPTIONAL] Recalibration table produced by generate_BQSR function.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param threads [OPTIONAL] Number of threads to split the work. Default 4
#' @param ram [OPTIONAL] RAM memory to asing to each thread. Default 4
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor [OPTIONAL] Task executor name. Default "recalCovariates"
#' @param task [OPTIONAL] Task name. Default "recalCovariates"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export


recal_covariates_gatk=function(bin_path="tools/gatk/gatk",before="",after="",
output_dir="",verbose=FALSE,threads=4,ram=4,mode="local",
executor=make_unique_id("recalCovariates"),task="recalCovariates",time="48:0:0",
update_time=60,wait=FALSE,hold=""){

  out_file_dir=set_dir(name=output_dir)
 
  if (before!="" & after==""){
    exec_code=paste0(bin_path," AnalyzeCovariates -bqsr ",before," -plots ",
    paste0(out_file_dir,"/",get_file_name(before),"_covariates_analysis_before.pdf"))
   
  }else{
    exec_code=paste0(bin_path," AnalyzeCovariates -before ",before," -after ",after,
      " -plots ",paste0(out_file_dir,"/",get_file_name(before),"_covariates_analysis.pdf"))

  }

  out_file_dir2=set_dir(dir=out_file_dir,name="batch")
  job=build_job(executor=executor,task=make_unique_id(task))
  if(mode=="batch"){
       exec_batch=build_job_exec(job=job,hold=hold,time=time,ram=ram,threads=threads,
       output_dir=out_file_dir2)
       exec_code=paste("echo 'source ~/.bashrc;",exec_code,"'|",exec_batch)
  }

  if(verbose){
      print(exec_code)
  }

  error=system(exec_code)
  
  if(error!=0){
    stop("gatk failed to run due to unknown error.
    Check std error for more information.")
  }

  if(wait&&mode=="batch"){
    job_validator(job=job,time=update_time,
    verbose=verbose,threads=threads)
  }

  return(job)
}





