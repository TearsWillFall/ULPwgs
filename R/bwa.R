
#' Read alignment using BWA
#'
#' This function aligns a sequence of reads to a reference genome
#' using a Burrows-Wheeler aligner. The reference genome has to be previously
#' been sorted and indexed. It generates a BAM file with the aligned sequence.
#'
#' @param file_R1 Path to the input file with the sequence.
#' @param file_R2 [Optional] Path to the input with the reverse read sequence.
#' @param bin_bwa Path to bwa executable. Default path tools/bwa/bwa.
#' @param bin_samtools Path to samtools executable. Default path tools/samtools/samtools.
#' @param ref_genome Path to input file with the reference genome sequence.
#' @param id_tag Read group identifier. Default NA.
#' @param pu_tag Platform unit identifier. Default NA.
#' @param sm_tag Sample name tag. If not given will use sample name as default.
#' @param pl_tag Platform/technology used to produce the read. Default "ILLUMINA". Options ["ILLUMINA","SOLID", "LS454", "HELICOS","PACBIO"]
#' @param lb_tag DNA preparation library identifier. 
#' @param sort Sort aligned file. Default TRUE.
#' @param coord_sort Sort BAM file by coordinate. Alternatively sort by name. Default TRUE.
#' @param index Generate index file if BAM sorted by coordinate. Default TRUE.
#' @param clean Clean intermediary files. Default TRUE.
#' @param ram RAM memory to use for sorting and indexing. Provided in GB.
#' @param threads Number of CPU cores to use. Default 3.
#' @param stats Generate BAM stats. Default all. Options ["all","flag","index",""]
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Task EXECUTOR ID. Default "alignment"
#' @param task_nam Task nam. Default "alignment"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @export

alignment_bwa=function(
  bin_bwa=build_default_binary_list()$alignment$bin_bwa,
  bin_samtools=build_default_binary_list()$alignment$bin_samtools,
  file_R1="",file_R2="",threads=3,ram=4,id_tag="NA",
  pu_tag="NA",pl_tag="ILLUMINA",lb_tag="NA",
  sm_tag="",sort=TRUE,coord_sort=TRUE,index=TRUE,
  clean=TRUE,stats="all",ref_genome="",output_dir=".",
  verbose=FALSE,batch_config=build_default_preprocess_config(),
  executor_id=make_unique_id("alignment"),
  task_name="alignment",mode="local",time="48:0:0",
  update_time=60,wait=FALSE,hold=NULL
){
    
    argg <- as.list(environment())

    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir,name="bwa_reports")
   
    input_files=file_R1
    if (!check_missing(file_R2)){
      sm_tag=ifelse(sm_tag=="",intersect_file_name(file_R1,file_R2),sm_tag)
      input_files=paste(file_R1,file_R2)
    }else{
      sm_tag=ifelse(sm_tag=="",get_file_name(file_R1),sm_tag)
    }
    out_file=paste0(out_file_dir,"/",sm_tag,".bam")
    GPU=paste0("\"@RG\\tID:",id_tag,"\\tPL:",pl_tag,"\\tPU:",pu_tag,"\\tLB:",
    lb_tag,"\\tSM:",sm_tag,"\"")

    exec_code=paste(bin_bwa,"mem -t", threads," -v 2 -R",GPU,"-M",ref_genome,
        input_files, "| ",paste0(bin_samtools)," view -h -b >",out_file)
    
    job=build_job(executor_id=executor_id,task_id=task_id)

    if(mode=="batch"){
      out_file_dir2=set_dir(dir=out_file_dir,name="batch")
      batch_code=build_job_exec(job=job,time=time,ram=ram,
      threads=threads,output_dir=out_file_dir2,hold=hold)
      exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
    }

    if(verbose){
       print_verbose(job=job,arg=argg,exec_code=exec_code)
    }

    error=execute_job(exec_code=exec_code)
    if(error!=0){
      stop("bwa failed to run due to unknown error.
      Check std error for more information.")
    }

    job_report=build_job_report(
      job_id=job,
      executor_id=executor_id, 
      task_id=task_id,
      input_args=argg,
      out_file_dir=out_file_dir,
      out_files=list(
        bam=out_file)
      )
    
    if(sort){
      job_report[["steps"]][["sort_and_index"]]=sort_and_index_bam_samtools(
        bin_samtools=bin_samtools,
        bam=out_file,output_dir=out_file_dir,batch_config=batch_config,
        executor_id=task_id,ram=ram,verbose=verbose,threads=threads,
        coord_sort=coord_sort,index=index,stats=stats,clean=clean,mode=mode,time=time,
        update_time=update_time,wait=FALSE,hold=job_report$job_id)
    
    }

    if(wait&&mode=="batch"){
        job_validator(job=c(
          job_report[["steps"]][["sort_and_index"]][["steps"]][["sort"]]$job_id,
          job_report[["steps"]][["sort_and_index"]][["steps"]][["index"]]$job_id,
          job_report[["steps"]][["sort_and_index"]][["steps"]][["stats"]][["steps"]][["index"]]$job_id,
          job_report[["steps"]][["sort_and_index"]][["steps"]][["stats"]][["steps"]][["flag"]]$job_id),
        time=update_time,verbose=verbose,threads=threads)
    }

    return(job_report)

  }

#' Index reference genome
#'
#' This function indexes a reference genome using the bwa aligner
#'
#' @param file Path to the input file with the reference genome in FASTA format.
#' @param bin_bwa Path to bwa executable. Default path tools/bwa/bwa.
#' @param executor_id Task EXECUTOR ID. Default "refIndex"
#' @param task_name Task name. Default "refIndex"
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @export


index_ref_bwa=function(
  bin_bwa=build_default_tool_binary_list()$bin_bwa,
  file="",threads=4,ram=4,verbose=FALSE,
  batch_config=build_default_preprocess_config(),
  executor_id=make_unique_id("refIndex"),
  task_name="refIndex",mode="local",time="48:0:0",
  update_time=60,wait=FALSE,hold=NULL
 ){

  argg <- as.list(environment())

  task_id=make_unique_id(task_name)
  exec_code=paste(bin_bwa,"index", file)
  
  job=build_job(executor_id=executor_id,task_id=task_id)
  if(mode=="batch"){

    out_file_dir2=set_dir(dir=".",name="batch")
    batch_code=build_job_exec(job=job,time=time,ram=ram,
    threads=threads,output_dir=out_file_dir2,hold=hold)
    exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
  }

  if(verbose){
       print_verbose(job=job,arg=argg,exec_code=exec_code)
    }

  error=execute_job(exec_code=exec_code)
  if(error!=0){
    stop("bwa failed to run due to unknown error.
    Check std error for more information.")
  }
 
  job_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    task_id=task_id,
    input_args=argg,
    exec_code=exec_code, 
    out_file_dir=out_file_dir,
    out_files=list(
      index=paste0(file,".fai")
      )
    )
  
  if(wait&&mode=="batch"){
      job_validator(job=job_report$job_id,
      time=update_time,verbose=verbose,threads=threads)
  }

  return(job_report)

}






#' Align Paired-End FASTQ Files using BWA
#'
#' This function aligns paired-end FASTQ files to a reference genome using the BWA mem algorithm.
#' It generates a BAM file and supports read group tagging for sample/library identification.
#'
#' The alignment is performed with the following options:
#' \itemize{
#'   \item Algorithm: BWA mem (suitable for reads > 70bp)
#'   \item Mark shorter split hits as secondary (-M flag)
#'   \item Verbose level: 2
#'   \item Pipe output directly to samtools view for BAM conversion
#' }
#'
#' Read group tags follow SAM format specification and are essential for multi-sample analysis.
#'
#' For more information on BWA:
#' https://github.com/lh3/bwa
#'
#' @param bin_bwa [REQUIRED] Path to BWA executable. Default: from build_default_binary_list().
#' @param bin_samtools [REQUIRED] Path to samtools executable. Default: from build_default_binary_list().
#' @param ref_genome [REQUIRED] Path to indexed reference genome. Default: HG19 from build_default_reference_list().
#' @param fastq [REQUIRED] Path to paired-end FASTQ files or data structure containing R1 and R2 file paths.
#' @param tags [OPTIONAL] Read group tags as a list with elements:
#'   \describe{
#'     \item{id_tag}{(ID) Read group identifier}
#'     \item{pu_tag}{(PU) Platform unit (flowcell.lane.barcode)}
#'     \item{pl_tag}{(PL) Platform (default: ILLUMINA)}
#'     \item{lb_tag}{(LB) Library name}
#'     \item{sm_tag}{(SM) Sample name}
#'   }
#' @param soft_clipping_supplementary [OPTIONAL] Enable soft-clipping of supplementary alignments 
#'   (-Y flag). Default: FALSE.
#' @param ... Additional arguments passed to internal functions for environment setup and job execution.
#'
#' @return Aligned BAM file. Output path tracked in job report.
#'
#' @seealso
#'   \link{sort_and_index_bam_samtools} for sorting and indexing the resulting BAM file
#'
#' @export

new_alignment_bwa=function(
  bin_bwa=build_default_binary_list()$alignment$bin_bwa,
  bin_samtools=build_default_binary_list()$alignment$bin_samtool,
  ref_genome=build_default_reference_list()$HG19$reference$genome,
  fastq=NULL,
  tags=list(
    id_tag="NA",
    pu_tag="NA",
    pl_tag="ILLUMINA",
    lb_tag="NA",
    sm_tag="NA"),
  soft_clipping_supplementary=FALSE,
  ...
){

  
  run_main=function(
    .env
  ){
    .this.env=environment()
    append_env(to=.this.env,from=.env)
    set_main(.env=.this.env)

    .main$out_files$bam=paste0(out_file_dir,"/",input_id,".bam")

    if(!is.null(tags)){
        tag_annot=paste0(
          " -R \"@RG\\tID:",tags$id_tag,
          "\\tPL:",tags$pl_tag,
          "\\tPU:",tags$pu_tag,
          "\\tLB:",tags$lb_tag,
          "\\tSM:",tags$sm_tag,"\""
          )
    }

    .main$exec_code=paste(
      bin_bwa," mem -t ",threads,
      " -v 2 ", ifelse(!is.null(tags),tags_annot,""),
      ifelse(soft_clipping_supplementary," -Y ",""),
      " -M ",ref_genome,
      input$fastq_r1,input$fastq_r2, " | ",
      bin_samtools, " view -bh >",  
      .main$out_files$bam
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

