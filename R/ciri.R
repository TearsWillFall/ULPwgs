#' Generate a quality control (QC) report from a fastaqc file
#'
#' This function takes a set of sequence files (fastq,SAM,BAM...) and
#' returns a report in HTML format.
#'
#'
#' @param file_R1 Path to the input file with the sequence.
#' @param file_R2 [Optional] Path to the input with the reverse read sequence.
#' @param bin_path Path to fastQC executable. Default path tools/CIRI/CIRI_Full_v2.1.1.jar.
#' @param threads Number of CPU cores to use. Default 3.
#' @param ref_genome Path to reference genome
#' @param db_annot Path to annotation GTF file
#' @param output_name File output name
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Executor ID. Default "fastQC"
#' @param task_name Name of the task. Default "fastQC"
#' @param ram RAM memory for batched job. Default 4
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param hold Job to hold on in batched mode.
#' @export


circ_rna=function(bin_path="tools/CIRI/CIRI_Full_v2.1.1.jar",file_R1="",file_R2="",
output_dir="",ref_genome="database/Homo_sapiens_GRCh38.fa",
db_annot="database/human_gencode_vch38.gtf",output_name="",verbose=FALSE,executor_id=make_unique_id("circRNA"),
task_name="circRNA",mode="local",threads=3,ram=4,time="48:0:0",update_time=60,wait=FALSE,hold=""){

  
  argg <- as.list(environment())

  task_id=make_unique_id(task_name)

  out_file_dir=set_dir(dir=output_dir,name="ciri_reports")


  exec_code=paste("java -jar ",bin_path,"Pipeline -1",file_R1, " -2",file_R2, " -d ",
  out_file_dir,"-t ",threads,"-a ",db_annot,"-r",ref_genome,"-o",output_name)


  job=build_job(executor_id=executor_id,task_id=task_id)

  if(mode=="batch"){
    out_file_dir2=set_dir(dir=out_file_dir,name="batch")
    batch_code=build_job_exec(job=job,time=time,ram=ram,threads=threads,
    output_dir=out_file_dir2,hold=hold)
    exec_code=paste0("echo 'source ~/.bashrc;",exec_code,"'|",batch_code)
  }

  if(verbose){
      print_verbose(exec_code=exec_code)
  }

  error=system(exec_code)
  if(error!=0){
    stop("circa failed to run due to unknown error.
    Check std error for more information.")
  }


  job_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    task_id=task_id,
    input_args = argg,
    out_file_dir=out_file_dir,
    out_files=list(
      file=NA)
    )

  if(wait&&mode=="batch"){
    job_validator(job=job_report$job_id,
    time=update_time,verbose=verbose,threads=threads)
  }

  return(job_report)
  
  }
