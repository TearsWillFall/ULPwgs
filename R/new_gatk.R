
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


new_markdups_gatk=function(
  sif_gatk=build_default_sif_list()$sif_gatk,
  bam=NULL,
  remove_duplicates=TRUE,
  ...
  ){


      run_main=function(
              .env
          ){
              .this.env=environment()
              append_env(to=.this.env,from=.env)
              

              set_main(.env=.this.env)
              
           
            
                add=""
                if(remove_duplicates){
                    add=" --remove-all-duplicates"
                }

                .main$out_files$report=paste0(out_file_dir,"/",input_id,".gatk_rmdup.txt")
                
                .main$exec_code=paste0(out_file_dir,"/",input_id,".gatk_rmdup.txt")
                    exec_code=paste0(" singularity exec -H ",getwd(),":/home ",
                    sif_gatk," /gatk/gatk MarkDuplicatesSpark -I ",bam, " -O ",paste0(out_file_dir,"/",input_id),
                    " -M ",.main$out_files$report," ", paste0(" --tmp-dir ",tmp_dir),
                    " --conf \'spark.executor.cores=",threads,"\'", add
                )
              run_job(
                .env=.this.env
              )

              .env$.main<-.main
          }

      .base.env=environment()
      list2env(list(...),envir=.base.env)
      set_env_vars(
          .env=.base.env,
          vars="bam"
      )

      launch(.env=.base.env)
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


new_recal_gatk=function(
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

new_generate_BQSR_gatk=function(
  region=NULL,
  sif_gatk=build_default_sif_list()$sif_gatk,
  bam=NULL,
  ref_genome=build_default_reference_list()$HG19$reference$genome,
  dbsnp=build_default_reference_list()$HG19$database$all_common,
  ...
){  

 
      run_main=function(
              .env
          ){
              .this.env=environment()
              append_env(to=.this.env,from=.env)
              

              set_main(.env=.this.env)
               
            if (is.null(region)){
                .main$out_files$recal_table=paste0(out_file_dir,"/",input_id,".recal.table")
            }else{
                if(is.list(region)){
                    region=input
                }else{
                    bam=input
                }
                reg=paste0(" -L ",region)
                .main$out_files$recal_table=paste0(out_file_dir,"/",get_file_name(bam),".",region,".recal.table")
            }

            ## Multiple vcf with snps can be given
            
            add=""
            if (!is.null(dbsnp)){
                add=paste(" --known-sites ",dbsnp,collapse=" ")
            }

            .main$exec_code=paste0(
                "singularity exec -H ",
                getwd(),":/home ",sif_gatk,
                " /gatk/gatk BaseRecalibrator -I ",bam,
                " -R ", ref_genome, dbsnp,
                reg," -O ",out_file,paste0(" --tmp-dir ",tmp_dir)
            )

            run_job(
                .env=.this.env
            )

            .env$.main<-.main

            .base.env=environment()
            list2env(list(...),envir=.base.env)
            set_env_vars(
                .env= .base.env,
                vars=list("bam","region")
            )

            launch(.env=.base.env)
    }   

}

