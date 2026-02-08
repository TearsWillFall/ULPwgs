
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



