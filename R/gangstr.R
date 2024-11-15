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

call_gangstr=function(
  bin_gangstr=build_default_tool_binary_list()$bin_gangstr,
  bam=NULL,
  regions=build_default_reference_list()$HG19$reference$chrx_tr,
  ref_genome=build_default_reference_list()$HG19$reference$genome,
  chromosomes=c(1:22,"X","Y"),
  ...
){  

      run_main=function(
              .env
          ){
              .this.env=environment()
              append_env(to=.this.env,from=.env)
        
              set_main(.env=.this.env)

              add=""
              if(!is.null(chromosomes)){
                add=paste(" --chrom ",paste0(chromosomes,collpase=","))
              }
              
              .main$exec_code=paste0(
                    bin_gangstr, 
                    " --bam ",bam,
                    " --ref ", ref_genome,
                    " --regions ", regions,
                    " --out ", paste0(out_file_dir,"/",ULPwgs::get_file_name(bam)),
                    " --bam-samps ", ULPwgs::get_file_name(bam),
                    add
             )

            run_job(
                .env=.this.env
            )
            
            .env$.main<-.main

        } 


            .base.env=environment()
            list2env(list(...),envir=.base.env)
            set_env_vars(
                .env= .base.env,
                vars=list("bam")
            )

            launch(.env=.base.env)  

}
