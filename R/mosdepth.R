#' Extract Structural Variants from BAM File Using LUMPY
#'
#' This function calls structural variants (SVs) from a BAM file using LUMPY Express. It generates a VCF file with SV calls, as well as BAM files containing split and discordant reads used for SV detection. Output files are written to disk and tracked in the environment.
#'
#' @param env_lumpy Path to the LUMPY conda environment. Default: from build_default_python_enviroment_list().
#' @param bam Path to the input BAM file. (Required)
#' @param ... Additional arguments passed to environment setup and job execution.
#'
#' @return No direct return value. Output VCF and BAM files are written to disk and tracked in the environment.
#' @export
#' 
coverage_wgs_mosdepth=function(
  bin_mosdepth=build_default_python_enviroment_list()$bin_mosdepth,
  bam=NULL,
  ...
  ){
     run_main=function(
    .env
  ){
    .this.env=environment()
    append_env(to=.this.env,from=.env)
    set_main(.env=.this.env)

    .main$out_files$mosdepth$dist=paste0(out_file_dir,"/",input_id,".mosdepth.global.dist.txt")
    .main$out_files$mosdepth$summary=paste0(out_file_dir,"/",input_id,".mosdepth.summary.txt")
    .main$exec_code=paste0(
        bin_mosdepth," -n --fast-mode -t ",threads," ", 
        out_file_dir,"/",input_id," ", input
    )
    
    run_job(.env=.this.env)
    .env$.main <- .main
  }

   .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="bam"
    )
    launch(.env=.base.env)
}
