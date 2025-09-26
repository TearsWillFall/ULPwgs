
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
extract_svs_lumpy=function(
  env_lumpy=build_default_python_enviroment_list()$env_lumpy,
  bam=NULL,
  ...
  ){
     run_main=function(
    .env
  ){
    .this.env=environment()
    append_env(to=.this.env,from=.env)
    set_main(.env=.this.env)

    .main$out_files$lumpy$vcf=paste0(out_file_dir,"/",input_id,".tmp/",input_id,".lumpy.SV.vcf")
    .main$out_files$lumpy$split_bam=paste0(out_file_dir,"/",input_id,".tmp/",input_id,".lumpy.split.bam")
    .main$out_files$lumpy$disc_bam=paste0(out_file_dir,"/",input_id,".tmp/",input_id,".lumpy.disc.bam")
    .main$exec_code=paste0(
      "conda activate ",env_lumpy,"; lumpyexpress -B ",input,
      " -T ",input_id,".tmp"," -k ",
      " -o ", .main$out_files$lumpt$vcf
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



