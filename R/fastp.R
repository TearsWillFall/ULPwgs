



#' Quality Control and Trimming of FASTQ Files using fastp
#'
#' This function performs quality control and adapter trimming on paired-end FASTQ files
#' using the fastp tool. It removes low-quality reads, trims adapters, and generates
#' comprehensive quality reports in both JSON and HTML formats.
#'
#' The function applies the following parameters:
#' \itemize{
#'   \item Quality threshold: Q20
#'   \item Minimum length: 75bp
#'   \item Sliding window size: 5bp
#'   \item Unqualified base limit: 40\%
#'   \item Poly-G tail removal enabled
#'   \item Per-base quality correction enabled
#' }
#'
#' For more information on fastp, visit:
#' https://github.com/OpenGene/fastp
#'
#' @param env_fastp [REQUIRED] Path to fastp conda environment.
#' @param fastq [REQUIRED] Path to input paired-end FASTQ files or data structure containing R1 and R2 file paths.
#' @param ... Additional arguments passed to internal functions.
#'
#' @return A list containing:
#'   \item{trim_fastq_r1}{Path to trimmed R1 FASTQ file}
#'   \item{trim_fastq_r2}{Path to trimmed R2 FASTQ file}
#'   \item{fastp_json}{Path to fastp JSON report}
#'   \item{fastp_html}{Path to fastp HTML report}
#'
#' @export
#' 

trim_umi_fastp=function(
  env_fastp=build_default_python_enviroment_list()$env_fastp,
  fastq=NULL,
  output_name="sample",
  ...
){
   run_main=function(
    .env
  ){
    .this.env=environment()
    append_env(to=.this.env,from=.env)
    set_main(.env=.this.env)

    .main$out_files$fastq_r1=paste0(out_file_dir,"/",input_id,".trimmed_R1.fastq")
    .main$out_files$fastq_r2=paste0(out_file_dir,"/",input_id,".trimmed_R2.fastq")
    .main$out_files$fastp_json=paste0(out_file_dir,"/",input_id,".fastp.json")
    .main$out_files$fastp_html=paste0(out_file_dir,"/",input_id,".fastp.html")

    .main$exec_code=paste(
        "conda activate ",env_fastp,
        "; fastp -i ", input$fastq_r1, "-I ", input$fastq_r2, 
        " -o ", .main$out_files$fastq_r1,
        " --out2 ", .main$out_files$fastq_r2,
        " -g -W 5 -q 20 -u 40 -x -3 -l 75 -c",
        " -j ",.main$out_files$fastp_json,
        " -h ",.main$out_files$fastp_html, 
        " -w ",threads
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







