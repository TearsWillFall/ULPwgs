

#' BAM Querying using BAMQL
#'
#' This function queries a BAM file using BAMQL and separates reads into accepted
#' and rejected BAM files based on the provided BAMQL query.
#'
#' @param sif_bamql [REQUIRED] Path to bamql sif file.
#' @param bam [REQUIRED] Path to BAM file.
#' @param query [REQUIRED] BAM query in BAMQL format.
#' @param accepted_id [OPTIONAL] Id to use for accepted reads.
#' @param rejected_id [OPTIONAL] Id to use for rejected reads.
#' @export

query_bam_bamql=function(
  sif_bamql=build_default_sif_list()$sif_bamql,
  bam=NULL,
  query=NULL,
  accepted_id="accepted",
  rejected_id="rejected",
  ...
){

  run_main=function(
    .env
  ){
    .this.env=environment()
    append_env(to=.this.env,from=.env)
    set_main(.env=.this.env)

    if(is.null(query)){
        stop("Provide a query")
    }

    .main$out_files$accepted_reads_bam=paste0(out_file_dir,"/",input_id,".",accepted_id,".bam")
    .main$out_files$rejected_reads_bam=paste0(out_file_dir,"/",input_id,".",rejected_id,".bam")
    .main$exec_code=paste(
      "singularity exec -H ",paste0(getwd(),":/home "),sif_bamql,
      " bamql -b -f ",input, 
      " -o ", .main$out_files$accepted_reads_bam, 
      " -O ", .main$out_files$rejected_reads_bam,
      ,"\'",query,"\'"
    )

    run_job(.env=.this.env)
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
