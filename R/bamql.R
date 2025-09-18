

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
  bin_samtools=build_default_tool_binary_list()$bin_samtools,
  sif_bamql=build_default_sif_list()$sif_bamql,
  bam=NULL,
  query=NULL,
  accepted_id="accepted",
  rejected_id="rejected",
  sort=TRUE,
  coord_sort=TRUE,
  index=TRUE,
  stats=TRUE,
  ...
){

  run_main=function(
    .env
  ){
    .this.env=environment()
    append_env(to=.this.env,from=.env)
    set_main(.env=.this.env)

    .main$steps[[fn_id]]<-.this.env
    .main.step=.main$steps[[fn_id]]

    if(is.null(query)){
        stop("Provide a query")
    }

    .main$out_files$accepted_reads_bam$unsorted=paste0(tmp_dir,"/",input_id,".",accepted_id,".",input_ext)
    .main$out_files$rejected_reads_bam$unsorted=paste0(tmp_dir,"/",input_id,".",rejected_id,".",input_ext)
    .main$exec_code=paste0(
      "singularity exec ",sif_bamql,
      " bamql -b -f ",input, 
      " -o ", .main$out_files$accepted_reads_bam$unsorted, 
      " -O ", .main$out_files$rejected_reads_bam$unsorted,
      " \'",query,"\'"
    )

    run_job(.env=.this.env)

    if(sort){
        .main.step$steps <-append(
          .main.step$steps ,new_sort_and_index_bam_samtools(
              bin_samtools = bin_samtools,
              bam=.main$out_files$accepted_reads_bam$unsorted,
              index=index,
              stats=stats,
              tmp_dir=tmp_dir,
              output_dir=paste0(out_file_dir,"/accepted"),
              env_dir=env_dir,
              batch_dir=batch_dir,
              ram=ram,
              verbose=verbose,
              threads=threads,
              coord_sort=coord_sort,
              err_msg=err_msg,
              clean=clean,
              executor_id=task_id,
              fn_id="accepted"
            )
        )
        .this.step=.main.step$steps$new_sort_and_index_bam_samtools.accepted
        .main.step$out_files$accepted_reads_bam$sorted=.this.step$out_files

        .main.step$steps <-append(
        .main.step$steps ,new_sort_and_index_bam_samtools(
            bin_samtools = bin_samtools,
            bam=.main$out_files$rejected_reads_bam$unsorted,
            index=index,
            stats=stats,
            tmp_dir=tmp_dir,
            output_dir=paste0(out_file_dir,"/rejected"),
            env_dir=env_dir,
            batch_dir=batch_dir,
            ram=ram,
            verbose=verbose,
            threads=threads,
            coord_sort=coord_sort,
            err_msg=err_msg,
            clean=clean,
            executor_id=task_id,
            fn_id="rejected"
          )
      )
      .this.step=.main.step$steps$new_sort_and_index_bam_samtools.rejected
      .main.step$out_files$rejected_reads_bam$sorted=.this.step$out_files
    }
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
