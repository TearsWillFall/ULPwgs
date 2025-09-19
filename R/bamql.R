

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

    .main$out_files$accepted_reads_bam$unsorted$bam=paste0(tmp_dir,"/",input_id,".",accepted_id,".",input_ext)
    .main$out_files$rejected_reads_bam$unsorted$bam=paste0(tmp_dir,"/",input_id,".",rejected_id,".",input_ext)
    .main$exec_code=paste0(
      "singularity exec ",sif_bamql,
      " bamql -b -f ",input, 
      " -o ", .main$out_files$accepted_reads_bam$unsorted$bam, 
      " -O ", .main$out_files$rejected_reads_bam$unsorted$bam,
      " \'",query,"\'"
    )

    run_job(.env=.this.env)

    
    .main.step<-.main$steps[[fn_id]]

    if(sort){
        .main.step$steps <-append(
          .main.step$steps ,new_sort_and_index_bam_samtools(
              bin_samtools = bin_samtools,
              bam=.main$out_files$accepted_reads_bam$unsorted$bam,
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
            bam=.main$out_files$rejected_reads_bam$unsorted$bam,
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


#' Query BAM File Using BAMQL and Calculate Fragment Size Metrics
#'
#' This function queries a BAM file using BAMQL, separates reads into accepted and rejected BAM files, and then calculates fragment size metrics for both sets using Picard. Output BAM files can be sorted, indexed, and statistics generated as needed.
#'
#' @param bin_picard Path to the Picard binary. Default: from build_default_tool_binary_list().
#' @param bin_samtools Path to the samtools binary. Default: from build_default_tool_binary_list().
#' @param sif_bamql Path to the BAMQL singularity image file. Default: from build_default_sif_list().
#' @param bam Path to the input BAM file. (Required)
#' @param query BAMQL query string to filter reads. (Required)
#' @param accepted_id Identifier for accepted reads output file. Default: "accepted".
#' @param rejected_id Identifier for rejected reads output file. Default: "rejected".
#' @param sort Logical; whether to sort output BAM files. Default: TRUE.
#' @param coord_sort Logical; whether to coordinate sort output BAM files. Default: TRUE.
#' @param index Logical; whether to index output BAM files. Default: TRUE.
#' @param stats Logical; whether to generate statistics for output BAM files. Default: TRUE.
#' @param deviations Numeric; number of standard deviations for fragment size metrics. Default: NULL.
#' @param min_width Numeric; minimum fragment width for metrics. Default: NULL.
#' @param width Numeric; fragment width for metrics. Default: NULL.
#' @param ... Additional arguments passed to environment setup and job execution.
#'
#' @return No direct return value. Output files and metrics are written to disk and tracked in the environment.
#' @export



query_and_fragment_size_bamql=function(
  bin_picard=build_default_tool_binary_list()$bin_picard,
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
  deviations=NULL,
  min_width=NULL,
  width=NULL,
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

    .main.step$steps=append(
        .main.step$steps,query_bam_bamql(
          bin_samtools=bin_samtools,
          sif_bamql=sif_bamql,
          bam=bam,
          query=query,
          accepted_id=accepted_id,
          rejected_id=rejected_id,
          sort=sort,
          coord_sort=coord_sort,
          index=index,
          stats=stats,
          tmp_dir=tmp_dir,
          output_dir=paste0(out_file_dir,"/bams"),
          env_dir=env_dir,
          batch_dir=batch_dir,
          ram=ram,
          verbose=verbose,
          threads=threads,
          coord_sort=coord_sort,
          err_msg=err_msg,
          clean=clean,
          executor_id=task_id
      ) 
    )
      .this.step=.main.step$steps$query_bam_bamql
      .main.step$out_files=.this.step$out_files
   
      .main.step$steps=append(
          .main.step$steps ,
          new_insertsize_metrics_bam_picard(
            bin_picard=bin_picard,
            bam=.main.step$out_files$accepted_reads_bam$sorted$srt_bam,
            deviations=deviations,
            min_width=min_width,
            width=width,
            tmp_dir=tmp_dir,
            output_dir=paste0(out_file_dir,"/insert_size/accepted"),
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
      .this.step=.main.step$steps$new_insertsize_metrics_bam_picard.accepted
      .main.step$out_files$fragments=.this.step$out_files

       .main.step$steps=append(
          .main.step$steps ,
          new_insertsize_metrics_bam_picard(
            bin_picard=bin_picard,
            bam=.main.step$out_files$rejected_reads_bam$sorted$srt_bam,
            deviations=deviations,
            min_width=min_width,
            width=width,
            tmp_dir=tmp_dir,
            output_dir=paste0(out_file_dir,"/insert_size/rejected"),
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
      .this.step=.main.step$steps$new_insertsize_metrics_bam_picard.rejected
      .main.step$out_files$fragments=.this.step$out_files
  
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



