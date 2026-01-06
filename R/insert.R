#' Extract insert size metrics from BAM file using samtools
#'
#' This function extracts the sequence motifs from the ends of fragments in a BAM file using samtools and awk. The motifs are counted, sorted, and written to an output file. Useful for analyzing fragment end composition (e.g., for nucleosome positioning or bias detection).
#'
#' @param bin_samtools Path to the samtools binary. Default: from build_default_tool_binary_list().
#' @param bam Path to the input BAM file. (Required)
#' @param region Optional genomic region to restrict the analysis (e.g., "chr1:1000-2000"). Default is NULL (whole BAM).
#' @param ... Additional arguments passed to environment setup and job execution.
#'
#' @return No direct return value. Output file with fragment end motif counts is written to disk and tracked in the environment.
#' @export

#' @export

insertsize_metrics_samtools=function(
    bin_samtools=build_default_tool_binary_list()$bin_samtools,
    bam=NULL,
    region=NULL,
    ...
){

    run_main=function(
        .env
    ){
        .this.env=environment()
        append_env(to=.this.env,from=.env)
   
        set_main(.env=.this.env)

        .main$out_files$insert_size=paste0(out_file_dir,"/",
            input_id,ifelse(is.null(region),
            "",paste0(".",region)),"_insert_size.txt")
    
        .main$exec_code=paste(
            "echo \'N insert_size region id\' > ",
            .main$out_files$insert_size,";",
            bin_samtools," view ",
            input," -@ ",
            threads,
            region,
            paste0(" | awk \'{print ",
            "sqrt($9^2)\" ", ifelse(is.null(region),"genome",region)," ",input_id),
            "\"}\' | sort -n | uniq -c"
        )

        .main$exec_code=paste0(.main$exec_code,">>",.main$out_files$insert_size)
  

    run_job(.env=.this.env)
    .env$.main <- .main
        
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
