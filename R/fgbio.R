

#' Extract UMI Tags from BAM file using fgbio
#'
#' This function extracts Unique Molecular Identifier (UMI) tags from BAM files.
#' The output is a BAM file with UMI information extracted and stored in the RX tag.
#' This tool is useful for deduplication and error correction in sequencing data with UMI barcodes.
#' 
#' For more information read:
#' https://fulcrumgenomics.github.io/fgbio/tools/latest/ExtractUmisFromBam.html
#'
#' @param env_fgbio [REQUIRED] Path to fgbio conda environment.
#' @param bam [REQUIRED] Path to input BAM file containing UMI sequences.
#' @export
#' 

extract_umi_fgbio=function(
  env_fgbio=build_default_python_enviroment_list()$env_fg_bio,
  bam=NULL,
  ...
){
   run_main=function(
    .env
  ){
    .this.env=environment()
    append_env(to=.this.env,from=.env)
    set_main(.env=.this.env)

    .main$out_files$bam=paste0(out_file_dir,"/",input_id,".umi.bam")

    .main$exec_code=paste(
      "conda activate ",env_fgbio,
      "; fgbio ExtractUmisFromBam -i",input,
      " -o ",.main$out_files$bam,
      " -r 3M3S+T 3M3S+T -t RX -a true"
    )

     run_job(.env=.this.env)

    .env$.main<-.main
  } 
    
   .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="bam"
    )

    launch(.env=.base.env)

}



#' Group Reads by UMI using fgbio
#'
#' This function groups sequencing reads by their Unique Molecular Identifier (UMI) tags
#' using fgbio's GroupReadsByUmi tool. It clusters reads with the same UMI and similar
#' sequences to identify derived sequences from the same original DNA molecule.
#'
#' Grouping strategy:
#' \itemize{
#'   \item Strategy: adjacency (connects UMIs that differ by up to --edits mismatches)
#'   \item Edit distance: 1 (allows for single nucleotide variations)
#'   \item UMI tag: RX
#'   \item Generates family size count statistics
#' }
#'
#' This is typically used after UMI extraction to group reads from the same original
#' molecule for downstream deduplication and consensus calling.
#'
#' For more information read:
#' https://fulcrumgenomics.github.io/fgbio/tools/latest/GroupReadsByUmi.html
#'
#' @param env_fgbio [REQUIRED] Path to fgbio conda environment.
#' @param bam [REQUIRED] Path to input BAM file with extracted UMI tags in the RX field.
#' @param ... Additional arguments passed to internal functions for environment setup and job execution.
#'
#' @return List of output files:
#'   \item{grouped_umi_bam}{BAM file with reads grouped by UMI}
#'   \item{family_size_counts}{Text file with UMI family size distribution statistics}
#'
#' @seealso
#'   \link{exctract_umi_fgbio} for extracting UMI tags from fastq/BAM files
#'
#' @export
#' 

group_by_umi_fgbio=function(
  env_fgbio=build_default_python_enviroment_list()$env_fg_bio,
  bam=NULL,
  ...
){
   run_main=function(
    .env
  ){
    .this.env=environment()
    append_env(to=.this.env,from=.env)
    set_main(.env=.this.env)

    .main$out_files$bam=paste0(out_file_dir,"/",input_id,".grouped.bam")
    .main$out_files$family_size_counts=paste0(out_file_dir,"/",input_id,".family_size_counts.txt")
    .main$exec_code=paste(
      "conda activate ",env_fgbio,
      "; fgbio GroupReadsByUmi --input=",input,
      " --output=",.main$out_files$bam,
      " --strategy=adjacency
      --edits=1
      -t RX
      -f ",.main$out_files$family_size_counts
    )

     run_job(.env=.this.env)

    .env$.main<-.main
  } 
    
   .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="bam"
    )

    launch(.env=.base.env)

}







#' Call Molecular Consensus Reads from UMI-Grouped BAM using fgbio
#'
#' This function calls consensus sequences from reads grouped by UMI tags using fgbio's
#' CallMolecularConsensusReads tool. It generates high-confidence consensus reads by combining
#' information from multiple reads derived from the same original DNA molecule.
#'
#' Consensus calling parameters:
#' \itemize{
#'   \item Error rate post-UMI: 40 (quality score, post-consensus error rate)
#'   \item Error rate pre-UMI: 45 (quality score, pre-consensus error rate)
#'   \item Minimum reads per consensus: 2
#'   \item Maximum reads per consensus: 50
#'   \item Minimum input base quality: 20
#'   \item Per-base consensus tags: disabled
#'   \item Read name prefix: "consensus"
#' }
#'
#' This is typically used after grouping reads by UMI to generate high-fidelity consensus
#' sequences for error correction and artifact removal in sequencing data.
#'
#' For more information read:
#' https://fulcrumgenomics.github.io/fgbio/tools/latest/CallMolecularConsensusReads.html
#'
#' @param env_fgbio [REQUIRED] Path to fgbio conda environment.
#' @param bam [REQUIRED] Path to input BAM file with UMI-grouped reads (from group_by_umi_fgbio).
#' @param ... Additional arguments passed to internal functions for environment setup and job execution.
#'
#' @return Consensus BAM file with molecular consensus reads. Output path tracked in job report.
#'
#' @seealso
#'   \link{group_by_umi_fgbio} for grouping reads by UMI before consensus calling
#'
#' @export
#' 

call_consensus_fgbio=function(
  env_fgbio=build_default_python_enviroment_list()$env_fg_bio,
  bam=NULL,
  ...
){
   run_main=function(
    .env
  ){
    .this.env=environment()
    append_env(to=.this.env,from=.env)
    set_main(.env=.this.env)

    .main$out_files$bam=paste0(out_file_dir,"/",input_id,".consensus.unmapped.bam")
   
    .main$exec_code=paste(
      "conda activate ",env_fgbio,
      "; fgbio CallMolecularConsensusReads --input=" ,input,
      " --output=",.main$out_files$bam,
      " --error-rate-post-umi 40
        --error-rate-pre-umi 45
        --output-per-base-tags false
        --min-reads 2
        --max-reads 50
        --min-input-base-quality 20
        --read-name-prefix=\'consensus\'"
    )

     run_job(.env=.this.env)

    .env$.main<-.main
  } 
    
   .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="bam"
    )

    launch(.env=.base.env)

}






