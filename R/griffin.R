#' Calculate GC Frequency Distribution from Reference Genome Using GRIFFIN
#'
#' This function computes the GC frequency distribution across mappable regions of a reference genome using GRIFFIN. It generates TSV files containing GC frequency data for specified genomic size ranges. The function uses a Snakemake workflow to process the mappable regions and produce frequency statistics at the specified bin sizes.
#'
#' @param env_griffin Path to the GRIFFIN conda environment. Default: from build_default_python_enviroment_list().
#' @param bin_griffin Path to the GRIFFIN scripts directory containing the Snakemake pipeline.
#' @param mappable_bed Path to the BED file containing mappable regions of the reference genome. Default: HG19 100bp mappable regions.
#' @param ref_genome Path to the reference genome FASTA file. Default: HG19 reference genome.
#' @param chrom_size Path to the chromosome size file. Default: HG19 chromosome sizes.
#' @param range Numeric vector specifying bin sizes (in bp) for which to calculate GC frequency. Default: c(1:501).
#' @param read_length Integer specifying the sequencing read length. Default: 100.
#' @param ... Additional arguments passed to environment setup and job execution.
#'
#' @return No direct return value. Output TSV files are written to disk and tracked in the environment under .main$out_files$gc_frequency.
#' @export


genome_GC_frequency_griffin=function(
  env_griffin=build_default_python_enviroment_list()$env_griffin,
  bin_griffin=build_default_tool_binary_list()$bin_griffin,
  sf_genome_gc=build_default_tool_snakemake_list()$sf_genome_gc,
  mappable_bed=build_default_reference_list()$HG19$reference$mappable_bed_100bp,
  ref_genome=build_default_reference_list()$HG19$reference$genome,
  chrom_size=build_default_reference_list()$HG19$reference$chrom_size,
  range=c(1:500),
  read_length=100,
  ...
){
     run_main=function(
    .env
    ){
      .this.env=environment()
      append_env(to=.this.env,from=.env)
      set_main(.env=.this.env)
      
      .main$out_files$gc_frequency=paste0(
        out_file_dir,"/",
        sub(".bed","",mappable_bed),".",
        input,".bp.GC_frequency.tsv")

      build_griffin_config_snakemake(.env=.main)
  
     
      
      .main$exec_code=paste0(
        paste0(
          "conda activate ",env_griffin,
          "; snakemake -s ",sf_genome_gc,
          " --configfile ",.main$out_files$config, 
          " --cores ",threads)
        )

      
      run_job(.env=.this.env)
      .env$.main <- .main
    }

   .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="range"
    )
    launch(.env=.base.env)
}

#' Build GRIFFIN Snakemake Configuration File
#'
#' This internal helper function generates a YAML configuration file required for the GRIFFIN Snakemake workflow. It sets up paths to input files (mappable regions, reference genome, chromosome sizes), output directory, size range parameters, and read length specifications.
#'
#' @param .env Environment object containing pipeline variables including out_file_dir, bin_griffin, mappable_bed, ref_genome, chrom_size, range, and read_length.
#'
#' @return No direct return value. A YAML configuration file is written to disk at .main$out_files$config.
#' @keywords internal
#'
build_griffin_config_snakemake=function(.env=NULL){

    .this.env=environment()
    append_env(to=.this.env,from=.env)

    out_file_dir=paste0(out_file_dir,"/results")
    out_file_dir_cfg=paste0(out_file_dir,"/config")
  
    out_files$config=paste0(out_file_dir_cfg,"/griffin_GC_frequency_config.",input,".yaml")
    dir.create(out_file_dir_cfg,showWarnings = FALSE,recursive = TRUE)

    cat(x=
        paste0(
        "griffin_scripts_dir: ", bin_griffin,"\n\n",
        "mappable_regions: ", mappable_bed,"\n\n",
        "reference_genome: ", ref_genome,"\n\n",
        "chrom_sizes: ", chrom_size,"\n\n",
        "out_dir: ",  out_file_dir,"\n\n",
        "size_range: ", "[",paste0(c(min(input),max(input)), collapse = ","),"]\n\n",
        "read_length: ", read_length
    ),
    file=out_files$config
    )

  append_env(to=.env,from=.this.env)

}



#' Build GRIFFIN Snakemake Configuration File
#'
#' This internal helper function generates a YAML configuration file required for the GRIFFIN Snakemake workflow. It sets up paths to input files (mappable regions, reference genome, chromosome sizes), output directory, size range parameters, and read length specifications.
#'
#' @param .env Environment object containing pipeline variables including out_file_dir, bin_griffin, mappable_bed, ref_genome, chrom_size, range, and read_length.
#'
#' @return No direct return value. A YAML configuration file is written to disk at .main$out_files$config.
#' @keywords internal
#'
build_griffin_gc_config_snakemake=function(.env=NULL){

    .this.env=environment()
    append_env(to=.this.env,from=.env)

    out_file_dir=paste0(out_file_dir,"/results")
    out_file_dir_cfg=paste0(out_file_dir,"/config")
  
    out_files$config=paste0(out_file_dir_cfg,"/griffin_GC_bias_config.",input_id,".yaml")
    dir.create(out_file_dir_cfg,showWarnings = FALSE,recursive = TRUE)

    cat(x=
        paste0(
        "encode_exclude: ../../Ref/encode_unified_GRCh38_exclusion_list.bed\n\n",
        "centromeres: ../../Ref/hg38_centromeres.bed\n\n",
        "gaps: ../../Ref/hg38_gaps.bed\n\n",
        "patches: ../../Ref/hg38_fix_patches.bed\n\n",
        "alternative_haplotypes: ../../Ref/hg38_alternative_haplotypes.bed\n\n",
        "mappability_bw: ../../Ref/k100.Umap.MultiTrackMappability.bw","\n\n",
        "mappability_correction: False\n\n",
        "map_quality: ", map_quality,"\n\n",
        "griffin_scripts_dir: ", bin_griffin,"\n\n",
        "mappable_regions: ", mappable_bed,"\n\n",
        "reference_genome: ", ref_genome,"\n\n",
        "chrom_sizes: ", chrom_size,"\n\n",
        "out_dir: ",  out_file_dir,"\n\n",
        "genome_GC_frequency: ",gc_frequency,"\n\n",
        "size_range: ", paste0(c(min(input),max(input)), collapse = " "),"\n\n",
        "samples:
          sample_name_1:", input,"\n\n"
    ),
    file=out_files$config
    )

  append_env(to=.env,from=.this.env)

}












#' Calculate GC Frequency Distribution from Reference Genome Using GRIFFIN
#'
#' This function computes the GC frequency distribution across mappable regions of a reference genome using GRIFFIN. It generates TSV files containing GC frequency data for specified genomic size ranges. The function uses a Snakemake workflow to process the mappable regions and produce frequency statistics at the specified bin sizes.
#'
#' @param env_griffin Path to the GRIFFIN conda environment. Default: from build_default_python_enviroment_list().
#' @param bin_griffin Path to the GRIFFIN scripts directory containing the Snakemake pipeline.
#' @param mappable_bed Path to the BED file containing mappable regions of the reference genome. Default: HG19 100bp mappable regions.
#' @param ref_genome Path to the reference genome FASTA file. Default: HG19 reference genome.
#' @param chrom_size Path to the chromosome size file. Default: HG19 chromosome sizes.
#' @param range Numeric vector specifying bin sizes (in bp) for which to calculate GC frequency. Default: c(1:501).
#' @param read_length Integer specifying the sequencing read length. Default: 100.
#' @param ... Additional arguments passed to environment setup and job execution.
#'
#' @return No direct return value. Output TSV files are written to disk and tracked in the environment under .main$out_files$gc_frequency.
#' @export


GC_bias_correction_griffin=function(
  env_griffin=build_default_python_enviroment_list()$env_griffin,
  bin_griffin=build_default_tool_binary_list()$bin_griffin,
  sf_gc_correction=build_default_tool_snakemake_list()$sf_gc_correction,
  mappable_bed=build_default_reference_list()$HG19$reference$mappable_bed_100bp,
  ref_genome=build_default_reference_list()$HG19$reference$genome,
  gc_frequency=build_default_reference_list()$HG19$gc_frequency,
  chrom_size=build_default_reference_list()$HG19$reference$chrom_size,
  range=c(15:500),
  map_quality=20,
  bam=NULL,
  ...
){
     run_main=function(
    .env
    ){
      .this.env=environment()
      append_env(to=.this.env,from=.env)
      set_main(.env=.this.env)
      
      .main$out_files$gc_frequency=paste0(
        out_file_dir,"/GC_bias/",input,".GC_bias.txt")

      build_griffin_gc_config_snakemake(.env=.main)
  
     
      
      .main$exec_code=paste0(
        paste0(
          "conda activate ",env_griffin,
          "; snakemake -s ",sf_genome_gc,
          " --configfile ",.main$out_files$config, 
          " --cores ",threads)
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

