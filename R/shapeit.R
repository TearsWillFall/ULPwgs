
#' Strelka wrapper for SNV variant calling
#'
#' This function wraps the STRELKA functions for variant calling
#' 
#' @param bin_shapeit Path to SHAPEIT binary file
#' @param haplotype Path to haplotype information from 1000G
#' @param legend Path to legend information from 1000G
#' @param map Path to map information from 1000G
#' @param vcf Path to VCF file
#' @export


phase_shapeit=function(
    bin_shapeit=build_default_tool_binary_list()$bin_shapeit,
    ref_panel=build_default_reference_list()$HG19$phasing$G1000$new$vcf,
    gmap=build_default_reference_list()$HG19$phasing$G1000$new$gmap,
    vcf=NULL,
    chr=NULL,
    scaffold=NULL,
    ...
){
    run_main=function(.env){

        .this.env=environment()
        append_env(to=.this.env,from=.env)
        set_main(.env=.this.env)

        .main$out_files$phased_vcf=paste0(out_file_dir,"/",input_id,".",chr,".phased.vcf")

        if(!is.null(ref_panel)){
            add=paste0(" --reference ",ref_panel[grepl(paste0("chr",chr,"."),ref_panel)]) 
        }

        if(!is.null(scaffold)){
            add=paste0(" --scaffold ",scaffold)
        }
        .main$exec_code=paste(
            bin_shapeit," --input ",input,
            " --map ",gmap[grepl(paste0("chr",chr,"."),gmap)],
            add,
            " --region ",chr,
            " --thread ",threads,
            " --output ",.main$out_files$phased_vcf,
            " --sequencing"
        )

        run_job(.env=.this.env)

        .env$.main <- .main
        
    }

    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="vcf"
    )
  
    launch(.env=.base.env)

}


