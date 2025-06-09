#' Generate MSI report for each sample using MSING
#' https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5904988/
#' 
#' @param env_msing Path to MSING environment.
#' @param bin_msing Path to MSING binary.
#' @param regions Path to BED with MSI regions.
#' @param baseline Path to noise background in reference samples.
#' @param ref_genome Path to the reference genome.
#' @param bam Path to BAM file.
#' @export


call_msing=function(
    env_msing=build_default_python_enviroment_list()$env_msing,
    bin_msing=build_default_tool_binary_list()$bin_msing,
    regions=build_default_reference_list()$HG19$panel$PCF_V3$msi$bed,
    baseline=build_default_reference_list()$HG19$panel$PCF_V3$msi$baseline,
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    bam=NULL,
    ...
){

    options(scipen = 999)
    run_main=function(
        .env
    ){
        .this.env=environment()
        append_env(to=.this.env,from=.env)
   
        set_main(.env=.this.env)

        add=""
        if(!is.null(regions)){
          add=paste("-T ",regions)
        }
        

        .main$out_files$msi_summary<-paste0(out_file_dir,"/",input_id,"/",input_id,".MSI_Analysis.txt")
        .main$out_files$msi<-paste0(out_file_dir,"/",input_id,"/",input_id,".msi.txt")
        .main$out_files$mpileup<-paste0(out_file_dir,"/",input_id,"/",input_id,".mpileup.txt")
        .main$exec_code=paste(
                ". ",
                bin_msing,
                env_msing,
                input,
                regions,
                baseline,
                ref_genome,
                out_file_dir
        )
        run_job(.env=.this.env)
        .env$.main <- .main
    }

    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env=.base.env,
      vars="bam"
    )

     launch(.env=.base.env)

}