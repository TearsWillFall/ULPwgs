#' Wrapper around assuite
#' 
#'
#' Wrapper around AmpliconArchitectSuite
#' 
#' @param sif_aa Path to AmpliconArchitect sif
#' @param bin_aasuite Path to AmpliconArchitect-suite 
#' @param bam Path to BAM file.
#' @param cns Path to BED file with seeds or segmentation
#' @param cnsize_min Minimum cnsize
#' @param cn_gain Minimum cn_gain
#' @param downsample Downsample BAM file
#' @param run_AA Run AmpliconArchitect
#' @param run_AC Run AmpliconClassifier
#' @param AA_runmode AmpliconArchitect runmode
#' @param AA_extendmode AmpliconArchitect extendedmode
#' @param AA_insert_sdevs AmpliconArchitect inster sdevs
#' @export
#' 
call_aasuite=function(
    sif_aa="/lustre/scratch/scratch/regmova/Singularity_Images/ampliconsuite-pipeline.sif",
    bin_aasuite="/lustre/scratch/scratch/regmova/tools/AmpliconSuite-pipeline/singularity/run_paa_singularity.py",
    cns=NULL,
    bam=NULL,
    cn_gain=4.5,
    cnsize_min=50000,
    downsample=10,
    run_AA=TRUE,
    run_AC=TRUE,
    AA_runmode="FULL",
    AA_extendmode="EXPLORE",
    AA_insert_sdevs=3,
    ...
    ){

        run_main=function(
        .env
        ){
        
            .this.env=environment()
            append_env(to=.this.env,from=.env)
            set_main(.env=.this.env)
            .main$exec_code=paste(
            bin_aasuite,
            " --sif  ",dirname(sif_aa), 
            " -s ",input_id, 
            " -t ", threads,
            " --cnv_bed ",cns, 
            " --bam ", bam, 
            " --cngain ",cn_gain,
            " -o ",paste0(out_file_dir,"/",input_id),
            " --cnsize_min ",cnsize_min,
            " --downsample ", downsample,
            " --AA_runmode ", AA_runmode,
            " --AA_extendmode ", AA_extendmode,
            " --AA_insert_sdevs ",AA_insert_sdevs,
            ifelse(run_AA," --run_AA "," "),
            ifelse(run_AC," --run_AC "," ")
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