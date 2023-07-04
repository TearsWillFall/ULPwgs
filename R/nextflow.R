#' Run ampliconArchitect suite on nextflow
#' 
#' @param sif Path to SIF file. Default NULL
#' @param reference_dir Path to reference directory for AA
#' @param nf_circdna Nextflow circDNA
#' @export

aa_nextflow=function(
    sif=NULL,
    reference_dir=build_default_reference_list()$HG19$aa,
    nf_circdna=build_default_nf_list()$nf_circdna,
    resume=TRUE,
    ...
){
    run_main=function(
        .env
    ){
        
            .this.env=environment()
            append_env(to=.this.env,from=.env)

            set_main(.env=.this.env)

            .main$exec_code=paste(
                set_nf_envir(),
                " nextflow run ",nf_circdna$name,
                " -r ",nf_circdna$version,
                " --input ",normalizePath(input),
                " --input_format BAM ",
                " --outdir ",out_file_dir, 
                " --genome GRCh37 ",
                " -profile singularity ",
                " --circle_identifier ampliconarchitect ",
                " --reference_build GRCh37 ",
                " --mosek_license_dir ",license_dir,
                " --aa_data_repo ", reference_dir,
                " --bam_sorted ", 
                " --skip_qc ",
                " --skip_multiqc ", 
                " --skip_markduplicates ", 
                " --max_cpus ", threads,
                " --max_memory ", ram, 
                ifelse(resume," -resume ","")
            )

            run_job(.env=.this.env)
            .main.step=.main$steps[[fn_id]]
    
            .env$.main <- .main   
    }
         
    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="sif"
    )

    launch(.env=.base.env)

}

