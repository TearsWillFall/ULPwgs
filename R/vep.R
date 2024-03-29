#' Create Read Orientation Model for Mutect2 Filter
#'
#' This function creates read orientation model for mutect2 filtering
#'
#' @param bin_vep [REQUIRED] Path to VEP binary.
#' @param bin_bgzip [REQUIRED] Path to bgzip binary.
#' @param bin_tabix [REQUIRED] Path to tabix binary.
#' @param vep_cache [REQUIRED] Path to vep cache location.
#' @param vcf [REQUIRED] Path to VCF file.
#' @param compress [OPTIONAL] Generate a compressed VCF
#' @param output_name [OPTIONAL] Name of output file
#' @param clean [OPTIONAL]Remove extra files.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param threads [OPTIONAL] Number of threads to split the work. Default 4
#' @param ram [OPTIONAL] RAM memory to asing to each thread. Default 4
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Task EXECUTOR ID. Default "recalCovariates"
#' @param task_name Task name. Default "recalCovariates"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export


annotate_vep=function(
    bin_vep=build_default_tool_binary_list()$bin_vep,
    bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
    bin_tabix=build_default_tool_binary_list()$bin_tabix,
    cache_vep=build_default_cache_list()$cache_vep,
    vcf=NULL,
    tabulate=TRUE,
    compress=TRUE,
    patient_id=NULL,
    tumour_id=NULL,
    normal_id=NULL,
    chromosomes=NULL,
    fmt="vcf",
    ...
){




    run_main=function(
      .env
    ){


      .this.env=environment()
      append_env(to=.this.env,from=.env)

      set_main(.env=.this.env)


      chr=""
      chrms=""
      if(!is.null(chromosomes)){
        chr=paste0(".",paste0(chromosomes,collapse="_"))
        chrms=paste0("--chr ",paste0(chromosomes,collapse=","))
      }


      .main$out_files$vep_vcf=paste0(
        out_file_dir,"/",input_id,".vep",chr,".vcf")

      .main$exec_code=paste(bin_vep," -format ",fmt,
        "-i",vcf,"-o", .main$out_files$vep_vcf,
        "--cache --offline --everything --force_overwrite --vcf --fork ",
        threads," --dir ",cache_vep,chrms
      )


      run_job(.this.env)
      .main.step=.main$steps[[fn]]


   

      if(compress){
          .main.step$steps<-append(
            .main.step$steps,
            compress_and_index_vcf_htslib(
              bin_bgzip=bin_bgzip,
              bin_tabix=bin_tabix,
              vcf=.main$out_files$vep_vcf,
              output_name=paste0(input_id,".vep",chr),
              compress=compress,
              index_format=index_format,
              bgzip_index=bgzip_index,
              output_dir=out_file_dir,
              tmp_dir=tmp_dir,
              env_dir=env_dir,
              batch_dir=batch_dir,
              err_msg=err_msg,
              verbose=verbose,
              threads=threads,
              ram=ram
            ) 
        )
        .this.step=.main.step$steps$compress_and_index_vcf_htslib
        .main.step$out_files <- append(
          .main.step$out_files,
          .this.step$out_files
          )
    }

       if(tabulate){
              .main.step$steps<-append(
                .main.step$steps,
                tabulate_vcf(
                  vcf=.main$out_files$vep_vcf,
                  output_dir=paste0(out_file_dir,"/tabulated"),
                  output_name=output_name,
                  patient_id=patient_id,
                  tumour_id=tumour_id,
                  normal_id=normal_id,
                  chromosomes=chromosomes,
                  tmp_dir=tmp_dir,
                  env_dir=env_dir,
                  batch_dir=batch_dir,
                  err_msg=err_msg,
                  threads=threads,
                  ram=ram,
                  verbose=verbose,
                  executor_id=task_id
              )
            )

          .this.step=.main.step$steps$tabulate_vcf
          .main.step$out_files <- append(
            .main.step$out_files,
            .this.step$out_files
          )
      }
    .env$.main<-.main

  }


  .base.env=environment()
  list2env(list(...),envir=.base.env)
  set_env_vars(
      .env=.base.env,
      vars="vcf"
  )

  launch(.env=.base.env)



}


#' Create Read Orientation Model for Mutect2 Filter
#'
#' This function creates read orientation model for mutect2 filtering
#'
#' @param bin_vep [REQUIRED] Path to VEP binary.
#' @param bin_bgzip [REQUIRED] Path to bgzip binary.
#' @param bin_tabix [REQUIRED] Path to tabix binary.
#' @param vep_cache [REQUIRED] Path to vep cache location.
#' @param vcf [REQUIRED] Path to VCF file.
#' @param compress [OPTIONAL] Generate a compressed VCF
#' @param output_name [OPTIONAL] Name of output file
#' @param clean [OPTIONAL]Remove extra files.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param threads [OPTIONAL] Number of threads to split the work. Default 4
#' @param ram [OPTIONAL] RAM memory to asing to each thread. Default 4
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Task EXECUTOR ID. Default "recalCovariates"
#' @param task_name Task name. Default "recalCovariates"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export


annotate_strelka_vep=function(
    bin_vep=build_default_tool_binary_list()$bin_vep,
    bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
    bin_tabix=build_default_tool_binary_list()$bin_tabix,
    cache_vep=build_default_cache_list()$cache_vep,
    vcf=NULL,
    patient_id=NULL,
    tumour_id=NULL,
    normal_id=NULL,
    fmt="vcf",
    tabulate=TRUE,
    type="snv",
    chromosomes=NULL,
    ...
){

    run_main=function(
      .env
    ){


      .this.env=environment()
      append_env(to=.this.env,from=.env)

      set_main(.env=.this.env)
      output_name=paste0(input_id,".",type)
      fn=paste0(fn,".",type)

      .main$steps[[fn]]<-.this.env
      .main.step=.main$steps[[fn]]
     

      .main$steps[[fn]]$steps<-append(
        .main$steps[[fn]]$steps,
        annotate_vep(
          bin_bgzip=bin_bgzip,
          bin_tabix=bin_tabix,
          bin_vep=bin_vep,
          cache_vep=cache_vep,
          vcf=input,
          chromosomes=chromosomes,
          tabulate=tabulate,
          patient_id=patient_id,
          tumour_id=tumour_id,
          normal_id=normal_id,
          index_format=index_format,
          bgzip_index=bgzip_index,
          output_dir=out_file_dir,
          output_name=output_name,
          tmp_dir=tmp_dir,
          env_dir=env_dir,
          batch_dir=batch_dir,
          err_msg=err_msg,
          verbose=verbose,
          threads=threads,
          ram=ram
        ) 
      )

     .this.step=.main$steps$annotate_vep
     .main.step$out_files[[type]] <- .this.step$out_files
     .env$.main<-.main

   }


    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
        .env=.base.env,
        vars="vcf"
    )

    launch(.env=.base.env)


}


