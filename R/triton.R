


run_triton=function(
  env_triton=build_default_python_enviroment_list()$env_triton,
  bin_triton=build_default_tool_binary_list()$bin_triton,
  ref_genome=build_default_reference_list()$HG19$reference$genome,
  nc_dict=build_default_reference_list()$OTHER$nc_fitting,
  bam_gc=NULL,
  bed=NULL,
  method="region",
  ...
){
     run_main=function(
    .env
    ){
      .this.env=environment()
      append_env(to=.this.env,from=.env)
      set_main(.env=.this.env)
      
      .main$out_files$skipped_sites=paste0(
        out_file_dir,"/",get_file_name(input$bam),"/",
        get_file_name(input$bam),"_SkippedSites.txt"
    )

    .main$out_files$triton_features=paste0(
        out_file_dir,"/",get_file_name(input$bam),"/",
        get_file_name(input$bam),"_TritonFeatures.tsv"
    )

    .main$out_files$triton_fragmentation_profiles=paste0(
        out_file_dir,"/",get_file_name(input$bam),"/",
        get_file_name(input$bam),"_TritonFragmentationProfiles.npz"
    )

    .main$out_files$triton_signal_profiles=paste0(
        out_file_dir,"/",get_file_name(input$bam),"/",
        get_file_name(input$bam),"_TritonSignalProfiles.npz"
    )

    

      .main$exec_code=paste0(
        paste0(
          "conda activate ",env_triton,
          "; ",bin_triton,
          " -n ",get_file_name(input$bam),
          " -i ",input$bam,
          " -b ",input$gc_bias,
          " -g ",ref_genome,
          " -r ",out_file_dir,
          " -m ", method,
          " -c ", threads,
          " -d ", nc_dict,
          " -a ", bed
          )
        )

       
      
      run_job(.env=.this.env)
      .env$.main <- .main
    }

   .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="bam_gc"
    )
    launch(.env=.base.env)
}

