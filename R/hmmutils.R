#' Generate a WIG file
#'
#' This function generates a WIG file.
#'
#'
#' @param bam Path to the BAM file .
#' @param bin_samtools Path to readCounter executable. Default path tools/samtools/samtools.
#' @param bin_readcount Path to readCounter executable. Default path tools/hmmcopy_utils/bin/readCounter.
#' @param output_dir Path to the output directory.
#' @param chrs String of chromosomes to include. c()
#' @param win Size of non overlaping windows. Default 500000.
#' @param format Output format [wig/seg] . Default wig
#' @param threads Number of threads to use. Default 3
#' @param verbose Enables progress messages. Default False.
#' @export


read_counter_region=function(
  bin_readcount=build_default_tool_binary_list()$bin_readcount,
  bam=NULL,
  region=NULL,
  win=500000,
  mapq=20,
  format="wig",
  ...
){  
    run_main=function(
        .env
    ){
        .this.env=environment()
        append_env(to=.this.env,from=.env)

        set_main(.env=.this.env)

        .main$out_files[[format]]=paste0(
            out_file_dir,"/",
            input_id,".",
            input,".",
            format
        )

        options(scipen=999)
        
        fmt=""
        if (format=="seg"){
            fmt="-s"
        }

        .main$exec_code=paste(
            bin_readcount,
            fmt,"--window", win,
            "--quality ", mapq ,
            "--chromosome ",input,
            bam,">", .main$out_files
        )

        run_job(.env=.this.env)
        .env$.main <- .main   
    }

    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="region"
    )

    launch(.env=.base.env)

  }


  #' Generate a WIG file
#'
#' This function generates a WIG file.
#'
#'
#' @param bam Path to the BAM file .
#' @param bin_samtools Path to readCounter executable. Default path tools/samtools/samtools.
#' @param bin_readcount Path to readCounter executable. Default path tools/hmmcopy_utils/bin/readCounter.
#' @param output_dir Path to the output directory.
#' @param chrs String of chromosomes to include. c()
#' @param win Size of non overlaping windows. Default 500000.
#' @param format Output format [wig/seg] . Default wig
#' @param threads Number of threads to use. Default 3
#' @param verbose Enables progress messages. Default False.
#' @export


read_counter=function(
  bin_samtools=build_default_tool_binary_list()$bin_samtools,
  bin_readcount=build_default_tool_binary_list()$bin_readcount,
  bam=NULL,
  region=NULL,
  patient_id=NULL,
  chromosomes=c(1:22,"X","Y"),
  win=500000,
  mapq=20,
  format="wig",
  ...
){  

    run_main=function(
        .env
    ){
        .this.env=environment()
        append_env(to=.this.env,from=.env)


        out_file_dir=set_dir(
            out_file_dir,
            name=paste0(patient_id,"/read_counts/",input_id)
        )

        set_main(.env=.this.env)

        .main$steps[[fn_id]]<-.this.env
        .main.step=.main$steps[[fn_id]]

        if(is.null(region)){
            .main.step$steps <-append(
                .main.step$steps,
                get_sq_bam(
                bin_samtools=bin_samtools,
                bam=input,
                output_name=input_id,
                output_dir=tmp_dir,
                header=TRUE,
                tmp_dir=tmp_dir,
                env_dir=env_dir,
                batch_dir=batch_dir,
                err_msg=err_msg,
                verbose=verbose,
                threads=threads,
                ram=ram,
                executor_id=task_id
                )
            )
            .this.step=.main.step$steps$get_sq_bam
            .main.step$out_files$region=.this.step$out_files$index_bed
            region=.main.step$out_files$region
        }

       
        if(length(region)==1){
            ### IF PATH READ AS BED
            if (file.exists(region)){
                region=read_bed(
                bed=region
                )$body
            }

            ### IF DATA.FRAME GENERATE GID
            if (is.data.frame(region)){
                region=region[region$chrom %in% chromosomes,]
                region=unlist(region$chrom)
            }
        }

            ###  ONLY RUN LOCALLY PARALLEL JOBS AT REGION LEVEL IF SAMPLE LEVEL LOCAL/BATCH MODE
            ###  LIMITATION OF MCLAPPLY

        run_mode="local_parallel"
        if(mode=="local_parallel"){
            run_mode="local"
        }

        .main.step$steps <-append(
        .main.step$steps,
        read_counter_region(
            bin_readcount=bin_readcount,
            bam=input,
            region=region,
            win=win,
            mapq=mapq,
            format=format,
            mode=run_mode,
            output_name=input_id,
            output_dir=paste0(out_file_dir,"/read_counts"),
            tmp_dir=tmp_dir,
            env_dir=env_dir,
            batch_dir=batch_dir,
            err_msg=err_msg,
            verbose=verbose,
            threads=threads,
            ram=ram,
            executor_id=task_id
            )
        )


        .this.step=.main.step$steps
        .main.step$out_files$region_counts=get_variable_env(env=.this.step)
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