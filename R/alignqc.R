

#' Generate quality control metrics for aligned sequences
#'
#' Comprehensive QC metrics generator for aligned BAM files supporting WGS, targeted sequencing (panel), and RNA-seq data.
#' Generates multiple metric types including MapQ distribution, insert size, artifact detection, and method-specific metrics.
#'
#' @description
#' This function generates comprehensive quality control metrics for an aligned sequence:
#' - **WGS metrics**: Generated when no target intervals are specified
#' - **Panel/Targeted metrics**: Generated when bait and target intervals are provided
#' - **RNA-seq metrics**: Generated when ribosomal intervals and reference flat files are specified
#'
#' Target and bait BEDs can be converted to interval format using Picard's BedToIntervalList:
#' `java -jar picard.jar BedToIntervalList I=targets.bed O=targets.interval_list SD=reference.fa`
#'
#' Off-target regions can be generated using bedtools complement:
#' `bedtools complement -i targets.bed -g reference.fa.fai > offtarget_regions.bed`
#'
#' @param bam Path to the input BAM file.
#' @param bin_samtools Path to samtools executable. Default: tools/samtools/samtools.
#' @param bin_picard Path to picard jar executable. Default: tools/picard/build/libs/picard.jar.
#' @param bin_bedtools Path to bedtools executable. Default: tools/bedtools2/bin/bedtools. Required for panel analysis.
#' @param ref_genome Path to reference genome FASTA file.
#' @param bi Bait capture target intervals (Picard interval format). For panel/targeted data.
#' @param ti Primary target intervals (Picard interval format). For panel/targeted data.
#' @param ri Ribosomal intervals file. Required for RNA-seq analysis.
#' @param ref_flat Reference gene model in flat format. Required for RNA-seq analysis.
#' @param mapq Minimum mapping quality threshold for metrics. Default: 0.
#' @param method Type of analysis method. Default: "tg". Options: ["wgs", "tg", "rna"].
#' @param ... Additional parameters passed to launch function, including:
#'   \itemize{
#'     \item{output_dir}{Path to output directory. Default: "."}
#'     \item{verbose}{Enable progress messages. Default: FALSE}
#'     \item{batch_config}{Batch configuration settings. Default: build_default_preprocess_config()}
#'     \item{threads}{Number of threads for parallel processing. Default: 3}
#'     \item{ram}{RAM in GB per thread. Default: 4}
#'     \item{executor_id}{Job executor ID. Default: unique ID with "alignQC" prefix}
#'     \item{task_name}{Task name for identification. Default: "alignQC"}
#'     \item{mode}{Execution mode. Options: ["local", "batch"]. Default: "local"}
#'     \item{time}{Max runtime for batch jobs. Default: "48:0:0"}
#'     \item{update_time}{Job status update interval in seconds. Default: 60}
#'     \item{wait}{Wait for batch job completion. Default: FALSE}
#'     \item{hold}{Job ID to hold until completion}
#'   }
#'
#' @return A job report list containing:
#'   \itemize{
#'     \item{job_id}{ID of the submitted job}
#'     \item{steps}{List of metric results for each QC step}
#'     \item{out_files}{Output file paths for generated metrics}
#'   }
#'
#' @details
#' The function executes the following QC steps in sequence:
#' 1. **mapq_qc**: Distribution of mapping quality scores
#' 2. **summary_qc**: Summary statistics using Picard CollectWgsMetrics
#' 3. **insert_size**: Fragment/insert size distribution
#' 4. **artifact_metrics**: Detection of sequencing artifacts (OxoG, FFPE)
#' 5. **oxog_metrics**: Oxidative guanine metrics
#' 6. **method-specific**: Additional metrics based on method (WGS/Panel/RNA)
#'
#' Note: BAM file must be sorted and indexed. Reference genome and intervals must use consistent coordinate systems.
#'
#' @export

new_metrics_alignqc=function(
  bin_samtools=build_default_tool_binary_list()$bin_samtools,
  bin_picard=build_default_tool_binary_list()$bin_picard,
  bin_bedtools=build_default_tool_binary_list()$bin_bedtools,
  ref_genome=build_default_reference_list()$HG19$reference$genome,
  bi=build_default_reference_list()$HG19$panel$PCF_V3$intervals$bi,
  ti=build_default_reference_list()$HG19$panel$PCF_V3$intervals$ti,
  bam=NULL,
  mapq=0,
  method="tg",
  ...
  ){
    # Main execution function that runs the QC pipeline
    # Defines the nested function that will be called by launch()
    
    run_main=function(
              .env
          ){
              # Initialize this environment and copy variables from parent
              .this.env=environment()
              append_env(to=.this.env,from=.env)
              # Initialize main job structure
              set_main(.env=.this.env)
            
            .main$steps[[fn_id]]<-.this.env
            .main.step=.main$steps[[fn_id]]

             # Record pipeline start time for elapsed time tracking
            start_time <- Sys.time()
            
      
            steps=c(
                    "mapq_qc",        # Step 1: Calculate MapQ distribution
                    "summary_qc",              # Step 2: Collect WGS summary metrics
                    "insert_size",         # Step 3: Analyze insert size distribution
                    "artifact",               # Step 4: Detect sequencing artifacts
                    "oxog"             # Step 5: Calculate OxoG metrics
                )

            # Add method-specific metrics based on analysis type
            if(!is.null(method)){
                    if(method=="targeted"){
                        steps=append(steps,"tg_qc")  # Add targeted/panel metrics
                    }else if (method=="wgs"){
                        steps=append(steps,"wgs_qc")  # Add whole genome metrics
                    }
            }
                
                # Total number of steps for progress tracking
                total_steps=length(steps)


            # Execute each QC step sequentially
            for(step in 1:total_steps){
        

                 # Log pipeline progress with current step number and name
                logger(paste("Running step", step, "of", total_steps, ":", steps[step]),start_time)

                tryCatch({
                    ### STEP 1
                    if(steps[step]=="mapq_qc"){
                                    .main.step$steps <-append(
                                    .main.step$steps,
                                    new_mapq_metrics_bam_samtools(
                                        bin_samtools=bin_samtools,
                                        bam=input,
                                        output_dir=paste0(out_file_dir,"/mapq"),
                                        output_name=input_id,
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

                    .this.step=.main.step$steps$new_mapq_metrics_bam_samtool
                    .main.step$out_files=append(.main.step$out_files,.this.step$out_files)


                        }

                    ### STEP 2

                    if(steps[step]=="summary_qc"){

                                    .main.step$steps <-append(
                                    .main.step$steps,
                                    new_summary_metrics_bam_picard(
                                        bin_samtools=bin_picard,
                                        bam=input,
                                        output_dir=paste0(out_file_dir,"/summary"),
                                        output_name=input_id,
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

                        .this.step=.main.step$steps$new_summary_metrics_bam_picard
                        .main.step$out_files=append(.main.step$out_files,.this.step$out_files)
                    }


                    ### STEP 3


                    if(steps[step]=="insert_size"){
                                    .main.step$steps <-append(
                                    .main.step$steps,
                                    new_insertsize_metrics_bam_picard(
                                        bin_picard=bin_picard,
                                        bam=input,
                                        output_dir=paste0(out_file_dir,"/insert_size"),
                                        output_name=input_id,
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

                        .this.step=.main.step$steps$new_insertsize_metrics_bam_picard
                        .main.step$out_files=append(.main.step$out_files,.this.step$out_files)
                    }

                    ### STEP 4
                    if(steps[step]=="artifact"){
                                    .main.step$steps <-append(
                                    .main.step$steps,
                                    artifact_metrics_bam_picard(
                                        bin_picard=bin_picard,
                                        bam=input,
                                        output_dir=paste0(out_file_dir,"/artifact"),
                                        output_name=input_id,
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
                        .this.step=.main.step$steps$artifact_metrics_bam_picard
                        .main.step$out_files=append(.main.step$out_files,.this.step$out_files)
                    }

                    ### STEP 5
                    if(steps[step]=="oxog"){
                                    .main.step$steps <-append(
                                    .main.step$steps,
                                    oxog_metrics_bam_picard(
                                        bin_picard=bin_picard,
                                        bam=input,
                                        output_dir=paste0(out_file_dir,"/artifact"),
                                        output_name=input_id,
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
                        .this.step=.main.step$steps$oxog_metrics_bam_picard
                        .main.step$out_files=append(.main.step$out_files,.this.step$out_files)
                    }


                    ### STEP 6
                    if(steps[step]=="tg_qc"){
                            .main.step$steps <-append(
                                    .main.step$steps,
                                    new_tg_summary_metrics_bam_picard(
                                        bin_picard=bin_picard,
                                        bam=input,
                                        bi=bi,
                                        ti=ti,
                                        output_dir=paste0(out_file_dir,"/tg_metrics"),
                                        output_name=input_id,
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
                        .this.step=.main.step$steps$new_tg_summary_metrics_bam_picard
                        .main.step$out_files=append(.main.step$out_files,.this.step$out_files)
                    }


                    ### STEP 7
                    if(steps[step]=="wgs_qc"){
                            .main.step$steps <-append(
                                    .main.step$steps,
                                    new_wgs_summary_metrics_bam_picard(
                                        bin_picard=bin_picard,
                                        bam=input,
                                        output_dir=paste0(out_file_dir,"/wgs_metrics"),
                                        output_name=input_id,
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
                        
                        .this.step=.main.step$steps$new_wgs_summary_metrics_bam_picard
                        .main.step$out_files=append(.main.step$out_files,.this.step$out_files)
                    }


                    # Log successful step completion
                    logger(paste("Completed step", step, "of", total_steps, ":", steps[step]),start_time)

                     }, error=function(e){
                    # Handle step execution errors with informative message
                    logger(paste("ERROR in step", step, ":", steps[step]),start_time)
                    stop(paste("Step '" , steps[step], "' failed. Error:", e$message,
                              "\nReview input files and parameters before retrying."))
                        })
                }

                                ### Remove input after processing
                if(clean){
                    file.remove(c(input,paste0(input,".bai")))
                }



                 # Log pipeline completion with total runtime
                total_elapsed <- as.numeric(difftime(Sys.time(), start_time, units="secs"))
                total_elapsed_str <- sprintf("%.1f", total_elapsed)
                logger(paste("AlignQC processing pipeline completed successfully."),start_time)
                logger(paste("Total steps executed:", total_steps, "| Total runtime:", total_elapsed_str, "seconds"),start_time)
                
                # Return main job structure to parent environment
                .env$.main <- .main
          }

    # Setup base environment and merge parameters passed via ...
    .base.env=environment()
    list2env(list(...),envir=.base.env)  # Add additional parameters
    # Configure environment variables for batch processing
    set_env_vars(
        .env= .base.env,
        vars="bam"
    )

    # Launch the pipeline with prepared environment
    launch(.env=.base.env)

    
  }

