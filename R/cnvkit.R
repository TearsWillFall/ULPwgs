


#'
#' This function functions calls Mutect2 for variant calling.
#' If a vector of tumour samples are provided these will be processed in multi-sample mode.
#' To run in tumour-normal mode suppply a single tumour and normal sample.
#' If no normal is supplied this will run in tumour only.
#' TO DO// Implement mitochondrial mode feature
#' 
#' For more information read:
#' https://cnvkit.readthedocs.io/en/stable/pipeline.html
#'
#'
#' @param sif_cnvkit [REQUIRED] Path to cnvkit sif file.
#' @param tumour [REQUIRED] Path to tumour BAM file.
#' @param normal [OPTIONAL] Path to normal BAM file.
#' @param ref_genome [REQUIRED] Path to reference genome fasta file.
#' @param rdata [OPTIONAL] Import R data information with list of BAM.
#' @param selected [OPTIONAL] Select BAM from list.
#' @param pon [OPTIONAL] Path to Panel of Normals. Default none
#' @param access [OPTIONAL] Path to reference genome accessibility. Default none
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
#' @param pon [OPTIONAL] Path to panel of normal.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param seg_method [OPTIONAL] Method to use to generate segmentation calls.Default  cbs. Options ["cbs","flasso","haar","none","hmm","hmm-tumor","hmm-germline"]
#' @param seq_method [OPTIONAL] Sequenced methods used. Default hybrid. Options ["hybrid","amplicon","wgs"]
#' @param trend_scatter [OPTIONAL] Show trendline. Default TRUE
#' @param range_scatter [OPTIONAL] Range to show in chr:start-end format. Default none
#' @param range_list_scatter [OPTIONAL] Path to BED file with range list. Default none
#' @param genes_scatter [OPTIONAL] List of genes to focus on. Default none.
#' @param margin_width_scatter [OPTIONAL] Default gene margin width. Default 1000000
#' @param plot_by_bin_scatter [OPTIONAL] Assume all bins are the same size. Default FALSE
#' @param trend_scatter [OPTIONAL] Show trendline. Default TRUE.
#' @param antitarget_symbol_scatter [OPTIONAL] Show antitarget probes with this symbol. Default "@"
#' @param segment_colour_scatter [OPTIONAL] Show antitarget probes with this symbol. Default "@"
#' @param y_max_scatter [OPTIONAL] Y max range. Default NULL.
#' @param y_min_scatter [OPTIONAL] Y min range. Default NULL.
#' @param cn_thr_diagram [OPTIONAL] Copy number threshold. Default 0.5
#' @param min_probes_diagram [OPTIONAL] Minimum number of probes to show a gene. Default 3
#' @param male_reference_diagram [OPTIONAL]  Assume inputs were normalized to a male reference. Default FALSE
#' @param gender_diagram [OPTIONAL] Assume sample gender. Default male
#' @param shift_diagram [OPTIONAL] Adjust the X and Y chromosomes according to sample sex. Default TRUE
#' @param normal_ichor [OPTIONAL] Path to normal samples. Default c(0.5,0.6,0.7,0.8,0.9)
#' @param ploidy_ichor [OPTIONAL] Initial tumour ploidy; can be more than one value if additional ploidy initializations are desired. Default: 2
#' @param lambda_ichor [OPTIONAL] Initial Student's t precision; must contain 4 values (e.g. c(1500,1500,1500,1500)); if not provided then will automatically use based on variance of data
#' @param scStates_ichor [OPTIONAL]  Subclonal states to consider. Default NULL
#' @param output_name_ichor [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
#' @param lambdaScaleHyperParam_ichor [OPTIONAL] Hyperparameter (scale) for Gamma prior on Student's-t precision. Default 3
#' @param maxCN_ichor [OPTIONAL] Total clonal states. Default 7.
#' @param estimateNormal_ichor [OPTIONAL] Estimate normal. Default TRUE.
#' @param estimateScPrevalence_ichor [OPTIONAL] Estimate subclonal prevalence. Default TRUE.
#' @param estimatePloidy_ichor [OPTIONAL] Estimate tumour ploidy. Default TRUE.
#' @param maxFracGenomeSubclone_ichor [OPTIONAL] Exclude solutions with subclonal genome fraction greater than this value. Default 0.5
#' @param maxFracCNASubclone_ichor  [OPTIONAL] Exclude solutions with fraction of subclonal events greater than this value. Default 0.7
#' @param minSegmentBins_ichor [OPTIONAL]  Minimum number of bins for largest segment threshold required to estimate tumor fraction; if below this threshold, then will be assigned zero tumor fraction
#' @param altFracThreshold_ichor [OPTIONAL] Minimum proportion of bins altered required to estimate tumor fraction; if below this threshold, then will be assigned zero tumor fraction. Default: [0.05]
#' @param includeHOMD_ichor [OPTIONAL] If FALSE, then exclude HOMD state. Useful when using large bins (e.g. 1Mb). Default FALSE.
#' @param txnE_ichor [OPTIONAL] Self-transition probability. Increase to decrease number of segments. Default: [0.9999999]
#' @param txnStrength_ichor [OPTIONAL] Transition pseudo-counts. Exponent should be the same as the number of decimal places of --txnE. Default: [1e+07]
#' @param plotFileType_ichor [OPTIONAL] File format for output plots. Default pdf
#' @param plotYLim_ichor [OPTIONAL] ylim to use for chromosome plots. Default: [c(-2,2)]
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


process_cnvkit=function(
    target=build_default_reference_list()$HG19$panel$PCF_V3$binned$target,
    antitarget=build_default_reference_list()$HG19$panel$PCF_V3$binned$antitarget,
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    access=build_default_reference_list()$HG19$reference$access_5k,
    sif_cnvkit=build_default_sif_list()$sif_cnvkit,
    pon=build_default_reference_list()$HG19$panel$PCF_V3$variant$pon_cn_male,
    tumour=NULL,
    seg_method="cbs",
    seq_method="hybrid",
    diagram=TRUE,
    scatter=TRUE,
    read_count=FALSE,
    min_mapq=0,
    gc=TRUE,
    edge=TRUE,
    rmask=TRUE,
    smooth=FALSE,
    drop_low_coverage=TRUE,
    drop_outliers=10,
    range_scatter=NULL,
    range_list_scatter=NULL,
    genes_scatter=NULL,
    margin_width_scatter=1000000,
    plot_by_bin_scatter=FALSE,
    trend_scatter=FALSE,
    antitarget_symbol_scatter="@",
    segment_colour_scatter="red",
    y_max_scatter=NULL,
    y_min_scatter=NULL,
    cn_thr_diagram=0.5,
    min_probes_diagram=3,
    male_reference_diagram=FALSE,
    gender_diagram="male",
    shift_diagram=TRUE,
    ...
  ){

  
    run_main=function(
      .env
    ){


      .this.env=environment()
      append_env(to=.this.env,from=.env)
      set_main(.env=.this.env)

      .main$steps[[fn]]<-.this.env
      .main.step=.main$steps[[fn]]

      .main.step$steps <-append(
        .main.step$steps,
        call_coverage_cnvkit(
          sif_cnvkit=sif_cnvkit,
          ref_genome=ref_genome,
          bam=input,
          bed=target,
          type="target",
          output_dir=out_file_dir,
          tmp_dir=tmp_dir,
          env_dir=env_dir,
          batch_dir=batch_dir,
          err_msg=err_msg,
          read_count=read_count,
          min_mapq=min_mapq,
          verbose=verbose,
          threads=threads,
          ram=ram,
          executor_id=task_id
      )
    )

      .this.step=.main.step$call_coverage_cnvkit.target
      .main.step$out_files$cnn=.this.step$out_files

      .main.step$steps <-append(
          .main.step$steps,
          call_coverage_cnvkit(
            sif_cnvkit=sif_cnvkit,
            ref_genome=ref_genome,
            bam=input,
            bed=target,
            type="antitarget",
            output_dir=out_file_dir,
            tmp_dir=tmp_dir,
            env_dir=env_dir,
            batch_dir=batch_dir,
            err_msg=err_msg,
            read_count=read_count,
            min_mapq=min_mapq,
            verbose=verbose,
            threads=threads,
            ram=ram,
            executor_id=task_id
        )
      )

      .this.step=.main.step$coverage_cnvkit.antitarget
      .main.step$out_files$cnn=append(.main.step$out_files$cnn,.this.step$out_files)


       .main.step$steps <-append(
          .main.step$steps,
          fix_cnvkit(
            sif_cnvkit=sif_cnvkit,
            pon=pon,
            target=.main.step$out_files$cnn$target,
            antitarget=.main.step$out_files$cnn$antitarget,
            output_name=input_id,
            output_dir=out_file_dir,
            tmp_dir=tmp_dir,
            env_dir=env_dir,
            batch_dir=batch_dir,
            verbose=verbose,
            err_msg=err_msg,
            gc=gc,
            edge=edge,
            rmask=rmask,
            threads=threads,
            ram=ram,
            executor_id=task_id
          )
       )

      .this.step=.main.step$fix_cnvkit
      .main.step$out_files=append(.main.step$out_files,.this.step$out_files)
          
      .main.step$steps <-append(
          .main.step$steps,
          segment_cnvkit(
            sif_cnvkit=sif_cnvkit,
            cnr=.main.step$out_files$cnr,
            seg_method=seg_method,
            output_name=input_id,
            output_dir=out_file_dir,
            tmp_dir=tmp_dir,
            env_dir=env_dir,
            batch_dir=batch_dir,
            err_msg=err_msg,
            smooth=smooth,
            drop_low_coverage=drop_low_coverage,
            drop_outliers=drop_outliers,
            verbose=verbose,
            threads=threads,
            ram=ram,
            executor_id=task_id
          )
        )

      .this.step=.main.step$segment_cnvkit
      .main.step$out_files=append(.main.step$out_files,.this.step$out_files)

    if(scatter){
          .main.step$steps <-append(
              .main.step$steps,
              scatter_cnvkit(
                sif_cnvkit=sif_cnvkit,
                cnr=.main.step$out_files$cnr,
                cns=.main.step$out_files$cns,
                output_name=input_id,
                output_dir=out_file_dir,
                tmp_dir=tmp_dir,
                env_dir=env_dir,
                batch_dir=batch_dir,
                err_msg=err_msg,
                range=range_scatter,
                range_list=range_list_scatter,
                genes=genes_scatter,
                margin_width=margin_width_scatter,
                plot_by_bin=plot_by_bin_scatter,
                trend=trend_scatter,
                antitarget_symbol=antitarget_symbol_scatter,
                segment_colour=segment_colour_scatter,
                title=input_id,
                y_max=y_max_scatter,
                y_min=y_max_scatter,
                verbose=verbose,
                threads=threads,ram=ram,
                executor_id=task_id
           )
        )
          .this.step=.main.step$scatter_cnvkit
          .main.step$out_files=append(.main.step$out_files,.this.step$out_files)



      }

    if(diagram){
      .main.step$steps <-append(
        .main.step$steps,
          diagram_cnvkit(
            sif_cnvkit=sif_cnvkit,
            cnr=jobs_report[["steps"]][["fixCNVkit"]]$out_files$cnr,
            cns=jobs_report[["steps"]][["segmentCNVkit"]]$out_files$cns,
            output_name=input_id,
            output_dir=out_file_dir,
            tmp_dir=tmp_dir,
            env_dir=env_dir,
            batch_dir=batch_dir,
            err_msg=err_msg,
            cn_thr=cn_thr_diagram,
            min_probes=min_probes_diagram,
            male_reference=male_reference_diagram,
            gender=gender_diagram,
            shift=shift_diagram,
            title=input_id,
            verbose=verbose,
            threads=threads,
            ram=ram,
            executor_id=task_id
        )
      )

        .this.step=.main.step$scatter_cnvkit
        .main.step$out_files=append(.main.step$out_files,.this.step$out_files)
    }
    

    .env.$.main<-.main

    }
 
 
  .base.env=environment()
  list2env(list(...),envir=.base.env)
  set_env_vars(
    .env= .base.env,
    vars="bam"
  )

  launch(.env=.base.env)

 
}












#' Wrapper around access command in CNVkit
#'
#' Create a BED file with non-mappable regions to exclude from the genome.
#' Additional regions can be supplied in a BED format using the exclude_regions argument
#' 
#' For more information read:
#' https://cnvkit.readthedocs.io/en/stable/pipeline.html
#'
#' @param sif_cnvkit [REQUIRED] Path to cnvkit sif file.
#' @param ref_genome [REQUIRED] Path to reference genome fasta file.
#' @param exclude_regions [REQUIRED] Additional regions to exclude in BED format. Default none.
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
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


access_cnvkit=function(
    sif_cnvkit=build_default_sif_list()$sif_cnvkit,
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    exclude_regions=NULL,
    gap_size=5000,
    ...
  ){


    run_main=function(
      .env
    ){


      .this.env=environment()
      append_env(to=.this.env,from=.env)
      set_main(.env=.this.env)

      .main$out_files$access=paste0(out_file_dir,"/",input_id,".",gap_size,".access.bed")

      if(exclude_regions!=""){
        exclude_regions=paste0(" -x ",exclude_regions)
      }


      .main$exec_code=paste("singularity exec -H ",paste0(getwd(),":/home "),sif_cnvkit,
        " cnvkit.py access -o ",.main$out_files$access,
        exclude_regions," -s ",gap_size, ref_genome
      )

      run_job(.env=.this.env)

      .env$.main<-.main

    }
   
      .base.env=environment()
      list2env(list(...),envir=.base.env)
      set_env_vars(
        .env= .base.env,
        vars="ref_genome"
      )

      launch(.env=.base.env)


  }









  

#' Wrapper around target function from CNVkit
#'
#' This function wraps around target function for CNVkit
#' This function generates an target BED file from an input target file. 
#' Additional parameters can be used to exclude regions and modify the average bin size
#' 
#' 
#' For more information read:
#' https://cnvkit.readthedocs.io/en/stable/pipeline.html
#' 
#' 
#' @param sif_cnvkit [REQUIRED] Path to cnvkit sif file.
#' @param ref_genome [REQUIRED] Path to reference genome.
#' @param bin_size_target [OPTIONAL] Size of bins for targets. Default 75
#' @param bin_size_antitarget [OPTIONAL] Size of bins for antitargets. Default 500000
#' @param min_bin_size_antitarget [OPTIONAL] Size of bins for antitargets. Default NULL
#' @param access [OPTIONAL] Path to reference genome accessibility. Default none
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
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


 create_pon_cnvkit=function(
    sif_cnvkit=build_default_sif_list()$sif_cnvkit,
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    access=build_default_reference_list()$HG19$reference$access_5k,
    normal=NULL,
    target=NULL,
    bin_size_target=75,
    bin_size_antitarget=500000,
    min_bin_size_antitarget=NULL,
    ...
  ){


    run_main=function(
      .env
    ){

        .this.env=environment()
        append_env(to=.this.env,from=.env)
        set_main(.env=.this.env)

        .main$out_files$pon=paste0(out_file_dir,"/",
          input_id,
          ".target_",bin_size_target,
          ".antitarget_",bin_size_antitarget,
          ".pon.cnn"
        )


        if(!is.null(target)){
          target=paste0(" -t ",target)
        }


        if(!is.null(access)){
          paste0(" -g ",access)
        }

        if(!is.null(min_bin_size_antitarget)){
          min_bin_size_antitarget=paste0(
            " --antitarget-min-size ",
            min_bin_size_antitarget
            )
        }

        .main$exec_code=paste("singularity exec -H /:/home ",sif_cnvkit,
          " cnvkit.py batch -p ",threads, " -n ",gsub(";"," ",input),
          " --output-reference ",.main$out_files$pon,
          " -f ", ref_genome,access,target,
          " --target-avg-size ",bin_size_target,
          " --antitarget-avg-size ",bin_size_antitarget,
          min_bin_size_antitarget
        )

        run_job(.env=.this.env)
        .env$.main<-.main
    }
    
      .base.env=environment()
      list2env(list(...),envir=.base.env)
      set_env_vars(
        .env= .base.env,
        vars="normal"
      )

      launch(.env=.base.env)


  }





#' Wrapper around reference function from CNVkit
#'
#' This function wraps around target function for CNVkit
#' This function generates an target BED file from an input target file. 
#' Additional parameters can be used to exclude regions and modify the average bin size
#' 
#' 
#' For more information read:
#' https://cnvkit.readthedocs.io/en/stable/pipeline.html
#'
#' @param sif_cnvkit [REQUIRED] Path to cnvkit sif file.
#' @param ref_genoma [REQUIRED] Path to reference genoms.
#' @param cnn [REQUIRED] Path to normal coverage profiles. Default none
#' @param gender [OPTIONAL] Sample gender. Default male.
#' @param target [OPTIONAL] Path to BED file with target regions. Default none.
#' @param antitarget [OPTIONAL] Path to BED file with target regions. Default none.
#' @param cluster [OPTIONAL] Calculate and store summary stats for clustered subsets of the normal samples with similar coverage profiles. Default FALSE
#' @param min_cluster_size [OPTIONAL] Minimum size to keep in reference profiles.
#' @param gc [OPTIONAL] Disble GC bias correction. Default FALSE.
#' @param edge [OPTIONAL] Disble edge correction. Default FALSE.
#' @param rmask [OPTIONAL] Disble repeat mask correction. Default FALSE.
#' @param male_reference [OPTIONAL] Adjust X chromosome to log2 0. Default FALSE.
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
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


  reference_cnvkit=function(
    sif_cnvkit=build_default_sif_list()$sif_cnvkit,
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    cnn=NULL,
    gender="male",
    gc=TRUE,
    edge=TRUE,
    rmask=TRUE,
    verbose=FALSE,
    cluster=FALSE,
    min_cluster_size=10,
    male_reference=FALSE,
    threads=1,ram=1,
    ...
  ){

  run_main=function(
    .env
  ){
    if(!is.null(gender)){
      gender=paste0(" -x ",gender)
    }

    add=""

    if(!gc){
      add=paste(add," --no-gc ")
    }

    if(!edge){
      add=paste(add," --no-edge ")
    }

    if(!rmask){
      add=paste(add," --no-rmask ")
    }


    if(male_reference){
      add=paste(add," -y ")
    }

    if(cluster){
      add=paste(add," -c ")
      min_cluster_size=paste0(" --min-cluster-size ",min_cluster_size) 
    }else{
      min_cluster_size=""
    }
  


    .main$out_file$reference=paste0(out_file_dir,"/",input_id,".cnn")

    .main$exec_code=paste("singularity exec -H /:/home ",sif_cnvkit,
      " cnvkit.py reference ", gender, min_cluster_size,
      " -f ",normalizePath(ref_genome)," -o ",
      .main$out_file$pon, add, gsub(";"," ",cnn)
    )

    run_job(.env=.this.env)
    .env$.main<-.main

  }



  .base.env=environment()
  list2env(list(...),envir=.base.env)
  set_env_vars(
    .env= .base.env,
    vars="cnn"
  )

  launch(.env=.base.env)

    
}







  

#' Wrapper around target function from CNVkit
#'
#' This function wraps around target function for CNVkit
#' This function generates an target BED file from an input target file. 
#' Additional parameters can be used to exclude regions and modify the average bin size
#' 
#' 
#' For more information read:
#' https://cnvkit.readthedocs.io/en/stable/pipeline.html
#'
#' @param sif_cnvkit [REQUIRED] Path to cnvkit sif file.
#' @param bed [REQUIRED] Path to input BED file with target regions. Default none
#' @param annotation [OPTIONAL] BED file with region annotation information. Default none
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param short_names [OPTIONAL] Use short annotation names. Default FALSE
#' @param bin_size [OPTIONAL] Average bin size. Default 100.
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


 create_target_cnvkit=function(
    sif_cnvkit=build_default_sif_list()$sif_cnvkit,
    bed=NULL,
    annotation=NULL,
    split=FALSE,
    short_names=FALSE,
    bin_size=100,
    ...
  ){
    

    run_main=function(
      .env
    ){

      .this.env=environment()
      append_env(to=.this.env,from=.env)
      set_main(.env=.this.env)

      
      if(!is.null(annotation)){
        annotation=paste0(" --annotate ",annotation)
      }

      add=""
      if(short_names){
        add=paste(add," --short-names ")
      }

      .main$out_files$tg_bed=paste0(
        out_file_dir,"/",input_id,".binned.targets.bed"
      )

      
      .main$exec_code=paste(
        "singularity exec -H /:/home ",sif_cnvkit,
        " cnvkit.py target --split -a ",bin_size," -o ",
        .main$out_files$tg_bed, add, input
      )

      run_job(.env=.this.env)


      .env$.main<-.main

    }


    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="bed"
    )

    launch(.env=.base.env)



  }




#' Wrapper around antitarget function from CNVkit
#'
#' This function wraps around antitarget function for CNVkit
#' This function generates an antitarget BED file from an input target file. 
#' Additional parameters can be used to exclude regions and modify the average bin size
#' 
#' 
#' For more information read:
#' https://cnvkit.readthedocs.io/en/stable/pipeline.html
#'
#' @param sif_cnvkit [REQUIRED] Path to cnvkit sif file.
#' @param bed [REQUIRED] Path to input BED file with target regions. Default none
#' @param access [OPTIONAL] Path to non-accessible regions to exclude. Default none
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param bin_size [OPTIONAL] Average bin size. Default 100.
#' @param min_bin_size [OPTIONAL] Minimum average bin size. Bins smaller that this will be dropped. Default NULL.
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

  create_antitarget_cnvkit=function(
    sif_cnvkit=build_default_sif_list()$sif_cnvkit,
    access=build_default_reference_list()$HG19$reference$access_5k,
    bed=NULL,
    bin_size=500000,
    min_bin_size=NULL,
    ...
  ){


    run_main=function(
      .env
    ){
      
        .this.env=environment()
        append_env(to=.this.env,from=.env)
        set_main(.env=.this.env)


          
        if(!is.null(min_bin_size)){
          add=paste0(" -m ",min_bin_size)
        }

        .main$out_files$atg_bed=paste0(
          out_file_dir,"/",input_id,".binned.antitargets.bed"
        )

       
        if(!is.null(access)){
          access=paste0(" -g ",access)
        }



        .main$exec_code=paste(
          "singularity exec -H   /:/home ",sif_cnvkit,
          " cnvkit.py antitarget -a ",
          bin_size,access,
          " -o ",.main$out_files$atg_bed, add, normalizePath(input)
        )


      run_job(.env=.this.env)
    
      .env$.main<-.main

    }

    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="bed"
    )

    launch(.env=.base.env)

   

  }









#' Create binned target and anitarget beds from target BED file
#'
#' This function wraps around target and antitarget functiosn for CNVkit
#' This function generates an target BED file from an input target file. 
#' Additional parameters can be used to exclude regions and modify the average bin size
#' 
#' 
#' For more information read:
#' https://cnvkit.readthedocs.io/en/stable/pipeline.html
#'
#' @param sif_cnvkit [REQUIRED] Path to cnvkit sif file.
#' @param access [OPTIONAL] Path to non-accessible regions to exclude. Default none
#' @param bed [REQUIRED] Path to input BED file with target regions. Default none
#' @param annotation [OPTIONAL] BED file with region annotation information. Default none
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param short_names [OPTIONAL] Use short annotation names. Default FALSE
#' @param bin_size_target [OPTIONAL] Average bin size for target. Default 120.
#' @param bin_size_antitarget [OPTIONAL] Average bin size for target. Default 500000.
#' @param min_bin_size_antitarget [OPTIONAL] Average bin size for target. Default null.
#' @param split [OPTIONAL] Average bin size for target. Default null.
#' @param output_dir [OPTIONAL] Split annotation names.
#' @param output_dir [OPTIONAL] Split annotation names.
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



bin_targets_cnvkit=function(
    sif_cnvkit=build_default_sif_list()$sif_cnvkit,
    access=build_default_reference_list()$HG19$reference$access_5k,
    bed=NULL,
    annotation=NULL,
    bin_size_antitarget=100000,
    bin_size_target=100,
    min_bin_size_antitarget=NULL,
    split=FALSE,
    short_names=FALSE,
    ...
){

   run_main=function(
    .env
   ){

      .this.env=environment()
      append_env(to=.this.env,from=.env)


      set_main(.env=.this.env)
    

      .main$steps[[fn]]<-.this.env
      .main.step=.main$steps[[fn]]
    

      .main.step$steps=append(
        .main.step$steps,
            create_target_cnvkit(
              sif_cnvkit=sif_cnvkit,
              bed=bed,
              annotation=annotation,
              output_name=output_name,
              output_dir=out_file_dir,
              tmp_dir=tmp_dir,
              batch_dir=batch_dir,
              env_dir=env_dir,
              split=split,
              short_names=short_names,
              bin_size=bin_size_target,
              verbose=verbose,
              threads=threads,ram=ram,
              executor_id=task_id
          )



      )

      .this.step=.main.step$steps$create_target_cnvkit
      .main.step$out_files<-.this.step$out_files

      .main.step$steps=append(
        .main.step$steps,
        create_antitarget_cnvkit(
          sif_cnvkit=sif_cnvkit,
          access=access,
          bed=bed,
          output_name=output_name,
          output_dir=out_file_dir,
          tmp_dir=tmp_dir,
          batch_dir=batch_dir,
          env_dir=env_dir,
          bin_size=bin_size_antitarget,
          min_bin_size=min_bin_size_antitarget,
          verbose=verbose,
          threads=threads,ram=ram,
          executor_id=task_id
        )

      )

      .this.step=.main.step$steps$create_antitarget_cnvkit
      .main.step$out_files<-append(
        .main.step$out_files,
        .this.step$out_files
      )

      .env$.main<-.main


   }
    
    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="bed"
    )
  
    launch(.env=.base.env)

}




  
#' Wrapper around autobin function from CNVkit
#'
#' This function wraps around target function for CNVkit
#' This function generates an target BED file from an input target file. 
#' Additional parameters can be used to exclude regions and modify the average bin size
#' 
#' 
#' For more information read:
#' https://cnvkit.readthedocs.io/en/stable/pipeline.html
#'
#' @param sif_cnvkit [REQUIRED] Path to cnvkit sif file.
#' @param ref_genoma [REQUIRED] Path to reference genoms.
#' @param bed [REQUIRED] Path to input BED file with target regions. Default none
#' @param bams [REQUIRED] Path to BAM files. Default none
#' @param annotation [OPTIONAL] BED file with region annotation information. Default none
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param short_names [OPTIONAL] Use short annotation names. Default FALSE
#' @param seq_method [OPTIONAL] Sequenced methods used. Default hybrid. Options ["hybrid","amplicon","wgs"]
#' @param seq_type [OPTIONAL] Use short annotation names. Default FALSE
#' @param min_bin_size_target [OPTIONAL] Mininimum target bin size. Default 20.
#' @param max_bin_size_target [OPTIONAL] Mininimum target bin size. Default 20000.
#' @param max_bin_size_antitarget [OPTIONAL] Mininimum target bin size. Default 500000.
#' @param min_bin_size_antitarget [OPTIONAL] Mininimum target bin size. Default 500.
#' @param bp_per_bin [OPTIONAL] Bases per bin. Default 100000.
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


autobin_cnvkit=function(
    sif_cnvkit=build_default_sif_list()$sif_cnvkit,
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    access=build_default_reference_list()$HG19$reference$access_5k,
    bam=NULL,
    bed=NULL,
    annotation=NULL,
    seq_method="hybrid",
    split=FALSE,
    short_names=FALSE,
    bp_per_bin=100000,
    min_bin_size_target=20,
    max_bin_size_target=20000,
    min_bin_size_antitarget=500,
    max_bin_size_antitarget=500000,
    ...
){


      run_main=function(
        .env
      ){

          .this.env=environment()
          append_env(to=.this.env,from=.env)
          set_main(.env=.this.env)
        
          if(!is.null(annotation)){
            annotation=paste0(" --annotate ",annotation)
          }

          if(!is.null(access)){
            access=paste0(" -g ",access)
          }

          add=""
          if(short_names){
            add=paste(add," --short-names ")
          }

          .main$exec_code=paste("singularity exec -H /:/home ",sif_cnvkit,
            " cnvkit.py autobin -t ",bed," -f ",ref_genome,
            " -b ",bp_per_bin, " -m ",seq_method,access,
            "  --target-max-size ",max_bin_size_target,
            "  --target-min-size ",min_bin_size_target,
            "  --antitarget-min-size ",min_bin_size_antitarget,
            "  --antitarget-max-size ",max_bin_size_antitarget,
            add,annotation,
            gsub(";"," ",input)
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







  #' Wrapper around autobin function from CNVkit
#'
#' This function wraps around target function for CNVkit
#' This function generates an target BED file from an input target file. 
#' Additional parameters can be used to exclude regions and modify the average bin size
#' 
#' 
#' For more information read:
#' https://cnvkit.readthedocs.io/en/stable/pipeline.html
#'
#' @param sif_cnvkit [REQUIRED] Path to cnvkit sif file.
#' @param ref_genoma [REQUIRED] Path to reference genoms.
#' @param bed [REQUIRED] Path to input BED file with target regions. Default none
#' @param bam [REQUIRED] Path to BAM files. Default none
#' @param read_count [OPTIONAL] Alternative method for coverage. Default FALSE
#' @param min_mapq [OPTIONAL] Minimum mapping quality to count a read for coverage. Default 0.
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
#' @param output_dir [OPTIONAL] Path to the output directory.
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


coverage_cnvkit=function(
  sif_cnvkit=build_default_sif_list()$sif_cnvkit,
  ref_genome=build_default_reference_list()$HG19$reference$genome,
  bam=NULL,
  bed=NULL,
  read_count=FALSE,
  min_mapq=0,
  ...
){
  
  run_main=function(
    .env
  ){

    .this.env=environment()
    append_env(to=.this.env,from=.env)
    set_main(.env=.this.env)

    .main$out_files$cnn=paste0(out_file_dir,"/",input_id,".cnn")

    add=""
    if(read_count){
      add=" -c "
    }

    .main$exec_code=paste(
      "singularity exec -H /:/home ",sif_cnvkit,
      " cnvkit.py coverage -p ",threads, "-q ",min_mapq,
      " -f ", normalizePath(ref_genome)," -o ",.main$out_files$cnn,
      add,normalizePath(input), bed
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





  #' Wrapper around autobin function from CNVkit
#'
#' This function wraps around target function for CNVkit
#' This function generates an target BED file from an input target file. 
#' Additional parameters can be used to exclude regions and modify the average bin size
#' 
#' 
#' For more information read:
#' https://cnvkit.readthedocs.io/en/stable/pipeline.html
#'
#' @param sif_cnvkit [REQUIRED] Path to cnvkit sif file.
#' @param ref_genoma [REQUIRED] Path to reference genoms.
#' @param bed [REQUIRED] Path to input BED file with target regions. Default none
#' @param bam [REQUIRED] Path to BAM files. Default none
#' @param read_count [OPTIONAL] Alternative method for coverage. Default FALSE
#' @param min_mapq [OPTIONAL] Minimum mapping quality to count a read for coverage. Default 0.
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
#' @param output_dir [OPTIONAL] Path to the output directory.
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


call_coverage_cnvkit=function(
  sif_cnvkit=build_default_sif_list()$sif_cnvkit,
  ref_genome=build_default_reference_list()$HG19$reference$genome,
  bam=NULL,
  type="target",
  bed=NULL,
  read_count=FALSE,
  min_mapq=0,
  ...
){
  
  run_main=function(
    .env
  ){

    .this.env=environment()
    append_env(to=.this.env,from=.env)
    set_main(.env=.this.env)


    output_name=paste0(input_id,".",type,"coverage")
    fn=paste0(fn,".",type)


    .main$steps[[fn]]<-.this.env
    .main.step=.main$steps[[fn]]
    
    if(type!="target" && type != "antitarget"){
      stop("Valid values for type are : target / antitarget")
    }


    .main.step=append(
      .main.step,
      coverage_cnvkit(
        sif_cnvkit=sif_cnvkit,
        ref_genome=ref_genome,
        bam=bam,
        bed=bed,
        output_name=output_name,
        output_dir=out_file_dir,
        tmp_dir=tmp_dir,
        batch_dir=batch_dir,
        env_dir=env_dir,
        err_msg=err_msg,
        read_count=read_count,
        min_mapq=min_mapq,
        verbose=verbose,
        threads=threads,
        ram=ram,
        executor_id=task_id
      )
    )

    .this.step=.main.step$steps$variants_by_filters_vcf
    .main.step$out_files[[type]]=.this.step$out_files$cnn


  }

  .base.env=environment()
  list2env(list(...),envir=.base.env)
  set_env_vars(
    .env= .base.env,
    vars="bam"
  )

  launch(.env=.base.env)
  

}




  #' Wrapper around autobin function from CNVkit
#'
#' This function wraps around target function for CNVkit
#' This function generates an target BED file from an input target file. 
#' Additional parameters can be used to exclude regions and modify the average bin size
#' 
#' 
#' For more information read:
#' https://cnvkit.readthedocs.io/en/stable/pipeline.html
#'
#' @param sif_cnvkit [REQUIRED] Path to cnvkit sif file.
#' @param ref_genoma [REQUIRED] Path to reference genoms.
#' @param bed [REQUIRED] Path to input BED file with target regions. Default none
#' @param bam [REQUIRED] Path to BAM files. Default none
#' @param read_count [OPTIONAL] Alternative method for coverage. Default FALSE
#' @param min_mapq [OPTIONAL] Minimum mapping quality to count a read for coverage. Default 0.
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
#' @param output_dir [OPTIONAL] Path to the output directory.
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


extract_coverage_cnvkit=function(
  sif_cnvkit=build_default_sif_list()$sif_cnvkit,
  ref_genome=build_default_reference_list()$HG19$reference$genome,
  bam=NULL,
  bed=NULL,
  type="target",
  read_count=FALSE,
  min_mapq=0,
  ...
){
  
  run_main=function(
    .env
  ){

    .this.env=environment()
    append_env(to=.this.env,from=.env)

    output_name=paste0(input_id,".",type)
    fn=paste0(fn,".",type)

    set_main(.env=.this.env)

          .main$steps[[fn]]<-.this.env
      .main.step=.main$steps[[fn]]

    
    .main.step$steps=append(
      .main.step$steps,
      create_antitarget_cnvkit(
        sif_cnvkit=sif_cnvkit,
        access=access,
        bed=bed,
        output_name=output_name,
        output_dir=out_file_dir,
        tmp_dir=tmp_dir,
        batch_dir=batch_dir,
        env_dir=env_dir,
        err_msg=err_msg,
        bin_size=bin_size_antitarget,
        min_bin_size=min_bin_size_antitarget,
        verbose=verbose,
        threads=threads,
        ram=ram,
        executor_id=task_id
      )

    )




  }

  .base.env=environment()
  list2env(list(...),envir=.base.env)
  set_env_vars(
    .env= .base.env,
    vars="bam"
  )

  launch(.env=.base.env)
  

}









#' Wrapper around reference function from CNVkit
#'
#' This function wraps around target function for CNVkit
#' This function generates an target BED file from an input target file. 
#' Additional parameters can be used to exclude regions and modify the average bin size
#' 
#' 
#' For more information read:
#' https://cnvkit.readthedocs.io/en/stable/pipeline.html
#'
#' @param sif_cnvkit [REQUIRED] Path to cnvkit sif file.oms.
#' @param seg_method [OPTIONAL] Method to use for segmentation. Default cbs. Options ["cbs","flasso","haar","none","hmm","hmm-tumor","hmm-germline"]
#' @param smooth [OPTIONAL] Smooth before CBS. Default TRUE
#' @param drop_low_coverage [OPTIONAL] Drop low coverage bins before segmentation. Default TRUE
#' @param drop_outliers [OPTIONAL] Drop outlier bins before segmentation. Default TRUE
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
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

segment_cnvkit=function(
    sif_cnvkit=build_default_sif_list()$sif_cnvkit,
    cnr=NULL,
    seg_method="cbs",
    smooth=FALSE,
    drop_low_coverage=FALSE,
    drop_outliers=10,
    ...
  ){

    

    run_main=function(
      .env
    ){


      .this.env=environment()
      append_env(to=.this.env,from=.env)
      set_main(.env=.this.env)


      add=""

      if(smooth){
        add=paste(add," --smooth-cbs ")
      }


      if(drop_low_coverage){
        add=paste(add," --drop-low-coverage ")
      }


      if(drop_outliers){
          drop_outliers=paste0("  --drop-outliers ",drop_outliers)
      }
    
      .main$out_files$cns=paste0(out_file_dir,"/",input_id,".cns")

      .main$exec_code=paste(
        "singularity exec -H /:/home ", sif_cnvkit,
        " cnvkit.py segment -p ",
        threads,drop_outliers," -o ",
        .main$out_files$cns,
        add,cnr
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









#' Wrapper around reference scatter from CNVkit
#'
#' This function wraps around target function for CNVkit
#' This function generates an target BED file from an input target file. 
#' Additional parameters can be used to exclude regions and modify the average bin size
#' 
#' 
#' For more information read:
#' https://cnvkit.readthedocs.io/en/stable/pipeline.html
#'
#' @param sif_cnvkit [REQUIRED] Path to cnvkit sif file.oms.
#' @param trend [OPTIONAL] Smooth before CBS. Default TRUE
#' @param range [OPTIONAL] Range to show in chr:start-end format. Default none
#' @param range_list [OPTIONAL] Path to BED file with range list. Default none
#' @param genes [OPTIONAL] List of genes to focus on. Default none.
#' @param margin_width [OPTIONAL] Default gene margin width. Default 1000000
#' @param plot_by_bin [OPTIONAL] Assume all bins are the same size. Default FALSE
#' @param trend [OPTIONAL] Show trendline. Default TRUE.
#' @param antitarget_symbol [OPTIONAL] Show antitarget probes with this symbol. Default "@"
#' @param segment_colour [OPTIONAL] Show antitarget probes with this symbol. Default "@"
#' @param title [OPTIONAL] Sample title. Default NULL.
#' @param y_max [OPTIONAL] Y max range. Default NULL.
#' @param y_min[OPTIONAL] Y min range. Default NULL.
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
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

  scatter_cnvkit=function(
    sif_cnvkit=build_default_sif_list()$sif_cnvkit,
    chrs=build_default_chr_list()$canonical,
    cnr=NULL,
    cns=NULL,
    range=NULL,
    range_list=NULL,
    genes=NULL,
    margin_width=1000000,
    plot_by_bin=FALSE,
    trend=TRUE,
    antitarget_symbol="@",
    segment_colour="red",
    title=NULL,
    y_max=NULL,
    y_min=NULL,
    ...
  ){




    run_main=function(
      .env
    ){

        
        .this.env=environment()
        append_env(to=.this.env,from=.env)
        set_main(.env=.this.env)


        if(!is.null(range)){
          range=paste0(" -c ",range)
        }

        if(!is.null(segment_colour)){
        segment_colour=paste0(" --segment-color ",segment_colour)
        }

        
        if(!is.null(genes)){
          genes=paste0(" --gene ",genes)
        }

        if(!is.null(range_list)){
          range_list=paste0(" -l ",range_list)
        }

        add=""

        if(!is.null(title)){
          title=paste0(" --title ",title)
        }
        if(!is.null(y_max)){
          y_max=paste0(" --y-max ",y_max)
        }


        if(!is.null(y_min)){
          y_min=paste0(" --y-min ",y_min)
        }
        

        if(trend){
          add=paste(add," --trend ")
        }


        if(plot_by_bin){
          add=paste(add," --by-bin ")
        }


        if(!is.null(margin_width)){
          margin_width=paste0(" --width ",margin_width)
        }


        .main$out_files$scatter=paste0(out_file_dir,"/",input_id,".scatter.pdf")


        cnr_tmp=paste0(cnr,".tmp")
        cns_tmp=""
        .main$exec_code=paste0("cat ",cnr," | head -n 1 > ",
          cnr_tmp," && cat ",cnr," | grep \"",
          paste0(paste0("^",chrs),collapse="\\|"),"\" >> ",cnr_tmp
        )
        if(!is.null(cns)){
          cns_tmp=paste0(cns,".tmp")
          cns_code=paste(" -s ",cns_tmp)
          .main$exec_code=paste0(
            .main$exec_code," && cat ",cns," |  head -n 1 > ",cns_tmp,
            " && ",paste0("cat ",cns," | grep \" ",
              paste0(paste0("^",unlist(chrs),collapse="\\|")," \" >> ",cns_tmp)
  
            )
          )
        }

        .main$exec_code=paste(
          .main$exec_code,";singularity exec -H /:/home ",
          sif_cnvkit," cnvkit.py scatter -o ",.main$out_files$scatter,
          cns_code,add,
          title,
          segment_colour,
          y_max,
          y_min,
          range,
          margin_width,
          range_list,
          cnr_tmp,
          "&& rm ",
          cnr_tmp,cns_tmp
        )

        .env$.main <- .main
      }

      .base.env=environment()
      list2env(list(...),envir=.base.env)
      set_env_vars(
          .env=.base.env,
          vars="cnr"
      )

      launch(.env=.base.env)

  }









#' Wrapper around reference function from CNVkit
#'
#' This function wraps around target function for CNVkit
#' This function generates an target BED file from an input target file. 
#' Additional parameters can be used to exclude regions and modify the average bin size
#' 
#' 
#' For more information read:
#' https://cnvkit.readthedocs.io/en/stable/pipeline.html
#'
#' @param sif_cnvkit [REQUIRED] Path to cnvkit sif file.oms.
#' @param target [REQUIRED] Path to CNN file for target regions. Default none.
#' @param antitarget [REQUIRED] Path to BED file with target regions. Default none.
#' @param gc [OPTIONAL] Disble GC bias correction. Default FALSE.
#' @param edge [OPTIONAL] Disble edge correction. Default FALSE.
#' @param rmask [OPTIONAL] Disble repeat mask correction. Default FALSE.
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
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


 fix_cnvkit=function(
    sif_cnvkit=build_default_sif_list()$sif_cnvkit,
    pon=build_default_reference_list()$HG19$panel$PCF_V3$variant$pon_cn_male,
    target=NULL,
    antitarget=NULL,
    output_name=NULL,
    gc=TRUE,
    edge=TRUE,
    rmask=TRUE,
    ...
  ){

    run_main=function(.env){

      .this.env=environment()
      append_env(to=.this.env,from=.env)
      set_main(.env=.this.env)

      .main$out_files$cnr=paste0(out_file_dir,"/",input_id,".cnr")

      add=""

      if(!gc){
        add=paste(add," --no-gc ")
      }

      if(!edge){
        add=paste(add," --no-edge ")
      }

      if(!rmask){
        add=paste(add," --no-rmask ")
      }


      .main$exec_code=paste(
        "singularity exec -H /:/home ",sif_cnvkit,
        " cnvkit.py fix -o ",.main$out_files$cnr, 
        add, target, antitarget, pon
      )

      run_job(.env=this.env)

      .env$.main<-.main

    }

    
    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
        .env=.base.env,
        vars="target"
    )

    launch(.env=.base.env)

    

   
  }
  


#' Wrapper around reference scatter from CNVkit
#'
#' This function wraps around target function for CNVkit
#' This function generates an target BED file from an input target file. 
#' Additional parameters can be used to exclude regions and modify the average bin size
#' 
#' 
#' For more information read:
#' https://cnvkit.readthedocs.io/en/stable/pipeline.html
#'
#' @param sif_cnvkit [REQUIRED] Path to cnvkit sif file.oms.
#' @param cnr [REQUIRED] Path to CNR as produced by fix_cnvkit command.
#' @param cns [OPTIONAL] Path to CNS as produced by segment_cnvkit command.
#' @param cn_thr [OPTIONAL] Copy number threshold. Default 0.5
#' @param min_probes [OPTIONAL] Minimum number of probes to show a gene. Default 3
#' @param male_reference [OPTIONAL]  Assume inputs were normalized to a male reference. Default FALSE
#' @param gender [OPTIONAL] Assume sample gender. Default male
#' @param shift [OPTIONAL] Adjust the X and Y chromosomes according to sample sex. Default TRUE
#' @param title [OPTIONAL] Sample ID. Default NULL.
#' @param output_name [OPTIONAL] Name for the output. If not given the name of the first tumour sample of the samples will be used.
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

  diagram_cnvkit=function(
    chrs=build_default_chr_list()$canonical,
    sif_cnvkit=build_default_sif_list()$sif_cnvkit,
    cnr=NULL,
    cns=NULL,
    cn_thr=0.5,
    min_probes=3,
    male_reference=FALSE,
    gender="male",
    shift=TRUE,
    title=NULL,
    ...
  ){

    
    run_main=function(
      .env
    ){


      if(!is.null(gender)){
        gender=paste0(" -x ",gender)
      }


      if(is.null(cn_thr)){

        cn_thr=paste0(" -t ",cn_thr)

      }

      if(is.null(min_probes)){

        min_probes=paste0(" -m ",min_probes)

      }

      if(!is.null(title)){
        title=paste0(" --title ",title)
      }


      add=""

      if(!shift){
        add=paste(add," --no-shift-xy ")
      }

      if(male_reference){
        add=paste(add," -y ")
      }
      
      tmp_cns=paste0(cns,".tmp")
      tmp_cnr=paste0(cnr,".tmp")

      .main$out_files$diagram=paste0(out_file_dir,"/",input_id,".diagram.pdf")
      cnr_tmp=paste0(cnr,".tmp")
      cns_tmp=""
      .main$exec_code=paste0("cat ",cnr," | head -n 1 > ",cnr_tmp," && cat ",cnr," | grep \"",
      paste0(paste0("^",chrs),collapse="\\|"),"\" >> ",cnr_tmp)
      if(is.null(cns)){
        cns_tmp=paste0(cns,".tmp")
        cns_code=paste(" -s ",cns_tmp)
        .main$exec_code=paste0(.main$exec_code," && cat ",
        cns," |  head -n 1 > ",cns_tmp," && ",paste0("cat ",cns," | grep \" ",
        paste0(paste0("^",unlist(chrs),collapse="\\|")," \" >> ",cns_tmp)))
      }
      
      .main$exec_code=paste(
        .main$exec_code,"; singularity exec -H /:/home ",sif_cnvkit,
        " cnvkit.py diagram -o ",.main$out_files$diagram,
        cns_code,gender,
        add,title,
        min_probes,
        cn_thr,
        cnr_tmp,
        "&& rm ",
        cnr_tmp,
        cns_tmp
      )

      run_job(.env=.this.env)
      .env$.main <- .main

    }


    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="cnr"
    )
  
    launch(.env=.base.env)



  }

