
#' Variant Calling using Mutect2
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


segmentation_cnvkit=function(
    rdata=NULL,selected=FALSE,
    sif_cnvkit=build_default_sif_list()$sif_cnvkit,
    tumor="",
    normal="",
    pon="",
    seg_method="cbs",
    seq_method="hybrid",
    baits=build_default_reference_list()$HG19$panel$PCF_V2$bed$bait,
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    access=build_default_reference_list()$HG19$reference$access_5k,
    pon_name="",
    output_dir=".",
    diagram=TRUE,
    scatter=TRUE,
    male=TRUE,
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=4,ram=4,mode="local",
    executor_id=make_unique_id("segmentationCNVkit"),
    task_name="segmentationCNVkit",time="48:0:0",
    update_time=60,
    wait=FALSE,hold=""
  ){


  if(!is.null(rdata)){
      load(rdata)
      if(!is.null(selected)){
          region=region_list[selected]
      }
  }

    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir,name="segmentation")
    job=build_job(executor_id=executor_id,task=task_id)
    

    gender=""
    if(male){
      gender=" -y "
    }



    seg_method=paste0("--segment-method ",seg_method)
    seq_method=paste0("--method ",seq_method)


    add=""
    if(scatter){
      add=paste(add," --scatter ")
    }

    if(diagram){
      add=paste(add," --diagram ")
    }


    if (pon_name!=""){
      pon_name=paste(" --output-reference ",pon_name)
    }

    if (access!=""){
      access=paste(" --access ",access)
    }



    if (fasta!=""){
      fasta=paste(" --fasta ",fasta)
    }

    if (baits!=""){
      baits=paste(  " --targets ",baits)
    }

    exec_code=paste("singularity exec -H ",paste0(getwd(),":/home "),sif_cnvkit,
    " cnvkit.py access -o ",out_file,exclude_regions," -s ",gap_size, ref_genome)



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
    output_name="access",
    output_dir=".",
    exclude_regions="",
    gap_size=5000,
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=1,mode="local",
    executor_id=make_unique_id("accessCNVkit"),
    task_name="accessCNVkit",time="48:0:0",
    update_time=60,
    wait=FALSE,hold=""
  ){


    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir)
    job=build_job(executor_id=executor_id,task_id=task_id)


    out_file=paste0(out_file_dir,"/",output_name,".bed")

    if(exclude_regions!=""){
      exclude_regions=paste0(" -x ",exclude_regions)
    }
      
    exec_code=paste("singularity exec -H ",paste0(getwd(),":/home "),sif_cnvkit,
    " cnvkit.py access -o ",out_file,exclude_regions," -s ",gap_size, ref_genome)


    if(mode=="batch"){
        out_file_dir2=set_dir(dir=out_file_dir,name="batch")
        batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
        threads=threads,output_dir=out_file_dir2)
        exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
    }

    if(verbose){
        print_verbose(job=job,arg=argg,exec_code=exec_code)
    }

    error=execute_job(exec_code=exec_code)
    
    if(error!=0){
      stop("cnvkit failed to run due to unknown error.
      Check std error for more information.")
    }

    job_report=build_job_report(
      job_id=job,
      executor_id=executor_id,
      exec_code=exec_code, 
      task_id=task_id,
      input_args = argg,
      out_file_dir=out_file_dir,
      out_files=list(
        access=out_file)
    )


    if(wait&&mode=="batch"){
      job_validator(job=job_report$job_id,time=update_time,
      verbose=verbose,threads=threads)
    }

    return(job_report)

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
    output_name="reference",
    normals="",
    output_dir=".",
    target="",
    bin_size_target=75,
    bin_size_antitarget=500000,
    min_bin_size_antitarget=NULL,
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=1,mode="local",
    executor_id=make_unique_id("createPonCNVkit"),
    task_name="createPonCNVkit",time="48:0:0",
    update_time=60,
    wait=FALSE,hold=""
  ){

    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir)
    job=build_job(executor_id=executor_id,task_id=task_id)


    out_file=paste0(out_file_dir,"/",output_name,".cnn")


    if(target!=""){
      target=paste0(" -t ",target)
    }


    if(access!=""){
      paste0(" -g ",access)
    }
    if(!is.null(min_bin_size_antitarget)){
      min_bin_size_antitarget=paste0(" --antitarget-min-size ",min_bin_size_antitarget)
    }

    exec_code=paste("singularity exec -H ",paste0(getwd(),":/home "),sif_cnvkit,
    " cnvkit.py batch -p ",threads, " -n ",paste0(normals,collapse=" "),
     " --output-reference ",out_file," -f ", ref_genome,access,target,
     " --target-avg-size ",bin_size_target,
     " --antitarget-avg-size ",bin_size_antitarget,min_bin_size_antitarget
     )


    if(mode=="batch"){
        out_file_dir2=set_dir(dir=out_file_dir,name="batch")
        batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
        threads=threads,output_dir=out_file_dir2)
        exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
    }

    if(verbose){
        print_verbose(job=job,arg=argg,exec_code=exec_code)
    }

    error=execute_job(exec_code=exec_code)
    
    if(error!=0){
      stop("cnvkit failed to run due to unknown error.
      Check std error for more information.")
    }

    job_report=build_job_report(
      job_id=job,
      executor_id=executor_id,
      exec_code=exec_code, 
      task_id=task_id,
      input_args = argg,
      out_file_dir=out_file_dir,
      out_files=list(
        access=out_file)
    )


    if(wait&&mode=="batch"){
      job_validator(job=job_report$job_id,time=update_time,
      verbose=verbose,threads=threads)
    }

    return(job_report)

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
    bed="",
    annotation="",
    output_name="",
    output_dir=".",
    split=FALSE,
    short_names=FALSE,
    bin_size=100,
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=1,mode="local",
    executor_id=make_unique_id("targetCNVkit"),
    task_name="targetCNVkit",time="48:0:0",
    update_time=60,
    wait=FALSE,hold=""
  ){

    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir)
    job=build_job(executor_id=executor_id,task_id=task_id)

      
    id=""
    if(output_name!=""){
      id=output_name
    }else{
      id=get_file_name(bed)
    }

    if(annotation!=""){
      annotation=paste0(" --annotate ",annotation)
    }

    add=""
    if(short_names){
      add=paste(add," --short-names ")
    }

    out_file=paste0(out_file_dir,"/",id,".binned.targets.bed")

    exec_code=paste("singularity exec -H ",paste0(getwd(),":/home "),sif_cnvkit,
    " cnvkit.py target --split -a ",bin_size," -o ",out_file, add, bed)


    if(mode=="batch"){
        out_file_dir2=set_dir(dir=out_file_dir,name="batch")
        batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
        threads=threads,output_dir=out_file_dir2)
        exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
    }

    if(verbose){
        print_verbose(job=job,arg=argg,exec_code=exec_code)
    }

    error=execute_job(exec_code=exec_code)
    
    if(error!=0){
      stop("cnvkit failed to run due to unknown error.
      Check std error for more information.")
    }

    job_report=build_job_report(
      job_id=job,
      executor_id=executor_id,
      exec_code=exec_code, 
      task_id=task_id,
      input_args = argg,
      out_file_dir=out_file_dir,
      out_files=list(
        target=out_file)
    )


    if(wait&&mode=="batch"){
      job_validator(job=job_report$job_id,time=update_time,
      verbose=verbose,threads=threads)
    }

    return(job_report)

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
    bed="",
    output_name="",
    output_dir=".",
    bin_size=500000,
    min_bin_size=NULL,
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=1,mode="local",
    executor_id=make_unique_id("antitargetCNVkit"),
    task_name="antitargetCNVkit",time="48:0:0",
    update_time=60,
    wait=FALSE,hold=""
  ){

    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir)
    job=build_job(executor_id=executor_id,task_id=task_id)


    id=""
    if(output_name!=""){
      id=output_name
    }else{
      id=get_file_name(bed)
    }

    add=""
    if(!is.null(min_bin_size)){
      add=paste0(" -m ",min_bin_size)
    }

    out_file=paste0(out_file_dir,"/",id,".binned.antitargets.bed")


    if(access!=""){
      access=paste0(" -g ",access)
    }


    exec_code=paste("singularity exec -H ",paste0(getwd(),":/home "),sif_cnvkit,
    " cnvkit.py antitarget -a ",bin_size,access," -o ",out_file, add, bed)

    if(mode=="batch"){
        out_file_dir2=set_dir(dir=out_file_dir,name="batch")
        batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
        threads=threads,output_dir=out_file_dir2)
        exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
    }

    if(verbose){
        print_verbose(job=job,arg=argg,exec_code=exec_code)
    }

    error=execute_job(exec_code=exec_code)
    
    if(error!=0){
      stop("cnvkit failed to run due to unknown error.
      Check std error for more information.")
    }

    job_report=build_job_report(
      job_id=job,
      executor_id=executor_id,
      exec_code=exec_code, 
      task_id=task_id,
      input_args = argg,
      out_file_dir=out_file_dir,
      out_files=list(
        antitarget=out_file)
    )


    if(wait&&mode=="batch"){
      job_validator(job=job_report$job_id,time=update_time,
      verbose=verbose,threads=threads)
    }

    return(job_report)

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
    bed="",
    output_name="",
    annotation="",
    output_dir=".",
    bin_size_antitarget=100000,
    bin_size_target=100,
    min_bin_size_antitarget=NULL,
    split=FALSE,
    short_names=FALSE,
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=1,mode="local",
    executor_id=make_unique_id("bintargetCNVkit"),
    task_name="bintargetCNVkit",time="48:0:0",
    update_time=60,
    wait=FALSE,hold=""
){

    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir)
    job=build_job(executor_id=executor_id,task_id=task_id)

    
    jobs_report=build_job_report(
      job_id=job,
      executor_id=executor_id,
      exec_code=list(), 
      task_id=task_id,
      input_args = argg,
      out_file_dir=out_file_dir,
      out_files=list()
    )


    jobs_report[["steps"]][["TargetsCNVkit"]]<-create_target_cnvkit(
      sif_cnvkit=sif_cnvkit,
      bed=bed,
      annotation=annotation,
      output_name=output_name,
      output_dir=out_file_dir,
      split=split,
      short_names=short_names,
      bin_size=bin_size_target,
      verbose=verbose,
      batch_config=batch_config,
      threads=1,ram=1,mode=mode,
      executor_id=task_id,
      time=time,
      hold=hold
  )


   jobs_report[["steps"]][["AntitargetsCNVkit"]]<-create_antitarget_cnvkit(
      sif_cnvkit=sif_cnvkit,
      access=access,
      bed=jobs_report[["steps"]][["TargetsCNVkit"]]$out_files$target,
      output_name=output_name,
      output_dir=out_file_dir,
      bin_size=bin_size_antitarget,
      min_bin_size=min_bin_size_antitarget,
      verbose=verbose,
      batch_config=batch_config,
      threads=1,ram=1,mode=mode,
      executor_id=task_id,
      time=time,
      hold=jobs_report[["steps"]][["TargetsCNVkit"]]$job_id
  )

  return(jobs_report)

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
    bed="",
    bams="",
    annotation="",
    output_name="",
    output_dir=".",
    seq_method="hybrid",
    split=FALSE,
    short_names=FALSE,
    bp_per_bin=100000,
    min_bin_size_target=20,
    max_bin_size_target=20000,
    min_bin_size_antitarget=500,
    max_bin_size_antitarget=500000,
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=1,mode="local",
    executor_id=make_unique_id("autobinCNVkit"),
    task_name="autobinCNVkit",time="48:0:0",
    update_time=60,
    wait=FALSE,hold=""
  ){

    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir)
    job=build_job(executor_id=executor_id,task_id=task_id)


    if(annotation!=""){
      annotation=paste0(" --annotate ",annotation)
    }

    if(access!=""){
      access=paste0(" -g ",access)
    }

    add=""
    if(short_names){
      add=paste(add," --short-names ")
    }

    exec_code=paste("singularity exec -H ",paste0(getwd(),":/home "),sif_cnvkit,
    " cnvkit.py autobin -t ",bed," -f ",ref_genome,
    " -b ",bp_per_bin, " -m ",seq_method,access,
    "  --target-max-size ",max_bin_size_target,
    "  --target-min-size ",min_bin_size_target,
    "  --antitarget-min-size ",min_bin_size_antitarget,
    "  --antitarget-max-size ",max_bin_size_antitarget,add,annotation,
     paste0(bams,collapse=" "))
de
    if(mode=="batch"){
        out_file_dir2=set_dir(dir=out_file_dir,name="batch")
        batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
        threads=threads,output_dir=out_file_dir2)
        exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
    }

    if(verbose){
        print_verbose(job=job,arg=argg,exec_code=exec_code)
    }

    error=execute_job(exec_code=exec_code)
    
    if(error!=0){
      stop("cnvkit failed to run due to unknown error.
      Check std error for more information.")
    }

    job_report=build_job_report(
      job_id=job,
      executor_id=executor_id,
      exec_code=exec_code, 
      task_id=task_id,
      input_args = argg,
      out_file_dir=out_file_dir,
      out_files=list(
        target=out_file)
    )


    if(wait&&mode=="batch"){
      job_validator(job=job_report$job_id,time=update_time,
      verbose=verbose,threads=threads)
    }

    return(job_report)

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
    rdata=NULL,
    selected=NULL,
    sif_cnvkit=build_default_sif_list()$sif_cnvkit,
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    bed="",
    bam="",
    output_name="",
    output_dir=".",
    read_count=FALSE,
    min_mapq=0,
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=1,mode="local",
    executor_id=make_unique_id("coverageCNVkit"),
    task_name="coverageCNVkit",time="48:0:0",
    update_time=60,
    wait=FALSE,hold=""
  ){

    if(!is.null(rdata)){
      load(rdata)
      if(!is.null(selected)){
        bam=bam_list[selected]
      }
    }

    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir)
    job=build_job(executor_id=executor_id,task_id=task_id)

         
    id=""
    if(output_name!=""){
      id=output_name
    }else{
      id=get_file_name(bam)
    }
  
    add=""
    if(read_count){
      add=paste(add," -c ")
    }

    out_file=paste0(out_file_dir,"/",id,".cnn")

    exec_code=paste("singularity exec -H ",paste0(getwd(),":/home "),sif_cnvkit,
    " cnvkit.py coverage -p ",threads, "-q ",min_mapq,
    " -f ",ref_genome," -o ",out_file, add, bam, bed)

    if(mode=="batch"){
        out_file_dir2=set_dir(dir=out_file_dir,name="batch")
        batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
        threads=threads,output_dir=out_file_dir2)
        exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
    }

    if(verbose){
        print_verbose(job=job,arg=argg,exec_code=exec_code)
    }

    error=execute_job(exec_code=exec_code)
    
    if(error!=0){
      stop("cnvkit failed to run due to unknown error.
      Check std error for more information.")
    }

    job_report=build_job_report(
      job_id=job,
      executor_id=executor_id,
      exec_code=exec_code, 
      task_id=task_id,
      input_args = argg,
      out_file_dir=out_file_dir,
      out_files=list(
        cnn=out_file)
    )


    if(wait&&mode=="batch"){
      job_validator(job=job_report$job_id,time=update_time,
      verbose=verbose,threads=threads)
    }

    return(job_report)

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


  parallel_sample_coverage_cnvkit=function(
    sif_cnvkit=build_default_sif_list()$sif_cnvkit,
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    bed="",
    bams="",
    output_dir=".",
    read_count=FALSE,
    min_mapq=0,
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=1,mode="local",
    executor_id=make_unique_id("coverageCNVkit"),
    task_name="coverageCNVkit",time="48:0:0",
    update_time=60,
    wait=FALSE,hold=""
  ){

    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir,name="cnn")
    tmp_dir=set_dir(dir=output_dir,name="tmp")
    job=build_job(executor_id=executor_id,task_id=task_id)

      jobs_report=build_job_report(
          job_id=job,
          executor_id=executor_id,
          exec_code=list(), 
          task_id=task_id,
          input_args=argg,
          out_file_dir=out_file_dir,
          out_files=list(
            )
      ) 

    bam_list=bams
    names(bam_list)=Vectorize(get_file_name)(bams)
    

    if(mode=="local"){
    jobs_report[["steps"]][["par_sample_coverage_cnvkit"]]<-
    parallel::mclapply(bam_list,FUN=function(bam){
      job_report <- coverage_cnvkit(
                sif_cnvkit=sif_cnvkit,
                ref_genome=ref_genome,
                bed=bed,
                bam=bam,
                output_name=get_file_name(bam),
                output_dir=out_file_dir,
                read_count=read_count,
                min_mapq=min_mapq,
                verbose=verbose,
                batch_config=batch_config,
                executor_id=task_id
              )
    },mc.cores=threads)
    
    }else if(mode=="batch"){

          rdata_file=paste0(tmp_dir,"/",job,".samples.RData")
          output_dir=out_file_dir
          save(bam_list,bed,sif_cnvkit,ref_genome,
          output_dir,read_count,min_mapq,output_dir,verbose,file = rdata_file)
          exec_code=paste0("Rscript -e \"ULPwgs::coverage_cnvkit(rdata=\\\"",
          rdata_file,"\\\",selected=$SGE_TASK_ID)\"")
          out_file_dir2=set_dir(dir=out_file_dir,name="batch")
          batch_code=build_job_exec(job=job,time=time,ram=ram,
          threads=1,output_dir=out_file_dir2,
          hold=hold,array=length(bam_list))
          exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)

          if(verbose){
              print_verbose(job=job,arg=argg,exec_code=exec_code)
          }
          error=execute_job(exec_code=exec_code)
          if(error!=0){
              stop("cnvkit failed to run due to unknown error.
              Check std error for more information.")
          }
         
         jobs_report[["steps"]][["par_sample_coverage_cnvkit"]]<- build_job_report(
              job_id=job,
              executor_id=executor_id,
              exec_code=exec_code, 
              task_id=task_id,
              input_args=argg,
              out_file_dir=out_file_dir,
              out_files=list(
                  cnn=paste0(out_file_dir,"/",names(bam_list),".cnn")
              )
        )
    }

      return(jobs_report)


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
    cnn="",
    output_name="reference",
    output_dir=".",
    gender="male",
    target="",
    antitarget="",
    gc=TRUE,
    edge=TRUE,
    rmask=TRUE,
    verbose=FALSE,
    cluster=FALSE,
    min_cluster_size=10,
    male_reference=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=1,mode="local",
    executor_id=make_unique_id("referenceCNVkit"),
    task_name="referenceCNVkit",time="48:0:0",
    update_time=60,
    wait=FALSE,hold=""
  ){

    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir)
    job=build_job(executor_id=executor_id,task_id=task_id)

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
  


    out_file=paste0(out_file_dir,"/",output_name,".cnn")

    exec_code=paste("singularity exec -H ",paste0(getwd(),":/home "),sif_cnvkit,
    " cnvkit.py reference ", gender, min_cluster_size,
    " -f ",ref_genome," -o ",out_file, add, paste0(cnn,collapse=" "))

    if(mode=="batch"){
        out_file_dir2=set_dir(dir=out_file_dir,name="batch")
        batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
        threads=threads,output_dir=out_file_dir2)
        exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
    }

    if(verbose){
        print_verbose(job=job,arg=argg,exec_code=exec_code)
    }

    error=execute_job(exec_code=exec_code)
    
    if(error!=0){
      stop("cnvkit failed to run due to unknown error.
      Check std error for more information.")
    }

    job_report=build_job_report(
      job_id=job,
      executor_id=executor_id,
      exec_code=exec_code, 
      task_id=task_id,
      input_args = argg,
      out_file_dir=out_file_dir,
      out_files=list(
        cnn=out_file)
    )


    if(wait&&mode=="batch"){
      job_validator(job=job_report$job_id,time=update_time,
      verbose=verbose,threads=threads)
    }

    return(job_report)

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
    rdata=NULL,
    selected=NULL,
    sif_cnvkit=build_default_sif_list()$sif_cnvkit,
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    pon=build_default_reference_list()$HG19$panel$PCF_V3$variant$pon_cn,
    target="",
    antitarget="",
    output_name="",
    output_dir=".",
    gc=TRUE,
    edge=TRUE,
    rmask=TRUE,
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=1,mode="local",
    executor_id=make_unique_id("fixCNVkit"),
    task_name="fixCNVkit",time="48:0:0",
    update_time=60,
    wait=FALSE,hold=""
  ){

    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir)
    job=build_job(executor_id=executor_id,task_id=task_id)

    if(!is.null(rdata)){
      load(rdata)
      if(!is.null(selected)){
        target=sample_list[selected,]$targets
        antitarget=sample_list[selected,]$antitargets
      }
    }

    id=""
    if(output_name!=""){
      id=output_name
    }else{
      id=get_file_name(target)
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


    out_file=paste0(out_file_dir,"/",id,".cnr")

    exec_code=paste("singularity exec -H ",paste0(getwd(),":/home "),sif_cnvkit,
    " cnvkit.py fix -o ",out_file, add, target, antitarget, pon)

    if(mode=="batch"){
        out_file_dir2=set_dir(dir=out_file_dir,name="batch")
        batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
        threads=threads,output_dir=out_file_dir2)
        exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
    }

    if(verbose){
        print_verbose(job=job,arg=argg,exec_code=exec_code)
    }

    error=execute_job(exec_code=exec_code)
    
    if(error!=0){
      stop("cnvkit failed to run due to unknown error.
      Check std error for more information.")
    }

    job_report=build_job_report(
      job_id=job,
      executor_id=executor_id,
      exec_code=exec_code, 
      task_id=task_id,
      input_args = argg,
      out_file_dir=out_file_dir,
      out_files=list(
        cnn=out_file)
    )


    if(wait&&mode=="batch"){
      job_validator(job=job_report$job_id,time=update_time,
      verbose=verbose,threads=threads)
    }

    return(job_report)

  }




#' Wrapper around parallel intergration for fix function from CNVkit
#'
#' This function wraps around fix function for CNVkit
#' This function correct sample coverage using a reference panel of normals files.
#' Additional parameters can be used to exclude regions and modify the average bin size.
#' 
#' 
#' For more information read:
#' https://cnvkit.readthedocs.io/en/stable/pipeline.html
#'
#' @param sif_cnvkit [REQUIRED] Path to cnvkit sif file.
#' @param ref_genoma [REQUIRED] Path to reference genome.
#' @param pon [REQUIRED] Path to panel of normals location.
#' @param targets [REQUIRED] Path to target CNN files.
#' @param antitargets [REQUIRED] Path to antitarget CNN files.
#' @param gc [OPTIONAL] Disble GC bias correction. Default FALSE.
#' @param edge [OPTIONAL] Disble edge correction. Default FALSE.
#' @param rmask [OPTIONAL] Disble repeat mask correction. Default FALSE.
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

  parallel_sample_fix_cnvkit=function(
    sif_cnvkit=build_default_sif_list()$sif_cnvkit,
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    pon=build_default_reference_list()$HG19$panel$PCF_V3$variant$pon_cn,
    targets="",
    antitargets="",
    output_name="",
    output_dir=".",
    gc=TRUE,
    edge=TRUE,
    rmask=TRUE,
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=1,mode="local",
    executor_id=make_unique_id("parFixCNVkit"),
    task_name="parFixCNVkit",time="48:0:0",
    update_time=60,
    wait=FALSE,hold=""
  ){

    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir,name="cnr")
    tmp_dir=set_dir(dir=output_dir,name="tmp")
    job=build_job(executor_id=executor_id,task_id=task_id)

      jobs_report=build_job_report(
          job_id=job,
          executor_id=executor_id,
          exec_code=list(), 
          task_id=task_id,
          input_args=argg,
          out_file_dir=out_file_dir,
          out_files=list(
            )
      ) 

    sample_list=data.frame(targets=targets,
    antitargets=antitargets,
    names=Vectorize(get_file_name)(targets))
    row.names(sample_list)= sample_list$names
    

    if(mode=="local"){
    jobs_report[["steps"]][["par_sample_fix_cnvkit"]]<-
    parallel::mclapply(sample_list,FUN=function(smpl){
      job_report <- fix_cnvkit(
                sif_cnvkit=sif_cnvkit,
                ref_genome=ref_genome,
                pon=pon,
                target=smpl$targets,
                antitarget=smpl$antitargets,
                output_name=smpl$names,
                output_dir=out_file_dir,
                gc=gc,
                edge=edge,
                rmask=rmask,
                verbose=verbose,
                batch_config=batch_config,
                executor_id=task_id
              )
    },mc.cores=threads)
    
    }else if(mode=="batch"){

          rdata_file=paste0(tmp_dir,"/",job,".samples.RData")
          output_dir=out_file_dir
          save(sample_list,bed,sif_cnvkit,ref_genome,
          gc,edge,rmask,output_dir,verbose,file = rdata_file)
          exec_code=paste0("Rscript -e \"ULPwgs::fix_cnvkit(rdata=\\\"",
          rdata_file,"\\\",selected=$SGE_TASK_ID)\"")
          out_file_dir2=set_dir(dir=out_file_dir,name="batch")
          batch_code=build_job_exec(job=job,time=time,ram=ram,
          threads=1,output_dir=out_file_dir2,
          hold=hold,array=nrow(sample_list))
          exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)

          if(verbose){
              print_verbose(job=job,arg=argg,exec_code=exec_code)
          }
          error=execute_job(exec_code=exec_code)
          if(error!=0){
              stop("gatk failed to run due to unknown error.
              Check std error for more information.")
          }
         
         jobs_report[["steps"]][["par_sample_fix_cnvkit"]]<- build_job_report(
              job_id=job,
              executor_id=executor_id,
              exec_code=exec_code, 
              task_id=task_id,
              input_args=argg,
              out_file_dir=out_file_dir,
              out_files=list(
                  cnn=paste0(out_file_dir,"/",sample_list$names,".cnr")
              )
        )
    }

      return(jobs_report)


  }




