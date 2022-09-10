
#' Variant Calling using Mutect2
#'
#' This function functions calls Mutect2 for variant calling.
#' If a vector of tumour samples are provided these will be processed in multi-sample mode.
#' To run in tumour-normal mode suppply a single tumour and normal sample.
#' If no normal is supplied this will run in tumour only.
#' TO DO// Implement mitochondrial mode feature
#' 
#' For more information read:
#' https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2
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


    exec_code=paste0("singularity exec -H ",getwd(),":/home ",sif_cnvkit,
  " /gatk/gatk  ",f1r2_list," -O ",out_file)


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
#' @param sif_cnvkit [REQUIRED] Path to cnvkit sif file.
#' @param ref_genome [REQUIRED] Path to reference genome fasta file.
#' @param rdata [OPTIONAL] Import R data information with list of BAM.
#' @param selected [OPTIONAL] Select BAM from list.
#' @param pon [OPTIONAL] Path to Panel of Normals. Default none
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
    output_name="reference",
    normals="",
    output_dir=".",
    exclude_regions="",
    target="",
    antitarget="",
    gap_size=5000,
    verbose=FALSE,
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


    out_file=paste0(out_file_dir,"/",output_name,".cnn")


    if(target!=""){
      target=paste0(" -t ",target)
    }

    if(antitarget!=""){
      antitarget=paste0(" -a ", antitarget)
    }

    if(access!=""){
      paste0(" -g ",acceess)
    }

    exec_code=paste("singularity exec -H ",paste0(getwd(),":/home "),sif_cnvkit,
    " cnvkit.py batch -n ",paste0(normals,collapse=" "), " --output-reference ",
    out_file," -f ", ref_genome,access,target,antitarget)


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

    out_file=paste0(out_file_dir,"/",id,".targets.bed")

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
    bed="",
    access="",
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
    if(!is.null(min_avg_size)){
      add=paste0(" -m ",min_avg_size)
    }

    out_file=paste0(out_file_dir,"/",id,".antitargets.bed")


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