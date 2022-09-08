








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
#' @param mnps [OPTIONAL] Report MNPs in vcf file.
#' @param contamination [OPTIONAL] Produce sample cross-contamination reports. Default TRUE.
#' @param orientation [OPTIONAL] Produce read orientation inforamtion. Default FALSE
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
    baits=build_default_reference_list()$HG19$panel$PCF_V2$bed$bait,
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    access=build_default_reference_list()$HG19$reference$access_5k,
    pon_name="",
    output_dir=".",
    diagram=TRUE,
    scatter=TRUE,
    male=TRUE,
    threads=3,
    verbose=FALSE
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


  if (pool_ref==""){
    if(verbose){
      print(paste(bin_path,"batch --output-dir ",out_file_dir,tumor,
      normal,baits,fasta,gender,pon_name,seg_method,output_dir," --p ",threads,add))
    }
    system(paste(bin_path,"batch --output-dir ",out_file_dir,tumor,
    normal,baits,fasta,gender,pon_name,seg_method,output_dir," --p ",threads,add))
  }else{
    pool_ref=paste(" -r ",pon)
    if(verbose){
      print(paste(bin_path,"batch ",tumor,pon,gender,ref_output,output_dir," --p ",threads,add))
    }
    system(paste(bin_path,"batch ",tumor,pon,gender,ref_output,output_dir," --p ",threads,add))
  }
}