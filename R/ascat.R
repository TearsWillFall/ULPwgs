
call_ascat=function(
    rdata=NULL,
    selected=NULL,
    output_name="",
    output_dir=".",
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=4,ram=4,mode="local",
    executor_id=make_unique_id("hybridIchorCNA"),
    task_name="hybridIchorCNA",time="48:0:0",
    update_time=60,
    wait=FALSE,hold=NULL
){

  job_reports=prepare_ascat()
  ascat.bc = ascat.loadData(Tumor_LogR_file = "Tumor_LogR.txt", Tumor_BAF_file = "Tumor_BAF.txt", Germline_LogR_file = "Germline_LogR.txt", Germline_BAF_file = "Germline_BAF.txt", gender = 'XX', genomeVersion = "hg19")
  ascat.plotRawData(ascat.bc, img.prefix = "Before_correction_")
  ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile = "GC_file.txt", replictimingfile = "RT_file.txt")
  ascat.plotRawData(ascat.bc, img.prefix = "After_correction_")
  ascat.bc = ascat.aspcf(ascat.bc)
  ascat.plotSegmentedData(ascat.bc)
  ascat.output = ascat.runAscat(ascat.bc, gamma=1, write_segments = T)
  QC = ascat.metrics(ascat.bc,ascat.output)
  save(ascat.bc, ascat.output, QC, file = 'ASCAT_objects.Rdata')

}




#' Prepare ASCAT data from PCF select BAM files
#' 
#' @param bin_allele_counter [REQUIRED] Path to allele_counter binary.
#' @param tumour [REQUIRED] Path to tumour BAM file.
#' @param normal [REQUIRED] Path to normal BAM file.
#' @param version [OPTIONAL] PCF Select panel version. Default V3.
#' @param ref_dataset [OPTIONAL] Reference dataset for panel allele and loci selection. Default battenberg
#' @param gender [OPTIONAL] Sample gender. Default XY.
#' @param patient_id [OPTIONAL] Patient identifier. Default NULL
#' @param genome_version [OPTIONAL] Genome reference version. Default HG19.
#' @param ascat_ref [OPTIONAL] List with default references.
#' @param batch_config [OPTIONAL] List with default references.
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


prepare_ascat=function(
  rdata=NULL,
  selected=NULL,
  bin_allele_counter=build_default_tool_binary_list()$bin_allele_counter,
  tumour=NULL,
  normal=NULL,
  patient_id=NULL,
  ref_dataset="battenberg",
  gender="XY",
  genome_version="HG19",
  panel_version="V3",
  output_dir=".",
  verbose=FALSE,
  ascat_ref=build_default_reference_list(),
  batch_config=build_default_preprocess_config(),
  threads=8,
  ram=4,mode="local",
  executor_id=make_unique_id("prepareASCAT"),
  task_name="prepareASCAT",time="48:0:0",
  update_time=60,
  wait=FALSE,hold=NULL
){


if(!is.null(rdata)){
    load(rdata)
    if(!is.null(selected)){
      region=region_list[selected]
    }
  }

  if(is.null(tumour)|is.null(normal)){
    stop("A pair of tumour and normal BAM files is required")
  }

  
  if(is.null(patient_id)){
    stop("Patient ID argument is required")
  }


  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir,name=paste0(patient_id,"/",get_file_name(tumour),"/",genome_version,"/",panel_version,"/",ref_dataset))
  job=build_job(executor_id=executor_id,task=task_id)
  

  out_file_tumour_log2=paste0(get_file_name(tumour),".LogR.txt")
  out_file_tumour_baf=paste0(get_file_name(tumour),".BAF.txt")
  out_file_normal_log2=paste0(get_file_name(normal),".LogR.txt")
  out_file_normal_baf=paste0(get_file_name(normal),".BAF.txt")
  
  
  allele_prefix=ascat_ref[[toupper(genome_version)]]$panel[[paste0("PCF_",panel_version)]]$ascat[[ref_dataset]]$alleleData
  loci_prefix=ascat_ref[[toupper(genome_version)]]$panel[[paste0("PCF_",panel_version)]]$ascat[[ref_dataset]]$alleleData

  

  exec_code=paste0("Rscript -e \"","setwd(\\\"",
            out_file_dir,"\\\");library(ASCAT);ASCAT::ascat.prepareHTS(tumourseqfile =\\\"",tumour,"\\\"",
            ",normalseqfile =\\\"",normal,"\\\"",
            ",tumourname = \\\"",get_file_name(tumour),"\\\"",
            ",normalname = \\\"",get_file_name(normal),"\\\"",
            ",allelecounter_exe =\\\"",bin_allele_counter,"\\\"",
            ",alleles.prefix=\\\"",allele_prefix,"\\\"",
            ",loci.prefix=\\\"",loci_prefix,"\\\"",
            ",gender=\\\"",gender,"\\\"",
            ",genomeVersion=\\\"",tolower(genome_version),"\\\"",
            ",nthreads=\\\"", threads,"\\\"",
            ",tumourLogR_file=\\\"", out_file_tumour_log2,"\\\"",
            ",tumourBAF_file=\\\"", out_file_tumour_baf,"\\\"",
            ",normalLogR_file=\\\"", out_file_normal_log2,"\\\"",
            ",normalBAF_file=\\\"", out_file_normal_baf,"\\\")\""
  )


 if(mode=="batch"){
       out_file_dir2=set_dir(dir=out_file_dir,name="batch")
       batch_code=build_job_exec(job=job,time=time,ram=ram,threads=threads,
       output_dir=out_file_dir2,hold=hold)
       exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
  }


  if(verbose){
  print_verbose(job=job,arg=argg,exec_code=exec_code)
  }

  error=execute_job(exec_code=exec_code)
  
  if(error!=0){
      stop("gatk failed to run due to unknown error.
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
    tumour_log2=out_file_tumour_log2,
    tumour_baf=out_file_tumour_baf,
    normal_log2=out_file_normal_log2,
    normal_baf=out_file_normal_baf)
  )

  return(job_report)


}








