#' Analyze a pair of tumour-normal samples using ASCAT
#' 
#' @param bin_allele_counter [REQUIRED] Path to allele_counter binary.
#' @param tumour [REQUIRED] Path to tumour BAM file.
#' @param normal [REQUIRED] Path to normal BAM file.
#' @param panel_version [OPTIONAL] PCF Select panel version. Default V3.
#' @param ref_dataset [OPTIONAL] Reference dataset for panel allele and loci selection. Default battenberg
#' @param gender [OPTIONAL] Sample gender. Default XY.
#' @param genome_version [OPTIONAL] Genome reference version. Default HG19.
#' @param ascat_ref [OPTIONAL] List with default references.
#' @param batch_config [OPTIONAL] List with default references.
#' @param gamma [OPTIONAL] Gamma parameter. Default 1.
#' @param penalty [OPTIONAL] Penalty parameter. Default 25.
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


call_ascat=function(
  rdata=NULL,
  selected=NULL,
  bin_allele_counter=build_default_tool_binary_list()$bin_allele_counter,
  ascat_ref=build_default_reference_list(),
  tumour=NULL,
  normal=NULL,
  ref_dataset="battenberg",
  gender="XY",
  genome_version="HG19",
  panel_version="V3",
  gamma=1,
  penalty=25,
  output_dir=".",
  verbose=FALSE,
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
      tumour=tumours_list[selected]

    }
  }


  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir,name=paste0(get_file_name(tumour),"/ascat_reports"))
  job=build_job(executor_id=executor_id,task=task_id)
  
  if(is.null(tumour)|is.null(normal)){
    stop("A pair of tumour and normal BAM files is required")
  }

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




  jobs_report[["steps"]][["prepare_ascat"]]=prepare_ascat(
    bin_allele_counter=bin_allele_counter,
    tumour=tumour,
    normal=normal,
    ref_dataset=ref_dataset,
    gender=gender,
    genome_version=genome_version,
    panel_version=panel_version,
    output_dir=paste0(out_file_dir,"/prepare_ascat"),
    verbose=verbose,
    ascat_ref=ascat_ref,
    batch_config=batch_config,
    threads=threads,
    ram=ram,mode=mode,
    executor_id=task_id,
    time=time,hold=hold
  )

  jobs_report[["steps"]][["process_ascat"]]=process_ascat(
    tumour_log2=normalizePath(jobs_report[["steps"]][["prepare_ascat"]]$out_files$tumour_log2),
    normal_log2=normalizePath(jobs_report[["steps"]][["prepare_ascat"]]$out_files$normal_log2),
    tumour_baf=normalizePath(jobs_report[["steps"]][["prepare_ascat"]]$out_files$tumour_baf),
    normal_baf=normalizePath(jobs_report[["steps"]][["prepare_ascat"]]$out_files$normal_baf),
    ref_dataset=ref_dataset,
    gender=gender,
    genome_version=genome_version,
    panel_version=panel_version,
    gamma=gamma,
    penalty=penalty,
    output_dir=paste0(out_file_dir,"/process_ascat"),
    verbose=verbose,
    ascat_ref=ascat_ref,
    batch_config=batch_config,
    threads=threads,
    ram=ram,mode=mode,
    executor_id=task_id,
    time=time,
    hold=hold
  )

  return(jobs_report)

}


#' Analyze in parallel multiple pair of tumour-normal samples using ASCAT
#' 
#' \\\\ TO DO FIX PARALLELIZATION IN LOCAL. Seems to be issues with running multiple files
#' 
#' @param bin_allele_counter [REQUIRED] Path to allele_counter binary.
#' @param tumour [REQUIRED] Path to tumour BAM file.
#' @param normal [REQUIRED] Path to normal BAM file.
#' @param patient_id [REQUIRED] Path to normal BAM file.
#' @param panel_version [OPTIONAL] PCF Select panel version. Default V3.
#' @param ref_dataset [OPTIONAL] Reference dataset for panel allele and loci selection. Default battenberg
#' @param gender [OPTIONAL] Sample gender. Default XY.
#' @param genome_version [OPTIONAL] Genome reference version. Default HG19.
#' @param ascat_ref [OPTIONAL] List with default references.
#' @param batch_config [OPTIONAL] List with default references.
#' @param gamma [OPTIONAL] Gamma parameter. Default 1.
#' @param penalty [OPTIONAL] Penalty parameter. Default 25.
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



parallel_samples_call_ascat=function(
  bin_allele_counter=build_default_tool_binary_list()$bin_allele_counter,
  tumour=NULL,
  normal=NULL,
  patient_id=NULL,
  ref_dataset="battenberg",
  gender="XY",
  genome_version="HG19",
  panel_version="V3",
  output_dir=".",
  gamma=1,
  penalty=25,
  verbose=FALSE,
  ascat_ref=build_default_reference_list(),
  batch_config=build_default_preprocess_config(),
  threads=8,
  ram=4,mode="local",
  executor_id=make_unique_id("parPrepareASCAT"),
  task_name="parPrepareASCAT",time="48:0:0",
  update_time=60,
  wait=FALSE,hold=NULL
){
  
  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir,name=patient_id)
  tmp_dir=set_dir(dir=out_file_dir,name="tmp")
 

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


    tumour_list=tumour
    names(tumour_list)=Vectorize(get_file_name)(tumour)

    if(mode=="local"){
      jobs_report[["steps"]][["par_sample_call_ascat"]]<-
      parallel::mclapply(tumour_list,FUN=function(tumour){
      job_report <- call_ascat(
              bin_allele_counter=bin_allele_counter,
              tumour=tumour,
              normal=normal,
              ref_dataset=ref_dataset,
              gender=gender,
              genome_version=genome_version,
              panel_version=panel_version,
              output_dir=out_file_dir,
              gamma=gamma,
              penalty=penalty,
              verbose=verbose,
              ascat_ref=ascat_ref,
              batch_config=ascat_ref,
              threads=4,
              ram=ram,mode=mode,
              executor_id=task_id,
              time=time,
              hold=hold
            )
      },mc.cores=4)
    }else if(mode=="batch"){

            rdata_file=paste0(tmp_dir,"/",job,".samples.RData")
            output_dir=out_file_dir
            executor_id=task_id
            save(
              tumour_list,
              bin_allele_counter,
              ref_dataset,
              gender,
              genome_version,
              output_dir,
              gamma,
              penalty,
              verbose,
              ascat_ref,
              file = rdata_file
            )
            exec_code=paste0("Rscript -e \"ULPwgs::call_cnvkit(rdata=\\\"",
            rdata_file,"\\\",selected=$SGE_TASK_ID)\"")
            out_file_dir2=set_dir(dir=out_file_dir,name="batch")
            batch_code=build_job_exec(job=job,time=time,ram=ram,
            threads=1,output_dir=out_file_dir2,
            hold=hold,array=length(tumour_list))
            exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)

            if(verbose){
                print_verbose(job=job,arg=argg,exec_code=exec_code)
            }
            error=execute_job(exec_code=exec_code)
            if(error!=0){
                stop("ascat failed to run due to unknown error.
                Check std error for more information.")
            }
            

        
            out_file_tumour_log2=paste0(names(tumour_list),".LogR.txt")
            out_file_tumour_baf=paste0(names(tumour_list),".BAF.txt")
            out_file_normal_log2=paste0(get_file_name(normal),".LogR.txt")
            out_file_normal_baf=paste0(get_file_name(normal),".BAF.txt")

            
            
            jobs_report[["steps"]][["par_sample_call_ascat"]]<- build_job_report(
                  job_id=job,
                  executor_id=executor_id,
                  exec_code=exec_code, 
                  task_id=task_id,
                  input_args=argg,
                  out_file_dir=out_file_dir,
                  out_files=list(
                    prepare_ascat=list(
                      data=list(
                        tumour_log2=paste0(out_file_dir,"/",
                          paste0(names(tumour_list),"/ascat_reports"),
                          "/",out_file_tumour_log2),
                        tumour_baf=paste0(out_file_dir,"/",
                          paste0(names(tumour_list),"/ascat_reports"),
                          "/",out_file_tumour_baf),
                        normal_log2=paste0(out_file_dir,"/",
                          paste0(names(tumour_list),"/ascat_reports"),
                          "/",out_file_normal_log2),
                        normal_baf=paste0(out_file_dir,"/",
                          paste0(names(tumour_list),"/ascat_reports"),
                          "/",out_file_normal_baf)
                      )
                  ),
                  process_ascat=list(
                      plots=list(
                        before_corr_normal=paste0(out_file_dir,"/",
                          paste0(names(tumour_list),"/ascat_reports"),
                          "/Before_correction_",names(tumour_list),".germline.png"),
                        after_corr_normal=paste0(out_file_dir,"/",
                          paste0(names(tumour_list),"/ascat_reports"),
                          "/After_correction_",names(tumour_list),".germline.png"),
                        before_corr_tumour=paste0(out_file_dir,"/",
                          paste0(names(tumour_list),"/ascat_reports"),
                          "/Before_correction_",names(tumour_list),".tumour.png"),
                        after_corr_tumour=paste0(out_file_dir,"/",
                          paste0(names(tumour_list),"/ascat_reports"),
                          "/After_correction_",names(tumour_list),".tumour.png"),
                        aspcf=paste0(out_file_dir,"/",
                          paste0(names(tumour_list),"/ascat_reports"),
                          "/",names(tumour_list),".ASPCF.png"),
                        ascat_profile=paste0(out_file_dir,"/",
                          paste0(names(tumour_list),"/ascat_reports"),
                          "/",names(tumour_list),".ASCATfprofile.png"),
                        raw_profile=paste0(out_file_dir,"/",
                          paste0(names(tumour_list),"/ascat_reports"),
                          "/",names(tumour_list),".rawprofile.png"),
                        sunrise=paste0(out_file_dir,"/",
                          paste0(names(tumour_list),"/ascat_reports"),
                          "/",names(tumour_list),".sunrise.png")
                      ),
                      data=list(
                        rdata=paste0(out_file_dir,"/",
                          paste0(names(tumour_list),"/ascat_reports"),
                          "/",rdata),
                        tumour_log2_pcfed=paste0(out_file_dir,"/",
                          paste0(names(tumour_list),"/ascat_reports"),
                          "/",names(tumour_list),".LogR.PCFed.txt"),
                        tumour_baf_pcfed=paste0(out_file_dir,"/",
                          paste0(names(tumour_list),"/ascat_reports"),
                          "/",names(tumour_list),".BAF.PCFed.txt"),
                        normal_log2_pcfed=paste0(out_file_dir,"/",
                          paste0(names(tumour_list),"/ascat_reports"),
                          "/",normal,".LogR.PCFed.txt"),
                        normal_baf_pcfed=paste0(out_file_dir,"/",
                          paste0(names(tumour_list),"/ascat_reports"),
                          "/",normal,".BAF.PCFed.txt"),
                        segments=paste0(out_file_dir,"/",
                          paste0(names(tumour_list),"/ascat_reports"),
                          "/",names(tumour_list),".segments.txt"),
                        segments_raw=paste0(out_file_dir,"/",
                          paste0(names(tumour_list),"/ascat_reports"),
                          "/",names(tumour_list),".segments_raw.txt")
                      )
                  )
                  
              )
            )
    }

  if(wait&&mode=="batch"){
    job_validator(job=unlist_lvl(jobs_report[["steps"]],var="job_id"),time=update_time,
    verbose=verbose,threads=threads)
  }

  return(jobs_report)



}






#' Prepare ASCAT data from PCF select BAM files
#' 
#' @param bin_allele_counter [REQUIRED] Path to allele_counter binary.
#' @param tumour [REQUIRED] Path to tumour BAM file.
#' @param normal [REQUIRED] Path to normal BAM file.
#' @param panel_version [OPTIONAL] PCF Select panel version. Default V3.
#' @param ref_dataset [OPTIONAL] Reference dataset for panel allele and loci selection. Default battenberg
#' @param gender [OPTIONAL] Sample gender. Default XY.
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
  bin_allele_counter=build_default_tool_binary_list()$bin_allele_counter,
  tumour=NULL,
  normal=NULL,
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

  

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir)
  job=build_job(executor_id=executor_id,task=task_id)
  
  if(is.null(tumour)|is.null(normal)){
    stop("A pair of tumour and normal BAM files is required")
  }
  

  out_file_tumour_log2=paste0(get_file_name(tumour),".LogR.txt")
  out_file_tumour_baf=paste0(get_file_name(tumour),".BAF.txt")
  out_file_normal_log2=paste0(get_file_name(normal),".LogR.txt")
  out_file_normal_baf=paste0(get_file_name(normal),".BAF.txt")
  
  
  allele_prefix=ascat_ref[[toupper(genome_version)]]$panel[[paste0("PCF_",panel_version)]]$ascat[[ref_dataset]]$allele
  loci_prefix=ascat_ref[[toupper(genome_version)]]$panel[[paste0("PCF_",panel_version)]]$ascat[[ref_dataset]]$loci

  

  exec_code=paste0("Rscript -e \"","options(warn = -1);setwd(\\\"",
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
      tumour_log2=paste0(out_file_dir,"/",out_file_tumour_log2),
      tumour_baf=paste0(out_file_dir,"/",out_file_tumour_baf),
      normal_log2=paste0(out_file_dir,"/",out_file_normal_log2),
      normal_baf=paste0(out_file_dir,"/",out_file_normal_baf)
      )
  )

  return(job_report)


}




#' Prepare ASCAT data from PCF select BAM files
#' 
#' @param tumour [REQUIRED] Path to tumour BAM file.
#' @param normal [REQUIRED] Path to normal BAM file.
#' @param version [OPTIONAL] PCF Select panel version. Default V3.
#' @param ref_dataset [OPTIONAL] Reference dataset for panel allele and loci selection. Default battenberg
#' @param gender [OPTIONAL] Sample gender. Default XY.
#' @param genome_version [OPTIONAL] Genome reference version. Default HG19.
#' @param ascat_ref [OPTIONAL] List with default references.
#' @param gamma [OPTIONAL] Gamma parameter. Default 1.
#' @param penalty [OPTIONAL] Penalty parameter. Default 25.
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


process_ascat=function(
  tumour_log2=NULL,
  normal_log2=NULL,
  tumour_baf=NULL,
  normal_baf=NULL,
  ref_dataset="battenberg",
  gender="XY",
  genome_version="HG19",
  gamma=1,
  penalty=25,
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


  if(is.null(tumour_log2)|is.null(normal_log2)){
    stop("A pair of tumour and normal LogR files are required. These can be prepared using the prepare_ascat function")
  }

  if(is.null(tumour_baf)|is.null(normal_baf)){
    stop("A pair of tumour and normal BAF files are required. These can be prepared using the prepare_ascat function")
  }


  
  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir)
  job=build_job(executor_id=executor_id,task=task_id)

  tumour_id=get_file_name(tumour_log2)
  normal_id=get_file_name(normal_log2)



  

  gc=ascat_ref[[toupper(genome_version)]]$panel[[paste0("PCF_",toupper(panel_version))]]$ascat[[ref_dataset]]$gc
  rt=ascat_ref[[toupper(genome_version)]]$panel[[paste0("PCF_",toupper(panel_version))]]$ascat[[ref_dataset]]$rt

  
  rdata=paste0(tumour_id,'.ASCAT_objects.Rdata')

  exec_code=paste0("Rscript -e \"","options(warn = -1);setwd(\\\"",
            out_file_dir,"\\\");library(ASCAT,quietly=TRUE);",
            "ascat.bc = ASCAT::ascat.loadData(Tumor_LogR_file =\\\"",tumour_log2,"\\\"",
            ",Tumor_BAF_file =\\\"",tumour_baf,"\\\"",
            ",Germline_BAF_file =\\\"",normal_baf,"\\\"",
            ",Germline_LogR_file =\\\"",normal_log2,"\\\"",
            ",gender=\\\"",toupper(gender),"\\\"",
            ",genomeVersion=\\\"",tolower(genome_version),"\\\",isTargetedSeq=T);",
            "ASCAT::ascat.plotRawData(ascat.bc,img.prefix = \\\"Before_correction_\\\");",
            "ascat.bc = ULPwgs::correctLog2_ascat(ascat.bc, GCcontentfile =\\\"",gc,"\\\"",
            ",replictimingfile =\\\"",rt,"\\\",threads=",threads,");",
            "ASCAT::ascat.plotRawData(ascat.bc,img.prefix = \\\"After_correction_\\\");",
            "ascat.bc = ASCAT::ascat.aspcf(ascat.bc,penalty=",penalty,");",
            "ascat.output = ASCAT::ascat.runAscat(ascat.bc, gamma=",gamma,", write_segments = T);",
            "ASCAT::ascat.plotSegmentedData(ascat.bc);",
            "QC = ASCAT::ascat.metrics(ascat.bc,ascat.output);",
            "save(ascat.bc, ascat.output, QC, file =\\\"",rdata,"\\\")\""
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
        plots=list(
          before_corr_normal=paste0(out_file_dir,"/Before_correction_",tumour_id,".germline.png"),
          after_corr_normal=paste0(out_file_dir,"/After_correction_",tumour_id,".germline.png"),
          before_corr_tumour=paste0(out_file_dir,"/Before_correction_",tumour_id,".tumour.png"),
          after_corr_tumour=paste0(out_file_dir,"/After_correction_",tumour_id,".tumour.png"),
          aspcf=paste0(out_file_dir,"/",tumour_id,".ASPCF.png"),
          ascat_profile=paste0(out_file_dir,"/",tumour_id,".ASCATfprofile.png"),
          raw_profile=paste0(out_file_dir,"/",tumour_id,".rawprofile.png"),
          sunrise=paste0(out_file_dir,"/",tumour_id,".sunrise.png")
        ),
        data=list(
          rdata=paste0(out_file_dir,"/",rdata),
          tumour_log2_pcfed=paste0(out_file_dir,"/",tumour_id,".LogR.PCFed.txt"),
          tumour_baf_pcfed=paste0(out_file_dir,"/",tumour_id,".BAF.PCFed.txt"),
          normal_log2_pcfed=paste0(out_file_dir,"/",normal_id,".LogR.PCFed.txt"),
          normal_baf_pcfed=paste0(out_file_dir,"/",normal_id,".BAF.PCFed.txt"),
          segments=paste0(out_file_dir,"/",tumour_id,".segments.txt"),
          segments_raw=paste0(out_file_dir,"/",tumour_id,".segments_raw.txt")
        )
      )
  )

  return(job_report)
}




#' @title ascat.correctLogR
#' @description Corrects logR of the tumour sample(s) with genomic GC content (replication timing is optional)
#' @param ASCATobj an ASCAT object
#' @param GCcontentfile File containing the GC content around every SNP for increasing window sizes
#' @param replictimingfile File containing replication timing at every SNP for various cell lines (optional)
#' @details Note that probes not present in the GC content file will be lost from the results
#' @return ASCAT object with corrected tumour logR
#'
#' @export
#' 
#' 

correctLog2_ascat= function(ASCATobj, GCcontentfile = NULL, replictimingfile = NULL,threads=4) {

  readCorrectionFile=function(correction_file,threads=4) {
  if (file.exists(correction_file) && file.info(correction_file)$size>0) {
    return(data.frame(data.table::fread(correction_file,sep='\t',showProgress=F,header=T,nThread = threads),
    row.names=1,check.names=F,stringsAsFactors=F))
  } else {
    return(NULL)
  }
  }




  if (is.null(GCcontentfile)) {
    stop("No GC content file given!")
  } else {
    GC_newlist=readCorrectionFile(GCcontentfile,threads=threads)
    stopifnot(!is.null(GC_newlist))
    ovl=intersect(rownames(ASCATobj$Tumor_LogR),rownames(GC_newlist))
    stopifnot(length(ovl)>nrow(ASCATobj$Tumor_LogR)/10)
    GC_newlist=GC_newlist[ovl,]
    if (!is.null(replictimingfile)) {
      replic_newlist=readCorrectionFile(replictimingfile,threads=threads)
      stopifnot(!is.null(replic_newlist))
      stopifnot(all(ovl %in% rownames(replic_newlist)))
      replic_newlist=replic_newlist[ovl,]
    } else {
      print.noquote("Warning: no replication timing file given, proceeding with GC correction only!")
    }
    
    SNPpos = ASCATobj$SNPpos[ovl,]
    Tumor_LogR = ASCATobj$Tumor_LogR[ovl,,drop=F]
    Tumor_BAF = ASCATobj$Tumor_BAF[ovl,,drop=F]
    
    chrs = intersect(ASCATobj$chrs,unique(SNPpos[,1]))
    
    Germline_LogR = NULL
    Germline_BAF = NULL
    if(!is.null(ASCATobj$Germline_LogR)) {
      Germline_LogR = ASCATobj$Germline_LogR[ovl,,drop=F]
      Germline_BAF = ASCATobj$Germline_BAF[ovl,,drop=F]
    }
    
    last = 0;
    ch = list();
    for (i in 1:length(ASCATobj$chrs)) {
      chrke = SNPpos[SNPpos[,1]==ASCATobj$chrs[i],]
      chrpos = chrke[,2]
      names(chrpos) = rownames(chrke)
      chrpos = sort(chrpos)
      ch[[i]] = (last+1):(last+length(chrpos))
      last = last+length(chrpos)
    }
    
    GC_correction_before=c()
    GC_correction_after=c()
    RT_correction_before=c()
    RT_correction_after=c()
    AUTOSOMES=!GC_newlist$Chr %in% ASCATobj$sexchromosomes
    
    for (s in 1:length(ASCATobj$samples)) {
      print.noquote(paste("Sample ", ASCATobj$samples[s], " (",s,"/",length(ASCATobj$samples),")",sep=""))
      Tumordata = Tumor_LogR[,s]
      names(Tumordata) = rownames(Tumor_LogR)
      
      # Calculate correlation (explicit weighting now automatically done)
      corr = abs(cor(GC_newlist[AUTOSOMES, 3:ncol(GC_newlist)], Tumordata[AUTOSOMES], use="complete.obs")[,1])
      
      index_1kb = grep(pattern = "X?1([06]00bp|kb)", x = names(corr))
      maxGCcol_insert = names(which.max(corr[1:index_1kb]))
      # if no replication timing data, expand large windows to 500kb
      if (!is.null(replictimingfile)) {
        index_max = grep(pattern = "X?10((00|24)00bp|0kb)", x = names(corr))
      } else {
        index_max = grep(pattern = "X?1M", x = names(corr)) - 1
      }
      # start large window sizes at 5kb rather than 2kb to avoid overly correlated expl variables
      maxGCcol_amplic = names(which.max(corr[(index_1kb+2):index_max]))
      
      cat("GC correlation: ",paste(names(corr),format(corr,digits=2), ";"),"\n")   
      cat("Short window size: ",maxGCcol_insert,"\n")
      cat("Long window size: ",maxGCcol_amplic,"\n")
      
      corrdata = data.frame(logr = Tumordata,
                            GC_insert = GC_newlist[,maxGCcol_insert],
                            GC_amplic = GC_newlist[,maxGCcol_amplic])
      
      if (!is.null(replictimingfile)) {
        corr_rep = abs(cor(replic_newlist[AUTOSOMES, 3:ncol(replic_newlist)], Tumordata[AUTOSOMES], use="complete.obs")[,1])
        maxreplic = names(which.max(corr_rep))
        
        cat("Replication timing correlation: ",paste(names(corr_rep),format(corr_rep,digits=2), ";"),"\n") 
        cat("Replication dataset: " ,maxreplic,"\n")
        
        # Multiple regression 
        corrdata$replic <- replic_newlist[, maxreplic]
        model = lm(logr ~ splines::ns(x = GC_insert, df = 5, intercept = T) + splines::ns(x = GC_amplic, df = 5, intercept = T) + splines::ns(x = replic, df = 5, intercept = T), y=F, model = F, data = corrdata, na.action="na.exclude")
        Tumor_LogR[,s] = residuals(model)
        
        RT_correction_before=c(RT_correction_before,paste0(maxreplic,"=",round(corr_rep[maxreplic],4)))
        corr_rep_after = abs(cor(replic_newlist[, which(colnames(replic_newlist)==maxreplic),drop=F], Tumor_LogR[,s], use="complete.obs")[,1])
        RT_correction_after=c(RT_correction_after,paste0(maxreplic,"=",round(corr_rep_after[maxreplic],4)))
      } else {
        model = lm(logr ~ splines::ns(x = GC_insert, df = 5, intercept = T) + splines::ns(x = GC_amplic, df = 5, intercept = T), y=F, model = F, data = corrdata, na.action="na.exclude")
        Tumor_LogR[,s] = residuals(model)
        RT_correction_before=c(RT_correction_before,NA)
        RT_correction_after=c(RT_correction_after,NA)
      }
      
      if (!ASCATobj$isTargetedSeq) {
        chr=split_genome(SNPpos)
      } else {
        chr=ch
      }
      
      GC_correction_before=c(GC_correction_before,paste0(maxGCcol_insert,"=",round(corr[maxGCcol_insert],4),' / ',maxGCcol_amplic,'=',round(corr[maxGCcol_amplic],4)))
      corr_after=abs(cor(GC_newlist[, which(colnames(GC_newlist) %in% c(maxGCcol_insert,maxGCcol_amplic))], Tumor_LogR[,s], use="complete.obs")[,1])
      GC_correction_after=c(GC_correction_after,paste0(maxGCcol_insert,"=",round(corr_after[maxGCcol_insert],4),' / ',maxGCcol_amplic,'=',round(corr_after[maxGCcol_amplic],4)))
    }
    
    names(GC_correction_before)=colnames(Tumor_LogR)
    names(GC_correction_after)=colnames(Tumor_LogR)
    names(RT_correction_before)=colnames(Tumor_LogR)
    names(RT_correction_after)=colnames(Tumor_LogR)
    
    ASCATobj$Tumor_LogR=Tumor_LogR
    ASCATobj$Tumor_BAF=Tumor_BAF
    ASCATobj$Germline_LogR=Germline_LogR
    ASCATobj$Germline_BAF=Germline_BAF
    ASCATobj$Tumor_LogR_segmented=NULL
    ASCATobj$Tumor_BAF_segmented=NULL
    ASCATobj$SNPpos=SNPpos
    ASCATobj$ch=ch
    ASCATobj$chr=chr
    ASCATobj$chrs=chrs
    ASCATobj$GC_correction_before=GC_correction_before
    ASCATobj$GC_correction_after=GC_correction_after
    ASCATobj$RT_correction_before=RT_correction_before
    ASCATobj$RT_correction_after=RT_correction_after
    return(ASCATobj)
  }
}







