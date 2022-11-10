#' Process ASCAT data using REFPHASE
#' 
#' @param tumour [REQUIRED] Path to tumour BAM file/s.
#' @param normal [REQUIRED] Path to normal BAM file/s.
#' @param homozygous_cutoff  [OPTIONAL] Cut-off to consider homozygous snps . Default 0.7. 
#' @param patient_id [REQUIRED] Path to normal BAM file.
#' @param center_baf [OPTIONAL] Re-center BAF values. Default TRUE
#' @param fit_log2 [OPTIONAL] Better fit ASCAT log2 ratios. Default TRUE
#' @param panel_version [OPTIONAL] PCF Select panel version. Default V3.
#' @param ref_dataset [OPTIONAL] Reference dataset for panel allele and loci selection. Default battenberg
#' @param gender [OPTIONAL] Sample gender. Default XY.
#' @param gamma [OPTIONAL] Gamma parameter. Default 1.
#' @param penalty [OPTIONAL] Penalty parameter. Default 25.
#' @param genome_version [OPTIONAL] Genome reference version. Default HG19.
#' @param ascat_ref [OPTIONAL] List with default references.
#' @param batch_config [OPTIONAL] List with default references.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param threads [OPTIONAL] Number of threads to split the work. Default 4
#' @param ram [OPTIONAL] RAM memory to asing to each thread. Default 4
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Task EXECUTOR ID. Default "processRefphase"
#' @param task_name Task name. Default "processRefphase"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export


call_refphase=function(
    bin_allele_counter=build_default_tool_binary_list()$bin_allele_counter,
    ascat_ref=build_default_reference_list(),
    ref_dataset="battenberg",
    gender="XY",
    genome_version="HG19",
    panel_version="V3",
    tumour=NULL,
    normal=NULL,
    patient_id=NULL,
    homozygous_cutoff=0.7,
    center_baf=TRUE,
    fit_log2=TRUE,
    gamma=1,
    penalty=25,
    output_dir=".",
    verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=8,
    ram=4,mode="local",
    executor_id=make_unique_id("callRefphase"),
    task_name="callRefphase",time="48:0:0",
    update_time=60,
    wait=FALSE,
    hold=NULL
){

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir)
  job=build_job(executor_id=executor_id,task=task_id)

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


  jobs_report[["steps"]][["parallel_call_ascat"]]<-parallel_samples_call_ascat(
    bin_allele_counter=bin_allele_counter,
    tumour=tumour,
    normal=normal,
    patient_id=patient_id,
    ref_dataset=ref_dataset,
    gender=gender,
    genome_version=genome_version,
    panel_version=panel_version,
    output_dir=out_file_dir,
    gamma=gamma,
    penalty=penalty,
    verbose=verbose,
    ascat_ref=ascat_ref,
    batch_config=batch_config,
    threads=threads,
    ram=ram,mode=mode,
    executor_id=task_id,
    time=time,
    hold=hold
  )
 

  jobs_report[["steps"]][["process_refphase"]]<-process_refphase(
      ascat_rdata=normalizePath(unlist(unlist_lvl(jobs_report[["steps"]][["parallel_call_ascat"]],var="rdata"))),
      patient_id=patient_id,
      homozygous_cutoff=homozygous_cutoff,
      center_baf=center_baf,
      fit_log2=fit_log2,
      output_dir=paste0(out_file_dir,"/",patient_id),
      verbose=verbose,
      batch_config=batch_config,
      threads=threads,
      ram=ram,mode=mode,
      executor_id=task_id,
      time=time,
      hold=unlist_lvl(jobs_report[["steps"]][["parallel_call_ascat"]],var="job_id")

  )

   return(jobs_report)
}


#' Process ASCAT data using REFPHASE
#' 
#' @param ascat_data [REQUIRED] Path to allele_counter binary.
#' @param homozygous_cutoff  [OPTIONAL] Cut-off to consider homozygous snps . Default 0.7. 
#' @param patient_id [REQUIRED] Path to normal BAM file.
#' @param center_baf [OPTIONAL] Re-center BAF values. Default TRUE
#' @param fit_log2 [OPTIONAL] Better fit ASCAT log2 ratios. Default TRUE
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param threads [OPTIONAL] Number of threads to split the work. Default 4
#' @param ram [OPTIONAL] RAM memory to asing to each thread. Default 4
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Task EXECUTOR ID. Default "processRefphase"
#' @param task_name Task name. Default "processRefphase"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export


process_refphase=function(
  ascat_rdata=NULL,
  patient_id=NULL,
  homozygous_cutoff=0.7,
  center_baf=TRUE,
  fit_log2=TRUE,
  output_dir=".",
  verbose=FALSE,
  batch_config=build_default_preprocess_config(),
  threads=8,
  ram=4,mode="local",
  executor_id=make_unique_id("processRefphase"),
  task_name="processRefphase",time="48:0:0",
  update_time=60,
  wait=FALSE,
  hold=NULL
){

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir,name="rephase_reports")
  job=build_job(executor_id=executor_id,task=task_id)

  rdata=paste0(out_file_dir,"/",patient_id,".refphase.RData")


  exec_code=paste0("Rscript -e \"",
    "options(warn = -1);",
    "library(ASCAT,quietly=TRUE);library(refphase,quietly=TRUE);",
    "refphase_input <-ULPwgs::read_and_load_ascat_refphase(ascat_rdata=unlist(strsplit(split=\\\",\\\",\\\"",
    paste(ascat_rdata,collapse = ","),"\\\")),homozygous_cutoff=",homozygous_cutoff,");",
    ifelse(center_baf,"refphase_input <- center_baf(refphase_input);",""),
    ifelse(fit_log2,"refphase_input <- fit_logr_to_ascat(refphase_input);",""),
    "results <- refphase(refphase_input);",
    "results$name <-\\\"",patient_id,"\\\";",
    "refphase::export(refobj=results,output_folder=\\\"",out_file_dir,"\\\");",
    "save(refphase_input, results, file =\\\"",rdata,"\\\")\""
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
      stop("refphase failed to run due to unknown error.
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
        data=list(
          rdata=paste0(out_file_dir,"/",rdata),
          phased_segments=paste0(out_file_dir,"/",patient_id,"_phased_segments.tsv"),
          phased_snps=paste0(out_file_dir,"/",patient_id,"_phased_snps.tsv"),
          consensus_events=paste0(out_file_dir,"/",patient_id,"_consensus_events.tsv"),
          genome_proportions=paste0(out_file_dir,"/",patient_id,"_genome_proportions.tsv")
        )
      )
  )

  return(job_report)
}



#' Read and load ASCAT style data from rdata files
#' 
#' @param ascat_data [REQUIRED] Path to allele_counter binary.
#' @param homozygous_cutoff  [OPTIONAL] Cut-off to consider homozygous snps . Default 0.7. 
#' @export


read_and_load_ascat_refphase=function(ascat_rdata=NULL,homozygous_cutoff = 0.7){
  samples=Vectorize(get_file_name)(ascat_rdata)

    ascat_input <- list()
    ascat_output <- list()
    lapply(ascat_rdata,FUN=function(x){
      load(x)
      sample=get_file_name(x)
      ascat_input[[sample]] <<- ascat.bc
      ascat_output[[sample]] <<- ascat.output
    })

  modified_load_ascat_refphase <- function(samples, ascat_input, ascat_output, het_only = FALSE, homozygous_cutoff = 0.7) {
  




        # Validates loaded ASCAT input and output
      validate_ascat_data <- function(samples, ascat_input, ascat_output) {
        if (any(length(samples) != length(ascat_input) | length(samples) != length(ascat_output))) {
          stop("samples does not have same length as ascat_input or ascat_output")
        }
        if (any(samples != names(ascat_input) | samples != names(ascat_output))) {
          stop("samples differs from names(ascat_input) or names(ascat_output)")
        }
      }

    #Format chromosomes to fit numerical
    #Formatting such that chromosom 1 is "1". Chromosome X and Y are labelled as "X" and "Y", respectively
    
    format_chromosomes <- function(chromosomes) {
    chromosomes <- gsub("(23)|x", "X", chromosomes)
    chromosomes <- gsub("(24)|y", "Y", chromosomes)
    chromosomes <- gsub("(c|C)hr(om)?([0-9XYxy]{1,2})", "\\3", chromosomes)

    stopifnot(gsub("[0-9]{1,2}|[XYxy]", "", chromosomes) == "")

    chromosomes
    }

    format_snps <- function(snps, samples, het_only = FALSE) {
    for (sample in samples) {
        # Mask non-hz SNPs with NA
        GenomicRanges::mcols(snps[[sample]])$baf[GenomicRanges::mcols(snps[[sample]])$germline_zygosity != "het"] <- NA
        # Only keep heterozygous positions
        if (isTRUE(het_only)) {
        snps[[sample]] <- snps[[sample]][GenomicRanges::mcols(snps[[sample]])$germline_zygosity == "het"]
        }
    }

    snps <- GenomicRanges::GRangesList(snps, compress = FALSE)
    snps <- intersect_all_snps(snps)

    names(snps) <- samples

    snps
    }

    format_segs <- function(segs, samples) {
    if (any(segs$start > segs$end)) {
        stop("Error: Start position in segmentation must be <= end position")
    }

    segs <- GenomicRanges::GRangesList(segs, compress = FALSE)
    names(segs) <- samples

    segs
    }

    intersect_all_snps <- function(snps) {

    intersected_snps <- IRanges::subsetByOverlaps(snps[[1]], snps[[2]])
    if (length(snps) >= 3) {
        for (sample in names(snps)[3:length(snps)]) {
        intersected_snps <- IRanges::subsetByOverlaps(intersected_snps, snps[[sample]])
        }
    }
    for (sample in names(snps)) {
        snps[[sample]] <- IRanges::subsetByOverlaps(snps[[sample]], intersected_snps)
    }

    snps
    }
    sample_data_from_ploidy_purity <- function(samples, ploidy, purity) {
        sample_data <- data.frame("sample_id" = samples,
                                    "ploidy" = unname(unlist(ploidy[samples])),
                                    "purity" = unname(unlist(purity[samples])))
        sample_data$segmentation <- ""
        sample_data$snps <- ""

        sample_data
    }


    validate_ascat_data(samples, ascat_input, ascat_output)
    data <- list("segs" = list(), "snps" = list(), "purity" = list(), "ploidy" = list())




  for (sample in samples) {
    cur_ascat_input <- ascat_input[[sample]]
    cur_ascat_output <- ascat_output[[sample]]

    # 1) SNPs
    cur_chroms <- format_chromosomes(cur_ascat_input$SNPpos$Chromosome)
    cur_snps <- GPos(
      seqnames = Rle(cur_chroms),
      pos = cur_ascat_input$SNPpos$Position,
      baf = cur_ascat_input$Tumor_BAF[[1]],
      logr = cur_ascat_input$Tumor_LogR[[1]],
      germline_zygosity = ifelse(is.na(cur_ascat_input$Germline_BAF[[1]]) |
        abs(0.5 - cur_ascat_input$Germline_BAF[[1]]) >= (homozygous_cutoff / 2), "hom", "het")
    )

    data$snps[[sample]] <- cur_snps

    # 2) Segments
    cur_segs <- cur_ascat_output$segments
    names(cur_segs) <- c("sample_id", "chrom", "start", "end", "cn_major", "cn_minor")
    cur_segs$chrom <- format_chromosomes(cur_segs$chrom)
    cur_segs$sample_id <- sample

    cur_segs <- GenomicRanges::makeGRangesFromDataFrame(cur_segs,
      keep.extra.columns = TRUE,
      ignore.strand = TRUE,
      seqnames.field = "chrom",
      start.field = "start",
      end.field = "end",
      starts.in.df.are.0based = FALSE
    )

    data$segs[[sample]] <- cur_segs

    # 3) Purity and Ploidy
    data$purity[[sample]] <- cur_ascat_output$aberrantcellfraction[[1]]
    data$ploidy[[sample]] <- cur_ascat_output$ploidy[[1]]
  }

  
  data$segs <- format_segs(data$segs, samples)
  data$snps <- format_snps(data$snps, samples, het_only = het_only)

  data$sample_data <- sample_data_from_ploidy_purity(samples, data$ploidy, data$purity)
  class(data) <- "RefphaseInput"

  data
  }

  refphase_input_data=modified_load_ascat_refphase(samples,ascat_input,ascat_output)
 
  return(refphase_input_data)

}
