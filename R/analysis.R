
#' Merge BAM files in directory
#'
#' This function takes a BAM file and merges it with others found within a
#' directory. This function is still Work In Progress (WIP).
#'
#' @param bam Path to the input bam file with the sequence.
#' @param bin_samtools Path to samtools binary.
#' @param bam_dir Path to directory with BAM files to merge.
#' @param threads Number of threads to use.Default 3.
#' @param output_name Output file name
#' @param verbose Enables progress messages. Default False.
#' @export


merge_bams_samtools=function(
bin_samtools=build_default_tool_binary_list()$bin_samtools,
bams="",output_name="",
verbose=TRUE,threads=3,ram=4,
executor=make_unique_id("mergeBAMs"),
task="mergeBAMs",mode="local",
time="48:0:0",update_time=60,wait=FALSE,hold=NULL){

    out_file_dir=set_dir(dir=output_dir,name="merged")
    exec_code=paste(bin_samtools,"merge ",paste0(out_file_dir,"/",output_name,".bam"), " --threads",
      threads,paste(bams,collapse=" "))
      
    if(verbose){
      print_verbose(exec_code=exec_code)
    }
    system(exec_code)
}



#' Concatenate BAM files in directory.
#'
#' This function takes a BAM file and merges it with others found within a
#' directory. The ouput file is a non-sorted bam file
#'
#' @param bam Path to the input bam file with the sequence.
#' @param bam_dir Path to directory with BAM files to merge.
#' @param threads Number of threads to use.Default 3.
#' @param output_name Output file name
#' @param verbose Enables progress messages. Default False.
#' @export

concatenate_bams=function(
  bin_samtools=build_default_tool_binary_list()$bin_samtools,
  bams="",output_name="",verbose=FALSE,
  batch_config=build_default_preprocess_config(),threads=3){

    exec_code=paste(bin_samtools,"cat -o",paste0(output_name,".bam"), " --threads",
      threads,paste(bams,collapse=" "))
    if(verbose){
      print_verbose(exec_code=exec_code)
    }
    system(exec_code)
}










#' Generate quality control metrics for aligned sequences
#'
#' This function generates quality control metrics for an aligned sequence. Works for WGS and Panel data.
#' WGS metrics will be generated if bait and target intervals are nor given, otherwise Panel metrics will be generateds.
#' Target and bait BEDs can be converted to interval format using Picards BedToIntervalList functions.
#' For example java -jar ~/tools/picard/build/libs/picard.jar BedToIntervalList I=PCFv2newchem_primary_targets.bed O=PCFv2newchem_primary_targets.interval_list SD=~/Scratch/RefGenome/hs37d5.fa
#' For more information about interval format check: https://gatk.broadinstitute.org/hc/en-us/articles/360036883931-BedToIntervalList-Picard-
#'
#' Off target BED can be created by generating the complementary regions from the target BED using bedtools complement function.
#' For example ~/tools/bedtools2/bin/bedtools complement -i PCFv2newchem_capture_targets.bed -g ~/Scratch/RefGenome/hs37d5.fa.fai > Probes_Off_target_regions.bed
#' Note: This BED has to be sorted using the same reference as the BAM files.
#' Using different reference to generate the BED from the one used to align the BAM may cause issues even when sorted due to scaffold chromosomes.
#'
#'
#' @param bam Path to the BAM file.
#' @param bin_samtools Path to samtools executable. Default path tools/samtools/samtools.
#' @param bin_picard Path to picard executable. Default path tools/picard/build/libs/picard.jar.
#' @param bin_bedtools Path to bedtools executable. Default path tools/bedtools2/bin/bedtools. Only required if analyzing panel data.
#' @param ref_genome Path to input file with the reference genome sequence.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param executor_id Task EXECUTOR ID. Default "gatherBAM"
#' @param task_name Task name. Default "gatherBAM"
#' @param thread Threads to use. Default 3.
#' @param ram RAM in GB to use. Default 4 Gb.
#' @param tmp_dir Path to tmp directory.
#' @param mapq Minimum MapQ for Picard Wgs metrics. Default 0.
#' @param ri Path to ribosomal intervals file. Only for RNAseq.
#' @param ref_flat Path to flat refrence. Only for RNAseq.
#' @param bi Bait capture target interval for panel data. Requires ti and off_target and on_tar arguments. Interval format.
#' @param ti Primary target intervals for panel data. Requires bi and off_tar and on_tar argmunets. Interval format.
#' @param mode Type of data to generate metrics for. Default tg. Options ["wgs","tg","rna"]
#' @export

metrics_alignqc=function(
  bin_samtools=build_default_tool_binary_list()$bin_samtools,
  bin_picard=build_default_tool_binary_list()$bin_picard,
  bin_bedtools=build_default_tool_binary_list()$bin_bedtools,
  ref_genome=build_default_reference_list()$HG19$reference$genome,
  bam="",output_dir=".",
  verbose=FALSE,batch_config=build_default_preprocess_config(),
  tmp_dir=".",mapq=0,
  bi=build_default_reference_list()$HG19$panel$PCF_V3$intervals$bi,
  ti=build_default_reference_list()$HG19$panel$PCF_V3$intervals$ti,
  ri=build_default_reference_list()$HG19$rnaseq$intervals$ri,
  ref_flat=build_default_reference_list()$HG19$rnaseq$reference$ref_flat,
  method="tg",
  mode="local",executor_id=make_unique_id("alignQC"),
  task_name="alignQC",time="48:0:0",
  threads=3,ram=4,update_time=60,wait=FALSE, hold=NULL){
    

    argg <- as.list(environment())
    task_id=make_unique_id(task_name)

    out_file_dir=set_dir(dir=output_dir,name="alignqc_report")
    
    job=build_job(executor_id=executor_id,task_id=task_id)

    ## Generate alignment metrics
    job_report=build_job_report(
      job_id=job,
      executor_id=executor_id, 
      task_id=task_id,
      input_args=argg,
      out_file_dir=out_file_dir,
      out_files=list(
        )
      )
  job_report[["steps"]][["mapq_metrics"]] <- mapq_metrics_bam_samtools(
   bin_samtools=bin_samtools,
   bam=bam,output_dir=out_file_dir,
   executor_id=task_id,
   verbose=verbose,time=time,mode=mode,
   threads=threads,ram=ram,
   update_time=update_time,
   batch_config=batch_config,
   wait=FALSE, hold=hold)

  job_report[["steps"]][["summary_metrics"]] <- summary_metrics_bam_picard(
    bin_picard=bin_picard,
    bam=bam,output_dir=out_file_dir,
    verbose=verbose,tmp_dir=tmp_dir,
    mode=mode,executor_id=task_id,
    threads=threads,
    batch_config=batch_config,
    ram=ram,update_time=update_time,
    wait=FALSE, hold=hold)
    
  job_report[["steps"]][["insertsize_metrics"]] <- insertsize_metrics_bam_picard(
  bin_picard=bin_picard,bam=bam,
  output_dir=out_file_dir,
   verbose=verbose,tmp_dir=tmp_dir,
   mode=mode,executor_id=task_id,
   threads=threads,ram=ram,
   batch_config=batch_config,
   update_time=update_time,
   wait=FALSE,hold=hold)


    ## Only call metrics for panel data if bait and target intervals are supplied, Otherwise WGS metrics.
    
    if (method=="tg"){
   
     job_report[["steps"]][["tg_summary_metrics"]] <- tg_summary_metrics_bam_picard(
      bin_picard=bin_picard,
      bam=bam,output_dir=out_file_dir,
      verbose=verbose,
      tmp_dir=tmp_dir,bi=bi,
      ti=ti,mode=mode,
      batch_config=batch_config,
      executor_id=task_id,
      threads=threads,ram=ram,
      update_time=update_time,
      wait=FALSE,hold=hold)
      # ## Picard doesn't output coverage stats for off-target regions therefore we have to estimate this manually.

      # ## I use this function to get the mean coverage on and off target.

      # ## For target regions

      # bed_coverage(bin_path=bin_path3,bam=bam,bed=on_tar,verbose=verbose,sorted=TRUE,
      #   mean=TRUE,hist=TRUE,fai=paste0(ref_genome,".fai"),suffix="on_Target",output_dir=out_file_dir)
      # bed_coverage(bin_path=bin_path3,bam=bam,bed=on_tar,verbose=verbose,sorted=TRUE,
      #   mean=TRUE,hist=FALSE,fai=paste0(ref_genome,".fai"),suffix="on_Target",output_dir=out_file_dir)

      # ## For off target regions

      # bed_coverage(bin_path=bin_path3,bam=bam,bed=off_tar,verbose=verbose,sorted=TRUE,
      #   mean=TRUE,hist=TRUE,fai=paste0(ref_genome,".fai"),suffix="off_Target",output_dir=out_file_dir)
      # bed_coverage(bin_path=bin_path3,bam=bam,bed=off_tar,verbose=verbose,sorted=TRUE,
      #   mean=TRUE,hist=FALSE,fai=paste0(ref_genome,".fai"),suffix="off_Target",output_dir=out_file_dir)

      # ## Generate violin and cummulative plots for target and off target regions

      # plot_coverage_panel(on_target=paste0(out_file_dir,"/coverage/",get_file_name(bam),".on_Target.Per_Region_Coverage.txt"),
      # off_target=paste0(out_file_dir,"/coverage/",get_file_name(bam),col=c(5,4),height=6,width=12,output_dir=out_file_dir))
      # plot_cumulative_cov(on_target=paste0(out_file,".on_Target.Histogram_Coverage.txt"),
      # off_target=paste0(out_file,".off_Target.Histogram_Coverage.txt"),height=6,width=12,output_dir=out_file_dir)

    }else if(method=="rna"){
        job_report[["steps"]][["rna_summary_metrics"]] <- rnaseq_summary_metrics_bam_picard(
        bin_picard=bin_picard,
        bam=bam,output_dir=out_file_dir,
        verbose=verbose,tmp_dir=tmp_dir,
        ri=ri,ref_flat=ref_flat,
        mode=mode,executor_id=task_id,
        threads=threads,ram=ram,
        batch_config=batch_config,
        update_time=update_time,
        wait=FALSE,hold=hold)
    }else if(method=="wgs"){
        job_report[["steps"]][["rna_summary_metrics"]] <- wgs_summary_metrics_bam_picard(
        bin_picard=bin_picard,
        bam=bam,
        output_dir=out_file_dir,
        verbose=verbose,
        ref_genome=ref_genome,
        batch_config=batch_config,
        tmp_dir=tmp_dir,mode=mode,
        executor_id=task_id,
        threads=threads,ram=ram,
        update_time=update_time,
        wait=FALSE,hold=hold)
    }

  }



#' Filter BAM to specific regions
#'
#' This function takes a BAM file and filter it to genomic coordinate or otherwise BED file with
#' multiple genomic coordinates. The output is a BAM file with reads for inputed genomic region/s
#'
#'
#' @param bin_samtools Path to readCounter executable. Default path tools/samtools/samtools.
#' @param bam Path to the BAM file .
#' @param position String of genomic position to filter. Ex chr6:1000-100000
#' @param output_name Name of the output file.
#' @param output_dir Path to output directory.
#' @param bed Size of non overlaping windows. Default 500000.
#' @param threads Number of threads to use. Default 1
#' @param verbose Enables progress messages. Default False.
#' @export

filter_bam=function(
  bin_samtools=build_default_tool_binary_list()$bin_samtools,
  bam="",position="",bed="",verbose=FALSE,
  batch_config=build_default_preprocess_config(),
  output_dir=".",threads=1){


  out_file_dir=set_dir(dir=output_dir,name="filtered")
  if (position!="" &bed!=""){
    cat("Position and bed arguments are mutually exclusive")
    quit()
  }

  if(bed!=""){
    bed=paste("-L",bed)
  }
  out_file=paste0(out_file_dir,"/",get_file_name(bam),".filtered.bam")

  exec_code=paste(bin_samtools,"view -b",bed,"-@",threads,bam,position,">",
    out_file
  )

  if (verbose){
    print_verbose(exec_code=exec_code)
  }
  system(exec_code)
}


#' Generate report for ULP-WGS samples
#'
#' This function generates a report that helps to segment the genome, predict large-scale
#' copy number alterations, and estimate tumor fraction in ULP-WGS samples.
#'
#' @param wig Path to the WIG file.
#' @param norm_wig Path to normal WIG file.
#' @param bin_ichor Path to ichorCNA executable. Default path tools/ichorCNA/scripts/runIchorCNA.R.
#' @param output_dir Path to the output directory.
#' @param bed Path to BED file with target regions.
#' @param sample_id String with sample name.
#' @param ploidy Initial tumour ploidy. Default 2,3
#' @param tumour_content Initial normal contamination. Default 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9
#' @param homozygous_del Include Homozygous deleteions. Default FALSE.
#' @param gc Path to GC-content WIG with . Default tools/ichorCNA/inst/extdata/gc_hg19_500kb.wig
#' @param map Path to mappability score WIG with GC content. Default tools/ichorCNA/inst/extdata/map_hg19_500kb.wig
#' @param centromere Path to file containing centromere locations. Default tools/ichorCNA/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt
#' @param normal_panel Path to file containing normal panel data. Default tools/ichorCNA/inst/extdata/HD_ULP_PoN_hg38_500kb_median_normAutosome_median.rds
#' @param subclonal_states Subclonal states to consider. Default NULL
#' @param male_tresh Minimum percentage of reads in chromosome Y to call male. Default 0.0001
#' @param libdir Path to dir containing ichorCNA libraries.
#' @param chrs Chromosomes to analyze. Default 'c(1:22,\"X\")'
#' @param chrsTrain Chromosomes to estimate parameters. Default 'c(1:22)'
#' @param verbose Enables progress messages. Default False.
#' @export

ichorCNA=function(bin_ichor=build_default_tool_binary_list()$bin_ichor,
sample_id="",
wig="",norm_wig="",bed="",ploidy="2,3",
tumour_content="0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9",
homozygous_del="False",subclonal_states=NULL,
gc="tools/ichorCNA/inst/extdata/gc_hg19_500kb.wig",
map="tools/ichorCNA/inst/extdata/map_hg19_500kb.wig",
centromere="tools/ichorCNA/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt",
normal_panel="tools/ichorCNA/inst/extdata/HD_ULP_PoN_500kb_median_normAutosome_mapScoreFiltered_median.rds",
output_dir=".",verbose=TRUE,libdir="tools/ichorCNA",male_tresh=0.0001,chrs="'c(1:22,\"X\")'",chrTrain="'c(1:22)'"){


    out_file_dir=set_dir(dir=output_dir,name="ichorcna")

   
    if(!sample_id==""){
      sample_name=get_file_name(wig)
    }else{
      sample_id=sample_name
    }


    if (norm_wig!=""){
      norm_wig=paste0("--NORMWIG ",norm_wig)
    }

    if (bed!=""){
      bed=paste0("--exons.bed ",bed)
    }
    if (!is.null(subclonal_states)){
      subclonal_states=paste("--scStates",paste0("'c(",subclonal_states,")'"))
    }else{
      subclonal_states=""
    }

    exec_code=paste("Rscript",bin_ichor,"--id",sample_id,"--WIG",wig,norm_wig,bed,
      "--ploidy",paste0("'c(",ploidy,")'"),"--normal",paste0("'c(",tumour_content,")'"),
      "--maxCN 7 --gcWig", gc,"--mapWig",map,"--centromere",centromere,subclonal_states,
      "--normalPanel",normal_panel,"--includeHOMD",homozygous_del,"--chrs",chrs,
      "--fracReadsInChrYForMale",male_tresh,"--chrTrain",chrTrain,
      "--estimateNormal True --estimatePloidy True --estimateScPrevalence True --outDir",out_file_dir,"--libdir",libdir)

    if(verbose){
      print_verbose(exec_code=exec_code)
    }
    system(exec_code)
  }


#' Generate panel of normals for ichorCNA
#'
#' This function generates a report that helps to segment the genome, predict large-scale
#' copy number alterations, and estimate tumor fraction in ULP-WGS samples.
#'
#' @param wigs Path to the WIG files.
#' @param wigs_dir Path to dir with WIG files.
#' @param bin_ichor_pon Path to ichorCNA executable. Default path tools/ichorCNA/scripts/createPanelOfNormals.
#' @param output_dir Path to the output directory.
#' @param bed Path to BED file with target regions.
#' @param output_name File output name.
#' @param gc Path to GC-content WIG with . Default tools/ichorCNA/inst/extdata/gc_hg19_500kb.wig
#' @param map Path to mappability score WIG with GC content. Default tools/ichorCNA/inst/extdata/map_hg19_500kb.wig
#' @param centromere Path to file containing centromere locations. Default tools/ichorCNA/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt
#' @param male_tresh Minimum percentage of reads in chromosome Y to call male. Default 0.0001
#' @param verbose Enables progress messages. Default False.
#' @export

panel_of_normals_ichorCNA=function(
  bin_ichor_pon=build_default_tool_binary_list()$bin_ichor_pon,
  wigs_dir="",wigs="",bed="",output_name="PoN_ichorCNA",
  gc="tools/ichorCNA/inst/extdata/gc_hg19_500kb.wig",
  map="tools/ichorCNA/inst/extdata/map_hg19_500kb.wig",
  centromere="tools/ichorCNA/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt",
  male_tresh=0.0001,verbose=TRUE){

    if (wigs!="" && wigs_dir!=""){
      cat("Arguments wigs and wigs_dir are mutually exclusive")
      quit()
    }

    if (bed!=""){
      wig_list=paste0("ichorCNAnormalList",ULPwgs::get_file_name(bed),".tmp")
    }else{
      wig_list=paste0("ichorCNAnormalList.tmp")
    }

    if (wigs!=""){
          writeLines(wigs,file(wig_list))
    }

    if (wigs_dir!=""){
      files=list.files(wigs_dir,pattern=".wig",recursive=TRUE,full.names=TRUE)
      writeLines(wigs,file(wig_list))
    }

    if (bed!=""){
      bed=paste0("--exons.bed ",bed)
    }

    exec_code=paste("Rscript",bin_ichor_pon,"--filelist",wig_list,bed," --gcWig", gc,
      "--mapWig",map,"--centromere",centromere,"--fracReadsInChrYForMale",male_tresh," --outfile",output_name)
    if(verbose){
      print_verbose(exec_code=exec_code)
    }
    system(exec_code)
    system(paste("rm",wig_list))

}
