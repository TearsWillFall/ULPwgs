
#' Merge BAM files in directory
#'
#' This function takes a BAM file and merges it with others found within a
#' directory. This function is still Work In Progress (WIP).
#'
#' @param bam Path to the input bam file with the sequence.
#' @param bam_dir Path to directory with BAM files to merge.
#' @param threads Number of threads to use.Default 3.
#' @param output_name Output file name
#' @param verbose Enables progress messages. Default False.
#' @export


merge_bams_samtools=function(bin_path="tools/samtools/samtools",bams="",output_name="",
verbose=TRUE,threads=3,ram=4,executor=make_unique_id("mergeBAMs"), task="mergeBAMs",mode="local",
time="48:0:0",update_time=60,wait=FALSE,hold=""){
    out_file_dir=set_dir(dir=output_dir,name="merged")
    exec_code=paste(bin_path,"merge ",paste0(out_file_dir,"/",output_name,".bam"), " --threads",
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

concatenate_bams=function(bin_path="tools/samtools/samtools",bams="",output_name="",
verbose=FALSE,threads=3){
    exec_code=paste(bin_path,"cat -o",paste0(output_name,".bam"), " --threads",
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
#' @param bin_path Path to samtools executable. Default path tools/samtools/samtools.
#' @param bin_path2 Path to picard executable. Default path tools/picard/build/libs/picard.jar.
#' @param bin_path3 Path to bedtools executable. Default path tools/bedtools2/bin/bedtools. Only required if analyzing panel data.
#' @param ref_genome Path to input file with the reference genome sequence.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param ram RAM in GB to use. Default 4 Gb.
#' @param tmp_dir Path to tmp directory.
#' @param mapq Minimum MapQ for Picard Wgs metrics. Default 0.
#' @param ri Path to ribosomal intervals file. Only for RNAseq.
#' @param ref_flat Path to flat refrence. Only for RNAseq.
#' @param bi Bait capture target interval for panel data. Requires ti and off_target and on_tar arguments. Interval format.
#' @param ti Primary target intervals for panel data. Requires bi and off_tar and on_tar argmunets. Interval format.
#' @param mode Type of data to generate metrics for. Default tg. Options ["wgs","tg","rna"]
#' @export

align_qc_metrics=function(bin_path="tools/samtools/samtools",
  bin_path2="tools/picard/build/libs/picard.jar",bin_path3="tools/bedtools2/bin/bedtools",
  bam="",output_dir="",ref_genome="",verbose=FALSE,tmp_dir=".",mapq=0,bi="",
  ti="",ri="",ref_flat="",method="tg",mode="local",executor=make_unique_id("alignQC"),
  task="alignQC",time="48:0:0",threads=4,ram=4,update_time=60,wait=FALSE, hold=""){
    
    out_file_dir=set_dir(dir=output_dir,name="alignqc_report")
    
    ref=""
    if (!ref_genome==""){
      ref=paste0(" R=",ref_genome)
    }

  

    ## Generate alignment metrics

   mapq_metrics_bam_samtools(bin_path=bin_path,bam=bam,output_dir=out_file_dir,
   verbose=verbose,time=time,mode=mode,
   threads=threads,ram=ram,update_time=update_time,wait=FALSE, hold=hold)

   summary_metrics_bam_picard(bin_path=bin_path2,bam=bam,output_dir=out_file_dir,
    verbose=verbose,tmp_dir=tmp_dir,mode=mode,executor=executor,threads=threads,
    ram=ram,update_time=update_time,wait=FALSE, hold=hold)
    
   insertsize_metrics_bam_picard(bin_path=bin_path2,bam=bam,output_dir=out_file_dir,
   verbose=verbose,tmp_dir=tmp_dir,mode=mode,executor=executor,threads=threads,ram=ram,
   update_time=update_time,wait=FALSE,hold=hold)


    ## Only call metrics for panel data if bait and target intervals are supplied, Otherwise WGS metrics.
    
    if (method=="tg"){
   
     tg_summary_metrics_bam_picard(bin_path=bin_path2,bam=bam,output_dir=out_file_dir,
      verbose=verbose,tmp_dir=tmp_dir,bi=bi,ti=ti,mode=mode,executor=executor,threads=threads,ram=ram,
      update_time=update_time,wait=FALSE,hold=hold)
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
        rnaseq_summary_metrics_bam_picard(bin_path=bin_path2,
        bam=bam,output_dir=out_file_dir,verbose=verbose,tmp_dir=tmp_dir,
        ri=ri,ref_flat=ref_flat,mode=mode,executor=executor,
        threads=threads,ram=ram,update_time=update_time,wait=FALSE,hold=hold)
    }else if(method=="wgs"){
        wgs_summary_metrics_bam_picard(bin_path=bin_path2,
        bam=bam,output_dir=out_file_dir,verbose=verbose,tmp_dir=tmp_dir,mode=mode,
        executor=executor,threads=threads,ram=ram,
        update_time=update_time,wait=FALSE,hold=hold)
    }

  }

#' Generate a WIG file
#'
#' This function generates a WIG file.
#'
#'
#' @param bam Path to the BAM file .
#' @param bin_path Path to readCounter executable. Default path tools/samtools/samtools.
#' @param bin_path2 Path to readCounter executable. Default path tools/hmmcopy_utils/bin/readCounter.
#' @param output_dir Path to the output directory.
#' @param chrs String of chromosomes to include. c()
#' @param win Size of non overlaping windows. Default 500000.
#' @param format Output format [wig/seg] . Default wig
#' @param threads Number of threads to use. Default 3
#' @param verbose Enables progress messages. Default False.
#' @export


read_counter=function(bin_path="tools/samtools/samtools",bin_path2="tools/hmmcopy_utils/bin/readCounter",
chrs=c(1:22,"X","Y"),win=500000, format="wig", bam="",output_dir="",verbose=FALSE,threads=3){

    win=format(win,scientific=F)


    out_file_dir=set_dir(dir=output_dir,name=paste0("read_counter/",format))
    out_file=paste0(out_file_dir,"/",get_file_name(bam))

    ## Check chr notation
    exec_code=paste(bin_path,"view",bam," | head -n 1 | awk -F \"\t\" '{cat $3}'")
    
    if (verbose){
      print_verbose(exec_code=exec_code)
    }
    chr=system(exec_code,intern=TRUE)

    fmt=""

    if (format=="seg"){
      fmt="-s"
    }

     if (threads>1){
        parallel::mclapply(seq(1,length(chrs)),FUN=function(x){
          exec_code=cat(paste(bin_path2,fmt,"--window", win,"--quality 20 --chromosome",
            paste0("chr",chrs[x],collapse=","), bam,">", paste0(out_file,".",x,".",format)))
          if (verbose){
            print_verbose(exec_code=exec_code)
          }
          system(exec_code)
        },mc.cores=threads
      )
      system(paste0("ls -v -d ",out_file_dir,"/* | xargs cat >",out_file,".",format))
      system(paste0("rm ",out_file,".*.",format))
      }else{
        exec_code=paste(bin_path2,fmt,"--window", win,"--quality 20 --chromosome",
          paste0("chr",chrs,collapse=","), bam,">", paste0(out_file,".",format))
        if (verbose){
          print_verbose(exec_code=exec_code)
        }
        system(exec_code)
      }
    
    if (grepl("chr",chr)){
     
      exec_code=paste("sed -i 's/chrom=chr/chrom=/g'",paste0(out_file,".",format))
      if (verbose){
          print_verbose(exec_code=exec_code)
      }
      system(exec_code)
    }
  }

#' Filter BAM to specific regions
#'
#' This function takes a BAM file and filter it to genomic coordinate or otherwise BED file with
#' multiple genomic coordinates. The output is a BAM file with reads for inputed genomic region/s
#'
#'
#' @param bin_path Path to readCounter executable. Default path tools/samtools/samtools.
#' @param bam Path to the BAM file .
#' @param position String of genomic position to filter. Ex chr6:1000-100000
#' @param output_name Name of the output file.
#' @param output_dir Path to output directory.
#' @param bed Size of non overlaping windows. Default 500000.
#' @param threads Number of threads to use. Default 1
#' @param verbose Enables progress messages. Default False.
#' @export

filter_bam=function(bin_path="tools/samtools/samtools",bam="",position="",bed="",
verbose=FALSE,output_dir="",threads=1){


  out_file_dir=set_dir(dir=output_dir,name="filtered")
  if (position!="" &bed!=""){
    cat("Position and bed arguments are mutually exclusive")
    quit()
  }

  if(bed!=""){
    bed=paste("-L",bed)
  }

  exec_code=paste(bin_path,"view -b",bed,"-@",threads,bam,position,">", paste0(out_file_dir,"/",get_file_name(bam),".filtered.bam"))

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
#' @param bin_path Path to ichorCNA executable. Default path tools/ichorCNA/scripts/runIchorCNA.R.
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

ichorCNA=function(bin_path="tools/ichorCNA/scripts/runIchorCNA.R",sample_id="",
wig="",norm_wig="",bed="",ploidy="2,3",tumour_content="0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9",
homozygous_del="False",subclonal_states=NULL,
gc="tools/ichorCNA/inst/extdata/gc_hg19_500kb.wig",
map="tools/ichorCNA/inst/extdata/map_hg19_500kb.wig",
centromere="tools/ichorCNA/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt",
normal_panel="tools/ichorCNA/inst/extdata/HD_ULP_PoN_500kb_median_normAutosome_mapScoreFiltered_median.rds",
output_dir="",verbose=TRUE,libdir="tools/ichorCNA",male_tresh=0.0001,chrs="'c(1:22,\"X\")'",chrTrain="'c(1:22)'"){


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

    exec_code=paste("Rscript",bin_path,"--id",sample_id,"--WIG",wig,norm_wig,bed,
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
#' @param bin_path Path to ichorCNA executable. Default path tools/ichorCNA/scripts/createPanelOfNormals.
#' @param output_dir Path to the output directory.
#' @param bed Path to BED file with target regions.
#' @param output_name File output name.
#' @param gc Path to GC-content WIG with . Default tools/ichorCNA/inst/extdata/gc_hg19_500kb.wig
#' @param map Path to mappability score WIG with GC content. Default tools/ichorCNA/inst/extdata/map_hg19_500kb.wig
#' @param centromere Path to file containing centromere locations. Default tools/ichorCNA/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt
#' @param male_tresh Minimum percentage of reads in chromosome Y to call male. Default 0.0001
#' @param verbose Enables progress messages. Default False.
#' @export

panel_of_normals_ichorCNA=function(bin_path="tools/ichorCNA/scripts/createPanelOfNormals.R",wigs_dir="",wigs="",bed="",
output_name="PoN_ichorCNA",gc="tools/ichorCNA/inst/extdata/gc_hg19_500kb.wig",
map="tools/ichorCNA/inst/extdata/map_hg19_500kb.wig",
centromere="tools/ichorCNA/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt",male_tresh=0.0001,verbose=TRUE){

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

    exec_code=paste("Rscript",bin_path,"--filelist",wig_list,bed," --gcWig", gc,
      "--mapWig",map,"--centromere",centromere,"--fracReadsInChrYForMale",male_tresh," --outfile",output_name)
    if(verbose){
      print_verbose(exec_code=exec_code)
    }
    system(exec_code)
    system(paste("rm",wig_list))

}
