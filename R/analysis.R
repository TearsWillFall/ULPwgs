#' Generate a quality control (QC) report from a fastaqc file
#'
#' This function takes a set of sequence files (fastq,SAM,BAM...) and
#' returns a report in HTML format.
#'
#'
#' @param file_R1 Path to the input file with the sequence.
#' @param file_R2 [Optional] Path to the input with the reverse read sequence.
#' @param bin_path Path to fastQC executable. Default path tools/FastQC/bin/fastqc.
#' @param n_cores Number of CPU cores to use. Default 3.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @export



fastqc=function (bin_path="tools/FastQC/bin/fastqc",file_R1="",file_R2="",n_cores=3,output_dir="",verbose=FALSE){

  sep="/"

  if(output_dir==""){
    sep=""
  }

  sample_name=get_sample_name(file_R1)

  if (!file_R2==""){
  sample_name=intersect_sample_name(file_path=file_R1,file_path2=file_R2)
  output_dir=paste0(output_dir,sep,sample_name,"_FastQC_reports")

  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  if(verbose){
    print(paste(bin_path,"-o ",output_dir,"-t ",n_cores,"--noextract",file_R1,file_R2))
  }
  system(paste(bin_path,"-o ",output_dir,"-t ",n_cores,"--noextract",file_R1,file_R2))

}else{
  output_dir=paste0(output_dir,sep,sample_name)

  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }
    if(verbose){
      print(paste(bin_path,"-o ",output_dir,"-t ",n_cores,"--noextract",file_R1))
    }
    system(paste(bin_path,"-o ",output_dir,"-t ",n_cores,"--noextract",file_R1))
  }
}



#' Adapter sequence trimmer
#'
#' This function takes a single/multiple fastq (if paired-end reads) files and
#' trims the adapter sequences found within them.
#'
#' @param file_R1 Path to the input file with the sequence.
#' @param file_R2 [Optional] Path to the input with the reverse read sequence.
#' @param xadapt [Optional] Adapter sequence/file.
#' @param yadapt [Optional] Adapter sequence/file.
#' @param bin_path Path to skewer executable. Default path tools/skewer/skewer.
#' @param n_cores Number of CPU cores to use. Default 3.
#' @param mean_quality Minimum mean quality of reads to be kept.Dedault 0
#' @param min_length Minimum length of reads to keep.Default 35
#' @param max_length Maximum length of reads to keep.Default NA.
#' @param verbose Enables progress messages. Default False.
#' @export


trimming=function(bin_path="tools/skewer/skewer",file_R1="",file_R2="",xadapt=NA,yadapt=NA,n_cores=3,output_dir="",verbose=FALSE,mean_quality=0,min_length=35,max_length=NA){

  sep="/"

  if(output_dir==""){
    sep=""
  }

  sample_name=get_sample_name(file_R1)
  if ((!is.na(xadapt)) & (!is.na(yadapt))){
    func=paste(bin_path,"-m tail -t",n_cores,"-x", xadapt,"-y", yadapt,"-Q",mean_quality,"-l",min_length)
  }
  else{
    func=paste(bin_path,"-m tail -t",n_cores,"-Q",mean_quality,"-l",min_length)
  }
  if(!is.na(max_length)){
    func=paste(func,"-L",max_length)
  }

  if (!file_R2==""){
  sample_name=intersect_sample_name(file_path=file_R1,file_path2=file_R2)

    if(verbose){
      print(paste(func,"-z -f sanger --quiet -o",sample_name,file_R1,file_R2))
    }
    system(paste(func,"-z -f sanger --quiet -o",sample_name,file_R1,file_R2))


  }else{


      if(verbose){
        print(paste(func,"-z -f sanger --quiet -o",sample_name,file_R1))
      }
      system(paste(func,"-z -f sanger --quiet -o",sample_name,file_R1))
    }
    output_dir=paste0(output_dir,sep,sample_name,"_trimmed")
    if(!dir.exists(output_dir)){
      dir.create(output_dir)
    }

    if(verbose){
      print(paste("mv",paste0(sample_name,"-trimmed*"),output_dir))
    }
    system(paste("mv",paste0(sample_name,"-trimmed*"),output_dir))


  }


#' Merge BAM files in directory
#'
#' This function takes a BAM file and merges it with others found within a
#' directory. This function is still Work In Progress (WIP).
#'
#' @param bam Path to the input bam file with the sequence.
#' @param bam_dir Path to directory with BAM files to merge.
#' @param verbose Enables progress messages. Default False.
#' @export

merge_bam=function(bin_path="tools/samtools/samtools",bam="",bam_dir="",verbose=FALSE){
    if(verbose){
      print(paste(bin_path,"merge",bam, paste0(bam_dir,"/*.bam")))
    }
    system(paste(bin_path,"merge",bam, paste0(bam_dir,"/*.bam")))
  }


#' Read alignment
#'
#' This function aligns a sequence of reads to a reference genome
#' using a Burrows-Wheeler aligner. The reference genome has to be previously
#' been sorted and indexed. It generates a BAM file with the aligned sequence.
#'
#' @param file_R1 Path to the input file with the sequence.
#' @param file_R2 [Optional] Path to the input with the reverse read sequence.
#' @param bin_path Path to bwa executable. Default path tools/bwa/bwa.
#' @param bin_path2 Path to samtools executable. Default path tools/samtools/samtools.
#' @param ref_genome Path to input file with the reference genome sequence.
#' @param n_cores Number of CPU cores to use. Default 3.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @export

alignment=function(bin_path="tools/bwa/bwa",bin_path2="tools/samtools/samtools",file_R1="",file_R2="",n_cores=3,ref_genome="",output_dir="",verbose=FALSE){

    sep="/"

    if(output_dir==""){
      sep=""
    }

    sample_name=get_sample_name(file_R1)
    GPU=paste0("\"@RG\\tID:",sample_name,"\\tPL:ILLUMINA\\tPU:NA\\tLB:",sample_name,"\\tSM:",sample_name,"\"")

    if (!file_R2==""){
    sample_name=intersect_sample_name(file_path=file_R1,file_path2=file_R2)
    output_dir=paste0(output_dir,sep,sample_name,"_BAM")
    if(!dir.exists(output_dir)){
      dir.create(output_dir)
    }

    GPU=paste0("\"@RG\\tID:",sample_name,"\\tPL:ILLUMINA\\tPU:NA\\tLB:",sample_name,"\\tSM:",sample_name,"\"")
    out_file=paste0(output_dir,"/",sample_name,".bam")
      if(verbose){
          print(paste(bin_path,"mem -t", n_cores," -v 2 -R",GPU,"-M",ref_genome, file_R1,file_R2, "|",bin_path2," view -h -b >",out_file))
      }
      system(paste(bin_path,"mem -t", n_cores," -v 2 -R",GPU,"-M",ref_genome, file_R1,file_R2, "| ",bin_path2," view -h -b >",out_file))

      }
    else{
      output_dir=paste0(output_dir,sep,sample_name,"_BAM")
      out_file=paste0(output_dir,"/",sample_name,".bam")
      if(!dir.exists(output_dir)){
        dir.create(output_dir)
      }

      if(verbose){
          print(paste(bin_path,"mem -t", n_cores," -v 2 -R",GPU,"-M",ref_genome, file_R1, "| ",paste0("./",bin_path2)," view -h -b >",out_file))
      }
      system(paste(bin_path,"mem -t", n_cores," -v 2 -R",GPU,"-M",ref_genome, file_R1, "| ",paste0("./",bin_path2)," view -h -b >",out_file))

    }
  }



#' Sort and index a sequence file
#'
#' This function sorts and indexes genomic sequence files.
#'
#' @param file Path to the input file with the sequence.
#' @param bin_path Path to bwa executable. Default path tools/samtools/samtools.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param threads Number of threads. Default 3
#' @param ram Ram memory to use per thread in GB. Default 1GB
#' @export


sort_and_index=function(bin_path="tools/samtools/samtools",file="",output_dir="",ram=1,verbose=FALSE,threads=3){

  sep="/"

  if(output_dir==""){
    sep=""
  }

  sample_name=get_sample_name(file)
  file_ext=get_file_extension(file)
  out_file=paste0(output_dir,sep,sample_name,"_SORTED.",toupper(file_ext))


  if (!dir.exists(out_file)){
      dir.create(out_file)
    }
  out_file=paste0(out_file,"/",sample_name)


    if (verbose){
      print("Sorting BAM file:")
      print(paste0(bin_path," sort ",file," -@ ",threads," -m ",ram,"G"," -o ",out_file,".SORTED.",file_ext))
    }
    system(paste0(bin_path," sort ",file," -@ ",threads," -m ",ram,"G"," -o ",out_file,".SORTED.",file_ext))
    file=paste0(out_file,".SORTED.",file_ext)

    if (verbose){
      print("Indexing sorted BAM file:")
      print(paste(bin_path," index",file," -@ ",threads))
    }
    system(paste(bin_path," index",file," -@ ",threads))
    if (verbose){
      print("Generating Flag stats:")
      print(paste0(bin_path," flagstat ",file," -@ ",threads," > ",paste0(out_file,".flagstat.txt")))
    }
    system(paste0(bin_path," flagstat ",file," -@ ",threads," > ",paste0(out_file,".flagstat.txt")))
    if (verbose){
      print("Generating Index stats:")
      print(paste0(bin_path," idxstats ",file," > ",paste0(out_file,".idxstats.txt")))
    }
    system(paste0(bin_path," idxstats ",file," > ",paste0(out_file,".idxstats.txt")))
  }




#' Remove duplicated reads
#'
#' This function removes duplicated reads (artifacts) found in aligned sequences.
#'
#' @param file Path to the input file with the aligned sequence.
#' @param bin_path Path to picard executable. Default path tools/picard/build/libs/picard.jar.
#' @param output_dir Path to the output directory.
#' @param tmp_dir Path to tmp directory.
#' @param verbose Enables progress messages. Default False.
#' @param hnd Maximum number of file handles. Default 1000.
#' @param ram RAM in GB to use. Default 4 Gb.
#' @export


remove_duplicates=function(bin_path="tools/picard/build/libs/picard.jar",file="",output_dir="",verbose=FALSE,hnd=1000,ram=4,tmp_dir=""){

    sep="/"

    if(output_dir==""){
      sep=""
    }

    sample_name=get_sample_name(file)
    file_ext=get_file_extension(file)

    out_file=paste0(output_dir,sep,sample_name,"_RMDUP.",toupper(file_ext))
    if (!dir.exists(out_file)){
        dir.create(out_file)
    }

    out_file=paste0(out_file,"/",sample_name)

    tmp=""
    if (!tmp_dir==""){
      tmp=paste0("TMP_DIR=",tmp_dir)
    }


      if(verbose){
        print(paste0("java -Xmx",ram,"g", " -Djava.io.tmpdir=",tmp_dir," -jar ",bin_path," MarkDuplicates I=",file, " O=",paste0(out_file,".RMDUP.",file_ext)," M=",paste0(out_file,".picard_rmdup.txt")," REMOVE_DUPLICATES=true AS=true VALIDATION_STRINGENCY=LENIENT ",paste0("MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=",hnd)," ",tmp))

      }
      system(paste0("java -Xmx",ram,"g", " -Djava.io.tmpdir=",tmp_dir," -jar ",bin_path," MarkDuplicates I=",file, " O=",paste0(out_file,".RMDUP.",file_ext)," M=",paste0(out_file,".picard_rmdup.txt")," REMOVE_DUPLICATES=true AS=true VALIDATION_STRINGENCY=LENIENT " ,paste0("MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=",hnd)," ",tmp))

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
#' @param off_tar Path to off target regions bed file. Has to be chromosome sorted. Requires bi and ti arguments.
#' @param on_tar Path to target region bed file. Has to be chromosome sorted. Requires bi and ti arguments.
#' @param bi Bait capture target interval for panel data. Requires ti and off_target and on_tar arguments. Interval format.
#' @param ti Primary target intervals for panel data. Requires bi and off_tar and on_tar argmunets. Interval format.
#' @export


qc_metrics=function(bin_path="tools/samtools/samtools",bin_path2="tools/picard/build/libs/picard.jar",bin_path3="tools/bedtools2/bin/bedtools",bam="",output_dir="",ref_genome="",verbose=FALSE,ram="",tmp_dir="",mapq=0,bi="",ti="",off_tar="",on_tar=""){
    sep="/"

    if(output_dir==""){
      sep=""
    }

    sample_name=get_sample_name(bam)

    out_file_dir=paste0(output_dir,sep,sample_name,"_alignQC_report")
    if (!dir.exists(out_file)){
        dir.create(out_file)
    }

    out_file=paste0(out_file_dir,"/",sample_name)

    ref=""
    if (!ref_genome==""){
      ref=paste0(" R=",ref_genome)
    }
    tmp=""
    if (!tmp_dir==""){
      tmp=paste0(" TMP_DIR=",tmp_dir)
    }

    ## Generate alignment metrics

    if (verbose){
      print("Generate MapQ distance map:")
      print(paste(bin_path,"view",bam," | awk -F", "'\\t'", "'{c[$5]++} END { for (i in c) printf(\"%s\\t%s\\n\",i,c[i]) }'"," | sort -t$'\\t' -k 1 -g >>", paste0(out_file,".mapq_dist.txt")))

    }
    system(paste(bin_path,"view",bam," | awk -F", "'\\t'", "'{c[$5]++} END { for (i in c) printf(\"%s\\t%s\\n\",i,c[i]) }'"," | sort -t$'\\t' -k 1 -g >>", paste0(out_file,".mapq_dist.txt")))

    if (verbose){
      print("Generate Alignment Metrics:")
      print(paste0("java -Xmx",ram,"g", " -Djava.io.tmpdir=",tmp_dir," -jar ",bin_path2," CollectAlignmentSummaryMetrics "," I=",bam," O=",paste0(out_file,".picard_summary.txt "),tmp))

    }
    system(paste0("java -Xmx",ram,"g", " -Djava.io.tmpdir=",tmp_dir," -jar ",bin_path2," CollectAlignmentSummaryMetrics "," I=",bam," O=",paste0(out_file,".picard_summary.txt "),tmp))

    if (verbose){
      print("Generate Insert Size Metrics:")
      print(paste0("java -Xmx",ram,"g", " -Djava.io.tmpdir=",tmp_dir," -jar ",bin_path2," CollectInsertSizeMetrics "," I=",bam," O=",paste0(out_file,".picard_insert_size.txt")," H=",paste0(out_file,".picard_insert_size.pdf "),tmp))

    }
    system(paste0("java -Xmx",ram,"g", " -Djava.io.tmpdir=",tmp_dir," -jar ",bin_path2, " CollectInsertSizeMetrics "," I=",bam," O=",paste0(out_file,".picard_insert_size.txt")," H=",paste0(out_file,".picard_insert_size.pdf "),tmp))


    ## Only call metrics for panel data if bait and target intervals are supplied, Otherwise WGS metrics.
    if (bi!="" & ti!=""){
      if (verbose){
        print("Generate Panel Metrics:")
        print(paste0("java -Xmx",ram,"g", " -Djava.io.tmpdir=",tmp_dir," -jar ",bin_path2," CollectHsMetrics BI=",bi," TI=",ti," I=",bam," THEORETICAL_SENSITIVITY_OUTPUT=",paste0(out_file,".TS.txt"),ref," O=",paste0(out_file,".CollectHSmetrics.txt "),tmp))

      }
      system(paste0("java -Xmx",ram,"g", " -Djava.io.tmpdir=",tmp_dir," -jar ",bin_path2," CollectHsMetrics BI=",bi," TI=",ti," I=",bam," THEORETICAL_SENSITIVITY_OUTPUT=",paste0(out_file,".TS.txt"),ref," O=",paste0(out_file,".CollectHSmetrics.txt "),tmp))


      ## Picard doesn't output coverage stats for off-target regions therefore we have to estimate this manually.

      ## I use this function to get the mean coverage on and off target.

      ## For target regions

      bed_coverage(bin_path=bin_path3,bam=bam,bed=on_tar,verbose=verbose,sorted=TRUE,mean=TRUE,hist=TRUE,fai=paste0(ref_genome,".fai"),suffix="on_Target",output_dir=out_file_dir)
      bed_coverage(bin_path=bin_path3,bam=bam,bed=on_tar,verbose=verbose,sorted=TRUE,mean=TRUE,hist=FALSE,fai=paste0(ref_genome,".fai"),suffix="on_Target",output_dir=out_file_dir)

      ## For off target regions

      bed_coverage(bin_path=bin_path3,bam=bam,bed=off_tar,verbose=verbose,sorted=TRUE,mean=TRUE,hist=TRUE,fai=paste0(ref_genome,".fai"),suffix="off_Target",output_dir=out_file_dir)
      bed_coverage(bin_path=bin_path3,bam=bam,bed=off_tar,verbose=verbose,sorted=TRUE,mean=TRUE,hist=FALSE,fai=paste0(ref_genome,".fai"),suffix="off_Target",output_dir=out_file_dir)

      ## Generate violin and cummulative plots for target and off target regions

      plot_coverage_panel(on_target=paste0(out_file,".on_Target.Per_Region_Coverage.txt"),off_target=paste0(out_file,".off_Target.Per_Region_Coverage.txt"),col=c(5,4),height=6,width=12,output_dir=out_file_dir)
      plot_cumulative_cov(on_target=paste0(out_file,".on_Target.Histogram_Coverage.txt"),off_target=paste0(out_file,".off_Target.Histogram_Coverage.txt"),height=6,width=12,output_dir=out_file_dir)

    }else{
          if (verbose){
            print(paste0("Generate WGS Metrics for minimum MAPq=",mapq,":"))
            print(paste0("java -Xmx",ram,"g", " -Djava.io.tmpdir=",tmp_dir," -jar ",bin_path2," CollectWgsMetrics MINIMUM_MAPPING_QUALITY=",mapq," ",ref," I=",bam," O=",paste0(out_file,".picard_wgs_q00.txt "),tmp))

          }
          system(paste0("java -Xmx",ram,"g", " -Djava.io.tmpdir=",tmp_dir," -jar ",bin_path2," CollectWgsMetrics MINIMUM_MAPPING_QUALITY=",mapq," ",ref," I=",bam," O=",paste0(out_file,".picard_wgs_q00.txt "),tmp))

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
#' @param chr String of chromosomes to include. c()
#' @param win Size of non overlaping windows. Default 500000.
#' @param verbose Enables progress messages. Default False.
#' @export


read_counter=function(bin_path="tools/samtools/samtools",bin_path2="tools/hmmcopy_utils/bin/readCounter",chr=c(1:22,"X","Y"),win=500000,bam="",output_dir="",verbose=FALSE){


    win=format(win,scientific=F)
    sep="/"

    if(output_dir==""){
      sep=""
    }

    sample_name=get_sample_name(bam)

    out_file=paste0(output_dir,sep,sample_name,"_WIG")
    if (!dir.exists(out_file)){
        dir.create(out_file)
    }

    out_file=paste0(out_file,"/",sample_name)

    if (verbose){
      print(paste(bin_path,"view",bam," | head -n 1 | awk -F \"\t\" '{print $3}'"))
    }
    chr=system(paste(bin_path,"view",bam," | head -n 1 | awk -F \"\t\" '{print $3}'"),intern=TRUE)

    if (grepl("chr",chr)){
      if (verbose){
        print(paste(bin_path2,"--window", win,"--quality 20 --chromosome",paste0("chr",chr,collapse=","), bam,">", paste0(out_file,".wig")))
      }
      system(paste(bin_path2,"--window", win,"--quality 20 --chromosome",paste0("chr",chr,collapse=","), bam,">" ,paste0(out_file,".wig")))

        if (verbose){
          print(paste("sed -i 's/chrom=chr/chrom=/g'",paste0(out_file,".wig")))
      }
      system(paste("sed -i 's/chrom=chr/chrom=/g'",paste0(out_file,".wig")))
    }
    else{
      if (verbose){
        print(paste(bin_path2,"--window", win,"--quality 20 --chromosome",paste0(chr,collapse=","), bam,">", paste0(out_file,".wig")))
      }
      system(paste(bin_path2,"--window", win,"--quality 20 --chromosome",paste0(chr,collapse=","), bam,">" ,paste0(out_file,".wig")))


    }
  }


#' Generate report for ULP-WGS samples
#'
#' This function generates a report that helps to segment the genome, predict large-scale
#' copy number alterations, and estimate tumor fraction in ULP-WGS samples.
#'
#' @param wig Path to the WIG file.
#' @param bin_path Path to ichorCNA executable. Default path tools/ichorCNA/scripts/runIchorCNA.R.
#' @param output_dir Path to the output directory.
#' @param sample_id String with sample name.
#' @param ploidy Initial tumour ploidy. Default 2,3
#' @param tumour_content Initial normal contamination. Default 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9
#' @param homozygous_del Include Homozygous deleteions. Default FALSE.
#' @param gc Path to GC-content WIG with . Default tools/ichorCNA/inst/extdata/gc_hg19_500kb.wig
#' @param map Path to mappability score WIG with GC content. Default tools/ichorCNA/inst/extdata/map_hg19_500kb.wig
#' @param centromere Path to file containing centromere locations. Default tools/ichorCNA/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt
#' @param normal_panel Path to file containing normal panel data. Default tools/ichorCNA/inst/extdata/HD_ULP_PoN_hg38_500kb_median_normAutosome_median.rds
#' @param libdir Path to dir containing ichorCNA libraries.
#' @param chrs Chromosomes to analyze. Default 'c(1:22,\"X\")'
#' @param chrsTrain Chromosomes to estimate parameters. Default 'c(1:22)'
#' @param verbose Enables progress messages. Default False.
#' @export

ichorCNA=function(bin_path="tools/ichorCNA/scripts/runIchorCNA.R",sample_id="",wig="",ploidy="2,3",tumour_content="0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9",homozygous_del="False",subclonal_states="NULL",gc="tools/ichorCNA/inst/extdata/gc_hg19_500kb.wig",map="tools/ichorCNA/inst/extdata/map_hg19_500kb.wig",centromere="tools/ichorCNA/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt",normal_panel="tools/ichorCNA/inst/extdata/HD_ULP_PoN_hg38_500kb_median_normAutosome_median.rds",output_dir="",verbose=TRUE,libdir="tools/ichorCNA",chrs="'c(1:22,\"X\")'",chrTrain="'c(1:22)'"){

    sep="/"

    if(output_dir==""){
      sep=""
    }
    sample_name=get_sample_name(wig)
    if(!sample_id==""){
      sample_name=sample_id
    }else{
      sample_id=sample_name
    }

    out_file=paste0(output_dir,sep,sample_name,"_CNAreport")
    if (!dir.exists(out_file)){
        dir.create(out_file)
    }
    if(verbose){
      print(paste("Rscript",bin_path,"--id",sample_id,"--WIG",wig,"--ploidy",paste0("'c(",ploidy,")'"),"--normal",paste0("'c(",tumour_content,")'"),"--maxCN 7 --gcWig", gc,"--mapWig",map,"--centromere",centromere,"--normalPanel",normal_panel,"--includeHOMD",homozygous_del,"--chrs",chrs,"--chrTrain",chrTrain,"--estimateNormal True --estimatePloidy True --estimateScPrevalence True --outDir",out_file,"--libdir",libdir))

    }
    system(paste("Rscript",bin_path,"--id",sample_id,"--WIG",wig,"--ploidy",paste0("'c(",ploidy,")'"),"--normal",paste0("'c(",tumour_content,")'"),"--maxCN 7 --gcWig", gc,"--mapWig",map,"--centromere",centromere,"--normalPanel",normal_panel,"--includeHOMD",homozygous_del,"--chrs",chrs,"--chrTrain",chrTrain,"--estimateNormal True --estimatePloidy True --estimateScPrevalence True --outDir",out_file,"--libdir",libdir))
  }
