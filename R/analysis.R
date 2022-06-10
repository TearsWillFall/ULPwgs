#' Generate a quality control (QC) report from a fastaqc file
#'
#' This function takes a set of sequence files (fastq,SAM,BAM...) and
#' returns a report in HTML format.
#'
#'
#' @param file_R1 Path to the input file with the sequence.
#' @param file_R2 [Optional] Path to the input with the reverse read sequence.
#' @param bin_path Path to fastQC executable. Default path tools/FastQC/bin/fastqc.
#' @param threads Number of CPU cores to use. Default 3.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @export


qc_fastqc=function (bin_path="tools/FastQC/bin/fastqc",file_R1="",file_R2="",threads=3,
output_dir="",verbose=FALSE){

  out_file_dir=set_dir(dir=output_dir,name="fastqc_reports")


  if (!file_R2==""){

    exec_code=paste(bin_path,"-o ", out_file_dir,"-t ",threads,"--noextract",file_R1,file_R2)


  }else{

    exec_code=paste(bin_path,"-o ", out_file_dir,"-t ",threads,"--noextract",file_R1)
  }

  if(verbose){
      print(exec_code)
  }

    system(exec_code)
  }




#' Adapter sequence trimmer using skewer
#'
#' Detects and removes adapter sequences from single and paired read data
#'
#' @param file_R1 Path to the input file with the sequence.
#' @param file_R2 [Optional] Path to the input with the reverse read sequence.
#' @param xadapt [Optional] Adapter sequence/file.
#' @param yadapt [Optional] Adapter sequence/file.
#' @param bin_path Path to skewer executable. Default path tools/skewer/skewer.
#' @param threads Number of CPU cores to use. Default 3.
#' @param mean_quality Minimum mean quality of reads to be kept.Dedault 0
#' @param min_length Minimum length of reads to keep.Default 35
#' @param max_length Maximum length of reads to keep.Default NA.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @export


trimming_skewer=function(bin_path="tools/skewer/skewer",file_R1="",file_R2="",xadapt=NA,
yadapt=NA,threads=3,output_dir="",verbose=FALSE,mean_quality=0,min_length=35,max_length=NA){

  out_file_dir=set_dir(dir=output_dir,name="skewer_reports")

  if ((!is.na(xadapt)) & (!is.na(yadapt))){
    func=paste(bin_path,"-m tail -t",threads,"-x", xadapt,"-y", yadapt,"-Q",mean_quality,"-l",min_length)
  }

  else{
    func=paste(bin_path,"-m tail -t",threads,"-Q",mean_quality,"-l",min_length)
  }
  if(!is.na(max_length)){
    func=paste(func,"-L",max_length)
  }

  if (!file_R2==""){
    exec_code=paste(func,"-z -f sanger --quiet -o",paste0(out_file_dir,"/",intersect_file_name(file_R1,file_R2)),file_R1,file_R2)
  }else{
    exec_code=paste(func,"-z -f sanger --quiet -o",paste0(out_file_dir,"/",get_file_name(file_R1)),file_R1)
  }

  if(verbose){
      print(exec_code)
    }
    system(exec_code)
}

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


merge_bams=function(bin_path="tools/samtools/samtools",bams="",output_name="",
  verbose=FALSE,threads=3){
    exec_code=paste(bin_path,"merge ",paste0(output_name,".bam"), " --threads",
      threads,paste(bams,collapse=" "))

    if(verbose){
      print(exec_code)
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
      print(exec_code)
    }
    system(exec_code)
}



#' Read alignment using BWA
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
#' @param id_tag Read group identifier. Default NA.
#' @param pu_tag Platform unit identifier. Default NA.
#' @param sm_tag Sample name tag. If not given will use sample name as default.
#' @param pl_tag Platform/technology used to produce the read. Default "ILLUMINA". Options ["ILLUMINA","SOLID", "LS454", "HELICOS","PACBIO"]
#' @param lb_tag DNA preparation library identifier. 
#' @param sort Sort aligned file. Default TRUE.
#' @param coord_sort Sort BAM file by coordinate. Alternatively sort by name. Default TRUE.
#' @param index Generate index file if BAM sorted by coordinate. Default TRUE.
#' @param ram RAM memory to use for sorting and indexing. Provided in GB.
#' @param threads Number of CPU cores to use. Default 3.
#' @param stats Generate BAM stats. Default all. Options ["all","flag","index",""
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @export

alignment_bwa=function(bin_path="tools/bwa/bwa",bin_path2="tools/samtools/samtools",
file_R1="",file_R2="",threads=3,ram=1,id_tag="NA",pu_tag="NA",pl_tag="ILLUMINA",lb_tag="NA",
sm_tag="",sort=TRUE,coord_sort=TRUE,index=TRUE,stats="all",ref_genome="",output_dir="",
verbose=FALSE){
    
    out_file_dir=set_dir(dir=output_dir,name="bwa_reports")
   
    input_files=file_R1
    if (!file_R2==""){
      sm_tag=ifelse(sm_tag=="",get_file_name(file_R1),sm_tag)
    }else{
      sm_tag=ifelse(sm_tag=="",intersect_file_name(file_R1,file_R2),sm_tag)
      input_files=paste(file_R1,file_R2)
    }
    out_file=paste0(out_file_dir,"/",sm_tag,".bam")
    GPU=paste0("\"@RG\\tID:",id_tag,"\\tPL:",pl_tag,"\\tPU:",pu_tag,"\\tLB:",
    lb_tag,"\\tSM:",sm_tag,"\"")

    exec_code=paste(bin_path,"mem -t", threads," -v 2 -R",GPU,"-M",ref_genome,
        input_files, "| ",paste0(bin_path2)," view -h -b >",out_file)
    if(verbose){
        print(exec_code)
    }
    system(exec_code)

    if(sort){
      sort_and_index_samtools(bin_path=bin_path2,bam=out_file,output_dir=out_file_dir,
      ram=ram,verbose=verbose,threads=threads,coord_sort=coord_sort,index=index)
      system(paste0("rm ",out_file))
    }

  }



#' Sort and index a sequence file
#'
#'
#' @param file Path to the input file with the sequence.
#' @param bin_path Path to bwa executable. Default path tools/samtools/samtools.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param threads Number of threads. Default 3
#' @param coord_sort Generate a coord sorted file. Otherwise queryname sorted. Default TRUE
#' @param ram Ram memory to use per thread in GB. Default 1GB
#' @param index Generate an index file for sorted BAM. Default TRUE
#' @param stats Generate BAM stats. Default all. Options ["all","flag","index",""]
#' @export

sort_and_index_samtools=function(bin_path="tools/samtools/samtools",bam="",output_dir="",
ram=1,verbose=FALSE,threads=3,coord_sort=TRUE,index=TRUE,stats="all"){

  out_file_dir=set_dir(dir=output_dir,name="sorted")

  bam_sort_samtools(bin_path=bin_path,bam=bam,output_dir=out_file_dir,ram=ram,
  verbose=verbose,threads=threads,coord_sort=coord_sort)

  bam=paste0(out_file_dir,"/",get_file_nanme(bam),".sorted.",get_file_ext(bam)(bam))

  if (coord_sort){

    if(index){
        index_samtools(bin_path=bin_path,bam=bam,verbose=verbose,threads=threads)
      if(stats){
        bam_stats_samtools(bin_path=bin_path,bam=bam,output_dir=out_file_dir,
        verbose=verbose,threads=threads,stats=stats)
      }
     
    }else{
        bam_stats_samtools(bin_path=bin_path,bam=bam,output_dir=out_file_dir,
        verbose=verbose,threads=threads,stats="flag")
    }
  }
}

#' Mark duplicated reads
#'
#' This function marks duplicated reads (artifacts) found in aligned sequences.
#'
#' @param file Path to the input file with the aligned sequence.
#' @param bin_path Path to picard executable. Default path tools/picard/build/libs/picard.jar.
#' @param output_dir Path to the output directory.
#' @param tmp_dir Path to tmp directory.
#' @param verbose Enables progress messages. Default False.
#' @param remove_duplicates Do not write duplicates to the output file. Default FALSE
#' @param hnd Maximum number of file handles. Default 1000.
#' @param ram RAM in GB to use. Default 4 Gb.
#' @export


mark_duplicates_picard=function(bin_path="tools/picard/build/libs/picard.jar",bam="",
output_dir="",verbose=FALSE,hnd=1000,ram=4,tmp_dir="",remove_duplicates=TRUE){


    out_file_dir=set_dir(dir=output_dir,name="markdups_reports")

    tmp=""
    if (!tmp_dir==""){
      tmp=paste0("TMP_DIR=",tmp_dir)
    }

    if(remove_duplicates){
      remove_duplicates=" REMOVE_DUPLICATES=true "
    }else{
      remove_duplicates=" REMOVE_DUPLICATES=false "
    }

    exec_code=paste0("java -Xmx",ram,"g", " -Djava.io.tmpdir=",tmp_dir," -jar ",
      bin_path," MarkDuplicates I=",file, " O=",paste0(out_file_dir,"/",
      get_file_name(bam),".rmdup.",get_file_ext(bam)),
      " M=",paste0(out_file_dir,"/",
      get_file_name(bam),".picard_rmdup.txt"),
      remove_duplicates, " AS=true VALIDATION_STRINGENCY=LENIENT ",
      paste0("MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=",hnd)," ",tmp)
    if(verbose){
      print(exec_code)

    }
    system(exec_code)

  }

#' Wrapper for MarkDuplicatesSpark from gatk
#'
#' This function removes duplicated reads (artifacts) found in aligned sequences and sorts the output bam.
#'
#' @param file Path to the input file with the aligned sequence.
#' @param bin_path Path to picard executable. Default path tools/picard/build/libs/picard.jar.
#' @param output_dir Path to the output directory.
#' @param tmp_dir Path to tmp directory.
#' @param threads Number of threads to split work.
#' @param verbose Enables progress messages. Default False.
#' @param remove_duplicates Remove all sequencing duplicates from BAM file. Default TRUE.
#' @export


mark_duplicates_gatk=function(bin_path="tools/gatk/gatk",bam="",output_dir="",
verbose=FALSE,tmp_dir="",threads=3,remove_duplicates=TRUE){

      out_file_dir=set_dir(dir=output_dir,name="markdups_reports")

      tmp=""
      if (!tmp_dir==""){
        tmp=paste0("--tmp-dir ",tmp_dir)
      }

      dups=""
      if(remove_duplicates){
          dups="--remove-all-duplicates"
      }
    
      if(verbose){
        print(paste0(bin_path," MarkDuplicatesSpark -I ",file, " -O ",
        paste0(out_file,".sorted.rmdup.",get_file_ext(bam))," -M ",paste0(out_file_dir,"/",get_file_name(bam),".gatk_rmdup.txt"),
        " ",tmp," --conf \'spark.executor.cores=",threads,"\'"), dups)
      }
        system(paste0(bin_path," MarkDuplicatesSpark -I ",file, " -O ",
        paste0(out_file,".sorted.rmdup.",get_file_ext(bam))," -M ",paste0(out_file_dir,"/",get_file_name(bam),".gatk_rmdup.txt"),
        " ",tmp," --conf \'spark.executor.cores=",threads,"\'"), dups)
    }



#' Function for base quality recalibration
#'
#' This function recalibrates the base quality of the reads in two steps process based on GATK best practices guides.
#'
#' @param bin_path [REQUIRED] Path to picard executable. Default path tools/samtools/samtools.
#' @param bin_path2 [REQUIRED] Path to picard executable. Default path tools/gatk/gatk.
#' @param bin_path3 [REQUIRED] Path to picard executable. Default path tools/picard/build/libs/picard.jar.
#' @param bam [REQUIRED]  Path to BAM file.
#' @param ref_genome [REQUIRED]  Path to reference genome.
#' @param snpdb [REQUIRED] Known variant database.Requires atleast 1.
#' @param threads [REQUIRED]Number of threads to split the work. Only relevant if region_bed file is given.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @export


recalibrate_bq_gatk=function(bin_path="tools/samtools/samtools",bin_path2="tools/gatk/gatk",
bin_path3="tools/picard/build/libs/picard.jar",bam="",ref_genome="",snpdb="",
threads=3,ram=4,output_dir="",verbose=FALSE){

  out_file_dir=set_dir(dir=output_dir,name="recal_reports/recal_before")
  out_file_dir2=set_dir(dir=output_dir,name="recal_reports/recal_after")
  out_file_dir3=set_dir(dir=output_dir,name="recal_reports/recal_tmp")
  out_file_dir4=set_dir(dir=output_dir,name="recal_reports")


  parallel_generate_BQSR_gatk(bin_path=bin_path,bin_path2=bin_path2,bam=bam,
    ref_genome=ref_genome,snpdb=snpdb,
    threads=threads,output_dir=out_file_dir,verbose=verbose)

  parallel_apply_BQSR_gatk(bin_path=bin_path,bin_path2=bin_path2,bin_path3=bin_path3,
    bam=bam,ref_genome=ref_genome,rec_table=paste0(out_file_dir,"/",get_file_name(bam),".recal.table"),
    output_dir=out_file_dir3,verbose=verbose,threads=threads)

  system(paste(paste0("mv ",out_file_dir3,"/*"),out_file_dir4))
  system(paste("rm -rf ",out_file_dir3))

  sort_and_index_samtools(bin_path=bin_path,bam=paste0(out_file_dir4,"/",
  get_file_name(bam),".recal.",get_file_ext(bam)),output_dir=out_file_dir4,
  ram=ram,verbose=verbose,threads=threads,coord_sort=TRUE,index=TRUE)

  parallel_generate_BQSR_gatk(bin_path=bin_path,bin_path2=bin_path2,
    bam=paste0(out_file_dir4,"/",get_file_name(bam),".sorted.recal.",get_file_ext(bam)),
    ref_genome=ref_genome,snpdb=snpdb,threads=threads,output_dir=out_file_dir2,verbose=verbose)

  recal_covariates_gatk(bin_path=bin_path2,before=paste0(out_file_dir,"/",get_file_name(bam),".recal.table"),
    after=paste0(out_file_dir2,"/",get_file_name(bam),".recal.table"),output_dir=out_file_dir4)


  system(paste0("rm ",paste0(out_file_dir4,"/",get_file_name(bam),".recal.",get_file_ext(bam))))
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
#' @param ri Path to ribosomal intervals file. Only for RNAseq.
#' @param ref_flat Path to flat refrence. Only for RNAseq.
#' @param bi Bait capture target interval for panel data. Requires ti and off_target and on_tar arguments. Interval format.
#' @param ti Primary target intervals for panel data. Requires bi and off_tar and on_tar argmunets. Interval format.
#' @param mode Type of data to generate metrics for. Default tg. Options ["wgs","tg","rna"]
#' @export

alignment_qc_metrics=function(bin_path="tools/samtools/samtools",
  bin_path2="tools/picard/build/libs/picard.jar",bin_path3="tools/bedtools2/bin/bedtools",
  bam="",output_dir="",ref_genome="",verbose=FALSE,ram=4,tmp_dir=".",mapq=0,bi="",
  ti="",off_tar="",on_tar="",ri="",ref_flat="",mode="tg"){
    
    out_file_dir=set_dir(dir=output_dir,name="alignqc_report")
    
    ref=""
    if (!ref_genome==""){
      ref=paste0(" R=",ref_genome)
    }

  

    ## Generate alignment metrics



    bam_metrics_mapq_samtools(bin_path=bin_path,bam=bam,output_dir=out_file_dir,
    verbose=verbose,threads=threads)

    bam_metrics_summary_samtools(bin_path=bin_path2,bam=bam,output_dir=out_file_dir,
    verbose=verbose,threads=threads,tmp_dir=tmp_dir,ram=ram)
    
    bam_metrics_insertsize_picard(bin_path=bin_path2,bam=bam,output_dir=out_file_dir,
    verbose=verbose,tmp_dir=tmp_dir,ram=ram)


    ## Only call metrics for panel data if bait and target intervals are supplied, Otherwise WGS metrics.
    
    if (mode=="tg"){
   
      bam_metrics_tg_summary_picard(bin_path=bin_path2,bam=bam,output_dir=out_file_dir,
      verbose=verbose,tmp_dir=tmp_dir,ram=ram,bi=bi,ti=ti)
      ## Picard doesn't output coverage stats for off-target regions therefore we have to estimate this manually.

      ## I use this function to get the mean coverage on and off target.

      ## For target regions

      bed_coverage(bin_path=bin_path3,bam=bam,bed=on_tar,verbose=verbose,sorted=TRUE,
        mean=TRUE,hist=TRUE,fai=paste0(ref_genome,".fai"),suffix="on_Target",output_dir=out_file_dir)
      bed_coverage(bin_path=bin_path3,bam=bam,bed=on_tar,verbose=verbose,sorted=TRUE,
        mean=TRUE,hist=FALSE,fai=paste0(ref_genome,".fai"),suffix="on_Target",output_dir=out_file_dir)

      ## For off target regions

      bed_coverage(bin_path=bin_path3,bam=bam,bed=off_tar,verbose=verbose,sorted=TRUE,
        mean=TRUE,hist=TRUE,fai=paste0(ref_genome,".fai"),suffix="off_Target",output_dir=out_file_dir)
      bed_coverage(bin_path=bin_path3,bam=bam,bed=off_tar,verbose=verbose,sorted=TRUE,
        mean=TRUE,hist=FALSE,fai=paste0(ref_genome,".fai"),suffix="off_Target",output_dir=out_file_dir)

      ## Generate violin and cummulative plots for target and off target regions

      plot_coverage_panel(on_target=paste0(out_file_dir,"/coverage/",get_file_name(bam),".on_Target.Per_Region_Coverage.txt"),
      off_target=paste0(out_file_dir,"/coverage/",get_file_name(bam),col=c(5,4),height=6,width=12,output_dir=out_file_dir))
      plot_cumulative_cov(on_target=paste0(out_file,".on_Target.Histogram_Coverage.txt"),
      off_target=paste0(out_file,".off_Target.Histogram_Coverage.txt"),height=6,width=12,output_dir=out_file_dir)

    }else if (mode=="rna"){
        bam_metrics_rnaseq_summary_picard(bin_path=bin_path2,
        bam=bam,output_dir=out_file_dir,verbose=verbose,tmp_dir=tmp_dir,
        ram=ram,ri=ri,ref_flat=ref_flat)
    }else if(mode=="wgs"){
        bam_metrics_wgs_summary_picard(bin_path=bin_path2,
        bam=bam,output_dir=out_file_dir,verbose=verbose,tmp_dir=tmp_dir,ram=ram)
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
    exec_code=paste(bin_path,"view",bam," | head -n 1 | awk -F \"\t\" '{print $3}'")
    
    if (verbose){
      print(exec_code)
    }
    chr=system(exec_code,intern=TRUE)

    fmt=""

    if (format=="seg"){
      fmt="-s"
    }

     if (threads>1){
        parallel::mclapply(1:length(chrs),FUN=function(x){
          exec_code=print(paste(bin_path2,fmt,"--window", win,"--quality 20 --chromosome",
            paste0("chr",chrs[x],collapse=","), bam,">", paste0(out_file,".",x,".",format)))
          if (verbose){
            print(exec_code)
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
          print(exec_code)
        }
        system(exec_code)
      }
    
    if (grepl("chr",chr)){
     
      exec_code=paste("sed -i 's/chrom=chr/chrom=/g'",paste0(out_file,".",format))
      if (verbose){
          print(exec_code)
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
    print("Position and bed arguments are mutually exclusive")
    quit()
  }

  if(bed!=""){
    bed=paste("-L",bed)
  }

  exec_code=paste(bin_path,"view -b",bed,"-@",threads,bam,position,">", paste0(out_file_dir,"/",get_file_name(bam),".filtered.bam"))

  if (verbose){
    print(exec_code)
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
      print(exec_code)
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
      print("Arguments wigs and wigs_dir are mutually exclusive")
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
      print(exec_code)
    }
    system(exec_code)
    system(paste("rm",wig_list))

}
