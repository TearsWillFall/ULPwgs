

options(scipen=999)

fastqc=function (bin="./tools/FastQC/bin/fastqc",fastqc_R1="",fastqc_R2="",n_cores=3,output_dir="",verbose=FALSE){

  sep="/"

  if(output_dir==""){
    sep=""
  }

  sample_name=get_sample_name(fastqc_R1)

  if (!fastqc_R2==""){
  tmp_name=get_sample_name(fastqc_R2)
  sample_name=sapply(sapply(c(0:(nchar(tmp_name)-1)),function (i) substr(tmp_name,1,nchar(tmp_name)-i)),function (x) grepl(x,fastqc_R1))
  sample_name=which(l)[2]
  output_dir=paste0(output_dir,sep,sample_name)

  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  if(verbose){
    print(paste(bin,"-o ",output_dir,"-t ",n_cores,"--noextract",fastqc_R1,fastqc_R2))
  }
  system(paste(bin,"-o ",output_dir,"-t ",n_cores,"--noextract",fastqc_R1,fastqc_R2))

}else{
  output_dir=paste0(output_dir,sep,sample_name)

  if(!dir.exists(output_dir)){
    dir.create(output_dir)
  }
    if(verbose){
      print(paste(bin,"-o ",output_dir,"-t ",n_cores,"--noextract",fastqc_R1))
    }
    system(paste(bin,"-o ",output_dir,"-t ",n_cores,"--noextract",fastqc_R1))
  }
}



trimming=function(path="./tools/skewer/skewer",fastqc_R1="",fastqc_R2="",xadapt=NA,yadapt=NA,n_cores=3,output_dir="",verbose=FALSE){

  sep="/"

  if(output_dir==""){
    sep=""
  }

  sample_name=get_sample_name(fastqc_R1)
  if ((!is.na(xadapt)) & (!is.na(yadapt))){
    func=paste(path,"-m tail -t",n_cores,"-x", xadapt,"-y", yadapt)
  }
  else{
    func=paste(path,"-m tail -t",n_cores)
  }
  if (!fastqc_R2==""){
  tmp_name=get_sample_name(fastqc_R2)
  sample_name=sapply(sapply(c(0:(nchar(tmp_name)-1)),function (i) substr(tmp_name,1,nchar(tmp_name)-i)),function (x) grepl(x,fastqc_R1))
  sample_name=which(l)[2]
  output_dir=paste0(output_dir,sep,sample_name)

    if(verbose){
      print(paste(func,"-z -l 35 -f sanger --quiet -o",output_dir,fastqc_R1,fastqc_R2))
    }
    system(paste(func,"-z -l 35 -f sanger --quiet -o",output_dir,fastqc_R1,fastqc_R2))


  }else{
    output_dir=paste0(output_dir,sep,sample_name)

      if(verbose){
        print(paste(func,"-z -l 35 -f sanger --quiet -o",output_dir,fastqc_R1))
      }
      system(paste(func,"-z -l 35 -f sanger --quiet -o",output_dir,fastqc_R1))
    }
  }


merge=function(path="./tools/samtools/samtools",bam="",bam_dir="",verbose=FALSE){
    if(verbose){
      print(paste(path,"merge",bam, paste0(bam_dir,"*.bam")))
    }
    system(paste(path,"merge",bam, paste0(bam_dir,"*.bam")))
  }


alignment=function(path="./tools/bwa/bwa",path2="./tools/samtools/samtools",fastqc_R1="",fastqc_R2="",n_cores=3,ref_genome="",output_dir="",verbose=F){

    sep="/"

    if(output_dir==""){
      sep=""
    }

    sample_name=get_sample_name(fastqc_R1)
    out_file=paste0(output_dir,sep,sample_name,".bam")
    GPU=paste0("\"@RG\\tID:",sample_name,"\\tPL:ILLUMINA\\tPU:NA\\tLB:",sample_name,"\\tSM:",sample_name,"\"")

    if (!fastqc_R2==""){
    tmp_name=get_sample_name(fastqc_R2)
    sample_name=sapply(sapply(c(0:(nchar(tmp_name)-1)),function (i) substr(tmp_name,1,nchar(tmp_name)-i)),function (x) grepl(x,fastqc_R1))
    sample_name=which(l)[2]
    GPU=paste0("\"@RG\\tID:",sample_name,"\\tPL:ILLUMINA\\tPU:NA\\tLB:",sample_name,"\\tSM:",sample_name,"\"")
    out_file=paste0(output_dir,sep,sample_name,".bam")
      if(verbose){
          print(paste(path,"mem -t", n_cores," -v 2 -R",GPU,"-M",ref_genome, fastqc_R1,fastqc_R2, "|",path2," view -h -b >",out_file))
      }
      system(paste(path,"mem -t", n_cores," -v 2 -R",GPU,"-M",ref_genome, fastqc_R1,fastqc_R2, "| ",path2," view -h -b >",out_file))

      }
    else{
      if(verbose){
          print(paste(path,"mem -t", n_cores," -v 2 -R",GPU,"-M",ref_genome, fastqc_R1, "| ",path2," view -h -b >",out_file))
      }
      system(paste(path,"mem -t", n_cores," -v 2 -R",GPU,"-M",ref_genome, fastqc_R1, "| ",path2," view -h -b >",out_file))

    }
  }

sort_and_index=function(path="./tools/samtools/samtools",file="",output_dir="",verbose=FALSE){

  sep="/"

  if(output_dir==""){
    sep=""
  }

  sample_name=get_sample_name(file)
  file_ext=get_file_extension(file)
  out_file=paste0(output_dir,sep,"SORTED.",file_ext)


  if (!dir.exists(out_file)){
      dir.create(out_file)
    }
  out_file=paste0(out_file,"/",sample_name)


    if (verbose){
      print("Sorting BAM file:")
      print(paste0(path," sort ",file," -o ",out_file,".SORTED.",file_ext))
    }
    system(paste0(path," sort ",file," -o ",out_file,".SORTED.",file_ext))
    file=paste0(out_file,".SORTED.",file_ext)

    if (verbose){
      print("Indexing sorted BAM file:")
      print(paste(path," index",file))
    }
    system(paste(path," index",file))
    if (verbose){
      print("Generating Flag stats:")
      print(paste0(path," flagstat ",file," > ",paste0(out_file,".flagstat.txt")))
    }
    system(paste0(path," flagstat ",file," > ",paste0(out_file,".flagstat.txt")))
    if (verbose){
      print("Generating Index stats:")
      print(paste0(path," idxstats ",file," > ",paste0(out_file,".idxstats.txt")))
    }
    system(paste0(path," idxstats ",file," > ",paste0(out_file,".idxstats.txt")))
  }

remove_duplicates=function(path="tools/picard/build/libs/picard.jar",file="",output_dir="",verbose=FALSE,ref_genome=""){

    sep="/"

    if(output_dir==""){
      sep=""
    }

    sample_name=get_sample_name(file)
    file_ext=get_file_extension(file)

    out_file=paste0(output_dir,sep,"RMDUP.",file_ext)
    if (!dir.exists(out_file)){
        dir.create(out_file)
    }

    out_file=paste0(out_file,"/",sample_name)

    if(verbose){
      print(paste0("java -jar ",path," MarkDuplicates I=",file, " O=",paste0(out_file,".RMDUP.",file_ext)," M=",paste0(out_file,".picard_rmdup.txt")," REMOVE_DUPLICATES=true AS=true VALIDATION_STRINGENCY=LENIENT"))

    }
    system(paste0("java -jar ",path," MarkDuplicates I=",file, " O=",paste0(out_file,".RMDUP.",file_ext)," M=",paste0(out_file,".picard_rmdup.txt")," REMOVE_DUPLICATES=true AS=true VALIDATION_STRINGENCY=LENIENT"))
  }


  bam="SORTED.bam/SRR11742820-trimmed.SORTED.bam"
  ref_genome="/home/osvaldas/Workshop_low_pass/ref/GRCh37.p13.genome.fa.gz"





qc_metrics=function(path="./tools/samtools/samtools",path2="tools/picard/build/libs/picard.jar",bam="",output_dir="",verbose=FALSE,ref_genome=""){
    sep="/"

    if(output_dir==""){
      sep=""
    }

    sample_name=get_sample_name(bam)

    out_file=paste0(output_dir,sep,"alignQC")
    if (!dir.exists(out_file)){
        dir.create(out_file)
    }

    out_file=paste0(out_file,"/",sample_name)

    if (verbose){
      print("Generate MapQ distance map:")
      print(paste( path,"view",bam," | awk -F", "'\\t'", "'{c[$5]++} END { for (i in c) printf(\"%s\\t%s\\n\",i,c[i]) }'"," | sort -t$'\\t' -k 1 -g >>", paste0(out_file,".mapq_dist.txt")))

    }
    system(paste(path,"view",bam," | awk -F", "'\\t'", "'{c[$5]++} END { for (i in c) printf(\"%s\\t%s\\n\",i,c[i]) }'"," | sort -t$'\\t' -k 1 -g >>", paste0(out_file,".mapq_dist.txt")))

    if (verbose){
      print("Generate Alignment Metrics:")
      print(paste0("java -jar ",path2," CollectAlignmentSummaryMetrics R=",ref_genome," I=",bam," O=",paste0(out_file,".picard_summary.txt")))

    }
    system(paste0("java -jar ",path2," CollectAlignmentSummaryMetrics R=",ref_genome," I=",bam," O=",paste0(out_file,".picard_summary.txt")))

    if (verbose){
      print("Generate Insert Size Metrics:")
      print(paste0("java -jar ",path2," CollectInsertSizeMetrics R=",ref_genome," I=",bam," O=",paste0(out_file,".picard_insert_size.txt")," H=",paste0(out_file,".picard_insert_size.pdf")))

    }
    system(paste0("java -jar ",path2, " CollectInsertSizeMetrics R=",ref_genome," I=",bam," O=",paste0(out_file,".picard_insert_size.txt")," H=",paste0(out_file,".picard_insert_size.pdf")))

    if (verbose){
      print("Generate WGS Metrics for minimum MAPq=0:")
      print(paste0("java -jar ",path2," CollectWgsMetrics MINIMUM_MAPPING_QUALITY=0 R=",ref_genome," I=",bam," O=",paste0(out_file,".picard_wgs_q00.txt")))

    }
    system(paste0("java -jar ",path2," CollectWgsMetrics MINIMUM_MAPPING_QUALITY=0 R=",ref_genome," I=",bam," O=",paste0(out_file,".picard_wgs_q00.txt")))

    if (verbose){
      print("Generate WGS Metrics for minimum MAPq=20:")
      print(paste0("java -jar ",path2," CollectWgsMetrics MINIMUM_MAPPING_QUALITY=20 R=",ref_genome," I=",bam," O=",paste0(out_file,".picard_wgs_q20.txt")))

    }
    system(paste0("java -jar ",path2, " CollectWgsMetrics MINIMUM_MAPPING_QUALITY=20 R=",ref_genome," I=",bam," O=",paste0(out_file,".picard_wgs_q20.txt")))

    if (verbose){
      print("Generate WGS Metrics for minimum MAPq=37:")
      print(paste0("java -jar ",path2, " CollectWgsMetrics MINIMUM_MAPPING_QUALITY=37 R=",ref_genome," I=",bam," O=",paste0(out_file,".picard_wgs_q37.txt")))

    }
    system(paste0("java -jar ",path2, " CollectWgsMetrics MINIMUM_MAPPING_QUALITY=37 R=",ref_genome," I=",bam," O=",paste0(out_file,".picard_wgs_q37.txt")))

  }

read_counter=function(path="./tools/hmmcopy_utils/bin/readCounter",verbose=FALSE,win=500000,bam="",output_dir="",chrm="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY"){

    sep="/"

    if(output_dir==""){
      sep=""
    }

    sample_name=get_sample_name(bam)

    out_file=paste0(output_dir,sep,"wig")
    if (!dir.exists(out_file)){
        dir.create(out_file)
    }

    out_file=paste0(out_file,"/",sample_name)

      if (verbose){
        print(paste(path,"--window", win,"--quality 20 --chromosome",paste0("'",chrm,"'"), bam,">", paste0(out_file,".wig")))
      }
      system(paste(path,"--window", win,"--quality 20 --chromosome",paste0("'",chrm,"'"), bam,">" ,paste0(out_file,".wig")))

      if (grepl("chr",chrm)){
        if (verbose){
          print(paste("sed -i 's/chrom=chr/chrom=/g'",paste0(out_file,".wig")))
      }
      system(paste("sed -i 's/chrom=chr/chrom=/g'",paste0(out_file,".wig")))
    }
  }

ichorCNA=function(path="tools/ichorCNA/scripts/runIchorCNA.R",sample_id="",wig="",PLOIDY="2,3",TUMOR_CONTENT="0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9",HOMOZYGOUS_DEL="False",SUBCLONAL_STATES="NULL",gc="tools/ichorCNA/inst/extdata/gc_hg19_500kb.wig",map="tools/ichorCNA/inst/extdata/map_hg19_500kb.wig",centromere="tools/ichorCNA/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt",output_dir="",verbose=TRUE){

    sep="/"

    if(output_dir==""){
      sep=""
    }
    if(sample_id==""){
      sample_id=get_sample_name(wig)
    }

    out_file=paste0(output_dir,sep,"report")
    if (!dir.exists(out_file)){
        dir.create(out_file)
    }
    if(verbose){
      print(paste("Rscript",path,"--id",sample_id,"--WIG",wig,"--ploidy",paste0("'c(",PLOIDY,")'"),"--normal",paste0("'c(",TUMOR_CONTENT,")'"),"--maxCN 7 --gcWig", gc,"--mapWig",map,"--centromere",centromere,"--includeHOMD",HOMOZYGOUS_DEL,"--chrs 'c(1:22,\"X\")' --chrTrain \'c(1:22)\' --estimateNormal True --estimatePloidy True --estimateScPrevalence True --outDir",out_file))

    }
    system(paste("Rscript",path,"--id",sample_id,"--WIG",wig,"--ploidy",paste0("'c(",PLOIDY,")'"),"--normal",paste0("'c(",TUMOR_CONTENT,")'"),"--maxCN 7 --gcWig", gc,"--mapWig",map,"--centromere",centromere,"--includeHOMD",HOMOZYGOUS_DEL,"--chrs 'c(1:22,\"X\")' --chrTrain \'c(1:22)\' --estimateNormal True --estimatePloidy True --estimateScPrevalence True --outDir",out_file))

  }
