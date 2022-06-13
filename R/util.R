#' Get the filename of input file
#' This function takes the absolute/relative path to a file and
#' returns its base name without the file extension suffix.
#'
#' @param bam_path Path to the input file
#' @return A string with the name of the file
#' @export


get_file_name=function(file_path=""){
  filename=unlist(strsplit(basename(file_path),"\\."))[1]
  return(filename)
}

#' Set the current dir
#' This function sets and/or creates the directory
#' for a function
#'
#' @param dir Path to the set the directory
#' @param name Name of the directory
#' @return A path to the new directory
#' @export


set_dir=function(dir="",name=""){
  sep="/"

  if(dir==""){
    sep=""
  }

  new_dir=paste0(dir,sep,name)

  if (!dir.exists(new_dir)){
      dir.create(new_dir,recursive=TRUE)
  }
  return(new_dir)
}





#' Get the extension of a file
#'
#' This function takes the absolute/relative path to a file and
#' returns the file extension suffix.
#'
#' @param file_path Path to the input file
#' @return A string with the extension of the file
#' @export

get_file_ext=function(file_path=""){
    ext = strsplit(basename(file_path), split="\\.")[[1]]
    ext = paste(ext[-1],collapse=".")
    return(ext)
}

#' Get sample name from two sample replicates
#'
#' This function takes the absolute/relative path to two files
#' and returns the longest common string among their basenames
#'
#' @param file_path Path to the input file
#' @param file_path2 Path to the second input file
#' @return A string with the longest common basename
#' @export

intersect_file_name=function(file_path="",file_path2=""){
  tmp_name=get_file_name(file_path2)
  sample_name=sapply(sapply(c(0:(nchar(tmp_name)-1)),
  function (i) substr(tmp_name,1,nchar(tmp_name)-i)),function (x) grepl(x,file_path))
  sample_name=names(which(sample_name)[1])
  sample_name=sub("(.*)[_.-].*","\\1",sample_name)
  return(sample_name)

}


#' Index reference genome
#'
#' This function indexes a reference genome
#'
#' @param file Path to the input file with the reference genome in FASTA format.
#' @param bin_path Path to bwa executable. Default path tools/bwa/bwa.
#' @param verbose Enables progress messages. Default False.
#' @export


index_ref_bwa=function(bin_path="tools/bwa/bwa",file="",verbose=FALSE){
exec_code=paste(bin_path,"index", file)

  if(verbose){
    print(exec_code)
  }
  error=system(exec_code)
  if(error!=0){
    stop("bwa failed to run due to unknown error.
    Check std error for more information.")
  }

}




#' Sort a BAM file
#'
#' This function sorts a genome sequence file (BAM/SAM)
#'
#' @param bam Path to the input file with the sequence.
#' @param bin_path Path to bwa executable. Default path tools/samtools/samtools.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param threads Number of threads. Default 3
#' @param coord_sort Generate a coord sorted file. Otherwise queryname sorted. Default TRUE
#' @param ram Ram memory to use per thread in GB. Default 1GB
#' @export

bam_sort_samtools=function(bin_path="tools/samtools/samtools",bam="",output_dir="",ram=1,
verbose=FALSE,threads=3,coord_sort=TRUE){

  out_file_dir=set_dir(dir=output_dir,name="sorted")

  sort_type=""
  
  if(!coord_sort){
    sort_type=" -n "
  }
  exec_code=paste0(bin_path," sort ",sort_type, bam," -@ ",threads," -m ",ram,"G"," -o ",
  paste0(out_file_dir,"/",get_file_name(bam),".sorted.",get_file_ext(bam)))
  if (verbose){
    print(exec_code)
  }
  error=system(exec_code)
  if(error!=0){
    stop("samtools failed to run due to unknown error.
    Check std error for more information.")
  }
}


#' Generate BAM file flag and index stats
#'
#'
#' @param bam Path to the input file with the sequence.
#' @param bin_path Path to bwa executable. Default path tools/samtools/samtools.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param threads Number of threads. Default 3
#' @param stats Generate BAM stats. Default all. Options ["all","flag","index",""]
#' @export


bam_stats_samtools=function(bin_path="tools/samtools/samtools",bam="",output_dir="",
verbose=FALSE,threads=3,stats="all"){


    out_file_dir=set_dir(dir=output_dir,name="stats")

    if(stats=="all"|stats=="flag"){
      bam_stats_flag(bin_path=bin_path,bam=bam,output_dir=out_file_dir,
      verbose=verbose,threads=threads)
    }
    if(stats=="all"|stats=="index"){
      bam_stats_index(bin_path=bin_path,bam=bam,output_dir=out_file_dir,
      verbose=verbose,threads=threads)
    }
}

#' Generate BAM file flagstats
#'
#'
#' @param bam Path to the input file with the sequence.
#' @param bin_path Path to bwa executable. Default path tools/samtools/samtools.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param threads Number of threads. Default 3
#' @export


bam_stats_flag=function(bin_path="tools/samtools/samtools",bam="",output_dir="",
verbose=FALSE,threads=3){

  out_file_dir=set_dir(dir=output_dir,name="flag")

  exec_code=paste0(bin_path," flagstat ",bam," -@ ",threads," > ",
    paste0(out_file_dir,"/",get_file_name(bam),".flagstat.txt"))
  
  if (verbose){
    print(exec_code)
  }


   error=system(exec_code)
  if(error!=0){
    stop("samtools failed to run due to unknown error.
    Check std error for more information.")
  }
}

#' Generate BAM file indexstats
#'
#'
#' @param bam Path to the input file with the sequence.
#' @param bin_path Path to bwa executable. Default path tools/samtools/samtools.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param threads Number of threads. Default 3
#' @export


bam_stats_index=function(bin_path="tools/samtools/samtools",bam="",output_dir="",
verbose=FALSE,threads=3){

  out_file_dir=set_dir(dir=output_dir,name="index")

  exec_code=paste0(bin_path," idxstats ",bam," > ",paste0(out_file_dir,"/",get_file_name(bam),".idxstats.txt"))
  if (verbose){
    print(exec_code)
  }
     error=system(exec_code)
  if(error!=0){
    stop("samtools failed to run due to unknown error.
    Check std error for more information.")
  }

}

#' Generate BAM MapQ metrics
#'
#'
#' @param bam Path to the input file with the sequence.
#' @param bin_path Path to bwa executable. Default path tools/samtools/samtools.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param threads Number of threads. Default 3
#' @export

bam_metrics_mapq_samtools=function(bin_path="tools/samtools/samtools",bam="",output_dir="",
verbose=FALSE,threads=3){


  out_file_dir=set_dir(dir=output_dir,name="mapq")


  exec_code=paste(bin_path,"view",bam," -@ ",threads, " | awk -F", "'\\t'",
    "'{c[$5]++} END { for (i in c) printf(\"%s\\t%s\\n\",i,c[i]) }'",
    " | sort -t$'\\t' -k 1 -g >>", paste0(out_file_dir,"/",get_file_name(bam),".mapq_dist.txt"))
  
  if (verbose){
    print(exec_code)
  }

    error=system(exec_code)
  if(error!=0){
    stop("samtools failed to run due to unknown error.
    Check std error for more information.")
  }
}


#' Generate BAM General Summary Metrics
#'
#'
#' @param bam Path to the input file with the sequence.
#' @param bin_path Path to bwa executable. Default path tools/samtools/samtools.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param ram RAM memory to use in GB. Default 4.
#' @param tmp_dir Path to TMP directory. Default .
#' @export

bam_metrics_summary_picard=function(bin_path="tools/picard/build/libs/picard.jar",bam="",output_dir="",
verbose=FALSE,threads=3,tmp_dir=".",ram=4){

  out_file_dir=set_dir(dir=output_dir,name="summary")

  tmp=""
  if (!tmp_dir==""){
    tmp=paste0(" TMP_DIR=",tmp_dir)
  }

  exec_code=paste0("java -Xmx",ram,"g", " -Djava.io.tmpdir=",tmp_dir,
        " -jar ",bin_path," CollectAlignmentSummaryMetrics ",
        "VALIDATION_STRINGENCY=SILENT I=",bam," O=",paste0(out_file_dir,"/",get_file_name(bam),".picard_summary.txt "),tmp)
  
  if (verbose){
        print(exec_code)

  }
  error=system(exec_code)
  if(error!=0){
    stop("samtools failed to run due to unknown error.
    Check std error for more information.")
  }

}



#' Generate BAM Insert Size Metrics
#'
#'
#' @param bam Path to the input file with the sequence.
#' @param bin_path Path to bwa executable. Default path tools/samtools/samtools.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param ram RAM memory to use in GB. Default 4.
#' @param tmp_dir Path to TMP directory. Default .
#' @export

bam_metrics_insertsize_picard=function(bin_path="tools/picard/build/libs/picard.jar",
bam="",output_dir="",verbose=FALSE,tmp_dir=".",ram=4){

  out_file_dir=set_dir(dir=output_dir,name="insertsize")


  tmp=""
  if (!tmp_dir==""){
    tmp=paste0(" TMP_DIR=",tmp_dir)
  }

  exec_code=paste0("java -Xmx",ram,"g", " -Djava.io.tmpdir=",tmp_dir," -jar ",
      bin_path," CollectInsertSizeMetrics ","VALIDATION_STRINGENCY=SILENT I=",
      bam," O=",paste0(out_file_dir,"/",get_file_name(bam),".picard_insert_size.txt")," H=",
      paste0(out_file_dir,"/",get_file_name(bam),".picard_insert_size.pdf "),tmp)

  if(verbose){
      print(exec_code)
  }
    error=system(exec_code)
  if(error!=0){
    stop("picard failed to run due to unknown error.
    Check std error for more information.")
  }

}




#' Generate BAM Summary for Targeted data
#'
#'
#' @param bam Path to the input file with the sequence.
#' @param bin_path Path to bwa executable. Default path tools/samtools/samtools.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param ram RAM memory to use in GB. Default 4.
#' @param tmp_dir Path to TMP directory. Default .
#' @param bi Bait capture target interval for panel data. Interval format.
#' @param ti Primary target intervals for panel data. Interval format.
#' @export

bam_metrics_tg_summary_picard=function(bin_path="tools/picard/build/libs/picard.jar",bam="",output_dir="",
verbose=FALSE,tmp_dir=".",ram=4,bi="",ti=""){

  out_file_dir=set_dir(dir=output_dir,name="summary")


  tmp=""
  if (!tmp_dir==""){
    tmp=paste0(" TMP_DIR=",tmp_dir)
  }

  exec_code=print(paste0("java -Xmx",ram,"g", " -Djava.io.tmpdir=",tmp_dir,
        " -jar ",bin_path," CollectHsMetrics VALIDATION_STRINGENCY=SILENT BI=",
        bi," TI=",ti," I=",bam," THEORETICAL_SENSITIVITY_OUTPUT=",
        paste0(out_file_dir,"/",get_file_name(bam),".picard_TS.txt"),ref," O=",
        paste0(out_file_dir,"/",get_file_name(bam),".picard_CollectHSmetrics.txt "),tmp))

   if (verbose){
        print(exec_code)
      }
       error=system(exec_code)
  if(error!=0){
    stop("picard failed to run due to unknown error.
    Check std error for more information.")
  }


}


#' Generate BAM Summary for RNAseq data
#'
#'
#' @param bam Path to the input file with the sequence.
#' @param bin_path Path to bwa executable. Default path tools/samtools/samtools.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param ram RAM memory to use in GB. Default 4.
#' @param tmp_dir Path to TMP directory. Default .
#' @param ri Ribosomal interval for RNAseq data. Interval format.
#' @param ref_flat Path to flat refrence. Interval format.
#' @export

bam_metrics_rnaseq_summary_picard=function(bin_path="tools/picard/build/libs/picard.jar",
bam="",output_dir="",verbose=FALSE,tmp_dir=".",ram=4,ri="",ref_flat=""){

  out_file_dir=set_dir(dir=output_dir,name="summary")

  tmp=""
  if (!tmp_dir==""){
    tmp=paste0(" TMP_DIR=",tmp_dir)
  }

  exec_code=paste0("java -Xmx",ram,"g", " -Djava.io.tmpdir=",tmp_dir," -jar ",bin_path,
        " CollectRnaSeqMetrics VALIDATION_STRINGENCY=SILENT STRAND_SPECIFICITY='NONE' REF_FLAT=",
         ref_flat, " RIBOSOMAL_INTERVALS=",ri,
         " I=",bam," O=",paste0(out_file_dir,"/",get_file_name(bam),".CollectRNAseqMetrics.txt "),tmp)

 if (verbose){

      print(exec_code)

  }
      error=system(exec_code)
  if(error!=0){
    stop("picard failed to run due to unknown error.
    Check std error for more information.")
  }
}



#' Generate BAM Summary for WGS data
#'
#'
#' @param bam Path to the input file with the sequence.
#' @param bin_path Path to bwa executable. Default path tools/samtools/samtools.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param ram RAM memory to use in GB. Default 4.
#' @param tmp_dir Path to TMP directory. Default .
#' @export

bam_metrics_wgs_summary_picard=function(bin_path="tools/picard/build/libs/picard.jar",
bam="",output_dir="",verbose=FALSE,tmp_dir=".",ram=4){

 out_file_dir=set_dir(dir=output_dir,name="summary")


  tmp=""
  if (!tmp_dir==""){
    tmp=paste0(" TMP_DIR=",tmp_dir)
  }

  exec_code=paste0("java -Xmx",ram,"g", " -Djava.io.tmpdir=",tmp_dir," -jar ",
            bin_path," CollectWgsMetrics VALIDATION_STRINGENCY=SILENT MINIMUM_MAPPING_QUALITY=",
            mapq," I=",bam," O=",paste0(out_file_dir,"/",get_file_name(bam),".picard_wgs_q00.txt "),tmp)

 if (verbose){
      print(exec_code)
 }
      error=system(exec_code)
  if(error!=0){
    stop("picard failed to run due to unknown error.
    Check std error for more information.")
  }
}



#' Index a BAM file
#'
#' This function indexes a genomic sequence file (BAM/SAM).
#'
#' @param bam Path to the input file with the sequence.
#' @param bin_path Path to bwa executable. Default path tools/samtools/samtools.
#' @param verbose Enables progress messages. Default False.
#' @param threads Number of threads. Default 3
#' @export

bam_index_samtools=function(bin_path="tools/samtools/samtools",bam="",verbose=FALSE,threads=3){

  exec_code=paste(bin_path," index",bam," -@ ",threads)
  if (verbose){
    print(exec_code)
  }
     error=system(exec_code)
  if(error!=0){
    stop("samtools failed to run due to unknown error.
    Check std error for more information.")
  }
}


#' Create BED with antitarget regions with padding
#'
#' This function indexes a genomic sequence file (BAM/SAM).
#'
#' @param bin_path Path to bwa executable. Default path tools/bedtools2/bin/bedtools.
#' @param bed Path to the input file with the sequence.
#' @param pad Pad distance. Default 10
#' @param output_name Output file name.
#' @param genome Path to genome fa.fai file
#' @param verbose Enables progress messages. Default False.
#' @export

complement_bed=function(bin_path="tools/bedtools2/bin/bedtools",bed="",pad=10,
output_name="Complement",genome="",verbose=FALSE){

  if (pad!=0){

    pad_bed(bin_path=bin_path,bed=bed,pad=pad,output_name=paste0(ULPwgs::get_file_name(bed),"_",pad),
    genome=genome,verbose=verbose)

    exec_code=paste0(bin_path," complement -i ",paste0(ULPwgs::get_file_name(bed),"_",pad,".bed"),
       " -g ", genome, " > ",paste0(output_name,".bed"))
    if(verbose){
      print(exec_code)
    }
         error=system(exec_code)
  if(error!=0){
    stop("bedtools failed to run due to unknown error.
    Check std error for more information.")
  }

  }else{
    exec_code=paste0(bin_path," complement -i ",bed, " -g ", genome, " > ",paste0(output_name,".bed"))
    if(verbose){
      print(exec_code)
    }
          error=system(exec_code)
  if(error!=0){
    stop("bedtools failed to run due to unknown error.
    Check std error for more information.")
  }
  }
}

#' Pad a BED file
#'
#' This function takes a BED file and pad each regions in both directions
#'
#' @param bed Path to the input file with the sequence.
#' @param bin_path Path to bwa executable. Default path tools/bedtools2/bin/bedtools.
#' @param pad Pad distance. Default 10
#' @param output_name Output file name.
#' @param genome Path to genome fa.fai file
#' @param verbose Enables progress messages. Default False.
#' @export

pad_bed=function(bin_path="tools/bedtools2/bin/bedtools",bed="",pad=10,
output_name="Padded",genome="",verbose=FALSE){

  exec_code=paste0(bin_path," slop -i ",bed, " -g ", genome," -b ",pad, " > ",
    paste0(output_name,".bed"))
  if(verbose){
    print(exec_code)
  }
  error=system(exec_code)
  if(error!=0){
    stop("bedtools failed to run due to unknown error.
    Check std error for more information.")
  }
}



#' Wrapper of BaseRecalibrator function of gatk
#'
#' Generates a recalibration table based on various covariates.
#' This function wraps around gatk BaseRecalibrator function.
#' For more information about this function: https://gatk.broadinstitute.org/hc/en-us/articles/360036898312-BaseRecalibrator
#'
#' @param bam [REQUIRED] Path to the BAM file.
#' @param bin_path [REQUIRED] Path to gatk executable. Default tools/gatk/gatk.
#' @param ref_genome [REQUIRED] Path to reference genome
#' @param dbsnp [REQUIRED] Path to known snp positions in VCF format. Multiple vcf can be supplied as a vector.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @export

generate_BQSR_gatk=function(region="",bin_path="tools/gatk/gatk",bam="",ref_genome="",
dbsnp="",output_dir="",verbose=FALSE){

  out_file_dir=set_dir(dir=output_dir)

  reg=""
  if (region==""){
      out_file=paste0(out_file_dir,get_file_name(bam),".recal.table")
  }else{
      reg=paste0(" -L ",region)
      out_file=paste0(out_file_dir,get_file_name(bam),".",region,".recal.table")
  }

  ## Multiple vcf with snps can be given

  if (dbsnp!=""){
    dbsnp=paste(" --known-sites ",dbsnp,collapse=" ")
  }

  exec_code=paste0(bin_path," BaseRecalibrator -I ",bam, " -R ", ref_genome,dbsnp,
  reg," -O ",out_file)

  if(verbose){
    print(exec_code)
  }
  error=system(exec_code)
  if(error!=0){
    stop("gatk failed to run due to unknown error.
    Check std error for more information.")
  }
}

#' Multiregion parallelization of generate_BQSR function
#'
#' Generates a recalibration table based on various covariates.
#' This function wraps around gatk BaseRecalibrator function.
#' For more information about this function: https://gatk.broadinstitute.org/hc/en-us/articles/360036898312-BaseRecalibrator
#'
#' @param bam [REQUIRED] Path to the BAM file.
#' @param bin_path [REQUIRED] Path to gatk executable. Default tools/samtools/samtools.
#' @param bin_path2 [REQUIRED] Path to gatk executable. Default tools/gatk/gatk.
#' @param ref_genome [REQUIRED] Path to reference genome
#' @param dbsnp [REQUIRED] Path to known snp positions in VCF format. Multiple vcf can be supplied as a vector.
#' @param threads Number of threads to split the work. Default 3
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param ram [OPTIONAL] If batch mode. RAM memory in GB per job. Default 1
#' @param update_time [OPTIONAL] If batch mode. Show job updates every update time. Default 60
#' @export
#' @import pbapply



parallel_generate_BQSR_gatk=function(bin_path="tools/samtools/samtools",
bin_path2="tools/gatk/gatk",bam="",ref_genome="",dbsnp="",threads=3,
output_dir="",verbose=FALSE,mode="local",time="48:0:0",ram=1,update_time=60){

  options(scipen = 999)

  out_file_dir=set_dir(dir=output_dir)

  dat=get_bam_reference_chr(bin_path=bin_path,bam=bam,verbose=verbose)
  dat$start=dat$start+1
  dat=dat %>% dplyr::mutate(Region=paste0(chr,":",start,"-",end))

  fun <- system.file("shell", "generate_BQSR_gatk.sh", package = "ULPwgs")
    
  

  task_id=sample(1:10000000,1)
  task_name="generate_BQSR"
  input_name=get_file_name(bam)

  job_name=paste0(c(task_id,task_name, input_name),collapse="_")

  parallel::mclapply(1:nrow(dat),FUN=function(x){
    tmp=dat[x,]
    if(mode=="local"){
        generate_BQSR_gatk(region=tmp$Region,
        bin_path=bin_path2,bam=bam,ref_genome=ref_genome,dbsnp=dbsnp,
        output_dir=out_file_dir,verbose=verbose)


      }else if (mode=="batch"){
        batch_name=paste0(c(tmp$chr,tmp$start,tmp$end),collapse="_")
        exec_code=paste("qsub -N ",paste0(c(job_name,batch_name),collapse="_"),paste0(" -l h_rt=",time),
        paste0(" -l mem=",ram,"G"), paste0(" -pe smp 2"), paste0(" -wd ."),
         fun, tmp$Region, bin_path,
         bam, ref_genome, dbsnp, out_file_dir,verbose)
        if(verbose){
          print(exec_code)
        }

        error=system(exec_code)
        if(error!=0){
          stop("gatk failed to run due to unknown error.
          Check std error for more information.")
        }

    }else{
      stop("Wrong Mode supplied. Available modes are ['local','batch']")
    }
  })

  if(mode=="batch"){
    batch_job_validator(job=job_name,time=update_time,verbose=verbose)
  }
   

  gather_BQSR_reports_gatk(bin_path=bin_path2,reports_dir=out_file_dir,output_dir=out_file_dir,
  output_name=get_file_name(bam),verbose=verbose)
  system(paste0("rm ",out_file_dir,"/*:*.recal.table"))
}




#' Wrapper around gatk GatherBQSRReports function
#'
#' This functions collects the Recalibration reports generated from scattered parallel_generate_BQSR output
#' This function wraps around gatk GatherBQSRReports function.
#' For more information about this function: https://gatk.broadinstitute.org/hc/en-us/articles/360036829851-GatherBQSRReports
#'
#' @param bin_path [REQUIRED] Path to gatk executable. Default tools/gatk/gatk.
#' @param reports_dir [REQUIRED] Path to the directory where reports are stored.
#' @param output_name [OPTIONAL] Name for the output report file.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @export

gather_BQSR_reports_gatk=function(bin_path="tools/gatk/gatk",reports_dir="",
output_name="Report",output_dir="",verbose=FALSE){

  out_file_dir=set_dir(dir=output_dir)
  files=list.files(reports_dir,full.names=TRUE,pattern=":")


  exec_code=paste0(bin_path," GatherBQSRReports ",paste(" -I ",files,collapse=" "),
    " -O ",paste0(out_file_dir,output_name,".recal.table"))
  if(verbose){
    print(exec_code)
  }
  error=system(exec_code)
  if(error!=0){
    stop("gatk failed to run due to unknown error.
    Check std error for more information.")
  }
}



#' Wrapper of applyBQSR function gatk
#'
#' Applies numerical corrections to each individual basecall based on the covariates analyzed before.
#' This function wraps around gatk applyBQSR  function.
#' For more information about this function: https://gatk.broadinstitute.org/hc/en-us/articles/360050814312-ApplyBQSR
#'
#' @param bam [REQUIRED] Path to the BAM file.
#' @param bin_path [REQUIRED] Path to gatk executable. Default tools/gatk/gatk.
#' @param ref_genome [REQUIRED] Path to reference genome
#' @param rec_table [REQUIRED] Path to covariates table generated by generate_BSQR.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @export



apply_BQSR_gatk=function(region="",bin_path="tools/gatk/gatk",bam="",ref_genome="",
rec_table="",output_dir="",verbose=FALSE){

  out_file_dir=set_dir(dir=output_dir)
 
  reg=""
  if (region==""){
      out_file=paste0(" ", out_file_dir,"/", get_file_name(bam),".recal.",get_file_ext(bam))
  }else{
      reg=paste0(" -L ",strsplit(region,"_")[[1]][2], " ")
      out_file=paste0(out_file_dir,"/", get_file_name(bam),".",region,".recal.",get_file_ext(bam))
  }
  exec_code=paste(bin_path," ApplyBQSR -I ",bam, " -R ", ref_genome,
    " --bqsr-recal-file ",rec_table, " -O ",out_file,reg)
  if(verbose){
    print(exec_code)
  }
  error=system(exec_code)
  if(error!=0){
    stop("gatk failed to run due to unknown error.
    Check std error for more information.")
  }
}


#' Multiregion parallelization of apply_BQSR function
#'
#' Recalibrates
#' Applies numerical corrections to each individual basecall based on the covariates analyzed before.
#' For more information about this function: https://gatk.broadinstitute.org/hc/en-us/articles/360050814312-ApplyBQSR
#'
#' @param bam [REQUIRED] Path to the BAM file.
#' @param bin_path [REQUIRED] Path to samtools executable. Default tools/samtools/samtools.
#' @param bin_path2 [REQUIRED] Path to gatk executable. Default tools/gatk/gatk.
#' @param bin_path3 [REQUIRED] Path to picard executable. Default tools/picard/build/libs/picard.jar
#' @param ref_genome [REQUIRED] Path to reference genome
#' @param rec_table [REQUIRED] Path to the recalibratio table.
#' @param threads [OPTIONAL] Number of threads to split the work. Default 4
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param ram [OPTIONAL] If batch mode. RAM memory in GB per job. Default 1
#' @param update_time [OPTIONAL] If batch mode. Show job updates every update time. Default 60
#' @export
#' @import pbapply

parallel_apply_BQSR_gatk=function(bin_path="tools/samtools/samtools",bin_path2="tools/gatk/gatk",
bin_path3="tools/picard/build/libs/picard.jar",bam="",ref_genome="",rec_table="",
output_dir="",verbose=FALSE,threads=4,mode="local",time="48:0:0",ram=4,update_time=60){

  options(scipen = 999)

  out_file_dir=set_dir(dir=output_dir)
  dat=get_bam_reference_chr(bin_path=bin_path,bam=bam,verbose=verbose)
  dat$start=dat$start+1
  dat$pos=1:nrow(dat)
  dat=dat %>% dplyr::mutate(Region=paste0(pos,"_",chr,":",start,"-",end))


  task_id=sample(1:10000000,1)
  task_name="apply_BQSR"
  input_name=get_file_name(bam)

  job_name=paste0(c(task_id,task_name, input_name),collapse="_")
 
  parallel::mclapply(1:nrow(dat),FUN=function(x){
    tmp=dat[x,]
    if(mode=="local"){
      apply_BQSR_gatk(region=tmp$Region,
      bin_path=bin_path2,bam=bam,ref_genome=ref_genome,
      rec_table=rec_table, output_dir=out_file_dir,verbose=verbose)
    }else if (mode=="batch"){
      batch_name=paste0(c(tmp$chr,tmp$start,tmp$end),collapse="_")
      exec_code=paste("qsub -N ",paste0(c(job_name,batch_name),collapse="_"),
          paste0(" -l h_rt=",time), paste0(" -l mem=",ram,"G"), paste0(" -pe smp 2"), paste0(" -wd ."),
          fun, tmp$Region, bin_path,
          bam, ref_genome, rec_table, out_file_dir,verbose)
          if(verbose){
            print(exec_code)
          }

          error=system(exec_code)
          if(error!=0){
            stop("gatk failed to run due to unknown error.
            Check std error for more information.")
          }

    }else{
      stop("Wrong Mode supplied. Available modes are ['local','batch']")
    } 

    
  },mc.cores=threads)
    
  if(mode=="batch"){
    batch_job_validator(job=job_name,time=update_time,verbose=verbose)
  }
    
  

  gather_bam_files(bin_path=bin_path3,bams_dir=out_file_dir,output_dir=out_file_dir,
  output_name=paste0(get_file_name(bam),".recal.sorted.rmdup.sorted"))
  system(paste0("rm ", out_file_dir,"/*:*.recal*.ba*"))

}


#' Validate job submited on batch on sungrid based cluster
#'
#' @param job Name of job or jobs.
#' @param time Time in seconds between checks. Default 10.
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @export


batch_job_validator=function(job="",time=10,verbose=FALSE){
  error=FALSE
  col_names=c("job_id","job_priority","job_name","user","status","start_time","cores")
  exec_code="qstat -xml | tr '\n' ' ' | sed 's#<job_list[^>]*>#\\n#g'   | sed 's#<[^>]*>##g' | grep \" \" | column -t"
  tryCatch({
      dat_info=read.table(text=system(exec_code,intern=TRUE))},error=function(x){
        return()
      }
  )

  names(dat_info)=col_names
  while(nrow(dat_info)!=0& !error){
    if(verbose){
          print(dat_info)
    }
    dat_info=read.table(text=system(exec_code,intern=TRUE))
    names(dat_info)=col_names
    if(any(!grepl("qw",dat_info$status)&!grepl("r",dat_info$status))){
      error=TRUE
    }
    Sys.sleep(time)
  }
  if(error){
    stop("One or more jobs failed")
  }
  return()
}



#' Wrapper around gatk GatherBamFiles function
#'
#' This functions collects the Recalibration reports generated from scattered parallel_apply_BQSR output
#' This function wraps around gatk GatherBamFiles function.
#' For more information about this function: https://gatk.broadinstitute.org/hc/en-us/articles/360037055512-GatherBamFiles-Picard-
#'
#' @param bin_path [REQUIRED] Path to gatk executable. Default tools/picard/build/libs/picard.jar.
#' @param bams_dir [REQUIRED] Path to the directory where BAM files are stored.
#' @param output_name [OPTIONAL] Name for the output file name.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @export

gather_bam_files=function(bin_path="tools/picard/build/libs/picard.jar",bams_dir="",
  output_name="File",output_dir="",verbose=FALSE){

  out_file_dir=set_dir(dir=output_dir)

  files=list.files(bams_dir,full.names=TRUE,pattern=":")
  files=files[grepl("bam$",files)]
  files=files[order(as.numeric(lapply(lapply(lapply(lapply(lapply(lapply(basename(files),
  FUN=strsplit,split="\\."),FUN="[[",index=1),FUN="[",index=2),FUN=strsplit,split="_"),FUN="[[",index=1),FUN="[",index=1)))]

  exec_code=paste0("java -jar ",bin_path," GatherBamFiles ",
    paste0(" I=",files,collapse=" ")," O=",paste0(out_file_dir,"/",output_name,".bam"))


  if(verbose){
    print(exec_code)
  }

  error=system(exec_code)




  if(error!=0){
    stop("picard failed to run due to unknown error.
    Check std error for more information.")
  }
}


#' Wrapper around samtools addreplacerg function
#'
#' This functions add/replaces RG tags lines in BAM files
#' This function wraps around samtools addreplacerg function
#'
#' @param bin_path [REQUIRED] Path to santools binary. Default tools/samtools/samtools.
#' @param bam [REQUIRED] Path to the BAM file/s.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param index [OPTIONAL] Generate an indexed file. Default False.
#' @param ID [REQUIRED] ID tag for RG tag line.
#' @param PL [OPTIONAL] PL tag for RG tag line.
#' @param PU [OPTIONAL] PU tag for RG tag line.
#' @param LB [OPTIONAL] LB tag for RG tag line.
#' @param SM [OPTIONAL] SM tag for RG tag line.
#' @param threads [OPTIONAL] Number of threads per jobs.
#' @param jobs [OPTIONAL] Number of jobs to run.
#' @export

replace_rg=function(bin_path="tools/samtools/samtools",bam="",output_dir="",
  verbose=FALSE,index=TRUE,ID="",PL="",PU="",LB="",SM="",threads=3,jobs=1){

  out_file_di=set_dir(dir=output_dir)

  if(ID!=""){
    ID=paste0("ID:",ID)
  }

  if(PL!=""){
    PL=paste0("PL:",PL)
  }

  if(PU!=""){
    PU=paste0("PU:",PU)
  }

  if(LB!=""){
    LB=paste0("LB:",LB)
  }

  if(SM!=""){
    SM=paste0("SM:",SM)
  }

  tag=c(ID,PL,PU,LB,SM)
  tag=tag[!tag==""]

  tag=paste0(" -r ",paste0(tag,collapse=" -r "))

  parallel::mclapply(1:length(bam),FUN=function(x){
    exec_code=paste(bin_path," addreplacerg ",tag," -o ",paste0(out_file_dir,"/",
      basename(sub("bam","rh.bam",bam[x]))), " -@ ",threads,bam[x])
  if(verbose){
      print(exec_code)
    }
    exec_code

    if(index){
      bam_index_samtools(bin_path=bin_path,bam=paste0(out_file_dir,"/",
      basename(sub("bam","rh.bam",bam[x]))),verbose=verbose,threads=threads)
    }
  },mc.cores=jobs)
}


#' Wrapper of AnalyzeCovariates function in gatk
#'
#' Generates a report of the recalibrated values.
#' This function wraps around gatk AnalyzeCovariates function.
#' For more information about this function: https://gatk.broadinstitute.org/hc/en-us/articles/360037066912-AnalyzeCovariates
#'
#' @param bin_path [REQUIRED] Path to gatk executable. Default tools/gatk/gatk.
#' @param before [REQUIRED] Recalibration table produced by generate_BQSR function.
#' @param after [OPTIONAL] Recalibration table produced by generate_BQSR function.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @export


recal_covariates=function(bin_path="tools/gatk/gatk",before="",after="",output_dir="",verbose=FALSE){

  out_file_dir=set_dir(name=output_dir)
 

  
  if (before!="" & after==""){
    exec_code=paste0(bin_path," AnalyzeCovariates -bqsr ",before," -plots ",
    paste0(out_file_dir,"/",get_file_name(before),"_covariates_analysis_before.pdf"))
    if(verbose){
      print(exec_code)
    }
  error=system(exec_code)
  if(error!=0){
    stop("gatk failed to run due to unknown error.
    Check std error for more information.")
  }
  }else{
    exec_code=paste0(bin_path," AnalyzeCovariates -before ",before," -after ",after,
      " -plots ",paste0(out_file_dir,"/",get_file_name(before),"_covariates_analysis.pdf"))
    if(verbose){
      print(exec_code)
    }
        error=system(exec_code)
  if(error!=0){
    stop("gatk failed to run due to unknown error.
    Check std error for more information.")
  }

  }


}


#' Estimate coverage for BED file
#'
#' This function estimates mean/per_base_coverage coverage for regions in bed file.
#' This function is intended for estimating off target mean coverage in panel data, however it can be used as a standalone.
#' It is highly recommended to use the sorted option and sort the BED and BAM files beforehand as it highly decreases RAM
#' consumption and increases computation speed. The easiest way of doing this is using sort -k1,1 -k2,2n on the bed file.
#' However, this may not work if the BED has been build based on different reference than the BAM. For example, hs37 vs hg19.
#' For more information: https://bedtools.readthedocs.io/en/latest/content/tools/coverage.html
#'
#'
#' @param bam Path to the input BAM file.
#' @param bin_path Path to bwa executable. Default tools/bedtools2/bin/bedtools.
#' @param bed Path to the input bed file.
#' @param sorted Are the input files sorted. Default TRUE
#' @param mean Estimate mean coverage per region. Default TRUE. FALSE produces coverage per base per target.
#' @param hist Enables progress messages. Default False.
#' @param fai Indexed genome to which sequece has been aligned. Default none. Only require if there is an issue with chromosome naming.
#' @param verbose Enables progress messages. Default False.
#' @param output_dir Output directory path. Default none.
#' @export


bed_coverage=function(bin_path="tools/bedtools2/bin/bedtools",bam="",bed="",
verbose=FALSE,sorted=TRUE,mean=TRUE,fai="",suffix="",output_dir="",hist=FALSE){
    
    out_file_dir=set_dir(dir=output_dir,name="coverage")

    srt=""
    if (sorted){
      srt="-sorted"
    }

    if (fai!=""){
      fai=paste("-g",fai)
    }
    if (suffix!=""){
      suffix=paste0(".",suffix)
    }
    mode=""
    if (hist){
      mode="-hist"

      ## Filter to reduce the size of the output as it produces the coverage per base stats too

      out_file=paste0("| grep \"all\"",">",out_file_dir,get_file_name(bam),suffix,".Histogram_Coverage.txt")
    }else{
      if (mean){
        mode="-mean"
        out_file=paste0(">",out_file_dir,get_file_name(bam),suffix,".Per_Region_Coverage.txt")
      }else{
        mode="-d"
        out_file=paste0(">",out_file_dir,get_filename(bam),suffix,".Per_Base_Coverage.txt")
      }
    }

    exec_code=paste(bin_path,"coverage -a",bed, "-b" ,bam,fai, mode,srt,out_file)
    if(verbose){
        print(exec_code)
    }
        error=system(exec_code)
  if(error!=0){
    stop("bedtools failed to run due to unknown error.
    Check std error for more information.")
  }

}



#' Function to collect chromosome data in bam
#'
#' This function takes a BAM file and collects the chr names from the bam
#' file.
#'
#' @param bin_path Path to samtools executable. Default path tools/samtools/samtools.
#' @param bam Path to directory with BAM files to merge.
#' @param verbose Enables progress messages. Default False.
#' @export

get_bam_reference_chr=function(bin_path="tools/samtools/samtools",bam="",verbose=FALSE){
  options(scipen = 999)

  exec_code=paste0(bin_path," view -H ",bam," | grep @SQ")
  if(verbose){
    print(exec_code)
  }
  SQ=read.table(text=system(exec_code,intern=TRUE),stringsAsFactors=FALSE)
  chr=SQ[,2]
  chr=unlist(lapply(chr,FUN=function(x){strsplit(x,":")[[1]][2]}))

  size=SQ[,3]
  size=unlist(lapply(size,FUN=function(x){strsplit(x,":")[[1]][2]}))

  ref_chr=data.frame(chr=chr,start=0,end=as.numeric(size))
  return(ref_chr)
}



#' Function to generate seq with trailing ner
#'
#' This functions genereates a sequence of number
#'
#'
#' @param from Start of sequence
#' @param to End of sequence
#' @param by Step size
#' @export



seqlast <- function (from, to, by)
{
  vec <- do.call(what = seq, args = list(from, to, by))
  if ( tail(vec, 1) != to ) {
    return(c(vec, to))
  } else {
    return(vec)
  }
}



#' Split chromosomes into bins
#'
#' This functions takes a BED file of chromosome positions for BAM file input
#'
#'
#' @param bin_path Path to samtools executable. Default path tools/samtools/samtools.
#' @param bam Path to directory with BAM files to merge.
#' @param verbose Enables progress messages. Default False.
#' @param bin_size Bin size. Default 400000000 pb
#' @export


bin_chromosomes <- function(bin_path="tools/samtools/samtools",bam="",verbose=FALSE,
bin_size=40000000){
  options(scipen = 999)
  chr=get_bam_reference_chr(bin_path=bin_path,bam=bam,verbose=verbose)
  bed=chr%>% dplyr::group_by(chr) %>%
  dplyr::summarise(start=seqlast(start,end,bin_size)) %>%
  dplyr::mutate(end=dplyr::lead(start)) %>% tidyr::drop_na()
  bed=bed[stringr::str_order(paste0(bed$chr,"_",bed$start), numeric = TRUE),]
  return(bed)
}






#' Get current script path
#'
#' This functions obtains the current script path
#'
#' @export




thisFile <- function() {
        cmdArgs <- commandArgs(trailingOnly = FALSE)
        needle <- "--file="
        match <- grep(needle, cmdArgs)
        if (length(match) > 0) {
                # Rscript
                return(normalizePath(sub(needle, "", cmdArgs[match])))
        } else {
                # 'source'd via R console
                return(normalizePath(sys.frames()[[1]]$ofile))
        }
}
