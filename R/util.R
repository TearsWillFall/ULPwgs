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
