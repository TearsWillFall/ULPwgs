#' Get the name of a sample from input file name
#'
#' This function takes the absolute/relative path to a file and
#' returns its base name without the file extension suffix.
#'
#' @param file_path Path to the input file
#' @return A string with the name of the file
#' @export


get_sample_name=function(file_path=""){
  sample_name=unlist(strsplit(basename(file_path),"\\."))[1]
  return(sample_name)
}


#' Get the extension of a file
#'
#' This function takes the absolute/relative path to a file and
#' returns the file extension suffix.
#'
#' @param file_path Path to the input file
#' @return A string with the extension of the file
#' @export
get_file_extension=function(file_path=""){
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

intersect_sample_name=function(file_path="",file_path2=""){
  tmp_name=get_sample_name(file_path2)
  sample_name=sapply(sapply(c(0:(nchar(tmp_name)-1)),function (i) substr(tmp_name,1,nchar(tmp_name)-i)),function (x) grepl(x,file_path))
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


index_ref=function(bin_path="tools/bwa/bwa",file="",verbose=FALSE){

    if(verbose){
        print(paste(bin_path,"index", file) )
    }
    system(print(paste(bin_path,"index", file) ) )


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
#' @param bam Path to the input bam file.
#' @param bin_path Path to bwa executable. Default tools/bedtools2/bin/bedtools.
#' @param bed Path to the input bed file.
#' @param sorted Are the input files sorted. Default TRUE
#' @param mean Estimate mean coverage per region. Default TRUE. FALSE produces coverage per base per target.
#' @param hist Enables progress messages. Default False.
#' @param verbose Enables progress messages. Default False.
#' @param output_dir Output directory path. Default none.
#' @export


bed_coverage=function(bin_path="tools/bedtools2/bin/bedtools",bam="",bed="",verbose=FALSE,sorted=TRUE,mean=TRUE,suffix="",output_dir="",hist=FALSE){
    sep="/"

    if(output_dir==""){
      sep=""
    }else{
      if (!dir.exists(output_dir)){
          dir.create(output_dir)
      }

    }

    sample_name=get_sample_name(bam)

    srt=""
    if (sorted){
      srt="-sorted"
    }


    if (suffix!=""){
      suffix=paste0(".",suffix)
    }
    mode=""
    if (hist){
      mode="-hist"

      ## Filter to reduce the size of the output as it produces the coverage per base stats too

      out_file=paste0("| grep \"all\"",">",output_dir,"/",sample_name,suffix,".Histogram_Coverage.txt")
    }else{
      if (mean){
        mode="-mean"
        out_file=paste0(">",output_dir,"/",sample_name,suffix,".Per_Region_Coverage.txt")
      }else{
        mode="-d"
        out_file=paste0(">",output_dir,"/",sample_name,suffix,".Per_Base_Coverage.txt")
      }
    }




    if(verbose){
        print(paste(bin_path,"coverage -a", "-b" ,bed, mode,srt,out_file))
    }
    system(paste(bin_path,"coverage  -a", "-b" ,bed, mode,srt,out_file))

    }
