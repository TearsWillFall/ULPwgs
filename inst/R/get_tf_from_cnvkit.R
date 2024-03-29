#!/usr/bin/env Rscript
options(tidyverse.quiet = TRUE)
options(dplyr.summarise.inform = FALSE)
library("optparse")
library("tidyverse")

option_list = list(
  make_option(c("-r", "--cnr"), type="character", default=NULL, 
              help="CNR data", metavar="character"),
  make_option(c("-s", "--cns"), type="character", default=NULL, 
              help="CNS data", metavar="character"),
  make_option(c("-t", "--tf"), type="character", default=NULL, 
              help="TF data", metavar="character"),
  make_option(c("-o", "--output_file"), type="character", default=NULL, 
              help="Output file", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);



import_cnvkit_data=function(
        cnr=NULL,
        cns=NULL,
        threads=1,
        chromosomes=c(1:22,"X","Y")
  ){
    jobs=length(cnr)
    dat=parallel::mclapply(1:jobs,FUN=function(x){
      lst=list(cnr=NA,cns=NA)
      cnr=read.table(cnr[[x]],header=TRUE,sep="\t")
      cnr=cnr[cnr$chromosome %in% chromosomes,]
      cnr$rid=paste0(cnr$chromosome,":",cnr$start,"-",cnr$end)
      cnr=GenomicRanges::makeGRangesFromDataFrame(cnr,keep.extra.columns=TRUE)
      lst$cnr=cnr
      cns=read.table(cns[[x]],header=TRUE,sep="\t")
      cns=cns[cns$chromosome %in% chromosomes,]
      cns$sid=paste0(cns$chromosome,":",cns$start,"-",cns$end)
      cns$size=cns$end-cns$start
      cns=GenomicRanges::makeGRangesFromDataFrame(cns,keep.extra.columns=TRUE)
      lst$cns=cns
      return(lst)
    },mc.cores=threads)
    names(dat)=Vectorize(ULPwgs::get_file_name)(cnr)
    return(dat)
  }

import_tf_data=function(
    tf=NULL,
    threads=1,
    chromosomes=c(1:22,"X","Y")
  ){
    jobs=length(tf)
    dat=parallel::mclapply(1:jobs,FUN=function(x){
      tfbs=ULPwgs::read_gpos(tf[[x]],sep=" ",rename=FALSE,sort=FALSE)$body
      tfbs=tfbs[tfbs$chrom %in% chromosomes,]
      tfbs=GenomicRanges::makeGRangesFromDataFrame(tfbs,start.field="pos",end.field="pos",keep.extra.columns=TRUE)
      return(tfbs)
    },mc.cores=threads)
    names(dat)=Vectorize(ULPwgs::get_file_name)(tf)
    return(dat)
}


get_tf_from_cnvkit=function(
    cnr=NULL,
    cns=NULL,
    tfbs=NULL,
    output_name=NULL,
    tf=NULL,
    id=NULL
){

    depth=median(cns$depth)
    overlaps=GenomicRanges::findOverlaps(tfbs,cnr)
    missing_tfbs=tfbs[-as.data.frame(overlaps)$queryHits]
    chromosomes=as.character(unique(GenomeInfoDb::seqnames(tfbs)))
    hit_tfbs=lapply(chromosomes,FUN=function(chr){
        tfbs_tmp=tfbs[GenomeInfoDb::seqnames(tfbs)==chr]
        cnr_tmp=cnr[GenomeInfoDb::seqnames(cnr)==chr]
        cns_tmp=cns[GenomeInfoDb::seqnames(cns)==chr]

        tmp_ranges=as.data.frame(IRanges::ranges(cnr_tmp))
        cnr_tmp$cnr_tfbs_pos=(tmp_ranges$start+tmp_ranges$end)/2
        cnr_tmp$cnr_tfbs_bin_size=tmp_ranges$width
   
        invisible(lapply(5:1,FUN=function(x){
            sol=dplyr::lag(cnr_tmp$rid,n=x)
            S4Vectors::mcols(cnr_tmp)[paste0("cnr_rid_left_",x)]<<-ifelse(
            is.na(sol),cnr_tmp$rid,sol)
        }))
        cnr_tmp$cnr_rid_central=cnr_tmp$rid
        
        invisible(lapply(1:5,FUN=function(x){
          sol=dplyr::lead(cnr_tmp$rid,n=x)
          S4Vectors::mcols(cnr_tmp)[paste0("cnr_rid_right_",x)]<<-ifelse(
             is.na(sol),cnr_tmp$rid,sol)
        }))

      
        invisible(lapply(5:1,FUN=function(x){
          sol=dplyr::lag(cnr_tmp$log2,n=x)
          S4Vectors::mcols(cnr_tmp)[paste0("cnr_tfbs_left_log2_",x)]<<-ifelse(
            is.na(sol),cnr_tmp$log2,sol)
        }))
        cnr_tmp$cnr_tfbs_central_log2=cnr_tmp$log2

        invisible(lapply(1:5,FUN=function(x){
          sol=dplyr::lead(cnr_tmp$log2,n=x)
          S4Vectors::mcols(cnr_tmp)[,paste0("cnr_tfbs_right_log2_",x)]<<-ifelse(
            is.na(sol),cnr_tmp$log2,sol
          )
        }))

        invisible(lapply(5:1,FUN=function(x){
            sol=dplyr::lag(cnr_tmp$depth,n=x)
            S4Vectors::mcols(cnr_tmp)[,paste0("cnr_tfbs_left_depth_",x)]<<-ifelse(
            is.na(sol),cnr_tmp$log2,sol
          )
        }))

        cnr_tmp$cnr_tfbs_central_depth=cnr_tmp$depth

        invisible(lapply(1:5,FUN=function(x){
          sol=dplyr::lead(cnr_tmp$depth,n=x)
          S4Vectors::mcols(cnr_tmp)[,paste0("cnr_tfbs_right_depth_",x)]<<-ifelse(
            is.na(sol),cnr_tmp$depth,sol
          )
         }))

        cnr_tmp=cnr_tmp[,grepl("cnr",names(S4Vectors::mcols(cnr_tmp)))]


        hit_tfbs=plyranges::join_overlap_left(
          tfbs_tmp,cnr_tmp,suffix=NULL)
        hit_tfbs$cnr_tfbs_pos_from_bin=as.data.frame(IRanges::ranges(hit_tfbs))$start-hit_tfbs$cnr_tfbs_pos

        tmp_ranges=as.data.frame(IRanges::ranges(cns_tmp))
        cns_tmp$cns_seg_pos=(tmp_ranges$start+tmp_ranges$end)/2
        cns_tmp$cns_seg_size=tmp_ranges$width
        cns_tmp$cns_seg_left_log2=dplyr::lag(cns_tmp$log2)
        cns_tmp$cns_seg_central_log2=cns_tmp$log2
        cns_tmp$cns_seg_right_log2=dplyr::lead(cns_tmp$log2)
        cns_tmp$cns_sid_left=dplyr::lag(cns_tmp$sid)
        cns_tmp$cns_sid_central=cns_tmp$sid
        cns_tmp$cns_sid_right=dplyr::lead(cns_tmp$sid)

        cns_tmp=cns_tmp[,grepl("cns",names(S4Vectors::mcols(cns_tmp)))]
        hit_tfbs=plyranges::join_overlap_left(
          hit_tfbs,cns_tmp,suffix=NULL
        )
        hit_tfbs$cnr_seg_pos_from_bin=as.data.frame(IRanges::ranges(hit_tfbs))$start-hit_tfbs$cns_seg_pos

          return(hit_tfbs)
        }
    )

    hit_tfbs=plyranges::bind_ranges(hit_tfbs)
    hit_tfbs$sample_depth=depth
    hit_tfbs$tf=tf
    hit_tfbs$id=id
    hit_tfbs=plyranges::mutate(
      plyranges::group_by(hit_tfbs,gid),
      BP=ifelse(plyranges::n()>1,1,0)
    )
    scores_tfbs= as.data.frame(hit_tfbs) %>% 
    dplyr::group_by(id,tf) %>%
    dplyr::summarise(
      dplyr::across(
        c(
          dplyr::starts_with("cnr_tfbs_"),
          dplyr::starts_with("cns_seg")),
            list(
                mean=~mean(.,na.rm=TRUE),
                median=~median(.,na.rm=TRUE),
                sd=~sd(.,na.rm=TRUE),
                iqr=~IQR(.,na.rm=TRUE),
                n=~sum(!is.na(.))
        ),
        .names = "{.col}.{.fn}"
      )
    )
    

    scores_tfbs=tidyr::pivot_longer(scores_tfbs,!id:tf)
    scores_tfbs=scores_tfbs %>% tidyr::separate_wider_delim(name,delim=".",names=c("var","type"))
    scores_tfbs=scores_tfbs %>% tidyr::pivot_wider(names_from=type,values_from=value)
    data.table::fwrite(as.data.frame(hit_tfbs),file=paste0(output_name,".hits.txt"),nThread=1)
    data.table::fwrite(as.data.frame(missing_tfbs),file=paste0(output_name,".miss.txt"),nThread=1)
    data.table::fwrite(scores_tfbs,file=paste0(output_name,".scores.txt"),nThread=1)
    return()
}

cnvkit_data=import_cnvkit_data(cnr=opt$cnr,cns=opt$cns)
tf_data=import_tf_data(tf=opt$tf)
output_name=opt$output_file
if(is.null(opt$output_file)){
    output_name=paste0(names(cnvkit_data),".",names(tf_data),".txt")
}

invisible(get_tf_from_cnvkit(
    cnr=cnvkit_data[[1]]$cnr,
    cns=cnvkit_data[[1]]$cns,
    tfbs=tf_data[[1]],
    output_name=output_name,
    tf=names(tf_data),
    id=names(cnvkit_data)
))
