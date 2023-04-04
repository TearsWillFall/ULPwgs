#' Wrapper around assuite
#' 
#'
#' Wrapper around AmpliconArchitectSuite
#' 
#' @param cns CNS file from cnvKit in data.frame format
#' @param chrom List of chromosomes to calculate ploidy on
#' @export
#' 

ploidy_from_cns=function(cns=NULL,chrom=c(1:22)){
    cns=cns[cns$chromosome %in% chrom,]
    cns$width=as.numeric(cns$end)-as.numeric(cns$start)
    cns$seg_weight=cns$width*cns$weight/cns$probes
    cns$weighted_log2=cns$seg_weight*cns$log2
    sol=data.frame(
        ploidy_all=(sum(cns$weighted_log2)/sum(cns$width))/2
    )
    return(sol)
}




#' Wrapper around assuite
#' 
#'
#' Wrapper around AmpliconArchitectSuite
#' 
#' @param cnr CNS file from cnvKit in data.frame format
#' @param chrom List of chromosomes to calculate ploidy on
#' @export
#' 



ploidy_from_cnr=function(cnr=NULL,chrom=c(1:22)){
    cnr=cnr[cnr$chromosome %in% chrom,]
    cnr$width=as.numeric(cnr$end)-as.numeric(cnr$start)
    cnr$bin_type=ifelse(cnr$gene=="Antitarget","Antitarget","Target")
    cnr$bin_depth=ifelse(
        cnr$bin_type=="Antitarget",
        sum(cnr[cnr$bin_type=="Target",]$depth)/
        sum(cnr[cnr$bin_type=="Antitarget",]$depth),
        1
    )

    cnr$bin_weight=cnr$width/cnr$bin_depth*cnr$weight
    cnr$weighted_log2=cnr$bin_weight*cnr$log2
    cnr_target=cnr %>% dplyr::filter(!grepl("Antitarget",gene))
    cnr_antitarget=cnr %>% dplyr::filter(grepl("Antitarget",gene))
    sol=data.frame(
        ploidy_all=(sum(cnr$weighted_log2)/sum(cnr$bin_weight))/2,
        ploidy_target=(sum(cnr_target$weighted_log2)/sum(cnr$bin_weight))/2,
        ploidy_antitarget=(sum(cnr_antitarget$weighted_log2)/sum(cnr$bin_weight))/2
    )
    return(sol)
}





