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
    cns$weighted_log2=cns$width*cns$log2*cns$weight/cns$probes
    cns_target=cns %>% dplyr::filter(!grepl("Antitarget",gene))
    cns_antitarget=cns %>% dplyr::filter(grepl("Antitarget",gene))
    sol=data.frame(
        ploidy_all=(sum(cns$width)/sum(cns$weighted_log2))/2,
        ploidy_target=(sum(cns_target$width)/sum(cns_target$weighted_log2))/2,
        ploidy_antitarget=(sum(cns_antitarget$width)/sum(cns_antitarget$weighted_log2))/2
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
        20000,
        1000
    )
    cnr$bin_weight=cnr$width/cnr$bin_depth*cnr$weight
    cnr$weighted_log2=cnr$bin_weight*cnr$log2
    cnr_target=cnr %>% dplyr::filter(!grepl("Antitarget",gene))
    cnr_antitarget=cnr %>% dplyr::filter(grepl("Antitarget",gene))
    sol=data.frame(
        ploidy_all=(sum(cnr$bin_weight)/sum(cnr$weighted_log2))/2,
        ploidy_target=(sum(cnr$bin_weight)/sum(cnr_target$weighted_log2))/2,
        ploidy_antitarget=(sum(cnr$bin_weight)/sum(cnr_antitarget$weighted_log2))/2
    )
    return(sol)
}





