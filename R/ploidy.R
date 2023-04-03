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
    cns$length=as.numeric(cns$end)-as.numeric(cns$start)
    cns$weighted_log2=cns$length*cns$log2
    cns_target=cns %>% dplyr::filter(!grepl("Antitarget",gene))
    cns_antitarget=cns %>% dplyr::filter(grepl("Antitarget",gene))
    sol=data.frame(
        ploidy_all=(sum(cns$length)/sum(cns$weighted_log2))/2,
        ploidy_target=(sum(cns_target$length)/sum(cns_target$weighted_log2))/2,
        ploidy_antitarget=(sum(cns_antitarget$length)/sum(cns_antitarget$weighted_log2))/2
    )
    return(sol)
}


