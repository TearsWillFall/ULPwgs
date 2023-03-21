#' Read reference genome fasta file
#' 
#' Read reference genome fasta file
#' 
#' fasta object structure:
#'      time: System time at read time
#'      fasta_origin: PATH to fasta origin file
#'      body: Fasta body
#'       
#' 
#' @param fasta PATH to FASTA file
#' @param fai PATH to FAI file. If not given assumes same directory as fasta file
#' @param sort Sort genomic positions alphanumerically. Default TRUE
#' @param threads Number of threads to use. Default 1
#' @return A fasta_object
#' @export



read_pileup=function(
    pileup=NULL,
    header=TRUE,
    sep="\t",
    rename=TRUE,
    sort=TRUE,
    threads=1
){
  options(scipen=999)
 
  sort_pileup=function(pileup_body=NULL){
        return(pileup_body %>% dplyr::arrange(
            gtools::mixedorder(gid))
        )
  }

    pileup_origin=NULL
    if(is.null(pileup)){
        stop("pileup argument is of type NULL")
    }else if(is.data.frame(pileup)){
        body=pileup
        origin_file_type="data.frame"
    }else if(file.exists(pileup)){
        if(grepl(".pileup$",pileup)){
            body=data.table::fread(
                file=pileup,sep=sep,header=header,
                colClasses="character",nThread=threads
            )
            origin_file_type="pileup"
            
        }else{
            stop("Not valid file format. Valid file formats are PILEUP.")
        }
   
    }else{
        stop("Not recognized PILEUP format")
    }


    if(sort){
        body=sort_pileup(body)
    }

    pileup_object=list(
        time=Sys.time(),
        origin_file_type=origin_file_type,
        pileup_origin=pileup_origin,
        body=body
    )

    return(pileup_object)

}