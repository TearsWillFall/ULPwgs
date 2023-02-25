#' Read a BED file
#' This function reads a BED file and stores it in a list format.
#' Assumes 0 based BED file
#' 
#' 
#' bed structure:
#'      time: System time at read time
#'      bed_origin: Path to BED origin file
#'      body: BED body
#'       
#' 
#' @param bed Path to the BED file
#' @param sep File separator 
#' @param header File contains header. Default TRUE
#' @param rename Rename BED column names to standard BED file format.
#' @return A list with time, body,bed_origin information
#' @export




read_bed=function(
    bed=NULL,sep="\t",
    header=TRUE,
    rename=TRUE,
    sort=TRUE
){
  options(scipen=999)
  col_names=c(
    "chrom","chromStart","chromEnd","name",
    "score","start","thickStart","thickEnd",
    "itemRgb","blockCount","blockSizes","blockStarts")
   
    sort_bed=function(bed=NULL){
        return(bed_body %>% dplyr::arrange(
            gtools::mixedorder(chrom),chromStart,chromEnd)
        )
    }

  if(is.null(bed)){
    stop("bed argument is of type NULL")
  }else if(is.data.frame(bed)){
    body=bed
  }else if(file.exists(bed)){
    if(grepl(".bed$",bed)){
        body=read.csv(
            file,sep=sep,header=header,
            colClasses="character",
            stringAsFactors=FALSE
        )
        if(!header|rename){
            names(body)=col_names[1:ncol(body)]
        }
    }else{
        stop("Not valid file format. Valid file formats are BED.")
    }
   
  }else{
    stop("Not recognized input BED format")
  }

  bed_object=list(
    time=Sys.time(),
    bed_origin=normalizePath(bed),
    body=ifelse(sort,sort_bed(body),body)
  )
  return(bed_object)
}







