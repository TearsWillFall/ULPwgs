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
    sort=FALSE,
    threads=1
){
  options(scipen=999)
  path=NA
  col_names=c(
    "chrom","chromStart","chromEnd","name",
    "score","start","thickStart","thickEnd",
    "itemRgb","blockCount","blockSizes","blockStarts")
   
    sort_bed=function(bed=NULL){
        return(
            bed %>% 
            dplyr::arrange(
              gtools::mixedorder(chrom),chromStart,chromEnd
            )
        )
    }

  if(is.null(bed)){
    stop("bed argument is of type NULL")
  }
  
  if(file.exists(bed)){
    path=bed

    bed=data.table::fread(
        file=bed,sep=sep,header=header,
        colClasses="character",
        nThread=threads
    )
    if(!header|rename){
        names(bed)=col_names[1:ncol(bed)]
    }

  }

  if(is.data.frame(bed)){
    body=bed
  }
  names(body)[1:3]<-c("chrom","chromStart","chromEnd")

  bed_object=list(
    time=Sys.time(),
    bed_origin=normalizePath(path),
    body=if(sort){sort_bed(body)}else{body}
  )

  return(bed_object)
}



