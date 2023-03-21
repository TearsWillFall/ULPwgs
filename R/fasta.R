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



read_fasta=function(
    fasta=NULL,
    fai=NULL,
    sort=TRUE,
    threads=1){
   
  
    fasta_origin=NULL
    ### Check if input is a file or data.frame
    if(is.null(fasta)){
        stop("fasta argument of type NULL")
    }else if (is.data.frame(fasta)){
        body=fasta
        origin_file_type="data.frame"
    }else if (file.exists(fasta)){
        ### Read VCF file and return body for region information
        if(grepl(".fa",fasta)){

              if(!is.null(fai)){
                fai=read_fai(fai) 
              }else{
                fai=read_fai(paste0(fasta,".fai"))
              }


            body=mclapply_os(X=seq(1,nrow(fai$body)),FUN=function(x){
                    info=fai$body[x,]
                    info$SEQ=system(paste0("sed '/^>",info$NAME," /,/^>/p'",fasta,"| head -n -1 | tail -n -2"),intern=TRUE)
                    return(info)
            },mc.cores=threads)

            body=dplyr::bind_rows(body)   
            origin_file_type="fasta"
            fasta_origin=normalizePath(fasta)

        }else{
           stop("Not valid file format. Valid file formats are FASTA") 
        }
    }else{
        stop("Not recognized FASTA format")
    }

    if(sort){
        body=body %>% arrange(body$OFFSET)
    }

    fasta_object=list(
        time=Sys.time(),
        origin_file_type=origin_file_type,
        fasta_origin=fasta_origin,
        body=body,
        fai=fai
    )
    return(fasta_object)
}

#' Read reference genome fasta index file
#' 
#' Read reference genome fasta index file
#' 
#' fai object structure:
#'      time: System time at read time
#'      fasta_origin: PATH to fasta origin file
#'      body: Fasta body
#'       
#'
#' @param fai PATH to FAI file.
#' @param sort Sort genomic positions alphanumerically. Default TRUE
#' @return A fasta_object
#' @export


read_fai=function(fai=NULL,sort=TRUE){
    body=read.csv(fai,header=FALSE,sep="\t")
    names(body)=c("NAME","LENGTH","OFFSET","LINEBASES","LINEWIDTH")
    
    if(sort){
        body=body %>% dplyr::arrange(OFFSET)
    }

    fai_object=list(
        time=Sys.time(),
        fai_origin=normalizePath(fai),
        body=body
    )
    return(fai_object)
}


#' Get reference base in FASTA at genomic position
#' 
#' Get reference base in FASTA at genomic position
#' 
#' gpos_object structure:
#'      time: System time at read time
#'      fasta_origin: PATH to fasta origin file
#'      body: Gpos body
#'       
#'
#' @param fasta PATH to FAI file.
#' @param fai PATH to FAI file.
#' @param gpos Genomic position.
#' @param sort Sort genomic positions alphanumerically. Default TRUE
#' @param threads Number of threads to use. Default 1
#' @return A fasta_object
#' @export

get_base_fasta=function(
      fasta=NULL,
      fai=NULL,
      gpos=NULL,
      sort=TRUE,
      header=TRUE,
      threads=1
){
    
      fasta=read_fasta(fasta=fasta,fai=fai,sort=sort,threads=threads)
      gpos=read_gpos(gpos=gpos,sort=sort,header=header)
      dat=mclapply_os(X=unique(gpos$body$chrom),FUN=function(x){
            tmp_gpos=gpos$body %>% dplyr::filter(chrom==x)
            tmp_fasta=fasta$body %>% dplyr::filter(NAME==x)
            this_fai=fasta$fai %>% dplyr::filter(NAME==x)
            this_chrom=lapply(1:nrow(tmp_gpos),FUN=function(y){
                this_gpos=tmp_gpos[y,]
                this_line=as.integer(this_gpos$pos/this_fai$LINEBASES)+1
                this_base=this_gpos$pos-(this_line*this_fai$LINEBASES)
                this_gpos$REF=substring(tmp_fasta$body[this_line,],this_base,this_base)
                return(this_gpos)
            })
            this_chrom=dplyr::bind_rows(this_chrom)
            this_chrom$OFFSET=this_fai$OFFSET
            return(this_chrom)
      },mc.cores=threads)

      dat=dplyr::bind_rows(dat) %>% dplyr::arrange(OFFSET,pos) %>% dplyr::select(-OFFSET)
      gpos_object=list(
            time=Sys.time(),
            fasta_origin=normalizePath(fasta),
            origin_file_type=gpos$origin_file_type,
            gpos_origin=gpos$gpos_origin,
            body=body
      )
      return(gpos_object)
}   

