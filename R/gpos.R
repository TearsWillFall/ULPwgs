#' Read file with Genomic Positions
#' 
#' Valid positional formats include data.frame with chrom and position columns,
#' as well as VCF and BED files.
#' 
#' Gpos structure:
#'      time: System time at read time
#'      origin_file_type: data.frame/vcf/bed
#'      gpos_origin: Path to gpos origin file
#'      body: Gpos body
#'       
#' 
#' @param pos Data.frame or Path to BED/VCF style files
#' @param sort Sort genomic positions alphanumerically. Default TRUE
#' @return A data.frame with genomic position information
#' @export


read_gpos=function(
    gpos=NULL,
    sort=TRUE
){  
    origin_file_type=NULL
    gpos_origin=NULL
    col_names=c("chrom","pos")


    bed_to_pos=function(
        bed_body=NULL
    ){
        chrom=bed_body[,1]
        pos=as.numeric(bed_body[,2])+1
        dat_pos=data.frame(chrom=chrom,pos=pos)
        return(dat_pos)
    }


    sort_gpos=function(gpos=NULL){
        return(gpos %>% dplyr::arrange(
            gtools::mixedorder(chrom),pos)
        )
    }


    ### Check if input is a file or data.frame
    if(is.null(gpos)){
        stop("Pos argument of type NULL")
    }else if (is.data.frame(gpos)){
        body=gpos
        gpos_origin="data.frame"
    }else if (file.exists(gpos)){
        ### Read VCF file and return body for region information
        if(grepl(".vcf",gpos)){
            origin_file_type="vcf"
            body=read_vcf(gpos)$body[,1:2]
        }else if(grepl(".bed",gpos)){
             origin_file_type="bed"
            ### Read BED file and transform into position based (VCF). 1 Based
            body=bed_to_pos(read_bed(gpos)$body)
        }else{
           stop("Not valid file format. Valid file formats are VCF/BED.") 
        }
        gpos_origin=normalizePath(gpos)
    }else{
        stop("Not recognized input Genomic Position format")
    }


    names(body)=col_names

    if(sort){
        body=sort_gpos(body)
    }

    gpos_object=list(
        time=Sys.time(),
        origin_file_type=origin_file_type,
        gpos_origin=gpos_origin,
        body=body
    )

    return(gpos_object)

}




#' Read file with Genomic Positions
#' 
#' Valid positional formats include data.frame with chrom and position columns,
#' as well as VCF and BED files.
#' 
#' Gpos structure:
#'      time: System time at read time
#'      origin_file_type: data.frame/vcf/bed
#'      gpos_origin: Path to gpos origin file
#'      body: Gpos body
#'       
#' 
#' @param env_cov Path to conda enviroment location for coverage
#' @param bam Path to BAM file
#' @param gpos Gpos structure or path to VCF/BED file with genomic location information
#' @return A data.frame with genomic position information
#' @export



get_coverage=function(
    env_cov=build_default_python_enviroment_list()$env_cov,
    bam=NULL,
    gpos=NULL
){
    
  run_main=function(
    .env
  ){
    .this.env=environment()
    append_env(to=.this.env,from=.env)
   
    set_main(.env=.this.env)
    
    gpos=read_gpos(gpos)

    fpath <- system.file("inst/python", paste0(fs,".py"), package=ns)

    .main$out_files$pileup=paste0(
      out_file_dir,"/",input_id,".pileup.txt"
    )

    tmp=write.table(x=gpos$body,
        file=paste0(tmp_dir,"/",input_id,".tmp"),
        sep="\t",col.names=FALSE,quote=FALSE
    )

    .main$exec_code=paste0(
        set_conda_envir(env_cov),
        fpath, 
        " -b ", bam, 
        " -g ",tmp,
        " -t ",threads,
        " -o ",.main$out_files,
        " -f True "
    )
  
    .main$exec_code=paste0(.main$exec_code," && rm",tmp) 
    
    run_job(.env=.this.env)

    .this.step=.main$steps[[fn]]
    
    .this.step$out_files$index_stats <- .main$out_file

    .env$.main <- .main

  }

  .base.env=environment()
  list2env(list(...),envir=.base.env)
  set_env_vars(
    .env=.base.env,
    vars="bam"
  )

  launch(.env=.base.env)

}
