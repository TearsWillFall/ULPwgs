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
#' @param gpos Data.frame or Path to BED/VCF style files
#' @param sort Sort genomic positions alphanumerically. Default TRUE
#' @return A data.frame with genomic position information
#' @export


read_gpos=function(
    gpos=NULL,
    sort=TRUE,
    header=TRUE,
    sep="\t",
    rename=TRUE,
    threads=1
){  
    origin_file_type=NULL
    gpos_origin=NULL

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
        if(!all(c("chrom","pos") %in% names(body))){
            stop("Columns chrom and/or pos missing in data.frame")
        }
        body[,c("ref","alt","gt","gid")[!c("ref","alt","gt","gid") %in% names(body)]]="."
        body=body %>% dplyr::select(chrom,pos,ref,alt,gt,gid)
        origin_file_type="data.frame"
    }else if (file.exists(gpos)){
        ### Read VCF file and return body for region information
        if(grepl(".vcf",gpos)){
            origin_file_type="vcf"
            tmp=read_vcf(gpos,sep=sep,threads=threads)
            if(!is.null(tmp$descriptors$FORMAT$GT)){
                body=tmp$body %>% unnest_vcf_body() %>% 
                dplyr::filter(FORMAT=="GT") %>% 
                dplyr::select(CHROM,POS,REF,ALT,VALUE)
            }else{
                body=tmp$body %>% dplyr::select(CHROM,POS,REF,ALT)
                body[,"gt"]="."
                body[,"gid"]=paste0(body$CHROM,":",body$POS,"_",body$REF,"-",body$ALT)
            }
        }else if(grepl(".bed",gpos)){
            origin_file_type="bed"
            body=bed_to_pos(
                read_bed(gpos,
                header=header,
                sep=sep,
                rename=rename,
                sort=sort)$body
            )
            body[,c("ref","alt","gt")]="."
            body[,"gid"]=paste0(body$chrom,":",body$pos)
        }else if (grepl(".pileup",gpos)){
            origin_file_type="pileup"
            body=read_pileup(pileup=gpos,header=header,sep=sep,rename=rename,sort=FALSE)$body[,c("chrom","pos","ref","alt","gt","gid")]
        }else if(grepl(".gpos",gpos)){
            body=data.table::fread(file=gpos,header=header,sep=sep,nThread=threads)[,c("chrom","pos")]
            body[,c("ref","alt","gt")]="."
            body[,"gid"]=paste0(body$chrom,":",body$pos)
        }else{
           stop("Not valid file format. Valid file formats are VCF/BED.") 
        }
        gpos_origin=normalizePath(gpos)
    }else{
        stop("Not recognized input Genomic Position format")
    }


    names(body)=c("chrom","pos","ref","alt","gt","gid")

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
    gpos=NULL,
    gt="het",
    write=TRUE,
    ...
){
    
  run_main=function(
    .env
  ){
    .this.env=environment()
    append_env(to=.this.env,from=.env)
   
    set_main(.env=.this.env)
    
    gpos=read_gpos(gpos,threads=threads) %>% only_snps_gpos()
    if(gt=="het"){
        gpos=gpos %>% only_het_gpos()
    }else if(gt=="hom")(
        gpos=gpos %>% only_hom_gpos()
    )

    fpath <- system.file("python", paste0(fn,".py"), package=ns)
    system(paste0("chmod 755 ",fpath))
    fpath<-paste0("python ",fpath)

    .main$out_files$pileup=paste0(
      out_file_dir,"/",input_id,".pileup"
    )
    tmp=paste0(tmp_dir,"/",input_id,".tmp")
    write.table(x=gpos$body,
        file=tmp,
        sep="\t",
        col.names=FALSE,
        row.names=FALSE,
        quote=FALSE
    )
    
    add=""
    if(write){
        add=paste0(" -o ",.main$out_files)
    }

    .main$exec_code=paste0(
        set_conda_envir(env_cov),
        fpath, 
        " -b ", bam, 
        " -g ",tmp,
        " -t ",threads,
        add,
        " -i ",input_id,
        " -f True "
    )
  
    .main$exec_code=paste0(.main$exec_code," && rm ",tmp) 
    
    run_job(.env=.this.env)

    .this.step=.main$steps[[fn_id]]
    
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



#' @export

only_het_gpos=function(gpos=NULL){
    gpos$body=gpos$body %>% dplyr::filter(grepl("0\\|1|1\\|0",gt))
    return(gpos)
}

#' @export

only_hom_gpos=function(gpos=NULL){
    gpos$body=gpos$body %>% dplyr::filter(grepl("1\\|1|0\\|0",gt))
    return(gpos)
}

#' @export

only_snps_gpos=function(gpos=NULL){
    gpos$body=gpos$body %>% dplyr::rowwise() %>% 
    dplyr::filter(nchar(ref)==1&nchar(alt)==1)
    return(gpos)
}