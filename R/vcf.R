#' Read a VCF file
#' This function read a VCF file and stores it in a list format
#' 
#' 
#'
#' @param vcf Path to the VCF file
#' @return A list with the header, body and column names of the VCF
#' @export


read_vcf=function(vcf="",sep="\t"){
  cols=c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")
  if(check_if_compressed(vcf)){
      read=system(paste0("gunzip -c ",vcf),intern=TRUE)
  }else{
      read=readLines(vcf)
  }
  header=read[grepl("^#",read)]
  body=read[!grepl("^#",read)]
  col_names=header[length(header)]
  raw_header=header[-length(header)]
  descriptors=extract_descriptors_vcf(raw_header)
  body=read.delim(text=body,stringsAsFactors = FALSE,sep=sep,
  header=FALSE)
  names(body)=read.table(text=sub("#","",col_names),stringsAsFactors = FALSE)
  samples=setdiff(names(body),cols)
  body=extract_body_vcf(body,samples)
  vcf_object=list(
    vcf_origin=normalizePath(vcf),
    samples=samples,
    descriptors=descriptors,
    body=body
  )
  return(vcf_object)
}

#' Extract VCF body
#' This function VCF body columns in a easy to access format
#' 
#'
#' @param vcf_body VCF body structure
#' @param samples Sample columns in VCF file
#' @export

extract_body_vcf=function(vcf_body,vcf_samples){
  extract_col_vcf=function(col,sep=";"){
    tmp=Vectorize(strsplit,USE.NAMES=FALSE)(col,split=sep)
    return(tmp)
  }
  vcf_body$FILTER=extract_col_vcf(vcf_body$FILTER)
  extract_info_vcf=function(col){
     info=extract_col_vcf(col,sep=";")
     extract_key_value=function(row){
      split=unlist(stringi::stri_split_fixed(row,pattern="=",n=2))
      hash=list(values=ifelse(is.na(split[2]),"",split[2]))
      names(hash)=split[1]
      return(hash)
     }
    info=lapply(info,FUN=function(row){
        tmp=unlist(lapply(row,FUN=extract_key_value),recursive=FALSE)
       
     })
    return(info)
  }

  vcf_body$INFO=extract_info_vcf(vcf_body$INFO)
  vcf_body$FORMAT=extract_col_vcf(vcf_body$FORMAT,sep=":")


  lapply(vcf_samples,FUN=function(sample){
     vcf_body[[sample]]<<-as.list(extract_col_vcf(vcf_body[[sample]],sep=":"))
    return()
  })

  vcf_body=vcf_body %>% tidyr::pivot_longer(names_to="SAMPLE",
  cols=vcf_samples,values_to="VALUE") %>% 
  tidyr::nest(SAMPLE=SAMPLE,FORMAT=FORMAT,VALUE=VALUE)

  return(vcf_body)
}


#' Add SNV AF information to Strelka VCF file
#' This function reads and modifies Strelka produced VCF files to 
#' add an entry for AF within the VCF file
#' 
#'
#' @param bin_bgzip Path to bgzip executable.
#' @param bin_tabix Path to TABIX executable.
#' @param vcf Path to VCF file
#' @param compress Compress VCF file. Default TRUE.
#' @param index Index VCF file. Default TRUE.
#' @param index_format VCF index format. Default tbi. Options [tbi,cbi].
#' @param bgzip_index Create BGZIP index for compressed file. Default FALSE
#' @param output_dir Path to the output directory.
#' @param clean Remove input VCF after completion. Default FALSE.
#' @param verbose Enables progress messages. Default False.
#' @param executor_id Task EXECUTOR ID. Default "gatherBQSR"
#' @param task_name Task name. Default "gatherBQSR"
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param threads Number of threads to split the work. Default 3
#' @param ram [OPTIONAL] If batch mode. RAM memory in GB per job. Default 1
#' @param update_time [OPTIONAL] If batch mode. Show job updates every update time. Default 60
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export


add_snv_af_strelka_vcf=function(
  bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
  bin_tabix=build_default_tool_binary_list()$bin_tabix,
  vcf,overwrite=FALSE,
  compress=TRUE,index=TRUE,
  index_format="tbi",bgzip_index=FALSE,
  clean=TRUE,verbose=FALSE, 
  executor_id=make_unique_id("addAFsnvsStrelka"),
  task_name="addAFsnvsStrelka",
  batch_config=build_default_preprocess_config(),
  mode="local",time="48:0:0",
  threads=1,ram=4,update_time=60,
  wait=FALSE,hold=NULL){
    
    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file=paste0(get_file_name(vcf),".snvs")
    out_file_dir=dirname(vcf)
    job=build_job(executor_id=executor_id,task_id=task_id)

    job_report=build_job_report(
      job_id=job,
      executor_id=executor_id,
      exec_code=list(),
      task_id=task_id,
      input_args = argg,
      out_file_dir=out_file_dir,
      out_files=list()
    )

    vcf_dat=read_vcf(vcf)
    vcf_dat$body=vcf_dat$body %>% unnest_vcf_body()
    vcf_dat$body=vcf_dat$body %>% dplyr::group_by_at(dplyr::vars(-VALUE,-FORMAT)) %>% 
    dplyr::group_modify(~dplyr::add_row(.x,FORMAT="AF"))
    ##Extract Tier1 read information for REF and ALT
    vcf_dat$body=vcf_dat$body %>% dplyr::mutate(
      UREF=strsplit(VALUE[FORMAT==paste0(REF,"U")],split=",")[[1]][1],
      UALT=strsplit(VALUE[FORMAT==paste0(ALT,"U")],split=",")[[1]][1]) %>%
      dplyr::mutate(VALUE=ifelse(FORMAT=="AF",
      as.numeric(UALT)/(as.numeric(UREF)+as.numeric(UALT)),VALUE))%>% dplyr::select(-c(UALT,UREF))
    vcf_dat$body=vcf_dat$body %>% 
    dplyr::mutate(VALUE=ifelse(is.na(VALUE),"",VALUE)) %>% nest_vcf_body()
   

    
    add_af_descriptor<-function(){
        list(Number="1",Type="Float",Description="\"Variant allelic frequency for tier 1 reads\"")
    }

    vcf_dat$descriptors$FORMAT[["AF"]]<-add_af_descriptor()
    
    if(!overwrite){
      out_file=paste0(out_file,".af")
    }

    job_report[["steps"]][["writeVCF"]]<-write_vcf(
      vcf=vcf_dat,
      output_name=out_file,
      output_dir=out_file_dir,
      bin_bgzip=bin_bgzip,
      bin_tabix=bin_tabix,
      compress=compress,
      index=index,index_format=index_format,
      bgzip_index=bgzip_index,
      clean=clean,verbose=verbose,mode=mode,
      batch_config=batch_config,
      time=time,threads=threads,
      ram=ram,hold=hold
    )
  
    return(job_report)
}



#' Add indel AF information to Strelka VCF file
#' This function reads and modifies Strelka produced VCF files to 
#' add an entry for AF within the VCF file
#' 
#' @param bin_bgzip Path to bgzip executable.
#' @param bin_tabix Path to TABIX executable.
#' @param vcf Path to VCF file
#' @param compress Compress VCF file. Default TRUE.
#' @param index Index VCF file. Default TRUE.
#' @param index_format VCF index format. Default tbi. Options [tbi,cbi].
#' @param bgzip_index Create BGZIP index for compressed file. Default FALSE
#' @param output_dir Path to the output directory.
#' @param clean Remove input VCF after completion. Default FALSE.
#' @param verbose Enables progress messages. Default False.
#' @param executor_id Task EXECUTOR ID. Default "gatherBQSR"
#' @param task_name Task name. Default "gatherBQSR"
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param threads Number of threads to split the work. Default 3
#' @param ram [OPTIONAL] If batch mode. RAM memory in GB per job. Default 1
#' @param update_time [OPTIONAL] If batch mode. Show job updates every update time. Default 60
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export


add_indel_af_strelka_vcf=function(
  bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
  bin_tabix=build_default_tool_binary_list()$bin_tabix,
  vcf,overwrite=FALSE,
  compress=TRUE,index=TRUE,
  index_format="tbi",bgzip_index=FALSE,
  clean=TRUE,verbose=FALSE, 
  executor_id=make_unique_id("addAFindelStrelka"),
  task_name="addAFindelStrelka",
  batch_config=build_default_preprocess_config(),
  mode="local",time="48:0:0",
  threads=1,ram=4,update_time=60,
  wait=FALSE,hold=NULL){
    
    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file=paste0(get_file_name(vcf),".indel")
    out_file_dir=dirname(vcf)
    job=build_job(executor_id=executor_id,task_id=task_id)


    job_report=build_job_report(
      job_id=job,
      executor_id=executor_id,
      exec_code=list(),
      task_id=task_id,
      input_args = argg,
      out_file_dir=out_file_dir,
      out_files=list()
    )


    vcf_dat=read_vcf(vcf)
    vcf_dat$body=vcf_dat$body %>% unnest_vcf_body()
    vcf_dat$body=vcf_dat$body %>% dplyr::group_by_at(dplyr::vars(-VALUE,-FORMAT)) %>% 
    dplyr::group_modify(~dplyr::add_row(.x,FORMAT="AF"))
    ##Extract Tier1 read information for REF and ALT
    vcf_dat$body=vcf_dat$body %>% dplyr::mutate(
      UREF=strsplit(VALUE[FORMAT=="TAR"],split=",")[[1]][1],
      UALT=strsplit(VALUE[FORMAT=="TIR"],split=",")[[1]][1]) %>%
      dplyr::mutate(VALUE=ifelse(FORMAT=="AF",
      as.numeric(UALT)/(as.numeric(UREF)+as.numeric(UALT)),VALUE))%>% dplyr::select(-c(UALT,UREF))
    vcf_dat$body=vcf_dat$body %>% 
    dplyr::mutate(VALUE=ifelse(is.na(VALUE),"",VALUE)) %>% nest_vcf_body()

    
    add_af_descriptor<-function(){
        list(Number="1",Type="Float",Description="\"Variant allelic frequency for tier 1 reads\"")
    }

    vcf_dat$descriptors$FORMAT[["AF"]]<-add_af_descriptor()
    
    if(!overwrite){
      out_file=paste0(out_file,".af")
    }

    job_report[["steps"]][["writeVCF"]]<-write_vcf(
      vcf=vcf_dat,
      output_name=out_file,
      output_dir=out_file_dir,
      bin_bgzip=bin_bgzip,
      bin_tabix=bin_tabix,
      compress=compress,
      index=index,index_format=index_format,
      bgzip_index=bgzip_index,
      clean=clean,verbose=verbose,mode=mode,
      batch_config=batch_config,
      time=time,threads=threads,
      ram=ram,hold=hold
    )
  
    return(job_report)
}




#' Add structural variant AF information to Strelka VCF file
#' This function reads and modifies Strelka produced VCF files to 
#' add an entry for AF within the VCF file
#' 
#' @param bin_bgzip Path to bgzip executable.
#' @param bin_tabix Path to TABIX executable.
#' @param vcf Path to VCF file
#' @param compress Compress VCF file. Default TRUE.
#' @param index Index VCF file. Default TRUE.
#' @param index_format VCF index format. Default tbi. Options [tbi,cbi].
#' @param bgzip_index Create BGZIP index for compressed file. Default FALSE
#' @param output_dir Path to the output directory.
#' @param clean Remove input VCF after completion. Default FALSE.
#' @param verbose Enables progress messages. Default False.
#' @param executor_id Task EXECUTOR ID. Default "gatherBQSR"
#' @param task_name Task name. Default "gatherBQSR"
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param threads Number of threads to split the work. Default 3
#' @param ram [OPTIONAL] If batch mode. RAM memory in GB per job. Default 1
#' @param update_time [OPTIONAL] If batch mode. Show job updates every update time. Default 60
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export



add_sv_af_strelka_vcf=function(
  bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
  bin_tabix=build_default_tool_binary_list()$bin_tabix,
  vcf,overwrite=FALSE,
  compress=TRUE,index=TRUE,
  index_format="tbi",bgzip_index=FALSE,
  clean=TRUE,verbose=FALSE, 
  executor_id=make_unique_id("addAFindelStrelka"),
  task_name="addAFindelStrelka",
  batch_config=build_default_preprocess_config(),
  mode="local",time="48:0:0",
  threads=1,ram=4,update_time=60,
  wait=FALSE,hold=NULL){
    
    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file=paste0(get_file_name(vcf),".indel")
    out_file_dir=dirname(vcf)
    job=build_job(executor_id=executor_id,task_id=task_id)


    job_report=build_job_report(
      job_id=job,
      executor_id=executor_id,
      exec_code=list(),
      task_id=task_id,
      input_args = argg,
      out_file_dir=out_file_dir,
      out_files=list()
    )


    vcf_dat=read_vcf(vcf)
    vcf_dat$body=vcf_dat$body %>% unnest_vcf_body()
    vcf_dat$body=vcf_dat$body %>% dplyr::group_by_at(dplyr::vars(-VALUE,-FORMAT)) %>% 
    dplyr::group_modify(~dplyr::add_row(.x,FORMAT="AFP")) %>% 
    dplyr::group_modify(~dplyr::add_row(.x,FORMAT="AFS"))
    ##Extract Tier1 read information for REF and ALT
    vcf_dat$body=vcf_dat$body %>% dplyr::mutate(
      PREF=tryCatch({strsplit(VALUE[FORMAT=="PR"],split=",")[[1]][1]},error=function(e){NA}),
      PALT=tryCatch({strsplit(VALUE[FORMAT=="PR"],split=",")[[1]][2]},error=function(e){NA}),
      SREF=tryCatch({strsplit(VALUE[FORMAT=="SR"],split=",")[[1]][1]},error=function(e){NA}),
      SALT=tryCatch({strsplit(VALUE[FORMAT=="SR"],split=",")[[1]][2]},error=function(e){NA})) %>%
      dplyr::mutate(VALUE=ifelse(FORMAT=="AFP",
      as.numeric(PALT)/(as.numeric(PREF)+as.numeric(PALT)),VALUE))%>% 
      dplyr::mutate(VALUE=ifelse(FORMAT=="AFS",
      as.numeric(SALT)/(as.numeric(SREF)+as.numeric(SALT)),VALUE)) %>% 
      dplyr::select(-c(PALT,PREF,SALT,SREF))
    vcf_dat$body=vcf_dat$body %>% 
    dplyr::mutate(VALUE=ifelse(is.na(VALUE),"",VALUE)) %>% nest_vcf_body()

    
    add_afp_descriptor<-function(){
        list(Number="1",Type="Float",Description="\"Variant allelic frequency by paired reads spanning the region.\"")
    }

    add_afs_descriptor<-function(){
        list(Number="1",Type="Float",Description="\"Variant allelic frequency by split-reads spanning the region.\"")
    }

    vcf_dat$descriptors$FORMAT[["AFP"]]<-add_afp_descriptor()
    vcf_dat$descriptors$FORMAT[["AFS"]]<-add_afs_descriptor()
    
    if(!overwrite){
      out_file=paste0(out_file,".af")
    }

    job_report[["steps"]][["writeVCF"]]<-write_vcf(
      vcf=vcf_dat,
      output_name=out_file,
      output_dir=out_file_dir,
      bin_bgzip=bin_bgzip,
      bin_tabix=bin_tabix,
      compress=compress,
      index=index,index_format=index_format,
      bgzip_index=bgzip_index,
      clean=clean,verbose=verbose,mode=mode,
      batch_config=batch_config,
      time=time,threads=threads,
      ram=ram,hold=hold
    )
  
    return(job_report)
}












#' Add AF information to Strelka VCF file
#' This function reads and modifies Strelka produced VCF files to 
#' add an entry for AF within the VCF file
#' 
#' @param bin_bgzip Path to bgzip executable.
#' @param bin_tabix Path to TABIX executable.
#' @param vcf_snv Path to VCF file with snv
#' @param vcf_indel Path to VCF file with snv
#' @param compress Compress VCF file. Default TRUE.
#' @param index Index VCF file. Default TRUE.
#' @param index_format VCF index format. Default tbi. Options [tbi,cbi].
#' @param bgzip_index Create BGZIP index for compressed file. Default FALSE
#' @param output_dir Path to the output directory.
#' @param clean Remove input VCF after completion. Default FALSE.
#' @param verbose Enables progress messages. Default False.
#' @param executor_id Task EXECUTOR ID. Default "gatherBQSR"
#' @param task_name Task name. Default "gatherBQSR"
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param threads Number of threads to split the work. Default 3
#' @param ram [OPTIONAL] If batch mode. RAM memory in GB per job. Default 1
#' @param update_time [OPTIONAL] If batch mode. Show job updates every update time. Default 60
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export


  add_af_strelka_vcf=function(
      bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
      bin_tabix=build_default_tool_binary_list()$bin_tabix,
      vcf_snv=NULL,
      vcf_indel=NULL,
      vcf_sv=NULL,
      verbose=FALSE, 
      executor_id=make_unique_id("addAFStrelka"),
      task_name="addAFStrelka",
      batch_config=build_default_preprocess_config(),
      mode="local",time="48:0:0",
      threads=1,ram=4,update_time=60,
      wait=FALSE,hold=NULL){

          argg <- as.list(environment())
          task_id=make_unique_id(task_name)
          job=build_job(executor_id=executor_id,task_id=task_id)

          jobs_report=build_job_report(
            job_id=job,
            executor_id=executor_id,
            exec_code=list(), 
            task_id=task_id,
            input_args=argg,
            out_files=list(
              )
            )
          if(!is.null(vcf_snv)){
              jobs_report[["steps"]][["annotateAFsnvStrelka"]]<-add_snv_af_strelka_vcf(
                  bin_bgzip=bin_bgzip,
                  bin_tabix=bin_tabix,
                  vcf=vcf_snv,
                  compress=TRUE,index=TRUE,
                  index_format="tbi",bgzip_index=FALSE,
                  clean=TRUE,verbose=verbose, 
                  executor_id=task_id,
                  batch_config=batch_config,
                  mode=mode,time=time,
                  threads=1,ram=1,
                  hold=hold
              )
          }
            
        

        if(!is.null(vcf_indel)){
          jobs_report[["steps"]][["annotateAFindelStrelka"]]<-add_indel_af_strelka_vcf(
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            vcf=vcf_indel,
            compress=TRUE,index=TRUE,
            index_format="tbi",
            bgzip_index=FALSE,
            clean=TRUE,verbose=verbose, 
            executor_id=task_id,
            batch_config=batch_config,
            mode=mode,time=time,
            threads=threads,ram=ram,
            hold=hold
        )
        }

     if(!is.null(vcf_sv)){
        
        jobs_report[["steps"]][["annotateAFsvStrelka"]]<-add_sv_af_strelka_vcf(
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            vcf=vcf_indel,
            compress=TRUE,index=TRUE,
            index_format="tbi",
            bgzip_index=FALSE,
            clean=TRUE,verbose=verbose, 
            executor_id=task_id,
            batch_config=batch_config,
            mode=mode,time=time,
            threads=threads,ram=ram,
            hold=hold
        )
     }


    jobs_report$out_files=list(
        vcf_snv=unlist_lvl( jobs_report[["steps"]][["annotateAFsnvStrelka"]],var="compressed_vcf"),
        vcf_indel=unlist_lvl(jobs_report[["steps"]][["annotateAFindelStrelka"]],var="compressed_vcf"),
        vcf_sv=unlist_lvl(jobs_report[["steps"]][["annotateAFsvStrelka"]],var="compressed_vcf")
    )

  return(jobs_report)

    }











#' Unnest VCF body
#' 
#' Unnest VCF body structure
#' 
#' By default only SAMPLE,FORMAT and VALUE columns will be unnested
#' To unnested also INFO columns argmument full has to be provided
#' 
#' 
#' 
#' @param vcf_body Nested VCF body data structure
#' @param full Fully unnest vcf body data structure. Default FALSE
#' @export


unnest_vcf_body=function(vcf_body,full=FALSE){
  vcf_body=vcf_body %>% dplyr::ungroup() %>% 
    tidyr::unnest(c(SAMPLE,FORMAT,VALUE)) %>%
    tidyr::unnest(c(FORMAT,VALUE))
  if(full){
   vcf_body=vcf_body %>% tidyr::unnest(FILTER)
  }
  return(vcf_body)
}



unnest_to_column_vcf=function(vcf_body,column="INFO"){
    vcf_body=vcf_body %>% tidyr::unnest_wider(var(column),names_sep=paste0(column,"_"))
    return(vcf_body)
}


#' Removes variants with from VCF file vased on the values in the FILTER column
#'
#' VCF datastructure is in list format and contains a header, a body and
#' the corresponding col_names
#'
#' @param bin_bgzip Path to bgzip executable.
#' @param bin_tabix Path to TABIX executable.
#' @param vcf Path to the input VCF file.
#' @param filters Filters to filter VCF by. Default PASS
#' @param exclusive Only keep variants with exactly these filters. Default TRUE.
#' @param output_name Output file name.
#' @param compress Compress VCF file. Default TRUE.
#' @param index Index VCF file. Default TRUE.
#' @param index_format VCF index format. Default tbi. Options [tbi,cbi].
#' @param bgzip_index Create BGZIP index for compressed file. Default FALSE
#' @param output_dir Path to the output directory.
#' @param clean Remove input VCF after completion. Default FALSE.
#' @param verbose Enables progress messages. Default False.
#' @param executor_id Task EXECUTOR ID. Default "gatherBQSR"
#' @param task_name Task name. Default "gatherBQSR"
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param threads Number of threads to split the work. Default 3
#' @param ram [OPTIONAL] If batch mode. RAM memory in GB per job. Default 1
#' @param update_time [OPTIONAL] If batch mode. Show job updates every update time. Default 60
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export
#' 

variants_by_filters_vcf=function(
  rdata=NULL,
  selected=NULL,
  bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
  bin_tabix=build_default_tool_binary_list()$bin_tabix,
  vcf="",filters="PASS",
  output_name="",
  exclusive=TRUE,
  compress=TRUE,index=TRUE,
  index_format="tbi",
  bgzip_index=FALSE,
  clean=TRUE,
  output_dir=".",verbose=FALSE,sep="\t",
  batch_config=build_default_preprocess_config(),
  threads=4,ram=4,mode="local",
  executor_id=make_unique_id("variantsByFiltersVCF"),
  task_name="variantsByFiltersVCF",time="48:0:0",
  update_time=60,wait=FALSE,hold=NULL
){

  if(!is.null(rdata)){
    load(rdata)
  if(!is.null(selected)){
      vcf=vcf_list[selected]
    }
  }


  id=""
  if(output_name!=""){ 
    id=output_name
  }else{
    id=paste0(get_file_name(vcf),".",paste0(filters,collapse="."))
  }


  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir)
  job=build_job(executor_id=executor_id,task_id=task_id)



  job_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    exec_code=list(),
    task_id=task_id,
    input_args = argg,
    out_file_dir=out_file_dir,
    out_files=list()
  )


  

  vcf=read_vcf(vcf=vcf,sep=sep)

  
  if(exclusive){
    vcf$body=vcf$body %>% dplyr::filter(lengths(FILTER)==length(filters))
  }

  vcf$body=vcf$body %>% dplyr::filter(eval(parse(text=paste0(lapply(
    filters,FUN=function(x){paste0("grepl(\"",x,"\",FILTER)")}),collapse="&"))))

  vcf$descriptors$variants_by_filters<-paste0(filters,collapse=";")
  
  job_report[["steps"]][["writeVCF"]]<-write_vcf(
    bin_bgzip=bin_bgzip,
    bin_tabix=bin_tabix,
    compress=compress,index=index,
    index_format=index_format,
    verbose=verbose,
    bgzip_index=bgzip_index,clean=clean,
    output_dir=out_file_dir,executor_id=task_id,
    mode=mode,time=time,threads=threads,ram=ram,
    hold=hold,
    vcf=vcf,output_name=id
  )
  return(job_report)
}





#' Removes variants from VCF file based on the values in the FILTER column across multiple VCFS in parallel
#'
#' VCF datastructure is in list format and contains a header, a body and
#' the corresponding col_names
#'
#' @param bin_bgzip Path to bgzip executable.
#' @param bin_tabix Path to TABIX executable.
#' @param vcf Path to the input VCF file.
#' @param filters Filters to filter VCF by. Default PASS
#' @param exclusive Only keep variants with exactly these filters. Default TRUE.
#' @param output_name Output file name.
#' @param compress Compress VCF file. Default TRUE.
#' @param index Index VCF file. Default TRUE.
#' @param index_format VCF index format. Default tbi. Options [tbi,cbi].
#' @param bgzip_index Create BGZIP index for compressed file. Default FALSE
#' @param output_dir Path to the output directory.
#' @param clean Remove input VCF after completion. Default FALSE.
#' @param verbose Enables progress messages. Default False.
#' @param executor_id Task EXECUTOR ID. Default "gatherBQSR"
#' @param task_name Task name. Default "gatherBQSR"
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param threads Number of threads to split the work. Default 3
#' @param ram [OPTIONAL] If batch mode. RAM memory in GB per job. Default 1
#' @param update_time [OPTIONAL] If batch mode. Show job updates every update time. Default 60
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export
#' 


parallel_vcfs_variants_by_filters_vcf=function(
  bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
  bin_tabix=build_default_tool_binary_list()$bin_tabix,
  vcf="",filters="PASS",
  exclusive=TRUE,
  compress=TRUE,index=TRUE,
  index_format="tbi",
  bgzip_index=FALSE,
  clean=TRUE,
  output_dir=".",verbose=FALSE,sep="\t",
  batch_config=build_default_preprocess_config(),
  threads=4,ram=4,mode="local",
  executor_id=make_unique_id("parVCFvariantsByFiltersVCF"),
  task_name="parVCFvariantsByFiltersVCF",time="48:0:0",
  update_time=60,wait=FALSE,hold=NULL
){

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir)
  tmp_dir=set_dir(dir=out_file_dir,name="tmp")
 

  job=build_job(executor_id=executor_id,task_id=task_id)


  jobs_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    exec_code=list(), 
    task_id=task_id,
    input_args=argg,
    out_file_dir=out_file_dir,
    out_files=list(
      )
    )



    vcf_list=vcf
    names(vcf_list)=Vectorize(get_file_name)(vcf)

    if(mode=="local"){
      jobs_report[["steps"]][["par_vcf_filter_variants"]]<-
      parallel::mclapply(vcf_list,FUN=function(vcf){
        job_report <- variants_by_filters_vcf(
              bin_bgzip=bin_bgzip,
              bin_tabix=bin_tabix,
              vcf=vcf,
              filters=filters,
              exclusive=exclusive,
              compress=compress,
              index=index,
              index_format=index_format,
              bgzip_index=bgzip_index,
              sep=sep,
              clean=clean,
              output_dir=out_file_dir,
              verbose=verbose,
              threads=threads,
              executor_id=task_id)
      },mc.cores=threads)
    }else if(mode=="batch"){
            rdata_file=paste0(tmp_dir,"/",job,".vcfs.RData")
            output_dir=out_file_dir
            executor_id=task_id
            save(vcf_list,
              bin_bgzip,
              bin_tabix,
              filters,
              exclusive,
              compress,
              index,
              sep,
              index_format,
              bgzip_index,
              executor_id,
              output_dir,
              clean,
              verbose,file = rdata_file)
            exec_code=paste0("Rscript -e \"ULPwgs::variants_by_filters_vcf(rdata=\\\"",
            rdata_file,"\\\",selected=$SGE_TASK_ID)\"")
            out_file_dir2=set_dir(dir=out_file_dir,name="batch")
            batch_code=build_job_exec(job=job,time=time,ram=ram,
            threads=1,output_dir=out_file_dir2,
            hold=hold,array=length(vcf_list))
            exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)

            if(verbose){
                print_verbose(job=job,arg=argg,exec_code=exec_code)
            }
            error=execute_job(exec_code=exec_code)
            if(error!=0){
                stop("vcf_process failed to run due to unknown error.
                Check std error for more information.")
            }
 
          jobs_report[["steps"]][["par_vcf_filter_variants"]]<- build_job_report(
                job_id=job,
                executor_id=executor_id,
                exec_code=exec_code, 
                task_id=task_id,
                input_args=argg,
                out_file_dir=out_file_dir,
                out_files=list(
                  extract_vcf=paste0(out_file_dir,"/",
                  names(vcf_list),".",paste0(filters,collapse="."),".vcf")
                )
          )
             }
  

  if(wait&&mode=="batch"){
    job_validator(job=unlist_lvl(jobs_report[["steps"]],var="job_id"),time=update_time,
    verbose=verbose,threads=threads)
  }

  return(jobs_report)

}




#' Extract PASS variants in VCF file
#'
#' VCF datastructure is in list format and contains a header, a body and
#' the corresponding col_names
#'
#' @param bin_bgzip Path to bgzip executable.
#' @param bin_tabix Path to TABIX executable.
#' @param vcf Path to the input VCF file.
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param executor_id Task EXECUTOR ID. Default "gatherBQSR"
#' @param task_name Task name. Default "gatherBQSR"
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param threads Number of threads to split the work. Default 3
#' @param ram [OPTIONAL] If batch mode. RAM memory in GB per job. Default 1
#' @param update_time [OPTIONAL] If batch mode. Show job updates every update time. Default 60
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export
#' 



extract_pass_variants_strelks_vcf=function(
      bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
      bin_tabix=build_default_tool_binary_list()$bin_tabix,
      vcf_snv=NULL,
      vcf_indel=NULL,
      vcf_sv=NULL,
      output_dir=".",verbose=FALSE,sep="\t",
      batch_config=build_default_preprocess_config(),
      threads=4,ram=4,mode="local",
      executor_id=make_unique_id("extractPASStrelka"),
      task_name="extractPASStrelka",time="48:0:0",
      update_time=60,wait=FALSE,hold=NULL
    ){

      argg <- as.list(environment())
      task_id=make_unique_id(task_name)
      out_file_dir=set_dir(dir=output_dir)
      job=build_job(executor_id=executor_id,task_id=task_id)

      jobs_report=build_job_report(
        job_id=job,
        executor_id=executor_id,
        exec_code=list(), 
        task_id=task_id,
        input_args=argg,
        out_file_dir=out_file_dir,
        out_files=list(
          )
        )

   if(!is.null(vcf_snv)){

    jobs_report[["steps"]][["extractPASSsnvVCF"]]<-
        variants_by_filters_vcf(
          bin_bgzip=bin_bgzip,
          bin_tabix=bin_tabix,
          vcf=vcf_snv,
          filters="PASS",
          output_name="somatic.snv",
          exclusive=TRUE,
          compress=TRUE,
          clean=TRUE,
          output_dir=out_file_dir,
          verbose=verbose,
          batch_config=batch_config,
          threads=1,ram=1,mode=mode,
          executor_id=task_id,
          time=time,
          hold=hold
        )

   }

    if(!is.null(vcf_indel)){

    jobs_report[["steps"]][["extractPASSindelVCF"]]<-
        variants_by_filters_vcf(
          bin_bgzip=bin_bgzip,
          bin_tabix=bin_tabix,
          vcf=vcf_indel,
          filters="PASS",
          output_name="somatic.indel",
          exclusive=TRUE,
          compress=TRUE,
          clean=TRUE,
          output_dir=out_file_dir,
          verbose=verbose,
          batch_config=batch_config,
          threads=1,ram=1,mode=mode,
          executor_id=task_id,
          time=time,
          hold=hold
       )
    }

  if(!is.null(vcf_sv)){

    jobs_report[["steps"]][["extractPASSindelVCF"]]<-
        variants_by_filters_vcf(
          bin_bgzip=bin_bgzip,
          bin_tabix=bin_tabix,
          vcf=vcf_sv,
          filters="PASS",
          output_name="somaticSV",
          exclusive=TRUE,
          compress=TRUE,
          clean=TRUE,
          output_dir=out_file_dir,
          verbose=verbose,
          batch_config=batch_config,
          threads=1,ram=1,mode=mode,
          executor_id=task_id,
          time=time,
          hold=hold
       )

  }

    jobs_report$out_files=list(
        vcf_snv=unlist_lvl( jobs_report[["steps"]][["extractPASSsnvVCF"]],var="compressed_vcf"),
        vcf_indel=unlist_lvl(jobs_report[["steps"]][["extractPASSindelVCF"]],var="compressed_vcf"),
        vcf_sv=unlist_lvl(jobs_report[["steps"]][["extractPASSsvVCF"]],var="compressed_vcf")
    )


      return(jobs_report)
  }







#' Tabulate INFO column information extracted from VCF file
#'
#' VCF datastructure with tabulated INFO column body
#'
#' @param vcf VCF datastructure with tabulated INFO column body

tabulate_info_vcf=function(vcf){
  vcf$body=vcf$body %>% unnest_wider(INFO)  
  return(vcf)
}


#' Return tabulated INFO column information to single INFO column
#'
#' VCF datastructure with tabulated INFO column body
#'
#' @param vcf VCF datastructure with tabulated INFO column body

untabulate_info_vcf=function(tab_vcf){
  columns=intersect(names(tab_vcf$descriptors$INFO),names(tab_vcf$body))
  tab_vcf$body=tab_vcf$body%>%
  pivot_longer(names_to="name_INFO",values_to="value_INFO",columns)
  names(tab_vcf$body$value_INFO)=tab_vcf$body$name_INFO
  tab_vcf$body=tab_vcf$body %>% na.omit(value_INFO)%>% 
  group_by(across(-c(name_INFO,value_INFO))) %>%
  summarise(INFO=list(as.list(value_INFO)))
  return(tab_vcf)
}




#' Return tabulated INFO column information to single INFO column
#'
#' VCF datastructure with tabulated INFO column body
#'
#' @param vcf 

extract_csq_info_vcf=function(vcf){
  cols_to_save=c(setdiff(names(vcf$body),"INFO"),"CSQ")
  cols_in_csq=unlist(strsplit(strsplit(
    split="Format: ",gsub("\"","",
    vcf$descriptors$INFO$CSQ$Description))[[1]][2],split="\\|"))
  vcf=tabulate_info_vcf(vcf)
  vcf$body=vcf$body %>% dplyr::select(cols_to_save)%>% rowwise() %>% mutate(CSQ=list(as.list(dplyr::bind_rows(lapply(unlist(stringr::str_split(CSQ,pattern=",")),
  FUN=function(x){lst=unlist(stringr::str_split(x,pattern="\\|"));
  names(lst)=cols_in_csq;lst}))))) %>% tidyr::unnest_wider(CSQ)
  return(vcf)
}







#' Generate tabulated table from annotate VCF
#'
#' VCF datastructure with tabulated INFO column body
#'
#' @param vcf Path to VCF file
#' @param output_dir Path to the output directory.
#' @param threads [OPTIONAL] Number of threads to split the work. Default 4
#' @param ram [OPTIONAL] RAM memory to asing to each thread. Default 4
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Task EXECUTOR ID. Default "recalCovariates"
#' @param task_name Task name. Default "recalCovariates"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export 


tabulate_vcf=function(
  rdata=NULL,
  selected=NULL,
  vcf=NULL,
  ns="ULPwgs",
  output_dir=".",
  mode="local",
  time="48:0:0",
  threads=1,
  ram=1,
  verbose=FALSE,
  batch_config=build_default_preprocess_config(),
  executor_id=make_unique_id("tabVCF"),
  task_name="tabVCF",
  hold=NULL
){


  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir)
  job=build_job(executor_id=executor_id,task_id=task_id)
  func_name=as.character(rlang::call_name(rlang::current_call()))
 

  
  envir=environment()

  main_tabulate_vcf=function(
      vcf=NULL,
      output_dir=".",
      output_name="",
      executor_id=make_unique_id("main"),
      task_name="main"
    ){

      argg <- as.list(environment())
      task_id=make_unique_id(task_name)
      out_file_dir=set_dir(dir=output_dir)
      
      job=build_job(executor_id=executor_id,task_id=task_id)


      id=""
      if(output_name!=""){ 
        id=output_name
      }else{
        id=get_file_name(vcf)
      }


      out_file=paste0(out_file_dir,"/",id,".tabulated.tsv")


      vcf=read_vcf(vcf)
      vcf=extract_csq_info_vcf(vcf)

      ##### Extract body information from VCF

      vcf_body=vcf$body %>% unnest(cols=Allele:TRANSCRIPTION_FACTORS)
      vcf_body=vcf_body %>% unnest_vcf_body() %>% 
      pivot_wider(values_from=VALUE,names_from=c(SAMPLE,FORMAT))

      #### Write to file 


      job_report=build_job_report(
        job_id=job,
        executor_id=executor_id,
        exec_code=list(), 
        task_id=task_id,
        input_args=argg,
        out_file_dir=out_file_dir,
        out_files=list(
          vfc_tab=out_file
          )
        )

      write.table(vcf_body,file=out_file,sep="\t",
      quote=FALSE,row.names=FALSE,col.names=TRUE)

      return(job_report)

  }

    if(!is.null(selected)){
        vcf=slist[selected]
        job_report=main_tabulate_vcf(
          vcf=vcf,
          output_dir=out_file_dir,
          output_name=get_file_name(vcf),
          executor_id=task_id
        )
        return(job_report)
        
    }else{

      slist=vcf
      names(slist)=Vectorize(get_file_name)(slist)

      jobs_report=run_job(
        envir=envir,
        slist=slist,
        job_id=job,
        output_dir=out_file_dir,
        nspace=ns,
        fun=func_name,
        mode=mode,
        hold=hold,
        time=time,
        verbose=verbose,
        threads=threads,
        ram=ram,
        batch_config=batch_config
      )
      return(jobs_report)
    }

      
}




filter_format_vcf=function(
  vcf,exprs=NULL,
  descriptor_id=NULL,
  exclude=TRUE,
  filter_name="filter",
  overwrite_previous_filters=FALSE,
  preserve=TRUE
){
  if(is.null(exprs)){
    stop("Expression is missing")
  }
  if(is.null(descriptor_id)){
    stop("Column descriptor ID is missing")
  }
  vcf$body=check_body_format_vcf(vcf$body)
  vcf$body=vcf$body %>% dplyr::group_by(across(-all_of(c("FORMAT","VALUE","SAMPLE")))) %>% 
    dplyr::mutate(FIL=any(eval(parse(text=paste0(
     "VALUE[across(all_of(column))==descriptor_id&SAMPLE==vcf$descriptors$tumor_sample]",
     exprs)))))

  if(preserve){
    vcf$body=vcf$body %>% mutate(FILTER=ifelse(ifelse(exclude,!FIL,FIL),
    purrr::list_merge(FILTER,filter_name),FILTER))
    vcf[["descriptors"]][["FILTER"]][[filter_name]][["Description"]]=paste0("\"",descriptor_id,exprs,"\"")
  }else{
    vcf$body=vcf$body %>% dplyr::filter(ifelse(exclude,!FIL,FIL))
  }
  vcf$body=vcf$body %>% dplyr::select(-FIL) %>% nest_vcf_body()
  return(vcf)
}





check_body_format_vcf=function(vcf_body){
  if(typeof(vcf_body$FORMAT)=="list"&
     typeof(vcf_body$VALUE)=="list"&
     typeof(vcf_body$SAMPLE)=="list"){
       vcf_body=vcf_body%>% unnest_vcf_body()
     }
    return(vcf_body)
}



#' Nest VCF body
#' 
#' Nest VCF body structure
#' 
#' By default only SAMPLE,FORMAT and VALUE columns will be unnested
#' To unnested also INFO columns argmument full has to be provided
#' 
#' 
#' 
#' @param vcf_body Nested VCF body data structure
#' @param full Fully unnest vcf body data structure. Default FALSE
#' @export


nest_vcf_body=function(vcf_body,full=FALSE){
  vcf_body=vcf_body %>% dplyr::ungroup() %>% 
    tidyr::nest(FORMAT=FORMAT,VALUE=VALUE) %>% 
    tidyr::nest(SAMPLE=SAMPLE,FORMAT=FORMAT,VALUE=VALUE) 
  if(full){
   vcf_body=vcf_body %>% tidyr::nest(FILTER=FILTER)
  }
  return(vcf_body)
}






#' Extract VCF header descriptors
#' This function extracts VCF header descriptors
#' 
#'
#' @param vcf_header VCF header structure
#' @export

extract_descriptors_vcf=function(vcf_header){
    extract_key_value=function(row){
      split=unlist(stringi::stri_split_fixed(row,pattern="=",n=2))
      id=split[1]
      if(grepl("<",split[2])&&grepl(">",split[2])){
           tmp=unlist(lapply(
                unlist(
                  stringi::stri_split_fixed(
                    gsub("<|>","",split[2]),
                    pattern=",",n=stringr::str_count(split[2],"="))
                ),FUN=extract_key_value),
                recursive=FALSE)
          
          values=list(tmp[-1])
          names(values)<-tmp[1]

        }else{
            values=split[2]
        }

      hash=list(
        values=values
      )

    names(hash)<- id
    return(hash)
  }
  rslt=unlist(lapply(gsub("##","",vcf_header),FUN=extract_key_value),recursive=FALSE)
  rslt=Map(function(x) Reduce(append, rslt[names(rslt)==x]), unique(names(rslt)))
  return(rslt)
}



#' Writes VCF file from a VCF data.structure
#'
#' VCF datastructure is in list format and contains a header, a body and
#' the corresponding col_names
#'
#' @param bin_bgzip Path to bgzip executable.
#' @param bin_tabix Path to TABIX executable.
#' @param vcf Path to the input VCF file.
#' @param output_name Output file name.
#' @param compress Compress VCF file. Default TRUE.
#' @param index Index VCF file. Default TRUE.
#' @param index_format VCF index format. Default tbi. Options [tbi,cbi].
#' @param bgzip_index Create BGZIP index for compressed file. Default FALSE
#' @param output_dir Path to the output directory.
#' @param clean Remove input VCF after completion. Default FALSE.
#' @param verbose Enables progress messages. Default False.
#' @param executor_id Task EXECUTOR ID. Default "gatherBQSR"
#' @param task_name Task name. Default "gatherBQSR"
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param threads Number of threads to split the work. Default 3
#' @param ram [OPTIONAL] If batch mode. RAM memory in GB per job. Default 1
#' @param update_time [OPTIONAL] If batch mode. Show job updates every update time. Default 60
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export
#' 

write_vcf=function(
  bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
  bin_tabix=build_default_tool_binary_list()$bin_tabix,
  vcf="",
  output_name="",
  compress=TRUE,index=TRUE,
  index_format="tbi",bgzip_index=FALSE,
  clean=TRUE,output_dir=".",
  verbose=FALSE, 
  executor_id=make_unique_id("writeVCF"),
  task_name="writeVCF",
  batch_config=build_default_preprocess_config(),
  mode="local",time="48:0:0",
  threads=1,ram=4,update_time=60,
  wait=FALSE,hold=NULL
){  
    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir)
    job=build_job(executor_id=executor_id,task_id=task_id)

    build_header_vcf=function(vcf){
      header=unlist(lapply(names(vcf$descriptors),
        FUN=function(col){
          val=vcf$descriptors[[col]];
        
          if(length(val)>1){
            ids=names(val)
            lapply(ids,FUN=function(id){
                paste0(paste0("##",col,"=<ID=",id),",",paste0(names(val[[id]]),"=",val[[id]],collapse=","),">")
            })
          }else{
            paste0("##",col,"=",val)      }
        }
      )) 
    }

    build_body_vcf=function(vcf){

      to_nest=c("FORMAT",vcf$samples)
      names(to_nest)=to_nest
      vcf_body=vcf$body %>% unnest_vcf_body() %>% 
      tidyr::pivot_wider(names_from=SAMPLE,values_from=VALUE) %>%
      dplyr::group_by(dplyr::across(-to_nest)) %>% 
      dplyr::summarise(dplyr::across(to_nest,list))

      vcf_body=vcf_body%>% dplyr::rowwise() %>%
       dplyr::mutate(across(c("FORMAT",all_of(vcf$samples)),
       function(x) paste0(unlist(x),collapse=":")))%>%
       dplyr::mutate(FILTER=paste0(FILTER,collapse=";"),
       INFO=paste0(paste0(names(unlist(INFO)),
       ifelse(unlist(INFO)!="","=",""),unlist(INFO)),collapse=";"))
     
      sort_vcf_body=function(vcf_body,vcf_descriptors){
          chrom=data.frame(CHROM=names(vcf_descriptors$contig),
          order=seq(1,length(vcf_descriptors$contig)))
          vcf_body=dplyr::left_join(vcf_body,chrom) %>% 
          dplyr::arrange(order,POS) %>% dplyr::select(-order)
          
          return(vcf_body)
       }
    
      return(sort_vcf_body(vcf_body,vcf$descriptors))
    }


    out_file=paste0(out_file_dir,"/",output_name,".vcf")

    job_report=build_job_report(
      job_id=job,
      executor_id=executor_id,
      exec_code=list(),
      task_id=task_id,
      input_args = argg,
      out_file_dir=out_file_dir,
      out_files=list(
        vcf=out_file)
    )


    if(output_name==""){
      stop("File output name can't be empty.")
    }
    
    ###Construct VCF header
    vcf_header=build_header_vcf(vcf)

    ###Construct VCF body
    if(nrow(vcf$body)>1){
      vcf_body=build_body_vcf(vcf=vcf)
      colnames=paste0("#",paste0(names(vcf_body),collapse="\t"))
    }else{
      vcf_body=""
      colnames=paste0("#",paste0(c("CHROM","POS","ID",
      "REF","ALT","QUAL","FILTER","INFO","FORMAT",
       vcf$samples),collapse="\t"))
    }
    
    

    write.table(
      x=vcf_header,
      file=out_file,
      sep="\t",
      quote=FALSE,
      row.names=FALSE,
      col.names=FALSE,
    )

    write.table(
      x=colnames,
      file=out_file,
      sep="\t",
      quote=FALSE,
      row.names=FALSE,
      col.names=FALSE,
      append=TRUE
    )

    write.table(
      x=vcf_body,
      file=out_file,
      sep="\t",
      quote=FALSE,
      row.names=FALSE,
      col.names=FALSE,
      append=TRUE
    )

    if(compress){
      job_report[["steps"]][["CompressAndIndexVCF"]]<-compress_and_index_vcf_htslib(
        bin_bgzip=bin_bgzip,
        bin_tabix=bin_tabix,
        vcf=out_file,compress=compress,
        index=index,index_format=index_format,
        bgzip_index=bgzip_index,
        output_dir=out_file_dir,output_name=output_name,
        clean=clean,verbose=verbose,mode=mode,
        batch_config=batch_config,
        time=time,threads=threads,
        ram=ram,hold=hold
    )

  }
  return(job_report)

}




