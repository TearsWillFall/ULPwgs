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



unnest_to_column=function(vcf_body,column="INFO"){
    vcf_body=vcf_body %>% tidyr::unnest_wider(var(column),names_sep=paste0(column,"_"))
    return(vcf_body)
}





vcf=dat



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

filter_format_vcf


tmp=vcf %>% filter_vcf(exprs=">0.6",exclude=FALSE,descriptor_id="AF",filter_name="af_above_0.6")


tmp

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




#' Add sv AF information to Strelka VCF file
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
  executor_id=make_unique_id("addAFsvStrelka"),
  task_name="addAFsvStrelka",
  batch_config=build_default_preprocess_config(),
  mode="local",time="48:0:0",
  threads=1,ram=4,update_time=60,
  wait=FALSE,hold=NULL){
    
    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file=paste0(get_file_name(vcf),".sv")
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

    vcf_dat$body=vcf_dat$body %>% unnest_vcf_body() %>% 
      dplyr::group_by_at(dplyr::vars(-VALUE,-FORMAT)) %>% 
      dplyr::group_modify(~dplyr::add_row(.x,FORMAT="AF")
    )
   
   
    ##Extract Tier1 read information for REF and ALT


    vcf_dat$body=vcf_dat$body %>% dplyr::mutate(
      UREFPR=tryCatch({strsplit(VALUE[FORMAT=="PR"],split=",")[[1]][1]},error=function(e){NA}),
      UALTPR=tryCatch({strsplit(VALUE[FORMAT=="PR"],split=",")[[1]][2]},error=function(e){NA}),
      UREFSR=tryCatch({strsplit(VALUE[FORMAT=="SR"],split=",")[[1]][1]},error=function(e){NA}),
      UALTSR=tryCatch({strsplit(VALUE[FORMAT=="SR"],split=",")[[1]][2]},error=function(e){NA})
    ) %>%
      dplyr::mutate(
        VALUE=ifelse(
          FORMAT=="AF",list(
            as.numeric(UALTPR)/(as.numeric(UREFPR)+as.numeric(UALTPR)),
            as.numeric(UALTSR)/(as.numeric(UREFSR)+as.numeric(UALTSR)))
        ,VALUE))%>% dplyr::select(-c(UALTPR,UREFPR,UALTSR,UREFSR))

    vcf_dat$body=vcf_dat$body %>% 
    dplyr::mutate(VALUE=ifelse(is.na(VALUE),"",VALUE))%>%
    nest_vcf_body()
    

    add_af_descriptor<-function(){
        list(
          Number="2",Type="Float",
          Description="\"Variant allelic frequency for supporting  reads and split-reads in the listed order\"")
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

    build_header_vcf=function(vcf_descriptors){
      header=unlist(lapply(names(vcf_descriptors),
        FUN=function(col){
          val=vcf_descriptors[[col]];
        
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

    build_body_vcf=function(vcf_body,vcf_samples){
      to_nest=c("FORMAT",vcf_samples)
      names(to_nest)=to_nest
      vcf_body=sv$body %>% unnest_vcf_body() %>% 
      tidyr::pivot_wider(names_from=SAMPLE,values_from=VALUE) %>%
      dplyr::group_by(dplyr::across(-to_nest)) %>% 
      dplyr::summarise(dplyr::across(to_nest,list))

      vcf_body=vcf_body%>% dplyr::rowwise() %>%
       dplyr::mutate(across(c("FORMAT",all_of(vcf_samples)),
       function(x) paste0(unlist(x),collapse=":")))%>%
       dplyr::mutate(FILTER=paste0(FILTER,collapse=","),
       INFO=paste0(paste0(names(unlist(INFO)),
       ifelse(unlist(INFO)!="","=",""),unlist(INFO)),collapse=","))
      return(vcf_body)
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
    vcf_header=build_header_vcf(vcf_descriptors=vcf$descriptors)

    ###Construct VCF body
    vcf_body=build_body_vcf(vcf_body=vcf$body,vcf_samples=vcf$samples)

    write.table(
      x=vcf_header,
      file=out_file,
      sep="\t",
      quote=FALSE,
      row.names=FALSE,
      col.names=FALSE,
    )

    write.table(
      x=paste0("#",paste0(names(vcf_body),collapse="\t")),
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




