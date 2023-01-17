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

  
  if(is.null(vcf)){
    stop("vcf arguments is of type NULL")
  }else if(check_if_compressed(vcf)){
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
  header=FALSE,colClasses="character")
  names(body)=read.table(text=sub("#","",col_names),
  stringsAsFactors = FALSE,colClasses="character")
  samples=setdiff(names(body),cols)
  body=extract_body_vcf(body,samples)
  body$POS=as.numeric(body$POS)
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
  vcf=NULL,
  ...
  ){
    
    run_main=function(.env){

        .this.env=environment()
        append_env(to=.this.env,from=.env)
    
        set_main(.env=.this.env)


        .main$steps[[fn]]<-.this.env
        .main.step=.main$steps[[fn]]

        vcf_dat=read_vcf(input)
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
        
      
      

        .main$steps[[fn]]$steps <- append(
          .main$steps[[fn]]$steps,
          write_vcf(
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            vcf=vcf_dat,
            output_name=paste0(input_id,".af"),
            output_dir=out_file_dir,
            tmp_dir=tmp_dir,
            env_dir=env_dir,
            batch_dir=batch_dir
          )
        )
        
        .this.step=.main.step$steps$write_vcf
        .main.step$out_files=.this.step$out_files
        .env$.main<-.main
  }

  
     .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="vcf"
    )
  
    launch(.env=.base.env)
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
  vcf=NULL,
  ...
  ){
    
    run_main=function(.env){

        .this.env=environment()
        append_env(to=.this.env,from=.env)
    
        set_main(.env=.this.env)

        .main$steps[[fn]]<-.this.env
        .main.step=.main$steps[[fn]]

        vcf_dat=read_vcf(input)
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
      
      

        .main$steps[[fn]]$steps <- append(
          .main$steps[[fn]]$steps,
          write_vcf(
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            vcf=vcf_dat,
            output_name=paste0(input_id,".af"),
            output_dir=out_file_dir,
            tmp_dir=tmp_dir,
            env_dir=env_dir,
            batch_dir=batch_dir
          )
        )
        
        .this.step=.main.step$steps$write_vcf
        .main.step$out_files=.this.step$out_files
        .env$.main<-.main
  }

  
     .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="vcf"
    )
  
    launch(.env=.base.env)
  
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
  vcf=NULL,
  ...
){
    
    run_main=function(.env){

        .this.env=environment()
        append_env(to=.this.env,from=.env)
    
        set_main(.env=.this.env)

        .main$steps[[fn]]<-.this.env
        .main.step=.main$steps[[fn]]

        vcf_dat=read_vcf(input)
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
        vcf_dat$body=vcf_dat$body %>% dplyr::filter(!is.na(VALUE)) %>% nest_vcf_body()

  
        add_afp_descriptor<-function(){
            list(
              Number="1",
              Type="Float",
              Description="\"Variant allelic frequency by paired reads spanning the region.\"")
        }

        add_afs_descriptor<-function(){
            list(
              Number="1",
              Type="Float",
              Description="\"Variant allelic frequency by split-reads spanning the region.\"")
        }

        vcf_dat$descriptors$FORMAT[["AFP"]]<-add_afp_descriptor()
        vcf_dat$descriptors$FORMAT[["AFS"]]<-add_afs_descriptor()
        
      
      

        .main.step$steps <- append(
          .main.step$steps,
          write_vcf(
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            vcf=vcf_dat,
            output_name=paste0(input_id,".af"),
            output_dir=out_file_dir,
            tmp_dir=tmp_dir,
            env_dir=env_dir,
            batch_dir=batch_dir
          )
        )
        
        .this.step=.main.step$steps$write_vcf
        .main.step$out_files=.this.step$out_files
        .env$.main<-.main

    }

   
    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="vcf"
    )
  
    launch(.env=.base.env)
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



add_af_strelka_vcf=function(
  bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
  bin_tabix=build_default_tool_binary_list()$bin_tabix,
  vcf=NULL,
  type="snv",
  ...
){
  
  run_main=function(
    .env
  ){


      .this.env=environment()
      append_env(to=.this.env,from=.env)
  
      set_main(.env=.this.env)

      output_name=paste0(input_id,".",type)
      fn=paste0(fn,".",type)


      .main$steps[[fn]]<-.this.env
      .main.step=.main$steps[[fn]]

      if(type=="snv"){


        .main$steps[[fn]]$steps<-append(
        .main$steps[[fn]]$steps,
          add_snv_af_strelka_vcf(
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            vcf=input,
            output_dir=out_file_dir,
            tmp_dir=tmp_dir,
            env_dir=env_dir,
            batch_dir=batch_dir,
            output_name=output_name,
            verbose=verbose,
            threads=threads,
            ram=ram,
            executor_id=task_id
          )
        )
        .this.step=.main.step$steps$add_snv_af_strelka_vcf
        .main.step$out_files[[type]]=.this.step$out_files
      }else if(type=="indel"){
        .main$steps[[fn]]$steps<-append(
        .main$steps[[fn]]$steps,
          add_indel_af_strelka_vcf(
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            vcf=input,
            output_dir=out_file_dir,
            tmp_dir=tmp_dir,
            env_dir=env_dir,
            batch_dir=batch_dir,
            output_name=output_name,
            verbose=verbose,
            threads=threads,
            ram=ram,
            executor_id=task_id
          )
        )
        .this.step=.main.step$steps$add_indel_af_strelka_vcf
        .main.step$out_files[[type]]=.this.step$out_files
      }else if(type=="sv"){
        .main$steps[[fn]]$steps<-append(
        .main$steps[[fn]]$steps,
          add_sv_af_strelka_vcf(
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            vcf=input,
            output_dir=out_file_dir,
            tmp_dir=tmp_dir,
            env_dir=env_dir,
            batch_dir=batch_dir,
            output_name=output_name,
            verbose=verbose,
            threads=threads,
            ram=ram,
            executor_id=task_id
          )
        )
      .this.step=.main.step$steps$add_sv_af_strelka_vcf
      .main.step$out_files[[type]]=.this.step$out_files 
      }else{
        stop("Wrong type argument")
      }

      .env$.main<-.main


  }
  .base.env=environment()
  list2env(list(...),envir=.base.env)
  set_env_vars(
    .env= .base.env,
    vars="vcf"
  )

  launch(.env=.base.env)

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
  bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
  bin_tabix=build_default_tool_binary_list()$bin_tabix,
  vcf=NULL,
  filters="PASS",
  exclusive=TRUE,
  ...
){

  run_main=function(.env){

    .this.env=environment()
    append_env(to=.this.env,from=.env)

    set_main(.env=.this.env)

     output_name=paste0(input_id,".",paste0(filters,collapse="."))

    .main$steps[[fn]]<-.this.env
    .main.step<-.main$steps[[fn]]
   

    vcf_dat=read_vcf(vcf=input)

    if(exclusive){
      vcf_dat$body=vcf$body %>% dplyr::filter(lengths(FILTER)==length(filters))
    }

    vcf_dat$body=vcf_dat$body %>% dplyr::filter(eval(parse(text=paste0(lapply(
      filters,FUN=function(x){paste0("grepl(\"",x,"\",FILTER)")}),collapse="&"))))

    vcf_dat$descriptors$variants_by_filters<-paste0(filters,collapse=";")

    .main.step$steps <- append(
      .main.step$steps,
      write_vcf(
        bin_bgzip=bin_bgzip,
        bin_tabix=bin_tabix,
        vcf=vcf_dat,
        output_name=output_name,
        output_dir=out_file_dir,
        tmp_dir=tmp_dir,
        env_dir=env_dir,
        batch_dir=batch_dir,
        err_msg=err_msg
      )
    )

    .this.step=.main.step$steps$write_vcf
    .main.step$out_files=append(
      .main.step$out_files,
      .this.step$out_files
    )
    .env$.main<-main

  }

     
  .base.env=environment()
  list2env(list(...),envir=.base.env)
  set_env_vars(
    .env= .base.env,
    vars="vcf"
  )
  launch(.env=.base.env)

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



extract_pass_variants_strelka_vcf=function(
      bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
      bin_tabix=build_default_tool_binary_list()$bin_tabix,
      vcf=NULL,
      type="snv",
      ...
    ){

    run_main=function(.env){


      .this.env=environment()
      append_env(to=.this.env,from=.env)
  
      set_main(.env=.this.env)

      output_name=paste0(input_id,".",type)
      fn=paste0(fn,".",type)


      .main$steps[[fn]]<-.this.env
      .main.step=.main$steps[[fn]]
 
     

     .main$steps[[fn]]$steps<-append(
        .main$steps[[fn]]$steps,
          variants_by_filters_vcf(
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            vcf=input,
            filters="PASS",
            output_name=output_name,
            output_dir=out_file_dir,
            tmp_dir=tmp_dir,
            env_dir=env_dir,
            batch_dir=batch_dir,
            verbose=verbose,
            threads=threads,
            ram=ram,
            executor_id=task_id
          )
        )
      .this.step=.main.step$steps$variants_by_filters_vcf
      .main.step$out_files[[type]]=.this.step$out_files

      .env$.main<-.main

    }

    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="vcf"
    )
  
    launch(.env=.base.env)

  }






#' Tabulate INFO column information extracted from VCF file
#'
#' VCF datastructure with tabulated INFO column body
#'
#' @param vcf VCF datastructure with tabulated INFO column body
#' @export

tabulate_info_vcf=function(vcf){
  vcf$body=vcf$body %>% tidyr::unnest_wider(INFO)  
  return(vcf)
}


#' Return tabulated INFO column information to single INFO column
#'
#' VCF datastructure with tabulated INFO column body
#'
#' @param vcf VCF datastructure with tabulated INFO column body
#' @export

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
#' @export

extract_csq_info_vcf=function(vcf){
  cols_to_save=c(setdiff(names(vcf$body),"INFO"),"CSQ")
  cols_in_csq=unlist(strsplit(strsplit(
    split="Format: ",gsub("\"","",
    vcf$descriptors$INFO$CSQ$Description))[[1]][2],split="\\|"))
  vcf=tabulate_info_vcf(vcf)
  vcf$body=vcf$body %>% dplyr::select(cols_to_save)%>% dplyr::rowwise() %>% 
  dplyr::mutate(CSQ=list(as.list(dplyr::bind_rows(lapply(unlist(stringr::str_split(CSQ,pattern=",")),
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
  vcf=NULL,
  ...
){


    run_main=function(
        .env
      ){


        .this.env=environment()
        append_env(to=.this.env,from=.env)

        set_main(.env=.this.env)

        .main$out_files$tab_vcf=paste0(out_file_dir,"/",input_id,".tabulated.tsv")


        ignore=TRUE
        #### Catch empty VCF 
        tryCatch({vcf=read_vcf(input);ignore=FALSE},error=function(error){
            warning("No variants detected in VCF. Ignoring input VCF.")
        })


        if(ignore){

            write.table("No variants detected in input VCF.",file=.main$out_files)
    
        }else{
            vcf=extract_csq_info_vcf(vcf)

            ##### Extract body information from VCF

            vcf_body=vcf$body %>% tidyr::unnest(cols=Allele:TRANSCRIPTION_FACTORS)
            vcf_body=vcf_body %>% unnest_vcf_body(full=TRUE) %>% 
            tidyr::pivot_wider(values_from=VALUE,names_from=c(SAMPLE,FORMAT))
            write.table(vcf_body,file=steps$out_files,sep="\t",quote=FALSE,
            row.names=FALSE,col.names=TRUE)
        }

        .main$steps[[fn]]<-.this.env
    

        .env$.main<- .main 

    }
  
        
    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="vcf"
    )
  
    launch(.env=.base.env)


      
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
  output_name=NULL,
  vcf=NULL,
  compress=TRUE,
  index=TRUE,
  ...
){  



  run_main=function(.env){


    .this.env=environment()
    append_env(to=.this.env,from=.env)

    set_main(.env=.this.env)


    .main$steps[[fn]] <- .this.env
    .main.step=.main$steps[[fn]]

    .main.step$out_files$vcf <- paste0(out_file_dir,"/",input_id,".vcf")

      build_header_vcf=function(vcf){
        header=unlist(lapply(names(vcf$descriptors),
          FUN=function(col){
            val=vcf$descriptors[[col]];
          
            if(length(val)>1){
              ids=names(val)
              lapply(ids,FUN=function(id){
                  paste0(paste0("##",col,"=<ID=",id),
                  ",",paste0(names(val[[id]]),"=",val[[id]],collapse=","),">")
              })
            }else{
              paste0("##",col,"=",val) }
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
            dplyr::arrange(order,POS) %>% 
            dplyr::select(-order)
            
            return(vcf_body)
          }
      
        return(sort_vcf_body(vcf_body,vcf$descriptors))
    }



    

      ###Construct VCF header
      vcf_header=build_header_vcf(vcf)

      ###Construct VCF body
      if(nrow(vcf$body)>1){
        vcf_body=build_body_vcf(vcf)
        colnames=paste0("#",paste0(names(vcf_body),collapse="\t"))
      }else{
        vcf_body=""
        colnames=paste0("#",paste0(c("CHROM","POS","ID",
        "REF","ALT","QUAL","FILTER","INFO","FORMAT",
        vcf$samples),collapse="\t"))
      }
      

      write.table(
        x=vcf_header,
        file=.main.step$out_files$vcf,
        sep="\t",
        quote=FALSE,
        row.names=FALSE,
        col.names=FALSE,
      )

      write.table(
        x=colnames,
        file=.main.step$out_files$vcf,
        sep="\t",
        quote=FALSE,
        row.names=FALSE,
        col.names=FALSE,
        append=TRUE
      )

      write.table(
        x=vcf_body,
        file=.main.step$out_files$vcf,
        sep="\t",
        quote=FALSE,
        row.names=FALSE,
        col.names=FALSE,
        append=TRUE
      )

    

      if(compress){
        .main.step$steps<-append(
          .main.step$steps,
          compress_and_index_vcf_htslib(
            bin_bgzip=bin_bgzip,
            bin_tabix=bin_tabix,
            vcf=.main.step$out_files,
            compress=compress,
            index=index,
            index_format=index_format,
            bgzip_index=bgzip_index,
            output_dir=out_file_dir,
            tmp_dir=tmp_dir,
            env_dir=env_dir,
            output_name=input_id,
            verbose=verbose,
            threads=threads,
            ram=ram
          ) 
        )
     .this.step=.main.step$steps$compress_and_index_vcf_htslib
     .main.step$out_files <- append(
      .main.step$out_files,
      .this.step$out_files
      )
    }
    .env$.main<-.main
  }

  .base.env=environment()
  list2env(list(...),envir=.base.env)
  set_env_vars(
    .env= .base.env,
    vars="output_name"
  )
 
  launch(.env=.base.env)

}




