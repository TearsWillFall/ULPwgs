#' Get the filename of input file
#' This function takes the absolute/relative path to a file and
#' returns its base name without the file extension suffix.
#'
#' @param bam_path Path to the input file
#' @return A string with the name of the file
#' @export


get_file_name=function(file_path=""){
  filename=unlist(strsplit(basename(file_path),"\\."))[1]
  return(filename)
}

#' Set the current dir
#' This function sets and/or creates the directory
#' for a function
#'
#' @param dir Path to the set the directory
#' @param name Name of the directory
#' @return A path to the new directory
#' @export


set_dir=function(dir="",name=""){
  sep="/"

  if(dir==""){
    sep=""
  }

  new_dir=paste0(dir,sep,name)

  if (!dir.exists(new_dir)){
      dir.create(new_dir,recursive=TRUE)
  }
  return(new_dir)
}


#' Set the current name
#' This function sets a new a name for a level
#'
#' @param current_name Path to the set the directory
#' @param name Name of the level
#' @return A path to the new directory
#' @export


set_name=function(current_name="",name=""){
  sep="_"

  if(current_name==""){
    sep=""
  }

  new_name=paste0(current_name,sep,name)

  return(new_name)
}


#' @export

unlist_lvl=function(named_list,var,recursive=FALSE){

  vars=names(named_list)
  lvl_found=any(vars %in% var)

  if(lvl_found){
      out=named_list[[var]]
      if(recursive){
        if(!is.null(names(named_list))){
            out=append(out,unlist(lapply(names(named_list),FUN=function(name){
            out=unlist_lvl(named_list=named_list[[name]],var)
        })))
      }
    }
  }else{
    out=unlist(lapply(names(named_list),FUN=function(name){
    out=unlist_lvl(named_list=named_list[[name]],var)
    }))
  }
  return(out)
}






#' @export

print_verbose=function(exec_code,arg,job,ws=1){
      rep(cat("    \n"),ws)
      cat(crayon::blue("Job:"))
      rep(cat("    \n"),ws)
      cat(paste0(crayon::red(job,"\n")))
      rep(cat("    \n"),ws)
      cat(crayon::blue("Arguments:"))
      rep(cat("    \n"),ws)
      lapply(names(arg),FUN=function(ag){
          rep(cat("    \n"),ws)
          cat(paste0(crayon::silver(ag),": ",arg[[ag]],"\n"))
      })
      rep(cat("    \n"),ws)
      cat(crayon::blue("Command:"))
      rep(cat("    \n"),ws)
      cat(paste0(crayon::green(exec_code,"\n")))
      rep(cat("    \n"),ws)
}


#' Check required columns in sample sheet
#' 
#'
#' @param sample_sheet Input sample sheet
#' @param req_cols Required columns in sample sheet
#' @export

check_req_cols=function(sample_sheet=build_default_sample_sheet(),
req_cols){
  cols_in_sheet=names(sample_sheet)[names(sample_sheet) %in% req_cols]
  miss_cols=cols_in_sheet[!req_cols %in% cols_in_sheet]
  len_miss_cols=length(miss_cols)
  if(len_miss_cols>0){
   stop(paste0("Error: The following columns are required and were not found in the sample sheet: ",
    paste(miss_cols,collapse=", "),
    "."))
  }

}



#' Check required variable values for a column
#' 
#'
#' @param sheet_col Sample sheet column
#' @param col_name Column name
#' @param types Allowed values types for columns
#' @export


check_req_types=function(sheet_col,col_name,types){
  types_in_col=unique(sheet_col[sheet_col %in% types])
  wrong_types=unique(sheet_col[!sheet_col %in% types])
  len_wrong_types=length(wrong_types)

  if(len_wrong_types>0){
   stop(paste0("Error: Column ",col_name," contains wrong value types in sample sheet: ",
   paste(wrong_types,collapse=", "),
    ". Allowed values for column ",col_name," are : ",paste(types,collapse=", ")))
  }
  return(TRUE)
}




#' Validate sample sheet input columns
#' 
#'
#' @param sample_sheet Input sample sheet
#' @param vars_list Default variable list
#' @export


validate_sample_sheet=function(sample_sheet=build_default_sample_sheet(),
  vars_list=build_default_variable_list(),opts_list=build_default_option_list()){
    req_cols=vars_list$variable[vars_list$required==TRUE]
    check_req_cols(req_cols=req_cols)
    req_type_cols=vars_list$variable[vars_list$needs_type_validation==TRUE]
    rtrn=lapply(req_type_cols,FUN=function(col){
          check_req_types(sheet_col=sample_sheet[,col,drop=TRUE],
           col_name=col,types=unlist(opts_list[col]))

    })
}

#' Parse all tool parameters from sample sheet
#' 
#'
#' @param sample_sheet Input sample sheet
#' @param config Default tool configure
#' @export



parse_tool_parameters=function(sample_sheet=build_default_sample_sheet(),
config=suppressWarnings(build_default_config()),steps_list=build_default_steps_list()){
  parameters=names(config)[names(config)!="name"&names(config)!="order"]
  parameters_in_sheet=names(sample_sheet)[names(sample_sheet) %in% parameters]
  missing_parameters=parameters[!parameters %in% parameters_in_sheet]
  len_missing_parameters=length(missing_parameters)
  if(len_missing_parameters>0){
    message(paste0("Warning: The following columns were not found for tool parameters in the sample sheet: ",
    paste(missing_parameters,collapse=", "),
    ". These parameters will be initiate with default configuration."))
  }
  default_config=config
  rslt=lapply(parameters_in_sheet,FUN=function(parameter){
    rslt=lapply(sample_sheet[,parameter],FUN=function(steps){
        step_list=unlist(strsplit(steps,";"))
          rslt=lapply(step_list,FUN=function(step){
            step_value_list=strsplit(step,"=")
            step_name=step_value_list[[1]][1]
            parameter_values=step_value_list[[1]][-1]
            if(parameter=="args"){
              parameter_values=paste0(parameter_values,collapse="=")
              parameter_values=gsub("\\{|\\}","",parameter_values)
              input_args=suppressWarnings(parse_args(args=parameter_values,step=step_name,steps_list=steps_list))
              default_args=suppressWarnings(parse_args(args=default_config[step_name,parameter],
              step=step_name,steps_list=steps_list))

              validated_args=validate_input_args(input_args,default_args)
              default_config[step_name,parameter]<<-parameter_values
            }else{
                   default_config[step_name,parameter]<<-parameter_values
            }
          })
      })
  
  })
  return(default_config)
}



#' Check if value is missing
#' 
#'
#' @param var Value of variable
#' @export

check_missing=function(var){
  return(kutils::isNA(var)|is.null(var))
}



#' Unnest job report
#' 
#'
#' @param job_report Job report to read
#' @param index Index level to read at
#' @export

read_job_report=function(job_report,index=1){
    lapply(job_report,FUN="[[",index=index)
}




#' Validate tool argument for step
#' 
#'
#' @param arg Tool argument list
#' @param step Pipeline step
#' @param step_list Pipeline step list
#' @export


validate_arg=function(step,arg,steps_list=build_default_steps_list()){
  
  val_arg=steps_list[[step]][["args"]][[arg]]

  if(is.null(val_arg)){
    stop(paste0(arg, " is an invalid argument for step: ",step))
  }

}

#' Parse tool arguments
#' 
#'
#' @param args Tool argument list
#' @param step Pipeline step
#' @param step_list Pipeline step list
#' @export


parse_args=function(args,step,steps_list=build_default_steps_list()){
  args_list=strsplit(args,"\\|")
  out=lapply(args_list[[1]],FUN=function(arg){
                  arg_value_list=strsplit(arg,"=")
                  validate_arg(step=step,arg=arg_value_list[[1]][1],steps_list=steps_list)
                  out=list(arg=arg_value_list[[1]][1],
                  value=arg_value_list[[1]][2])
  })
  out=dplyr::bind_rows(out)
  out$step=step
  row.names(out)=out$arg
  return(out)
}

#' Validate tool args inputs
#' 
#'
#' @param input Input tool argument values
#' @param default Default tool argument values
#' @export


validate_input_args=function(input,default){
  undetermined=input[!input$arg %in% default$arg,]
  n_und=nrow(undetermined)
  if(n_und>0){
    message(paste0("Warning: Ignoring input argument",
    ifelse(n_und==1,": ","s: ")),
    paste0(undetermined$arg,collapse=", "),". ",ifelse(n_und==1,
    "This argument is not valid ","These arguments are not valid "),
    "for ",unique(undetermined$step), " step.")
  }
  determined=input[input$arg %in% default$arg,]
  default[default$arg %in%determined$arg,"value"]=determined$value
  return(default)
}




#' Find sequencing instrument name from sequecing information
#' 
#' Find sequencing instrument name using instrument_id and/or flowcell_id for a sequencing sample
#' 
#' @param instrument_id Dataframe with matched data for instrument_id and instrument name
#' @param flowcell_id Dataframe with matched data for flowcell_id, flowcell_type and instrument name
#' @param seq_info Sequencing info from a sample
#' @return A dataframe with found matching instrument names
#' @export


find_instrument=function(instrument_id=build_instrument_id(),
  flowcell_id=build_flowcell_id(),seq_info){
  instrument=lapply(instrument_id$pattern,FUN=function(instrument_pattern){
    grepl(instrument_pattern,seq_info$instrument_id,perl=TRUE)
    })
  instrument_id_found=instrument_id[unlist(instrument),]

  flowcell=lapply(flowcell_id$pattern,FUN=function(flowcell_pattern){
    grepl(flowcell_pattern,seq_info$flowcell_id,perl=TRUE)
    })
  
  flowcell_id_found=flowcell_id[unlist(flowcell),]

  found=list(instrument_by_flowcell_id=flowcell_id_found$instrument,
  flowcell_type=flowcell_id_found$flowcell,instrument_by_intrument_id=instrument_id_found$instrument)
  
  return(found)
}



#' Infer sequencing information from sequencing file
#' 
#' Infer sequencing information for:
#' instrument
#' run
#' flowcell
#' lane 
#'
#' @param bin_samtools Path to samtools binary. Default "tools/samtools/samtools".
#' @param file_path Path to the input file
#' @return A string with the extension of the file
#' @export


infer_sequencing_info=function(bin_samtools=build_default_tool_binary_list()$bin_samtools,file_path){
    header=extract_read_header(bin_samtools=bin_samtools,file_path=file_path)
    seq_info=parse_read_header(header)
    instrument=find_instrument(seq_info=seq_info)
    seq_info$flowcell_type=instrument$flowcell_type
    seq_info$instrument_by_flowcell_id=instrument$instrument_by_flowcell_id
    seq_info$instrument_by_insttrument_id=instrument$instrument_by_instrument_id
    return(seq_info)
}

#' Extract the a header for a random read from a sequencing file
#'
#' Extract the header from the top most read from any sequencing file
#' Acepted formats: fasta | bam | sam 
#' Compressed fasta also accepted
#'
#' @param bin_samtools Path to samtools binary. Default "tools/samtools/samtools".
#' @param file_path Path to the input file
#' @return A string with the extension of the file
#' @export

extract_read_header=function(
  bin_samtools=build_default_tool_binary_list()$bin_samtools,file_path){
  

  file_ext=get_file_ext(file_path)
    
  if(!file.exists(file_path)){
    stop("ERROR: ",file_path," could not be found")
  }
  file.exists(file_path)
  if(grepl("f*q",file_ext)){
    if(check_if_compressed(file_path)){
      read=system(paste0("gunzip -c ",file_path,"| head -n 1 "),intern=TRUE)
    }else{
      read=system(paste0("cat ",file_path,"| head -n 1 "),intern=TRUE)
    }
    
  }else if(grepl("bam$",file_ext)){
    read=system(paste0(bin_samtools," view | head -n 1 | awk '{print $1}' "),intern=TRUE)
  }else{
   stop("Error: Could not figure out file format for file : ",file_path)
  }
}

#' Parse read header
#'
#' Parse read header info to extract ids for:
#' instrument
#' run
#' flowcell
#' lane 
#'
#' @param header Sequencing read header
#' @return A list with sequencing information
#' @export


parse_read_header=function(header){
  info_list=strsplit(header,":")[[1]]
  instrument_id=sub("@","",info_list[1])
  run_id=info_list[2]
  flowcell_id=info_list[3]
  lane_id=info_list[4]
  return(list(
    instrument_id=instrument_id,
    run_id=run_id,
    flowcell_id=flowcell_id,
    lane_id=lane_id
    )
    )
}


#' Check if sample is compressed
#'
#' Check sequencing info for each read group for each sample
#'
#' @param file_path Path to the input file
#' @return A string with the extension of the file
#' @export


check_if_compressed=function(file_path){
  rslt=system(paste0("file ",file_path),intern=TRUE)
  return(grepl("compr",rslt))
}



#' Check sequencing information for samples
#'
#' Check sequencing info for each read group for each sample
#'
#' @param sample_sheet Path to the input file
#' @param vars_list List with variables
#' @return A string with the extension of the file
#' @export


seq_info_check=function(sample_sheet=build_default_sample_sheet(),
vars_list=build_default_variable_list()){
  vars=vars_list$variable
  seq_info=lapply(seq(1,nrow(sample_sheet)),FUN=function(x){
    R1_seq_info=infer_sequencing_info(file_path=sample_sheet[x,]$R1)
    R1_seq_info=append(R1_seq_info,sample_sheet[x,vars[vars %in% names(sample_sheet)]])
    R1_seq_info$read_group="R1"
    R1_seq_info$path=sample_sheet[x,]$R1
    R2_seq_info=infer_sequencing_info(file_path=sample_sheet[x,]$R2)
    R2_seq_info=append(R2_seq_info,sample_sheet[x,vars[vars %in% names(sample_sheet)]])
    R2_seq_info$read_group="R2"
    R2_seq_info$path=sample_sheet[x,]$R2
    seq_info=dplyr::bind_rows(R1_seq_info,R2_seq_info)
    seq_info=seq_info %>% 
      tidyr::pivot_longer(cols=!c(vars[vars %in% names(sample_sheet)],
      read_group,path),names_to="platform",values_to="value")
    seq_info=seq_info %>% dplyr::group_by(platform) %>% 
      dplyr::mutate(validate=value[read_group=="R1"]==value[read_group=="R2"])
  }) %>% dplyr::bind_rows() %>% dplyr::ungroup()%>% tidyr::pivot_wider(values_from=value,names_from=platform)
  if(any(seq_info$validate==FALSE)){
      stop("Miss-matching R1 and R2 information found for : \n",seq_info[seq_info$validate==FALSE,])
  }
  return(seq_info)
}



#' Check parameter configuration in sample sheet
#'
#' Check parameter configuration in sample sheet
#'
#' 

#' @param sample_sheet Dataframe with sample information
#' @param config Default tool config
#' @param vars_list List with variables
#' @return A string with the extension of the file
#' @export


parameter_config_check=function(sample_sheet=build_default_sample_sheet(),
config=suppressWarnings(build_default_config()),vars_list=build_default_variable_list(),
steps_list=build_default_steps_list()){
  vars=vars_list$variable
  tool_configs=lapply(seq(1,nrow(sample_sheet)),FUN=function(x){
    tool_config=parse_tool_parameters(sample_sheet=sample_sheet[x,],config=config,steps_list=steps_list)
    tool_config=append(tool_config,sample_sheet[x,c(vars[vars %in% names(sample_sheet)],"R1","R2")])
    return(tool_config)
  })
  tool_config=dplyr::bind_rows(tool_configs) %>% tidyr::pivot_longer(cols=c(R1,R2),
  names_to="read_group",values_to="path")
  return(tool_config)
}


#' Parse HS metrics and WGS metrics picard info
#'
#' 
#'
#' 

#' @param summary Path to Picard summary file
#' @param output_dir Path to output directory
#' @param output_name File output name
#' @param extract Metric to extract
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param ram RAM memory to use in GB. Default 4.
#' @param executor [OPTIONAL] Task executor name. Default "recalCovariates"
#' @param task [OPTIONAL] Task name. Default "recalCovariates"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export

parse_picard_metrics=function(summary="",output_dir="",output_name="",
extract="stats",verbose=FALSE,threads=1,ram=4,
mode="local",executor_id=make_unique_id("parsePicardMetrics"),
task_name="parsePicardMetrics",time="48:0:0",
update_time=60,wait=FALSE,hold=""){

  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir,name="parse_summary")

  metrics=read.table(summary,sep="\t",quote="\\",nrows=1,header=TRUE)
  metrics$SAMPLE=output_name
  histogram=read.table(summary,sep="\t",quote="\\",skip=8,header=TRUE)
  histogram$SAMPLE=output_name
  job=build_job(executor_id = executor_id,task_id=task_id)

  if(extract=="stats"|extract=="all"){
    write.table(metrics,file=paste0(out_file_dir,"/metrics_stats.txt"),quote=FALSE,row.names=FALSE,col.names=TRUE)
  }

  if(extract=="histogram"|extract=="all"){
    write.table(histogram,file=paste0(out_file_dir,"/metrics_histogram.txt"),quote=FALSE,row.names=FALSE,col.names=TRUE)
  }


}


#' Get the extension of a file
#'
#' This function takes the absolute/relative path to a file and
#' returns the file extension suffix.
#'
#' @param file_path Path to the input file
#' @return A string with the extension of the file
#' @export

get_file_ext=function(file_path=""){
    ext = strsplit(basename(file_path), split="\\.")[[1]]
    ext = paste(ext[-1],collapse=".")
    return(ext)
}

#' Get sample name from two sample replicates
#'
#' This function takes the absolute/relative path to two files
#' and returns the longest common string among their basenames
#'
#' @param file_path Path to the input file
#' @param file_path2 Path to the second input file
#' @return A string with the longest common basename
#' @export

intersect_file_name=function(file_path="",file_path2=""){
  tmp_name=get_file_name(file_path2)
  sample_name=sapply(sapply(c(0:(nchar(tmp_name)-1)),
  function (i) substr(tmp_name,1,nchar(tmp_name)-i)),function (x) grepl(x,file_path))
  sample_name=names(which(sample_name)[1])
  sample_name=sub("(.*)[_.-].*","\\1",sample_name)
  return(sample_name)

}




#' Create BED with antitarget regions with padding
#'
#' This function indexes a genomic sequence file (BAM/SAM).
#'
#' @param bin_bedtools Path to bedtools executable. Default path tools/bedtools2/bin/bedtools.
#' @param bed Path to the input file with the sequence.
#' @param pad Pad distance. Default 10
#' @param output_name Output file name.
#' @param genome Path to genome fa.fai file
#' @param verbose Enables progress messages. Default False.
#' @param threads Number of threads . Default 4
#' @param ram RAM memory. Default 4
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Task EXECUTOR ID. Default "mardupsGATK"
#' @param task_name Task name. Default "mardupsGATK"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param verbose [OPTIONAL] Enables progress messages. Default False.#
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] Hold job until job is finished. Job ID. 
#' @export

complement_bed=function(
  bin_betools=build_default_tool_binary_list()$bin_bedtools,
  bed="",pad=10,output_name="Complement",genome="",verbose=FALSE,
  threads=3,ram=1,coord_sort=TRUE,mode="local",
  executor_id=make_unique_id("complementBED"),clean=TRUE,
  task_name="complementBED",time="48:0:0",update_time=60,wait=FALSE,hold=""
  ){


  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir)

  out_file=paste0(out_file_dir,output_name,".bed")


  job_report=build_job_report(
    job_id=job, 
    executor_id=executor_id, 
    task_id=task_id,
    input_args=argg,
    out_file_dir=out_file_dir,
    out_files=list(
      complement=out_file
    )
  )
  if (pad!=0){
    job_report[["steps"]][["pad_bed"]]=pad_bed(
      bin_bedtools=bin_bedtools,bed=bed,pad=pad,
      output_name=paste0(ULPwgs::get_file_name(bed),"_",pad),
      output_dir=out_file_dir,genome=genome,verbose=verbose)



    exec_code=paste0(bin_samtools," complement -i ",
    job_report[["steps"]][["pad_bed"]]$out_file$pad,
       " -g ", genome, " > ",out_file)
    hold=job_report[["steps"]][["pad_bed"]]$job_id



  }else{
    exec_code=paste0(bin_bedtools," complement -i ",bed, " -g ", genome, " > ",out_file)
  }

  if(mode=="batch"){
        out_file_dir2=set_dir(dir=out_file_dir,name="batch")
        exec_batch=build_job_exec(job=job,time=time,ram=ram,threads=threads,
        output_dir=out_file_dir2,hold=hold)
        exec_code=paste("echo 'source ~/.bashrc;",exec_code,"'|",exec_batch)
  }

  if(verbose){
        print_verbose(job=job,arg=argg,exec_code=exec_code)
  }
      error=system(exec_code)
  if(error!=0){
    stop("bedtools failed to run due to unknown error.
    Check std error for more information.")
  }

  if(wait&&mode=="batch"){
        job_validator(job=job_report$job_id,
        time=update_time,verbose=verbose,threads=threads)
  } 

    return(job_report)
}

#' Pad a BED file
#'
#' This function takes a BED file and pad each regions in both directions
#'
#' @param bed Path to the input file with the sequence.
#' @param bin_bedtools Path to bedtools executable. Default path tools/bedtools2/bin/bedtools.
#' @param pad Pad distance. Default 10
#' @param output_name Output file name.
#' @param output_dir Path to output directory
#' @param genome Path to genome fa.fai file
#' @param verbose Enables progress messages. Default False.
#' @param threads Number of threads . Default 4
#' @param ram RAM memory. Default 4
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Task EXECUTOR ID. Default "mardupsGATK"
#' @param task_name Task name. Default "mardupsGATK"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param verbose [OPTIONAL] Enables progress messages. Default False.#
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] Hold job until job is finished. Job ID. 
#' @export

pad_bed=function(bin_bedtools=build_default_tool_binary_list()$bin_bedtools,bed="",pad=10,
  output_name="Padded",output_dir="",genome="",verbose=FALSE,threads=3,ram=1,
  coord_sort=TRUE,mode="local",executor_id=make_unique_id("padBED"),task_name="padBED",
  time="48:0:0",update_time=60,wait=FALSE,hold=""
){


  argg <- as.list(environment())
  task_id=make_unique_id(task_name)
  out_file_dir=set_dir(dir=output_dir)


  out_file=paste0(out_file_dir,output_name,".bed")
  exec_code=paste0(bin_bedtools," slop -i ",bed, " -g ", genome," -b ",pad, " > ",
    out_file)

  job=build_job(executor_id=executor_id,task_id=task_id)
 
  if(mode=="batch"){
        out_file_dir2=set_dir(dir=out_file_dir,name="batch")
        exec_batch=build_job_exec(job=job,time=time,ram=ram,threads=threads,
        output_dir=out_file_dir2,hold=hold)
        exec_code=paste("echo 'source ~/.bashrc;",exec_code,"'|",exec_batch)
    }
   if(verbose){
         print_verbose(job=job,arg=argg,exec_code=exec_code)
    }
  error=system(exec_code)
  if(error!=0){
    stop("bedtools failed to run due to unknown error.
    Check std error for more information.")
  }

    job_report=build_job_report(
    job_id=job, 
    executor_id=executor_id, 
    task_id=task_id,
    input_args=argg,
    out_file_dir=out_file_dir,
    out_files=list(
      pad=out_file
    )
  )
  return(job_report)
}




#' @export

add_bl=function(){
  break_line="|----"
  return(break_line)
}


#' @export

add_l=function(){
  line="----"
  return(line)
}

#' @export

add_fill=function(fill="\t",n=2){
  paste0(rep(fill,n),collapse="")
}

#' @export

add_nesting_level=function(n=1){
  nest=paste0(" ",add_fill(" ",n=4))

  if(n>0){
    paste0(add_fill(" ",n=4),paste0(base::rep(nest,n),collapse=""),add_bl())
  }else{
    add_fill(" ",n=9)
  }
}

#' @export

add_nesting_ws=function(nesting="",nest="|    ",n=1){
    paste0(rep(paste0(nesting,nest,"\n"),n),collapse="")
}


#' @export


add_nest=function(){
  return("|    ")
}

#' @export

break_nest=function(count,info,nesting){
  if(count==1&length(info)>1){
          nesting=paste0(nesting,add_nest())
  }else{
          nesting=paste0(nesting,"     ")
  }
  return(nesting)
}

#' @export

add_arrow=function(nesting="",n=2,bold=FALSE){
  if(bold){
    txt=paste0(paste0(rep(paste0(nesting,"        ",crayon::bold("      ..   "),"\n"),n),collapse=""),
    paste0(nesting,"        ",crayon::bold("      \\/   "),"\n"))
  }else{
    txt=paste0(paste0(rep(paste0(nesting,"        ","      ..   ","\n"),n),collapse=""),
    paste0(nesting,"        ","      \\/   ","\n"))
  }
  
  return(txt)
}






#' Wrapper around samtools addreplacerg function
#'
#' This functions add/replaces RG tags lines in BAM files
#' This function wraps around samtools addreplacerg function
#'
#' @param bin_samtools [REQUIRED] Path to santools binary. Default tools/samtools/samtools.
#' @param bam [REQUIRED] Path to the BAM file/s.
#' @param output_dir [OPTIONAL] Path to the output directory.
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param index [OPTIONAL] Generate an indexed file. Default False.
#' @param ID [REQUIRED] ID tag for RG tag line.
#' @param PL [OPTIONAL] PL tag for RG tag line.
#' @param PU [OPTIONAL] PU tag for RG tag line.
#' @param LB [OPTIONAL] LB tag for RG tag line.
#' @param SM [OPTIONAL] SM tag for RG tag line.
#' @param threads [OPTIONAL] Number of threads per jobs.
#' @param jobs [OPTIONAL] Number of jobs to run.
#' @export

replace_rg=function(
  bin_samtools=build_default_tool_binary_list()$bin_samtools,bam="",output_dir="",
  verbose=FALSE,index=TRUE,ID="",PL="",PU="",LB="",SM="",threads=3,jobs=1){

  out_file_di=set_dir(dir=output_dir)

  if(ID!=""){
    ID=paste0("ID:",ID)
  }

  if(PL!=""){
    PL=paste0("PL:",PL)
  }

  if(PU!=""){
    PU=paste0("PU:",PU)
  }

  if(LB!=""){
    LB=paste0("LB:",LB)
  }

  if(SM!=""){
    SM=paste0("SM:",SM)
  }

  tag=c(ID,PL,PU,LB,SM)
  tag=tag[!tag==""]

  tag=paste0(" -r ",paste0(tag,collapse=" -r "))

  parallel::mclapply(seq(1,length(bam)),FUN=function(x){
    exec_code=paste(bin_samtools," addreplacerg ",tag," -o ",paste0(out_file_dir,"/",
      basename(sub("bam","rh.bam",bam[x]))), " -@ ",threads,bam[x])
  if(verbose){
      print(exec_code)
    }
    exec_code

    if(index){
      bam_index_samtools(bin_samtools=bin_samtools,bam=paste0(out_file_dir,"/",
      basename(sub("bam","rh.bam",bam[x]))),verbose=verbose,threads=threads)
    }
  },mc.cores=jobs)
}




#' Estimate coverage for BED file
#'
#' This function estimates mean/per_base_coverage coverage for regions in bed file.
#' This function is intended for estimating off target mean coverage in panel data, however it can be used as a standalone.
#' It is highly recommended to use the sorted option and sort the BED and BAM files beforehand as it highly decreases RAM
#' consumption and increases computation speed. The easiest way of doing this is using sort -k1,1 -k2,2n on the bed file.
#' However, this may not work if the BED has been build based on different reference than the BAM. For example, hs37 vs hg19.
#' For more information: https://bedtools.readthedocs.io/en/latest/content/tools/coverage.html
#'
#' @param bin_bedtools Path to bwa executable. Default tools/bedtools2/bin/bedtools.
#' @param bam Path to the input BAM file.
#' @param bed Path to the input bed file.
#' @param sorted Are the input files sorted. Default TRUE
#' @param mean Estimate mean coverage per region. Default TRUE. FALSE produces coverage per base per target.
#' @param hist Enables progress messages. Default False.
#' @param fai Indexed genome to which sequece has been aligned. Default none. Only require if there is an issue with chromosome naming.
#' @param verbose Enables progress messages. Default False.
#' @param output_dir Output directory path. Default none.
#' @export


bed_coverage=function(
  bin_bedtools=build_default_tool_binary_list()$bin_bedtools,bam="",bed="",
  verbose=FALSE,sorted=TRUE,mean=TRUE,fai="",suffix="",output_dir="",hist=FALSE
){
    
    out_file_dir=set_dir(dir=output_dir,name="coverage")

    srt=""
    if (sorted){
      srt="-sorted"
    }

    if (fai!=""){
      fai=paste("-g",fai)
    }
    if (suffix!=""){
      suffix=paste0(".",suffix)
    }
    mode=""
    if (hist){
      mode="-hist"

      ## Filter to reduce the size of the output as it produces the coverage per base stats too

      out_file=paste0("| grep \"all\"",">",out_file_dir,get_file_name(bam),suffix,".Histogram_Coverage.txt")
    }else{
      if (mean){
        mode="-mean"
        out_file=paste0(">",out_file_dir,get_file_name(bam),suffix,".Per_Region_Coverage.txt")
      }else{
        mode="-d"
        out_file=paste0(">",out_file_dir,get_filename(bam),suffix,".Per_Base_Coverage.txt")
      }
    }

    exec_code=paste(bin_bedtools,"coverage -a",bed, "-b" ,bam,fai, mode,srt,out_file)
    if(verbose){
        print(exec_code)
    }
        error=system(exec_code)
  if(error!=0){
    stop("bedtools failed to run due to unknown error.
    Check std error for more information.")
  }

}



#' Function to collect chromosome data in bam
#'
#' This function takes a BAM file and collects the chr names from the bam
#' file.
#'
#' @param bin_samtools Path to samtools executable. Default path tools/samtools/samtools.
#' @param bam Path to directory with BAM files to merge.
#' @param verbose Enables progress messages. Default False.
#' @param threads Number of threads . Default 4
#' @param ram RAM memory. Default 4
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Task EXECUTOR ID. Default "mardupsGATK"
#' @param task_name Task name. Default "mardupsGATK"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param verbose [OPTIONAL] Enables progress messages. Default False.#
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] Hold job until job is finished. Job ID. 
#' @export

get_bam_reference_chr=function(
    bin_samtools=build_default_tool_binary_list()$bin_samtools,
    bam="",output_name="chrReference",
    output_dir="",verbose=FALSE,
    executor_id=make_unique_id("getBAMchr"),task_name="getBAMchr",
    mode="local",time="48:0:0",
    threads=4,ram=4,update_time=60,wait=FALSE,hold=""
  ){
 
    options(scipen = 999)
    
    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir)

    job=build_job(executor_id=executor_id,task_id=task_id)

    
    out_file=paste0(out_file_dir,output_name,".txt")
    exec_code=paste0(bin_samtools," view -H ",bam,
    " | grep @SQ| awk -F  \"\\t|:\" \'{print $3\"\\t\"0\"\\t\"$5}\' |  awk \'BEGIN{print \"chr\\tstart\\tend\"}1\' >",out_file)
    
    job=build_job(executor_id = executor_id,task_id=task_id)
    
    if(mode=="batch"){
    
        out_file_dir2=set_dir(dir=out_file_dir,name="batch")
        exec_batch=build_job_exec(job=job,
        time=time,ram=ram,threads=threads,
        output_dir=out_file_dir2,hold=hold)
        exec_code=paste("echo 'source ~/.bashrc;",exec_code,"'|",exec_batch)
    }


    if(verbose){
      print_verbose(job=job,arg=argg,exec_code=exec_code)
    }


  error=system(exec_code)
  if(error!=0){
    stop("samtools failed to run due to unknown error.
    Check std error for more information.")
  }

  job_report=build_job_report(
    job_id=job,
    executor_id=executor_id,
    task_id=task_id, 
    input_args = argg,
    out_file_dir=out_file_dir,
      out_files=list(
        ref=out_file)
  )

    if(wait&&mode=="batch"){
    batch_validator(job=job_report$job_id,
    time=update_time,verbose=verbose,threads=threads)
  }


  return(job_report)
}




#' Function to collect chromosome data in fasta fai
#'
#' This function takes a FAI file and collects the chr names from the fai
#' file.
#'
#' @param fasta Path to directory with fasta file.
#' @param output_name Output name
#' @param verbose Enables progress messages. Default False.
#' @param threads Number of threads . Default 4
#' @param ram RAM memory. Default 4
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Task EXECUTOR ID. Default "mardupsGATK"
#' @param task_name Task name. Default "mardupsGATK"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param verbose [OPTIONAL] Enables progress messages. Default False.#
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] Hold job until job is finished. Job ID. 
#' @export

get_fai_reference_chr=function(
    fasta="",output_name="chrRef",output_dir="",
    verbose=FALSE,executor_id=make_unique_id("getFAIchr"),
    task_name="getFAIrchr",
    mode="local",time="48:0:0",
    threads=4,ram=4,update_time=60,wait=FALSE,hold=""
  ){
 
    options(scipen = 999)
    
    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir)

    job=build_job(executor_id=executor_id,task_id=task_id)

    
    out_file=paste0(out_file_dir,output_name,".txt")
    exec_code=paste("cat",paste0(fasta,".fai"),
    " | awk \'{print $1\"\\t\"0\"\\t\"$2}\' |  awk \'BEGIN{print \"chr\\tstart\\tend\"}1\' >",out_file)
    
    job=build_job(executor_id = executor_id,task_id=task_id)
    
    if(mode=="batch"){
    
        out_file_dir2=set_dir(dir=out_file_dir,name="batch")
        exec_batch=build_job_exec(job=job,
        time=time,ram=ram,threads=threads,
        output_dir=out_file_dir2,hold=hold)
        exec_code=paste("echo 'source ~/.bashrc;",exec_code,"'|",exec_batch)
    }


    if(verbose){
      print_verbose(job=job,arg=argg,exec_code=exec_code)
    }


    error=system(exec_code)
    if(error!=0){
      stop("samtools failed to run due to unknown error.
      Check std error for more information.")
    }

    job_report=build_job_report(
      job_id=job,
      executor_id=executor_id,
      task_id=task_id, 
      input_args = argg,
      out_file_dir=out_file_dir,
        out_files=list(
          ref=out_file)
    )

      if(wait&&mode=="batch"){
      batch_validator(job=job_report$job_id,
      time=update_time,verbose=verbose,threads=threads)
    }


  return(job_report)
}




#' Function to generate seq with trailing ner
#'
#' This functions genereates a sequence of number
#'
#'
#' @param from Start of sequence
#' @param to End of sequence
#' @param by Step size
#' @export



seqlast <- function (from, to, by)
{
  vec <- do.call(what = seq, args = list(from, to, by))
  if ( tail(vec, 1) != to ) {
    return(c(vec, to))
  } else {
    return(vec)
  }
}



#' Split chromosomes into bins
#'
#' This functions takes a BED file of chromosome positions for BAM file input
#'
#'
#' @param bin_samtools Path to samtools executable. Default path tools/samtools/samtools.
#' @param bam Path to directory with BAM files to merge.
#' @param verbose Enables progress messages. Default False.
#' @param bin_size Bin size. Default 400000000 pb
#' @export


bin_chromosomes <- function(bin_samtools=build_default_tool_binary_list()$bin_samtools,bam="",verbose=FALSE,
bin_size=40000000){
  options(scipen = 999)
  chr=get_bam_reference_chr(bin_samtools=bin_samtools,bam=bam,verbose=verbose)
  bed=chr%>% dplyr::group_by(chr) %>%
  dplyr::summarise(start=seqlast(start,end,bin_size)) %>%
  dplyr::mutate(end=dplyr::lead(start)) %>% tidyr::drop_na()
  bed=bed[stringr::str_order(paste0(bed$chr,"_",bed$start), numeric = TRUE),]
  return(bed)
}






#' Get current script path
#'
#' This functions obtains the current script path
#'
#' @export




thisFile <- function() {
        cmdArgs <- commandArgs(trailingOnly = FALSE)
        needle <- "--file="
        match <- grep(needle, cmdArgs)
        if (length(match) > 0) {
                # Rscript
                return(normalizePath(sub(needle, "", cmdArgs[match])))
        } else {
                # 'source'd via R console
                return(normalizePath(sys.frames()[[1]]$ofile))
        }
}


#' @export


make_triangles <- function(x, y, point = "up") {
  x <- as.integer((x))
  y <- as.integer((y))

  if (point == "up") {
    newx <- sapply(x, function(x) {
      c(x - 0.5, x - 0.5, x + 0.5)
    }, simplify = FALSE)
    newy <- sapply(y, function(y) {
      c(y - 0.5, y + 0.5, y + 0.5)
    }, simplify = FALSE)
  } else if (point == "down") {
    newx <- sapply(x, function(x) {
      c(x - 0.5, x + 0.5, x + 0.5)
    }, simplify = FALSE)
    newy <- sapply(y, function(y) {
      c(y - 0.5, y - 0.5, y + 0.5)
    }, simplify = FALSE)
  }
  data.frame(x = unlist(newx), y = unlist(newy))
}