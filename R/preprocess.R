

#' Preprocess sequencing data
#' 
#'
#' @param sample_sheet Input sample sheet
#' @param config Default tool configure
#' @param vars Default variables
#' @param executor_id Task EXECUTOR ID . Default "preprocessSEQ"
#' @param task_name Task name . Default "preprocessSEQ"
#' @param output_dir Path to output directory
#' @param merge_level Which level to merge samples
#' @param nest_ws Nesting white-space separator
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Executor ID. Default "fastQC"
#' @param task_name Name of the task. Default "fastQC"
#' @param ram RAM memory for batched job. Default 4
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param hold Job to hold on in batched mode.
#' @export


preprocess_seq=function(
    sample_sheet=build_default_sample_sheet(),
    vars_list=build_default_variable_list(),
    config=suppressWarnings(build_default_config()),
    opts_list=build_default_option_list(),
    pmts_list=build_default_parameter_list(),
    steps_list=build_default_steps_list(),
    bin_list=build_default_binary_list(),
    ref_list=build_default_reference_list(),
    merge_level="library",
    nest_ws=1,nesting="",
    executor_id=make_unique_id("preprocessSEQ"), 
    task_name="preprocessSEQ",output_dir=".",
    ram=1,
    header=TRUE,
    sep="",
    threads=1,
    mode="local",
    batch_config=build_default_preprocess_config(),
    time="48:0:0",
    update_time=60,wait=FALSE,hold=NULL,
    verbose=FALSE
){
          
    argg <- as.list(environment())

  
    task_id=make_unique_id(task_name)

    job=build_job(executor_id=executor_id,task=task_id)
   

    
    out_file_dir=set_dir(dir=output_dir)

    job_report=build_job_report(
        job_id=job,
        executor_id=executor_id, 
        task_id=task_id,
        input_args=argg,
        out_file_dir=out_file_dir,
        out_files=list()
      )
    if(!is.na(sample_sheet)){
        
        if(!is.data.frame(sample_sheet)){
                sample_sheet=read.csv(sample_sheet,header=header,sep=sep,
                stringsAsFactors=FALSE)
                if(!header){
                    names(sample_sheet)=c("project_id","patient_id","sample_id","sequencing_type","method_type",
                    "method_version","reference","library_id","R1","R2","step","threads",
                    "ram","batch_config","time","mode","verbose","args")
                }
        }
    }
    sample_sheet=sample_sheet %>% dplyr::arrange(dplyr::across(vars_list[vars_list$required,]$variable))
    
 
    validate_sample_sheet(
    sample_sheet=sample_sheet,
    vars_list=vars_list,opts_list=opts_list)
    
    seq_info=seq_info_check(sample_sheet=sample_sheet,vars_list=vars_list)

    seq_info=suppressMessages(dplyr::left_join(seq_info,
    parameter_config_check(sample_sheet=sample_sheet,config=config,
    vars_list=vars_list,steps_list=steps_list)))
    job=build_job(executor_id=executor_id,task_id=task_id)

    for_id(seq_info=seq_info,output_dir=out_file_dir,
    vars_list=vars_list,nesting=nesting,
    merge_level=merge_level,executor_id = task_id,
    pmts_list=pmts_list,bin_list=bin_list,
    ref_list=ref_list,print_tree=TRUE,ram=ram,threads=threads,
    mode=mode,batch_config=batch_config,
    verbose=verbose,hold=hold,wait=wait,update_time=update_time)

    for_id(seq_info=seq_info,output_dir=output_dir,
    vars_list=vars_list,nesting=nesting,
    merge_level=merge_level,executor_id = task_id,
    pmts_list=pmts_list,bin_list=bin_list,
    ref_list=ref_list,print_tree=FALSE,
    ram=ram,threads=threads,
    mode=mode,batch_config=batch_config,
    verbose=verbose,hold=hold,wait=wait,update_time=update_time)
    cat("COMPLETE!\n")

}


#' For each variable in sample sheet 
#' 
#' Run through each variable in sample sheet and report information
#'
#' @param seq_info Sample sequencing information
#' @param var_list List with variables
#' @param pmts_list List with parameters
#' @param bin_list List with binaries
#' @param ref_list List with references
#' @param batch_config Batch method configuration
#' @param nesting Starting nesting level
#' @param nest_ws Nesting ws to provide. Default 1
#' @param merge_level Level to merge samples. Default library
#' @param executor_id Task EXECUTOR ID . Default "preprocessSEQ"
#' @param task_name Task name . Default "preprocessSEQ"
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Executor ID. Default "fastQC"
#' @param task_name Name of the task. Default "fastQC"
#' @param ram RAM memory for batched job. Default 4
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param hold Job to hold on in batched mode.
#' @export

for_id=function(
    seq_info,output_dir=".",name="",
    vars_list=build_default_variable_list(),
    pmts_list=build_default_parameter_list(),
    bin_list=build_default_binary_list(),
    ref_list=build_default_reference_list(),
    nesting="",merge_level="library",
    nest_ws=1,print_tree=FALSE,
    executor_id=make_unique_id("loopVariables"),
    task_name="loopVariables",
    ram=1,
    threads=1,
    mode="local",
    batch_config=build_default_preprocess_config(),
    time="48:0:0",
    update_time=60,wait=FALSE,hold=NULL,
    verbose=FALSE
   ){              
                
                process_variable=function(
                    ct,seq_info,var,var_text,
                    vars_list_left,info,
                    output_dir=".",name="",
                    pmts_list=build_default_parameter_list(),
                    bin_list=build_default_binary_list(),
                    ref_list=build_default_reference_list(),
                    print_tree=FALSE,
                    nesting="",
                    merge_level="library",
                    nest_ws=1,
                    executor_id=make_unique_id("loopSteps"),
                    task_name="loopSteps",
                    ram=1,
                    threads=1,
                    mode="local",
                    batch_config=build_default_preprocess_config(),
                    clean=TRUE,time="48:0:0",
                    update_time=60,wait=FALSE,hold=NULL,
                    verbose=FALSE){
                                        argg <- as.list(environment())
                                        task_id=make_unique_id(task_name)
                                        merge=FALSE
                                        id=info[ct]
                                        ## Filter sequencing info for id
                                        seq_info_id=seq_info[seq_info[,var,drop=TRUE]==id,]
                                        out_file_dir=set_dir(dir=output_dir,name=id)
                                        new_name=set_name(current_name=name,name=id)
                                        if(print_tree){
                                                if(var!="project_id"){
                                                    cat(add_nesting_ws(nesting,n=nest_ws))
                                            }
                                            if(grepl(merge_level,var)&nrow(seq_info_id %>% dplyr::distinct(path))>2){
                                                merge=TRUE
                                            }
                        
                                            merge_txt=""
                                            if(merge){
                                                merge_txt=crayon::bold(" <<<<===== INFO::SAMPLES WILL BE MERGED AT THIS LEVEL")
                                                seq_info_id[seq_info_id$name=="merge_bam",]$step="TRUE"
                                                samples=seq_info_id %>% dplyr::distinct(path)
                                                samples$last=FALSE
                                                samples[seq(nrow(samples)-1,nrow(samples)),]$last=TRUE
                                                seq_info_id=dplyr::left_join(seq_info_id,samples,by="path")
                                                seq_info_id[seq_info_id$last!=TRUE&seq_info_id$order>5,]$step="FALSE"
                                            }
                                
                                            instrument_name=""
                                            if(var=="flowcell_id"){
                                                instrument_name=paste0("   Platform: ",unique(seq_info_id$instrument_by_flowcell_id))
                                            }

                                            cat(paste0(nesting,"|----",crayon::blue(var_text),crayon::red(id),
                                            crayon::silver(instrument_name),merge_txt,"\n"))
                        
                                            nesting=break_nest(count=ct,info=info,nesting=nesting)


                                        }

                                        
                                        ## Call recursively if variables
                                    
                                        if(length(vars_list_left$variable)>0){
                                            for_id(seq_info=seq_info_id,output_dir=out_file_dir,vars_list=vars_list_left,
                                            nesting=nesting,nest_ws=nest_ws,name=new_name,ref_list=ref_list,print_tree=print_tree,
                                              ram=ram,threads=threads,mode=mode,batch_config=batch_config,verbose=verbose,hold=hold,wait=wait,update_time=update_time)
                                        }else{
                                            tool_config_id=seq_info_id %>% dplyr::select(-c(read_group,path)) %>%  
                                            dplyr::distinct() %>% dplyr::filter(step==TRUE)

                                            seq_info_id=seq_info_id %>% dplyr::select(-c("order",pmts_list$parameter)) %>%  
                                            dplyr::distinct()

                                            out_file_dir_tmp=set_dir(dir=out_file_dir,name="tmp")
                                            out_file_dir_job_report=set_dir(dir=out_file_dir,name="job_report")
                                            

                                            seq_info_R1=seq_info_id[seq_info_id$read_group=="R1",]
                                            seq_info_R2=seq_info_id[seq_info_id$read_group=="R2",]

                                            file_R1=seq_info_R1$path
                                            file_R2=seq_info_R2$path
                                     
                                            

                                            if(print_tree){
                                                cat(add_nesting_ws(nesting,n=nest_ws))
                                                # cat(paste0(nesting,"|----",crayon::green(paste0("R1: ",seq_info_R1$path)),"\n"))
                                                # cat(add_nesting_ws(nesting,n=nest_ws))
                                                # cat(paste0(nesting,"|----",crayon::green(paste0("R2: ",seq_info_R2$path)),"\n")) 
                                                bold_text=FALSE
                                                lapply(seq(1,nrow(tool_config_id)),FUN=function(step){
                                                    
                                                    cat(add_nesting_ws(nesting=nesting,nest="        "))

                                                    if(tool_config_id[step,]$name=="merge_bam"){
                                                            bold_text<<-TRUE
                                                    }
                                                    lapply(seq(1,length(pmts_list$parameter)),FUN=function(pmt){
                                                        
                                                        space="       |"
                                                        if(pmt==ceiling(length(pmts_list$parameter)/2)){
                                                            space=paste0("STEP ",step," |")
                                                        }
                                                    
                                                        txt=paste0(space,pmts_list$text[pmt],
                                                            tool_config_id[step,pmts_list$parameter[pmt]],"\n")

                                                        if(bold_text){
                                                            cat(paste0(nesting,crayon::bold(txt)))
                                                        }else{
                                                            cat(paste0(nesting,txt))
                                                        }

                                                    })
                                

                                                    if(step!=nrow(tool_config_id)){
                                                            if(bold_text){
                                                                cat(add_arrow(nesting=nesting,n=2,bold=TRUE))
                                                            }else{
                                                                cat(add_arrow(nesting=nesting,n=2))
                                                            }
                                                    }else{
                                                        cat(add_nesting_ws(nesting=nesting,nest="        "))
                                                    }
                                
                                                })
                                            }else{
                                
                                                rdata_file=paste0(out_file_dir,"/",new_name,".RData")
                                                save(list=ls(),file = rdata_file)
                                                job=build_job(executor_id=executor_id,task=task_id)
                                                exec_code=paste0("Rscript -e \"ULPwgs::process_sample(rdata=\\\"",rdata_file,"\\\")\"")
                                                if(mode=="batch"){
                                                    out_file_dir2=set_dir(dir=out_file_dir,name="batch")
                                                    batch_code=build_job_exec(job=job,time=time,ram=ram,threads=threads,
                                                    output_dir=out_file_dir2,hold=hold)
                                                    exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
                                                }
                                                if(verbose){
                                                    print_verbose(job=job,exec_code=exec_code)
                                                }
            
                                                error=execute_job(exec_code=exec_code)
                                                if(error!=0){
                                                    stop("Process sample failed to run due to unknown error.
                                                    Check std error for more information.")
                                                }
                                    
                                            }
                                    }
                }

                task_id=make_unique_id(task_name)
                var=vars_list$variable[1]
                var_text=vars_list$text[1]
                vars_list_left=vars_list[-1,]

                info=unique(seq_info[,var,drop=TRUE])



                scroll=seq(1,length(info))
                
            
                lapply(X=scroll,
                FUN=process_variable,
                output_dir=output_dir,
                seq_info=seq_info,
                name=name,
                var=var,
                info=info,
                var_text=var_text,
                vars_list_left=vars_list_left,
                pmts_list=pmts_list,
                bin_list=bin_list,
                ref_list=ref_list,
                nesting=nesting,
                nest_ws=nest_ws,
                print_tree=print_tree,
                merge_level=merge_level,
                executor_id=task_id,
                ram=ram,threads=threads,
                mode=mode,
                batch_config=batch_config,
                verbose=verbose,hold=hold,
                wait=wait,
                update_time=update_time)
                

}




#' For each variable in sample sheet 
#' 
#' Run through each variable in sample sheet and report information
#'
#' @param rdata Path to RData file
#' @export
process_sample=function(rdata=""){
      
                load(file=rdata)
                report=list()
                hold=NULL
                bam=""
                cat("\t\n")
                cat(crayon::magenta(paste0("Processing sample: ",new_name,"\n")))
                cat("\t\n")
                lapply(seq(1,nrow(tool_config_id)),FUN=function(step){
                    if(tool_config_id[step,]$name=="pre_fastqc"){
                        cat("\t\n")
                        cat(crayon::bold("pre_fastqc: \n"))
                        cat("\t\n")
                            report[[new_name]][["steps"]][["pre_fastqc"]]<<-qc_fastqc(
                                bin_fastqc=bin_list$pre_fastqc$bin_fastqc,
                                file_R1=file_R1,
                                file_R2=file_R2,
                                output_dir=paste0(out_file_dir,"/fastqc_reports/pre_trim"),
                                executor_id=task_id,
                                verbose=tool_config_id[step,]$verbose,
                                mode=tool_config_id[step,]$mode,
                                batch_config=tool_config_id[step,]$batch_config,
                                threads=tool_config_id[step,]$threads,
                                ram=tool_config_id[step,]$ram,
                                time=tool_config_id[step,]$time,
                                update_time=60,wait=FALSE,hold=hold)
                    }

            if(tool_config_id[step,]$name=="trimming"){
                    cat("\t\n")
                    cat(crayon::bold("trimming: \n"))
                    cat("\t\n")

                    args=suppressWarnings(parse_args(tool_config_id[step,]$args,step="trimming"))

                    report[[new_name]][["steps"]][["trimming"]]<<-trimming_skewer(
                        bin_skewer=bin_list$trimming$bin_skewer,
                        file_R1=file_R1,
                        file_R2=file_R2,
                        output_dir=out_file_dir,
                        xadapt=args["xadapt",]$value,
                        yadapt=args["yadapt",]$value,
                        mean_quality=args["mean_quality",]$value,
                        min_length=args["min_length",]$value,
                        max_length=args["max_length",]$value,
                        threads=tool_config_id[step,]$threads,
                        output_name=new_name,
                        ram=tool_config_id[step,]$ram,
                        batch_config=tool_config_id[step,]$batch_config,
                        verbose=tool_config_id[step,]$verbose,
                        mode=tool_config_id[step,]$mode,
                        time=tool_config_id[step,]$time,
                        executor_id=task_id,
                        update_time=60,wait=FALSE,hold=hold)

                file_R1 <<- report[[new_name]][["steps"]][["trimming"]]$out_files$r1
                file_R2 <<- report[[new_name]][["steps"]][["trimming"]]$out_files$r2
                hold <<- unlist_lvl(report[[new_name]][["steps"]][["trimming"]],var="job_id",recursive=TRUE)
            }

            if(tool_config_id[step,]$name=="post_fastqc"){
                    cat("\t\n")
                    cat(crayon::bold("post_fastqc: \n"))
                    cat("\t\n")
                report[[new_name]][["steps"]][["post_fastqc"]]<<-qc_fastqc(
                    bin_fastqc=bin_list$pre_fastqc$bin_fastqc,
                    file_R1=file_R1,
                    file_R2=file_R2,
                    output_dir=paste0(out_file_dir,"/fastqc_reports/post_trim"),
                    executor_id=task_id,
                    batch_config=tool_config_id[step,]$batch_config,
                    verbose=tool_config_id[step,]$verbose,
                    mode=tool_config_id[step,]$mode,
                    threads=tool_config_id[step,]$threads,
                    ram=tool_config_id[step,]$ram,
                    time=tool_config_id[step,]$time,
                    update_time=60,wait=FALSE,
                    hold=hold)
            }

            if(tool_config_id[step,]$name=="alignment"){
                cat("\t\n")
                cat(crayon::bold("alignment \n"))
                cat("\t\n")

                args=suppressWarnings(parse_args(tool_config_id[step,]$args,step="alignment"))

                report[[new_name]][["steps"]][["alignment"]]<<-alignment_bwa(
                    bin_bwa=bin_list$alignment$bin_bwa,
                    bin_samtools=bin_list$alignment$bin_samtools,
                    file_R1=file_R1,
                    file_R2=file_R2,
                    output_dir=out_file_dir,
                    id_tag=paste0(seq_info_R1$flowcell_id,".",seq_info_R1$lane_id),
                    pu_tag=paste0(seq_info_R1$flowcell_id,".",seq_info_R1$lane_id,".",
                    seq_info_R1$library_id),
                    pl_tag=seq_info_R1$instrument_by_flowcell_id,
                    lb_tag=seq_info_R1$library_id,
                    sm_tag=new_name,
                    threads=tool_config_id[step,]$threads,
                    ram=tool_config_id[step,]$ram,
                    ref_genome=ref_list[[seq_info_R1$reference]]$reference$genome,
                    coord_sort=as.logical(args["coord_sort",]$value),
                    stats=args["stats",]$value,
                    clean=as.logical(args["clean",]$value),
                    batch_config=tool_config_id[step,]$batch_config,
                    verbose=tool_config_id[step,]$verbose,
                    mode=tool_config_id[step,]$mode,
                    time=tool_config_id[step,]$time,
                    executor_id=task_id,
                    update_time=60,
                    wait=FALSE,
                    hold= hold)
                
                bam <<- report[[new_name]][["steps"]][["alignment"]][["steps"]][["sort_and_index"]][["steps"]][["sort"]]$out_files$bam
                hold <<- unlist_lvl(report[[new_name]][["steps"]][["alignment"]],var="job_id",recursive=TRUE)
            }



              if(tool_config_id[step,]$name=="pre_alignqc"){
                    cat("\t\n")
                    cat(crayon::bold("pre_alignqc: \n"))
                    cat("\t\n")

                    bi=""
                    ti=""
                    ri=""
                    ref_flat=""
                    if(seq_info_R1$method_type=="CAPTURE"|seq_info_R1$method_type=="EXOME"){
                        method="tg"
                        bi=ref_list[[seq_info_R1$reference]][["panel"]][[seq_info_R1$method_version]]$intervals$bi
                        ti=ref_list[[seq_info_R1$reference]][["panel"]][[seq_info_R1$method_version]]$intervals$ti
                    } else if(seq_info_R1$method_type=="WGS"){
                        method="wgs"
                    } else if(seq_info_R1$method_type=="RNASEQ"){
                        method="rna"
                        ri=ref_list[[seq_info_R1$reference]][["rnaseq"]][[seq_info_R1$method_version]]$intervals$ri
                        ref_flat=ref_list[[seq_info_R1$reference]][["rnaseq"]][[seq_info_R1$method_version]]$reference$ref_flat
                    }
                
                    args=suppressWarnings(parse_args(tool_config_id[step,]$args,step="pre_alignqc"))

                    report[[new_name]][["steps"]][["metrics_pre_alignqc"]] <<- metrics_alignqc(
                        bin_samtools=bin_list$pre_alignqc$bin_samtools,
                        bin_picard=bin_list$pre_alignqc$bin_picard,
                        bin_bedtools=bin_list$pre_alignqc$bin_bedtools,
                        bam=bam,
                        output_dir=paste0(out_file_dir,"/alignqc_reports/pre_alignqc"),
                        ref_genome=ref_list[[seq_info_R1$reference]]$reference$genome,
                        verbose=tool_config_id[step,]$verbose,
                        tmp_dir=out_file_dir,
                        bi=bi,
                        ti=ti,
                        ri=ri,
                        ref_flat=ref_flat,
                        method=method,
                        executor_id=task_id,
                        batch_config=tool_config_id[step,]$batch_config,
                        mode=tool_config_id[step,]$mode,
                        time=tool_config_id[step,]$time,
                        threads=tool_config_id[step,]$threads,
                        ram=tool_config_id[step,]$ram,
                        update_time=60,wait=FALSE,
                        hold=hold
                    )
                }  

            if(!merge){

                if(tool_config_id[step,]$name=="markdups"){
                    cat("\t\n")
                    cat(crayon::bold("markdups: \n"))
                    cat("\t\n")
                    
                    args=suppressWarnings(parse_args(tool_config_id[step,]$args,step="markdups"))

                    report[[new_name]][["steps"]][["markdups"]]<<-markdups_gatk(
                        sif_gatk=bin_list$markdups$bin_gatk,
                        bam=bam,
                        output_dir=out_file_dir,
                        remove_duplicates=as.logical(args["remove_duplicates",]$value),
                        batch_config=tool_config_id[step,]$batch_config,
                        mode=tool_config_id[step,]$mode,
                        verbose=tool_config_id[step,]$verbose,
                        threads=tool_config_id[step,]$threads,
                        ram=tool_config_id[step,]$ram,
                        time=tool_config_id[step,]$time,
                        tmp=out_file_dir_tmp,
                        executor_id=task_id,
                        update_time=60,wait=FALSE,
                        hold=hold)
                    bam <<- report[[new_name]][["steps"]][["markdups"]]$out_files$bam
                    hold <<- unlist_lvl(report[[new_name]][["steps"]][["markdups"]],var="job_id",recursive=TRUE)
                    }
            
                if(tool_config_id[step,]$name=="recalibrate"){
                    cat("\t\n")
                    cat(crayon::bold("recalibrate: \n"))
                    cat("\t\n")
                
                    args=suppressWarnings(parse_args(tool_config_id[step,]$args,step="recalibrate"))
                    
                        report[[new_name]][["steps"]][["recalibrate"]] <<- recal_gatk(
                            bin_samtools=bin_list$recalibrate$bin_samtools,
                            sif_gatk=bin_list$recalibrate$bin_gatk,
                            bin_picard=bin_list$recalibrate$bin_picard,
                            bam=bam,tmp_dir=out_file_dir_tmp,
                            output_dir=out_file_dir,
                            ref_genome=ref_list[[seq_info_R1$reference]]$reference$genome,
                            dbsnp=ref_list[[seq_info_R1$reference]]$database$all_common,
                            clean=as.logical(args["clean",]$value),
                            verbose=tool_config_id[step,]$verbose,
                            batch_config=tool_config_id[step,]$batch_config,
                            threads=tool_config_id[step,]$threads,
                            ram=tool_config_id[step,]$ram,
                            mode=tool_config_id[step,]$mode,
                            time=tool_config_id[step,]$time,
                            executor_id=task_id,
                            update_time=60,wait=FALSE,
                            hold=hold)
                        bam <<- report[[new_name]][["steps"]][["recalibrate"]][["steps"]][["sort_and_index"]][["steps"]][["sort"]]$out_files$bam
                        hold <<- unlist_lvl(report[[new_name]][["steps"]][["recalibrate"]],var="job_id",recursive = TRUE)
                    
                    }


                if(tool_config_id[step,]$name=="post_alignqc"){
                    cat("\t\n")
                    cat(crayon::bold("post_alignqc: \n"))
                    cat("\t\n")

                    bi=""
                    ti=""
                    ri=""
                    ref_flat=""
                    if(seq_info_R1$method_type=="CAPTURE"|seq_info_R1$method_type=="EXOME"){
                        method="tg"
                        bi=ref_list[[seq_info_R1$reference]][["panel"]][[seq_info_R1$method_version]]$intervals$bi
                        ti=ref_list[[seq_info_R1$reference]][["panel"]][[seq_info_R1$method_version]]$intervals$ti
                    } else if(seq_info_R1$method_type=="WGS"){
                        method="wgs"
                    } else if(seq_info_R1$method_type=="RNASEQ"){
                        method="rna"
                        ri=ref_list[[seq_info_R1$reference]][["rnaseq"]][[seq_info_R1$method_version]]$intervals$ri
                        ref_flat=ref_list[[seq_info_R1$reference]][["rnaseq"]][[seq_info_R1$method_version]]$reference$ref_flat
                    }
                
                    args=suppressWarnings(parse_args(tool_config_id[step,]$args,step="post_alignqc"))

                    report[[new_name]][["steps"]][["metrics_post_alignqc"]] <<- metrics_alignqc(
                        bin_samtools=bin_list$post_alignqc$bin_samtools,
                        bin_picard=bin_list$post_alignqc$bin_picard,
                        bin_bedtools=bin_list$post_alignqc$bin_bedtools,
                        bam=bam,
                        output_dir=paste0(out_file_dir,"/alignqc_reports/post_alignqc"),
                        ref_genome=ref_list[[seq_info_R1$reference]]$reference$genome,
                        verbose=tool_config_id[step,]$verbose,
                        tmp_dir=out_file_dir,
                        bi=bi,
                        ti=ti,
                        ri=ri,
                        ref_flat=ref_flat,
                        method=method,
                        executor_id=task_id,
                        batch_config=tool_config_id[step,]$batch_config,
                        mode=tool_config_id[step,]$mode,
                        time=tool_config_id[step,]$time,
                        threads=tool_config_id[step,]$threads,
                        ram=tool_config_id[step,]$ram,
                        update_time=60,wait=FALSE,
                        hold=hold
                    )
                }
            }    
                    
        })

        save(report,file=paste0(out_file_dir_job_report,"/",new_name,".job_report.RData"))
       
}


        
       



