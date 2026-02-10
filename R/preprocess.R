

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
#' Process UMI-tagged sequencing data
#'
#' Processes sequencing data that include Unique Molecular Identifiers (UMIs).
#' The function runs a multi-step pipeline that extracts UMIs, trims adapters,
#' aligns reads, groups reads by UMI, collapses UMIs to consensus sequences and
#' remaps consensus reads. It is intended to reduce PCR and sequencing errors
#' by collapsing reads originating from the same original molecule.
#'
#' @details
#' The pipeline performs the following logical stages (implemented as ordered
#' steps in the function): extraction of UMIs, read trimming, mapping of
#' trimmed reads, tagging/merging to preserve UMI information, grouping by UMI,
#' consensus calling (with fgbio), conversion between BAM/FASTQ where necessary,
#' remapping consensus reads and final tagging/merging.
#'
#' Each step logs progress and is wrapped in error handling so failures report
#' which pipeline step failed and abort early.
#'
#' @param sif_gatk Path to the GATK Singularity image (used for GATK-based helpers).
#'   Defaults to \code{build_default_sif_list()$sif_gatk}.
#' @param env_fgbio Python environment identifier for fgbio tools. Defaults to
#'   \code{build_default_python_enviroment_list()$env_fg_bio}.
#' @param env_fastp Python environment identifier for fastp. Defaults to
#'   \code{build_default_python_enviroment_list()$env_fastp}.
#' @param bin_bwa Path or identifier for the BWA binary. Defaults to
#'   \code{build_default_binary_list()$alignment$bin_bwa}.
#' @param bin_samtools Path or identifier for the samtools binary. Defaults to
#'   \code{build_default_binary_list()$alignment$bin_samtool}.
#' @param ref_genome Path to the reference genome FASTA used for alignments.
#'   Defaults to \code{build_default_reference_list()$HG19$reference$genome}.
#' @param fastq Named list or object with FASTQ file paths. Expected keys are
#'   \code{fastq_r1} and \code{fastq_r2} (or a single fastq value for single-end).
#' @param project_id Optional project identifier used when constructing output
#'   directory structure.
#' @param patient_id Optional patient identifier used when constructing output
#'   directory structure and read-group tags.
#' @param sample_id Optional sample identifier; used as \code{input_id} when
#'   building output filenames and read-group tags.
#' @param sequencing_type method_type,method_version,reference Optional strings
#'   describing sequencing assay details (used in output path creation).
#' @param library_id,run_id,flowcell_id,lane_id Optional run/library identifiers.
#'   If not provided the function attempts to infer them from the FASTQ files
#'   using \code{new_check_seq_info}.
#' @param ... Additional arguments forwarded to internal helper functions (for
#'   example \code{tmp_dir}, \code{batch_dir}, \code{threads}, \code{ram},
#'   \code{verbose}, \code{executor_id} and \code{batch_config}).
#'
#' @return Invisibly returns the internal `.main` object containing pipeline
#'   steps and collected output paths; the function primarily writes files and
#'   job reports to the specified output directories.
#'
#' @examples
#' \dontrun{
#' fastq <- list(fastq_r1 = "sample_R1.fastq.gz", fastq_r2 = "sample_R2.fastq.gz")
#' preprocess_umi(fastq = fastq, patient_id = "P001", sample_id = "S001",
#'                output_dir = "./results", threads = 8, ram = 16)
#' }
#'
#' @export

preprocess_umi=function(
    sif_gatk=build_default_sif_list()$sif_gatk,
    env_fgbio=build_default_python_enviroment_list()$env_fgbio,
    env_fastp=build_default_python_enviroment_list()$env_fastp,
    bin_bwa=build_default_binary_list()$alignment$bin_bwa,
    bin_picard=build_default_tool_binary_list()$bin_picard,
    bin_bedtools=build_default_tool_binary_list()$bin_bedtools,
    bin_samtools=build_default_binary_list()$alignment$bin_samtool,
    bin_fastqc=build_default_tool_binary_list()$bin_fastqc,
    ref_genome=build_default_reference_list()$HG19$reference$genome,
    dbsnp=build_default_reference_list()$HG19$database$all_common,
    bi=build_default_reference_list()$HG19$panel$PCF_V3$intervals$bi,
    ti=build_default_reference_list()$HG19$panel$PCF_V3$intervals$ti,
    fastq_r1=NULL,
    fastq_r2=NULL,
    project_id=NULL,
    patient_id=NULL,
    sample_id=NULL,
    sequencing_type=NULL,
    method_type=NULL,
    method_version=NULL,
    reference=NULL,
    library_id=NULL,
    run_id=NULL,
    flowcell_id=NULL,
    lane_id=NULL,
    chromosomes=NULL,
    clean_tmp=TRUE,
    ...
){
    
      run_main=function(
            .env
      ){

        # Copy environment variables from parent scope to current environment
        .this.env=environment()
        append_env(to=.this.env,from=.env)

        # Validate all required metadata parameters are provided
        # These parameters define the sample hierarchy and processing context
        for(id in c("project_id",
                    "patient_id",
                    "sample_id",
                    "sequencing_type",
                    "method_type",
                    "method_version",
                    "reference",
                    "library_id")){
                if(is.null(get(id))){
                        stop("Variable ",id, " required to continue. Please assign a value")}
        }
        
        # Organize input FASTQ files into named list structure
        # Supports both paired-end (R1/R2) and single-end (fastq_r1 only) sequencing
        fastq=list(fastq_r1=fastq_r1,fastq_r2=fastq_r2)
        # Use sample_id as the primary identifier for output file naming
        input_id=sample_id
        
        # Extract sequencing metadata from FASTQ file headers
        # Infers instrument, run, flowcell, and lane information when available
        info=new_check_seq_info(fastq=fastq)
        
        # Auto-populate run metadata from FASTQ headers if not explicitly provided
        # This enables flexible parameter specification - required fields can be auto-detected
        for(id in c("library_id","run_id","flowcell_id","lane_id")){
            assign(id,ifelse(is.null(get(id)),
                    (info %>% dplyr::filter(name==id))$r1,get(id))
            )
        }
        
        # Construct hierarchical output directory structure based on sample metadata
        # Format: project/patient/sample/sequencing_type/method_type/method_version/reference/library/run/flowcell/lane
        
        project_structure=paste0(
                project_id,"/",
                patient_id,"/",
                sample_id,"/",
                sequencing_type,"/",
                method_type,"/",
                method_version,"/",
                reference,"/",
                library_id,"/",
                run_id,"/",
                flowcell_id,"/",
                lane_id
        )
        
        # Set env directory
        set_env_dirs(.this.env,name=project_structure)

        # Initialize main job structure and set primary execution environment
        set_main(.env=.this.env)

        # Store this function's environment in main steps registry
        .main$steps[[fn_id]]<-.this.env
        # Reference the current processing step for appending results
        .main.step=.main$steps[[fn_id]]

            # Record pipeline start time for elapsed time tracking
            start_time <- Sys.time()
        

            # Define UMI processing pipeline steps in logical execution order
            # Each step name directly corresponds to conditional processing blocks below
            # Steps are designed to handle: raw reads → UMI extraction → trimming → mapping → deduplication → consensus → remapping
            steps=c(
                "pre_trim_fastqc",
                "raw_fastq_to_bam",        # Step 1: Convert raw FASTQ to unmapped BAM format
                "extract_umi",              # Step 2: Extract molecular barcodes (UMI) from reads
                "raw_bam_to_fastq",         # Step 3: Convert UMI-tagged BAM back to FASTQ for processing
                "trim_adapt",
                "post_trim_fastqc",               # Step 4: Trim sequencing adapters and low-quality bases with fastp
                "map_trimmed",              # Step 5: Align trimmed reads to reference genome with BWA
                "tag_trimmed",              # Step 6: Merge mapped/unmapped BAM files and tag with attributes
                "filter_paired",            # Step 7: Filter for properly paired reads (flag -f 2)
                "sort_filtered",            # Step 8: Sort and index filtered BAM file by coordinate
                "pre_dedup_qc",             # Step 9: Generate QC metrics before deduplication
                "group_umi",                # Step 10: Group reads by UMI/molecular barcode for deduplication
                "collapse_consensus",       # Step 11: Generate consensus sequences from UMI-grouped reads
                "consensus_bam_to_fastq",           # Step 16: Perform BQSR on consensus BAM
                "post_dedup_qc",
                "post_dedup_fastqc",   # Step 12: Convert consensus BAM to FASTQ for remapping
                "remap_consensus",          # Step 13: Realign consensus sequences to reference genome
                "tag_consensus",            # Step 14: Merge and tag final consensus BAM with read group info
                "index_consensus",          # Step 15: Index tagged consensus BAM file
                "recal_bam"             # Step 17: Generate QC metrics after deduplication
            )

            # Append cleaning step if required
            
            if(clean_tmp){
                steps=append(steps,"clean_tmp")
            }
        
            # Total number of pipeline steps for progress reporting and loop control
            total_steps=length(steps)

            # Execute each pipeline step in order with progress logging
            for(step in 1:total_steps){
                
                # Log pipeline progress with current step number and name
                logger(paste("Running step", step, "of", total_steps, ":", steps[step]),start_time)
                
                # Wrap step execution in error handling to enable graceful failure reporting
                # If tryCatch catches error, it logs the step and error message, then aborts
                tryCatch({



                
                ### STEP 1: Convert raw FASTQ to unmapped BAM format
                if(steps[step]=="pre_trim_fastqc"){
                    
                    .main.step$steps <-append(
                        .main.step$steps,
                        new_qc_fastqc(
                                bin_fastqc=bin_fastqc,
                                fastq=list(fastq),
                                output_dir=paste0(out_file_dir,"/fastqc/pre_trim"),
                                output_name=paste0(input_id),
                                tmp_dir=tmp_dir,
                                env_dir=env_dir,
                                batch_dir=batch_dir,
                                err_msg=err_msg,
                                verbose=verbose,
                                threads=threads,
                                fn_id="pre_trim",
                                ram=ram,
                                executor_id=task_id
                        )
                    )

                  
                    .this.step=.main.step$steps$new_qc_fastqc.pre_trim
                    .main.step$out_files$fastqc$pre_trim=.this.step$out_files
                }


                ### STEP 1: Convert raw FASTQ to unmapped BAM format
                if(steps[step]=="raw_fastq_to_bam"){
                    
                    .main.step$steps <-append(
                        .main.step$steps,
                        fastq_to_sam_gatk(
                                sif_gatk=sif_gatk,
                                fastq=list(fastq),
                                tags=list(
                                    id_tag=patient_id,
                                    pu_tag="TPU",
                                    pl_tag="ILLUMINA",
                                    lb_tag=library_id,
                                    sm_tag=input_id
                                ),
                                output_dir=tmp_dir,
                                output_name=paste0(input_id,".unmapped"),
                                tmp_dir=tmp_dir,
                                env_dir=env_dir,
                                batch_dir=batch_dir,
                                err_msg=err_msg,
                                verbose=verbose,
                                threads=threads,
                                fn_id="raw",
                                ram=ram,
                                executor_id=task_id
                        )
                    )

                    .this.step=.main.step$steps$fastq_to_sam_gatk.raw
                    .main.step$out_files$raw$bam$unmapped=.this.step$out_files
                }

                if(steps[step]=="extract_umi"){

               

                    ### STEP 2: Extract molecular barcodes (UMI) from reads
                    .main.step$steps <-append(
                    .main.step$steps,
                    extract_umi_fgbio(
                        env_fgbio = env_fgbio,
                        bam=.main.step$out_files$raw$bam$unmapped,
                        output_dir=tmp_dir,
                        output_name=paste0(input_id,".unmapped"),
                        tmp_dir=tmp_dir,
                        env_dir=env_dir,
                        batch_dir=batch_dir,
                        err_msg=err_msg,
                        verbose=verbose,
                        threads=threads,
                        ram=ram,
                        executor_id=task_id
                        )
                    )

                    .this.step=.main.step$steps$extract_umi_fgbio
                    .main.step$out_files$raw$bam$unmapped$umi=.this.step$out_files
                }

                ### STEP 3: Convert UMI-tagged BAM back to FASTQ for processing
                if(steps[step]=="raw_bam_to_fastq"){

                    .main.step$steps <-append(
                        .main.step$steps,
                        sam_to_fastq_gatk(
                                sif_gatk=sif_gatk,
                                bam=.main.step$out_files$raw$bam$unmapped$umi,
                                output_dir=tmp_dir,
                                output_name=paste0(input_id,".unmapped.umi"),
                                tmp_dir=tmp_dir,
                                env_dir=env_dir,
                                batch_dir=batch_dir,
                                err_msg=err_msg,
                                verbose=verbose,
                                threads=threads,
                                fn_id="raw",
                                ram=ram,
                                executor_id=task_id
                        )
                    )

                    .this.step=.main.step$steps$sam_to_fastq_gatk.raw
                    .main.step$out_files$raw$fastq$untrimmed=.this.step$out_files
                }

                ### STEP 4: Trim sequencing adapters and low-quality bases with fastp
                if(steps[step]=="trim_adapt"){
                
                    .main.step$steps <-append(
                            .main.step$steps,
                        trim_umi_fastp(
                                env_fastp=env_fastp,
                                fastq=list(.main.step$out_files$raw$fastq$untrimmed),
                                output_dir=paste0(out_file_dir,"/fastp"),
                                output_name=paste0(input_id,".unmapped.umi"),
                                tmp_dir=tmp_dir,
                                env_dir=env_dir,
                                batch_dir=batch_dir,
                                err_msg=err_msg,
                                verbose=verbose,
                                threads=threads,
                                ram=ram,
                                executor_id=task_id
                        )
                    )

                    .this.step=.main.step$steps$trim_umi_fastp
                    .main.step$out_files$raw$fastq$trimmed=.this.step$out_files
                }

                           ### STEP 1: Convert raw FASTQ to unmapped BAM format
                if(steps[step]=="post_trim_fastqc"){
                    
                    .main.step$steps <-append(
                        .main.step$steps,
                        new_qc_fastqc(
                                bin_fastqc=bin_fastqc,
                                fastq=list(.main.step$out_files$raw$fastq$trimmed),
                                output_dir=paste0(out_file_dir,"/fastqc/post_trim"),
                                output_name=paste0(input_id),
                                tmp_dir=tmp_dir,
                                env_dir=env_dir,
                                batch_dir=batch_dir,
                                err_msg=err_msg,
                                verbose=verbose,
                                threads=threads,
                                fn_id="post_trim",
                                ram=ram,
                                executor_id=task_id
                        )
                    )

                    .this.step=.main.step$steps$new_qc_fastqc.post_trim
                    .main.step$out_files$fastqc$post_trim=.this.step$out_files
                }






                ### STEP 5: Align trimmed reads to reference genome with BWA
                if(steps[step]=="map_trimmed"){
                
                    .main.step$steps <-append(
                            .main.step$steps,
                            new_alignment_bwa(
                                    bin_bwa=bin_bwa,
                                    bin_samtools=bin_samtools,
                                    ref_genome=ref_genome,
                                    fastq=list(.main.step$out_files$raw$fastq$trimmed),
                                    tags= list(
                                        id_tag=patient_id,
                                        pu_tag="TPU",
                                        pl_tag="ILLUMINA",
                                        lb_tag=library_id,
                                        sm_tag=input_id
                                    ),
                                    output_dir=tmp_dir,
                                    output_name=paste0(input_id,".mapped.umi"),
                                    tmp_dir=tmp_dir,
                                    env_dir=env_dir,
                                    batch_dir=batch_dir,
                                    err_msg=err_msg,
                                    verbose=verbose,
                                    threads=threads,
                                    ram=ram,
                                    fn_id="raw",
                                    executor_id=task_id
                            )
                    )

                    .this.step=.main.step$steps$new_alignment_bwa.raw
                    .main.step$out_files$raw$bam$mapped$untagged=.this.step$out_files$bam

                }

                ### STEP 6: Merge mapped/unmapped BAM files and tag with attributes
                if(steps[step]=="tag_trimmed"){
                
                    .main.step$steps <-append(
                            .main.step$steps,
                            merge_bam_umi_gatk(
                                    sif_gatk=sif_gatk,
                                    bin_samtools = bin_samtools,
                                    ref_genome=ref_genome,
                                    bam=list(list(
                                        mapped=.main.step$out_files$raw$bam$mapped$untagged,
                                        unmapped=.main.step$out_files$raw$bam$unmapped$umi)),
                                    attributes=c("XO","NM","MD"),
                                    sort_order="queryname",
                                    aligned_reads_only=TRUE,
                                    add_mate_cigar=FALSE,
                                    output_dir=tmp_dir,
                                    output_name=paste0(input_id,".mapped.umi"),
                                    tmp_dir=tmp_dir,
                                    env_dir=env_dir,
                                    batch_dir=batch_dir,
                                    err_msg=err_msg,
                                    verbose=verbose,
                                    threads=threads,
                                    ram=ram,
                                    fn_id="raw",
                                    executor_id=task_id
                            )
                    )

                    .this.step=.main.step$steps$merge_bam_umi_gatk.raw
                    .main.step$out_files$raw$bam$mapped$tagged$raw=.this.step$out_files$bam

                }


                ### STEP 7: Filter for properly paired reads (flag -f 2)
                if(steps[step]=="filter_paired"){

                    .main.step$steps <-append(
                        .main.step$steps,
                        filter_samtools(
                                bin_samtools=bin_samtools,
                                bam=.main.step$out_files$raw$bam$mapped$tagged$raw,
                                flag=2,
                                chromosomes=chromosomes,
                                output_dir=tmp_dir,
                                output_name=paste0(input_id,".mapped.umi.tagged"),
                                tmp_dir=tmp_dir,
                                env_dir=env_dir,
                                batch_dir=batch_dir,
                                err_msg=err_msg,
                                verbose=verbose,
                                threads=threads,
                                ram=ram,
                                executor_id=task_id
                        )
                    )
                    
                    .this.step=.main.step$steps$filter_samtools
                    .main.step$out_files$raw$bam$mapped$tagged$filtered$ungrouped$unsorted=.this.step$out_files$bam

                }

                ### STEP 8: Sort and index filtered BAM file by coordinate
                if(steps[step]=="sort_filtered"){
                     .main.step$steps <-append(
                        .main.step$steps,
                            new_sort_and_index_bam_samtools(
                                    bin_samtools=bin_samtools,
                                    bam=.main.step$out_files$raw$bam$mapped$tagged$filtered$ungrouped$unsorted,
                                    sort=TRUE,
                                    index=TRUE,
                                    coord_sort=TRUE,
                                    stats=TRUE,
                                    output_dir=paste0(out_file_dir,"/raw"),
                                    output_name=paste0(input_id,".mapped.umi.tagged.filtered"),
                                    tmp_dir=tmp_dir,
                                    env_dir=env_dir,
                                    batch_dir=batch_dir,
                                    err_msg=err_msg,
                                    verbose=verbose,
                                    threads=threads,
                                    ram=ram,
                                    fn_id="pre",
                                    executor_id=task_id
                            )
                     )

                    .this.step=.main.step$steps$new_sort_and_index_bam_samtools.pre
                    .main.step$out_files$raw$bam$mapped$tagged$filtered$ungrouped$sorted=.this.step$out_files

                }


                ### STEP 9: Generate QC metrics before deduplication
                if(steps[step]=="pre_dedup_qc"){
                     .main.step$steps <-append(
                        .main.step$steps,
                            new_metrics_alignqc(
                                    bin_samtools=bin_samtools,
                                    bin_picard=bin_picard,
                                    bin_bedtools=bin_bedtools,
                                    ref_genome=ref_genome,
                                    bi=bi,
                                    ti=ti,
                                    bam=.main.step$out_files$raw$bam$mapped$tagged$filtered$ungrouped$sorted$srt_bam,
                                    mapq=0,
                                    method=tolower(method_type),
                                    output_dir=paste0(out_file_dir,"/alignqc/pre_dedup"),
                                    output_name=paste0(input_id),
                                    tmp_dir=tmp_dir,
                                    env_dir=env_dir,
                                    batch_dir=batch_dir,
                                    err_msg=err_msg,
                                    verbose=verbose,
                                    threads=threads,
                                    fn_id="raw",
                                    ram=ram,
                                    executor_id=task_id
                            )
                     )

                    .this.step=.main.step$steps$new_metrics_alignqc.raw
                    .main.step$out_files$raw$alignqc=.this.step$out_files

                }

                ### STEP 10: Group reads by UMI/molecular barcode for deduplication
                if(steps[step]=="group_umi"){
                    
                    .main.step$steps <-append(
                        .main.step$steps,
                        group_by_umi_fgbio(
                                env_fgbio=env_fgbio,
                                bam=.main.step$out_files$raw$bam$mapped$tagged$filtered$ungrouped$unsorted,
                                output_dir=tmp_dir,
                                output_name=paste0(input_id,".mapped.umi.tagged.filtered"),
                                tmp_dir=tmp_dir,
                                env_dir=env_dir,
                                batch_dir=batch_dir,
                                err_msg=err_msg,
                                verbose=verbose,
                                threads=threads,
                                ram=ram,
                                executor_id=task_id
                        )
                    )
                    
                    .this.step=.main.step$steps$group_by_umi_fgbio
                    .main.step$out_files$raw$bam$mapped$tagged$filtered$grouped=.this.step$out_files
                }
                ### STEP 11: Generate consensus sequences from UMI-grouped reads
                if(steps[step]=="collapse_consensus"){
                    
                    .main.step$steps <-append(
                    .main.step$steps,
                        call_consensus_fgbio(
                                env_fgbio=env_fgbio,
                                bin_samtools = bin_samtools,
                                bam=.main.step$out_files$raw$bam$mapped$tagged$filtered$grouped$bam,
                                tags= list(
                                        id_tag=patient_id,
                                        pu_tag="TPU",
                                        pl_tag="ILLUMINA",
                                        lb_tag=library_id,
                                        sm_tag=input_id
                                    ),
                                output_dir=tmp_dir,
                                output_name=input_id,
                                tmp_dir=tmp_dir,
                                env_dir=env_dir,
                                batch_dir=batch_dir,
                                err_msg=err_msg,
                                verbose=verbose,
                                threads=threads,
                                ram=ram,
                                executor_id=task_id
                        )
                    )
                    
                    .this.step=.main.step$steps$call_consensus_fgbio
                    .main.step$out_files$consensus$bam$unmapped=.this.step$out_files
                }


                ### STEP 12: Convert consensus BAM to FASTQ for remapping
                if(steps[step]=="consensus_bam_to_fastq"){
                
                    .main.step$steps <-append(
                        .main.step$steps,
                    sam_to_fastq_gatk(
                                sif_gatk=sif_gatk,
                                bam= .main.step$out_files$consensus$bam$unmapped$bam,
                                output_dir=tmp_dir,
                                output_name=paste0(input_id,".consensus"),
                                tmp_dir=tmp_dir,
                                env_dir=env_dir,
                                batch_dir=batch_dir,
                                err_msg=err_msg,
                                verbose=verbose,
                                threads=threads,
                                fn_id="consensus",
                                ram=ram,
                                executor_id=task_id
                        )
                    )

                    .this.step=.main.step$steps$sam_to_fastq_gatk.consensus
                    .main.step$out_files$consensus$fastq=.this.step$out_files
                }


                if(steps[step]=="post_dedup_fastqc"){
                    
                    .main.step$steps <-append(
                        .main.step$steps,
                        new_qc_fastqc(
                                bin_fastqc=bin_fastqc,
                                fastq=list(.main.step$out_files$consensus$fastq),
                                output_dir=paste0(out_file_dir,"/fastqc/post_dedup"),
                                output_name=paste0(input_id),
                                tmp_dir=tmp_dir,
                                env_dir=env_dir,
                                batch_dir=batch_dir,
                                err_msg=err_msg,
                                verbose=verbose,
                                threads=threads,
                                fn_id="post_dedup",
                                ram=ram,
                                executor_id=task_id
                        )
                    )

                  
                    .this.step=.main.step$steps$new_qc_fastqc.pre_trim
                    .main.step$out_files$fastqc$pre_trim=.this.step$out_files
                }

                ### STEP 13: Realign consensus sequences to reference genome
                if(steps[step]=="remap_consensus"){

                    .main.step$steps <-append(
                    .main.step$steps,
                        new_alignment_bwa(
                                bin_bwa=bin_bwa,
                                bin_samtools=bin_samtools,
                                ref_genome=ref_genome,
                                fastq=list(.main.step$out_files$consensus$fastq),
                                tags=list(
                                    id_tag=patient_id,
                                    pu_tag="TPU",
                                    pl_tag="ILLUMINA",
                                    lb_tag=library_id,
                                    sm_tag=input_id
                                ),
                                output_dir=tmp_dir,
                                output_name=paste0(input_id,".consensus.mapped.untagged"),
                                tmp_dir=tmp_dir,
                                env_dir=env_dir,
                                batch_dir=batch_dir,
                                err_msg=err_msg,
                                verbose=verbose,
                                threads=threads,
                                ram=ram,
                                fn_id="consensus",
                                executor_id=task_id
                                )
                        )

                    .this.step=.main.step$steps$new_alignment_bwa.consensus
                    .main.step$out_files$consensus$bam$mapped$untagged=.this.step$out_files
                }
        

                ### STEP 14: Merge and tag final consensus BAM with read group info
                if(steps[step]=="tag_consensus"){

                    .main.step$steps <-append(
                    .main.step$steps,
                    merge_bam_umi_gatk(
                            sif_gatk=sif_gatk,
                            bin_samtools = bin_samtools,
                            ref_genome=ref_genome,
                            bam=list(list(
                                mapped=.main.step$out_files$consensus$bam$mapped$untagged,
                                unmapped=.main.step$out_files$consensus$bam$unmapped)),
                            attributes=c("X0","RX"),
                            sort_order="coordinate",
                            aligned_reads_only=FALSE,
                            add_mate_cigar=TRUE,
                            output_dir=tmp_dir,
                            output_name=paste0(input_id,".consensus.mapped.tagged"),
                            tmp_dir=tmp_dir,
                            env_dir=env_dir,
                            batch_dir=batch_dir,
                            err_msg=err_msg,
                            verbose=verbose,
                            threads=threads,
                            ram=ram,
                            fn_id="consensus",
                            executor_id=task_id
                            )
                    )

                    .this.step=.main.step$steps$merge_bam_umi_gatk.consensus
                    .main.step$out_files$consensus$bam$mapped$tagged$bam=.this.step$out_files
                }

                

                ### STEP 15: Index tagged consensus BAM file
                if(steps[step]=="index_consensus"){
                     .main.step$steps <-append(
                        .main.step$steps,
                            new_sort_and_index_bam_samtools(
                                    bin_samtools=bin_samtools,
                                    bam=.main.step$out_files$consensus$bam$mapped$tagged$bam,
                                    sort=FALSE,
                                    index=TRUE,
                                    stats=FALSE,
                                    tmp_dir=tmp_dir,
                                    env_dir=env_dir,
                                    batch_dir=batch_dir,
                                    err_msg=err_msg,
                                    verbose=verbose,
                                    threads=threads,
                                    ram=ram,
                                    fn_id="post",
                                    executor_id=task_id
                            )
                     )

                    .this.step=.main.step$steps$new_sort_and_index_bam_samtools.post
                    .main.step$out_files$consensus$bam$mapped$tagged$index=.this.step$out_files

                }


                ### This step will fail if there is not enough reads to recalibrate
                ### STEP 16: Perform BQSR on consensus BAM
                if(steps[step]=="recal_bam"){
                     .main.step$steps <-append(
                        .main.step$steps,
                            new_recal_gatk(
                                    bin_samtools=bin_samtools,
                                    sif_gatk=sif_gatk,
                                    bin_picard=bin_picard,
                                    ref_genome=ref_genome,
                                    dbsnp=dbsnp,
                                    chromosomes=chromosomes,
                                    bam=.main.step$out_files$consensus$bam$mapped$tagged$bam,
                                    output_dir=paste0(out_file_dir),
                                    output_name=input_id,
                                    tmp_dir=tmp_dir,
                                    env_dir=env_dir,
                                    batch_dir=batch_dir,
                                    err_msg=err_msg,
                                    verbose=verbose,
                                    threads=threads,
                                    ram=ram,
                                    executor_id=task_id
                            )
                     )

                    .this.step=.main.step$steps$new_recal_gatk
                    .main.step$out_files$recal=.this.step$out_files$recal

                }


                ### STEP 17: Generate QC metrics after deduplication
                if(steps[step]=="post_dedup_qc"){
                     .main.step$steps <-append(
                        .main.step$steps,
                            new_metrics_alignqc(
                                    bin_samtools=bin_samtools,
                                    bin_picard=bin_picard,
                                    bin_bedtools=bin_bedtools,
                                    ref_genome=ref_genome,
                                    bi=bi,
                                    ti=ti,
                                    bam=.main.step$out_files$recal$sorted$srt_bam,
                                    mapq=0,
                                    method=tolower(method_type),
                                    output_dir=paste0(out_file_dir,"/alignqc/post_dedup"),
                                    output_name=paste0(input_id),
                                    tmp_dir=tmp_dir,
                                    env_dir=env_dir,
                                    batch_dir=batch_dir,
                                    err_msg=err_msg,
                                    verbose=verbose,
                                    threads=threads,
                                    fn_id="consensus",
                                    ram=ram,
                                    executor_id=task_id
                            )
                     )

                    .this.step=.main.step$steps$new_metrics_alignqc.raw
                    .main.step$out_files$raw$alignqc=.this.step$out_files
                }


                


                if(steps[step]=="clean_tmp"){
                    unlink(tmp_dir,recursive = TRUE,force=TRUE)
                }


                    # Log successful step completion
                    logger(paste("Completed step", step, "of", total_steps, ":", steps[step]),start_time)
                    
                }, error=function(e){
                    # Handle step execution errors with informative message
                    logger(paste("ERROR in step", step, ":", steps[step]),start_time)
                    stop(paste("Step '" , steps[step], "' failed. Error:", e$message,
                              "\nReview input files and parameters before retrying."))
                })
        
        }
            
        # Log pipeline completion with total runtime
        total_elapsed <- as.numeric(difftime(Sys.time(), start_time, units="secs"))
        total_elapsed_str <- sprintf("%.1f", total_elapsed)
        logger(paste("UMI processing pipeline completed successfully."),start_time)
        logger(paste("Total steps executed:", total_steps, "| Total runtime:", total_elapsed_str, "seconds"),start_time)
          
        # Return main object to parent environment for job tracking and output reporting
        .env$.main <- .main

    }
    
    # Setup execution environment by copying parent scope variables
    .base.env=environment()
    # Merge additional parameters passed via ... into execution environment
    list2env(list(...),envir=.base.env)
    # Configure environment variables required for batch processing and temp directories
    set_env_vars(
        .env= .base.env,
        vars="fastq_r1"
    )

    # Launch the UMI processing pipeline with fully prepared environment
    launch(.env=.base.env)
        

}


        
       



