

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
#' @export


preprocess_seq=function(sample_sheet=build_default_sample_sheet(),
    executor_id=make_unique_id("preprocessSEQ"), 
    vars_list=build_default_variable_list(),
    config=build_default_config(),
    opts_list=build_default_option_list(),
    pmts_list=build_default_parameter_list(),
    steps_list=build_default_steps_list(),
    bin_list=build_default_binaries_list(),
    task_name="preprocessSEQ",output_dir="",
    merge_level="library",nest_ws=1,
    nesting=""){

    task_id=make_unique_id(task_name)
    validate_sample_sheet(sample_sheet=sample_sheet,
    vars_list=vars_list,opts_list=opts_list)
    
    seq_info=seq_info_check(sample_sheet=sample_sheet,vars_list=vars_list)

    seq_info=dplyr::left_join(seq_info,
    parameter_config_check(sample_sheet=sample_sheet,config=config,
    vars_list=vars_list,steps_list=steps_list))
    job=build_job(executor_id=executor_id,task_id=task_id)
    for_id(seq_info=seq_info,output_dir=output_dir,
    vars_list=vars_list,nesting=nesting,merge_level=merge_level,
    pmts_list=pmts_list,bin_list=bin_list,print_tree=TRUE)

}




#' For each variable in sample sheet 
#' 
#' Run through each variable in sample sheet and report information
#'
#' @param seq_info Sample sequencing information
#' @param var_list List with variables
#' @param pmts_list List with parameters
#' @param nesting Starting nesting level
#' @param nest_ws Nesting ws to provide. Default 1
#' @param merge_level Level to merge samples. Default library
#' @param executor_id Task EXECUTOR ID . Default "preprocessSEQ"
#' @param task_name Task name . Default "preprocessSEQ"
#' @param output_dir Path to output directory. Default none
#' @export



for_id=function(seq_info,output_dir="",
             vars_list=build_default_variable_list(),
             pmts_list=build_default_parameter_list(),
             bin_list=build_default_binaries_list(),
             nesting="",merge_level="library",
             nest_ws=1,print_tree=TRUE){  
                var=vars_list$variable[1]
                var_text=vars_list$text[1]
                vars_list_left=vars_list[-1,]
                info=unique(seq_info[,var,drop=TRUE])
                count=1
                merge=FALSE
                rs=lapply(X=info,FUN=function(id){
                

                        ## Filter sequencing info for id
                        seq_info_id=seq_info[seq_info[,var,drop=TRUE]==id,]
                        out_file_dir=set_dir(dir=output_dir,name=id)
                
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
                                seq_info_id[seq_info_id$last!=TRUE&seq_info_id$order>4,]$step="FALSE"
                            }
                
                            instrument_name=""
                            if(var=="flowcell_id"){
                                instrument_name=paste0("   Platform: ",unique(seq_info_id$instrument_by_flowcell_id))
                            }

                            cat(paste0(nesting,"|----",crayon::blue(var_text),crayon::red(id),
                            crayon::silver(instrument_name),merge_txt,"\n"))
                        

                            nesting=break_nest(count=count,info=info,nesting=nesting)


                        }

                        

                        ## Call recursively if variables
                        if(length(vars_list_left$variable)>0){
                            for_id(seq_info=seq_info_id,output_dir=out_file_dir,vars_list=vars_list_left,
                            nesting=nesting,nest_ws=nest_ws)
                        }else{
                            tool_config_id=seq_info_id %>% dplyr::select(-c(read_group,path)) %>%  
                            dplyr::distinct() %>% dplyr::filter(step==TRUE)

                            seq_info_id=seq_info_id %>% dplyr::select(-c("order",pmts_list$parameter)) %>%  
                            dplyr::distinct()

                            seq_info_R1=seq_info_id[seq_info_id$read_group=="R1",]
                            seq_info_R2=seq_info_id[seq_info_id$read_group=="R2",]
                            
                            if(print_tree){
                                cat(add_nesting_ws(nesting,n=nest_ws))
                                cat(paste0(nesting,"|----",crayon::green(paste0("R1: ",seq_info_R1$path)),"\n"))
                                cat(add_nesting_ws(nesting,n=nest_ws))
                                cat(paste0(nesting,"|----",crayon::green(paste0("R2: ",seq_info_R2$path)),"\n")) 
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
                             }
                            }
    count<<-count+1
    })
}



            if(tool_config_id[step,]$name=="pre_fastqc"){
                                sink(paste0(out_file_dir,"fastqc.out"))
                                job_report=qc_fastqc(bin_path=bin_list$pre_fastqc$bin_fastqc,
                                    file_R1=seq_info_R1$path,
                                    file_R2=seq_info_R2$path,
                                    output_dir=paste0(out_file_dir,"/fastqc_reports/pre_trim"),
                                    executor_id=task_id,
                                    verbose=tool_config_id[step,]$verbose,
                                    mode=tool_config_id[step,]$mode,
                                    threads=tool_config_id[step,]$threads,
                                    ram=tool_config_id[step,]$ram,
                                    time=tool_config_id[step,]$time,
                                    update_time=60,wait=FALSE,hold="")
                            }   



#   ## Run only preselected steps

#             tool_config=tool_config %>% filter(step==TRUE)
#             out_file_dir=set_dir(dir=out_file_dir_main,name=paste0(lane_id))

#             ## Check which steps are selected in sample sheet

#             if(grepl("pre_fastqc",rownames(parameters))){

                job_report=qc_fastqc(bin_path=bin_fastqc,
                file_R1=sub_sub_sample_info$file[1],
                file_R2=sub_sub_sample_info$file[2],
                output_dir=paste0(out_file_dir,"/fastqc_reports/pre_trim"),
                executor_id=task_id,
                verbose=tool_parameters["pre_fastqc","verbose"],
                mode=tool_parameters["pre_fastqc","mode"],
                threads=tool_parameters["pre_fastqc","threads"],
                ram=tool_parameters["pre_fastqc","ram"],
                time=tool_parameters["pre_fastqc","time"],
                update_time=60,wait=FALSE,hold="")
        
#             }
            
#             if(grepl("trimming",rownames(parameters))){

#                 job_report=trimming_skewer(bin_path=bin_skewer,
#                 file_R1=sub_sub_sample_info$file[1],
#                 file_R2=sub_sub_sample_info$file[2],
#                 output_dir=out_file_dir,
#                 xadapt=tool_parameters["trimming","xadapt"],
#                 yadapt=tool_parameters["trimming","yadapt"],
#                 threads=tool_parameters["trimming","threads"],
#                 ram=tool_parameters["trimming","ram"],
#                 verbose=tool_parameters["trimming","verbose"],
#                 mean_quality=tool_parameters["trimming","mean_quality"],
#                 min_length=tool_parameters["trimming","min_length"],
#                 max_length=tool_parameters["trimming","min_length"],
#                 mode=tool_parameters["trimming","mode"],
#                 time=tool_parameters["trimming","time"],
#                 executor_id=task_id,
#                 update_time=60,wait=FALSE)

#             }

        
        
#                 if(grepl("alignment",rownames(parameters))){
#                     job_report=alignment_bwa(bin_path=bin_bwa,
#                         bin_path2=bin_samtools,
#                         file_R1=sub_sub_sample_info$file[1],
#                         file_R2=sub_sub_sample_info$file[2],
#                         output_dir=out_file_dir,
#                         id_tag=sample_parameters["tags","id"],
#                         pu_tag=sample_parameters["tags","pu"],
#                         pl_tag=sample_parameters["tags","platform"],
#                         lb_tag=sample_parameters["tags","library"],
#                         sm_tag=sample_parameters["tags","sample"],
#                         threads=tool_parameters["alignment","threads"],
#                         ram=tool_parameters["alignment","ram"],
#                         ref_genome=tool_parameters["alignment","ref_genome"],
#                         coord_sorted=tool_parameters["alignment","coord_sorted"],
#                         stats=tool_parameters["alignment","stats"],
#                         verbose=tool_parameters["alignment","verbose"],
#                         mode=tool_parameters["alignment","mode"],
#                         time=tool_parameters["alignment","time"],
#                         executor_id=task_id,
#                         update_time=60,
#                         wait=FALSE,
#                         hold=job_report$job_id)
#                 }
#                 return(job_report)
        
            
#             if(n_lanes>1){
#                 job_report=merge_bams_samtools(
#                     bin_path=bin_samtools,
#                     bams=job_report$out_files$bam,
#                     output_name=sample_parameters["merge_bams","verbose"],
#                     verbose=tool_parameters["merge_bams","verbose"],
#                     threads=tool_parameters["merge_bams","threads"],
#                     ram=tool_parameters["merge_bams","time"],
#                     mode=tool_parameters["merge_bams","mode"],
#                     time=tool_parameters["merge_bams","time"],
#                     executor_id=task_id,
#                     update_time=60,wait=FALSE,hold=job_report$job_id)
#             }

#             if(grepl("markdups",rownames(parameters))){
#                 job_report=markdups_gatk(
#                     bin_path=bin_markdups_gatk,
#                     bam=job_report$out_files$bam,
#                     output_dir=out_file_dir_main,
#                     verbose=tool_parameters["markdups","verbose"],
#                     tmp_dir=tool_parameters["markdups","tmp_dir"],
#                     threads=tool_parameters["markdups","threads"],
#                     ram=tool_parameters["markdups","ram"],
#                     remove_duplicates=tool_parameters["markdups","remove_duplicates"],
#                     mode=tool_parameters["markdups","mode"],
#                     time=tool_parameters["markdups","time"],
#                     executor_id=task_id,
#                     update_time=60,wait=FALSE,hold=job_report$job_id)
#             }

#             if(grepl("recalibrate",rownames(parameters))){
#                 job_report=recal_gatk(
#                     bin_path=bin_samtools,
#                     bin_path2=bin_gatk,
#                     bin_path3=bin_picard,
#                     bam=bam,
#                     output_dir=out_file_dir,
#                     ref_genome=tool_parameters["recal_gatk","ref_genome"],
#                     dbsnp=tool_parameters["recal_gatk","dbsnp"],
#                     clean=tool_parameters["recal_gatk","clean"],
#                     verbose=tool_parameters["recal_gatk","verbose"],
#                     threads=tool_parameters["recal_gatk","threads"],
#                     ram=tool_parameters["recal_gatk","ram"],
#                     mode=tool_parameters["recal_gatk","mode"],
#                     time=tool_parameters["recal_gatk","time"],
#                     executor_id=task_id,
#                     update_time=60,wait=FALSE,hold=job)
#             }

#             if(grepl("alignqc",rownames(parameters))){
#                 align_qc_metrics(bin_path=bin_samtools,
#                 bin_path2=bin_picard,bin_path3=bin_bedtools,
#                 bam=bam,output_dir="",ref_genome="",verbose=FALSE,tmp_dir=".",mapq=0,bi="",
#                 ti="",ri="",ref_flat="",method="tg",mode="local",executor=make_unique_id("alignQC"),
#                 task="alignQC",time="48:0:0",threads=4,ram=4,update_time=60,wait=FALSE, hold="")
#             }




#   if(grepl("post_fastqc",rownames(parameters))){

#                             job_report=qc_fastqc(bin_path=bin_fastqc,
#                             file_R1=seq_info_R1$path,
#                             file_R2=seq_info_R2$path,
#                             output_dir=paste0(out_file_dir,"/fastqc_reports/post_trim"),
#                             executor_id=task_id,
#                             verbose=parameters["pre_fastqc","verbose"],
#                             mode=parameters["pre_fastqc","mode"],
#                             threads=parameters["pre_fastqc","threads"],
#                             ram=parameters["pre_fastqc","ram"],
#                             time=parameters["pre_fastqc","time"],
#                             update_time=60,wait=FALSE,hold=job_report$job_id)
#                         }