

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
    vars_list=build_default_variable_list(),config=build_default_config(),
    pmts_list=build_default_parameter_list(),
    task_name="preprocessSEQ",output_dir="",merge_level="library",nest_ws=1){
    task_id=make_unique_id(task_name)
    validate_sample_sheet(sample_sheet=sample_sheet,vars_list=vars)
    sample_info=list()
    sample_info$seq_info=seq_info_check(sample_sheet=sample_sheet,vars_list=vars)
    sample_info$tool_config=parameter_config_check(sample_sheet=sample_sheet,
    config=config,vars_list=vars)
    job=build_job(executor_id=executor_id,task_id=task_id)
    for_id(seq_info=sample_info$seq_info,tool_config=sample_info$tool_config,output_dir=output_dir,vars_list=vars_list,
    nesting=nesting,merge_level=merge_level,pmts_list=pmts_list)

}




#' For each variable in sample sheet 
#' 
#' Run through each variable in sample sheet and report information
#'
#' @param seq_info Sample sequencing information
#' @param tool_config Tool config for each sample
#' @param var_list List with variables
#' @param pmts_list List with parameters
#' @param nesting Starting nesting level
#' @param nest_ws Nesting ws to provide. Default 1
#' @param merge_level Level to merge samples. Default library
#' @param executor_id Task EXECUTOR ID . Default "preprocessSEQ"
#' @param task_name Task name . Default "preprocessSEQ"
#' @param output_dir Path to output directory. Default none
#' @export



for_id=function(seq_info,tool_config, output_dir="",
             vars_list=build_default_variable_list(),
             pmts_list=build_default_parameter_list(),
             nesting="",merge_level="library",nest_ws=1){  
                
             var=vars_list$variable[1]
             var_text=vars_list$text[1]
             vars_list_left=vars_list[-1,]
             info=unique(seq_info[,var,drop=TRUE])
             count=1
             rs=lapply(X=info,FUN=function(id){
            
                    if(var!="project_id"){
                          add_nesting_ws(nesting,n=nest_ws)
                    }
                    merge=""
                    if(grepl(merge_level,var)){
                        merge=crayon::bold(" <<<<===== INFO::SAMPLES WILL BE MERGED AT THIS LEVEL")
                    }


                    
                    seq_info_id=seq_info[seq_info[,var,drop=TRUE]==id,]
                    tool_config_id=tool_config[tool_config[,var,drop=TRUE]==id,]
            
                    out_file_dir=set_dir(dir=output_dir,name=id)
                   
                    instrument_name=""
                    if(var=="flowcell_id"){
                        instrument_name=paste0("   Platform: ",unique(seq_info_id$instrument_by_flowcell_id))
                    }

                    cat(paste0(nesting,"|----",crayon::blue(var_text),crayon::red(id),crayon::silver(instrument_name),merge,"\n"))
                   

                    nesting=break_nest(count=count,info=info,nesting=nesting)

               
                    if(length(vars_left$variable)>0){
                        for_id(seq_info=seq_info_id,output_dir=out_file_dir,vars=vars_list_left,nesting=nesting,nest_ws=nest_ws)
                    }else{
                        seq_info_R1=seq_info_id[seq_info_id$read_group=="R1",]
                        seq_info_R2=seq_info_id[seq_info_id$read_group=="R2",]
                        tool_config_id=tool_config_id %>% dplyr::distinct(-c(R1,R2))
                        add_nesting_ws(nesting,n=nest_ws)
                        cat(paste0(nesting,"|----",crayon::green(paste0("R1: ",seq_info_R1$path)),"\n"))
                        add_nesting_ws(nesting,n=nest_ws)
                        cat(paste0(nesting,"|----",crayon::green(paste0("R2: ",seq_info_R2$path)),"\n"))    
                        lapply(seq(1,nrow(tool_config_id)),FUN=function(step){
                            lapply(seq(1,nrow(pmts_list$parameter)),FUN=function(pmt){
                                cat(paste0(nesting,"|----    ",paste0(pmts_list[pmt,]$text,tool_config[step,pmt]),"\n"))
                            })
                        })
                     
                    }
                    count<<-count+1
            })
}

















#   ## Run only preselected steps

#             tool_config=tool_config %>% filter(step==TRUE)
#             out_file_dir=set_dir(dir=out_file_dir_main,name=paste0(lane_id))

#             ## Check which steps are selected in sample sheet

#             if(grepl("pre_fastqc",rownames(parameters))){

#                 job_report=qc_fastqc(bin_path=bin_fastqc,
#                 file_R1=sub_sub_sample_info$file[1],
#                 file_R2=sub_sub_sample_info$file[2],
#                 output_dir=paste0(out_file_dir,"/fastqc_reports/pre_trim"),
#                 executor_id=task_id,
#                 verbose=tool_parameters["pre_fastqc","verbose"],
#                 mode=tool_parameters["pre_fastqc","mode"],
#                 threads=tool_parameters["pre_fastqc","threads"],
#                 ram=tool_parameters["pre_fastqc","ram"],
#                 time=tool_parameters["pre_fastqc","time"],
#                 update_time=60,wait=FALSE,hold="")
        
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