

#' Preprocess sequencing data
#' 
#'
#' @param sample_sheet Input sample sheet
#' @param config Default tool configure
#' @param executor_id Task EXECUTOR ID . Default "preprocessSEQ"
#' @param task_name Task name . Default "preprocessSEQ"
#' @param output_dir Path to output directory
#' @export


preprocess_seq=function(sample_sheet=build_default_sample_sheet(),
    executor_id=make_unique_id("preprocessSEQ"), 
    variables=build_default_variable_list(),config=build_default_config(),
    task_name="preprocessSEQ",output_dir=""){
        
    task_id=make_unique_id(task_name)
    sample_info=list()
    sample_info$seq_info=seq_info_check(sample_info=sample_sheet)
    sample_info$tool_config=parameter_config_check(sample_sheet=sample_sheet,config=config)
    seq_info=sample_info$seq_info %>% dplyr::filter(validate==TRUE)

    job=build_job(executor_id=executor_id,task_id=task_id)


 
    ## Go through each patient
    smt=lapply(unique(seq_info$project_id),FUN=function(project_id){
        seq_info_project=seq_info %>% dplyr::filter(project_id==project_id)
        out_file_dir_project=set_dir(dir=output_dir,name=project_id)
        ## Go through each patient
        lapply(unique(seq_info_project$patient_id),FUN=function(patient_id){
            seq_info_patient=seq_info_project %>% dplyr::filter(project_id==project_id,patient_id==patient_id)
            out_file_dir_patient=set_dir(dir=out_file_dir_project,name=patient_id)
            ## Go through each sample
            lapply(unique(seq_info_patient$sample_id),FUN=function(sample_id){
                seq_info_sample=seq_info_patient %>% dplyr::filter(project_id==project_id,patient_id==patient_id,sample_id==sample_id)
                out_file_dir_sample=set_dir(dir=out_file_dir_patient,name=sample_id)
                    ## Go through each method
                    lapply(unique(seq_info_sample$method_id),FUN=function(method_id){
                        seq_info_method=seq_info_sample %>% dplyr::filter(project_id==project_id,patient_id==patient_id,
                        sample_id==sample_id,method_id==method_id)
                        out_file_dir_method=set_dir(dir=out_file_dir_sample,name=method_id)
                     
                        ## Go through each flowcell ID
                        lapply(unique(seq_info_method$flowcell_id),FUN=function(flowcell_id){
                            
                            seq_info_flowcell=seq_info_method %>% dplyr::filter(project_id==project_id,patient_id==patient_id,
                            sample_id==sample_id,method_id==method_id,flowcell_id==flowcell_id)
                            out_file_dir_flowcell=set_dir(dir=out_file_dir_method,name=flowcell_id)
                         
                            ## Go through each lane
                            lapply(unique(seq_info_flowcell$lane_id),FUN=function(lane_id){
                                
                                seq_info_lane=seq_info_flowcell %>% dplyr::filter(project_id==project_id,patient_id==patient_id,
                                sample_id==sample_id,method_id==method_id,flowcell_id==flowcell_id,lane_id==lane_id)
                                out_file_dir_lane=set_dir(dir=out_file_dir_flowcell,name=lane_id)
                        
                                ## Go through each library
                                lapply(unique(seq_info_lane$library_id),FUN=function(library_id){
                                        
                                        seq_info_library=seq_info_lane %>% dplyr::filter(project_id==project_id,patient_id==patient_id,
                                        sample_id==sample_id,method_id==method_id,flowcell_id==flowcell_id,lane_id==lane_id,library_id=library_id)
                                        out_file_dir_library=set_dir(dir=out_file_dir_lane,name=library_id)

                                        cat(paste0("Project ID: ",project_id,"\n"))
                                        cat(paste0(add_fill(n=1),add_bl(),"Patient ID: ",patient_id,"\n"))                                        
                                        cat(paste0(add_fill(n=2),add_bl(),"Sample ID: ",sample_id,"\n"))
                                        cat(paste0(add_fill(n=3),add_bl(),"Method ID: ",method_id,"\n"))
                                        cat(paste0(add_fill(n=4),add_bl(),"Flowcell ID: ",flowcell_id,"\n"))
                                        cat(paste0(add_fill(n=5),add_bl(),"Lane ID: ",lane_id,"\n"))
                                        cat(paste0(add_fill(n=6),add_bl(),"Library ID: ",library_id,"\n"))

                               
                                        # seq_info_R1=seq_info %>% dplyr::filter(project_id==project_id,
                                        # patient_id==patient_id,sample_id==sample_id,method_id==method_id,
                                        # flowcell_id==flowcell_id,lane_id==lane_id,library_id==library_id,
                                        # read_group=="R1")
                
                                        # seq_info_R2=seq_info %>% dplyr::filter(project_id==project_id,
                                        # patient_id==patient_id,sample_id==sample_id,method_id==method_id,
                                        # flowcell_id==flowcell_id,lane_id==lane_id,library_id==library_id,
                                        # read_group=="R2")
                                        
                                        # cat(paste0(add_fill(n=7),"|----R1: ",seq_info_R1$path,"\n"))
                                        # cat(paste0(add_fill(n=7),"|      \n"))
                                        # cat(paste0(add_fill(n=7),"|----R2: ",seq_info_R2$path,"\n"))
                                })
                            })
                        })
                    })
                })
            })
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

#                 if(grepl("post_fastqc",rownames(parameters))){

#                     job_report=qc_fastqc(bin_path=bin_fastqc,
#                     file_R1=sub_sub_sample_info$file[1],
#                     file_R2=sub_sub_sample_info$file[2],
#                     output_dir=paste0(out_file_dir,"/fastqc_reports/post_trim"),
#                     executor_id=task_id,
#                     verbose=parameters["pre_fastqc","verbose"],
#                     mode=parameters["pre_fastqc","mode"],
#                     threads=parameters["pre_fastqc","threads"],
#                     ram=parameters["pre_fastqc","ram"],
#                     time=parameters["pre_fastqc","time"],
#                     update_time=60,wait=FALSE,hold=job_report$job_id)
#                 }
        
        
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