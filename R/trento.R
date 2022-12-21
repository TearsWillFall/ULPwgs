#' Generate a quality control (QC) report from a fastaqc file
#'
#' This function takes a set of sequence files (fastq,SAM,BAM...) and
#' returns a report in HTML format.
#'
#'
#' @param sif_path Path to singularity image file. 
#' @param fastq_dir Path to fastq directory. 
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Executor ID. Default "preprocess_trento"
#' @param task_name Name of the task. Default "preprocess_trento"
#' @param threads Number of CPU cores to use. Default 3.
#' @param ram RAM memory for batched job. Default 4
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param hold Job to hold on in batched mode.
#' @export



preprocess_seq_trento=function(
    sif_preprocess=build_default_sif_list()$sif_preprocess,
    fastq_dir="", threads=3,ram=4,output_dir=".",verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    executor_id=make_unique_id("preprocess_trento"),tmp_dir=".",
    task_name="reprocess_trento",mode="local",time="48:0:0",
    update_time=60,wait=FALSE,hold=NULL){

    argg <- as.list(environment())

    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir,name="preprocess")
    
 
    exec_code=paste(" singularity run --app preProcess -H ",paste0(getwd(),":/home"),
    sif_preprocess, " -i ", fastq_dir ," -o ", out_file_dir," -t ",tmp_dir,
    " -n " , threads, " -m ", ram, 
    "-f .+_R1_.+[.]f.+[.]gz,.+_1[.].+[.]gz,R1[.].+[.]gz -r .+_R2_.+[.]f.+[.]gz,.+_2[.].+[.]gz,R2[.].+[.]gz")
    
    
    out_file_dir2=set_dir(dir=out_file_dir,name="batch")
    job=build_job(executor_id=executor_id,task_id=make_unique_id(task_id))

    if(mode=="batch"){
        batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,threads=threads,
        output_dir=out_file_dir2)
        exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
    }

    if(verbose){
        print_verbose(job=job,arg=argg,exec_code=exec_code)
    }


    job_report=build_job_report(
        job_id=job,
        executor_id=executor_id,
        exec_code=exec_code, 
        task_id=task_id,
        input_args=argg,
        out_file_dir=out_file_dir,
        out_files=list(
            NA
        )
    )

    error=execute_job(exec_code=exec_code)
    if(error!=0){
        stop("Preprocess failed to run due to unknown error.
        Check std error for more information.")
    }
    
    if(wait&&mode=="batch"){
        batch_validator(job=job_report$job_id,
        time=update_time,verbose=verbose,threads=threads)
    }

    return(job_report)


}



#' Process multiple tumour samples using CLONET
#'
#' This function identifies a set of BAM files as tumour and normal
#' and processes them using the CLONET pipeline in parallel
#' If sample sheet is provided data has to be supplied in the following format:
#' patient_id   tumour  normal  version
#' 
#' Header can be ommitted if data is given in the order above and header argument is set to FALSE.
#' Sample sheet separator can be set using the sep argument.
#' 
#' 
#' @param sample_sheet Path to sample sheet or data.frame. 
#' @param bam_dir Path to bam directory. Only if sample sheet is not provided. 
#' @param normal_id Normal sample identifier. Only if sample sheet is not provided. 
#' @param patient_id Patient id. Only if sample sheet is not provided. 
#' @param tc Pre-computed tumour content. Default NULL.
#' @param ploidy Pre-computed ploidy. Default NULL.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Executor ID. Default "preprocess_trento"
#' @param task_name Name of the task. Default "preprocess_trento"
#' @param threads Number of CPU cores to use. Default 3.
#' @param ram RAM memory for batched job. Default 4
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param output_dir Path to the output directory.
#' @param tmp_dir Path to temporary file directory.
#' @param verbose Enables progress messages. Default False.
#' @param hold Job to hold on in batched mode.
#' @export

multisample_clonet_trento=function(
    sample_sheet=NA,bam_dir="",normal_id="",patient_id="",version="V3",
    tc=NULL,ploidy=NULL,
    tmp_dir=".",header=TRUE,sep="",threads=3,
    ram=4,output_dir=".",verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    executor_id=make_unique_id("multi_clonet"),
    task_name="multi_clonet",mode="local",time="48:0:0",
    update_time=60,wait=FALSE,hold=NULL
){

        argg <- as.list(environment())

        task_id=make_unique_id(task_name)
        out_file_dir=set_dir(dir=output_dir)

        job=build_job(executor_id=executor_id,task_id=task_id)

        job_report=build_job_report(
                job_id=job,
                executor_id=executor_id,
                exec_code=list(), 
                task_id=task_id,
                input_args=argg,
                out_file_dir=out_file_dir,
                out_files=list(
                )
        )



    columns=c("patient_id","tumour","normal","version","tc","ploidy","verbose",
    "batch_config","threads","ram","time","mode","hold")


    if(!is.na(sample_sheet)){
      
        if(!is.data.frame(sample_sheet)){
                file_info=read.csv(sample_sheet,header=header,sep=sep,stringsAsFactors=FALSE)
                if(!header){
                    names(file_info)=columns
                }
        }else{
                file_info=sample_sheet
        }
        

        job_report[["steps"]][["clonet"]]=parallel::mclapply(seq(1,nrow(file_info)),FUN=function(x){
            
            lapply(columns,FUN=function(col){
                if(is.null(file_info[[col]])){
                    file_info[[col]]<<-get(col)
                }

        
            })
        
            job_report<- clonet_trento(
                tumour=file_info[x,]$tumour,
                normal=file_info[x,]$normal,
                patient_id=file_info[x,]$patient_id,
                version=file_info[x,]$version,
                threads=file_info[x,]$threads,
                tc=file_info[x,]$tc,
                ploidy=file_info[x,]$ploidy,
                ram=file_info[x,]$ram,
                output_dir=paste0(out_file_dir,file_info[x,]$patient_id),
                verbose=file_info[x,]$verbose,
                executor_id=task_id,
                mode=file_info[x,]$mode,
                time=file_info[x,]$time,
                hold=file_info[[x,"hold"]])
            },mc.cores=ifelse(mode=="local",1,3))

    }else{
        bam_dir_path=system(paste("realpath",bam_dir),intern=TRUE)
        bam_files=system(paste0("find ",bam_dir_path,"| grep bam$"),intern=TRUE)
        t_files=bam_files[!grepl(normal_id,bam_files)]
        normal=bam_files[grepl(normal_id,bam_files)]

        job_report[["steps"]][["clonet"]]=parallel::mclapply(t_files,FUN=function(tumour){
        job_report<-clonet_trento(
                        tumour=tumour,normal=normal,
                        patient_id=patient_id,
                        version=version,
                        threads=threads,
                        tc=tc,
                        ploidy=ploidy,
                        ram=ram,output_dir=paste0(out_file_dir,patient_id),verbose=verbose,
                        executor_id=task_id,mode=mode,time=time,
                        hold=hold
                    )
        },mc.cores=ifelse(mode=="local",1,3))

    }


    if(wait&&mode=="batch"){
        job_validator(job=unlist_level(named_list=job_report[["steps"]][["clonet"]],var="job_id"),
        time=update_time,verbose=verbose,threads=threads)
    }

    return(job_report)
    }





#' Process a pair of tumour-normal samples using CLONET
#'
#' This function takes a pair of tumour and 
#' normal BAMS and applies the CLONET pipeline
#'
#'
#' @param sif_path Path to singularity image file
#' @param version PCF Select panel version to use
#' @param tumour Path to tumour BAM file 
#' @param normal Path to normal BAM file
#' @param patient_id Patient id. 
#' @param tc Pre-computed tumour content. Default NULL.
#' @param ploidy Pre-computed ploidy. Default NULL.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Executor ID. Default "clonet"
#' @param task_name Name of the task. Default "clonet"
#' @param threads Number of CPU cores to use. Default 3.
#' @param ram RAM memory for batched job. Default 4
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param output_dir Path to the output directory.
#' @param tmp_dir Path to temporary file directory.
#' @param verbose Enables progress messages. Default False.
#' @param hold Job to hold on in batched mode.
#' @export


clonet_trento=function(
    sif_clonet=build_default_sif_list()$sif_clonet$V3,
    version="V3",tumour="",normal="",
    patient_id="",tc=NULL,ploidy=NULL,tmp_dir=".",threads=3,
    ram=4,output_dir=".",verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    executor_id=make_unique_id("clonet"),
    task_name="clonet",mode="local",time="48:0:0",
    update_time=60,wait=FALSE,hold=NULL
){

    argg <- as.list(environment())

    if(!is.na(version)){
        sif_clonet=build_default_sif_list()$sif_clonet[version]
        if(is.null(sif_clonet)){
            stop(paste0(ver, " is not a valid PCF Select panel version"))
        }
    }
  
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir,name="clonet_reports")
    out_file_dir_tmp=set_dir(dir=output_dir,name="clonet_tmp")
    out_file_dir=set_dir(dir=out_file_dir,name=get_file_name(tumour))
    out_file_dir_tmp=set_dir(dir=out_file_dir_tmp,name=get_file_name(tumour))
    

    if(!is.null(tc)){
        tc=paste(" -c ",tc)
    }

    if(!is.null(ploidy)){
        ploidy=paste(" -p ",ploidy)
    }
    
    file_info=data.frame(Patient=patient_id,Tumour=normalizePath(tumour),Normal=normalizePath(normal))

    sample_sheet=paste0(out_file_dir_tmp,"/",patient_id,"_",get_file_name(tumour),"_tmp.txt")
    write.table(file_info,file=sample_sheet,quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")


    exec_code=paste(" singularity run -H ",paste0(normalizePath(getwd()),":/home"), " --app pcfs ",
    sif_clonet, " -s ", normalizePath(sample_sheet) ," -o ", out_file_dir, " -t ", 
    out_file_dir_tmp, " -n " , threads, tc, ploidy)
    

 
    out_file_dir2=set_dir(dir=out_file_dir,name="batch")
    job=build_job(executor_id=executor_id,task_id=make_unique_id(task_id))

    if(mode=="batch"){
        batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,threads=threads,
        output_dir=out_file_dir2)
        exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
    }

    if(verbose){
        print_verbose(job=job,arg=argg,exec_code=exec_code)
    }

    job_report=build_job_report(
        job_id=job,
        executor_id=executor_id,
        exec_code=exec_code, 
        task_id=task_id,
        input_args=argg,
        out_file_dir=out_file_dir,
        out_files=list(
            NA
        )
    )

    error=execute_job(exec_code=exec_code)
    if(error!=0){
        stop("Clonet failed to run due to unknown error.
        Check std error for more information.")
    }
    
    if(wait&&mode=="batch"){
        batch_validator(job=job_report$job_id,
        time=update_time,verbose=verbose,threads=threads)
    }

    return(job_report)
}


#' View CLONET output
#'
#' Allows to visualize CLONET output files
#'
#'
#' @param method Select view mode. Options [beta_log2, ai, snps, stats]
#' @param clonet_dir Path to clonet directory of directories.
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Executor ID. Default "preprocess_trento"
#' @param task_name Name of the task. Default "preprocess_trento"
#' @param threads Number of CPU cores to use. Default 3.
#' @param ram RAM memory for batched job. Default 4
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param output_dir Path to the output directory.
#' @param tmp_dir Path to temporary file directory.
#' @param verbose Enables progress messages. Default False.
#' @param hold Job to hold on in batched mode.
#' @export


clonet_view_trento=function(method="log2_beta", clonet_dir="",threads=3,
    ram=4,output_dir=".",verbose=FALSE,sample_labels=NA,
    batch_config=build_default_preprocess_config(),
    executor_id=make_unique_id("clonet_view"),
    task_name="clonet_view",mode="local",time="48:0:0",
    update_time=60,wait=FALSE,hold=NULL,cn_list=build_default_cn_list(),
    clonet_dirs=build_default_clonet_dir_list()
){

    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir)
    validation=unlist(lapply(clonet_dir,FUN=check_clonet_output,clonet_dirs=clonet_dirs))
    clonet_dir=clonet_dir[validation]
    plt_data=list()
    
    ## Reac CN data
    cn_input_files=paste0(clonet_dir,"/",clonet_dirs$cn_snv_calls$CN_SNVs_calls.csv)
    cn_data=lapply(cn_input_files,FUN=function(x){
        tryCatch({
               cn_data=read.table(x,sep=",",header=TRUE)
        },error=function(e){
            return()
            ##warning(paste0("Could not find file ",x))
            
        })
    })


    ### Read TC data
    tc_input_files=paste0(clonet_dir,"/",clonet_dirs$tcEstimation$tc_estimations_CLONETv2.tsv)
    tc_data=lapply(tc_input_files,FUN=function(x){
        tryCatch({
               tc_data=read.table(x,sep=",",header=TRUE)
        },error=function(e){
            return()
            ##warning(paste0("Could not find file ",x))

        })
    })

    
    ### Gene Annotations
    annot_input_file <- system.file("extdata", "gene_annotations_hg19.bed", package="ULPwgs")
    annot_data=read.delim( annot_input_file,header=TRUE)
    panel=annot_data %>% dplyr::group_by(PANEL_VERSION) %>% dplyr::summarise(N=dplyr::n())
    panel$diff=abs(panel$N-length(unique(cn_data$gene)))
    

  
    cn_data=dplyr::bind_rows(cn_data)
    cn_data=dplyr::left_join(cn_data,cn_list,by=c("cn.call.corr"="cn"))
    cn_data=dplyr::left_join(cn_data,annot_data %>% 
    dplyr::filter(PANEL_VERSION==panel[which.min(panel$diff),]$PANEL_VERSION),
    by=c("gene"="hgnc_symbol"))
    if(!is.na(sample_labels)){
        label_annotation=data.frame(s_name=names(sample_labels),id=sample_labels)
        cn_data=dplyr::left_join(cn_data,label_annotation,by=c("sample"="id"))
        cn_data=cn_data %>% dplyr::arrange(s_name,sample)
        cn_data$g_order=as.numeric(as.factor(cn_data$gene))
      
        cn_data$s_order=as.numeric(as.factor(paste0(cn_data$s_name,"_",cn_data$sample)))
    }else{
        cn_data=cn_data %>% dplyr::arrange(sample,gene)
        cn_data$s_name=cn_data$sample
        cn_data$s_order=as.numeric(as.factor(cn_data$sample))
        cn_data$g_order=as.numeric(as.factor(cn_data$gene))
    }

   
    tc_data=dplyr::bind_rows(tc_data)


    plt_data[["cn_data"]]=cn_data
    plt_data[["tc_data"]]=tc_data


    if(method=="log2_beta"){
        clonet_log2_beta(plt_data=plt_data[["cn_data"]])
    } else if(method=="ai"){
        clonet_ai(plt_data=plt_data[["cn_data"]])
    }else if(method=="cn"){
        clonet_cn(plt_data=plt_data[["cn_data"]])
    }else if(method=="snp"){
        clonet_snp(plt_data=plt_data)
    }else if(method=="stats"){
        clonet_stats(plt_data=plt_data[["cn_data"]])
    }else if(method=="pathways"){
        clonet_pathways(plt_data=plt_data[["cn_data"]])
    }else if(method=="tc"){
        clonet_tc(plt_data=plt_data[["cn_data"]])
    }
}




#' @export
clonet_log2_beta=function(plt_data){
    server_log2_beta=function(input,output,session){

        boxes=list()
        lapply(unique(plt_data$s_name),FUN=function(id){
            tmp_plt_data=plt_data %>% dplyr::filter(s_name==id)
            output[[paste0(id,"_plot")]]<- shiny::renderPlot({
                plot_log2_beta(tmp_plt_data,
                gene_tg=any(grepl("Target",input[[paste0(id,"_gene_type")]])),
                gene_ctrl=any(grepl("Control",input[[paste0(id,"_gene_type")]])),
                gene_other=any(grepl("Other",input[[paste0(id,"_gene_type")]])),
                log2_limit=input[[paste0(id,"_log2_limit")]],
                gene_lbl_evi=input[[paste0(id,"_gene_lbl_evi")]],
                gene_lbl_beta_low=input[[paste0(id,"_gene_lbl_beta_low")]],
                gene_lbl_beta_high=input[[paste0(id,"_gene_lbl_beta_high")]],
                gene_lbl=input[[paste0(id,"_gene_lbl")]],
                gene_lbl_size=input[[paste0(id,"_gene_lbl_size")]]
                )
            })
        
            boxes[[id]] <<- shinydashboardPlus::box(
            width = 12, 
            id =paste0(id,"_box"),
            footer=paste0("Tumour Content=",unique(tmp_plt_data$tc),
            "; Ploidy=",unique(tmp_plt_data$ploidy)),
            sidebar =shinydashboardPlus::boxSidebar(
                id=paste0(id,"_sb"),
                width=40,
             
                shiny::radioButtons(paste0(id,"_gene_lbl"), "Labels:", c(
                "Genes" = 1,
                "SNPS" = 2, "Informative SNPS" = 3, "No labels" = 4
                ),
                selected = 4
                ),
                shiny::sliderInput(paste0(id,"_gene_lbl_evi"), "AI Evidence:",
                min = 0, max = 1, value = 0.2, step = 0.1, ticks = FALSE
                ),
                shiny::sliderInput(paste0(id,"_gene_lbl_beta_low"), "Beta Low Labels:",
                min = 0, max = 1, value = 0.1, step = 0.1, ticks = FALSE
                ),
                shiny::sliderInput(paste0(id,"_gene_lbl_beta_high"), "Beta High Labels:",
                min = 0, max = 1, value = 1, step = 0.1, ticks = FALSE
                ),
                shiny::sliderInput(paste0(id,"_log2_limit"), "Log2_Limit:",
                min = 1, max = 10, value = 2, step = 0.1, ticks = FALSE
                ),
                shiny::sliderInput(paste0(id,"_gene_lbl_size"), "Label Size:",
                min = 0, max = 10, value = 2, step = 0.1, ticks = FALSE
                ),
                shinyWidgets::awesomeCheckbox(
                    inputId = paste0(id,"_mdl_mrg"),
                    label = "Show Marginal Plots",
                    value = TRUE
                ),
                shinyWidgets::awesomeCheckboxGroup(
                    inputId = paste0(id,"_gene_type"),
                    label = "Genes:",
                    choices = c("Target", "Control", "Other"),
                    selected = c("Target")
                )
        ),
        shiny::plotOutput(paste0(id,"_plot")),
        title =id, collapsible = TRUE,
        collapsed = FALSE, solidHeader = TRUE
    )
    })
    
    output[["UI"]] <- shiny::renderUI({
        fluidRow(boxes)
    })
    session$onSessionEnded(function() {
        stopApp()
    })
}
    shiny::shinyApp(ui = build_ui, server = server_log2_beta)
}




#' @export
clonet_pathways=function(plt_data){

    server_pathways=function(input,output,session){
        boxes=list()
        summ_plt_data=cn_data %>%dplyr::filter(Class!="") %>% 
        dplyr::group_by(s_name,s_order,Class,cn_class,class_col,tc,ploidy) %>% 
        dplyr::summarise(N=n()) %>% dplyr::group_by(s_name,s_order,Class) %>% 
        dplyr::mutate(Freq=N/sum(N))%>% dplyr::ungroup() %>% 
        tidyr::complete(tidyr::nesting(sample,tc,ploidy),
        Class,tidyr::nesting(cn_class,class_col),fill=list(Freq=0,N=0))

        lapply(unique(plt_data$sample),FUN=function(id){
            tmp_plt_data=summ_plt_data %>% dplyr::filter(sample==id)
            output[[paste0(id,"_plot")]]<- shiny::renderPlot({
                plot_pathways(summ_plt_data)
            })

            boxes[[id]] <<- shinydashboardPlus::box(
            width = 12, 
            id =paste0(id,"_box"),
            footer=paste0("Tumour Content=",unique(tmp_plt_data$tc),
            "; Ploidy=",unique(tmp_plt_data$ploidy)),
            sidebar=shinydashboardPlus::boxSidebar(
                id=paste0(id,"_sb"),
                width=26
        ),
        shiny::plotOutput(paste0(id,"_plot")),
        title =id, collapsible = TRUE,
        collapsed = FALSE, solidHeader = TRUE
    )
    })
    
    output[["UI"]] <- shiny::renderUI({
        fluidRow(boxes)
    })
    session$onSessionEnded(function() {
        stopApp()
    })
}
    shiny::shinyApp(ui = build_ui, server = server_pathways)
}



#' @export
clonet_ai=function(plt_data){

    server_ai=function(input,output,session){
            output[[paste0("ai_plot")]]<- shiny::renderPlot({
        
               plot_ai(
                    plt_data,
                    ai_limit=input[["ai_gene_lbl_evi"]],
                    gene_tg=any(grepl("Target",input[["ai_gene_type"]])),
                    gene_ctrl=any(grepl("Control",input[["ai_gene_type"]])),
                    gene_other=any(grepl("Other",input[["ai_gene_type"]])),
                    min_log2_limit=input[["ai_gene_min_log2"]],
                    max_log2_limit = input[["ai_gene_max_log2"]]
                )
         
            })

            my_box <- shinydashboardPlus::box(
                width = 12,
                id="ai_box",
                footer=paste0("Min TC=",min(plt_data$tc,na.rm=TRUE),"; ",
                "Max TC=",max(plt_data$tc,na.rm=TRUE),"; ",
                "Min Ploidy=",min(plt_data$ploidy,na.rm=TRUE),"; ",
                "Max Ploidy=",max(plt_data$ploidy,na.rm=TRUE)),
                sidebar =shinydashboardPlus::boxSidebar(
                    id="ai_sb",
                    width=26,
                    shiny::sliderInput("ai_gene_lbl_evi", "AI Evidence:",
                    min = 0, max = 1, value = 0.2, step = 0.1, ticks = FALSE
                    ),
                    shiny::sliderInput("ai_gene_min_log2", "Min Log Limit:",
                    min = -1, max = 0, value = -0.3, step = 0.1, ticks = FALSE
                    ),
                    shiny::sliderInput("ai_gene_max_log2", "Max Log Limit:",
                    min = 0, max = 4, value = 0.5, step = 0.1, ticks = FALSE
                    ),
                    shinyWidgets::awesomeCheckboxGroup(
                        inputId = "ai_gene_type",
                        label = "Genes:",
                        choices = c("Target", "Control", "Other"),
                        selected = c("Target")
                    )
            ),
            shiny::plotOutput("ai_plot",height=length(unique(plt_data$gene))*20),
            title ="Allelic Imbalance", collapsible = TRUE,
            collapsed = FALSE, solidHeader = TRUE
     )
        output[["UI"]] <- shiny::renderUI({
            fluidRow(my_box)
        })
    session$onSessionEnded(function() {
        stopApp()
    })
    }

    shiny::shinyApp(ui = build_ui, server = server_ai)
}

#' @export
clonet_cn=function(plt_data){

    server_cn=function(input,output,session){

        boxes=list()
        lapply(unique(plt_data$s_name),FUN=function(id){
            tmp_plt_data=plt_data %>% dplyr::filter(s_name==id)
            output[[paste0(id,"_plot")]]<- shiny::renderPlot({
                
                plot_cn(tmp_plt_data,
                gene_tg=any(grepl("Target",input[[paste0(id,"_gene_type")]])),
                gene_ctrl=any(grepl("Control",input[[paste0(id,"_gene_type")]])),
                gene_other=any(grepl("Other",input[[paste0(id,"_gene_type")]])),
                cn_limit=input[[paste0(id,"_cn_limit")]],
                gene_lbl_evi=input[[paste0(id,"_gene_lbl_evi")]],
                gene_lbl_beta_low=input[[paste0(id,"_gene_lbl_beta_low")]],
                gene_lbl_beta_high=input[[paste0(id,"_gene_lbl_beta_high")]],
                gene_lbl=input[[paste0(id,"_gene_lbl")]],
                gene_lbl_size=input[[paste0(id,"_gene_lbl_size")]]
                )
            })
    
            boxes[[id]] <<- shinydashboardPlus::box(
            width = 12, 
            id =paste0(id,"_box"),
            footer=paste0("Tumour Content=",unique(tmp_plt_data$tc),
            "; Ploidy=",unique(tmp_plt_data$ploidy)),
            sidebar =shinydashboardPlus::boxSidebar(
                
                id=paste0(id,"_sb"),
                width=26,
             
                    shiny::radioButtons(paste0(id,"_gene_lbl"), "Labels:", c(
                "Genes" = 1,
                "SNPS" = 2, "Informative SNPS" = 3, "No labels" = 4
                ),
                selected = 4
                ),
                shiny::sliderInput(paste0(id,"_gene_lbl_evi"), "AI Evidence:",
                min = 0, max = 1, value = 0.2, step = 0.1, ticks = FALSE
                ),
                shiny::sliderInput(paste0(id,"_gene_lbl_beta_low"), "Beta Low Labels:",
                min = 0, max = 1, value = 0.1, step = 0.1, ticks = FALSE
                ),
                shiny::sliderInput(paste0(id,"_gene_lbl_beta_high"), "Beta High Labels:",
                min = 0, max = 1, value = 1, step = 0.1, ticks = FALSE
                ),
                shiny::sliderInput(paste0(id,"_cn_limit"), "CN_Limit:",
                min = 1, max = 10, value = 3, step = 1, ticks = FALSE
                ),
                shiny::sliderInput(paste0(id,"_gene_lbl_size"), "Label Size:",
                min = 0, max = 10, value = 2, step = 0.1, ticks = FALSE
                ),
                shinyWidgets::awesomeCheckbox(
                    inputId = paste0(id,"_mdl_mrg"),
                    label = "Show Marginal Plots",
                    value = TRUE
                ),
                shinyWidgets::awesomeCheckboxGroup(
                    inputId = paste0(id,"_gene_type"),
                    label = "Genes:",
                    choices = c("Target", "Control", "Other"),
                    selected = c("Target")
                )
        ),
        shiny::plotOutput(paste0(id,"_plot"),width="500px"),
        title =id, collapsible = TRUE,
        collapsed = FALSE, solidHeader = TRUE
    )
    })
    
    output[["UI"]] <- shiny::renderUI({
        fluidRow(boxes)
    })
    session$onSessionEnded(function() {
        stopApp()
    })
}
    shiny::shinyApp(ui = build_ui, server = server_cn)
}



#' @export
clonet_tc=function(plt_data){

    server_tc=function(input,output,session){

        boxes=list()
        lapply(unique(plt_data$s_name),FUN=function(id){
            tmp_plt_data=plt_data %>% dplyr::filter(s_name==id)
            output[[paste0(id,"_plot")]]<- shiny::renderPlot({
                plot_tc(tmp_plt_data,
                gene_tg=any(grepl("Target",input[[paste0(id,"_gene_type")]])),
                gene_ctrl=any(grepl("Control",input[[paste0(id,"_gene_type")]])),
                gene_other=any(grepl("Other",input[[paste0(id,"_gene_type")]])),
                gene_ai_limit=input[[paste0(id,"_gene_ai_limit")]],
                gene_log2_loss=input[[paste0(id,"_gene_log2_loss")]]
                )
            })

            output[[paste0(id,"_table")]]<- DT::renderDataTable({
                tc_table(tmp_plt_data,
                gene_tg=any(grepl("Target",input[[paste0(id,"_gene_type")]])),
                gene_ctrl=any(grepl("Control",input[[paste0(id,"_gene_type")]])),
                gene_other=any(grepl("Other",input[[paste0(id,"_gene_type")]])),
                gene_ai_limit=input[[paste0(id,"_gene_ai_limit")]],
                gene_log2_loss=input[[paste0(id,"_gene_log2_loss")]]
                )
            },extensions = c('Buttons','FixedHeader',
            'KeyTable','ColReorder','FixedColumns'),

                options = list(scrollX = TRUE,dom = 'Bfrtip', 
                colReorder = TRUE,fixedHeader = TRUE,
                keys = TRUE,
                dom = 't',
                fixedColumns = list(leftColumns = 3),
                buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))



            boxes[[id]] <<- shinydashboardPlus::box(
            width = 12, 
            id =paste0(id,"_box"),
            footer=paste0("TC=",unique(tmp_plt_data$tc),";TC_Manual=",
            unique(tmp_plt_data$tc_manual),
            "; Ploidy=",unique(tmp_plt_data$ploidy)),
            sidebar =shinydashboardPlus::boxSidebar(
                id=paste0(id,"_sb"),
                width=26,
                shiny::sliderInput(paste0(id,"_gene_ai_limit"), "AI Evidence:",
                min = 0, max = 1, value = 0, step = 0.01, ticks = FALSE
                ),
                shiny::sliderInput(paste0(id,"_gene_log2_loss"), "Log2 Loss Limit:",
                min = -1, max = 0, value = -0.05, step = 0.05, ticks = FALSE
                ),
                shinyWidgets::awesomeCheckboxGroup(
                    inputId = paste0(id,"_gene_type"),
                    label = "Genes:",
                    choices = c("Target", "Control", "Other"),
                    selected = c("Target")
                )
        ),
        shiny::plotOutput(paste0(id,"_plot")),
        DT::dataTableOutput(paste0(id,"_table"))
        ,
        title =id, collapsible = TRUE,
        collapsed = FALSE, solidHeader = TRUE
     )
    })
    
    output[["UI"]] <- shiny::renderUI({
        fluidRow(boxes)
    })
    session$onSessionEnded(function() {
        stopApp()
    })
}
    shiny::shinyApp(ui = build_ui, server = server_tc)
}






#' @export
clonet_stats=function(plt_data){
    server_stats=function(input,output,session){
            output[[paste0("stats_plot")]]<- shiny::renderPlot({
                
               plot_stats(
                    plt_data,
                    gene_tg=any(grepl("Target",input[["stats_gene_type"]])),
                    gene_ctrl=any(grepl("Control",input[["stats_gene_type"]])),
                    gene_other=any(grepl("Other",input[["stats_gene_type"]]))
                )
         
            })

            my_box <- shinydashboardPlus::box(
                width = 12,
                id="stats_box",
                footer=paste0("Min TC=",min(plt_data$tc,na.rm=TRUE),"; ",
                "Max TC=",max(plt_data$tc,na.rm=TRUE),"; ",
                "Min Ploidy=",min(plt_data$ploidy,na.rm=TRUE),"; ",
                "Max Ploidy=",max(plt_data$ploidy,na.rm=TRUE)),
                sidebar =shinydashboardPlus::boxSidebar(
                    id="stats_sb",
                    width=26,
                    shiny::sliderInput("stats_gene_lbl_evi", "AI Evidence:",
                    min = 0, max = 1, value = 0.2, step = 0.1, ticks = FALSE
                    ),
                    shinyWidgets::awesomeCheckboxGroup(
                        inputId = "stats_gene_type",
                        label = "Genes:",
                        choices = c("Target", "Control", "Other"),
                        selected = c("Target")
                    )
            ),
            shiny::plotOutput("stats_plot"),
            title ="Aberration Statistics", collapsible = TRUE,
            collapsed = FALSE, solidHeader = TRUE
     )
        output[["UI"]] <- shiny::renderUI({
            fluidRow(my_box)
        })
    session$onSessionEnded(function() {
        stopApp()
    })
    }

    shiny::shinyApp(ui = build_ui, server = server_stats)
}


#' @export
clonet_snp=function(plt_data){

    server_snp=function(input,output,session){

        boxes=list()
        lapply(unique(plt_data$sample),FUN=function(id){
            tmp_plt_data=plt_data %>% dplyr::filter(sample==id)
            output[[paste0(id,"_plot")]]<- shiny::renderPlot({
                
                plot_cn(tmp_plt_data,
                gene_tg=any(grepl("Target",input[[paste0(id,"_gene_type")]])),
                gene_ctrl=any(grepl("Control",input[[paste0(id,"_gene_type")]])),
                gene_other=any(grepl("Other",input[[paste0(id,"_gene_type")]])),
                cn_limit=input[[paste0(id,"_cn_limit")]],
                gene_lbl_evi=input[[paste0(id,"_gene_lbl_evi")]],
                gene_lbl_beta_low=input[[paste0(id,"_gene_lbl_beta_low")]],
                gene_lbl_beta_high=input[[paste0(id,"_gene_lbl_beta_high")]],
                gene_lbl=input[[paste0(id,"_gene_lbl")]],
                gene_lbl_size=input[[paste0(id,"_gene_lbl_size")]]
                )
            })
    
            boxes[[id]] <<- shinydashboardPlus::box(
            width = 12, 
            id =paste0(id,"_box"),
            footer=paste0("Tumour Content=",unique(tmp_plt_data$tc),
            "; Ploidy=",unique(tmp_plt_data$ploidy)),
            sidebar =shinydashboardPlus::boxSidebar(
                
                id=paste0(id,"_sb"),
                width=26,
             
                    shiny::radioButtons(paste0(id,"_gene_lbl"), "Labels:", c(
                "Genes" = 1,
                "SNPS" = 2, "Informative SNPS" = 3, "No labels" = 4
                ),
                selected = 4
                ),
                shiny::sliderInput(paste0(id,"_gene_lbl_evi"), "AI Evidence:",
                min = 0, max = 1, value = 0.2, step = 0.1, ticks = FALSE
                ),
                shiny::sliderInput(paste0(id,"_gene_lbl_beta_low"), "Beta Low Labels:",
                min = 0, max = 1, value = 0.1, step = 0.1, ticks = FALSE
                ),
                shiny::sliderInput(paste0(id,"_gene_lbl_beta_high"), "Beta High Labels:",
                min = 0, max = 1, value = 1, step = 0.1, ticks = FALSE
                ),
                shiny::sliderInput(paste0(id,"_cn_limit"), "CN_Limit:",
                min = 1, max = 10, value = 3, step = 1, ticks = FALSE
                ),
                shiny::sliderInput(paste0(id,"_gene_lbl_size"), "Label Size:",
                min = 0, max = 10, value = 2, step = 0.1, ticks = FALSE
                ),
                shinyWidgets::awesomeCheckbox(
                    inputId = paste0(id,"_mdl_mrg"),
                    label = "Show Marginal Plots",
                    value = TRUE
                ),
                shinyWidgets::awesomeCheckboxGroup(
                    inputId = paste0(id,"_gene_type"),
                    label = "Genes:",
                    choices = c("Target", "Control", "Other"),
                    selected = c("Target")
                )
        ),
        shiny::plotOutput(paste0(id,"_plot")),
        title =id, collapsible = TRUE,
        collapsed = FALSE, solidHeader = TRUE
    )
    })
    
    output[["UI"]] <- shiny::renderUI({
        fluidRow(boxes)
    })
    session$onSessionEnded(function() {
        stopApp()
    })
}
    shiny::shinyApp(ui = build_ui, server = server_cn)
}




#' @export


check_clonet_output=function(clonet_dir="",
clonet_dirs=build_default_clonet_dir_list()){

    validation=all(unlist(lapply(unlist(clonet_dirs),FUN=function(file){

        search_dir=paste0(clonet_dir,"/",file)

        file.exists(search_dir)
    })))

    return(validation)
    
}



