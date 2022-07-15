#' Build default tool configuration
#' 
#' @param steps Build configuration for steps
#' @param steps_list List with steps
#' 
#' @export


build_default_config=function(steps=names(build_default_steps_list()),
steps_list=build_default_steps_list()){
    config_list=lapply(steps,FUN=build_step,steps_list=steps_list)
    config_list=dplyr::bind_rows(config_list)
    row.names(config_list)=config_list$name
    return(config_list)
}




#' Build default tool configuration
#' 
#'
#' @param step Step to build for configuration
#' @param steps_list Step list
#' @export


build_step=function(step,steps_list=build_default_steps_list()){
    step_info=steps_list[[step]]
    parameter=names(step_info[["args"]])
    if(!is.null(parameter)){
         step_info$args=paste0(paste0(parameter,"=",
         step_info[["args"]]),collapse="|")
    }else{

        step_info$args=""
    }
   
    step_info$name=step
    return(data.frame(step_info,stringsAsFactors = FALSE))
}


#' Build default tool configuration
#' 
#'
#' @param steps List with pipeline step information
#' @export


build_default_steps_list=function(
    steps=list(
        pre_fastqc=list(
            order=1,
            step=TRUE,
            threads=6,
            ram=8,
            mode="local",
            verbose=TRUE,
            time="48:0:0",
            args=list()),
        trimming=list(
            order=2,
            step=TRUE,
            threads=6,
            ram=8,
            mode="local",
            verbose=TRUE,
            time="48:0:0",
            args=list(
                mean_quality=0,
                min_length=18,
                max_length="",
                xadapt="",
                yadapt=""
            )),

        post_fastqc=list(
            order=3,
            step=TRUE,
            threads=6,
            ram=8,
            mode="local",
            verbose=TRUE,
            time="48:0:0",
            args=list()),
        
        alignment=list(
            order=4,
            step=TRUE,
            threads=12,
            ram=16,
            mode="local",
            verbose=TRUE,
            time="48:0:0",
            args=list(
                clean=TRUE,
                stats="all",
                coord_sort=FALSE
            )),
        merge_bam=list(
            order=5,
            step=FALSE,
            threads=6,
            ram=16,
            mode="local",
            verbose=TRUE,
            time="48:0:0",
            args=list(
            )),

        markdups=list(
            order=6,
            step=TRUE,
            threads=12,
            ram=16,
            mode="local",
            verbose=TRUE,
            time="48:0:0",
            args=list(
                remove_duplicates=TRUE)
            ),
         recalibrate=list(
            order=7,
            step=TRUE,
            threads=2,
            ram=8,
            mode="local",
            verbose=TRUE,
            time="48:0:0",
            args=list(
                clean=TRUE
            )
            ),
        alignqc=list(
            order=8,
            step=TRUE,
            threads=2,
            ram=8,
            mode="local",
            verbose=TRUE,
            time="48:0:0",
            args=list()))
        ){
    return(steps)
}





#' Build default tool configuration
#' 
#'
#' @param config List with name,threads,ram,step,mode,verbose,time,args config
#' @export


build_default_parameter_list=function(
    parameters=list(
        parameter=c(
        "name",
        "step",
        "mode",
        "time",
        "threads",
        "ram",
        "args")
        ,
        text=c(
            "Step Name: ",
            "Run: ",
            "Run Mode: ",
            "Run Time: ",
            "Threads: ",
            "Ram: ",
            "Tool Parameters : "
        )
    )
){
    return(parameters)
    
}






#' Build default sample sheet
#' 
#'
#' @param sample_sheet List with project_id,patient_id,sample_id,R1,R2,genome,step,threads,ram,time,mode,verbose configuration
#' @export


build_default_sample_sheet=function(
    sample_sheet=list(
        project_id="TRAILS",
        patient_id=c("TR001","TR001","TR001","TR002","TR002"),
        sample_id=c("ID1","ID1","ID2","ID1","ID2"),
        sequencing_type="BULK",
        method_type="CAPTURE",
        method_version="PCF_V3",
        reference="HG19",
        library_id="LB1",
        R1=c("test/test_data/multi_lane/P1_S1_L1_R1.fq.gz",
        "test/test_data/multi_lane/P1_S1_L2_R1.fq.gz",
        "test/test_data/multi_lane/P1_S2_L1_R1.fq.gz",
        "test/test_data/multi_lane/P2_S1_L1_R1.fq.gz",
        "test/test_data/multi_lane/P2_S2_L1_R1.fq.gz"),
        R2=c("test/test_data/multi_lane/P1_S1_L1_R2.fq.gz",
        "test/test_data/multi_lane/P1_S1_L2_R2.fq.gz",
        "test/test_data/multi_lane/P1_S2_L1_R2.fq.gz",
        "test/test_data/multi_lane/P2_S1_L1_R2.fq.gz",
        "test/test_data/multi_lane/P2_S2_L1_R2.fq.gz"),
        step="pre_fastqc=TRUE;trimming=TRUE;post_fastqc=TRUE;alignment=TRUE;markdups=TRUE;recalibrate=TRUE;alignqc=TRUE",
        threads="pre_fastqc=6;trimming=6;post_fastqc=6;alignment=12;merge_bam=6;markdups=12;recalibrate=2;alignqc=6",
        ram="pre_fastqc=8;trimming=8;post_fastqc=8;alignment=16;merge_bam=16;markdups=16;recalibrate=8;alignqc=8",
        time="pre_fastqc=48:0:0;trimming=48:0:0;post_fastqc=48:0:0;alignment=48:0:0;merge_bam=48:0:0;markdups=48:0:0;recalibrate=48:0:0;alignqc=48:0:0",
        mode="pre_fastqc=local;trimming=local;post_fastqc=local;alignment=local;merge_bam=local;markdups=local;recalibrate=local;alignqc=local",
        verbose="pre_fastqc=TRUE;trimming=TRUE;post_fastqc=TRUE;alignment=TRUE;merge_bam=TRUE;markdups=TRUE;recalibrate=TRUE;alignqc=TRUE",
        args="trimming={mean_quality=0|min_length=18|max_length=|xadapt=|yadapt=};alignment={clean=TRUE|stats=all|coord_sort=FALSE};markdups={remove_duplicates=TRUE};recalibrate={clean=TRUE}"
    )
){

    sample_sheet=data.frame(sample_sheet,stringsAsFactors = FALSE)
    return(sample_sheet)
}

#' Build default binaries config
#' 
#'
#' @param binaries List with tool_names and paths 
#' @export



build_default_binary_list=function(
    binaries=
        list(
                pre_fastqc=list(
                    bin_fastqc="~/Scratch/tools/FastQC/bin/fastqc"
                ),
                trimming=list(
                    bin_skewer="~/Scratch/tools/skewer/skewer"
                ),
                post_fastq=list(
                    bin_fastqc="~/Scratch/tools/FastQC/bin/fastqc"
                ),
                alignment=list(
                    bin_bwa="~/Scratch/tools/bwa/bwa",
                    bin_samtools="~/Scratch/tools/samtools/samtools"
                ),
                merge_bams=list(
                    bin_samtools="~/Scratch/tools/samtools/samtools"
                ),
                markdups=list(
                    bin_gatk="~/Scratch/tools/gatk/gatk"
                ),
                recal_gatk=list(
                    bin_samtools="~/Scratch/tools/samtools/samtools",
                    bin_gatk="~/Scratch/tools/gatk/gatk",
                    bin_picard="~/Scratch/tools/picard/build/libs/picard.jar"
                ),
                alignqc=list(  
                    bin_samtools="~/Scratch/tools/samtools/samtools",
                    bin_picard="~/Scratch/tools/picard/build/libs/picard.jar",
                    bin_bedtools="~/Scratch/tools/bedtools2/bin/bedtools"
                )
            )
        ){
         return(binaries)
}

#' Build default binaries config
#' 
#'
#' @param variables List with variable names
#' @export


build_default_variable_list=function(
    variables=list(
        name=c(
            "project",
            "patient",
            "sample",
            "sequencing",
            "method",
            "version",
            "reference",
            "library",
            "run",
            "flowcell",
            "lane"),
        variable=c(
            "project_id",
            "patient_id",
            "sample_id",
            "sequencing_type",
            "method_type",
            "method_version",
            "reference",
            "library_id",
            "run_id",
            "flowcell_id",
            "lane_id"),

        required=c(
            TRUE,
            TRUE,
            TRUE,
            TRUE,
            TRUE,
            TRUE,
            TRUE,
            TRUE,
            FALSE,
            FALSE,
            FALSE),

        needs_type_validation=c(
            FALSE,
            FALSE,
            FALSE,
            TRUE,
            TRUE,
            TRUE,
            TRUE,
            FALSE,
            FALSE,
            FALSE,
            FALSE
        ),
        text=c(
            "Project ID: ",
            "Patient ID: ",
            "Sample ID: ",
            "Sequencing Type: ",
            "Method Type: ",
            "Method Version: ",
            "Reference Genome: ",
            "Library ID: ",
            "Run ID: ",
            "Flowcell ID: ",
            "Lane ID: "
            )
))(
    return(data.frame(variables,stringsAsFactors = FALSE))
)


#' Build default binaries config
#' 
#'
#' @param options List with options for each variable
#' @export

build_default_option_list=function(
    options=list(
                project_id="",
                patient_id="",
                sample_id="",
                sequencing_type=c("BULK","SINGLE_CELL"),
                method_type=c("WGS","CAPTURE","EXOME","RNASEQ"),
                method_version=list(
                    WGS=c("HC","LP","ULP"),
                    CAPTURE=c("PCF_V1","PCF_V2","PCF_V3"),
                    EXOME=c("ALL_EXOME"),
                    RNASEQ=list("ALL_TRANSCRIPTOME")
                ),
                reference=c("HG19","HG38"),
                library_id="",
                flowcell_id="",
                lane_id="")
){
    return(options)
}



#' Build default references
#' 
#'
#' @param references List with reference files
#' @export

build_default_reference_list=function(
    references=list(
            HG19=list(
                reference=list(
                    genome="references/hg19/reference/hs37d5.fa"
                ),
                database=list(
                    all_snps="references/hg19/database/00-All.vcf.gz",
                    all_common="references/hg19/database/00-common_all.vcf.gz"
                )
            ),
            HG38=list(
                reference=list(
                    genome="references/hg38/reference/ucsc.hg38.fa"
                    
                ),
                database=list(
                    all_snps="references/hg38/database/00-All.vcf.gz",
                    all_common="references/hg38/database/00-common_all.vcf.gz"
                )
            )
        )
    ){
        return(references)
    }


