#' Build default tool configuration
#' 
#'
#' @param config List with name,threads,ram,step,mode,verbose,time,args config
#' @export


build_default_config=function(
    config=list(name=c("pre_fastqc","trimming",
    "post_fastqc","alignment","merge_bam","markdups","recalibrate","alignqc"),
    threads=c(6,6,6,12,6,12,2,6),
    ram=c(8,8,8,16,16,16,8,8),
    step=TRUE,
    mode="local",
    verbose="TRUE",
    time="48:0:0",
    args=c("",
    "mean_quality=0|min_length=18|max_length=''|xadapt=''|yadapt=''",
    "",
    "clean=TRUE|stats='all'|coord_sorted=FALSE",
    "",
    "remove_duplicates=TRUE",
    "",
    ""))){
    config=data.frame(config,stringsAsFactors = FALSE)
    row.names(config)=config$name
    return(config)
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
        method_id="TG",
        library_id=1,
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
        genome="hg19",
        step="pre_fastqc=TRUE;trimming=TRUE;post_fastqc=TRUE;alignment=TRUE;merge_bam=TRUE;markdups=TRUE;recalibrate=TRUE;alignqc=TRUE",
        threads="pre_fastqc=6;trimming=6;post_fastqc=6;alignment=12;merge_bam=6;markdups=12;recalibrate=2;alignqc=6",
        ram="pre_fastqc=8;trimming=8;post_fastqc=8;alignment=16;merge_bam=16;markdups=16;recalibrate=8;alignqc=8",
        time="pre_fastqc=48:0:0;trimming=48:0:0;post_fastqc=48:0:0;alignment=48:0:0;merge_bam=48:0:0;markdups=48:0:0;recalibrate=48:0:0;alignqc=48:0:0",
        mode="pre_fastqc=local;trimming=local;post_fastqc=local;alignment=local;merge_bam=local;markdups=local;recalibrate=local;alignqc=local",
        verbose="pre_fastqc=TRUE;trimming=TRUE;post_fastqc=TRUE;alignment=TRUE;merge_bam=TRUE;markdups=TRUE;recalibrate=TRUE;alignqc=TRUE",
        args="trimming={mean_quality=0|min_length=18|max_length=''|xadapt=''|yadapt=''};alignment={clean=TRUE|stats='all'|coord_sorted=FALSE};markdups={remove_duplicates=TRUE}"
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



build_default_binaries_config=function(binaries=list(tool=c("fastqc","skewer","bwa",
"samtools","gatk","picard","bedtools"),path=c("tools/FastQC/bin/fastqc",
"tools/skewer/skewer","tools/bwa/bwa","tools/samtools/samtools",
"tools/gatk/gatk","tools/picard/build/libs/picard.jar",
"tools/bedtools2/bin/bedtools"))){
    binaries=data.frame(binaries,stringsAsFactors = FALSE)
    row.names(binaries)=binaries$name
    return(binaries)
}