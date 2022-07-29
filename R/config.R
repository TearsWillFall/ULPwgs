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
        R1=c(
        "/lustre/scratch/scratch/PCF/test/test_data/multi_lane/P1_S1_L1_R1.fq.gz",
        "/lustre/scratch/scratch/PCF/test/test_data/multi_lane/P1_S1_L2_R1.fq.gz",
        "/lustre/scratch/scratch/PCF/test/test_data/multi_lane/P1_S2_L1_R1.fq.gz",
        "/lustre/scratch/scratch/PCF/test/test_data/multi_lane/P2_S1_L1_R1.fq.gz",
        "/lustre/scratch/scratch/PCF/test/test_data/multi_lane/P2_S2_L1_R1.fq.gz"),
        R2=c(
        "/lustre/scratch/scratch/PCF/test/test_data/multi_lane/P1_S1_L1_R2.fq.gz",
        "/lustre/scratch/scratch/PCF/test/test_data/multi_lane/P1_S1_L2_R2.fq.gz",
        "/lustre/scratch/scratch/PCF/test/test_data/multi_lane/P1_S2_L1_R2.fq.gz",
        "/lustre/scratch/scratch/PCF/test/test_data/multi_lane/P2_S1_L1_R2.fq.gz",
        "/lustre/scratch/scratch/PCF/test/test_data/multi_lane/P2_S2_L1_R2.fq.gz"),
        step=build_default_steps(),
        threads=build_default_step_threads(),
        ram=build_default_step_ram(),
        batch_config=build_default_step_batch(),
        time=build_default_step_time(),
        mode=build_default_step_mode(),
        verbose=build_default_step_verbose(),
        args="trimming={mean_quality=0|min_length=18|max_length=|xadapt=|yadapt=};alignment={clean=TRUE|stats=all|coord_sort=FALSE};markdups={remove_duplicates=FALSE};recalibrate={clean=TRUE}"
    )
){

    sample_sheet=data.frame(sample_sheet,stringsAsFactors = FALSE)
    return(sample_sheet)
}


#' @export

build_default_step_time=function(
    time=list(
        pre_fastqc="48:0:0",
        trimming="48:0:0",
        post_fastqc="48:0:0",
        alignment="48:0:0",
        merge_bam="48:0:0",
        markdups="48:0:0",
        recalibrate="48:0:0",
        alignqc="48:0:0"
    )
){
    paste0(names(time),"=",time,collapse=";")    
}

#' @export
build_default_step_threads=function(
   threads=list(
        pre_fastqc=6,
        trimming=6,
        post_fastqc=6,
        alignment=12,
        merge_bam=12,
        markdups=12,
        recalibrate=12,
        alignqc=6
    )
){
    paste0(names(threads),"=",threads,collapse=";")    
}

#' @export
build_default_steps=function(
    steps=list(
        pre_fastqc=TRUE,
        trimming=TRUE,
        post_fastqc=TRUE,
        alignment=TRUE,
        merge_bam=TRUE,
        markdups=TRUE,
        recalibrate=TRUE,
        alignqc=TRUE
    )
){
    paste0(names(steps),"=",steps,collapse=";")    
}

#' @export

build_default_step_mode=function(
    mode=list(
        pre_fastqc="local",
        trimming="local",
        post_fastqc="local",
        alignment="local",
        merge_bam="local",
        markdups="local",
        recalibrate="local",
        alignqc="local"
    )
){
    paste0(names(mode),"=",mode,collapse=";")    
}

#' @export

build_default_step_verbose=function(
    verbose=list(
        pre_fastqc=TRUE,
        trimming=TRUE,
        post_fastqc=TRUE,
        alignment=TRUE,
        merge_bam=TRUE,
        markdups=TRUE,
        recalibrate=TRUE,
        alignqc=TRUE
    )
){
    paste0(names(verbose),"=",verbose,collapse=";")    
}

#' @export

build_default_step_batch=function(
     config=list(
        pre_fastqc=build_default_preprocess_config(),
        trimming=build_default_preprocess_config(),
        post_fastqc=build_default_preprocess_config(),
        alignment=build_default_preprocess_config(),
        merge_bam=build_default_preprocess_config(),
        markdups=build_default_preprocess_config(),
        recalibrate=build_default_preprocess_config(),
        alignqc=build_default_preprocess_config()

    )
){
       paste0(names(config),"=",config,collapse=";")    
}

#' @export
build_default_step_ram=function(
    ram=list(
        pre_fastqc=6,
        trimming=6,
        post_fastqc=6,
        alignment=6,
        merge_bam=6,
        markdups=8,
        recalibrate=8,
        alignqc=6
    )
){
    paste0(names(ram),"=",ram,collapse=";")    
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
                    bin_fastqc=build_default_tool_binary_list()$bin_fastqc
                ),
                trimming=list(
                    bin_skewer=build_default_tool_binary_list()$bin_skewer
                ),
                post_fastq=list(
                    bin_fastqc=build_default_tool_binary_list()$bin_fastqc
                ),
                alignment=list(
                    bin_bwa=build_default_tool_binary_list()$bin_bwa,
                    bin_samtools=build_default_tool_binary_list()$bin_samtools
                ),
                merge_bams=list(
                    bin_samtools=build_default_tool_binary_list()$bin_samtools
                ),
                markdups=list(
                    bin_gatk=build_default_tool_binary_list()$bin_gatk
                ),
                recalibrate=list(
                    bin_samtools=build_default_tool_binary_list()$bin_samtools,
                    bin_gatk=build_default_tool_binary_list()$bin_gatk,
                    bin_picard=build_default_tool_binary_list()$bin_picard
                ),
                alignqc=list(  
                    bin_samtools=build_default_tool_binary_list()$bin_samtools,
                    bin_picard=build_default_tool_binary_list()$bin_picard,
                    bin_bedtools=build_default_tool_binary_list()$bin_bedtools
                )
            )
        ){
         return(binaries)
}




#' Build default tools binaries config
#' 
#'
#' @param binaries List with tool_names and paths 
#' @export


build_default_tool_binary_list=function(
    binaries=
        list(
                bin_fastqc="/lustre/scratch/scratch/regmova/tools/FastQC/bin/fastqc",
                bin_skewer="/lustre/scratch/scratch/regmova/tools/skewer/skewer",
                bin_fastqc="/lustre/scratch/scratch/regmova/tools/FastQC/bin/fastqc",
                bin_bwa="/lustre/scratch/scratch/regmova/tools/bwa/bwa",
                bin_ichor_pon="/lustre/scratch/scratch/regmova/tools/ichorCNA/scripts/createPanelOfNormals.R",
                bin_ichor="/lustre/scratch/scratch/regmova/tools/ichorCNA/scripts/runIchorCNA.R",
                bin_samtools="/lustre/scratch/scratch/regmova/tools/samtools/samtools",
                bin_gatk="/lustre/scratch/scratch/regmova/tools/gatk/gatk",
                bin_readcount="/lustre/scratch/scratch/regmova/tools/hmmcopy_utils/bin/readCounter",
                bin_ciri="/lustre/scratch/scratch/regmova/tools/CIRI/CIRI_Full_v2.1.1.jar",
                bin_picard="/lustre/scratch/scratch/regmova/tools/picard/build/libs/picard.jar",
                bin_samtools="/lustre/scratch/scratch/regmova/tools/samtools/samtools",
                bin_picard="/lustre/scratch/scratch/regmova/tools/picard/build/libs/picard.jar",
                bin_bedtools="/lustre/scratch/scratch/regmova/tools/bedtools2/bin/bedtools",
                bin_ciri_quant="/lustre/scratch/scratch/regmova/tools/CIRIquant/bin/CIRIquant"
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
                    genome="/lustre/scratch/scratch/PCF/references/hg19/reference/hs37d5.fa"
                ),
                database=list(
                    all_snps="/lustre/scratch/scratch/PCF/references/hg19/database/00-All.vcf.gz",
                    all_common="/lustre/scratch/scratch/PCF/references/hg19/database/00-common_all.vcf.gz"
                ),
                panel=list(
                    PCF_V2=list(
                        intervals=list(
                            bi="/lustre/scratch/scratch/PCF/references/hg19/panel/v2/cap_tg.interval_list",
                            ti="/lustre/scratch/scratch/PCF/references/hg19/panel/v2/prim_tg.interval_list"
                        ),
                        bed=list(
                            bait="/lustre/scratch/scratch/PCF/references/hg19/panel/v2/cap_tg.bed",
                            target="/lustre/scratch/scratch/PCF/references/hg19/panel/v2/prim_tg.bed",
                            antitarget="/lustre/scratch/scratch/PCF/references/hg19/panel/v2/off_tg.bed"
                        )
                    ),
                    PCF_V3=list(
                        intervals=list(
                            bi="/lustre/scratch/scratch/PCF/references/hg19/panel/v3/cap_tg.interval_list",
                            ti="/lustre/scratch/scratch/PCF/references/hg19/panel/v3/prim_tg.interval_list"
                        ),
                        bed=list(
                            bait="/lustre/scratch/scratch/PCF/references/hg19/panel/v3/cap_tg.bed",
                            target="/lustre/scratch/scratch/PCF/references/hg19/panel/v3/prim_tg.bed",
                            antitarget="/lustre/scratch/scratch/PCF/references/hg19/panel/v3/off_tg.bed"
                        )
                
                    )

                ),
                rnaseq=list(
                    intervals=list(
                        ri="/lustre/scratch/scratch/PCF/references/hg19/rnaseq//intervals/rRNA.interval_list"
                    ),
                    reference=list(
                        ref_flat="/lustre/scratch/scratch/PCF/references/hg19/rnaseq/reference/refFlat.txt"
                    )
                )
            ),
            HG38=list(
                reference=list(
                    genome="/lustre/scratch/scratch/PCF/references/hg38/reference/ucsc.hg38.fa"
                    
                ),
                database=list(
                    all_snps="/lustre/scratch/scratch/PCF/references/hg38/database/00-All.vcf.gz",
                    all_common="/lustre/scratch/scratch/PCF/references/hg38/database/00-common_all.vcf.gz"
                )
            )
        )
    ){
        return(references)
    }



#' Build default references
#' 
#'
#' @param clonet_dirs List with reference files
#' @export

build_default_clonet_dir_list=function(
    clonet_dirs=list(
    beta_computation=list(
        segmentation_amplicons_log2r.tsv="beta_computation/segmentation_amplicons_log2r.tsv",
        segmentation.seg="beta_computation/segmentation.seg"
    ),
    cn_snv_calls=list(
        CN_SNVs_calls.csv="cn_snv_calls/CN_SNVs_calls.csv",
        SNVs_calls.csv="cn_snv_calls/SNVs_calls.csv",
        SNVs_calls_corrected.csv="cn_snv_calls/SNVs_calls_corrected.csv"
    ),
    focalTables=list(
        segmentation_focal.seg="focalTables/segmentation_focal.seg",
        bed_with_rc.Rdata="focalTables/bed_with_rc.Rdata"
    ),
    pacbam=list(

    ),
    peak_correction=list(
        germline_distribution_shifts.Rdata="peak_correction/germline_distribution_shifts.Rdata",
        peak_shifts.tsv="peak_correction/peak_shifts.tsv"
    ),
    tcEstimation=list(
        tc_estimations_CLONETv2.tsv="tcEstimation/tc_estimations_CLONETv2.tsv"
    ),
    abemus=list(
        table_mutations.tsv="abemus/table_mutations.tsv",
        table_mutations_nocommonSNPs.tsv="abemus/table_mutations_nocommonSNPs.tsv",
        tabindex_optimalR.tsv="abemus/tabindex_optimalR.tsv",
        samples_info_file_rpa.tsv="abemus/samples_info_file_rpa.tsv"
    ),
    allelicImbalance=list(
        ai_log2_table.chrX.csv="allelicImbalance/ai_log2_table.chrX.csv",
        ai_log2_table.csv="allelicImbalance/ai_log2_table.csv"
    ),
    annovar=list(
    ))
){
    return(clonet_dirs)
}



#' Build default references
#' 
#'
#' @param cn_list List cn events
#' @export

build_default_cn_list=function(
    cn_list=list(
        cn=c(
            "Unb.Gain",
            "wt",
            "CNNL",
            "NoImb_WTGrayZoneCN",
            "Likely Loss",
            "HemiDel",
            "Imb_WTGrayZoneCN",
            "Bal.Gain",
            "Likely Gain",
            "Uncertain Bal.Gain",
            "HomoDel",
            "Uncertain HomoDel",
            "Deletion on chrX"
        ),
       cn_class=c(
            "Gain",
            "Wt",
            "Wt",
            "Wt",
            "Loss",
            "Loss",
            "Wt",
            "Gain",
            "Gain",
            "Gain",
            "Loss",
            "Loss",
            "Loss"
       ),
       class_col=c(
            "#FF0000",
            "#CDCDCD",
            "#CDCDCD",
            "#CDCDCD",
            "#00B0F0",
            "#00B0F0",
            "#CDCDCD",
            "#FF0000",
            "#FF0000",
            "#FF0000",
            "#00B0F0",
            "#00B0F0",
            "#00B0F0"
       ),
        col=c(
            "#ED7D31",
            "#CDCDCD",
            "#FFE699",
            "#CDCDCD",
            "#00B0F0",
            "#94F6F6",
            "#CDCDCD",
            "#C00000",
            "#FF0000",
            "#C00000",
            "#0070C0",
            "#0070C0",
            "#305496"
        ),
        alpha=c(
            1,
            0.75,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            0.75,
            1,
            0.75,
            0.75
        ),
       order=c(
          3,
          0,
          10,
          1,
          -2,
          -3,
          -1,
          5,
          2,
          4,
          -5,
          -4,
          -10
       )
    )
){
    return(data.frame(cn_list,stringsAsFactors = FALSE))
}




#' Build default myriad module list
#' 
#'
#' @param modules List modules
#' @export

#' @export

build_default_myriad_module_list=function(
    modules=list(
        core=list(
            "rcps-core/1.0.0"
        ),
        compilers=list(
            mpi="compilers mpi",
            intel="compilers/intel/2018/update3",
            openblas=list(
                v0214="openblas/0.2.14/gnu-4.9.2",
                v0370="openblas/0.3.7-serial/gnu-4.9.2"
            )
        ),
        python3=list(
            v3910="python/3.9.10",
            recommended="python3/recommended"
        ),
        r=list(
            recommended="r/recommended"
        )
    )
){
    return(modules)
}



#' Build default myriad module list
#' 
#'
#' @param modules List of modules
#' @export

#' @export

build_default_preprocess_config=function(
   modules=c(
    build_default_myriad_module_list()$r,
    build_default_myriad_module_list()$r)
){
    return(paste0(lapply(modules,FUN=myriad_module,mode="load"),collapse=";"))
}

