#' Build default myriad module list
#' 
#'
#' @param modules List modules
#' @export

build_default_myriad_module_list=function(
    modules=list(
        core=list(
            "rcps-core/1.0.0"
        ),
        compilers=list(
            mpi="compilers mpi",
            gcc="compilers mpi gcc-libs",
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
        ),
        java=list(
            default="java/1.8.0_92",
            v17="java/temurin-17/17.0.2+8"
        )

    )
){
    return(modules)
}


#' Load/Unload modules in myriad
#'
#' This functions creates a constructor for loading/unloading modules in myriad cluster
#'
#'
#' @param mode Path to samtools executable. Default path tools/samtools/samtools.
#' @param module Module to load.
#' @param force Force unload module. Default FALSE
#' @export

myriad_module=function(mode="load",module="",force=FALSE){
    
    tag=""
    if(force){
      tag=" -f "
    }
    if(mode=="load"){
      mdl=paste0("module load ",tag,module)
    }else if (mode=="unload"){
      mdl=paste0("module unload ",tag,module)
    }
    
    return(mdl)
}



#' Build default config for preprocessing step
#' 
#'
#' @param modules List of modules
#' @export

build_default_preprocess_config=function(
   modules=list(
    r=myriad_module(mode="load",module=build_default_myriad_module_list()$r))
){
    return(paste0(modules,collapse=";"))
}





#' Build default tool configuration
#' 
#' @param steps Build configuration for steps
#' @param steps_list List with steps
#' 
#' @export

build_default_config=function(
    steps=names(build_default_steps_list()),
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


build_step=function(
    step,steps_list=build_default_steps_list(),sep="||"){
    step_info=steps_list[[step]]
    parameter=names(step_info[["args"]])
    if(!is.null(parameter)){
         step_info$args=paste0(paste0(parameter,"=",
         step_info[["args"]]),collapse=sep)
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
            step=build_default_step_list()$pre_fastqc,
            threads=build_default_step_threads_list()$pre_fastqc,
            ram=build_default_step_ram_list()$pre_fastqc,
            batch_config=build_default_preprocess_config(),
            mode=build_default_step_mode_list()$pre_fastqc,
            verbose=build_default_step_verbose_list()$pre_fastqc,
            time=build_default_step_time_list()$pre_fastqc,
            args=list()),
        trimming=list(
            order=2,
            step=build_default_step_list()$trimming,
            threads=build_default_step_threads_list()$trimming,
            ram=build_default_step_ram_list()$trimming,
            batch_config=build_default_preprocess_config(),
            mode=build_default_step_mode_list()$trimming,
            verbose=build_default_step_verbose_list()$trimming,
            time=build_default_step_time_list()$trimming,
            args=list(
                mean_quality=0,
                min_length=18,
                max_length="",
                xadapt="",
                yadapt=""
            )),

        post_fastqc=list(
            order=3,
           step=build_default_step_list()$post_fastqc,
            threads=build_default_step_threads_list()$post_fastqc,
            ram=build_default_step_ram_list()$post_fastqc,
            batch_config=build_default_preprocess_config(),
            mode=build_default_step_mode_list()$post_fastqc,
            verbose=build_default_step_verbose_list()$post_fastqc,
            time=build_default_step_time_list()$post_fastqc,
            args=list()),
        
        alignment=list(
            order=4,
            step=build_default_step_list()$alignment,
            threads=build_default_step_threads_list()$alignment,
            ram=build_default_step_ram_list()$alignment,
            batch_config=build_default_preprocess_config(),
            mode=build_default_step_mode_list()$alignment,
            verbose=build_default_step_verbose_list()$alignment,
            time=build_default_step_time_list()$alignment,
            args=list(
                clean=TRUE,
                stats="all",
                coord_sort=FALSE
            )),

        pre_alignqc=list(
            order=5,
            step=build_default_step_list()$pre_alignqc,
            threads=build_default_step_threads_list()$pre_alignqc,
            ram=build_default_step_ram_list()$pre_alignqc,
            batch_config=build_default_preprocess_config(),
            mode=build_default_step_mode_list()$pre_alignqc,
            verbose=build_default_step_verbose_list()$pre_alignqc,
            time=build_default_step_time_list()$pre_alignqc,
            args=list()),
        merge_bam=list(
            order=6,
            step=build_default_step_list()$merge_bam,
            threads=build_default_step_threads_list()$merge_bam,
            ram=build_default_step_ram_list()$merge_bam,
            batch_config=build_default_preprocess_config(),
            mode=build_default_step_mode_list()$merge_bam,
            verbose=build_default_step_verbose_list()$merge_bam,
            time=build_default_step_time_list()$merge_bam,
            args=list(
            )),

        markdups=list(
            order=7,
            step=build_default_step_list()$markdups,
            threads=build_default_step_threads_list()$markdups,
            ram=build_default_step_ram_list()$markdups,
            batch_config=build_default_preprocess_config(),
            mode=build_default_step_mode_list()$markdups,
            verbose=build_default_step_verbose_list()$markdups,
            time=build_default_step_time_list()$markdups,
            args=list(
                remove_duplicates=FALSE)
            ),
         recalibrate=list(
            order=8,
            step=build_default_step_list()$recalibrate,
            threads=build_default_step_threads_list()$recalibrate,
            ram=build_default_step_ram_list()$recalibrate,
            batch_config=build_default_preprocess_config(),
            mode=build_default_step_mode_list()$recalibrate,
            verbose=build_default_step_verbose_list()$recalibrate,
            time=build_default_step_time_list()$recalibrate,
            args=list(
                clean=TRUE
            )
            ),
        post_alignqc=list(
            order=9,
            step=build_default_step_list()$post_alignqc,
            threads=build_default_step_threads_list()$post_alignqc,
            ram=build_default_step_ram_list()$post_alignqc,
            batch_config=build_default_preprocess_config(),
            mode=build_default_step_mode_list()$post_alignqc,
            verbose=build_default_step_verbose_list()$post_alignqc,
            time=build_default_step_time_list()$post_alignqc,
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
        "batch_config",
        "args"),
        text=c(
            "Step Name: ",
            "Run: ",
            "Run Mode: ",
            "Run Time: ",
            "Threads: ",
            "Ram: ",
            "Batch Config: ",
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
        "/myriadfs/home/regmova/Scratch/PCF/test/test_data/multi_lane/P1_S1_L1_R1.fq.gz",
        "/myriadfs/home/regmova/Scratch/PCF/test/test_data/multi_lane/P1_S1_L2_R1.fq.gz",
        "/myriadfs/home/regmova/Scratch/PCF/test/test_data/multi_lane/P1_S2_L1_R1.fq.gz",
        "/myriadfs/home/regmova/Scratch/PCF/test/test_data/multi_lane/P2_S1_L1_R1.fq.gz",
        "/myriadfs/home/regmova/Scratch/PCF/test/test_data/multi_lane/P2_S2_L1_R1.fq.gz"),
        R2=c(
        "/myriadfs/home/regmova/Scratch/PCF/test/test_data/multi_lane/P1_S1_L1_R2.fq.gz",
        "/myriadfs/home/regmova/Scratch/PCF/test/test_data/multi_lane/P1_S1_L2_R2.fq.gz",
        "/myriadfs/home/regmova/Scratch/PCF/test/test_data/multi_lane/P1_S2_L1_R2.fq.gz",
        "/myriadfs/home/regmova/Scratch/PCF/test/test_data/multi_lane/P2_S1_L1_R2.fq.gz",
        "/myriadfs/home/regmova/Scratch/PCF/test/test_data/multi_lane/P2_S2_L1_R2.fq.gz"),
        step=collapse_step_list(build_default_step_list()),
        threads=collapse_step_list(build_default_step_threads_list()),
        ram=collapse_step_list(build_default_step_ram_list()),
        batch_config=collapse_step_list(build_default_step_batch_list()),
        time=collapse_step_list(build_default_step_time_list()),
        mode=collapse_step_list(build_default_step_mode_list()),
        verbose=collapse_step_list(build_default_step_verbose_list()),
        args=build_default_args()
    )
){

    sample_sheet=data.frame(sample_sheet,stringsAsFactors = FALSE)
    return(sample_sheet)
}




#' @export

build_default_args=function(
    steps_list=build_default_steps_list(),int_sep=";",
    ext_sep="||"){
    lapply(names(steps_list),FUN=function(x){
        if(length(steps_list[[x]]$args)>0){
            args=collapse_step_list(steps_list[[x]]$args,sep=int_sep)
        }else{
            return()
        }
        paste0(x,"={",args,"}")
    })%>% purrr::discard(is.null) %>% paste(collapse=ext_sep)

}

#' @export

build_default_step_time_list=function(
    time=list(
        pre_fastqc="48:0:0",
        trimming="48:0:0",
        post_fastqc="48:0:0",
        alignment="48:0:0",
        merge_bam="48:0:0",
        pre_alignqc="48:0:0",
        markdups="48:0:0",
        recalibrate="48:0:0",
        post_alignqc="48:0:0"
    )
){
    return(time)  
}

#' @export
build_default_step_threads_list=function(
   threads=list(
        pre_fastqc=6,
        trimming=6,
        post_fastqc=6,
        alignment=8,
        merge_bam=4,
        pre_alignqc=6,
        markdups=8,
        recalibrate=6,
        post_alignqc=6
    )
){
    return(threads)   
}

#' @export
build_default_step_list=function(
    steps=list(
        pre_fastqc=TRUE,
        trimming=TRUE,
        post_fastqc=TRUE,
        alignment=TRUE,
        merge_bam=TRUE,
        pre_alignqc=TRUE,
        markdups=TRUE,
        recalibrate=TRUE,
        post_alignqc=TRUE
    )
){
   return(steps)    
}

#' @export

build_default_step_mode_list=function(
    mode=list(
        pre_fastqc="local",
        trimming="local",
        post_fastqc="local",
        alignment="local",
        merge_bam="local",
        pre_alignqc="local",
        markdups="local",
        recalibrate="local",
        post_alignqc="local"
    )
){
    return(mode)
}

#' @export

build_default_step_verbose_list=function(
    verbose=list(
        pre_fastqc=TRUE,
        trimming=TRUE,
        post_fastqc=TRUE,
        alignment=TRUE,
        merge_bam=TRUE,
        pre_alignqc=TRUE,
        markdups=TRUE,
        recalibrate=TRUE,
        post_alignqc=TRUE
    )
){
   return(verbose)   
}

#' @export

build_default_step_batch_list=function(
     config=list(
        pre_fastqc=build_default_preprocess_config(),
        trimming=build_default_preprocess_config(),
        post_fastqc=build_default_preprocess_config(),
        alignment=build_default_preprocess_config(),
        merge_bam=build_default_preprocess_config(),
        pre_alignqc=build_default_preprocess_config(),
        markdups=build_default_preprocess_config(),
        recalibrate=build_default_preprocess_config(),
        post_alignqc=build_default_preprocess_config()

    )
){
       return(config)
}

#' @export
build_default_step_ram_list=function(
    ram=list(
        pre_fastqc=6,
        trimming=6,
        post_fastqc=6,
        alignment=6,
        merge_bam=6,
        pre_alignqc=6,
        markdups=8,
        recalibrate=8,
        post_alignqc=6
    )
){
    return(ram)  
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
                pre_alignqc=list(  
                    bin_samtools=build_default_tool_binary_list()$bin_samtools,
                    bin_picard=build_default_tool_binary_list()$bin_picard,
                    bin_bedtools=build_default_tool_binary_list()$bin_bedtools
                ),
                markdups=list(
                    bin_gatk=build_default_sif_list()$sif_gatk
                ),
                recalibrate=list(
                    bin_samtools=build_default_tool_binary_list()$bin_samtools,
                    bin_gatk=build_default_sif_list()$sif_gatk,
                    bin_picard=build_default_tool_binary_list()$bin_picard
                ),
                post_alignqc=list(  
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
                bin_fastqc="/myriadfs/home/regmova/Scratch/tools/FastQC/bin/fastqc",
                bin_skewer="/myriadfs/home/regmova/Scratch/tools/skewer/skewer",
                bin_fastqc="/myriadfs/home/regmova/Scratch/tools/FastQC/bin/fastqc",
                bin_bcftools="/myriadfs/home/regmova/Scratch/tools/bcftools/bcftools",
                bin_bgzip="/myriadfs/home/regmova/Scratch/tools/htslib/bin/bgzip",
                bin_tabix="/myriadfs/home/regmova/Scratch/tools/htslib/bin/tabix",
                bin_htsfile="/myriadfs/home/regmova/Scratch/tools/htslib/bin/htsfile",
                bin_nextflow="/myriadfs/home/regmova/bin/nextflow",
                bin_aasuite="/myriadfs/home/regmova/Scratch/tools/AmpliconSuite-pipeline/singularity/run_paa_singularity.py",
                bin_bwa="/myriadfs/home/regmova/Scratch/tools/bwa/bwa",
                bin_shapeit="/myriadfs/home/regmova/Scratch/tools/shapeit4/bin/shapeit",
                bin_allele_counter="/myriadfs/home/regmova/Scratch/tools/alleleCount/bin/alleleCounter",
                bin_strelka=list(
                    somatic="/myriadfs/home/regmova/Scratch/tools/strelka-2.9.10/bin/configureStrelkaSomaticWorkflow.py",
                    germline="/myriadfs/home/regmova/Scratch/tools/strelka-2.9.10/bin/configureStrelkaGermlineWorkflow.py"
                ),
                bin_ucsc=list(
                    bigBedToBed="/myriadfs/home/regmova/Scratch/tools/OtherTools/bigBedToBed",
                    bedToBigBed="/myriadfs/home/regmova/Scratch/tools/OtherTools/bedToBigBed",
                    liftOver="/myriadfs/home/regmova/Scratch/tools/OtherTools/liftOver"
                ),
                bin_sshpass="/myriadfs/home/regmova/Scratch/tools/sshpass/sshpass",
                bin_vep="/myriadfs/home/regmova/Scratch/tools/ensembl-vep/vep",    
                bin_manta="/myriadfs/home/regmova/Scratch/tools/manta-1.6.0/bin/configManta.py",
                bin_ichor_pon="/myriadfs/home/regmova/Scratch/tools/ichorCNA/scripts/createPanelOfNormals.R",
                bin_ichor="/myriadfs/home/regmova/Scratch/tools/ichorCNA/scripts/runIchorCNA.R",
                bin_samtools="/myriadfs/home/regmova/Scratch/tools/samtools/samtools",
                bin_gatk="/myriadfs/home/regmova/Scratch/tools/gatk/gatk",
                bin_readcount="/myriadfs/home/regmova/Scratch/tools/hmmcopy_utils/bin/readCounter",
                bin_ciri="/myriadfs/home/regmova/Scratch/tools/CIRI/CIRI_Full_v2.1.1.jar",
                bin_picard="/myriadfs/home/regmova/Scratch/tools/picard/build/libs/picard.jar",
                bin_samtools="/myriadfs/home/regmova/Scratch/tools/samtools/samtools",
                bin_picard="/myriadfs/home/regmova/Scratch/tools/picard/build/libs/picard.jar",
                bin_bedtools="/myriadfs/home/regmova/Scratch/tools/bedtools2/bin/bedtools",
                bin_ciri_quant="/myriadfs/home/regmova/Scratch/tools/CIRIquant/bin/CIRIquant",
                bin_msing="/myriadfs/home/regmova/Scratch/tools/msings/scripts/run_msings.sh",
                bin_gangstr="/myriadfs/home/regmova/Scratch/tools/GangSTR/bin/GangSTR"
        )
    ){
         return(binaries)
}






#' Build default nf-core pipeline list
#' 
#'
#' @param pipelines List with tool_names and paths 
#' @export


build_default_nf_list=function(
    pipelines=
        list(
            nf_circdna=list(
                name="nf-core/circdna",
                version="1.0.3"
            ),
            nf_sarek="nf-core/sarek",
            nf_rnaseq="nf-core/rnaseq"
        )
    ){
         return(pipelines)
}


#' Build default license list
#' 
#'
#' @param licenses List with tool_names and paths 
#' @export


build_default_license_list=function(
    licenses=list(
        dir="/myriadfs/home/regmova/lic",
        licenses=
        list(
            mosek="/myriadfs/home/regmova/lic/mosek.lic",
            gurobi="/myriadfs/home/regmova/lic/gurobi.lic"
        )
    )
    ){
    return(licenses)
}







#' Build default python enviroment list
#' 
#'
#' @param enviroments List with tool_names and paths 
#' @export


build_default_python_enviroment_list=function(
    enviroments=
        list(
                env_circlemap="/myriadfs/home/regmova/miniconda3/envs/circle-map",
                env_hatchet="/myriadfs/home/regmova/miniconda3/envs/hatchet",
                env_cov="/myriadfs/home/regmova/miniconda3/envs/cov",
                env_medicc2="/myriadfs/home/regmova/miniconda3/envs/medicc_env",
                env_msing="/myriadfs/home/regmova/Scratch/tools/msings/msings-env/bin/activate"
        )
    ){
         return(enviroments)
}




#' Build chromosome list
#' 
#'
#' @param chromosomes List with tool_names and paths 
#' @export


build_default_chr_list=function(
    chromosomes=list(
        canonical=c(seq(1,22),"X","Y")
    )
){
 return(chromosomes)
}
build_default_chr_list()






#' Build default sif list
#' 
#'
#' @param sifs List of sif files
#' @export


build_default_sif_list=function(
    sifs=
        list(
            sif_clonet=list(
                V2="/myriadfs/home/regmova/Scratch/Singularity_Images/pcf_select_v2_14_11_2022.sif",
                V3="/myriadfs/home/regmova/Scratch/Singularity_Images/pcf_select_v3_14_11_2022.sif"
            ),
            sif_preprocess="/myriadfs/home/regmova/Scratch/Singularity_Images/preProcess_latest.sif",
            sif_aa="/myriadfs/home/regmova/Scratch/Singularity_Images/ampliconsuite-pipeline.sif",
            sif_gatk="/myriadfs/home/regmova/Scratch/Singularity_Images/gatk_latest.sif",
            sif_cnvkit="/myriadfs/home/regmova/Scratch/Singularity_Images/cnvkit_latest.sif"
        )
    ){
         return(sifs)
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


#' Build default cache list
#' 
#'
#' @param cache List with options for each variable
#' @export

build_default_cache_list=function(
    cache=list(
        cache_vep="/myriadfs/home/regmova/Scratch/tools/ensembl-vep/.vep"
    )
){
    return(cache)
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
                    chrx_tr="/myriadfs/home/regmova/Scratch/PCF/references/hg19/reference/chrX_TR_database.bed",
                    genome="/myriadfs/home/regmova/Scratch/PCF/references/hg19/reference/hs37d5.fa",
                    access_5k="/myriadfs/home/regmova/Scratch/PCF/references/hg19/reference/hs37d5.access-5k-mappable.bed",
                    access_10k="/myriadfs/home/regmova/Scratch/PCF/references/hg19/reference/hs37d5.access-10k-mappable.bed"
                ),
                phasing=list(
                    G1000=list(
                        old=list(
                            haplotype=paste0("/myriadfs/home/regmova/Scratch/PCF/references/hg19/phasing/G1000/old/1000GP_Phase3/1000GP_Phase3_chr",c(1:22,"X"),".hap.gz"),
                            legend=paste0("/myriadfs/home/regmova/Scratch/PCF/references/hg19/phasing/G1000/old/1000GP_Phase3/1000GP_Phase3_chr",c(1:22,"X"),".legend.gz"),
                            gmap=paste0("/myriadfs/home/regmova/Scratch/PCF/references/hg19/phasing/G1000/old/1000GP_Phase3/genetic_map_chr",c(1:22,"X"),"_combined_b37.txt"),
                            sample="/myriadfs/home/regmova/Scratch/PCF/references/hg19/phasing/G1000/old/1000GP_Phase3/1000GP_Phase3.sample"
                        ),
                        new=list(
                            vcf=paste0("/myriadfs/home/regmova/Scratch/PCF/references/hg19/phasing/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr",c(1:22,"X"),".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"),
                            gmap=paste0("/myriadfs/home/regmova/Scratch/PCF/references/hg19/phasing/G1000/new/chr",c(1:22,"X"),".b37.gmap.gz")
                        )
                    )
                ),
                aa=list(
                    dir="/myriadfs/home/regmova/aa"
                ),
                database=list(
                    all_snps="/myriadfs/home/regmova/Scratch/PCF/references/hg19/database/00-All.vcf.gz",
                    all_common="/myriadfs/home/regmova/Scratch/PCF/references/hg19/database/00-common_all.vcf.gz"
                ),
                annotation=list(
                    reflat="/myriadfs/home/regmova/Scratch/PCF/references/hg19/reference/refFlat.txt",
                    genes="/myriadfs/home/regmova/Scratch/PCF/references/hg19/reference/genes.txt"
                ),
                wgs=list(
                    binned=list(
                            target_1k="/myriadfs/home/regmova/Scratch/PCF/references/hg19/wgs/reference_1k/hs37d5.access-5k-mappable.target.bed",
                            target_5k="/myriadfs/home/regmova/Scratch/PCF/references/hg19/wgs/reference_5k/hs37d5.access-5k-mappable.target.bed",
                            target_50k="/myriadfs/home/regmova/Scratch/PCF/references/hg19/wgs/reference_50k/hs37d5.access-5k-mappable.target.bed"
                    ),
                    variant=list(
                        pon_cn_1k="/myriadfs/home/regmova/Scratch/PCF/references/hg19/wgs/reference_1k/pon_trails_wgs.target_1000.pon.cnn",
                        pon_cn_5k="/myriadfs/home/regmova/Scratch/PCF/references/hg19/wgs/reference_5k/pon_trails_wgs.target_5000.pon.cnn",
                        pon_cn_50k="/myriadfs/home/regmova/Scratch/PCF/references/hg19/wgs/reference_50k/pon_trails_wgs.target_50000.pon.cnn"
                    ),
                    ichorcna=list(
                        pon_500k="/myriadfs/home/regmova/Scratch/tools/ichorCNA/inst/extdata/HD_ULP_PoN_500kb_median_normAutosome_mapScoreFiltered_median.rds",
                        gc_500k="/myriadfs/home/regmova/Scratch/tools/ichorCNA/inst/extdata/gc_hg19_500kb.wig",
                        map_500k="/myriadfs/home/regmova/Scratch/tools/ichorCNA/inst/extdata/map_hg19_500kb.wig",
                        centromere="/myriadfs/home/regmova/Scratch/tools/ichorCNA/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt"
                    )
                ),
                panel=list(
                    PCF_V2=list(
                        intervals=list(
                            bi="/myriadfs/home/regmova/Scratch/PCF/references/hg19/panel/v2/cap_tg.interval_list",
                            ti="/myriadfs/home/regmova/Scratch/PCF/references/hg19/panel/v2/prim_tg.interval_list"
                        ),
                        bed=list(
                            bait="/myriadfs/home/regmova/Scratch/PCF/references/hg19/panel/v2/cap_tg.bed",
                            target="/myriadfs/home/regmova/Scratch/PCF/references/hg19/panel/v2/prim_tg.bed"
                        ),
                        variant=list(
                            pon_cn="/myriadfs/home/regmova/Scratch/PCF/references/hg19/panel/v2/pon_seg/pon_seg.cnn",
                            pon_cn_male="/myriadfs/home/regmova/Scratch/PCF/references/hg19/panel/v2/pon_seg/pon_seg_male.cnn",
                            pon_muts="/myriadfs/home/regmova/Scratch/PCF/references/hg19/panel/v2/pon_var/pon_gatkPoN.vcf.gz"
                        ),
                        binned=list(
                            target="/myriadfs/home/regmova/Scratch/PCF/references/hg19/panel/v2/cap_tg.binned.targets.bed",
                            antitarget="/myriadfs/home/regmova/Scratch/PCF/references/hg19/panel/v2/cap_tg.binned.antitargets.bed"
                        ),
                         ascat=list(
                            G1000=list(
                                allele="/myriadfs/home/regmova/Scratch/PCF/references/hg19/panel/v2/ascat/G1000/alleleData_chr",
                                loci="/myriadfs/home/regmova/Scratch/PCF/references/hg19/panel/v2/ascat/G1000/loci_chr",
                                gc="/myriadfs/home/regmova/Scratch/PCF/references/hg19/panel/v2/ascat/G1000/gc.txt",
                                rt="/myriadfs/home/regmova/Scratch/PCF/references/hg19/panel/v2/ascat/G1000/rt.txt"
                            ),
                            battenberg=list(
                                allele="/myriadfs/home/regmova/Scratch/PCF/references/hg19/panel/v2/ascat/battenberg/alleleData_chr",
                                loci="/myriadfs/home/regmova/Scratch/PCF/references/hg19/panel/v2/ascat/battenberg/loci_chr",
                                gc="/myriadfs/home/regmova/Scratch/PCF/references/hg19/panel/v2/ascat/battenberg/gc.txt",
                                rt="/myriadfs/home/regmova/Scratch/PCF/references/hg19/panel/v2/ascat/battenberg/rt.txt"
                            )
                        ),
                          annotation=list(
                            genes="/myriadfs/home/regmova/Scratch/PCF/references/hg19/panel/v2/genes.bed")
                    ),
                    PCF_V3=list(
                        intervals=list(
                            bi="/myriadfs/home/regmova/Scratch/PCF/references/hg19/panel/v3/cap_tg.interval_list",
                            ti="/myriadfs/home/regmova/Scratch/PCF/references/hg19/panel/v3/prim_tg.interval_list"
                        ),
                        bed=list(
                            bait="/myriadfs/home/regmova/Scratch/PCF/references/hg19/panel/v3/cap_tg.bed",
                            target="/myriadfs/home/regmova/Scratch/PCF/references/hg19/panel/v3/prim_tg.bed"
                        ),
                        variant=list(
                            pon_cn="/myriadfs/home/regmova/Scratch/PCF/references/hg19/panel/v3/pon_seg/pon_seg.cnn",
                            pon_cn_male="/myriadfs/home/regmova/Scratch/PCF/references/hg19/panel/v3/pon_seg/pon_seg_male.cnn",
                            pon_muts="/myriadfs/home/regmova/Scratch/PCF/references/hg19/panel/v3/pon_var/pon_gatkPoN.vcf.gz"
                        ),
                        binned=list(
                            target="/myriadfs/home/regmova/Scratch/PCF/references/hg19/panel/v3/cap_tg.binned.targets.bed",
                            antitarget="/myriadfs/home/regmova/Scratch/PCF/references/hg19/panel/v3/cap_tg.binned.antitargets.bed"
                        ),
                        ascat=list(
                            G1000=list(
                                allele="/myriadfs/home/regmova/Scratch/PCF/references/hg19/panel/v3/ascat/G1000/alleleData_chr",
                                loci="/myriadfs/home/regmova/Scratch/PCF/references/hg19/panel/v3/ascat/G1000/loci_chr",
                                gc="/myriadfs/home/regmova/Scratch/PCF/references/hg19/panel/v3/ascat/G1000/gc.txt",
                                rt="/myriadfs/home/regmova/Scratch/PCF/references/hg19/panel/v3/ascat/G1000/rt.txt"
                            ),
                            battenberg=list(
                                allele="/myriadfs/home/regmova/Scratch/PCF/references/hg19/panel/v3/ascat/battenberg/alleleData_chr",
                                loci="/myriadfs/home/regmova/Scratch/PCF/references/hg19/panel/v3/ascat/battenberg/loci_chr",
                                gc="/myriadfs/home/regmova/Scratch/PCF/references/hg19/panel/v3/ascat/battenberg/gc.txt",
                                rt="/myriadfs/home/regmova/Scratch/PCF/references/hg19/panel/v3/ascat/battenberg/rt.txt"
                            )
                        ),
                        annotation=list(
                            genes="/myriadfs/home/regmova/Scratch/PCF/references/hg19/panel/v3/genes.bed"),
                        msi=list(
                            bed="/myriadfs/home/regmova/Scratch/PCF/references/hg19/panel/v3/msi/regions.bed",
                            baseline="/myriadfs/home/regmova/Scratch/PCF/references/hg19/panel/v3/msi/baseline.txt"
                        )
                        )

                ),
                rnaseq=list(
                    intervals=list(
                        ri="/myriadfs/home/regmova/Scratch/PCF/references/hg19/rnaseq//intervals/rRNA.interval_list"
                    ),
                    reference=list(
                        ref_flat="/myriadfs/home/regmova/Scratch/PCF/references/hg19/rnaseq/reference/refFlat.txt"
                    )
                ),
                variant=list(
                    germ_reference="/myriadfs/home/regmova/Scratch/PCF/references/hg19/variant/af-only-gnomad.raw.sites.vcf",
                    biallelic_reference="/myriadfs/home/regmova/Scratch/PCF/references/hg19/variant/small_exac_common_3.vcf",
                    hapmap_reference="/myriadfs/home/regmova/Scratch/PCF/references/hg19/variant/hapmap_3.3.b37.vcf",
                    mills_reference="/myriadfs/home/regmova/Scratch/PCF/references/hg19/variant/Mills_and_1000G_gold_standard.indels.b37.vcf"
                )
            ),
            HG19_ALT=list(
                reference=list(
                    genome="/myriadfs/home/regmova/Scratch/PCF/references/hg19_alt/reference/human_g1k_v37.fasta"
                )


            ),
            HG38=list(
                reference=list(
                    genome="/myriadfs/home/regmova/Scratch/PCF/references/hg38/reference/ucsc.hg38.fa"
                    
                ),
                panel=list(
                    PCF_V3=list(
                        intervals=list(
                            bi="/myriadfs/home/regmova/Scratch/PCF/references/hg19_alt/panel/v3/cap_tg.interval_list",
                            ti="/myriadfs/home/regmova/Scratch/PCF/references/hg19_alt/panel/v3/prim_tg.interval_list"
                        )
                    )
                ),

                database=list(
                    all_snps="/myriadfs/home/regmova/Scratch/PCF/references/hg38/database/00-All.vcf.gz",
                    all_common="/myriadfs/home/regmova/Scratch/PCF/references/hg38/database/00-common_all.vcf.gz"
                )
            )
        )
    ){
        return(references)
}



#' Build default references
#' 
#'
#' @param configs List with reference files
#' @export


build_default_hatchet_config=function(
    configs=list(
        download_panel=list(
            run=FALSE,
            config=list(
                ref_panel = "1000GP_Phase3",
                refpaneldir = "/myriadfs/home/regmova/Scratch/TRIALS/HATCHET/reference/panel"
            )
        ),
        genotype_snps=list(
            run=TRUE,
            config=list(
                mincov = 8,
                maxcov = 300,
                reference_version = "hg19",
                chr_notation = "False"
            )
        ),
        count_alleles=list(
            run=TRUE,
            config=list(
                mincov = 8,
                maxcov = 300
            )
        ),
        count_reads=list(
            run=TRUE,
            config=list(
            )
        ),
        phase_snps=list(
            run=TRUE,
            config=list(
                refpaneldir = "/myriadfs/home/regmova/Scratch/TRIALS/HATCHET/reference/panel"
            )
        ),
        combine_counts=list(
            run=TRUE,
            config=list(
                msr = 3000,
                mtr = 5000
            )
        ),
        cluster_bins=list(
            run=TRUE,
            config=list(
                diploidbaf = 0.10,
                minK = 2,
                maxK = 30,
                tau=0.000001,
                transmat = "\"diag\"",
                decoding = "\"map\"",
                covar = "\"diag\"",
                restarts=10
            )
        ),
        loc_clust=list(
            run=TRUE,
            config=list(
            )
        ),
        plot_bins=list(
            run=TRUE,
            config=list(
            )
        ),
        compute_cn=list(
            run=TRUE,
            config=list(
                solver = "cpp",
                clones = "2,8",
                ghostprop = 0.3,
                tolerancerdr = 0.08,
                tolerancebaf = 0.04,
                seeds = 400,
                minprop = 0.03
            )
        ),
        plot_cn=list(
            run=TRUE,
            config=list(
            )
        )
    )
){
    return(configs)
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








