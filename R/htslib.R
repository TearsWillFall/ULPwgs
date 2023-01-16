
#' Compress VCF using BGZIP from HTSLIB
#'
#'
#' @param bin_bgzip Path to bgzip executable.
#' @param vcf Path to the input VCF file.
#' @param output_name Output file name.
#' @param output_dir Path to the output directory.
#' @param bgzip_index Create BGZIP index for compressed file. Default FALSE.
#' @param clean Remove input VCF after completion. Default FALSE.
#' @param verbose Enables progress messages. Default False.
#' @param batch_config Default config to use in BATCH mode.
#' @param ram RAM memory to use in GB. Default 4.
#' @param tmp_dir Path to TMP directory. Default .
#' @param executor_id [OPTIONAL] Task executor name. Default "compressVCF
#' @param task_name [OPTIONAL] Task name. Default "compressVCF"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export

compress_vcf_htslib=function(
    bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
    vcf=NULL,
    ...
){
   
    build_main=function(.env){

        .this.env=environment()
        append_env(to=.this.env,from=.env)
        set_main(.env=.this.env)

        .main$out_files$bgzip_vcf=paste0(out_file_dir,"/",input_id,".vcf.gz")

        idx=""
        
        if(bgzip_index){
            idx=" -i "
            .main$out_files$bgzip_vcf_idx=paste0(.main$out_files$bgzip_vcf,".gzi")
        }

        .main$exec_code=paste(
            bin_bgzip,idx," -f -@ ",threads,
            " -c",input,">",.main$out_files$bgzip_vcf)

        run_job(.env=.this.env)
        .main.step=.main$steps[[fn]]
    


        if(index){
            .main$steps[[fn]]$steps<-append(
                .main$steps[[fn]]$steps,
                    index_vcf_htslib(
                    bin_tabix=bin_tabix,
                    vcf=.main.step$out_files$bgzip_vcf,
                    index_format=index_format,
                    verbose=verbose,
                    threads=threads,ram=ram,
                    executor_id=task_id
                )
            )

            .this.step=.main.step$steps$index_vcf_htslib
            .main.step$out_files=append(.main.step$out_files,.this.step$out_files) 

        }

        .env$.main <- .main
        
    }

    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="vcf"
    )
  
    launch(.env=.base.env)
    
}

#' Uncompress VCF using BGZIP from HTSLIB
#'
#'
#' @param bin_bgzip Path to bgzip executable.
#' @param vcf Path to the input VCF file.
#' @param output_name Output file name.
#' @param output_dir Path to the output directory.
#' @param clean Remove input VCF after completion. Default FALSE.
#' @param verbose Enables progress messages. Default False.
#' @param batch_config Default config to use in BATCH mode.
#' @param ram RAM memory to use in GB. Default 4.
#' @param tmp_dir Path to TMP directory. Default .
#' @param executor_id [OPTIONAL] Task executor name. Default "compressVCF
#' @param task_name [OPTIONAL] Task name. Default "compressVCF"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export



uncompress_vcf_htslib=function(
    bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
    vcf=NULL,
    ...
){

    build_main=function(.env){

        .this.env=environment()
        append_env(to=.this.env,from=.env)
        set_main(.env=.this.env)

        .main$out_files$vcf=paste0(out_file_dir,"/",input_id,".vcf")

        .main$exec_code=paste(bin_bgzip," -f -d -@ ",threads,
            " -c ",vcf," > ",.main$out_files$vcf
        )


        run_job(.env=.this.env)

        .env$.main <- .main

    }

    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="vcf"
    )
  
    launch(.env=.base.env)
  

}


#' Index VCF using TABIX from HTSLIB
#'
#'
#' @param bin_tabix Path to TABIX executable.
#' @param vcf Path to the input VCF file.
#' @param index_format Indexing method. Default tbi. Options [tbi,csi]
#' @param verbose Enables progress messages. Default False.
#' @param batch_config Default config to use in BATCH mode.
#' @param ram RAM memory to use in GB. Default 4.
#' @param tmp_dir Path to TMP directory. Default .
#' @param executor_id [OPTIONAL] Task executor name. Default "indexVCF"
#' @param task_name [OPTIONAL] Task name. Default "indexVCF"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export



index_vcf_htslib=function(
    bin_tabix=build_default_tool_binary_list()$bin_tabix,
    vcf=NULL,
    ...
){
    

    build_main=function(.env){
        

        .this.env=environment()
        append_env(to=.this.env,from=.env)
        set_main(.env=.this.env)

        if(index_format=="csi"){
            fmt=" -C "
            .main$out_files$bgzip_vcf_idx=paste0(vcf,".csi")
        }else if(index_format=="tbi"){
            fmt=""
            .main$out_files$vcf_idx_tbi=paste0(vcf,".tbi")
        }else{
            stop("Non valid index format provided. Valid formats are tbi/csi.")
        }
        
        .main$exec_code=paste0(bin_tabix, fmt," -f ",vcf)


        run_job(.env=.this.env)

        .env$.main<-.main
    }

    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="vcf"
    )
  
    launch(.env=.base.env)


}



#' Compress VCF using BGZIP from HTSLIB
#'
#'
#' @param bin_bgzip Path to bgzip executable.
#' @param bin_tabix Path to TABIX executable.
#' @param vcf Path to the input VCF file.
#' @param compress Compress VCF file. Default TRUE.
#' @param index Index VCF file. Default TRUE.
#' @param index_format VCF index format. Default tbi. Options [tbi,cbi].
#' @param bgzip_index Create BGZIP index for compressed file. Default FALSE
#' @param output_name Output file name.
#' @param output_dir Path to the output directory.
#' @param clean Remove input VCF after completion. Default FALSE.
#' @param verbose Enables progress messages. Default False.
#' @param batch_config Default config to use in BATCH mode.
#' @param ram RAM memory to use in GB. Default 4.
#' @param tmp_dir Path to TMP directory. Default .
#' @param executor_id [OPTIONAL] Task executor name. Default "compressVCF
#' @param task_name [OPTIONAL] Task name. Default "compressVCF"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] HOld job until job is finished. Job ID. 
#' @export


compress_and_index_vcf_htslib=function(
    bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
    bin_tabix=build_default_tool_binary_list()$bin_tabix,
    vcf=NULL,
    compress=TRUE,
    ...
){


    build_main=function(.env){

        .this.env=environment()
        append_env(to=.this.env,from=.env)
        set_main(.env=.this.env)
        .main$steps[[fn]]<-.this.env
        .main.step=.main$steps[[fn]]

        if(compress){
            .main$steps[[fn]]$steps<- append(
                .main$steps[[fn]]$steps,
                compress_vcf_htslib(
                    bin_bgzip=bin_bgzip,
                    vcf=vcf,
                    index=index,
                    index_format=index_format,
                    bgzip_index=bgzip_index,
                    output_dir=out_file_dir,
                    clean=clean,
                    verbose=verbose,
                    output_name=output_name,
                    threads=threads,ram=ram,
                    executor_id=task_id
                )
            ) 
        }else{

                if(index){
                .main$steps[[fn]]$steps<-append(
                    .main$steps[[fn]]$steps,
                        index_vcf_htslib(
                        bin_tabix=bin_tabix,
                        vcf=vcf,
                        index_format=index_format,
                        verbose=verbose,
                        threads=threads,ram=ram,
                        executor_id=task_id
                    )
                )
            }

            
            .this.step=.main.step$steps$index_vcf_htslib
            .main.step$out_files=append(.main.step$out_files,.this.step$out_files) 

            

        }

        

        .env$.main<-.main


    }


    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="vcf"
    )
  
    launch(.env=.base.env)

   

 

}



