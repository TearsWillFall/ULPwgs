
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
    bin_tabix=build_default_tool_binary_list()$bin_tabix,
    vcf=NULL,
    index=TRUE,
    ...
){
   
    run_main=function(.env){

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

        .main.step=.main$steps[[fn_id]]
    
        if(index){
            .main.step$steps<-append(
                .main.step,
                    index_vcf_htslib(
                        bin_tabix=bin_tabix,
                        vcf=.main.step$out_files$bgzip_vcf,
                        index_format=index_format,
                        output_dir=out_file_dir,
                        tmp_dir=tmp_dir,
                        env_dir=env_dir,
                        batch_dir=batch_dir,
                        err_msg=err_msg,
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

#' Compress BED using BGZIP from HTSLIB
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

compress_bed_htslib=function(
    bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
    bin_tabix=build_default_tool_binary_list()$bin_tabix,
    bed=NULL,
    index=TRUE,
    ...
){
   
    run_main=function(.env){

        .this.env=environment()
        append_env(to=.this.env,from=.env)
        set_main(.env=.this.env)

        .main$out_files$bgzip_bed=paste0(out_file_dir,"/",input_id,".bed.gz")

        idx=""
        
        if(bgzip_index){
            idx=" -i "
            .main$out_files$bgzip_bed_idx=paste0(.main$out_files$bgzip_bed,".gzi")
        }

        .main$exec_code=paste(
            bin_bgzip,idx," -f -@ ",threads,
            " -c",input,">",.main$out_files$bgzip_bed)

        run_job(.env=.this.env)

        .main.step=.main$steps[[fn_id]]
    
        if(index){
            .main.step$steps<-append(
                .main.step,
                    index_bed_htslib(
                        bin_tabix=bin_tabix,
                        bed=.main.step$out_files$bgzip_bed,
                        index_format=index_format,
                        output_dir=out_file_dir,
                        tmp_dir=tmp_dir,
                        env_dir=env_dir,
                        batch_dir=batch_dir,
                        err_msg=err_msg,
                        verbose=verbose,
                        threads=threads,ram=ram,
                        executor_id=task_id
                )
            )

            .this.step=.main.step$steps$index_bed_htslib
            .main.step$out_files=append(.main.step$out_files,.this.step$out_files) 

        }

        .env$.main <- .main
        
    }

    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="bed"
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

    run_main=function(.env){

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



#' Uncompress BED using BGZIP from HTSLIB
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



uncompress_bed_htslib=function(
    bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
    bed=NULL,
    ...
){

    run_main=function(.env){

        .this.env=environment()
        append_env(to=.this.env,from=.env)
        set_main(.env=.this.env)

        .main$out_files$bed=paste0(out_file_dir,"/",input_id,".bed")

        .main$exec_code=paste(bin_bgzip," -f -d -@ ",threads,
            " -c ",bed," > ",.main$out_files$bed
        )


        run_job(.env=.this.env)

        .env$.main <- .main

    }

    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="bed"
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
    
    run_main=function(.env){
        

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




#' Index BED using TABIX from HTSLIB
#'
#'
#' @param bin_tabix Path to TABIX executable.
#' @param bed Path to the input VCF file.
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



index_bed_htslib=function(
    bin_tabix=build_default_tool_binary_list()$bin_tabix,
    bed=NULL,
    ...
){
    
    run_main=function(.env){
        

        .this.env=environment()
        append_env(to=.this.env,from=.env)
        set_main(.env=.this.env)

        if(index_format=="csi"){
            fmt=" -C "
            .main$out_files$bgzip_bed_idx=paste0(bed,".csi")
        }else if(index_format=="tbi"){
            fmt=""
            .main$out_files$bed_idx_tbi=paste0(bed,".tbi")
        }else{
            stop("Non valid index format provided. Valid formats are tbi/csi.")
        }
        
        .main$exec_code=paste0(bin_tabix, fmt," -f ",bed)


        run_job(.env=.this.env)

        .env$.main<-.main
    }

    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="bed"
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
    index=TRUE,
    ...
){


    run_main=function(.env){

        .this.env=environment()
        append_env(to=.this.env,from=.env)
        set_main(.env=.this.env)
        .main$steps[[fn_id]]<-.this.env
        .main.step=.main$steps[[fn_id]]

        if(compress){
            .main.step$steps<- append(
                .main$steps[[fn_id]]$steps,
                compress_vcf_htslib(
                    bin_bgzip=bin_bgzip,
                    vcf=vcf,
                    index=index,
                    index_format=index_format,
                    bgzip_index=bgzip_index,
                    output_dir=out_file_dir,
                    tmp_dir=tmp_dir,
                    env_dir=env_dir,
                    batch_dir=batch_dir,
                    err_msg=err_msg,
                    clean=clean,
                    verbose=verbose,
                    output_name=output_name,
                    threads=threads,
                    ram=ram,
                    executor_id=task_id
                )
            ) 
              .this.step=.main.step$steps$compress_vcf_htslib
              .main.step$out_files=append(.main.step$out_files,.this.step$out_files)
        }else{

                if(index){
                .main.step$steps<-append(
                    .main.step$steps,
                        index_vcf_htslib(
                        bin_tabix=bin_tabix,
                        vcf=vcf,
                        output_dir=out_file_dir,
                        tmp_dir=tmp_dir,
                        env_dir=env_dir,
                        batch_dir=batch_dir,
                        err_msg=err_msg,
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





#' Compress BED using BGZIP from HTSLIB
#'
#'
#' @param bin_bgzip Path to bgzip executable.
#' @param bin_tabix Path to TABIX executable.
#' @param bed Path to the input VCF file.
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


compress_and_index_bed_htslib=function(
    bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
    bin_tabix=build_default_tool_binary_list()$bin_tabix,
    bed=NULL,
    compress=TRUE,
    index=TRUE,
    ...
){


    run_main=function(.env){

        .this.env=environment()
        append_env(to=.this.env,from=.env)
        set_main(.env=.this.env)
        .main$steps[[fn_id]]<-.this.env
        .main.step=.main$steps[[fn_id]]

        if(compress){
            .main.step$steps<- append(
                .main$steps[[fn_id]]$steps,
                compress_bed_htslib(
                    bin_bgzip=bin_bgzip,
                    bed=bed,
                    index=index,
                    index_format=index_format,
                    bgzip_index=bgzip_index,
                    output_dir=out_file_dir,
                    tmp_dir=tmp_dir,
                    env_dir=env_dir,
                    batch_dir=batch_dir,
                    err_msg=err_msg,
                    clean=clean,
                    verbose=verbose,
                    output_name=output_name,
                    threads=threads,
                    ram=ram,
                    executor_id=task_id
                )
            ) 
              .this.step=.main.step$steps$compress_bed_htslib
              .main.step$out_files=append(.main.step$out_files,.this.step$out_files)
        }else{

                if(index){
                .main.step$steps<-append(
                    .main.step$steps,
                        index_bed_htslib(
                        bin_tabix=bin_tabix,
                        bed=bed,
                        output_dir=out_file_dir,
                        tmp_dir=tmp_dir,
                        env_dir=env_dir,
                        batch_dir=batch_dir,
                        err_msg=err_msg,
                        index_format=index_format,
                        verbose=verbose,
                        threads=threads,ram=ram,
                        executor_id=task_id
                    )
                )
            }

            
            .this.step=.main.step$steps$index_bed_htslib
            .main.step$out_files=append(.main.step$out_files,.this.step$out_files)
        }

        .env$.main<-.main


    }


    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="bed"
    )
  
    launch(.env=.base.env)

   

 

}
