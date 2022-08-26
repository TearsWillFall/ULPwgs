
compress_vcf_htslib=function(
    bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
    vcf="",output_name="",bgzip_index=FALSE,
    output_dir=".",clean=FALSE,verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=4,mode="local",
    executor_id=make_unique_id("compressVCF"),
    task_name="compressVCF",time="48:0:0",
    update_time=60,wait=FALSE,hold=""
){
    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir)
    job=build_job(executor_id=executor_id,task_id=task_id)


    id=""
    if(output_name!=""){
        id=output_name
    }else{
        id=get_file_name(vcf)
    }

     out_file=paste0(out_file_dir,"/",id,".vcf.gz")

    idx=""
    out_index=""
    if(bgzip_index){
        idx=" -i "
        out_index=paste0(out_file,".gzi")
    }
    
    exec_code=paste(bin_bgzip, idx," -f -@ ",threads," -c",vcf,">",out_file)


    if(clean){
        exec_code=paste(exec_code," && rm",vcf)
    }

    if(mode=="batch"){
        out_file_dir2=set_dir(dir=out_file_dir,name="batch")
        batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
        threads=threads,output_dir=out_file_dir2)
        exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
    }

    if(verbose){
        print_verbose(job=job,arg=argg,exec_code=exec_code)
    }

    error=execute_job(exec_code=exec_code)
    
    if(error!=0){
        stop("htslib failed to run due to unknown error.
        Check std error for more information.")
    }

    job_report=build_job_report(
        job_id=job,
        executor_id=executor_id,
        exec_code=exec_code, 
        task_id=task_id,
        input_args = argg,
        out_file_dir=out_file_dir,
        out_files=list(
            compressed_vcf=out_file,
            compression_idx=out_index)
    )


    if(wait&&mode=="batch"){
        job_validator(job=job_report$job_id,time=update_time,
        verbose=verbose,threads=threads)
    }

    return(job_report)



}


uncompress_vcf_htslib=function(
    bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
    vcf="",output_dir=".",output_name="",
    clean=FALSE,verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=4,mode="local",
    executor_id=make_unique_id("uncompressVCF"),
    task_name="uncompressVCF",time="48:0:0",
    update_time=60,wait=FALSE,hold=""
){
    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir)
    job=build_job(executor_id=executor_id,task_id=task_id)


    id=""
    if(output_name!=""){
        id=output_name
    }else{
        id=get_file_name(vcf)
    }

    out_file=paste0(out_file_dir,"/",id,".vcf")

    exec_code=paste(bin_bgzip," -f -d -@ ",threads," -c ",vcf," > ",out_file)


    if(clean){
        exec_code=paste(exec_code," && rm",vcf)
    }

    if(mode=="batch"){
        out_file_dir2=set_dir(dir=out_file_dir,name="batch")
        batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
        threads=threads,output_dir=out_file_dir2)
        exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
    }

    if(verbose){
        print_verbose(job=job,arg=argg,exec_code=exec_code)
    }

    error=execute_job(exec_code=exec_code)
    
    if(error!=0){
        stop("htslib failed to run due to unknown error.
        Check std error for more information.")
    }

    job_report=build_job_report(
        job_id=job,
        executor_id=executor_id,
        exec_code=exec_code, 
        task_id=task_id,
        input_args = argg,
        out_file_dir=out_file_dir,
        out_files=list(
            uncompressed_vcf=out_file)
    )


    if(wait&&mode=="batch"){
        job_validator(job=job_report$job_id,time=update_time,
        verbose=verbose,threads=threads)
    }

    return(job_report)


}



index_vcf_htslib=function(
    bin_tabix=build_default_tool_binary_list()$bin_tabix,
    vcf="",index_format="tbi",verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=4,mode="local",
    executor_id=make_unique_id("indexVCF"),
    task_name="indexVCF",time="48:0:0",
    update_time=60,wait=FALSE,hold=""
){
    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir)
    job=build_job(executor_id=executor_id,task_id=task_id)

    
    if(index_format=="csi"){
        fmt=" -C "
        out_file=paste0(vcf,".csi")
    }else if(index_format=="tbi"){
        fmt=""
        out_file=paste0(vcf,".tbi")
    }
    
    exec_code=paste0(bin_bgzip, fmt," -f ",vcf)

    if(mode=="batch"){
        out_file_dir2=set_dir(dir=out_file_dir,name="batch")
        batch_code=build_job_exec(job=job,hold=hold,time=time,ram=ram,
        threads=threads,output_dir=out_file_dir2)
        exec_code=paste0("echo '. $HOME/.bashrc;",batch_config,";",exec_code,"'|",batch_code)
    }

    if(verbose){
        print_verbose(job=job,arg=argg,exec_code=exec_code)
    }

    error=execute_job(exec_code=exec_code)
    
    if(error!=0){
        stop("htslib failed to run due to unknown error.
        Check std error for more information.")
    }

    job_report=build_job_report(
        job_id=job,
        executor_id=executor_id,
        exec_code=exec_code, 
        task_id=task_id,
        input_args = argg,
        out_file_dir=out_file_dir,
        out_files=list(
            vcf_idx=paste0(vcf,".tbi"))
    )


    if(wait&&mode=="batch"){
        job_validator(job=job_report$job_id,time=update_time,
        verbose=verbose,threads=threads)
    }

    return(job_report)



}

compress_and_index_vcf_htslib=function(
    bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
    bin_tabix=build_default_tool_binary_list()$bin_tabix,
    vcf="",compress=TRUE,index=TRUE,index_format="tbi",
    bgzip_index=FALSE,output_dir=".",clean=FALSE,verbose=FALSE,
    batch_config=build_default_preprocess_config(),
    threads=1,ram=4,mode="local",
    executor_id=make_unique_id("CompressAndIndexVCF"),
    task_name="CompressAndIndexVCF",time="48:0:0",
    update_time=60,wait=FALSE,hold=""

){
    argg <- as.list(environment())
    task_id=make_unique_id(task_name)
    out_file_dir=set_dir(dir=output_dir)
    job=build_job(executor_id=executor_id,task=task_id)
    
    job_report=build_job_report(
        job_id=job,
        executor_id=executor_id, 
        exec_code=list(), 
        task_id=task_id,
        input_args=argg,
        out_file_dir=out_file_dir,
        out_files=list()
    )


    if(compress){
        job_report[["step"]][["compressVCF"]]<-compress_vcf_htslib(
            bin_bgzip=bin_bgzip,
            vcf=vcf,bgzip_index=bgzip_index,output_dir=out_file_dir,
            clean=clean,verbose=verbose,
            batch_config=batch_config,
            threads=threads,ram=ram,mode=mode,
            executor_id=task_id,
            time=time,
            hold=hold
        )
        vcf=job_report[["step"]][["compressVCF"]]$out_files$compressed_vcf
    }

    if(index){
        job_report[["step"]][["indexVCF"]]<-index_vcf_htslib(
            bin_tabix=bin_tabix,
            vcf=job_report[["step"]][["compressVCF"]]$out_files$compressed_vcf,
            index_format=index_format,verbose=verbose,
            batch_config=build_config,
            threads=threads,ram=ram,
            mode=mode,
            executor_id=task_id,
            time=time,
            hold=hold
        )
    }



    if(wait&&mode=="batch"){
        job_validator(job=job_report$job_id,time=update_time,
        verbose=verbose,threads=threads)
    }

    return(job_report)

}



