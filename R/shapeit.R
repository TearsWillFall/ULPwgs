
#' Strelka wrapper for SNV variant calling
#'
#' This function wraps the STRELKA functions for variant calling
#' 
#' @param bin_shapeit Path to SHAPEIT binary file
#' @param haplotype Path to haplotype information from 1000G
#' @param legend Path to legend information from 1000G
#' @param map Path to map information from 1000G
#' @param vcf Path to VCF file
#' @export


phase_shapeit=function(
    bin_bcftools=build_default_tool_binary_list()$bin_bcftools,
    bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
    bin_tabix=build_default_tool_binary_list()$bin_tabix,
    bin_shapeit=build_default_tool_binary_list()$bin_shapeit,
    ref_panel=build_default_reference_list()$HG19$phasing$G1000$new$vcf,
    gmap=build_default_reference_list()$HG19$phasing$G1000$new$gmap,
    vcf=NULL,
    chr=c(1:22,"X"),
    scaffold=NULL,
    output_name=NULL,
    ...
){
    




     run_main=function(.env){

            .this.env=environment()
            append_env(to=.this.env,from=.env)
            set_main(.env=.this.env)
            
            out_file_dir=paste0(out_file_dir,"/",input_id,"/shapeit")

            .main$steps[[fn_id]]<-.this.env
            .main.step=.main$steps[[fn_id]]

        
            vcf=read_vcf(vcf=input,threads=threads)
            vcf$body=vcf$body %>% 
                only_gt_vcf() %>% 
                no_gt_to_pass_vcf() %>% 
                only_monoallel_vcf()
                
            vcf_tmp=vcf

            mclapply_os(chr,
                FUN=function(x){

                    this.chr=as.character(x)
                    
                    vcf_tmp$body=vcf$body %>% 
                        only_chr_vcf(chr=x)

                    .main.step$steps <- append(
                        .main.step$steps,
                        write_vcf(
                            bin_bgzip=bin_bgzip,
                            bin_tabix=bin_tabix,
                            vcf=vcf_tmp,
                            output_name=paste0(input_id,".",x,".filtered"),
                            output_dir=tmp_dir,
                            tmp_dir=tmp_dir,
                            env_dir=env_dir,
                            batch_dir=batch_dir,
                            err_msg=err_msg,
                            threads=2,
                            fn_id=x,
                            ram=ram,
                            executor=task_id
                            )
                        )
                .this.step=.main.step$steps[[paste0("write_vcf.",this.chr)]]
                .main.step$out_files$split_vcf[[this.chr]]=.this.step$out_files
                
      
                .main.step$steps<- append(
                .main.step$steps,
                    phase_chr_shapeit(
                        bin_bcftools=bin_bcftools,
                        bin_bgzip=bin_bgzip,
                        bin_tabix=bin_tabix,
                        bin_shapeit=bin_shapeit,
                        ref_panel=ref_panel,
                        gmap=gmap,
                        vcf=.main.step$out_files$split_vcf[[this.chr]]$bgzip_vcf,
                        chr=this.chr,
                        output_name=input_id,
                        scaffold=scaffold,
                        output_dir=out_file_dir,
                        tmp_dir=tmp_dir,
                        env_dir=env_dir,
                        fn_id=x,
                        batch_dir=batch_dir,
                        ram=ram,
                        verbose=verbose,
                        threads=2,
                        err_msg=err_msg,
                        executor_id=task_id
                    )
                )
            .this.step=.main.step$steps[[paste0("phase_chr_shapeit.",this.chr)]]
            .main.step$out_files$phased_vcf[[this.chr]]=.this.step$out_files
            return()
            },mc.cores=floor(threads/2))
           
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




#' Strelka wrapper for SNV variant calling
#'
#' This function wraps the STRELKA functions for variant calling
#' 
#' @param bin_shapeit Path to SHAPEIT binary file
#' @param haplotype Path to haplotype information from 1000G
#' @param legend Path to legend information from 1000G
#' @param map Path to map information from 1000G
#' @param vcf Path to VCF file
#' @export



phase_chr_shapeit=function(
            bin_bcftools=build_default_tool_binary_list()$bin_bcftools,
            bin_bgzip=build_default_tool_binary_list()$bin_bgzip,
            bin_tabix=build_default_tool_binary_list()$bin_tabix,
            bin_shapeit=build_default_tool_binary_list()$bin_shapeit,
            ref_panel=build_default_reference_list()$HG19$phasing$G1000$new$vcf,
            gmap=build_default_reference_list()$HG19$phasing$G1000$new$gmap,
            vcf=NULL,
            chr=c(1:22,"X"),
            scaffold=NULL,
            output_name=NULL,
            ...
        ){

                    run_main=function(.env){

                        .this.env=environment()
                        append_env(to=.this.env,from=.env)
                        set_main(.env=.this.env)

                        output_name=paste0(input_id,".",chr)

                        .main$steps[[fn_id]]<-.this.env
                        .main.step=.main$steps[[fn_id]]

                        .main.step$out_files$phased_vcf=paste0(out_file_dir,"/",input_id,".phased.vcf")


                        if(!is.null(ref_panel)){
                            add=paste0(" --reference ",ref_panel[grepl(paste0("chr",input,"."),ref_panel)]) 
                        }

                        if(!is.null(scaffold)){
                            add=paste0(" --scaffold ",scaffold)
                        }

                        .main$exec_code=paste(
                            bin_shapeit," --input ", vcf,
                            " --map ",gmap[grepl(paste0("chr",input,"."),gmap)],
                            add,
                            " --region ",input,
                            " --thread ",threads,
                            " --output ",.main.step$out_files$phased_vcf,
                            " --sequencing"
                        )

                        run_job(.env=.this.env)

                        .env$.main <- .main
                        
                    }

                    .base.env=environment()
                    list2env(list(...),envir=.base.env)
                    set_env_vars(
                    .env= .base.env,
                    vars="chr"
                    )
                
                    launch(.env=.base.env)

}

