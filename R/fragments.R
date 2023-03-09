
#' Get coverage at TFBS
#' This function calculates coverages based around gpos of TFBS
#' 
#' 

#' @param chrom Chromosome
#' @param pos Position
#' @param bam Path to BAM file.
#' @param region Number of bases around the genomic position
#' @export




get_coverage_tfbs=function(
    chrom=NULL,
    pos=NULL,
    bam=NULL,
    region=1000,
    ...
){



    run_main=function(
        .env
    ){
        .this.env=environment()
        append_env(to=.this.env,from=.env)
   
        set_main(.env=.this.env)

        .main$steps[[fn_id]]<-.this.env
        .main.step=.main$steps[[fn_id]]

        tfbs_id=paste0(chrom,"_",pos)
        tfbs=data.frame(chrom=chrom,pos=(pos-region):(pos+region))
        
        .main.step$steps<-append(
            .main.step$steps,
            get_coverage(
                bam=bam,
                gpos=tfbs,
                output_dir=out_file_dir,
                tmp_dir=tmp_dir,
                env_dir=env_dir,
                batch_dir=batch_dir,
                err_msg=err_msg,
                verbose=verbose,
                output_name=paste0(input_id,".",tfbs_id),
                threads=threads,
                ram=ram,
                executor_id=task_id
            )
        )
        .this.step=.main.step$steps$get_coverage
        .main.step$out_files$tfbs_coverage[[tfbs_id]]=.this.step$out_files$pileup

        .env$.main<-.main
    }


    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="bam"
    )
  
    launch(.env=.base.env)
}


#' Get coverage at TFBS
#' This function calculates coverages based around gpos of TFBS
#' 
#' 

#' @param gpos File with genomic position or data.frame
#' @param bam Path to BAM file.
#' @param region Number of bases around the genomic position
#' @export






evaluate_tf=function(
    gpos=NULL,
    bam=NULL,
    region=1000,
    ...
){

     run_main=function(
        .env
    ){
        .this.env=environment()
        append_env(to=.this.env,from=.env)
   
        set_main(.env=.this.env)

        .main$steps[[fn_id]]<-.this.env
        .main.step=.main$steps[[fn_id]]

        tf=ULPwgs::get_file_name(gpos)
        gpos=read_gpos(gpos=gpos,threads=threads) %>% 
        dplyr::distinct() %>% 
        dplyr::filter(pos>=1000)

        mclapply_os(1:nrow(gpos),function(pos){
                 id=paste0(gpos[pos,],collapse="_")
                .main.step$steps<-append(
                    .main.step$steps,
                    get_coverage_tfbs(
                            gpos=gpos[pos,],
                            bam=bam,
                            region=region,
                            output_dir=paste0(out_file_dir,"/",tf),
                            tmp_dir=tmp_dir,
                            env_dir=env_dir,
                            batch_dir=batch_dir,
                            output_name=paste0(input_id,".",tf),
                            err_msg=err_msg,
                            verbose=verbose,
                            threads=1,
                            ram=ram,
                            executor_id=task_id,
                            fn_id=id
                        )
                    )
                    .this.step=.main.step$steps[[paste0("get_coverage_tfbs.",id)]]
                    .main.step$out_files<-append(.main.step$out_files,this.step$out_files)
                    return()
        },mc.cores=threads)

        .env$.main<-.main
    }


    .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
        .env=.base.env,
        vars="bam"
    )

    launch(.env=.base.env)

}



