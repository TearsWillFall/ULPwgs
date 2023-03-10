
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
    gpos=NULL,
    bam=NULL,
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

     
        .main.step$steps<-append(
            .main.step$steps,
            get_coverage(
                bam=bam,
                gpos=gpos,
                output_dir=out_file_dir,
                tmp_dir=tmp_dir,
                env_dir=env_dir,
                batch_dir=batch_dir,
                err_msg=err_msg,
                verbose=verbose,
                output_name=input_id,
                threads=threads,
                gt="",
                ram=ram,
                executor_id=task_id
            )
        )
        .this.step=.main.step$steps$get_coverage
        .main.step$out_files$tfbs_coverage[[tf]]=.this.step$out_files$pileup

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


#' Get across TF
#' This function calculates coverages for TFBS of a single TF
#' 
#' 

#' @param gpos File with genomic position or data.frame
#' @param bam Path to BAM file.
#' @param region Number of bases around the genomic position
#' @param sep Separator
#' @param header Header in GPOS file. Defaults TRUE
#' @export


evaluate_tf=function(
    gpos=NULL,
    bam=NULL,
    region=1000,
    sep="\t",
    header=TRUE,
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
        gpos=read_gpos(gpos=gpos,threads=threads,header=header,sep=sep)

        gpos$body=gpos$body %>% 
            dplyr::distinct() %>% 
            dplyr::filter(pos>=1000)

        tfbs=gpos$body %>% rowwise() %>%
        mutate(tfbs_id=paste0(chrom,"_",pos),from=pos-1000,to=pos+1000) %>% 
        group_by(tfbs_id) %>%
        mutate(pos = map2(from, to, seq)) %>%
        unnest()%>%
        select(chrom,pos)

        tfbs[,c("ref","alt","gt")]<-"."
        
        .main.step$steps<-append(
            .main.step$steps,
            get_coverage_tfbs(
                    gpos=tfbs,
                    bam=bam,
                    output_dir=paste0(out_file_dir,"/",tf),
                    tmp_dir=tmp_dir,
                    env_dir=env_dir,
                    batch_dir=batch_dir,
                    output_name=paste0(input_id,".",tf),
                    err_msg=err_msg,
                    verbose=verbose,
                    threads=threads,
                    ram=ram,
                    executor_id=task_id,
                    fn_id=tf
                )
            )
            .this.step=.main.step$steps[[paste0("get_coverage_tfbs.",tf)]]
            .main.step$out_files<-append(.main.step$out_files,.this.step$out_files)
       
       
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



