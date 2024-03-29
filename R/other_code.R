
# get_tf_from_cnvkit=function(
#   cnvkit_data=NULL,
#   tf=NULL,
#   threads=1
# ){
#   jobs=length(tf)
#   hit_tfbs=parallel::mclapply(1:jobs,FUN=function(x){
#     cnr=cnvkit_data$cnr
#     cns=cnvkit_data$cns
#     depth=median(cns$depth)
#     lst=list(missing_tfbs=NA,hit_tfbs=NA)
#     tfbs=tf[[x]]
#     overlaps=GenomicRanges::findOverlaps(tfbs,cnr)
#     missing_tfbs=tfbs[-as.data.frame(overlaps)$queryHits]
#     lst$missing_tfbs=missing_tfbs
#     chromosomes=as.character(unique(GenomeInfoDb::seqnames(tfbs)))
#     chunks=length(chromosomes)
#     hit_tfbs=lapply(chromosomes,FUN=function(chr){
#         tfbs_tmp=tfbs[GenomeInfoDb::seqnames(tfbs)==chr]
#         cnr_tmp=cnr[GenomeInfoDb::seqnames(cnr)==chr]
#         cns_tmp=cns[GenomeInfoDb::seqnames(cns)==chr]

#         tmp_ranges=as.data.frame(IRanges::ranges(cnr_tmp))
#         cnr_tmp$cnr_pos=(tmp_ranges$start+tmp_ranges$end)/2
#         cnr_tmp$cnr_bin_size=tmp_ranges$width

#         invisible(lapply(5:1,FUN=function(x){
#             sol=dplyr::lag(cnr_tmp$rid,n=x)
#             S4Vectors::mcols(cnr_tmp)[paste0("cnr_rid_left_",x)]<<-ifelse(
#             is.na(sol),cnr_tmp$rid,sol)
#         }))
#         cnr_tmp$cnr_rid_central=cnr_tmp$rid
        
#         invisible(lapply(1:5,FUN=function(x){
#           sol=dplyr::lead(cnr_tmp$rid,n=x)
#           S4Vectors::mcols(cnr_tmp)[paste0("cnr_rid_right_",x)]<<-ifelse(
#              is.na(sol),cnr_tmp$rid,sol)
#         }))

      
#         invisible(lapply(5:1,FUN=function(x){
#           sol=dplyr::lag(cnr_tmp$log2,n=x)
#           S4Vectors::mcols(cnr_tmp)[paste0("cnr_tfbs_left_log2_",x)]<<-ifelse(
#             is.na(sol),cnr_tmp$log2,sol)
#         }))
#         cnr_tmp$cnr_tfbs_central_log2=cnr_tmp$log2

#         invisible(lapply(1:5,FUN=function(x){
#           sol=dplyr::lead(cnr_tmp$log2,n=x)
#           S4Vectors::mcols(cnr_tmp)[,paste0("cnr_tfbs_right_log2_",x)]<<-ifelse(
#             is.na(sol),cnr_tmp$log2,sol
#           )
#         }))

#         invisible(lapply(5:1,FUN=function(x){
#             sol=dplyr::lag(cnr_tmp$depth,n=x)
#             S4Vectors::mcols(cnr_tmp)[,paste0("cnr_tfbs_left_depth_",x)]<<-ifelse(
#             is.na(sol),cnr_tmp$log2,sol
#           )
#         }))

#         cnr_tmp$cnr_tfbs_central_depth=cnr_tmp$depth

#         invisible(lapply(1:5,FUN=function(x){
#           sol=dplyr::lead(cnr_tmp$depth,n=x)
#           S4Vectors::mcols(cnr_tmp)[,paste0("cnr_tfbs_right_depth_",x)]<<-ifelse(
#             is.na(sol),cnr_tmp$depth,sol
#           )
#          }))

#         cnr_tmp=cnr_tmp[,grepl("cnr",names(S4Vectors::mcols(cnr_tmp)))]


#         hit_tfbs=plyranges::join_overlap_left(
#           tfbs_tmp,cnr_tmp,suffix=NULL)

#         tmp_ranges=as.data.frame(IRanges::ranges(cns_tmp))
#         cns_tmp$cns_pos=(tmp_ranges$start+tmp_ranges$end)/2
#         cns_tmp$cns_seg_size=tmp_ranges$width
#         cns_tmp$cns_seg_left_log2=dplyr::lag(cns_tmp$log2)
#         cns_tmp$cns_seg_central_log2=cns_tmp$log2
#         cns_tmp$cns_seg_right_log2=dplyr::lead(cns_tmp$log2)
#         cns_tmp$cns_sid_left=dplyr::lag(cns_tmp$sid)
#         cns_tmp$cns_sid_central=cns_tmp$sid
#         cns_tmp$cns_sid_right=dplyr::lead(cns_tmp$sid)

#         cns_tmp=cns_tmp[,grepl("cns",names(S4Vectors::mcols(cns_tmp)))]
#         hit_tfbs=plyranges::join_overlap_left(
#           hit_tfbs,cns_tmp,suffix=NULL
#         )
#           return(hit_tfbs)
#         }
#     )
#     hit_tfbs=plyranges::bind_ranges(hit_tfbs)
#     hit_tfbs$sample_depth=depth
#     lst$hit_tfbs=hit_tfbs
#     return(lst)
#     },mc.cores=ifelse(jobs==1,1,threads))
#   names(hit_tfbs)=names(tf)
#   return(hit_tfbs)
# }



# get_tf_from_sample=function(
#   cnvkit_data=NULL,
#   tf=NULL,
#   threads=1
# ){
#   jobs=length(cnvkit_data)
#   dat=parallel::mclapply(1:jobs,function(x){
#     sol=get_tf_from_cnvkit(
#         cnvkit_data=cnvkit_data[[x]],
#         tf=tf,
#         threads=ifelse(jobs==1,threads,1))
#   },mc.cores=threads)
#   names(dat)=names(cnvkit_data)
#   return(dat)
# }




# import_cnvkit_data=function(
#     cnr=NULL,
#     cns=NULL,
#     threads=1,
#     chromosomes=c(1:22,"X","Y")
# ){
#   jobs=length(cnr)
#   dat=parallel::mclapply(1:jobs,FUN=function(x){
#     lst=list(cnr=NA,cns=NA)
#     cnr=read.table(cnr[[x]],header=TRUE,sep="\t")
#     cnr=cnr[cnr$chromosome %in% chromosomes,]
#     cnr$rid=paste0(cnr$chromosome,":",cnr$start,"-",cnr$end)
#     cnr=GenomicRanges::makeGRangesFromDataFrame(cnr,keep.extra.columns=TRUE)
#     lst$cnr=cnr
#     cns=read.table(cns[[x]],header=TRUE,sep="\t")
#     cns=cns[cns$chromosome %in% chromosomes,]
#     cns$sid=paste0(cns$chromosome,":",cns$start,"-",cns$end)
#     cns$size=cns$end-cns$start
#     cns=GenomicRanges::makeGRangesFromDataFrame(cns,keep.extra.columns=TRUE)
#     lst$cns=cns
#     return(lst)
#   },mc.cores=threads)
#   names(dat)=Vectorize(ULPwgs::get_file_name)(cnr)
#   return(dat)
# }


# import_tf_data=function(
#     tf=NULL,
#     threads=1,
#     chromosomes=c(1:22,"X","Y")
# ){
#   jobs=length(tf)
#   dat=parallel::mclapply(1:jobs,FUN=function(x){
#     tfbs=ULPwgs::read_gpos(tf[[x]],sep=" ",rename=FALSE,sort=FALSE)$body
#     tfbs=tfbs[tfbs$chrom %in% chromosomes,]
#     tfbs=GenomicRanges::makeGRangesFromDataFrame(tfbs,start.field="pos",end.field="pos",keep.extra.columns=TRUE)
#     return(tfbs)
#   },mc.cores=threads)
#   names(dat)=Vectorize(ULPwgs::get_file_name)(tf)
#   return(dat)
# }

# transform_tf_data=function(tf_data=NULL){
#   sol=lapply(names(tf_data),FUN=function(x){dat=tf_data[[x]]$GRHL2$hit_tfbs;dat$id=names(tf_data[x]);dat})
#   sol=plyranges::bind_ranges(sol)
#   sol=as.data.frame(sol)
#   return(sol)
# }


# sol_data_long=sol_data %>% select(starts_with("cnr_tfbs"),id,seqnames,gid) %>% pivot_longer(cols=!id:gid)
# sol_data_long=sol_data_long[grepl("log2",sol_data_long$name),]
# sol_data_long$pos=ifelse(grepl("left",sol_data_long$name),"left",ifelse(grepl("right",sol_data_long$name),"right","central"))
# sol_data_long=sol_data_long %>% group_by(id,gid) %>% mutate(dist_to_tfbs=-5:5,abs_dist_to_tfbs=abs(-5:5))
# sol_data_long=sol_data_long %>% group_by(id,gid,abs_dist_to_tfbs) %>% mutate(abs_value=median(value,na.rm=TRUE))




#' Filter BAM file by size using samtools
#'
#'
#' @param bam Path to the input file with the sequence.
#' @param bin_samtools Path to samtools executable. Default path tools/samtools/samtools.
#' @param mapq Mapping quality of the read to analyze. Default 60.
#' @param mapq Flags of the reads to read. Default c(99, 147, 83, 163)
#' @param region Genomic region to search
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param threads Number of threads. Default 3
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Job EXECUTOR ID. Default "mardupsGATK"
#' @param task_name Name of the task. Default "mardupsGATK"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] Hold job until job is finished. Job ID. 
#' @export


get_tf_from_cnvkit=function(
  cnr=NULL,
  cns=NULL,
  tf=NULL,
  ...
){

   run_main=function(
    .env
  ){

    .this.env=environment()
    append_env(to=.this.env,from=.env)
    set_main(.env=.this.env)

    fpath <- system.file("R", paste0(fn,".R"), package=ns)
    system(paste0("chmod 755 ",fpath))
    fpath<-paste0("Rscript ",fpath)

    tf_name=get_file_name(input)
    .main$out_files[[tf_name]]$hits<-paste0(
      out_file_dir,"/",input_id,".", tf_name,".hits.txt"
    )
    .main$out_files[[tf_name]]$miss<-paste0(
      out_file_dir,"/",input_id,".", tf_name,".miss.txt"
    )
    .main$exec_code<-paste0(
        fpath,
        " -r ",cnr,
        " -s ",cns,
        " -t ",input,
        " -o ",paste0(out_file_dir,"/",input_id,".",tf_name)
    )
    run_job(.env=.this.env)
    .env$.main<-.main

  }

  .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="tf"
    )

    launch(.env=.base.env)
}




#' Filter BAM file by size using samtools
#'
#'
#' @param bam Path to the input file with the sequence.
#' @param bin_samtools Path to samtools executable. Default path tools/samtools/samtools.
#' @param mapq Mapping quality of the read to analyze. Default 60.
#' @param mapq Flags of the reads to read. Default c(99, 147, 83, 163)
#' @param region Genomic region to search
#' @param output_dir Path to the output directory.
#' @param verbose Enables progress messages. Default False.
#' @param threads Number of threads. Default 3
#' @param mode [REQUIRED] Where to parallelize. Default local. Options ["local","batch"]
#' @param executor_id Job EXECUTOR ID. Default "mardupsGATK"
#' @param task_name Name of the task. Default "mardupsGATK"
#' @param time [OPTIONAL] If batch mode. Max run time per job. Default "48:0:0"
#' @param update_time [OPTIONAL] If batch mode. Job update time in seconds. Default 60.
#' @param wait [OPTIONAL] If batch mode wait for batch to finish. Default FALSE
#' @param hold [OPTIONAL] Hold job until job is finished. Job ID. 
#' @export


get_tf_from_sample=function(
  cnvkit_data=NULL,
  tf_data=NULL,
  patient_id=NULL,
  ...
){

   run_main=function(
    .env
  ){

    .this.env=environment()
    append_env(to=.this.env,from=.env)
    out_file_dir=set_dir(
        out_file_dir,
        name=paste0(patient_id,"/tf_reports/",get_file_name(input[[1]]$cns))
    )

    set_main(.env=.this.env)
    
    run_mode="local_parallel"
    if(mode=="local_parallel"){
      run_mode="local"
    }
   

    .main$steps[[fn_id]]<-.this.env
    .main.step=.main$steps[[fn_id]]


    .main.step$steps <-append(
    .main.step$steps,
      get_tf_from_cnvkit(
        cnr=input[[1]]$cnr,
        cns=input[[1]]$cns,
        tf=tf_data,
        output_name=get_file_name(input[[1]]$cns),
        output_dir=out_file_dir,
        tmp_dir=tmp_dir,
        env_dir=env_dir,
        batch_dir=batch_dir,
        err_msg=err_msg,
        mode=run_mode,
        verbose=verbose,
        threads=threads,
        ram=ram,
        executor_id=task_id
      )
    )
    .this.step=.main.step$steps$get_tf_from_cnvkit
    .main.step$out_files=append(.main.step$out_files,.this.step$out_files)
    .env$.main<-.main
  }

  .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
      .env= .base.env,
      vars="cnvkit_data"
    )

    launch(.env=.base.env)
}


# options(scipen=999)
# sol=lapply(-1,
#   FUN=function(y){
#     lapply(unique(dat$id)[918:1056],
#       FUN=function(x){
#         dat_tmp=dat %>% filter(id!=x);
#         ovlps=findOverlaps(dat_tmp,dat %>% filter(id==x),maxgap=y);
#         sol=dat_tmp[as.data.frame(ovlps)$queryHits];
#         sol$tf_hit=x;sol$dist=y;data.table::fwrite(file=paste0("overlaps/tf_",x,".size_",y,".overlaps"),as.data.frame(sol));return()})})


# sol_wider=dat %>% 
# as.data.frame()%>% 
# group_by(seqnames,start,end,id,dist)%>% 
# summarise(N=length(unique(tf_hit))) %>% 
# arrange(dist,desc(N)) 
# dat_only_dt=dplyr::left_join(as.data.frame(dat_only),sol_wider)
# dat_only_dt=dat_only_dt %>% mutate_if(is.numeric,dplyr::coalesce,0)
# dat_only_dt_longer=dat_only_dt %>% pivot_longer(`-1`:`1000000`)