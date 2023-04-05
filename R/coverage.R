
#' @export

generate_pga=function(
    cnr=NULL,
    gain=seq(0,1,by=0.01),
    loss=seq(-1,0,by=0.01),
    chrom=c(1:22),
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

        .main$out_files$pga=paste0(out_file_dir,"/",input_id,".pga")
        all_info=calculate_pga(cnr=input,id=input_id,gain=gain,loss=loss,chrom=chrom,threads=threads)

        write.table(
            all_info,
            file=.main$out_files$pga,
            sep="\t",
            col.names=TRUE,
            row.names=FALSE,
            quote=FALSE
        )

        .env$.main <- .main
    }

     .base.env=environment()
    list2env(list(...),envir=.base.env)
    set_env_vars(
        .env= .base.env,
        vars="cnr"
    )

    launch(.env=.base.env)



}

#' @export


calculate_pga=function(
    cnr=NULL,
    id=NULL,
    gain=seq(0,1,by=0.1),
    loss=seq(-1,0,by=0.1),
    chrom=c(1:22),
    threads=1
){


    dat=read.table(cnr,sep="\t",header=TRUE)
    dat=dat[dat$chromosome %in% chrom,]
    ploidy=ploidy_from_cnr(cnr=dat,chrom=chrom)
    all_info=mclapply_os(loss,FUN=function(x){
        info=lapply(gain,FUN=function(y){
            dat_tmp=dat
           
    
            dat_tmp$width=dat_tmp$end-dat_tmp$start
            dat_tmp$bin_type=ifelse(
                grepl("Antitarget",dat_tmp$gene),
                "Antitarget","Target"
            )
            dat_tmp$log2=dat_tmp$log2
            dat_tmp$TYPE=ifelse(dat_tmp$log2>=y,"GAIN",
                ifelse(dat_tmp$log2<=x,"LOSS","WT"
            ))


            dat_rslt=dat_tmp%>% 
                dplyr::group_by(TYPE) %>% 
                dplyr::summarise(
                    N=dplyr::n(),
                    N_target=length(log2[bin_type=="Target"]),
                    N_antitarget=length(log2[bin_type=="Antitarget"]),
                    G_size=sum(width),
                    G_size_target=sum(width[bin_type=="Target"]),
                    G_size_antitarget=sum(width[bin_type=="Antitarget"]),
                    Weight=median(weight),
                    Weight_target=median(as.numeric(weight[bin_type=="Target"]),na.rm=TRUE),
                    Weight_antitarget=median(as.numeric(weight[bin_type=="Antitarget"]),na.rm=TRUE)
                ) %>% 
                dplyr::ungroup() %>% 
                dplyr::mutate(
                    FRACTION=N/sum(N),
                    FRACTION_target=N_target/sum(N_target),
                    FRACTION_antitarget=N_antitarget/sum(N_antitarget),
                    FRACTION_genome=G_size/sum(G_size),
                    FRACTION_genome_target=G_size_target/sum(G_size_target),
                    FRACTION_genome_antitarget=G_size_antitarget/sum(G_size_antitarget),
                    T_N=sum(N),
                    T_N_target=sum(N_target),
                    T_N_antitarget=sum(N_antitarget),
                    T_G_size=sum(G_size),
                    T_G_size_target=sum(G_size_target),
                    T_G_size_antitarget=sum(G_size_antitarget)
                )

            dat_rslt$gain=y
            dat_rslt$loss=x
            return(dat_rslt)
        }

        )
        info=dplyr::bind_rows(info)
        return(info)
        
        },mc.cores=threads

    )
    all_info=dplyr::bind_rows(all_info)
    all_info$ploidy_all=ploidy$ploidy_all
    all_info$ploidy_target=ploidy$ploidy_target
    all_info$ploidy_antitarget=ploidy$ploidy_antitarget
    
    all_info=all_info %>% 
        dplyr::arrange(gain,dplyr::desc(loss)) %>% 
        dplyr::filter(TYPE!="WT")
    all_info$id=id

    return(all_info)

}



#' @export

extract_pga=function(tumour=NULL,normal=NULL){

    tumour=read.table(tumour,sep="\t",header=TRUE)
    normal=read.table(normal,sep="\t",header=TRUE)
    tumour=dplyr::left_join(
        tumour,normal,
        by=c("TYPE"="TYPE","loss"="loss","gain"="gain"),
        suffix=c(".tumour",".normal")
    )
    tumour=tumour %>%
        dplyr::rowwise()%>%
        dplyr::mutate(
            FRACTION=FRACTION.tumour-FRACTION.normal,
            FRACTION_target=FRACTION_target.tumour-FRACTION_target.normal,
            FRACTION_antitarget=FRACTION_antitarget.tumour-FRACTION_antitarget.normal,
            FRACTION_genome=FRACTION_genome.tumour-FRACTION_genome.normal,
            FRACTION_genome_target=FRACTION_genome_target.tumour-FRACTION_genome_target.normal,
            FRACTION_genome_antitarget=FRACTION_genome_antitarget.tumour-FRACTION_genome_antitarget.normal
    ) %>% 
    dplyr::group_by(TYPE) %>%
    dplyr::mutate(
        corrFRACTION=FRACTION[gain==0&loss==0],
        corrFRACTION_target=FRACTION_target[gain==0&loss==0],
        corrFRACTION_antitarget=FRACTION_antitarget[gain==0&loss==0],
        corrFRACTION_genome=FRACTION_genome[gain==0&loss==0],
        corrFRACTION_genome_target=FRACTION_genome_target[gain==0&loss==0],
        corrFRACTION_genome_antitarget=FRACTION_genome_antitarget[gain==0&loss==0]
    ) %>% dplyr::ungroup() %>%
        dplyr::mutate(
            fcorrFRACTION=FRACTION-corrFRACTION,
            fcorrFRACTION_target=FRACTION_target-corrFRACTION_target,
            fcorrFRACTION_antitarget=FRACTION_antitarget-corrFRACTION_antitarget,
            fcorrFRACTION_genome=FRACTION_genome-corrFRACTION_genome,
            fcorrFRACTION_genome_target=FRACTION_genome_target-corrFRACTION_genome_target,
            fcorrFRACTION_genome_antitarget=FRACTION_genome_antitarget-corrFRACTION_genome_antitarget
        ) %>%dplyr::group_by(id.tumour,loss,gain) %>%
         dplyr::mutate(
            sle=gain-loss,
            sumfcorrFRACTION=sum(fcorrFRACTION),
            sumfcorrFRACTION_target=sum(fcorrFRACTION_target),
            sumfcorrFRACTION_antitarget=sum(fcorrFRACTION_antitarget),
            sumfcorrFRACTION_genome=sum(fcorrFRACTION_genome),
            sumfcorrFRACTION_genome_target=sum(fcorrFRACTION_genome_target),
            sumfcorrFRACTION_genome_antitarget=sum(fcorrFRACTION_genome_antitarget)
        )









    bin_all=tumour %>%dplyr::ungroup() %>%
        dplyr::filter(sumfcorrFRACTION==max(sumfcorrFRACTION)) %>%
        dplyr::filter(sle==min(sle)) %>%
        dplyr::select(
            TYPE,
            loss,
            gain,
            T_N.tumour,
            T_N_target.tumour,
            T_N_antitarget.tumour,
            T_G_size.tumour,
            T_G_size_target.tumour,
            T_G_size_antitarget.tumour,
             id.tumour,
            FRACTION,
            corrFRACTION,
            fcorrFRACTION,
            sumfcorrFRACTION
        ) %>% 
        dplyr::mutate(method="bin",source="all") %>%
        dplyr::rename(PGA=sumfcorrFRACTION)
    bin_target=tumour %>%dplyr::ungroup() %>%
        dplyr::filter(sumfcorrFRACTION_target==max(sumfcorrFRACTION_target)) %>%
        dplyr::filter(sle==min(sle)) %>%
        dplyr::select(
            TYPE,
            loss,
            gain,
            T_N.tumour,
            T_N_target.tumour,
            T_N_antitarget.tumour,
            T_G_size.tumour,
            T_G_size_target.tumour,
            T_G_size_antitarget.tumour,
             id.tumour,
            FRACTION_target,
            corrFRACTION_target,
            fcorrFRACTION_target,
            sumfcorrFRACTION_target
        ) %>% 
        dplyr::mutate(method="bin",source="target") %>%
        dplyr::rename(
            PGA=sumfcorrFRACTION_target,
            FRACTION=FRACTION_target,
            corrFRACTION=corrFRACTION_target,
            fcorrFRACTION=fcorrFRACTION_target
        )
    bin_antitarget=tumour %>%dplyr::ungroup() %>%
        dplyr::filter(sumfcorrFRACTION_antitarget==max(sumfcorrFRACTION_antitarget)) %>%
        dplyr::filter(sle==min(sle))%>%
        dplyr::select(
            TYPE,
            loss,
            gain,
            T_N.tumour,
            T_N_target.tumour,
            T_N_antitarget.tumour,
            T_G_size.tumour,
            T_G_size_target.tumour,
            T_G_size_antitarget.tumour,
             id.tumour,
            FRACTION_antitarget,
            corrFRACTION_antitarget,
            fcorrFRACTION_antitarget,
            sumfcorrFRACTION_antitarget
        ) %>% 
        dplyr::mutate(method="bin",source="antitarget") %>%
        dplyr::rename(
            PGA=sumfcorrFRACTION_antitarget,
            FRACTION=FRACTION_antitarget,
            corrFRACTION=corrFRACTION_antitarget,
            fcorrFRACTION=fcorrFRACTION_antitarget
        )
    base_all=tumour %>%dplyr::ungroup() %>%
        dplyr::filter(sumfcorrFRACTION_genome==max(sumfcorrFRACTION_genome)) %>%
          dplyr::filter(sle==min(sle)) %>%
        dplyr::select(
            TYPE,
            loss,
            gain,
            T_N.tumour,
            T_N_target.tumour,
            T_N_antitarget.tumour,
            T_G_size.tumour,
            T_G_size_target.tumour,
            T_G_size_antitarget.tumour,
            id.tumour,
            FRACTION_genome,
            corrFRACTION_genome,
            fcorrFRACTION_genome,
            sumfcorrFRACTION_genome
        ) %>% 
        dplyr::mutate(method="base",source="all") %>%
        dplyr::rename(
            PGA=sumfcorrFRACTION_genome,
            FRACTION=FRACTION_genome,
            corrFRACTION=corrFRACTION_genome,
            fcorrFRACTION=fcorrFRACTION_genome
        )
    base_target=tumour %>%dplyr::ungroup() %>%
        dplyr::filter(sumfcorrFRACTION_genome_target==max(sumfcorrFRACTION_genome_target)) %>%
        dplyr::filter(sle==min(sle)) %>%
        dplyr::select(
            TYPE,
            loss,
            gain,
            T_N.tumour,
            T_N_target.tumour,
            T_N_antitarget.tumour,
            T_G_size.tumour,
            T_G_size_target.tumour,
            T_G_size_antitarget.tumour,
            id.tumour,
            FRACTION_genome_target,
            corrFRACTION_genome_target,
            fcorrFRACTION_genome_target,
            sumfcorrFRACTION_genome_target
        ) %>% 
        dplyr::mutate(method="base",source="target") %>%
        dplyr::rename(
            PGA=sumfcorrFRACTION_genome_target,
            FRACTION=FRACTION_genome_target,
            corrFRACTION=corrFRACTION_genome_target,
            fcorrFRACTION=fcorrFRACTION_genome_target
        )
    base_antitarget=tumour %>%dplyr::ungroup() %>%
        dplyr::filter(sumfcorrFRACTION_genome_antitarget==max(sumfcorrFRACTION_genome_antitarget)) %>%
        dplyr::filter(sle==min(sle)) %>%
        dplyr::select(
            TYPE,
            loss,
            gain,
            T_N.tumour,
            T_N_target.tumour,
            T_N_antitarget.tumour,
            T_G_size.tumour,
            T_G_size_target.tumour,
            T_G_size_antitarget.tumour,
            id.tumour,
            FRACTION_genome_antitarget,
            corrFRACTION_genome_antitarget,
            fcorrFRACTION_genome_antitarget,
            sumfcorrFRACTION_genome_antitarget
        ) %>% 
        dplyr::mutate(method="base",source="antitarget") %>%
        dplyr::rename(
            PGA=sumfcorrFRACTION_genome_antitarget,
            FRACTION=FRACTION_genome_antitarget,
            corrFRACTION=corrFRACTION_genome_antitarget,
            fcorrFRACTION=fcorrFRACTION_genome_antitarget
        )
    dat_summ=dplyr::bind_rows(bin_all,bin_target,bin_antitarget,base_all,base_target,base_antitarget)
    return(dat_summ)

}