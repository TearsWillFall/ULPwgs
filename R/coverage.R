
#' @export

generate_pga=function(
    cnr=NULL,
    gain=seq(0,1,by=0.01),
    loss=seq(-1,0,by=0.01),
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
        dat=read.table(input,sep="\t",header=TRUE)
        
    
        all_info=mclapply_os(loss,FUN=function(x){
            info=lapply(gain,FUN=function(y){
                dat_tmp=dat
                dat_tmp$width=dat_tmp$end-dat_tmp$start
                dat_tmp$TYPE=ifelse(dat_tmp$log2>=y,"GAIN",
                    ifelse(dat_tmp$log2<=x,"LOSS","WT"
                ))
        
                dat_rslt=dat_tmp%>% 
                    dplyr::group_by(TYPE) %>% 
                    dplyr::summarise(
                        N=dplyr::n(),
                        N_target=length(log2[gene!="Antitarget"]),
                        N_antitarget=length(log2[gene=="Antitarget"]),
                        G_size=sum(width),
                        G_size_target=sum(width[gene!="Antitarget"]),
                        G_size_antitarget=sum(width[gene=="Antitarget"]),
                        Weight=median(weight),
                        Weight_target=median(as.numeric(weight[gene!="Antitarget"]),na.rm=TRUE),
                        Weight_antitarget=median(as.numeric(weight[gene=="Antitarget"]),na.rm=TRUE)
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
        
        all_info=all_info %>% 
            dplyr::arrange(gain,dplyr::desc(loss)) %>% 
            dplyr::filter(TYPE!="WT")
        all_info$id=input_id
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
    )

    bin_all=tumour %>% dplyr::group_by(TYPE) %>% 
        dplyr::mutate(BASE=FRACTION[gain==0&loss==0]) %>% 
        dplyr::filter(FRACTION==max(FRACTION)) %>% 
        dplyr::select(TYPE,loss,gain,BASE,FRACTION) %>% 
        dplyr::ungroup()%>%
        dplyr::rename_with(~sub("_.*","",.x),starts_with("FRACTION")) %>%
        dplyr::mutate(PGA=sum(FRACTION),method="bin",source="all")
    bin_target=tumour %>% dplyr::group_by(TYPE)%>% 
        dplyr::mutate(BASE=FRACTION_target[gain==0&loss==0]) %>% 
        dplyr::filter(FRACTION_target==max(FRACTION_target)) %>% 
        dplyr::select(TYPE,loss,gain,BASE,FRACTION_target)%>% 
        dplyr::ungroup()%>%
        dplyr::rename_with(~sub("_.*","",.x),starts_with("FRACTION"))
        dplyr::mutate(PGA=sum(FRACTION),method="bin",source="target")
    bin_antitarget=tumour %>% dplyr::group_by(TYPE)%>% 
        dplyr::mutate(BASE=FRACTION_antitarget[gain==0&loss==0]) %>% 
        dplyr::filter(FRACTION_antitarget==max(FRACTION_antitarget)) %>% 
        dplyr::select(TYPE,loss,gain,BASE,FRACTION_antitarget)%>% 
        dplyr::rename_with(~sub("_.*","",.x),starts_with("FRACTION")) %>%
        dplyr::ungroup()%>%
        dplyr::mutate(PGA=sum(FRACTION),method="bin",source="antitarget")
    base_all=tumour %>% dplyr::group_by(TYPE)%>% 
        dplyr::mutate(BASE=FRACTION_genome[gain==0&loss==0]) %>% 
        dplyr::filter(FRACTION_genome==max(FRACTION_genome))%>% 
        dplyr::select(TYPE,loss,gain,BASE,FRACTION_genome)%>% 
        dplyr::ungroup()%>%
        dplyr::rename_with(~sub("_.*","",.x),starts_with("FRACTION")) %>%
        dplyr::mutate(PGA=sum(FRACTION),method="base",source="all")
    base_target=tumour %>% dplyr::group_by(TYPE)%>% 
        dplyr::mutate(BASE=FRACTION_genome_target[gain==0&loss==0]) %>% 
        dplyr::filter(FRACTION_genome_target==max(FRACTION_genome_target))%>% 
        dplyr::select(TYPE,loss,gain,BASE,FRACTION_genome_target)%>% 
        dplyr::ungroup()%>%
        dplyr::rename_with(~sub("_.*","",.x),starts_with("FRACTION")) %>%
        dplyr::mutate(PGA=sum(FRACTION),method="base",source="target")
    base_antitarget=tumour %>% 
        dplyr::group_by(TYPE)%>% 
        dplyr::mutate(BASE=FRACTION_genome_antitarget[gain==0&loss==0]) %>% 
        dplyr::filter(FRACTION_genome_antitarget==max(FRACTION_genome_antitarget))%>% 
        dplyr::select(TYPE,loss,gain,BASE,FRACTION_genome_antitarget)%>% 
        dplyr::ungroup()%>%
        dplyr::rename_with(~sub("_.*","",.x),starts_with("FRACTION")) %>% 
        dplyr::mutate(PGA=sum(FRACTION),method="base",source="antitarget")

    dat_summ=dplyr::bind_rows(bin_all,bin_target,bin_antitarget,base_all,base_target,base_antitarget)
    return(dat_summ)

}