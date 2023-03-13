
#' @export

fraction_genome_altered=function(
    cnr=NULL,
    gain=seq(0,1,by=0.01),
    loss=seq(-1,0,by=0.01),
    threads=1
){

    dat=read.table(cnr,sep="\t",header=TRUE)
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
    
    return(all_info %>% dplyr::arrange(gain,dplyr::desc(loss)) %>% dplyr::filter(TYPE!="WT"))

}