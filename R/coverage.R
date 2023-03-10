
#' @export

fraction_genome_altered=function(
    cnr=NULL,
    gain=seq(0,0.3,by=0.05),
    loss=seq(-0.6,0,by=0.05),
    threads=1
){

    dat=read.table(cnr,sep="\t",header=TRUE)
    all_info=mclapply_os(loss,FUN=function(x){
        info=lapply(gain,FUN=function(y){
            dat_tmp=dat
            dat_tmp$TYPE=ifelse(dat_tmp$log2>=y,"GAIN",
                ifelse(dat_tmp$log2<=x,"LOSS","WT"
            ))
    
            dat_rslt=dat_tmp%>% 
                dplyr::group_by(TYPE) %>% 
                dplyr::summarise(
                    N=dplyr::n(),
                    N_target=length(log2[gene!="Antitarget"]),
                    N_antitarget=length(log2[gene=="Antitarget"]),
                    weight=median(weight),
                    weight_target=median(weight[gene!="Antitarget"]),
                    weight_antitarget=median(weight[gene=="Antitarget"])
                ) %>% 
                dplyr::ungroup() %>% 
                dplyr::mutate(
                    FRACTION=N/sum(N),
                    FRACTION_target=N_target/sum(N_target),
                    FRACTION_antitarget=N_antitarget/sum(N_antitarget)
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
    
    return(all_info %>% dplyr::arrange(gain,dplyr::desc(loss)))

}