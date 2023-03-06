#' @export

plot_phased=function(
    normal=NULL,
    tumour=NULL,
    normal_cov=10,
    normal_af=0.6,
    save=TRUE,
    format="png"
){
    
    library(patchwork)
    out_file=paste0(get_file_name(normal),".",format)
    normal=read.table(normal,sep="\t",header=TRUE,stringsAsFactors=FALSE)
    tumours=lapply(tumour,read.table,sep="\t",header=TRUE,stringsAsFactors=FALSE)
    tumours=dplyr::bind_rows(tumours)
    tumours_wider=tumours %>% tidyr::pivot_wider(id_cols=c(chrom,pos),names_from=id,values_from=af)
    tumours_wider=dplyr::left_join(normal %>% dplyr::select(chrom,pos,ref,alt,gt,af,depth),tumours_wider)
    tumours_long=tumours_wider %>% tidyr::pivot_longer(cols=!chrom:depth)
    tumours_long$gt_col=ifelse(tumours_long$gt=="1|0","blue","yellow")
    tumours_long_filt=tumours_long %>% 
    dplyr::filter(depth>=normal_cov,af<=normal_af&(1-normal_af)>=af)%>% 
    dplyr::arrange(pos)

    tumours_wider_cov=tumours %>% tidyr::pivot_wider(id_cols=c(chrom,pos),names_from=id,values_from=depth)
    tumours_wider_cov=dplyr::left_join(normal %>% dplyr::select(chrom,pos,ref,alt,gt,af,depth),tumours_wider_cov)
    tumours_long_cov=tumours_wider_cov %>% tidyr::pivot_longer(cols=!chrom:depth)
    tumours_long_cov$gt_col=ifelse(tumours_long_cov$gt=="1|0","blue","yellow")
    tumours_long_cov_filt=tumours_long_cov %>% 
    dplyr::filter(depth>=normal_cov,af<=normal_af&(1-normal_af)>=af) %>% 
    dplyr::arrange(pos)

    p1<-ggplot(tumours_long_filt,aes(x=as.numeric(as.factor(pos)),y=value,col=gt_col))
    p1<-p1+geom_hline(aes(yintercept=0.5),linetype="longdash")+
    geom_hline(aes(yintercept=0.25),alpha=0.5,linetype="longdash")+
    geom_hline(aes(yintercept=0.75),alpha=0.5,linetype="longdash")+
    facet_grid(name~"")
    p1<-p1+geom_point(size=0.1)+geom_smooth(se=FALSE)+
    scale_colour_identity()+
    theme_bw()+scale_y_continuous(limits=c(0,1),expand=c(0,0))

    p2<-ggplot(tumours_long_cov_filt,aes(x=as.numeric(as.factor(pos)),y=log2(value/depth),col=gt_col))
    p2<-p2+geom_hline(aes(yintercept=0.5),linetype="longdash")+
    geom_hline(aes(yintercept=0.25),alpha=0.5,linetype="longdash")+
    geom_hline(aes(yintercept=0.75),alpha=0.5,linetype="longdash")+
    facet_grid(name~"")
    p2<-p2+geom_point(size=0.1)+geom_smooth(se=FALSE)+
    scale_colour_identity()+
    theme_bw()+
    facet_grid(name~"")


    if(save){
         ggsave(
             out_file,
            plot=p1|p2,height=length(tumour)*600,
            width=2*1200,units="px"
        )
    }
   
}   