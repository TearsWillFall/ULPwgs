#' @export

plot_phased=function(
    normal=NULL,
    tumour=NULL,
    normal_cov=10,
    het_af=0.6,
    save=TRUE,
    format="png",
    method="single",
    plot_type="bar",
    threads=1
){
    
    library(patchwork)

    max_af_het=het_af
    min_af_het=1-het_af

    normals=lapply(normal,read.table,sep="\t",header=TRUE,stringsAsFactors=FALSE)
    tumours=lapply(tumour,read.table,sep="\t",header=TRUE,stringsAsFactors=FALSE)
    normals=dplyr::bind_rows(normals) %>% arrange(gtools::mixedorder(chromosome),pos)
    tumours=dplyr::bind_rows(tumours) %>% arrange(gtools::mixedorder(chromosome),pos)
    tumours_wider=tumours %>% tidyr::pivot_wider(id_cols=c(chrom,pos),names_from=id,values_from=af)
    tumours_wider=dplyr::left_join(normals %>% dplyr::select(chrom,pos,ref,alt,gt,af,depth),tumours_wider)
    tumours_long=tumours_wider %>% tidyr::pivot_longer(cols=!chrom:depth)
    tumours_long$gt_col=ifelse(tumours_long$gt=="1|0","blue","yellow")
    tumours_long_filt=tumours_long %>% dplyr::filter(depth>=normal_cov) %>%
    dplyr::filter(af<=max_af_het&min_af_het<=af)%>% 
    dplyr::arrange(pos)

    tumours_wider_cov=tumours %>% tidyr::pivot_wider(id_cols=c(chrom,pos),names_from=id,values_from=depth)
    tumours_wider_cov=dplyr::left_join(normals %>% dplyr::select(chrom,pos,ref,alt,gt,af,depth),tumours_wider_cov)
    tumours_long_cov=tumours_wider_cov %>% tidyr::pivot_longer(cols=!chrom:depth)
    tumours_long_cov$gt_col=ifelse(tumours_long_cov$gt=="1|0","blue","yellow")
    tumours_long_cov_filt=tumours_long_cov %>% dplyr::filter(depth>=normal_cov) %>%
    dplyr::filter(af<=max_af_het&min_af_het<=af) %>% 
    dplyr::arrange(pos)



    if(method=="single"){
        mclapply_os(unique(tumours_long_filt$name),FUN=function(id){


            p1<-ggplot(tumours_long_filt %>% dplyr::filter(name==id),
            aes(x=as.numeric(as.factor(pos)),
            y=value,col=gt_col))
            p1<-p1+geom_hline(aes(yintercept=0.5),linetype="longdash")+
            geom_hline(aes(yintercept=0.25),alpha=0.5,linetype="longdash")+
            geom_hline(aes(yintercept=0.75),alpha=0.5,linetype="longdash")+
            facet_grid(""~reorder(chromosome,gtools::mixedorder(chromosome)),space="free",scale="free_x")
        

            p2<-ggplot(tumours_long_cov_filt %>% dplyr::filter(name==id),
            aes(x=as.numeric(as.factor(pos)),xend=as.numeric(as.factor(pos)),
            y=log2(value/depth),yend=0.5,col=gt_col))
            p2<-p2+geom_hline(aes(yintercept=0),linetype="longdash")+
            geom_hline(yintercept=c(1:3),alpha=0.5,linetype="longdash")+
            geom_hline(yintercept=-c(1:3),alpha=0.5,linetype="longdash")+
            facet_grid(""~reorder(chromosome,gtools::mixedorder(chromosome)),space="free",scale="free_x")

            if(plot_type=="point"){

                p1<-p1+geom_point(size=0.1)+geom_smooth(se=FALSE)+
                scale_colour_identity()+
                theme_bw()+scale_y_continuous(limits=c(0,1))+
                facet_grid(""~reorder(chromosome,gtools::mixedorder(chromosome)),space="free",scale="free_x")+
                   facet_grid(name~"")+ theme(
                    axis.title.x=element_blank(),
                    panel.spacing = unit(0, "lines"),
                    axis.text.x=element_blank(),
                    panel.border = element_rect(colour = "black", fill=NA, size=0.1),
                    axis.ticks.x=element_blank())

                p2<-p2+geom_point(size=0.1)+geom_smooth(se=FALSE)+
                scale_colour_identity()+
                theme_bw()+
                facet_grid(""~reorder(chromosome,gtools::mixedorder(chromosome)),space="free",scale="free_x")+
                   facet_grid(name~"")+ theme(
                    axis.title.x=element_blank(),
                    panel.spacing = unit(0, "lines"),
                    axis.text.x=element_blank(),
                    panel.border = element_rect(colour = "black", fill=NA, size=0.1),
                    axis.ticks.x=element_blank())

            }else if(plot_type=="segment"){
                p1<-p1+geom_segment(aes(xend=as.numeric(as.factor(pos)),yend=0.5))+geom_smooth(se=FALSE)+
                scale_colour_identity()+
                theme_bw()+
                facet_grid(""~reorder(chromosome,gtools::mixedorder(chromosome)),space="free",scale="free_x")+
                   facet_grid(name~"")+ theme(
                    axis.title.x=element_blank(),
                    panel.spacing = unit(0, "lines"),
                    axis.text.x=element_blank(),
                    panel.border = element_rect(colour = "black", fill=NA, size=0.1),
                    axis.ticks.x=element_blank())

                p2<-p2+geom_segment(aes(xend=as.numeric(as.factor(pos)),yend=0.5))+geom_smooth(se=FALSE)+
                scale_colour_identity()+
                theme_bw()+
                facet_grid(""~reorder(chromosome,gtools::mixedorder(chromosome)),space="free",scale="free_x")+
                   facet_grid(name~"")+ theme(
                    axis.title.x=element_blank(),
                    panel.spacing = unit(0, "lines"),
                    axis.text.x=element_blank(),
                    panel.border = element_rect(colour = "black", fill=NA, size=0.1),
                    axis.ticks.x=element_blank())
            }
          
            out_file=paste0(get_file_name(id),".",format)

            ggsave(
                out_file,
                plot=p1|p2,height=600,
                width=2*1200,units="px"
            )

        },mc.cores=threads)
           
    }

    if(method=="panel"){


        p1<-ggplot(tumours_long_filt,
        aes(x=as.numeric(as.factor(pos)),y=value,col=gt_col))
        p1<-p1+geom_hline(aes(yintercept=0.5),linetype="longdash")+
        geom_hline(aes(yintercept=0.25),alpha=0.5,linetype="longdash")+
        geom_hline(aes(yintercept=0.75),alpha=0.5,linetype="longdash")+
        facet_grid(name~reorder(chromosome,gtools::mixedorder(chromosome)),space="free",scale="free_x")
    

        p2<-ggplot(tumours_long_cov_filt,
        aes(x=as.numeric(as.factor(pos)),y=log2(value/depth),col=gt_col))
        p2<-p2+geom_hline(aes(yintercept=0.5),linetype="longdash")+
        geom_hline(aes(yintercept=0.25),alpha=0.5,linetype="longdash")+
        geom_hline(aes(yintercept=0.75),alpha=0.5,linetype="longdash")+
        facet_grid(name~reorder(chromosome,gtools::mixedorder(chromosome)),space="free",scale="free_x")+
      
          if(plot_type=="point"){

                p1<-p1+geom_point(size=0.1)+geom_smooth(se=FALSE)+
                scale_colour_identity()+
                theme_bw()+scale_y_continuous(limits=c(0,1))+
                facet_grid(name~reorder(chromosome,gtools::mixedorder(chromosome)),space="free",scale="free_x")+ theme(
                    axis.title.x=element_blank(),
                    panel.spacing = unit(0, "lines"),
                    axis.text.x=element_blank(),
                    panel.border = element_rect(colour = "black", fill=NA, size=0.1),
                    axis.ticks.x=element_blank())

                p2<-p2+geom_point(size=0.1)+geom_smooth(se=FALSE)+
                scale_colour_identity()+
                theme_bw()+facet_grid(name~reorder(chromosome,gtools::mixedorder(chromosome)),space="free",scale="free_x")+theme_classic() + theme(
                    axis.title.x=element_blank(),
                    panel.spacing = unit(0, "lines"),
                    axis.text.x=element_blank(),
                    panel.border = element_rect(colour = "black", fill=NA, size=0.1),
                    axis.ticks.x=element_blank())


            }else if(plot_type=="segment"){

                p1<-p1+geom_segment(aes(xend=as.numeric(as.factor(pos)),yend=0.5))+geom_smooth(se=FALSE)+
                scale_colour_identity()+
                theme_bw()+
                facet_grid(name~reorder(chromosome,gtools::mixedorder(chromosome)),space="free",scale="free_x")+
                   facet_grid(name~"")+ theme(
                    axis.title.x=element_blank(),
                    panel.spacing = unit(0, "lines"),
                    axis.text.x=element_blank(),
                    panel.border = element_rect(colour = "black", fill=NA, size=0.1),
                    axis.ticks.x=element_blank())

                p2<-p2+geom_segment(aes(xend=as.numeric(as.factor(pos)),yend=0.5))+geom_smooth(se=FALSE)+
                scale_colour_identity()+
                theme_bw()+
                facet_grid(name~reorder(chromosome,gtools::mixedorder(chromosome)),space="free",scale="free_x")+
                   facet_grid(name~"")+ theme(
                    axis.title.x=element_blank(),
                    panel.spacing = unit(0, "lines"),
                    axis.text.x=element_blank(),
                    panel.border = element_rect(colour = "black", fill=NA, size=0.1),
                    axis.ticks.x=element_blank())

            }



        out_file=paste0(get_file_name(normal[1]),".",format)

       if(save){
         ggsave(
             out_file,
            plot=p1|p2,height=length(unique(tumours$id))*600,
            width=2*1200,units="px",limitsize = FALSE
        )
       }

    }




   
   
}  


#' @export

extract_het_pileup=function(
    tumour=NULL,
    normal=NULL,
    normal_cov=20,
    het_af=0.6,
    min_ai=0.1,
    threads=1,
    best_cov=1000
){

    max_af_het=het_af
    min_af_het=1-het_af

    normals=ULPwgs::mclapply_os(normal,function(x){read.table(x,sep="\t",header=TRUE,stringsAsFactors=FALSE)},mc.cores=threads)
    tumours=ULPwgs::mclapply_os(tumour,function(x){read.table(x,sep="\t",header=TRUE,stringsAsFactors=FALSE)},mc.cores=threads)
    normals=dplyr::bind_rows(normals)
    normals$gid=paste0(normals$chr,":",normals$pos,"_",normals$ref,"_",normals$alt)
    tumours=dplyr::bind_rows(tumours)
    tumours$gid=paste0(tumours$chr,":",tumours$pos,"_",tumours$ref,"_",tumours$alt)
    tumours$id=sub("TRAILS_TR..._","",sub("_BULK.*","",tumours$id))
    n_normals=length(unique(normals$id))
    tumours=dplyr::left_join(
        normals %>% 
        dplyr::filter(af>=min_af_het&af<=max_af_het,depth>=normal_cov) %>% 
        dplyr::group_by(gid) %>% 
        dplyr::mutate(N=n())%>%
        dplyr::filter(N==n_normals)%>%
        dplyr::ungroup()%>%
        dplyr::select(gid) %>%
        dplyr::distinct(),
        tumours,
        by=c("gid"),
        multiple = "all"
        )
    tumours=tumours %>%
        dplyr::group_by(id) %>%
        dplyr::mutate(
            depth_id=mean(depth)) %>%
        dplyr::ungroup()%>%
        dplyr::mutate(
            weight=depth_id/best_cov)%>%
        dplyr::group_by(id) %>%
        dplyr::mutate(
            BAF=0.5-af,
            aBAF=abs(0.5-af)
        )

    segmentation=ULPwgs::mclapply_os(unique(tumours$id),function(y){
        segmentation=lapply(unique(tumours$chrom),FUN=function(x){
            tumours_tmp=tumours %>% 
                dplyr::ungroup() %>% 
                dplyr::filter(chrom==x,id==y) %>% 
                dplyr::arrange(pos)
            seg=segment(CNA(tumours_tmp$aBAF, tumours_tmp$chrom, tumours_tmp$pos))
            tmp=cbind(seg$output,seg$segRows)
            tmp=tmp[-c(1:2)]
            merged=lapply(1:nrow(tmp),FUN=function(z){
                cbind(tumours_tmp[tmp[z,]$startRow:tmp[z,]$endRow,],tmp[z,])

        })
        
        return(merged)
    
    })
    segmentation=dplyr::bind_rows(segmentation)

    },mc.cores=threads)

    segmentation=dplyr::bind_rows(segmentation)
    

    
    segmentation=segmentation %>% 
        dplyr::rowwise()%>%
        dplyr::mutate(
            ai=ifelse(seg.mean>min_ai,1,0))%>%
        dplyr::group_by(gid) %>%
        dplyr::mutate(maxBAF=BAF[which.max(seg.mean)],
                    maxID=id[which.max(seg.mean)],
                    refAI=ai[which.max(seg.mean)]
        ) %>% dplyr::ungroup()

    ref=segmentation %>% dplyr::filter(id==maxID)
    ref$id="Reference"
    segmentation=rbind(segmentation,ref) 
  

    segmentation=segmentation %>% dplyr::rowwise()%>%
        dplyr::mutate(
            class=ifelse(maxBAF>0,"yellow",
                ifelse((maxBAF<0),"blue","grey")),
            obs_class=ifelse(ai==1&BAF>0,"yellow",
                ifelse(ai==1&BAF<0,"blue","grey"))
        )%>% dplyr::mutate(
            switch=ifelse(
                class=="yellow"&obs_class=="blue",1,
                ifelse(class=="blue"&obs_class=="yellow",1,0)),
            ref=ifelse(id==maxID,1,0)
        ) %>% dplyr::group_by(id) %>% 
        dplyr::mutate(
            n_snps_with_switch=sum(switch),
            f_switch_id=sum(switch)/n(),
            n_snps_with_ai=sum(ai),
            f_switch_id_ai=ifelse(sum(ai)==0,0,sum(switch)/sum(ai))
        
            )%>% 
        dplyr::group_by(gid) %>% 
        dplyr::mutate(
            n_samples_with_switch=sum(switch),
            f_switch_gid=sum(switch)/n(),
            n_samples_with_ai=sum(ai),
            f_switch_gid_ai=ifelse(sum(ai)==0,0,sum(switch)/sum(ai))
        )

    segmentation= segmentation  %>% 
    dplyr::ungroup() %>% 
    dplyr::arrange(id,chrom,pos) %>% 
    dplyr::mutate(gid=forcats::fct_inorder(gid))
}

