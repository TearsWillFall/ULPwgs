#' Generate coverage plots for panel on and/or off target regions
#'
#' This function takes a tab separated file/s containing coverage information of on-target and/or off-target
#' regions and generates a plot of the coverage distribution.
#'
#'
#' @param on_target [REQUIRED] Path to file with on-target region coverage.
#' @param off_target [OPTIONAL] Path to file with off-target region coverage.
#' @param col [OPTIONAL] Column/s where coverage information is located. Default 4.
#' @param height [OPTIONAL] Plot height in inches. Default 6.
#' @param width [OPTIONAL] Plot width in inches. Default 12.
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param output_dir [OPTIONAL] Output directory path.
#' @import ggplot2
#' @export




plot_coverage_panel=function(on_target="",off_target="",col=c(5,4),
height=6,width=12,verbose=FALSE,output_dir="."){
  
  out_file_dir=set_dir(dir=output_dir,name="plots")

  dat1=read.table(on_target, stringsAsFactors=FALSE)

  dat1$type="On_Target"
  dat1=dat1[,c(col[1],ncol(dat1))]
  names(dat1)=c("Coverage","Type")
  dat1$Coverage=as.numeric(dat1$Coverage)

  dat=dat1
  if (off_target!=""){
    dat2=read.table(off_target, stringsAsFactors=FALSE)

    dat2$type="Off_Target"
    dat2=dat2[,c(col[2],ncol(dat2))]
    names(dat2)=c("Coverage","Type")
    dat2$Coverage=as.numeric(dat2$Coverage)
    bind_dat=dplyr::bind_rows(dat1,dat2)
    dat=bind_dat
  }

  p=ggplot(dat,aes(x=Type,y=Coverage))+geom_violin(aes(fill=Type),alpha=0.5)+geom_boxplot(width=0.1) +
   stat_summary(fun=median, geom="text", show.legend = FALSE,
               vjust=-1.2,hjust=-1.4, aes( label=round(..y.., digits=1)))+theme_classic()

  ggsave(paste0(out_file_dir,"/",get_file_name(on_target),".Region_Coverage.png"),width=width,height=height)
}



#' Generate cumulative coverage plots for panel on and/or off target regions
#'
#' This function takes a tab separated file/s containing coverage information of on-target and/or off-target
#' regions and generates a plot of the coverage distribution.
#'
#'
#' @param on_target [REQUIRED] Path to file with on-target region coverage.
#' @param off_target [OPTIONAL] Path to file with off-target region coverage.
#' @param col [OPTIONAL] Column/s where coverage information is located. Default columns 2 and 5.
#' @param height [OPTIONAL] Plot height in inches. Default 6.
#' @param width [OPTIONAL] Plot width in inches. Default 12.
#' @param verbose [OPTIONAL] Enables progress messages. Default False.
#' @param output_dir [OPTIONAL] Output directory path.
#' @import ggplot2
#' @export


plot_cumulative_cov=function(on_target="",off_target="",col=list(c(2,5),c(2,5)),
height=6,width=12,verbose=FALSE,output_dir="."){
  
  out_file_dir=set_dir(dir=output_dir,name="plots")

  dat1=read.table(on_target, stringsAsFactors=FALSE)
  dat1$type="On_Target"
  dat1=dat1[,c(unlist(col[1]),ncol(dat1))]
  names(dat1)=c("Depth","Fraction","Type")
  dat1$Depth=as.numeric(dat1$Depth)
  dat1$Fraction_targets_above_depth=1-cumsum(dat1$Fraction)
  dat1=dat1[dat1$Depth<4000,]

  dat=dat1
  if (off_target!=""){
    dat2=read.table(off_target, stringsAsFactors=FALSE)
    dat2$type="Off_Target"
    dat2=dat2[,c(unlist(col[2]),ncol(dat2))]
    names(dat2)=c("Depth","Fraction","Type")
    dat2$Depth=as.numeric(dat2$Depth)
    dat2$Fraction_targets_above_depth=1-cumsum(dat2$Fraction)
    dat2=dat2[dat2$Depth<20,]
    bind_dat=dplyr::bind_rows(dat1,dat2)
    dat=bind_dat
  }

  p=ggplot(dat,aes(y=Fraction_targets_above_depth,x=Depth))+geom_hline(aes(yintercept=0.9),linetype="dotted",
  alpha=0.75)+geom_hline(aes(yintercept=0.5),
  linetype="dotted",alpha=0.75)+geom_line(aes(col=Type),size=2) +
  theme_classic()+facet_wrap(Type~.,scales="free")
  
  ggsave(paste0(out_file_dir,"/",get_file_name(on_target),
  ".Cumulative_Region_Coverage.png"),width=width,height=height)
}

#' @export


plot_log2_beta=function(
    plt_data,gene_tg=TRUE,gene_ctrl=FALSE,gene_other=FALSE,
    log2_limit=2, gene_lbl=1, gene_lbl_evi=0.2,
    gene_lbl_beta_low=0.2, gene_lbl_beta_high=1,
    gene_lbl_size=2){


        filt <- c()
        if (gene_tg) {
          filt <- c(filt, "target")
        }

        if (gene_ctrl) {
          filt <- c(filt, "control")
        }
        if (gene_other) {
          filt <- c(filt, "other")
        }
        plt_data <- plt_data[plt_data$gene_type %in% filt, ]
        p <- ggplot(plt_data) +
        geom_vline(aes(xintercept = 0), linetype = "longdash", alpha = 0.25) +
        geom_vline(aes(xintercept = -unique(log2(ploidy/2))),
        linetype = "longdash",
        col = "red", alpha = 0.5
        ) +
        geom_point(aes(all_log2, beta,
        fill = chr,
        shape = ifelse(evidence >= 0.2, 24, 21)
        ),
        col = "black"
        ) +
        scale_shape_identity() +
        theme_bw() +
        scale_y_continuous(name = "BETA", limits = c(0, 1)) +
        scale_x_continuous(name = "ALL_LOG2", limits = c(
        -log2_limit,
        log2_limit
        )) +
        theme(legend.position = "none")

    if (gene_lbl == 1) {
        p <- p + ggrepel::geom_label_repel(
        data = plt_data %>% dplyr::filter(
            evidence >= gene_lbl_evi, beta >= gene_lbl_beta_low,
            beta <= gene_lbl_beta_high
        ),
        aes(all_log2, beta, label = gene, fill = col), max.overlaps = Inf,size=gene_lbl_size
        ) + scale_fill_identity()
    } else if (gene_lbl == 2) {
        p <- p + ggrepel::geom_label_repel(
        data = plt_data %>% dplyr::filter(
            evidence >= gene_lbl_evi, beta >= gene_lbl_beta_low,
            beta <= gene_lbl_beta_high
        ),
        aes(all_log2, beta, label = snps, fill = col), max.overlaps = Inf,size=gene_lbl_size
        ) + scale_fill_identity()
    } else if (gene_lbl == 3) {
        p <- p + ggrepel::geom_label_repel(
        data = plt_data %>% dplyr::filter(
            evidence >= gene_lbl_evi, beta >= gene_lbl_beta_low,
            beta <= gene_lbl_beta_high
        ),
        aes(all_log2, beta, label = informative_snps, fill = col), max.overlaps = Inf,
        size=gene_lbll_size
        ) + scale_fill_identity()
    }
    print(p)

}

#' @export


plot_cn=function(plt_data, cn_limit=3, gene_tg=TRUE, gene_ctrl=FALSE, gene_other=FALSE,
gene_lbl=1, gene_lbl_evi=0.2, gene_lbl_beta_low=0.2, gene_lbl_beta_high=1, gene_lbl_size=2){
  plt_data <- plt_data %>%
    dplyr::filter(!is.na(cnA), !is.na(cnB)) %>%
    dplyr::mutate(
      cnA = ifelse(cnA < 0, 0, ifelse(cnA > cn_limit, cn_limit, cnA)),
      cnB = ifelse(cnB < 0, 0, ifelse(cnB > cn_limit, cn_limit, cnB))
    )

  filt <- c()
  if (gene_tg) {
    filt <- c(filt, "target")
  }
  if (gene_ctrl) {
    filt <- c(filt, "control")
  }
  if (gene_other) {
    filt <- c(filt, "other")
  }
  plt_data <- plt_data[plt_data$gene_type %in% filt, ]

  p <- ggplot() +
    geom_abline(aes(intercept = 0, slope = 1), alpha = 0.4) +
    geom_vline(aes(xintercept = 1:(cn_limit - 1)), linetype = "longdash", alpha = 0.4) +
    geom_hline(aes(yintercept = 1:(cn_limit - 1)), linetype = "longdash", alpha = 0.4) +
    geom_point(data = plt_data, aes(cnA, cnB,
      fill = chr,
      shape = ifelse(evidence >= 0.2, 24, 21)
    ), col = "black") +
    scale_y_continuous(limits = c(0, cn_limit), expand = c(0.01, 0.01)) +
    scale_x_continuous(limits = c(0, cn_limit), expand = c(0.01, 0.01)) +
    scale_shape_identity() +
    theme_bw() +
    theme(legend.position = "none")
 if (gene_lbl == 1) {
        p <- p + ggrepel::geom_label_repel(
        data = plt_data %>% dplyr::filter(
            evidence >= gene_lbl_evi, beta >= gene_lbl_beta_low,
            beta <= gene_lbl_beta_high
        ),
        aes(cnA, cnB, label = gene, fill = col), max.overlaps = Inf,size=gene_lbl_size
        ) + scale_fill_identity()
    } else if (gene_lbl == 2) {
        p <- p + ggrepel::geom_label_repel(
        data = plt_data %>% dpylr::filter(
            evidence >= gene_lbl_evi, beta >= gene_lbl_beta_low,
            beta <= gene_lbl_beta_high
        ),
        aes(cnA, cnB, label = snps, fill = col), max.overlaps = Inf,size=gene_lbl_size
        ) + scale_fill_identity()
    } else if (gene_lbl == 3) {
        p <- p + ggrepel::geom_label_repel(
        data = plt_data %>% dplyr::filter(
            evidence >= gene_lbl_evi, beta >= gene_lbl_beta_low,
            beta <= gene_lbl_beta_high
        ),
        aes(cnA, cnB, label = informative_snps, fill = col), max.overlaps = Inf,
        size=gene_lbll_size
        ) + scale_fill_identity()
    }
    print(p)
}


#' @import patchwork
#' @export

plot_ai=function(plt_data,ai_limit=0.2,min_log2_limit=-0.3,
max_log2_limit=0.5,gene_tg=TRUE,gene_ctrl=FALSE,gene_other=FALSE){


  filt <- c()
  if (gene_tg) {
    filt <- c(filt, "target")
  }
  if (gene_ctrl) {
    filt <- c(filt, "control")
  }
  if (gene_other) {
    filt <- c(filt, "other")
  }
  plt_data <- plt_data[plt_data$gene_type %in% filt, ]
  
  map_triangles=function(plt_data){
    newcoord_up <- make_triangles(plt_data$s_order,
      plt_data$g_order)

    newcoord_down <- make_triangles(plt_data$s_order,
    plt_data$g_order,point="down")

    newcoord_down <- newcoord_down %>% dplyr::select(xdown = x, ydown = y)

    repdata=purrr::map_df(1:nrow(plt_data), function(i) plt_data[rep(i, 3), ])
    newdata <- dplyr::bind_cols(repdata, newcoord_up, newcoord_down)
    return(newdata)
  }
  


  plt_data=plt_data %>% dplyr::mutate(cn_t=ifelse(chr=="X",round(2**all_log2),
  2*round(2**all_log2)),
  beta=ifelse(evidence>=ai_limit,beta,1),log2.c=ifelse(min_log2_limit>all_log2,
  min_log2_limit,ifelse(max_log2_limit<all_log2,max_log2_limit,all_log2))) %>%
  dplyr::mutate(log2.c=ifelse(log2_cutoff_pass,log2.c,0))


  to_plot=list()
 
  to_plot[["X"]]<-map_triangles(plt_data %>% dplyr::filter(chr=="X") %>% 
  dplyr::arrange(s_order,gene) %>% dplyr::mutate(g_order=as.numeric(as.factor(gene))))

  to_plot[["autosome"]]<- map_triangles(plt_data%>% dplyr::filter(chr!="X") %>% 
  dplyr::arrange(s_order,gene) %>%dplyr::mutate(g_order=as.numeric(as.factor(gene))))



  tc_pl_plot=function(plt_data){
    tc_data=plt_data %>% dplyr::group_by(s_name,s_order) %>% 
    dplyr::distinct(tc,ploidy)%>% dplyr::ungroup() %>%
    dplyr::mutate(nc=1-tc) %>% 
    tidyr::pivot_longer(cols = -c(s_name,s_order,ploidy)) %>% 
    dplyr::mutate(col=ifelse(name=="tc","red","grey"),total=1,n_samples=max(s_order))

    df.grobs <- tc_data %>%
        dplyr::group_by(s_name,s_order,ploidy, total,n_samples) %>%
        dplyr::do(subplots = ggplot(., aes(1, value, fill = col)) +
          geom_col(position = "fill", colour = "black") +
          coord_polar(theta = "y") +
          theme_void() +
          theme(legend.position = "none") +
          scale_fill_identity()) %>%
        dplyr::mutate(subgrobs = list(annotation_custom(ggplot2::ggplotGrob(subplots),
          x = s_order - total/2 , y = 2.5*n_samples - n_samples / 2,
          xmax = s_order + total/2 , ymax = 2.5*n_samples + n_samples / 2
        )))



    p <- df.grobs %>% {
        p <- ggplot(.) +
          scale_fill_gradient2(low = "blue", high = "red", mid = "grey", midpoint = 2)+
          geom_point(aes(x = s_order, y = 2.5*n_samples)) +
          theme_void() +
          ylab("TC/PL") +
          theme(
            plot.margin = unit(c(0, 0, 0, 0), "pt"),
            axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            axis.line.x = element_blank(), axis.title.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(), legend.position = "none"
          )
      p <- p+ geom_tile(aes(x = s_order, y = 2.5*n_samples,fill=ploidy), col = "black")
      p<- p + .$subgrobs
      p<- p + scale_x_continuous(expand = c(0, 0))+scale_y_continuous(expand=c(0,0))
      
      
    }
    print(p&coord_equal())
  
  }

main_plot=function(data,type="autosome",show_cn=TRUE){
    tmp=data[[type]] %>% dplyr::arrange(g_order)
    p=ggplot(data[[type]]) +
    geom_polygon(aes(x = x, y = y, fill =(1-beta)/tc, 
    group = interaction(s_order, gene)),col="black",size=0.1) +
    scale_fill_gradient(low = "grey", high = "#fbff00", limits = c(0, 1))+
    ggnewscale::new_scale_fill() +
    geom_polygon(aes(x=xdown, y = ydown, fill = log2.c,
    group = interaction(s_order, gene)),col="black",size=0.1)+
    scale_fill_gradient2(low="blue",mid="grey",high="red")+
    theme(
        plot.margin = unit(c(0, 0, 0, 0), "pt"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_text(), axis.ticks.y = element_blank(),
        axis.title.y = element_text(angle = 90, vjust = 0.5), legend.position = "none"
    )+ylab(ifelse(type=="autosome","GENES AUTOSOMES","GENES X"))+
    scale_y_continuous(breaks = seq_along(unique(tmp$gene)), 
                      labels = unique(tmp$gene),expand=c(0,0))+
    scale_x_continuous(expand=c(0,0))

    if(show_cn){
      p<-p+geom_text(data=data[[type]],aes(s_order,as.numeric(as.factor(gene)),
      label=cn_t),col="white",size=3)
    }
    print(p&coord_equal())
  }

  

  name_plot=function(plt_data){
   p <- ggplot(plt_data %>% dplyr::distinct(s_name,s_order),aes(x = s_order, y = 2.5*max(s_order)))+ geom_tile(fill=NA, col = NA)
   p<-p+theme_void()+scale_x_continuous(expand = c(0,0))
   p<- p+geom_text(aes(x = s_order, y = 1.25*max(s_order),label=s_name),angle=45, hjust=0)

  print(p&coord_equal(clip="off"))
  }
  
  plts=list()
  plts[["names"]]=name_plot(plt_data)
  plts[["tc_pl"]]=tc_pl_plot(plt_data)
  plts[["autosome"]]=main_plot(to_plot,show_cn=TRUE)
  plts[["X"]]=main_plot(to_plot,type="X",show_cn=TRUE)



  p=patchwork::wrap_plots(plts,ncol=1)

  print(p)
}



plot_stats=function(plt_data,gene_tg=TRUE, gene_ctrl=FALSE, gene_other=FALSE){

    filt <- c()
    if (gene_tg) {
      filt <- c(filt, "target")
    }
    if (gene_ctrl) {
      filt <- c(filt, "control")
    }
    if (gene_other) {
      filt <- c(filt, "other")
    }
    plt_data <- plt_data[plt_data$gene_type %in% filt, ]

    tc_pl_plot=function(plt_data){
    
    tc_data=plt_data %>% dplyr::group_by(s_name,s_order) %>% 
    dplyr::distinct(tc,ploidy)%>% dplyr::ungroup() %>%
    dplyr::mutate(nc=1-tc) %>% 
    tidyr::pivot_longer(cols = -c(s_name,s_order,ploidy)) %>% 
    dplyr::mutate(col=ifelse(name=="tc","red","grey"),total=1,n_samples=max(s_order))

    df.grobs <- tc_data %>%
        dplyr::group_by(s_name,s_order,ploidy, total,n_samples) %>%
        dplyr::do(subplots = ggplot(., aes(1, value, fill = col)) +
          geom_col(position = "fill", colour = "black") +
          coord_polar(theta = "y") +
          theme_void() +
          theme(legend.position = "none") +
          scale_fill_identity()) %>%
        dplyr::mutate(subgrobs = list(annotation_custom(ggplot2::ggplotGrob(subplots),
          y = s_order - total/2 , x = 1 - total / 2,
          ymax = s_order + total/2 , xmin= 1 + total / 2
        )))



    p <- df.grobs %>% {
        p <- ggplot(.) +
          scale_fill_gradient2(low = "blue", high = "red", mid = "grey", midpoint = 2)+
          geom_point(aes(y = s_order, x = 1)) +
          theme_void() +
          ylab("TC/PL") +
          theme(
            plot.margin = unit(c(0, 0, 0, 0), "pt"),
            axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            axis.line.x = element_blank(), axis.title.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(), legend.position = "none"
          )
      p <- p+ geom_tile(aes(y = s_order, x= 1,fill=ploidy), col = "black")
      p<- p + .$subgrobs
      p<- p + scale_x_continuous(expand = c(0, 0))+scale_y_continuous(expand=c(0,0))
      
      
    }
    print(p&coord_equal())
  
  }

  main_plot=function(plt_data){

    summ=plt_data %>%
    dplyr::group_by(s_name,cn.call.corr,order,alpha,col,s_order) %>% 
    dplyr::summarise(N=dplyr::n(),genes=paste0(gene,collapse=";")) %>% 
    dplyr::arrange(desc(s_order),desc(N)) %>%
    dplyr::group_by(s_name,s_order) %>% 
    dplyr::mutate(Freq=N/sum(N),col=col) %>% dplyr::ungroup()
  

    p1<-ggplot(summ)+geom_bar(position="stack",stat="identity",aes(y=reorder(paste0(s_name,"_",s_order),s_order),
        x=Freq,fill=col,alpha=alpha,group = order),col="black")+
        scale_fill_identity()+
        theme_classic()+scale_x_continuous(expand=c(0,0))+
        xlab("Aberration Frequency %")+ylab("Samples")+scale_alpha_identity()
}
  plts=list()
  plts[["main"]] <-main_plot(plt_data)
  plts[["tc_pl"]] <-tc_pl_plot(plt_data)
  patchwork::wrap_plots(plts,ncol=2)

}




plot_pathways=function(plt_data){


      main_plot <- function(plt_data){

        p <- ggplot(plt_data,stat="identity",position="stack",
        aes(x=Class,y=Freq))+geom_point(shape=21,aes(fill=class_col))+geom_polygon(fill=NA,aes(
        group=cn_class,col=class_col))+ggiraphExtra::coord_radar()+theme_bw()+scale_y_continuous(limits=c(0,1))+
        theme(axis.text.y =element_blank(),axis.ticks.y =element_blank())+
        scale_color_identity()+scale_fill_identity()
        
        print(p)
      }
      p <- main_plot(plt_data)

      print(p)

  }
    
plot_tc=function(plt_data,gene_ai_limit=0,
gene_log2_loss=-0.05,gene_tg=TRUE,gene_ctrl=FALSE,
gene_other=FALSE){


  filt <- c()
  if (gene_tg) {
    filt <- c(filt, "target")
  }
  if (gene_ctrl) {
    filt <- c(filt, "control")
  }
  if (gene_other) {
    filt <- c(filt, "other")
  }
  plt_data <- plt_data[plt_data$gene_type %in% filt, ]

  plt_data <- plt_data %>% dplyr::group_by(sample) %>% dplyr::summarise(beta=mean(
  as.numeric(beta[(all_log2+log2(ploidy/2)) < gene_log2_loss & evidence>gene_ai_limit]),na.rm=TRUE),
  genes=paste0(gene[(all_log2+log2(ploidy/2)) < gene_log2_loss & evidence>gene_ai_limit ],collapse=";")) %>% 
  dplyr::mutate(tc=1-round(beta/(2-beta), digits = 2)) %>% dplyr::mutate(nc=1-tc)

  main_plot=function(plt_data){
    tc_data=plt_data[1,] %>%
    dplyr::distinct(tc,nc,genes)%>% dplyr::ungroup() %>%
    tidyr::pivot_longer(cols = -c(genes)) %>% 
    dplyr::mutate(col=ifelse(name=="tc","red","grey"),total=1)

    p<-ggplot(tc_data, aes(1, value, fill = col)) +
          geom_col(position = "fill", colour = "black",size=1.5) +
          coord_polar(theta = "y") +
          theme_void() +
          theme(legend.position = "none") +
          scale_fill_identity()+
          geom_text(position = position_stack(vjust = 0.5),
          aes(label=paste0(value*100,"%")),size=8)

    print(p)
  
  }

  p=main_plot(plt_data)
  p<-p+ggtitle("TC Manual")

  print(p)
}


      

