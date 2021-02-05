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




plot_coverage_panel=function(on_target="",off_target="",col=c(4,4),height=6,width=12,verbose=FALSE,output_dir=""){
  sep="/"

  if(output_dir==""){
    sep=""
  }

  sample_name=get_sample_name(on_target)
  dat1=read.table(on_target, stringsAsFactors=FALSE)
  if(is.character(dat1[1,1])){
    dat1=dat1[-1,]
  }
  dat1$type="On_Target"
  dat1=dat1[,c(col[1],ncol(dat1))]
  names(dat1)=c("Coverage","Type")
  dat1$Coverage=as.numeric(dat1$Coverage)

  dat=dat1
  if (off_target!=""){
    dat2=read.table(off_target, stringsAsFactors=FALSE)
    if(is.character(dat2[1,1])){
      dat2=dat2[-1,]
    }
    dat2$type="Off_Target"
    dat2=dat2[,c(col[2],ncol(dat2))]
    names(dat2)=c("Coverage","Type")
    dat2$Coverage=as.numeric(dat2$Coverage)
    bind_dat=dplyr::bind_rows(dat1,dat2)
    dat=bind_dat
  }

  p=ggplot(dat,aes(x=Type,y=Coverage))+geom_violin(aes(fill=Type),alpha=0.5)+geom_boxplot(width=0.1) + stat_summary(fun=median, geom="text", show.legend = FALSE,
               vjust=-1.2,hjust=-1.4, aes( label=round(..y.., digits=1)))+theme_classic()

  out_file=paste0(output_dir,sep,paste0(sample_name,".Region_Coverage.png"))
  ggsave(out_file,width=width,height=height)
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


plot_cumulative_cov=function(on_target="",off_target="",col=list(c(2,5),c(2,5)),height=6,width=12,verbose=FALSE,output_dir=""){
  sep="/"

  if(output_dir==""){
    sep=""
  }

  sample_name=get_sample_name(on_target)
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

  p=ggplot(dat,aes(y=Fraction_targets_above_depth,x=Depth))+geom_hline(aes(yintercept=0.9),linetype="dotted",alpha=0.75)+geom_hline(aes(yintercept=0.5),linetype="dotted",alpha=0.75)+geom_line(aes(col=Type),size=2) +theme_classic()+facet_wrap(Type~.,scales="free")
  out_file=paste0(output_dir,sep,paste0(sample_name,".Cumulative_Region_Coverage.png"))
  ggsave(out_file,width=width,height=height)
}
