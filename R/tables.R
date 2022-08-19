tc_table<-function(plt_data,gene_ai_limit=0,
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

  plt_data <- plt_data %>% dplyr::group_by(sample) %>% 
  dplyr::filter((all_log2+log2(ploidy/2)) < gene_log2_loss & evidence>gene_ai_limit) %>%
  dplyr::mutate(mean_beta=mean(as.numeric(beta),na.rm=TRUE))
 return(plt_data)

}