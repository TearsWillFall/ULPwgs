get_coverage_tfbs=function(gpos=NULL,bam=NULL,region=1000,central=150,threads=1){

    gpos=read_gpos(gpos=gpos,threads=threads) %>% 
    dplyr::distinct() %>% 
    dplyr::filter(start>1000)
  
    mclapply_os(1:nrow(gpos),FUN=function(pos){
           this.gpos=gpos[pos,]
           tfbs=data.frame(chrom=this.gpos$chrom,pos=(this.gpos$pos-1000):(this.gpos$pos+1000))
           get_coverage(bam=bam,gpos=tfbs,threads=1)
        }
    )
   


}