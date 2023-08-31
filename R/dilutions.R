dilute_tumour=function(
    bin_samtools=build_default_tool_binary_list()$bin_samtools,
    tumour=NULL,
    normal=NULL,
    read_fraction=1,
    tc_start=NULL,
    pl_start=NULL,
    tc_end=NULL,
    pl_end=NULL,
    ...
){

    if(read_fraction<1){
        subsample_samtools(bam=tumour,read_fraction=read_fraction)
    }


}