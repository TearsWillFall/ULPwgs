install.packages("roxygen2")
devtools::install_github("TearsWillFall/ULPwgs")

getwd()
setwd("ULPwgs")
devtools::load_all()
devtools::document()


?install_required_tools
install_required_tools()

fastqc(fastqc_R1="/home/osvaldas/Workshop_low_pass/fastq/SRR11742820.fastq")
trimming(fastqc_R1="/home/osvaldas/Workshop_low_pass/fastq/SRR11742820.fastq",verbose=TRUE)
fastqc(fastqc_R1="SRR11742820-trimmed.fastq.gz")
alignment(fastqc_R1="SRR11742820-trimmed.fastq.gz",
ref_genome="/home/osvaldas/Workshop_low_pass/ref/GRCh37.p13.genome.fa")
sort_and_index(file="SRR11742820-trimmed.bam",verbose=TRUE)
remove_duplicates(file="SORTED.bam/SRR11742820-trimmed.SORTED.bam",verbose=TRUE,ref_genome="/home/osvaldas/Workshop_low_pass/ref/GRCh37.p13.genome.fa")
sort_and_index(file="RMDUP.SORTED.bam/SRR11742820-trimmed.RMDUP.SORTED.bam")
qc_metrics(bam="SORTED.RMDUP.SORTED.bam/SRR11742820-trimmed.SORTED.RMDUP.SORTED.bam",
ref_genome="/home/osvaldas/Workshop_low_pass/ref/GRCh37.p13.genome.fa")
read_counter(verbose=TRUE,bam="SORTED.RMDUP.SORTED.bam/SRR11742820-trimmed.SORTED.RMDUP.SORTED.bam")
ichorCNA(wig="wig/SRR11742820-trimmed.wig")
