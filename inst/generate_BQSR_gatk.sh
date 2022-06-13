#!/bin/bash -l

RScript -e "generate_BQSR_gatk(region=$1,bin_path=$2,bam=$3,ref_genome=$4,
snpdb=$5,output_dir=$6,verbose=$7)"