# ULPwgs

ULPwgs is a convenient tool that integrates existing tools to help to perform analysis of ultra-low-pass Whole Genome Sequencing data (ULP-WGS) in R.

## 1. Description

ULPwgs simplifies the process of analyzing ultra-low-pass Whole Genome Sequencing data (ULP-WGS) by giving a simplified installation process for commonly used bioinformatics tools, as well as providing a series of functions that allow calling these tools in the coding language of R.

**Tools:**
* [FastQC](https://github.com/s-andrews/FastQC):  A quality control analysis tool for high throughput sequencing data 
* [skewer](https://github.com/relipmoc/skewer): Tool for processing next-generation sequencing (NGS) paired-end sequences
* [bwa](https://github.com/lh3/bwa): Burrow-Wheeler Aligner for short-read alignment
* [htslib](https://github.com/samtools/htslib): C library for high-throughput sequencing data formats 
* [samtools](https://github.com/samtools/samtools/):Tools (written in C using htslib) for manipulating next-generation sequencing data
* [picard](https://github.com/broadinstitute/picard): A set of command line tools (in Java) for manipulating high-throughput sequencing (HTS) data and formats such as SAM/BAM/CRAM and VCF
* [hmmcopy_utils](https://github.com/shahcompbio/hmmcopy_utils): Tools for extracting read counts and gc and mappability statistics in preparation for running HMMCopy
* [ichorCNA](https://github.com/broadinstitute/ichorCNA): Estimating tumor fraction in cell-free DNA from ultra-low-pass whole genome sequencing

## 2. System Requirements
Tested on Ubuntu and Arch Linux.

In order to be able to download and compile the source files of all the required tools the following programs are required:
* GNU make. 
* gcc
* ant
* cmake
* autoconf

## 3. Instructions

Before installing the package
install.packages("devtools")
If you are in are in a fresh Ubuntu it's very likely sudo apt install build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev

