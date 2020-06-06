# ULPwgs

ULPwgs is a convenient tool that integrates existing tools to help to perform analysis of ultra-low-pass Whole Genome Sequencing data (ULP-WGS) in R.

## 1. Description

ULPwgs simplifies the process of analyzing ultra-low-pass Whole Genome Sequencing data (ULP-WGS) by giving a simplified installation process for commonly used bioinformatics tools, as well as providing a series of functions that allow calling these tools in R coding language.

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


Tested on fresh Ubuntu and Arch Linux installations. Should be working on other Linux distros too. 

In order to be able to download and compile the source files of all the required tools the following programs are required:
* make
* gcc
* ant
* cmake
* autoconf

These tools can and should be installed using the terminal with the following commands:

* **For Ubuntu:**

  ```
  sudo apt install make gcc ant cmake autoconf
  ```

* **For Arch Linux:**

  ```
  sudo pacman -S ant make cmake gcc autoconf
  ```
  
Additional dependencies may need to be installed to succesfully install `devtools` package in R:
  
* **For Ubuntu:**

  ```
  sudo apt install build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev
  ```
  
* **For Arch Linux:**

  ```
  sudo pacman -S build-essential libcurl-gnutls libxml2-dev openssl
  ```


## 3. Instructions

In order to install `ULPwgs` package we will be using R `devtools`:

```
install.packages("devtools")
devtools::install_github("TearsWillFall/ULPwgs")
```
If `devtools` package installation fails check System Requirements section, as you may be missing a dependency.

Once `ULPwgs` package is installed we can use the function `install_required_tools()` to download and install all the tools required for the bioinformatic analysis. This will create a directory named `tools` in the current working directory.

```
ULPwgs::install_required_tools()
```


