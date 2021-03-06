# ULPwgs

ULPwgs is a convenient tool that integrates existing tools to help to perform analysis of Ultra-Low-Pass Whole Genome Sequencing data (ULP-WGS) in R.

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
* [ichorCNA](https://github.com/broadinstitute/ichorCNA): Estimating tumor fraction in cell-free DNA from Ultra-Low-Pass Whole Genome Sequencing

## 2. System Requirements


Tested on fresh Ubuntu and Arch Linux installations. Should be working on other Linux distros too as long as equivalent packages are provided. 

In order to be able to download and compile the source files of all the required tools the following programs are required:
* git
* make
* gcc
* ant
* cmake
* autoconf

These tools can and should be installed using the terminal with the following commands:

* **For Ubuntu:**

  ```
  sudo apt install git make gcc ant cmake autoconf
  ```

* **For Arch Linux:**

  ```
  sudo pacman -S git ant make cmake gcc autoconf
  ```
  
Additional dependencies may be needed to succesfully install `devtools` package in R:
  
* **For Ubuntu:**

  ```
  sudo apt install build-essential libcurl4-gnutls-dev libxml2-dev openssl
  ```
  
* **For Arch Linux:**

  ```
  sudo pacman -S build-essential libcurl-gnutls libxml2-dev openssl
  ```


## 3. Installation Instructions

In order to install `ULPwgs` package we will be using R `devtools`:

```
install.packages("devtools")
devtools::install_github("TearsWillFall/ULPwgs")
```
If `devtools` package installation fails check System Requirements section, as you may be missing a dependency.

Once the `ULPwgs` package is installed we can use the function `install_required_tools()` to download and set up all the tools required for the bioinformatic process. This will create a directory named `tools` in the current working directory with all the tools. **Note: All functions within this package call scripts from the `tools` directory, therefore if this directory is moved, deleted or the current working directory is changed, these function will fail.**

```
ULPwgs::install_required_tools()
```

Or, alternatively.

```
library("ULPwgs")
install_required_tools()
```
## 4. Usage
![Bioinformatic Workflow example of ULPwgs](https://github.com/TearsWillFall/ULPwgs/blob/master/Graph.png?raw=true)
