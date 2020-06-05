#' Installs all the required tools for the analysis.
#'
#' This function downloads and compiles the source files of all the
#' tools needed for the ULP-WGS analysis. However, it DOESN'T provide
#' all the libraries and dependencies necessary for their succesful installation,
#' therefore, they have to be installed beforehand as described in the README.md
#' or at https://github.com/TearsWillFall/ULPwgs.

install_required_tools=function(){
  dir.create("tools")

  setwd("./tools")

  urls=c("https://github.com/s-andrews/FastQC.git",
  "https://github.com/relipmoc/skewer.git","https://github.com/samtools/samtools.git",
  "https://github.com/lh3/bwa.git","https://github.com/samtools/htslib.git",
  "https://github.com/broadinstitute/picard.git","https://github.com/shahcompbio/hmmcopy_utils.git",
  "https://github.com/broadinstitute/ichorCNA.git")

  sapply(urls,function (x) system(paste("git clone",x)))
  setwd("./FastQC")
  system("ant")
  system("chmod 755 bin/fastqc")
  setwd("..")
  setwd("./skewer")
  system("make")
  setwd("..")
  setwd("./bwa")
  system("make")
  setwd("..")
  setwd("./htslib")
  system("autoheader")
  system("autoconf -Wno-syntax")
  system("./configure")
  system("make")
  setwd("..")
  setwd("./samtools")
  system("autoheader")
  system("autoconf -Wno-syntax")
  system("./configure")
  system("make")
  setwd("..")
  setwd("./picard")
  system("./gradlew shadowJar")
  setwd("..")
  setwd("./hmmcopy_utils/")
  system(paste("cmake",getwd()))
  system("make")
  setwd("..")
  setwd("..")
}
