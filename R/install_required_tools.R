#' Installs all the required tools for the analysis.
#'
#' This function downloads and compiles the source files of all the
#' tools needed for the ULP-WGS analysis. However, it DOESN'T provide
#' all the libraries and dependencies necessary for their succesful installation,
#' therefore, they have to be installed beforehand as described in the README.md
#' or at https://github.com/TearsWillFall/ULPwgs.
#' @export


install_required_tools=function(){
  BTools::install_tools(whitelist=c("FastQC","bwa","skewer",
  "samtools","picard","hmmcopy_utils","ichorCNA"))
}
