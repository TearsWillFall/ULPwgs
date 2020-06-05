#' Get the name of a sample from input file name
#'
#' This function takes the absolute/relative path to a file and
#' returns its base name without the file extension suffix.
#'
#' @param file_path Path to the input file
#' @return A string with the name of the file
#' @export
get_sample_name=function(file_path=""){
  sample_name=unlist(strsplit(basename(file_path),"\\."))[1]
  return(sample_name)
}

getwd()
#' Get the extension of a file
#'
#' This function takes the absolute/relative path to a file and
#' returns the file extension suffix.
#'
#' @param file_path Path to the input file
#' @return A string with the extension of the file
#' @export
get_file_extension=function(file_path=""){
    ext = strsplit(basename(file_path), split="\\.")[[1]]
    ext = paste(ext[-1],collapse=".")
    return(ext)
}
