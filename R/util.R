get_sample_name=function(file_path=""){
  sample_name=unlist(strsplit(basename(file_path),"\\."))[1]
  return(sample_name)
}


get_file_extension=function(file_path=""){
    ext = strsplit(basename(file_path), split="\\.")[[1]]
    ext = paste(ext[-1],collapse=".")
    return(ext)
}
