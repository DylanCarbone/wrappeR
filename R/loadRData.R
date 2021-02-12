# function to deal with different file formats
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
loadRfile <- function(fileName){
  ext_file <- tolower(tools::file_ext(fileName))
  if(ext_file == 'rds'){
    ret <- readRDS(fileName)
  } else if(ext_file == 'rdata'){
    ret <- loadRData(fileName)
  } else {
    ret <- NULL
  }
}

