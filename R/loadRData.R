#' \code{loadRFile} - function to deal with different file formats
#' 
#' @description This function checks file format (rds or rdata) and loads it
#'
#' @param fileName name of the file to load
#' 
#' @return file loaded



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

