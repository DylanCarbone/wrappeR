#' \code{extractMetadata} - Extract metadata from master table created with createMetadata
#' 
#' @description This function extract more metadata from input data and model output files 
#'
#' @param meta table created by createMetadata
#' 
#' @return Dataframe of metadata
#' 
#' @import dplyr
#'         
#' @export


extractMetadata <- function(meta){
  library(dplyr)
  meta_complete <- NULL
  for(r in 1:dim(meta)[1]){
    meta_tmp <- meta[r,]
    if(meta_tmp$data_type == "occmod_outputs") {
      tryCatch({out <- loadRfile(list.files(meta_tmp$data_location, 
                                            full.names = TRUE)[1])
      meta_tmp$min_year <- out$min_year
      meta_tmp$max_year <- out$max_year
      
      meta_tmp$n_species <- length(unique(sub("_.*", "", 
                                              list.files(meta_tmp$data_location))))
      
      meta_tmp$regions <- paste(out$regions, collapse = "; ") 
      meta_tmp$sparta_v <- ifelse(!is.null(attr(out, "metadata")), 
                                  attr(out, "metadata")$analysis$session.info$otherPkgs$sparta$Version,
                                  NA)},
      error=function(e) NULL)
      
    } else {
      tryCatch({taxa_data <- loadRfile(list.files(meta_tmp$data_location, 
                                                  full.names = TRUE))
      meta_tmp$n_species <- ifelse(class(taxa_data) == "data.frame",
                                   length(unique(taxa_data$CONCEPT)),
                                   dim(taxa_data$spp_vis[-1])[2])},
      error=function(e) NULL)
    }
    # meta_complete <- rbind(meta_complete, meta_tmp)
    meta_complete <- bind_rows(meta_complete, meta_tmp)
  }
  return(meta_complete)
}


