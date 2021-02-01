#' \code{createRoster} - Specify how to filter the occupancy model outputs.
#' 
#' @description This function can be used to specify the filters that will be applied to 
#'              the occupancy model outputs using \code{applyFilters}. It works by creating a list
#'              of 1-row dataframes with all the information needed for \code{applyFilters}. This list 
#'              is then applied to the applyFilters function later. Arguments should be provided
#'              as vectors of equal length, with each element in the vectors corresponding to one 
#'              call to \code{applyFilters}. 
#'              
#' @param index Numeric. Index of the number of taxonomic groups to 
#'              \code{applyFilters} across.
#'
#' @param modPath A character string or vector of strings. Location(s) of the 
#'                occupancy model outputs. 
#'
#' @param metaPath A character string or vector of strings. Location(s) of the 
#'                 occupancy model metadata.
#'   
#' @param ver A character string or vector of strings. Which set of occupancy 
#'            model outputs to use? Can be manually specified e.g. Charlie's 
#'            are "2017_Charlie"; or to source the most recent model versions 
#'            use "most_recent" (default), which uses model metadata stored at 
#'            "/data-s3/most_recent_meta" to identify the most recent model 
#'            version per taxonomic group.
#'       
#' @param group A character string or vector of strings. Taxonomic group(s), e.g. "Ants"
#' 
#' @param indicator A character string or vector of strings. Whether or not to 
#'                  subset species and, if so, based on what. Options are: 
#'                  "priority" for priority species; "pollinators" for 
#'                  pollinators; and all to return all species in the group.
#'                  
#' @param region A character string or vector of strings. One of "UK", "GB", 
#'               "ENGLAND", "WALES", "SCOTLAND", or "NORTHERN.IRELAND" per 
#'               taxonomic group.
#'
#' @param nSamps Numeric or numeric vector. Number of samples to extract from 
#'               each species' posterior distribution.
#'               
#' @param minObs Numeric or numeric vector. 
#'               Threshold number of observation below which a species is 
#'               dropped from the sample.
#'               
#' @param write Logical or logical vector. If TRUE then the outputs are 
#'              written as a .rdata file to outPath.
#'              
#' @param outPath A character string or vector of strings. Location to store 
#'                the outputs if write = TRUE. 
#'                
#' @param speciesToKeep A character vector of strings. the names of species to
#' include, this is used in combination with 'indicator'. ONLY species on both
#' lists will be included on the output.
#'  	  
#' @param clipBy A character string or vector of strings. One of "species" or 
#'               "group" indicating whether to clip outputs by the first and 
#'               last years of data for each species or for the whole group, 
#'               respectively.
#' 	  
#' @return A list of 1-row dataframes containing all arguments needed for \code{applyFilters}. \code{applyFilters}
#'         can then be applied to this list to filter models outputs for different taxonomic groups, which may 
#'         come from different rounds (e.g. Charlie's or later), for different regions, etc. 
#'         
#' @export
#' 

createRoster <- function(index,
                         modPath = "/data-s3/occmods/", 
                         metaPath,
                         ver = "most_recent",
                         group, 
                         indicator, 
                         region,
                         nSamps = 1000,
                         minObs = 50,
                         write,
                         outPath,
                         speciesToKeep = NA,
                         clipBy = "species",
                         t0,
                         tn) {
  
  if("most_recent" %in% ver) {
    
    # find metadata for most recent models
    mr_files <- list.files("/data-s3/most_recent_meta")
    mr_files_ver <- gsub("metadata_", "", mr_files)
    mr_files_ver <- as.numeric(gsub(".csv", "", mr_files_ver))
    
    # load metadata for most recent models
    mr <- read.csv(paste0("/data-s3/most_recent_meta/", mr_files[which.max(mr_files_ver)]),
                   stringsAsFactors = FALSE)
    
    # small data frame of ver and group
    tdf <- data.frame(ver = ver, group = group)
    
    # subset metadata to matching taxonomic groups and most recent models
    mr <- mr[mr$taxa %in% tdf$group & mr$most_recent == TRUE, ]
    
    # replace version with most recent model name
    ver <- ifelse(tdf$ver == "most_recent", mr$dataset_name, tdf$ver)
    
  }
 
  if (all(region %in% c("GB", "UK", "ENGLAND", "SCOTLAND", "WALES", "NORTHERN.IRELAND")) == FALSE) {
    
    warning("Warning: not all regions match either GB, UK, ENGLAND, SCOTLAND, WALES, or NORTHERN.IRELAND is this what you expected?")

  }
  
  df <- data.frame(index = index,
                   modPath = modPath,
                   datPath = paste0("/data-s3/occmods/",
                                    group, "/", ver, "/"),
                   metaPath = metaPath,
                   ver = ver, 
                   group = group, 
                   indicator = indicator,
                   region = region, 
                   nSamps = as.numeric(nSamps), 
                   minObs = minObs, 
                   write = write, 
                   outPath = outPath,
                   speciesToKeep = ifelse(test = is.na(speciesToKeep), 
                                          yes = NA, 
                                          no = as.character(paste(speciesToKeep, collapse = ','))),
                   clipBy = clipBy,
                   t0 = t0,
                   tn = tn)
  
  roster <- split(df, seq(nrow(df)))
  
}
