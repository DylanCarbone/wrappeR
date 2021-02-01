#' \code{applySamp} - Reads a series of occupancy model and summarises the outputs

#' @param roster list
#' @param parallel Logical
#' @param sample Logical. Should the model sample from the posterior distribution or just get the a parameters instead?
#' @export
#' 

applySamp <- function(roster, parallel = TRUE, sample = TRUE) {
  
  data(speciesInfo)
  
  if (roster$indicator == "priority") {
    
    keepInds <- which(!is.na(speciesInfo[, roster$region])) 
    
    ## use both latin names and concept codes to screen for priority species 
    
    keep <- c(as.character(speciesInfo$Species[keepInds]), 
              as.character(speciesInfo$concept[keepInds]))
    
    keep <- keep[-which(is.na(keep))]
    
  } else if (roster$indicator == "pollinators") {
    
    keep <- sampSubset("pollinators",
                       inPath = roster$metaPath)
    
  } else {
    
    modFilePath <- file.path(roster$modPath, roster$group, "occmod_outputs", roster$ver)
    modFiles <- list.files(modFilePath)
    
    # read the suffix of the first model files
    filetype <- strsplit(modFiles[1], "\\.")[[1]][2] 
    if(!filetype %in% c("rdata", "rds")) stop("Model files must be either .rds or .rdata")
    
    # strip out the files
    modFiles <- modFiles[grepl(paste0(filetype, "$"), modFiles)] # dollar sign ensures the filetype suffix is at end of name
    
    # retain the species names
    keep <- gsub(pattern = paste0("\\.", filetype), repl = "", modFiles)
  }
  
  first_spp <- keep[[1]]
  
  if (substr(first_spp, (nchar(first_spp) + 1) - 2, nchar(first_spp)) %in% c("_1", "_2", "_3")) {
    
    keep <- gsub("(.*)_\\w+", "\\1", keep) # remove all after last underscore (e.g., chain "_1")
    keep <- gsub("(.*)_\\w+", "\\1", keep) # remove all after last underscore (e.g., iteration "_2000")
    
    keep <- unique(keep) # unique species names
    
  }
  
  ## select which species to drop based on scheme advice etc. These are removed by stackFilter
  
  drop <- which(!is.na(speciesInfo$Reason_not_included) & speciesInfo$Reason_not_included != "Didn't meet criteria")
  
  drop <- c(as.character(speciesInfo$Species[drop]), 
            as.character(speciesInfo$concept[drop]))
  
  if(sample == TRUE)
    out <- tempSampPost(indata = paste0(roster$modPath, roster$group, "/occmod_outputs/", roster$ver, "/"),
                        keep = keep,
                        output_path = NULL,
                        REGION_IN_Q = paste0("psi.fs.r_", roster$region),
                        sample_n = roster$nSamps,
                        group_name = roster$group,
                        combined_output = TRUE,
                        #max_year_model = 2018,
                        #min_year_model = 1970,
                        write = FALSE,
                        minObs = roster$minObs,
                        t0 = roster$t0,
                        tn = roster$tn,
                        parallel = parallel,
                        filetype = filetype)
  else
    out <- getA(indata = paste0(roster$modPath, roster$group, "/occmod_outputs/", roster$ver, "/"),
                        keep = keep,
                        REGION_IN_Q = paste0("a_", roster$region),
                        group_name = roster$group,
                        combined_output = TRUE,
                        write = FALSE,
                        minObs = roster$minObs,
                        t0 = roster$t0,
                        tn = roster$tn,
                        parallel = parallel)
  
  samp_post <- out[[1]]
  
  samp_post$species <- tolower(samp_post$species)
  
  meta <- out[[2]]
  
  meta[ ,1] <- tolower(meta[, 1])
  
  return(list(samp_post = samp_post, 
              meta = meta,
              indicator = roster$indicator,
              group_name = roster$group,
              region = roster$region,
              clipBy = roster$clipBy,
              minObs = roster$minObs))
}
