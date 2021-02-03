#' \code{applyClipping} - Apply filters to occupancy model outputs.
#' 
#' @description This function can be used to subset the occupancy model
#'              outputs based on number of records, priority/ pollinator status
#'              and region. It works at the level of the taxonomic group so
#'              must be applied across multiple groups if needed.
#'
#' @param roster String. A dataframe with columns: datPath, modPath, ver, indicator, region,
#'               nSamps, minObs, write, outPath, clipBy, group (see \code{createRoster}). 
#'
#' @param meta metadata produced by applySamp
#' @param clipBy String. Must be either Species or Group               
#' @param parallel Boolean. Should the operatoin run in parallel? If so then use n.cores-1.
#' 	  
#' @return A dataframe with processed model outputs to be passed to \code{calcMSI}.
#'         
#' @export


applyClipping <- function(data, parallel = TRUE, clipBy = "group") {
  # parallel argument currently not used
  
  data(speciesInfo)
  
  if (clipBy == "group" | data$clipBy != "species") { 
    # check that this arrangement makes sense
    # I think it's ok, given only two possibilities. Group overrides......
    data$meta[,3] <- min(data$meta[,3])
    data$meta[,4] <- max(data$meta[,4])
    
    }

  ## select which species to drop based on scheme advice etc. These are removed by stackFilter
  
  drop <- which(!is.na(speciesInfo$Reason_not_included) & speciesInfo$Reason_not_included != "Didn't meet criteria")
  
  drop <- c(as.character(speciesInfo$Species[drop]), 
            as.character(speciesInfo$concept[drop]))
  
  stacked_samps <- tempStackFilter(input = "memory",
                                   dat = data$samp_post,
                                   indata = NULL,
                                   output_path = NULL, 
                                   group_name = paste0(data$indicator, data$group), 
                                   metadata = data$meta, 
                                   region = data$region,
                                   minObs = data$minObs, 
                                   maxStartGap = 0, 
                                   maxEndGap = 0,
                                   maxMiddleGap = 10, 
                                   keepSpecies = NULL, 
                                   removeSpecies = drop,
                                   ClipFirst = TRUE, 
                                   ClipLast = TRUE)
  
  if (roster$write == TRUE) {
   save(stacked_samps, file = paste0(roster$outPath, roster$group, "_", roster$indicator,
                                     "_", roster$region, ".rdata"))
  }
  
  return(stacked_samps)
}
