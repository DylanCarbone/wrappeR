#' \code{tempStackFilter} - Version of StackFilter from TrendSummaries, edited
#' @export

tempStackFilter <- function (input = "memory", dat = NULL, indata = NULL, output_path, group_name, metadata, region = NULL, 
                             minObs = NULL, maxStartGap = NULL, maxEndGap = NULL, maxMiddleGap = NULL, 
                             keepSpecies = NULL, removeSpecies = NULL, ClipFirst = FALSE, 
                             ClipLast = FALSE, write = FALSE) 
{
  if (input == "file") {
    load(indata)
  } else {
    stacked_samps <- dat
  }
  
  spp_list <- stacked_samps$species
  if (any(!is.null(c(minObs, maxStartGap, maxEndGap, maxMiddleGap)))) {
    if (is.null(metadata)) {
      stop("Metadata is required for this set of filters but not supplied")
    }
    if (length(setdiff(spp_list, metadata$Species)) > 0) {
      stop("Not all species have metadata available. Please ensure all species are included in the metadata and that all names match")
    }
  }
  if (!is.null(keepSpecies)) {
    if (length(setdiff(keepSpecies, spp_list)) > 0) {
      warning("Not all keepSpecies are contained within the indata")
    }
    if (!is.null(removeSpecies)) {
      if (length(setdiff(keepSpecies, removeSpecies) < 
                 length(keepSpecies))) {
        stop("Some species are listed under both keepSpecies and removeSpecies. Please ensure these vectors do not overlap")
      }
    }
  }
  if (!is.null(keepSpecies)) {
    stacked_samps <- droplevels(subset(stacked_samps, stacked_samps$species %in% 
                                         keepSpecies))
  }
  if (!is.null(removeSpecies)) {
    stacked_samps <- droplevels(subset(stacked_samps, !(stacked_samps$species %in% 
                                                          removeSpecies)))
  }
  if (is.null(region)) {
    if (!is.null(minObs)) {
      goodspecies <- droplevels(subset(metadata, metadata$n_obs >= 
                                         minObs))
      stacked_samps <- droplevels(subset(stacked_samps, 
                                         stacked_samps$species %in% goodspecies$Species))
    }
    if (!is.null(maxStartGap)) {
      goodspecies <- droplevels(subset(metadata, metadata$gap_start <= 
                                         maxStartGap))
      stacked_samps <- droplevels(subset(stacked_samps, 
                                         stacked_samps$species %in% goodspecies$Species))
    }
    if (!is.null(maxEndGap)) {
      goodspecies <- droplevels(subset(metadata, metadata$gap_end <= 
                                         maxEndGap))
      stacked_samps <- droplevels(subset(stacked_samps, 
                                         stacked_samps$species %in% goodspecies$Species))
    }
    if (!is.null(maxMiddleGap)) {
      goodspecies <- droplevels(subset(metadata, metadata$gap_middle <= 
                                         maxMiddleGap))
      stacked_samps <- droplevels(subset(stacked_samps, 
                                         stacked_samps$species %in% goodspecies$Species))
    }
  }
  else {
    if (!is.null(minObs)) {
      goodspecies <- droplevels(subset(metadata, metadata[, 
                                                          paste0("n_obs_r_", region)] >= minObs))
      stacked_samps <- droplevels(subset(stacked_samps, 
                                         stacked_samps$species %in% goodspecies$Species))
    }
    if (!is.null(maxStartGap)) {
      goodspecies <- droplevels(subset(metadata, metadata[, 
                                                          paste0("gap_start_r_", region)] <= maxStartGap))
      stacked_samps <- droplevels(subset(stacked_samps, 
                                         stacked_samps$species %in% goodspecies$Species))
    }
    if (!is.null(maxEndGap)) {
      goodspecies <- droplevels(subset(metadata, metadata[, 
                                                          paste0("gap_end_r_", region)] <= maxEndGap))
      stacked_samps <- droplevels(subset(stacked_samps, 
                                         stacked_samps$species %in% goodspecies$Species))
    }
    if (!is.null(maxMiddleGap)) {
      goodspecies <- droplevels(subset(metadata, metadata[, 
                                                          paste0("gap_middle_r_", region)] <= maxMiddleGap))
      stacked_samps <- droplevels(subset(stacked_samps, 
                                         stacked_samps$species %in% goodspecies$Species))
    }
  }
  if (ClipFirst == TRUE) {
    if (is.null(region)) {
      for (i in spp_list) {
        min_year_data <- metadata[metadata$Species == 
                                    i, "min_year_data"]
        min_year_model <- metadata[metadata$Species == 
                                     i, "min_year_model"]
        if (min_year_data > min_year_model) {
          stacked_samps[stacked_samps$species == i, names(stacked_samps) %in% 
                          paste0("year_", min_year_model:(min_year_data - 
                                                            1))] <- NA
        }
      }
    }
    else {
      for (i in spp_list) { 
        # this counter-intuitively loops over every iteration, not every species
        # although the syntax is clunky and ugly, it does turn years prior to first observation to NA
        min_year_data <- metadata[metadata$Species == 
                                    i, paste0("min_year_data_r_", region)]
        min_year_model <- metadata[metadata$Species == 
                                     i, paste0("min_year_model_r_", region)]
        if (min_year_data > min_year_model) {
          stacked_samps[stacked_samps$species == i, names(stacked_samps) %in% 
                          paste0("year_", min_year_model:(min_year_data - 
                                                            1))] <- NA
        }
      }
    }
  }
  if (ClipLast == TRUE) {
    if (is.null(region)) {
      for (i in spp_list) {
        max_year_data <- metadata[metadata$Species == 
                                    i, "max_year_data"]
        max_year_model <- metadata[metadata$Species == 
                                     i, "max_year_model"]
        if (max_year_data < max_year_model) {
          stacked_samps[stacked_samps$species == i, names(stacked_samps) %in% 
                          paste0("year_", (max_year_data + 1):max_year_model)] <- NA
        }
      }
    }
    else {
      for (i in spp_list) {
        max_year_data <- metadata[metadata$Species == 
                                    i, paste0("max_year_data_r_", region)]
        max_year_model <- metadata[metadata$Species == 
                                     i, paste0("max_year_model_r_", region)]
        if (max_year_data < max_year_model) {
          stacked_samps[stacked_samps$species == i, names(stacked_samps) %in% 
                          paste0("year_", (max_year_data + 1):max_year_model)] <- NA
        }
      }
    }
  }
  
  if (write == TRUE) {
    save(stacked_samps, file = paste0(output_path, group_name, 
                                      "stacked_Filter", ".rdata")) 
  }
  
  return(stacked_samps)
}


