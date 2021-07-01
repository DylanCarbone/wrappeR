#' \code{createMetadata} - Create master metadata table from objectstore data
#' repository
#' 
#' @description This function creates the metadata from a hierarchical data
#'  structure. It assumes a structure of:
#'  file_location/<taxa>/<data_type>/<dataset_name>. 
#'  It assumes that all dataset names contain the year of entry (in the
#'  format 2019 or 2019.2 for the second dataset in 2019). It also assumes
#'  that there is an equivalent folder at:
#'  file_location/<taxa>/input_data/<dataset_name>
#'  Which contains the input data.
#'  A lot of the parameters are hardcoded, and as such it is unlikely to return
#'  perfect metadata if the structure of outputs is altered.
#'  This takes about an hour to run as of 12 May 2021. You can pass it the
#'  previous metadata as a parameter, so it will not run metadata collection
#'  on any old datasets.
#'
#' @param file_location Path to root directory of data outputs on object store
#' @param oldMetadata Dataframe object previously created by this function. All
#'  rows of data in this object will not be run. Note that if data has changed
#'  in any of the datasets contained within the old metadata dataframe, these
#'  will not be re-run. Defaults to NULL (metadata will be calculated for all
#'  datasets in the object store)
#' 
#' @return Dataframe of metadata
#' 
#' @import pbapply
#' @import stringr
#' @import dplyr
#'         
#' @export

createMetadata <- function(file_location, oldMetadata = NULL){
  metadata <- suppressWarnings({pblapply(list.files(file_location), function(group){
    cat('Creating metadata for',group,'\n')
    data_types <- list.files(file.path(file_location, group))
    # We don't want the input data yet, we'll link that later
    data_types <- data_types[data_types != 'input_data']
    
    # Create a placeholder empty dataframe
    ds_placeholder <- data.frame(group = group,
                                 data_type = 'None',
                                 dataset_name = 'None',
                                 data_location = 'None',
                                 input_file = 'None',
                                 year = 'None',
                                 most_recent_year = 'None',
                                 most_recent = FALSE,
                                 min_year  = 'NA',
                                 max_year = 'NA',
                                 n_species = 'NA',
                                 n_species_input = 'NA',
                                 regions = 'NA',
                                 regions_aggs = 'NA',
                                 sparta_v = 'NA',
                                 provenance = 'NA',
                                 user = 'NA',
                                 date = 'NA',
                                 stringsAsFactors = FALSE)
    if(length(data_types)==0){
      # No output data types, just return placeholder
      ds <- ds_placeholder
    } else {
      # There's some data types. So get the metadata
      ds <- lapply(data_types, FUN = function(data_type){
        datasets <- list.files(file.path(file_location, group, data_type))
        if(length(datasets)==0){
          # We have no datasets in this data_type folder. Return appropriate dataframe
          df <- ds_placeholder
          df$data_type <- data_type
          return(df)
        } else {
          alreadyRun <- NULL
          if(!is.null(oldMetadata)){
            datasetsToRun <- !(file.path(file_location, group, data_type, datasets) %in%
                                 oldMetadata$data_location)
            alreadyRun <- oldMetadata[oldMetadata$data_location %in%
                                        file.path(file_location, group, data_type, datasets),]
            if(!any(datasetsToRun)){
              # All our datasets are in metadata already. Just return the previous metadata
              return(alreadyRun)
            } else {
              # We have new datasets. Subset the datasets to the new ones.
              datasets <- datasets[datasetsToRun]
            }
          }
          
          # Find the year of the datasets
          # Later datasets in the same year are listed as XXXX.1 or XXXX.2 etc
          strregex <- '[0-9]{4}([\\.0-9]+)?'
          datasets_with_years <- grepl(strregex, datasets)
          years <- lapply(str_extract(string = datasets,
                                      pattern = strregex), as.numeric) %>% unlist()
          # Find the maximum year, but first we have to drop NAs
          years_num <- years[!is.na(years)]
          if(length(years_num)>0){
            most_recent <- max(years_num)
          } else {
            most_recent <- 'Unknown'
          }
          # For all years where the regex failed, return "Unknown"
          years[is.na(years)] <- 'Unknown'
          
          df <- lapply(1:length(datasets), FUN = function(i){
            # Build data location
            data_location <- file.path(file_location, group, data_type, datasets[i])
            
            # Find input data
            input_path <- file.path(file_location, group, 'input_data', datasets[i])
            if(!dir.exists(input_path)){
              input_file <- 'Not Available'
              n_species_input <- 'Unknown'
            } else {
              input_file <- list.files(input_path, full.names = TRUE)
              if(length(input_file)>1){
                # There's more than 1. Perhaps one is the 'visitData' formatted file. Let's check.
                tmp_input <- input_file[!grepl('visitData',basename(input_file))]
                if(length(tmp_input)==1){
                  input_file <- tmp_input
                } else {
                  # We still have more than one. Maybe one is a source csv?
                  tmp_input <- input_file[!grepl('csv$',basename(input_file))]
                  if(length(tmp_input)==1){
                    input_file <- tmp_input
                  } else {
                    # Oh well, we still have two. Return 'Multiple Available'
                    input_file <- 'Multiple Available' 
                    n_species_input <- 'Unknown'
                  }
                }
              } else if(length(input_file)==0){
                input_file <- 'Not Available'
                n_species_input <- 'Unknown'
              }
            }
            
            if(!(input_file %in% c('Multiple Available', 'Not Available'))){
              speciesListInput <- suppressWarnings({loadRfile(input_file)})
              if('spp_vis' %in% names(speciesListInput)){
                # Visit Data format. Find the names
                speciesList <- as.character(names(speciesListInput$spp_vis)[-1])
              } else {
                # Probably standard format. Find names in first column (hopefully species)
                speciesList <- as.character(unique(data.frame(speciesListInput)[,1]))
              }
              n_species_input <- length(speciesList)
            } else {
              speciesListInput <- NULL
            }
            
            # Pull in some key metadata
            ourfiles <- list.files(data_location, full.names = TRUE)
            ourfiles <- ourfiles[!file.info(ourfiles)$isdir]
            if(length(ourfiles)==0){
              # There's no data. Return what we can.
              dq <- ds_placeholder
              dq$dataset_name <- datasets[i]
              dq$data_type <- data_type
              dq$data_location <- data_location
              dq$input_file <- input_file
              dq$n_species_input <- n_species_input
              dq$year <- as.character(years[i])
              dq$most_recent_year <- as.character(most_recent)
              dq$most_recent <- years[i]==most_recent
            } else {
              # We have data.
              # If the files are daisychained, we need the first item in a daisychain Check for that
              firstRun <-
                str_extract(ourfiles, pattern = '(?<=_)[0-9]+(?=_[1-3]{1}\\.[RrDdSsAaTt]+$)') %>%
                as.numeric()
              if(any(!is.na(firstRun))){
                # We have daisychains
                minDaisy <- min(firstRun[!is.na(firstRun)])
                refFile <- ourfiles[which(firstRun == minDaisy)[1]]
              } else {
                #No daisies. Just return the first file
                refFile <- ourfiles[1]
              }
              min_year <- max_year <- regions <- regions_aggs <- sparta_v <-
                provenance <- user <- submit_date <- n_species <- n_species_input <- 'Unknown'
              if(length(ourfiles)==0){
                # We have no outputs. Just set n_species to 0 so we know we tried to find data
                n_species <- 0
              } else {
                if(data_type == "occmod_outputs"){
                  tryCatch({
                    out <- suppressWarnings({loadRfile(refFile)})
                    min_year <- out$min_year
                    max_year <- out$max_year
                    
                    # This bit of code finds the input species, then checks each one off
                    # against the input species list.
                    # It's a bit slow, but only needs to be run once per metadata output,
                    # and is probably more trustworthy than the hacky regex option
                    if(!is.null(speciesListInput)){
                      n_species <- lapply(speciesList, FUN = function(species){
                        any(grepl(tolower(species),tolower(basename(ourfiles))))
                      }) %>% unlist() %>% sum()
                    } else {
                      # We don't have the relevant input data, so instead use a hacky bit of regex
                      n_species <- length(unique(sub("_[0-9]+_[1-3]{1}\\.[RrDdSsAaTt]+$",
                                                     "",
                                                     basename(ourfiles))))
                    }
                    
                    regions <- paste(out$regions, collapse = "; ")
                    regions_aggs <- paste(names(out$region_aggs), collapse = "; ")
                    if(!is.null(attr(out, "metadata"))){
                      metadata_analysis <- attr(out, "metadata")$analysis
                      if('session_info' %in% names(metadata_analysis)){
                        sparta_v <- metadata_analysis$session_info[[2]][['sparta']]
                      } else {
                        sparta_v <- metadata_analysis$session.info$otherPkgs$sparta$Version
                      }
                      params <- lapply(c('provenance','user','date'), FUN = function(parameter){
                        output_param <- metadata_analysis[[parameter]]
                        if(is.null(output_param)) output_param <- 'Unknown'
                        as.character(output_param)
                      })
                      provenance <- params[[1]]
                      user <- params[[2]]
                      submit_date <- params[[3]]
                    } else {
                      sparta_v <- provenance <- user <- submit_date <- 'Unknown'
                    }
                  }, error=function(e){
                    min_year <- max_year <- n_species <- n_species_input <- regions <- regions_aggs <-
                      sparta_v <- provenance <- user <- submit_date <- 'Unknown'
                  })
                } else {
                  tryCatch({
                    out <- suppressWarnings({loadRfile(refFile)})
                    n_species <- length(unique(sub("_[0-9]+_[1-3]{1}\\.[RrDdSsAaTt]+$",
                                                   "",
                                                   basename(ourfiles))))
                    if(class(out)=='data.frame'){
                      if('CONCEPT' %in% names(out)){
                        n_species <- length(unique(out$CONCEPT))
                      }
                    } else {
                      if('spp_vis' %in% names(out)){
                        n_species <- dim(out$spp_vis[-1])[2]
                      }
                    }
                    
                    if('min_year' %in% names(out)){
                      min_year <- out$min_year
                    }
                    if('max_year' %in% names(out)){
                      max_year <- out$max_year
                    }
                    if('regions' %in% names(out)){
                      regions <- paste(out$regions, collapse = "; ")
                    }
                    if('regions_aggs' %in% names(out)){
                      regions_aggs <- paste(names(out$regions_aggs), collapse = "; ")
                    }
                  }, error=function(e){
                    min_year <- max_year <- n_species <- n_species_input <- regions <- regions_aggs <-
                      sparta_v <- provenance <- user <- submit_date <- 'Unknown'
                  })
                }
              }
              dq <-
                data.frame(group = group,
                           data_type = data_type,
                           dataset_name = datasets[i],
                           data_location = data_location,
                           input_file = input_file,
                           year = as.character(years[i]),
                           most_recent_year = as.character(most_recent),
                           most_recent = (years[i]==most_recent),
                           min_year = as.character(min_year),
                           max_year = as.character(max_year),
                           n_species = as.character(n_species),
                           n_species_input = as.character(n_species_input),
                           regions = regions,
                           regions_aggs = regions_aggs,
                           sparta_v = sparta_v,
                           provenance = provenance,
                           user = user,
                           date = submit_date,
                           stringsAsFactors = FALSE)
            }
            # dq is the dataframe for a specific dataset within a data type within a group
            # i.e. will have a unique 'group + data_type + dataset_name' combination
            return(dq)
          }) %>% bind_rows()
          df <- bind_rows(alreadyRun, df)
          # df is the dataframe for all datasets with a data type within a group
          # i.e. will be the dataframe for a unique 'group + data_type' combination
          return(df)
        }
      }) %>% bind_rows()
    }
    # ds is the dataframe for all datasets within all data types within a group
    # i.e. will be the dataframe for a 'group'
    return(ds)
  })}) %>% bind_rows()
  # metadata combines all groups, their sub-'data_types' and sub-'datasets'
  return(metadata)
}