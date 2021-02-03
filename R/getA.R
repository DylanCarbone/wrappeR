#' \code{getA} - Extracts the "a" parameter from occupancy model outputs.
#' 
#' @description "a" is the occupancy on the logit scale
#' Currently this only works for Regions. For the whole domain some recoding and calculation would be required.
#' This code has been copied from tempSampPost(). There is potential redundancy that could be streamlined at a later date
#' The data extracted this way are what we need for the bma method
#'
#' @importFrom parallel mclapply
#' @importFrom parallel detectCores 
#' @importFrom stats aggregate
#' @importFrom stats sd
#' @importFrom utils read.csv
#' @importFrom utils str
#' @importFrom utils write.csv
#' 
#' @export


getA <- function(indata = "../data/model_runs/", 
                         keep,
                         #output_path = "../data/sampled_posterior_1000/",
                         REGION_IN_Q = "a", 
                         group_name = "",
                         combined_output = TRUE,
                         max_year_model = NULL, 
                         min_year_model = NULL,
                         write = FALSE,
                         minObs = NULL,
                         t0, 
                         tn,
                         parallel = TRUE,
                         n.cores = NULL){
  
 if(parallel & is.null(n.cores)) n.cores <- parallel::detectCores() - 1
  
  ### set up species list we want to loop though ###
  spp.list <- list.files(indata, 
                         pattern = ".rdata") # species for which we have models
  
  spp.list <- gsub(".rdata", "", spp.list)
  spp.list <- spp.list[tolower(spp.list) %in% tolower(keep)]
  
  # loop through species
  
  readModel <- function(species, minObs) { 
    out <- NULL
    #raw_occ <- NULL
    load(paste0(indata, species, ".rdata"))
    
    nRec <- out$species_observations

    if(nRec >= minObs) {
      # Get the a parameter 
      out1 <- data.frame(
                    year = out$min_year:out$max_year,
                    logit_occ = out$BUGSoutput$mean[REGION_IN_Q],
                    logit_occ_sd = out$BUGSoutput$sd[REGION_IN_Q]
                    )
      out1$species <- species
      
      # metadata: from here is unchanged from get
      dat <- out$model$data()
      dat <- data.frame(year = dat$Year,
                        rec = dat$y)
      
      first <- min(dat$year[dat$rec == 1]) + (t0 - 1)
      last <- max(dat$year[dat$rec == 1]) + (t0 - 1)
      
      firstMod <- t0
      lastMod <- tn
      
      yrs <- sort(unique(dat$year[dat$rec == 1]), decreasing = F)
      
      gaps <- NULL
      
      if (length(yrs) > 1) {
        
        for (i in (1:length(yrs) - 1)) {
          gaps <- c(gaps, yrs[i+1] - yrs[i])
        }
      }
      
      if (!is.null(gaps)) {
        
        gap <- max(gaps)
        
      } else {
        gap <- 1
      } 

      out2 <- data.frame(species, nRec, first, last, gap, firstMod, lastMod)
      return(list(out1, out2))
    } else return(NULL)
  }
  
  if(parallel) outputs <- parallel::mclapply(spp.list, mc.cores = n.cores,
                          readModel, minObs=minObs)
  else outputs <- lapply(spp.list, 
                           readModel, minObs=minObs)
  
  
  if(parallel) a_post <- parallel::mclapply(outputs, mc.cores = n.cores,
                                   function(x)  x[[1]])
  else a_post <- lapply(outputs, function(x)  x[[1]])
  
  a_post <- do.call("rbind", a_post)
  
  if(parallel) meta <- parallel::mclapply(outputs, mc.cores = n.cores,
                              function(x) x[[2]])
  else meta <- lapply(outputs, function(x) x[[2]])
  
  meta <- do.call("rbind", meta)
  
  meta <- data.frame(Species = meta$species,
                     n_obs = meta$nRec,
                     min_year_data = meta$first,
                     max_year_data = meta$last,
                     min_year_model = meta$firstMod,
                     max_year_model = meta$lastMod,
                     gap_start = 0,
                     gap_end = 0,
                     gap_middle = meta$gap)
  
  colnames(meta) <- paste0(colnames(meta), "_r_", gsub("a_", "", REGION_IN_Q))

  return(list(a_post, meta))

}
