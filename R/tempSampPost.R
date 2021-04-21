#' \code{tempSampPost} - from TrendSummaries, edited
#' @importFrom parallel mclapply
#' @importFrom parallel detectCores
#' @importFrom stats aggregate
#' @importFrom stats sd
#' @importFrom utils read.csv
#' @importFrom utils str
#' @importFrom utils write.csv
#' @export

tempSampPost <- function(indata = "../data/model_runs/", 
                         keep, # species to keep
                         keep_iter, # species to keep with iteration numbers
                         output_path = "../data/sampled_posterior_1000/",
                         region,
                         sample_n = 999,
                         tolerance = 0, # number of iterations above or below sample_n to be acceptable
                         group_name = "",
                         combined_output = TRUE,
                         max_year_model = NULL, 
                         min_year_model = NULL,
                         write = FALSE,
                         minObs = NULL,
                         scaleObs = "global", # scale at which to evaluate the number of records
                         t0, 
                         tn,
                         parallel = TRUE,
                         n.cores = NULL,
                         filetype = "rdata"){
  
  if(parallel & is.null(n.cores)) n.cores <- parallel::detectCores() - 1
  
  REGION_IN_Q <- paste0("psi.fs.r_", region)
  
  ### set up species list we want to loop though ###
  
  # extract minimum iteration number for chained models
  if(!is.null(keep_iter)) {
    
    # function to find minimum iteration for JASMIN models - Tom August
    findMinIteration <- function(list_of_file_names){
      
      if(length(list_of_file_names) < 1) stop('Error: list_of_file_names is empty')
      if(!is.character(list_of_file_names)) stop('Error: list_of_file_names must be a character')
      
      # remove the last number and file extension
      # find '_' followed by a signal number and a '.' and remove
      # that and everything that follows
      list_of_file_names <- gsub('_[[:digit:]]{1}$', '', list_of_file_names)
      
      # Extract the iterations number
      iterations <- regmatches(list_of_file_names, regexpr('[[:digit:]]+$', list_of_file_names))
      
      # Get minimum
      return(min(as.numeric(iterations)))
      
    }
    
    min_iter <- findMinIteration(keep_iter)
    
  }
  
  samp_post <- NULL # create the stacked variable, will be used if combined_output is TRUE.
  
  # load_rdata function
  # loads an RData file, and assigns it to an object name
  load_rdata <- function(fileName) {
    load(fileName)
    get(ls()[ls() != "fileName"])
  }
  
  # loop through species
  
  combineSamps <- function(species, minObs) { 
    # NJBI this function refers to several global variables, e.g. tn - not good practice
    
    # set up defaults
    out_dat <- NULL
    out_meta <- NULL
    raw_occ <- NULL
    nRec_glob <- NA
    nRec_reg <- NA
    nRec <- NA
    rot <- NULL
    
    if(!is.null(keep_iter)) {
      
      # chained models
      try(out_dat <- load_rdata(paste0(indata, species, "_20000_1.rdata")))  # where the first part of the model is stored for JASMIN models
      try(out_meta <- load_rdata(paste0(indata, species, "_", min_iter, "_1.rdata"))) # where metadata is stored for JASMIN models

      } else {
      
        # non-chained models
        
        if(filetype == "rds") {
      
          try(out_dat <- readRDS(paste0(indata, species, ".rds")))
          out_meta <- out_dat
          
        }
    
        else if(filetype == "rdata") {
          
          try(out_dat <- load_rdata(paste0(indata, species, ".rdata")))
          out_meta <- out_dat
          
        }
      
      }
    
    if(!is.null(out_dat$model) & !is.null(out_meta)) { # there is a model object to read from with metadata
      
      # global scale evaluation
      # HOW COULD WE MAKE THIS TEMPORALLY EXPLICIT?
      nRec_glob <- out_meta$species_observations # total number of observations for species
      
      # regional scale metadata
      if(scaleObs != "global") {
        
        dat <- out_meta$model$data() # retrieve input data
        
        # region vs region_aggs
        if(region %in% out_meta$regions) {
          
          # explicit region - region
          
          # sites within selected region
          region_site <- dat[[paste0("r_", region)]][dat$Site]
          
          } else {
            
            # aggregate region - region_aggs
            
            # this gives you the names of the regions that make up the region_agg
            region_aggs <- unlist(out_meta$region_aggs[region])
            
            # sites within selected aggregate region (i.e., sites across all nested regions)
            region_site <- rowSums(sapply(region_aggs, function(x) dat[[paste0("r_", x)]][dat$Site]))
            
          }
        
        dat <- data.frame(year = dat$Year, # year
                          rec = dat$y, # records
                          region_site = region_site) 
          
        # subset to region and temporal window
        dat <- dat[dat$region_site == 1 & dat$year >= (t0 - (out_meta$min_year - 1)) & dat$year <= (tn - (out_meta$min_year - 1)), ]
          
        nRec_reg <- sum(dat$rec) # number of observations within region within time window t0 - tn
        
      }
    }
    
    if(scaleObs == "global") 
      # global number of observations
      nRec <- nRec_glob # CAUTION - not temporally explicit
    else
      # regional number of observations
      nRec <- nRec_reg # CAUTION - temporally explicit
    
    print(paste0("load: ", species, ", ", scaleObs, " records: ", nRec))
    
    if(nRec >= minObs & # there are enough observations globally or in region
       region %in% c(out_meta$regions, names(out_meta$region_aggs)) & # the region or region_agg is listed for the species
       !is.null(out_dat$model) # there is a model object to read from
       ) { # the conditions are met
      
      if(!is.null(keep_iter)) {
        
        # chained models
        out_dat1 <- NULL
        out_dat2 <- NULL
        out_dat3 <- NULL
        
        try(out_dat1 <- load_rdata(paste0(indata, species, "_20000_1.rdata"))) # where occupancy data is stored for JASMIN models 
        raw_occ1 <- data.frame(out_dat1$BUGSoutput$sims.list[REGION_IN_Q])
        try(out_dat2 <- load_rdata(paste0(indata, species, "_20000_2.rdata"))) # where occupancy data is stored for JASMIN models 
        raw_occ2 <- data.frame(out_dat2$BUGSoutput$sims.list[REGION_IN_Q])
        try(out_dat3 <- load_rdata(paste0(indata, species, "_20000_3.rdata"))) # where occupancy data is stored for JASMIN models 
        raw_occ3 <- data.frame(out_dat3$BUGSoutput$sims.list[REGION_IN_Q])
        
        if(!is.null(out_dat1) & !is.null(out_dat2) & !is.null(out_dat3)) # if all models loaded correctly
          raw_occ <- rbind(raw_occ1, raw_occ2, raw_occ3)
        
      } else {
        
        # non-chained models
        raw_occ <- data.frame(out_dat$BUGSoutput$sims.list[REGION_IN_Q])
      
      }
      
      if(!is.null(raw_occ)) {
        
        # check whether the number of sims is enough to sample 
        # first calculate the difference between n.sims and sample_n.
        # positive numbers indicate we have more than we need
        if(!is.null(keep_iter)) {
        
          # chained models - sims from three chains
          diff <- (out_dat$BUGSoutput$n.sims * 3) - sample_n
        
        } else {
        
          diff <- out_dat$BUGSoutput$n.sims - sample_n

        }
      
        if(diff > tolerance) {
        
          # we have more sims in the model than we want, so we need to sample them
          raw_occ <- raw_occ[sample(1:nrow(raw_occ), sample_n), ]
        
        } else 
        
          if(abs(diff) <= tolerance) {
            # The number of sims is very close to the target, so no need to sample
            print(paste("no sampling required: n.sims =", out_dat$BUGSoutput$n.sims))
          
          } else
            stop("Error: Not enough iterations stored. Choose a smaller value of sample_n")
      
        colnames(raw_occ) <- paste("year_", out_meta$min_year:out_meta$max_year, sep = "")

        raw_occ$iteration <- 1:sample_n
        raw_occ$species <- species
      
        if(combined_output != TRUE) {
          write.csv(raw_occ, file = paste(output_path, gsub(".rdata", "" ,i), "_sample_", sample_n, "_post_", REGION_IN_Q, ".csv", sep = ""), row.names = FALSE)
        } 
      
        out1 <- raw_occ
      
        dat <- out_meta$model$data()
        dat <- data.frame(year = dat$Year,
                          rec = dat$y)
      
        first <- min(dat$year[dat$rec == 1]) + (t0 - 1)
        last <- max(dat$year[dat$rec == 1]) + (t0 - 1)
      
        firstMod <- t0
      
        lastMod <- tn
      
        yrs <- sort(unique(dat$year[dat$rec == 1]), decreasing = FALSE)
      
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
      
        # metadata data frame
        out2 <- data.frame(species, nRec_glob, nRec_reg, first, last, gap, firstMod, lastMod)
        
        # add rules of thumb metrics
        rot <- as.data.frame(attr(out_meta, "metadata")$analysis$spp_Metrics)
        
        # if model doesn't have rule of thumb data
        if(is.null(rot)) rot <- data.frame(median = NA, P90 = NA, visits_median = NA, visits_P90 = NA, prop_list_one = NA, prop_repeats_grp = NA, prop_abs = NA)
        
        # EqualWt and HighSpec decision trees (see https://www.biorxiv.org/content/10.1101/813626v1.full)
        rot$EqualWt <- ifelse(rot$prop_abs >= 0.990, rot$P90 >= 3.1, rot$P90 >= 6.7)
        rot$HighSpec <- ifelse(rot$prop_abs >= 0.958, rot$P90 >= 9.5, rot$P90 >= 29)
        
        # join rules of thumb data
        out2 <- cbind(out2, rot)
        
        print(paste("Sampled:", species))
      
        return(list(out1, out2))
      
      } else {
      
        print(paste("Error loading model:", species)) # this can happen for JASMIN models if one of the models in the chain doesn't load properly
      
        return(NULL)
      
      }
      
    } else {
      
      # informative messages
      if(!is.na(nRec) & !nRec >= minObs) 
        print(paste("Dropped (lack of observations):", species)) 
      else if(!is.na(nRec) & !region %in% c(out_meta$regions, names(out_meta$region_aggs)))
        print(paste("Dropped (region or region_aggs not present for species):", species))
      else print(paste("Error loading model:", species))
      
      return(NULL)
    }
  }
  
  if(parallel) outputs <- parallel::mclapply(keep, mc.cores = n.cores,
                                             combineSamps, minObs = minObs)
  else outputs <- lapply(keep, 
                         combineSamps, minObs = minObs)
  
  
  if(parallel) samp_post <- parallel::mclapply(outputs, mc.cores = n.cores,
                                               function(x)  y <- x[[1]])
  else samp_post <- lapply(outputs, 
                           function(x)  y <- x[[1]])
  
  samp_post <- do.call("rbind", samp_post)
  
  if(parallel) meta <- parallel::mclapply(outputs, mc.cores = n.cores,
                                          function(x) y <- x[[2]])
  else meta <- lapply(outputs, 
                      function(x) y <- x[[2]])
  
  meta <- do.call("rbind", meta)
  
  meta <- data.frame(Species = meta$species,
                     n_obs_global = meta$nRec_glob,
                     n_obs_regional = meta$nRec_reg,
                     min_year_data = meta$first,
                     max_year_data = meta$last,
                     min_year_model = meta$firstMod,
                     max_year_model = meta$lastMod,
                     gap_start = 0,
                     gap_end = 0,
                     gap_middle = meta$gap)
  
  colnames(meta) <- paste0(colnames(meta), "_r_", gsub("psi.fs.r_", "", REGION_IN_Q))
  
  if (write == TRUE) {
    save(samp_post, file = paste(output_path, group_name, "_all_spp_sample_", sample_n, "_post_", REGION_IN_Q, ".rdata", sep = ""))
  }
  
  return(list(samp_post, meta))
}
