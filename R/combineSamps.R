


combineSamps <- function(species, minObs, region, sample_n) { 
  # NJBI this function refers to several global variables, e.g. tn - not good practice
  
  # set up defaults
  out_dat <- NULL
  out_meta <- NULL
  raw_occ <- NULL
  nRec_glob <- NA
  nRec_reg <- NA
  nRec <- NA
  gaps <- NULL
  rot <- NULL
  REGION_IN_Q <- paste0("psi.fs.r_", region)
  
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
    
    # retrieve input data
    dat <- out_meta$model$data() 
    
    # non-temporally explicit observation dataframe
    dat <- data.frame(year = dat$Year, # year
                      rec = dat$y) # records
    
    if(scaleObs == "global") {
      
      # subset to temporal window
      dat_glob <- dat[dat$year >= (t0 - (out_meta$min_year - 1)) & dat$year <= (tn - (out_meta$min_year - 1)), ]
      
      # number of global observations within time window t0 - tn
      nRec_glob <- sum(dat_glob$rec)
      
    } else { # regional scale metadata
      
      # region vs region_aggs
      if(region %in% out_meta$regions) { # explicit region - region
        
        # sites within selected region
        region_site <- dat[[paste0("r_", region)]][dat$Site]
        
      } else { # aggregate region - region_aggs
        
        # this gives you the names of the regions that make up the region_agg
        region_aggs <- unlist(out_meta$region_aggs[region])
        
        # sites within selected aggregate region (i.e., sites across all nested regions)
        region_site <- rowSums(sapply(region_aggs, function(x) dat[[paste0("r_", x)]][dat$Site]))
        
      }
      
      dat_reg <- data.frame(year = dat$year, # year
                        rec = dat$rec, # records
                        region_site = region_site) # sites included within region
      
      # subset to temporal window
      dat_reg <- dat_reg[dat_reg$region_site == 1 & dat_reg$year >= (t0 - (out_meta$min_year - 1)) & dat_reg$year <= (tn - (out_meta$min_year - 1)), ]
      
      # number of observations within region within time window t0 - tn
      nRec_reg <- sum(dat_reg$rec)
      
    }
    
  }
  
  if(scaleObs == "global") 
    # global number of observations
    nRec <- nRec_glob
  else
    # regional number of observations
    nRec <- nRec_reg
  
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
      
      if(scaleObs == "global") 
        datm <- dat_glob # temporally explicit global scale metadata
      else
        datm <- dat_reg # temporally explicit regional scale metadata
      
      first <- min(datm$year[datm$rec == 1]) + (t0 - 1)
      last <- max(datm$year[datm$rec == 1]) + (t0 - 1)
      
      firstMod <- t0
      
      lastMod <- tn
      
      yrs <- sort(unique(datm$year[datm$rec == 1]), decreasing = FALSE)
      
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
      if(!is.null((attr(out_meta, "metadata")$analysis$spp_Metrics)))
        rot <- as.data.frame(attr(out_meta, "metadata")$analysis$spp_Metrics)
      else # if model doesn't have rule of thumb data
        rot <- data.frame(median = NA, P90 = NA, visits_median = NA, visits_P90 = NA, prop_list_one = NA, prop_repeats_grp = NA, prop_abs = NA)
      
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
