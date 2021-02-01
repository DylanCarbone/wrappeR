sampSubset <- function(subset, inPath) {
  
  spp.list <- read.csv(paste0(inPath, "speciesInfo.csv"), stringsAsFactors = FALSE)
  
  prio.list <- read.csv(paste0(inPath, "PrioritySpeciesNames_v2.csv"), 
                        stringsAsFactors = FALSE,
                        fileEncoding = "latin1")
  
  poll.list <- read.csv(paste0(inPath, "pollinators.csv"), stringsAsFactors = FALSE)
   
  poll2Drop.list <- read.csv(paste0(inPath, "poll2Drop.csv"))
  
  if (subset == "priority") {
    
    x <- spp.list$Species[which(tolower(spp.list$Species) %in% tolower(prio.list$NBN_Name) |
                                  tolower(spp.list$Species) %in% tolower(prio.list$DesigName) |
                                  tolower(spp.list$Species) %in% tolower(prio.list$PREFERRED_NBN_NAME) |
                                  tolower(spp.list$concept) %in% tolower(prio.list$MatchName))]
    
    y <- spp.list$concept[which(tolower(spp.list$Species) %in% tolower(prio.list$NBN_Name) |
                                  tolower(spp.list$Species) %in% tolower(prio.list$DesigName) |
                                  tolower(spp.list$Species) %in% tolower(prio.list$PREFERRED_NBN_NAME) |
                                  tolower(spp.list$concept) %in% tolower(prio.list$MatchName))]
    
    spp <- c(x,y)
    
  } else if (subset == "pollinators") {
    
    x <- spp.list$Species[which(tolower(spp.list$Species) %in% tolower(poll.list$species))]
    
    y <- spp.list$concept[which(tolower(spp.list$Species) %in% tolower(poll.list$species))]
    
    z <- read.csv(paste0(inPath, "hovConc.csv"))
    
    z <- unique(z)
    
    z <- z[,1]
    
    spp <- c(x,y, as.character(z))
    
  } else if (subset == "poll2Drop") {
    
    spp <- poll2Drop.list$species
    
  }
  
  return(spp)
  
}

extractMeta <- function(inPath, group, outPath, write, region) {
  
  file <- list.files(inPath,
                     full.names=T,
                     pattern = group)
  
  if (length(file) > 0) {
    
    load(file)
    
    colnames(taxa_data) <- toupper(colnames(taxa_data))
    
    if (group == "Hoverflies") {
      
      taxa_data <- taxa_data[-which(taxa_data$YEAR < 1980), ]
      
    }
    
    getMeta <- function(spp) {
      
      nRec <- length(taxa_data$CONCEPT[taxa_data$CONCEPT == spp])
      
      first <- min(taxa_data$YEAR[taxa_data$CONCEPT == spp])
      
      last <- max(taxa_data$YEAR[taxa_data$CONCEPT == spp])
      
      firstMod <- min(taxa_data$YEAR)
      
      lastMod <- max(taxa_data$YEAR)
      
      yrs <- sort(unique(taxa_data$YEAR[taxa_data$CONCEPT == spp]), decreasing = F)
      
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
      
      return(data.frame(spp, nRec, first, last, gap, firstMod, lastMod))
      
    }
    
    spp <- unique(taxa_data$CONCEPT)
    
    names(spp) <- spp
    
    taxa_meta <- purrr::map_df(.x = spp, .f = getMeta)
    
    out <-   data.frame(Species = taxa_meta$spp,
                        n_obs = taxa_meta$nRec,
                        min_year_data = taxa_meta$first,
                        max_year_data = taxa_meta$last,
                        min_year_model = taxa_meta$firstMod,
                        max_year_model = taxa_meta$lastMod,
                        gap_start = 0,
                        gap_end = 0,
                        gap_middle = taxa_meta$gap)
    
    colnames(out) <- paste0(colnames(out), "_r_", region)
    
    if (write == TRUE) {
      
      save(out, file=paste0(outPath, group, ".rdata"))
      
    }
    
    return(out)
    
  }
  
}



sampArray <- function(dat, startYear, endYear) {
  
  #combined.df <- dat
  
  #combined.df <- combined.df[,-ncol(combined.df)]
  
  #combined.df$iteration <- as.numeric(combined.df$iteration)
  
  #arr <- simplify2array(by(combined.df, combined.df$iteration, as.matrix))
  #print(str(arr))

  arr <- reshape2::acast(melt(dat, id=c("iteration", "species")), 
               species ~ variable ~iteration,
               value.var = "value")
  
  start <- (startYear - startYear) + 1
  
  end <- (endYear - startYear) + 1
  
  arr <- arr[,start:end,]
  
  #dimnames(arr)[[1]] <- 1:length(dimnames(arr)[[1]])
  
  #dimnames(arr)[[2]] <- 1:length(dimnames(arr)[[2]])
  
  return(arr)
  
}
