#' @importFrom parallel mclapply
#' @importFrom parallel detectCores
#' @importFrom stats aggregate
#' @importFrom stats sd
#' @importFrom utils read.csv
#' @importFrom utils str
#' @importFrom utils write.csv
#' @export

tempSampPost <- function(indata = "../data/model_runs/", 
                         keep,
                         output_path = "../data/sampled_posterior_1000/",
                         REGION_IN_Q = "psi.fs",
                         sample_n = 1000,
                         tolerance = 3, # number of iterations above or below sample_n to be acceptable
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
  
  samp_post <- NULL # create the stacked variable, will be used if combined_output is TRUE.
  
  # loop through species
  
  combineSamps <- function(species, minObs) { 
    # NJBI this function refers to several global variables, e.g. tn - not good practice
    #print(species)
    out <- NULL
    raw_occ <- NULL
    load(paste0(indata, species, ".rdata"))
    
    nRec <- out$species_observations
    print(paste(species, nRec))
    
    if(nRec >= minObs) {
      raw_occ <- data.frame(out$BUGSoutput$sims.list[REGION_IN_Q])
      
      # check whether the number of sims is enough to sample 
      # first calculate the difference between n.sims and sample_n.
      # positive numbers indicate we have more than we need
      diff <- out$BUGSoutput$n.sims - sample_n
      if(diff > tolerance){
        # we have more sims in the model than we want, so we need to sample them
        raw_occ <- raw_occ[sample(1:nrow(raw_occ), sample_n),]
      } else 
        if(abs(diff) <= tolerance){
          # The number of sims is very close to the target, so no need to sample
          print(paste0("no sampling required: n.sims=", n.sims))
        } else
          stop("Not enough iterations stored. Choose a smaller value of sample_n")
      
      colnames(raw_occ) <- paste("year_", out$min_year:out$max_year, sep = "")
      raw_occ$iteration <- 1:sample_n
      raw_occ$species <- species
      
      if(combined_output != TRUE) {
        write.csv(raw_occ, file = paste(output_path, gsub(".rdata", "" ,i), "_sample_", sample_n, "_post_", REGION_IN_Q, ".csv", sep = ""), row.names = FALSE)
      } 
      
      out1 <- raw_occ
      
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
                    combineSamps, minObs=minObs)
  else outputs <- lapply(spp.list, 
                           combineSamps, minObs=minObs)
  
  
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
                     n_obs = meta$nRec,
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
