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
  
  ### set up species list we want to loop though ###
  
  # default
  iter <- NULL
  
  # extract minimum and maximum iteration number for chained models
  if(!is.null(keep_iter)) {
    
    # function to find minimum and maximum iterations for JASMIN models - Tom August
    findIteration <- function(list_of_file_names){
      
      if(length(list_of_file_names) < 1) stop('Error: list_of_file_names is empty')
      if(!is.character(list_of_file_names)) stop('Error: list_of_file_names must be a character')
      
      # remove the last number and file extension
      # find '_' followed by a signal number and a '.' and remove
      # that and everything that follows
      list_of_file_names <- gsub('_[[:digit:]]{1}$', '', list_of_file_names)
      
      # Extract the iterations number
      iterations <- regmatches(list_of_file_names, regexpr('[[:digit:]]+$', list_of_file_names))
      
      # Get minimum and maximum
      return(c(min(as.numeric(iterations)), max(as.numeric(iterations))))
      
    }
    
    iter <- findIteration(keep_iter)
    
  }
  
  samp_post <- NULL # create the stacked variable, will be used if combined_output is TRUE.
  
  # loop through species

  if(parallel) outputs <- parallel::mclapply(keep, mc.cores = n.cores,
                                             combineSamps, indata = indata, keep_iter = keep_iter, region = region, sample_n = sample_n, tolerance = tolerance, combined_output = combined_output, minObs = minObs, scaleObs = scaleObs, t0 = t0, tn = tn, filetype = filetype, iter = iter)
  else outputs <- lapply(keep, 
                         combineSamps, indata = indata, keep_iter = keep_iter, region = region, sample_n = sample_n, tolerance = tolerance, combined_output = combined_output, minObs = minObs, scaleObs = scaleObs, t0 = t0, tn = tn, filetype = filetype, iter = iter)
  
  
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
                     gap_middle = meta$gap,
                     rot_median = meta$median,
                     rot_P90 = meta$P90,
                     rot_visits_median = meta$visits_median,
                     rot_visits_P90 = meta$visits_P90,
                     rot_prop_list_one = meta$prop_list_one,
                     rot_prop_repeats_grp = meta$prop_repeats_grp,
                     rot_prop_abs = meta$prop_abs,
                     rot_EqualWt = meta$EqualWt,
                     rot_HighSpec = meta$HighSpec)
  
  colnames(meta) <- paste0(colnames(meta), "_r_", region)
  
  if (write == TRUE) {
    save(samp_post, file = paste(output_path, group_name, "_all_spp_sample_", sample_n, "_post_", region, ".rdata", sep = ""))
  }
  
  return(list(samp_post, meta))
}
