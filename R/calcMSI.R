#' \code{calcMSI} - Calculate multispecies indicator.
#' 
#' @description This function takes the output of \code{applyFilters} and produces a multispecies indicator
#'              using the chosen method.
#'
#' @param dat String. Object returned by \code{applyFilters}.
#' 
#' @param method String. Which indicator method to use. One of "lambda" or "bma". 
#' 
#' @param write Logical. Whether or not to write the outputs to  file.
#'                  
#' @param outPath String. Where to store the outputs if write = TRUE. 
#'
#' @param plotLabel String. Title for the indicator plot.
#' 
#' @param bmaInd String. Use "prime" if all species data are available for all species in all years, and 
#'               "reg" otherwise. See Freeman et al. 2020 for details. 
#' 
#' @param ... String. Title for the indicator plot.Additional arguments to be passed to \code{BRCindicators::bma}
#' 	  
#' @return A list with elements indicator, short term assessment, long term assessment and plot.
#'         
#' @export
#' 

calcMSI <- function(dat, 
                    method, 
                    write, 
                    outPath,
                    plotLabel,
                    bmaInd = NULL,
                    ...) {
  
  if (!method %in% c("lambda","bma")) stop("Method must be one of lambda or bma")
  
  if (method == "lambda") {
    
    # lambda method
    
    colnames(dat)[1:(ncol(dat) -2)] <- gsub("year_", "", colnames(dat[1:(ncol(dat) -2)]))
    
    minYr <- as.numeric(colnames(dat)[1])
    
    maxYr <- as.numeric(colnames(dat)[ncol(dat) - 2])
    
    arr <- sampArray(dat = dat,
                     startYear = minYr,
                     endYear = maxYr)
    
    ind <- BRCindicators::lambda_indicator( 
      input = arr, 
      index = 100, 
      threshold_sd = Inf, 
      threshold_Rhat = Inf, 
      threshold_yrs = 10, 
      upperQuantile = 0.975, 
      lowerQuantile = 0.025 
    )
    
    summary <- ind$summary
    
    lt <- BRCindicators::trend_assessment(ind, 
                                          start_year = minYr,  
                                          end_year = maxYr) 
    
    st <- BRCindicators::trend_assessment(ind, 
                                          start_year = maxYr - 5,  
                                          end_year = maxYr) 
    
    final <- BRCindicators::trend_assessment(ind, 
                                             start_year = maxYr,  
                                             end_year = maxYr) 
    
  } else {
    
    # bma method
    
    # fudge factor for occupancy as log of 0 is undefined
    fudgeOcc <- function(x, fudgeFac = 0.0001) {
      x[x == 0] <- fudgeFac
      x[x == 1] <- 1 - fudgeFac
      return(x)
    }
    
    # year column names
    year_cols <- colnames(dat)[grep("year", colnames(dat))]
    
    # select necessary columns
    dat <- dat[, c(year_cols, "species")]
    
    # log transform occupancy with adjustment for 1 and 0
    dat[, year_cols] <- apply(dat[, year_cols], 2, function(x) {log(fudgeOcc(x))})
    
    means <- aggregate(. ~ species, data = dat, mean, na.action = na.pass)
    
    sds <- aggregate(. ~ species, data = dat, sd, na.action = na.pass)
    
    getSumStats <- function(stat) {
      
      out <- reshape2::melt(stat, 
                  id.vars = "species",
                  variable.name = "year",
                  value.name = "index")
      
      return(out) 
      
    }
    
    means <- getSumStats(means)
    
    se <- getSumStats(sds)[,3]
    
    inDat <- data.frame(means, 
                        se = se)
    
    # prevent Error in node tau.obs when se = 0
    inDat$se[inDat$se == 0] <- 0.0001
    
    inDat$year <- as.numeric(gsub("year_", "", inDat$year))
    
    ind <- BRCindicators::bma(data = inDat,
               seFromData = TRUE,
               m.scale = "loge", # default to loge scale, based on log transformation of occupancy above
               ...)
    
    if (bmaInd != "prime") {
      
      summary <- data.frame(indicator = ind$Index.M,
                            lower = ind$lowerCI.M,
                            upper = ind$upperCI.M,
                            Species_Number = length(unique(inDat$species)))
      
    } else {
      
      summary <- data.frame(indicator = ind$Index.Mprime,
                            lower = ind$lowerCI.Mprime,
                            upper = ind$upperCI.Mprime,
                            Species_Number = length(unique(inDat$species)))
      
    }
    
    lt <- BRCindicators::trend_assessment(ind, 
                                          method = method,
                                          start_year = min(ind$year),  
                                          end_year = max(ind$year)) 
    
    st <- BRCindicators::trend_assessment(ind, 
                                          method = method,
                                          start_year = (max(ind$year) - 5),  
                                          end_year = max(ind$year)) 
    
    final <- BRCindicators::trend_assessment(ind, 
                                             method = method,
                                             start_year = max(ind$year),  
                                             end_year = max(ind$year)) 
    
  }

  out <- list(summary, ind, st, lt, final)
  
  names(out) <- c("Summary", "MetaData", "st", "lt", "final")

  return(out)
  
}