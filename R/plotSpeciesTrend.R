#' \code{plotSpeciesTrend} - Plot one species' trend, with 95% credible intervals
#' @param dat String. The samp_post element of the \code{applySamp} output
#' @param pspecies Name of the species to plot
#' @export
#' 

plotSpeciesTrend <- function(dat, species) {
  
  s1 <- dat[dat$species == species, ]
  
  s1Mean <- as.numeric(colMeans(s1[,1:(ncol(s1) - 2)]))
  
  quants <- sapply(1:(ncol(s1) - 2),
                   function(x) {
                     quantile(s1[,x], probs = c(0.025, 0.975))
                   })
  
  lower <- quants[1,]
  upper <- quants[2,]
  years <- gsub("year_", "", colnames(dat)[1:(ncol(dat) - 2)])
  
  plotDat <- data.frame(mean = s1Mean,
                        lower = lower,
                        upper = upper,
                        years = as.numeric(years))
  
  p <- ggplot2::ggplot(data = plotDat, aes(x = years, y = mean)) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(alpha = 0.5, aes(ymax = upper,
                                          ymin = lower)) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(species) +
    ggplot2::labs(x = "",
                  y = "occupancy")
  
  print(p)
  return(p)
  
}