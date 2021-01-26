#' speciesInfo
#'
#' Vaious information about species included in the UK occupancy indicators. Includes e.g. latin names and BRC concept codes
#' where available, and whether or not the species is on any of the four country's priority species lists.
#'
#' @docType data
#'
#' @usage data(speciesInfo)
#'
#' @format An object of class \code{"data frame"}. Columns for species latin name, taxonomic group, 
#'         source of latin name, whether or not it should be included in the final dataset, 
#'         reason for exclusion if applicable, additional details, BRC concept code (if available, otherwise gives latin name Charlie had),
#'         first year of records (will become outdated), last year with records (will become outdated), 
#'         whether or not the species is on the English priority list, Northern Irish priority list, 
#'         Scottish list and Welsh list.  
#'         .
#'
#' @keywords datasets
#'
#' @references JNCC 2020
#'
#' @source \href{https://jncc.gov.uk/our-work/conservation-designations-for-uk-taxa/}{QTL Archive}
#'
#' @examples
#' data(speciesInfo)
#' 
"speciesInfo"