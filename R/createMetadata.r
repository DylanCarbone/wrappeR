createMetadata <- function(file_location){
  library(dplyr)
  library(pbapply)
  metadata <- lapply(list.files(file_location), function(taxa){
    cat('Creating metadata for',taxa,'\n')
    data_types <- list.files(file.path(file_location, taxa))
    ds <- pblapply(data_types, FUN = function(data_type){
      datasets <- list.files(file.path(file_location, taxa, data_type))
      years <- str_extract(string = datasets, pattern = '[0-9]+') %>% as.numeric()
      data.frame(taxa = taxa,
                 data_type = data_type,
                 data_location = file.path(file_location, taxa, data_type, datasets),
                 dataset_name = datasets,
                 most_recent = (years==max(years)),
                 stringsAsFactors = FALSE)
    }) %>% bind_rows()
  }) %>% bind_rows()
  return(metadata)
}
