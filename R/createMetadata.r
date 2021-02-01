library(dplyr)
library(pbapply)
objst <- '/data-s3/occmods'
list.files(objst)
metadata <- lapply(list.files(objst)[1:3], function(taxa){
  cat('Creating metadata for',taxa,'\n')
  data_types <- list.files(file.path(objst, taxa))
  ds <- pblapply(data_types, FUN = function(data_type){
    datasets <- list.files(file.path(objst, taxa, data_type))
    years <- str_extract(string = datasets, pattern = '[0-9]+') %>% as.numeric()
    data.frame(taxa = taxa,
               data_type = data_type,
               data_location = file.path(objst, taxa, data_type, datasets),
               dataset_name = datasets,
               most_recent = (years==max(years)),
               stringsAsFactors = FALSE)
  }) %>% bind_rows()
}) %>% bind_rows()

saveRDS(metadata, 'metadata.rds')
