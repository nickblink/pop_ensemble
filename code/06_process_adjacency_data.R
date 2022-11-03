current_path <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))
setwd('../')

D <- data.table::fread('data/county_adjacency.txt', 
                       col.names = c('name','GEOID','nbr_name','nbr_GEOID'))

for(i in 1:nrow(D)){
  nn_original <- UNFINISHED. USED THE ALREADY PROCESSED GITHUB VERSION BUT NEED TO QC THAT
}