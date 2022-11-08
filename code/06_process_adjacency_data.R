current_path <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))
setwd('../')


D <- read.csv('data/countyadj.csv')

# read in saved wp data
wp <- read.csv('data/merged_wp_census_data_270922.csv')

# read in state abbrevs
abb <- read.csv('data/state_abbreviations.csv')

# match this ish up 
count <- 0
for(i in 1:nrow(D)){
  nn <- D$X[i]
  st <- substr(nn, nchar(nn) - 1, nchar(nn))
  ind <- which(abb$Code == st)
  
  D$X[i] <- gsub(substr(nn, nchar(nn) - 2, nchar(nn)), paste0(', ',abb$State[ind]), nn)
  
  # if(!(test %in% wp$NAME)){
  #   print(nn)
  # }
}

setdiff(wp$NAME, D$X)
setdiff(D$X, wp$NAME)
# Ok so there's some difference that I'm ignoring now for the sake of getting this running. But for the final work this should be looked into more

### Reformat to match the names in each
wp <- wp[wp$NAME %in% D$X,]

D2 <- D
rownames(D2) <- D2$X
D2 <- D2[,-1]
ind <- match(wp$NAME, rownames(D2))
D2 <- D2[ind, ind]

identical(rownames(D2), wp$NAME)

ind <- which(D2[,ncol(D2)] == 1)
rownames(D2)[ind]
# looks good to me from an online map

### Write the results
write.csv(wp, row.names = F, file = 'data/merged_wp_census_data2_081122.csv')

write.csv(D2, file = 'data/countyadj2.csv')

# Looking at the other adjacency matrix data from online (more annoying structure, but from the census, so perhaps better)
{
D2 <- data.table::fread('data/county_adjacency.txt', 
                       col.names = c('name','GEOID','nbr_name','nbr_GEOID'))

length(unique(D2$GEOID))
length(unique(D2$name))
# so we should merge by GEOID

D2[D2$name == 'Todd County, MN',]
}
