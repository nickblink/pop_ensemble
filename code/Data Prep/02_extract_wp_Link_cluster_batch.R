# sbatch --array=47-100 -J WP_processing run_wp_extraction.sh

library(dplyr)
library(sf)
library(sp)
library(raster)
library(ggplot2)
library(tigris)
options(tigris_use_cache = FALSE)

setwd('/n/holyscratch01/nethery_lab/Lab/nlink')
# setwd('C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Population Estimation/')

# get input data
inputs <- commandArgs(trailingOnly = TRUE)
print(inputs)
print(class(inputs))
i <- as.integer(inputs[[1]])

data_file <- 'data/usa_ppp_2019.tif'
data_file_UN <- 'data/usa_ppp_2019_UNadj.tif'

#### Make counties list ####
# get states list and filter to only contiguous/continental US
states <- tigris::states() 
states <- states[!(states$NAME %in% c("Alaska", "American Samoa", "Commonwealth of the Northern Mariana Islands", "Guam", "Hawaii", "Puerto Rico", "United States Virgin Islands")),]

# get counties and filter to only contiguous US
counties <- tigris::counties(year = 2019) 
counties <- counties[counties$STATEFP %in% states$STATEFP,]
counties <- as(counties, 'Spatial')
# counties1 <- spTransform(counties, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

counties <- sp::spTransform(counties, CRS("+proj=longlat +datum=WGS84 +no_defs"))

if(i == 1){
  save(counties, file = 'data/wp_county_aggregate_2019/counties_2019.RData')
}

#### Sum spatial points in batches ####
wp <- raster(data_file)
wp_UN <- raster(data_file_UN)

num_x <- 10
num_y <- 10
# US boundaries
xmin<-(-135)
ymin<-20
xmax<-(-60)
ymax<-50

# make grid 
x_along <- seq(from = xmin, to = xmax, length.out = num_x+1)
y_along <- seq(from = ymin, to = ymax, length.out = num_y+1)

# for(i in 35:100){
#   t0 <- Sys.time()
x_ind <- (i-1) %% num_x + 1
y_ind <- floor(i/num_y-.0001) + 1

# create the cropbox
xmin_i <- x_along[x_ind]
xmax_i <- x_along[x_ind + 1]
ymin_i <- y_along[y_ind]
ymax_i <- y_along[y_ind + 1]
cropbox<-extent(xmin_i, xmax_i, ymin_i, ymax_i)

# crop the data
wp_crop <- crop(wp, cropbox)
wp_crop_UN <- crop(wp_UN, cropbox)

# convert raster to spatial points
print('converting to spatial')
wp_pts <- rasterToPoints(wp_crop, spatial = T)
print('converting to spatial-UN')
wp_pts_UN <- rasterToPoints(wp_crop_UN, spatial = T)

print('summing across counties')
sum_overct <- sp::over(counties, wp_pts, fn=sum)
sum_overct_UN <- sp::over(counties, wp_pts_UN, fn=sum)

# counties$WP_1km_estimate <- sum_overct[,1]
# counties$WP_1km_estimate_UN <- sum_overct_UN[,1]

save(wp_pts, wp_pts_UN,
  sum_overct, sum_overct_UN, 
  #counties, 
  file = sprintf('data/wp_county_aggregate_2019/WP_2019_%s.RData', i))
# }


#### Aggregate the results ####

if(FALSE){
  # initialize the lists
  sum_list <- list()
  sum_list_UN <- list()
  
  dd <- dir('data/wp_county_aggregate_2019/')
  ff <- sapply(dd, function(nn){
    #nn <- sprintf('WP_2019_%s.RData', i)
    return(stringr::str_match(nn, '2019_\\s*(.*?)\\s*.R')[2])
  })
  setdiff(1:100, as.integer(ff))
  # good
  
  # load all the results into the lists
  for(i in 1:100){
    print(i)
    load(sprintf('data/wp_county_aggregate_2019/WP_2019_%s.RData', i))
    sum_list[[i]] <- sum_overct
    sum_list_UN[[i]] <- sum_overct_UN
  }
  
  # sum across the lists
  sum_all <- dplyr::bind_cols(sum_list) %>%
    rowSums(na.rm = T)
  
  sum_all_UN <- dplyr::bind_cols(sum_list_UN) %>%
    rowSums(na.rm = T)
}

load('data/wp_county_aggregate_2019/counties_2019.RData')


counties$WP <- sum_all
counties$WP_UN <- sum_all_UN
# save(counties, file = 'data/wp_2019.RData')
