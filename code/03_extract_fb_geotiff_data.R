####################################################
## THIS CODE PROCESSES RAW WORLDPOP DATA FROM URL ##
## NEEDS TO BE RUN ON CANNON IN PARALLEL          ##
####################################################

## at command line run the following for interactive R job
## module load gcc/7.1.0-fasrc01 R/3.3.3-fasrc01 udunits/2.2.26-fasrc01 gdal/2.3.0-fasrc01 proj/5.0.1-fasrc01 geos/3.6.2-fasrc01
## export R_LIBS_USER=$HOME/apps/R:$R_LIBS_USER
## srun -p test --pty --mem 10000 -t 0-02:00 /bin/bash
## R --quiet


library(rgdal)
library(sf)
library(sp)
library(raster)
library(ggplot2)
library(USAboundaries)
library(tigris)
library(maptools)
options(tigris_use_cache = FALSE)

current_path <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))
setwd('..')

link <- "/Users/liyanran/Desktop/Research/Rachel/pop_ensemble/code/population_usa_2019-07-01.vrt"
test <- raster(link)


# pop_csv1 <- read.csv("population_usa_2019-07-01_part_1_of_6.csv")
# pop_csv2 <- read.csv("population_usa_2019-07-01_part_2_of_6.csv")
# pop_csv3 <- read.csv("population_usa_2019-07-01_part_3_of_6.csv")
# pop_csv4 <- read.csv("population_usa_2019-07-01_part_4_of_6.csv")
# pop_csv5 <- read.csv("population_usa_2019-07-01_part_5_of_6.csv")
# pop_csv6 <- read.csv("population_usa_2019-07-01_part_6_of_6.csv")

summary(pop_csv)
summary(pop_csv2)
summary(pop_csv3)
summary(pop_csv4)
summary(pop_csv5)
summary(pop_csv6)

link <- "/Users/liyanran/Desktop/Research/Rachel/pop_ensemble/code/meta_tif/population_usa18_-90_2019-07-01.tif"
test <- raster(link)
fb_pts1 <- rasterToPoints(test, spatial = T)


setwd("/Users/liyanran/Desktop/Research/Rachel/pop_ensemble/code/meta_tif")
list <- list.files()
data <- data.frame()
# 
# allrasters <- Reduce(merge, lapply(list, raster))
# allrasters <- lapply(allrasters, rasterToPoints)

#to check the index numbers of all imported raster list elements
allrasters

#call single raster element
allrasters[[1]]
xx1 <- rasterToPoints(raster(list[1]), spatial = T)
xx2 <- rasterToPoints(raster(list[2]), spatial = T)
colnames(xx1@data) <- "fb"
colnames(xx2@data) <- "fb"
xx3 <- spRbind(xx1, xx2)


for(i in list[-1]){
  path <- raster(i)
  points <- rasterToPoints(path, spatial = T)
  colnames(points@data) <- "fb"
  xx1 <- spRbind(xx1, points)
}

# get states list and filter to only contiguous/continental US
states <- tigris::states() 
states <- states[!(states$NAME %in% c("Alaska", "American Samoa", "Commonwealth of the Northern Mariana Islands", "Guam", "Hawaii", "Puerto Rico", "United States Virgin Islands")),]

# get counties and filter to only contiguous US
counties <- tigris::counties() 
counties <- counties[counties$STATEFP %in% states$STATEFP,]
counties <- as(counties, 'Spatial')
counties <- spTransform(counties, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

sum_overct <- sp::over(counties, xx1, fn=sum)

counties$WP_1km_estimate <- sum_overct[,1]


