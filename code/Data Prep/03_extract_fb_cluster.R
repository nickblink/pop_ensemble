####################################################
## THIS CODE PROCESSES RAW WORLDPOP DATA FROM URL ##
## NEEDS TO BE RUN ON CANNON IN PARALLEL          ##
####################################################

## at command line run the following for interactive R job
## module load gcc/7.1.0-fasrc01 R/3.3.3-fasrc01 udunits/2.2.26-fasrc01 gdal/2.3.0-fasrc01 proj/5.0.1-fasrc01 geos/3.6.2-fasrc01
## export R_LIBS_USER=$HOME/apps/R:$R_LIBS_USER
## srun -p test --pty --mem 10000 -t 0-02:00 /bin/bash
## R --quiet




library(raster)
library(rgdal)
library(sf)
library(sp)
library(raster)
library(ggplot2)
library(USAboundaries)
library(tigris)
library(maptools)
options(tigris_use_cache = FALSE)


setwd("/n/home02/liyr8/pop_ensemble/meta_tif")
list <- list.files()
xx1 <- rasterToPoints(raster(list[1]), spatial = T)
colnames(xx1@data) <- "fb"
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

counties$fb_estimate <- sum_overct[,1]
save(xx1, sum_overct, counties, file = 'xx1.RData')
