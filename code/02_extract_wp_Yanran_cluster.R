####################################################
## THIS CODE PROCESSES RAW WORLDPOP DATA FROM URL ##
## NEEDS TO BE RUN ON CANNON IN PARALLEL          ##
####################################################

## at command line run the following for interactive R job
## module load gcc/7.1.0-fasrc01 R/3.3.3-fasrc01 udunits/2.2.26-fasrc01 gdal/2.3.0-fasrc01 proj/5.0.1-fasrc01 geos/3.6.2-fasrc01
## export R_LIBS_USER=$HOME/apps/R:$R_LIBS_USER
## srun -p test --pty --mem 10000 -t 0-02:00 /bin/bash
## R --quiet

## read command line arguments ##
#args<-commandArgs(TRUE)
#for (i in 1:length(args)) { eval (parse (text = args[[i]] )) }

#library(raster)
#library(rgdal)
library(sf)
library(sp)
library(raster)
library(ggplot2)
#library(USAboundaries)
library(tigris)
options(tigris_use_cache = FALSE)

current_path <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))
setwd('..')

## years of data:2010 ##

if (file.exists('us_wp.tif')){
}else{
  
  ## MA lat/long boundaries ##
  # xmin<-(-73.508142)
  # ymin<-41.237964
  # xmax<-(-69.928393)
  # ymax<-42.886589
  
  ## US lat/long boundaries ##
  xmin<-(-124.736342) # Westernmost -179.1479
  ymin<- 24.521208 #Southernmost 18.91042
  xmax<-(-66.945392) #Easternmost 179.7779
  ymax<-49.382808 #Northernmost 71.39042
  
  
  # xmin<-(-179.1479)
  # ymin<- 18.91042
  # xmax<-179.7779
  # ymax<-71.39042
  # 
  
  cropbox<-extent(xmin,xmax,ymin,ymax)
  
  ## read in raster ##
  #rname<- "https://data.worldpop.org/GIS/Population/Global_2000_2020/2010/USA/usa_ppp_2010.tif"
  #download.file(url=rname,destfile='/n/home02/liyr8/pop_ensemble/usa_ppp_2010.tif')
  
  wp<-raster('/n/home02/liyr8/pop_ensemble/usa_ppp_2010_1km_Aggregated.tif')
  
  ## view attributes ##
  wp
  
  ## crop to MA boundaries ##
  us_pop<-crop(wp,cropbox)
  
  ## export MA raster ##
  writeRaster(us_pop, 'us_wp.tif', overwrite=TRUE)
  
  ## delete big raster file ##
  #file.remove('/n/home02/liyr8/pop_ensemble/usa_ppp_2010.tif')
  
}


##### Aggregate the WP data into counties #####
wp <- raster('data/usa_ppp_2010_1km_Aggregated.tif')

wp_pts <- rasterToPoints(wp, spatial = T)

# get states list and filter to only contiguous/continental US
states <- tigris::states() 
states <- states[!(states$NAME %in% c("Alaska", "American Samoa", "Commonwealth of the Northern Mariana Islands", "Guam", "Hawaii", "Puerto Rico", "United States Virgin Islands")),]

# get counties and filter to only contiguous US
counties <- tigris::counties() 
counties <- counties[counties$STATEFP %in% states$STATEFP,]
counties <- as(counties, 'Spatial')
counties <- spTransform(counties, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

sum_overct <- sp::over(counties, wp_pts, fn=sum)

counties$WP_1km_estimate <- sum_overct[,1]

# save(wp_pts, sum_overct, counties, file = 'data/wp_1km_county_aggregate_2010.RData')
