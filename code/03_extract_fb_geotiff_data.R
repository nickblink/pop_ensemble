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
options(tigris_use_cache = FALSE)

current_path <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))
setwd('..')

link <- "/Users/liyanran/Desktop/Research/Rachel/pop_ensemble/code/population_usa_2019-07-01.vrt"
test <- raster(link)


#pop_csv <- read.csv("population_usa_2019-07-01_part_1_of_6.csv")


link <- "/Users/liyanran/Desktop/Research/Rachel/pop_ensemble/code/population_usa18_-90_2019-07-01.tif"
test <- raster(link)
fb_pts1 <- rasterToPoints(test, spatial = T)






