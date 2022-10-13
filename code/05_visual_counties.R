library(maptools)
library(spdep)
library(MASS)
library(msm)
library(tigris)
library(ggplot2)
library(sf)
library(tidyr)
library(purrr)
library(plyr)
library(sys)


current_path <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))
load("../data/merged_fb_census_data_280922.RData")

## extract CT shapefile ##
us_shp<-counties(state = NULL, cb = FALSE, resolution = "500k", year = 2010)


## merge with pop and rate data ##
us_shp<-merge(us_shp,df,by.x='GEOID10',by.y='GEOID',all.x=T)

## put in sf format to use nice mapping packages ##
us_shp<-st_as_sf(us_shp)

## remove CTs with zero population (these are the islands off the MA coast) ##
us_shp<-us_shp[which(us_shp$acs>0 | us_shp$pep>0 | us_shp$worldpop >0 |us_shp$fb >0),]



## get quantiles of expected counts for maps ##
ma_q<-c(0,10,100,1000,10000,100000,1000000,10000000)

us_shp$exp_factor<-cut(us_shp$acs,breaks=ma_q,labels=c('<1','[1,5)','[5,10)','[10,20)','[20,30)','[30,50)','50+'),right=F)

## race-stratified expected counts, MA ##
us_maps<-ggplot(us_shp, aes(fill = exp_factor)) +
  geom_sf(colour=NA) +
  scale_fill_brewer(palette = "YlOrRd")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_line(colour = 'transparent'),
        panel.grid.minor = element_line(colour = 'transparent'),
        legend.title=element_blank(),text=element_text(size=18))
us_maps
