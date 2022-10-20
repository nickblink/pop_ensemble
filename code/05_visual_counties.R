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

# df$GEOID = as.character(df$GEOID)
# 
# write.csv(df,"../data/merged_fb_census_data_280922.csv", row.names = FALSE)


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

us_shp$exp_acs<-cut(us_shp$acs,breaks=ma_q,labels=c('<10','[10,100)','[100,1000)','[1000,10000)','[10000,100000)','[100000,1000000)','100000+'),right=F)
us_shp$exp_pep<-cut(us_shp$pep,breaks=ma_q,labels=c('<10','[10,100)','[100,1000)','[1000,10000)','[10000,100000)','[100000,1000000)','100000+'),right=F)
us_shp$exp_worldpop<-cut(us_shp$worldpop,breaks=ma_q,labels=c('<10','[10,100)','[100,1000)','[1000,10000)','[10000,100000)','[100000,1000000)','100000+'),right=F)
us_shp$exp_fb<-cut(us_shp$fb,breaks=ma_q,labels=c('<10','[10,100)','[100,1000)','[1000,10000)','[10000,100000)','[100000,1000000)','100000+'),right=F)

## race-stratified expected counts, MA ##
us_acs<-ggplot(us_shp, aes(fill = exp_acs)) +
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
us_acs
pdf('../pic/us_acs.pdf',width=7,height=6)


## race-stratified expected counts, MA ##
us_pep<-ggplot(us_shp, aes(fill = exp_pep)) +
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
us_pep
pdf('../pic/us_pep.pdf',width=7,height=6)

## race-stratified expected counts, MA ##
us_worldpop<-ggplot(us_shp, aes(fill = exp_worldpop)) +
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
us_worldpop
pdf('../pic/us_worldpop.pdf',width=7,height=6)

## race-stratified expected counts, MA ##
us_fb<-ggplot(us_shp, aes(fill = exp_fb)) +
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
us_fb
pdf('../pic/us_fb.pdf',width=7,height=6)



