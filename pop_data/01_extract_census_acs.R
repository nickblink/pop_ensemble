################################################################################
## extract ACS and census denominator data for ages <65 only using tidycensus ##
## then format and merge together for analyses                                ##
################################################################################

library(tidycensus)
library(tigris)
library(sp)
library(stringr)
library(haven)
library(reshape2)
library(rstudioapi)

current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path))

#################
## 1. ACS DATA ##
#################

## variable names ##
v12 <- load_variables(2012, "acs1", cache = TRUE)

ma_acs<- get_acs(geography = "county",
                 variables = 'B01003_001', # Total Population
                 year=2010)
ma_acs<-as.data.frame(ma_acs)
ma_acs<-ma_acs[order(ma_acs$GEOID),]


####################
## 2. CENSUS DATA ##
####################

## census 2010 data ##
## variable names ##
v10 <- load_variables(2010, "sf1")

##extract using tidycensus ##
ma_ce<-get_decennial(geography = "county",
                     variables='P001001',
                     year=2010,sumfile='sf1')

## organize ##
ma_ce<-as.data.frame(ma_ce)
ma_ce<-ma_ce[order(ma_ce$GEOID),]

###########################
## 3. PEP estimated DATA ##
###########################

ma_pep <- 
  get_estimates(geography = "county", 
                product = "population", 
                year = 2010)
ma_pep<-as.data.frame(ma_pep)


##################################
## 4. merge ACS and census data ##
##################################

## read each processed dataset ##
# load('acs_data2206.RData')
# 
# load('ce_data2206.RData')
# 
# ## merge them by fips code, race, age, and sex ##
adat<-merge(ma_acs, ma_ce,by=c('GEOID','NAME','agecat','sex'))
# 
# save(adat,file='merged_denom_cov_data2206.RData')



# ## merge in the ice for race from the census (P003002-P003003)/P001001 ##
# ice<-get_decennial(geography = "tract",
#                    variables=c('P001001','P003002','P003003'),
#                    state = "MA",year=2010,sumfile='sf1',output = 'wide')
# ice<-as.data.frame(ice)
# ice<-data.frame('GEOID'=ice$GEOID,'ce_ice_racewb'=(ice$P003002-ice$P003003)/ice$P001001)
# 
# ma_ce<-merge(ma_ce,ice,by='GEOID')
# 
# ma_ce<-ma_ce[order(ma_ce$GEOID),]
# 
# save(ma_ce,file='ce_data2206.RData')
# 
# rm(list=ls())



