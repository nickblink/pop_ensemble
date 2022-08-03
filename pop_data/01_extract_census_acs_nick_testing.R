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
library(dplyr)

current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path))

#################
## 1. ACS DATA ##
#################

## variable names ##
v12 <- load_variables(2012, "acs1", cache = TRUE)

ma_acs<- get_acs(geography = "county",
                 variables = c('B01003_001'), # Total Population
                 year=2010)
ma_acs<-as.data.frame(ma_acs)
ma_acs<-ma_acs[order(ma_acs$GEOID),]

a <- ma_acs %>%
  dplyr::filter(variable == 'B01003_001')
b <- ma_acs %>%
  dplyr::filter(variable == 'B01001_001')

c <- merge(a, b, by = 'GEOID')
sum(c$estimate.x == c$estimate.y)

# so these two columns are the same. Good to know

####################
## Comparing ACS 1, 3, and 5 
####################

acs1 <- get_acs(geography = "county",
                 variables = c('B01003_001', 'B01001_001'), # Total Population
                 year=2010,
                survey = 'acs1')
# only for geographies with populations 65k or greater



## 3 year ACS estimates are no longer available, so nvm

acs1 <- get_acs(geography = "county",
                variables = c('B01003_001'), # Total Population
                year=2010,
                survey = 'acs1')

acs5 <- get_acs(geography = "county",
                variables = c('B01003_001'), # Total Population
                year=2010,
                survey = 'acs5')

nrow(acs1)/nrow(acs5)
# around 1/4


## trying different geographies

# to get tract or block group you need to specify the state, so I'd need to cycle through all states to get all the tracts and block groups.
acs_tract <- get_acs(geography = "tract",
                     state = 'MA',
                variables = c('B01003_001'), # Total Population
                year=2016,
                survey = 'acs1')

# block group is not available for ACS before 2013
acs_block_group <- get_acs(geography = "block group",
                     state = 'MA',
                     variables = c('B01003_001'), # Total Population
                     year=2016,
                     survey = 'acs5')


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

ma_pep <- get_estimates(
  geography = "county",
  product = "population",
  time_series = TRUE) %>% filter(DATE == 3) %>% filter(variable == "POP")

# choose date == 3 to map to July 1, 2010 based off this info: https://www.census.gov/data/developers/data-sets/popest-popproj/popest/popest-vars/2019.html

ma_pep<-as.data.frame(ma_pep)
ma_pep<-ma_pep[order(ma_pep$GEOID),]

##################################
## 4. merge ACS and census data ##
##################################

## read each processed dataset ##
# load('acs_data2206.RData')
# 
# load('ce_data2206.RData')
# 
# ## merge them by fips code, race, age, and sex ##
# Not by race, age, and sex, right?

ma_acs <- ma_acs %>%
  select(GEOID, NAME, acs = estimate)
ma_ce <- ma_ce %>%
  select(GEOID, NAME, census = value)
ma_pep <- ma_pep %>%
  select(GEOID, NAME, pep = value)

adat<-merge(ma_acs, ma_ce, by=c('GEOID','NAME'), all = TRUE) %>%
  merge(ma_pep, by=c('GEOID','NAME'), all = TRUE)




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



