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

ma_pep <- get_estimates(
  geography = "county",
  product = "population",
  time_series = TRUE) %>% filter(DATE == 3) %>% filter(variable == "POP")
ma_pep<-as.data.frame(ma_pep)
ma_pep<-ma_pep[order(ma_pep$GEOID),]

##################################
## 4. merge ACS and census data ##
##################################

## read each processed dataset ##
# load('acs_data2206.RData')
# 
# load('ce_data2206.RData')

ma_acs <- ma_acs %>%
  select(GEOID, NAME, acs = estimate)
ma_ce <- ma_ce %>%
  select(GEOID, NAME, census = value)
ma_pep <- ma_pep %>%
  select(GEOID, NAME, pep = value)

## Fixing mismatched names
ma_pep$NAME <- gsub('Petersburg Borough, Alaska', 'Petersburg Census Area, Alaska', ma_pep$NAME)
ma_pep$NAME <- gsub('LaSalle Parish, Louisiana', 'La Salle Parish, Louisiana', ma_pep$NAME)
ma_pep$NAME[grep('Ana County, New Mexico', ma_pep$NAME)] <- 'Dona Ana County, New Mexico'

ma_acs$NAME[grep('Ana County, New Mexico', ma_acs$NAME)] <- 'Dona Ana County, New Mexico'
# this county was renamed in 2015
ind <- grep('Shannon County, South Dakota', ma_acs$NAME)
ma_acs$GEOID[ind] <- '46102'
ma_acs$NAME[ind] <- 'Oglala Lakota County, South Dakota'

ind <- grep('Shannon County, South Dakota', ma_ce$NAME)
ma_ce$GEOID[ind] <- '46102'
ma_ce$NAME[ind] <- 'Oglala Lakota County, South Dakota'

# this county was also renamed in 2015
ind <- grep('Wade Hampton Census Area, Alaska', ma_acs$NAME)
ma_acs$GEOID[ind] <- '02158'
ma_acs$NAME[ind] <- 'Kusilvak Census Area, Alaska'

ind <- grep('Wade Hampton Census Area, Alaska', ma_ce$NAME)
ma_ce$GEOID[ind] <- '02158'
ma_ce$NAME[ind] <- 'Kusilvak Census Area, Alaska'

adat<-merge(ma_acs, ma_ce, by=c('GEOID','NAME'), all = TRUE) %>%
  merge(ma_pep, by=c('GEOID','NAME'), all = TRUE) %>%
  na.omit() # omits one row with "bedford city"




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



