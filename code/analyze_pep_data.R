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

pep <- get_estimates(geography = "county",
  product = "population",
  time_series = TRUE) %>% 
  dplyr::filter(variable == "POP") %>%
  rename(pep = value)

res <- data.frame(DATE = 1:12, cor = NA, prop_bias = NA, MAPE = NA)
for(i in 1:12){
  pepyear <- pep %>% dplyr::filter(DATE == i)
  df <- merge(ma_ce, pepyear, by = c('NAME','GEOID'))
  res[i,2] <-  cor(df$pep, df$value)
  res[i,3] <- sum(df$pep)/sum(df$value)
  res[i,4] <- mean(abs(df$pep - df$value)/df$value)
}

res

# ok so it's almost certainly based on the census
