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

##### Trying the PEP API #####
library(httr)
library(jsonlite)
library(censusapi)


url <- 'https://api.census.gov/data/2000/pep/int_population'

api_call <- httr::GET(url)

api_char <- rawToChar(api_call$content)

api_JSON <- jsonlite::fromJSON(api_char, flatten = T)

# hmm not sure this gives me anything though. I tried all the links and nothing is showing up

df <- listCensusApis()


grep('population_estimat', df$title, value = T, ignore.case = T)

df %>% 
  filter(title == '2000 Population Estimates - 2000-2010 Intercensal Estimates: Population')

tt <- getCensus(name = 'pep/int_population', key = api_census_key, vars = 'test')

makeVarlist(name = 'pep/int_population')
# oh boy

listCensusMetadata(name = 'pep/int_population')
# not available