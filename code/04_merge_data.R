#####
# Merge the census data products, worldpop, and facebook data into one data frame for analysis

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

# Pull in census data products 
load('../data/census_products_data220822.RData')

# Pull in Worldpop data
load('../data/wp_1km_county_aggregate_2010.RData')
wp_df <- data.frame(GEOID = counties$GEOID,
                    worldpop = sum_overct[,1])

# Pull in Facebook data (not done)

# Merge
df <- merge(adat, wp_df, by = 'GEOID')

# checking NAs
apply(df, 2, function(xx) {sum(is.na(xx))})

# checking correlation
cor(df[,c('acs','census','pep','worldpop')])

# checking means 
colMeans(df[,c('acs','census','pep','worldpop')])

# save results
# save(df, file = '../data/merged_wp_census_data_270922.RData')
