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
load('../data/merged_wp_census_data_270922.RData')

# Pull in Facebook data
load('../data/xx1.RData')
fb_df <- data.frame(GEOID = counties$GEOID,
                    fb = sum_overct[,1])

# Pull in Facebook data (not done)

# Merge
df <- merge(df, fb_df, by = 'GEOID')

# checking NAs
apply(df, 2, function(xx) {sum(is.na(xx))})

# checking correlation
cor(df[,c('acs','census','pep','worldpop', 'fb')])

# checking means 
colMeans(df[,c('acs','census','pep','worldpop', 'fb')])

# save results
# save(df, file = '../data/merged_fb_census_data_280922.RData')
