# ACS 2019, PEP 2019, WP 2019, FB 2019, Census 2020 

library(tidycensus)
library(tidyverse)
library(tigris)
library(sp)
library(stringr)
library(haven)
library(reshape2)
library(rstudioapi)
library(dplyr)

current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path))

#### Census 2020 ####
load_variables(2020, "pl")
census <- get_decennial(
  geography = "county",
  variables = "P1_001N", # Total population
  year = 2020
)

#### ACS 2019 ####
ACS <- get_acs(geography = "county",
          variables = 'B01003_001', # Total Population
          year=2019)

#### PEP 2019 ####
PEP <- get_estimates(geography = "county",
              product = "population",
              time_series = TRUE, 
              vintage = 2019) %>%
  filter(DATE == 12, variable == 'POP')

#### WP 2019 ####
load("C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Population Estimation/Data/wp_2019.RData")

#### FB 2019 ####

#### Bring together ####

# Change value names to be each method output.
census <- census %>%
  dplyr::select(GEOID, NAME, census = value)
ACS <- ACS %>% 
  dplyr::select(GEOID, NAME, acs = estimate)
PEP <- PEP %>%
  dplyr::select(GEOID, NAME, pep = value)
WP <- counties %>%
  as.data.frame() %>%
  dplyr::select(GEOID, NAME, wp = WP, wp_un = WP_UN)

# fix naming of counties to match each other

# merge!
df <- merge(ACS, census, by=c('GEOID','NAME'), all = TRUE) %>%
  merge(PEP, by=c('GEOID','NAME'), all = TRUE) %>%
  merge(WP, by=c('GEOID'), all = TRUE)

# tt <- sapply(1:3222, function(i){
#   grepl(df$NAME.y[i], df$NAME.x[i])
# }) # good
df$NAME <- df$NAME.x
df$NAME.y <- df$NAME.x <- NULL

#### Filter to mainland US ####
tg <- tigris::fips_codes
tg$fips <- paste0(tg$state_code, tg$county_code)

# look at unique state codes
unique(tg[,c('state','state_code')])

# remove non-mainland states
tg <- tg %>% 
  filter(!(state_code %in% c('02','15', '60','66','69','72','74','78')))

df <- df %>%
  filter(GEOID %in% tg$fips)

df <- tg %>% 
  dplyr::select(state = state_name, GEOID = fips) %>%
  merge(df, by = 'GEOID')

# look into places still missing.
NA_rows <- apply(df, 1, function(xx) any(is.na(xx)))
df[NA_rows,]
# good

# gucci goo. Pausing here.
# save(df, file = '../../data/census_ACS_PEP_WP_08292024.RData')

cor(df[,3:7])
colSums(df[,3:7])
# ok looking good.
#### American Indian data ####
# Set options for tidycensus
# options(tigris_use_cache = TRUE)  # Cache spatial data
# options(timeout = 300)            # Increase timeout limit

v19 <- load_variables(2019, "acs1", cache = TRUE)
ind <- grep('B02001', v19$name)
v19[ind,]
v19$label[v19$name == 'B02001_004']

# ACS American Indian data.
ACS_AIAN <- get_acs(
  geography = "county",
  variables = "B02001_004", # Code for American Indian and Alaska Native alone
  year = 2019#,
  #survey = "acs5" # Using 5-year estimates for more reliable county-level data
) %>%
  select(GEOID,
         ACS = estimate)

# PEP American Indian data
PEP_raw <- get_estimates(
  geography = "county",
  product = "characteristics",
  breakdown = "RACE",
  breakdown_labels = TRUE,
  year = 2019,
  state = NULL  # NULL gets all states
)

# Clean and format the PEP data
PEP_AIAN <- PEP_raw %>%
  # Filter for AIAN alone
  filter(RACE == "American Indian and Alaska Native alone") %>%
  # Create clean names
  separate(NAME, into = c("county_name", "state"), sep = ", ") %>%
  # Select and rename columns
  select(
    GEOID,
    PEP = value
  ) %>%
  # Sort by population in descending order
  arrange(desc(PEP))

# Census data
# Get 2020 Census data for American Indian and Alaska Native population
census_AIAN <- get_decennial(
  geography = "county",
  variables = "P1_005N",  # American Indian and Alaska Native alone
  year = 2020,
  sumfile = "pl"  # PL 94-171 Redistricting Data
) %>%
  select(GEOID,
         census = value)
# ARGH small counts are bad.

# load in WP. Scale based off full population values.
load('../../data/census_ACS_PEP_WP_cleaned_08292024.RData')
WP_AIAN <- merge(df, ACS_AIAN, by = 'GEOID')
WP_AIAN$WP <- WP_AIAN$wp/WP_AIAN$acs*WP_AIAN$ACS 
WP_AIAN <- WP_AIAN %>%
  select(GEOID, state, NAME, WP)

# merge!
df_AIAN <- merge(WP_AIAN, census_AIAN, by=c('GEOID')) %>%
  merge(PEP_AIAN, by=c('GEOID')) %>%
  merge(ACS_AIAN, by=c('GEOID')) %>%
  # select names to exactly match names in other dataset.
  select(GEOID, state, acs = ACS, census, pep = PEP, wp = WP, NAME)

# save(df_AIAN, adjacency, file = '../../data/census_ACS_PEP_WP_AIAN_10282024.RData')
