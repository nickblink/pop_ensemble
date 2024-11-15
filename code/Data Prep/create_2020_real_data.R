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

#### Creating density data ####
## Full pop
load('../../data/census_ACS_PEP_WP_08292024.RData')

counties <- counties(cb = TRUE, year = 2019)

# Extract county name, FIPS code, and land area
land_area <- counties %>%
  sf::st_set_geometry(NULL) %>% # Remove geometry to keep only data
  select(GEOID, ALAND) %>%
  mutate(ALAND_sqmi = ALAND * 0.0000003861)

# If using later data (like 2022), the 8 counties of CT are excluded because they don't report their population numbers directly. These have to be gathered from online if they are going to be used.
df = merge(df, land_area, by = 'GEOID')
df$acs_density <- df$acs/df$ALAND_sqmi
df$pep_density <- df$pep/df$ALAND_sqmi
# these are highly correlated. They will need to be logged.

# save(df, file = '../../data/census_ACS_PEP_WP_wDensity_11152024.RData')

### Now with AIAN
load('../../data/census_ACS_PEP_WP_AIAN_10282024.RData')

counties <- counties(cb = TRUE, year = 2019)

# Extract county name, FIPS code, and land area
land_area <- counties %>%
  sf::st_set_geometry(NULL) %>% # Remove geometry to keep only data
  select(GEOID, ALAND) %>%
  mutate(ALAND_sqmi = ALAND * 0.0000003861)

df_AIAN = merge(df_AIAN, land_area, by = 'GEOID')
df_AIAN$acs_density <- df_AIAN$acs/df_AIAN$ALAND_sqmi
df_AIAN$pep_density <- df_AIAN$pep/df_AIAN$ALAND_sqmi

# Get the proportions.
df2 <- merge(df %>% select(GEOID, acs_full = acs, pep_full = pep), 
              df_AIAN %>% select(GEOID, acs, pep))
df2$acs_AIAN_proportion <- df2$acs/df2$acs_full
df2$pep_AIAN_proportion <- df2$pep/df2$pep_full

df_AIAN <- merge(df_AIAN, df2 %>% select(GEOID, acs_AIAN_proportion, pep_AIAN_proportion))

# accidentally overwrote the old data anyway.
#save(df_AIAN, adjacency, file = '../../data/census_ACS_PEP_WP_AIAN_wDensity_11152024.RData')

