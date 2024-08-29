# ACS 2019, PEP 2019, WP 2019, FB 2019, Census 2020 

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