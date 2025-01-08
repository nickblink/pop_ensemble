# ACS 2019, PEP 2019, WP 2019, FB 2019, Census 2020 

library(tidycensus)
library(tidyverse)
library(tigris)
library(sp)
library(stringr)
library(haven)
library(reshape2)
library(rstudioapi)

current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path))

#### Checking AIAN PEP Vintage vs. Year ####
load('../../data/census_ACS_PEP_WP_AIAN_10282024.RData')

# PEP American Indian data <- now with VINTAGE!!
PEP_raw <- get_estimates(
  geography = "county",
  product = "characteristics",
  breakdown = "RACE",
  breakdown_labels = TRUE,
  vintage = 2019,
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
    PEP_vintage = value
  ) %>%
  # Sort by population in descending order
  arrange(desc(PEP_vintage))

test <- merge(df_AIAN %>% select(GEOID, census, acs, wp, pep),
              PEP_AIAN, by = 'GEOID')

lm.fit <- lm(census ~ acs + wp + pep, data = df_AIAN)

### Testing things
# PEP American Indian data <- now with VINTAGE!!
PEP_raw <- get_estimates(
  geography = "county",
  product = "characteristics",
  breakdown = "RACE",
  breakdown_labels = TRUE,
  vintage = 2024,
  year = 2019,
  state = 'MA'  # NULL gets all states
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
    PEP_vintage = value
  ) %>%
  # Sort by population in descending order
  arrange(desc(PEP_vintage))



#
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
load('../../data/census_ACS_PEP_WP_cleaned_08292024.RData')

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

# save(df, adjacency, file = '../../data/census_ACS_PEP_WP_wDensity_11152024.RData')

### Now with AIAN
load('../../data/census_ACS_PEP_WP_wDensity_11152024.RData')
df_full <- df
load('../../data/census_ACS_PEP_WP_AIAN_cleaned_10282024.RData')

counties <- counties(cb = TRUE, year = 2019)

# Extract county name, FIPS code, and land area
land_area <- counties %>%
  sf::st_set_geometry(NULL) %>% # Remove geometry to keep only data
  select(GEOID, ALAND) %>%
  mutate(ALAND_sqmi = ALAND * 0.0000003861)

df = merge(df, land_area, by = 'GEOID')
df$acs_density <- df$acs/df$ALAND_sqmi
df$pep_density <- df$pep/df$ALAND_sqmi

# Get the proportions.
df2 <- merge(df_full %>% select(GEOID, acs_full = acs, pep_full = pep, acs_density_full = acs_density, pep_density_full = pep_density), 
              df %>% select(GEOID, acs, pep))
df2$acs_AIAN_proportion <- df2$acs/df2$acs_full
df2$pep_AIAN_proportion <- df2$pep/df2$pep_full

df <- merge(df, df2 %>% select(GEOID, acs_AIAN_proportion, pep_AIAN_proportion, acs_density_full, pep_density_full))

# accidentally overwrote the old data anyway.
# save(df, adjacency, file = '../../data/census_ACS_PEP_WP_AIAN_wDensity_11152024.RData')


#### Creating AIAN death data ####
# Pull in COVID data.
cvd <- read.table('C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Population Estimation/Data/Wonder CDC Underlying Cause of Death, 2018-2022, Single Race Non Hispanic.txt', 
                  header = T, 
                  sep = '\t', )
cvd$County.Code <- as.character(cvd$County.Code)
# adding a 0 to the four digit codes that are supposed to have one.
cvd$County.Code <- ifelse(nchar(cvd$County.Code) == 4,
                          paste0('0', cvd$County.Code), 
                          cvd$County.Code)

cvd <- cvd %>% 
  filter(Year.Code == 2020,
         Single.Race.6 == 'American Indian or Alaska Native') %>%
  select(GEOID = County.Code, AIAN_deaths = Deaths)

# save(cvd, file = '../../data/AmericanIndian_COVID_Deaths_2020.RData')


#### AIAN 2018 data ####
load('../../data/census_ACS_PEP_WP_AIAN_wDensity_11152024.RData')

# ACS 2018 data.
ACS_AIAN <- get_acs(
  geography = "county",
  variables = "B02001_004", # Code for American Indian and Alaska Native alone
  year = 2018#,
  #survey = "acs5" # Using 5-year estimates for more reliable county-level data
) %>%
  select(GEOID,
         acs_2018 = estimate)

# PEP 2018 data.
PEP_raw <- get_estimates(
  geography = "county",
  product = "characteristics",
  breakdown = "RACE",
  breakdown_labels = TRUE,
  year = 2018,
  state = NULL  # NULL gets all states
)


# PEP 2018 data.
PEP_raw1 <- get_estimates(
  geography = "county",
  product = "characteristics",
  breakdown = "RACE",
  breakdown_labels = TRUE,
  year = 2015,
  vintage = 2019,
  state = NULL  # NULL gets all states
)
# PEP 2018 data.
PEP_raw2 <- get_estimates(
  geography = "county",
  product = "characteristics",
  breakdown = "RACE",
  breakdown_labels = TRUE,
  year = 2015,
  vintage = 2015,
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
    pep_2018 = value
  ) %>%
  # Sort by population in descending order
  arrange(desc(pep_2018))

df2 <- df %>%
  merge(ACS_AIAN, by = 'GEOID') %>%
  merge(PEP_AIAN, by = 'GEOID')

df2 <- df2 %>%
  mutate(pep_diff = pep - pep_2018,
         acs_diff = acs - acs_2018)

df <- df2
save(df, adjacency, file = '../../data/census_ACS_PEP_WP_AIAN_wDensity_and2018_01022024.RData')

#
#### PEP AND ACS 2018 data ####
load('../../data/census_ACS_PEP_WP_wDensity_11152024.RData')

# comparing vintage estimates
{
PEP_2018 <- get_estimates(geography = "county",
                     product = "population",
                     time_series = TRUE, 
                     vintage = 2018) %>%
  filter(DATE == 11, variable == 'POP')

PEP_2019 <- get_estimates(geography = "county",
                          product = "population",
                          time_series = TRUE, 
                          vintage = 2019) %>%
  filter(DATE == 11, variable == 'POP')

test = merge(PEP_2018[,c('GEOID','value')],
             PEP_2019[,c('GEOID','value')],
             by = 'GEOID')
# Cor is high.
cor(test$value.x, test$value.y)

# MAPE is 0.3%. Ok
mean(abs(test$value.x - test$value.y)/test$value.x)
}

# getting 2018 values
PEP_2018 <- get_estimates(geography = "county",
                          product = "population",
                          time_series = TRUE, 
                          vintage = 2018) %>%
  filter(DATE == 11, variable == 'POP') %>%
  select(GEOID, pep_2018 = value)

ACS_2018 <- get_acs(geography = "county",
               variables = 'B01003_001', # Total Population
               year=2018) %>%
  select(GEOID, acs_2018 = estimate)

df2 <- df %>% 
  merge(PEP_2018, by = 'GEOID') %>%
  merge(ACS_2018, by = 'GEOID')

df2 <- df2 %>%
  mutate(pep_diff = pep - pep_2018,
         acs_diff = acs - acs_2018)

df <- df2
# save(df, adjacency, file = '../../data/census_ACS_PEP_WP_wDensity_and2018_01022024.RData')

# regression analysis of the variables.
lm.fit.1 <- lm(census ~ pep + acs + wp, data = df2)
summary(lm.fit.1)
# ok so ACS does seem to actually be useful here. Interesting.

lm.fit.2 <- lm(census ~ pep, data = df2)
summary(lm.fit.2)

lm.fit.3 <- lm(census ~ pep_2018 + acs_2018 + wp, data = df2)
summary(lm.fit.3)

lm.fit.4 <- lm(census ~ pep_2018 + acs_2018, data = df2)
summary(lm.fit.4)
# interesting. This makes sense though - it's sort of getting at a trend in the data that is showing the difference in PEP and ACS to show the trends over the years.

lm.fit.5 <- lm(census ~ pep + acs + pep_2018 + acs_2018 + +wp, data = df2)
summary(lm.fit.5)

lm.fit.6 <- lm(census ~ acs + wp, data = df2)
summary(lm.fit.6)
# hm not great.

###

#### PCA on the data? ####
load('../../data/census_ACS_PEP_WP_wDensity_and2018_01022024.RData')

X <- df %>%
  select(acs, pep, wp, acs_2018, pep_2018)

#pca_result <- prcomp(X, center = TRUE, scale. = TRUE)
pca_result <- prcomp(X, center = F, scale. = F)

# Summary of PCA
summary(pca_result)

head(pca_result$x)

pca_result$rotation
# ah that makes sense. The first PC is the common direction of all of them (towards the magnitude of census), and then the following PCs differentiate them.
#
plot(pca_result, type = "l", main = "Scree Plot")
biplot(pca_result, scale = 0)

PC_cols <- pca_result$x
colnames(PC_cols) <- colnames(PC_cols) %>% tolower()

# make all columns positive on average.
PC_cols_pos <- PC_cols
for(j in 1:ncol(PC_cols)){
  if(sum(PC_cols[,j]) < 0){
    PC_cols_pos[,j] <- -PC_cols_pos[,j]
  }
}

df <- cbind(df, PC_cols_pos)

# save(df, adjacency, file = '../../data/census_ACS_PEP_WP_wDensity_and2018_01022024.RData')


#
#### Regression Analysis AIAN ####
load('../../data/census_ACS_PEP_WP_AIAN_wDensity_11152024.RData')
lm.fit.1 <- lm(census ~ pep + acs + wp, data = df_AIAN)
summary(lm.fit.1)
# interesting. 

lm.fit.2 <- lm(census ~ pep + acs, data = df_AIAN)
summary(lm.fit.2)
# So PEP is really dominating here. Why? 

lm.fit.3 <- lm(log1p(census) ~ log1p(pep) + log1p(acs) + log1p(wp), data = df_AIAN)
summary(lm.fit.3)
# ok. Still not sure what to make of that. Darn darn.

lm.fit.5 <- lm(census ~ acs + wp, data = df_AIAN)
summary(lm.fit.5)
# interesting. WP could be negative?
