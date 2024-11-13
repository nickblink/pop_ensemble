
df <- read.csv('C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Population Estimation/Data/Monthly_COVID-19_Death_Rates_per_100_000_Population_by_Age_Group__Race_and_Ethnicity__Sex__and_Region_with_Double_Stratification_20241113.csv')
head(df)
length(unique(df$jurisdiction_residence))

df <- read.csv('C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Wastewater Surveillance/data/Provisional_COVID-19_Deaths_by_County__and_Race_and_Hispanic_Origin_20241113.csv')
length(unique(df$FIPS.Code)) # ~1/3 of counties?
# a lot of White. Not a lot of Am In and Al Native. Hispanic?

df <- read.csv('C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Wastewater Surveillance/data/Healthdata_Weekly_United_States_COVID-19_Cases_and_Deaths_by_County_-_ARCHIVED.csv')
df$date <- as.Date(df$date, format = '%m/%d/%Y')
head(df)
length(unique(df$fips_code))
# this looks good. Dates go from 1/22/2020 - 5/10/2023

df <- read.table('C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Population Estimation/Data/Wonder CDC COVID Underlying Cause of Death, 2018-2022, Single Race.txt', header = T, sep = '\t')

df2 <- read.table('C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Population Estimation/Data/Wonder CDC Underlying Cause of Death, 2018-2022, Single Race Non Hispanic.txt', header = T, sep = '\t')
tmp <- df2 %>% filter(Single.Race.6 == 'American Indian or Alaska Native')

# testing WW data
df <- read.table('C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Wastewater Surveillance/data/Wonder CDC Underlying Cause of Death, 2018-2022, Single Race.txt', header = T, sep = '\t')
sum(df$Deaths) # seems reasonable.
length(unique(df$County)) # also seems reasonable. But not all?
table(table(df$County)) # is this because of low counts or what? Judging from "Adair County, IA", NO. This is probably a lack of reporting.
