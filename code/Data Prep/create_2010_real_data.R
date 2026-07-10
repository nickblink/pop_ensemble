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

#### Census 2010 ####
load_variables(2020, "pl")
census_2010 <- get_decennial(
  geography = "county",
  variables = "P001001",  # Total population for 2010
  year = 2010
)
# Error right now - Idk why.

#### ACS 2009 ####
ACS <- get_acs(geography = "county",
          variables = 'B01003_001', # Total Population
          year=2009)
# These should be vintage according to Claude.

#### PEP 2009 ####
pep_vintage <- read_csv("../../data/pep_2000_2010_vintage_precensus.csv") %>%
  filter(SUMLEV == "050")  # County level only
# https://www.census.gov/programs-surveys/popest/technical-documentation/research/evaluation-estimates/2010-evaluation-estimates.html

# Sanity check: compare precensus vs postcensus 2009 estimates
pep_postcensus <- read_csv("../../data/pep_2000_2010_postcensus.csv") %>%
  filter(SUMLEV == 50)
# https://www.census.gov/programs-surveys/popest/data/data-sets.2009.List_221933341.html#list-tab-List_221933341

# Compare totals
cat("Precensus 2009 total:", sum(pep_vintage$POPESTIMATE072009), "\n")
cat("Postcensus 2009 total:", sum(pep_postcensus$POPESTIMATE2009), "\n")

# Compare totals
cat("Precensus 2010 total:", sum(pep_vintage$POPESTIMATE072010), "\n")
cat("Postcensus 2010 total:", sum(pep_postcensus$POPESTIMATE2010), "\n")
