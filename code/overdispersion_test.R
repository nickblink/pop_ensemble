library(tidycensus)
library(tigris)
library(sp)
library(stringr)
library(haven)
library(reshape2)
library(rstudioapi)
library(dplyr)
library(ggplot2)

current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path))

load('../data/tmp_merged_census_products_08132024.RData')

df$census_sqrt <- sqrt(df$census)
df$acs_resid <- df$acs - df$census
df$pep_resid <- df$pep - df$census

ggplot(df, aes(x = census_sqrt, y = abs(acs_resid))) +
  geom_point() +  # Plot the data points
  geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Add a line of best fit
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") + 
  scale_x_log10() +  # Apply log scale to x-axis
  scale_y_log10() +  # Apply log scale to y-axis
  labs(title = "ACS residual comparison",
       x = "Poisson sd: sqrt(census)",
       y = "Empirical sd: |acs-census|") +
  theme_minimal()


ggplot(df, aes(x = census_sqrt, y = abs(pep_resid))) +
  geom_point() +  # Plot the data points
  geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Add a line of best fit
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") + 
  scale_x_log10() +  # Apply log scale to x-axis
  scale_y_log10() +  # Apply log scale to y-axis
  labs(title = "PEP residual comparison",
       x = "Poisson sd: sqrt(census)",
       y = "Empirical sd: |PEP-census|") +
  theme_minimal()

get_dispersion <- function(poisson_model){
  dispersion_statistic <- sum(residuals(poisson_model, type = "pearson")^2) / poisson_model$df.residual
 dispersion_statistic 
}

# test dispersion of ACS
acs.fit <- glm(census ~ log(acs), family = poisson(), data = df)
#acs.fit <- glm(census ~ log(acs), family = quasipoisson(), data = df)
get_dispersion(acs.fit)

# test dispersion of PEP
pep.fit <- glm(census ~ log(pep), family = poisson(), data = df)
get_dispersion(pep.fit)

# test dispersion of ACS and PEP together
both.fit <- glm(census ~ log(acs) + log(pep), family = poisson(), data = df)
get_dispersion(both.fit)



