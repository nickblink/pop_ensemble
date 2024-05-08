library(dplyr)
library(rstan)
# library(ggplot2)

rstan_options(auto_write = TRUE)

# load extra functions
source('code/extra_functions_CAR.R')

## pull in the data
D2010 = read.csv('data/merged_wp_census_data2_081122.csv')
county_adj = read.csv('data/countyadj2.csv', row.names = 1)

## subset data by state
NY_lst <- subset_data_by_state(D2010, county_adj, 'New York', 'NY')

# parameters for simulations and MCMC fitting
models = c('X1','X2','X3')
n.sample = 1000
burnin = 500

#### 4/25/2024: Testing the CV blocking ####
# 5-fold validation.
system.time({
  res_lst <- multiple_sims(NY_lst, models, N_sims = 1, n.sample = n.sample, burnin = burnin, family = 'normal', use_softmax = T, # (shared params)
                           variances = c(10^2, 10^2, 10^2),  means = c(100,100,100), rho = 0.3, tau2 = 1, sigma2 = 10^2, # (DGP params)
                           sigma2_prior_shape = 50, sigma2_prior_rate = 0.5, tau2_prior_shape = 1, tau2_prior_rate = 1, num_y_samples = 3, stan_path = 'code/CAR_leroux_sparse_normal.stan') # (stan params)
})

save(res_lst, file = 'results/onerun_results_04302024.RData')