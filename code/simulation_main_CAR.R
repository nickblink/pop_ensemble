library(dplyr)
library(rstan)

rstan_options(auto_write = TRUE)

if(file.exists('C:/Users/Admin-Dell')){
  setwd("C:/Users/Admin-Dell/Documents/github_projects/pop_ensemble/")
}else{
  setwd("C:/Users/nickl/Documents/github_projects/pop_ensemble/")
}
source('code/extra_functions_CAR.R')


# pull in the data
D2010 = read.csv('data/merged_wp_census_data2_081122.csv')
county_adj = read.csv('data/countyadj2.csv', row.names = 1)
models = c('pep', 'worldpop')


# subset data by state
lst <- subset_data_by_state(D2010, county_adj, 'New York', 'NY')


# Now I want to be able to simulate all sorts of data:
# - different CAR precision matrices (could be an input? Like Leroux or BYM)
# - different # of models
# - simulated values of models
# - being able to scale down the values in the data
# - adding some noise to the models (no need right now)

# simulate the data
phi_true, u_true, data = simulate_data(data_NY[:],
                                       adj_NY[:], 
                                       sim_numbers = False,
                                       scale_down = 1,
                                       poisson_noise = False,
                                       pivot = -1, 
                                       one_model = False, 
                                       models = models)


# fit the model


# analysis
# pull the gradients, log likelihood values
# plot the gradients and the log posterior, splitting likelihood and prior

# show the ESS and rhat

# plot density of phi samples and density of phi averages

# plot density of u values

# plot u values fitted vs. true

# plot predictions vs. true


### (later) plot the chloroploth maps