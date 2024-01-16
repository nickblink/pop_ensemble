library(dplyr)
library(rstan)

rstan_options(auto_write = TRUE)

## set working directory
if(file.exists('C:/Users/Admin-Dell')){
  setwd("C:/Users/Admin-Dell/Documents/github_projects/pop_ensemble/")
}else{
  setwd("C:/Users/nickl/Documents/github_projects/pop_ensemble/")
}
source('code/extra_functions_CAR.R')


## pull in the data
D2010 = read.csv('data/merged_wp_census_data2_081122.csv')
county_adj = read.csv('data/countyadj2.csv', row.names = 1)
models = c('pep', 'worldpop')


## subset data by state
NY_lst <- subset_data_by_state(D2010, county_adj, 'New York', 'NY')

# Now I want to be able to simulate all sorts of data:
# - different CAR precision matrices (could be an input? Like Leroux or BYM)
# - different # of models
# - simulated values of models
# - being able to scale down the values in the data
# - adding some noise to the models (no need right now)


## simulate the data
# data2 <- simulate_models(data = NY_lst$data, models = c('acs','pep'), means = c(100, 200), variances = c(10^2, 10^2))
data_lst <- simulate_data(NY_lst$data, NY_lst$adjacency, models = models, precision_type = 'Leroux', tau2 = 1, rho = 0.3)


## run MAP estimation?


## fit the model
# function to fit:
### Fits the CAR model using rstan. Prepares the data for rstan and runs it.
# data: input data with output column y and covariate columns according to models
# adjacency: the adjacency matrix for the data
# models: the models for the ensemble.
# precision_type: Cressie or Leroux, for the type of precision matrix
# pivot: what pivot index to use for data creation (-1 indicates no pivot)

run_stan_CAR <- functon(data, adjacency, models = c('acs','pep','worldpop'), precision_type = 'Cressie', pivot = -1, run_MAP = F, seed = 10){
  
}

# create the posterior

# create stan data prep

# analysis
# pull the gradients, log likelihood values
# plot the gradients and the log posterior, splitting likelihood and prior

# show the ESS and rhat

# plot density of phi samples and density of phi averages

# plot density of u values

# plot u values fitted vs. true

# plot predictions vs. true


### (later) plot the chloroploth maps