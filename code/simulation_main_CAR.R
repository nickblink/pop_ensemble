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
# data: input data with output column y and covariate columns according to models.
# adjacency: the adjacency matrix for the data.
# models: the models for the ensemble.
# precision_type: Cressie or Leroux, for the type of precision matrix.
# n.sample: the number of iterations to run the rstan code.
# burnin: the number of burnin iterations to run the rstan code.
# seed: a seed for reproducability
run_stan_CAR <- function(data, adjacency, models = c('acs','pep','worldpop'), precision_type = 'Leroux', n.sample = 10000, burnin = 5000, seed = 10){
  # error checking for precision matrix type
  if(precision_type != 'Leroux'){stop('only have Leroux precision coded')}
  
  # prep the data
  stan_data <- prep_stan_data_leroux_sparse(data, adjacency, models)
  
  # fit the stan model
  stan_fit <- stan(file = "code/CAR_leroux_sparse.stan",
                   data = stan_data, 
                   iter = n.sample, 
                   warmup = burnin,
                   chains = 1, 
                   init = '0',
                   cores = 1,
                   seed = seed)
  
  # extract important info
  stan_out <- extract(stan_fit)
  stan_summary = summary(stan_fit, pars = c('tau2','rho', 'phi'))$summary
  stan_lst <- list(stan_out = stan_out, 
                   stan_summary = stan_summary)
  return(stan_lst)
}

test <- run_stan_CAR(data_lst$data, data_lst$adjacency, models = models)

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