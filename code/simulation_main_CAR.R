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

## fit the model
stan_fit <- run_stan_CAR(data_lst$data, data_lst$adjacency, models = models)


# debug and analyze the model
shinystan::launch_shinystan(stan_fit)

stan_out <- extract(stan_fit)

tt <- stan_out$mu
which(apply(tt, 1, function(xx) any(is.na(xx))))
# [94,55] is one place

ss <- stan_out$u
ss[94, 55, ]
# NaN 0. Ok ok 

vv <- stan_out$phi
vv[94, 55, ]
# 709.5587 -252.6754
# ok in general I am getting crazy big values for phi.
# so no wonder exponentiating these gets crazy big.

# what is tau2?
tau <- stan_out$tau2
colMeans(tau)
# 97193.93 100957.00
# ok obviously something is wrong.

aa <- stan_out$log_detQ
plot(density(aa))
# so on the log scale not much variation. But that's huge on the not log scale, right? No, it's super super small. Because of the tau2, amiright?

bb <- stan_out$ldet_vec
# ok so the tau2 term dominates

cc <- stan_out$lp__
# ah so it's friggin huge.

stan_summary = summary(stan_fit, pars = c('tau2','rho', 'phi'))$summary
stan_lst <- list(stan_fit = stan_fit,
                 stan_out = stan_out,
                 stan_summary = stan_summary)

stan_summary = summary(stan_list$stan_fit, pars = c('tau2','rho', 'phi', 'u'))$summary


### Trying with smaller values
data_lst2 <- simulate_data(NY_lst$data, NY_lst$adjacency, models = models, precision_type = 'Leroux', tau2 = 1, rho = 0.3, scale_down = 1000)

## fit the model
stan_fit2 <- run_stan_CAR(data_lst2$data, data_lst2$adjacency, models = models)

### Trying with different y values
data2 <- simulate_models(data = NY_lst$data, models = c('acs','pep'), means = c(100, 200), variances = c(10^2, 10^2))
data_lst3 <- simulate_data(data2, NY_lst$adjacency, models = models, precision_type = 'Leroux', tau2 = 1, rho = 0.3)

stan_fit3 <- run_stan_CAR(data_lst3$data, data_lst3$adjacency, models = models)


# create the posterior
process_results <- function(data, models, stan_list){
  ## get the ESS of phi, u, and others
  
  ## get the posterior split into likelihood and log posterior
  
  ## compare the true phi values with the estimated phi values
  
  ## compare the true u values with the estimated u values
  
  ## compare the true outcomes with the estimated outcomes (compared to just using one model in the outcomes)
}

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