library(dplyr)
library(rstan)
library(ggplot2)

rstan_options(auto_write = TRUE)

# set working directory
if(file.exists('C:/Users/Admin-Dell')){
  setwd("C:/Users/Admin-Dell/Documents/github_projects/pop_ensemble/")
}else{
  setwd("C:/Users/nickl/Documents/github_projects/pop_ensemble/")
}

# load extra functions
source('code/extra_functions_CAR.R')

#### Running with 2 and 3 models ####

## pull in the data
D2010 = read.csv('data/merged_wp_census_data2_081122.csv')
county_adj = read.csv('data/countyadj2.csv', row.names = 1)
models = c('acs','pep')

## subset data by state
NY_lst <- subset_data_by_state(D2010, county_adj, 'New York', 'NY')

## Running with simulated model and y values - 2 models
data <- simulate_models(data = NY_lst$data, models = c('acs','pep'), means = c(100, 200), variances = c(10^2, 10^2))

data_lst <- simulate_data(data, NY_lst$adjacency, models = models, precision_type = 'Leroux', tau2 = 1, rho = 0.3)

stan_fit <- run_stan_CAR(data_lst$data, data_lst$adjacency, models = models, n.sample = 10000, burnin = 5000)

# save(stan_fit, data_lst, file = 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Population Estimation/Results/results_sim_01232024.RData')

pp <- process_results(data_lst, models, stan_fit)

# ggsave(pp, file = 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Population Estimation/Figures/ACS100_PEP200_01242024.png', height = 10, width = 5)


## Running with simulated model and y values - 3 models
models = c('acs','pep','worldpop')
data <- simulate_models(data = NY_lst$data, models = models, means = c(100, 200, 300), variances = c(10^2, 10^2, 10^2))

data_lst <- simulate_data(data, NY_lst$adjacency, models = models, precision_type = 'Leroux', tau2 = 1, rho = 0.3)

stan_fit <- run_stan_CAR(data_lst$data, data_lst$adjacency, models = models, n.sample = 10000, burnin = 5000)

pp <- process_results(data_lst, models, stan_fit)

# ggsave(pp, file = 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Population Estimation/Figures/ACS100_PEP200_WP300_01242024.png', height = 10, width = 7)

#
#### Running repeatedly - 2 models ####
res_lst <- list()
N_iters = 10
for(i in 1:N_iters){
  data <- simulate_models(data = NY_lst$data, models = models, means = c(100, 200), variances = c(10^2, 10^2), seed = i)
  
  data_lst <- simulate_data(data, NY_lst$adjacency, models = models, precision_type = 'Leroux', tau2 = 1, rho = 0.3, seed = i)
  
  stan_fit <- run_stan_CAR(data_lst$data, data_lst$adjacency, models = models, n.sample = 10000, burnin = 5000, seed = i)
  
  res_lst[[i]] <- list(data_list = data_lst, stan_fit = stan_fit)
}

# save(res_lst, file = 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Population Estimation/Results/results_2models_10runs_01252024.RData')

plot_list <- list()
for(i in 1:N_iters){
  pp <- process_results(res_lst[[i]]$data_list, models, res_lst[[i]]$stan_fit, ESS = F, likelihoods = F, phi_estimates = F, u_estimates = F, y_estimates = F)
  plot_list[[i]] <- pp
}

tt <- plot_grid(plotlist = plot_list, ncol = 2)
# ggsave(plot = tt, filename = 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Population Estimation/Figures/spatial_params_2models_10runs_01252024.png', width = 10, height = 15)

# 
#### Running repeatedly - 3 models ####
models = c('acs','pep','worldpop')
res_lst <- list()
N_iters = 10
for(i in 1:N_iters){
  data <- simulate_models(data = NY_lst$data, models = models, means = c(100, 200, 300), variances = c(10^2, 10^2, 10^2), seed = i)
  
  data_lst <- simulate_data(data, NY_lst$adjacency, models = models, precision_type = 'Leroux', tau2 = 1, rho = 0.3, seed = i)
  
  stan_fit <- run_stan_CAR(data_lst$data, data_lst$adjacency, models = models, n.sample = 10000, burnin = 5000, seed = i)
  
  res_lst[[i]] <- list(data_list = data_lst, stan_fit = stan_fit)
}

# save(res_lst, file = 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Population Estimation/Results/results_3models_10runs_01252024.RData')

plot_list <- list()
for(i in 1:N_iters){
  pp <- process_results(res_lst[[i]]$data_list, models, res_lst[[i]]$stan_fit, ESS = F, likelihoods = F, phi_estimates = F, u_estimates = F, y_estimates = F)
  plot_list[[i]] <- pp
}

tt <- plot_grid(plotlist = plot_list, ncol = 2)
#ggsave(plot = tt, filename = 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Population Estimation/Figures/spatial_params_3models_10runs_01252024.png', width = 10, height = 15)

# 
#### Testing precision matrix determinant calculation ####
# determinant from two models
stan_out <- extract(stan_fit)

# get the 2-D matrix of the determinant
stan_det <- stan_out$log_detQ

ind <- 20
tau2 <- stan_out$tau2[ind,]
rho <- stan_out$rho[ind,]
Q <- generate_precision_mat(W = data_lst$adjacency, type = 'Leroux', tau2 = tau2[1], rho = rho[1]) 

ldets <- c()
for(i in 1:10){
  tau2 <- stan_out$tau2[i,1]
  rho <- stan_out$rho[i,1]
  Q <- generate_precision_mat(W = data_lst$adjacency, type = 'Leroux', tau2 = tau2, rho = rho) 
  ldets[i] <- log(det(as.matrix(Q)))
}

stan_det[1:10,1]
ldets
stan_det[1:10,1] - ldets
# PHEW THIS IS GOOD

#### Testing data generation tau2 ####
data_lst <- simulate_data(data, NY_lst$adjacency, models = models, precision_type = 'Leroux', tau2 = c(1,10), rho = 0, seed = 18)

head(data_lst$phi_true)
var(data_lst$phi_true$acs)
var(data_lst$phi_true$pep)
# ok good

#### Checking the MVN sampling ####
Q <- generate_precision_mat(W = data_lst$adjacency, type = 'Leroux', tau2 = 1, rho = 0.9) 
par(mfrow = c(3,3))

for(i in 1:9){
  test <- sample_MVN_from_precision(n = 1, Q = Q)
  Sigma <- solve(Q)
  test2 <- MASS::mvrnorm(n = 1, mu = rep(0,nrow(Sigma)), Sigma = Sigma)
  plot(density(test))
  lines(density(test2), col = 'red')
  print(cor(test, test2))
}

seed = 23
# cholesky
data_lst1 <- simulate_data(data, NY_lst$adjacency, models = models, precision_type = 'Leroux', tau2 = 1, rho = 0.3, seed = seed)

stan_fit1 <- run_stan_CAR(data_lst1$data, data_lst1$adjacency, models = models, n.sample = 10000, burnin = 5000)

pp1 <- process_results(data_lst1, models, stan_fit1)

# MVRnorm
data_lst2 <- simulate_data(data, NY_lst$adjacency, models = models, precision_type = 'Leroux', tau2 = 1, rho = 0.3, cholesky = F, seed = seed)

stan_fit2 <- run_stan_CAR(data_lst2$data, data_lst2$adjacency, models = models, n.sample = 10000, burnin = 5000)

pp2 <- process_results(data_lst2, models, stan_fit2)

tt = plot_grid(pp1, pp2, nrow = 1)
tt

#
#### previous runs ####
{
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
}


### (later) plot the chloroploth maps