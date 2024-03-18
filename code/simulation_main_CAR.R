library(dplyr)
library(rstan)
library(ggplot2)

rstan_options(auto_write = TRUE)

# set working directory
if(file.exists('C:/Users/Admin-Dell')){
  root_dir = 'C:/Users/Admin-Dell'
}else{
  root_dir = 'C:/Users/nickl'
}
setwd(sprintf('%s/Documents/github_projects/pop_ensemble/', root_dir))

# load extra functions
source('code/extra_functions_CAR.R')

## pull in the data
D2010 = read.csv('data/merged_wp_census_data2_081122.csv')
county_adj = read.csv('data/countyadj2.csv', row.names = 1)

## subset data by state
NY_lst <- subset_data_by_state(D2010, county_adj, 'New York', 'NY')

# parameters for simulations and MCMC fitting
models = c('X1','X2','X3')
n.sample = 10000
burnin = 5000

#### 3/18/2024: Normal, different means, direct weights, 1 run ####
# run the simulations
system.time({
  res_lst <- multiple_sims(NY_lst, models, variances = c(10^2, 10^2, 10^2), means = c(100,200,300), N_sims = 1, rho = 0.3, tau2 = 0.01, tau2_fixed = F, family = 'normal', sigma2 = 10^2, direct_weights = T, n.sample = n.sample, burnin = burnin, sigma2_prior_shape = 5, sigma2_prior_rate = 0.05, stan_path = 'code/CAR_leroux_sparse_normal.stan')
}) # 

pp <- process_results(res_lst$sim_list[[1]]$data_list, res_lst$models, res_lst$sim_list[[1]]$stan_fit, tau2_estimates = T, likelihoods = F, sigma2_estimates = T)

# save things
{
  warnings = warnings()
  save(res_lst, warnings, file = sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Results/03182024_normal_3models_difmeans_direct_weights_tau001_sigma5005.RData', root_dir))
  
  ggsave(pp, filename = sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Figures/03182024_normal_3models_difmeans_direct_weights_tau001_sigma5005(single sim).png', root_dir), height = 10, width = 5)
}

#
#### 3/18/2024: Normal, softmax, one run ####

# run the simulations
system.time({
  res_lst <- multiple_sims(NY_lst, models, variances = c(10^2, 10^2, 10^2), means = c(100,100,100), N_sims = 5, rho = 0.3, tau2 = 1, tau2_fixed = F, family = 'normal', sigma2 = 10^2, direct_weights = F, n.sample = n.sample, burnin = burnin, sigma2_prior_shape = 5, sigma2_prior_rate = 0.05, stan_path = 'code/CAR_leroux_sparse_normal.stan')
}) # ~3 minutes

pp <- process_results(res_lst$sim_list[[1]]$data_list, res_lst$models, res_lst$sim_list[[1]]$stan_fit, tau2_estimates = T, likelihoods = F, sigma2_estimates = T)

# save things
{
  warnings = warnings()
  save(res_lst, warnings, file = sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Results/03182024_normal_3models_mean100_softmax_tau1_sigma5005.RData', root_dir))
  
  ggsave(pp, filename = sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Figures/03182024_normal_3models_mean100_softmax_tau1_sigma5005(single sim).png', root_dir), height = 10, width = 5)
}


#
#### 3/18/2024: Normal, direct weights, 1 run ####
# run the simulations
system.time({
  res_lst <- multiple_sims(NY_lst, models, variances = c(10^2, 10^2, 10^2), means = c(100,100,100), N_sims = 1, rho = 0.3, tau2 = 0.01, tau2_fixed = F, family = 'normal', sigma2 = 10^2, direct_weights = T, n.sample = n.sample, burnin = burnin, sigma2_prior_shape = 5, sigma2_prior_rate = 0.05, stan_path = 'code/CAR_leroux_sparse_normal.stan')
}) # ~3 minutes

pp <- process_results(res_lst$sim_list[[1]]$data_list, res_lst$models, res_lst$sim_list[[1]]$stan_fit, tau2_estimates = T, likelihoods = F, sigma2_estimates = T)

# save things
{
  warnings = warnings()
  save(res_lst, warnings, file = sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Results/03182024_normal_3models_mean100_direct_weights_tau001_sigma5005.RData', root_dir))
  
  ggsave(pp, filename = sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Figures/03182024_normal_3models_mean100_direct_weights_tau001_sigma5005(single sim).png', root_dir), height = 10, width = 5)
}

#
#### 3/18/2024: (New) Poisson, direct weights ####
system.time({
  res_lst <- multiple_sims(NY_lst, models, variances = c(10^2, 10^2, 10^2), means = c(100,100,100), N_sims = 1, rho = 0.3, tau2 = 0.01, tau2_fixed = F, family = 'poisson', direct_weights = T, n.sample = n.sample, burnin = burnin, stan_path = 'code/CAR_leroux_sparse_poisson.stan')
}) 

pp <- process_results(res_lst$sim_list[[1]]$data_list, res_lst$models, res_lst$sim_list[[1]]$stan_fit, tau2_estimates = T, likelihoods = F, sigma2_estimates = F)

# save things
{
  warnings = warnings()
  save(res_lst, warnings, file = sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Results/03182024_poisson_3models_mean100_direct_weights_tau001.RData', root_dir))
  
  ggsave(pp, filename = sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Figures/03182024_poisson_3models_mean100_direct_weights_tau001(single sim).png', root_dir), height = 10, width = 5)
}

#

#### 3/18/2024: (New) Poisson, softmax ####
system.time({
  res_lst <- multiple_sims(NY_lst, models, variances = c(10^2, 10^2, 10^2), means = c(100,100,100), N_sims = 1, rho = 0.3, tau2 = 1, tau2_fixed = F, family = 'poisson', direct_weights = F, n.sample = n.sample, burnin = burnin, stan_path = 'code/CAR_leroux_sparse_poisson.stan')
}) 

pp <- process_results(res_lst$sim_list[[1]]$data_list, res_lst$models, res_lst$sim_list[[1]]$stan_fit, tau2_estimates = T, likelihoods = F, sigma2_estimates = F)

# save things
{
  warnings = warnings()
  save(res_lst, warnings, file = sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Results/03182024_poisson_3models_mean100_softmax_tau1.RData', root_dir))
  
  ggsave(pp, filename = sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Figures/03182024_poisson_3models_mean100_softmax_tau1(single sim).png', root_dir), height = 10, width = 5)
}

#
#### Fixing rho ####
system.time({
  res_lst <- multiple_sims(NY_lst, models, variances = c(10^2, 10^2, 10^2), means = c(100,100,100), N_sims = 1, rho = 0.99, tau2 = 0.1, tau2_fixed = F, family = 'normal', sigma2 = 10^2, direct_weights = T, n.sample = n.sample, burnin = burnin, fix_rho_value = 0.99, sigma2_prior_shape = 1000, sigma2_prior_rate = 10, stan_path = 'code/CAR_leroux_sparse_normal.stan')
}) # 4-8 minutes
# with fix_rho == 1 and setting rho_used = 1, I get "Chain 1: Rejecting initial value: Chain 1:   Log probability evaluates to log(0), i.e. negative infinity."
# with fix_rho == 1 and setting rho_used = 0.9, I get no errors. Ok so it has to do with rho = 1. Even with rho = 0.99, it works

panel_plot <- make_panel_plot(res_lst)                                                                
panel_plot

# save things
{
  warnings = warnings()
  save(res_lst, warnings, file = sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Results/03012024_normal_3models_mean100_sigma_prior100010_rho_fixed099.RData', root_dir))
  
  ggsave(panel_plot, filename = sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Figures/03012024_normal_3models_mean100_sigma_prior100010_rho_fixed099.png', root_dir), height = 10, width = 20)
}

#
#### Same mean, stronger sigma priors ####

# run the simulations
system.time({
  res_lst <- multiple_sims(NY_lst, models, variances = c(10^2, 10^2, 10^2), means = c(100,100,100), N_sims = 5, rho = 0.3, tau2 = 1, tau2_fixed = F, family = 'normal', sigma2 = 10^2, direct_weights = F, n.sample = n.sample, burnin = burnin, sigma2_prior_shape = 5, sigma2_prior_rate = 0.05, stan_path = 'code/CAR_leroux_sparse_normal.stan')
}) # ~3 minutes

panel_plot <- make_panel_plot(res_lst)

# save things
{
  warnings = warnings()
  save(res_lst, warnings, file = sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Results/03012024_normal_3models_mean100_sigma_prior505_tau1.RData', root_dir))
  
  ggsave(panel_plot, filename = sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Figures/03012024_normal_3models_mean100_sigma_prior505_tau1.png', root_dir), height = 10, width = 20)
}

### run the simulations Strong prior
system.time({
  res_lst <- multiple_sims(NY_lst, models, variances = c(10^2, 10^2, 10^2), means = c(100,100,100), N_sims = 5, rho = 0.3, tau2 = 1, tau2_fixed = F, family = 'normal', sigma2 = 10^2, direct_weights = F, n.sample = n.sample, burnin = burnin, sigma2_prior_shape = 1000, sigma2_prior_rate = 10, stan_path = 'code/CAR_leroux_sparse_normal.stan')
}) # ~3 minutes

panel_plot <- make_panel_plot(res_lst)

# save things
{
  warnings = warnings()
  save(res_lst, warnings, file = sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Results/03012024_normal_3models_mean100_sigma_prior100010_tau1.RData', root_dir))
  
  ggsave(panel_plot, filename = sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Figures/03012024_normal_3models_mean100_sigma_prior100010_tau1.png', root_dir), height = 10, width = 20)
}

# 
#### No softmax/direct weights simulation - same mean, stronger sigma priors ####

# Gamma(0.001, 0.001)
system.time({
  res_lst <- multiple_sims(NY_lst, models, variances = c(10^2, 10^2, 10^2), means = c(100,100,100), N_sims = 5, rho = 0.3, tau2 = 1, tau2_fixed = F, family = 'normal', sigma2 = 10^2, direct_weights = T, n.sample = n.sample, burnin = burnin, sigma2_prior_shape = 0.001, sigma2_prior_rate = 0.001, stan_path = 'code/CAR_leroux_sparse_normal.stan')
}) # ~30 minutes for 10 simulations

panel_plot <- make_panel_plot(res_lst)

# save things
{
  warnings = warnings()
  save(res_lst, warnings, file = sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Results/03012024_normal_3models_mean100_direct_weights_sigma_prior001001_tau1.RData', root_dir))
  
  ggsave(panel_plot, filename = sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Figures/03012024_normal_3models_mean100_direct_weights_sigma_prior001001_tau1.png', root_dir), height = 10, width = 20)
}




# Gamma (5, 0.05)
system.time({
  res_lst <- multiple_sims(NY_lst, models, variances = c(10^2, 10^2, 10^2), means = c(100,100,100), N_sims = 5, rho = 0.3, tau2 = 1, tau2_fixed = F, family = 'normal', sigma2 = 10^2, direct_weights = T, n.sample = n.sample, burnin = burnin, sigma2_prior_shape = 5, sigma2_prior_rate = 0.05, fix_rho_value = -1, stan_path = 'code/CAR_leroux_sparse_normal.stan')
}) # ~3 minutes

panel_plot <- make_panel_plot(res_lst)

# save things
{
  warnings = warnings()
  save(res_lst, warnings, file = sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Results/02262024_normal_3models_mean100_direct_weights_sigma_prior505_tau1.RData', root_dir))
  
  ggsave(panel_plot, filename = sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Figures/02262024_normal_3models_mean100_direct_weights_sigma_prior505_tau1.png', root_dir), height = 10, width = 20)
}

# Gamma(1000, 10)
system.time({
  res_lst <- multiple_sims(NY_lst, models, variances = c(10^2, 10^2, 10^2), means = c(100,100,100), N_sims = 5, rho = 0.3, tau2 = 1, tau2_fixed = F, family = 'normal', sigma2 = 10^2, direct_weights = T, n.sample = n.sample, burnin = burnin, sigma2_prior_shape = 1000, sigma2_prior_rate = 10, stan_path = 'code/CAR_leroux_sparse_normal.stan')
}) # ~3 minutes

panel_plot <- make_panel_plot(res_lst)

# save things
{
  warnings = warnings()
  save(res_lst, warnings, file = sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Results/03012024_normal_3models_mean100_direct_weights_sigma_prior100010_tau1.RData', root_dir))
  
  ggsave(panel_plot, filename = sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Figures/03012024_normal_3models_mean100_direct_weights_sigma_prior100010_tau1.png', root_dir), height = 10, width = 20)
}


# 
#### No softmax/direct weights simulation - same mean ####

# run the simulations
system.time({
  res_lst <- multiple_sims(NY_lst, models, variances = c(10^2, 10^2, 10^2), means = c(100,100,100), N_sims = 5, rho = 0.3, tau2 = 0.01, tau2_fixed = F, family = 'normal', sigma2 = 10^2, direct_weights = T, n.sample = n.sample, burnin = burnin, stan_path = 'code/CAR_leroux_sparse_normal_noSoftmax.stan')
}) # ~23 minutes

panel_plot <- make_panel_plot(res_lst)

ggsave(panel_plot, filename = sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Figures/02202024_normal_3models_mean100_direct_weights.png', root_dir), height = 10, width = 20)

# 
#### No softmax/direct weights simulation - separated means ####

# run the simulations
system.time({
  res_lst <- multiple_sims(NY_lst, models, variances = c(10^2, 10^2, 10^2), means = c(100,200,300), N_sims = 5, rho = 0.3, tau2 = 0.01, tau2_fixed = F, family = 'normal', sigma2 = 10^2, direct_weights = T, n.sample = n.sample, burnin = burnin, stan_path = 'code/CAR_leroux_sparse_normal_noSoftmax.stan')
}) # ~23 minutes

panel_plot <- make_panel_plot(res_lst)

ggsave(panel_plot, filename = sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Figures/02202024_normal_3models_difmeans_direct_weights.png', root_dir), height = 10, width = 20)

# 
#### MVN model simulation - negative correlation ####
# run the simulations
system.time({
  res_lst <- multiple_sims(NY_lst, models, variances = c(10^2, 10^2, 10^2), means = c(100,100,100), N_sims = 5, rho = 0.3, tau2 = 1, tau2_fixed = F, family = 'normal', sigma2 = 10^2, n.sample = n.sample, burnin = burnin, M_MVN_alpha = -50, stan_path = 'code/CAR_leroux_sparse_normal.stan')
})
# 78m, though maybe due to comp sleeping?

# make dah plot!
panel_plot <- make_panel_plot(res_lst)

# save things
{
  warnings = warnings()
  save(res_lst, warnings, file = sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Results/02122024_normal_3models_mean100_MVNneg_5runs_results.RData', root_dir))
  
  ggsave(panel_plot, filename = sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Figures/02122024_normal_3models_mean100_MVNneg_panel_plot.png', root_dir), height = 10, width = 20)
}

#### MVN model simulation - positive correlation ####
# run the simulations
system.time({
  res_lst <- multiple_sims(NY_lst, models, variances = c(10^2, 10^2, 10^2), means = c(100,100,100), N_sims = 5, rho = 0.3, tau2 = 1, tau2_fixed = F, family = 'normal', sigma2 = 10^2, n.sample = n.sample, burnin = burnin, M_MVN_alpha = 50, stan_path = 'code/CAR_leroux_sparse_normal.stan')
})
# 78m, though maybe due to comp sleeping?

# make dah plot!
panel_plot <- make_panel_plot(res_lst)

# save things
{
  warnings = warnings()
  save(res_lst, warnings, file = sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Results/02122024_normal_3models_mean100_MVNpos_5runs_results.RData', root_dir))
  
  ggsave(panel_plot, filename = sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Figures/02122024_normal_3models_mean100_MVNpos_panel_plot.png', root_dir), height = 10, width = 20)
}

#
#### 3 models, CAR + normal variance (BYM) ####

# run the simulations
system.time({
  res_lst <- multiple_sims(NY_lst, models, variances = c(10^2, 10^2, 10^2), means = c(100,100,100), N_sims = 5, rho = 0.3, tau2 = 1, tau2_fixed = F, family = 'normal', sigma2 = 10^2, n.sample = n.sample, burnin = burnin, M_CAR_tau2 = 100, M_CAR_rho = 0.4, M_BYM_variance = T, stan_path = 'code/CAR_leroux_sparse_normal.stan')
})
# 78m, though maybe due to comp sleeping?

# make dah plot!
panel_plot <- make_panel_plot(res_lst)

# save things
{
  warnings = warnings()
  save(res_lst, warnings, file = sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Results/02122024_normal_3models_mean100_BYMvar_5runs_results.RData', root_dir))
  
  ggsave(panel_plot, filename = sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Figures/02122024_normal_3models_mean100_BYMvar_panel_plot.png', root_dir), height = 10, width = 20)
}

#
#### 3 models + CAR variance ####
# run the simulations
system.time({
  res_lst <- multiple_sims(NY_lst, models, variances = c(10^2, 10^2, 10^2), means = c(100,100,100), N_sims = 5, rho = 0.3, tau2 = 1, tau2_fixed = F, family = 'normal', sigma2 = 10^2, n.sample = n.sample, burnin = burnin, M_CAR_tau2 = 100, M_CAR_rho = 0.4, stan_path = 'code/CAR_leroux_sparse_normal.stan')
})
# 78m, though maybe due to comp sleeping?

# make dah plot!
panel_plot <- make_panel_plot(res_lst)

# save things
{
  warnings = warnings()
  save(res_lst, warnings, file = sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Results/02122024_normal_3models_mean100_CARvar_5runs_results.RData', root_dir))
  
  ggsave(panel_plot, filename = sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Figures/02122024_normal_3models_mean100_CARvar_panel_plot.png', root_dir), height = 10, width = 20)
}

#
#### 3 models, same mean ####


# run the simulations
system.time({
  res_lst <- multiple_sims(NY_lst, models, variances = c(10^2, 10^2, 10^2), means = c(100,100,100), N_sims = 5, rho = 0.3, tau2 = 1, tau2_fixed = F, family = 'normal', sigma2 = 10^2, n.sample = n.sample, burnin = burnin, stan_path = 'code/CAR_leroux_sparse_normal.stan')
})
# 78m, though maybe due to comp sleeping?

# make dah plot!
panel_plot <- make_panel_plot(res_lst)

# save things
{
  swarnings = warnings()
  save(res_lst, warnings, file = sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Results/02122024_normal_3models_mean100_5runs_results.RData', root_dir))
  
  ggsave(panel_plot, filename = sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Figures/02122024_normal_3models_mean100_panel_plot.png', root_dir), height = 10, width = 20)
}

#
#### Running with normal distribution, 3 models ####
models = c('M1','M2','M3')

# run the simulations
system.time({
  res_lst <- multiple_sims(NY_lst, models, variances = c(10^2, 10^2, 10^2), means = c(100,200,300), N_sims = 10, rho = 0.3, tau2 = 1, tau2_fixed = F, family = 'normal', sigma2 = 10^2, stan_path = 'code/CAR_leroux_sparse_normal.stan')
})
# 78m, though maybe due to comp sleeping?

# plot a single simulation results
pp <- process_results(res_lst[[1]]$data_list, models, res_lst[[1]]$stan_fit, tau2_estimates = T, likelihoods = T, sigma2_estimates = T)

# plot multiple simulation parameter estimates and u estimates
tt1 <- plot_multiple_sims(res_lst, models, sigma2_estimates = T)
tt2 <- plot_multiple_sims(res_lst, models, u_estimates = T, rho_estimates = F, tau2_estimates = F)

# save things
{
  #save(res_lst, file = sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Results/02052024_normal_3models_10runs_results.RData', root_dir))
  
  ggsave(pp, filename = sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Figures/02052024_normal_3models_1run.png', root_dir), height = 10, width = 5)
  
  ggsave(plot = tt1, filename = sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Figures/02052024_normal_3models_10runs_spatial_params.png', root_dir), width = 10, height = 15)
  
  ggsave(plot = tt2, filename = sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Figures/02052024_normal_3models_10runs_u_estimates.png', root_dir), width = 10, height = 15)
}

#
#### Running with normal distribution, 2 models ####
models = c('M1','M2')

# run the simulations
system.time({
  res_lst <- multiple_sims(NY_lst, models, variances = c(10^2, 10^2), means = c(100,200), N_sims = 10, rho = 0.3, tau2 = 1, tau2_fixed = F, family = 'normal', sigma2 = 10^2, stan_path = 'code/CAR_leroux_sparse_normal.stan')
})
# 18m

# plot a single simulation results
pp <- process_results(res_lst[[1]]$data_list, models, res_lst[[1]]$stan_fit, tau2_estimates = T, likelihoods = T, sigma2_estimates = T)

# plot multiple simulation parameter estimates and u estimates
tt1 <- plot_multiple_sims(res_lst, models, sigma2_estimates = T)
tt2 <- plot_multiple_sims(res_lst, models, u_estimates = T, rho_estimates = F, tau2_estimates = F)

# save things
{
  save(res_lst, file = 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Population Estimation/Results/02052024_normal_2models_10runs_results.RData')
  
  ggsave(pp, file = 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Population Estimation/Figures/02052024_normal_2models_1run.png', height = 10, width = 5)
  
  ggsave(plot = tt1, filename = 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Population Estimation/Figures/02052024_normal_2models_10runs_spatial_params.png', width = 10, height = 15)
  
  ggsave(plot = tt2, filename = 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Population Estimation/Figures/02052024_normal_2models_10runs_u_estimates.png', width = 10, height = 15)
}
#
#### Running with tau2 fixed, 2 models ####
models = c('M1','M2')

# run the simulations
system.time({
  res_lst <- multiple_sims(NY_lst, models, variances = c(10^2, 10^2), means = c(100,200), N_sims = 10, rho = 0.3, tau2_fixed = T, stan_path = 'code/CAR_leroux_sparse_fixtau2.stan')
})
# 4-5 minutes

# plot a single simulation results
pp <- process_results(res_lst[[1]]$data_list, models, res_lst[[1]]$stan_fit, tau2_estimates = F)

# plot overall spatial param estimates
tt <- plot_multiple_sims(res_lst, models, tau2_estimates = F)

# save things
{
  save(res_lst, file = 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Population Estimation/Results/results_23models_10runs_fixtau2_01312024.RData')
  
  ggsave(pp, file = 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Population Estimation/Figures/2models_fixtau2_01312024.png', height = 10, width = 5)
  
  ggsave(plot = tt, filename = 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Population Estimation/Figures/spatial_params_2models_10runs_fixtau2_01312024.png', width = 10, height = 15)
}


#
#### Running with tau2 fixed, 3 models ####
models = c('M1','M2','M3')

# run the simulations
system.time({
  res_lst <- multiple_sims(NY_lst, models, variances = c(10^2, 10^2, 10^2), means = c(100,200,300), N_sims = 10, rho = 0.3, tau2_fixed = T, stan_path = 'code/CAR_leroux_sparse_fixtau2.stan')
})
# 7-8 minutes

# plot a single simulation results
pp <- process_results(res_lst[[1]]$data_list, models, res_lst[[1]]$stan_fit, tau2_estimates = F)

# plot overall spatial param estimates
tt <- plot_multiple_sims(res_lst, models, tau2_estimates = F)

# save things
{
  save(res_lst, file = 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Population Estimation/Results/results_3models_10runs_fixtau2_01312024.RData')
  
  ggsave(pp, file = 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Population Estimation/Figures/3models_fixtau2_01312024.png', height = 10, width = 5)
  
  ggsave(plot = tt, filename = 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Population Estimation/Figures/spatial_params_3models_10runs_fixtau2_01312024.png', width = 10, height = 15)
}

#
#### Run with high n.sample ####
models = c('M1','M2','M3')
res_lst <- multiple_sims(NY_lst, models, variances = c(10^2, 10^2, 10^2), means = c(100,200, 300), N_sims = 1, rho = 0.3, n.sample = 50000, burnin = 10000)
# warnings? Some divergent transitions, some BFMI low.
# but no ESS issues, because I sampled enough.

#
#### Running with 3 models and rho = 0.7 ####
models = c('M1','M2','M3')

# run the simulations
system.time({
  res_lst <- multiple_sims(NY_lst, models, variances = c(10^2, 10^2, 10^2), means = c(100,200, 300), N_sims = 10, rho = 0.7)
})
# 7 minutes

# plot a single simulation results
pp <- process_results(res_lst[[1]]$data_list, models, res_lst[[1]]$stan_fit)

# plot overall spatial param estimates
tt <- plot_multiple_sims(res_lst, models)

# save things
{
  save(res_lst, file = 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Population Estimation/Results/results_3models_10runs_rho07_01312024.RData')
  
  ggsave(pp, file = 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Population Estimation/Figures/M1100_M2200_M3300_rho07_01312024.png', height = 10, width = 5)
  
  ggsave(plot = tt, filename = 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Population Estimation/Figures/spatial_params_3models_10runs_rho07_01312024.png', width = 10, height = 15)
}

#
#### Running with 2 models and rho = 0.7 ####
models = c('M1','M2')

# run the simulations
system.time({
  res_lst <- multiple_sims(NY_lst, models, variances = c(10^2, 10^2), means = c(100,200), N_sims = 10, rho = 0.7)
})

# plot a single simulation results
pp <- process_results(res_lst[[1]]$data_list, models, res_lst[[1]]$stan_fit)

# plot overall spatial param estimates
tt <- plot_multiple_sims(res_lst, models)

# save things
{
  save(res_lst, file = 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Population Estimation/Results/results_2models_10runs_rho07_01312024.RData')
  
  ggsave(pp, file = 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Population Estimation/Figures/M1100_M2200_rho07_01312024.png', height = 10, width = 5)
  
  ggsave(plot = tt, filename = 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Population Estimation/Figures/spatial_params_2models_10runs_rho07_01312024.png', width = 10, height = 15)
}

#
#### Running with 2 and 3 models ####

## pull in the data
D2010 = read.csv('data/merged_wp_census_data2_081122.csv')
county_adj = read.csv('data/countyadj2.csv', row.names = 1)
models = c('M1','M2')

## subset data by state
NY_lst <- subset_data_by_state(D2010, county_adj, 'New York', 'NY')

## Running with simulated model and y values - 2 models
data <- simulate_models(data = NY_lst$data, models = c('M1','M2'), means = c(100, 200), variances = c(10^2, 10^2))

data_lst <- simulate_y(data, NY_lst$adjacency, models = models, precision_type = 'Leroux', tau2 = 1, rho = 0.3)

stan_fit <- run_stan_CAR(data_lst$data, data_lst$adjacency, models = models, n.sample = 10000, burnin = 5000)

# save(stan_fit, data_lst, file = 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Population Estimation/Results/results_sim_01232024.RData')

pp <- process_results(data_lst, models, stan_fit)

# ggsave(pp, file = 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Population Estimation/Figures/ACS100_PEP200_01242024.png', height = 10, width = 5)


## Running with simulated model and y values - 3 models
models = c('M1','M2','M3')
data <- simulate_models(data = NY_lst$data, models = models, means = c(100, 200, 300), variances = c(10^2, 10^2, 10^2))

data_lst <- simulate_y(data, NY_lst$adjacency, models = models, precision_type = 'Leroux', tau2 = 1, rho = 0.3)

stan_fit <- run_stan_CAR(data_lst$data, data_lst$adjacency, models = models, n.sample = 10000, burnin = 5000)

pp <- process_results(data_lst, models, stan_fit)

# ggsave(pp, file = 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Population Estimation/Figures/ACS100_PEP200_WP300_01242024.png', height = 10, width = 7)

#
#### Running repeatedly - 2 models ####
models = c('M1','M2')

res_lst <- multiple_sims(NY_lst, models, variances = c(10^2, 10^2), means = c(100,200), N_sims = 10)

tt <- plot_multiple_sims(res_lst, models)

# save(res_lst, file = 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Population Estimation/Results/results_2models_10runs_01252024.RData')

# ggsave(plot = tt, filename = 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Population Estimation/Figures/spatial_params_2models_10runs_01252024.png', width = 10, height = 15)

# 
#### Running repeatedly - 3 models ####
models = c('M1', 'M2', 'M3')

res_lst <- multiple_sims(NY_lst, models = models, variances = c(10^2, 10^2, 10^2), means = c(100,200, 300), N_sims = 3)

tt <- plot_multiple_sims(res_lst, models, ncol = 1)

# save(res_lst, file = 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Population Estimation/Results/results_3models_10runs_01252024.RData')

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
data_lst <- simulate_y(data, NY_lst$adjacency, models = models, precision_type = 'Leroux', tau2 = c(1,10), rho = 0, seed = 18)

head(data_lst$phi_true)
var(data_lst$phi_true$M1)
var(data_lst$phi_true$M2)
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
data_lst1 <- simulate_y(data, NY_lst$adjacency, models = models, precision_type = 'Leroux', tau2 = 1, rho = 0.3, seed = seed)

stan_fit1 <- run_stan_CAR(data_lst1$data, data_lst1$adjacency, models = models, n.sample = 10000, burnin = 5000)

pp1 <- process_results(data_lst1, models, stan_fit1)

# MVRnorm
data_lst2 <- simulate_y(data, NY_lst$adjacency, models = models, precision_type = 'Leroux', tau2 = 1, rho = 0.3, cholesky = F, seed = seed)

stan_fit2 <- run_stan_CAR(data_lst2$data, data_lst2$adjacency, models = models, n.sample = 10000, burnin = 5000)

pp2 <- process_results(data_lst2, models, stan_fit2)

tt = plot_grid(pp1, pp2, nrow = 1)
tt

#
