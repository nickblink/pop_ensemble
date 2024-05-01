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

# Take in the inputs
inputs = c('N_sims=1:dataset=NY:N_models=3:n.sample=1000:burnin=500:family=normal:use_softmax=T:variances=c(100,100,100):means=c(100,100,100):rho=0.3:tau2=1:sigma2=100:sigma2_prior_shape=50:sigma2_prior_rate=0.5:tau2_prior_shape=1:tau2_prior_rate=1:num_y_samples=3:stan_path=code/CAR_leroux_sparse_normal.stan:CV_blocks=5\r','3')
inputs <- commandArgs(trailingOnly = TRUE)
params <- list()

inputs[[1]] <- gsub('\r', '', inputs[[1]])
params[['job_id']] <- as.integer(inputs[[2]])
for(str in strsplit(inputs[[1]],':')[[1]]){
  tmp = strsplit(str, '=')[[1]]
  print(tmp)
  nn = tmp[1]
  val = tolower(tmp[2])
  if(nn %in% c('N_sims', 'n.sample', 'burnin', 'rho', 'tau2', 'sigma2', 'sigma2_prior_shape', 'sigma2_prior_rate', 'tau2_prior_shape', 'tau2_prior_rate', 'num_y_samples', 'CV_blocks')){
    val = as.numeric(val)
  }else if(nn == 'N_models'){
    models <- paste0('X', 1:as.numeric(val))
    params[['models']] <- models
    next
  }else if(nn %in% c('means', 'variances')){
    val <- as.numeric(strsplit(gsub('c\\(|\\)', '', val), ',')[[1]])
  }else if(val %in% c('t','f')){
    val = as.logical(ifelse(val == 't', T, F))
  }
  params[[nn]] = val
}

params[['seed_start']] = (params[['job_id']] - 1)*params[['N_sims']]

# pull in the data
D2010 = read.csv('data/merged_wp_census_data2_081122.csv')
county_adj = read.csv('data/countyadj2.csv', row.names = 1)

# subset data by state
if(params[['dataset']] == 'ny'){
  NY_lst <- subset_data_by_state(D2010, county_adj, 'New York', 'NY')
  
  # store the data in params
  params[['raw_data']] <- NY_lst
}else{
  stop('only taking NY data right now')
}

# run the code!
system.time({
  res_lst <- do.call(multiple_sims, params)
})


# RUN WITH PARALLEL POWER

# SAVE THE RESULTS
### Need to create the output folder for this run.