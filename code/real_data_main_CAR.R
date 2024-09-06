library(dplyr)
library(rstan)
library(ggplot2)
library(doRNG)
library(doParallel)

rstan_options(auto_write = F)

# set working directory for home or laptop
if(file.exists('C:/Users/Admin-Dell')){
  root_dir = 'C:/Users/Admin-Dell'
  setwd(sprintf('%s/Documents/github_projects/pop_ensemble/', root_dir))
}else if(file.exists('C:/Users/nickl')){
  root_dir = 'C:/Users/nickl'
  setwd(sprintf('%s/Documents/github_projects/pop_ensemble/', root_dir))
}

# load extra functions
source('code/extra_functions_CAR.R')

inputs = c('dataset=MA:n.sample=100:burnin=50:family=negbin:use_softmax=T:models=acs,pep,wp:outcome=census:fixed_rho=-1:fixed_tau2=-1:sigma2_prior_shape=0.001:sigma2_prior_rate=0.001:tau2_prior_shape=1:tau2_prior_rate=1:theta_prior_shape=0.001:theta_prior_rate=0.001:stan_path=code/CAR_leroux_sparse_negbin_alpha.stan:CV_blocks=5:return_quantiles=F:parallel=F:alpha_variance_prior=-1:output_path_addition=test')

# cluster inputs
inputs <- commandArgs(trailingOnly = TRUE)

print(inputs)

# create params list from inputs
{
params <- list()
inputs <- gsub('\r', '', inputs)
for(str in strsplit(inputs,':')[[1]]){
  tmp = strsplit(str, '=')[[1]]
  nn = tmp[1]
  if(nn %in% c('stan_path', 'output_path')){
    val = tmp[2]
  }else{
    val = tolower(tmp[2])
  }
  
  if(nn %in% c('n.sample', 'burnin', 'fixed_rho',  'fixed_tau2', 'sigma2_prior_shape', 'sigma2_prior_rate', 'tau2_prior_shape', 'tau2_prior_rate', 'CV_blocks','theta_prior_shape','theta_prior_rate', 'alpha_variance_prior')){
    val = as.numeric(val)
  }else if(nn == 'N_models'){
    models <- paste0('X', 1:as.numeric(val))
    params[['models']] <- models
    next
  }else if(val %in% c('t','f')){
    val = as.logical(ifelse(val == 't', T, F))
  }else if(nn == 'models'){
    val = as.character(strsplit(val, ',')[[1]])
  }
  params[[nn]] = val
}

# make results file
{
  date <- gsub('-','_', Sys.Date())
  
  
  if(!is.null(params$output_path_addition)){
    results_file <- sprintf('results/real_data_results/real_data_fit_%s_ID%s_%s.RData', params$output_path_addition, sample(100000, 1), date)
  }else{
    results_file <- sprintf('results/real_data_results/real_data_fit_%s_%s.RData', sample(100000, 1), date)
  }
  
  # if there is already a results file, change name to not over-write it.
  if(file.exists(results_file)){
    iter = 0
    while(file.exists(results_file)){
      iter = iter + 1
      results_file <- sprintf('results/real_data_results/real_data_fit_%s(%i).RData', date, iter)
    }
  }
  params$results_file <- results_file
}

print(params)

# pull in the data
load('data/census_ACS_PEP_WP_cleaned_08292024.RData')

# subset data by state
if(params[['dataset']] == 'all'){
  # params[['raw_data']]
  # NY_lst <- subset_data_by_state(df, adjacency, 'New York', 'NY')
  
  # store the data in params
  params[['raw_data']] <- list(data = df, adjacency = adjacency)
}else{
  params[['raw_data']] <- subset_data_by_state(df, adjacency, abbrev = params[['dataset']])
}
}

# fit the real data!
res <- do.call(fit_model_real, params)
# res <- fit_model_real(params)

# save the results
save(res, params, file = params$results_file)

