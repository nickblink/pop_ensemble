library(dplyr)
library(rstan)
library(ggplot2)
library(doRNG)
library(doParallel)

# set working directory for home or laptop
home_dir = T
if(file.exists('C:/Users/Admin-Dell')){
  root_dir = 'C:/Users/Admin-Dell'
  setwd(sprintf('%s/Documents/github_projects/pop_ensemble/', root_dir))
  rstan_options(auto_write = T)
}else if(file.exists('C:/Users/nickl')){
  root_dir = 'C:/Users/nickl'
  setwd(sprintf('%s/Documents/github_projects/pop_ensemble/', root_dir))
  rstan_options(auto_write = T)
}else{
  home_dir = F
  rstan_options(auto_write = F)
}

# load extra functions
source('code/extra_functions_CAR.R')

if(home_dir){
  inputs = 'dataset=fullsubset:models=acs,pep,wp:n.sample=100:burnin=50:outcome=census:family=negbin:use_softmax=T:fixed_rho=-1:fixed_tau2=-1:sigma2_prior_shape=50:sigma2_prior_rate=0.5:tau2_prior_shape=1:tau2_prior_rate=1:theta_gamma_prior=0:stan_path=code/CAR_leroux_sparse_negbin_alpha_FE_NC.stan:CV_blocks=10:return_quantiles=F:output_path_addition=full_softmax_pepdensity_alpha0001_noncentered:chains_cores=10:alpha_variance_prior=1e-04:preprocess_scale=F:fixed_effects=pep_density:phi_noncentered=1'
}else{
  # cluster inputs
  inputs <- commandArgs(trailingOnly = TRUE)
}

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

  if(val %in% c('NULL', 'null')){
    val <- NULL
  }else if(nn %in% c('n.sample', 'burnin', 'fixed_rho',  'fixed_tau2', 'sigma2_prior_shape', 'sigma2_prior_rate', 'tau2_prior_shape', 'tau2_prior_rate', 'CV_blocks','theta_prior_shape','theta_prior_rate', 'theta_gamma_prior', 'alpha_variance_prior', 'chains_cores','phi_noncentered')){
    val = as.numeric(val)
  }else if(nn == 'N_models'){
    models <- paste0('X', 1:as.numeric(val))
    params[['models']] <- models
    next
  }else if(val %in% c('t','f')){
    val = as.logical(ifelse(val == 't', T, F))
  # }else if(val %in% c('NULL', 'null')){
  #   val <- NULL
  }else if(nn == 'models'){
    val = as.character(strsplit(val, ',')[[1]])
  }
  # else if(val %in% c('NULL', 'null')){
  #   val <- NULL
  # }
  
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

## pull in the data
# American Indian Alaska Native data
if(params[['dataset']] == 'aian'){
  load('data/census_ACS_PEP_WP_AIAN_wDensity_and2017_01022025.RData')
  params[['raw_data']] <- list(data = df, adjacency = adjacency)

  # American Indian Alaska Native subset in the SW
}else if(params[['dataset']] %in% c('aiansub', 'aiansubset','aian2')){
  load('data/census_ACS_PEP_WP_AIANsubset_wDensity_and2018_05212025.RData')
  params[['raw_data']] <- list(data = df, adjacency = adjacency)

  # full data subset in the SW
}else if(params[['dataset']] %in% c('fullpopsub', 'fullsub','fullsubset')){
  load('data/census_ACS_PEP_WP_fullpopsubset_wDensity_and2018_06202025.RData')
  params[['raw_data']] <- list(data = df, adjacency = adjacency)

  # full data elsewhere
}else{
  load('data/census_ACS_PEP_WP_wDensity_and2017_01022025.RData')
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

}

print(colnames(params$raw_data$data))

# fit the real data!
res <- do.call(fit_model_real, params)

if(F){
  tt <- res$sim_list
  tt$stan_summary$summary -> ss
  grep('zeta', rownames(ss)) -> ind
  head(ss[ind,])
  grep('beta', rownames(ss)) -> ind2
  head(ss[ind2,])
}

print('model run')

# save the results
save(res, params, file = params$results_file)

