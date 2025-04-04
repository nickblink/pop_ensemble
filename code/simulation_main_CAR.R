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

inputs = c('R=2:dataset=NY:N_models=3:n.sample=100:burnin=50:family=negbin:use_softmax=F:variances=100,100,100:means=100,100,100:rho=0.99:fixed_rho=-1:tau2=0.01:fixed_tau2=-1:tau2_prior_shape=1:tau2_prior_rate=1:theta_true=100:num_y_samples=3:stan_path=code/CAR_leroux_sparse_negbin_test.stan:CV_blocks=-1:return_quantiles=F:parallel=F:output_path=simulation_200GBmem_TEST_rho_1_theta_100_direct_est_negbin_3models_CV5_ID32183_2024_08_16','3')
#inputs = c('R=40:dataset=NY:N_models=3:n.sample=2000:burnin=1000:family=negbin:use_softmax=T:variances=100,100,100:means=100,100,100:rho=0.99:fixed_rho=-1:tau2=1:fixed_tau2=-1:tau2_prior_shape=1:tau2_prior_rate=1:theta_true=100:num_y_samples=3:stan_path=code/CAR_leroux_sparse_negbin_test.stan:CV_blocks=-1:return_quantiles=F:parallel=T:output_path=simulation_rho_099_theta_100_softmax_invchi_prior_2025_04_01','1')

# cluster inputs
inputs <- commandArgs(trailingOnly = TRUE)

# create params list from inputs
{
params <- list()
inputs[[1]] <- gsub('\r', '', inputs[[1]])
params[['job_id']] <- as.integer(inputs[[2]])
for(str in strsplit(inputs[[1]],':')[[1]]){
  tmp = strsplit(str, '=')[[1]]
  nn = tmp[1]
  if(nn %in% c('stan_path', 'output_path')){
    val = tmp[2]
  }else{
    val = tolower(tmp[2])
  }
  
  if(nn %in% c('R', 'n.sample', 'burnin', 'rho', 'fixed_rho', 'tau2', 'fixed_tau2', 'sigma2', 'sigma2_prior_shape', 'sigma2_prior_rate', 'theta_true', 'tau2_prior_shape', 'tau2_prior_rate', 'num_y_samples', 'CV_blocks','theta_prior_shape','theta_prior_rate')){
    val = as.numeric(val)
  }else if(nn == 'N_models'){
    models <- paste0('X', 1:as.numeric(val))
    params[['models']] <- models
    next
  }else if(nn %in% c('means', 'variances')){
    val <- as.numeric(strsplit(val, ',')[[1]])
  }else if(val %in% c('t','f')){
    val = as.logical(ifelse(val == 't', T, F))
  }
  params[[nn]] = val
}

params[['seed_start']] = (params[['job_id']] - 1)*params[['R']]

print(params)

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
}

# set up the output folder
{
  date <- gsub('-','_', Sys.Date())
  if(!('output_path' %in% names(params))){
    params[['output_path']] <- sprintf('results/simulation_%s', date)
  }else{
    if(!grepl('results', params[['output_path']])){
      params[['output_path']] <- sprintf('results/%s', params[['output_path']])
    }
  }
  
  if(!file.exists(params[['output_path']])){
    dir.create(params[['output_path']], recursive = T)
  }
  results_file <- sprintf('%s/sim_results_%i.RData', params[['output_path']],  params[['job_id']])
  
  # if there is already a results file, change name to not over-write it.
  if(file.exists(results_file)){
    iter = 0
    while(file.exists(results_file)){
      iter = iter + 1
      results_file <- sprintf('%s/sim_results_%i(%i).RData', params[['output_path']], params[['job_id']], iter)
    }
  }
}

print('set up the output folder')

# run the model
if(params[['parallel']]){
  # register the cores
  registerDoParallel(cores = params[['R']])
  
  # set number of sims within one node to be 1.
  params[['N_sims']] <- 1
  
  # one run of the code.
  one_run <- function(params, i){
    library(dplyr)
    library(rstan)
    # update the seed start.
    params[['seed_start']] <- params[['seed_start']] + i
    
    # run model and return!
    tmp <- do.call(multiple_sims, params)
    return(tmp)
  }

  print('calling one_run')
  system.time({
    res_lst <- foreach(i = 1:params[['R']]) %dorng% {
      tryCatch({
        message(paste("Starting iteration", i))
        one_run(params, i)
      }, error = function(e) {
        message(paste("Error in iteration", i, ":", e$message))
        message(capture.output(traceback(), type = "message"))
        NULL  # Or use NA, or a named list with error info
      })
    }
  })
  
}else{
  params[['N_sims']] <- params[['R']]
  system.time({
    res_lst <- do.call(multiple_sims, params)
  })
}

# save the results
save(res_lst, params, file = results_file)

if(F){
  tt <- res_lst$sim_list[[1]]
  sf <- tt$stan_fit
  head(sf[,1:5])
}
