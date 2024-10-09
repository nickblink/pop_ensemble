library(dplyr)

## This script makes bash commands for given simulations
bash_command <- function(R=40, dataset='NY', N_models=3, n.sample=10000, burnin=5000, family='normal',use_softmax=F,variances=c(100,100,100), means=c(100,100,100), rho = 0.3,  fixed_rho = -1, tau2 = .01, fixed_tau2 = -1, sigma2 = 100, sigma2_prior_shape = 50, sigma2_prior_rate = 0.5, tau2_prior_shape = 1, tau2_prior_rate=1, num_y_samples=3, theta = 10, theta_prior_shape = 0.001, theta_prior_rate = 0.001, stan_path='code/CAR_leroux_sparse_normal.stan', CV_blocks = 5, return_quantiles = T,parallel = T, output_path_addition = NULL, array_length = 5, alpha_variance_prior=-1, chains_cores=1){
  
  if(!grepl(family, stan_path)){
    warning('family not in stan path - is that correct?')
  }
  
  if(!use_softmax & tau2 > 0.1){
    print('Are you sure you want tau2 so high?')
  }
  
  # make the output path
  if(is.null(output_path_addition)){
    output_path <- sprintf('simulation_%s_%s_%smodels_CV%s_ID%s_%s', ifelse(use_softmax, 'softmax', 'direct_est'), family, N_models, CV_blocks, sample(1e6, size = 1), gsub('-','_',Sys.Date()))
    output_path <- gsub('\\.','',output_path)
  }else{
    output_path <- sprintf('simulation_%s_%s_%s_%smodels_CV%s_ID%s_%s', output_path_addition, ifelse(use_softmax, 'softmax', 'direct_est'), family, N_models, CV_blocks, sample(1e6, size = 1), gsub('-','_',Sys.Date()))
    output_path <- gsub('\\.','',output_path)
  }
  
  job_name = sprintf('%s_%s_%s_%smodels_CV%s', round(runif(1)*1000), ifelse(use_softmax, 'softmax', 'directest'), family, N_models, CV_blocks)
  
  params <- list(R=R, dataset=dataset, N_models=N_models, n.sample=n.sample, burnin=burnin, family=family,use_softmax=use_softmax,variances=variances, means=means, rho = rho, fixed_rho = fixed_rho, tau2 = tau2, fixed_tau2 = fixed_tau2, sigma2 = sigma2, sigma2_prior_shape = sigma2_prior_shape, sigma2_prior_rate = sigma2_prior_rate, tau2_prior_shape = tau2_prior_shape, tau2_prior_rate=tau2_prior_rate, theta = theta, theta_prior_shape = theta_prior_shape, theta_prior_rate = theta_prior_rate, num_y_samples=num_y_samples, stan_path=stan_path, CV_blocks = CV_blocks, return_quantiles = return_quantiles, parallel = parallel, output_path = output_path, chains_cores = chains_cores, alpha_variance_prior = alpha_variance_prior)
  
  param_str = paste(paste(names(params), params, sep = '='), collapse=':') %>%
    gsub(' |c\\(|\\)','',.) %>%
    gsub('TRUE','T',.) %>%
    gsub('FALSE','F',.)

  if(array_length > 1){
    command_str = sprintf('sbatch --array=1-%s -J %s run_sim_pop_est.bash %s', array_length, job_name, param_str)
  }else{
    command_str = sprintf('sbatch --array=1-%s -J %s run_sim_pop_est.bash %s', array_length, job_name, param_str)
  }
  
  
  return(command_str)

}

## This script makes bash commands for given simulations
bash_command_real <- function(dataset='all', models='acs,pep,wp', n.sample=2000, burnin=1000, outcome = 'census', family='negbin', use_softmax=F, fixed_rho = -1, fixed_tau2 = -1, sigma2_prior_shape = 50, sigma2_prior_rate = 0.5, tau2_prior_shape = 1, tau2_prior_rate=1, theta_prior_shape = 0.001, theta_prior_rate = 0.001, stan_path='code/CAR_leroux_sparse_negbin_alpha.stan', CV_blocks = 5, return_quantiles = F, output_path_addition = NULL, alpha_variance_prior=-1, chains_cores=10, preprocess_scale = F){
  
  if(!grepl(family, stan_path)){
    warning('family not in stan path - is that correct?')
  }
  
  job_name = sprintf('%s_%s_%smodels_CV%s', round(runif(1)*1000), ifelse(use_softmax, 'softmax', 'directest'), family, CV_blocks)
  
  params <- list(dataset=dataset, models=models, n.sample=n.sample, burnin=burnin, outcome = outcome, family=family,use_softmax=use_softmax, fixed_rho = fixed_rho, fixed_tau2 = fixed_tau2,  sigma2_prior_shape = sigma2_prior_shape, sigma2_prior_rate = sigma2_prior_rate, tau2_prior_shape = tau2_prior_shape, tau2_prior_rate=tau2_prior_rate,  theta_prior_shape = theta_prior_shape, theta_prior_rate = theta_prior_rate, stan_path=stan_path, CV_blocks = CV_blocks, return_quantiles = return_quantiles, output_path_addition = output_path_addition, chains_cores = chains_cores, alpha_variance_prior = alpha_variance_prior, preprocess_scale = preprocess_scale)
  
  param_str = paste(paste(names(params), params, sep = '='), collapse=':') %>%
    gsub(' |c\\(|\\)','',.) %>%
    gsub('TRUE','T',.) %>%
    gsub('FALSE','F',.)

    command_str = sprintf('sbatch -J %s run_real_data_pop_est.bash %s', job_name, param_str)
  
  return(command_str)
}

bash_wrapper <- function(bash_file = NULL, theta_vec = NULL, rho_vec = NULL, output_path_addition = NULL, ...){
  
  # get the bash commands
  if(!is.null(theta_vec) | !is.null(rho_vec)){
    params <- expand.grid(if(is.null(theta_vec)) theta else theta_vec,
                          if(is.null(rho_vec)) rho else rho_vec)

    if(is.null(output_path_addition)){
      cmds <- lapply(1:nrow(params), function(ii) {
        bash_command(theta = params[ii,1], 
                     rho = params[ii,2], 
                     output_path_addition = sprintf('rho_%s_theta_%s',as.character(params[ii,2]), as.character(params[ii,1])), ...)})
    }else{
      cmds <- lapply(1:nrow(params), function(ii) {
        bash_command(theta = params[ii,1], 
                     rho = params[ii,2], 
                     output_path_addition = sprintf('%s_rho_%s_theta_%s', output_path_addition, as.character(params[ii,2]), as.character(params[ii,1])), ...)})
    }
  }
  
  # write the commands
  if(!is.null(bash_file)){
    # if it already exists, update it
    if(file.exists(bash_file)){
      lapply(cmds, write, bash_file, append = T, sep = '', )
      close(bash_file)
    # if it doesn't exist, create it
    }else{
      out_file <- file(bash_file, open='wb')
      lapply(cmds, write, out_file, append = T, sep = '')
      close(out_file)
    }
    
    # check that there are no repeats in commands
    test <- read.table(bash_file)
    col <- lapply(test[,ncol(test)], function(str){
      tmp <- strsplit(str, ':')[[1]]
      tmp <- tmp[-grep('output_path', tmp)]
      paste(tmp, collapse = ':')
    })
    
    if(length(unique(col)) != length(col)){
      stop('there are repeating simulation commands')
    }
  }
  
  return(cmds)
}

#### Real bash - NY with higher MCMC ####
bash_command_real(use_softmax = T, preprocess_scale = T, alpha_variance_prior = 0.01,  output_path_addition = 'softmax_preprocess_alpha', dataset = 'NY', n.sample = 20000, burnin = 5000)

#### Real bash commands ####
bash_command_real(use_softmax = F)
bash_command_real(use_softmax = T, output_path_addition = 'softmax_parallel')
bash_command_real(use_softmax = F, preprocess_scale = T,  output_path_addition = 'directest_preprocess')
bash_command_real(use_softmax = T, preprocess_scale = T,  output_path_addition = 'softmax_preprocess')
bash_command_real(use_softmax = T, alpha_variance_prior = 0.01,  output_path_addition = 'softmax_alpha')
bash_command_real(use_softmax = T, preprocess_scale = T, alpha_variance_prior = 0.01,  output_path_addition = 'softmax_preprocess_alpha')


# 
#### Testing NB ####
bash_command(CV_blocks = 10, tau2 = 0.01, rho = 0.3, theta = 10, output_path_addition = 'theta10_rho03')

bash_wrapper(theta_vec = c(10,100), rho_vec = c(0.3, 0.99), family = 'negbin', stan_path = 'code/CAR_leroux_sparse_negbin.stan', bash_file = 'code/bash_commands/negbin_bash_08232024.txt')

bash_wrapper(theta_vec = c(10,100), rho_vec = c(0.3, 0.99), tau2 = 1, use_softmax = T, family = 'negbin', stan_path = 'code/CAR_leroux_sparse_negbin.stan', bash_file = 'code/bash_commands/negbin_bash_08232024.txt')

#### Testing a different sigma2 ####
bash_command(CV_blocks = 10, tau2 = 0.01, rho = 0.3, sigma2 = 25, sigma2_prior_shape = 1, sigma2_prior_rate = .01, output_path_addition = 'sigma2eq25_rho03')

bash_command(CV_blocks = 10, tau2 = 1, use_softmax = T, rho = 0.3, sigma2 = 25, sigma2_prior_shape = 1, sigma2_prior_rate = .01, output_path_addition = 'sigma2eq25_rho03')

bash_command(CV_blocks = 10, tau2 = 0.01, rho = 0.99, sigma2 = 25, sigma2_prior_shape = 1, sigma2_prior_rate = .01, output_path_addition = 'sigma2eq25_rho099')

bash_command(CV_blocks = 10, tau2 = 1, use_softmax = T, rho = 0.99, sigma2 = 25, sigma2_prior_shape = 1, sigma2_prior_rate = .01, output_path_addition = 'sigma2eq25_rho099')

bash_command(CV_blocks = 10, tau2 = 0.01, rho = 0.3, sigma2 = 4, sigma2_prior_shape = 1, sigma2_prior_rate = .01, output_path_addition = 'sigma2eq4_rho03')

bash_command(CV_blocks = 10, tau2 = 1, use_softmax = T, rho = 0.3, sigma2 = 4, sigma2_prior_shape = 1, sigma2_prior_rate = .01, output_path_addition = 'sigma2eq4_rho03')

bash_command(CV_blocks = 10, tau2 = 0.01, rho = 0.99, sigma2 = 4, sigma2_prior_shape = 1, sigma2_prior_rate = .01, output_path_addition = 'sigma2eq4_rho099')

bash_command(CV_blocks = 10, tau2 = 1, use_softmax = T, rho = 0.99, sigma2 = 4, sigma2_prior_shape = 1, sigma2_prior_rate = .01, output_path_addition = 'sigma2eq4_rho099')

#### Testing fixed rho ####
bash_command(CV_blocks = 10, tau2 = 0.01, rho = 0.3, fixed_rho = 0.99, output_path_addition = 'fixedrho099_rho03')

bash_command(CV_blocks = 10, tau2 = 0.01, rho = 0.7, fixed_rho = 0.99, output_path_addition = 'fixedrho099_rho07')

bash_command(CV_blocks = 10, tau2 = 0.01, rho = 0.99, fixed_rho = 0.99, output_path_addition = 'fixedrho099_rho099')

bash_command(CV_blocks = 10, use_softmax = T, tau2 = 1, rho = 0.3, fixed_rho = 0.99, output_path_addition = 'fixedrho099_rho03')

bash_command(CV_blocks = 10, use_softmax = T, tau2 = 1,  rho = 0.7, fixed_rho = 0.99, output_path_addition = 'fixedrho099_rho07')

bash_command(CV_blocks = 10, use_softmax = T, tau2 = 1, rho = 0.99, fixed_rho = 0.99, output_path_addition = 'fixedrho099_rho099')

#####
