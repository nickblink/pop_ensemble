library(dplyr)

## This script makes bash commands for given simulations
bash_command <- function(R=40, dataset='NY', N_models=3, n.sample=10000, burnin=5000, family='negbin',use_softmax=T,variances=c(100,100,100), means=c(100,100,100), rho = 0.3,  fixed_rho = -1, tau2 = 1, fixed_tau2 = -1, tau2_prior_shape = 1, tau2_prior_rate=1, theta_true = 10, theta_gamma_prior = 0, num_y_samples=3,  stan_path='code/CAR_leroux_sparse_negbin.stan', CV_blocks = 5, return_quantiles = F,parallel = T, output_path_addition = NULL, array_length = 5, alpha_variance_prior=-1, chains_cores=1, phi_noncentered = NULL){

  if(!grepl(family, stan_path)){
    warning('family not in stan path - is that correct?')
  }
  
  if(!use_softmax & tau2 > 0.1){
    print('Are you sure you want tau2 so high?')
  }
  
  if(use_softmax & alpha_variance_prior > 0.0001){
    print('Are you sure you want the alpha prior variance so high?')
  }
  
  if(!use_softmax & alpha_variance_prior > 0){
    stop('dont use alpha variance prior with direct estimate.')
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
  
  params <- list(R=R, dataset=dataset, N_models=N_models, n.sample=n.sample, burnin=burnin, family=family,use_softmax=use_softmax,variances=variances, means=means, rho = rho, fixed_rho = fixed_rho, tau2 = tau2, fixed_tau2 = fixed_tau2, tau2_prior_shape = tau2_prior_shape, tau2_prior_rate=tau2_prior_rate, theta_true = theta_true, theta_gamma_prior = theta_gamma_prior, num_y_samples=num_y_samples, stan_path=stan_path, CV_blocks = CV_blocks, return_quantiles = return_quantiles, parallel = parallel, output_path = output_path, chains_cores = chains_cores, alpha_variance_prior = alpha_variance_prior, phi_noncentered = phi_noncentered)
  
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
bash_command_real <- function(dataset='all', models='acs,pep,wp', n.sample=2000, burnin=1000, outcome = 'census', family='negbin', use_softmax=F, fixed_rho = -1, fixed_tau2 = -1, sigma2_prior_shape = 50, sigma2_prior_rate = 0.5, tau2_prior_shape = 1, tau2_prior_rate = 1, theta_gamma_prior = 0, stan_path='code/CAR_leroux_sparse_negbin_alpha_FE.stan', CV_blocks = 5, return_quantiles = F, output_path_addition = NULL, alpha_variance_prior=-1, chains_cores=10, preprocess_scale = F, fixed_effects = NULL, phi_noncentered = NULL){
  
  ## Parameter error checking.
  if(!grepl(family, stan_path)){
    warning('family not in stan path - is that correct?')
  }
  
  if(!is.null(fixed_effects)){
    if(!(fixed_effects %in% c('intercept', 'pep_density', 'acs_density', 'pep_fulldensity', 'acs_fulldensity', 'pep_density_proportion', 'pep_fulldensity_proportion'))){
      stop('Incorrect fixed effects')
    }
  }
  
  if(!(tolower(dataset) %in% c('all','aian'))){
    stop('Incorrect dataset')
  }
  
  # if(!use_softmax & (tau2_prior_shape > 0.1 | tau2_prior_rate > 0.1)){
  #   print('Are you sure you want tau2 so high?')
  # } ## Gamma(1,1) is fine for low values of tau2.
  
  job_name = sprintf('%s_%s_%smodels_CV%s', round(runif(1)*1000), ifelse(use_softmax, 'softmax', 'directest'), family, ifelse(is.null(CV_blocks), 'none', CV_blocks))
  
  params <- list(dataset=dataset, models=models, n.sample=n.sample, burnin=burnin, outcome = outcome, family=family,use_softmax=use_softmax, fixed_rho = fixed_rho, fixed_tau2 = fixed_tau2,  sigma2_prior_shape = sigma2_prior_shape, sigma2_prior_rate = sigma2_prior_rate, tau2_prior_shape = tau2_prior_shape, tau2_prior_rate = tau2_prior_rate, theta_gamma_prior = theta_gamma_prior, stan_path=stan_path, CV_blocks = CV_blocks, return_quantiles = return_quantiles, output_path_addition = output_path_addition, chains_cores = chains_cores, alpha_variance_prior = alpha_variance_prior, preprocess_scale = preprocess_scale, fixed_effects = fixed_effects, phi_noncentered = phi_noncentered)
  
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
    params <- expand.grid(if(is.null(theta_vec)) theta_true else theta_vec,
                          if(is.null(rho_vec)) rho else rho_vec)

    if(is.null(output_path_addition)){
      cmds <- lapply(1:nrow(params), function(ii) {
        bash_command(theta_true = params[ii,1], 
                     rho = params[ii,2], 
                     output_path_addition = sprintf('rho_%s_theta_%s',as.character(params[ii,2]), as.character(params[ii,1])), ...)})
    }else{
      cmds <- lapply(1:nrow(params), function(ii) {
        bash_command(theta_true = params[ii,1], 
                     rho = params[ii,2], 
                     output_path_addition = sprintf('%s_rho_%s_theta_%s', output_path_addition, as.character(params[ii,2]), as.character(params[ii,1])), ...)})
    }
  }else{
    cmds <- bash_command(...)
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

bash_wrapper_real <- function(bash_file = NULL, ...){
  cmds <- bash_command_real(...)
  
  if(!is.null(bash_file)){
    # if it already exists, update it
    if(file.exists(bash_file)){
      lapply(cmds, write, bash_file, append = T, sep = '')
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

#### 4/9/2025: Updated bash commands real - testing noncentering ####
bash_wrapper_real(dataset = 'all', 
                  use_softmax = T, 
                  alpha_variance_prior = .0001, 
                  fixed_effects = 'pep_density', 
                  phi_noncentered = 0, 
                  stan_path = 'code/CAR_leroux_sparse_negbin_alpha_FE_NC.stan',
                  output_path_addition = 'full_softmax_pepdensity_alpha0001_centered', 
                  bash_file = 'code/bash_commands/real_data_noncentering_04092025.txt')

bash_wrapper_real(dataset = 'all', 
                  use_softmax = T, 
                  alpha_variance_prior = .0001, 
                  fixed_effects = 'pep_density', 
                  phi_noncentered = 1, 
                  stan_path = 'code/CAR_leroux_sparse_negbin_alpha_FE_NC.stan',
                  output_path_addition = 'full_softmax_pepdensity_alpha0001_noncentered', 
                  bash_file = 'code/bash_commands/real_data_noncentering_04092025.txt')

bash_wrapper_real(dataset = 'AIAN', 
                  use_softmax = T, 
                  alpha_variance_prior = .0001, 
                  fixed_effects = 'pep_density', 
                  phi_noncentered = 0, 
                  stan_path = 'code/CAR_leroux_sparse_negbin_alpha_FE_NC.stan',
                  output_path_addition = 'full_softmax_pepdensity_alpha0001_centered', 
                  bash_file = 'code/bash_commands/real_data_noncentering_04092025.txt')

bash_wrapper_real(dataset = 'AIAN', 
                  use_softmax = T, 
                  alpha_variance_prior = .0001, 
                  fixed_effects = 'pep_density', 
                  phi_noncentered = 1, 
                  stan_path = 'code/CAR_leroux_sparse_negbin_alpha_FE_NC.stan',
                  output_path_addition = 'full_softmax_pepdensity_alpha0001_noncentered', 
                  bash_file = 'code/bash_commands/real_data_noncentering_04092025.txt')



## Simulation commands
bash_wrapper(use_softmax = T,
             alpha_variance_prior = .0001, 
             theta_vec = 100, 
             rho_vec = c(0.3, 0.99),
             CV_blocks = 10, 
             stan_path = 'code/CAR_leroux_sparse_negbin_NC.stan',
             phi_noncentered = 0, 
             output_path_addition = 'centered',
             bash_file = 'code/bash_commands/simulation_noncentering_04092025.txt')

bash_wrapper(use_softmax = T,
             alpha_variance_prior = .0001, 
             theta_vec = 100, 
             rho_vec = c(0.3, 0.99),
             CV_blocks = 10, 
             stan_path = 'code/CAR_leroux_sparse_negbin_NC.stan',
             phi_noncentered = 1, 
             output_path_addition = 'noncentered', 
             bash_file = 'code/bash_commands/simulation_noncentering_04092025.txt')

#
#### 3/12/2025: Simulation bash again! Creating round 1 sim commands for paper ####
bash_wrapper(use_softmax = T, tau2 = 1, CV_blocks = 10, theta_vec = c(100), rho_vec = c(0.3, 0.99), family = 'negbin', stan_path = 'code/CAR_leroux_sparse_negbin.stan', bash_file = 'code/bash_commands/simulation_bash_03122025.txt')

bash_wrapper(use_softmax = F, tau2 = .01, CV_blocks = 10, theta_vec = c(100), rho_vec = c(0.3, 0.99), family = 'negbin', stan_path = 'code/CAR_leroux_sparse_negbin.stan', bash_file = 'code/bash_commands/simulation_bash_03122025.txt')

#
#### 1/13/2025 2017 and 2018 prediction models ####

## Directest intercept only
{
# (1.1) Full pop,  directest, interceptonly: 2018 only
bash_wrapper_real(dataset = 'all', bash_file = 'code/bash_commands/real_data_laggedyear_01132025.txt', use_softmax = F, fixed_effects = 'intercept', output_path_addition = 'directest_interceptonly_2018only', models = 'acs_2018,pep_2018,wp')

# (1.2) Full pop,  directest, interceptonly: 2018 and 2019
bash_wrapper_real(dataset = 'all', bash_file = 'code/bash_commands/real_data_laggedyear_01132025.txt', use_softmax = F, fixed_effects = 'intercept', output_path_addition = 'directest_interceptonly_20182019', models = 'acs_2018,pep_2018,acs,pep,wp')

# (1.3) Full pop,  directest, interceptonly: 2017 and 2018
bash_wrapper_real(dataset = 'all', bash_file = 'code/bash_commands/real_data_laggedyear_01132025.txt', use_softmax = F, fixed_effects = 'intercept', output_path_addition = 'directest_interceptonly_20172018', models = 'acs_2018,pep_2018,acs_2017,pep_2017 ,wp')

# (1.4) AIAN,  directest, interceptonly: 2018 only
bash_wrapper_real(dataset = 'AIAN', bash_file = 'code/bash_commands/real_data_laggedyear_01132025.txt', use_softmax = F, fixed_effects = 'intercept', output_path_addition = 'aian_directest_interceptonly_2018only', models = 'acs_2018,pep_2018,wp')

# (1.5) AIAN,  directest, interceptonly: 2018 and 2019
bash_wrapper_real(dataset = 'AIAN', bash_file = 'code/bash_commands/real_data_laggedyear_01132025.txt', use_softmax = F, fixed_effects = 'intercept', output_path_addition = 'aian_directest_interceptonly_20182019', models = 'acs_2018,pep_2018,acs,pep,wp')

# (1.6) Full pop,  directest, interceptonly: 2017 and 2018
bash_wrapper_real(dataset = 'AIAN', bash_file = 'code/bash_commands/real_data_laggedyear_01132025.txt', use_softmax = F, fixed_effects = 'intercept', output_path_addition = 'aian_directest_interceptonly_20172018', models = 'acs_2018,pep_2018,acs_2017,pep_2017 ,wp')
}

## Softmax alpha density 
{
  
  # (2.1) Full pop,  softmax_alpha, density: 2018 only
  bash_wrapper_real(dataset = 'all', bash_file = 'code/bash_commands/real_data_laggedyear_01132025.txt', use_softmax = T,  alpha_variance_prior = .01, fixed_effects = 'pep_density', output_path_addition = 'softmax_alpha_density_2018only', models = 'acs_2018,pep_2018,wp')
  
  # (2.2) Full pop,  softmax_alpha, density: 2018 and 2019
  bash_wrapper_real(dataset = 'all', bash_file = 'code/bash_commands/real_data_laggedyear_01132025.txt', use_softmax = T,  alpha_variance_prior = .01, fixed_effects = 'pep_density', output_path_addition = 'softmax_alpha_density_20182019', models = 'acs_2018,pep_2018,acs,pep,wp')
  
  # (2.3) Full pop,  softmax_alpha, density: 2017 and 2018
  bash_wrapper_real(dataset = 'all', bash_file = 'code/bash_commands/real_data_laggedyear_01132025.txt', use_softmax = T,  alpha_variance_prior = .01, fixed_effects = 'pep_density', output_path_addition = 'softmax_alpha_density_20172018', models = 'acs_2018,pep_2018,acs_2017,pep_2017 ,wp')
  
  # (2.4) AIAN,  softmax_alpha, pep_fulldensity: 2018 only
  bash_wrapper_real(dataset = 'AIAN', bash_file = 'code/bash_commands/real_data_laggedyear_01132025.txt', use_softmax = T,  alpha_variance_prior = .01, fixed_effects = 'pep_fulldensity', output_path_addition = 'aian_softmax_alpha_pepfulldensity_2018only', models = 'acs_2018,pep_2018,wp')
  
  # (2.5) AIAN,  softmax_alpha, pep_fulldensity: 2018 and 2019
  bash_wrapper_real(dataset = 'AIAN', bash_file = 'code/bash_commands/real_data_laggedyear_01132025.txt', use_softmax = T,  alpha_variance_prior = .01, fixed_effects = 'pep_fulldensity', output_path_addition = 'aian_softmax_alpha_pepfulldensity_20182019', models = 'acs_2018,pep_2018,acs,pep,wp')
  
  # (2.6) Full pop,  softmax_alpha, pep_fulldensity: 2017 and 2018
  bash_wrapper_real(dataset = 'AIAN', bash_file = 'code/bash_commands/real_data_laggedyear_01132025.txt', use_softmax = T,  alpha_variance_prior = .01, fixed_effects = 'pep_fulldensity', output_path_addition = 'aian_softmax_alpha_pepfulldensity_20172018', models = 'acs_2018,pep_2018,acs_2017,pep_2017 ,wp')
  
}

#
#### 1/13/2025 Re-doing main models with CV = 10 ####
## 6 models to run.
## Make sure to increase total run time!

# Full pop, CV10, directest, interceptonly
bash_wrapper_real(dataset = 'all', bash_file = 'code/bash_commands/real_data_CV10models_01132025.txt', use_softmax = F, fixed_effects = 'intercept', CV_blocks = 10, output_path_addition = 'directest_CV10_interceptonly')

# Full pop, CV10, softmax, alpha, pep_density (double check param name)
bash_wrapper_real(dataset = 'all', bash_file = 'code/bash_commands/real_data_CV10models_01132025.txt', use_softmax = T, alpha_variance_prior = .01, fixed_effects = 'pep_density', CV_blocks = 10, output_path_addition = 'softmax_alpha_CV10_density')

# Full pop, CV10, softmax, density
bash_wrapper_real(dataset = 'all', bash_file = 'code/bash_commands/real_data_CV10models_01132025.txt', use_softmax = T, fixed_effects = 'pep_density', CV_blocks = 10,output_path_addition = 'softmax_CV10_density')

# AIAN, CV10, directest, interceptonly
bash_wrapper_real(dataset = 'aian', bash_file = 'code/bash_commands/real_data_CV10models_01132025.txt', use_softmax = F, fixed_effects = 'intercept', CV_blocks = 10, output_path_addition = 'AIAN_directest_CV10_interceptonly')

# AIAN,  CV10, softmax, alpha, pep_fulldensity
bash_wrapper_real(dataset = 'aian', bash_file = 'code/bash_commands/real_data_CV10models_01132025.txt', use_softmax = T, alpha_variance_prior = .01, fixed_effects = 'pep_fulldensity', CV_blocks = 10, output_path_addition = 'AIAN_softmax_alpha_CV10_pepfulldensity')

# AIAN, CV10, softmax, pep_fulldensity
bash_wrapper_real(dataset = 'aian', bash_file = 'code/bash_commands/real_data_CV10models_01132025.txt', use_softmax = T, fixed_effects = 'pep_fulldensity', CV_blocks = 10, output_path_addition = 'AIAN_softmax_CV10_pepfulldensity')


#
#### 5 models runs - take 2, with no preprocess, and with PCA ####
bash_wrapper_real(bash_file = 'code/bash_commands/real_data_01052025.txt', use_softmax = F, models = 'acs,pep,wp,acs_2018,pep_2018', fixed_effects = 'intercept', preprocess_scale = F, output_path_addition = 'directest_interceptonly_5models')

bash_wrapper_real(bash_file = 'code/bash_commands/real_data_01052025.txt', use_softmax = F, models = c('acs,pep,wp,acs_diff,pep_diff'), fixed_effects = 'intercept', preprocess_scale = F, output_path_addition = 'directest_interceptonly_5modelsDiff')

bash_wrapper_real(bash_file = 'code/bash_commands/real_data_01052025.txt', use_softmax = F, models = c('PC1,PC2,PC3,PC4,PC5'), fixed_effects = 'intercept', preprocess_scale = F, output_path_addition = 'directest_interceptonly_5modelsPCA')

bash_wrapper_real(bash_file = 'code/bash_commands/real_data_01052025.txt', use_softmax = T, models = c('acs,pep,wp,acs_2018,pep_2018'), fixed_effects = 'intercept', preprocess_scale = F, alpha_variance_prior = 1, output_path_addition = 'softmax_interceptonly_5models')

bash_wrapper_real(bash_file = 'code/bash_commands/real_data_01052025.txt', use_softmax = T, models = c('acs,pep,wp,acs_diff,pep_diff'), fixed_effects = 'intercept', preprocess_scale = F, alpha_variance_prior = 1, output_path_addition = 'softmax_interceptonly_5modelsDiff')

bash_wrapper_real(bash_file = 'code/bash_commands/real_data_01052025.txt', use_softmax = T, models = 'PC1,PC2,PC3,PC4,PC5', fixed_effects = 'intercept', preprocess_scale = F, alpha_variance_prior = 1, output_path_addition = 'softmax_interceptonly_5modelsPCA')

#
#### TESTING ####
bash_command_real(use_softmax = F, models = c('acs,pep,wp,acs_2018,pep_2018'), fixed_effects = 'intercept', preprocess_scale = T, output_path_addition = 'directest_preprocess_interceptonly_5models', CV_blocks = NULL, n.sample = 100, burnin = 50, chains_cores = 1)


#### Real bash - Full pop with 5 models ####
bash_wrapper_real(bash_file = 'code/bash_commands/real_data_01032025.txt', use_softmax = F, models = c('acs,pep,wp,acs_2018,pep_2018'), fixed_effects = 'intercept', preprocess_scale = T, output_path_addition = 'directest_preprocess_interceptonly_5models')

bash_wrapper_real(bash_file = 'code/bash_commands/real_data_01032025.txt', use_softmax = T, models = c('acs,pep,wp,acs_2018,pep_2018'), fixed_effects = 'intercept', preprocess_scale = T, output_path_addition = 'softmax_preprocess_interceptonly_5models')

bash_wrapper_real(bash_file = 'code/bash_commands/real_data_01032025.txt', use_softmax = F, models = c('acs,pep,wp,acs_diff,pep_diff'), fixed_effects = 'intercept', preprocess_scale = T, output_path_addition = 'directest_preprocess_interceptonly_5modelsDiff')

bash_wrapper_real(bash_file = 'code/bash_commands/real_data_01032025.txt', use_softmax = T, models = c('acs,pep,wp,acs_diff,pep_diff'), fixed_effects = 'intercept', preprocess_scale = T, output_path_addition = 'softmax_preprocess_interceptonly_5modelsDiff')

bash_wrapper_real(bash_file = 'code/bash_commands/real_data_01032025.txt', dataset = 'AIAN', use_softmax = F, models = c('acs,pep,wp,acs_2018,pep_2018'), fixed_effects = 'intercept', preprocess_scale = T, output_path_addition = 'AIAN_directest_preprocess_interceptonly_5models')

bash_wrapper_real(bash_file = 'code/bash_commands/real_data_01032025.txt', dataset = 'AIAN', use_softmax = T, models = c('acs,pep,wp,acs_2018,pep_2018'), fixed_effects = 'intercept', preprocess_scale = T, output_path_addition = 'AIAN_softmax_preprocess_interceptonly_5models')

bash_wrapper_real(bash_file = 'code/bash_commands/real_data_01032025.txt', dataset = 'AIAN', use_softmax = F, models = c('acs,pep,wp,acs_diff,pep_diff'), fixed_effects = 'intercept', preprocess_scale = T, output_path_addition = 'AIAN_directest_preprocess_interceptonly_5modelsDiff')

bash_wrapper_real(bash_file = 'code/bash_commands/real_data_01032025.txt', dataset = 'AIAN',use_softmax = T, models = c('acs,pep,wp,acs_diff,pep_diff'), fixed_effects = 'intercept', preprocess_scale = T, output_path_addition = 'AIAN_softmax_preprocess_interceptonly_5modelsDiff')

#
#### Testing the theta multiplier ####
bash_wrapper_real(use_softmax = T, 
                  preprocess_scale = F, 
                  theta_multiplier = 1000,
                  theta_prior_shape = 0.001, 
                  theta_prior_rate = 0.001,
                  fixed_effects = 'intercept',
                  output_path_addition = 'softmax_interceptonly_theta1000_001')

bash_wrapper_real(use_softmax = T, 
                  preprocess_scale = F, 
                  theta_multiplier = 1000,
                  theta_prior_shape = 1, 
                  theta_prior_rate = 1,
                  fixed_effects = 'intercept',
                  output_path_addition = 'softmax_interceptonly_theta1000_1')

#### Real bash - 8 models using different fixed effects ####
#(1) Full pop, softmax, preprocess, intercept + density.
bash_1118 <- 'code/bash_commands/real_data_bash_11182024.txt'

bash_wrapper_real(use_softmax = T, 
                  preprocess_scale = T, 
                  fixed_effects = 'pep_density',
                  output_path_addition = 'softmax_preprocess_density', 
                  bash_file = bash_1118)

# (2) Full pop, softmax, alpha, intercept + density.
bash_wrapper_real(use_softmax = T, 
                  alpha_variance_prior = .01,
                  fixed_effects = 'pep_density',
                  output_path_addition = 'softmax_alpha_density',
                  bash_file = bash_1118)

# (3) Full pop, direct est, no prep, intercept + density.
bash_wrapper_real(use_softmax = F, 
                  fixed_effects = 'pep_density',
                  output_path_addition = 'directest_density',
                  bash_file = bash_1118)

# (4) AIAN, softmax, preprocess, intercept.
bash_wrapper_real(dataset = 'AIAN',
                  use_softmax = T, 
                  preprocess_scale = T, 
                  fixed_effects = 'intercept',
                  output_path_addition = 'AIAN_softmax_interceptonly',
                  bash_file = bash_1118)

# (5) AIAN, softmax, preprocess, AIAN density.
bash_wrapper_real(dataset = 'AIAN',
                  use_softmax = T, 
                  preprocess_scale = T, 
                  fixed_effects = 'pep_density',
                  output_path_addition = 'AIAN_softmax_AIANdensity',
                  bash_file = bash_1118)

# (6) AIAN, softmax, preprocess, fullpop density.
bash_wrapper_real(dataset = 'AIAN',
                  use_softmax = T, 
                  preprocess_scale = T, 
                  fixed_effects = 'pep_fulldensity',
                  output_path_addition = 'AIAN_softmax_FULLdensity',
                  bash_file = bash_1118)

# (7) AIAN, softmax, preprocess, AIAN density and AIAN proportion
bash_wrapper_real(dataset = 'AIAN',
                  use_softmax = T, 
                  preprocess_scale = T, 
                  fixed_effects = 'pep_density_proportion',
                  output_path_addition = 'AIAN_softmax_AIANdensityANDprop',
                  bash_file = bash_1118)

# (8) AIAN, softmax, preprocess, full density and AIAN proportion
bash_wrapper_real(dataset = 'AIAN',
                  use_softmax = T, 
                  preprocess_scale = T, 
                  fixed_effects = 'pep_fulldensity_proportion',
                  output_path_addition = 'AIAN_softmax_FULLdensityANDprop',
                  bash_file = bash_1118)

#
#### Real bash - Full pop with intercept ####

## running the top three versions with the intercept - woohoo!
bash_command_real(use_softmax = F, fixed_effects = 'intercept', output_path_addition = 'directest_interceptonly')
bash_command_real(use_softmax = T, fixed_effects = 'intercept', output_path_addition = 'softmax_interceptonly')
bash_command_real(use_softmax = T, fixed_effects = 'intercept', alpha_variance_prior = 0.01,  output_path_addition = 'softmax_alpha_interceptonly')
bash_command_real(use_softmax = T, fixed_effects = 'intercept', preprocess_scale = T,  output_path_addition = 'softmax_preprocess_interceptonly')
bash_command_real(use_softmax = T, fixed_effects = 'intercept', alpha_variance_prior = 0.01,  output_path_addition = 'softmax_alpha_interceptonly')
bash_command_real(use_softmax = T, fixed_effects = 'intercept', preprocess_scale = T, alpha_variance_prior = 0.01,  output_path_addition = 'softmax_preprocess_alpha_interceptonly')

#### Real bash - AIAN ####
bash_command_real(dataset = 'AIAN', use_softmax = T, output_path_addition = 'softmax_AIAN')

bash_command_real(dataset = 'AIAN', use_softmax = T, preprocess_scale = T,  output_path_addition = 'softmax_preprocess_AIAN')

bash_command_real(dataset = 'AIAN', use_softmax = T, alpha_variance_prior = 0.01,  output_path_addition = 'softmax_alpha_AIAN')

bash_command_real(dataset = 'AIAN', use_softmax = T, preprocess_scale = T, alpha_variance_prior = 0.01,  output_path_addition = 'softmax_preprocess_alpha_AIAN')

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
