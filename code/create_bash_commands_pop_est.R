## This script makes bash commands for given simulations

bash_command <- function(R=40, dataset='NY', N_models=3, n.sample=10000, burnin=5000, family='normal',use_softmax=F,variances=c(100,100,100), means=c(100,100,100), rho = 0.3, tau2 = 1, sigma2 = 100, sigma2_prior_shape = 50, sigma2_prior_rate = 0.5, tau2_prior_shape = 1, tau2_prior_rate=1, num_y_samples=3, stan_path='code/CAR_leroux_sparse_normal.stan', CV_blocks = 5, return_quantiles = T,parallel = T, output_path = NULL, array_length = 5){
  
  # make the output path
  if(is.null(output_path)){
    output_path <- sprintf('simulation_%s_%s_%smodels_CV%s_ID%s_%s', ifelse(use_softmax, 'softmax', 'direct_est'), family, N_models, CV_blocks, sample(1e6, size = 1), gsub('-','_',Sys.Date()))
    output_path <- gsub('\\.','',output_path)
  }
  
  job_name = sprintf('%s_%s_%smodels_CV%s', ifelse(use_softmax, 'softmax', 'directest'), family, N_models, CV_blocks)
  
  params <- list(R=R, dataset=dataset, N_models=N_models, n.sample=n.sample, burnin=burnin, family=family,use_softmax=F,variances=variances, means=means, rho = rho, tau2 = tau2, sigma2 = sigma2, sigma2_prior_shape = sigma2_prior_shape, sigma2_prior_rate = sigma2_prior_rate, tau2_prior_shape = tau2_prior_shape, tau2_prior_rate=tau2_prior_rate, num_y_samples=num_y_samples, stan_path=stan_path, CV_blocks = CV_blocks, return_quantiles = return_quantiles, parallel = parallel, output_path = output_path)
  
  param_str = paste(paste(names(params), params, sep = '='), collapse=':') %>%
    gsub(' |c\\(|\\)','',.) %>%
    gsub('TRUE','T',.) %>%
    gsub('FALSE','F',.)

  command_str = sprintf('sbatch --array=1-%s -J %s run_sim_pop_est.bash %s', array_length, job_name, param_str)
  
  return(command_str)

}

bash_command(CV_blocks = 5)

bash_command(CV_blocks = 10)

