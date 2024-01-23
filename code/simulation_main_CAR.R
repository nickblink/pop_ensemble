library(dplyr)
library(rstan)
library(ggplot2)

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
models = c('acs','pep')

## subset data by state
NY_lst <- subset_data_by_state(D2010, county_adj, 'New York', 'NY')

## Trying with different y values
data <- simulate_models(data = NY_lst$data, models = c('acs','pep'), means = c(100, 200), variances = c(10^2, 10^2))

data_lst <- simulate_data(data, NY_lst$adjacency, models = models, precision_type = 'Leroux', tau2 = 1, rho = 0.3)

stan_fit <- run_stan_CAR(data_lst$data, data_lst$adjacency, models = models, n.sample = 10000, burnin = 5000)

save(stan_fit, data_lst, file = 'C:/Users/Admin-Dell/Dropbox/Academic/HSPH/Research/Population Estimation/Results/results_sim_01232024.RData')
# ('C:/Users/nickl/Dropbox/Academic/HSPH/Research/Population Estimation/Results/results_sim_01232024.RData')

# process the results
process_results <- function(data_list, models, stan_fit){
  N = nrow(data_list$data)
  
  ## grab the results
  stan_summary = summary(stan_fit, pars = c('tau2','rho', 'phi', 'u','y_exp','lp__'))$summary
  stan_out <- extract(stan_fit)
  
  ## get the convergence parameters
  # ESS of tau2 and rho
  ESS_spatial <- data.frame(stan_summary[1:(2*length(models)), c('n_eff', 'Rhat')])
  print(ESS_spatial)
  
  # ESS of phi
  ind = grep('phi', rownames(stan_summary))
  # ind1 = grep('phi\\[[0-9]{1,2},1\\]', rownames(stan_summary))
  # ind2 = grep('phi\\[[0-9]{1,2},2\\]', rownames(stan_summary))
  print(sprintf('median ESS for phi is %s and median rhat is %s', round(median(stan_summary[ind, 'n_eff']), 1), round(median(stan_summary[ind, 'Rhat']), 3)))
  x = data.frame(n_eff = stan_summary[ind, 'n_eff'])
  p_ESS <- ggplot(data = x, aes(x = n_eff)) + 
    geom_density() +
    scale_x_continuous(trans='log2') + 
    ggtitle('ESS of phis')
  
  ## get the posterior split into likelihood and log posterior (meh)
  
  ## get the spatial parameter estimates
  tau2 <- NULL
  rho <- NULL
  for(i in 1:length(models)){
    tau2 <- rbind(tau2, 
                  data.frame(value = stan_out$tau2[,i], 
                             model = models[i]))
    rho <- rbind(rho, 
                  data.frame(value = stan_out$rho[,i], 
                             model = models[i]))
  }
  
  # get the true spatial params
  true_vals <- data.frame(model = models, 
                          tau2 = data_lst$tau2, 
                          rho = data_lst$rho)
  
  # plot the tau2 param
  p_tau2 <- ggplot(data = tau2, aes(x = model, y = value)) + 
    geom_boxplot() + 
    geom_point(data = true_vals, aes(x = model, y = tau2, col = 'red')) + 
    ggtitle('tau2 estimates') + 
    theme(legend.position = 'none')
  
  # plot the rho param 
  p_rho <- ggplot(data = rho, aes(x = model, y = value)) + 
    geom_boxplot() + 
    geom_point(data = true_vals, aes(x = model, y = rho, col = 'red')) + 
    ggtitle('rho estimates') + 
    theme(legend.position = 'none')
  
  ## compare the true phi values with the estimated phi values
  {
    phi_est <- as.data.frame(matrix(0, nrow = N, ncol = length(models)))
    colnames(phi_est) <- models
    for(i in 1:length(models)){
      ind = grep(sprintf('phi\\[[0-9]{1,2},%s\\]', i), rownames(stan_summary))
      phi_est[,i] <- stan_summary[ind,'50%']
    }
    phi_est$index = 1:N
    
    # convert estimates to long
    phi_est_long <- tidyr::gather(phi_est, key = model, value = phi_median_est, -index)
    
    # convert true vals to long
    phi_true_long <- tidyr::gather(data_list$phi_true, key = model, value = phi_true, -index)
    
    # merge them
    phi_mat <- merge(phi_est_long, phi_true_long, by = c('index','model'))
  
    # plot 'em
    p_phi <- ggplot(phi_mat, aes(x = phi_true, phi_median_est)) + 
      geom_point() + 
      geom_smooth(method='lm') + 
      geom_abline(slope = 1, intercept = 0, col = 'red') + 
      facet_wrap(~model) + 
      xlab('true phi values') + 
      ylab('median estimated phi values') + 
      ggtitle('phi estimates')
  }
  
  ## compare the true u values with the estimated u values
  {
    u_est <- as.data.frame(matrix(0, nrow = N, ncol = length(models)))
    colnames(u_est) <- models
    for(i in 1:length(models)){
      ind = grep(sprintf('u\\[[0-9]{1,2},%s\\]', i), rownames(stan_summary))
      u_est[,i] <- stan_summary[ind,'50%']
    }
    u_est$index = 1:N
    
    # convert estimates to long
    u_est_long <- tidyr::gather(u_est, key = model, value = u_median_est, -index)
    
    # convert true vals to long
    u_true_long <- tidyr::gather(data_list$u_true, key = model, value = u_true, -index)
    
    # merge them
    u_mat <- merge(u_est_long, u_true_long, by = c('index','model'))
    
    # plot 'em
    p_u <- ggplot(u_mat, aes(x = u_true, u_median_est)) + 
      geom_point() + 
      geom_smooth(method='lm') + 
      geom_abline(slope = 1, intercept = 0, col = 'red') + 
      facet_wrap(~model) + 
      xlab('true u values') + 
      ylab('median estimated u values') + 
      ggtitle('u estimates')
  }
  
  ## compare the true outcomes with the estimated outcomes (compared to just using one model in the outcomes)
  ind = grep('y_exp', rownames(stan_summary))
  data_list$data$y_predicted <- stan_summary[ind,'50%']
  p_y <- ggplot(data_list$data, aes(x = y, y_predicted)) + 
    geom_point() +
    geom_smooth(method='lm') + 
    geom_abline(slope = 1, intercept = 0, col = 'red') + 
    xlab('observed y') + 
    ylab('median estimated y') + 
    ggtitle('y estimation')
  
  full_plot <- cowplot::plot_grid(p_ESS,
                                  cowplot::plot_grid(p_rho, p_tau2),
                                  p_phi,
                                  p_u,
                                  p_y, 
                                  ncol = 1)
  
  return(full_plot)
}

pp <- process_results(data_lst, models, stan_fit)


# Now I want to be able to simulate all sorts of data:
# - different CAR precision matrices (could be an input? Like Leroux or BYM)
# - different # of models
# - simulated values of models
# - being able to scale down the values in the data
# - adding some noise to the models (no need right now)


# previous runs
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