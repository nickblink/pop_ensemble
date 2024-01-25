library(cowplot)

### Subset full dataset and adjacency by a specific state
subset_data_by_state <- function(data, adjacency, state, abbrev = NULL){
  # get the indices corresponding to the state
  ind1 <- grep(state, data$NAME)
  ind2 <- grep(state, rownames(adjacency))
  
  # check that the indices match between data and adjacency matrix
  if(!identical(ind1, ind2)){
    stop('indices of data names and adjacency names do not match')
  }
  
  # check abbreviated adjacency column names
  if(!is.null(abbrev)){
    ind3 <- grep(abbrev, colnames(adjacency))
    if(!identical(ind1, ind3)){
      stop('indices of data names and adjacency columnss do not match')
    }
  }
  
  # subset the data
  data_subset <- data[ind1, ]
  adjacency_subset <- adjacency[ind1, ind1]
  
  # return it!
  return(list(data = data_subset, adjacency = adjacency_subset))
}


### simulate the numbers in the data according to a normal distribution
# data: the input data
# models: the models to simulate data for
# means: the means of the normals for each model
# variances: the variances of the normals for each model
simulate_models <- function(data, models, means, variances, seed = 10){
  set.seed(seed)
  # check that the lengths of the inputs match
  if(length(models) != length(means) | length(models) != length(variances)){
    stop('length of models, means, and variances, not matching up')
  }
  
  # sample data for each model
  for(i in 1:length(models)){
    m = models[i]
    data[,m] <- rnorm(nrow(data), means[i], sd = sqrt(variances[i]))
  }
  
  return(data)
}
  

### Generate the precision matrix based. Can generate an Cressie, Cressie, or Leroux precision matrix
# type: "Cressie" or "Leroux", for the type of precision matrix
# rho: the spatial correlation parameter, ranging from 0 to 1
# tau2: the variance parameter
generate_precision_mat <- function(W, type, tau2, rho){
  D = diag(rowSums(W))
  if(type == 'Cressie'){
    Q = (D - rho*W)
  }else if(type == 'Leroux'){
    Q = (rho*(D - W) + (1-rho)*diag(rep(1,nrow(D))))
  }else{
    stop('please input a proper precision matrix type')
  }
  Q = Q/tau2
  return(Q)
}


### Sampling a multivariate normal from a precision matrix using cholesky decomposition
sample_MVN_from_precision <- function(n = 1, mu=rep(0, nrow(Q)), Q){
  p <- length(mu)
  Z <- matrix(rnorm(p*n), p, n)
  U <- chol(Q) # By default R's chol fxn returns upper cholesky factor
  X <- backsolve(U, Z) # more efficient and stable than acctually inverting
  X <- sweep(X, 1, mu, FUN=`+`)
  return(X)
}


### Simulate the data!
# data: the input data
# adjacency: the adjacency matrix 
# models: the models to use in the ensemble
# scale_down: a factor to scale down the covariate values
# pivot: what pivot index to use for data creation (-1 indicates no pivot)
# precision_type: "Cressie" or "Leroux", determining the CAR precision matrix.
# tau2: the CAR variance parameter
# rho: the CAR spatial correlation parameter
# seed: random seed to initialize function with
simulate_data <- function(data, adjacency, models = c('acs','pep','worldpop'), scale_down = 1, pivot = -1, precision_type = 'Cressie', tau2 = 1, rho = 0.3, seed = 10, cholesky = T){
  # set seed for reproducability 
  set.seed(seed)
  
  # transform to be for each model
  if(length(tau2) == 1 & length(models) > 1){
    tau2 = rep(tau2, length(models))
  }
  if(length(rho) == 1 & length(models) > 1){
    rho = rep(rho, length(models))
  }
  
  # scale down the data size
  data[,models] <- data[,models]/scale_down
  
  # initialize phi_true
  phi_true <- matrix(0, nrow = nrow(data), ncol = length(models))
  
  # cycle through models to generate phi values
  for(i in 1:length(models)){
    # create the precision matrix
    Q <- generate_precision_mat(W = adjacency, type = precision_type, tau2 = tau2[i], rho = rho[i])
    
    if(pivot == -1){
      # sample from MVN based on the precision matrix
      if(cholesky){
        phi_true[,i] <- sample_MVN_from_precision(n = 1, Q = Q)
      }else{
        Sigma <- solve(Q)
        phi_true[,i] <- MASS::mvrnorm(n = 1, mu = rep(0,nrow(Sigma)), Sigma = Sigma)
      }
    }else{
      stop('havent coded for pivot data generation')
    }
  }
  colnames(phi_true) <- models
    
  # get exponentiated values and sum across models
  exp_phi = exp(phi_true)
  exp_phi_rows = rowSums(exp_phi)
  
  # get model weights and calculate the mean estimate
  u_true <- exp_phi/exp_phi_rows
  
  # get the expected census values
  data$y_expected <- rowSums(u_true*data[,models])
  
  # simulate the y values
  data$y <- rpois(n = nrow(data), lambda = data$y_expected)
  
  data$index = 1:nrow(data)
  phi_true <- as.data.frame(phi_true) %>%
    mutate(index = 1:nrow(phi_true))
  u_true <- as.data.frame(u_true) %>%
    mutate(index = 1:nrow(u_true))
  return(list(data = data, adjacency = adjacency, phi_true = phi_true, u_true = u_true, tau2 = tau2, rho = rho))
}


### prep the data for fitting the stan model
# data: the observed data with the outcome "y" and the covariates as models. 
# W: the adjacency matrix
# models: a vector of the models to use in the ensemble
prep_stan_data_leroux_sparse <- function(data, W, models){
  # checking columns
  if(!('y' %in% colnames(data))){
    stop('need y as a column in data')
  }
  if(any(!(models %in% colnames(data)))){
    stop('models input do not exist in data')
  }
  
  N = nrow(data)
  
  W_n = sum(W)/2
  # W2 = make_district_W2_matrix(df)
  
  # get eigenvalues for determinant calculation
  lambda = eigen(diag(rowSums(W)) - W - diag(1, nrow(W)))$values
  
  # make the design matrix
  X = as.matrix(data[,models])
  
  # get the outcome
  y = data$y
  
  # # comparing missingness.
  # if(!identical(as.integer(rownames(X_obs)), which(!is.na(df$y)))){
  #   stop('mismatch of model matrix and df missing rows')
  # }
  
  # missingness data
  N_miss = sum(is.na(y))
  N_obs = sum(!is.na(y))
  ind_miss = which(is.na(y))
  ind_obs = which(!is.na(y))
  y_obs = y[ind_obs]
  
  # make the stan data frame
  stan_data <- list(
    M = length(models), # number of models
    N = N, # number of observations
    N_miss = N_miss, # number of missing y points
    N_obs = N_obs,
    ind_miss = ind_miss, # indices of missing y points
    ind_obs = ind_obs,
    X = X, # design matrix
    y_obs = y_obs, # outcome variable 
    W = W,
    W_n = W_n,
    I = diag(1.0, N),
    lambda = lambda)
  
  return(stan_data)
}


### Fits the CAR model using rstan. Prepares the data for rstan and runs it.
# data: input data with output column y and covariate columns according to models.
# adjacency: the adjacency matrix for the data.
# models: the models for the ensemble.
# precision_type: Cressie or Leroux, for the type of precision matrix.
# n.sample: the number of iterations to run the rstan code.
# burnin: the number of burnin iterations to run the rstan code.
# seed: a seed for reproducability
run_stan_CAR <- function(data, adjacency, models = c('acs','pep','worldpop'), precision_type = 'Leroux', n.sample = 10000, burnin = 5000, seed = 10){
  # error checking for precision matrix type
  if(precision_type != 'Leroux'){stop('only have Leroux precision coded')}
  
  # prep the data
  stan_data <- prep_stan_data_leroux_sparse(data, adjacency, models)
  
  # fit the stan model
  stan_fit <- stan(file = "code/CAR_leroux_sparse.stan",
                   data = stan_data, 
                   iter = n.sample, 
                   warmup = burnin,
                   chains = 1, 
                   init = '0',
                   cores = 1,
                   seed = seed)
  
  # # extract important info
  # stan_out <- extract(stan_fit)
  # stan_summary = summary(stan_fit, pars = c('tau2','rho', 'phi'))$summary
  # stan_lst <- list(stan_fit = stan_fit,
  #                  stan_out = stan_out, 
  #                  stan_summary = stan_summary)
  #return(stan_lst)
  return(stan_fit)
}


### process the results. Prints ESS for spatial params and returns a plot of various results for parameter estimates
# data_list: List containing data used for the fit, the true phi values, and the true u values.
# models: Vector of models used for fitting.
# stan_fit: The result of the stan fit on this data.
process_results <- function(data_list, models, stan_fit, ESS = T, likelihoods = T, spatial_params = T, phi_estimates = T, u_estimates = T, y_estimates = T){
  N = nrow(data_list$data)
  
  plot_list = NULL
  
  ## grab the results
  stan_summary = summary(stan_fit, pars = c('tau2','rho', 'phi', 'u','y_exp','lp__'))$summary
  stan_out <- extract(stan_fit)
  
  ## get the convergence parameters
  if(ESS){
    # ESS of tau2 and rho
    ESS_spatial <- data.frame(stan_summary[1:(2*length(models)), c('n_eff', 'Rhat')])
    colnames(ESS_spatial)[1] <- 'ESS'
    p_ESS_spatial <- gridExtra::tableGrob(round(ESS_spatial, 3))
    
    # ESS of phi
    ind = grep('phi', rownames(stan_summary))
    # ind1 = grep('phi\\[[0-9]{1,2},1\\]', rownames(stan_summary))
    # ind2 = grep('phi\\[[0-9]{1,2},2\\]', rownames(stan_summary))
    print(sprintf('median ESS for phi is %s and median rhat is %s', round(median(stan_summary[ind, 'n_eff']), 1), round(median(stan_summary[ind, 'Rhat']), 3)))
    x = data.frame(n_eff = stan_summary[ind, 'n_eff'])
    p_ESS_phi <- ggplot(data = x, aes(x = n_eff)) + 
      geom_density() +
      scale_x_continuous(trans='log2') + 
      ggtitle('ESS of phis')
    
    plot_list = append(plot_list, list(plot_grid(p_ESS_spatial, p_ESS_phi)))
  }
  
  ## get the posterior split into likelihood and log posterior
  if(likelihoods){
    sq = round(seq(1, length(stan_out$lp__), length.out = 200))
    lkl <- data.frame(i = sq, 
                      log_post = stan_out$lp__[sq],
                      log_likelihood = stan_out$log_likelihood[sq])
    lkl$log_prior = lkl$log_post - lkl$log_likelihood
    p1 <- ggplot(lkl, aes(x = i, y = log_post)) + 
      geom_point() + 
      geom_smooth()
    p2 <- ggplot(lkl, aes(x = i, y = log_likelihood)) + 
      geom_point() + 
      geom_smooth()
    p3 <- ggplot(lkl, aes(x = i, y = log_prior)) + 
      geom_point() + 
      geom_smooth()
    lkl_plot <- plot_grid(p1, p2, p3, nrow = 1)
    plot_list = c(plot_list, list(lkl_plot))
  }
  
  ## get the spatial parameter estimates
  if(spatial_params){
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
                            tau2 = data_list$tau2, 
                            rho = data_list$rho)
    
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
    
    plot_list <- append(plot_list, list(plot_grid(p_rho, p_tau2)))
  }
  
  ## compare the true phi values with the estimated phi values
  if(phi_estimates){
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
    
    plot_list <- append(plot_list, list(p_phi))
  }
  
  ## compare the true u values with the estimated u values
  if(u_estimates){
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
    
    plot_list <- append(plot_list, list(p_u))
  }
  
  ## compare the true outcomes with the estimated outcomes (compared to just using one model in the outcomes)
  if(y_estimates){
    ind = grep('y_exp', rownames(stan_summary))
    data_list$data$y_predicted <- stan_summary[ind,'50%']
    p_y <- ggplot(data_list$data, aes(x = y, y_predicted)) + 
      geom_point() +
      geom_smooth(method='lm') + 
      geom_abline(slope = 1, intercept = 0, col = 'red') + 
      xlab('observed y') + 
      ylab('median estimated y') + 
      ggtitle('y estimation')
    
    plot_list <- append(plot_list, list(p_y))
  }
  
  ## Combine all the plots!
  full_plot <- plot_grid(plotlist = plot_list,
                         ncol = 1)
  
  return(full_plot)
}


  