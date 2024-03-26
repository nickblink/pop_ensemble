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
# seed: random seed
# adjacency: adjacency matrix. Only relevant for CAR distribution simulation.
# M_CAR_rho: rho for the CAR distribution of each model. If null, CAR is not added to each model.
# M_CAR_tau2: tau2 for the CAR distribution of each model. If null CAR is not added to each model.
# M_BYM_variance: if running CAR, to add in normal variance as well (== T) or just CAR variance (== F).
# precision_type: precision type for the CAR model. Only relevant if running CAR simulation.
# M_MVN_alpha: value of multivariate normal covariance level (same covariance for all off-diagonals).

simulate_models <- function(data, models, means, variances, seed = 10, adjacency = NULL, M_CAR_rho = NULL, M_CAR_tau2 = NULL, M_BYM_variance = F, precision_type = 'Leroux', M_MVN_alpha = NULL, ...){
  # set random seed.
  set.seed(seed)
  
  run_CAR = F
  # check if running CAR
  if(is.null(M_CAR_rho) + is.null(M_CAR_tau2) == 1){
    stop('either M_CAR_rho and M_CAR_tau2 should both be NULL or both be values.')
  }else if(is.null(M_CAR_rho) + is.null(M_CAR_tau2) == 0){
    print('running CAR distribution simulation for model values.')
    run_CAR = T
    if(M_CAR_rho < 0 | M_CAR_rho > 1 | M_CAR_tau2 < 0){
      stop('Improper CAR values.')
    }
    if(!is.null(M_MVN_alpha)){
      stop('cant run M_MVN simulation and CAR simulation right now.')
    }
  }
  
  # check that the lengths of the inputs match.
  if(length(models) != length(means) | length(models) != length(variances)){
    stop('length of models, means, and variances, not matching up')
  }
  
  # sample from MVN value
  if(!is.null(M_MVN_alpha)){
    # create the covariance matrix
    Sigma = diag(variances)
    Sigma[Sigma == 0] <- M_MVN_alpha
    
    # sample the data
    data[,models] = MASS::mvrnorm(n = nrow(data), mu = means, Sigma = Sigma)
  }else{
    # sample data independently for each model.
    for(i in 1:length(models)){
      m = models[i]
      if(run_CAR){
        # create the precision matrix
        Q <- generate_precision_mat(W = adjacency, type = precision_type, tau2 = M_CAR_tau2, rho = M_CAR_rho)
        
        # sample phi values
        phi <- sample_MVN_from_precision(n = 1, Q = Q)
        
        # sample model values according to BYM or just CAR 
        if(M_BYM_variance){
          data[,m] <- rnorm(nrow(data), means[i], sd = sqrt(variances[i])) + phi
        }else{
          data[,m] <- means[i] + phi
        }
        
      }else{
        data[,m] <- rnorm(nrow(data), means[i], sd = sqrt(variances[i]))
      }
      
    }
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


### Simulate the outcome y
# data: the input data.
# adjacency: the adjacency matrix.
# models: the models to use in the ensemble.
# scale_down: a factor to scale down the covariate values.
# pivot: what pivot index to use for data creation (-1 indicates no pivot).
# precision_type: "Cressie" or "Leroux", determining the CAR precision matrix.
# tau2: the CAR variance parameter.
# rho: the CAR spatial correlation parameter.
# seed: random seed to initialize function with.
# cholesky: whether to simulate phi's from cholesky decomposition.
# family: Poisson or Normal/Gaussian, for the type of outcome family.
# sigma2: variance of the normal distribution.
# use_softmax: If true, softmax is used. If false, the phi values are directly used as weights (centered at 1/M).

simulate_y <- function(data, adjacency, models = c('M1','M2','M3'), scale_down = 1, pivot = -1, precision_type = 'Leroux', tau2 = 1, rho = 0.3, seed = 10, cholesky = T, family = 'poisson', sigma2 = NULL, use_softmax = NULL, num_y_samples = 1, use_pivot = F, ...){
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
    
  if(use_softmax){
    # get exponentiated values and sum across models
    exp_phi = exp(phi_true)
    exp_phi_rows = rowSums(exp_phi)
    
    # get model weights and calculate the mean estimate
    u_true <- exp_phi/exp_phi_rows
  }else{
    # directly use phi values to create model weights
    if(use_pivot){
      u_true = matrix(NA, 
                      nrow = nrow(data), 
                      ncol = length(models))
      colnames(u_true) <- models
      u_true[,1:(length(models) - 1)] = 1/length(models) + phi_true[,1:(length(models) - 1)]
      u_true[,length(models)] = 1 - rowSums(u_true[,1:(length(models) - 1)])
    }else{
      u_true = 1/length(models) + phi_true
    }
  }
  
  # get the expected census values
  data$y_expected <- rowSums(u_true*data[,models])
  
  # simulate the y values
  for(i in 1:num_y_samples){
    if(tolower(family) == 'poisson'){
      # first Poisson sample
      data$y <- rpois(n = nrow(data), lambda = data$y_expected)
      
      # sample repeated y instances
      if(i >= 2){
        data[,paste0('y',i)] <- rpois(n = nrow(data), lambda = data$y_expected)
      }
    }else if(tolower(family) %in% c('normal','gaussian')){
      if(is.null(sigma2)){
        stop('please put in a sigma2 value if simulating from normal')
      }
      # first Normal sample
      data$y <- rnorm(n = nrow(data), mean = data$y_expected, sd = sqrt(sigma2))
      
      # sample repeated y instances
      if(i >= 2){
        data[,paste0('y',i)] <- rnorm(n = nrow(data), mean = data$y_expected, sd = sqrt(sigma2))
      }
    }else{
      stop('please input a proper family')
    }
  }
  
  # create phi and u data frames
  data$index = 1:nrow(data)
  phi_true <- as.data.frame(phi_true) %>%
    mutate(index = 1:nrow(phi_true))
  u_true <- as.data.frame(u_true) %>%
    mutate(index = 1:nrow(u_true))
  
  # create the list to return
  res_list = list(data = data, adjacency = adjacency, phi_true = phi_true, u_true = u_true, tau2 = tau2, rho = rho)
  if(!is.null(sigma2)){
    res_list[['sigma2']] <- sigma2
  }
  
  return(res_list)
}


### prep the data for fitting the stan model
# data: the observed data with the outcome "y" and the covariates as models. 
# W: the adjacency matrix.
# models: a vector of the models to use in the ensemble.
# sigma2_prior_shape: Shape of the gamma distribution prior.
# sigma2_prior_rate: rate of the gamma distribution prior.
prep_stan_data_leroux_sparse <- function(data, W, models, use_softmax = F, use_pivot = F, use_normal = T, sigma2_prior_shape = 1, sigma2_prior_rate = 10, tau2_prior_shape = 1, tau2_prior_rate = 1, fix_rho_value = - 1, fix_tau2_value = -1, ...){

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
    lambda = lambda,
    use_softmax = as.integer(use_softmax),
    use_pivot = as.integer(use_pivot),
    use_normal = as.integer(use_normal),
    sigma2_prior_shape = sigma2_prior_shape,
    sigma2_prior_rate = sigma2_prior_rate,
    tau2_prior_shape = tau2_prior_shape,
    tau2_prior_rate = tau2_prior_rate,
    rho_value = fix_rho_value,
    tau2_value = fix_tau2_value)

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
run_stan_CAR <- function(data, adjacency, models = c('M1','M2','M3'), precision_type = 'Leroux', n.sample = 10000, burnin = 5000, seed = 10, stan_m = NULL, stan_path = "code/CAR_leroux_sparse.stan", tau2 = NULL, use_softmax = NULL, use_normal = T, use_pivot = F, init_vals = '0', ...){
  
  # error checking for precision matrix type
  if(precision_type != 'Leroux'){stop('only have Leroux precision coded')}
  
  # prep the data
  stan_data <- prep_stan_data_leroux_sparse(data, adjacency, models, use_softmax = use_softmax, use_normal = use_normal, use_pivot = use_pivot, ...)
  
  # # update fixed tau2 value to be of M dimensions.
  # if(!is.null(tau2)){
  #   if(length(tau2) == 1 & length(models) > 1){
  #     tau2 = rep(tau2, length(models))
  #   }
  # 
  #   # add in tau2 fixed
  #   stan_data$tau2 = tau2
  # }
  
  # create the stan model if not done already
  if(is.null(stan_m)){
    stan_m <- stan_model(stan_path)
  }

  # fit the stan model
  stan_fit <- rstan::sampling(object = stan_m,
                   data = stan_data, 
                   iter = n.sample, 
                   warmup = burnin,
                   chains = 1, 
                   init = init_vals,
                   cores = 1,
                   seed = seed)

  if(!use_softmax & tau2 > 0.1){
    print('Are you sure you want tau2 so high?')
  }
  
  # return the results!
  return(stan_fit)
}


### Run multiple simulation runs. This calls the data creation functions and run_stan_CAR.
# raw_data: data list containing "data" and "adjacency".
# models: list of input models to put in function.
# stan_path: path to stan code.
# means: means of the input models' data creation.
# variances: variances of the input models' data creation.
# family: family of y distribution for simulation.

### Optional arguments
# precision_type: Leroux or Cressie.
# tau2: scalar or vector of CAR variance parameter.
# rho: scalar or vector of spatial correlation parameter.
# n.sample: number of stan chain samples.
# burnin: length of burnin period for stan.
# sigma2: sigma2 value of y distribution.
multiple_sims <- function(raw_data, models, means, variances, family = 'poisson', N_sims = 10, stan_path = "code/CAR_leroux_sparse_poisson.stan", init_vals = '0', family_name_check = T, use_softmax = F, ...){
  
  ### Parameter error checks
  {
    # Checking the name and family match
    if(!grepl(family, stan_path) & family_name_check){
      stop('The code does not have the family name in it')
    }
    
    # checking that the rho in the DGP and the model fit are equal (if the rho is fixed in the model fit)
    if('rho' %in% names(list(...)) & 'fix_rho_value' %in% names(list(...))){
      if(list(...)$fix_rho_value > 1){
        stop('rho value cannot be greater than 1')
      }
      if(list(...)$fix_rho_value > 0 & list(...)$fix_rho_value != list(...)$rho){
        print('WARNING: rho in DGP and rho in model not equal')
      }
    }
    
    # checking that the rho in the DGP and the model fit are equal (if the rho is fixed in the model fit)
    if('tau2' %in% names(list(...)) & 'fix_tau2_value' %in% names(list(...))){
      if(list(...)$fix_tau2_value > 0 & list(...)$fix_tau2_value != list(...)$tau2){
        print('WARNING: tau2 in DGP and tau2 in model not equal')
      }
    }
  
  }
  # capture the arguments
  arguments <- match.call()
  
  # set use_normal variable
  if(tolower(family) %in% c('normal','gaussian')){
    use_normal = T
  }else{
    use_normal = F
  }
  
  # compile the stan program
  m <- stan_model(stan_path)

  # initialize results
  sim_lst <- list()
  
  # cycle through N simulations
  for(i in 1:N_sims){
    # simulate input model values
    data <- simulate_models(data = raw_data$data, adjacency = raw_data$adjacency, models = models, seed = i, means = means, variances = variances, ...)
    
    # simulate y values from input models
    data_lst <- simulate_y(data, raw_data$adjacency, models = models, seed = i, family = family, use_softmax = use_softmax, ...)
    
    # update the initialization to start at the true values.
    if(tolower(init_vals) == 't' | tolower(init_vals) == 'truth'){
      init_list = list(phi = as.matrix(data_lst$phi_true[,-ncol(data_lst$phi_true)]))
      init_vals <- function(){init_list}
    }
    
   
    # fit the Bayesian model
    stan_fit <- run_stan_CAR(data_lst$data, data_lst$adjacency, models = models, seed = i, stan_m = m, use_softmax = use_softmax, ...)
    
    # store results
    sim_lst[[i]] <- list(data_list = data_lst, stan_fit = stan_fit)
  }
  
  # store the final set of the results 
  res_lst <- list(sim_list = sim_lst, arguments = arguments, models = models)
  
  # return the results
  return(res_lst)
}


### process the results. Prints ESS for spatial params and returns a plot of various results for parameter estimates
# data_list: List containing data used for the fit, the true phi values, and the true u values.
# models: Vector of models used for fitting.
# stan_fit: The result of the stan fit on this data.
# ESS: Whether to include the ESS of the parameters.
# likelihoods: whether to include the prior, likelihood, and posterior
# <XX>_estimates: whether to include the <XX> estimates
process_results <- function(data_list, models, stan_fit, ESS = T, likelihoods = T, rho_estimates = T, tau2_estimates = T, sigma2_estimates = F, phi_estimates = T, u_estimates = T, y_estimates = T){
  N = nrow(data_list$data)
  
  plot_list = NULL
  
  ## grab the results
  #stan_summary = summary(stan_fit, pars = c('tau2','rho', 'phi', 'u','y_exp','lp__'))$summary
  stan_summary = summary(stan_fit)$summary
  stan_out <- extract(stan_fit)
  
  ## get the convergence parameters
  if(ESS){
    # ESS of tau2 and rho
    ESS_spatial <- data.frame(stan_summary[1:(2*length(models)), c('n_eff', 'Rhat')])
    colnames(ESS_spatial)[1] <- 'ESS'
    p_ESS_spatial <- gridExtra::tableGrob(round(ESS_spatial, 3))
    
    # ESS of phi
    ind = grep('^phi', rownames(stan_summary))
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
      geom_smooth(method = 'loess', formula = y ~ x)
    p2 <- ggplot(lkl, aes(x = i, y = log_likelihood)) + 
      geom_point() + 
      geom_smooth(method = 'loess', formula = y ~ x)
    p3 <- ggplot(lkl, aes(x = i, y = log_prior)) + 
      geom_point() + 
      geom_smooth(method = 'loess', formula = y ~ x)
    lkl_plot <- plot_grid(p1, p2, p3, nrow = 1)
    plot_list = c(plot_list, list(lkl_plot))
  }
  
  ## get the spatial parameter estimates
  if(rho_estimates){
    rho <- NULL
    for(i in 1:length(models)){
      rho <- rbind(rho, 
                   data.frame(value = stan_out$rho[,i], 
                              model = models[i]))
    }
    
    # get the true spatial params
    true_vals <- data.frame(model = models, 
                            rho = data_list$rho)
    
    # plot the rho param 
    p_rho <- ggplot(data = rho, aes(x = model, y = value)) + 
      geom_boxplot() + 
      geom_point(data = true_vals, aes(x = model, y = rho, col = 'red')) + 
      ggtitle('rho estimates') + 
      theme(legend.position = 'none')
  }
  
  if(tau2_estimates){
    tau2 <- NULL
    for(i in 1:length(models)){
      tau2 <- rbind(tau2, 
                    data.frame(value = stan_out$tau2[,i], 
                               model = models[i]))
    }
    
    # get the true spatial params
    true_vals <- data.frame(model = models, 
                            tau2 = data_list$tau2)
    
    # plot the tau2 param
    p_tau2 <- ggplot(data = tau2, aes(x = model, y = value)) + 
      geom_boxplot() + 
      geom_point(data = true_vals, aes(x = model, y = tau2, col = 'red')) + 
      ggtitle('tau2 estimates') + 
      theme(legend.position = 'none')
  }
  
  if(sigma2_estimates){
    true_val = data.frame(val = data_list$sigma2)
    df <- data.frame(estimates = stan_out$sigma2)
    
    p_sigma2 <- ggplot(data = df, aes(y = estimates)) + 
      geom_boxplot() + 
      geom_point(data = true_val, aes(y = val, x = 0, col = 'red')) +
      ggtitle('sigma2 estimates') + 
      theme(legend.position = 'none')
  }
  
  # combine rho and tau2
  if(rho_estimates + tau2_estimates + sigma2_estimates > 0){
    # create the hyperparam plot list
    hyperparam_plot = NULL
    if(rho_estimates){
      hyperparam_plot <- append(hyperparam_plot, list(p_rho))
    }
    if(tau2_estimates){
      hyperparam_plot <- append(hyperparam_plot, list(p_tau2))
    }
    if(sigma2_estimates){
      hyperparam_plot <- append(hyperparam_plot, list(p_sigma2))
    }

    # add the plot back to the overall one.
    plot_list <- append(plot_list, list(plot_grid(plotlist = hyperparam_plot, nrow = 1)))
  }
  
  ## compare the true phi values with the estimated phi values
  if(phi_estimates){
    phi_est <- as.data.frame(matrix(0, nrow = N, ncol = length(models)))
    colnames(phi_est) <- models
    for(i in 1:length(models)){
      ind = grep(sprintf('^phi\\[[0-9]{1,2},%s\\]', i), rownames(stan_summary))
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
      geom_smooth(method='lm', formula = y ~ x) + 
      geom_abline(slope = 1, intercept = 0, col = 'red') + 
      facet_wrap(~model) + 
      xlab('true phi values') + 
      ylab('median est') + 
      ggtitle('phi estimates')
    
    plot_list <- append(plot_list, list(p_phi))
  }
  
  ## compare the true u values with the estimated u values
  if(u_estimates){
    u_est <- as.data.frame(matrix(0, nrow = N, ncol = length(models)))
    colnames(u_est) <- models
    for(i in 1:length(models)){
      ind = grep(sprintf('^u\\[[0-9]{1,2},%s\\]', i), rownames(stan_summary))
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
      geom_smooth(method='lm', formula = y ~ x) + 
      geom_abline(slope = 1, intercept = 0, col = 'red') + 
      facet_wrap(~model) + 
      xlab('true u values') + 
      ylab('median est') + 
      ggtitle('u estimates')
    
    plot_list <- append(plot_list, list(p_u))
  }
  
  ## compare the true outcomes with the estimated outcomes (compared to just using one model in the outcomes)
  if(y_estimates){
    # first plot
    ind = grep('y_exp', rownames(stan_summary))
    data_list$data$y_predicted <- stan_summary[ind,'50%']
    p_y <- ggplot(data_list$data, aes(x = y, y_predicted)) + 
      geom_point() +
      geom_smooth(method='lm', formula = y ~ x) + 
      geom_abline(slope = 1, intercept = 0, col = 'red') + 
      xlab('observed y') + 
      ylab('median est') + 
      ggtitle('y estimation')
    
    if('y2' %in% colnames(data_list$data) & 'y3' %in% colnames(data_list$data)){
      
      p_y2 <- ggplot(data_list$data, aes(x = y2, y_predicted)) + 
        geom_point() +
        geom_smooth(method='lm', formula = y ~ x) + 
        geom_abline(slope = 1, intercept = 0, col = 'red') + 
        xlab('observed y') + 
        ylab('median est') + 
        ggtitle('y OOS est')
      
      p_y3 <- ggplot(data_list$data, aes(x = y3, y_predicted)) + 
        geom_point() +
        geom_smooth(method='lm', formula = y ~ x) + 
        geom_abline(slope = 1, intercept = 0, col = 'red') + 
        xlab('observed y') + 
        ylab('median est') + 
        ggtitle('y OOS est')
      
      plot_list <- append(plot_list, list(plot_grid(p_y, p_y2, p_y3, nrow = 1)))
    }else{
      plot_list <- append(plot_list, list(p_y))
    }
    
  }
  
  ## Combine all the plots!
  full_plot <- plot_grid(plotlist = plot_list,
                         ncol = 1)
  
  return(full_plot)
}


### Plot the results for multiple simulations. Currently just plots the spatial parameters, but that will likely be adjusted.
# sim_lst: simulation list from multiple_sims, containing a list with "data_list" and "stan_fit" for each simulation run.
# models: models used in the fitting.
# ncol: number of columns in output plot.
plot_multiple_sims <- function(sim_lst, models, ncol = 2, ESS = F, likelihoods = F, rho_estimates = T, tau2_estimates = T, sigma2_estimates = F, phi_estimates = F, u_estimates = F, y_estimates = F){
  
  # initialize the plot list
  plot_list <- list()
  
  # cycle through each simulation
  for(i in 1:length(sim_lst)){
    # create the param estimate plots
    pp <- process_results(sim_lst[[i]]$data_list, models, sim_lst[[i]]$stan_fit, rho_estimates = rho_estimates, tau2_estimates = tau2_estimates, sigma2_estimates = sigma2_estimates, ESS = ESS, likelihoods = likelihoods, phi_estimates = phi_estimates, u_estimates = u_estimates, y_estimates = y_estimates)
    
    # store the plot
    plot_list[[i]] <- pp
  }
  
  # make the final plot
  final_plot <- plot_grid(plotlist = plot_list, ncol = ncol)
  
  return(final_plot)
}


### Make the panel plot that extracts different parameters from multiply simulations and plots them all together. This calls plot_multiple_sims
# res_lst: results list from running multiple sims
make_panel_plot <- function(res_lst){
  # plot a single simulation results
  pp <- process_results(res_lst$sim_list[[1]]$data_list, res_lst$models, res_lst$sim_list[[1]]$stan_fit, tau2_estimates = T, likelihoods = T, sigma2_estimates = T)
  
  # plot multiple simulation parameter estimates, u estimates, and y estimates
  tt1 <- plot_multiple_sims(res_lst$sim_list, res_lst$models, sigma2_estimates = T, ncol = 1)
  tt2 <- plot_multiple_sims(res_lst$sim_list, res_lst$models, u_estimates = T, rho_estimates = F, tau2_estimates = F, ncol = 1)
  tt3 <- plot_multiple_sims(res_lst$sim_list, res_lst$models, y_estimates = T, u_estimates = F, rho_estimates = F, tau2_estimates = F, ncol = 1)
  
  # combine them all
  full_plot <- plot_grid(pp, tt1, tt2, tt3, nrow = 1, rel_widths = c(3,3,3,2))#, labels = c('One Sim','Params','Weights','Y'))

  return(full_plot)
}
