library(dplyr)
library(rstan)

### Print the most recent files in a directory (makes post processing easier).
# l: Number of files to show.
recent_files <- function(l = 5){
  dd <- dir()
  file_df <- data.frame(file = dd,
                        time = unlist(lapply(dd, function(xx) {return(as.character(file.info(xx)$ctime))}))) %>% 
    mutate(time = as.Date(time)) %>%
    arrange(desc(time))
  print(head(file_df, l))
}


### Compare the parameters between two simulation runs.
# folder1: First folder to compare.
# folder2: Second folder to compare.
compare_parameters <- function(folder1, folder2){
  # load the parameters
  load(dir(folder1, full.names = T)[1])
  params1 <- params
  load(dir(folder2, full.names = T)[1])
  params2 <- params
  rm(params)
  
  sd_12 <- setdiff(names(params1), names(params2))
  sd_21 <- setdiff(names(params2), names(params1))
  if(length(sd_12) > 0){
    print(sprintf('%s is in folder 1 parameters but not folder 2', sd_12))
  }
  if(length(sd_21) > 0){
    print(sprintf('%s is in folder 2 parameters but not folder 1', sd_21))
  }
  
  params_to_compare <- setdiff(intersect(names(params1), names(params2)), c('raw_data', 'output_path'))
  for(pp in params_to_compare){
    if(!identical(params1[[pp]], params2[[pp]])){
      print(sprintf('%s is different in folder 1 (%s) than in folder 2 (%s)', pp, params1[[pp]], params2[[pp]]))
    }
  }
}


### Print out the median and confidence intervals for a given vector.
# vec: A vector or 1-column data frame.
# prob_vec: The probability values to take the quantiles at.
make_medianCI_string <- function(vec, prob_vec = c(0.05, 0.95)){
  # if dataframe, take first column
  if(is.data.frame(vec)){
    if(ncol(vec) > 1){
      stop('please only put in a vector or single column df to make_medianCI_string')
    }
    vec <- vec[,1]
  }
  median_val <- median(vec)
  quantiles <- quantile(vec, probs = prob_vec)
  result_string <- sprintf("%.2f (%.2f, %.2f)", 
                           median_val, quantiles[1], quantiles[2])
  return(result_string)
}

### Convert vector of outputs to a matrix
# vec: Vector where values are arranged by column.
# n_models: number of models to split the vector by (e.g. number of rows of output matrix.)
vec_to_mat <- function(vec, n_models = 3){
  n <- length(vec)
  n_rows <- n/n_models
  
  # checking that the ordering is correcting (that it's by column, NOT row).
  tmp <- names(vec)[2]
  if(substr(tmp, nchar(tmp) - 1, nchar(tmp) - 1) != '1'){
    stop('ordering should be by column, but index doesnt match that')
  }
  
  # checking that the number of models is correct.
  tmp <- names(vec)[n]
  if(as.numeric(substr(tmp, nchar(tmp) - 1, nchar(tmp) - 1)) != n_models){
    stop('incorrect number of models')
  }
  
  mat <- matrix(NA, nrow = n_rows, ncol = n_models)
  for(j in 1:n_models){
    ind <- ((j - 1)*n_rows + 1):(j*n_rows)
    mat[,j] <- vec[ind]
  }
  
  return(mat)
}


### Make an adjacency list out of an adjacency matrix.
# adj_mat: An adjacency matrix of 0's and 1's
make_adjacency_list <- function(adj_mat){
  adj_list <- list()
  
  # cycle through each row
  for(i in 1:nrow(adj_mat)){
    adj_list[[rownames(adj_mat)[i]]] <- which(adj_mat[i,] == 1)
  }
  
  return(adj_list)
}


### Make K non-neighboring folds of data from an adjacency matrix.
# adj_mat: An adjacency matrix of 0's and 1's
# K: Number of folds. If equal to 0, perform LOOCV. Otherwise, this should be at least 5.
make_data_folds <- function(adj_mat, K){
  
  # LOOCV
  if(K == 0){
    groups <- 1:nrow(adj_mat)
  }else if(K < 5){
    stop('too few folds.')
  }else{
    # make the adjacency list
    adj_list <- make_adjacency_list(adj_mat)
    
    # make the map
    groups <- tmaptools::map_coloring(adj_list, ncols = K)
  }
  
  print(table(groups))
  return(groups)
}


### Subset full dataset and adjacency by a specific state
subset_data_by_state <- function(data, adjacency, state, abbrev = NULL){
  df_states <- read.csv('data/state_abbreviations.csv')
  if(!is.null(abbrev)){
    ind <- which(df_states$Code == toupper(abbrev))
    state <- df_states$State[ind]
  }
  
  # get the indices corresponding to the state
  ind1 <- grep(state, data$NAME)
  ind2 <- grep(state, rownames(adjacency))
  
  # check that the indices match between data and adjacency matrix
  if(!identical(ind1, ind2)){
    stop('indices of data names and adjacency names do not match')
  }
  
  # check abbreviated adjacency column names
  # if(!is.null(abbrev)){
  #   ind3 <- grep(abbrev, colnames(adjacency))
  #   browser()
  #   if(!identical(ind1, ind3)){
  #     stop('indices of data names and adjacency columnss do not match')
  #   }
  # }
  
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
simulate_X <- function(data, models, means, variances, seed = 10, adjacency = NULL, M_CAR_rho = NULL, M_CAR_tau2 = NULL, M_BYM_variance = F, precision_type = 'Leroux', M_MVN_alpha = NULL, ...){
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
# theta: Overdispersion parameter for the negative binomial distribution.
# use_softmax: If true, softmax is used. If false, the phi values are directly used as weights (centered at 1/M).
simulate_phi_u_y <- function(data, adjacency, models = c('M1','M2','M3'), scale_down = 1, pivot = -1, precision_type = 'Leroux', tau2 = 1, rho = 0.3, seed = 10, cholesky = T, family = 'poisson', sigma2 = NULL, theta_true = NULL, use_softmax = NULL, num_y_samples = 1, use_pivot = F, resample_phi = T, ...){
  # set seed for reproducability 
  # set.seed(seed)
  
  if(resample_phi){
    # setting a seed based on input
    set.seed(seed)
  }else{
    # setting the seed so that u and phi are the same each time.
    set.seed(1)
  }
  
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
  
  # setting the seed so that y is now sampled differently.
  set.seed(seed)
  
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
    }else if(tolower(family) %in% c('nb','negbin','negative_binomial')){
      data$y <- MASS::rnegbin(n = nrow(data), mu = data$y_expected, theta = theta_true)
      
      # sample repeated y instances
      if(i >= 2){
        data[,paste0('y',i)] <- MASS::rnegbin(n = nrow(data), mu = data$y_expected, theta = theta_true)
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
  }else if(!is.null(theta_true)){
    res_list[['theta_true']] <- theta_true
  }
  
  return(res_list)
}


### Get quantiles from a stan fit
# stan_fit: A stan fit object.
# quantiles: what quantiles to return.
get_stan_quantiles <- function(stan_fit, quantiles = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975)){
  # convert to a data frame.
  tmp <- as.matrix(stan_fit)
  
  # get the quantiles.
  quants <- apply(tmp, 2, function(x) {quantile(x, probs = quantiles)})
  
  # rename the rows to be in [0,1]
  rownames(quants) <- as.numeric(gsub('\\%','',rownames(quants)))/100
  
  # return!
  return(quants)
}


### Get the MAP from the stan posterior distribution ###
# stan_fit: A stan fit object.
get_stan_MAP <- function(stan_fit, inc_warmup = T){
  stan_array <- as.matrix(stan_fit)
  max_ind <- which.max(stan_array[,'lp__'])
  max_vec <- stan_array[max_ind,]
  # u_max <- extract(stan_fit, pars = 'u')[[1]][max_ind,,]
  return(max_vec)
}


### prep the data for fitting the stan model
# data: the observed data with the outcome "y" and the covariates as models. 
# W: the adjacency matrix.
# models: a vector of the models to use in the ensemble.
# sigma2_prior_shape: Shape of the gamma distribution prior.
# sigma2_prior_rate: rate of the gamma distribution prior.
prep_stan_data_leroux_sparse <- function(data, W, models, outcome = 'y', use_softmax = F, use_pivot = F, use_normal = T, sigma2_prior_shape = 1, sigma2_prior_rate = 10, theta_prior_shape = NULL, theta_prior_rate = NULL, tau2_prior_shape = 1, tau2_prior_rate = 1, fixed_rho = - 1, fixed_tau2 = -1, theta_gamma_prior = 0, alpha_variance_prior = NULL, family = NULL, rho = NULL, tau2 = NULL, Z = NULL, ...){
  # checking columns
  if(!(outcome %in% colnames(data))){
    stop(sprintf('need outcome %s as a column in data', outcome))
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
  y = data[,outcome]
  
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
    ind_miss = as.array(ind_miss), # indices of missing y points
    ind_obs = ind_obs,
    X = X, # design matrix
    y_obs = y_obs, # outcome variable 
    W = W,
    W_n = W_n,
    I = diag(1.0, N),
    lambda = lambda,
    use_softmax = as.integer(use_softmax),
    use_pivot = as.integer(use_pivot),
    theta_gamma_prior = theta_gamma_prior,
    tau2_prior_shape = tau2_prior_shape,
    tau2_prior_rate = tau2_prior_rate,
    fixed_rho = fixed_rho,
    fixed_tau2 = fixed_tau2)
  
  if(tolower(family) == 'normal'){
    stan_data <- c(stan_data,
                   list(sigma2_prior_shape = sigma2_prior_shape,
                        sigma2_prior_rate = sigma2_prior_rate))
  }else if(tolower(family) == 'negbin'){
    if(!is.null(theta_prior_shape)){
      if(!is.na(theta_prior_shape)){
        stan_data <- c(stan_data,
                       list(theta_prior_shape = theta_prior_shape,
                            theta_prior_rate = theta_prior_rate))
      }
    }
    
  }else{
    stop('please input a proper family')
  }
  
  if(!is.null(alpha_variance_prior)){
    stan_data <- c(stan_data, list(alpha_variance_prior = alpha_variance_prior))
  }
  
  if(!is.null(Z)){
    stan_data <- c(stan_data, list(num_vars = ncol(Z), Z = Z))
  }
  
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
run_stan_CAR <- function(data, adjacency, models = c('M1','M2','M3'), precision_type = 'Leroux', n.sample = 10000, burnin = 5000, seed = 10, stan_m = NULL, stan_path = NULL, tau2 = NULL, use_softmax = NULL, use_normal = T, use_pivot = F, init_vals = '0', family = family, chains_cores = 1, Z = NULL, ...){

  # error checking for precision matrix type.
  if(precision_type != 'Leroux'){stop('only have Leroux precision coded')}
  
  # checking about softmax and tau2 values.
  if(!is.null(tau2)){
    if(!use_softmax & tau2 > 0.1){
      print('Are you sure you want tau2 so high?')
    }
  }
  
  # prep the data.
  stan_data <- prep_stan_data_leroux_sparse(data, adjacency, models, use_softmax = use_softmax, use_normal = use_normal, use_pivot = use_pivot, family = family, Z = Z, ...)
  
  # create the stan model if not done already.
  if(is.null(stan_m)){
    stan_m <- rstan::stan_model(stan_path)
  }

  # fit the stan model.
  stan_fit <- rstan::sampling(object = stan_m,
                   data = stan_data, 
                   iter = n.sample, 
                   warmup = burnin,
                   init = init_vals,
                   chains = chains_cores, 
                   cores = chains_cores,
                   seed = seed, 
                   thin = 10,
                   show_messages = F,
                   verbose = F)

  # return the results!
  return(list(stan_fit, stan_data))
}


### Run multiple simulation runs. This calls the data creation functions and run_stan_CAR.
# raw_data: data list containing "data" and "adjacency".
# models: list of input models to put in function.
# stan_path: path to stan code.
# means: means of the input models' data creation.
# variances: variances of the input models' data creation.
# family: family of y distribution for simulation.
# CV_blocks: Number of blocks for running cross-validation. If -1, only running the full model on the data.
# seed_start: Value to shift the seed starting by.
# return_quantiles: If true, only return quantiles of simulation results. If false, return the full simulation results.
## Optional arguments
# precision_type: Leroux or Cressie.
# tau2: scalar or vector of CAR variance parameter.
# rho: scalar or vector of spatial correlation parameter.
# n.sample: number of stan chain samples.
# burnin: length of burnin period for stan.
# sigma2: sigma2 value of y distribution.
multiple_sims <- function(raw_data, models, means, variances, family = 'poisson', N_sims = 10, stan_path = "code/CAR_leroux_sparse_poisson.stan", init_vals = '0', family_name_check = T, use_softmax = F, CV_blocks = -1, seed_start = 0, return_quantiles = T, ...){
  print(getwd())
  rstan_options(auto_write = F)
  
  
  ### Parameter error checks
  {
    # checking stan path exists 
    if(!file.exists(stan_path)){
      print(getwd())
      print(sprintf('stan_path = %s', stan_path))
      stop('stan path does not exist.')
    }
    
    # Checking the name and family match
    if(!grepl(family, stan_path) & family_name_check){
      stop('The code does not have the family name in it')
    }
    
    # checking that the rho in the DGP and the model fit are equal (if the rho is fixed in the model fit)
    if('rho' %in% names(list(...)) & 'fixed_rho' %in% names(list(...))){
      if(list(...)$fixed_rho > 1){
        stop('rho value cannot be greater than 1')
      }
      if(list(...)$fixed_rho > 0 & list(...)$fixed_rho != list(...)$rho){
        print('WARNING: rho in DGP and rho in model not equal')
      }
    }
    
    # checking that the rho in the DGP and the model fit are equal (if the rho is fixed in the model fit)
    if('tau2' %in% names(list(...)) & 'fixed_tau2' %in% names(list(...))){
      if(list(...)$fixed_tau2 > 0 & list(...)$fixed_tau2 != list(...)$tau2){
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
  print('check 1')
  # compile the stan program
  m <- rstan::stan_model(stan_path)
  print('check 1.5')

  # initialize results
  sim_lst <- list()
  
  # cycle through N simulations
  for(i in 1:N_sims){
    # set the seed value
    seed_val = i + seed_start
    
    print('check 2')
    # simulate input model values
    data <- simulate_X(data = raw_data$data, adjacency = raw_data$adjacency, models = models, seed = seed_val, means = means, variances = variances, ...)
    
    print('check 3')
    # simulate y values from input models
    data_lst <- simulate_phi_u_y(data, raw_data$adjacency, models = models, seed = seed_val, family = family, use_softmax = use_softmax, ...)
    
    print('check 3.5')
    # update the initialization to start at the true values.
    if(tolower(init_vals) == 't' | tolower(init_vals) == 'truth' | tolower(init_vals) == 'true'){
      init_list = list(phi = as.matrix(data_lst$phi_true[,-ncol(data_lst$phi_true)]),
                       rho = data_lst$rho,
                       tau2 = data_lst$tau2,
                       sigma2 = data_lst$sigma2,
                       theta = data_lst$theta_true)
      init_vals <- function(){init_list}
    }
    
    print('check 4')
    # run cross-validation.
    if(CV_blocks >= 0){
      # make the folds
      folds = make_data_folds(data_lst$adjacency, K = CV_blocks)
      
      # make the data frame to store block predictions.
      block_y_pred <- list()
      
      # cycle through the blocks.
      nrow(data_lst$data)
      for(k in 1:ifelse(CV_blocks == 0, nrow(data_lst$data), CV_blocks)){
        print(sprintf('------------%s------------', k))
        # pull out the data
        block_data <- data_lst$data
        
        # get indices of current block
        ind <- which(folds == k)
        
        # set y values to 0 of current block
        block_data$y[ind] <- NA
        
        print('check 5')
        # run the model!
        tmp_lst <- run_stan_CAR(block_data, data_lst$adjacency, models = models, seed = seed_val, stan_m = m, use_softmax = use_softmax, init_vals = init_vals, family = family, ...)
        tmp_stan_fit <- tmp_lst[[1]]
        print('check 6')
        
        # store the outcome values:
        tmp_y_pred <- t(extract(tmp_stan_fit, pars = 'y_pred')[[1]])
        for(i_block in ind){
          block_y_pred[[i_block]] <- tmp_y_pred[i_block,]
        }
      }
      
      CV_pred <- do.call('rbind', block_y_pred)
    }
   
    # fit the Bayesian model on the full data
    stan_fit_lst <- run_stan_CAR(data_lst$data, data_lst$adjacency, models = models, seed = seed_val, stan_m = m, use_softmax = use_softmax, init_vals = init_vals, family = family, ...)
    stan_fit <- stan_fit_lst[[1]]

    # get the MAP posterior values.    
    stan_MAP <- get_stan_MAP(stan_fit)

    # store results
    if(return_quantiles){
      stan_quants <- get_stan_quantiles(stan_fit)
      tmp_lst <- list(data_list = data_lst, stan_fit = stan_quants, stan_MAP = stan_MAP, stan_data = stan_fit_lst[[2]])
    }else{
      stan_quants <- get_stan_quantiles(stan_fit)
      tmp_lst <- list(data_list = data_lst, stan_fit = stan_fit, stan_quants = stan_quants, stan_MAP = stan_MAP, stan_data = stan_fit_lst[[2]])
    }
    
    
    if(CV_blocks >= 0){
      if(return_quantiles){
        CV_quants <- apply(CV_pred, 1, function(x) {quantile(x, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975))})
        rownames(CV_quants) <- as.numeric(gsub('\\%','',rownames(CV_quants)))/100
        tmp_lst[['CV_pred']] <- CV_quants
      }else{
        tmp_lst[['CV_pred']] <- CV_pred
      }
    }
    
    # store the seed for checking later.
    tmp_lst[['seed_val']] <- seed_val
    
    # store the list!
    sim_lst[[i]] <- tmp_lst
  }
  
  # store the final set of the results 
  res_lst <- list(sim_list = sim_lst, arguments = arguments, models = models)
  
  # return the results
  return(res_lst)
}


### Fit the real data
# raw_data: data list containing "data" and "adjacency".
# models: list of input models to put in function.
# stan_path: path to stan code.
# family: family of y distribution for simulation.
# CV_blocks: Number of blocks for running cross-validation. If null, only running the full model on the data.
# seed_start: Value to shift the seed starting by.
# return_quantiles: If true, only return quantiles of simulation results. If false, return the full simulation results.
# alpha_variance_prior: Prior on the alpha variance to estimate. If NULL or < 0, alpha not estimated.
# preprocess_scale: 
## Optional arguments
# n.sample: number of stan chain samples.
# burnin: length of burnin period for stan.
fit_model_real <- function(raw_data, models=c('acs','pep','wp'), family = 'poisson', stan_path = "code/CAR_leroux_sparse_poisson.stan", init_vals = '0', family_name_check = T, use_softmax = F, CV_blocks = -1, seed_start = 0, return_quantiles = T, alpha_variance_prior=NULL, preprocess_scale = F, fixed_effects = 'none', ...){
  
  print(sprintf('stan_path = %s', stan_path))
  
  # capture the arguments
  arguments <- match.call()
  
  ### Parameter error checks
  {
    # checking stan path exists 
    if(!file.exists(stan_path)){
      stop('stan path does not exist.')
    }
    
    # Checking the name and family match
    if(!grepl(family, stan_path) & family_name_check){
      stop('The code does not have the family name in it')
    }
    
    if(!is.null(alpha_variance_prior)){
      if(alpha_variance_prior > 0 & use_softmax == F){
        stop('cant use alpha for direct estimate. Only for softmax.')
      }
    }
    
    if(!is.null(fixed_effects)){
      if(!grepl('FE', stan_path)){
        stop('If using fixed effects, need to call a stan script with FE.')
      }
      if(fixed_effects == 'intercept'){
        Z = matrix(1, nrow = nrow(raw_data$data), ncol = 1)
      }else if(fixed_effects == 'pep_density'){
        Z = cbind(rep(1, nrow(raw_data$data)),
                  log(.01 + raw_data$data$pep_density))
      }else if(fixed_effects == 'acs_density'){
        Z = cbind(rep(1, nrow(raw_data$data)),
                  log(.01 + raw_data$data$acs_density))
      }else if(fixed_effects == 'pep_fulldensity'){
        Z = cbind(rep(1, nrow(raw_data$data)),
                  log(.01 + raw_data$data$pep_density_full))
      }else if(fixed_effects == 'acs_fulldensity'){
        Z = cbind(rep(1, nrow(raw_data$data)),
                  log(.01 + raw_data$data$acs_density_full))
      }else if(fixed_effects == 'pep_density_proportion'){
        Z = cbind(rep(1, nrow(raw_data$data)),
                  log(.01 + raw_data$data$pep_density),
                  log(.01 + raw_data$data$pep_AIAN_proportion))
      }else if(fixed_effects == 'pep_fulldensity_proportion'){
        Z = cbind(rep(1, nrow(raw_data$data)),
                  log(.01 + raw_data$data$pep_density_full),
                  log(.01 + raw_data$data$pep_AIAN_proportion))
      }else{
        stop('not a valid fixed effects situation.')
        # load(...)
      }

      # normalize the columns. 
      if(ncol(Z) > 1){
        for(p in 2:ncol(Z)){
          Z[,p] <- (Z[,p] - mean(Z[,p]))/sd(Z[,p])
        }
      }
      
    }
  }
  
  # set use_normal variable.
  if(tolower(family) %in% c('normal','gaussian')){
    use_normal = T
  }else{
    use_normal = F
  }
  
  # pre-scale the input models to be centered at the census.
  if(preprocess_scale){
    census_sum <- sum(raw_data$data$census)
    for(m in models){
      raw_data$data[,m] <- raw_data$data[,m]*census_sum/sum(raw_data$data[,m]) 
    }
  }
  
  print('check 1')
  # compile the stan program
  m <- rstan::stan_model(stan_path)
  print('check 1.5')
  
  # run cross-validation.
  if(CV_blocks >= 0){
    # make the folds
    folds = make_data_folds(raw_data$adjacency, K = CV_blocks)
    
    # make the data frame to store block predictions.
    block_y_pred <- list()
    
    # cycle through the blocks.
    for(k in 1:ifelse(CV_blocks == 0, nrow(raw_data$data), CV_blocks)){
      print(sprintf('------------%s------------', k))
      # pull out the data
      block_data <- raw_data$data
      
      # get indices of current block
      ind <- which(folds == k)
      
      # set y values to 0 of current block
      block_data$census[ind] <- NA
      
      print('check 5')
      # run the model!
      print(sprintf('use softmax = %s', use_softmax))
      tmp_lst <- run_stan_CAR(block_data,
                                   raw_data$adjacency,
                                   models = models,
                                   stan_m = m,
                                   use_softmax = use_softmax, 
                                   init_vals = init_vals, 
                                   family = family,
                                   alpha_variance_prior = alpha_variance_prior, 
                                   Z = Z, ...)
      tmp_stan_fit <- tmp_lst[[1]]
      print('check 6')
      
      # store the outcome values:
      tmp_y_pred <- t(extract(tmp_stan_fit, pars = 'y_pred')[[1]])
      for(i_block in ind){
        block_y_pred[[i_block]] <- tmp_y_pred[i_block,]
      }
    }
    
    CV_pred <- do.call('rbind', block_y_pred)
  }
  
  # fit the Bayesian model on the full data
  stan_lst <- run_stan_CAR(raw_data$data,
                           raw_data$adjacency, 
                           models = models, 
                           stan_m = m, 
                           use_softmax = use_softmax,
                           init_vals = init_vals,
                           family = family,
                           alpha_variance_prior = alpha_variance_prior, 
                           Z = Z, ...)
  stan_fit <- stan_lst[[1]]
  
  # get the MAP posterior values.    
  stan_MAP <- get_stan_MAP(stan_fit)
  
  # store results
  tmp_lst <- list(data_list = raw_data, 
                  stan_fit = stan_fit,  
                  stan_MAP = stan_MAP, 
                  stan_summary = summary(stan_fit), 
                  stan_out = extract(stan_fit), 
                  stan_data = stan_lst[[2]])
  
  if(return_quantiles){
    tmp_lst[['stan_quants']] <- stan_quants
  }
  
  if(CV_blocks >= 0){
    if(return_quantiles){
      CV_quants <- apply(CV_pred, 1, function(x) {quantile(x, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975))})
      rownames(CV_quants) <- as.numeric(gsub('\\%','',rownames(CV_quants)))/100
      tmp_lst[['CV_pred']] <- CV_quants
    }else{
      tmp_lst[['CV_pred']] <- CV_pred
    }
  }
  
  # store the final set of the results 
  res_lst <- list(sim_list = tmp_lst, arguments = arguments, models = models)
  
  # return the results
  return(res_lst)
} # fit model real


### Function to generate a list of results across simulations 
# folder: name of folder containing results files.
# root: directory where this folder is located.
# debug_mode: pauses the code right after loading results.
generate_metrics_list <- function(folder = NULL, root = NULL, hmc_diag = F, debug_mode = F){
  
  # get the root if necessary
  if(is.null(root)){
    if(file.exists('C:/Users/Admin-Dell')){
      root_dir = 'C:/Users/Admin-Dell'
    }else{
      root_dir = 'C:/Users/nickl'
    }
    root <- sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Results/simulated_results',
                    root_dir)
  }
  
  # get the names of the files.
  file_names <- dir(sprintf('%s/%s', root, folder), full.names = T)
  
  # print # of files and return if none found.
  if(length(file_names) == 0){
    stop('no files found')
  }else{
    print(sprintf('found %s results files', length(file_names)))
  }
  
  # initialize values.
  metrics_lst <- list()
  iter <- 0
  
  # cycle through each results file and pull in results.
  for(f in file_names){
    load(f)
    if(iter == 0){
      param_printer = params
      param_printer[['raw_data']] = NULL
      print(param_printer)
    }
    
    if(debug_mode){
      browser()
    }
    
    # cycle through simulations within this file.
    for(i in 1:length(res_lst)){
      iter <- iter + 1
      
      # check length for error.
      if(length(res_lst[[i]]$sim_list) > 1){
        browser()
      }else{
        # extract the simulation results.
        tmp <- res_lst[[i]]$sim_list[[1]]  
      }
      
      # extract the quantiles
      if('stan_quants' %in% names(tmp)){
        stan_quants <- tmp$stan_quants
      }else{
        stan_quants <- tmp$stan_fit
      }
      
      # extract the medians.
      medians <- stan_quants['0.5',]
      rho_medianX <- median(medians[grep('rho_', names(medians))])

      # printing for error checking.
      # print(sprintf('seed start = %s: sum(u) = %s: sum(y) = %s',  res_lst[[i]]$arguments$seed_start, sum(tmp$data_list$u_true), sum(tmp$data_list$data$y)))
      
      # pull out the y predictions
      ind_y_pred <- grep('y_pred', names(medians))
      median_y_pred <- medians[ind_y_pred]
      y_pred_025 <- stan_quants['0.025',ind_y_pred]
      y_pred_05 <- stan_quants['0.05',ind_y_pred]
      y_pred_95 <- stan_quants['0.95',ind_y_pred]
      y_pred_975 <- stan_quants['0.975',ind_y_pred]
      int_widths_95 <- y_pred_975 - y_pred_025
      
      # pull out the y values
      y <- tmp$data_list$data$y
      y2 <- tmp$data_list$data$y2
      
      # get true u and phi values
      phi_true <- tmp$data_list$phi_true
      u_true <- tmp$data_list$u_true
      phi_true_flat <- as.vector(as.matrix(phi_true[,-ncol(phi_true)]))
      u_true_flat <- as.vector(as.matrix(u_true[,-ncol(u_true)]))
      u_MAP_vec <- tmp$stan_MAP[grep('^u\\[', names(tmp$stan_MAP))]
      u_MAP_mat <- matrix(u_MAP_vec, ncol = 3, byrow = F)

      # get estimates u and phi values
      ind_phi <- grep('^phi\\[', colnames(stan_quants))
      ind_u <- grep('^u\\[', colnames(stan_quants))
      phi_est_05 <- stan_quants['0.05', ind_phi] 
      phi_est_95 <- stan_quants['0.95', ind_phi]
      median_phi <- medians[ind_phi]
      u_est_05 <- stan_quants['0.05', ind_u]
      u_est_95 <- stan_quants['0.95', ind_u]
      median_u_mat <- medians[ind_u] %>% 
        vec_to_mat(., n_models = 3)

      metrics_lst[[iter]] <- list(mean_y = mean(y),
                                  median_y = median(y),
                                  RMSE_train = sqrt(mean((median_y_pred - y)^2)),
                                  RMSE_general = sqrt(mean((median_y_pred - y2)^2)),
                                  RMSE_CV = sqrt(mean((tmp$CV_pred['0.5',] - y)^2)),
                                  MAE_train = mean(abs(y - median_y_pred)),
                                  MAE_CV = mean(abs(y - tmp$CV_pred['0.5',])),
                                  MAPE_train = 100*mean(abs(median_y_pred - y)/y),
                                  MAPE_CV = 100*mean(abs(tmp$CV_pred['0.5',] - y)/y),
                                  CP_90_train = (y >= y_pred_05 & y <= y_pred_95),
                                  CP_95_train = (y >= y_pred_025 & y <= y_pred_975),
                                  CP_90_general = (y2 >= y_pred_05 & y2 <= y_pred_95),
                                  CP_90_CV = (y >= tmp$CV_pred['0.05',] & y <= tmp$CV_pred['0.95',]),
                                  CP_95_CV = (y >= tmp$CV_pred['0.025',] & y <= tmp$CV_pred['0.975',]),
                                  CP_90_phi = (phi_true_flat >= phi_est_05 & phi_true_flat <= phi_est_95),
                                  CP_90_u = (u_true_flat >= u_est_05 & u_true_flat <= u_est_95),
                                  median_int_width_train = median(int_widths_95),
                                  median_int_width_CV = median(tmp$CV_pred['0.975',] - tmp$CV_pred['0.025',]),
                                  rank_equal = sapply(1:nrow(median_u_mat), function(xx){
                                    res <- all(rank(median_u_mat[xx,]) == rank(u_true[xx,-ncol(u_true)]))
                                    res
                                  }),
                                  MAP_rank_equal = sapply(1:nrow(u_MAP_mat),  function(xx){
                                    res <- all(rank(u_MAP_mat[xx,]) == rank(u_true[xx,-ncol(u_true)]))
                                    res
                                  }),
                                  median_rhoX = rho_medianX)
      
      if(hmc_diag){
        fit <- tmp$stan_fit
        
        if (!is.null(fit)) {
          sampler_params <- tryCatch(get_sampler_params(fit, inc_warmup = FALSE), error = function(e) NULL)
          
          if (!is.null(sampler_params) && length(sampler_params) > 0) {
            n_divergent <- sum(sapply(sampler_params, function(chain) sum(chain[, "divergent__"])))
            max_treedepth_hit <- any(sapply(sampler_params, function(chain) any(chain[, "treedepth__"] >= 10)))
            bfmi_low <- any(sapply(sampler_params, function(chain) {
              mean_energy <- mean(chain[, "energy__"])
              var_energy <- var(chain[, "energy__"])
              bfmi <- mean_energy^2 / var_energy
              bfmi < 0.3
            }))
            
            metrics_lst[[iter]] <- c(metrics_lst[[iter]],
                                     n_divergent = n_divergent,
                                     max_treedepth_hit = max_treedepth_hit,
                                     bfmi_low = bfmi_low)
          } else {
            metrics_lst[[iter]] <- c(metrics_lst[[iter]],
                                     n_divergent = NA_integer_,
                                     max_treedepth_hit = NA,
                                     bfmi_low = NA)
            
            warning(paste("No sampler params in", f, "- likely due to 0 samples"))
          }
        }
      }
    }
  }
  
  return_lst <- list(metrics_list = metrics_lst, 
                     single_sim = tmp)
  return(return_lst)
}


### process the results. Prints ESS for spatial params and returns a plot of various results for parameter estimates.
# data_list: List containing data used for the fit, the true phi values, and the true u values.
# stan_fit: The result of the stan fit on this data.
# stan_fit_quantiles: True if the stan fit input is just in the quantiles of the estimates.
# models: Vector of models used for fitting.
# ESS: Whether to include the ESS of the parameters.
# likelihoods: whether to include the prior, likelihood, and posterior.
# <XX>_estimates: whether to include the <XX> estimates.
# metrics_values: whether to include the RMSE values and coverage probabilities.
process_results <- function(data_list, stan_fit, stan_fit_quantiles = F, models = NULL, CV_pred = NULL, ESS = T, likelihoods = T, rho_estimates = T, tau2_estimates = T, sigma2_estimates = F, phi_estimates = T, u_estimates = T, y_estimates = T, metrics_values = T){
  
  # checking stan fit type.
  if(any(class(stan_fit) == 'matrix') & stan_fit_quantiles == F){
    print('setting stan_fit_quantiles to true.')
    stan_fit_quantiles = T
  }
  
  # extracting N and models
  N = nrow(data_list$data)
  if(is.null(models)){
    models = grep('^X', colnames(data_list$data), value = T)
  }
  
  plot_list = NULL
  
  ## grab the results
  #stan_summary = summary(stan_fit, pars = c('tau2','rho', 'phi', 'u','y_exp','lp__'))$summary
  if(!stan_fit_quantiles){
    stan_summary = summary(stan_fit)$summary
    stan_out <- extract(stan_fit)
  }else{
    if(ESS | likelihoods){
      stop('cant return ESS or likelihoods when only quantiles of MCMC results are input.')
    }
  }
  
  ## get the convergence parameters
  if(ESS){
    # ESS of tau2 and rho
    ESS_spatial <- data.frame(stan_summary[1:(2*length(models)), c('n_eff', 'Rhat')])
    colnames(ESS_spatial)[1] <- 'ESS'
    p_ESS_spatial <- gridExtra::tableGrob(round(ESS_spatial, 3))
    
    # ESS of phi
    ind = grep('^phi', rownames(stan_summary))
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
    print('rho est')
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
    print('tau2 est')
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
    print('sigma2 est')
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
    print('phi est')
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
    print('u est')
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
    print('y est')
    
    # first plot
    if(stan_fit_quantiles){
      ind = grep('y_pred', colnames(stan_fit))
      data_list$data$y_predicted <- stan_fit['0.5', ind]
    }else{
      ind = grep('y_pred', rownames(stan_summary))
      data_list$data$y_predicted <- stan_summary[ind, '50%']
    }
    
    p_y <- ggplot(data_list$data, aes(x = y, y_predicted)) + 
      geom_point() +
      geom_smooth(method='lm', formula = y ~ x) + 
      geom_abline(slope = 1, intercept = 0, col = 'red') + 
      xlab('observed y') + 
      ylab('median est') + 
      ggtitle('Training Set')
    
    p_y_plots <- list(p_y)
    if('y2' %in% colnames(data_list$data)){
      
      p_y2 <- ggplot(data_list$data, aes(x = y2, y_predicted)) + 
        geom_point() +
        geom_smooth(method='lm', formula = y ~ x) + 
        geom_abline(slope = 1, intercept = 0, col = 'red') + 
        xlab('observed y') + 
        ylab('median est') + 
        ggtitle('y OOS est')
      
      p_y_plots <- append(p_y_plots, list(p_y2))
    }
    # get the CV RMSE, if CV was run.
    if(!is.null(CV_pred)){
      if(!stan_fit_quantiles){
        data_list$data$y_predicted_CV = apply(CV_pred, 1, median)
      }else{
        data_list$data$y_predicted_CV <- CV_pred['0.5',]
      }
      
      p_y3 <- ggplot(data_list$data, aes(x = y, y_predicted_CV)) + 
        geom_point() +
        geom_smooth(method='lm', formula = y ~ x) + 
        geom_abline(slope = 1, intercept = 0, col = 'red') + 
        xlab('observed y') + 
        ylab('median est') + 
        ggtitle('Cross-Validation')
      
      p_y_plots <- append(p_y_plots, list(p_y3))
    }
    # add the plot back to the overall one.
    plot_list <- append(plot_list, list(plot_grid(plotlist = p_y_plots, nrow = 1)))
  }
  
  if(metrics_values){
    print('RMSE est')
    # initialize data frame:
    metrics_df <- data.frame(dataset = as.character(NA), RMSE = as.numeric(NA), CP.95 = as.numeric(NA))
    
    # get the median, lower, and upper predictions from this set.
    ind = grep('y_pred', rownames(stan_summary))
    median_y_pred <- stan_summary[ind,'50%']
    lower_y_pred <- stan_summary[ind,'2.5%']
    upper_y_pred <- stan_summary[ind,'97.5%']
    
    # get the training set RMSE and CP
    y = data_list$data$y
    RMSE_train = sqrt(mean((median_y_pred - y)^2))
    CP_train = mean(y >= lower_y_pred & y <= upper_y_pred)
    metrics_df[1,1] <- 'train'
    metrics_df[1,2:3] <- c(RMSE_train, CP_train)
    
    # get the generalization RMSE in a new set.
    y2 = data_list$data$y2
    RMSE_general = sqrt(mean((median_y_pred - y2)^2))
    CP_general = mean(y2 >= lower_y_pred & y2 <= upper_y_pred)
    metrics_df[2,1] <- 'generalize'
    metrics_df[2,2:3] <- c(RMSE_general, CP_general)
    
    # get the CV RMSE, if CV was run.
    if(!is.null(CV_pred)){
      CV_quants = t(apply(CV_pred, 1, function(xx){
        quantile(xx, probs = c(0.025, 0.5, 0.975))
      }))
      
      RMSE_CV = sqrt(mean((CV_quants[,2] - y)^2))
      CP_CV = mean(y >= CV_quants[,1] & y <= CV_quants[,3])
      metrics_df[3,1] <- 'train-CV'
      metrics_df[3,2:3] <- c(RMSE_CV, CP_CV)
    }
    
    # rounding for display
    metrics_df[,2:3] <- round(metrics_df[,2:3], 3)
    p_metrics <- gridExtra::tableGrob(metrics_df)
    
    plot_list <- append(plot_list, list(p_metrics))
  }
  
  ## Combine all the plots!
  full_plot <- plot_grid(plotlist = plot_list,
                         ncol = 1)
  
  return(full_plot)
}

just_metrics <- function(data_list, stan_fit, stan_fit_quantiles = F, stan_summary, models, CV_pred = NULL){
  res <- plot_real_results(data_list = data_list,
                    stan_fit = stan_fit,
                    stan_summary = stan_summary,
                    models = models,
                    CV_pred = CV_pred,
                    ESS = F,
                    rhats = F,
                    alpha_estimates = F, 
                    tau2_estimates = F, 
                    rho_estimates = F,
                    theta_estimates = F, 
                    phi_estimates = F,
                    pairwise_phi_estimates = F,
                    u_estimates = F,
                    y_estimates = F,
                    beta_estimates = F,
                    metrics_values = T,
                    return_table = T)
  return(res)
}

### plot_real_results. Prints ESS for spatial params and returns a plot of various results for parameter estimates.
# data_list: List containing data used for the fit, the true phi values, and the true u values.
# stan_fit: The result of the stan fit on this data.
# stan_fit_quantiles: True if the stan fit input is just in the quantiles of the estimates.
# models: Vector of models used for fitting.
# ESS: Whether to include the ESS of the parameters.
# likelihoods: whether to include the prior, likelihood, and posterior.
# <XX>_estimates: whether to include the <XX> estimates.
# metrics_values: whether to include the RMSE values and coverage probabilities.
plot_real_results <- function(data_list, stan_fit, stan_fit_quantiles = F, stan_summary = NULL, models = c('acs','pep','wp'), thin = NULL, CV_pred = NULL, ESS = T, rhats = T, rho_estimates = T, tau2_estimates = T, sigma2_estimates = F, theta_estimates = T, alpha_estimates = F, phi_estimates = F, pairwise_phi_estimates = T, u_estimates = T, y_estimates = T, metrics_values = T, beta_estimates = T, beta_varnames = NULL, return_table = F){
  # check the divergences
  check_divergences(stan_fit)
  
  # checking stan fit type.
  if(any(class(stan_fit) == 'matrix') & stan_fit_quantiles == F){
    print('setting stan_fit_quantiles to true.')
    stan_fit_quantiles = T
  }
  
  # extracting N and models
  N = nrow(data_list$data)
  if(is.null(models)){
    models = grep('^X', colnames(data_list$data), value = T)
  }
  
  if(!is.null(thin)){
    stop('code not finished - it appears too slow to do')
    draws_array <- as.array(stan_fit)
    thinned_draws <- draws_array[seq(1, dim(draws_array)[1], by = thin)]
    
    thin_indices <- seq(1, fit@sim$iter, by = 10)
    
    # Modify the number of iterations and the list of samples
    fit@sim$iter <- length(thin_indices)
    fit@sim$samples <- lapply(fit@sim$samples, function(chain) {
      lapply(chain, function(param) param[thin_indices])
    })
    
    # Thin the samples for each chain directly without looping through each parameter
    fit@sim$samples <- lapply(fit@sim$samples, function(chain) {
      # Use vectorized indexing to thin all parameters at once in each chain
      lapply(chain, `[`, thin_indices)
    })
    
    # If your model also involves warmup draws, ensure to adjust the `warmup` parameter if needed
    fit@sim$warmup2 <- fit@sim$warmup / 10
  }
  
  plot_list = NULL
  rel_heights <- c()
  
  ## grab the results
  #stan_summary = summary(stan_fit, pars = c('tau2','rho', 'phi', 'u','y_exp','lp__'))$summary
  if(is.null(stan_summary)){
    stan_summary = summary(stan_fit)$summary
  }else{
    print('using input stan summary')
  }
  stan_out <- extract(stan_fit)
  
  # if(!stan_fit_quantiles){
  #   stan_summary = summary(stan_fit)$summary
  #   stan_out <- extract(stan_fit)
  # }else{
  #   if(ESS | likelihoods){
  #     stop('cant return ESS or likelihoods when only quantiles of MCMC results are input.')
  #   }
  # }
  
  ## get the effective sample size convergence parameters
  if(ESS){
    print('ESS est')
    # ESS of tau2 and rho
    ESS_spatial <- data.frame(stan_summary[1:(2*length(models) + 1), c('n_eff', 'Rhat')])
    colnames(ESS_spatial)[1] <- 'ESS'
    p_ESS_spatial <- gridExtra::tableGrob(round(ESS_spatial, 3))
    
    # ESS of u
    ind <- grep('^u',rownames(stan_summary))
    print(sprintf('median ESS for u is %s and median rhat is %s', round(median(stan_summary[ind, 'n_eff']), 1), round(median(stan_summary[ind, 'Rhat']), 3)))
    
    n_eff_df <- NULL
    for(i in 1:length(models)){
      ind = grep(sprintf('^u\\[[0-9]*,%s\\]', i), rownames(stan_summary))
      n_eff_df <- rbind(n_eff_df, 
                        data.frame(ESS = stan_summary[ind, 'n_eff'],
                                   model = models[i]))
      
    }
    
    # x = data.frame(n_eff = stan_summary[ind, 'n_eff'])
    # p_ESS_phi <- ggplot(data = x, aes(x = n_eff)) + 
    #   geom_density() +
    #   scale_x_continuous(trans='log2') + 
    #   ggtitle('ESS of phis')
    
    p_ESS_u <- ggplot(data = n_eff_df, aes(x = ESS)) + geom_density(aes(colour = model)) + 
      scale_x_continuous(trans='log2') + 
      ggtitle('ESS of u values')
    
    plot_list = append(plot_list, list(plot_grid(p_ESS_spatial, p_ESS_u)))
    rel_heights <- c(rel_heights, 1)
  }
  
  ## get the rhat convergence parameters.
  if(rhats){
    ## u rhats
    ind <- grep('^u',rownames(stan_summary))
    
    n_rhat_df <- NULL
    for(i in 1:length(models)){
      ind = grep(sprintf('^u\\[[0-9]*,%s\\]', i), rownames(stan_summary))
      n_rhat_df <- rbind(n_rhat_df, 
                        data.frame(rhat = stan_summary[ind, 'Rhat'],
                                   model = models[i]))
      
    }
    
    p_rhat_u <- ggplot(data = n_rhat_df, aes(x = rhat)) + geom_density(aes(colour = model)) + 
      scale_x_continuous(trans='log10') +
      scale_y_continuous(trans='sqrt') + 
      ggtitle('rhat of u values')

    ## y rhats
    ind = grep('y_pred', rownames(stan_summary))
    y_pred_rhat <- stan_summary[ind, 'Rhat',drop = F]
    p_rhat_y <- ggplot(data = y_pred_rhat, aes(x = Rhat)) + geom_density() + 
      scale_x_continuous(trans='log10') + 
      scale_y_continuous(trans='sqrt') + 
      ggtitle('rhat of y values')
    
    plot_list = append(plot_list, list(plot_grid(p_rhat_u, p_rhat_y)))
    rel_heights <- c(rel_heights, 1)
  }
  
  ## get the spatial parameter estimates
  if(rho_estimates){
    print('rho est')
    rho <- NULL
    for(i in 1:length(models)){
      rho <- rbind(rho, 
                   data.frame(value = stan_out$rho[,i], 
                              model = models[i]))
    }
    
    # plot the rho param 
    p_rho <- ggplot(data = rho, aes(x = model, y = value)) + 
      geom_boxplot() + 
      ggtitle('rho estimates') + 
      theme(legend.position = 'none')
  }
  
  if(tau2_estimates){
    print('tau2 est')
    tau2 <- NULL
    for(i in 1:length(models)){
      tau2 <- rbind(tau2, 
                    data.frame(value = stan_out$tau2[,i], 
                               model = models[i]))
    }
    
    # plot the tau2 param
    p_tau2 <- ggplot(data = tau2, aes(x = model, y = value)) + 
      geom_boxplot() + 
      ggtitle('tau2 estimates') + 
      theme(legend.position = 'none')
  }
  
  if(sigma2_estimates){
    print('sigma2 est')
    df <- data.frame(estimates = stan_out$sigma2)
    
    p_sigma2 <- ggplot(data = df, aes(y = estimates)) + 
      geom_boxplot() + 
      ggtitle('sigma2 estimates') + 
      theme(legend.position = 'none')
  }
  
  if(theta_estimates){
    print('theta est')
    df <- data.frame(estimates = stan_out$theta)
    
    p_theta <- ggplot(data = df, aes(y = estimates)) + 
      geom_boxplot() + 
      ggtitle('theta estimates') + 
      theme(legend.position = 'none')
  }
  
  if(alpha_estimates){
    print('alpha est')
    alpha <- NULL
    for(i in 1:length(models)){
      alpha <- rbind(alpha, 
                     data.frame(value = stan_out$alpha[,i], 
                                model = models[i]))
    }
    
    # plot the alpha param
    p_alpha <- ggplot(data = alpha, aes(x = model, y = value)) + 
      geom_boxplot() + 
      ggtitle('alpha estimates') + 
      theme(legend.position = 'none') 

    # plot_list <- append(plot_list, list(p_alpha))
    # rel_heights <- c(rel_heights, 1)
  }
  
  # combine rho and tau2
  if(rho_estimates + tau2_estimates + sigma2_estimates + theta_estimates + alpha_estimates > 0){
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
    if(theta_estimates){
      hyperparam_plot <- append(hyperparam_plot, list(p_theta))
    }
    if(alpha_estimates){
      hyperparam_plot <- append(hyperparam_plot, list(p_alpha))
    }
    
    # add the plot back to the overall one.
    plot_list <- append(plot_list, list(plot_grid(plotlist = hyperparam_plot, nrow = 1)))
    rel_heights <- c(rel_heights, 1)
  }

  ## compare the true phi values with the estimated phi values
  if(phi_estimates){
    print('phi est')

    phi_est <- as.data.frame(matrix(0, nrow = N, ncol = length(models)))
    colnames(phi_est) <- models
    for(i in 1:length(models)){
      ind = grep(sprintf('^phi\\[[0-9]*,%s\\]', i), rownames(stan_summary))
      phi_est[,i] <- stan_summary[ind,'50%']
    }
    phi_est$index = 1:N
    
    # convert estimates to long
    phi_est_long <- tidyr::gather(phi_est, key = model, value = phi_median_est, -index)
    
    # plot 'em
    p_phi <- ggplot(phi_est_long, aes(x = phi_median_est)) + 
      geom_density() + 
      facet_wrap(~model) + 
      xlab('estimate') + 
      ylab('density') + 
      theme(axis.text.y = element_blank()) +
      ggtitle('phi estimates')
    
    plot_list <- append(plot_list, list(p_phi))
    rel_heights <- c(rel_heights, 1)
  }
  
  ## pairwise phi estimates
  if(pairwise_phi_estimates){
    phi_est <- as.data.frame(matrix(0, nrow = N, ncol = length(models)))
    colnames(phi_est) <- models
    for(i in 1:length(models)){
      ind = grep(sprintf('^phi\\[[0-9]*,%s\\]', i), rownames(stan_summary))
      phi_est[,i] <- stan_summary[ind,'50%']
    }
    p_pairs_phi <- GGally::ggpairs(phi_est,
                    upper = list(continuous = 'density'),
                    diag = list(continuous = "barDiag"),
                    lower = list(continuous = 'cor'),title = 'phi pairwise')
    
    plot_list <- append(plot_list, list(GGally::ggmatrix_gtable(p_pairs_phi)))
    rel_heights <- c(rel_heights, 2)
  }
  
  ## compare the true u values with the estimated u values
  if(u_estimates){
    print('u est')
    u_est <- as.data.frame(matrix(0, nrow = N, ncol = length(models)))
    colnames(u_est) <- models
    for(i in 1:length(models)){
      ind = grep(sprintf('^u\\[[0-9]*,%s\\]', i), rownames(stan_summary))
      u_est[,i] <- stan_summary[ind,'50%']
    }
    u_est$index = 1:N
    
    # convert estimates to long
    u_est_long <- tidyr::gather(u_est, key = model, value = u_median_est, -index)
    
    # plot 'em
    p_u <- ggplot(u_est_long, aes(x = u_median_est)) + 
      geom_density() + 
      facet_wrap(~model, scales = 'free') + 
      xlab('estimate') + 
      ylab('density') + 
      theme(axis.text.y = element_blank()) +
      ggtitle('u estimates')
    
    plot_list <- append(plot_list, list(p_u))
    rel_heights <- c(rel_heights, 1)
  }
  
  ## compare the observed outcomes with the estimated outcomes 
  if(y_estimates){
    print('y est')
    # first plot
    
    if(stan_fit_quantiles){
      ind = grep('y_pred', colnames(stan_fit))
      data_list$data$y_predicted <- stan_fit['0.5', ind]
    }else{
      ind = grep('y_pred', rownames(stan_summary))
      data_list$data$y_predicted <- stan_summary[ind, '50%']
    }
    
    ind <- sample(nrow(data_list$data), 100)
    
    p_y <- ggplot(data_list$data[ind,], aes(x = census, y_predicted)) + 
      geom_point() +
      geom_smooth(method='lm', formula = y ~ x) + 
      geom_abline(slope = 1, intercept = 0, col = 'red') +
      scale_x_continuous(trans='log2') +
      scale_y_continuous(trans='log2') +
      xlab('observed y') + 
      ylab('median est') + 
      ggtitle('y estimation')
    
    # get the CV RMSE, if CV was run.
    if(!is.null(CV_pred)){
      if(!stan_fit_quantiles){
        data_list$data$y_predicted_CV = apply(CV_pred, 1, median)
      }else{
        data_list$data$y_predicted_CV <- CV_pred['0.5',]
      }
      
      p_y2 <- ggplot(data_list$data[ind,], aes(x = census, y_predicted_CV)) +
        geom_point() +
        geom_smooth(method='lm', formula = y ~ x) + 
        geom_abline(slope = 1, intercept = 0, col = 'red') +
        scale_x_continuous(trans='log2') +
        scale_y_continuous(trans='log2') +
        xlab('observed y') + 
        ylab('CV est') + 
        ggtitle('y CV est')
      
      plot_list <- append(plot_list, list(plot_grid(p_y, p_y2)))
      rel_heights <- c(rel_heights, 1)
    }else{
      plot_list <- append(plot_list, list(p_y))
      rel_heights <- c(rel_heights, 1)
    }
  }

  ## Plot the beta values
  if(beta_estimates){
    print('beta est')
    beta <- NULL
    
    # get number of variables and their names.
    n_vars <- length(stan_out$beta[1,,1])
    if(is.null(beta_varnames)){
      beta_varnames <- paste0('var',1:n_vars)
    }
    
    # cycle through models and variables.
    for(m in 1:length(models)){
      for(ivar in 1:n_vars){
        beta <- rbind(beta, 
                      data.frame(value = stan_out$beta[,ivar,m], 
                                 model = models[m],
                                 var = beta_varnames[ivar]))
      }
    }
    
    # plot the beta param
    p_beta <- ggplot(data = beta, aes(x = var, y = value, fill = model)) + 
      geom_boxplot(position = position_dodge(width = 0.75)) + 
      geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
      ggtitle('Beta Estimates') + 
      theme(legend.position = 'right') + 
      scale_fill_brewer(palette = "Set2")
    
    # update the most recent CP_RMSE plot list and rel heights
    plot_list <- append(plot_list, list(p_beta))
    rel_heights <- c(rel_heights, 1)
  }
  
  ## get the RMSE CP values!
  if(metrics_values){
    print('RMSE est')
    # initialize data frame:
    metrics_df <- data.frame(dataset = as.character(NA), 
                             RMSE = as.numeric(NA), 
                             rRMSE = as.numeric(NA),
                             logRMSE = as.numeric(NA),
                             MAPE = as.numeric(NA),
                             MAE = as.numeric(NA),
                             CP.95 = as.numeric(NA),
                             med_int = as.numeric(NA),
                             mean_p_int = as.numeric(NA))
    
    # get the median, lower, and upper predictions from this set.
    ind = grep('y_pred', rownames(stan_summary))
    median_y_pred <- stan_summary[ind,'50%']
    lower_y_pred <- stan_summary[ind,'2.5%']
    upper_y_pred <- stan_summary[ind,'97.5%']
    int_widths <- upper_y_pred - lower_y_pred
    
    # get the training set y.
    y <- data_list$data$census
    nonzero_ind <- which(y > 0)
    
    # Calcute the RMSE, relative-RMSE, log RMSE, MAPE, and CV.
    RMSE_train <- sqrt(mean((median_y_pred - y)^2))
    rRMSE_train <- sqrt(mean(((median_y_pred[nonzero_ind] - y[nonzero_ind])/y[nonzero_ind])^2))
    logRMSE_train <-  sqrt(mean((log1p(median_y_pred[nonzero_ind]) - log1p(y[nonzero_ind]))^2))
    MAPE_train <- 100*mean(abs((median_y_pred[nonzero_ind] - y[nonzero_ind])/y[nonzero_ind]))
    MAE_train <- mean(abs(y - median_y_pred))
    CP_train <- mean(y >= lower_y_pred & y <= upper_y_pred)
    med_int <- median(int_widths)
    mean_p_int <- mean(int_widths/y)
    
    # store the results (separate rows because of separate character types).
    metrics_df[1,1] <- 'train'
    metrics_df[1,2:9] <- c(RMSE_train, rRMSE_train, logRMSE_train, MAPE_train, MAE_train, CP_train, med_int, mean_p_int)

    # get the CV RMSE, if CV was run.
    if(!is.null(CV_pred)){
      CV_quants = t(apply(CV_pred, 1, function(xx){
        quantile(xx, probs = c(0.025, 0.5, 0.975))
      }))
      int_widths_CV <- CV_quants[,3] - CV_quants[,1]
      
      # Calcute the RMSE, relative-RMSE, log RMSE, MAPE, and CV.
      RMSE_CV = sqrt(mean((CV_quants[,2] - y)^2))
      rRMSE_CV <- sqrt(mean(((CV_quants[nonzero_ind,2] - y[nonzero_ind])/y[nonzero_ind])^2))
      logRMSE_CV <-  sqrt(mean((log1p(CV_quants[nonzero_ind,2]) - log1p(y[nonzero_ind]))^2))
      MAPE_CV <- 100*mean(abs((CV_quants[nonzero_ind,2] - y[nonzero_ind])/y[nonzero_ind]))
      MAE_CV <- mean(abs(CV_quants[,2] - y))
      CP_CV = mean(y >= CV_quants[,1] & y <= CV_quants[,3])
      med_int_CV <- median(int_widths_CV)
      mean_p_int_CV <- mean(int_widths_CV/y)
      
      # store the results.
      metrics_df[2,1] <- 'CV'
      metrics_df[2,2:9] <- c(RMSE_CV, rRMSE_CV, logRMSE_CV, MAPE_CV, MAE_CV, CP_CV, med_int_CV, mean_p_int_CV)
    }

    # rounding for display
    metrics_df[,2:9] <- round(metrics_df[,2:9], 3)
    rownames(metrics_df) <- NULL
    p_metrics <- gridExtra::tableGrob(metrics_df)
    
    plot_list <- append(plot_list, list(p_metrics))
    rel_heights <- c(rel_heights, 1)
  }

  ## Combine all the plots!
  full_plot <- plot_grid(plotlist = plot_list,
                         rel_heights = rel_heights,
                         ncol = 1)
  
  if(return_table){
    return(list(plot = full_plot, 
                metrics = metrics_df))
  }else{
    return(full_plot)  
  }
}


### Plot the results for multiple simulations. Currently just plots the spatial parameters, but that will likely be adjusted.
# sim_lst: simulation list from s, containing a list with "data_list" and "stan_fit" for each simulation run.
# models: models used in the fitting.
# ncol: number of columns in output plot.
plot_s_estimates <- function(sim_lst, models, ncol = 2, ESS = F, likelihoods = F, rho_estimates = T, tau2_estimates = T, sigma2_estimates = F, phi_estimates = F, u_estimates = F, y_estimates = F){
  
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


### Plots the metrics across a set of simulations.
# input_list: Results from the function "generate_metrics_list"
plot_metrics <- function(input_lst, single_sim_res = NULL, include_MAP_rank = F, make_table = T){
  # get correct list of metrics and singe simulation results
  if('metrics_list' %in% names(input_lst)){
    print('extracting metrics AND single simulation results.')
    metrics_lst <- input_lst$metrics_list
    single_sim_res <- input_lst$single_sim
  }else{
    metrics_lst <- input_lst
  }
  
  # number of locations in dataset.
  n_loc <- length(metrics_lst[[1]]$CP_90_train)
  
  # initialize plot list
  plot_list <- NULL
  phi_u_plot <- NULL
  
  ## phi and u coverage plots
  {
    # phi coverage plot
    phi_CP <- data.frame(CP = rowMeans(sapply(metrics_lst, function(xx) xx[['CP_90_phi']])))
    
    # is this averaged over location or simulation?
    p_phi_CP <- ggplot(data = phi_CP, aes(y = CP)) + 
      geom_boxplot() + 
      ggtitle('phi 90% coverage')
    
    phi_u_plot <- append(phi_u_plot, list(p_phi_CP))
    
    # u coverage plot
    u_CP <- data.frame(CP = rowMeans(sapply(metrics_lst, function(xx) xx[['CP_90_u']])))
    
    p_u_CP <- ggplot(data = u_CP, aes(y = CP)) + 
      geom_boxplot() + 
      ggtitle('u 90% coverage')
    
    phi_u_plot <- append(phi_u_plot, list(p_u_CP))
  }
  
  ## u rank plot
  u_rank <- data.frame(rank = rowMeans(sapply(metrics_lst, function(xx) xx[['rank_equal']])))
  
  p_u_rank <- ggplot(data = u_rank, aes(y = rank)) + 
    geom_boxplot() + 
    ggtitle('u-rank scores')
  
  phi_u_plot <- append(phi_u_plot, list(p_u_rank))
  
  if(include_MAP_rank){
    u_MAP_rank <- data.frame(rank = rowMeans(sapply(metrics_lst, function(xx) xx[['MAP_rank_equal']])))
    
    p_MAP_rank <- ggplot(data = u_MAP_rank, aes(y = rank)) + 
      geom_boxplot() + 
      ggtitle('MAP u-rank scores')
    phi_u_plot <- append(phi_u_plot, list(p_MAP_rank))
  }
  
  plot_list <- append(plot_list, list(plot_grid(plotlist = phi_u_plot, nrow = 1)))
  
  
  ## RMSE plot
  {
    RMSE_df <- NULL
    for(rmse in c('RMSE_train', 'RMSE_general', 'RMSE_CV')){
      tmp <- data.frame(val = sapply(metrics_lst, function(xx) xx[[rmse]]),
                        source = rmse)
      RMSE_df <- rbind(RMSE_df, tmp)
    }
    
    RMSE_df$source <- factor(RMSE_df$source, levels = c('RMSE_train', 'RMSE_general', 'RMSE_CV'))
    
    p_RMSE <- ggplot(data = RMSE_df, aes(x = source, y = val)) + 
      geom_boxplot() + 
      ggtitle('RMSE values')
    
    plot_list <- append(plot_list, list(p_RMSE))
  }
  
  ## y coverage plot
  {
    y_CP_df <- NULL
    for(CP in c('CP_90_train', 'CP_90_general', 'CP_90_CV')){
      tmp <- data.frame(val = rowMeans(sapply(metrics_lst, function(xx) xx[[CP]])),
                        source = CP)
      y_CP_df <- rbind(y_CP_df, tmp)
    }
    
    y_CP_df$source <- factor(y_CP_df$source, levels = c('CP_90_train', 'CP_90_general', 'CP_90_CV'))
    
    p_y_CP <- ggplot(data = y_CP_df, aes(x = source, y = val)) + 
      geom_boxplot() + 
      ggtitle('y prediction 90% coverage')
    
    plot_list <- append(plot_list, list(p_y_CP))
  }
  
  # single sim plot
  if(!is.null(single_sim_res)){
    # plot the y predictions for a single simulation.
    y <- single_sim_res$data_list$data$y
    y2 <- single_sim_res$data_list$data$y2
    
    p_yfit <- process_results(single_sim_res$data_list, 
                          CV_pred = single_sim_res$CV_pred, 
                          stan_fit = single_sim_res$stan_fit, 
                          ESS = F, likelihoods = F, rho_estimates = F, tau2_estimates = F, sigma2_estimates = F, phi_estimates = F, u_estimates = F, metrics_values = F, 
                          y_estimates = T)
    
    plot_list <- append(plot_list, list(p_yfit))
  }
  
  if(make_table){
    res_tbl <- data.frame(metric = c('u-rank','u-coverage-90','RMSE','y-coverage-90'), 
                          train = c(make_medianCI_string(u_rank),
                                    make_medianCI_string(u_CP),
                                    make_medianCI_string(RMSE_df %>% filter(source == 'RMSE_train') %>% select(val)),
                                    make_medianCI_string(y_CP_df %>% filter(source == 'CP_90_train') %>% select(val))),
                          CV = c(NA,
                                 NA,
                                 make_medianCI_string(RMSE_df %>% filter(source == 'RMSE_CV') %>% select(val)),
                                 make_medianCI_string(y_CP_df %>% filter(source == 'CP_90_CV') %>% select(val))))
  }
  print(res_tbl)
  print(xtable::xtable(res_tbl))
  
  ## full plot
  full_plot <- cowplot::plot_grid(plotlist = plot_list,
                                  ncol = 1)
  
  return(full_plot)
}


### Make the panel plot that extracts different parameters from multiple simulations and plots them all together. This calls plot_s
# res_lst: results list from running multiple sims
make_panel_plot <- function(res_lst){
  # plot a single simulation results
  pp <- process_results(res_lst$sim_list[[1]]$data_list, res_lst$models, res_lst$sim_list[[1]]$stan_fit, tau2_estimates = T, likelihoods = T, sigma2_estimates = T)
  
  # plot multiple simulation parameter estimates, u estimates, and y estimates
  tt1 <- plot_s(res_lst$sim_list, res_lst$models, sigma2_estimates = T, ncol = 1)
  tt2 <- plot_s(res_lst$sim_list, res_lst$models, u_estimates = T, rho_estimates = F, tau2_estimates = F, ncol = 1)
  tt3 <- plot_s(res_lst$sim_list, res_lst$models, y_estimates = T, u_estimates = F, rho_estimates = F, tau2_estimates = F, ncol = 1)
  
  # combine them all
  full_plot <- plot_grid(pp, tt1, tt2, tt3, nrow = 1, rel_widths = c(3,3,3,2))#, labels = c('One Sim','Params','Weights','Y'))

  return(full_plot)
}


### Make a chloropleth map of a given result.
# df: data frame to plot. 
# val: the name of the column in df to plot.
# counties: The counties sf data frame.
# states: The states sf data frame.
make_chloropleth_plot <- function(df, val, counties = NULL, states = NULL){
  
  # get county shapefiles.
  if(is.null(counties)){
    counties <- counties(cb = TRUE, year = 2020, class = "sf")
  }
  
  # Load state boundaries.
  if(is.null(states)){
    states <- states(cb = TRUE, year = 2020, class = "sf")
  }
  
  # merge in counties data.
  data <- left_join(counties, df, by = 'GEOID')
  
  # get the outcome of interest.
  data$outcome = data[,val,drop=T] %>% as.numeric()
  
  p1 <- ggplot(data = data) +
    geom_sf(aes(fill = outcome), color = NA) +
    geom_sf(data = states, fill = NA, color = "black", size = 0.5) + 
    scale_fill_viridis_c(option = "plasma", na.value = "white", trans = 'log',
                         #breaks = pretty_breaks(n = 5),  # Use pretty breaks
                         labels = function(x) round(x, digits = 0)) +  # Color scale
    theme_minimal() +
    coord_sf(xlim = c(-130, -65), ylim = c(24, 50)) +
    labs(title = sprintf("County %s", val),
         fill = "Value") +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank()) 
  
  return(p1)
}


### Plot the posterior parameter correlation from a stan fit and input parameters.
plot_param_correlation <- function(fit, pars = c('rho_estimated','tau2_estimated','theta')) {
  # Convert to draws_matrix
  draws <- as_draws_matrix(fit)
  
  # Match parameter names
  param_names <- colnames(draws)
  matching_pars <- grep(paste0("^(", paste(pars, collapse = "|"), ")"), param_names, value = TRUE)
  
  if (length(matching_pars) == 0) stop("No matching parameters found in fit.")
  
  # Subset to selected parameters
  draws_subset <- draws[, matching_pars]
  
  # Compute correlation matrix
  cor_mat <- cor(draws_subset)
  
  # Melt to long format
  cor_df <- melt(cor_mat)
  cor_df$value_label <- sprintf("%.2f", cor_df$value)  # format to 2 decimal places
  
  # Plot
  ggplot(cor_df, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = value_label), size = 3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                         limit = c(-1, 1), name = "Correlation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Posterior Correlation Matrix", x = NULL, y = NULL)
}


### Plot divergent pairs in results for a given set of parameters.
plot_divergent_pairs <- function(stan_fit, parameter_names) {
  library(posterior)
  library(GGally)
  library(ggplot2)
  library(rstan)
  
  # Convert to posterior draws
  draws_df <- as_draws_df(stan_fit)
  
  # Extract sampler diagnostics
  sampler_params <- get_sampler_params(stan_fit, inc_warmup = FALSE)
  divergent_vec <- unlist(lapply(sampler_params, function(x) x[, "divergent__"]))
  
  # Safely create a factor with correct labels
  if (all(divergent_vec == 0)) {
    draws_df$divergent <- factor(divergent_vec, levels = 0, labels = "No Divergence")
  } else if (all(divergent_vec == 1)) {
    draws_df$divergent <- factor(divergent_vec, levels = 1, labels = "Divergence")
  } else {
    draws_df$divergent <- factor(divergent_vec, levels = c(0, 1), labels = c("No Divergence", "Divergence"))
  }
  
  # Subset to selected parameters + divergence
  plot_df <- draws_df[, c(parameter_names, "divergent")]
  
  # Clean column names for ggplot
  colnames(plot_df) <- gsub("_estimated", "", colnames(plot_df))
  colnames(plot_df) <- gsub("\\[", "_", gsub("\\]", "", colnames(plot_df)))
  
  # Make plot
  p <- ggpairs(
    plot_df,
    columns = 1:length(parameter_names),
    aes(color = divergent, alpha = 0.4),
    upper = list(continuous = wrap("points", size = 0.5)),
    lower = list(continuous = wrap("points", size = 0.5)),
    diag = list(continuous = wrap("densityDiag", alpha = 0.6))
  ) +
    scale_color_manual(values = c("black", "green")) +
    theme_minimal()
  
  return(p)
}

