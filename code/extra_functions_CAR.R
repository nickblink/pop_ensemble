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
  # set.seed(seed)
  
  # setting the seed so that u and phi are the same each time.
  set.seed(1)
  
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


### Get quantiles from a stan fit
# stan_fit: An stan fit object.
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
                   seed = seed,
                   show_messages = F,
                   verbose = F)

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
# CV_blocks: Number of blocks for running cross-validation. If null, only running the full model on the data.
# seed_start: Value to shift the seed starting by.
# return_quantiles: If true, only return quantiles of simulation results. If false, return the full simulation results.
## Optional arguments
# precision_type: Leroux or Cressie.
# tau2: scalar or vector of CAR variance parameter.
# rho: scalar or vector of spatial correlation parameter.
# n.sample: number of stan chain samples.
# burnin: length of burnin period for stan.
# sigma2: sigma2 value of y distribution.
multiple_sims <- function(raw_data, models, means, variances, family = 'poisson', N_sims = 10, stan_path = "code/CAR_leroux_sparse_poisson.stan", init_vals = '0', family_name_check = T, use_softmax = F, CV_blocks = NULL, seed_start = 0, return_quantiles = T, ...){
  
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
    # set the seed value
    seed_val = i + seed_start
    
    # simulate input model values
    data <- simulate_models(data = raw_data$data, adjacency = raw_data$adjacency, models = models, seed = seed_val, means = means, variances = variances, ...)
    
    # simulate y values from input models
    data_lst <- simulate_y(data, raw_data$adjacency, models = models, seed = seed_val, family = family, use_softmax = use_softmax, ...)
    
    # update the initialization to start at the true values.
    if(tolower(init_vals) == 't' | tolower(init_vals) == 'truth'){
      init_list = list(phi = as.matrix(data_lst$phi_true[,-ncol(data_lst$phi_true)]))
      init_vals <- function(){init_list}
    }
    
    if(!is.null(CV_blocks)){
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
        
        # run the model!
        tmp_stan_fit <- run_stan_CAR(block_data, data_lst$adjacency, models = models, seed = seed_val, stan_m = m, use_softmax = use_softmax, ...)
        # store the outcome values:
        tmp_y_pred <- t(extract(tmp_stan_fit, pars = 'y_pred')[[1]])
        for(i_block in ind){
          block_y_pred[[i_block]] <- tmp_y_pred[i_block,]
        }
      }
      
      CV_pred <- do.call('rbind', block_y_pred)
    }
   
    # fit the Bayesian model
    stan_fit <- run_stan_CAR(data_lst$data, data_lst$adjacency, models = models, seed = seed_val, stan_m = m, use_softmax = use_softmax, ...)
    
    # store results
    if(return_quantiles){
      stan_quants <- get_stan_quantiles(stan_fit)
      tmp_lst <- list(data_list = data_lst, stan_fit = stan_quants)
    }else{
      tmp_lst <- list(data_list = data_lst, stan_fit = stan_fit)
    }
    
    
    if(!is.null(CV_blocks)){
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


### Function to generate a list of results across simulations 
# folder: name of folder containing results files.
# root: directory where this folder is located.
# debug_mode: pauses the code right after loading results.
generate_metrics_list <- function(folder = NULL, root = NULL, debug_mode = F){
  
  # get the root if necessary
  if(is.null(root)){
    if(file.exists('C:/Users/Admin-Dell')){
      root_dir = 'C:/Users/Admin-Dell'
    }else{
      root_dir = 'C:/Users/nickl'
    }
    root <- sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Results',
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
      
      # extract the medians.
      medians <- tmp$stan_fit['0.5',]
      
      # printing for error checking.
      # print(sprintf('seed start = %s: sum(u) = %s: sum(y) = %s',  res_lst[[i]]$arguments$seed_start, sum(tmp$data_list$u_true), sum(tmp$data_list$data$y)))
      
      # pull out the y predictions
      ind_y_pred <- grep('y_pred', names(medians))
      median_y_pred <- medians[ind_y_pred]
      y_pred_05 <- tmp$stan_fit['0.05',ind_y_pred]
      y_pred_95 <- tmp$stan_fit['0.95',ind_y_pred]
      
      # pull out the y values
      y <- tmp$data_list$data$y
      y2 <- tmp$data_list$data$y2
      
      # get true u and phi values
      phi_true <- tmp$data_list$phi_true
      u_true <- tmp$data_list$u_true
      phi_true_flat <- as.vector(as.matrix(phi_true[,-ncol(phi_true)]))
      u_true_flat <- as.vector(as.matrix(u_true[,-ncol(u_true)]))
      
      # get estimates u and phi values
      ind_phi <- grep('^phi\\[', colnames(tmp$stan_fit))
      ind_u <- grep('^u\\[', colnames(tmp$stan_fit))
      phi_est_05 <- tmp$stan_fit['0.05', ind_phi] 
      phi_est_95 <- tmp$stan_fit['0.95', ind_phi]
      median_phi <- medians[ind_phi]
      u_est_05 <- tmp$stan_fit['0.05', ind_u]
      u_est_95 <- tmp$stan_fit['0.95', ind_u]
      median_u_mat <- medians[ind_u] %>% 
        vec_to_mat(., n_models = 3)

      metrics_lst[[iter]] <- list(RMSE_train = sqrt(mean((median_y_pred - y)^2)),
                                  RMSE_general = sqrt(mean((median_y_pred - y2)^2)),
                                  RMSE_CV = sqrt(mean((tmp$CV_pred['0.5',] - y)^2)),
                                  CP_90_train = (y >= y_pred_05 & y <= y_pred_95),
                                  CP_90_general = (y2 >= y_pred_05 & y2 <= y_pred_95),
                                  CP_90_CV = (y >= tmp$CV_pred['0.05',] & y <= tmp$CV_pred['0.95',]),
                                  CP_90_phi = (phi_true_flat >= phi_est_05 & phi_true_flat <= phi_est_95),
                                  CP_90_u = (u_true_flat >= u_est_05 & u_true_flat <= u_est_95),
                                  rank_equal = sapply(1:nrow(median_u_mat), function(xx){
                                    res <- all(rank(median_u_mat[xx,]) == rank(u_true[xx,-ncol(u_true)]))
                                    res
                                  })
      )
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
# RMSE_CP_values: whether to include the RMSE values and coverage probabilities.
process_results <- function(data_list, stan_fit, stan_fit_quantiles = F, models = NULL, CV_pred = NULL, ESS = T, likelihoods = T, rho_estimates = T, tau2_estimates = T, sigma2_estimates = F, phi_estimates = T, u_estimates = T, y_estimates = T, RMSE_CP_values = T){
  
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
      ggtitle('y estimation')
    
    if('y2' %in% colnames(data_list$data)){
      
      p_y2 <- ggplot(data_list$data, aes(x = y2, y_predicted)) + 
        geom_point() +
        geom_smooth(method='lm', formula = y ~ x) + 
        geom_abline(slope = 1, intercept = 0, col = 'red') + 
        xlab('observed y') + 
        ylab('median est') + 
        ggtitle('y OOS est')
      
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
          ylab('CV est') + 
          ggtitle('y CV est')
        
        plot_list <- append(plot_list, list(plot_grid(p_y, p_y2, p_y3, nrow = 1)))
      }else if('y3' %in% colnames(data_list$data)){
        p_y3 <- ggplot(data_list$data, aes(x = y3, y_predicted)) + 
          geom_point() +
          geom_smooth(method='lm', formula = y ~ x) + 
          geom_abline(slope = 1, intercept = 0, col = 'red') + 
          xlab('observed y') + 
          ylab('median est') + 
          ggtitle('y OOS est')
        
        plot_list <- append(plot_list, list(plot_grid(p_y, p_y2, p_y3, nrow = 1)))
      }else{
        plot_list <- append(plot_list, list(plot_grid(p_y, p_y2, nrow = 1)))
      }
      
    }else{
      plot_list <- append(plot_list, list(p_y))
    }
    
  }
  
  if(RMSE_CP_values){
    # initialize data frame:
    RMSE_CP_df <- data.frame(dataset = as.character(NA), RMSE = as.numeric(NA), CP.95 = as.numeric(NA))
    
    # get the median, lower, and upper predictions from this set.
    ind = grep('y_pred', rownames(stan_summary))
    median_y_pred <- stan_summary[ind,'50%']
    lower_y_pred <- stan_summary[ind,'2.5%']
    upper_y_pred <- stan_summary[ind,'97.5%']
    
    # get the training set RMSE and CP
    y = data_list$data$y
    RMSE_train = sqrt(mean((median_y_pred - y)^2))
    CP_train = mean(y >= lower_y_pred & y <= upper_y_pred)
    RMSE_CP_df[1,1] <- 'train'
    RMSE_CP_df[1,2:3] <- c(RMSE_train, CP_train)
    
    # get the generalization RMSE in a new set.
    y2 = data_list$data$y2
    RMSE_general = sqrt(mean((median_y_pred - y2)^2))
    CP_general = mean(y2 >= lower_y_pred & y2 <= upper_y_pred)
    RMSE_CP_df[2,1] <- 'generalize'
    RMSE_CP_df[2,2:3] <- c(RMSE_general, CP_general)
    
    # get the CV RMSE, if CV was run.
    if(!is.null(CV_pred)){
      CV_quants = t(apply(CV_pred, 1, function(xx){
        quantile(xx, probs = c(0.025, 0.5, 0.975))
      }))
      
      RMSE_CV = sqrt(mean((CV_quants[,2] - y)^2))
      CP_CV = mean(y >= CV_quants[,1] & y <= CV_quants[,3])
      RMSE_CP_df[3,1] <- 'train-CV'
      RMSE_CP_df[3,2:3] <- c(RMSE_CV, CP_CV)
    }
    
    # rounding for display
    RMSE_CP_df[,2:3] <- round(RMSE_CP_df[,2:3], 3)
    p_RMSE_CP <- gridExtra::tableGrob(RMSE_CP_df)
    
    plot_list <- append(plot_list, list(p_RMSE_CP))
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
plot_multiple_sims_estimates <- function(sim_lst, models, ncol = 2, ESS = F, likelihoods = F, rho_estimates = T, tau2_estimates = T, sigma2_estimates = F, phi_estimates = F, u_estimates = F, y_estimates = F){
  
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
plot_metrics <- function(input_lst, single_sim_res = NULL){
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
                          ESS = F, likelihoods = F, rho_estimates = F, tau2_estimates = F, sigma2_estimates = F, phi_estimates = F, u_estimates = F, RMSE_CP_values = F, 
                          y_estimates = T)
    
    plot_list <- append(plot_list, list(p_yfit))
  }
  
  ## full plot
  full_plot <- cowplot::plot_grid(plotlist = plot_list,
                                  ncol = 1)
  
  return(full_plot)
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
