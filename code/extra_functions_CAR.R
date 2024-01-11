
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
simulate_models <- function(data, models, means, variances){
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
  
### Generate the precision matrix based. Can generate an ICAR, Cressie, or Leroux precision matrix
# type: "ICAR" or "Leroux", for the type of precision matrix
# rho: the spatial correlation parameter, ranging from 0 to 1
# tau2: the variance parameter
generate_precision_mat <- function(W, type, tau2, rho){
  D = diag(rowSums(W))
  if(type == 'ICAR'){
    Q = (D - rho*W)
  }else if(type == 'Leroux'){
    Q = (rho*(D - W) + (1-rho)*diag(rep(1,nrow(D))))
  }else{
    stop('please input a proper precision matrix type')
  }
  return(Q)
}

### Sampling a multivariate normal from a precision matrix
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
# precision_type: "ICAR" or "Leroux", determining the CAR precision matrix.
# tau2: the CAR variance parameter
# rho: the CAR spatial correlation parameter
simulate_data <- function(data, adjacency, models = c('acs','pep','worldpop'), scale_down = 1, pivot = -1, precision_type = 'ICAR', tau2 = 1, rho = 0.3){
  
  # scale down the data size
  data[,models] <- data[,models]/scale_down
  
  # create the precision matrix
  Q <- generate_precision_mat(W = adjacency, type = precision_type, tau2 = tau2, rho = rho)
  
  if(pivot == -1){
    # sample from MVN based on the precision matrix
    phi_true <- sample_MVN_from_precision(n = length(models), Q = Q)
  }else{
    stop('havent coded for pivot data generation')
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
  
  return(list(data = data, phi_true = phi_true, u_true = u_true))
}









  