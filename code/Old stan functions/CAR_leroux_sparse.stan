functions {
  /**
  * Return the log probability of a proper conditional autoregressive (CAR) prior 
  * with a sparse representation for the adjacency matrix
  *
  * @param phi CAR prior phi for a given model
  * @param tau2 variance parameter for the CAR prior (real)
  * @param alpha Dependence (usually spatial) parameter for the CAR prior (real)
  * @param W_sparse Sparse representation of adjacency matrix (int array)
  * @param n Length of phi (int)
  * @param W_n Number of adjacent pairs (int)
  * @param D_sparse Number of neighbors for each location (vector)
  * @param log_det_Q The precomputed log determinant of the precision matrix
  *
  * @return Log probability density of CAR prior up to additive constant
  */
  real sparse_car_lpdf(vector phi, real tau2, real rho, int[,] W_sparse, 
  vector D_sparse, real log_det_Q, int n, int W_n) {
      row_vector[n] phit_D; // phi' * D
      row_vector[n] phit_W; // phi' * W
    
      phit_D = (phi .* D_sparse)';
      phit_W = rep_row_vector(0, n);
      for (i in 1:W_n) {
        phit_W[W_sparse[i, 1]] = phit_W[W_sparse[i, 1]] + phi[W_sparse[i, 2]];
        phit_W[W_sparse[i, 2]] = phit_W[W_sparse[i, 2]] + phi[W_sparse[i, 1]];
      }
    
      return 0.5 * (log_det_Q - (1/tau2) * (rho * (phit_D * phi) - rho * (phit_W * phi) + (1 - rho)*dot_self(phi)));
  }
}data {
  int<lower=0> M; // number of models
  int<lower=0> N;  // number of observations
  int<lower=0> N_miss; // number of missing y points
  int<lower=0> N_obs; // number of observed y points
  int<lower=0, upper=N> ind_miss[N_miss]; // indices of missing y points
  int<lower=0, upper=N> ind_obs[N_obs]; // indices of observed y points
  matrix[N,M] X; // design matrix of ensemble models
  vector[N_obs] y_obs;  // output
  matrix<lower=0, upper = 1>[N, N] W; //adjacency matrix
  int W_n; // Number of adjacency pairs
  matrix[N,N] I; // Identity matrix
  vector[N] lambda; // the eigenvalues of the D - W - I matrix
  int<lower=0, upper=1> use_softmax; // 0 - no softmax, 1 - use softmax on phi.
  int<lower=0, upper=1> use_normal; // 0 - Poisson, 1 - Normal
  real<lower=0> sigma2_prior_shape; // prior shape for sigma2
  real<lower=0> sigma2_prior_rate; // prior rate for sigma2
  real<upper=1> rho_value; // the fixed rho value. If < 0, then rho is estimated.
}
transformed data {
  int W_sparse[W_n, 2];   // adjacency pairs
  int<lower=0, upper=1> estimate_rho;
  vector[N] D_sparse;     // diagonal of D (number of neigbors for each site)
  vector[M] v_ones = rep_vector(1, M); // vector for computing row sums
  int<lower=0> y_obs_int[use_normal ? 0 : N_obs];
  
  { // generate sparse representation for W
  int counter;
  counter = 1;
  // loop over upper triangular part of W to identify neighbor pairs
    for (i in 1:(N - 1)) {
      for (j in (i + 1):N) {
        if (W[i, j] == 1) {
          W_sparse[counter, 1] = i;
          W_sparse[counter, 2] = j;
          counter = counter + 1;
        }
      }
    }
  }
  for (i in 1:N) D_sparse[i] = sum(W[i]); // Compute the sparse representation of D
  
  if(rho_value >= 0){
    estimate_rho = 0;
  }else{
    estimate_rho = 1;
  }
  
  // create integer valued outcomes for Poisson
  if(use_normal == 0){
    for(i in 1:N_obs){
	  y_obs_int[i] = 0;
	  while(y_obs_int[i] <= y_obs[i]){
	    y_obs_int[i] = y_obs_int[i] + 1;
	  }
	}
  }
}
parameters {
  real<lower=0> sigma2[use_normal ? 1 : 0]; // y variance of outcome
  real<lower=0> tau2[M]; // CAR variance parameter for each model
  real<lower=0, upper=1> rho_estimated[estimate_rho ? M : 0]; // spatial correlation for each model (set to size 0 if rho is fixed)
  matrix[N, M] phi; // CAR parameter: number of observations x number of models
}
transformed parameters {
  // variable declarations
  matrix[use_softmax ? N : 0, use_softmax ? M : 0] exp_phi;
  matrix[use_softmax ? N : 0, use_softmax ? M : 0] exp_phi_sum;
  matrix[N, M] u;
  vector[N] mu;
  vector[N_obs] observed_est;
  real log_detQ[M];
  matrix[N + 1, M] ldet_vec;
  real<lower = 0> sigma[use_normal ? 1 : 0];
  real rho[M];
  
  // variable calculations
  if(use_softmax == 1){
    exp_phi = exp(phi);
    for(m in 1:M){
      exp_phi_sum[1:N,m] = (exp_phi * v_ones);
    }
    u = exp_phi ./ exp_phi_sum;
  }else{
    u = 1.0/M + phi;
  }
  mu = (X .* u)*v_ones;
  observed_est = mu[ind_obs];
  
  // store the rho used
  if(estimate_rho == 0){
    for(m in 1:M){
	  rho[m] = rho_value;
    }
  }else{
    rho = rho_estimated;
  }
  
  // calculate the log determinants
  for (m in 1:M){
	ldet_vec[N + 1,m] = -N*log(tau2[m]);
	for (i in 1:N){
		ldet_vec[i,m] = log1p(rho[m]*lambda[i]);
	}
	log_detQ[m] = sum(ldet_vec[1:N+1,m]);
  }
  
  // transform the y Gaussian variance to standard deviation.
  if(use_normal == 1){
    sigma = sqrt(sigma2);
  }
}
model {
  if(use_normal == 1){
    sigma2 ~ gamma(sigma2_prior_shape, sigma2_prior_rate); // sigma2 prior
    y_obs ~ normal(observed_est, sigma); // normal likelihood
  }else{
    y_obs_int ~ poisson(observed_est); // Poisson likelihood
  }

  // CAR prior
  for(m in 1:M){
	phi[1:N, m] ~ sparse_car(tau2[m], rho[m], W_sparse, D_sparse, log_detQ[m], N, W_n);
  }
  // gamma prior on tau2 
  tau2 ~ gamma(1, 5);
}
generated quantities {
  if(use_normal == 1){
    vector[N] y_exp = mu;
	real y_pred[N] = normal_rng(mu, sigma);
	real log_likelihood = normal_lpdf(y_obs | observed_est, sigma);
  }else{
    vector[N] y_exp = exp(mu);
	int y_pred[N] = poisson_rng(y_exp);
	real log_likelihood = poisson_lpmf(y_obs_int | observed_est);
  }
}
