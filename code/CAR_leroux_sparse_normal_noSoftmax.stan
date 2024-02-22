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
}
transformed data {
  int W_sparse[W_n, 2];   // adjacency pairs
  vector[N] D_sparse;     // diagonal of D (number of neigbors for each site)
  vector[M] v_ones = rep_vector(1, M); // vector for computing row sums
  
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
}
parameters {
  real<lower=0> sigma2; // y variance of outcome
  real<lower=0> tau2[M]; // CAR variance parameter for each model
  real<lower=0, upper=1> rho[M]; // spatial correlation for each model
  matrix[N, M] phi; // CAR parameter: number of observations x number of models
}
transformed parameters {
  // variable declarations
  matrix[N, M] exp_phi;
  matrix[N, M] exp_phi_sum;
  matrix[N, M] u;
  vector[N] mu;
  vector[N_obs] observed_est;
  real log_detQ[M];
  matrix[N + 1, M] ldet_vec;
  real<lower = 0> sigma;
  
  // variable calculations
  u = 1.0/M + phi;
  mu = (X .* u)*v_ones;
  observed_est = mu[ind_obs];
  
  // calculate the log determinants
  for (m in 1:M){
	ldet_vec[N + 1,m] = -N*log(tau2[m]);
	for (i in 1:N){
		ldet_vec[i,m] = log1p(rho[m]*lambda[i]);
	}
	log_detQ[m] = sum(ldet_vec[1:N+1,m]);
  }
  
  // transform the y Gaussian variance to standard deviation.
  sigma = sqrt(sigma2);
}
model {
  // prior on sigma2
  sigma2 ~ gamma(1, 10);
  // likelihood
  // y_obs ~ poisson(observed_est);
  y_obs ~ normal(observed_est, sigma);
  // CAR prior
  for(m in 1:M){
	phi[1:N, m] ~ sparse_car(tau2[m], rho[m], W_sparse, D_sparse, log_detQ[m], N, W_n);
  }
  // gamma prior on tau2 
  tau2 ~ gamma(1, 5);
}
generated quantities {
  vector[N] y_exp = mu;
  int y_pred[N] = poisson_rng(mu);
  real log_likelihood = normal_lpdf(y_obs | observed_est, sigma);
}
