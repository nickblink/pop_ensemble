functions {
  real partial_sparse_car_lpdf(array[] phi_slice, int start, int end, real tau2, real rho, 
                               int[,] W_sparse, vector D_sparse, real log_det_Q, int W_n) {
    row_vector[end - start + 1] phit_D;
    row_vector[end - start + 1] phit_W;
    
    phit_D = (phi_slice .* D_sparse[start:end])';
    phit_W = rep_row_vector(0, end - start + 1);
    for (i in 1:W_n) {
      if (W_sparse[i, 1] >= start && W_sparse[i, 1] <= end) {
        phit_W[W_sparse[i, 1] - start + 1] += phi_slice[W_sparse[i, 2] - start + 1];
      }
      if (W_sparse[i, 2] >= start && W_sparse[i, 2] <= end) {
        phit_W[W_sparse[i, 2] - start + 1] += phi_slice[W_sparse[i, 1] - start + 1];
      }
    }

    return 0.5 * (log_det_Q - (1/tau2) * (rho * (phit_D * phi_slice) - rho * (phit_W * phi_slice) + 
                                          (1 - rho) * dot_self(phi_slice)));
  }
  real partial_sparse_car_sum(array[] phi_slice, int start, int end, 
                                   real tau2, real rho, int[,] W_sparse, 
                                   vector D_sparse, real log_det_Q, int W_n) {
    return partial_sparse_car_lpdf(phi_slice | start, end, tau2, rho, W_sparse, 
                                   D_sparse, log_det_Q, W_n);
  }
}
data {
  int<lower=0> M; // number of models
  int<lower=0> N;  // number of observations
  int<lower=0> N_miss; // number of missing y points
  int<lower=0> N_obs; // number of observed y points
  int<lower=0, upper=N> ind_miss[N_miss]; // indices of missing y points
  int<lower=0, upper=N> ind_obs[N_obs]; // indices of observed y points
  matrix[N,M] X; // design matrix of ensemble models
  int<lower=0> y_obs[N_obs];  // output
  matrix<lower=0, upper = 1>[N, N] W; //adjacency matrix
  int W_n; // Number of adjacency pairs
  matrix[N,N] I; // Identity matrix
  vector[N] lambda; // the eigenvalues of the D - W - I matrix
  int<lower=0, upper=1> use_softmax; // 0 - no softmax, 1 - use softmax on phi.
  int<lower=0, upper=1> use_pivot; // 0 - no direct pivot, 1 - use pivot in last X value.
  real<lower=0> theta_prior_shape; // prior shape for theta
  real<lower=0> theta_prior_rate; // prior rate for theta
  real<lower=0> tau2_prior_shape; // prior shape for tau2
  real<lower=0> tau2_prior_rate; // prior rate for tau2
  real<upper=1> fixed_rho; // the fixed rho value. If < 0, then rho is estimated.
  real fixed_tau2; // the fixed tau2 value. If < 0, then tau2 is estimated
  real alpha_variance_prior; // the prior variance for alpha. If < 0, then alpha is not used.
  //int<lower=1> grainsize; //number of 
}
transformed data {
  int W_sparse[W_n, 2];   // adjacency pairs
  int<lower=0, upper=1> estimate_rho; // whether to estimate rho
  int<lower=0, upper=1> estimate_tau2; // whether to estimate tau2
  int<lower=0, upper=1> use_alpha; // whether to use alpha
  vector[N] D_sparse;     // diagonal of D (number of neigbors for each site)
  vector[M] v_ones = rep_vector(1, M); // vector for computing row sums
  int M_phi; // num models estimated for (can be different if using pivot).
  int grainsize = 10; 
  
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
  
  if(fixed_rho >= 0){
    estimate_rho = 0;
  }else{
    estimate_rho = 1;
  }
  
  if(fixed_tau2 >= 0){
    estimate_tau2 = 0;
  }else{
    estimate_tau2 = 1;
  }
  
  if(alpha_variance_prior > 0){
    use_alpha = 1;
  }else{
    use_alpha = 0;
  }
  
  if(use_pivot == 1){
    M_phi = M - 1;
  }else{
    M_phi = M;
  }
}

parameters {
  real<lower=0> theta; // y dispersion parameter.
  //real<lower=0> tau2[M]; // CAR variance parameter for each model.
  real<lower=0> tau2_estimated[estimate_tau2 ? M : 0]; // 
  real<lower=0, upper=1> rho_estimated[estimate_rho ? M : 0]; // spatial correlation for each model (set to size 0 if rho is fixed).
  real alpha[use_alpha ? M : 0]; // parameter for scaling models.
  matrix[N, M_phi] phi; // CAR parameter: number of observations x number of models.
}
transformed parameters {
  // variable declarations
  matrix[use_softmax ? N : 0, use_softmax ? M : 0] exp_phi;
  matrix[use_softmax ? N : 0, use_softmax ? M : 0] exp_phi_sum;
  matrix[N, M] u;
  vector[N] mu;
  vector[N] mu_bounded;
  vector[N_obs] observed_est;
  real log_detQ[M];
  matrix[N + 1, M] ldet_vec;
  real rho[M];
  real tau2[M];
  real phi_array_1[N];//
  real phi_array_2[N];//
  real phi_array_3[N];//
  
  // variable calculations
  if(use_softmax == 1){
    exp_phi = exp(phi);
    for(m in 1:M){
      exp_phi_sum[1:N,m] = (exp_phi * v_ones);
    }
    u = exp_phi ./ exp_phi_sum;
  }else{
    if(use_pivot == 1){
	  u[1:N,1:(M-1)] = 1.0/M + phi;
	  u[1:N,M] = 1.0 - u[1:N,1:(M-1)] * rep_vector(1, M-1);
	}else{
	  u = 1.0/M + phi;
	}
  }
  
  // scale the model weights.
  if(use_alpha == 1){
    for(m in 1:M){
	  u[1:N,m] = u[1:N,m]*alpha[m];
	}
  }
  
  // calculate mu.
  mu = (X .* u)*v_ones;
  
  // bound mu by 0.
  for (n in 1:N){
    if(mu[n] <= 0){
	  mu_bounded[n] = 0.01;
	}else{
	  mu_bounded[n] = mu[n];
	}
  }
  
  // get the observed predictions.
  observed_est = mu_bounded[ind_obs];
  
  // store the rho used
  if(estimate_rho == 0){
    for(m in 1:M){
	  rho[m] = fixed_rho;
    }
  }else{
    rho = rho_estimated;
  }
  
  // store the tau2 used
  if(estimate_tau2 == 0){
    for(m in 1:M){
	  tau2[m] = fixed_tau2;
    }
  }else{
    tau2 = tau2_estimated;
  }
  
  // calculate the log determinants
  for (m in 1:M){
	ldet_vec[N + 1,m] = -N*log(tau2[m]);
	for (i in 1:N){
		ldet_vec[i,m] = log1p(rho[m]*lambda[i]);
	}
	log_detQ[m] = sum(ldet_vec[1:N+1,m]);
  }
  
  for (i in 1:N){
	phi_array_1[i] = phi[i,1];
  }
  for (i in 1:N){
	phi_array_2[i] = phi[i,2];
  }
  for (i in 1:N){
	phi_array_3[i] = phi[i,3];
  }
}
model {
  theta ~ gamma(theta_prior_shape, theta_prior_rate); // prior on theta
  y_obs ~ neg_binomial_2(observed_est, theta); // likelihood

  // Parallelized CAR prior calculation
  target += reduce_sum(partial_sparse_car_sum, phi_array_1, grainsize, tau2[1], rho[1], 
  W_sparse, D_sparse, log_detQ[1], W_n);
  target += reduce_sum(partial_sparse_car_sum, phi_array_2, grainsize, tau2[2], rho[2], 
  W_sparse, D_sparse, log_detQ[2], W_n);
  target += reduce_sum(partial_sparse_car_sum, phi_array_3, grainsize, tau2[3], rho[3], 
  W_sparse, D_sparse, log_detQ[3], W_n);
  
  if (estimate_tau2 == 1) {
    tau2_estimated ~ gamma(tau2_prior_shape, tau2_prior_rate);
  }
  
  if (use_alpha == 1) {
    alpha ~ normal(1, alpha_variance_prior);
  }
}
generated quantities {
  real phi_dprior = 0;
  vector[N] y_exp = mu;
  real y_pred[N] = neg_binomial_2_rng(mu_bounded, theta);
}
