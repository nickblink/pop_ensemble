#!/usr/bin/env python
# coding: utf-8

# # 0_b High-level Wrapper Functions

# ## 0_b_1 Generate 1D Data
# * generate_data_1d

# In[62]:


OutcomeDistribution = Callable[[np.ndarray, np.ndarray], Any]

def generate_data_1d(N_train: int, 
                     N_base: int = 50, 
                     N_test: int = 1000,
                     y_dist: OutcomeDistribution = np.random.lognormal, 
                     seed: int = 0):
  """Generates 1D data for training and testingensemble model.
  
  Args:
    N_train: Number of extra training data points for ensemble model.
    N_base: Number of data points to train base model.
    N_test: Number of testing locations.
    y_dist: The data generating distribution of y that takes into two 
      vector-valued arguments (e.g., mean and variance).
    seed: Random seed
  
  Returns: 
    The (X, Y) data for base model, ensemble training, and ensemble testing. X is the analogue of location for an 
    air pollution application; Y is the analogue for concentration. 
  """
  np.random.seed(seed)
  tf.random.set_seed(seed)
  
  # Define data ranges
  train_range = (-3 * np.pi, 3 * np.pi)
  test_range = (-5 * np.pi, 5 * np.pi)

  # Build inputs X
  process_X_fn = lambda X: np.sort(X).astype(np.float32)

  X_base = process_X_fn(
      np.random.uniform(*train_range, size=(N_base, 1)))  # Data for training base model.
  X_train_extra = process_X_fn(
      np.random.uniform(*train_range, size=(N_train, 1)))  # Additional training data for ensemble.
  X_train = np.concatenate([X_base, X_train_extra])
  X_test = process_X_fn(
      np.random.uniform(*test_range, size=(N_test, 1)))  # Testing data covering whole feature space.

  # True functions for moments.
  f_mean = lambda x: np.sin(x)
  f_stddev = lambda x: tf.exp(np.cos(x) - 1.0).numpy()

  # Compute loc and scale as functions of input X
  loc_base = f_mean(X_base)
  loc_train = f_mean(X_train)
  loc_test = f_mean(X_test)

  std_base = f_stddev(X_base)
  std_train = f_stddev(X_train)
  std_test = f_stddev(X_test)

  # Sample outputs Y and convert to float32 for TF compatibility.
  Y_base = y_dist(loc_base, std_base).astype(np.float32)
  Y_train = y_dist(loc_train, std_train).astype(np.float32)
  Y_test = y_dist(loc_test, std_test).astype(np.float32)

  # Compute true mean and variance.
  N_test_sample = 1000
  Y_test_sample = np.stack(
      [y_dist(loc_test, std_test) for _ in range(N_test_sample)])
  mean_test = np.mean(Y_test_sample, axis=0)

  return X_base, X_train, X_test, Y_base, Y_train, Y_test, mean_test


# ## 0_b_2 Generate 2D Data
# * generate_data_2d
# * bird 
# * townsend 
# * rosenbrock
# * goldstein

# In[19]:


def generate_data_2d(N_train: int, 
                     N_base: int = 100, 
                     N_test: int = 5000,
                     y_dist: OutcomeDistribution = np.random.lognormal, 
                     seed: int = 0,
                     explicit_skewness: bool = False):
  """Generates 2D data for training ensemble model.
  
  Args:
    N_train: Number of extra training data points for ensemble model.
    N_base: Number of data points to train base model.
    N_test: Number of testing locations.
    y_dist: The data generating distribution of y that takes into two 
      vector-valued arguments (e.g., mean and variance).
    seed: Random seed.
    explicit_skewness: Whether to explicitly introduce skewness into the data-generating process.
  
  Returns: 
    The (X, Y) data for base model, ensemble training, and ensemble testing.
  """
  np.random.seed(seed)
  tf.random.set_seed(seed)

  # Data ranges.
  train_range = (-np.pi, np.pi)
  test_range = (-1.25 * np.pi, 1.25 * np.pi)

  # Build inputs X
  sample_X = lambda low, high, n: np.random.uniform(low, high, size=(n, 2)).astype(np.float32)

  X_base = sample_X(*train_range, N_base)  # Data for training base model.
  X_train_extra = sample_X(*train_range, N_train)  # Additional training data for ensemble.
  X_train = np.concatenate([X_base, X_train_extra])
  X_test = sample_X(*test_range, N_test)  # Testing data covering whole feature space.

  f_mean = lambda x, mplr=1.5: bird(x[:,0:1]*mplr, x[:,1:2]*mplr)/5.  # mplr 1.25 to 1.5
  # f_std = lambda x, mplr=1.0: townsend(x[:,0:1]*mplr, x[:,1:2]*mplr)/5.
  f_std = lambda x, mplr=0.5: rosenbrock(x[:,0:1]*mplr, x[:,1:2]*mplr)  # mplr 0.25 to 0.5
  f_skew = lambda x, mplr=0.5: np.clip(
      townsend(x[:,0:1]*mplr, x[:,1:2]*mplr), a_min=1e-5, a_max=np.inf)

  # Generate mean, std and skewness surfaces.
  loc_base, loc_train, loc_test = map(f_mean, [X_base, X_train, X_test])
  std_base, std_train, std_test = map(f_std, [X_base, X_train, X_test])
  skew_base, skew_train, skew_test = map(f_skew, [X_base, X_train, X_test])

  # DEBUG: make std very small.
  std_base = 1e-3
  std_train = 1e-3
  std_test = 1e-3

  dist_args_base = (loc_base, std_base)
  dist_args_train = (loc_train, std_train)
  dist_args_test = (loc_test, std_test)

  if explicit_skewness:
    dist_args_base += (skew_base,)
    dist_args_train += (skew_train,)
    dist_args_test += (skew_test,)

  # Sample outputs Y and convert to float32 for TF compatibility.
  Y_base = y_dist(*dist_args_base).astype(np.float32)
  Y_train = y_dist(*dist_args_train).astype(np.float32)
  Y_test = y_dist(*dist_args_test).astype(np.float32)

  # Compute true mean and variance.
  N_test_sample = 1000
  Y_test_sample = np.stack(
      [y_dist(*dist_args_test) for _ in range(N_test_sample)])
  mean_test = np.mean(Y_test_sample, axis=0)
  # std_test = np.std(Y_test_sample, axis=0)

  return X_base, X_train, X_test, Y_base, Y_train, Y_test, mean_test

# Helper functions.
bird = lambda x, y: (  # Mishra's Bird function.
    np.sin(y) * np.exp((1.-np.cos(x))**2) + 
    np.cos(x) * np.exp((1.-np.sin(y))**2) + (x-y)**2)/100. - 1.
townsend = lambda x, y: ((  # (Modified) Townsend function.
    -(np.cos(x - 1.) * y)**2 - x * np.sin(3 * x + y)) + 10.)
rosenbrock = lambda x, y: (  # (Modified) Rosenbrock function.
    100. * (y - x**2)**2 + (1 - x)**2)/25.
goldstein = lambda x, y: (
    1. + (x + y + 1.)**2 * (19. - 14.*x + 3.*x**2 - 14.*y + 6*x*y + 3. * y**2)
    ) * (30. + (2.*x - 3.*y)**2 * (18. - 32.*x + 12.*x**2 + 48.*y - 36.*x*y + 27.*y**2))


# ## 0_b_3 Run Posterior Interfence
# * run_posterior_inference

# In[21]:


def run_posterior_inference(model_dist: tfd.Distribution, 
                            Y: tf.Tensor, 
                            mcmc_config: Dict[str, Any], 
                            map_config: Optional[Dict[str, Any]] = None, 
                            model_config: Optional[Dict[str, Any]] = None,
                            initialize_from_map: bool = True):
  """Wrapper function for running MCMC with MAP initialization."""
  # Defines posterior log likelihood function, and also a 
  # randomly-sampled initial state from model prior.
  nchain = mcmc_config['nchain']
  init_state, target_log_prob_fn = prepare_mcmc(model_dist, Y, nchain=nchain)  
  
  if initialize_from_map:
    # Initializes at MAP, shape (num_chains, param_shape_0, param_shape_1).
    print('Running MAP:', end='\t')
    init_state = run_map(target_log_prob_fn=target_log_prob_fn, 
                         gp_config=model_config,
                         **map_config)

    init_state = tf.stack([init_state] * mcmc_nchain, axis=0)

  # Run MCMC, shape (param_shape_0, param_shape_1, num_chains).
  print('Running MCMC:', end='\t')
  gp_w_samples, _ = run_mcmc(init_state=init_state,
                             target_log_prob_fn=target_log_prob_fn,
                             **mcmc_config)  
  
  return gp_w_samples


# ## 0_b_4 Run Base Models
# * run_base_models

# In[22]:


def run_base_models(X_base, X_train, X_test, 
                    Y_base, Y_train, Y_test, kernels=None, 
                    num_train_steps=100, debug_mode=False):
  if kernels is None:
    kernels = [gpf.kernels.Matern52(lengthscales=0.5), 
               gpf.kernels.Polynomial(degree=3.), 
               gpf.kernels.ArcCosine(weight_variances=1., bias_variance=1.),
               gpf.kernels.Periodic(gpf.kernels.Exponential(lengthscales=2.))]

  n_models = len(kernels)
  kernel_names = [kernel.name for kernel in kernels]  

  models = [
    get_base_prediction(X_base, Y_base, X_train, kernel=k,
                        num_train_steps=num_train_steps, 
                        debug_mode=debug_mode) for k in kernels]
  base_preds_train = tf.stack([
    get_base_prediction(X_base, Y_base, X_train, model=m, debug_mode=debug_mode) for m in models], axis=-1)
  base_preds_test = tf.stack([
    get_base_prediction(X_base, Y_base, X_test, model=m, debug_mode=debug_mode) for m in models], axis=-1)
  
  return base_preds_train, base_preds_test, kernel_names


# ## 0_b_5 Run BMA Model
# * run_bma_model

# In[23]:


def run_bma_model(X_train, X_test, Y_train,
                  base_preds_train, base_preds_test,
                  gp_lengthscale=1.,
                  gp_l2_regularizer=0.1,
                  y_noise_std=0.1,
                  map_step_size=0.1,
                  map_num_steps=10_000,
                  mcmc_step_size=0.1,
                  mcmc_num_steps=10_000,
                  mcmc_nchain=10,
                  mcmc_burnin=2_500,
                  mcmc_initialize_from_map=False,
                  n_samples_eval=1000,
                  n_samples_train=100,
                  n_samples_test=200,
                  return_mcmc_examples=True,
                  seed=0, 
                  debug_mode=False):
  # Assemble model configs.
  model_config = DEFAULT_GP_CONFIG.copy()
  map_config = DEFAULT_MAP_CONFIG.copy()
  mcmc_config = DEFAULT_MCMC_CONFIG.copy()

  model_config['lengthscale']=gp_lengthscale  
  model_config['l2_regularizer']=gp_l2_regularizer
  model_config['y_noise_std']=y_noise_std

  map_config['learning_rate']=map_step_size
  map_config['num_steps']=map_num_steps

  mcmc_config['nchain']=mcmc_nchain
  mcmc_config['burnin']=mcmc_burnin
  mcmc_config['step_size']=mcmc_step_size
  mcmc_config['num_steps']=mcmc_num_steps
  
  # Model prior.
  bma_prior, gp_config = bma_dist(X_train, 
                                  base_preds_train, 
                                  **model_config)

  model_config.update(gp_config)  
  if debug_mode:
    print(model_config)

  # Run MCMC estimation.
  weight_samples = run_posterior_inference(
      model_dist=bma_prior, 
      model_config=model_config,
      Y=Y_train, 
      map_config=map_config,
      mcmc_config=mcmc_config,
      initialize_from_map=mcmc_initialize_from_map)  
  
  del bma_prior
  # Get posterior sample for all model parameters.
  bma_joint_samples = make_bma_samples(X_test, None, 
                                       base_preds_test, 
                                       bma_weight_samples=weight_samples[0],
                                       bma_model_config=model_config, 
                                       n_samples=n_samples_eval, 
                                       seed=seed,
                                       y_samples_only=False,
                                       debug_mode=debug_mode)

  if return_mcmc_examples:
    # MCMC data for BNE training, shapes (num_samples * num_data, ...)
    means_train_mcmc, X_train_mcmc, Y_train_mcmc = make_bma_samples(
        X_train, Y_train, 
        base_preds_train, 
        bma_weight_samples=weight_samples[0],
        bma_model_config=model_config,
        n_samples=n_samples_train,
        seed=seed, 
        prepare_mcmc_training=True,
        debug_mode=debug_mode)

    # MCMC data for BNE testing, shape (num_samples, num_data, num_output).
    means_test_mcmc = make_bma_samples(X_test, None, 
                                       base_preds_test, 
                                       bma_weight_samples=weight_samples[0],
                                       bma_model_config=model_config,
                                       n_samples=n_samples_test,
                                       seed=seed,
                                       debug_mode=debug_mode)
    
    return (bma_joint_samples, X_train_mcmc, Y_train_mcmc, 
            means_train_mcmc, means_test_mcmc)
    
  return bma_joint_samples


# ## 0_b_6 Run BNE Model
# * run_bne_model

# In[24]:


BNE_MOMENT_TYPES = ('none', 'mean', 'variance', 'skewness')
def run_bne_model(X_train: tf.Tensor,
                  Y_train: tf.Tensor,
                  X_test: tf.Tensor, 
                  base_model_samples_train: tf.Tensor,
                  base_model_samples_test: tf.Tensor,
                  moment_mode: str = 'none',
                  gp_lengthscale: float = 1., 
                  gp_l2_regularizer: float = 10., 
                  variance_prior_mean: float = 0., 
                  skewness_prior_mean: float = 0.,
                  map_step_size: float = 5e-3,
                  map_num_steps: float = 10_000,
                  mcmc_step_size: float = 1e-2,
                  mcmc_num_steps: float = 10_000,
                  mcmc_burnin: int = 2_500,
                  mcmc_nchain: int = 10,
                  mcmc_initialize_from_map: bool = False,
                  seed: int = 0,
                  debug_mode: bool = False):
  """Runs the full BNE model end-to-end.

  This function performs end-to-end estimation of BNE model. It supports four 
  modes as determined by `moment_mode`:

    * 'none': Estimate a constant-variance model y ~ N(m(x), scale).
    * 'mean': Estimate a residual process model y ~ N(m(x) + r(x), scale).
    * 'variance': Estimate a heterogeneous residual process model 
        y ~ N(m(x) + r(x), scale(x)).
    * 'skewness': Estimate a heterogeneous residual process model with 
        skewness y ~ EMG(m(x) + r(x), scale(x), rate(x)).

  Args:
    X_train: The training features (num_train, num_dim).
    Y_train: The training responses (num_train, 1).
    X_test: The testing features (num_test, num_dim).
    base_model_samples_train: The posterior samples m(x) from the base ensemble 
      model at train locations, shape (num_train, 1).
    base_model_samples_test: The posterior samples m(x) from the base ensemble 
      model at test locations, shape (num_test, 1).
    moment_mode: The highest moment for the BNE model to estimate. 
      Must be one of ('none', 'mean', 'variance', 'skewness').
    gp_lengthscale:  length scale for model Gaussian processes.
    gp_l2_regularizer: l2_regularizer for model Gaussian processes.
    bne_variance_prior_mean: Prior mean to variance parameter.  
    bne_skewness_prior_mean: Prior mean to skewness parameter.
    map_step_size: Step size for MAP optimization.
    map_num_steps: Number of train steps for MAP optimization. 
    mcmc_step_size: Step size for HMC sampling
    mcmc_num_steps: Number of train steps for HMC sampling.
    seed: The seed for generating posterior samples.

  Returns:
    A dictionary with model parameters and its MCMC samples. 
  """
  # Prepares model and estimation configs.
  gp_config = DEFAULT_GP_CONFIG.copy()
  model_config = DEFAULT_BNE_CONFIG.copy()
  map_config = DEFAULT_MAP_CONFIG.copy()
  mcmc_config = DEFAULT_MCMC_CONFIG.copy()

  gp_config['lengthscale'] = gp_lengthscale
  gp_config['l2_regularizer'] = gp_l2_regularizer

  model_config['estimate_mean'] = moment_mode in ('mean', 'variance', 'skewness')
  model_config['estimate_variance'] = moment_mode in ('variance', 'skewness')
  model_config['estimate_skewness'] = moment_mode in ('skewness',)
  model_config['variance_prior_mean'] = variance_prior_mean
  model_config['skewness_prior_mean'] = skewness_prior_mean
  model_config.update(gp_config)

  map_config['learning_rate'] = map_step_size
  map_config['num_steps']=map_num_steps

  mcmc_config['nchain'] = mcmc_nchain
  mcmc_config['burnin'] = mcmc_burnin
  mcmc_config['step_size'] = mcmc_step_size
  mcmc_config['num_steps'] = mcmc_num_steps

  # Constructs prior distribution.
  bne_prior, bne_gp_config = bne_model_dist(X_train,
                                            mean_preds=base_model_samples_train,
                                            **model_config)

  model_config.update(bne_gp_config)

  # Estimates GP weight posterior using MCMC.
  gp_weight_samples = run_posterior_inference(
      model_dist=bne_prior,
      model_config=bne_gp_config,
      Y=Y_train,
      map_config=map_config,
      mcmc_config=mcmc_config,
      initialize_from_map=mcmc_initialize_from_map)
  
  # Generates the posterior sample for all model parameters. 
  joint_samples = make_bne_samples(X_test,
                                   mean_preds=base_model_samples_test,
                                   bne_model_config=model_config,
                                   bne_weight_samples=gp_weight_samples[0],
                                   seed=seed,
                                   debug_mode=debug_mode)
  del bne_prior  
  return joint_samples

