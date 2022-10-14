## libraries
from typing import Any, Callable, Dict, List, Optional, Union, Tuple

import os
import gc
import time
import pickle
import functools

import multiprocessing as mp

# from google.colab import files
# from google.colab import 

import numpy as np
import tensorflow as tf

import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns
import edward2 as ed
import tensorflow_probability as tfp

tfd = tfp.distributions
tfb = tfp.bijectors

dtype = tf.float32
import gpflow as gpf
import logging

logging.getLogger('tensorflow').setLevel(logging.ERROR)  # suppress pfor warnings
# Verify versions.
print(f'TensorFlow version: {tf.__version__}. Expected: 2.7.0')
print(f'TensorFlow Probability version: {tfp.__version__}. Expected: 0.15.0')
tf.test.gpu_device_name()

os.getcwd()

## Default Configs
# GP configs.
y_noise_std = 0.1  # @param
hidden_units = 128  # @param
lengthscale=1.  # @param
l2_regularizer=0.1  # @param

DEFAULT_GP_CONFIG = dict(lengthscale=lengthscale,
                         l2_regularizer=l2_regularizer, 
                         hidden_units=hidden_units, 
                         y_noise_std=y_noise_std)


# BNE model configs.
estimate_mean = "True" # @param ["True", "False"]
estimate_variance = "False" # @param ["True", "False"]
estimate_skewness = "False" # @param ["True", "False"]
variance_prior_mean=0. # @param
skewness_prior_mean=0. # @param

estimate_mean = eval(estimate_mean)
estimate_variance = eval(estimate_variance)
estimate_skewness = eval(estimate_skewness)

DEFAULT_BNE_CONFIG = dict(estimate_mean=estimate_mean,
                          estimate_variance=estimate_variance,
                          estimate_skewness=estimate_skewness,
                          variance_prior_mean=variance_prior_mean,
                          skewness_prior_mean=skewness_prior_mean)

# MAP configs.
map_step_size=0.1 # @param
map_num_steps=10_000 # @param

DEFAULT_MAP_CONFIG = dict(learning_rate=map_step_size,
                          num_steps=map_num_steps)

# MCMC configs.
mcmc_step_size=0.1 # @param
mcmc_sample_size=500 # @param
mcmc_num_steps=10_000 # @param
mcmc_burnin=2_500 # @param
mcmc_nchain=10 # @param
mcmc_seed=0 # @param

DEFAULT_MCMC_CONFIG = dict(step_size=mcmc_step_size, 
                           num_steps=mcmc_sample_size, 
                           burnin=mcmc_burnin, 
                           nchain=mcmc_nchain, 
                           seed=mcmc_seed)


## High-level Wrapper Functions
# @title Wrapper: generate_data_1d
OutcomeDistribution = Callable[[np.ndarray, np.ndarray], Any]

def generate_data_1d(N_train: int, 
                     N_base: int = 50, 
                     N_test: int = 1000,
                     y_dist: OutcomeDistribution = np.random.lognormal, 
                     seed: int = 0):
  """Generates 1D data for training ensemble model.
  
  Args:
    N_train: Number of extra training data points for ensemble model.
    N_base: Number of data points to train base model.
    N_test: Number of testing locations.
    y_dist: The data generating distribution of y that takes into two 
      vector-valued arguments (e.g., mean and variance).
    seed: Random seed
  
  Returns: 
    The (X, Y) data for base model, ensemble training, and ensemble testing.
  """
  np.random.seed(seed)
  tf.random.set_seed(seed)
  
  # Data ranges.
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
  X_test = np.sort(X_test, axis=0)  # Sort testing data for visualization.

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

  # @title Wrapper: generate_data_2d
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
    explicit_skewness: Whether to explicitly introduce skewness.
  
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





# @title Wrapper: run_posterior_inference
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

# @title Wrapper: run_base_models
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

# @title Wrapper: run_bma_model
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

# @title Wrapper: run_bne_model
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
  
  # Generates the e for all model parameters. 
  joint_samples = make_bne_samples(X_test,
                                   mean_preds=base_model_samples_test,
                                   bne_model_config=model_config,
                                   bne_weight_samples=gp_weight_samples[0],
                                   seed=seed,
                                   debug_mode=debug_mode)
  del bne_prior  
  return joint_samples

## Core Model Functions

# @title Distribution: ExponentiallyModifiedGaussian
from tensorflow_probability.python.bijectors import identity as identity_bijector
from tensorflow_probability.python.bijectors import softplus as softplus_bijector
from tensorflow_probability.python.distributions import distribution
from tensorflow_probability.python.distributions import exponential as exponential_lib
from tensorflow_probability.python.distributions import normal as normal_lib
from tensorflow_probability.python.internal import assert_util
from tensorflow_probability.python.internal import dtype_util
from tensorflow_probability.python.internal import parameter_properties
from tensorflow_probability.python.internal import reparameterization
from tensorflow_probability.python.internal import samplers
from tensorflow_probability.python.internal import special_math
from tensorflow_probability.python.internal import tensor_util
from tensorflow_probability.python.math import generic as tfp_math


class ExponentiallyModifiedGaussian(
    distribution.AutoCompositeTensorDistribution):
  """Exponentially modified Gaussian distribution.
  #### Mathematical details
  The exponentially modified Gaussian distribution is the sum of a normal
  distribution and an exponential distribution.
  ```none
  X ~ Normal(loc, scale)
  Y ~ Exponential(rate)
  Z = X + Y
  ```
  is equivalent to
  ```none
  Z ~ ExponentiallyModifiedGaussian(loc, scale, rate)
  ```
  #### Examples
  ```python
  tfd = tfp.distributions
  # Define a single scalar ExponentiallyModifiedGaussian distribution
  dist = tfd.ExponentiallyModifiedGaussian(loc=0., scale=1., rate=3.)
  # Evaluate the pdf at 1, returing a scalar.
  dist.prob(1.)
  ```
  """

  def __init__(self,
               loc,
               scale,
               rate,
               validate_args=False,
               allow_nan_stats=True,
               name='ExponentiallyModifiedGaussian'):
    """Construct an exponentially-modified Gaussian distribution.
    The Gaussian distribution has mean `loc` and stddev `scale`,
    and Exponential distribution has rate parameter `rate`.
    The parameters `loc`, `scale`, and `rate` must be shaped in a way that
    supports broadcasting (e.g. `loc + scale + rate` is a valid operation).
    Args:
      loc: Floating-point `Tensor`; the means of the distribution(s).
      scale: Floating-point `Tensor`; the stddevs of the distribution(s). Must
        contain only positive values.
      rate: Floating-point `Tensor`; the rate parameter for the exponential
        distribution.
      validate_args: Python `bool`, default `False`. When `True` distribution
        parameters are checked for validity despite possibly degrading runtime
        performance. When `False` invalid inputs may silently render incorrect
        outputs.
      allow_nan_stats: Python `bool`, default `True`. When `True`, statistics
        (e.g., mean, mode, variance) use the value "`NaN`" to indicate the
        result is undefined. When `False`, an exception is raised if one or more
        of the statistic's batch members are undefined.
      name: Python `str` name prefixed to Ops created by this class.
    Raises:
      TypeError: if `loc`, `scale`, and `rate` are not all the same `dtype`.
    """
    parameters = dict(locals())
    with tf.name_scope(name) as name:
      dtype = dtype_util.common_dtype([loc, scale, rate], dtype_hint=tf.float32)
      self._loc = tensor_util.convert_nonref_to_tensor(
          loc, dtype=dtype, name='loc')
      self._scale = tensor_util.convert_nonref_to_tensor(
          scale, dtype=dtype, name='scale')
      self._rate = tensor_util.convert_nonref_to_tensor(
          rate, dtype=dtype, name='rate')
      super(ExponentiallyModifiedGaussian, self).__init__(
          dtype=dtype,
          reparameterization_type=reparameterization.FULLY_REPARAMETERIZED,
          validate_args=validate_args,
          allow_nan_stats=allow_nan_stats,
          parameters=parameters,
          name=name)

  @staticmethod
  def _param_shapes(sample_shape):
    return dict(
        zip(('loc', 'scale', 'rate'),
            ([tf.convert_to_tensor(sample_shape, dtype=tf.int32)] * 3)))

  @classmethod
  def _parameter_properties(cls, dtype, num_classes=None):
    return dict(
        loc=parameter_properties.ParameterProperties(),
        scale=parameter_properties.ParameterProperties(
            default_constraining_bijector_fn=(
                lambda: softplus_bijector.Softplus(low=dtype_util.eps(dtype)))),
        rate=parameter_properties.ParameterProperties(
            default_constraining_bijector_fn=(
                lambda: softplus_bijector.Softplus(low=dtype_util.eps(dtype)))))

  @property
  def loc(self):
    """Distribution parameter for the mean of the normal distribution."""
    return self._loc

  @property
  def scale(self):
    """Distribution parameter for standard deviation of the normal distribution."""
    return self._scale

  @property
  def rate(self):
    """Distribution parameter for rate parameter of exponential distribution."""
    return self._rate

  def _event_shape_tensor(self):
    return tf.constant([], dtype=tf.int32)

  def _event_shape(self):
    return tf.TensorShape([])

  def _sample_n(self, n, seed=None):
    normal_seed, exp_seed = samplers.split_seed(seed, salt='emg_sample')
    # need to make sure component distributions are broadcast appropriately
    # for correct generation of samples
    loc = tf.convert_to_tensor(self.loc)
    rate = tf.convert_to_tensor(self.rate)
    scale = tf.convert_to_tensor(self.scale)
    batch_shape = self._batch_shape_tensor(loc=loc, scale=scale, rate=rate)
    loc_broadcast = tf.broadcast_to(loc, batch_shape)
    rate_broadcast = tf.broadcast_to(rate, batch_shape)
    normal_dist = normal_lib.Normal(loc=loc_broadcast, scale=scale)
    exp_dist = exponential_lib.Exponential(rate_broadcast)
    x = normal_dist.sample(n, normal_seed)
    y = exp_dist.sample(n, exp_seed)
    return x + y

  def _log_prob(self, x):
    loc = tf.convert_to_tensor(self.loc)
    rate = tf.convert_to_tensor(self.rate)
    scale = tf.convert_to_tensor(self.scale)
    two = dtype_util.as_numpy_dtype(x.dtype)(2.)
    z = (x - loc) / scale
    w = rate * scale
    return (tf.math.log(rate) + w / two * (w - 2 * z) +
            special_math.log_ndtr(z - w))

  def _log_cdf(self, x):
    rate = tf.convert_to_tensor(self.rate)
    scale = tf.convert_to_tensor(self.scale)
    x_centralized = x - self.loc
    u = rate * x_centralized
    v = rate * scale
    vsquared = tf.square(v)
    return tfp_math.log_sub_exp(
        special_math.log_ndtr(x_centralized / scale),
        -u + vsquared / 2. + special_math.log_ndtr((u - vsquared) / v))

  def _mean(self):
    return tf.broadcast_to(
        self.loc + 1 / self.rate, self._batch_shape_tensor())

  def _variance(self):
    return tf.broadcast_to(
        tf.square(self.scale) + 1 / tf.square(self.rate),
        self._batch_shape_tensor())

  def _parameter_control_dependencies(self, is_init):
    assertions = []

    if is_init:
      try:
        self._batch_shape()
      except ValueError:
        raise ValueError(
            'Arguments `loc`, `scale`, and `rate` must have compatible shapes; '
            'loc.shape={}, scale.shape={}, rate.shape={}.'.format(
                self.loc.shape, self.scale.shape, self.rate.shape))
      # We don't bother checking the shapes in the dynamic case because
      # all member functions access both arguments anyway.

    if not self.validate_args:
      assert not assertions  # Should never happen.
      return []

    if is_init != tensor_util.is_ref(self.scale):
      assertions.append(assert_util.assert_positive(
          self.scale, message='Argument `scale` must be positive.'))

    if is_init != tensor_util.is_ref(self.rate):
      assertions.append(assert_util.assert_positive(
          self.rate, message='Argument `rate` must be positive.'))

    return assertions

  def _default_event_space_bijector(self):
    return identity_bijector.Identity(validate_args=self.validate_args)

# @title Model: Random Feature.

def make_random_feature(X: tf.Tensor, 
                        lengthscale: float = 1., 
                        hidden_units: int = 512, 
                        seed: int = 0, 
                        return_Wb: bool = False, 
                        **unused_kwargs) -> Union[tf.Tensor, Tuple[tf.Tensor]]:
  """Makes random feature for scalable RBF kernel approximation.

  Given data point x with shape (P, 1), the random feature is computed as 
  
              Phi(x) = sqrt(2/M) * cos(Wx + b), where

  W and b are untrainable random weights (shape [M, P]) and biases 
  (shape [M, 1]) initialized from Gaussian and uniform distribution, 
  respectively. Here M is the number of hidden units.
  
  Args:
    X: The input data of shape (batch_size, feature_dim).
    lengthscale: The length scale of RBF kernel.
    hidden_units: The number of hidden units used to approximate the exact 
      Gaussian process posterior. Higher number leads to more exact 
      approximation, although it doesn't necessarily leads to better prediction
      performance.
    seed: The random seed used to generate the random features.
    return_Wb: Whether to return the weights and bias (W, b) of the random 
      feature mapping cos().
    unused_kwargs: Unused keyword arguments for signiture consistency with other
      prior functions.

  Returns:
    Random feature Phi(X) with shape (batch_size, hidden_units). Or 
    the random weights and biases if return_Wb=True.
  """
  del unused_kwargs

  # Configure shapes.
  feature_dim = tf.shape(X)[1]
  W_shape = (feature_dim, hidden_units)
  b_shape = (hidden_units,)

  # Apply lengthscale.
  X = tf.convert_to_tensor(X, dtype=dtype)
  lengthscale = tf.convert_to_tensor(lengthscale, dtype=dtype)
  
  X = X / lengthscale

  # Sample random features.
  tf.random.set_seed(seed)
  W_dist = ed.initializers.OrthogonalRandomFeatures(stddev=1.)
  b_dist = tf.initializers.RandomUniform(0, 2*np.pi)
  W = W_dist(W_shape).numpy()
  b = b_dist(b_shape).numpy()

  multiplier = tf.sqrt(2./hidden_units)
  random_feature = multiplier * tf.math.cos(tf.matmul(X, W) + b)
  
  if return_Wb:
    return W, b
  return random_feature

# @title Model: Gaussian Process.
def rfgp_dist(inputs: tf.Tensor,
              units: int = 1,
              hidden_units: int = 128,               
              lengthscale: float = 1., 
              l2_regularizer: float = 1., 
              y_noise_std: float = 0.1,
              seed: int = 0,
              posterior_mode: bool = False,
              posterior_sample: Optional[tf.Tensor] = None,
              return_joint_dist: bool = True,
              verbose: bool = False
              ) -> Tuple[Union[tfd.Distribution, List[tfd.Distribution]], 
                         Dict[str, Any]]:
  """Specifies a random-feature-based Gaussian process (RFGP) distribution.
  
  Given random feature mapping Phi(x), a Gaussian process prior 
  f ~ GaussianProcess can be expressed as:

      f(x) ~ matmul(Phi(x), W);   W i.i.d. ~ Normal(0, s)

  where W is a i.i.d. Normal prior with shape (num_hidden_units, units), and 
  its standard deviation s serves as a regularizer to control the smoothness of
  f(x) (similar to the l2 penalty in Ridge regression).

  This function also allow user to specify a Gaussian process regression model
  (by setting y_noise_std>0), that is, it will specify a prior model like below:

      y | f ~ N(f(x), y_noise_std);   f ~ GaussianProcess (as specified above)

  Args:
    inputs: Input data of size (batch_size, feature_dim).
    units: Number of outputs from the Gaussian process.
    hidden_units: The number of hidden units used to approximate the exact 
      Gaussian process posterior. Higher number leads to more exact 
      approximation, although it doesn't necessarily leads to better prediction
      performance.
    lengthscale: The length scale of RBF kernel.
    l2_regularizer: The L2 regularizer to Gaussian process prior, which is 
      effectively the inverse standard deviation of the prior distribution 
      of beta.
    y_noise_std: The standard deviation of the output y given f. If -1 then 
      y | f is a deterministic distribution (i.e., y = f).
    seed: The random seed used to generate the random features.
    posterior_mode: Whether to specify the posterior distribution of the 
      Gaussian process instead of the prior distribution. If `True`, then 
      the distribution of `beta` is specified as a sampler from the posterior 
      sample.
    posterior_sample: The posterior sample of `beta` with dimension 
      (num_mcmc_sample, hidden_units, units) or 
      (hidden_units, units, num_mcmc_sample). Must be provided if `posterior_mode=True`.
    return_joint_dist: If True, return a tfd.Distribution object which 
      represents the joint distribution of the model random variables [W, y], 
      which can be used for MCMC sampling. If False, return a list of random
      variables which can be used as a building block for more complex models,
      and also add the model random features Phi(X) to the model_config.

  Returns:
    joint_dist: The random variables in the Gaussian process prior [W, y].
    model_config: A dictionary of keyword arguments to make_random_feature.
  """
  feature_dim = tf.shape(inputs)[1]
  # Latent random features.
  random_features = make_random_feature(inputs, 
                                        lengthscale=lengthscale, 
                                        hidden_units=hidden_units, 
                                        seed=seed)
    
  # Parameter Distributions.
  W_dist = tfd.Normal(loc=tf.zeros(shape=(hidden_units, units)), 
                      scale=1./l2_regularizer)
  if posterior_mode:
    # If in posterior mode, defines W_dist as the Empirical distribution 
    # constructed from the posterior samples.
    if (posterior_sample.shape[1] == hidden_units and 
        posterior_sample.shape[2] == units):
      # Convert the posterior_sample shape from  
      # (num_sample, hidden_units, output_units) to
      # (hidden_units, output_units, num_sample).
      posterior_sample = tf.transpose(posterior_sample, [1, 2, 0])

    sample_hidden_dim = posterior_sample.shape[0]
    sample_output_dim = posterior_sample.shape[1]
    num_sample = posterior_sample.shape[2]

    if (sample_hidden_dim != hidden_units or 
        sample_output_dim != units or 
        len(posterior_sample.shape) != 3):
      raise ValueError(
          "Expect posterior_sample shape to be "
          f"[num_sample, hidden_units({hidden_units}),units({units})] or "
          f"[hidden_units({hidden_units}),units({units}), num_sample]. Got "
          f"{posterior_sample.shape}.")

    if verbose:
      print(f"Constructing posterior from {num_sample} samples.")
    W_dist = tfd.Empirical(posterior_sample, event_ndims=0)
  
  # Outcome Distributions.
  if y_noise_std > 0.:
    y_dist = lambda gp_weights: tfd.Normal(
        loc=tf.linalg.matmul(random_features, gp_weights), scale=y_noise_std)
  else:
    y_dist = lambda gp_weights: tfd.Deterministic(
        loc=tf.linalg.matmul(random_features, gp_weights))
        
  joint_dist = dict(gp_weights=W_dist, y=y_dist)

  if return_joint_dist:
    joint_dist = tfd.JointDistributionNamedAutoBatched(joint_dist)

  # Model config.
  model_config = dict(units=units,
                      hidden_units=hidden_units,                             
                      lengthscale=lengthscale, 
                      l2_regularizer=l2_regularizer, 
                      y_noise_std=y_noise_std,
                      seed=seed)
  
  if not return_joint_dist:
    model_config['random_features'] = random_features

  return joint_dist, model_config

# @title Model: Bayesian Model Averaging.
def bma_dist(inputs: tf.Tensor, 
             base_model_preds: tf.Tensor, 
             y_noise_std: float = 0.1,
             posterior_mode: bool = False,
             posterior_sample: Optional[tf.Tensor] = None,
             sample_intermediate_variables: bool = False,
             return_joint_dist: bool = True,
             debug_mode: bool = False,
             **gp_kwargs):
  """Specifies an adaptive Bayesian model averaging (BMA) model.
  
  This function specifies an adaptive ensemble model with model weights 
  parameterized by Gaussian processes:

      y | f ~ N(m(x), y_noise_std)  where  m = sum_k p_k(x) * m_k(x),       
      
  where m_k(x) are base model predictions, and p_k's are adaptive ensemble 
  weights specified as:

      [p_1(x), ..., p_K(x)] = softmax([f_1(x), ..., f_K(x)])
      [f_1, ..., f_K]  ~ i.i.d.  GaussianProcess

  Args:
    inputs: Input data of size (batch_size, feature_dim).
    base_model_preds: Base model predictions of shape (batch_size, num_base_models).
    y_noise_std: The standard deviation of the output distribution y | f.
    posterior_mode: Whether to specify the posterior distribution of the 
      Gaussian process instead of the prior distribution. If `True`, then 
      the distribution of `beta` is specified as a sampler from the posterior 
      sample.
    posterior_sample: The posterior sample of `beta` with dimension 
      (hidden_units, units, num_mcmc_sample). Must be provided if 
      `posterior_mode=True`.
    sample_intermediate_variables: Whether to generate posterior samples of
      intermediate variables for model understanding when `posterior_mode=True`.
    return_joint_dist: If True, return a tfd.Distribution object which 
      represents the joint distribution of the model random variables [W, y], 
      which can be used for MCMC sampling. If False, return a list of random
      variables which can be used as a building block for more complex models,
      and also add the model random features Phi(X) to the model_config.

  Returns:
    joint_dist: The random variables in the Gaussian process prior [W, y].
    model_config: A dictionary of keyword arguments to make_random_feature.
  """
  num_model_data, num_model = base_model_preds.shape
  num_input_data, num_input_features = inputs.shape

  if num_model_data != num_input_data:
    raise ValueError(
        f'Number of data points in input data ({num_input_data}) and in model '
        f'predictions ({num_model_data}) should be equal.')

  # Specifies Gaussian process priors for ensemble weights.
  gp_kwargs.pop('y_noise_std', None)
  gp_kwargs.pop('units', None)
  gp_dists, gp_config = rfgp_dist(
      inputs, 
      units=num_model, 
      posterior_mode=posterior_mode,
      posterior_sample=posterior_sample,
      return_joint_dist=False,
      verbose=debug_mode,
      **gp_kwargs)

  gp_features = gp_config.pop('random_features')
  gp_weight_dist = gp_dists["gp_weights"]  

  # Specifies Bayesian model averaging model.
  joint_dist = dict()

  joint_dist["gp_weights"] = gp_weight_dist
  if posterior_mode and sample_intermediate_variables:
    # Generates posterior of all intermediate variables.
    joint_dist["gps"] = lambda gp_weights: tfd.Deterministic(
        loc=build_ensemble_weight_logits(gp_weights, gp_features))
    joint_dist["ensemble_weights"] = lambda gps: tfd.Deterministic(
        loc=build_ensemble_weights(gps))
    joint_dist["y"] = lambda ensemble_weights: tfd.Normal(
        loc=ensemble_prediction(ensemble_weights, base_model_preds), 
        scale=y_noise_std)
  else:
    # Use collapsed joint distribution for easy MCMC sampling.
    joint_dist["y"] = lambda gp_weights: tfd.Normal(
        loc=ensemble_mean(gp_weights, base_model_preds, gp_features), 
        scale=y_noise_std)
  
  if return_joint_dist:
    joint_dist = tfd.JointDistributionNamedAutoBatched(joint_dist) 
  
  # Store model config.
  model_config = gp_config
  model_config['y_noise_std']=y_noise_std

  return joint_dist, model_config
  

# Utility functions.
def build_ensemble_weight_logits(gp_weights, gp_features):
  """Builds Gaussian process prediction from random-feature weights."""
  return tf.linalg.matmul(gp_features, gp_weights) 

def build_ensemble_weights(logits):
  """Builds ensemble weights from Gaussian process prediction."""
  return tf.keras.activations.softmax(logits)

def ensemble_prediction(weights, base_model_preds):
  """Builds final ensemble prediction from ensemble weights and base models."""
  return tf.reduce_sum(base_model_preds * weights, axis=-1, keepdims=True)

def ensemble_mean(gp_weights, base_model_preds, gp_features):
  """Computes final ensemble prediction directly from random-feature weights"""
  logits = build_ensemble_weight_logits(gp_weights, gp_features)
  ensemble_weights = build_ensemble_weights(logits)
  return ensemble_prediction(ensemble_weights, base_model_preds)

# @title Model: Bayesian Nonparametric Ensemble
def bne_model_dist(inputs: tf.Tensor, 
                   mean_preds: Union[float, tf.Tensor],
                   estimate_mean=True,
                   estimate_variance=True,
                   estimate_skewness=True,
                   variance_prior_mean=0.,
                   skewness_prior_mean=0.,
                   posterior_mode: bool = False,
                   posterior_sample: Optional[tf.Tensor] = None,
                   return_joint_dist: bool = True,
                   sample_intermediate_variables: bool = False,
                   debug_mode: bool = False,
                   **gp_kwargs):
  """BNE model with spatially adaptive variances and skewness.

  Given the original prediction m(x) from a Bayesian model, the Bayesian 
  nonparametric ensemble (BNE) improves the original model's uncertainty 
  quality by deploying a semi-parametric distribution estimator that flexibly
  captures the higher moments of the data distribution (i.e., the aleatoric 
  uncertainty) while also quantifying model's epistemic uncertainty (due to the 
  lack of data).

  Specifically, we assume:

  y ~ ExponentiallyModifiedGaussian(m(x) + loc(x), scale(x), rate(x)) 

  where loc, scale, rate are (transformed) Gaussian processes that adaptively 
  estimates higher moments of the data distribution:

    loc         ~ GaussianProcess(0, k_loc)
    log(scale)  ~ GaussianProcess(0, k_scale)
    log(rate)   ~ GaussianProcess(0, k_rate)   
    
  
  Args:
    inputs: Input data of size (batch_size, feature_dim).
    mean_preds: Original mean predictions. Either a scalar or an array of shape 
      (batch_size, 1).
    y_noise_std: The standard deviation of the output distribution y | f.
    estimate_mean: Whether to add residual process to model.
    estimate_variance: Whether to estimate spatio-temporally adaptive variance 
      of the outcome distribution using Gaussian process.
    estimate_skewness: Whether to estimate spatio-temporally adaptive skewness 
      of the outcome distribution using Gaussian process.
    skewness_prior_mean: The prior mean for the Gaussian process logits for the 
      skewness parameter. If user which the model to not estimate the skewness,
      can set this to a very large negative value. e.g., -25.0.
    posterior_mode: Whether to specify the posterior distribution of the 
      Gaussian process instead of the prior distribution. If `True`, then 
      the distribution of `beta` is specified as a sampler from the posterior 
      sample.
    posterior_sample: The posterior sample of `beta` with dimension 
      (hidden_units, units, num_mcmc_sample). Must be provided if 
      `posterior_mode=True`.
    return_joint_dist: If True, return a tfd.Distribution object which 
      represents the joint distribution of the model random variables [W, y], 
      which can be used for MCMC sampling. If False, return a list of random
      variables which can be used as a building block for more complex models,
      and also add the model random features Phi(X) to the model_config.
    sample_intermediate_variables: Whether to generate posterior samples of
      intermediate variables for model understanding when `posterior_mode=True`.
    **gp_kwargs: Optional keyword arguments to the rfgp_dist() which generates 
      the Gaussian processes.

  Returns:
    joint_dist: The random variables in the Gaussian process prior [W, y].
    model_config: A dictionary of keyword arguments to make_random_feature.  
  """
  num_input_data, num_input_features = inputs.shape

  # Specify GP priors.
  gp_kwargs.pop('y_noise_std', None)
  gp_kwargs.pop('units', None)

  bne_units = 2 + int(estimate_skewness)
  bne_gp_dists, bne_gp_config = rfgp_dist(
      inputs, 
      # Perform inference for mean, optionlly variance and skewness.      
      units=bne_units,  
      y_noise_std=-1.,  
      posterior_mode=posterior_mode,
      posterior_sample=posterior_sample,
      return_joint_dist=False,
      verbose=debug_mode,
      **gp_kwargs)

  gp_weight_dist = bne_gp_dists["gp_weights"]
  gp_features = bne_gp_config.pop('random_features')

  if debug_mode:
    print(f'estimate_mean: {estimate_mean}')
    print(f'estimate_variance: {estimate_variance}')
    print(f'estimate_skewness: {estimate_skewness}')

  # Specify mean distribution.
  if estimate_mean:
    def mean_from_gps(gps):
      return gps[:, 0, tf.newaxis]
  else:
    def mean_from_gps(gps):
      return tf.zeros_like(gps)[:, 0, tf.newaxis]

  def mean_dist_from_gps(gps):
    return tfd.Deterministic(loc=mean_from_gps(gps))

  # Specify variance distribution.  
  if estimate_variance:
    # Model variance as a log Gaussian process.
    def scale_from_gps(gps, gp_weights):
      return tf.exp(gps[:, 1, tf.newaxis] + variance_prior_mean)
  else:
    # Model variance as a scalar (by taking an element from GP weights).
    def scale_from_gps(gps, gp_weights):
      return tf.exp(gp_weights[0, 1] + variance_prior_mean)
      
  def scale_dist_from_gps(gps, gp_weights):
    return tfd.Deterministic(loc=scale_from_gps(gps, gp_weights))

  # Specify skewness distribution.
  if estimate_skewness:
    # Estimate mean, variance and skewness using EMG distribution.
    def rate_from_gps(gps):
      return tf.exp(gps[:, 2, tf.newaxis] + skewness_prior_mean)

    def y_dist_from_moments(mean, scale, rate):
      return ExponentiallyModifiedGaussian(
          loc=mean_preds + mean, scale=scale, rate=rate)

    def y_dist_from_weights(gp_weights):
      gps = tf.matmul(gp_features, gp_weights)
      mean = mean_from_gps(gps)
      scale = scale_from_gps(gps, gp_weights)
      rate = rate_from_gps(gps)
      return ExponentiallyModifiedGaussian(
          loc=mean_preds + mean, scale=scale, rate=rate)
  else:
    # Estimate mean and variance using Gaussian distribution.
    def y_dist_from_moments(mean, scale):
      return tfd.Normal(loc=mean_preds + mean, scale=scale)

    def y_dist_from_weights(gp_weights):
      gps = tf.matmul(gp_features, gp_weights)
      mean = mean_from_gps(gps)
      scale = scale_from_gps(gps, gp_weights)
      return tfd.Normal(loc=mean_preds + mean, scale=scale)
  
  # Specify full joint distributions.
  joint_dist = dict()
  if posterior_mode and sample_intermediate_variables:
    joint_dist['gp_weights'] = gp_weight_dist
    joint_dist['gps'] = lambda gp_weights: tfd.Deterministic(
        loc=tf.matmul(gp_features, gp_weights))
    joint_dist['mean'] = mean_dist_from_gps
    joint_dist['scale'] = scale_dist_from_gps
    if estimate_skewness:
      joint_dist['rate'] = lambda gps: tfd.Deterministic(
          loc=rate_from_gps(gps))
      joint_dist['skewness'] = lambda rate: tfd.Exponential(rate=rate)
    joint_dist['y'] = y_dist_from_moments
  else: 
    joint_dist['gp_weights'] = gp_weight_dist
    joint_dist['y'] = y_dist_from_weights

  if return_joint_dist:
    joint_dist = tfd.JointDistributionNamedAutoBatched(joint_dist)
  
  return joint_dist, bne_gp_config

# @title MCMC: Low-level sampler.

@tf.function(autograph=False)
def run_chain(init_state: List[tf.Tensor], 
              step_size: float, 
              target_log_prob_fn: Callable[..., tf.Tensor], 
              num_steps: int = 500, 
              burnin: int = 100, 
              seed: int = 0,
              kernel_type: str = "hmc",
              step_adaptor_type: str = "simple"
              ) -> Union[List[tf.Tensor], Tuple[tf.Tensor]]:
  """Low-level function that runs MCMC sampling for a given model posterior.
  
  Args:
    init_state: The initial state for the MCMC sampler.
    step_size: The step size of a Hamiltonian Monte Carlo step.
    target_log_prob_fn: The log likelihood function for model posterior.
    num_steps: The number of total MCMC samples to return.
    burnin: The length of the burn-in period for MCMC warmup.
    seed: The random seed for MCMC sampling.
    kernel_type: Type of MCMC kernel to use, either ('hmc', 'nuts').
    step_adaptor_type: Type of MCMC kernel to use, one of 
      ('simple', 'dual_averaging').

  Returns:
    chain_state: Posterior sample from all MCMC chains.
    sampler stat: Sampling statistics, currently (step_size, acceptance ratio).
  """
  if kernel_type not in ('hmc', 'nuts'):
    raise ValueError(
        f"kernel_type {kernel_type} must be one of ('hmc', 'nuts').")

  if step_adaptor_type not in ('simple', 'dual_averaging'):
    raise ValueError(
        f"step_adaptor_type {step_adaptor_type} must be one of "
        "('simple', 'dual_averaging').")

  def trace_fn(_, pkr): 
    if kernel_type is 'hmc':
      step_size = pkr.inner_results.accepted_results.step_size
    else:
      step_size = pkr.inner_results.step_size

    return (step_size, pkr.inner_results.log_accept_ratio)

  if kernel_type is 'hmc':
    kernel = tfp.mcmc.HamiltonianMonteCarlo(
        target_log_prob_fn=target_log_prob_fn,
        num_leapfrog_steps=5,
        step_size=step_size)
    step_adaptation_kwargs = dict()
  else:
    kernel = tfp.mcmc.NoUTurnSampler(
        target_log_prob_fn=target_log_prob_fn,
        step_size=step_size)
    step_adaptation_kwargs = dict(
        step_size_setter_fn=lambda pkr, new_step_size: pkr._replace(
            step_size=new_step_size),
        step_size_getter_fn=lambda pkr: pkr.step_size,
        log_accept_prob_getter_fn=lambda pkr: pkr.log_accept_ratio,)

  if step_adaptor_type is 'simple':
    kernel = tfp.mcmc.SimpleStepSizeAdaptation(
      inner_kernel=kernel, 
      num_adaptation_steps=burnin)
  else:
    kernel = tfp.mcmc.DualAveragingStepSizeAdaptation(
      inner_kernel=kernel,
      num_adaptation_steps=burnin,
      target_accept_prob=0.75,
      **step_adaptation_kwargs)

  # Execute sampling.
  chain_state, sampler_stat = tfp.mcmc.sample_chain(
      num_results=num_steps,
      num_burnin_steps=burnin,
      current_state=init_state,
      kernel=kernel,
      trace_fn=trace_fn,
      seed=seed)
    
  return chain_state, sampler_stat

def mix_chain_samples(samples: Union[tf.Tensor, List[tf.Tensor]]):
  """Mix MCMC samplers from different chains.
  
    Given a posterior sample with shape [num_sample, num_chain, ...], 
    collapse the samples from different chains by reshaping it as 
    [..., num_sample * num_chain].

  Args:
    samples: The posterior sample from multiple chains.

  Returns:
    The collapsed posterior samples.
  """
  if not isinstance(samples, list):
    samples = [samples]

  mixed_samples = []
  for sample in samples:
    sample_shape = list(tf.shape(sample).numpy())
    sample = tf.reshape(sample, [-1] + sample_shape[2:])
    sample = tf.transpose(sample, [1, 2, 0])

    mixed_samples.append(sample)

  return mixed_samples

# @title MCMC: Main API.
def run_mcmc(init_state: Optional[List[tf.Tensor]] = None,
             target_log_prob_fn: Optional[Callable[..., tf.Tensor]] = None, 
             model_dist: Optional[List[tfd.Distribution]] = None,      
             y: Optional[tf.Tensor] = None,
             sample_size: int = 500, 
             nchain: int = 10,             
             num_steps: int = 500, 
             burnin: int = 100, 
             step_size: float = .1, 
             seed: int = 0, 
             debug_mode: bool = False,
             **mcmc_kwargs):
  """Executes MCMC training for a given model posterior.
  
  Args:
    target_log_prob_fn: The log likelihood function of modle posterior.
      If not provided, then a default set of (init_state, target_log_prob_fn)
      will be generated by `prepare_mcmc`.    
    init_state: The initial states to the MCMC sampler, a list of tf.tensors 
      with shape (num_chains, num_variables). If not provided, then a default 
      set of (init_state, target_log_prob_fn) will be generated by 
      `prepare_mcmc`.
    model_dist: The model posterior distribution to be used by `prepare_mcmc`. 
      Must be provided if init_state or target_log_prob_fn is None.
    y: The output variable with shape (batch_size, ) to be used by 
      `prepare_mcmc`.  Must be provided if init_state or target_log_prob_fn is
       None.
    sample_size: The number of the final MCMC samples to return after thinning
      the gathered MCMC samples.
    n_chain: The number of MCMC chain in sampling.
    num_steps: The number of total MCMC samples to generate.
    burnin: The length of the burn-in period for MCMC warmup.
    seed: The random seed for MCMC sampling.
    debug_mode: If True. also return the original unmixed samples.
    **mcmc_kwargs: Additional keyword arguments to pass to the low-level MCMC
      function.

  Return:
    mixed_samples: A list of posterior samples with shapes [sample_size, ...]
      for each variable in the model posterior. 
    sampler_stat: diagnostic statistics of the MCMC chain, which contains 
      the step size and the proposal acceptance of each HMC step.
  """
  # Prepares initial states and model log likelihoods for MCMC.
  if init_state is None or target_log_prob_fn is None:
    # By default, sample first parameter of a two-parameter model (W, y).
    init_state, target_log_prob_fn = prepare_mcmc(
        model_dist, y, nchain=nchain)
    
  # Perform MCMC.
  chain_samples, sampler_stat = run_chain(
      init_state=init_state, 
      step_size=step_size,
      target_log_prob_fn=target_log_prob_fn,
      num_steps=num_steps, 
      burnin=burnin, 
      seed=seed,
      **mcmc_kwargs)
  # Clear tf.function cache.
  try:
    run_chain._stateful_fn._function_cache.clear()
  except:
    run_chain._stateful_fn._function_cache.primary.clear()

  # Thining.
  sample_size_per_chain = int(sample_size / nchain)
  sample_ids = np.linspace(
      0, num_steps-1, sample_size_per_chain).astype(int)
  chain_samples_thinned = chain_samples.numpy()[sample_ids]

  # Mix examples from different chains, 
  # Shape [param_dim_1, param_dim_2, num_mcmc_samples].
  mixed_samples = mix_chain_samples(chain_samples_thinned)

  # Check acceptance probability.
  p_accept = tf.math.exp(tfp.math.reduce_logmeanexp(
    tf.minimum(sampler_stat[-1], 0.)))
  print(f'Acceptance Ratio: {p_accept}')
  
  if debug_mode:
    return mixed_samples, chain_samples, sampler_stat
  return mixed_samples, sampler_stat

# @title MAP: Main API.
def run_map(target_log_prob_fn, gp_config,
            learning_rate=0.1, 
            num_steps=20_000, print_every=1_000, seed=None):
  """Executes MAP estimation using Adam optimizer."""
  # Prepares input variable.
  W_shape = (gp_config['hidden_units'], gp_config['units'])
  W_init = tf.initializers.HeNormal(seed=seed)(shape=W_shape)
  W_var = tf.Variable(initial_value=W_init)

  # Prepares loss function (negative log-likelihood).
  @tf.function
  def loss():
    nll = -target_log_prob_fn(W_var)
    return tf.reduce_mean(nll)

  # Runs optimization loop using Adam optimizer.
  opt = tf.keras.optimizers.Adam(learning_rate=learning_rate)

  for iter in range(num_steps):
    if iter % print_every == 0:
      print(f'{loss().numpy()}...', end='')
    _ = opt.minimize(loss, [W_var])
  print('Done.')

  return tf.constant(W_var)

## Utility Functions
# @title Visualization: plot_distribution
def plot_distribution(X, Y, loc=None, scale=None, X_train=None, Y_train=None, 
                      loc_color="black", scale_color="silver",
                      alpha_k=(1., 0.6, 0.3), interval_only=False):
  if isinstance(X, tuple):
    X, _ = X
  if isinstance(X_train, tuple):
    X_train, _ = X_train    
  
  x = X.squeeze()
  
  if not interval_only:
    plt.figure(figsize=(15, 5))

  if loc is not None and scale is not None:
    for k in (3, 2, 1):
      lb = (loc - k * scale).squeeze()
      ub = (loc + k * scale).squeeze()
      plt.fill_between(x, lb, ub, color=scale_color, alpha=alpha_k[3-k])
    plt.plot(x, lb, color=scale_color, alpha=0.5)
    plt.plot(x, ub, color=scale_color, alpha=0.5)
    plt.plot(X, loc, color=loc_color)

  if not interval_only:
    plt.scatter(X, Y, color="gray", alpha=0.8)
    if X_train is not None:
      plt.scatter(X_train, Y_train, color="black", alpha=0.6)
      
def gpr_2d_visual(pred_mean, pred_cov,
                  X_train, y_train, X_test, y_test,
                  title="", save_addr="", fontsize=12):
#     if save_addr:
#         pathlib.Path(save_addr).parent.mkdir(parents=True, exist_ok=True)
#         plt.ioff()

    # prediction surface
    n_reshape = int(np.sqrt(pred_mean.size))
    pred_mean_plot = pred_mean.reshape(n_reshape, n_reshape)
    X_valid = X_test.reshape(n_reshape, n_reshape, 2)
    x_grid, y_grid = X_valid[:, :, 0], X_valid[:, :, 1]

    ax = plt.axes(projection='3d')
    if isinstance(X_train, np.ndarray):
        ax.scatter(X_train[:, 0], X_train[:, 1], y_train, c="black")
    ax.plot_surface(X=x_grid, Y=y_grid, Z=pred_mean_plot, cmap='inferno')
    ax.set_zlim(np.min(y_test), np.max(y_test))

#     # optionally, compute RMSE
#     if pred_mean.size == y_test.size:
#         rmse = metric_util.rmse(y_test, pred_mean)
#         title = "{}, RMSE={:.4f}".format(title, rmse)

    plt.title(title, fontsize=fontsize)

# @title Visualization: plot_base_models
def plot_base_models(base_preds_test, kernel_names,
                     X_test, Y_test, X_train, Y_train, X_base, Y_base, 
                     ax=plt):
  for base_pred in tf.transpose(base_preds_test):
    ax.plot(X_test, base_pred)

  ax.scatter(X_train, Y_train, c='b', s=10)
  ax.scatter(X_base, Y_base, c='r', s=10)
  ax.scatter(X_test, Y_test, c='k', alpha=0.1)

  ax.legend(kernel_names)

# @title Data: get_base_prediction.
def get_base_prediction(X: np.ndarray, Y: np.ndarray, 
                        X_test: Optional[np.ndarray] = None, 
                        kernel: Optional[gpf.kernels.Kernel] = None,                         
                        model: Optional[gpf.models.BayesianModel] = None, 
                        num_train_steps: int = 200,
                        debug_mode: bool = False) -> Union[gpf.models.GPR, np.ndarray]:
  """Utility functions to train a base model.
  
  Args:
    X: Training data features.
    Y: Training data outcomes.
    X_test: Testing data to evaluate the model on. Default to the training data
      X.
    model: A trained GPFlow model.
    kernel: A GPFlow kernel that can be used to define a GPFlow Gaussian process.
      Cannot be None if `model = None`.
    num_train_steps: The number of gradient steps to optimize the Gaussian 
      process parameters.

  Returns:
    A trained GPFlow model if `model=None` and kernel is provided.
    Model prediction on X_test if model is provided.  
  """
  # Convert all data to float64 for GPFlow compatibility.
  X = X.astype(np.float64)
  Y = Y.astype(np.float64)
  X_test = X if X_test is None else X_test.astype(np.float64)

  if model is None:
    # Define model and train with `num_train_steps` number of steps.
    if debug_mode:
      print(f'Training GP with kernel `{kernel.name}`.')
    model = gpf.models.GPR(data=(X, Y), kernel=kernel)
    opt = tf.optimizers.Adam(1e-3)
    for _ in range(num_train_steps):
      opt.minimize(model.training_loss, model.trainable_variables)
    return model
  else:
    # Return prediction.
    if debug_mode:
      print(f'Predicting GP with kernel "{model.kernel.name}".')
    prediction = model.predict_f(X_test)[0]
    prediction = prediction.numpy().flatten().astype(np.float32)
    return prediction

# @title Data: process_mcmc_data

def process_mcmc_data(mcmc_sample, X, Y, op_type='flatten', debug_mode=False):
  """Transforms data to / from a MCMC-compatible format.
  
  Given MCMC samples from a pre-trained model (shape (num_train, ..., num_mcmc)) 
  and original data (X, Y) (shape [num_train, ...]), transforms them into a MCMC
  compatible format. Namely, the mcmc_sample will be reshaped to 
  (batch_size * num_mcmc, ...), and (X, Y) will be repeated at stacked to 
  (batch_size * num_mcmc, ...) as well.

  Alternatively, this function can also 'deflatten' a mcmc_sample with shape
  (num_mcmc_new, num_train * num_mcmc, ...) to 
  (num_mcmc_new * num_mcmc, num_train, ...).

  Args:
    mcmc_sample: MCMC samples from a pre-trained model. 
      If op_type='flatten', shape (num_train, ..., num_mcmc).
      If op_type='deflatten', shape (num_mcmc_new, num_train * num_mcmc, ...).
    X: Training data features. Shape (num_train, feature_dim) 
    Y: Training data response, shape (num_train, ...)
    op_type: Type of operation to ('flatten', 'deflatten')
    
  Returns:
    If op_type='flatten':

    sample_processed: shape (num_train * num_mcmc, ...)
    X_processed: shape (num_train * num_mcmc, feature_dim)
    Y_processed: shape (num_train * num_mcmc, ...)

    If op_type='deflatten':

    Reshaped MCMC samples from (num_mcmc_new, num_train * num_mcmc, ...) to
    (num_mcmc_new * num_mcmc, num_train, ...).
  """
  mcmc_sample = tf.convert_to_tensor(mcmc_sample, dtype=tf.float32)
  X = tf.convert_to_tensor(X, dtype=tf.float32)
  Y = tf.convert_to_tensor(Y, dtype=tf.float32)

  if op_type == 'flatten':
    return flatten_mcmc_data(mcmc_sample, X, Y, debug_mode)
  if op_type == 'deflatten':
    return deflatten_mcmc_sample(mcmc_sample, X, Y, debug_mode)


def flatten_mcmc_data(mcmc_sample, X, Y, debug_mode=False):
  """Merge input data to shapes (num_train*num_mcmc, ...).

  Args:
    mcmc_sample: shape (num_train, ..., num_mcmc).
    X: shape (num_train, feature_dim)
    Y: shape (num_train, ...)
  
  Returns:
    sample_processed: shape (num_train*num_mcmc, ...)
    X_processed: shape (num_train*num_mcmc, feature_dim)
    Y_processed: shape (num_train*num_mcmc, ...)
  """
  num_train = mcmc_sample.shape.as_list()[0]
  num_mcmc =  mcmc_sample.shape.as_list()[-1]
  output_dims = mcmc_sample.shape.as_list()[1]
  num_train_X = X.shape.as_list()[0]
  num_train_Y = Y.shape.as_list()[0]
  output_dims_Y = Y.shape.as_list()[-1]

  assert num_train == num_train_X, f'num_train={num_train}, num_train_X={num_train_X}'
  assert num_train == num_train_Y, f'num_train={num_train}, num_train_Y={num_train_Y}'
  assert output_dims == output_dims_Y, f'output_dims={output_dims}, output_dims_Y={output_dims_Y}'

  if debug_mode:
    print(f'Process data with num_mcmc={num_mcmc}, num_train={num_train}, '
          f'output_dims={output_dims}')
  
  # Process MCMC sample, row-major stacking (contiguous predictions per X).
  mcmc_sample_ = tf.transpose(mcmc_sample, [0, 2, 1])  # (num_train, num_mcmc, ...).
  mcmc_sample_ = tf.reshape(mcmc_sample_, [-1, output_dims])  # (num_train * num_mcmc, ...).

  # Process X, row-major stacking from (num_train, orig_dims).
  X_ = tf.repeat(X, num_mcmc, axis=0)  # (num_train * num_mcmc, feature_dims)
  Y_ = tf.repeat(Y, num_mcmc, axis=0)  # (num_train * num_mcmc, output_dims)

  return mcmc_sample_, X_, Y_


def deflatten_mcmc_sample(mcmc_sample_new, X, unused_Y, debug_mode=False):
  """Reshapes new MCMC sample to (num_mcmc, num_train, ...).
  
  Given data with shape (num_train * num_mcmc, ...), the model 
  produces a sample of (num_mcmc_new, num_train*num_mcmc, ...). This function
  averages over the num_mcmc_new dimension and returns a tensor with shape 
  (num_mcmc * num_mcmc_new, num_train, ...).

  Args:
    mcmc_sample_new: New MCMC samples of shape 
      (num_mcmc_new, num_train * num_mcmc, ...).
    X: shape (num_train, feature_dim)
    unused_Y: shape (num_train, ...)
  
  Returns:
    Reshaped MCMC samples of shape (num_mcmc_new*num_mcmc, num_train, ...).
  """
  del unused_Y
  num_mcmc_new, num_total = mcmc_sample_new.shape.as_list()[:2]
  num_train = X.shape.as_list()[0]
  num_mcmc = int(num_total / num_train)

  if debug_mode:
    print(f'Process data with num_total={num_total}, num_train={num_train}, '
          f'num_mcmc_new={num_mcmc_new}, num_mcmc={num_mcmc}')

  # De-flatten back to original dimension, (num_mcmc_new, num_train, num_mcmc, ...).
  mcmc_sample_ = tf.reshape(mcmc_sample_new, [num_mcmc_new, num_train, num_mcmc, -1])

  # Change axis to original, (num_mcmc_new, num_mcmc, num_train, ...).
  mcmc_sample_ = tf.transpose(mcmc_sample_, [0, 2, 1, 3])

  # Collapse over MCMC dimension.
  return tf.reshape(mcmc_sample_, [num_mcmc_new * num_mcmc, num_train, -1])

# @title Model: prepare_mcmc
def prepare_mcmc(
    model_prior, y, param_name='gp_weights', nchain=10, seed=0):
  init_state = model_prior.sample(nchain)[param_name]

  def target_log_prob_fn(gp_w):
    # conditional_sample = model_prior.sample(batch_size)
    # Only perform inference on gp_w.
    conditional_sample = dict()
    conditional_sample[param_name] = gp_w
    conditional_sample['y'] = y

    log_probs = model_prior.log_prob(conditional_sample)
    return tf.reduce_mean(log_probs)

  return init_state, target_log_prob_fn

# @title Model: make_bma_samples.
def make_bma_samples(X, Y, base_model_preds, 
                     bma_weight_samples, 
                     bma_model_config, n_samples=50, seed=0, 
                     y_samples_only=True,
                     prepare_mcmc_training=False,
                     debug_mode=False):
  # Defines posterior distribution.
  bma_posterior, _ = bma_dist(X, base_model_preds,
                              posterior_mode=True,
                              posterior_sample=bma_weight_samples,
                              sample_intermediate_variables=not y_samples_only,
                              debug_mode=debug_mode,
                              **bma_model_config)  

  # Samples from posterior. Shape (num_mcmc, num_data, num_output).
  bma_samples = bma_posterior.sample(n_samples, seed=seed)

  if y_samples_only:
    # Only returns samples for y, shape (num_mcmc, num_data, num_output).
    bma_samples = bma_samples["y"]

    if prepare_mcmc_training:
      # Reshapes to (num_data, num_output, num_mcmc).
      bma_samples = tf.transpose(bma_samples, [1, 2, 0])

      # Reshape to (num_data * num_mcmc, ...).
      return process_mcmc_data(bma_samples, X, Y)

  del bma_posterior
  return bma_samples

# @title Model: make_bne_samples
def make_bne_samples(X: tf.Tensor, 
                     mean_preds: tf.Tensor,
                     bne_weight_samples: tf.Tensor, 
                     bne_model_config: Dict[str, Any], 
                     seed: int = 0, 
                     debug_mode: bool = False):
  """Generates posterior predictive samples for BNE.
  
  Args: 
    X: Input features, shape (num_data, num_input_dim).
    mean_preds: Predictive samples from a basic ensemble model, 
      shape (num_pred_samples, num_data, 1).
    bne_weight_samples: Posterior samples for random feature weights of the BNE
      model, shape (num_hidden_units, num_param_dim, num_bne_samples). The 
      `num_bne_samples` here is the number of posterior sample for the BNE's 
      internal parameters and does not need to equal to `num_pred_samples` from
       `mean_preds`.
    bne_model_config: Keyword arguments to `bne_model_dist`.
    seed: Random seed for posterior sampler.
    debug_mode: Whether to print out intermediate debugging information.

  Returns:
    bne_samples: A dictionary of BNE posterior predictions, with fields:
     * 'bne_gp_w': random feature weights, 
        (num_pred_samples, num_hidden_units, num_params)
     * 'bne_gp': GP posteriors evaluated at X, shape 
        (num_pred_samples, num_params)
     * 'mean': the mean parameter for the residual process,
       (num_pred_samples, 1)
     * 'scale': the scale (variance) parameter 
        (num_pred_samples, 1)
     * 'rate': the rate parameter for modeling skewness
        (num_pred_samples, 1). 
        Available only when `estimate_skewness=True` in bne_model_dist.
     * 'skewness': the skewness parameter computed as 
        tfd.Exponential(rate) (num_pred_samples, 1). Available only when 
        `estimate_skewness=True` in bne_model_dist.
     * 'mean_original': The predictive distribution from the original ensemble
        model, (num_pred_samples, 1).
     * 'y': The predictive distribution of the final model ensmeble, 
        (num_pred_samples, 1).
  """
  # Identify the number of predictive samples to create.
  num_pred_samples = mean_preds.shape[0]

  if debug_mode:
    print(
        f'Generate {num_pred_samples} predictive samples.'
        f'based on mean_preds with shape {mean_preds.shape}')
  
  # Construct the residual posterior. 
  bne_posterior, _ = bne_model_dist(
      inputs=X,
      mean_preds=0.,  # Only generate residual sample by fixing mean at zero.
      posterior_mode=True,
      posterior_sample=bne_weight_samples,
      sample_intermediate_variables=True,
      debug_mode=debug_mode,
      **bne_model_config)

  if debug_mode:
    print(f'BNE Posterior Model Graph: {bne_posterior.resolve_graph()}')

  # Sample from residual posterior.
  bne_samples = bne_posterior.sample(num_pred_samples, seed=seed)
  
  # Construct full model prediction by adding mean. Shape (num_sample, num_data, 1)
  bne_samples['mean_original'] = mean_preds
  bne_samples['resid'] = bne_samples['y']  
  bne_samples['y'] = bne_samples['resid'] + mean_preds  

  del bne_posterior  
  return bne_samples

## Simulation Functions
# @title Simulation: get_data
def get_data(seed, N_train, y_dist, data_gen_fn, num_train_steps=100):
  """Generates data and trains base models."""
  (X_base, X_train, X_test, 
   Y_base, Y_train, Y_test, 
   mean_test) = data_gen_fn(N_train, N_base, N_test, 
                            y_dist=y_dist, seed=seed)

  base_preds_train, base_preds_test, kernel_names = run_base_models(
      X_base, X_train, X_test, Y_base, Y_train, Y_test, 
      num_train_steps=num_train_steps)
  
  data_dict = dict(X_base=X_base,
                   X_train=X_train,
                   X_test=X_test,
                   Y_base=Y_base, 
                   Y_train=Y_train, 
                   Y_test=Y_test, 
                   mean_test=mean_test, 
                   base_preds_train=base_preds_train, 
                   base_preds_test=base_preds_test, 
                   base_model_names=kernel_names)
  return data_dict

# @title Simulation: get_bma_result
def get_bma_result(data_dict, bma_config):
  """Trains Adaptive Bayesian model averaging."""
  (bma_joint_samples, X_train_mcmc, Y_train_mcmc, 
   means_train_mcmc, means_test_mcmc) = run_bma_model(
       X_train=data_dict["X_train"], 
       X_test=data_dict["X_test"], 
       Y_train=data_dict["Y_train"], 
       base_preds_train=data_dict["base_preds_train"], 
       base_preds_test=data_dict["base_preds_test"], 
       return_mcmc_examples=True,
       **bma_config)

  data_dict['X_train_mcmc'] = X_train_mcmc
  data_dict['Y_train_mcmc'] = Y_train_mcmc
  data_dict['means_train_mcmc'] = means_train_mcmc
  data_dict['means_test_mcmc'] = means_test_mcmc
  data_dict['bma_mean_samples'] = bma_joint_samples['y']

  return data_dict

# @title Simulation: get_bne_result
def get_bne_result(data_dict, moment_mode, bne_config):
  """Trains Bayesian nonparametric ensemble."""
  mode_to_name_map = {'none': 'bma', 'mean': 'bae', 
                      'variance': 'bne_var', 'skewness': 'bne_skew'}
  model_name = mode_to_name_map[moment_mode]

  joint_samples = run_bne_model(X_train=data_dict['X_train_mcmc'], 
                                Y_train=data_dict['Y_train_mcmc'], 
                                X_test=data_dict['X_test'],
                                base_model_samples_train=data_dict['means_train_mcmc'],
                                base_model_samples_test=data_dict['means_test_mcmc'],
                                moment_mode=moment_mode,
                                **bne_config)
  
  data_dict[f'{model_name}_samples'] = joint_samples['y']
  return data_dict


















