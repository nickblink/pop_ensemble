## libraries
from typing import Any, Callable, Dict, List, Optional, Union, Tuple

import os
import gc
import time
import pickle
import functools
import scipy
#import pyreadr

import multiprocessing as mp

# from google.colab import files
# from google.colab import 

import numpy as np
import tensorflow as tf

import pandas as pd

import plotly.express as px

import matplotlib.pyplot as plt
import seaborn as sns
import edward2 as ed
import tensorflow_probability as tfp

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

tfd = tfp.distributions
tfb = tfp.bijectors

dtype = tf.float64
import gpflow as gpf
import logging

from sklearn.model_selection import KFold 
from sklearn.linear_model import LinearRegression

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

# simulate the data
def simulate_data(data_OG, adjacency, pivot = -1, sim_numbers = False, scale_down = 1, poisson_noise = False, one_model = False, models = ['acs', 'pep', 'worldpop']):
    """Simulated data for the CAR model. 
    
    Args:
        data: The input data file.
        adjacency: the adjacency matrix for the data.
        pivot: The column to be used as the pivot for the softmax ensemble weights. -1 indicates no pivot.
        sim_numbers: whether to simulate data values for the models. If false, the true values are used.
        scale_down: value to scale down the models data by.
        poisson_noise: Whether to add noise to input models by simulating a poisson from each
        one_model: whether to only have one model determine the output. Worldpop is the default if this is chosen.
    """
    # make copy of day so that the original data is not corruptedb
    data = data_OG[:]
    
    if sim_numbers:
        #data['acs'] = np.random.normal(80.0, 10.0, data.shape[0])
        #data['pep'] = np.random.normal(100.0, 10.0, data.shape[0])
        #data['worldpop'] = np.random.normal(120.0, 10.0, data.shape[0])
        data['pep'] = np.random.normal(100.0, 10.0, data.shape[0])
        data['worldpop'] = np.random.normal(data['pep'], 10.0, data.shape[0])
    
    # scale down the values of the data 
    data[models] = data[models][:]/scale_down
    
    if poisson_noise:
        for m in models:
            data[m] = np.random.poisson(data[m])
        
    if one_model:
        data['census'] = np.random.poisson(data['worldpop'].to_numpy())
        
        return _, _, data
    else:
        tau2 = 1
        rho = 0.3
        print('fixing tau2 and rho')
        nchain = 1
        
        Q = (1/tau2)*(np.diag(adjacency.sum(axis=1)) - rho*adjacency)
        Q = tf.constant(Q, dtype = tf.float64)

        if(pivot == -1):
            phi_true = tf.constant(np.array([mv_normal_sample(precision_matrix = Q, 
                                                          num_models = len(models)) for i in range(nchain)]),
                               dtype = tf.float64)

        elif(pivot in range(len(models))):
            nm = len(models) - 1 
            phi_np = np.array([mv_normal_sample(precision_matrix = Q,
                                          num_models = nm) for i in range(nchain)])
            # fix the added dimension if one is dropped
            if nm == 1:
                phi_np = phi_np[:,:,np.newaxis]
            phi_np = np.insert(phi_np, pivot, 0., axis = 2)
            phi_true = tf.constant(phi_np, dtype = tf.float64)

        else:
            raise Exception('Pivot needs to be -1, 0, 1, or 2')

        # get exponentiated values and sum across models
        exp_phi = tf.math.exp(phi_true)
        exp_phi_rows = tf.reduce_sum(exp_phi, 2)

        # get model weights and calculate mean estimate
        u_true = exp_phi/exp_phi_rows[...,None]

        tmp = data[models].values*u_true
        n = tf.reduce_sum(tmp, axis = 2)
        
        data['census_exp'] = n[0]

        #data['census'] = np.random.poisson(n)[0]
        data.loc[:,'census'] = np.random.poisson(n)[0]
        
        # add in the u values
        u_df = pd.DataFrame(u_true[0].numpy(), columns = ['u_' + m for m in models])
        data = pd.concat([data.reset_index(drop = True)[:], u_df], axis = 1)

    return phi_true, u_true, data

# subset data by state
def subset_data_by_state(data, adjacency, state, abbrev = None):
    str_vals = data['NAME'].str.find(state)
    indices = [i for i in range(len(data)) if str_vals[i] > -1]
    
    str_vals2 = adjacency.index.str.find(state)
    indices2 = [i for i in range(len(adjacency)) if str_vals2[i] > -1]
    
    if not indices == indices2:
        raise Exception('indices of the data names and adjacency names do not match')
        
    if abbrev != None:
        str_vals3 = adjacency.columns.str.find(abbrev)
        indices3 = [i for i in range(len(adjacency)) if str_vals3[i] > -1]
        if not indices == indices3:
            missing = [i for i in range(len(indices)) if indices[i] != indices3[i]]
            #print(missing)
            #print(indices[29:31])
            #print(indices3[29:31])
            raise Exception('indices of the data names and adjacency column names do not match')

    data2 = data.iloc[indices, :]
    adjacency2 = adjacency.iloc[indices, indices]
    return data2, adjacency2

# Convert phi values to u ensemble weights by computing the softmax
# If pivot is non-negative, first introduce the pivot into the data

def phi_to_u(phi, pivot = -1):
    # get the index of last dimension
    dim_n = len(phi.shape) - 1
    
    # add in the pivot
    if pivot > -1:
        phi = tf.constant(np.insert(phi, pivot, 0., axis = dim_n))
    
    # exponentiate values
    exp_phi = tf.math.exp(phi)
    
    # get sums of exponentiated phi values at each index
    exp_phi_rows = tf.reduce_sum(exp_phi, dim_n)
    
    # get the model weights by taking the softmax!
    u = exp_phi/exp_phi_rows[...,None]
    
    return u
    
    
def get_log_prob_from_results(res_dict):
    """ Gets the target log probability function from a set of results
    """
    _, target_log_prob_fn = prepare_mcmc_CAR(data = res_dict['data'],
    adjacency = res_dict['adjacency'], 
    nchain = res_dict['mcmc_config']['nchain'],
    pivot = res_dict['pivot_fit'],
    models = res_dict['models'],
    run_MAP = False)
    
    return(target_log_prob_fn)
    
##### CAR functions

# This function was taken from online
# Generate samples from a multi-variate normal distribution with provided precision matrix WITHOUT inverting
def mv_normal_sample(mu=0, precision_matrix=None, num_models=1):

    # Precision matrix must be a square matrix
    assert precision_matrix.shape[0] == precision_matrix.shape[1], 'Precision matrix must be a square matrix'

    dim = precision_matrix.shape[0]

    chol_U = scipy.linalg.cholesky(precision_matrix, lower=False)

    # Create num_models iid standard normal vectors
    z_vector_matrix = np.random.normal(loc=0, scale=1, size=[num_models, dim])

    # Sample from the MV normal with precision matrix by solving the Cholesky decomp for each normal vector
    samples = np.squeeze(np.array(
        [scipy.linalg.solve_triangular(a=chol_U, b=z_vector_matrix[i, :], unit_diagonal=False) + mu for i in
         range(num_models)]))

    return (np.transpose(samples))
    
# @title MAP: Main API.
def run_map_CAR(target_log_prob_fn, 
            init_state,
            learning_rate=0.1, 
            num_steps=20_000, 
            print_every=1_000, 
            seed=None,
            return_all = False):
  """Executes MAP estimation using Adam optimizer."""
  # Prepares input variable.
  phi = tf.Variable(initial_value=init_state)

  # Prepares loss function (negative log-likelihood).
  @tf.function
  def loss():
    nll = -target_log_prob_fn(phi)
    return tf.reduce_mean(nll)

  # Runs optimization loop using Adam optimizer.
  opt = tf.keras.optimizers.Adam(learning_rate=learning_rate)

  if return_all:
    loss_list = []
    for iter in range(num_steps):
      if iter % print_every == 0:
        print(iter)
        print(f'{loss().numpy()}...', end='')
        _ = opt.minimize(loss, [phi])
      loss_list.append(target_log_prob_fn(phi).numpy())
    print('Done.')
    
    return loss_list
   
  else:
    for iter in range(num_steps):
      if iter % print_every == 0:
        print(iter)
        print(f'{loss().numpy()}...', end='')
        _ = opt.minimize(loss, [phi])
    print('Done.')

    return tf.constant(phi)

#### Makes the MCMC sampler (with the log probability function)
def prepare_mcmc_CAR(data, 
                     adjacency, 
                     nchain, 
                     pivot, 
                     models = ['acs', 'pep', 'worldpop'],                      
                     run_MAP = True,
                     map_config: Optional[Dict[str, Any]] = None):
  """prepares the initial state and log prob function"""
  tau2 = 1
  rho = 0.3
  print('fixing tau2 and rho')
  print('when adding in tau2 and rho, need to update the likelihood function!')

  Q = (1/tau2)*(np.diag(adjacency.sum(axis=1)) - rho*adjacency)
  Q = tf.constant(Q, dtype = tf.float64)

  # define log likelihood function
  def target_log_prob_fn_CAR(phi, debug_return = False):
    #Q = (1/tau2)*(np.diag(adjacency.sum(axis=1)) - rho*adjacency)
    #Q = tf.constant(Q, dtype = tf.float64)
        
    ll = tf.Variable(0., dtype = tf.float64)
    A = tf.Variable(0., dtype = tf.float64)
    B = tf.Variable(0., dtype = tf.float64)
    C = tf.Variable(0., dtype = tf.float64)
    
    for chain in range(phi.shape[0]):
        # (1) Prob of the CAR random effect values
        A = A - 0.5*tf.reduce_mean(tf.linalg.diag_part(
            tf.linalg.matmul(phi[chain,:,:],tf.linalg.matmul(Q, phi[chain,:,:]), transpose_a = True))) 
    ll = ll + A
    
    # add in determinant values
    # these are multiplied by the number of chains because they are included in the likelihood for each
    log_det = tf.constant(np.linalg.slogdet(Q)[1], dtype = tf.float64)
    B = 0.5*phi.shape[0]*len(models)*log_det
    #log_det = tf.linalg.logdet(Q)[1], dtype = tf.float64
    ll = ll + B 
    
    if(pivot == -1):
        # get exponentiated values and sum across models
        exp_phi = tf.math.exp(phi)
        exp_phi_rows = tf.reduce_sum(exp_phi, 2)
    elif(pivot in range(3)):
        phi_np = phi.numpy()
        phi_np = np.insert(phi_np, pivot, 0., axis = 2)
        exp_phi = tf.math.exp(tf.constant(phi_np))
        exp_phi_rows = tf.reduce_sum(exp_phi, 2)
    else:
        raise Exception('Pivot needs to be -1, 0, 1, or 2')
    
    # get model weights and calculate mean estimate
    u = exp_phi/exp_phi_rows[...,None]
      
    tmp = data[models].values*u
    n = tf.reduce_sum(tmp, axis = 2)
    
    # update the log likelihood 
    C = tf.reduce_sum([np.sum(data['census']*np.log(n[chain,:]) - n[chain,:]) for chain in range(phi.shape[0])])
    ll = ll + C
    
    if debug_return:
        dct = {'phi': A, 
               'det': B,
               'Poisson': C,
               'total': ll}
        return(dct)
    else:
        return(ll)  

  if(pivot == -1):
    nm = len(models)
  else:
    nm = len(models) - 1
  init_state = tf.constant(np.array([mv_normal_sample(precision_matrix = Q, num_models = nm) for i in range(nchain)]), dtype = tf.float64)
  
  # adding in an extra dimension
  if nm == 1:
     init_state = init_state[:,:,np.newaxis]

  if run_MAP:
    print('running MAP')
    init_state = run_map_CAR(target_log_prob_fn_CAR, init_state)

  return init_state, target_log_prob_fn_CAR


# @tf.function(experimental_compile=True)
def run_mcmc_CAR(init_state: Optional[List[tf.Tensor]] = None,
             target_log_prob_fn: Optional[Callable[..., tf.Tensor]] = None, 
             data: Optional[pd.DataFrame] = None,      
             adjacency: Optional[pd.DataFrame] = None,
             #models: Optional[List] = None,
             #data: Optional = None,
             #adjacency: Optional = None,             
             pivot: int = -1,
             models: Optional = None,             
             sample_size: int = 500, 
             nchain: int = 10,             
             num_steps: int = 500, 
             burnin: int = 100, 
             num_adaptation_steps: int = None,
             step_size: float = .1, 
             run_MAP = True,
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
    init_state, target_log_prob_fn = prepare_mcmc_CAR(data, adjacency, nchain, pivot, models, run_MAP)
  else:
    nchain = init_state.shape[0]
    
  # Perform MCMC.
  chain_samples, sampler_stat = run_chain_CAR(
      init_state=init_state, 
      step_size=step_size,
      target_log_prob_fn=target_log_prob_fn,
      num_steps=num_steps, 
      burnin=burnin, 
      num_adaptation_steps=num_adaptation_steps,
      seed=seed,
      **mcmc_kwargs)
  # Clear tf.function cache.
  try:
    try:
      run_chain._stateful_fn._function_cache.clear()
    except:
      run_chain._stateful_fn._function_cache.primary.clear()
  except:
    print('no cache clearing')

  # Thinning.
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

def run_chain_CAR(init_state: List[tf.Tensor], 
              step_size: float, 
              target_log_prob_fn: Callable[..., tf.Tensor], 
              num_steps: int = 500, 
              burnin: int = 100, 
              seed: int = 0,
              kernel_type: str = "hmc",
              step_adaptor_type: str = "simple",
              num_adaptation_steps: int = None
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
  print('kernel type is ' + kernel_type)
  
  if kernel_type not in ('hmc', 'nuts'):
    raise ValueError(
        f"kernel_type {kernel_type} must be one of ('hmc', 'nuts').")

  if step_adaptor_type not in ('simple', 'dual_averaging', 'none'):
    raise ValueError(
        f"step_adaptor_type {step_adaptor_type} must be one of "
        "('simple', 'dual_averaging', 'none').")

  if num_adaptation_steps is None:
    num_adaptation_steps = int(burnin  * 0.8)

  def trace_fn(_, pkr): 
    if step_adaptor_type == 'none':
      step_size = 0.1
    elif kernel_type == 'hmc':
      step_size = pkr.inner_results.accepted_results.step_size
    else:
      step_size = pkr.inner_results.step_size

    return (step_size, pkr.inner_results.log_accept_ratio)

  if kernel_type == 'hmc':
    kernel = tfp.mcmc.HamiltonianMonteCarlo(
        target_log_prob_fn=target_log_prob_fn,
        num_leapfrog_steps=3,
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

  if step_adaptor_type == 'simple':
    print('simple step size')
    kernel = tfp.mcmc.SimpleStepSizeAdaptation(
      inner_kernel=kernel, 
      num_adaptation_steps= num_adaptation_steps,
      #num_adaptation_steps = 1000,
      target_accept_prob=0.7)
  elif step_adaptor_type == 'dual_averaging':
    print('dual averaging step size')
    kernel = tfp.mcmc.DualAveragingStepSizeAdaptation(
      inner_kernel=kernel,
      num_adaptation_steps= num_adaptation_steps,
      target_accept_prob=0.7,
      **step_adaptation_kwargs)
  else:
    print('no step adaptor. Step size = ' +  str(step_size))

  print(num_steps)
  print(trace_fn)
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


def pull_gradient(phis, log_prob_fn, skip_val = 100, max_iter = None):
    """ Given a set of phis and the log probability function, pull the log likelihood values and gradients
    Args:
        phis: A set of phi values from one simulation run
        log_prob_fn: The log probability function
        skip_val: Skipping between likelihood and gradient calculation
        max_iter: whether to stop after a certain point
        
    Returns:
        the iterations calculated, the likelihoods, and the gradients
    """
    if max_iter is None:
        max_iter = phis.shape[0]
    iter = 0
    iter_counts = []
    posteriors = []
    phi_likelihoods = []
    det_likelihoods = []
    pois_likelihoods = []
    gradients = list()
    while iter < max_iter:
        # store the iter counts
        iter_counts.append(iter)
        
        # get the phi values at this chain sample
        phi = phis[iter,:,:,:]
        
        # get the likelihoods
        tmp = log_prob_fn(phi, debug_return = True)
        phi_likelihoods.append(tmp['phi'].numpy())
        det_likelihoods.append(tmp['det'].numpy())
        pois_likelihoods.append(tmp['Poisson'].numpy())
        posteriors.append(tmp['total'].numpy())
        
        # get the gradients
        with tf.GradientTape() as g:
          g.watch(phi)
          y = log_prob_fn(phi)
        gradients.append(g.gradient(y, phi))
        
        iter = iter + skip_val
        
    return iter_counts, posteriors, gradients, phi_likelihoods, det_likelihoods, pois_likelihoods



def pull_gradient_wrapper(phis, log_prob_fn, skip_val = 100, max_iter = None, step_sizes = None):
    """ A wrapper function to call "pull gradient"
    Args:
        phis: A set of phi values from one simulation run
        log_prob_fn: The log probability function
        skip_val: Skipping between likelihood and gradient calculation
        max_iter: whether to stop after a certain point
        step_sizes: Input series of step sizes

    Returns:
        A data frame with the iterations, log likelihoods, and average phi gradients across all chains
    """
    res = pull_gradient(phis, log_prob_fn, skip_val = skip_val, max_iter = max_iter)
    
    chain_abs_grads = []
    for i in range(len(res[2])):
        tt = res[2][i].numpy()
        chain_abs_grads.append(np.mean(abs(tt), axis = (1,2)))
    
    # organize the gradients by chain and then average across
    chain_abs_grads = np.array(chain_abs_grads)
    all_abs_grads = np.mean(chain_abs_grads, axis = 1)
    
    # get the step sizes at the iterations selected (not the most efficient python way)
    if step_sizes is not None:
        step_subset = []
        for i in range(len(res[0])):
            step_subset.append(step_sizes[res[0][i]])
    
    # only including total mean gradients rather than individual chains
    res_df = pd.DataFrame(np.transpose(np.array([res[0], res[1], all_abs_grads, step_subset, res[3], res[4], res[5]])), columns = ['iter', 'log_posterior', 'mean_abs_grad', 'step_size', 'logL_CAR_prior', 'logL_det', 'log_likelihood'])

    return(res_df)
    