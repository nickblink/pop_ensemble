#!/usr/bin/env python
# coding: utf-8

# # Default Configs

# In[1]

import pdb
import os
# os.getcwd()
# set cwd to the pop_ensemble directory
os.chdir('C:/Users/nickl/Documents/pop_ensemble')

exec(open("code/BNE_code/bne_0_a_load_modules.py").read())
exec(open("code/BNE_code/bne_0_b_wrapper_functions.py").read())
exec(open("code/BNE_code/bne_0_c_core_model_functions.py").read())
exec(open("code/BNE_code/bne_0_d_utility_functions.py").read())
exec(open("code/BNE_code/bne_0_e_simulation_functions.py").read())


# In[46]:


# GP configs.
y_noise_std = 0.1  # @param
hidden_units = 128  # @param
lengthscale=1.  # @param
l2_regularizer=0.1  # @param

DEFAULT_GP_CONFIG = dict(lengthscale=lengthscale,
                         l2_regularizer=l2_regularizer, 
                         hidden_units=hidden_units, 
                         y_noise_std=y_noise_std)


# In[47]:


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


# In[48]:


# MAP configs.
map_step_size=0.1 # @param
map_num_steps=10_000 # @param

DEFAULT_MAP_CONFIG = dict(learning_rate=map_step_size,
                          num_steps=map_num_steps)


# In[49]:


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


# # Full Example: Time Series with BMA, BAE, and BNE

# ## Data

# In[60]:


N_train = 50  # @param
N_base = 50 # @param  
N_test = 1000  # @param
seed = 0  # @param


# In[67]:


# STR: added mean_test as output
(X_base, X_train, X_test, Y_base, Y_train, Y_test, mean_test) = generate_data_1d(N_train, N_base, N_test, 
                                             y_dist=np.random.lognormal,
                                             seed=seed)


# In[68]:


plot_distribution(X_test, Y_test, loc=None, scale=None,  
                  X_train=X_train, Y_train=Y_train)
plt.ylim(-1, 17.5)
plt.show()


# ## Base Ensemble Models

# In[69]:


# Specify base models.
kernels = [gpf.kernels.Matern52(lengthscales=0.5), 
           gpf.kernels.Polynomial(degree=3.), 
           gpf.kernels.ArcCosine(weight_variances=1., bias_variance=1.),
           gpf.kernels.Periodic(gpf.kernels.Exponential(lengthscales=2.))]

n_models = len(kernels)
kernel_names = [kernel.name for kernel in kernels]


# In[70]:

models = [
  get_base_prediction(X_base, Y_base, X_train, kernel=k) for k in kernels]

base_preds_train = tf.stack([
  get_base_prediction(X_base, Y_base, X_train, model=m) for m in models], axis=-1)

base_preds_test = tf.stack([
  get_base_prediction(X_base, Y_base, X_test, model=m) for m in models], axis=-1)


# In[102]:


# something is wrong with the plots

plt.figure(figsize=(12, 6))

plot_base_models(base_preds_test, kernel_names,
                 X_test, Y_test, X_train, Y_train, X_base, Y_base)

plt.ylim([-1, 15.])
plt.title("Base Model Predictions")
plt.show()


# ## Bayesian Model Averaging

# A Bayesian ensemble model where ensemble weights $w_k's$ are parameterized by Gaussian process priors:
# 
# $y \sim N(\mu(x), \sigma^2)$ 
# 
# $\mu(x) = \sum_{k=1}^K w_k(x) * m_k(x) \quad$  where $\{m_k\}_{k=1}^K$ are base model predictions.
# 
# $w(x) = softmax(f(x)) \qquad\;\;\;$ where $w=[w_1, \dots, w_K]$ and $f=[f_1, \dots, f_K]$
# 
# $f \stackrel{i.i.d.}{\sim} GaussianProcess(0, k)$
# 
# 

# In[103]:


# Model configs.
y_noise_std = 0.1  # @param
lengthscale=1.  # @param
l2_regularizer=0.1  # @param


# In[104]:


# MCMC configs.
map_step_size=0.1 # @param
map_num_steps=10_000 # @param

mcmc_step_size=0.1 # @param
mcmc_num_steps=10_000 # @param


# In[105]:


# Posterior configs.
bma_n_samples_train = 100 # @param
bma_n_samples_test = 200 # @param
bma_n_samples_eval = 1000  # @param

bma_seed = 0  # @param
bne_seed = 0 # @param


# In[106]:


# Assemble into configs.
bma_model_config = DEFAULT_GP_CONFIG.copy()
map_config = DEFAULT_MAP_CONFIG.copy()
mcmc_config = DEFAULT_MCMC_CONFIG.copy()

bma_model_config.update(dict(lengthscale=lengthscale,
                             l2_regularizer=l2_regularizer,
                             y_noise_std=y_noise_std))

map_config.update(dict(learning_rate=map_step_size,
                       num_steps=map_num_steps))

mcmc_config.update(dict(step_size=mcmc_step_size, 
                        num_steps=mcmc_num_steps))


# ### Build Model

# In[107]:


bma_prior, bma_gp_config = bma_dist(X_train, 
                                    base_preds_train, 
                                    **bma_model_config)

bma_model_config.update(bma_gp_config)


# In[108]:


# Check if the model graph is specified correctly.
bma_prior.resolve_graph()
# what is this supposed to look like?

# ### Run MCMC

# In[109]:

# list of dim [1][128]. The second number has something to do with hidden
# GP layers, which I don't really get here, but it's some approximation of the GP.
bma_gp_w_samples = run_posterior_inference(model_dist=bma_prior, 
                                           model_config=bma_model_config,
                                           Y=Y_train, 
                                           map_config=map_config,
                                           mcmc_config=mcmc_config)


# In[110]:


bma_joint_samples = make_bma_samples(X_test, Y_test, base_preds_test, 
                                     bma_weight_samples=bma_gp_w_samples[0],
                                     bma_model_config=bma_model_config, 
                                     n_samples=bma_n_samples_eval, 
                                     seed=bma_seed,
                                     y_samples_only=False)

bma_ensemble_weights = bma_joint_samples['ensemble_weights']
bma_y_samples = bma_joint_samples['y']


# ### Inspect Posterior

# In[111]:


# Plot has issue 
# Plot posterior prediction.
plt.figure(figsize=(12, 6))

plot_base_models(base_preds_test, kernel_names,
                 X_test, Y_test, X_train, Y_train, X_base, Y_base)

#for sample in bma_y_samples:
for i in range(100):
    sample = bma_y_samples[i]
    #pdb.set_trace()
    plt.plot(X_test, sample, alpha=0.01, c='k')

plt.ylim([-1, 15])
plt.title('Posterior Predictive, Bayesian Model Averaging')
plt.show()


# In[112]:


# Plot posterior ensemble weights for different models.
fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(20, 25))

for model_id, model_name in enumerate(kernel_names):
  ax = axes[model_id]
  color = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red'][model_id]

  plot_base_models(base_preds_test, kernel_names,
                   X_test, Y_test, X_train, Y_train, X_base, Y_base,
                   ax=ax)

  #samples = bma_ensemble_weights[:, :, model_id]
  samples = bma_ensemble_weights[:100, :, model_id]
  sample_mean = tf.reduce_mean(samples, axis=0)

  # Plot posterior predictive samples of ensemble weights.
  # Multiply by 10 so it visually more stands out.
  for sample in samples:
    ax.plot(X_test, sample * 10, alpha=0.005, c=color)
  
  # Plot posterior mean of ensemble weights.
  ax.plot(X_test, sample_mean * 10, c=color, 
          linestyle='--', linewidth=3)

  ax.set_ylim([-1., 15.])
  ax.set_title(f'Ensemble Weight Prediction, {model_name}')

plt.show()


# which is doing the best? I can't tell what the results are supposed to look
# like anyway

# ### Prepare Data for BAE/BNE

# In[113]:


# Construct data from BMA samples, shapes (num_samples * num_data, ...)
means_train_mcmc, X_train_mcmc, Y_train_mcmc = make_bma_samples(
    X_train, Y_train, base_preds_train, 
    bma_weight_samples=bma_gp_w_samples[0],
    bma_model_config=bma_model_config,
    n_samples=bma_n_samples_train,
    seed=bma_seed, 
    prepare_mcmc_training=True)

# Mean samples based on test data, shape (num_samples, num_data, num_output).
# It is used to generate final examples in `make_bne_samples()`.
means_test_mcmc = make_bma_samples(
    X_test, Y_test, base_preds_test, 
    bma_weight_samples=bma_gp_w_samples[0],
    bma_model_config=bma_model_config,
    n_samples=bma_n_samples_test,
    seed=bma_seed)


# ## Bayesian Additive Ensemble

# Given $\mu(x)$ the posterior of a Bayesian ensemble model, the Bayesian Additive Ensemble is defined as:    
# 
# $y \sim N(\mu(x) + r(x), \sigma^2)$
# 
# $r \sim GaussianProcess(0, k)$
# 
# The additive ensemble $r(x)$ services two purposes: 
# 
# 1. Mitigates systematic bias in model prediction; 
# 2. Quantifies the model's epistemic uncertainty.

# In[114]:


# BNE GP Configs.
lengthscale = 1. # @param
l2_regularizer = 10. # @param

# BNE model configs. 
# If estimate_mean=False, only estimates a constant variance on top of the 
# original model.
estimate_mean = "True" # @param ["True", "False"]
variance_prior_mean=0. # @param


# In[115]:


# MAP and MCMC configs
map_step_size=0.1 # @param
map_num_steps=10_000 # @param

mcmc_step_size=1e-2 # @param
mcmc_num_steps=10_000 # @param


# In[116]:


bne_gp_config = DEFAULT_GP_CONFIG.copy()
bne_model_config = DEFAULT_BNE_CONFIG.copy()

map_config = DEFAULT_MAP_CONFIG.copy()
mcmc_config = DEFAULT_MCMC_CONFIG.copy()


bne_gp_config.update(dict(lengthscale=lengthscale, 
                          l2_regularizer=l2_regularizer))
bne_model_config.update(dict(estimate_mean=eval(estimate_mean),
                             variance_prior_mean=variance_prior_mean,
                             **bne_gp_config))

map_config.update(dict(learning_rate=map_step_size,
                       num_steps=map_num_steps))
mcmc_config.update(dict(step_size=mcmc_step_size, 
                        num_steps=mcmc_num_steps))


# ### Define Model

# In[117]:


# Construct posterior sampler.
bne_prior, bne_gp_config = bne_model_dist(
    inputs=X_train_mcmc,
    mean_preds=means_train_mcmc,
    **bne_model_config)

bne_model_config.update(bne_gp_config)
print(f'prior model graph: {bne_prior.resolve_graph()}')


# ### Run MCMC

# In[118]:


# Estimates GP weight posterior using MCMC.
bne_gp_w_samples = run_posterior_inference(model_dist=bne_prior,
                                           model_config=bne_gp_config,
                                           Y=Y_train_mcmc,
                                           map_config=map_config,
                                           mcmc_config=mcmc_config,
                                           initialize_from_map=True)


# In[ ]:


# Generates the posterior sample for all model parameters. 
bne_joint_samples = make_bne_samples(X_test,
                                     mean_preds=means_test_mcmc,
                                     bne_model_config=bne_model_config,
                                     bne_weight_samples=bne_gp_w_samples[0],
                                     seed=bne_seed)


# ### Inspect Posterior

# In[ ]:


# Plot raw samples.
plt.figure(figsize=(12, 6))
plot_base_models(base_preds_test, kernel_names,
                 X_test, Y_test, X_train, Y_train, X_base, Y_base)

for sample in bne_joint_samples['y']:
  plt.plot(X_test, sample, alpha=0.01, c='k', linewidth=2)
  plt.ylim([-1, 15])

plt.title('Posterior Predictive, Bayesian Additive Ensemble')
plt.show()


# ## Bayesian Nonparametric Ensemble (Variance Only)

# So far, we are only estimating the mean-component of the model, i.e., we are assuming: 
# 
# $y \sim Gaussian(m(x), \sigma^2); \quad m(x) = GP(0, k)$.
# 
# By doing so, the model is implicitly assuming the distribution of $y$ is always a symmetric Gaussian distribution with constant mean across space and time. As a result, our model can only quantify model uncertainty (due to lack of data) via the GP prior, but cannot flexibly capture the data uncertainty that is inherent to the empirical distribution of y.
# 
# To resolve this, we extend the ensemble's outcome distribution $y | f$ by also estimating the higher moments of the data distribution (e.g., variance, skewness, etc) using flexible estimators. Specifically, we specify the outcome distribution family to the [maximum-entropy distribution](https://en.wikipedia.org/wiki/Principle_of_maximum_entropy) given the known moments, so the predictive distribution is [minimax](https://arxiv.org/pdf/math/0410076.pdf) and still statistically efficient to estimate.

# For example, when we want to estimate the first two moments (mean and variance) of the distribution, this leads to a Gaussian distribution with spatio-temporally adaptive variance $\sigma(x)^2$:
# 
# $$y \sim Gaussian(m(x), \sigma(x)^2); \quad \mbox{where} \quad m \sim GP(0, k_m), \sigma \sim GP(0, k_\sigma)$$
# 
# and when we want to estimate the first three moments (mean and variance) of the distribution, this leads to a [Exponentially-modifed Gaussian](https://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution) (EMG) distribution with spatio-temporally adaptive variance $\sigma(x)^2$ and skewness $\lambda(x)$:
# 
# $$y \sim EMG(m(x), \sigma(x)^2, \lambda(x)); \quad \mbox{where} \quad m \sim GP(0, k_m), \sigma \sim GP(0, k_\sigma), \lambda \sim GP(0, k_\lambda)$$
# 
# 

# In[ ]:


# BNE GP Configs.
lengthscale = 1. # @param
l2_regularizer = 10. # @param

# BNE model configs.
variance_prior_mean=0. # @param


# In[ ]:


# MAP and MCMC configs

map_step_size=5e-3 # @param
map_num_steps=10_000 # @param

mcmc_step_size=1e-2 # @param
mcmc_num_steps=10_000 # @param


# In[ ]:


bne_gp_config = DEFAULT_GP_CONFIG.copy()
bne_model_config = DEFAULT_BNE_CONFIG.copy()

map_config = DEFAULT_MAP_CONFIG.copy()
mcmc_config = DEFAULT_MCMC_CONFIG.copy()


bne_gp_config.update(dict(lengthscale=lengthscale, 
                          l2_regularizer=l2_regularizer))
bne_model_config.update(dict(estimate_variance=True,
                             variance_prior_mean=variance_prior_mean,
                             **bne_gp_config))

map_config.update(dict(learning_rate=map_step_size,
                       num_steps=map_num_steps))
mcmc_config.update(dict(step_size=mcmc_step_size, 
                        num_steps=mcmc_num_steps))


# ### Define Model

# In[ ]:


# Construct posterior sampler.
bne_prior, bne_gp_config = bne_model_dist(
    inputs=X_train_mcmc,
    mean_preds=means_train_mcmc,
    **bne_model_config)

bne_model_config.update(bne_gp_config)
print(f'prior model graph: {bne_prior.resolve_graph()}')


# ### Run MCMC

# In[ ]:


# Estimates GP weight posterior using MCMC.
bne_gp_w_samples = run_posterior_inference(model_dist=bne_prior,
                                           model_config=bne_gp_config,
                                           Y=Y_train_mcmc,
                                           map_config=map_config,
                                           mcmc_config=mcmc_config,
                                           initialize_from_map=True)


# In[ ]:


# Generates the posterior sample for all model parameters. 
bne_joint_samples = make_bne_samples(X_test,
                                     mean_preds=means_test_mcmc,
                                     bne_model_config=bne_model_config,
                                     bne_weight_samples=bne_gp_w_samples[0],
                                     seed=bne_seed)


# ### Inspect Posterior

# In[ ]:


# Plot raw samples.
plt.figure(figsize=(12, 6))
plot_base_models(base_preds_test, kernel_names,
                 X_test, Y_test, X_train, Y_train, X_base, Y_base)

for sample in bne_joint_samples['y']:
  plt.plot(X_test, sample, alpha=0.01, c='k', linewidth=2)
  plt.ylim(-1, 15)

plt.title('Posterior Predictive, Bayesian Nonparametric Ensemble')
plt.show()


# ## Bayesian Nonparametric Ensemble (Variance + Skewness)

# In[ ]:


# BNE GP Configs.
lengthscale = 1. # @param
l2_regularizer = 10. # @param

# BNE model configs.
variance_prior_mean=0. # @param
skewness_prior_mean=0. # @param


# In[ ]:


# MAP and MCMC configs
map_step_size=5e-3 # @param
map_num_steps=10_000 # @param

mcmc_step_size=1e-2 # @param
mcmc_num_steps=10_000 # @param


# In[ ]:


bne_gp_config = DEFAULT_GP_CONFIG.copy()
bne_model_config = DEFAULT_BNE_CONFIG.copy()

map_config = DEFAULT_MAP_CONFIG.copy()
mcmc_config = DEFAULT_MCMC_CONFIG.copy()


bne_gp_config.update(dict(lengthscale=lengthscale, 
                          l2_regularizer=l2_regularizer))
bne_model_config.update(dict(estimate_variance=True,
                             estimate_skewness=True,
                             variance_prior_mean=variance_prior_mean,
                             skewness_prior_mean=skewness_prior_mean,
                             **bne_gp_config))

map_config.update(dict(learning_rate=map_step_size,
                       num_steps=map_num_steps))
mcmc_config.update(dict(step_size=mcmc_step_size, 
                        num_steps=mcmc_num_steps))


# ### Define Model

# In[ ]:


# Construct prior distribution.
bne_prior, bne_gp_config = bne_model_dist(
    inputs=X_train_mcmc,
    mean_preds=means_train_mcmc,
    **bne_model_config)

bne_model_config.update(bne_gp_config)
print(f'prior model graph: {bne_prior.resolve_graph()}')


# ### Run MCMC

# In[ ]:


# Estimates GP weight posterior using MCMC.
bne_gp_w_samples = run_posterior_inference(model_dist=bne_prior,
                                           model_config=bne_gp_config,
                                           Y=Y_train_mcmc,
                                           map_config=map_config,
                                           mcmc_config=mcmc_config,
                                           initialize_from_map=True)


# In[ ]:


# Generates the posterior sample for all model parameters. 
bne_joint_samples = make_bne_samples(X_test,
                                     mean_preds=means_test_mcmc,
                                     bne_model_config=bne_model_config,
                                     bne_weight_samples=bne_gp_w_samples[0],
                                     seed=bne_seed)


# ### Inspect Posterior

# In[ ]:


# Plot raw samples.
plt.figure(figsize=(12, 6))
plot_base_models(base_preds_test, kernel_names,
                 X_test, Y_test, X_train, Y_train, X_base, Y_base)

for sample in bne_joint_samples['y']:
  plt.plot(X_test, sample, alpha=0.01, c='k', linewidth=2)
  plt.ylim(-1, 15)

plt.title('Posterior Predictive, Bayesian Nonparametric Ensemble')
plt.show()

