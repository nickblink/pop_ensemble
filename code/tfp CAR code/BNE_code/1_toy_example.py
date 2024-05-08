#!/usr/bin/env python
# coding: utf-8

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

# In[22]:



# # Default Configs

# In[2]:


# GP configs.
y_noise_std = 0.1  # @param
hidden_units = 128  # @param
lengthscale=1.  # @param
l2_regularizer=0.1  # @param

DEFAULT_GP_CONFIG = dict(lengthscale=lengthscale,
                         l2_regularizer=l2_regularizer, 
                         hidden_units=hidden_units, 
                         y_noise_std=y_noise_std)


# In[3]:


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


# In[9]:


# MAP configs.
map_step_size=0.1 # @param
map_num_steps=10_000 # @param

DEFAULT_MAP_CONFIG = dict(learning_rate=map_step_size,
                          num_steps=map_num_steps)


# In[8]:


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


# # Toy Example: Gaussian Process Regression

# In[7]:


n_train = 10  # @param
n_test = 1000  # @param

feature_dim=1  # @param
hidden_dim=512  # @param

data_seed=0  # @param
model_seed=0  # @param
mcmc_seed=0  # @param


# ## Data

# In[11]:


np.random.seed(data_seed)
X = np.sort(np.random.uniform(-np.pi, np.pi, size=(n_train, feature_dim))).astype(np.float32)
X_test = np.linspace(-2*np.pi, 2*np.pi, n_test).reshape((n_test, feature_dim)).astype(np.float32)

Y = np.cos(X)
Y_test = np.cos(X_test)


# In[14]:


X


# In[17]:


plt.plot(X_test, Y_test, c='k')
plt.scatter(X, Y, c='r')

plt.ylim([-5., 5.])


# ## Build Model

# In[18]:

# rfgp in bne_0_2_core_model_functions.
# creates and RBF Gaussian Process


gp_prior, gp_config = rfgp_dist(
    X, lengthscale=1.0, y_noise_std=0.1, l2_regularizer=0.1)


# In[55]:


gp_config


# In[56]:


gp_prior.log_prob(gp_prior.sample())


# ## Run MCMC

# In[57]:


samples, _ = run_mcmc(model_dist=gp_prior, y=Y,
                      step_size=.1, num_steps=500, 
                      burnin=100, nchain=4, seed=0)


# ## Examine Posterior.

# In[58]:


# Construct posterior.
posterior_sample = samples[0]  # Shape [num_samples, param_dim]

gp_dist, _ = rfgp_dist(
    X_test, 
    posterior_mode=True, 
    posterior_sample=posterior_sample,
    **gp_config)


# In[59]:


# Visualize posterior prediction.
y_samples = gp_dist.sample(1000)["y"]

[plt.plot(X_test, sample, alpha=0.01, c='k') for sample in y_samples.numpy()]

plt.plot(X_test, Y_test, c='k')
plt.scatter(X, Y, c='r')

plt.ylim([-5., 5.])
plt.show()

# So this shows how outside of the training window, the predictions go crazy. 
# I still don't get what sort of model fitting is really done here.
# We have a prior on the spread of Y from a GP, but that's not gonna do well
# for points far from the observed ones anyway. And then do we run some kind
# of feature regression or is it just a GP? Not sure this is super important to understand.