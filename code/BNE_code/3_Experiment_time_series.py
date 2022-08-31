#!/usr/bin/env python
# coding: utf-8

# # Default Configs

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


# # Experiment I: 1D Time Series

# ## Data Configs

# Use `DATA_PATH` to specify the Google Drive folder for storing the results. (make sure to create it in Google Drive beforehand). 

# In[ ]:


DRIVE_PATH = '/content/drive/'   # @param ['/content/drive/']
DATA_PATH = '/My Drive/BNE_simu'   # @param ['/My Drive/BNE_simu']

FULL_DATA_PATH = f'{DRIVE_PATH}{DATA_PATH[1:]}'
print(f'FULL_DATA_PATH: {FULL_DATA_PATH}')


# In[ ]:


N_train_grid = (10, 25, 50, 100, 250, 500, 750, 1000)   # @param
N_base = 50 # @param  
N_test = 1000  # @param
y_dist_name = "lognormal" # @param ["lognormal", "normal", "cauchy"]

num_runs = 100  # @param
num_run_groups = 20   # @param

seed_groups = [range(num_runs)[i:i+num_runs//num_run_groups] for 
               i in range(0, num_runs, num_runs//num_run_groups)]


dist_dict = {'lognormal': np.random.lognormal,
             'normal': np.random.normal,
             'cauchy': lambda m, s: np.random.standard_cauchy(size=(m.shape[0], 1)) * s + m}
y_dist = dist_dict[y_dist_name]


# ## Model Configs

# In[ ]:


# Optimization configs. 
# Consider reduce below parameters / set to `False` if MCMC is taking too long:
# mcmc_num_steps, mcmc_burnin, mcmc_nchain, mcmc_initialize_from_map.
map_step_size=5e-3 # @param
map_num_steps=10_000 # @param

mcmc_step_size=1e-3 # @param
mcmc_num_steps=500 # @param

mcmc_nchain=1 # @param
mcmc_burnin=50 # @param
mcmc_initialize_from_map=False # @param


# In[ ]:


# BMA parameters.
y_noise_std = 0.1 # @param
bma_gp_lengthscale = 1. # @param
bma_gp_l2_regularizer = 0.1 # @param

bma_n_samples_eval = 100 # @param
bma_n_samples_train = 50 # @param
bma_n_samples_test = 100 # @param
bma_seed = 0 # @param


# In[ ]:


# BNE parameters.
bne_gp_lengthscale = 1. # @param
bne_gp_l2_regularizer = 10. # @param
bne_variance_prior_mean = 0. # @param
bne_skewness_prior_mean = 0. # @param
bne_seed = 0 # @param


# In[ ]:


bma_config=dict(gp_lengthscale=bma_gp_lengthscale,
                gp_l2_regularizer=bma_gp_l2_regularizer,
                y_noise_std=y_noise_std,
                map_step_size=map_step_size,
                map_num_steps=map_num_steps,
                mcmc_step_size=mcmc_step_size,
                mcmc_num_steps=mcmc_num_steps,
                mcmc_initialize_from_map=mcmc_initialize_from_map,
                n_samples_eval=bma_n_samples_eval, 
                n_samples_train=bma_n_samples_train, 
                n_samples_test=bma_n_samples_test, 
                seed=bma_seed)

bne_config = dict(gp_lengthscale=bne_gp_lengthscale,
                  gp_l2_regularizer=bne_gp_l2_regularizer,
                  variance_prior_mean=bne_variance_prior_mean,
                  skewness_prior_mean=bne_skewness_prior_mean,
                  map_step_size=map_step_size,
                  map_num_steps=map_num_steps,
                  mcmc_step_size=mcmc_step_size,
                  mcmc_num_steps=mcmc_num_steps,
                  mcmc_nchain=mcmc_nchain,
                  mcmc_burnin=mcmc_burnin,
                  mcmc_initialize_from_map=mcmc_initialize_from_map,
                  seed=bne_seed)


# ## Run Pipeline

# In[ ]:


# Configure saving.
drive.mount(DRIVE_PATH, force_remount=True)

tf.io.gfile.makedirs(FULL_DATA_PATH)


# In[ ]:


for group_id, seeds in enumerate(seed_groups):
  # if group_id > 1:
  #   continue

  for N_train in N_train_grid:
    print(f'=======Run simulation with n={N_train}, group={group_id}=======')
    # Run pipeline.
    data_dicts = run_pipeline(seeds, N_train, group_id, 
                              data_gen_fn=generate_data_1d)
    save_to_drive(data_dicts, file_name=f'result_{N_train}_{group_id}')
    del data_dicts


# ## Evaluation

# 
# * BAE/BNE models always has better RMSE and Calibration Error.
# * BNE-skewness has best coverage probability for 95\% in small sample, and remains competitive has sample size grow.
# * BMA has decent coverage probability for $95\%$ only when averaged across the whole dataset but not locally. Therefore may lead to under- or over-estimation of uncertainty depending on the location in the feature space.

# In[ ]:


# Compute metrics for all N_train and models. 
metric_rows = []

for N_train in N_train_grid:
  data_dicts = {}
  for group_id in range(len(seed_groups)):
    data_dict_group = load_from_drive(f'result_{N_train}_{group_id}')
    data_dicts.update(data_dict_group)

  for model_name in ('bma', 'bae', 'bne_var', 'bne_skew'):
    metrics = [compute_metrics(data, model_name) for data in data_dicts.values()]
    metrics = np.stack(metrics)
    metric_means = np.mean(metrics, axis=0)
    metric_stds = np.std(metrics, axis=0)

    metric_row = dict(n=N_train, 
                      model_name=model_name, 
                      # Metric means.
                      nll_ind=metric_means[0], 
                      mse_ind=metric_means[1], 
                      ece_ind=metric_means[2], 
                      cov_prob_95_ind=metric_means[3],
                      cov_prob_90_ind=metric_means[4],
                      cov_prob_85_ind=metric_means[5],
                      cov_prob_80_ind=metric_means[6],
                      nll_all=metric_means[7], 
                      mse_all=metric_means[8], 
                      ece_all=metric_means[9], 
                      cov_prob_95_all=metric_means[10],
                      cov_prob_90_all=metric_means[11],
                      cov_prob_85_all=metric_means[12],
                      cov_prob_80_all=metric_means[13],
                      # Metric STDs.
                      nll_ind_std=metric_stds[0], 
                      mse_ind_std=metric_stds[1], 
                      ece_ind_std=metric_stds[2], 
                      cov_prob_95_ind_std=metric_stds[3],
                      cov_prob_90_ind_std=metric_stds[4],
                      cov_prob_85_ind_std=metric_stds[5],
                      cov_prob_80_ind_std=metric_stds[6],
                      nll_all_std=metric_stds[7], 
                      mse_all_std=metric_stds[8], 
                      ece_all_std=metric_stds[9], 
                      cov_prob_95_all_std=metric_stds[10],
                      cov_prob_90_all_std=metric_stds[11],
                      cov_prob_85_all_std=metric_stds[12],
                      cov_prob_80_all_std=metric_stds[13])
    
    metric_rows.append(metric_row)


# In[ ]:


metric_pd = pd.DataFrame(metric_rows)
metric_pd


# ## Visualization

# In[ ]:


# @title plot_1d_result
def plot_1d_result(data_dict, model_name='bma'):
  pred_samples = data_dict[f'{model_name}_samples']

  base_preds_test = data_dict['base_preds_test']
  base_model_names = data_dict['base_model_names']
  X_train = data_dict['X_train']
  X_test = data_dict['X_test']
  X_base = data_dict['X_base']
  
  Y_train = data_dict['Y_train']
  Y_test = data_dict['Y_test']
  Y_base = data_dict['Y_base']

  # Plot base models and data.
  plt.figure(figsize=(12, 6))
  plot_base_models(base_preds_test, base_model_names,
                   X_test, Y_test, X_train, Y_train, X_base, Y_base)

  # Plot predictive samples.
  for sample in pred_samples:
    plt.plot(X_test, sample, alpha=0.01, c='k', linewidth=2)
    plt.ylim(-1, 15)

  plt.title(f'Posterior Predictive, {model_name.upper()}')
  plt.show()


# In[ ]:


N_train = "1000"  # @param [10, 25, 50, 100, 250, 500, 750, 1000]
group_id = 0  # param


# In[ ]:


# Load data from disk and plot it's first run.
data_dicts = load_from_drive(f'result_{N_train}_{group_id}')
data_dict = data_dicts[0]


# In[ ]:


plot_result(data_dicts[0], model_name='bma')


# In[ ]:


plot_result(data_dicts[0], model_name='bae')


# In[ ]:


plot_result(data_dicts[0], model_name='bne_var')


# In[ ]:


plot_result(data_dicts[0], model_name='bne_skew')


# # Experiment II: 2D Spatial Field

# In[ ]:





# In[ ]:




