#!/usr/bin/env python
# coding: utf-8

# # 0_e Simulation Functions

# ## 0_e_1 Simulation: Get Data 
# * get_data

# In[40]:


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


# ## 0_e_2 Simulation: Get BMA Result 
# * get_bma_result

# In[41]:


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


# ## 0_e_3 Simulation: Get BNE Result 
# * get_bne_result

# In[42]:


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


# ## 0_e_4 Simulation: Save/Load from Drive
# * save_to_drive 
# * load_from_drive

# In[43]:


def save_to_drive(data, file_name):
  path_name = os.path.join(FULL_DATA_PATH, f'{file_name}.pkl')
  with tf.io.gfile.GFile(path_name, 'wb') as f:
    pickle.dump(data, f)

def load_from_drive(file_name):  
  path_name = os.path.join(FULL_DATA_PATH, f'{file_name}.pkl')
  with tf.io.gfile.GFile(path_name, 'rb') as f:
    data = pickle.load(f)
  return data


# ## 0_e_5 Simulation: Compute Metrics
# * compute_metrics

# In[44]:


def compute_metrics(data_dict, model_name, q_true=None, ind_ids=None, num_sample=None):
  if q_true is None:
    q_true = np.array(
        [0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25,
         0.75, 0.8, 0.85, 0.9, 0.925, 0.95, 0.975])

  if ind_ids is None:
    # Find IDs of in-domain test data via range comparison 
    # between X_train and X_test.
    X_train_min = np.min(data_dict['X_train'], axis=0)
    X_train_max = np.max(data_dict['X_train'], axis=0)

    test_ids_greater_than_min = np.all(
        data_dict['X_test'] > X_train_min, axis=-1)
    test_ids_less_than_max = np.all(
        data_dict['X_test'] < X_train_max, axis=-1)

    ind_ids = np.where(
        np.logical_and(test_ids_greater_than_min, test_ids_less_than_max))[0]

  samples = data_dict[f'{model_name}_samples']
  means_true = data_dict['mean_test']
  y_test = data_dict['Y_test']

  if num_sample is not None:
    samples = samples[:num_sample]

  means_pred = np.mean(samples, axis=0)
  stds_pred = np.std(samples, axis=0)
  quantile_pred = np.quantile(samples, q=q_true, axis=0)

  # Compute in-domain metrics.
  nll_ind = np.mean(
      ((means_pred[ind_ids] - means_true[ind_ids])/stds_pred[ind_ids])**2)
  mse_ind = np.mean((means_pred[ind_ids] - means_true[ind_ids])**2)  

  q_pred_ind = np.mean(y_test[ind_ids] < quantile_pred[:, ind_ids], axis=(1, 2))
  ece_ind = np.mean((q_pred_ind - q_true)**2)
  cov_prob_95_ind = q_pred_ind[-1] - q_pred_ind[0]
  cov_prob_90_ind = q_pred_ind[-2] - q_pred_ind[1]
  cov_prob_85_ind = q_pred_ind[-3] - q_pred_ind[2]
  cov_prob_80_ind = q_pred_ind[-4] - q_pred_ind[3]

  # Compute all-domain (ind + ood) metrics.
  nll_all = np.mean(((means_pred - means_true)/stds_pred)**2)
  mse_all = np.mean((means_pred - means_true)**2)

  q_pred_all = np.mean(y_test < quantile_pred, axis=(1, 2))
  ece_all = np.mean((q_pred_all - q_true)**2)
  cov_prob_95_all = q_pred_all[-1] - q_pred_all[0]
  cov_prob_90_all = q_pred_all[-2] - q_pred_all[1]
  cov_prob_85_all = q_pred_all[-3] - q_pred_all[2]
  cov_prob_80_all = q_pred_all[-4] - q_pred_all[3]

  return (nll_ind, mse_ind, ece_ind, cov_prob_95_ind, cov_prob_90_ind, cov_prob_85_ind, cov_prob_80_ind,
          nll_all, mse_all, ece_all, cov_prob_95_all, cov_prob_90_all, cov_prob_85_all, cov_prob_80_all)


# ## 0_e_6 Simulation: Run Pipeline
# * run_pipeline

# In[45]:


def run_pipeline(seeds, N_train, group_id, data_gen_fn, base_train_steps=200):
  # Data Generation
  data_dicts = {}
  print('Data:', end='', flush=True)
  t0 = time.time()
  for seed in seeds:
    print(f'Run {seed+1}...', end='', flush=True)
    data_dicts[seed] = get_data(seed, N_train, y_dist, data_gen_fn, 
                                num_train_steps=base_train_steps)
  print(f'Time: {(time.time()-t0)/60.:.4f} min.')

  # BMA-mean.
  print('BMA-mean:', flush=True)
  t0 = time.time()
  for seed in seeds:
    print(f'Run {seed+1}: ', end='', flush=True)
    data_dicts[seed] = get_bma_result(data_dicts[seed], bma_config=bma_config) 
  print(f'Time: {(time.time()-t0)/60.:.4f} min.', flush=True)
  tf.keras.backend.clear_session()
  gc.collect()

  # BMA.
  print('BMA:', flush=True)
  # Inhere BMA MCMC configs.
  bma_var_config = bne_config.copy()
  bma_var_config['mcmc_initialize_from_map'] = bma_config['mcmc_initialize_from_map']
  t0 = time.time()
  for seed in seeds:
    print(f'Run {seed+1}: ', end='', flush=True)
    data_dicts[seed] = get_bne_result(data_dicts[seed], moment_mode='none', 
                                      bne_config=bma_var_config) 
  print(f'Time: {(time.time()-t0)/60.:.4f} min.', flush=True)
  tf.keras.backend.clear_session()
  gc.collect()

  # BAE.
  print('BAE:', flush=True)
  t0 = time.time()
  for seed in seeds:
    print(f'Run {seed+1}: ', end='', flush=True)
    data_dicts[seed] = get_bne_result(data_dicts[seed], moment_mode='mean', 
                                      bne_config=bne_config) 
  print(f'Time: {(time.time()-t0)/60.:.4f} min.', flush=True)
  tf.keras.backend.clear_session()
  gc.collect()

  # BNE-Variance.
  print('BNE-Variance:', flush=True)
  t0 = time.time()
  for seed in seeds:
    print(f'Run {seed+1}: ', end='', flush=True)
    data_dicts[seed] = get_bne_result(data_dicts[seed], moment_mode='variance', 
                                      bne_config=bne_config) 
  print(f'Time: {(time.time()-t0)/60.:.4f} min.', flush=True)
  tf.keras.backend.clear_session()
  gc.collect()

  # BNE-Skewness.
  print('BNE-Skewness:', flush=True)
  t0 = time.time()
  for seed in seeds:
    print(f'Run {seed+1}: ', end='', flush=True)
    data_dicts[seed] = get_bne_result(data_dicts[seed], moment_mode='skewness', 
                                      bne_config=bne_config) 
  print(f'Time: {(time.time()-t0)/60.:.4f} min.', flush=True)
  tf.keras.backend.clear_session()
  gc.collect()

  return data_dicts

