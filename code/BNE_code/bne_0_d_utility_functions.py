#!/usr/bin/env python
# coding: utf-8

# # 0_d Utility Functions

# ## 0_d_1 Visualization: Plot Description
# * plot_description

# In[33]:


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
      


# ## 0_d_2 Visualization: Plot Base Models
# * plot_base_models

# In[34]:


def plot_base_models(base_preds_test, kernel_names,
                     X_test, Y_test, X_train, Y_train, X_base, Y_base, 
                     ax=plt):
  for base_pred in tf.transpose(base_preds_test):
    ax.plot(X_test, base_pred)

  ax.scatter(X_train, Y_train, c='b', s=10)
  ax.scatter(X_base, Y_base, c='r', s=10)
  ax.scatter(X_test, Y_test, c='k', alpha=0.1)

  ax.legend(kernel_names)


# ## 0_d_3 Data: Get Base Model Predictions
# * get_base_predictions

# In[35]:


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


# ## 0_d_4 Data: Process MCMC Data
# * process_mcmc_data
# * flatten_mcmc_data
# * deflatten_mcmc_sample

# In[36]:


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


# ## 0_d_5 Model: Prepare MCMC
# * prepare_mcmc

# In[37]:


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


# ## 0_d_6 Model: Make BMA Samples
# * make_bma_samples

# In[38]:


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


# ## 0_d_7 Model: Make BNESamples
# * make_bne_samples

# In[39]:


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

