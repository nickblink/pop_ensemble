#!/usr/bin/env python
# coding: utf-8

# # 0_c Core Model Functions

# ## 0_c_1 Exponentially Modified Gaussian
# * ExponentiallyModifiedGaussian

# In[25]:


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


# ## 1_c_2 Model: Random Feature
# * make_random_feature

# In[26]:


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


# ## 1_c_3 Model: Gaussian Process
# * rfgp_dist

# In[27]:


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


# ## 1_c_4 Model: Bayesian Model Averaging
# * bma_dist
# * build_ensemble_weight_logits
# * build_ensemble_weights
# * ensemble_prediction
# * ensemble_mean

# In[28]:


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


# ## 1_c_5 Model: Bayesian Nonparametric Ensemble
# * bne_model_dist

# In[29]:


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


# ## 1_c_6 MCMC: Low-level Sampler
# * run_chain
# * mix_chain_samples

# In[30]:


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


# ## 1_c_7 MCMC: Main API
# * run_mcmc

# In[31]:


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


# ## 1_c_8 MAP: Main API
# * run_map

# In[32]:


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

