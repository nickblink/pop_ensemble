#!/usr/bin/env python
# coding: utf-8

# # 0_a Load Modules 

# ## 0_a_1 Load Modules

# In[7]:


# import modules not pre-installed in colab
get_ipython().system(' pip install edward2 ')
get_ipython().system(' pip install gpflow')


# In[8]:


# install required packages that are pre-installd in colab 
get_ipython().system(' pip install numpy')
get_ipython().system(' pip install tensorflow ')
get_ipython().system(' pip install pandas')
get_ipython().system(' pip install matplotlib.pyplot')
get_ipython().system(' pip install seaborn')


# In[9]:


from typing import Any, Callable, Dict, List, Optional, Union, Tuple

import os
import gc
import time
import pickle
import functools

import multiprocessing as mp

#from google.colab import files
#from google.colab import drive


# In[10]:


import numpy as np
import tensorflow as tf

import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns


# In[11]:


import edward2 as ed
import tensorflow_probability as tfp

tfd = tfp.distributions
tfb = tfp.bijectors

dtype = tf.float32


# In[12]:


import gpflow as gpf


# In[13]:


import logging

logging.getLogger('tensorflow').setLevel(logging.ERROR)  # suppress pfor warnings


# ### import modules for exponentially-modified Gaussian

# In[ ]:


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


# In[14]:


# Verify versions.
print(f'TensorFlow version: {tf.__version__}. Expected: 2.7.0')
print(f'TensorFlow Probability version: {tfp.__version__}. Expected: 0.15.0')


# If possible (e.g., using Window, Unix, or Google Colab), enable GPU support. For Mac computers, gpu is not available and the user should anticipate longer run times. 

# In[16]:


tf.test.gpu_device_name()

