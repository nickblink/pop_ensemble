"""GPR Mixutre using MCMC, mean-field VI, and augmented VI.


#### References

[1]:    Michalis Titsias. Variational Learning of Inducing Variables in
        Sparse Gaussian Processes. _12th Proceedings of AISTATS, PMLR 5:567-574_, 2009.
        http://proceedings.mlr.press/v5/titsias09a.html
[2]:    Ching-An Cheng and Byron Boots. Variational Inference for Gaussian
        Process Models with Linear Complexity. _Advances in NIPS 30_, 2017.
        http://papers.nips.cc/paper/7103-variational-inference-for-gaussian-process-models-with-linear-complexity
[3]:    Marton Havasi, José Miguel Hernández-Lobato and Juan José Murillo-Fuentes.
        Deep Gaussian Processes with Decoupled Inducing Inputs.
        _arXiv preprint arXiv:1801.02939_, 2018.
        https://arxiv.org/pdf/1801.02939.pdf
"""

import os
import time

from importlib import reload

import numpy as np

import tensorflow as tf
import tensorflow_probability as tfp
from tensorflow_probability import edward2 as ed
import gpflowSlim as gpf

# sys.path.extend([os.getcwd()])

from calibre.model import gaussian_process as gp
from calibre.model import gp_regression

import calibre.util.data as data_util
import calibre.util.visual as visual_util
import calibre.util.gp_flow as gpflow_util

from calibre.util.data import simple_sin_curve_1d, simple_cos_curve_1d, generate_1d_data
from calibre.util.inference import make_value_setter

import matplotlib.pyplot as plt

tfd = tfp.distributions

_SAVE_ADDR_PREFIX = "./result/gpr/"

"""""""""""""""""""""""""""""""""
# 1. Generate data
"""""""""""""""""""""""""""""""""

N_train = 300
N_test = 300
N_valid = 500

X_train, y_train = data_util.generate_1d_data_multimodal(
    N=N_train, f_list=[simple_sin_curve_1d, simple_cos_curve_1d],
    noise_sd=0.01, seed=1000,
    uniform_x=True, uniform_x_range=[-0.5, 6.5])
X_test, y_test = data_util.generate_1d_data_multimodal(
    N=N_test, f_list=[simple_sin_curve_1d, simple_cos_curve_1d],
    noise_sd=0.01, seed=2000,
    uniform_x=True, uniform_x_range=[-0.5, 6.5])

X_train = np.expand_dims(X_train, 1).astype(np.float32)
y_train = y_train.astype(np.float32)
X_test = np.expand_dims(X_test, 1).astype(np.float32)
y_test = y_test.astype(np.float32)

std_y_train = np.std(y_train)

X_valid = np.concatenate((np.linspace(-0.5, 6.5, N_valid),
                          X_test.squeeze()))
X_valid = np.expand_dims(np.sort(X_valid), 1).astype(np.float32)
y_valid = np.concatenate([simple_sin_curve_1d(X_valid),
                          simple_cos_curve_1d(X_valid)])
X_valid = np.concatenate([X_valid, X_valid])

N, D = X_train.shape

#
plt.plot(np.linspace(-0.5, 6.5, 100),
         simple_sin_curve_1d(np.linspace(-0.5, 6.5, 100)), c='black')
plt.plot(np.linspace(-0.5, 6.5, 100),
         simple_cos_curve_1d(np.linspace(-0.5, 6.5, 100)), c='black')
plt.plot(X_train.squeeze(), y_train.squeeze(),
         'o', c='red', markeredgecolor='black')
plt.plot(X_test.squeeze(), y_test.squeeze(),
         'o', c='blue', markeredgecolor='black')
plt.close()

"""""""""""""""""""""""""""""""""
# 2. MCMC, simple mixture
"""""""""""""""""""""""""""""""""
"""2.1. sampler basic config"""
num_results = 10000
num_burnin_steps = 10000
ls_val = 1.
n_mix_val = 20

# define mcmc computation graph
mcmc_graph = tf.Graph()
with mcmc_graph.as_default():
    # build likelihood by explicitly
    log_joint = ed.make_log_joint_fn(gp_regression.model_mixture)


    def target_log_prob_fn(mix_prob, gp_f, sigma):
        """Unnormalized target density as a function of states."""
        return log_joint(X_train, y=y_train, ls=ls_val, n_mix=n_mix_val,
                         mix_prob=mix_prob, gp_f=gp_f, sigma=sigma)


    # set up state container
    initial_state = [
        tf.constant([1 / n_mix_val] * n_mix_val, name='init_mix_prob'),
        tf.random_normal([n_mix_val, N], stddev=0.01, name='init_gp_func'),
        tf.constant([0.1] * n_mix_val, name='init_sigma'),
    ]

    # set up HMC transition kernel
    step_size = tf.get_variable(
        name='step_size',
        initializer=1.,
        use_resource=True,  # For TFE compatibility.
        trainable=False)

    unconstraining_bijectors = [
        tfp.bijectors.SoftmaxCentered(),
        tfp.bijectors.Identity(),
        tfp.bijectors.Identity(),
    ]

    hmc = tfp.mcmc.HamiltonianMonteCarlo(
        target_log_prob_fn=target_log_prob_fn,
        num_leapfrog_steps=3,
        step_size=step_size,
        step_size_update_fn=tfp.mcmc.make_simple_step_size_update_policy(),
    )

    mcmc_kernel = tfp.mcmc.TransformedTransitionKernel(
        inner_kernel=hmc, bijector=unconstraining_bijectors
    )

    # set up main sampler
    state, kernel_results = tfp.mcmc.sample_chain(
        num_results=num_results,
        num_burnin_steps=num_burnin_steps,
        current_state=initial_state,
        kernel=mcmc_kernel,
        parallel_iterations=10
    )

    mix_prob_sample, gpf_sample, sigma_sample, = state

    # set up init op
    init_op = tf.global_variables_initializer()

    mcmc_graph.finalize()

""" 2.2. execute sampling"""
with tf.Session(graph=mcmc_graph) as sess:
    init_op.run()
    [mix_prob_sample_val,
     f_samples_full_val,
     sigma_sample_val,
     is_accepted_,
     ] = sess.run(
        [
            mix_prob_sample,
            gpf_sample,
            sigma_sample,
            kernel_results.inner_results.is_accepted,
        ])
    print('Acceptance Rate: {}'.format(np.mean(is_accepted_)))
    sess.close()

""" 2.3. prediction and visualization"""
# prediction by sample from posterior
f_sample_indices = np.asarray([
    np.argmax(np.random.multinomial(1, pvals=mix_prob))
    for mix_prob in mix_prob_sample_val
])

f_samples_val = np.asarray([f_samples_full_val[i, f_sample_indices[i]]
                            for i in range(num_results)])

f_valid_val = gp.sample_posterior_full(X_new=X_valid, X=X_train,
                                       f_sample=f_samples_val.T,
                                       ls=ls_val, kernel_func=gp.rbf)

# visualize
mu = np.mean(f_valid_val, axis=1)
cov = np.var(f_valid_val, axis=1)

visual_util.gpr_1d_visual(mu, cov, pred_quantiles=[],
                          X_train=X_train, y_train=y_train,
                          X_test=X_valid, y_test=y_valid,
                          title="GP Mixture Posterior, {} Mixtures, Hamilton MC".format(n_mix_val),
                          save_addr=os.path.join(_SAVE_ADDR_PREFIX,
                                                 "gpr_hmc_mixture.png"))

visual_util.gpr_1d_visual(pred_mean=None, pred_cov=None, pred_quantiles=[],
                          pred_samples=list(f_valid_val.T)[:10000],
                          X_train=X_train, y_train=y_train,
                          X_test=X_valid, y_test=y_valid,
                          title="GP Mixture Posterior Samples, {} Mixtures, Hamilton MC".format(n_mix_val),
                          save_addr=os.path.join(_SAVE_ADDR_PREFIX,
                                                 "gpr_hmc_mixture_sample.png")
                          )

"""""""""""""""""""""""""""""""""
# 2. MCMC, adaptive mixture
"""""""""""""""""""""""""""""""""
"""2.1. sampler basic config"""
num_results = 5000
num_burnin_steps = 5000
ls_val = 1.
n_mix_val = 5

# define mcmc computation graph
mcmc_graph = tf.Graph()
with mcmc_graph.as_default():
    # build likelihood by explicitly
    log_joint = ed.make_log_joint_fn(gp_regression.model_mixture)


    def target_log_prob_fn(mix_member, gp_weight, gp_comp, sigma):
        """Unnormalized target density as a function of states."""
        return log_joint(X_train, y=y_train, ls=ls_val, n_mix=n_mix_val,
                         gp_weight=gp_weight, gp_comp=gp_comp,
                         mix_member=mix_member, sigma=sigma)


    # set up state container
    initial_state = [
        tf.constant([1 / n_mix_val] * n_mix_val, name='init_mix_prob'),
        tf.random_normal([n_mix_val, N], stddev=0.01, name='init_gp_func'),
        tf.constant([0.1] * n_mix_val, name='init_sigma'),
    ]

    # set up HMC transition kernel
    step_size = tf.get_variable(
        name='step_size',
        initializer=1.,
        use_resource=True,  # For TFE compatibility.
        trainable=False)

    unconstraining_bijectors = [
        tfp.bijectors.SoftmaxCentered(),
        tfp.bijectors.Identity(),
        tfp.bijectors.Identity(),
    ]

    hmc = tfp.mcmc.HamiltonianMonteCarlo(
        target_log_prob_fn=target_log_prob_fn,
        num_leapfrog_steps=3,
        step_size=step_size,
        step_size_update_fn=tfp.mcmc.make_simple_step_size_update_policy(),
    )

    mcmc_kernel = tfp.mcmc.TransformedTransitionKernel(
        inner_kernel=hmc, bijector=unconstraining_bijectors
    )

    # set up main sampler
    state, kernel_results = tfp.mcmc.sample_chain(
        num_results=num_results,
        num_burnin_steps=num_burnin_steps,
        current_state=initial_state,
        kernel=mcmc_kernel,
        parallel_iterations=10
    )

    mix_prob_sample, gpf_sample, sigma_sample, = state

    # set up init op
    init_op = tf.global_variables_initializer()

    mcmc_graph.finalize()

""" 2.2. execute sampling"""
with tf.Session(graph=mcmc_graph) as sess:
    init_op.run()
    [mix_prob_sample_val,
     f_samples_full_val,
     sigma_sample_val,
     is_accepted_,
     ] = sess.run(
        [
            mix_prob_sample,
            gpf_sample,
            sigma_sample,
            kernel_results.inner_results.is_accepted,
        ])
    print('Acceptance Rate: {}'.format(np.mean(is_accepted_)))
    sess.close()

""" 2.3. prediction and visualization"""
# prediction by sample from posterior
f_sample_indices = np.asarray([
    np.argmax(np.random.multinomial(1, pvals=mix_prob))
    for mix_prob in mix_prob_sample_val
])

f_samples_val = np.asarray([f_samples_full_val[i, f_sample_indices[i]]
                            for i in range(num_results)])

f_valid_val = gp.sample_posterior_full(X_new=X_valid, X=X_train,
                                       f_sample=f_samples_val.T,
                                       ls=ls_val, kernel_func=gp.rbf)

# visualize
mu = np.mean(f_valid_val, axis=1)
cov = np.var(f_valid_val, axis=1)

visual_util.gpr_1d_visual(mu, cov, pred_quantiles=[],
                          X_train=X_train, y_train=y_train,
                          X_test=X_valid, y_test=y_valid,
                          title="GP Mixture Posterior, {} Mixtures, Hamilton MC".format(n_mix_val),
                          save_addr=os.path.join(_SAVE_ADDR_PREFIX,
                                                 "gpr_hmc_mixture.png"))

visual_util.gpr_1d_visual(pred_mean=None, pred_cov=None, pred_quantiles=[],
                          pred_samples=list(f_valid_val.T)[:10000],
                          X_train=X_train, y_train=y_train,
                          X_test=X_valid, y_test=y_valid,
                          title="GP Mixture Posterior Samples, {} Mixtures, Hamilton MC".format(n_mix_val),
                          save_addr=os.path.join(_SAVE_ADDR_PREFIX,
                                                 "gpr_hmc_mixture_sample.png")
                          )

"""""""""""""""""""""""""""""""""
# 3. Mean-field VI
"""""""""""""""""""""""""""""""""

""" 3.1. Set up the computational graph """
mfvi_graph = tf.Graph()
ls_val = 0.1

with mfvi_graph.as_default():
    # sample from variational family
    q_f, q_sig, qf_mean, qf_sdev = gp_regression.variational_mfvi(X=X_train)

    # compute the expected predictive log-likelihood
    with ed.tape() as model_tape:
        with ed.interception(make_value_setter(gp_f=q_f, sigma=q_sig)):
            _, _, y = gp_regression.model(X=X_train, ls=ls_val)

    log_likelihood = y.distribution.log_prob(y_train)

    # compute the KL divergence
    kl = 0.
    for rv_name, variational_rv in [("gp_f", q_f), ("sigma", q_sig)]:
        kl += tf.reduce_sum(
            variational_rv.distribution.kl_divergence(
                model_tape[rv_name].distribution)
        )

    # define loss op: ELBO = E_q(p(x|z)) + KL(q || p)
    elbo = tf.reduce_mean(log_likelihood - kl)

    # define optimizer
    optimizer = tf.train.AdamOptimizer(5e-3)
    train_op = optimizer.minimize(-elbo)

    # define init op
    init_op = tf.global_variables_initializer()

    mfvi_graph.finalize()

""" 3.2. execute optimization """
max_steps = 20000  # number of training iterations

with tf.Session(graph=mfvi_graph) as sess:
    start_time = time.time()

    sess.run(init_op)
    for step in range(max_steps):
        start_time = time.time()
        _, elbo_value = sess.run([train_op, elbo])
        if step % 1000 == 0:
            duration = time.time() - start_time
            print("Step: {:>3d} Loss: {:.3f} ({:.3f} sec)".format(
                step, elbo_value, duration))
    qf_mean_val, qf_sdev_val = sess.run([qf_mean, qf_sdev])

    sess.close()

""" 3.3. prediction & visualization """
with tf.Session() as sess:
    f_samples = gp_regression.variational_mfvi_sample(n_sample=10000,
                                                      qf_mean=qf_mean_val,
                                                      qf_vcov=qf_sdev_val)
    f_samples_val = sess.run(f_samples)

# still use exact posterior predictive
f_test_val = gp.sample_posterior_full(X_new=X_test, X=X_train,
                                      f_sample=f_samples_val.T,
                                      ls=ls_val, kernel_func=gp.rbf)

# visualize
mu = np.mean(f_test_val, axis=1)
cov = np.var(f_test_val, axis=1)

visual_util.gpr_1d_visual(mu, cov,
                          X_train=X_train, y_train=y_train,
                          X_test=X_test, y_test=y_test,
                          title="RBF, Mean-field VI",
                          save_addr="./result/gpr/gpr_mfvi.png")

"""""""""""""""""""""""""""""""""
# 4. Sparse GP (Structured VI)
"""""""""""""""""""""""""""""""""

""" 3.1. Set up the computational graph """
sgp_graph = tf.Graph()
ls_val = 0.1
X_induce = np.expand_dims(np.linspace(np.min(X_train),
                                      np.max(X_train), 10), 1).astype(np.float32)

with sgp_graph.as_default():
    # sample from variational family

    (q_f, q_sig, qf_mean, qf_cov,
     Sigma_pre, S, Kxx, Kxz,
     Kzz, Kzz_inv, Kxz_Kzz_inv) = gp_regression.variational_sgpr(X=X_train,
                                                                 Z=X_induce,
                                                                 ls=ls_val)

    # compute the expected predictive log-likelihood
    with ed.tape() as model_tape:
        with ed.interception(make_value_setter(gp_f=q_f, sigma=q_sig)):
            _, _, y = gp_regression.model(X=X_train, ls=ls_val)

    log_likelihood = y.distribution.log_prob(y_train)

    # compute the KL divergence
    kl = 0.
    for rv_name, variational_rv in [("gp_f", q_f), ("sigma", q_sig)]:
        kl += tf.reduce_sum(
            variational_rv.distribution.kl_divergence(
                model_tape[rv_name].distribution)
        )

    # define loss op: ELBO = E_q(p(x|z)) + KL(q || p)
    elbo = tf.reduce_mean(log_likelihood - kl)

    # define optimizer
    optimizer = tf.train.AdamOptimizer(5e-3)
    train_op = optimizer.minimize(-elbo)

    # define init op
    init_op = tf.global_variables_initializer()

    sgp_graph.finalize()

""" 3.2. execute optimization """
max_steps = 50000  # number of training iterations

with tf.Session(graph=sgp_graph) as sess:
    start_time = time.time()

    sess.run(init_op)
    for step in range(max_steps):
        start_time = time.time()
        (Sigma_pre_val, S_val, Kxx_val, Kxz_val,
         Kzz_val, Kzz_inv_val, Kxz_Kzz_inv_val) = sess.run([
            Sigma_pre, S, Kxx, Kxz, Kzz, Kzz_inv, Kxz_Kzz_inv])
        _, elbo_value = sess.run([train_op, elbo])
        if step % 1000 == 0:
            duration = time.time() - start_time
            print("Step: {:>3d} Loss: {:.3f} ({:.3f} sec)".format(
                step, elbo_value, duration))
    qf_mean_val, qf_cov_val = sess.run([qf_mean, qf_cov])

    sess.close()

""" 3.3. prediction & visualization """
with tf.Session() as sess:
    f_samples = gp_regression.variational_sgpr_sample(n_sample=10000,
                                                      qf_mean=qf_mean_val,
                                                      qf_vcov=qf_cov_val)
    f_samples_val = sess.run(f_samples)

# still use exact posterior predictive
f_test_val = gp.sample_posterior_full(X_new=X_test, X=X_train,
                                      f_sample=f_samples_val.T,
                                      ls=ls_val, kernel_func=gp.rbf)

# visualize
mu = np.mean(f_test_val, axis=1)
cov = np.var(f_test_val, axis=1)

visual_util.gpr_1d_visual(mu, cov,
                          X_train=X_train, y_train=y_train,
                          X_test=X_test, y_test=y_test,
                          X_induce=X_induce,
                          title="RBF, Structured VI (Sparse GP)",
                          save_addr="./result/gpr/gpr_sgp.png")

"""""""""""""""""""""""""""""""""
# 5. Decoupled Sparse GP
"""""""""""""""""""""""""""""""""
"""Implements the De-coupled Sparse GP VI method by [2].

Select two sets of inducing points Z_m, Z_s 
    for mean and covariance separately, then:
    
Consequently, q(F) becomes 
    q(F|a, B) ~ MVN(F| Mu =     Kxz_m a
                       Sigma =  Kxx - Kxz_s (B^{-1} + K_zz_s)^{-1} Kxz_s^T)

In practice, parametrize B = LL^T using the Cholesky decomposition, then:
 
    (B^{-1} + K_zz_s)^{-1} = L H^{-1} L^T, where
    H = I + L^T K_zz_s L

"""

"""""""""""""""""""""""""""""""""""""""""""""
# Appendix. MAP inference using GPflow-slim
"""""""""""""""""""""""""""""""""""""""""""""
# define computation graph
tf.reset_default_graph()
kern_func = gpf.kernels.RBF
kern_pars = {'ARD': False, 'lengthscales': 1.}

mu, var, par_val, m, k = gpflow_util.fit_gpflow(X_train, y_train, X_test, y_test,
                                                kern_func=kern_func, **kern_pars)
# visualization
visual_util.gpr_1d_visual(mu, var,
                          X_train, y_train, X_test, y_test,
                          title="RBF, MAP, GPflow Implementation",
                          save_addr="./result/gpr/gpr_gpflow.png")
