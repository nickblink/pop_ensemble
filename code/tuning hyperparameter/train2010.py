#!/usr/bin/env python
# coding: utf-8

# ### Import Libraries

from wrapper_functions import *


os.getcwd()
os.listdir('../../data')



import pyreadr




#training2010 = pd.read_csv('../data/merged_wp_census_data_280922.csv')
training2010 = pyreadr.read_r('../../data/merged_fb_census_data_280922.RData') # also works for Rds

print(training2010.keys()) # let's check what objects we got
type(training2010)


# In[5]:


#training2010 = pyreadr.read_r('../data/merged_fb_census_data_280922.RData') # also works for Rds
# ^I can't get this library to install properly (Nick)
#training2010 = pd.read_csv('../data/merged_wp_census_data_280922.csv')

# objects
#print(training2010.keys()) # let's check what objects we got
training2010 = training2010["df"] # extract the pandas data frame for object df1

coordinate = np.asarray(training2010[["lon", "lat"]].values.tolist()).astype(np.float32)
training2010=training2010.fillna(0)


# In[6]:


models = ['acs', 'pep', 'worldpop','fb']
# base_preds_train = tf.stack([training2010[m] for m in models], axis=-1).astype(np.float32)
# base_preds_test = tf.stack([training2010[m] for m in models], axis=-1).astype(np.float32)


# standardize
X_train1 = np.asarray(training2010[["lon", "lat"]].values.tolist()).astype(np.float32)
X_test1 = np.asarray(training2010[["lon", "lat"]].values.tolist()).astype(np.float32)
X_valid = np.concatenate((X_train1, X_test1), axis=0)
X_centr = np.nanmean(X_valid, axis=0)
X_scale = np.nanmax(X_valid, axis=0) - np.nanmin(X_valid, axis=0)

X_train1 = (X_train1 - X_centr) / X_scale
X_test1 = (X_test1 - X_centr) / X_scale

Y_train = np.expand_dims(training2010["census"], 1).astype(np.float32)

Y_test = np.expand_dims(training2010["census"], 1).astype(np.float32)

base_preds_train = tf.stack([training2010[m].astype(np.float32) for m in models], axis=-1)
base_preds_test = tf.stack([training2010[m].astype(np.float32) for m in models], axis=-1)
print(base_preds_train.shape, base_preds_test.shape)


print("2010 center and scale: ", X_centr, X_scale)

import pandas as pd


# GP configs.
y_noise_std = 0.1  # @param
hidden_units = 1024  # @param
lengthscale=0.1  # @param
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

# Posterior configs.
bma_n_samples_train = 100 # @param
bma_n_samples_test = 200 # @param
bma_n_samples_eval = 1000  # @param

bma_seed = 0  # @param
bne_seed = 0 # @param

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


# ### Model Configs

# In[13]:


# Optimization configs. 
# Consider reduce below parameters / set to `False` if MCMC is taking too long:
# mcmc_num_steps, mcmc_burnin, mcmc_nchain, mcmc_initialize_from_map.
map_step_size=5e-4   # @param
map_num_steps=10_000  # @param

mcmc_step_size=1e-4 # @param
mcmc_num_steps=1000 # @param

mcmc_nchain=1 # @param
mcmc_burnin=100 # @param
bne_mcmc_initialize_from_map="True" # @param ["False", "True"]

bne_mcmc_initialize_from_map = eval(bne_mcmc_initialize_from_map)




bma_config=dict(gp_lengthscale=bma_gp_lengthscale,
                gp_l2_regularizer=bma_gp_l2_regularizer,
                y_noise_std=y_noise_std,
                map_step_size=map_step_size,
                map_num_steps=map_num_steps,
                mcmc_step_size=mcmc_step_size,
                mcmc_num_steps=mcmc_num_steps,
                mcmc_initialize_from_map=False,
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
                  mcmc_initialize_from_map=bne_mcmc_initialize_from_map,
                  seed=bne_seed)


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



# bma_prior, bma_gp_config = bma_dist(X_train1, 
#                                     base_preds_train, 
#                                     **bma_model_config)

# bma_model_config.update(bma_gp_config)

# # Check if the model graph is specified correctly.
# bma_prior.resolve_graph()


# bma_prior



print(map_config,mcmc_config,bma_model_config)


ref_model = LinearRegression()
def rmse(y_obs, y_pred):
    return np.sqrt(np.mean((y_obs - y_pred) ** 2))



rmse_lr = []
rmse_bma = []
rmse_gam = []

kf = KFold(n_splits=10, shuffle=True, random_state=bma_seed)

for train_index, test_index in kf.split(training2010):

    X_tr, X_te = X_train1[train_index], X_train1[test_index]
    Y_tr, Y_te = Y_train[train_index], Y_train[test_index]

    base_preds_tr, base_preds_te = base_preds_train.numpy()[train_index], base_preds_train.numpy()[test_index]
    #print(X_tr.shape, X_te.shape, Y_tr.shape, Y_te.shape, base_preds_tr.shape, base_preds_te.shape)

    # Ref: linear regression
    ref_model.fit(X_tr, Y_tr)
    Y_pred = ref_model.predict(X_te)
    rmse_lr.append(rmse(Y_te, Y_pred))
    #print(rmse_lr)

    # #GMA
    # ens_feature = np.asarray(list(base_preds_tr))
    # term_list = [s(dim_index) for dim_index in range(ens_feature.shape[1])]
    # term_list += [te(*list(ens_feature.shape[1] + np.array(range(X_tr.shape[1]))))]
    # gam_feature_terms = TermList(*term_list)

    # gam_X_tr = np.concatenate([X_tr, base_preds_tr], axis=1)
    # gam_X_te = np.concatenate([X_te, base_preds_te], axis=1)
    # ref_gam = LinearGAM(gam_feature_terms)
    # #ref_gam = LinearGAM(te(0, 1, 2) + te(0, 1, 3) + te(0,1,4))
    # gam = ref_gam.fit(gam_X_tr, Y_tr)
    # Y_pred = gam.predict(gam_X_te)
    # rmse_gam.append(rmse(Y_te, Y_pred))

    #BMA
    bma_prior, bma_gp_config = bma_dist(X_tr,
                                    base_preds_tr,
                                    **bma_model_config)

    bma_model_config.update(bma_gp_config)


    bma_gp_w_samples = run_posterior_inference(model_dist=bma_prior,
                                           model_config=bma_model_config,
                                           Y=Y_tr,
                                           map_config=map_config,
                                           mcmc_config=mcmc_config)


    bma_joint_samples = make_bma_samples(X_te, None, base_preds_te,
                                     bma_weight_samples=bma_gp_w_samples[0],
                                     bma_model_config=bma_model_config,
                                     n_samples=bma_n_samples_eval,
                                     seed=bne_seed,
                                     y_samples_only=False)

    y_pred = bma_joint_samples['y']
    y_pred = tf.reduce_mean(y_pred, axis=0)

    rmse_bma.append(rmse(Y_te, y_pred))

print('RMSE of LR:', np.mean(rmse_lr)),  print('RMSE of BMA:', np.mean(rmse_bma))




# bma_ensemble_weights = bma_joint_samples['ensemble_weights']



# ensemble_weights_val = tf.reduce_mean(bma_ensemble_weights, axis=0)
# #coordinate = np.asarray(training2010[["lon", "lat"]].values.tolist()).astype(np.float32)
# weights_dict = {
#     "acs": ensemble_weights_val[:, 0],
#     "pep": ensemble_weights_val[:,1],
#     "worldpop": ensemble_weights_val[:, 2],
#     "fb": ensemble_weights_val[:, 3]
# }

# color_norm_weights = make_color_norm(
#     list(weights_dict.values())[1],   
#     method="percentile")


# # In[63]:


# for model_name in models:
#     output = pd.DataFrame(np.column_stack([training2010[["GEOID"]], weights_dict[model_name]]))
#     output = output.set_axis(['GEOID', model_name], axis=1)
#     output[model_name] = output[model_name].astype(float)
#     fig = px.choropleth_mapbox(output, geojson=counties, locations='GEOID', color=model_name,
#                            color_continuous_scale="Viridis",
#                            #range_color=(0.05,0.07),
#                            mapbox_style="carto-positron",
#                            #featureidkey="properties.MWS_ID",
#                            zoom=3, center = {"lat": 37.0902, "lon": -95.7129},
#                            opacity=0.5,
#                            labels={'unemp':'unemployment rate'}
#                           )
#     fig.update_layout(margin={"r":0,"t":0,"l":0,"b":0})
#     fig.show()


# # ! conda install geopandas==0.3.0
# # ! conda install pyshp==1.2.10
# # ! conda install shapely==1.6.3
