library(dplyr)
library(rstan)
library(ggplot2)
library(cowplot)

# set working directory
if(file.exists('C:/Users/Admin-Dell')){
  root_dir = 'C:/Users/Admin-Dell'
}else{
  root_dir = 'C:/Users/nickl'
}
root_results <- sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Results/',
                        root_dir)

root_git <- sprintf('%s/Documents/github_projects/pop_ensemble/', 
                                  root_dir)

setwd(root_git)

# load extra functions
source('code/extra_functions_CAR.R')

#### Intercept only results pt 2 ####
setwd(root_results)
load('real_data/real_data_fit_directest_interceptonly_ID70259_2024_11_12.RData')
p1 <- plot_real_results(data_list = res$sim_list$data_list,
                        stan_fit = res$sim_list$stan_fit,
                        stan_summary = res$sim_list$stan_summary$summary,
                        CV_pred = res$sim_list$CV_pred,
                        alpha_estimates = F)
ggsave(plot = p1, filename = '../Figures/11132024_fullpop_directest_interceptonly_real_data.png', height = 12, width = 7)
tmp <- res$sim_list$stan_summary$summary
tmp[grep('beta',rownames(tmp)),]
# median -0.3, 0.6, 0.3

load('real_data/real_data_fit_softmax_alpha_interceptonly_ID24617_2024_11_12.RData')
p1 <- plot_real_results(data_list = res$sim_list$data_list,
                        stan_fit = res$sim_list$stan_fit,
                        stan_summary = res$sim_list$stan_summary$summary,
                        CV_pred = res$sim_list$CV_pred,
                        alpha_estimates = F)
ggsave(plot = p1, filename = '../Figures/11132024_fullpop_softmax_alpha_interceptonly_real_data.png', height = 12, width = 7)
tmp <- res$sim_list$stan_summary$summary
tmp[grep('beta',rownames(tmp)),]
# median -1.3, 2.9, -1.6




#
#### Intercept only results ####
setwd(root_results)
load('../real_data/real_data_fit_softmax_preprocess_alpha_interceptonly_ID45884_2024_11_12.RData')

tt <- res$sim_list$stan_fit
ss <- extract(tt)

p1 <- plot_real_results(data_list = res$sim_list$data_list,
                        stan_fit = res$sim_list$stan_fit,
                        stan_summary = res$sim_list$stan_summary$summary,
                        CV_pred = res$sim_list$CV_pred,
                        alpha_estimates = F)

#### Investigating US American Indian results ####
# (1) Look at the distribution of alpha for PEP.
# (2) Look at a single prediction from a single high alpha value - is it correct?

setwd(root_results)
load('../real_data/real_data_fit_softmax_preprocess_alpha_aian_ID97760_2024_10_28.RData')

tt <- params$raw_data$data
head(tt)
colMeans(tt[,3:6])
# Ah yes the preprocessing occurs within the fitting function, so the preprocessed data is not saved.

tt <- res$sim_list$stan_fit
alpha <- extract(tt, 'alpha')[[1]]
colMeans(alpha)
apply(alpha, 2, median)
# what? Wtf is being plotted? Ah it was tau2. Dingus.

#
#### Results on US American Indian - four dif. models ####
setwd(root_results)
setwd('../')

# softmax
load('real_data/real_data_fit_softmax_aian_ID82309_2024_10_28.RData')
p1 <- plot_real_results(res$sim_list$data_list, res$sim_list$stan_fit, CV_pred = res$sim_list$CV_pred, alpha_estimates = F)
ggsave(plot = p1, filename = '../Figures/10312024_AmericanIndianAN_softmax_real_data.png', height = 12, width = 7)


# softmax alpha
load('real_data/real_data_fit_softmax_alpha_aian_ID7045_2024_10_28.RData')
p1 <- plot_real_results(res$sim_list$data_list, res$sim_list$stan_fit, CV_pred = res$sim_list$CV_pred, alpha_estimates = T)
ggsave(plot = p1, filename = '../Figures/10312024_AmericanIndianAN_softmax_alpha_real_data.png', height = 12, width = 7)

# softmax preprocess
load('real_data/real_data_fit_softmax_preprocess_aian_ID41695_2024_10_28.RData')
p1 <- plot_real_results(res$sim_list$data_list, res$sim_list$stan_fit, CV_pred = res$sim_list$CV_pred, alpha_estimates = F)
ggsave(plot = p1, filename = '../Figures/10312024_AmericanIndianAN_softmax_preprocess_real_data.png', height = 12, width = 7)

# softmax preprocess alpha
load('real_data/real_data_fit_softmax_preprocess_alpha_aian_ID97760_2024_10_28.RData')
p1 <- plot_real_results(res$sim_list$data_list, res$sim_list$stan_fit, CV_pred = res$sim_list$CV_pred, alpha_estimates = T)
ggsave(plot = p1, filename = '../Figures/10312024_AmericanIndianAN_softmax_preprocess_alpha_real_data.png', height = 12, width = 7)

#
#### Results on TX with high sample number ####
setwd(root_results)
setwd('../')

# TX, centering X, estimating alpha.
load('real_data/real_data_fit_softmax_preprocess_alpha_ID37706_2024_10_08.RData')
p1 <- plot_real_results(res$sim_list$data_list, res$sim_list$stan_fit, CV_pred = res$sim_list$CV_pred, alpha_estimates = T)
ggsave(plot = p1, filename = '../Figures/10082024_TXn20k_softmax_centeringX_alpha_real_data.png', height = 12, width = 7)
#
#### Revisiting simulations - making plots for paper ####
# Firstly, I need to redo these sims anyway because the results I am showing are for fixed u across simulations. Oh well. I am not going to redo them just yet, though.
setwd(root_results)
files <- grep('08_23', dir(root_results), value = T)
comparison <- 'simulation_softmax_normal_3models_CV0_ID769482_2024_05_16/'
compare_parameters(comparison, files[1])
files <- files[c(4,8)]

for(i in 1:length(files)){
  out_name <- sprintf('../Figures/08232024_forPaper_%s.png', stringr::str_match(files[i], 'simulation_(.*?)_negbin')[2])
  generate_metrics_list(files[i]) %>%
    plot_metrics(include_MAP_rank = T) #%>%
    #ggsave(., filename = out_name, height = 8, width = 6)
}


if(T){
  tt <- generate_metrics_list(files[i])
  df_list <- tt$single_sim$data_list
  df_list$data$y2 <- NULL
  # plot the y predictions for a single simulation.
  y <- tt$single_sim$data_list$data$y
  
  p_yfit <- process_results(df_list, 
                            CV_pred = tt$single_sim$CV_pred, 
                            stan_fit = tt$single_sim$stan_fit, 
                            ESS = F, likelihoods = F, rho_estimates = F, tau2_estimates = F, sigma2_estimates = F, phi_estimates = F, u_estimates = F, RMSE_CP_values = F, 
                            y_estimates = T)
  
  ggsave(plot = p_yfit, filename = '../../Figures/10072024_singlerun_plots.pdf',height = 2, width = 8)
}

# 
#### Results with pre-processing and whatnot. ####
# REDOING ON 11/3/2024 with updated plotting function.
setwd(root_results)
setwd('../')

# softmax, centering the X values, estimating alpha.
load('real_data/real_data_fit_softmax_preprocess_alpha_ID57152_2024_09_16.RData')
p1 <- plot_real_results(res$sim_list$data_list, res$sim_list$stan_fit, CV_pred = res$sim_list$CV_pred, alpha_estimates = T)
ggsave(plot = p1, filename = '../Figures/11032024_softmax_centeringX_alpha_real_data.png', height = 12, width = 7)

# softmax, centering the X values, no alpha.
load('real_data/real_data_fit_softmax_preprocess_ID26230_2024_09_16.RData')
p1 <- plot_real_results(res$sim_list$data_list, res$sim_list$stan_fit, CV_pred = res$sim_list$CV_pred)
ggsave(plot = p1, filename = '../Figures/11032024_softmax_centeringX_NOalpha_real_data.png', height = 12, width = 7)

# softmax, no centering, alpha.
load('real_data/real_data_fit_softmax_alpha_ID31198_2024_09_16.RData')
p1 <- plot_real_results(res$sim_list$data_list, res$sim_list$stan_fit, CV_pred = res$sim_list$CV_pred, alpha_estimates = T)
ggsave(plot = p1, filename = '../Figures/11032024_softmax_NOcenteringX_alpha_real_data.png', height = 12, width = 7)

# softmax, no centering, no alpha.
load('real_data/real_data_fit_softmax_parallel_ID99987_2024_09_16.RData')
p1 <- plot_real_results(res$sim_list$data_list, res$sim_list$stan_fit, CV_pred = res$sim_list$CV_pred)
ggsave(plot = p1, filename = '../Figures/11032024_softmax_NOcenteringX_NOalpha_real_data.png', height = 12, width = 7)

# direct estimate, centering the X values, no alpha.
load('real_data/real_data_fit_directest_preprocess_ID84612_2024_09_16.RData')
p1 <- plot_real_results(res$sim_list$data_list, res$sim_list$stan_fit, CV_pred = res$sim_list$CV_pred)
ggsave(plot = p1, filename = '../Figures/11032024_directest_centeringX_NOalpha_real_data.png', height = 12, width = 7)

# direct estimate, no centering, no alpha.
load('real_data/real_data_fit_null_ID82460_2024_09_16.RData')
p1 <- plot_real_results(res$sim_list$data_list, res$sim_list$stan_fit, CV_pred = res$sim_list$CV_pred)
ggsave(plot = p1, filename = '../Figures/11032024_directest_NOcenteringX_NOalpha_real_data.png', height = 12, width = 7)


#
#### Loading test parallelized models ####
setwd(root_results)
load('real_data_fit_test_ID77075_2024_09_12.RData')
tt <- res$sim_list
ss <- tt$CV_pred
# ok so it's really that easy! Jeez.

p1 <- plot_real_results(res$sim_list$data_list, res$sim_list$stan_fit, CV_pred = res$sim_list$CV_pred)

#
#### Results 9/9/2024 - NB real data - softmax ####
setwd(root_results)
load('real_data_fit_softmax_ID15462_2024_09_09.RData')

p1 <- plot_real_results(res$sim_list$data_list, res$sim_list$stan_fit, CV_pred = res$sim_list$CV_pred)

# ggsave(plot = p1, filename = '../Figures/09092024_softmax_real_data.png', height = 12, width = 7)

# with alpha
load('real_data_fit_softmax_alpha_ID42629_2024_09_09.RData')

p1 <- plot_real_results(res$sim_list$data_list, res$sim_list$stan_fit, CV_pred = res$sim_list$CV_pred, alpha_estimates = T)

# ggsave(plot = p1, filename = '../Figures/09092024_softmax_alpha_real_data.png', height = 12, width = 7)


#
#### Results 9/9/2024 - NB real data - direct est ####
# burnin = 2k. n.sample = 4k
setwd(root_results)
load('real_data_fit_direct_ID63941_2024_09_06.RData')

p1 <- plot_real_results(res$sim_list$data_list, res$sim_list$stan_fit, CV_pred = res$sim_list$CV_pred)

# ggsave(plot = p1, filename = '../Figures/09092024_direct_est_real_data.png', height = 12, width = 7)

# why are the phi's and u's not equal? Dingus it's because I add 1/3
stan_fit <- res$sim_list$stan_fit
stan_summary = summary(stan_fit)$summary
stan_out <- extract(stan_fit)

#
#### Results 9/6/2024 - NB real data - softmax? ####
setwd(root_results)
load('real_data_fit_softmax_ID17016_2024_09_05.RData')

plot_real_results(res$sim_list$data_list, res$sim_list$stan_fit, CV_pred = res$sim_list$CV_pred)

#
#### Results 9/3/2024 - Negative binomial real data ####
load('results/real_data_results/real_data_fit_2024_09_03(1).RData')

plot_real_results(res$sim_list$data_list, res$sim_list$stan_fit, CV_pred = res$sim_list$CV_pred)

#
#### Results 8/30/2024 - Negative binomial with resampled phi ####
setwd(root_results)
setwd('simulated_results/')
files <- grep('08_29', dir(), value = T)

out_name <- sprintf('../Figures/08292024_%s_resampled_phi.png', stringr::str_match(files, 'simulation_(.*?)_negbin')[2])
generate_metrics_list(files) %>%
  plot_metrics(include_MAP_rank = F) #%>%
  #ggsave(., filename = out_name, height = 8, width = 6)

# ^ NEED TO CHECK ABOUT THE METRICS OF COVERAGE 

#
#### Results 8/26/2024 - Negative binomial ####
setwd(root_results)
files <- grep('08_23', dir(root_results), value = T)
comparison <- 'simulation_softmax_normal_3models_CV0_ID769482_2024_05_16/'
compare_parameters(comparison, files[1])

# making all the plots!
for(i in 1:length(files)){
  out_name <- sprintf('../Figures/08232024_%s.png', stringr::str_match(files[i], 'simulation_(.*?)_negbin')[2])
  generate_metrics_list(files[i]) %>%
    plot_metrics(include_MAP_rank = T) %>%
    ggsave(., filename = out_name, height = 8, width = 6)
}


#
#### Results 8/07/2024 - Lower sigma2 values ####
setwd(root_results)
files <- grep('08_06', dir(root_results), value = T)
comparison <- 'simulation_softmax_normal_3models_CV0_ID769482_2024_05_16/'
compare_parameters(comparison, files[1])
compare_parameters(files[1], files[2])

# making all the plots!
for(i in 1:length(files)){
  out_name <- sprintf('../Figures/08062024_%s.png', stringr::str_match(files[i], 'simulation_(.*?)_normal')[2])
  generate_metrics_list(files[i]) %>%
    plot_metrics(include_MAP_rank = T) %>%
    ggsave(., filename = out_name, height = 8, width = 6)
  
}

load("simulation_sigma2eq25_rho03_direct_est_normal_3models_CV10_ID203001_2024_08_06/sim_results_1.RData")
tmp <- res_lst[[1]]$sim_list[[1]]
process_results(tmp$data_list, tmp$stan_fit, ESS = F, likelihoods = F, stan_fit_quantiles = T)
# ugh this is a pain in the butt.

for(i in 1:length(files)){
  load(dir(files[i], full.names = T)[1])
  tmp <- res_lst[[1]]$sim_list[[1]]
  out_name <- sprintf('../Figures/08062024_%s.png', stringr::str_match(files[i], 'simulation_(.*?)_normal')[2])
  print(out_name)
  print(tmp$stan_fit[,1])
}

#
#### Results 6/05/2024 - *ICAR* with different true rho values ####
setwd(root_results)
recent_files(l = 6)

direct_rho03 <- 'simulation_fixedrho099_rho03_direct_est_normal_3models_CV10_ID986502_2024_06_04/'
comparison <- 'simulation_softmax_normal_3models_CV0_ID769482_2024_05_16/'
compare_parameters(direct_rho03, comparison)

direct_rho07 <- 'simulation_fixedrho099_rho07_direct_est_normal_3models_CV10_ID890724_2024_06_04/'
direct_rho099 <- 'simulation_fixedrho099_rho099_direct_est_normal_3models_CV10_ID436361_2024_06_04/'
softmax_rho03 <- 'simulation_fixedrho099_rho03_softmax_normal_3models_CV10_ID723759_2024_06_04/'
softmax_rho07 <- 'simulation_fixedrho099_rho07_softmax_normal_3models_CV10_ID342875_2024_06_04/'
softmax_rho099 <- 'simulation_fixedrho099_rho099_softmax_normal_3models_CV10_ID987269_2024_06_04/'

generate_metrics_list(direct_rho03) %>%
  plot_metrics(include_MAP_rank = T) %>%
  ggsave(., filename = '../Figures/06052024_direct_est_rho03_fixedrho099.png', height = 8, width = 6)

generate_metrics_list(direct_rho07) %>%
  plot_metrics(include_MAP_rank = T) %>%
  ggsave(., filename = '../Figures/06052024_direct_est_rho07_fixedrho099.png', height = 8, width = 6)

generate_metrics_list(direct_rho099) %>%
  plot_metrics(include_MAP_rank = T) %>%
  ggsave(., filename = '../Figures/06052024_direct_est_rho099_fixedrho099.png', height = 8, width = 6)

generate_metrics_list(softmax_rho03) %>%
  plot_metrics(include_MAP_rank = T) %>%
  ggsave(., filename = '../Figures/06052024_softmax_rho03_fixedrho099.png', height = 8, width = 6)

generate_metrics_list(softmax_rho07) %>%
  plot_metrics(include_MAP_rank = T) %>%
  ggsave(., filename = '../Figures/06052024_softmax_rho07_fixedrho099.png', height = 8, width = 6)

generate_metrics_list(softmax_rho099) %>%
  plot_metrics(include_MAP_rank = T) %>%
  ggsave(., filename = '../Figures/06052024_softmax_rho099_fixedrho099.png', height = 8, width = 6)


#
#### Results 5/20/2024 - means = 0,0,0 ####
setwd(root_results)
recent_files()

mean0 <- 'simulation_mean0_softmax_normal_3models_CV20_ID612954_2024_05_20/'
comparison <- 'simulation_softmax_normal_3models_CV0_ID769482_2024_05_16/'
compare_parameters(mean0, comparison)

p_mean0 <- generate_metrics_list(mean0) %>%
  plot_metrics()
ggsave(plot = p_mean0, filename = '../Figures/05202024_softmax_CV20_means0.png', height = 8, width = 6)

#
#### Results 5/20/2024 - rho = 0.99 ####
setwd(root_results)
recent_files()

folder_softmax_rho099 <- 'simulation_rho099_softmax_normal_3models_CV0_ID956416_2024_05_17/'
folder_DE_rho099 <- 'simulation_rho099_direct_est_normal_3models_CV0_ID310314_2024_05_17/'
folder_compare <- 'simulation_rho07_direct_est_normal_3models_CV20_ID832432_2024_05_16/'
compare_parameters(folder_softmax_rho099, folder_compare)
compare_parameters(folder_softmax_rho099, folder_DE_rho099)
# good

metrics_softmax_rho099 <- generate_metrics_list(folder_softmax_rho099)
p_softmax_rho099 <- plot_metrics(metrics_softmax_rho099)
ggsave(plot = p_softmax_rho099, filename = '../Figures/05202024_softmax_CV0_rho099.png', height = 8, width = 6)

metrics_DE_rho099 <- generate_metrics_list(folder_DE_rho099)
p_DE_rho099 <- plot_metrics(metrics_DE_rho099)
ggsave(plot = p_DE_rho099, filename = '../Figures/05202024_direct_est_CV0_rho099.png', height = 8, width = 6)

#
#### Results 5/17/2024 ####
setwd(root_results)
recent_files()

folder_DE_CV0 <- 'simulation_direct_est_normal_3models_CV0_ID689308_2024_05_16/'
folder_softmax_CV0 <- 'simulation_softmax_normal_3models_CV0_ID769482_2024_05_16/'
folder_DE_rho07 <- 'simulation_rho07_direct_est_normal_3models_CV20_ID832432_2024_05_16/'
folder_softmax_rho07 <- 'simulation_rho07_softmax_normal_3models_CV20_ID431503_2024_05_16/'

compare_parameters(folder_DE_rho07, folder_DE_CV0)
compare_parameters(folder_softmax_rho07, folder_DE_CV0)
# good.

metrics_DE_CV0 <- generate_metrics_list(folder_DE_CV0)
p_DE_CV0 <- plot_metrics(metrics_DE_CV0)
ggsave(plot = p_DE_CV0, filename = '../Figures/05162024_direct_est_CV0.png', height = 8, width = 6)

metrics_softmax_CV0 <- generate_metrics_list(folder_softmax_CV0)
p_softmax_CV0 <- plot_metrics(metrics_softmax_CV0)
ggsave(plot = p_softmax_CV0, filename = '../Figures/05162024_softmax_CV0.png', height = 8, width = 6)

metrics_DE_rho07 <- generate_metrics_list(folder_DE_rho07)
p_DE_rho07 <- plot_metrics(metrics_DE_rho07)
ggsave(plot = p_DE_rho07, filename = '../Figures/05162024_direct_est_rho07.png', height = 8, width = 6)

metrics_softmax_rho07 <- generate_metrics_list(folder_softmax_rho07)
p_softmax_rho07 <- plot_metrics(metrics_softmax_rho07)
ggsave(plot = p_softmax_rho07, filename = '../Figures/05162024_softmax_rho07.png', height = 8, width = 6)

#
#### Results 5/15/2024 ####
setwd(root_results)
folder_DE_K10 <- 'simulation_direct_est_normal_3models_CV10_ID660099_2024_05_14/'
metrics_DE_K10 <- generate_metrics_list(folder_DE_K10)
p_DE_K10 <- plot_metrics(metrics_DE_K10)
ggsave(plot = p_DE_K10, filename = '../Figures/05152024_direct_est_CV10.png', height = 8, width = 6)

folder_DE_K20 <- 'simulation_direct_est_normal_3models_CV20_ID365609_2024_05_14/'
metrics_DE_K20 <- generate_metrics_list(folder_DE_K20, debug_mode = T)
p_DE_K20 <- plot_metrics(metrics_DE_K20)
ggsave(plot = p_DE_K20, filename = '../Figures/05152024_direct_est_CV20.png', height = 8, width = 6)

folder_DE_K0 <- 'simulation_direct_est_normal_3models_CV0_ID862997_2024_05_14/'
metrics_DE_K0 <- generate_metrics_list(folder_DE_K0, debug_mode = T)
p_DE_K0 <- plot_metrics(metrics_DE_K0)
ggsave(plot = p_DE_K0, filename = '../Figures/05152024_direct_est_CV0.png', height = 8, width = 6)

folder_softmax_K10 <- 'simulation_softmax_normal_3models_CV10_ID958515_2024_05_15/'
metrics_softmax_K10 <- generate_metrics_list(folder_softmax_K10)
p_softmax_K10 <- plot_metrics(metrics_softmax_K10)
ggsave(plot = p_softmax_K10, filename = '../Figures/05152024_softmax_CV10.png', height = 8, width = 6)

folder_softmax_K20 <- 'simulation_softmax_normal_3models_CV20_ID999614_2024_05_15/'
metrics_softmax_K20 <- generate_metrics_list(folder_softmax_K20)
p_softmax_K20 <- plot_metrics(metrics_softmax_K20)
ggsave(plot = p_softmax_K20, filename = '../Figures/05152024_softmax_CV20.png', height = 8, width = 6)

#folder_softmax_K0 <- 'simulation_softmax_normal_3models_CV0_ID900288_2024_05_14/'
metrics_softmax_K0 <- generate_metrics_list(folder_softmax_K0)
p_softmax_K0 <- plot_metrics(metrics_softmax_K0)
ggsave(plot = p_softmax_K0, filename = '../Figures/05152024_softmax_CV0.png', height = 8, width = 6)

#
#### Results 5/13/2024 ####
folder_K10 <- 'simulation_direct_est_normal_3models_CV10_ID512833_2024_05_10'
metrics_K10 <- generate_metrics_list(folder_K10)
p_K10 <- plot_metrics(metrics_K10)

folder_K20 <- 'simulation_direct_est_normal_3models_CV20_ID915154_2024_05_10'
metrics_K20 <- generate_metrics_list(folder_K20)
p_K20 <- plot_metrics(metrics_K20)

folder_LOOCV <- 'simulation_direct_est_normal_3models_CV0_ID723970_2024_05_10'
metrics_LOOCV <- generate_metrics_list(folder_LOOCV)
p_LOOCV <- plot_metrics(metrics_LOOCV)

#### Checking if u's and phi's stay the same across simulations. ####
setwd(root_results)
load('softmax_tau21_tauprior11_CV5_05082024/sim_results_1.RData')

df1 = res_lst$sim_list[[1]]$data_list$data
u1 = res_lst$sim_list[[1]]$data_list$u_true
params1 = params

load('softmax_tau21_tauprior11_CV5_05082024/sim_results_2.RData')

df2 = res_lst$sim_list[[1]]$data_list$data
u2 = res_lst$sim_list[[1]]$data_list$u_true
params2 = params

sum(df1$y)
sum(df2$y)
# dif y's

sum(df1$y2)
sum(df2$y2)
# dif y2's

sum(u1[,1])
sum(u2[,1])
# same u's






