library(dplyr)
library(rstan)
library(ggplot2)

rstan_options(auto_write = TRUE)

# set working directory
if(file.exists('C:/Users/Admin-Dell')){
  root_dir = 'C:/Users/Admin-Dell'
}else{
  root_dir = 'C:/Users/nickl'
}
root_results <- sprintf('%s/Dropbox/Academic/HSPH/Research/Population Estimation/Results',
                        root_dir)

root_git <- sprintf('%s/Documents/github_projects/pop_ensemble/', 
                                  root_dir)

setwd(root_git)

# load extra functions
source('code/extra_functions_CAR.R')

#### Master plotter ####

plot_all <- function(folder){
  metrics_lst <- generate_metrics_list(folder)
  metrics_plot <- plot_metrics(metrics_lst)
  
  
}

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

#### Testing metric list generation function. ####
folder <- 'softmax_tau21_tauprior11_CV5_05082024'
metrics_lst <- generate_metrics_list(folder)
plot_metrics(metrics_lst)

#
#### REDO THIS PLZ Checking if u's and phi's stay the same across simulations. ####
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






