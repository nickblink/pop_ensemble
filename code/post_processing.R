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

#### Testing metric list generation function. ####
folder <- 'softmax_tau21_tauprior11_CV5_05022024_v3'
metrics_lst <- generate_metrics_list(folder)


#### REDO THIS PLZ Checking if u's and phi's stay the same across simulations. ####
load('softmax_tau21_tauprior11_CV5_05022024_v3/sim_results_1.RData')

df1 = res_lst$sim_list[[1]]$data_list$data

load('softmax_tau21_tauprior11_CV5_05022024_v3/sim_results_2.RData')

df2 = res_lst$sim_list[[1]]$data_list$data

sum(df1$y)
sum(df2$y)
# dif y's

load('softmax_tau21_tauprior11_CV5_05022024_v3/sim_results_1.RData')

u1 = res_lst$sim_list[[1]]$data_list$u_true

load('softmax_tau21_tauprior11_CV5_05022024_v3/sim_results_2.RData')

u2 = res_lst$sim_list[[1]]$data_list$u_true

sum(u1[,1])
sum(u2[,1])
# dif u's. Change that!

# Code should take in a directory:
# For each file in the directory:
#    Load in the results
#    Compute the necessary metrics (RMSE, indicators for coverage)
#    

u1 = res_lst$sim_list[[1]]$data_list$u_true
u2 = res_lst$sim_list[[2]]$data_list$u_true

sum(u1[,1])
sum(u2[,1])

u1 = res_lst$sim_list[[1]]$data_list$data$y2
u2 = res_lst$sim_list[[2]]$data_list$data$y2

sum(u1)
sum(u2)


