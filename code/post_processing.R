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

res_file <- 'softmax_tau21_tauprior11_CV5_05022024_v3'
file_names <- dir(sprintf('%s/%s', root_results, res_file), full.names = T)

for(f in file_names){
  load(f)
  for(i in 1:length(res_lst$sim_list)){
    tmp <- res_lst$sim_list[[i]]  
    medians <- tmp$stan_fit['0.5',]
    
    # pull out the y predictions
    ind_y_pred <- grep('y_pred', names(medians))
    median_y_pred <- medians[ind_y_pred]
    y_pred_05 <- tmp$stan_fit['0.05',ind_y_pred]
    y_pred_95 <- tmp$stan_fit['0.95',ind_y_pred]
    
    y <- tmp$data_list$data$y
    y2 <- tmp$data_list$data$y2
    RMSE_train = sqrt(mean((median_y_pred - y)^2))
    RMSE_general = sqrt(mean((median_y_pred - y2)^2))
    CP_90_train <- (y >= y_pred_05 & y <= y_pred_95)
    CP_90_general <- (y2 >= y_pred_05 & y2 <= y_pred_95)
    
  }
  
}

###
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


