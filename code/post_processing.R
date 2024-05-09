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
folder <- 'softmax_tau21_tauprior11_CV5_05082024'
metrics_lst <- generate_metrics_list(folder)

### Plots the metrics across a set of simulations.
plot_metrics <- function(metrics_lst){
  # number of locations in dataset.
  n_loc <- length(metrics_lst[[1]]$CP_90_train)
  
  ## RMSE plot
  {
    RMSE_df <- NULL
    for(rmse in c('RMSE_train', 'RMSE_general', 'RMSE_CV')){
      tmp <- data.frame(val = sapply(metrics_lst, function(xx) xx[[rmse]]),
                        source = rmse)
      RMSE_df <- rbind(RMSE_df, tmp)
    }
    
    RMSE_df$source <- factor(RMSE_df$source, levels = c('RMSE_train', 'RMSE_general', 'RMSE_CV'))
    
    p_RMSE <- ggplot(data = RMSE_df, aes(x = source, y = val)) + 
      geom_boxplot()
  }
  
  ## y coverage plot
  {
    CP_90_train <- rowMeans(sapply(metrics_lst, function(xx) xx[['CP_90_train']]))
    y_CP_df <- data.frame(val = CP_90_train,
                          source = 'CP_train')
    
    CP_90_general <- rowMeans(sapply(metrics_lst, function(xx) xx[['CP_90_general']]))
    y_CP_df <- rbind(y_CP_df, data.frame(val = CP_90_general,
                                         source = 'CP_general'))
    y_CP_df$source <- factor(y_CP_df$source, levels = c('CP_train', 'CP_general'))
    
    p_y_CP <- ggplot(data = y_CP_df, aes(x = source, y = val)) + 
      geom_boxplot() + 
      ggtitle('y prediction 90% coverage')
  }
  
  ## phi and u coverage plots
  {
    # phi coverage plot
    phi_CP <- data.frame(CP = rowMeans(sapply(metrics_lst, function(xx) xx[['CP_90_phi']])))
    
    p_phi_CP <- ggplot(data = phi_CP, aes(y = CP)) + 
      geom_boxplot() + 
      ggtitle('phi 90% coverage')
    
    # u coverage plot
    u_CP <- data.frame(CP = rowMeans(sapply(metrics_lst, function(xx) xx[['CP_90_u']])))
    
    p_u_CP <- ggplot(data = u_CP, aes(y = CP)) + 
      geom_boxplot() + 
      ggtitle('u 90% coverage')
  }
  
  ## u rank plot
  u_rank <- data.frame(rank = rowMeans(sapply(metrics_lst, function(xx) xx[['rank_equal']])))
  
  p_u_rank <- ggplot(data = u_rank, aes(y = rank)) + 
    geom_boxplot() + 
    ggtitle('u-rank scores')
}


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






