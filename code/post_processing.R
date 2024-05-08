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

metrics_lst <- list()
iter <- 0
for(f in file_names){
  load(f)
  for(i in 1:length(res_lst$sim_list)){
    iter <- iter + 1
    tmp <- res_lst$sim_list[[i]]  
    medians <- tmp$stan_fit['0.5',]
    
    # pull out the y predictions
    ind_y_pred <- grep('y_pred', names(medians))
    median_y_pred <- medians[ind_y_pred]
    y_pred_05 <- tmp$stan_fit['0.05',ind_y_pred]
    y_pred_95 <- tmp$stan_fit['0.95',ind_y_pred]
    
    # pull out the y values
    y <- tmp$data_list$data$y
    y2 <- tmp$data_list$data$y2
    
    # # calculate the y metrics
    # RMSE_train = sqrt(mean((median_y_pred - y)^2))
    # RMSE_general = sqrt(mean((median_y_pred - y2)^2))
    # CP_90_train <- (y >= y_pred_05 & y <= y_pred_95)
    # CP_90_general <- (y2 >= y_pred_05 & y2 <= y_pred_95)
    
    # get true u and phi values
    phi_true <- tmp$data_list$phi_true
    u_true <- tmp$data_list$u_true
    phi_true_flat <- as.vector(as.matrix(phi_true[,-ncol(phi_true)]))
    u_true_flat <- as.vector(as.matrix(u_true[,-ncol(u_true)]))
    
    # get estimates u and phi values
    ind_phi <- grep('^phi\\[', colnames(tmp$stan_fit))
    ind_u <- grep('^u\\[', colnames(tmp$stan_fit))
    phi_est_05 <- tmp$stan_fit['0.05', ind_phi] 
    phi_est_95 <- tmp$stan_fit['0.95', ind_phi]
    median_phi <- medians[ind_phi]
    u_est_05 <- tmp$stan_fit['0.05', ind_u]
    u_est_95 <- tmp$stan_fit['0.95', ind_u]
    median_u_mat <- medians[ind_u] %>% 
      vec_to_mat(., n_models = 3)
    
    # # get coverage metrics
    # CP_90_phi <- (phi_true_flat >= phi_est_05 & phi_true_flat <= phi_est_95)
    # CP_90_u <- (u_true_flat >= u_est_05 & u_true_flat <= u_est_95)
    # 
    # # get the ranking metrics
    # rank_equal <- sapply(1:nrow(median_u_mat), function(xx){
    #   res <- all(rank(median_u_mat[xx,]) == rank(u_true[xx,-ncol(u_true)]))
    #   res
    # })
    
    metrics_lst[[iter]] <- list(RMSE_train = sqrt(mean((median_y_pred - y)^2)),
                                RMSE_general = sqrt(mean((median_y_pred - y2)^2)),
                                CP_90_train = (y >= y_pred_05 & y <= y_pred_95),
                                CP_90_general = (y2 >= y_pred_05 & y2 <= y_pred_95),
                                CP_90_phi = (phi_true_flat >= phi_est_05 & phi_true_flat <= phi_est_95),
                                CP_90_u = (u_true_flat >= u_est_05 & u_true_flat <= u_est_95),
                                rank_equal = sapply(1:nrow(median_u_mat), function(xx){
                                  res <- all(rank(median_u_mat[xx,]) == rank(u_true[xx,-ncol(u_true)]))
                                  res
                                })
    )
  }
}

### Convert vector of outputs to a matrix
# vec: Vector where values are arranged by column.
# n_models: number of models to split the vector by (e.g. number of rows of output matrix.)
vec_to_mat <- function(vec, n_models = 3){
  n <- length(vec)
  n_rows <- n/n_models
  
  # checking that the ordering is correcting (that it's by column, NOT row).
  tmp <- names(vec)[2]
  if(substr(tmp, nchar(tmp) - 1, nchar(tmp) - 1) != '1'){
    stop('ordering should be by column, but index doesnt match that')
  }
  
  # checking that the number of models is correct.
  tmp <- names(vec)[n]
  if(as.numeric(substr(tmp, nchar(tmp) - 1, nchar(tmp) - 1)) != n_models){
    stop('incorrect number of models')
  }
  
  mat <- matrix(NA, nrow = n_rows, ncol = n_models)
  for(j in 1:n_models){
    ind <- ((j - 1)*n_rows + 1):(j*n_rows)
    mat[,j] <- vec[ind]
  }
  
  return(mat)
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


