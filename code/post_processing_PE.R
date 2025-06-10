library(dplyr)
library(tidyverse)
library(rstan)
library(ggplot2)
library(cowplot)
library(reshape2)
library(GGally)
library(posterior)

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

make_table_line <- function(metric, cols = c('dataset','MAPE', 'MAE', 'CP.95', 'med_int')) {
  res <- apply(metric, 1, function(xx) {
    # Make a copy to preserve original values
    xx_fmt <- xx
    
    # Apply rounding
    if ('MAPE' %in% names(xx_fmt)) xx_fmt['MAPE'] <- formatC(as.numeric(xx['MAPE']), format = 'f', digits = 2)
    if ('MAE' %in% names(xx_fmt)) xx_fmt['MAE'] <- as.character(round(as.numeric(xx['MAE'])))
    if ('CP.95' %in% names(xx_fmt)) xx_fmt['CP.95'] <- formatC(as.numeric(xx['CP.95']), format = 'f', digits = 3)
    if ('med_int' %in% names(xx_fmt)) xx_fmt['med_int'] <- as.character(round(as.numeric(xx['med_int'])))
    
    paste(xx_fmt[cols], collapse = ' & ')
  })
  
  return(res)
}

#
#### 6/10/2025: Getting full pop HMC diagnostics and results AND plotting SM ####
setwd(root_results)
files <- grep('04_22', dir('real_data', full.names = T), value = T)
files <- files[c(2,4)]

## Inspecting the convergence diagnostics
{
  # Initialize a data frame to hold diagnostics
  diag_summary <- tibble(
    file = character(),
    dataset = character(),
    softmax = logical(),
    preprocess = logical(),
    alpha = logical(),
    effects = character(),
    n_divergent = integer(),
    mean_divergent = numeric(),
    max_treedepth_hit = logical(),
    bfmi_low = numeric(),
    phi_noncentered = logical()
  )
  
  for (f in files) {
    load(f)
    fit <- res$sim_list$stan_fit
    
    if (!is.null(fit)) {
      sampler_params <- tryCatch(get_sampler_params(fit, inc_warmup = FALSE), error = function(e) NULL)
      
      if (!is.null(sampler_params) && length(sampler_params) > 0) {
        n_divergent <- sum(sapply(sampler_params, function(chain) sum(chain[, "divergent__"])))
        mean_divergent <- mean(sapply(sampler_params, function(chain) mean(chain[, "divergent__"])))
        #max_treedepth_hit <- any(sapply(sampler_params, function(chain) any(chain[, "treedepth__"] >= 10)))
        mean_treedepth_hit <- mean(sapply(sampler_params, function(chain) mean(chain[, "treedepth__"] >= 10)))
        bfmi_by_chain <- sapply(sampler_params, function(chain) {
          energy <- chain[, "energy__"]
          numer <- sum(diff(energy)^2) / (length(energy) - 1)
          denom <- var(energy)
          bfmi <- numer / denom
          bfmi
        })
        
        bfmi_low <- mean(bfmi_by_chain < 0.3)
        
        diag_summary <- add_row(diag_summary,
                                file = basename(f),
                                dataset = params$dataset, 
                                softmax = ifelse(params$use_softmax, T, F),
                                preprocess = ifelse(params$preprocess_scale, T, F),
                                alpha = ifelse(params$alpha_variance_prior == -1, F,T),
                                effects = params$fixed_effects,
                                n_divergent = n_divergent,
                                mean_divergent = mean_divergent,
                                max_treedepth_hit = mean_treedepth_hit,
                                bfmi_low = bfmi_low,
                                phi_noncentered = (params$phi_noncentered == 1)
        )
      } else {
        diag_summary <- add_row(diag_summary,
                                file = basename(f),
                                dataset = NA,
                                softmax = NA,
                                preprocess = NA,
                                alpha = NA,
                                effects = NA,
                                n_divergent = NA_integer_,
                                mean_divergent = NA,
                                max_treedepth_hit = NA,
                                bfmi_low = NA,
                                phi_noncentere = NA
        )
        warning(paste("No sampler params in", f, "- likely due to 0 samples"))
      }
    } else {
      warning(paste("Failed to load stan_fit in", f))
    }
  }
  
  # save(diag_summary, file = 'processed_results/real_data_fullpop_06102025.RData')
  
  # tt <- diag_summary[,c(2, 10, 7, 8, 9)]
}

load("real_data/real_data_fit_full_softmax_pepdensity_alpha0001_noncentered_ID73323_2025_04_22.RData")
res_SM <- just_metrics(
  data_list = res$sim_list$data_list,
  stan_fit = res$sim_list$stan_fit,
  stan_summary = res$sim_list$stan_summary$summary,
  models = params$models,
  CV_pred = res$sim_list$CV_pred
)
make_table_line(res_SM$metrics)

load("real_data/real_data_fit_full_directest_pepdensity_alpha0001_noncentered_ID44508_2025_04_22.RData")
res_DE <- just_metrics(
  data_list = res$sim_list$data_list,
  stan_fit = res$sim_list$stan_fit,
  stan_summary = res$sim_list$stan_summary$summary,
  models = params$models,
  CV_pred = res$sim_list$CV_pred
)
make_table_line(res_DE$metrics)

### Plot weights all together. (couldnt get this version to work)
{
  load("real_data/real_data_fit_full_softmax_pepdensity_alpha0001_noncentered_ID73323_2025_04_22.RData")
  SM_summary <- res$sim_list$stan_summary$summary
  
  load("real_data/real_data_fit_full_directest_pepdensity_alpha0001_noncentered_ID44508_2025_04_22.RData")
  DE_summary <- res$sim_list$stan_summary$summary
  
  N = nrow(res$sim_list$data_list$data)
  models = c('acs','pep','wp')
  
  u_est_SM <- as.data.frame(matrix(0, nrow = N, ncol = 3))
  colnames(u_est_SM) <- models
  for(i in 1:length(models)){
    ind = grep(sprintf('^u\\[[0-9]*,%s\\]', i), rownames(SM_summary))
    u_est_SM[,i] <- SM_summary[ind,'50%']
  }
  u_est_SM$index = 1:N
  # convert estimates to long
  u_est_SM_long <- tidyr::gather(u_est_SM, key = model, value = u_median_est, -index) %>%
    mutate(model = toupper(as.character(model))) # Capitalize entire model names
  
  u_est_DE <- as.data.frame(matrix(0, nrow = N, ncol = 3))
  colnames(u_est_DE) <- models 
  for(i in 1:length(models)){
    ind = grep(sprintf('^u\\[[0-9]*,%s\\]', i), rownames(SM_summary))
    u_est_DE[,i] <- DE_summary[ind,'50%']
  }
  u_est_DE$index = 1:N
  # convert estimates to long
  u_est_DE_long <- tidyr::gather(u_est_DE, key = model, value = u_median_est, -index) %>%
    mutate(model = toupper(as.character(model))) # Capitalize entire model names
  
  u_est_long <- bind_rows(
    u_est_SM_long %>% mutate(set = "SM"),
    u_est_DE_long %>% mutate(set = "DE")
  )
  
  medians_df <- u_est_long %>%
    group_by(set, model) %>%
    summarize(u_median_true = median(u_median_est, na.rm = TRUE), .groups = "drop")
}

### Plot weights original version
{
  load("real_data/real_data_fit_full_softmax_pepdensity_alpha0001_noncentered_ID73323_2025_04_22.RData")
  pweights_SM <- plot_real_results(data_list = res$sim_list$data_list,
                          stan_fit = res$sim_list$stan_fit,
                          stan_summary = res$sim_list$stan_summary$summary,
                          models = params$models,
                          CV_pred = res$sim_list$CV_pred, rhats = F,
                          alpha_estimates = F,
                          ESS = F, rho_estimates = F, tau2_estimates = F, 
                          sigma2_estimates = F, theta_estimates = F, phi_estimates = F,
                          pairwise_phi_estimates = F, y_estimates = F, metrics_values = F, beta_estimates = F)
  ggsave(pweights_SM, file = '../Figures/06102025_fullpop_u_estimates_softmax.png', width = 6, height = 2)
  
  load("real_data/real_data_fit_full_directest_pepdensity_alpha0001_noncentered_ID44508_2025_04_22.RData")
  pweights_DE <- plot_real_results(data_list = res$sim_list$data_list,
                                   stan_fit = res$sim_list$stan_fit,
                                   stan_summary = res$sim_list$stan_summary$summary,
                                   models = params$models,
                                   CV_pred = res$sim_list$CV_pred, rhats = F,
                                   alpha_estimates = F,
                                   ESS = F, rho_estimates = F, tau2_estimates = F, 
                                   sigma2_estimates = F, theta_estimates = F, phi_estimates = F,
                                   pairwise_phi_estimates = F, y_estimates = F, metrics_values = F, beta_estimates = F)
  
  ggsave(pweights_DE, file = '../Figures/06102025_fullpop_u_estimates_directest.png', width = 6, height = 2)
  
  p1 <- pweights_SM + labs(tag = "Softmax")
  p2 <- pweights_DE + labs(tag = "Direct Estimate")
  p1/p2
  ggsave(file = '../Figures/06102025_fullpop_u_estimates_combined.png', width = 6, height = 4)
  
}

### Plot Chloropleth.
{
  load("real_data/real_data_fit_full_softmax_pepdensity_alpha0001_noncentered_ID73323_2025_04_22.RData")
  p_SM <- plot_weights_map(
    data_list = res$sim_list$data_list,
    stan_summary = res$sim_list$stan_summary$summary,
    facet = TRUE,
    show_state_abbr = F
  )
  
  load("real_data/real_data_fit_full_directest_pepdensity_alpha0001_noncentered_ID44508_2025_04_22.RData")
  p_DE <- plot_weights_map(
    data_list = res$sim_list$data_list,
    stan_summary = res$sim_list$stan_summary$summary,
    facet = TRUE,
    show_state_abbr = F
  )
  
  p1 <- p_SM + theme(plot.margin = margin(0, 0, 0, 0)) + labs(tag = "Softmax") + theme(
    plot.tag.position = c(0.5, 3),  # middle-top
    plot.tag = element_text(hjust = 0, vjust = -2)
  )
  p2 <- p_DE + theme(plot.margin = margin(0, 0, 0, 0)) + labs(tag = "Direct Estimate") + theme(
    plot.tag.position = c(0.5, 3),  # middle-top
    plot.tag = element_text(hjust = 0, vjust = -2)
  )
  p1 | p2
  
  ggsave('../Figures/06102025_fullpop_combined_chloropleth.png', height = 10, width = 8)
}

#
#### 5/29/2025: Making AIAN u-weight and chloropleth plots ####
setwd(root_results)
setwd('real_data/')
load('real_data_fit_aiansubset_directest_intercept_noncentered_ID32362_2025_05_22.RData')

### Direct est first.
### Weights.
{
  pweights_DE <- plot_real_results(data_list = res$sim_list$data_list,
                          stan_fit = res$sim_list$stan_fit,
                          stan_summary = res$sim_list$stan_summary$summary,
                          models = params$models,
                          CV_pred = res$sim_list$CV_pred, rhats = F,
                          alpha_estimates = F,
                          ESS = F, rho_estimates = F, tau2_estimates = F, 
                          sigma2_estimates = F, theta_estimates = F, phi_estimates = F,
                          pairwise_phi_estimates = F, y_estimates = F, metrics_values = F, beta_estimates = F)
  # inspect plot. How is it?
  # ggsave(pweights_DE, file = '../../Figures/05252025_aiansubset_u_estimates_directest.png', width = 6, height = 3)
}

### Chloropleth.
{
  p_DE <- plot_weights_map(
    data_list = res$sim_list$data_list,
    stan_summary = res$sim_list$stan_summary$summary,
    facet = TRUE,
    xlim = c(-125, -100),
    ylim = c(31, 42)
  )
  
  #ggsave('../../Figures/05292025_aiansubset_directest_chloropleth.png', height = 12, width = 7)
}

### Softmax 
load('real_data_fit_aiansubset_softmax_pepdensity_alpha0001_noncentered_ID88451_2025_05_22.RData')
### Weights.
{
  pweights_SM <- plot_real_results(data_list = res$sim_list$data_list,
                          stan_fit = res$sim_list$stan_fit,
                          stan_summary = res$sim_list$stan_summary$summary,
                          models = params$models,
                          CV_pred = res$sim_list$CV_pred, rhats = F,
                          alpha_estimates = F,
                          ESS = F, rho_estimates = F, tau2_estimates = F, 
                          sigma2_estimates = F, theta_estimates = F, phi_estimates = F,
                          pairwise_phi_estimates = F, y_estimates = F, metrics_values = F, beta_estimates = F)
  # inspect plot. How is it?
  # ggsave(pweights_SM, file = '../../Figures/05252025_aiansubset_u_estimates_softmax.png', width = 6, height = 3)
}

### Chloropleth.
{
  p_SM <- plot_weights_map(
    data_list = res$sim_list$data_list,
    stan_summary = res$sim_list$stan_summary$summary,
    facet = TRUE,
    xlim = c(-125, -100),
    ylim = c(31, 42)
  )
  
 # ggsave('../../Figures/05292025_aiansubset_softmax_chloropleth.png', height = 12, width = 7)
}

### Combining density plots
{
  p1 <- pweights_SM + labs(tag = "Softmax")
  p2 <- pweights_DE + labs(tag = "Direct Estimate")
  p1/p2
  ggsave(file = '../../Figures/06102025_aian_u_estimates_combined.png', width = 6, height = 4)
}

### Combining chloropleths
{
  p1 <- p_SM + theme(plot.margin = margin(0, 0, 0, 0)) + labs(tag = "Softmax") + theme(
    plot.tag.position = c(0.5, 3),  # middle-top
    plot.tag = element_text(hjust = 0, vjust = -2)
  )
  p2 <- p_DE + theme(plot.margin = margin(0, 0, 0, 0)) + labs(tag = "Direct Estimate") + theme(
    plot.tag.position = c(0.5, 3),  # middle-top
    plot.tag = element_text(hjust = 0, vjust = -2)
  )
  p1 | p2
  
  ggsave('../../Figures/06102025_aian_combined_chloropleth.png', height = 10, width = 8)
}


### temporary testing for weights map function
setwd(root_results)
setwd('real_data')
load('Older fits/real_data_fit_directest_cv10_interceptonly_ID81515_2025_01_15.RData')

plot_weights_map(
  data_list = res$sim_list$data_list,
  stan_summary = res$sim_list$stan_summary$summary,
  show_state_abbr = F,
  facet = T)



#
#### 5/29/2025: Getting AIAN subset HMC diagnostics and results ####
setwd(root_results)
files <- grep('05_22', dir('real_data', full.names = T), value = T)

## Inspecting the convergence diagnostics
{
  # Initialize a data frame to hold diagnostics
  diag_summary <- tibble(
    file = character(),
    dataset = character(),
    softmax = logical(),
    preprocess = logical(),
    alpha = logical(),
    effects = character(),
    n_divergent = integer(),
    mean_divergent = numeric(),
    max_treedepth_hit = logical(),
    bfmi_low = numeric(),
    phi_noncentered = logical()
  )
  
  for (f in files) {
    load(f)
    fit <- res$sim_list$stan_fit
    
    if (!is.null(fit)) {
      sampler_params <- tryCatch(get_sampler_params(fit, inc_warmup = FALSE), error = function(e) NULL)
      
      if (!is.null(sampler_params) && length(sampler_params) > 0) {
        n_divergent <- sum(sapply(sampler_params, function(chain) sum(chain[, "divergent__"])))
        mean_divergent <- mean(sapply(sampler_params, function(chain) mean(chain[, "divergent__"])))
        #max_treedepth_hit <- any(sapply(sampler_params, function(chain) any(chain[, "treedepth__"] >= 10)))
        mean_treedepth_hit <- mean(sapply(sampler_params, function(chain) mean(chain[, "treedepth__"] >= 10)))
        bfmi_by_chain <- sapply(sampler_params, function(chain) {
          energy <- chain[, "energy__"]
          numer <- sum(diff(energy)^2) / (length(energy) - 1)
          denom <- var(energy)
          bfmi <- numer / denom
          bfmi
        })
        
        bfmi_low <- mean(bfmi_by_chain < 0.3)
        
        diag_summary <- add_row(diag_summary,
                                file = basename(f),
                                dataset = params$dataset, 
                                softmax = ifelse(params$use_softmax, T, F),
                                preprocess = ifelse(params$preprocess_scale, T, F),
                                alpha = ifelse(params$alpha_variance_prior == -1, F,T),
                                effects = params$fixed_effects,
                                n_divergent = n_divergent,
                                mean_divergent = mean_divergent,
                                max_treedepth_hit = mean_treedepth_hit,
                                bfmi_low = bfmi_low,
                                phi_noncentered = (params$phi_noncentered == 1)
        )
      } else {
        diag_summary <- add_row(diag_summary,
                                file = basename(f),
                                dataset = NA,
                                softmax = NA,
                                preprocess = NA,
                                alpha = NA,
                                effects = NA,
                                n_divergent = NA_integer_,
                                mean_divergent = NA,
                                max_treedepth_hit = NA,
                                bfmi_low = NA,
                                phi_noncentere = NA
        )
        warning(paste("No sampler params in", f, "- likely due to 0 samples"))
      }
    } else {
      warning(paste("Failed to load stan_fit in", f))
    }
  }
  
  # save(diag_summary, file = 'processed_results/real_data_aiansubset_05222025.RData')
  
  # tt <- diag_summary[,c(2, 10, 7, 8, 9)]
}

### Make a plot
load('real_data/real_data_fit_aiansubset_softmax_pepdensity_alpha0001_noncentered_ID88451_2025_05_22.RData')
p1 <- plot_real_results(data_list = res$sim_list$data_list,
                        stan_fit = res$sim_list$stan_fit,
                        stan_summary = res$sim_list$stan_summary$summary,
                        models = params$models,
                        CV_pred = res$sim_list$CV_pred,
                        alpha_estimates = T,
                        ESS = T, rho_estimates = T, tau2_estimates = T, 
                        sigma2_estimates = F, theta_estimates = T, phi_estimates = F,
                        pairwise_phi_estimates = T, y_estimates = F, metrics_values = T, beta_estimates = T)

res_SM <- just_metrics(
  data_list = res$sim_list$data_list,
  stan_fit = res$sim_list$stan_fit,
  stan_summary = res$sim_list$stan_summary$summary,
  models = params$models,
  CV_pred = res$sim_list$CV_pred
)
make_table_line(res_SM$metrics)

# ggsave(p1, file = '../Figures/05292025_real_data_fit_aiansubset_noncentered_softmax.png', height = 12, width = 7)

## Direct est
load('real_data/real_data_fit_aiansubset_directest_intercept_noncentered_ID32362_2025_05_22.RData')
p1 <- plot_real_results(data_list = res$sim_list$data_list,
                        stan_fit = res$sim_list$stan_fit,
                        stan_summary = res$sim_list$stan_summary$summary,
                        models = params$models,
                        CV_pred = res$sim_list$CV_pred,
                        alpha_estimates = F,
                        ESS = T, rho_estimates = T, tau2_estimates = T, 
                        sigma2_estimates = F, theta_estimates = T, phi_estimates = F,
                        pairwise_phi_estimates = T, y_estimates = F, metrics_values = T, beta_estimates = T)

res_SM <- just_metrics(
  data_list = res$sim_list$data_list,
  stan_fit = res$sim_list$stan_fit,
  stan_summary = res$sim_list$stan_summary$summary,
  models = params$models,
  CV_pred = res$sim_list$CV_pred
)
make_table_line(res_SM$metrics)

# ggsave(p1, file = '../Figures/05292025_real_data_fit_aiansubset_noncentered_directest.png', height = 12, width = 7)
#
#### 5/28/2025: Update simulation results figures and tables ####
setwd(root_results)
setwd('simulated_results/')
files1 <- grep('2025_05_22', dir(), value = T)
# because I repeated analysis
files <- files1[c(1,2,4,5)]

results_list <- lapply(files, function(f){
  generate_metrics_list(f)})

warning('Hardcoding names of results - make sure they match')
names(results_list) <- c('DE_rho03', 'SM_rho03', 'DE_rho099', 'SM_rho099')

### Making the table metrics 
{
  # Function to compute median and 95% quantiles
  summary_stats <- function(vec) {
    q_values <- quantile(vec, probs = c(0.025, 0.975), na.rm = TRUE)
    
    return(c(
      median = median(vec, na.rm = TRUE),
      Q2.5 = unname(q_values[1]),  # Ensure it doesn't inherit unwanted names
      Q97.5 = unname(q_values[2])
    ))
  }
  
  all_sim_results <- lapply(results_list, function(xx){
    metrics_list <- xx$metrics_list
    
    # Extract numeric metrics
    MAPE_train <- sapply(metrics_list, function(x) x$MAPE_train)
    MAPE_CV <- sapply(metrics_list, function(x) x$MAPE_CV)
    MAE_train <- sapply(metrics_list, function(x) x$MAE_train)
    MAE_CV <- sapply(metrics_list, function(x) x$MAE_CV)
    median_int_width_train <- sapply(metrics_list, function(x) x$median_int_width_train)
    median_int_width_CV <- sapply(metrics_list, function(x) x$median_int_width_CV)
    
    # Compute interval widths across all.
    median_int_width_train_across_all <- sapply(metrics_list, function(x) as.integer(x$int_widths_train, na.rm = TRUE)) %>% median()
    median_int_width_CV_across_all <- sapply(metrics_list, function(x) as.integer(x$int_widths_CV, na.rm = TRUE)) %>% median()
    
    # Compute proportion of TRUE values for Boolean vectors
    CP_95_train <- sapply(metrics_list, function(x) mean(x$CP_95_train, na.rm = TRUE))
    CP_90_train <- sapply(metrics_list, function(x) mean(x$CP_90_train, na.rm = TRUE))
    CP_95_CV <- sapply(metrics_list, function(x) mean(x$CP_95_CV, na.rm = TRUE))
    CP_90_CV <- sapply(metrics_list, function(x) mean(x$CP_90_CV, na.rm = TRUE))
    
    # get all CPs to join together
    CP_95_train_across_all <- sapply(metrics_list, function(x) as.integer(x$CP_95_train, na.rm = TRUE)) %>% mean()
    CP_90_train_across_all <- sapply(metrics_list, function(x) as.integer(x$CP_90_train, na.rm = TRUE)) %>% mean()
    CP_95_CV_across_all <- sapply(metrics_list, function(x) as.integer(x$CP_95_CV, na.rm = TRUE)) %>% mean()
    CP_90_CV_across_all <- sapply(metrics_list, function(x) as.integer(x$CP_90_CV, na.rm = TRUE)) %>% mean()
    
    # Compute summary statistics for all extracted values
    final_results <- list(
      MAPE_train = summary_stats(MAPE_train),
      MAPE_CV = summary_stats(MAPE_CV),
      MAE_train = summary_stats(MAE_train),
      MAE_CV = summary_stats(MAE_CV),
      median_int_width_train = summary_stats(median_int_width_train),
      median_int_width_CV = summary_stats(median_int_width_CV),
      median_int_width_train_across_all = median_int_width_train_across_all,
      median_int_width_CV_across_all = median_int_width_CV_across_all,
      CP_95_train = summary_stats(CP_95_train),
      CP_90_train = summary_stats(CP_90_train),
      CP_95_CV = summary_stats(CP_95_CV),
      CP_90_CV = summary_stats(CP_90_CV),
      CP_95_train_across_all = CP_95_train_across_all,
      CP_90_train_across_all = CP_90_train_across_all,
      CP_95_CV_across_all = CP_95_CV_across_all,
      CP_90_CV_across_all = CP_90_CV_across_all
    )
    
    final_results
  })
  
  latex_rows <- generate_latex_values(data_list = all_sim_results, coverage = 95, scale_across_all = T, digits = 2)
  cat(paste(latex_rows, collapse = " \\\\\n"))
}


#
#### 5/28/2025: Getting simulated metrics and diagnostics (again) ####
setwd(root_results)
setwd('simulated_results/')

grep('5_22', dir(), value = T)

# figure out diff in softmax params.
compare_parameters(folder1 = "simulation_rho_03_theta_100_softmax_negbin_3models_CV10_ID823748_2025_05_22",
                   folder2 = "simulation_rho_03_theta_100_softmax_negbin_3models_CV10_ID835785_2025_05_22/")
# ok so these are the same?

# figure out diff in softmax params.
compare_parameters(folder1 = "simulation_rho_099_theta_100_softmax_negbin_3models_CV10_ID290522_2025_05_22/",
                   folder2 = "simulation_rho_099_theta_100_softmax_negbin_3models_CV10_ID663244_2025_05_22/")
# ok so also the same? I guess I just repeated a simulation run.

res_DE_rho03 <- generate_metrics_list(folder = "simulation_rho_03_theta_100_direct_est_negbin_3models_CV10_ID28206_2025_05_22", hmc_diag = T)

res_SM_rho03 <- generate_metrics_list(folder = "simulation_rho_03_theta_100_softmax_negbin_3models_CV10_ID823748_2025_05_22", hmc_diag = T)

res_DE_rho099 <- generate_metrics_list(folder = "simulation_rho_099_theta_100_direct_est_negbin_3models_CV10_ID641351_2025_05_22/", hmc_diag = T)

res_SM_rho099 <- generate_metrics_list(folder = "simulation_rho_099_theta_100_softmax_negbin_3models_CV10_ID290522_2025_05_22/", hmc_diag = T)

res_all <- list(res_DE_rho03 = res_DE_rho03, 
                res_SM_rho03 = res_SM_rho03, 
                res_DE_rho099 = res_DE_rho099, 
                res_SM_rho099 = res_SM_rho099)

# save(res_all, file = '../processed_results/simulation_results_05222025.RData')

for(i in 1:length(res_all)){
  tmp <- res_all[[i]]
  print(sprintf('---%s---', names(res_all)[[i]]))
  bfmi <- sapply(tmp$metrics_list, function(xx) xx$bfmi_low) %>% sum()
  trees <- sapply(tmp$metrics_list, function(xx) xx$max_treedepth_hit) %>% sum()
  divs <- sapply(tmp$metrics_list, function(xx) xx$n_divergent) %>% sum()
  print(sprintf('low bfmis: %s, tree depths hit: %s, number divergences:%s', bfmi, trees, divs))
}
# out of 900 samples per sim X 200 sims. = 180000
21/(900*200)*100
1981/(900*200) # still just 1.1%


### Getting u-rank scores 
{
  # pull out the u-rank of each file
  u_rank_scores <- lapply(res_all, function(xx){
    rank <- colMeans(sapply(xx$metrics_list, function(yy) yy[['rank_equal']]))
    rank
  })
  
  # Define LaTeX-style expressions for x-axis labels
  group_labels <- c(
    "res_SM_rho03" = expression(SM ~ rho == 0.3),
    "res_SM_rho099" = expression(SM ~ rho == 0.99),
    "res_DE_rho03" = expression(DE ~ rho == 0.3),
    "res_DE_rho099" = expression(DE ~ rho == 0.99)
  )
  
  # Convert list to a tidy data frame
  df <- u_rank_scores %>%
    enframe(name = "Group", value = "Values") %>%
    unnest(Values)
  df$Group <- factor(df$Group, levels = names(group_labels))  # Ensure correct ordering
  
  # Compute the 2.5%, 50% (median), and 97.5% quantiles for each group
  df_summary <- df %>%
    group_by(Group) %>%
    summarise(
      lowest = min(Values),
      lower = quantile(Values, 0.025),   # 2.5th percentile
      middle = quantile(Values, 0.50),  # Median (50th percentile)
      upper = quantile(Values, 0.975),   # 97.5th percentile
      highest = max(Values),
      .groups = "drop"
    )
  
  # Create the customized boxplot
  ggplot(df, aes(x = Group, y = Values)) +
    # Use geom_segment() for whiskers
    geom_segment(data = df_summary, aes(x = Group, xend = Group, y = lowest, yend = highest), color = "black") +
    # Use geom_crossbar() to create the box
    geom_crossbar(data = df_summary, aes(x = Group, ymin = lower, y = middle, ymax = upper), fill = "white", color = "black") +
    # Add a horizontal reference line at y = 1/6
    geom_hline(yintercept = 1/6, linetype = "dashed", color = "red") +
    # Clean theme
    theme_minimal() +
    labs(x = NULL, y = "u-rank within simulation run", title = NULL) +
    scale_x_discrete(labels = group_labels) +  # Apply LaTeX-style labels
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          legend.position = "none")  # Remove legend
  
  ggsave(filename = '../../Figures/05282025_urank_boxplot.png', height = 5, width = 7)
}


#
#### 4/18/2025: Getting simulated metrics and diagnostics ####
setwd(root_results)
setwd('simulated_results/')

res_C_rho03 <- generate_metrics_list(folder = 'simulation_centered_rho_03_theta_100_softmax_negbin_3models_CV10_ID777773_2025_04_18/', hmc_diag = T)

res_C_rho099 <- generate_metrics_list(folder = 'simulation_centered_rho_099_theta_100_softmax_negbin_3models_CV10_ID855930_2025_04_18/', hmc_diag = T)

res_NC_rho03 <- generate_metrics_list(folder = 'simulation_noncentered_rho_03_theta_100_softmax_negbin_3models_CV10_ID329578_2025_04_18/', hmc_diag = T)

res_NC_rho099 <- generate_metrics_list(folder = 'simulation_noncentered_rho_099_theta_100_softmax_negbin_3models_CV10_ID609610_2025_04_18/', hmc_diag = T)

res_all <- list(center_03 = res_C_rho03, 
                center_099 = res_C_rho099, 
                noncenter_03 =res_NC_rho03, 
                noncenter_099 = res_NC_rho099)

for(i in 1:length(res_all)){
  tmp <- res_all[[i]]
  print(sprintf('---%s---', names(res_all)[[i]]))
  bfmi <- sapply(tmp$metrics_list, function(xx) xx$bfmi_low) %>% sum()
  trees <- sapply(tmp$metrics_list, function(xx) xx$max_treedepth_hit) %>% sum()
  divs <- sapply(tmp$metrics_list, function(xx) xx$n_divergent) %>% sum()
  print(sprintf('low bfmis: %s, tree depths hit: %s, number divergences:%s', bfmi, trees, divs))
}

## I didn't do any cross-validation but I can still get u-rank:

### Getting u-rank scores 
{
  # pull out the u-rank of each file
  u_rank_scores <- lapply(res_all, function(xx){
    rank <- colMeans(sapply(xx$metrics_list, function(yy) yy[['rank_equal']]))
    rank
  })
  
  # Define LaTeX-style expressions for x-axis labels
  group_labels <- c(
    "center_03" = expression(C ~ rho == 0.3),
    "center_099" = expression(C ~ rho == 0.99),
    "noncenter_03" = expression(NC ~ rho == 0.3),
    "noncenter_099" = expression(NC ~ rho == 0.99)
  )
  
  # Convert list to a tidy data frame
  df <- u_rank_scores %>%
    enframe(name = "Group", value = "Values") %>%
    unnest(Values)
  df$Group <- factor(df$Group, levels = names(group_labels))  # Ensure correct ordering
  
  # Compute the 2.5%, 50% (median), and 97.5% quantiles for each group
  df_summary <- df %>%
    group_by(Group) %>%
    summarise(
      lowest = min(Values),
      lower = quantile(Values, 0.025),   # 2.5th percentile
      middle = quantile(Values, 0.50),  # Median (50th percentile)
      upper = quantile(Values, 0.975),   # 97.5th percentile
      highest = max(Values),
      .groups = "drop"
    )
  
  # Create the customized boxplot
  ggplot(df, aes(x = Group, y = Values)) +
    # Use geom_segment() for whiskers
    geom_segment(data = df_summary, aes(x = Group, xend = Group, y = lowest, yend = highest), color = "black") +
    # Use geom_crossbar() to create the box
    geom_crossbar(data = df_summary, aes(x = Group, ymin = lower, y = middle, ymax = upper), fill = "white", color = "black") +
    # Add a horizontal reference line at y = 1/6
    geom_hline(yintercept = 1/6, linetype = "dashed", color = "red") +
    # Clean theme
    theme_minimal() +
    labs(x = NULL, y = "u-rank within simulation run", title = NULL) +
    scale_x_discrete(labels = group_labels) +  # Apply LaTeX-style labels
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          legend.position = "none")  # Remove legend
  
  ggsave(filename = '../../Figures/03172025_urank_boxplot.png', height = 5, width = 7)
}



#
#### 4/18/2025: Getting real non-centered HMC diagnostics and plots ####
setwd(root_results)
files <- grep('04_09', dir('real_data', full.names = T), value = T)

## Inspecting the convergence diagnostics
{
  # Initialize a data frame to hold diagnostics
  diag_summary <- tibble(
    file = character(),
    dataset = character(),
    softmax = logical(),
    preprocess = logical(),
    alpha = logical(),
    effects = character(),
    n_divergent = integer(),
    max_treedepth_hit = logical(),
    bfmi_low = numeric(),
    phi_noncentered = logical()
  )
  
  for (f in files) {
    load(f)
    fit <- res$sim_list$stan_fit
    
    if (!is.null(fit)) {
      sampler_params <- tryCatch(get_sampler_params(fit, inc_warmup = FALSE), error = function(e) NULL)
      
      if (!is.null(sampler_params) && length(sampler_params) > 0) {
        n_divergent <- sum(sapply(sampler_params, function(chain) sum(chain[, "divergent__"])))
        #max_treedepth_hit <- any(sapply(sampler_params, function(chain) any(chain[, "treedepth__"] >= 10)))
        mean_treedepth_hit <- mean(sapply(sampler_params, function(chain) mean(chain[, "treedepth__"] >= 10)))
        bfmi_by_chain <- sapply(sampler_params, function(chain) {
          energy <- chain[, "energy__"]
          numer <- sum(diff(energy)^2) / (length(energy) - 1)
          denom <- var(energy)
          bfmi <- numer / denom
          bfmi
        })
        
        bfmi_low <- mean(bfmi_by_chain < 0.3)
        
        diag_summary <- add_row(diag_summary,
                                file = basename(f),
                                dataset = ifelse(params$dataset == 'all', 'fullpop', 'AIAN'), 
                                softmax = ifelse(params$use_softmax, T, F),
                                preprocess = ifelse(params$preprocess_scale, T, F),
                                alpha = ifelse(params$alpha_variance_prior == -1, F,T),
                                effects = params$fixed_effects,
                                n_divergent = n_divergent,
                                max_treedepth_hit = mean_treedepth_hit,
                                bfmi_low = bfmi_low,
                                phi_noncentered = (params$phi_noncentered == 1)
        )
      } else {
        diag_summary <- add_row(diag_summary,
                                file = basename(f),
                                dataset = NA,
                                softmax = NA,
                                preprocess = NA,
                                alpha = NA,
                                effects = NA,
                                n_divergent = NA_integer_,
                                max_treedepth_hit = NA,
                                bfmi_low = NA,
                                phi_noncentere = NA
        )
        warning(paste("No sampler params in", f, "- likely due to 0 samples"))
      }
    } else {
      warning(paste("Failed to load stan_fit in", f))
    }
  }
  
  # save(diag_summary, file = 'processed_results/real_data_softmax_noncentering.RData')
  
  tt <- diag_summary[,c(2, 10, 7, 8, 9)]
}

## inspecting an example correlation of phi and tau2.
{
  load("real_data/real_data_fit_full_softmax_pepdensity_alpha0001_noncentered_ID94010_2025_04_09.RData")
  fit1 <- res$sim_list$stan_fit
  
  samples <- rstan::extract(fit1)
  phi_1 <- samples$phi_estimated[, , 1]
  tau2_1 <- samples$tau2[, 1] 
  correlations1 <- apply(abs(phi_1), 2, function(x) cor(x, tau2_1))
  plot(density(correlations1), xlim = c(-0.5,1))
  
  load("real_data/real_data_fit_full_softmax_pepdensity_alpha0001_centered_ID60768_2025_04_09.RData")
  fit2 <- res$sim_list$stan_fit
  
  samples <- rstan::extract(fit2)
  phi_1 <- samples$phi_estimated[, , 1]
  tau2_1 <- samples$tau2[, 1] 
  correlations2 <- apply(abs(phi_1), 2, function(x) cor(x, tau2_1))
  lines(density(correlations2), col = 'red')
}

## inspecting parameter correlation plots 
{
  #load("real_data/real_data_fit_full_softmax_pepdensity_alpha0001_noncentered_ID94010_2025_04_09.RData")
  #load("real_data/real_data_fit_full_softmax_pepdensity_alpha0001_centered_ID60768_2025_04_09.RData")
  #load("real_data/real_data_fit_aian_softmax_pepdensity_alpha0001_noncentered_ID37322_2025_04_09.RData")
  load('real_data/real_data_fit_aian_softmax_pepdensity_alpha0001_centered_ID74511_2025_04_09.RData')
  fit <- res$sim_list$stan_fit
  check_hmc_diagnostics(fit)
  # ok so something is very, very off here. Why? How could we have 1000 divergences and no E-BFMI?
  
  # Cross plots of parameters
  target_pars <- c('theta','alpha[1]','alpha[2]','alpha[3]')
  plot_divergent_pairs(fit, target_pars)
  plot_divergent_pairs(fit, c("tau2_estimated[1]", "tau2_estimated[2]", "tau2_estimated[3]", 'rho[1]', 'rho[2]', 'rho[3]'))
}

## Plotting the results for example fits:
{
  load('real_data/real_data_fit_full_softmax_pepdensity_alpha0001_centered_ID60768_2025_04_09.RData')
  p1 <- plot_real_results(data_list = res$sim_list$data_list,
                          stan_fit = res$sim_list$stan_fit,
                          stan_summary = res$sim_list$stan_summary$summary,
                          models = params$models,
                          CV_pred = NULL,
                          alpha_estimates = T,
                          ESS = T, rho_estimates = T, tau2_estimates = T, 
                          sigma2_estimates = F, theta_estimates = T, phi_estimates = F,
                          pairwise_phi_estimates = T, y_estimates = F, metrics_values = T, beta_estimates = T)
  
  ggsave(p1, file = '../Figures/04182025_real_data_fit_full_centered_softmax.png', height = 12, width = 7)
  
  load('real_data/real_data_fit_full_softmax_pepdensity_alpha0001_noncentered_ID94010_2025_04_09.RData')
  p1 <- plot_real_results(data_list = res$sim_list$data_list,
                          stan_fit = res$sim_list$stan_fit,
                          stan_summary = res$sim_list$stan_summary$summary,
                          models = params$models,
                          CV_pred = NULL,
                          alpha_estimates = T,
                          ESS = T, rho_estimates = T, tau2_estimates = T, 
                          sigma2_estimates = F, theta_estimates = T, phi_estimates = F,
                          pairwise_phi_estimates = T, y_estimates = F, metrics_values = T, beta_estimates = T)
  
  ggsave(p1, file = '../Figures/04182025_real_data_fit_full_noncentered_softmax.png', height = 12, width = 7)
}


#
#### 4/7/2025: Getting real data HMC diagnostics ####
setwd(root_results)

# get the results files.
files <- grep('04_04', dir('real_data', full.names = T), value = T)

# Initialize a data frame to hold diagnostics
diag_summary <- tibble(
  file = character(),
  dataset = character(),
  softmax = logical(),
  preprocess = logical(),
  alpha = logical(),
  effects = character(),
  n_divergent = integer(),
  max_treedepth_hit = logical(),
  bfmi_low = logical()
)

for (f in files) {
  load(f)
  fit <- res$sim_list$stan_fit
  
  if (!is.null(fit)) {
    sampler_params <- tryCatch(get_sampler_params(fit, inc_warmup = FALSE), error = function(e) NULL)
    
    if (!is.null(sampler_params) && length(sampler_params) > 0) {
      n_divergent <- sum(sapply(sampler_params, function(chain) sum(chain[, "divergent__"])))
      #max_treedepth_hit <- any(sapply(sampler_params, function(chain) any(chain[, "treedepth__"] >= 10)))
      mean_treedepth_hit <- mean(sapply(sampler_params, function(chain) mean(chain[, "treedepth__"] >= 10)))
      bfmi_low <- any(sapply(sampler_params, function(chain) {
        mean_energy <- mean(chain[, "energy__"])
        var_energy <- var(chain[, "energy__"])
        bfmi <- mean_energy^2 / var_energy
        bfmi < 0.3
      }))
      
      diag_summary <- add_row(diag_summary,
                              file = basename(f),
                              dataset = ifelse(params$dataset == 'all', 'fullpop', 'AIAN'), 
                              softmax = ifelse(params$use_softmax, T, F),
                              preprocess = ifelse(params$preprocess_scale, T, F),
                              alpha = ifelse(params$alpha_variance_prior == -1, F,T),
                              effects = params$fixed_effects,
                              n_divergent = n_divergent,
                              max_treedepth_hit = mean_treedepth_hit,
                              bfmi_low = bfmi_low
      )
    } else {
      diag_summary <- add_row(diag_summary,
                              file = basename(f),
                              dataset = NA,
                              softmax = NA,
                              preprocess = NA,
                              alpha = NA,
                              effects = NA,
                              n_divergent = NA_integer_,
                              max_treedepth_hit = NA,
                              bfmi_low = NA
      )
      warning(paste("No sampler params in", f, "- likely due to 0 samples"))
    }
  } else {
    warning(paste("Failed to load stan_fit in", f))
  }
}

## Inspecting the fit for file 3
load(files[3])
fit <- res$sim_list$stan_fit

# Cross plots of parameters
target_pars <- c('theta','alpha[1]','alpha[2]','alpha[3]')
plot_divergent_pairs(fit, target_pars)
plot_divergent_pairs(fit, c("tau2_estimated[1]", "tau2_estimated[2]", "tau2_estimated[3]", 'rho[1]', 'rho[2]', 'rho[3]'))


# 

#
#### 4/4/2025: Exploring AIAN real data ####
load('data/census_ACS_PEP_WP_AIAN_wDensity_11152024.RData')

{library(ggplot2)
  library(patchwork)
  
  # PEP log-log
  pep_log <- ggplot(df |> filter(census > 0, pep > 0), aes(x = census, y = pep)) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = "Census (log scale)", y = "PEP (log scale)", title = "Census vs PEP") +
    theme_minimal()
  
  # ACS log-log
  acs_log <- ggplot(df |> filter(census > 0, acs > 0), aes(x = census, y = acs)) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = "Census (log scale)", y = "ACS (log scale)", title = "Census vs ACS") +
    theme_minimal()
  
  # WP log-log
  wp_log <- ggplot(df |> filter(census > 0, wp > 0), aes(x = census, y = wp)) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = "Census (log scale)", y = "WP (log scale)", title = "Census vs WP") +
    theme_minimal()
  
  # PEP zoomed (linear)
  pep_zoom <- ggplot(df, aes(x = census, y = pep)) +
    geom_point() +
    coord_cartesian(xlim = c(0, 10), ylim = c(0, 10)) +
    labs(x = "Census", y = "PEP", title = "Census vs PEP (Zoomed In)") +
    theme_minimal()
  
  # ACS zoomed
  acs_zoom <- ggplot(df, aes(x = census, y = acs)) +
    geom_point() +
    coord_cartesian(xlim = c(0, 10), ylim = c(0, 10)) +
    labs(x = "Census", y = "ACS", title = "Census vs ACS (Zoomed In)") +
    theme_minimal()
  
  # WP zoomed
  wp_zoom <- ggplot(df, aes(x = census, y = wp)) +
    geom_point() +
    coord_cartesian(xlim = c(0, 10), ylim = c(0, 10)) +
    labs(x = "Census", y = "WP", title = "Census vs WP (Zoomed In)") +
    theme_minimal()
  
  # Combine plots into a 2x3 grid
  (pep_log | acs_log | wp_log) /
    (pep_zoom | acs_zoom | wp_zoom)
}


#
#### 4/4/2025: Getting HMC diagnostics from previous real data runs ####
setwd(root_results)

# get the results files.
files <- grep('11_12|11_18|11_19', dir('real_data', full.names = T), value = T)

# Initialize a data frame to hold diagnostics
diag_summary <- tibble(
  file = character(),
  dataset = character(),
  softmax = logical(),
  preprocess = logical(),
  alpha = logical(),
  effects = character(),
  n_divergent = integer(),
  max_treedepth_hit = logical(),
  bfmi_low = logical()
)

for (f in files) {
  load(f)
  fit <- res$sim_list$stan_fit
  
  if (!is.null(fit)) {
    sampler_params <- tryCatch(get_sampler_params(fit, inc_warmup = FALSE), error = function(e) NULL)
    
    if (!is.null(sampler_params) && length(sampler_params) > 0) {
      n_divergent <- sum(sapply(sampler_params, function(chain) sum(chain[, "divergent__"])))
      max_treedepth_hit <- any(sapply(sampler_params, function(chain) any(chain[, "treedepth__"] >= 10)))
      bfmi_low <- any(sapply(sampler_params, function(chain) {
        mean_energy <- mean(chain[, "energy__"])
        var_energy <- var(chain[, "energy__"])
        bfmi <- mean_energy^2 / var_energy
        bfmi < 0.3
      }))
      
      diag_summary <- add_row(diag_summary,
                              file = basename(f),
                              dataset = ifelse(params$dataset == 'all', 'fullpop', 'AIAN'), 
                              softmax = ifelse(params$use_softmax, T, F),
                              preprocess = ifelse(params$preprocess_scale, T, F),
                              alpha = ifelse(params$alpha_variance_prior == -1, F,T),
                              effects = params$fixed_effects,
                              n_divergent = n_divergent,
                              max_treedepth_hit = max_treedepth_hit,
                              bfmi_low = bfmi_low
      )
    } else {
      diag_summary <- add_row(diag_summary,
                              file = basename(f),
                              dataset = NA,
                              softmax = NA,
                              preprocess = NA,
                              alpha = NA,
                              effects = NA,
                              n_divergent = NA_integer_,
                              max_treedepth_hit = NA,
                              bfmi_low = NA
      )
      warning(paste("No sampler params in", f, "- likely due to 0 samples"))
    }
  } else {
    warning(paste("Failed to load stan_fit in", f))
  }
}

#
#### 4/4/2025: Getting simulation run metrics AND HMC diagnostics ####
setwd(root_results)
setwd('simulated_results/')

res_gamma <- generate_metrics_list(folder = 'simulation_rho_099_theta_100_softmax_negbin_3models_CV10_ID29995_2025_04_01', hmc_diag = T)

res_invqui <- generate_metrics_list(folder = 'simulation_rho_099_theta_100_softmax_negbin_3models_CV10_invchi_2025_04_01', hmc_diag = T)

sapply(res_gamma$metrics_list, function(xx) xx$bfmi_low) %>% sum()

sapply(res_invqui$metrics_list, function(xx) xx$bfmi_low) %>% sum()
# ok so not much of the three hmc diagnostic parameters. 

## Look at parameter correlations for three runs. Start with SM rho = 0.99 and gamma prior 
# Find simulation runs with low, medium, and high u-rank.
{
  results <- res_gamma
  u_rank_scores <- colMeans(sapply(results$metrics_list, function(yy) yy[['rank_equal']]))
  
  which.min(u_rank_scores)  # 125
  which(u_rank_scores == median(u_rank_scores)) # 19
  which.max(u_rank_scores) # 37
  
  setwd('simulation_rho_099_theta_100_softmax_negbin_3models_CV10_ID29995_2025_04_01')
  load('sim_results_1.RData')
  good <- res_lst[[37]]
  mid <- res_lst[[19]]  
  
  load('sim_results_4.RData')
  bad <- res_lst[[5]]
  
  fit_good <- good$sim_list[[1]]$stan_fit
  fit_mid <- mid$sim_list[[1]]$stan_fit
  fit_bad <- bad$sim_list[[1]]$stan_fit
}


target_pars <- c("theta",
                 "tau2_estimated[1]", "tau2_estimated[2]", "tau2_estimated[3]", 'rho[1]', 'rho[2]', 'rho[3]')

plot_divergent_pairs(fit_good, target_pars)
plot_divergent_pairs(fit_mid, target_pars)
plot_divergent_pairs(fit_bad, target_pars)

## Now for invchi
{
  results <- res_invqui
  u_rank_scores <- colMeans(sapply(results$metrics_list, function(yy) yy[['rank_equal']]))
  
  which.min(u_rank_scores)  # 125
  which(u_rank_scores == median(u_rank_scores)) # 20
  which.max(u_rank_scores) # 37
  
  setwd('simulation_rho_099_theta_100_softmax_negbin_3models_CV10_invchi_2025_04_01')
  load('sim_results_1.RData')
  good <- res_lst[[37]]
  mid <- res_lst[[20]]  
  
  load('sim_results_4.RData')
  bad <- res_lst[[5]]
  
  fit_good <- good$sim_list[[1]]$stan_fit
  fit_mid <- mid$sim_list[[1]]$stan_fit
  fit_bad <- bad$sim_list[[1]]$stan_fit
}

plot_divergent_pairs(fit_good, target_pars)
plot_divergent_pairs(fit_mid, target_pars)
plot_divergent_pairs(fit_bad, target_pars)



##
#### 3/28/2025: Debugging AIAN SM alpha density results ####
library(bayesplot)
library(posterior)
setwd(root_results)
setwd('real_data')
#
load('real_data_fit_aian_softmax_alpha_cv10_pepfulldensity_ID17611_2025_01_15.RData')
fit_AIAN <- res$sim_list$stan_fit
load('real_data_fit_softmax_alpha_cv10_density_ID77695_2025_01_15.RData')
fit_full <- res$sim_list$stan_fit
load('real_data_fit_softmax_alpha_density_invchi_ID90998_2025_04_03.RData')
fit_full_invchi <- res$sim_list$stan_fit

# posterior::variables(as_draws_array(fit_AIAN))
pars_vec <- c('theta',"rho[1]", "rho[2]", "rho[3]",
              "tau2[1]", "tau2[2]", "tau2[3]")

## look at the parameter pairings to see if there are strong correlations between parameters.
AIAN_array <- as.array(fit_AIAN)
AIAN_thinned <- AIAN_array[seq(1, dim(AIAN_array)[1], by = 5), , ]
mcmc_parcoord(AIAN_thinned, pars = pars_vec, transformations = function(x) {(max(x) - x)/(max(x) - min(x))})
plot_param_correlation(fit_AIAN, pars = c('rho_estimated','tau2_estimated','theta', 'beta'))

full_array <- as.array(fit_full)
full_thinned <- full_array[seq(1, dim(full_array)[1], by = 5), , ]
mcmc_parcoord(full_thinned, pars = pars_vec, transformations = function(x) {(max(x) - x)/(max(x) - min(x))})
plot_param_correlation(fit_full, pars = c('rho_estimated','tau2_estimated','theta', 'beta'))

full_array_invchi <- as.array(fit_full_invchi)
full_thinned_invchi <- full_array_invchi[seq(1, dim(full_array_invchi)[1], by = 5), , ]
mcmc_parcoord(full_thinned_invchi, pars = pars_vec, transformations = function(x) {(max(x) - x)/(max(x) - min(x))})
plot_param_correlation(fit_full_invchi, pars = c('alpha','theta', 'beta'))

## Look at the diagnostics of the model fit.
check_hmc_diagnostics(fit_full)
check_hmc_diagnostics(fit_AIAN)
check_hmc_diagnostics(fit_full_invchi)

pars_vec <- c('theta', 'alpha[1]', 'alpha[2]', 'alpha[3]')
plot_divergent_pairs(fit_full, pars_vec)
plot_divergent_pairs(fit_AIAN, pars_vec)
plot_divergent_pairs(fit_full_invchi, pars_vec)

# Now for divergent values.
draws_df <- as_draws_df(fit_AIAN) 

# Extract sampler parameters per chain
sampler_params <- get_sampler_params(fit_AIAN, inc_warmup = FALSE)

# Combine all chains into one long vector
divergent_vec <- unlist(lapply(sampler_params, function(x) x[, "divergent__"]))

draws_df$divergent <- as.logical(divergent_vec)

# Step 3: Choose your 10 parameters
target_pars <- c("theta",
                 "tau2_estimated[1]", "tau2_estimated[2]", "tau2_estimated[3]",
                 "rho_estimated[1]", "rho_estimated[2]", "rho_estimated[3]",
                 "alpha[1]", "alpha[2]", "alpha[3]")
# target_pars <- c("theta", 
#                  "tau2_estimated[1]", "tau2_estimated[2]", "tau2_estimated[3]")

# Step 4: Subset and rename for plotting
plot_df <- draws_df[, c(target_pars, "divergent")]

# Optional: rename columns to make them more readable
colnames(plot_df) <- gsub("_estimated", "", colnames(plot_df))
colnames(plot_df) <- gsub("\\[", "_", gsub("\\]", "", colnames(plot_df)))

# Step 5: Plot
ggpairs(plot_df, 
        columns = 1:length(target_pars), 
        aes(color = divergent, alpha = 0.3),
        upper = list(continuous = wrap("points", size = 0.5)),
        lower = list(continuous = wrap("points", size = 0.5))) +
  scale_color_manual(values = c("black", "green"))

#shinystan::launch_shinystan(fit_AIAN)

#
#### 3/25/2025: Inspecting why higher instability in higher rho simulation ####
setwd(root_results)
setwd('simulated_results/')
files <- grep('2025_03_13', dir(), value = T)
# I want to look at this hypothesis: The low u-rank results correspond to poor rho estimation. Let's start with SM. Meh. Let's do it all!

# for each file, I want to pull out the median u-rank in each simulation and the median rho estimate in each simulation.
results_list <- lapply(files, function(f){
  generate_metrics_list(f)})

tt <- results_list[[1]]$metrics_list

# pull out the u-rank of each file
u_rank_scores <- lapply(results_list, function(xx){
  rank <- colMeans(sapply(xx$metrics_list, function(yy) yy[['rank_equal']]))
  rank
})

rho_medians <- lapply(results_list, function(xx){
  rank <- sapply(xx$metrics_list, function(yy) yy[['median_rhoX']])
  rank
})

par(mfrow = c(2,2))
for(i in 1:4){
  plot(rho_medians[[i]], u_rank_scores[[i]])
  print(files[i])
  print(cor(rho_medians[[i]], u_rank_scores[[i]]))
}
# so there doesnt appear to be a strong association, does there? Maybe. There is a correlation of 0.4 between u-rank and rho median value for SM rho = 0.99, but virtually none for DE rho = 0.99. Idk what to make of that.

#
#### 3/17/2025: Make simulation results figures and tables ####
setwd(root_results)
setwd('simulated_results/')
files <- grep('2025_03_13', dir(), value = T)

results_list <- lapply(files, function(f){
  generate_metrics_list(f)})

warning('Hardcoding names of results - make sure they match')
names(results_list) <- c('DE_rho03', 'SM_rho03', 'DE_rho099', 'SM_rho099')

### Getting u-rank scores 
{
  # pull out the u-rank of each file
  u_rank_scores <- lapply(results_list, function(xx){
    rank <- colMeans(sapply(xx$metrics_list, function(yy) yy[['rank_equal']]))
    rank
  })
  
  # Define LaTeX-style expressions for x-axis labels
  group_labels <- c(
    "SM_rho03" = expression(SM ~ rho == 0.3),
    "SM_rho099" = expression(SM ~ rho == 0.99),
    "DE_rho03" = expression(DE ~ rho == 0.3),
    "DE_rho099" = expression(DE ~ rho == 0.99)
  )
  
  # Convert list to a tidy data frame
  df <- u_rank_scores %>%
    enframe(name = "Group", value = "Values") %>%
    unnest(Values)
  df$Group <- factor(df$Group, levels = names(group_labels))  # Ensure correct ordering
  
  # Compute the 2.5%, 50% (median), and 97.5% quantiles for each group
  df_summary <- df %>%
    group_by(Group) %>%
    summarise(
      lowest = min(Values),
      lower = quantile(Values, 0.025),   # 2.5th percentile
      middle = quantile(Values, 0.50),  # Median (50th percentile)
      upper = quantile(Values, 0.975),   # 97.5th percentile
      highest = max(Values),
      .groups = "drop"
    )
  
  # Create the customized boxplot
  ggplot(df, aes(x = Group, y = Values)) +
    # Use geom_segment() for whiskers
    geom_segment(data = df_summary, aes(x = Group, xend = Group, y = lowest, yend = highest), color = "black") +
    # Use geom_crossbar() to create the box
    geom_crossbar(data = df_summary, aes(x = Group, ymin = lower, y = middle, ymax = upper), fill = "white", color = "black") +
    # Add a horizontal reference line at y = 1/6
    geom_hline(yintercept = 1/6, linetype = "dashed", color = "red") +
    # Clean theme
    theme_minimal() +
    labs(x = NULL, y = "u-rank within simulation run", title = NULL) +
    scale_x_discrete(labels = group_labels) +  # Apply LaTeX-style labels
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          legend.position = "none")  # Remove legend
  
  ggsave(filename = '../../Figures/03172025_urank_boxplot.png', height = 5, width = 7)
}

### Making the table metrics 
{
  # Function to compute median and 95% quantiles
  summary_stats <- function(vec) {
    q_values <- quantile(vec, probs = c(0.025, 0.975), na.rm = TRUE)
    
    return(c(
      median = median(vec, na.rm = TRUE),
      Q2.5 = unname(q_values[1]),  # Ensure it doesn't inherit unwanted names
      Q97.5 = unname(q_values[2])
    ))
  }
  
  all_sim_results <- lapply(results_list, function(xx){
    metrics_list <- xx$metrics_list
    
    # Extract numeric metrics
    MAPE_train <- sapply(metrics_list, function(x) x$MAPE_train)
    MAPE_CV <- sapply(metrics_list, function(x) x$MAPE_CV)
    MAE_train <- sapply(metrics_list, function(x) x$MAE_train)
    MAE_CV <- sapply(metrics_list, function(x) x$MAE_CV)
    median_int_width_train <- sapply(metrics_list, function(x) x$median_int_width_train)
    median_int_width_CV <- sapply(metrics_list, function(x) x$median_int_width_CV)
    
    # Compute proportion of TRUE values for Boolean vectors
    CP_95_train <- sapply(metrics_list, function(x) mean(x$CP_95_train, na.rm = TRUE))
    CP_90_train <- sapply(metrics_list, function(x) mean(x$CP_90_train, na.rm = TRUE))
    CP_95_CV <- sapply(metrics_list, function(x) mean(x$CP_95_CV, na.rm = TRUE))
    CP_90_CV <- sapply(metrics_list, function(x) mean(x$CP_90_CV, na.rm = TRUE))
    
    # Compute summary statistics for all extracted values
    final_results <- list(
      MAPE_train = summary_stats(MAPE_train),
      MAPE_CV = summary_stats(MAPE_CV),
      MAE_train = summary_stats(MAE_train),
      MAE_CV = summary_stats(MAE_CV),
      median_int_width_train = summary_stats(median_int_width_train),
      median_int_width_CV = summary_stats(median_int_width_CV),
      CP_95_train = summary_stats(CP_95_train),
      CP_90_train = summary_stats(CP_90_train),
      CP_95_CV = summary_stats(CP_95_CV),
      CP_90_CV = summary_stats(CP_90_CV)
    )
    
    final_results
  })
  
  generate_latex_values <- function(data_list, coverage = 95, digits = 2) {
    # Validate input
    if (!coverage %in% c(90, 95)) stop("Coverage must be 90 or 95.")
    
    # Define custom order for models: SM should come before DE
    model_order <- c("SM", "DE")
    
    # Define custom order for dataset type: Train should come before CV
    dataset_order <- c("train", "CV")
    
    # Initialize a list to store row data before sorting
    rows_list <- list()
    
    # Check for names.
    if(is.null(names(data_list))){
      stop('data list inputted needs to be named.')
    }
    
    # Iterate through each simulation run
    for (sim_name in names(data_list)) {
      sim_data <- data_list[[sim_name]]
      
      # Extract the prefix (e.g., "DE" or "SM") and rho value
      split_name <- strsplit(sim_name, "_rho")[[1]]
      model_name <- split_name[1]  # "DE" or "SM"
      rho_value <- ifelse(split_name[2] == '03', 0.3,
                          ifelse(split_name[2] == '099', 0.99, NA))  # Convert rho to numeric for sorting
      
      # Format the first column for LaTeX: "DE, $\rho=0.3$"
      latex_name <- paste0(model_name, ", $\\rho=", rho_value, "$")
      
      # Iterate over Train and CV
      for (type in c("train", "CV")) {
        # Determine variable names
        MAPE_var <- paste0("MAPE_", type)
        MAE_var <- paste0("MAE_", type)
        CP_var <- paste0("CP_", coverage, "_", type)
        width_var <- paste0("median_int_width_", type)
        
        # Extract and round values
        MAPE <- round(sim_data[[MAPE_var]], digits)
        MAE <- round(sim_data[[MAE_var]], digits)
        CP <- round(sim_data[[CP_var]], digits)
        width <- round(sim_data[[width_var]], digits)
        
        # Format as "median (Q2.5, Q97.5)"
        MAPE_str <- paste0(MAPE["median"], " (", MAPE["Q2.5"], ", ", MAPE["Q97.5"], ")")
        MAE_str <- paste0(MAE["median"], " (", MAE["Q2.5"], ", ", MAE["Q97.5"], ")")
        CP_str <- paste0(CP["median"], " (", CP["Q2.5"], ", ", CP["Q97.5"], ")")
        width_str <- paste0(width["median"], " (", width["Q2.5"], ", ", width["Q97.5"], ")")
        
        # Store row data in a list for sorting
        rows_list <- append(rows_list, list(
          data.frame(model_name = model_name, rho = rho_value, type = type, 
                     row_text = paste(latex_name, type, MAPE_str, MAE_str, CP_str, width_str, sep = " & "))
        ))
      }
    }
    
    # Convert list to data frame
    rows_df <- do.call(rbind, rows_list)
    
    # Convert factors to enforce sorting order
    rows_df$model_name <- factor(rows_df$model_name, levels = model_order)  # SM first, then DE
    rows_df$type <- factor(rows_df$type, levels = dataset_order)  # Train first, then CV
    
    # Order rows: SM first, then ascending rho, then Train before CV
    rows_df <- rows_df[order(rows_df$model_name, rows_df$rho, rows_df$type), ]
    
    # Extract ordered LaTeX rows
    output_lines <- rows_df$row_text
    
    return(output_lines)
  }
  
  latex_rows <- generate_latex_values(data_list = all_sim_results, coverage = 95, digits = 2)
  cat(paste(latex_rows, collapse = " \\\\\n"))
}
#
#### 1/21/2025: Make chloropleth plots of the results ####
setwd(root_results)
setwd('real_data/')
load('real_data_fit_directest_cv10_interceptonly_ID81515_2025_01_15.RData')

# stan_out <- extract(res$sim_list$stan_fit)

data_list = res$sim_list$data_list
N = nrow(data_list$data)
models <- params$models
stan_summary <- res$sim_list$stan_summary$summary
u_est <- as.data.frame(matrix(0, nrow = N, ncol = length(models)))
colnames(u_est) <- models
for(i in 1:length(models)){
  ind = grep(sprintf('^u\\[[0-9]*,%s\\]', i), rownames(stan_summary))
  u_est[,i] <- stan_summary[ind,'50%']
}
u_est$index = 1:N

# Bring back in county fips codes.
u_est$GEOID <- data_list$data$GEOID

# get county shapefiles.
counties <- tigris::counties(cb = TRUE, year = 2020, class = "sf")

# Load state boundaries.
states <- tigris::states(cb = TRUE, year = 2020, class = "sf")

merged_df <- merge(counties, u_est, by = 'GEOID')
merged_df$outcome <- merged_df$pep

p1 <- ggplot(data = merged_df) +
  geom_sf(aes(fill = outcome), color = NA) +
  geom_sf(data = states, fill = NA, color = "black", size = 0.5) + 
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0.97, # Centering the scale around 1
    limits = c(0.83, 1.07), # Adjusting to your data range
    na.value = "white"
  ) +
  theme_minimal() +
  coord_sf(xlim = c(-130, -65), ylim = c(24, 50)) +
  labs(fill = "PEP Weight") +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

ggsave(p1, file = '../../Figures/01152025_pep_estimates_chloropleth_directest_cv10_interceptonly.png',
       width = 10, height = 5)

## Other code for density plotting
if(F){
  # convert estimates to long
  u_est_long <- tidyr::gather(u_est, key = model, value = u_median_est, -index)
  
  # plot 'em
  p_u <- ggplot(u_est_long, aes(x = u_median_est)) + 
    geom_density() + 
    facet_wrap(~model, scales = 'free') + 
    xlab('estimate') + 
    ylab('density') + 
    theme(axis.text.y = element_blank()) +
    ggtitle('u estimates')
}


#
#### 1/21/2025: Make weight density plots ####
setwd(root_results)
setwd('real_data/')
#load('real_data_fit_directest_cv10_interceptonly_ID81515_2025_01_15.RData')
load('real_data_fit_softmax_alpha_cv10_density_ID77695_2025_01_15')

p1 <- plot_real_results(data_list = res$sim_list$data_list,
                        stan_fit = res$sim_list$stan_fit,
                        stan_summary = res$sim_list$stan_summary$summary,
                        models = params$models,
                        CV_pred = res$sim_list$CV_pred, rhats = F,
                        alpha_estimates = F,
                        ESS = F, rho_estimates = F, tau2_estimates = F, 
                        sigma2_estimates = F, theta_estimates = F, phi_estimates = F,
                        pairwise_phi_estimates = F, y_estimates = F, metrics_values = F, beta_estimates = F)
# inspect plot. How is it?
# ggsave(p1, file = '../../Figures/01152025_u_estimates_directest_cv10_interceptonly.png', width = 6, height = 3)
ggsave(p1, file = '../../Figures/01152025_u_estimates_softmax_alpha_density_cv10.png', width = 6, height = 3)

#
#### 1/20/2025: All the recent results! ####

setwd(root_results)
setwd('real_data/')
recent_files(2)
files <- grep('2025_01_15', dir(), value = T)
ind_aian <- grep('aian', files)
files_full <- files[-ind_aian]
files_aian <- files[ind_aian]

### Full model
plots <- list()
metrics <- list()
file_names <- c()

for (i in seq_along(files_full)) {
  print(i)
  print(files_full[i])
  tryCatch({
    load(files_full[i])
    stan_fit <-  res$sim_list$stan_fit
    if(dim(stan_fit)[1]*dim(stan_fit)[2] > 10000){
      print('skipping because dim(stan_fit) = ')
      print(dim(stan_fit))
      next
    }
    res <- just_metrics(
      data_list = res$sim_list$data_list,
      stan_fit = res$sim_list$stan_fit,
      stan_summary = res$sim_list$stan_summary$summary,
      models = params$models,
      CV_pred = res$sim_list$CV_pred
    )
    plots <- append(plots, list(res$plot))
    metrics <- append(metrics, list(res$metrics))
    file_names <- c(file_names, files_full[i])
  }, error = function(e) {
    cat("Error in file:", files_full[i], "\n")
    cat("Error message:", e$message, "\n")
  })
}
# real_data_fit_softmax_alpha_cv10_density_ID77695_2025_01_15.RData had 17/1000 divergence. IGNORING, though this could be an issue.
# real_data_fit_softmax_alpha_density_20172018_ID59308_2025_01_15.RData had 176/1000 divergence!!
# real_data_fit_softmax_alpha_density_2018only_ID54737_2025_01_15.RData had 63/1000 divergence.

names(metrics) <- file_names
CV_MAE <- sapply(metrics, function(tmp){ tmp[2,'MAE']})
CV_MAPE <- sapply(metrics, function(tmp){ tmp[2,'MAPE']})
MAE_ord <- order(CV_MAE)
metrics[MAE_ord]
# save(plots, metrics, files_full, file = '../01202025_full_population_metrics_CVand5models.RData')
load('../01202025_full_population_metrics_CVand5models.RData')


### AIAN model
plots_aian <- list()
metrics_aian <- list()
file_names_aian <- c()

for (i in seq_along(files_aian)) {
  print(i)
  print(files_aian[i])
  tryCatch({
    load(files_aian[i])
    stan_fit <-  res$sim_list$stan_fit
    if(dim(stan_fit)[1]*dim(stan_fit)[2] > 10000){
      print('skipping because dim(stan_fit) = ')
      print(dim(stan_fit))
      next
    }
    res <- just_metrics(
      data_list = res$sim_list$data_list,
      stan_fit = res$sim_list$stan_fit,
      stan_summary = res$sim_list$stan_summary$summary,
      models = params$models,
      CV_pred = res$sim_list$CV_pred
    )
    plots_aian <- append(plots_aian, list(res$plot))
    metrics_aian <- append(metrics_aian, list(res$metrics))
    file_names_aian <- c(file_names_aian, files_aian[i])
  }, error = function(e) {
    cat("Error in file:", files_aian[i], "\n")
    cat("Error message:", e$message, "\n")
  })
}
# real_data_fit_aian_softmax_alpha_pepfulldensity_2018only_ID88033_2025_01_15.RData had 278/1000 divergences!
# real_data_fit_aian_softmax_alpha_pepfulldensity_20182019_ID48754_2025_01_15.RData had 374/1000 divergences!
# real_data_fit_aian_softmax_alpha_pepfulldensity_20172018_ID60604_2025_01_15.RData had 366/1000 divergences!
# real_data_fit_aian_softmax_alpha_cv10_pepfulldensity_ID17611_2025_01_15.RData had 300/1000 divergences!
# the first four had 2-8 divergences.

names(metrics_aian) <- file_names_aian
CV_MAE <- sapply(metrics_aian, function(tmp){ tmp[2,'MAE']})
CV_MAPE <- sapply(metrics_aian, function(tmp){ tmp[2,'MAPE']})
MAE_ord_aian <- order(CV_MAE)
metrics_aian[MAE_ord_aian]
save(plots_aian, metrics_aian, files_full, file = '../01202025_aianmetrics_CVand5models.RData')



#
#### 1/13/2025: Getting metrics for all results - full pop and AIAN ####
setwd(root_results)
setwd('real_data/')

files_full <- dir()[-grep('aian', dir(), ignore.case = T)]

plotz <- list()
metricz <- list()
file_names <- c()

for (i in seq_along(files_full)) {
  print(i)
  print(files_full[i])
  tryCatch({
    load(files_full[i])
    if(params$dataset == 'aian'){next} # didnt originally have this but it seems important.
    stan_fit <-  res$sim_list$stan_fit
    if(dim(stan_fit)[1]*dim(stan_fit)[2] > 10000){
      print('skipping because dim(stan_fit) = ')
      print(dim(stan_fit))
      next
    }
    res <- just_metrics(
      data_list = res$sim_list$data_list,
      stan_fit = res$sim_list$stan_fit,
      stan_summary = res$sim_list$stan_summary$summary,
      models = params$models,
      CV_pred = res$sim_list$CV_pred
    )
    plotz <- append(plotz, list(res$plot))
    metricz <- append(metricz, list(res$metrics))
    file_names <- c(file_names, files_full[i])
  }, error = function(e) {
    cat("Error in file:", files_full[i], "\n")
    cat("Error message:", e$message, "\n")
  })
}
# real_data_fit_softmax_alpha_interceptonly_ID24617_2024_11_12.RData had 2.6% divergence. N bueno
# save the raw results. 
save(plotz, metricz, file_names, file = '../01092025_full_population_metrics.RData')
# save the plot! Extra long, of course.
p1 <- plot_grid(plotlist = plotz, ncol = 1, labels = file_names)

# get order of CV MAE
names(metricz) <- file_names
CV_MAE <- sapply(metricz, function(tmp){ tmp[2,'MAE']})
CV_MAPE <- sapply(metricz, function(tmp){ tmp[2,'MAPE']})
MAE_ord <- order(CV_MAE)

# look at them all.
metricz[MAE_ord]

# now do this for aian
{
files_aian <- dir()[grep('aian', dir(), ignore.case = T)]

plotz_aian <- list()
metrics_aian <- list()
file_names_aian <- c()

for (i in seq_along(files_aian)) {
  print(i)
  print(files_aian[i])
  tryCatch({
    load(files_aian[i])
    res <- just_metrics(
      data_list = res$sim_list$data_list,
      stan_fit = res$sim_list$stan_fit,
      stan_summary = res$sim_list$stan_summary$summary,
      models = params$models,
      CV_pred = res$sim_list$CV_pred
    )
    plotz_aian <- append(plotz_aian, list(res$plot))
    metrics_aian <- append(metrics_aian, list(res$metrics))
    file_names_aian <- c(file_names_aian, files_aian[i])
  }, error = function(e) {
    cat("Error in file:", files_aian[i], "\n")
    cat("Error message:", e$message, "\n")
  })
}
# "real_data_fit_softmax_preprocess_alpha_aian_ID97760_2024_10_28.RData" had 8/10k divergences.

# "real_data_fit_aian_directest_preprocess_interceptonly_5models_ID55316_2025_01_03.RData" had 11/1000 divergences (maybe concerning)

# save it!
save(plotz_aian, metrics_aian, file_names_aian, file = '../01092025_AIAN_metrics.RData')
}

outcome_aian <- data.frame()
for(i in seq_along(metrics_aian)){
  extracted_text <- stringr::str_extract(file_names_aian[i], "(?<=fit_).*?(?=_ID)")
  print(extracted_text)
  tmp_df <- data.frame(model = rep(extracted_text, 2)) %>%
    cbind(metrics_aian[[i]])
  outcome_aian <- rbind(outcome_aian, tmp_df)
}

names(metrics_aian) <- file_names_aian
CV_MAE <- sapply(metrics_aian, function(tmp){ tmp[2,'MAE']})
CV_MAPE <- sapply(metrics_aian, function(tmp){ tmp[2,'MAPE']})
MAE_ord <- order(CV_MAE)

# look at them all.
metrics_aian[MAE_ord]


#

#
#### Results with ACS and PEP PCs - direct est ####
setwd(root_results)
setwd('real_data/')
recent_files(2)

file <- recent_files(1)$file

for(f in file){
  load(f)
  #divergence_check(res)  
  p1 <- plot_real_results(data_list = res$sim_list$data_list,
                          stan_fit = res$sim_list$stan_fit,
                          stan_summary = res$sim_list$stan_summary$summary,
                          models = params$models,
                          CV_pred = res$sim_list$CV_pred,
                          alpha_estimates = F)
  out_name <- sprintf('../../Figures/%s_%s_real_data.png',
                      format(Sys.Date(), "%m%d%Y"), 
                      sub(".*fit_(.*?)_ID.*", "\\1", f))
  print(out_name)
  ggsave(plot = p1, filename = out_name, height = 12, width = 7)
}

#
#### Results with ACS and PEP 2018 and 2019 ####
setwd(root_results)
setwd('real_data/')
recent_files(2)

files <- recent_files(l = 2) %>%
  pull(file)

# save images of results
for(f in files){
  load(f)
  #divergence_check(res)  
  p1 <- plot_real_results(data_list = res$sim_list$data_list,
                          stan_fit = res$sim_list$stan_fit,
                          stan_summary = res$sim_list$stan_summary$summary,
                          models = params$models,
                          CV_pred = res$sim_list$CV_pred,
                          alpha_estimates = F)
  out_name <- sprintf('../../Figures/01052025_%s_real_data.png', sub(".*fit_(.*?)_ID.*", "\\1", f))
  print(out_name)
  #ggsave(plot = p1, filename = out_name, height = 12, width = 7)
}

# doing some error checking
{
load('real_data_fit_aian_directest_preprocess_interceptonly_5models_ID55316_2025_01_03.RData')

#load('real_data_fit_softmax_preprocess_density_ID10231_2024_11_19.RData')

fit <- res$sim_list$stan_fit
stan_diag(fit)

stan_diag(fit, info = 'sample')

stan_diag(fit, info = 'stepsize')

stan_diag(fit, info = 'treedepth')

stan_diag(fit, info = 'divergence')

check_divergences(fit)
}
#
#### 1/3/2025: Theta results part 2 ####
setwd(root_results)
setwd('real_data/')

# pull file names
recent_files(l = 6)
files <- recent_files(l = 3) %>%
  pull(file)

# save images of results
for(f in files){
  load(f)
  p1 <- plot_real_results(data_list = res$sim_list$data_list,
                          stan_fit = res$sim_list$stan_fit,
                          stan_summary = res$sim_list$stan_summary$summary,
                          CV_pred = res$sim_list$CV_pred,
                          alpha_estimates = F)
  out_name <- sprintf('../../Figures/01022025_AIAN_%s_real_data.png', sub(".*fit_(.*?)_ID.*", "\\1", f))
  ggsave(plot = p1, filename = out_name, height = 12, width = 7)
}

# Done! Now inspect the results!

#
#### 1/2/2025: Fixed tau2 results and theta (?) results from last month ####
setwd(root_results)
files <- grep('11_22|11_26', dir('real_data', full.names = T), value = T)

for(f in files){
  load(f)
  p1 <- plot_real_results(data_list = res$sim_list$data_list,
                          stan_fit = res$sim_list$stan_fit,
                          stan_summary = res$sim_list$stan_summary$summary,
                          CV_pred = res$sim_list$CV_pred,
                          alpha_estimates = F)
  out_name <- sprintf('../Figures/01022025_AIAN_%s_real_data.png', sub(".*fit_(.*?)_ID.*", "\\1", f))
  ggsave(plot = p1, filename = out_name, height = 12, width = 7)
}

#
#### 11/21/2024: AIAN get results and make chloropleth maps ####
setwd(root_results)

## Organizing the results
{
files <- grep('aian', dir('real_data', full.names = T), value = T)
# doing file 5 because it's not preprocessed.
load(files[5])
params$preprocess_scale

df <- res$sim_list$data_list$data %>%
  select(GEOID, state, census, acs, pep, wp)

for(f in files){
  load(f)
  outcome_name = sprintf('AIAN_%s_%s_%s_%s',
                      ifelse(params$use_softmax, 'softmax', 'directest'),
                      ifelse(params$preprocess_scale, 'centering', 'NOcentering'),
                      ifelse(params$alpha_variance_prior == -1, 'NOalpha','alpha'),
                      ifelse(!is.null(params$fixed_effects), params$fixed_effects, 'noFE'))
  print(f)
  print(outcome_name)
  print('----------------')
  
  # check that the GEOID is equal.
  if(any(res$sim_list$data_list$data$GEOID != df$GEOID)){
    stop('error - GEOIDs not equal.')
  }
  
  # get the median!
  if(!is.null(res$sim_list$stan_summary$summary)){
    tmp <- summary(res$sim_list$stan_fit)$summary
    out_col <- tmp[grep('y_pred', rownames(tmp)), '50%']
  }else{
    samples <- extract(res$sim_list$stan_fit, pars = "y_pred", permuted = TRUE)$y_pred
    out_col <- apply(samples,2, median)
  }
  df[,outcome_name] <- out_col
}

save(df, file = 'AIAN_medians_11212024.RData')
}

## save names
aian_names <- colnames(df[,7:14])
colnames(df)[7:14] <- paste0('model', 1:8)

## Getting the correlations and means.
cor_matrix <- cor(df[,3:14])

ggcorrplot::ggcorrplot(cor_matrix, lab = TRUE, lab_size = 3, colors = c("blue", "white", "red"))

colMeans(df[,3:14])

## Merging in mortality.
setwd(root_git)
load('data/AmericanIndian_COVID_Deaths_2020.RData')
df2 <- merge(cvd, df)

df2 %>% 
  group_by(state) %>% 
  summarize(n = n(), n_d = sum(AIAN_deaths)) %>% 
  arrange(n)

## Make chloropleth plots.
AIAN_cvd <- df2 %>%
  filter(state %in% c('Arizona','New Mexico','Colorado','Utah','Oklahoma'))

# Create new columns for mortality rates
for (i in 4:15) {
  col_name <- paste0("mortality_rate_", names(AIAN_cvd)[i])
  AIAN_cvd[[col_name]] <- AIAN_cvd$AIAN_deaths / AIAN_cvd[[i]]
}

# Inspect the updated dataset
head(AIAN_cvd)

# load the counties
counties <- tigris::counties(cb = TRUE, year = 2020, class = "sf")

# Get state boundaries for the year 2020
states <- tigris::states(cb = TRUE, year = 2020, class = "sf") %>%
  filter(NAME %in% AIAN_cvd$state)


# Merge the shapefile with the mortality rates data
merged_data <- merge(counties, AIAN_cvd, by = 'GEOID')  # Adjust FIPS column names if different

# Create a list to store the plots
plots <- list()

# Generate plots for each mortality rate column
for (i in 4:15) {
  col_name <- paste0("mortality_rate_", names(AIAN_cvd)[i])
  
  # Create the plot
  p <- ggplot(data = merged_data) +
    geom_sf(aes(fill = !!sym(col_name))) +
    geom_sf(data = states, fill = NA, color = "black") +
    scale_fill_viridis_c() +
    ggtitle(paste("Mortality Rate:", col_name)) +
    theme_minimal()
  
  # Add the plot to the list
  plots[[col_name]] <- p
}

library(patchwork)

# Define global limits for the color scale
range_values <- range(AIAN_cvd[ , paste0("mortality_rate_", names(AIAN_cvd)[4:15])], na.rm = TRUE)

# Update plots with consistent color scale
for (i in 4:15) {
  col_name <- paste0("mortality_rate_", names(AIAN_cvd)[i])
  plots[[col_name]] <- plots[[col_name]] +
    scale_fill_viridis_c(limits = range_values)
}

# Combine all plots into a grid
combined_plot <- wrap_plots(plots, ncol = 3)
print(combined_plot)

# save it
setwd(root_results)
ggsave("../Figures/mortality_chloropleth_11212024.png", combined_plot, width = 15, height = 10)

#
#### Fixed effects with covariates results ####
setwd(root_results)

# get the results files.
files <- grep('11_12|11_18|11_19', dir('real_data', full.names = T), value = T)

# cycle through and save results.
for(f in files){
  load(f)
  plot_name = sprintf('../Figures/%s_%s_%s_%s_%s_11192024.png',
                      ifelse(params$dataset == 'all', 'fullpop', 'AIAN'),
                      ifelse(params$use_softmax, 'softmax', 'directest'),
                      ifelse(params$preprocess_scale, 'centering', 'NOcentering'),
                      ifelse(params$alpha_variance_prior == -1, 'NOalpha','alpha'),
                      params$fixed_effects)
  print(f)
  print(plot_name)
  print('----------------')
  p1 <- plot_real_results(data_list = res$sim_list$data_list, 
                          stan_fit = res$sim_list$stan_fit,
                          stan_summary = res$sim_list$stan_summary$summary,
                          CV_pred = res$sim_list$CV_pred,
                          alpha_estimates = (params$alpha_variance_prior != -1))
  ggsave(plot = p1, filename = plot_name, height = 14, width = 7)
  
}

#
#### COVID rates ####
# Using our model 2020 estimates, census 2020 estimates, or PEP 2019 estimates, compute the COVID rates in all US counties for AIAN population in 2020.

# Pull in covid data
load('data/AmericanIndian_COVID_Deaths_2020.RData')

make_chloropleth_plot(cvd, 'AIAN_deaths')

#
#### Chloropleth plots! ####
setwd(root_results)
load('real_data/real_data_fit_directest_interceptonly_ID70259_2024_11_12.RData')

df<- res$sim_list$data_list$data
val <- 'census'

data <- merge(df, counties, by = 'GEOID')


  #
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






