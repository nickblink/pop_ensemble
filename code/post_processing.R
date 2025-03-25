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

#### 3/13/2025: Debugging AIAN SM alpha density results ####
setwd(root_results)
setwd('real_data')
load('real_data_fit_aian_softmax_alpha_cv10_pepfulldensity_ID17611_2025_01_15.RData')
fit <- res$sim_list$stan_fit

check_hmc_diagnostics(fit)

library(bayesplot)
mcmc_parcoord(as.array(fit), pars = c("rho", "tau2", "theta"))

# look at the pairs() to see if there are strong correlations between parameters.

# Try running with a higher adapt delta?

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
make_table_line <- function(metric, cols = c('dataset','MAPE', 'MAE', 'CP.95', 'med_int')){
  res <- apply(metric, 1, function(xx){
    paste(xx[cols], collapse = ' & ')
  })
  return(res)
}


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






