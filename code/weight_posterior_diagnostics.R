### Posterior diagnostics of SWEEP ensemble weights (u) for the four real-data fits
### used in the paper (Figures 2 and 4). Written for the AOAS revision, reviewer
### comment 1.2 (dispersion of SM vs DE weights), Nov 2026.
###
### Sections:
###   1. Posterior uncertainty of the weights (posterior SD, 95% CI width) vs. the
###      across-county spread of posterior medians (what Figures 2 and 4 display).
###   2. Within-county posterior correlation between product weights.
###   3. Convergence (Rhat / n_eff of u), prediction uncertainty (mu), tau2.
###   4. Per-chain weights for one fit, plus a chain x product grid plot of the
###      county-level within-chain median weights.
###   5. Full-population SM only: per-chain log posterior, and whether the
###      PEP-led and ACS-led chain modes give different predictions / MAPE.
###
### Usage (from the repo root):
###   Rscript code/weight_posterior_diagnostics.R "<folder with the .RData fits>" [section ...]
###   e.g. Rscript code/weight_posterior_diagnostics.R "C:/.../Pop Est Results (tmp)" 4
### With no sections listed, all five run. Full-population fits are 4-5 GB each.

suppressMessages({
  library(rstan)
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
fit_dir <- if (length(args) >= 1) args[1] else
  'C:/Users/Admin-Dell/Dropbox/Academic/Research/Pop Est Results (tmp)'
sections <- if (length(args) >= 2) as.integer(args[-1]) else 1:5
out_dir <- 'pic/ensemble_weights'
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Fits behind the paper's Figures 2 and 4 (see post_processing_PE.R, ~line 845).
fit_ids <- c(full_SM = 'ID73323', full_DE = 'ID44508',
             aian_SM = 'ID88451', aian_DE = 'ID32362')
fit_files <- list.files(fit_dir, pattern = '[.]RData$', full.names = TRUE)

load_fit <- function(nm) {
  f <- grep(fit_ids[[nm]], fit_files, value = TRUE)
  if (length(f) != 1) stop(sprintf('Expected one file for %s in %s', nm, fit_dir))
  env <- new.env(); load(f, envir = env)
  env$res
}

options(width = 200, digits = 3)
tab_uncert <- list(); tab_cor <- list(); tab_conv <- list()

if (any(sections %in% 1:3)) {
  for (nm in names(fit_ids)) {
    res <- load_fit(nm)
    models <- toupper(res$models)
    S <- res$sim_list$stan_summary$summary

    ## 1. posterior uncertainty vs across-county spread
    if (1 %in% sections) {
      for (j in seq_along(models)) {
        ind <- grep(sprintf('^u[[][0-9]+,%d[]]$', j), rownames(S))
        med <- S[ind, '50%']; sdv <- S[ind, 'sd']; wid <- S[ind, '97.5%'] - S[ind, '2.5%']
        tab_uncert[[length(tab_uncert) + 1]] <- data.frame(
          fit = nm, product = models[j], n_counties = length(ind),
          mean_post_sd = mean(sdv), median_post_sd = median(sdv),
          mean_ci95_width = mean(wid), median_ci95_width = median(wid),
          sd_across_counties_of_medians = sd(med))
      }
    }

    ## 2. within-county posterior correlation between weights
    if (2 %in% sections) {
      u <- rstan::extract(res$sim_list$stan_fit, pars = 'u')$u  # draws x N x M
      prs <- combn(seq_along(models), 2)
      for (k in seq_len(ncol(prs))) {
        a <- prs[1, k]; b <- prs[2, k]
        r <- sapply(seq_len(dim(u)[2]), function(i) suppressWarnings(cor(u[, i, a], u[, i, b])))
        tab_cor[[length(tab_cor) + 1]] <- data.frame(
          fit = nm, pair = paste(models[a], models[b], sep = '-'),
          mean_within_county_cor = mean(r, na.rm = TRUE),
          median_within_county_cor = median(r, na.rm = TRUE))
      }
      rm(u)
    }

    ## 3. convergence, prediction uncertainty, tau2
    if (3 %in% sections) {
      iu <- grep('^u[[]', rownames(S)); imu <- grep('^mu[[]', rownames(S))
      it <- grep('^tau2[[]', rownames(S))
      tab_conv[[nm]] <- data.frame(
        fit = nm,
        u_Rhat_median = median(S[iu, 'Rhat']), u_Rhat_max = max(S[iu, 'Rhat']),
        u_frac_Rhat_gt_1.05 = mean(S[iu, 'Rhat'] > 1.05),
        u_neff_median = median(S[iu, 'n_eff']),
        mu_relSD_median = median(S[imu, 'sd'] / S[imu, 'mean']),
        mu_relCIwidth_median = median((S[imu, '97.5%'] - S[imu, '2.5%']) / S[imu, 'mean']),
        tau2_means = paste(signif(S[it, 'mean'], 3), collapse = ','))
    }
    cat('processed', nm, '\n')
    rm(res); gc(verbose = FALSE)
  }
  if (length(tab_uncert)) { cat('\n=== 1. Posterior uncertainty of weights ===\n'); print(do.call(rbind, tab_uncert), row.names = FALSE) }
  if (length(tab_cor))    { cat('\n=== 2. Within-county posterior correlation ===\n'); print(do.call(rbind, tab_cor), row.names = FALSE) }
  if (length(tab_conv))   { cat('\n=== 3. Convergence / prediction uncertainty ===\n'); print(do.call(rbind, tab_conv), row.names = FALSE) }
}

## 4. per-chain weights: chain x product grid plot
plot_chain_weights <- function(nm, file_stub) {
  res <- load_fit(nm)
  models <- toupper(res$models)
  u <- rstan::extract(res$sim_list$stan_fit, pars = 'u', permuted = FALSE)  # iter x chain x param
  pn <- dimnames(u)[[3]]
  n_chain <- dim(u)[2]

  # county-level within-chain posterior median, matching how Figures 2/4 summarize weights
  df <- do.call(rbind, lapply(seq_along(models), function(j) {
    idx <- grep(sprintf('^u[[][0-9]+,%d[]]$', j), pn)
    do.call(rbind, lapply(seq_len(n_chain), function(ch) {
      data.frame(chain = ch, product = models[j], county = seq_along(idx),
                 weight = apply(u[, ch, idx], 2, median),
                 within_sd = apply(u[, ch, idx], 2, sd))
    }))
  }))

  # label each chain by which product it gives the highest average weight
  lead <- df %>% group_by(chain, product) %>% summarize(m = mean(weight), .groups = 'drop') %>%
    group_by(chain) %>% slice_max(m, n = 1) %>% ungroup()
  chain_lab <- setNames(sprintf('Chain %d\n(%s-led)', lead$chain, lead$product), lead$chain)
  df$chain_f <- factor(chain_lab[as.character(df$chain)], levels = chain_lab[order(lead$chain)])
  df$product <- factor(df$product, levels = c('ACS', 'PEP', 'WP'))

  summ <- df %>% group_by(chain_f, product) %>%
    summarize(med = median(weight), sd_w = mean(within_sd), .groups = 'drop')

  cat('\n=== 4. Per-chain summary for', nm, '===\n')
  print(as.data.frame(summ %>% mutate(chain_f = gsub('\n', ' ', chain_f))), row.names = FALSE)

  # categorical slots 1-3 of the default validated palette, fixed order
  pal <- c(ACS = '#2a78d6', PEP = '#eb6834', WP = '#1baf7a')
  # put the median label on the side away from the line
  summ <- summ %>% mutate(lab_x = ifelse(med > 0.5, 0.02, 0.98), lab_h = ifelse(med > 0.5, 0, 1))

  p <- ggplot(df, aes(x = weight, fill = product)) +
    # each panel scaled to its own peak so the WP spike doesn't flatten ACS/PEP
    geom_histogram(aes(y = after_stat(ncount)), binwidth = 0.02, boundary = 0, color = NA) +
    geom_vline(data = summ, aes(xintercept = med), linetype = 'dashed', linewidth = 0.4,
               color = '#3d3d3a') +
    geom_text(data = summ, aes(x = lab_x, y = 1, hjust = lab_h, label = sprintf('median %.2f', med)),
              vjust = 1, size = 2.6, color = '#3d3d3a', inherit.aes = FALSE) +
    facet_grid(chain_f ~ product) +
    scale_fill_manual(values = pal, guide = 'none') +
    scale_x_continuous(breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(breaks = NULL) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1.05)) +
    labs(x = 'County-level posterior median weight (within chain)',
         y = 'Share of counties (each panel scaled to its peak)',
         title = sprintf('Ensemble weights by MCMC chain: %s', nm),
         subtitle = 'Each panel: distribution across counties of the within-chain posterior median weight; dashed line = median across counties') +
    theme_minimal(base_size = 9) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          strip.text.y = element_text(angle = 0, hjust = 0),
          strip.text.x = element_text(face = 'bold', size = 10),
          panel.spacing.y = unit(0.25, 'lines'),
          plot.subtitle = element_text(size = 7.5, color = '#5f5e5a'))

  f_out <- file.path(out_dir, sprintf('%s_%s_chain_weights.png', format(Sys.Date(), '%m%d%Y'), file_stub))
  ggsave(f_out, p, width = 7.5, height = 11, dpi = 200, bg = 'white')
  cat('saved', f_out, '\n')
  # paper version: title/subtitle carried by the LaTeX caption instead
  f_paper <- sub('[.]png$', '_paper.png', f_out)
  ggsave(f_paper, p + labs(title = NULL, subtitle = NULL), width = 7.5, height = 10.5, dpi = 200, bg = 'white')
  cat('saved', f_paper, '\n')
  invisible(df)
}

if (4 %in% sections) {
  plot_chain_weights('full_SM', 'fullpop_SM')
}

## 5. Do the chain modes give different predictions? (full-population SM)
##    Compares county-level posterior mean of mu between PEP-led and ACS-led chains,
##    and each mode's in-sample MAPE against the census.
if (5 %in% sections) {
  res <- load_fit('full_SM')
  fit <- res$sim_list$stan_fit
  u  <- rstan::extract(fit, pars = 'u',  permuted = FALSE)
  mu <- rstan::extract(fit, pars = 'mu', permuted = FALSE)
  pep_j <- which(toupper(res$models) == 'PEP')
  pep_idx <- grep(sprintf('^u[[][0-9]+,%d[]]$', pep_j), dimnames(u)[[3]])
  chain_pep <- sapply(seq_len(dim(u)[2]), function(ch) mean(u[, ch, pep_idx]))
  pep_led <- chain_pep > 0.5
  lp <- rstan::extract(fit, pars = 'lp__', permuted = FALSE)[, , 1]
  cat('Mean lp__ by chain:', paste(sprintf('%d:%.0f', seq_along(chain_pep), colMeans(lp)), collapse = ' '), '\n')
  mu_mode <- function(chs) colMeans(apply(mu[, chs, , drop = FALSE], 3, c))
  mu_pep <- mu_mode(which(pep_led)); mu_acs <- mu_mode(which(!pep_led))
  dat <- res$sim_list$data_list$data
  y_col <- intersect(c('census', 'y'), colnames(dat))[1]
  y <- dat[[y_col]]
  cat('outcome column used:', y_col, ' (n =', length(y), ', NA =', sum(is.na(y)), ')\n')

  # negative binomial log likelihood per draw (same likelihood as the Stan model)
  theta <- rstan::extract(fit, pars = 'theta', permuted = FALSE)[, , 1]  # iter x chain
  obs <- which(!is.na(y))
  ll <- sapply(seq_len(dim(mu)[2]), function(ch) sapply(seq_len(dim(mu)[1]), function(it)
    sum(dnbinom(y[obs], size = theta[it, ch], mu = pmax(mu[it, ch, obs], 0.01), log = TRUE))))
  cat('Mean log likelihood by chain:', paste(sprintf('%d:%.0f', seq_along(chain_pep), colMeans(ll)), collapse = ' '), '\n')
  cat(sprintf('Mean log likelihood: PEP-led %.1f, ACS-led %.1f | mean lp__: PEP-led %.1f, ACS-led %.1f\n',
              mean(ll[, pep_led]), mean(ll[, !pep_led]), mean(lp[, pep_led]), mean(lp[, !pep_led])))
  cat(sprintf('Mean theta: PEP-led %.1f, ACS-led %.1f\n', mean(theta[, pep_led]), mean(theta[, !pep_led])))

  rel_diff <- abs(mu_pep - mu_acs) / ((mu_pep + mu_acs) / 2)
  cat('\n=== 5. Predictions by chain mode (full_SM) ===\n')
  cat('PEP-led chains:', which(pep_led), ' ACS-led chains:', which(!pep_led), '\n')
  cat(sprintf('Relative difference in county predicted mean between modes: median %.4f, 90th pct %.4f, max %.4f\n',
              median(rel_diff), quantile(rel_diff, 0.9), max(rel_diff)))
  if (!is.null(y)) {
    cat(sprintf('In-sample MAPE vs census: PEP-led %.2f%%, ACS-led %.2f%%\n',
                100 * mean(abs(mu_pep - y) / y), 100 * mean(abs(mu_acs - y) / y)))
  }
}
