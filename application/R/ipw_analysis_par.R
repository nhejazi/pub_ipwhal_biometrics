if (grepl("savio2", Sys.info()["nodename"])) {
  .libPaths("/global/scratch/nhejazi/R")
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
}

# read in command line arguments
args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE,
                             defaults = list(selector_type = "dcar"))

# reference for logging
print(args)

library(here)
library(tidyverse)
library(data.table)
library(origami)
library(hal9001)
library(future)
library(future.apply)
devtools::load_all(here("..", "uhalipw"))
source(here("..", "simulations", "R", "selectors.R"))

# parallelization options
options(future.globals.maxSize = 1e10)
plan(multiprocess, workers = 20)

# constants used throughout analysis
v_folds <- 5
ci_level <- 0.95
ci_mult <- abs(stats::qnorm(p = (1 - ci_level) / 2))
set.seed(281681)

# read in data and modify as per Hernan & Robins
nhefs <- read_csv(here("data", "nhefs.csv"))
nhefs$cens <- ifelse(is.na(nhefs$wt82), 1, 0)

# provisionally ignore subjects with missing values for weight in 1982
nhefs_nmv <- nhefs[which(!is.na(nhefs$wt82)), ]

###############################################################################
# Estimation of IP weights via logistic regression
###############################################################################
## NOTE: code for this analysis was adapted from Hernan & Robins (2020)
## https://remlapmot.github.io/cibookex-r/ip-weighting-and-marginal-structural-models.html#program-12.2

## 1) fit logistic regression model of smoking cessation
ip_glm <- glm(qsmk ~ sex + race + age + as.factor(education) + smokeintensity +
              smokeyrs + as.factor(exercise) + as.factor(active) + wt71,
              family = binomial(), data = nhefs_nmv)
### we'll also use the logistic model from Hernan & Robins - why have ^2 terms?
ip_glm_robins <- glm(
  qsmk ~ sex + race + age + I(age ^ 2) +
    as.factor(education) + smokeintensity +
    I(smokeintensity ^ 2) + smokeyrs + I(smokeyrs ^ 2) +
    as.factor(exercise) + as.factor(active) + wt71 + I(wt71 ^ 2),
  family = binomial(),
  data = nhefs_nmv
)

## 2) compute predicted probability of smoking cessation
ipw_classical <- lapply(list(ip_glm, ip_glm_robins), function(ip_glm) {
# predict for the given IPW model
  pred_qsmk_obs <- ifelse(nhefs_nmv$qsmk == 0,
                          1 - predict(ip_glm, type = "response"),
                          predict(ip_glm, type = "response"))
  nhefs_nmv$ip_glm <- 1 / pred_qsmk_obs

  ## 3) construct IPW estimator based on parametric IP weights
  tsm_tx_ipw_glm <- mean(with(nhefs_nmv, (qsmk * ip_glm) * wt82_71))
  eif_tx_ipw_glm <- with(nhefs_nmv, (qsmk * ip_glm) * wt82_71) - tsm_tx_ipw_glm
  var_tx_ipw_glm <- var(eif_tx_ipw_glm) / nrow(nhefs_nmv)
  psi_par_ipw_glm <- tsm_tx_ipw_glm - mean(nhefs_nmv$wt82_71)
  eif_mean_y <- with(nhefs_nmv, wt82_71 - mean(wt82_71))
  var_par_ipw_glm <- var(eif_tx_ipw_glm - eif_mean_y) / nrow(nhefs_nmv)

  ## 4) make summary table
  ipw_glm_summary <- list(qsmk = c("yes", "par"),
                          psi_ipw = c(tsm_tx_ipw_glm, psi_par_ipw_glm),
                          var_ipw = c(var_tx_ipw_glm, var_par_ipw_glm)) %>%
    as_tibble() %>%
    mutate(
      ci_lwr = psi_ipw - ci_mult * sqrt(var_ipw),
      ci_upr = psi_ipw + ci_mult * sqrt(var_ipw)
    )
  return(ipw_glm_summary)
}) %>% set_names(c("main_terms", "quadratic_robins"))

###############################################################################
# Estimation of IP weights via undersmoothed HAL
###############################################################################

## 1) set up simulation parameters and dataset
Mncv_max <- 500
lambda_seq <- exp(seq(-0.5, -20, length = 1e4))
gtrunc_seq <- round(seq(0.01, 0.2, length.out = 10), 2)

### setting up data
nhefs_nmv <- as.data.table(nhefs_nmv)
data_tx_tsm <- copy(nhefs_nmv)[, qsmk := 1]
folds <- make_folds(nhefs_nmv, V = v_folds)
folds_by_ids <- folds2foldvec(folds)
names_w <- c("sex", "race", "age", "education", "smokeintensity", "smokeyrs",
             "exercise", "active", "wt71")
names_noty <- c("qsmk", names_w)

## 2) run CV-HAL to get CV-estimate of L1 norm and sequence
gn_cvhal_fit <- fit_hal(X = as.matrix(nhefs_nmv[, ..names_w]),
                        Y = nhefs_nmv$qsmk,
                        max_degree = 4,
                        fit_type = "glmnet",
                        foldid = folds_by_ids,
                        use_min = TRUE,
                        family = "binomial",
                        cv_select = TRUE,
                        standardize = FALSE,
                        yolo = FALSE)
gn_cv_pred <- predict(gn_cvhal_fit, new_data = nhefs_nmv)
L1_cv <- sum(abs(gn_cvhal_fit$coef))
lambda_cv <- gn_cvhal_fit$lambda_star
L1_seq <- seq(L1_cv, Mncv_max * L1_cv, length.out = 5000)

## 2) run under-smoothed HAL for the sequence of L1 norms
gn_ipw_hal <- fit_gn_ipw(data_in = nhefs_nmv, folds = folds,
                         x_names = names_w, y_names = "qsmk",
                         L1_seq = L1_seq, lambda_seq = lambda_seq)

## 3) need Q_n --- but this doesn't need to be undersmoothed: use CV-HAL
Qn_cvhal_fit <- fit_hal(X = as.matrix(nhefs_nmv[, ..names_noty]),
                        Y = nhefs_nmv$wt82_71,
                        max_degree = 4,
                        fit_type = "glmnet",
                        foldid = folds_by_ids,
                        use_min = TRUE,
                        family = "gaussian",
                        cv_select = TRUE,
                        standardize = FALSE,
                        yolo = FALSE)
Qn_tx_tsm_pred <- predict(Qn_cvhal_fit,
                          new_data = as.matrix(data_tx_tsm[, ..names_noty]))

## 4) build IPW estimators and compute variance based on influence function
ipw_tx_gtrunc <- lapply(gtrunc_seq, function(gtrunc) {
  ipw <- build_ipw(gn_uhal = gn_ipw_hal$gn_cv,
                   a = nhefs_nmv$qsmk,
                   y = nhefs_nmv$wt82_71,
                   Qn_tsm = Qn_tx_tsm_pred,
                   bound_gn = gtrunc,
                   est_type = "ht")
  return(ipw)
}) %>% set_names(paste0("g_delta_", gtrunc_seq))


## 5) select IPW estimator and gn as minimizer of D_CAR/score criterion
### compute estimator + EIF for truncated estimator variant
if (args$selector_type == "dcar") {
  tx_est_trunc <- dcar_selector(ipw_gtrunc = ipw_tx_gtrunc,
                                L1_seq = L1_seq,
                                est_type = "usm",
                                gtrunc_type = "joint",
                                use_eff_bound = FALSE)
  print(tx_est_trunc$dcar_mean)
} else if (args$selector_type == "score") {
  tx_est_trunc <- score_selector(ipw_gtrunc = ipw_tx_gtrunc,
                                 gn_usm_fit = gn_ipw_hal,
                                 A_obs = as.numeric(nhefs_nmv$qsmk),
                                 L1_seq = L1_seq,
                                 est_type = "usm_mintrunc")
  print(tx_est_trunc$score_crit)
}

ipw_hal_trunc_summary <- list(
  qsmk = c("yes", "par"),
  psi_ipw = c(tx_est_trunc$psi, tx_est_trunc$psi - mean(nhefs_nmv$wt82_71)),
  var_ipw = c(tx_est_trunc$var_ic,
              var(tx_est_trunc$ic - with(nhefs_nmv, wt82_71 - mean(wt82_71))) /
                nrow(nhefs_nmv))
) %>%
as_tibble() %>%
mutate(
  ci_lwr = psi_ipw - ci_mult * sqrt(var_ipw),
  ci_upr = psi_ipw + ci_mult * sqrt(var_ipw)
)

### compute estimator + EIF for minimally truncated estimator variant
tx_est_mintrunc <- dcar_selector(ipw = ipw_tx_gtrunc,
                                 L1_seq = L1_seq,
                                 est_type = "usm_mintrunc",
                                 gtrunc_type = "joint",
                                 use_eff_bound = FALSE)
print(tx_est_mintrunc$dcar_mean)
ipw_hal_mintrunc_summary <- list(
  qsmk = c("yes", "ate"),
  psi_ipw = c(tx_est_mintrunc$psi,
              tx_est_mintrunc$psi - mean(nhefs_nmv$wt82_71)),
  var_ipw = c(tx_est_mintrunc$var_ic,
              var(tx_est_mintrunc$ic -
                  with(nhefs_nmv, wt82_71 - mean(wt82_71))) /
              nrow(nhefs_nmv))
) %>%
as_tibble() %>%
mutate(
  ci_lwr = psi_ipw - ci_mult * sqrt(var_ipw),
  ci_upr = psi_ipw + ci_mult * sqrt(var_ipw)
)


## 6) for reference, also compute IPW and IC for gn selected by global CV
ipw_gcv_tx <- build_ipw(gn_uhal = as.matrix(gn_cv_pred),
                        a = nhefs_nmv$qsmk,
                        y = nhefs_nmv$wt82_71,
                        Qn_tsm = Qn_tx_tsm_pred,
                        est_type = "ht")
ipw_gcv_summary <- list(
  qsmk = c("yes", "ate"),
  psi_ipw = c(ipw_gcv_tx$psi, ipw_gcv_tx$psi - mean(nhefs_nmv$wt82_71)),
  var_ipw = c(ipw_gcv_tx$var_ic,
              var(ipw_gcv_tx$var_ic -
                  with(nhefs_nmv, wt82_71 - mean(wt82_71))) /
              nrow(nhefs_nmv))
) %>%
as_tibble() %>%
mutate(
  ci_lwr = psi_ipw - ci_mult * sqrt(var_ipw),
  ci_upr = psi_ipw + ci_mult * sqrt(var_ipw)
)


## 7) construct output from all procedures for tables and visualizations
ipw_results <- list(hal_trunc = ipw_hal_trunc_summary,
                    hal_mintrunc = ipw_hal_mintrunc_summary,
                    hal_gcv = ipw_gcv_summary,
                    glm_main = ipw_classical[[1]],
                    glm_quad = ipw_classical[[2]]) %>%
  bind_rows(.id = "gn_estim") %>%
  arrange(qsmk)
saveRDS(ipw_results, here("data", paste0("nhefs_ipw_", args$selector_type,
                                         "_results.rds")))
