if (grepl("savio2", Sys.info()["nodename"])) {
  .libPaths("/global/scratch/nhejazi/R")
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
}

# set virtual environment and load packages
renv::activate(project = here::here())
library(here)
library(tidyverse)
library(data.table)
library(origami)
library(hal9001)
library(cobalt)
devtools::load_all(here("..", "uhalipw"))
source(here("..", "simulations", "R", "selectors.R"))
set.seed(281681)

# load results from Ashkan
#load(here("data", "nhefs_ipwhal_ashkan.rda"))
#ps_hal_gcv <- dat$preds.hal.CV
#ps_hal_dcar <- dat$preds.hal.DCAR

# constants used throughout analysis
v_folds <- 5
ci_level <- 0.95
ci_mult <- abs(stats::qnorm(p = (1 - ci_level) / 2))

# read in data and modify as per Hernan & Robins
nhefs <- read_csv(here("data", "nhefs.csv"))
nhefs$cens <- ifelse(is.na(nhefs$wt82), 1, 0)

# provisionally ignore subjects with missing values for weight in 1982
nhefs_nmv <- nhefs[which(!is.na(nhefs$wt82)), ]

###############################################################################
# Estimation of IP weights via logistic regression
# NOTE: code for this analysis was adapted from Hernan & Robins (2020)
###############################################################################
## 0) fit the unadjusted logistic regression model of smoking cessation
ip_unadj <- glm(qsmk ~ 1, family = binomial(), data = nhefs_nmv)

## 1) fit logistic regression model of smoking cessation
ip_glm <- glm(qsmk ~ sex + race + age + as.factor(education) + smokeintensity +
              smokeyrs + as.factor(exercise) + as.factor(active) + wt71,
              family = binomial(), data = nhefs_nmv)

## 2) we'll also use the logistic model from Hernan & Robins
ip_glm_robins <- glm(
  qsmk ~ sex + race + age + I(age ^ 2) +
    as.factor(education) + smokeintensity +
    I(smokeintensity ^ 2) + smokeyrs + I(smokeyrs ^ 2) +
    as.factor(exercise) + as.factor(active) + wt71 + I(wt71 ^ 2),
  family = binomial(),
  data = nhefs_nmv
)

## 3) compute predicted probability of smoking cessation
ipw_classical <- lapply(list(ip_unadj, ip_glm, ip_glm_robins),
                        function(ip_glm) {
  # predict for the given IPW model
  pred_qsmk_obs <- ifelse(nhefs_nmv$qsmk == 0,
                          1 - predict(ip_glm, type = "response"),
                          predict(ip_glm, type = "response"))
  nhefs_nmv$ip_glm <- 1 / pred_qsmk_obs

  ## 3) construct IPW estimator based on parametric IP weights
  tsm_tx_ipw_glm <- mean(with(nhefs_nmv, (qsmk * ip_glm) * wt82_71))
  eif_tx_ipw_glm <- with(nhefs_nmv, (qsmk * ip_glm) * wt82_71) - tsm_tx_ipw_glm
  var_tx_ipw_glm <- var(eif_tx_ipw_glm) / nrow(nhefs_nmv)
  tsm_ct_ipw_glm <- mean(with(nhefs_nmv, ((1 - qsmk) * ip_glm) * wt82_71))
  eif_ct_ipw_glm <- with(nhefs_nmv, ((1 - qsmk) * ip_glm) * wt82_71) -
    tsm_ct_ipw_glm
  var_ct_ipw_glm <- var(eif_ct_ipw_glm) / nrow(nhefs_nmv)
  psi_ate_ipw_glm <- tsm_tx_ipw_glm - tsm_ct_ipw_glm
  var_ate_ipw_glm <- var(eif_tx_ipw_glm - eif_ct_ipw_glm) / nrow(nhefs_nmv)

  ## 4) make summary table
  ipw_glm_summary <- list(qsmk = c("no", "yes", "ate"),
                          psi_ipw = c(tsm_ct_ipw_glm, tsm_tx_ipw_glm,
                                      psi_ate_ipw_glm),
                          var_ipw = c(var_ct_ipw_glm, var_tx_ipw_glm,
                                      var_ate_ipw_glm)) %>%
    as_tibble() %>%
    mutate(
      ci_lwr = psi_ipw - ci_mult * sqrt(var_ipw),
      ci_upr = psi_ipw + ci_mult * sqrt(var_ipw)
    )
  return(ipw_glm_summary)
}) %>% set_names(c("unadj", "main_terms", "quadratic_robins"))

## 4) create balance table of IP weights
nhefs_nmv <- nhefs_nmv %>%
  mutate(
    ps_unadj = ip_unadj$fitted.values %>% unname(),
    ps_mainterms = ip_glm$fitted.values %>% unname(),
    ps_robins = ip_glm_robins$fitted.values %>% unname(),
    ipw_unadj = case_when(qsmk == 1 ~ 1 / ps_unadj,
                          qsmk == 0 ~ 1 / (1 - ps_unadj)),
    ipw_mainterms = case_when(qsmk == 1 ~ 1 / ps_mainterms,
                              qsmk == 0 ~ 1 / (1 - ps_mainterms)),
    ipw_robins = case_when(qsmk == 1 ~ 1 / ps_robins,
                           qsmk == 0 ~ 1 / (1 - ps_robins))
  )
ps_covars <- nhefs_nmv %>%
  select(sex, race, age, education, smokeintensity, smokeyrs,
         exercise, active, wt71)

###############################################################################
# Estimation of IP weights via undersmoothed HAL
###############################################################################

## 1) set up analysis parameters and dataset
Mncv_max <- 500
lambda_seq <- exp(seq(-0.5, -20, length = 1e4))
gtrunc_seq <- round(seq(0.01, 0.2, length.out = 20), 2)

### setting up data
nhefs_nmv <- as.data.table(nhefs_nmv)
data_tx_tsm <- copy(nhefs_nmv)[, qsmk := 1]
data_ct_tsm <- copy(nhefs_nmv)[, qsmk := 0]
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
Qn_ct_tsm_pred <- predict(Qn_cvhal_fit,
                          new_data = as.matrix(data_ct_tsm[, ..names_noty]))

## 4) build IPW estimators and compute variance based on influence function
ipw_tx_gtrunc <- lapply(gtrunc_seq, function(gtrunc) {
  ipw <- build_ipw(gn_uhal = gn_ipw_hal$gn_cv,
                   a = nhefs_nmv$qsmk,
                   y = nhefs_nmv$wt82_71,
                   Qn_tsm = Qn_tx_tsm_pred,
                   bound_gn = gtrunc,
                   est_type = "ht")
  return(ipw)
  }) %>%
  set_names(paste0("g_delta_", gtrunc_seq))

ipw_ct_gtrunc <- lapply(gtrunc_seq, function(gtrunc) {
  ipw <- build_ipw(gn_uhal = (1 - gn_ipw_hal$gn_cv),
                   a = (1 - nhefs_nmv$qsmk),
                   y = nhefs_nmv$wt82_71,
                   Qn_tsm = Qn_ct_tsm_pred,
                   bound_gn = gtrunc,
                   est_type = "ht")
  }) %>%
  set_names(paste0("g_delta_", gtrunc_seq))


## 5) select IPW estimator and gn as minimizer of D_CAR criterion
### compute estimator + EIF for truncated estimator variant
tx_est_trunc <- dcar_selector(ipw_gtrunc = ipw_tx_gtrunc,
                              L1_seq = L1_seq,
                              est_type = "usm",
                              gtrunc_type = "joint",
                              use_eff_bound = FALSE)
print(tx_est_trunc$dcar_mean)
ct_est_trunc <- dcar_selector(ipw_gtrunc = ipw_ct_gtrunc,
                              L1_seq = L1_seq,
                              est_type = "usm",
                              gtrunc_type = "joint",
                              use_eff_bound = FALSE)
print(ct_est_trunc$dcar_mean)

ipw_hal_trunc_summary <- list(
    qsmk = c("no", "yes", "ate"),
    psi_ipw = c(ct_est_trunc$psi, tx_est_trunc$psi,
                tx_est_trunc$psi - ct_est_trunc$psi),
    var_ipw = c(ct_est_trunc$var_ic, tx_est_trunc$var_ic,
                var(tx_est_trunc$ic - ct_est_trunc$ic)  / nrow(nhefs_nmv))
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
ct_est_mintrunc <- dcar_selector(ipw = ipw_ct_gtrunc,
                                 L1_seq = L1_seq,
                                 est_type = "usm_mintrunc",
                                 gtrunc_type = "joint",
                                 use_eff_bound = FALSE)
print(ct_est_mintrunc$dcar_mean)
ipw_hal_mintrunc_summary <- list(
    qsmk = c("no", "yes", "ate"),
    psi_ipw = c(ct_est_mintrunc$psi, tx_est_mintrunc$psi,
                tx_est_mintrunc$psi - ct_est_mintrunc$psi),
    var_ipw = c(ct_est_mintrunc$var_ic, tx_est_mintrunc$var_ic,
                var(tx_est_mintrunc$ic - ct_est_mintrunc$ic) / nrow(nhefs_nmv))
  ) %>%
  as_tibble() %>%
  mutate(
    ci_lwr = psi_ipw - ci_mult * sqrt(var_ipw),
    ci_upr = psi_ipw + ci_mult * sqrt(var_ipw)
  )


## 6) for reference, also compute IPW and IC for gn selected by global CV
ipw_gcv_tx <- build_ipw(gn_uhal = as.matrix(gn_ipw_hal$gn_cv[, 1]),
                        a = nhefs_nmv$qsmk,
                        y = nhefs_nmv$wt82_71,
                        Qn_tsm = Qn_tx_tsm_pred,
                        est_type = "ht")
ipw_gcv_ct <- build_ipw(gn_uhal = as.matrix(1 - gn_ipw_hal$gn_cv[, 1]),
                        a = (1 - nhefs_nmv$qsmk),
                        y = nhefs_nmv$wt82_71,
                        Qn_tsm = Qn_ct_tsm_pred,
                        est_type = "ht")
ipw_gcv_summary <- list(
    qsmk = c("no", "yes", "ate"),
    psi_ipw = c(ipw_gcv_ct$psi, ipw_gcv_tx$psi,
                ipw_gcv_tx$psi - ipw_gcv_ct$psi),
    var_ipw = c(ipw_gcv_ct$var_ic, ipw_gcv_tx$var_ic,
                var(ipw_gcv_tx$ic - ipw_gcv_ct$ic) / nrow(nhefs_nmv))
  ) %>%
  as_tibble() %>%
  mutate(
    ci_lwr = psi_ipw - ci_mult * sqrt(var_ipw),
    ci_upr = psi_ipw + ci_mult * sqrt(var_ipw)
  )

## add HAL-based IP weights to NHEFS data
nhefs_nmv <- nhefs_nmv %>%
  mutate(
    ps1_hal_gcv = as.numeric(ipw_gcv_tx$gn), #dat$preds.hal.CV,
    ps0_hal_gcv = 1 - ps1_hal_gcv,
    ps1_hal_mintrunc = as.numeric(tx_est_mintrunc$gn), #dat$preds.hal.DCAR,
    ps0_hal_mintrunc = 1 - ps1_hal_mintrunc,
    ps1_hal_trunc = as.numeric(tx_est_trunc$gn),
    ps0_hal_trunc = 1 - ps1_hal_trunc,
    ipw_hal_gcv = case_when(qsmk == 1 ~ 1 / ps1_hal_gcv,
                            qsmk == 0 ~ 1 / ps0_hal_gcv),
    ipw_hal_mintrunc = case_when(qsmk == 1 ~ 1 / ps1_hal_mintrunc,
                                 qsmk == 0 ~ 1 / ps0_hal_mintrunc),
    ipw_hal_trunc = case_when(qsmk == 1 ~ 1 / ps1_hal_trunc,
                              qsmk == 0 ~ 1 / ps0_hal_trunc)
  )

###############################################################################
# Results from GLM-based and HAL-based IPW analyses
###############################################################################

## 7) create plots of the IPW solution path in lambda
ipw_ate_psi_path <- ipw_tx_gtrunc$g_delta_0.01$psi -
  ipw_ct_gtrunc$g_delta_0.01$psi
ipw_ate_dcar_path <- ipw_tx_gtrunc$g_delta_0.01$dcar -
  ipw_ct_gtrunc$g_delta_0.01$dcar
ipw_ate_eif_path <- ipw_tx_gtrunc$g_delta_0.01$ic -
  ipw_ct_gtrunc$g_delta_0.01$ic
ipw_ate_path_plot_data <- as_tibble(list(
  l1_norm = L1_seq,
  psi = ipw_ate_psi_path,
  dcar_pn = abs(colMeans(ipw_ate_dcar_path)),
  eif_pn = colMeans(ipw_ate_eif_path)
))

## 8) create balance table for IPW variants
ipw_balance_table <- bal.tab(
  x = ps_covars,
  treat = "qsmk", data = nhefs_nmv,
  weights = c("ipw_unadj", "ipw_mainterms", "ipw_robins",
              "ipw_hal_gcv", "ipw_hal_mintrunc"),
  estimand = "ATE", int = TRUE, poly = 3
)

## 9) construct output from all procedures for tables and visualizations
ipw_results <- list(hal_trunc = ipw_hal_trunc_summary,
                    hal_mintrunc = ipw_hal_mintrunc_summary,
                    hal_gcv = ipw_gcv_summary,
                    glm_unadj = ipw_classical[[1]],
                    glm_main = ipw_classical[[2]],
                    glm_quad = ipw_classical[[3]]) %>%
  bind_rows(.id = "gn_estim") %>%
  arrange(qsmk)

## 10) create output list and save results object
nhefs_out <- list(results = ipw_results, balance = ipw_balance_table,
                  path_plot = ipw_ate_path_plot_data)
saveRDS(nhefs_out, here("data", "nhefs_ipw_results.rds"))
