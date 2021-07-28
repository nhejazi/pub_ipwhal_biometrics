# This script is a variant of 02_est_npipw_hal.R, butchered so as to return the
# IPW estimates across a grid of L1-norms for use in better assessing the
# capabilities of the plateau selector based on Lepski's method

# use custom package library on Savio cluster
if (grepl("savio2", Sys.info()["nodename"])) {
  .libPaths("/global/scratch/nhejazi/R")
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
}

# read in command line arguments
args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE,
                             defaults = list(dgp = "1a",
                                             Mncv_max = 50,
                                             trunc_max = 0.2,
                                             trunc_len = 10,
                                             ipw = "ht",
                                             Q_reg = "Qn_hal",
                                             g_trunc = "profile",
                                             dcar_compat = TRUE,
                                             selector = "dcar"))

# reference for logging
print(args)

# packages, housekeeping
library(here)
library(data.table)
library(tidyverse)
library(origami)
library(hal9001)
library(foreach)
library(future)
library(doFuture)
library(doRNG)

# load helper package and scripts
source(here("R", "01_dgp.R"))
devtools::load_all(here("..", "uhalipw"))
source(here("R", "selectors.R"))
source(here("R", "utils.R"))

# parallelization
options(future.globals.maxSize = 1e32)
registerDoFuture()
plan(multiprocess, workers = 24)
set.seed(72439)

# fixed simulation parameters (not passed in via command line args)
n_iter <- 1000
v_folds <- 5
tsm_contrast <- 1                                # estimating Q(1,W) or Q(0,W)
n_obs <- 1000
lambda_seq <- exp(seq(-0.5, -20, length = 1e4))
gtrunc_seq <- round(seq(0.01, args$trunc_max, length.out = args$trunc_len), 2)

# truth for data-generating process
sim_truth <- get_truth(n_samp = 1e7, tsm_contrast = tsm_contrast,
                       dgp = args$dgp)

# parallelized loop for number of simulations
sim_iter <- foreach(this_iter = seq_len(n_iter),
                    .options.multicore = list(preschedule = FALSE),
                    .errorhandling = "stop") %dorng% {
  ## 0) generate data for this simulation
  dgp <- make_simple_data(n_obs, dgp = args$dgp)
  data_o <- dgp$data_obs
  folds <- make_folds(data_o, V = v_folds)
  folds_by_ids <- folds2foldvec(folds)
  names_w <- str_subset(colnames(data_o), "W")
  names_noty <- colnames(data_o)[!str_detect(colnames(data_o), "Y")]

  ### also, compute true g_0 and Q_0 from structural equations
  g0 <- with(data_o, dgp$dgp_funs$g0(W1 = W1, W2 = W2))
  Q0 <- with(data_o, dgp$dgp_funs$Q0(A = tsm_contrast, W1 = W1, W2 = W2))

  ### make counterfactual (intervened) data for prediction downstream
  data_tsm <- copy(data_o)[, A := tsm_contrast]

  ## 1) first, run CV-HAL to get CV-estimate of L1 norm and sequence
  gn_cvhal_fit <- fit_hal(X = as.matrix(data_o[, ..names_w]),
                           Y = data_o$A,
                           max_degree = NULL,
                           fit_type = "glmnet",
                           foldid = folds_by_ids,
                           use_min = FALSE,
                           family = "binomial",
                           cv_select = TRUE,
                           standardize = FALSE,
                           yolo = FALSE)
  gn_cv_pred <- predict(gn_cvhal_fit, new_data = data_o)
  L1_cv <- sum(abs(gn_cvhal_fit$coef))
  lambda_cv <- gn_cvhal_fit$lambda_star
  L1_seq <- seq(L1_cv, args$Mncv_max * L1_cv, length.out = 1000)

  ## 2) run under-smoothed HAL for the sequence of L1 norms
  ## NOTE: something really odd happens here (at least in some cases), where
  #        the first predicted set of propensity scores (w/ L1 norm matching
  #        matching that from global CV) does not match gn_cv_pred well?!
  gn_ipw_hal <- fit_gn_ipw(data_in = data_o, folds = folds,
                           x_names = names_w, y_names = "A",
                           L1_seq = L1_seq, lambda_seq = lambda_seq)

  ## 3) need Q_n --- but this doesn't need to be undersmoothed: use CV-HAL,
  ##    or use truth based on the data-generating mechanism, or just an LM.
  ### fit and predict using CV-HAL without undersmoothing
  if (args$Q_reg == "Qn_hal") {
    Qn_cvhal_fit <- fit_hal(X = as.matrix(data_o[, ..names_noty]),
                            Y = data_o$Y,
                            max_degree = NULL,
                            fit_type = "glmnet",
                            foldid = folds_by_ids,
                            use_min = FALSE,
                            family = "gaussian",
                            cv_select = TRUE,
                            standardize = FALSE,
                            yolo = FALSE)
    Qn_tsm_pred <- predict(Qn_cvhal_fit,
                           new_data = as.matrix(data_tsm[, ..names_noty]))
  } else if (args$Q_reg == "Qn_lm") {
    ### fit and predict with a simple linear model too
    Qn_cvlm_fit <- cross_validate(cv_fun = cv_lm, folds = folds,
                                  obs_data = data_o, tsm_data = data_tsm,
                                  reg_form = "Y ~ .^2", .combine = FALSE)
    fold_pred_idx <- order(do.call(c, lapply(folds, `[[`, "validation_set")))
    Qn_tsm_pred <- do.call(c, Qn_cvlm_fit$predictions)[fold_pred_idx]
  } else if (args$Q_reg == "Q0") {
    Qn_tsm_pred <- Q0
  }

  ## 4) build IPW estimator and compute variance based on influence function
  ipw_gtrunc <- lapply(gtrunc_seq, function(gtrunc) {
    # use an incompatible weight in D_CAR to easen tension between the score
    # and weights, i.e., the incompatible selector still uses (A - gn_uhal)
    # but switches out the term (Qn / gn_uhal) for (Qn / gn_dcar)
    if (!args$dcar_compat) {
      ipw <- build_ipw(gn_uhal = gn_ipw_hal$gn_cv,
                       a = data_o$A,
                       y = data_o$Y,
                       Qn_tsm = Qn_tsm_pred,
                       gn_dcar = gn_cv_pred,
                       bound_gn = gtrunc,
                       est_type = args$ipw)
    } else {
      ipw <- build_ipw(gn_uhal = gn_ipw_hal$gn_cv,
                       a = data_o$A,
                       y = data_o$Y,
                       Qn_tsm = Qn_tsm_pred,
                       bound_gn = gtrunc,
                       est_type = args$ipw)
    }
    return(ipw)
  }) %>%
  set_names(paste0("g_delta_", gtrunc_seq))

  ### for reference, also compute IPW and IC for gn selected by global CV
  ipw_gcv <- build_ipw(gn_uhal = as.matrix(gn_cv_pred),
                       a = data_o$A,
                       y = data_o$Y,
                       Qn_tsm = Qn_tsm_pred,
                       est_type = args$ipw)
  psi_ipw_gcv <- ipw_gcv$psi
  ic_ipw_gcv <- as.numeric(ipw_gcv$ic)
  var_ic_ipw_gcv <- var(ic_ipw_gcv - as.numeric(ipw_gcv$dcar)) / n_obs
  gn_ipw_gcv <- as.numeric(ipw_gcv$gn)
  dcar_mean_gcv <- mean(as.numeric(ipw_gcv$dcar))
  dvar_eff <- sqrt(sim_truth$effic_bound / n_obs)
  dcar_crit_gcv <- dcar_mean_gcv / dvar_eff

  # find regularization parameter and L1 norm selected by undersmoothing
  ipw_across_plateau <- lapply(ipw_gtrunc, function(ipw) {
    lambda_map <- lambda_seq[as.numeric(str_remove(colnames(ipw$dcar), "s"))]
    var_ic <- ipw$var_ic
    psi_ipw <- ipw$psi
    out <- list(psi = psi_ipw, var_ic = var_ic, lambda_seq = lambda_map,
                L1_seq = L1_seq) %>%
      as_tibble()
    return(out)
  }) %>%
  set_names(gtrunc_seq) %>%
  bind_rows(.id = "g_trunc") %>%
  mutate(g_trunc = as.numeric(g_trunc))
  return(ipw_across_plateau)
}

# re-organize output into separate frames
out <- bind_rows(sim_iter, .id = "iter") %>%
  mutate(
    lambda_seq = round(lambda_seq, 4),
    L1_seq = round(L1_seq, 2),
  )

# save final output
saveRDS(object = out,
        file = here("data", "plateau", "test_mc_var.rds"))


###############################################################################

# assess results locally
if (!grepl("savio2", Sys.info()["nodename"])) {
  # load data
  library(here)
  test_data_file <- "test_mc_var.rds"
  test_data <- readRDS(here("data", "plateau", test_data_file))

  # make summary plot
  data_to_plot <- test_data %>%
    mutate(L1_round = round(L1_seq, 1)) %>%
    dplyr::filter(L1_round < 300) %>%
    group_by(L1_round, g_trunc) %>%
    summarize(
      mc_var = mean(var(psi, na.rm = TRUE))
    )

  p_plateau <- data_to_plot %>%
    ggplot(aes(x = L1_round, y = mc_var)) +
    geom_point(alpha = 0.5) +
    #geom_line() +
    theme_bw()
    #facet_wrap(vars(g_trunc))
  ggsave(filename = "~/plateau_test.pdf", plot = p_plateau)
}
