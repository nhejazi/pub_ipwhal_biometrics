# COMMAND LINE ARGUMENTS
# 1) dgp: which DGP to use for experiment (see file 01_dgp.R for details)
# 2) Mncv_max: multiple of the L1 norm up to which the sequence ought to go
# 3) ipw: "ht" for Horvitz-Thompson (classical) or "hajek" for stabilized
# 4) Q_reg: "Q0" for using the true outcome mechanism from the selected DGP,
#           "Qn_hal" for using an estimated outcome mechanism based on HAL, or
#           "Qn_iden" for ignoring the outcome mechanism and using just 1's
# 5) g_trunc: "profile" for sequential minimization over {lambda, delta} or
#             "joint" for joint minimization over {lambda, delta}

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
n_iter <- 500
v_folds <- 5
tsm_contrast <- 1                                # estimating Q(1,W) or Q(0,W)
n_samp <- cumsum(rep(sqrt(100), 4))^2            # sample sizes at root-n scale
lambda_seq <- exp(seq(-0.5, -20, length = 1e4))
gtrunc_seq <- round(seq(0.01, args$trunc_max, length.out = args$trunc_len), 2)

# truth for data-generating process
sim_truth <- get_truth(n_samp = 1e7, tsm_contrast = tsm_contrast,
                       dgp = args$dgp)

# loop over sample sizes
sim_over_nsamp <- lapply(n_samp, function(n_obs) {
  # parallelized loop for number of simulations
  sim_iter <- foreach(this_iter = seq_len(n_iter),
                      .options.multicore = list(preschedule = FALSE),
                      .errorhandling = "remove") %dorng% {
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

    ### also for reference, compute the unadjusted estimator and IC
    est_unadj <- data_o[A == tsm_contrast, mean(Y)]
    ic_unadj <- (data_o[, A] / data_o[, mean(A)]) * (data_o[, Y] - est_unadj)
    var_ic_unadj <- var(ic_unadj) / n_obs

    ## 5) select IPW estimator and gn as minimizer of given criterion
    ### NOTE: allow options for D_CAR selector and score-based selector
    if (args$selector == "dcar") {
      # compute estimator + EIF for truncated estimator variant
      est_trunc <- dcar_selector(ipw_gtrunc = ipw_gtrunc,
                                 L1_seq = L1_seq,
                                 est_type = "usm",
                                 gtrunc_type = args$g_trunc)
      crit_val_trunc <- est_trunc$dcar_mean

      # compute estimator + EIF for minimally truncated estimator variant
      est_mintrunc <- dcar_selector(ipw = ipw_gtrunc,
                                    L1_seq = L1_seq,
                                    est_type = "usm_mintrunc",
                                    gtrunc_type = args$g_trunc)
      crit_val_mintrunc <- est_mintrunc$dcar_mean

    # call helper function for score-based selector
    } else if (args$selector == "score") {
      # compute estimator + EIF for truncated estimator variant
      est_trunc <- score_selector(ipw_gtrunc = ipw_gtrunc,
                                  gn_usm_fit =  gn_ipw_hal,
                                  A_obs = data_o$A,
                                  L1_seq = L1_seq,
                                  est_type = "usm")
      crit_val_trunc <- est_trunc$score_crit

      # compute estimator + EIF for minimally truncated estimator variant
      est_mintrunc <- score_selector(ipw = ipw_gtrunc,
                                     gn_usm_fit =  gn_ipw_hal,
                                     A_obs = data_o$A,
                                     L1_seq = L1_seq,
                                     est_type = "usm_mintrunc")
      crit_val_mintrunc <- est_mintrunc$score_crit
    } else if (args$selector == "plateau") {
      # compute estimator + EIF for truncated estimator variant
      est_trunc <- plateau_selector(ipw_gtrunc = ipw_gtrunc,
                                    gn_usm_fit =  gn_ipw_hal,
                                    L1_seq = L1_seq,
                                    est_type = "usm")
      crit_val_trunc <- est_trunc$plateau_crit

      # compute estimator + EIF for minimally truncated estimator variant
      est_mintrunc <- plateau_selector(ipw = ipw_gtrunc,
                                       gn_usm_fit =  gn_ipw_hal,
                                       L1_seq = L1_seq,
                                       est_type = "usm_mintrunc")
      crit_val_mintrunc <- est_mintrunc$plateau_crit
    }

    ### remove HAL fit objects and predictions since now unnecessary
    if (args$Q_reg == "Qn_hal") {
      rm(Qn_cvhal_fit, gn_cvhal_fit, gn_ipw_hal, ipw_gtrunc); gc()
    }

    ## 6) construct output for single iteration
    sim_out <- list(# first, undersmoothed and truncated estimator
                    psi_usm_trunc = est_trunc$psi,
                    var_usm_trunc = est_trunc$var_ic,
                    ic_usm_trunc = est_trunc$ic,
                    gn_usm_trunc = est_trunc$gn,
                    Mn_usm_trunc = est_trunc$Mn,
                    lambda_usm_trunc = est_trunc$lambda,
                    gtrunc_usm_trunc = est_trunc$gtrunc,
                    crit_usm_trunc = crit_val_trunc,
                    # next, undersmoothed estimator w/ minimal truncation
                    psi_usm_mintrunc = est_mintrunc$psi,
                    var_usm_mintrunc = est_mintrunc$var_ic,
                    ic_usm_mintrunc = est_mintrunc$ic,
                    gn_usm_mintrunc = est_mintrunc$gn,
                    Mn_usm_mintrunc = est_mintrunc$Mn,
                    lambda_usm_mintrunc = est_mintrunc$lambda,
                    gtrunc_usm_mintrunc = est_mintrunc$gtrunc,
                    crit_usm_mintrunc = crit_val_mintrunc,
                    # last, estimator selected by global CV
                    psi_cv = psi_ipw_gcv,
                    var_cv = var_ic_ipw_gcv,
                    ic_cv = ic_ipw_gcv,
                    gn_cv = gn_ipw_gcv,
                    Mn_cv = L1_cv,
                    lambda_cv = lambda_cv,
                    crit_cv = dcar_mean_gcv,
                    # finally, the unadjusted estimator
                    psi_unadj = est_unadj,
                    var_unadj = var_ic_unadj,
                    ic_unadj = ic_unadj)
    return(sim_out)
  }

  # re-organize output into separate frames
  out <- list(# first, undersmoothed, truncated estimator
              psi_usm_trunc = sapply(sim_iter, `[[`, "psi_usm_trunc"),
              var_usm_trunc = sapply(sim_iter, `[[`, "var_usm_trunc"),
              ic_usm_trunc = sapply(sim_iter, `[[`, "ic_usm_trunc"),
              gn_usm_trunc = sapply(sim_iter, `[[`, "gn_usm_trunc"),
              Mn_usm_trunc = sapply(sim_iter, `[[`, "Mn_usm_trunc"),
              lambda_usm_trunc = sapply(sim_iter, `[[`, "lambda_usm_trunc"),
              gtrunc_usm_trunc = sapply(sim_iter, `[[`, "gtrunc_usm_trunc"),
              crit_usm_trunc = sapply(sim_iter, `[[`, "crit_usm_trunc"),
              # next, undersmoothed estimator w/ minimal truncation
              psi_usm_mintrunc = sapply(sim_iter, `[[`, "psi_usm_mintrunc"),
              var_usm_mintrunc = sapply(sim_iter, `[[`, "var_usm_mintrunc"),
              ic_usm_mintrunc = sapply(sim_iter, `[[`, "ic_usm_mintrunc"),
              gn_usm_mintrunc = sapply(sim_iter, `[[`, "gn_usm_mintrunc"),
              Mn_usm_mintrunc = sapply(sim_iter, `[[`, "Mn_usm_mintrunc"),
              lambda_usm_mintrunc = sapply(sim_iter, `[[`,
                                           "lambda_usm_mintrunc"),
              gtrunc_usm_mintrunc = sapply(sim_iter, `[[`,
                                           "gtrunc_usm_mintrunc"),
              crit_usm_mintrunc = sapply(sim_iter, `[[`, "crit_usm_mintrunc"),
              # last, estimator selected by global CV
              psi_gcv = sapply(sim_iter, `[[`, "psi_cv"),
              var_gcv = sapply(sim_iter, `[[`, "var_cv"),
              ic_gcv = sapply(sim_iter, `[[`, "ic_cv"),
              gn_gcv = sapply(sim_iter, `[[`, "gn_cv"),
              Mn_gcv = sapply(sim_iter, `[[`, "Mn_cv"),
              lambda_gcv = sapply(sim_iter, `[[`, "lambda_cv"),
              crit_gcv = sapply(sim_iter, `[[`, "crit_cv"),
              # finally, the unadjusted estimator
              psi_unadj = sapply(sim_iter, `[[`, "psi_unadj"),
              var_unadj = sapply(sim_iter, `[[`, "var_unadj"),
              ic_unadj = sapply(sim_iter, `[[`, "ic_unadj")
             )
  return(out)
})

# save final output
names(sim_over_nsamp) <- paste0("n", n_samp)
timestamp <- str_replace_all(Sys.time(), " ", "_")
saveRDS(object = sim_over_nsamp,
        file = here("data", args$selector,
                    paste0("dgp_", args$dgp,
                           "_ipw_", args$ipw, "_",
                            args$Q_reg, "_",
                            "trunc_", args$g_trunc, "_",
                            args$Mncv_max, "Mncv_",
                            timestamp, ".rds")))
