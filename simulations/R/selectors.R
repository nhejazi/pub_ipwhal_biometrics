# Selectors for an IPW estimator based on undersmooth HAL for gn

# D_CAR-based selector minimizes the unsolved component of the IC
dcar_selector <- function(ipw_gtrunc, L1_seq, est_type, gtrunc_type,
                          use_eff_bound = TRUE) {
  # compute D_CAR component for all allowed truncations of estimator
  dcar_emp_gtrunc <- sapply(ipw_gtrunc, function(ipw) {
    dcar_emp <- apply(ipw$dcar, 2, function(dcar) abs(mean(dcar)))
    return(dcar_emp)
  })

  # method to choose optimal {delta, lambda} combination
  if (gtrunc_type == "profile") {
    # choose optimal undersmoothing lambda based on no truncation, then
    # subsequently choose a possible delta that minimizes still further
    opt_lambda_idx <- opt_lambda_mintrunc_idx <-
      unname(which.min(dcar_emp_gtrunc[, 1]))
    opt_gtrunc_idx <- unname(which.min(dcar_emp_gtrunc[opt_lambda_idx, ]))
  } else if (gtrunc_type == "joint") {
    # choose optimal undersmoothing lambda and possible truncation delta
    # based on joint minimization to catch any interaction between these
    dcar_emp_minidx <- arrayInd(which.min(dcar_emp_gtrunc),
                                dim(dcar_emp_gtrunc))
    opt_gtrunc_idx <- dcar_emp_minidx[2]
    opt_lambda_idx <- dcar_emp_minidx[1]
    opt_lambda_mintrunc_idx <- unname(which.min(dcar_emp_gtrunc[, 1]))
  }

  # set truncation index to "1" when considering minimally truncated estimator
  if (est_type == "usm_mintrunc") {
    opt_gtrunc_idx <- 1
    opt_lambda_idx <- opt_lambda_mintrunc_idx
  }

  # get D_CAR and estimate for estimator variant under consideration
  dcar_emp <- dcar_emp_gtrunc[opt_lambda_idx, opt_gtrunc_idx]
  ipw <- ipw_gtrunc[[opt_gtrunc_idx]]

  # get IPW estimate that minimizes criterion for D_{CAR}
  psi_ipw <- ipw$psi[opt_lambda_idx]
  ic_ipw <- ipw$ic[, opt_lambda_idx]
  var_ic_ipw <- var(ic_ipw - ipw$dcar[, opt_lambda_idx]) / nrow(ipw$ic)
  #var_ic_ipw <- var(ic_ipw) / nrow(ipw$ic)
  gn_ipw <- ipw$gn[, opt_lambda_idx]

  # compute variance of D_star or D_{CAR} and criterion for sufficiency
  # N.B., instead of IC variance, could simply use the efficiency bound
  if (use_eff_bound) {
    dvar_eff <- sqrt(sim_truth$effic_bound / nrow(ipw$ic))
    dcrit_ratio <- dcar_emp / dvar_eff
  } else {
    dcrit_ratio <- NA
  }

  # find regularization parameter and L1 norm selected by undersmoothing
  lambda_mapped <- lambda_seq[as.numeric(str_remove(colnames(ipw$dcar), "s"))]
  lambda_selected <- lambda_mapped[opt_lambda_idx]
  Mn_selected <- L1_seq[opt_lambda_idx]
  gtrunc_selected <- gtrunc_seq[opt_gtrunc_idx]

  # output
  out <- list(psi = psi_ipw,
              var_ic = var_ic_ipw,
              ic = ic_ipw,
              gn = gn_ipw,
              Mn = Mn_selected,
              lambda = lambda_selected,
              gtrunc = gtrunc_selected,
              dcar_mean = dcar_emp,
              dcar_crit = dcrit_ratio)
  return(out)
}

###############################################################################

## global selector satisfies an optimality criterion that avoids the IC

# from Mark:

# set the L1 norm C so that
# |min_{s,j} P_n phi_{s,j}(A-G_{n,C}) |<= delta sigma_n /(n^{1/2}log n) C^{-1}
# i.e. multiply the cut-off on right hand side in method you implemented by
# delta * sigma_n, where delta is lower bound for g and sigma_n^2 is variance
# of D_{CAR}.
# This might be worth trying since you already ran the other one. If g can get
# small this will undersmooth more.

score_selector <- function(ipw_gtrunc, gn_usm_fit, A_obs, L1_seq, est_type) {
  # NOTE: want to find the lambda/L1-norm whose solutions of this criterion is
  #       closest to the given bound, i.e., many of the choices (possibly) will
  #       satisfy the criterion, but there will be a maximizer.
  # So, the procedure is to find the choice of lambda that has a minimum emp
  # mean over the product of residuals (A - gn) and basis functions that is
  # closest (i.e., maximal) to the bound given by the criterion

  # get sample size and coerce HAL basis to matrix to avoid repeated casting
  n_obs <- length(A_obs)
  hal_basis_mat <- as.matrix(gn_usm_fit$gn_hal_basis)

  # apply global selector across all truncation values in grid
  score_crit_emp_gtrunc <- lapply(ipw_gtrunc, function(ipw) {
    # compute right-hand side of criterion
    delta_gn <- unname(apply(ipw$gn, 2, min))
    sigma_dcar <- sqrt(unname(apply(ipw$dcar, 2, var) / n_obs))
    vdl_cutoffs <- (delta_gn * sigma_dcar) / (sqrt(n_obs) * log(n_obs)) *
      (L1_seq^(-1))

    # compute selector based on global criterion
    score_crit_lambda <- apply(ipw$gn, 2, function(gn_lambda) {
      # should be mean of each column * residual
      emp_mean_basis_resids <- colMeans(hal_basis_mat * (A_obs - gn_lambda))

      # find minimum of residuals
      nz_basis_resids <- emp_mean_basis_resids[emp_mean_basis_resids > 0]
      min_basis_resids <- min(abs(nz_basis_resids))
      return(min_basis_resids)
    })

    # find product of residuals for given lambda closest to criterion
    crit_pass <- score_crit_lambda <= vdl_cutoffs
    score_crit_min <- which.min(abs(score_crit_lambda[crit_pass] -
                                    vdl_cutoffs[crit_pass]))
    score_norm <- min(which(names(score_crit_min) == names(score_crit_lambda)))

    # select estimate, variance, etc., to return
    psi_score <- ipw$psi[score_norm]
    var_score <- ipw$var_ic[score_norm]
    dcar_score <- ipw$dcar[, score_norm]
    gn_score <- ipw$gn[, score_norm]
    ic_score <- ipw$ic[, score_norm]

    # find lambda and L1-norm corresponding to choice by Lepski's method
    lambda_map <- lambda_seq[as.numeric(str_remove(colnames(ipw$dcar), "s"))]
    lambda_score <- lambda_map[score_norm]
    Mn_score <- L1_seq[score_norm]

    # return as simple list
    out <- list(psi_selected = psi_score,
                var_selected = var_score,
                ic_selected = ic_score,
                gn_selected = gn_score,
                Mn_selected = Mn_score,
                lambda_selected = lambda_score,
                dcar_emp = mean(dcar_score)
               )
    return(out)
  })

  # pick best truncation based on minimization of D_CAR
  score_dcar <- sapply(score_crit_emp_gtrunc, function(ipw_score) {
    return(ipw_score$dcar_emp)
  })
  opt_gtrunc_idx <- which.min(abs(score_dcar))

  # set truncation index to "1" when considering minimally truncated estimator
  if (est_type == "usm_mintrunc") {
    opt_gtrunc_idx <- 1
  }

  # pick best truncation and return components
  ipw_score <- score_crit_emp_gtrunc[[opt_gtrunc_idx]]

  # output
  out <- list(psi = ipw_score$psi_selected,
              var_ic = ipw_score$var_selected,
              ic = ipw_score$ic_selected,
              gn = ipw_score$gn_selected,
              Mn = ipw_score$Mn_selected,
              lambda = ipw_score$lambda_selected,
              gtrunc = gtrunc_seq[opt_gtrunc_idx],
              score_crit = ipw_score$dcar_emp)
  return(out)
}

###############################################################################

# plateau selector based on Lepski's method

# comments from Mark:
# This says that we should increase C till either upper or lower bound of
# confidence interval, ie.
# psi_{n,C}-1.96 sigma_n(C) or psi_{n,C}+1.96 sigma_n(C) reaches a plateau.
# Since d/dC sigma_n(C)>0, if psi_{n,C} is increasing, then we look for
# d/dC psi_{n,C}=1.96 d/dC sigma_n(C), while if psi_{n,C} is decreasing, then
# we look for d/dC psi_{n,C}=-1.96 d/dC sigma_n(C). It could be that for a
# while psi_{n,C} is increasing, and then goes down etc, so depending on what
# it is doing we look at the relevant upper or lower bound for a plateau.

plateau_selector <- function(ipw_gtrunc, gn_usm_fit, L1_seq, est_type,
                             ci_level = 0.95) {
  # make confidence interval multipliers
  ci_norm_mult <- abs(stats::qnorm(p = (1 - ci_level) / 2))

  # apply global selector across all truncation values in grid
  lepski_gtrunc <- lapply(ipw_gtrunc, function(ipw) {
    # find difference in estimate and CIs across sequence of lambdas
    psi_diff <- diff(ipw$psi)
    ci_psi_diff <- sign(psi_diff) * ci_norm_mult * diff(sqrt(ipw$var_ic))

    # Lepski's method looks for the first point when the change in the
    # parameter estimate is ~less than change in the confidence interval width
    valid_norms <- which(psi_diff <= ci_psi_diff)
    lepski_norm <- which.min(abs(psi_diff - ci_psi_diff))

    # select estimate, variance, etc., to return
    psi_lepski <- ipw$psi[lepski_norm]
    var_lepski <- ipw$var_ic[lepski_norm]
    dcar_lepski <- ipw$dcar[, lepski_norm]
    gn_lepski <- ipw$gn[, lepski_norm]
    ic_lepski <- ipw$ic[, lepski_norm]

    # find lambda and L1-norm corresponding to choice by Lepski's method
    lambda_map <- lambda_seq[as.numeric(str_remove(colnames(ipw$dcar), "s"))]
    lambda_lepski <- lambda_map[lepski_norm]
    Mn_lepski <- L1_seq[lepski_norm]

    # return as simple list
    out <- list(psi_selected = psi_lepski,
                var_selected = var_lepski,
                ic_selected = ic_lepski,
                gn_selected = gn_lepski,
                Mn_selected = Mn_lepski,
                lambda_selected = lambda_lepski,
                dcar_emp = mean(dcar_lepski)
               )
    return(out)
  })

  # pick best truncation based on minimization of D_CAR
  lepski_dcar <- sapply(lepski_gtrunc, function(ipw_lepski) {
    return(ipw_lepski$dcar_emp)
  })
  opt_gtrunc_idx <- which.min(abs(lepski_dcar))

  # set truncation index to "1" when considering minimally truncated estimator
  if (est_type == "usm_mintrunc") {
    opt_gtrunc_idx <- 1
  }

  # pick best truncation and return components
  ipw_lepski <- lepski_gtrunc[[opt_gtrunc_idx]]

  # output
  out <- list(psi = ipw_lepski$psi_selected,
              var_ic = ipw_lepski$var_selected,
              ic = ipw_lepski$ic_selected,
              gn = ipw_lepski$gn_selected,
              Mn = ipw_lepski$Mn_selected,
              lambda = ipw_lepski$lambda_selected,
              gtrunc = gtrunc_seq[opt_gtrunc_idx],
              plateau_crit = ipw_lepski$dcar_emp)
  return(out)
}
