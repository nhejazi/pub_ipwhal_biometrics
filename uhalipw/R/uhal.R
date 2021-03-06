utils::globalVariables("..x_names")

#' Undersmoothed Highly Adaptive Lasso for Cross-Validation
#'
#' @param fold Object specifying cross-validation folds as generated by a call
#'  to \code{\link[origami]{make_folds}}.
#' @param data_in An input data set, usually a \code{data.table}.
#' @param outcome_type A \code{character} specifying the error family of the
#'  regression model, passed directly to \code{\link[hal9001]{fit_hal}}.
#' @param x_names A \code{character} giving the name(s) of the column(s) in the
#'  input dataset corresponding to the regressors/predictors.
#' @param y_names A \code{character} giving the name of the column in the input
#'  dataset corresponding to the outcome.
#' @param L1_seq A \code{numeric} sequence of L1-norm values to be matched to
#'  the corresponding values of the regularization parameter lambda.
#' @param lambda_seq A \code{numeric} sequence of values of the regularization
#'  parameter, to be used in fitting HAL regressions and matched up to values
#'  of the L1 norm.
#' @param basis_list A \code{list} of basis functions produced by a prior call
#'  to \code{\link[hal9001]{fit_hal}}, used to ensure that basis functions are
#'  preserved across a sequence of undersmoothed HAL regressions.
#'
#' @import data.table
#' @importFrom origami training validation
#' @importFrom hal9001 fit_hal make_design_matrix
#' @importFrom stats plogis
#'
#' @export
cv_hal_usmooth <- function(fold, data_in, outcome_type, x_names, y_names,
                           L1_seq = seq(0.2, 2000, length.out = 10000),
                           lambda_seq = exp(seq(-0.5, -20, length = 1e4)),
                           basis_list = NULL) {
  ## 0) cross-validation via origami
  train_data <- origami::training(data_in)
  valid_data <- origami::validation(data_in)
  x_train <- as.matrix(train_data[, ..x_names])
  y_train <- as.numeric(train_data[, get(y_names)])
  x_valid <- as.matrix(valid_data[, ..x_names])

  ## 1) fit HAL over sequence of lambda values
  hal_fit <- hal9001::fit_hal(
    X = x_train, Y = y_train,
    max_degree = NULL,
    fit_type = "glmnet",
    family = outcome_type,
    return_lasso = TRUE,
    return_x_basis = TRUE,
    basis_list = basis_list,
    lambda = lambda_seq,
    cv_select = FALSE,
    standardize = FALSE,
    yolo = FALSE
  )
  beta_hat <- hal_fit$glmnet_lasso$beta
  alpha_hat <- hal_fit$glmnet_lasso$a0
  coef_mat <- rbind(alpha_hat, beta_hat)

  ### compute observed L1 norm and map to input sequence
  L1_emp <- apply(beta_hat, 2, function(x) sum(abs(x)))
  L1_ind <- sapply(L1_seq, function(x) min(which(L1_emp > x)))
  L1_ind[L1_ind == Inf] <- length(L1_emp)

  ## 2) predictions on validation data for each value of lambda
  valid_x_basis <- hal9001::make_design_matrix(x_valid, hal_fit$basis_list)
  valid_x_basis_clean <- valid_x_basis[, as.numeric(names(hal_fit$copy_map))]
  pred_mat <- cbind(rep(1, nrow(x_valid)), valid_x_basis_clean)
  preds_valid <- as.matrix(pred_mat %*% coef_mat)

  ## 3) re-scale predictions if propensity score regression
  if (outcome_type == "binomial") {
    preds_valid <- apply(preds_valid, 2, stats::plogis)
  }

  ## 4) output
  out <- list(
    hat_valid = as.matrix(preds_valid[, L1_ind]),
    basis_list = hal_fit$basis_list
  )
  return(out)
}

#' Fit propensity score g_n(A|W) for IPW estimate of treatment-specific mean
#'
#' @param data_in An input data set, usually a \code{data.table}.
#' @param folds A \code{list} of \code{Fold} objects specifying the full set of
#'  cross-validation folds, generated by \code{\link[origami]{make_folds}}.
#' @param x_names A \code{character} giving the name(s) of the column(s) in the
#'  input dataset corresponding to the regressors/predictors. For the treatment
#'  mechanism, this should be only the baseline covariates W.
#' @param y_names A \code{character} giving the name of the column in the input
#'  dataset corresponding to the outcome. For the treatment mechanism, this
#'  should be only the treatment A.
#' @param L1_seq A \code{numeric} sequence of L1-norm values to be matched to
#'  the corresponding values of the regularization parameter lambda.
#' @param lambda_seq A \code{numeric} sequence of values of the regularization
#'  parameter, to be used in fitting HAL regressions and matched up to values
#'  of the L1 norm.
#'
#' @import data.table
#' @importFrom origami cross_validate
#' @importFrom hal9001 fit_hal
#'
#' @export
fit_gn_ipw <- function(data_in, folds, x_names, y_names, L1_seq, lambda_seq) {
  # fit HAL on first training-validation split to get set of basis functions
  gn_cv_fit_hal_initial <- cv_hal_usmooth(
    fold = folds[[1]],
    data_in = data_in,
    outcome_type = "binomial",
    y_names = y_names,
    x_names = x_names,
    L1_seq = L1_seq,
    lambda_seq = lambda_seq,
    basis_list = NULL
  )

  # fit a cross-validated under-smoothed HAL for the propensity score
  # NOTE: use set of basis functions discovered in initial fit on first split
  gn_cv_fit_hal <- origami::cross_validate(
    cv_fun = cv_hal_usmooth,
    folds = folds[-1],
    data = data_in,
    outcome_type = "binomial",
    y_names = y_names,
    x_names = x_names,
    L1_seq = L1_seq,
    lambda_seq = lambda_seq,
    basis_list = gn_cv_fit_hal_initial$basis_list,
    use_future = FALSE,
    .combine = FALSE
  )

  # fit CV-HAL on full data to get reduced basis functions across all fits
  gn_hal_fit_full <- hal9001::fit_hal(
    X = as.matrix(data_in[, ..x_names]),
    Y = as.numeric(data_in[, get(y_names)]),
    max_degree = NULL,
    fit_type = "glmnet",
    family = "binomial",
    return_lasso = FALSE,
    return_x_basis = TRUE,
    cv_select = TRUE,
    standardize = FALSE,
    yolo = FALSE
  )

  # get indices of observations in validation folds
  idx_folds <- do.call(c, lapply(folds, `[[`, "validation_set"))

  # get predictions of propensity score for each L1-norm
  gn_cv <- do.call(rbind, c(list(gn_cv_fit_hal_initial$hat_valid),
                            gn_cv_fit_hal$hat_valid))[idx_folds, ]

  # find CV-risk (negative log-likelihood loss) of propensity score regression
  cv_nloglik <- apply(gn_cv, 2, function(g_hat) {
    nll(obs = data_in[, get(y_names)], probs = g_hat)
  })
  ghat_cv_select <- gn_cv[, which.min(cv_nloglik)]
  L1_cv <- L1_seq[which.min(cv_nloglik)]

  # output
  return(list(
    gn_cv = gn_cv,                           # predictions for each L1-norm
    gn_hal_basis = gn_hal_fit_full$x_basis,  # basis functions from CV-fit
    gn_cv_select = ghat_cv_select,           # predictions for CV-selector
    gn_select_idx = which.min(cv_nloglik),   # index of L1-norm of CV-selector
    L1_cv_select = L1_cv                     # L1-norm minimizer of CV-NLL
  ))
}

#' Outcome regression Q_n(A,W) for treatment-specific mean
#'
#' @param data_in An input data set, usually a \code{data.table}.
#' @param folds A \code{list} of \code{Fold} objects specifying the full set of
#'  cross-validation folds, generated by \code{\link[origami]{make_folds}}.
#' @param x_names A \code{character} giving the name(s) of the column(s) in the
#'  input dataset corresponding to the regressors/predictors. For the outcome
#'  mechanism, this should be the baseline covariates W and treatment A.
#' @param y_names A \code{character} giving the name of the column in the input
#'  dataset corresponding to the outcome. For the outcome regression, this
#'  should be the outcome Y.
#' @param L1_seq A \code{numeric} sequence of L1-norm values to be matched to
#'  the corresponding values of the regularization parameter lambda.
#' @param lambda_seq A \code{numeric} sequence of values of the regularization
#'  parameter, to be used in fitting HAL regressions and matched up to values
#'  of the L1 norm.
#'
#' @import data.table
#' @importFrom origami cross_validate
#'
#' @keywords internal
fit_Qn_ipw <- function(data_in, folds, x_names, y_names, L1_seq, lambda_seq) {
  # fit HAL on first training-validation split to get set of basis functions
  Qn_cv_fit_hal_initial <- cv_hal_usmooth(
    fold = folds[[1]],
    data_in = data_in,
    outcome_type = "gaussian",
    y_names = y_names,
    x_names = x_names,
    L1_seq = L1_seq,
    lambda_seq = lambda_seq,
    basis_list = NULL
  )

  # fit a cross-validated under-smoothed HAL for the outcome regression
  # NOTE: use set of basis functions discovered in initial fit on first split
  Qn_cv_fit_hal <- origami::cross_validate(
    cv_fun = cv_hal_usmooth,
    folds = folds[-1],
    data = data_in,
    outcome_type = "gaussian",
    y_names = y_names,
    x_names = x_names,
    L1_seq = L1_seq,
    lambda_seq = lambda_seq,
    basis_list = Qn_cv_fit_hal_initial$basis_list,
    use_future = FALSE,
    .combine = FALSE
  )

  # get indices of observations in validation folds
  idx_folds <- do.call(c, lapply(folds, `[[`, "validation_set"))
  Qn_cv <- rbind(
    Qn_cv_fit_hal_initial$hat_valid,
    do.call(rbind, Qn_cv_fit_hal$hat_valid)
  )[idx_folds, ]

  # find CV-MSE of outcome regression
  cv_mse <- apply(Qn_cv, 2, function(Q_hat) {
    mse(obs = data_in[, get(y_names)], estim = Q_hat)
  })
  Qhat_cv_select <- Qn_cv[, which.min(cv_mse)]
  L1_cv <- L1_seq[which.min(cv_mse)]

  # output
  return(list(
    Qn_cv = Qn_cv,
    Qn_cv_select = Qhat_cv_select,
    Qn_select_idx = which.min(cv_mse),
    L1_cv_select = L1_cv
  ))
}
