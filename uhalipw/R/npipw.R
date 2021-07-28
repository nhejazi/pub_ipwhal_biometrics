#' Construct undersmoothed IPW estimator for the treatment-specific mean
#'
#' While arguments are obvious, note that, by default, the IPW estimator
#' is computed for the treatment-specific mean E[E[Y | A = 1, W]] under
#' treatment do(A = 1), while the complementary "control"-specific mean
#' parameter E[E[Y | A = 0, W]] can be easily computed simply by setting
#' gn_uhal := 1 - gn_uhal, a := 1 - a, and Qn_tsm := Q(A = 0, W).
#'
#' @param gn_uhal Undersmoothed estimates of the propensity score gn.
#' @param a Observed values of the treatment.
#' @param y Observed values of the outcome.
#' @param Qn_tsm Counterfactual outcomes, predicted Qn or true Q0, used in the
#'  score equation D_CAR.
#' @param est_type Type of IPW estimator, Hajek's stabilized IPW estimator or
#'  the classical Horvitz-Thompson IPW estimator. Hajek's stabilized estimator
#'  has some advantages in cases where the propensity score may suffer from
#'  near-violations of the assumption of positivity. The default is to use the
#'  classical Horvitz-Thompson IPW estimator.
#' @param bound_gn A possible lower bound for propensity score truncation.
#' @param gn_dcar Estimates of the propensity score used in the re-weighting
#'  term (i.e., Qn / gn) of the EIF. The default is \code{NULL}, in which case
#'  \code{gn} is used in the re-weighting term as well as in the score.
#' @param gn_eif Estimates of the propensity score used in the re-weighting
#'  term (i.e., A / gn) of the EIF. The default is \code{NULL}, in which case
#'  \code{gn} is used in the re-weighting term as well as in the IPW estimator.
#'
#' @importFrom stats var
#'
#' @export
build_ipw <- function(gn_uhal,
                      a,
                      y,
                      Qn_tsm,
                      est_type = c("ht", "hajek"),
                      bound_gn = NULL,
                      gn_dcar = NULL,
                      gn_eif = NULL) {

  # set default to be Horvitz-Thompson (non-stabilized) IPW, not Hajek's
  est_type <- match.arg(est_type)

  # set estimated gn to respect bounds from true propensity score g0
  if (!is.null(bound_gn)) {
    gn_uhal[gn_uhal < min(bound_gn)] <- min(bound_gn)
  }

  # compute appropriate IPW estimator
  psi_ipw <- unname(apply(gn_uhal, 2, function(gn) {
    psi_est <- ipw_est(a = a, y = y, gn = gn, est_type = est_type)
    return(psi_est)
  }))

  # compute projection onto D_CAR for re-weighted estimator
  dcar_ipw <- apply(gn_uhal, 2, function(gn) {
    if (!is.null(gn_dcar)) {
      dcar_est <- ipw_dcar(a = a, y = y, gn = gn, gn_dcar = gn_dcar,
                           Qn_tsm = Qn_tsm, est_type = est_type)
    } else {
      dcar_est <- ipw_dcar(a = a, y = y, gn = gn, Qn_tsm = Qn_tsm,
                           est_type = est_type)
    }
    return(dcar_est)
  })

  # compute influence function for re-weighted estimator and scaled variance
  eif_ipw <- apply(gn_uhal, 2, function(gn) {
    if (!is.null(gn_eif)) {
      eif_est <- ipw_eif(a = a, y = y, gn = gn, gn_eif = gn_eif,
                         est_type = est_type)
    } else {
      eif_est <- ipw_eif(a = a, y = y, gn = gn, est_type = est_type)
    }
    return(eif_est)
  })

  # estimate variance directly from the EIF
  var_eif_ipw <- unname(apply(eif_ipw, 2, stats::var) / length(y))

  # output
  return(list(
    psi = psi_ipw,
    dcar = dcar_ipw,
    ic = eif_ipw,
    var_ic = var_eif_ipw,
    gn = gn_uhal
  ))
}

###############################################################################

#' Compute IPW Estimators
#'
#' @param a Observed values of the treatment.
#' @param y Observed values of the outcome.
#' @param gn_uhal Undersmoothed estimates of the propensity score gn.
#' @param est_type Type of IPW estimator, Hajek's stabilized IPW estimator or
#'  the classical Horvitz-Thompson IPW estimator. Hajek's stabilized estimator
#'  has some advantages in cases where the propensity score may suffer from
#'  near-violations of the assumption of positivity. The default is to use the
#'  classical Horvitz-Thompson IPW estimator.
#'
#' @keywords internal
ipw_est <- function(a, y, gn, est_type = c("ht", "hajek")) {
  # the Horvitz-Thompson IPW estimator with _non-stabilized_ weights
  if (est_type == "ht") {
    psi <- mean((a / gn) * y)

  # Hajek's IPW estimator with _stabilized_ weights
  } else if (est_type == "hajek") {
    psi <- mean((a / gn) * y) / mean(a / gn)
  }
  return(psi)
}

###############################################################################

#' Compute the Efficient Influence Function of IPW Estimators
#'
#' @param a Observed values of the treatment.
#' @param y Observed values of the outcome.
#' @param gn Estimates of the propensity score gn used in computing the IPW
#'  estimator. These same estimates may be reused in computing the EIF or be
#'  replaced by the alternative \code{gn_eif}.
#' @param gn_eif Estimates of the propensity score used in the re-weighting
#'  term (i.e., A / gn) of the EIF. The default is \code{NULL}, in which case
#'  \code{gn} is used in the re-weighting term as well as in the IPW estimator.
#' @param est_type Type of IPW estimator, Hajek's stabilized IPW estimator or
#'  the classical Horvitz-Thompson IPW estimator. Hajek's stabilized estimator
#'  has some advantages in cases where the propensity score may suffer from
#'  near-violations of the assumption of positivity. The default is to use the
#'  classical Horvitz-Thompson IPW estimator.
#'
#' @keywords internal
ipw_eif <- function(a, y, gn, gn_eif = NULL, est_type = c("ht", "hajek")) {
  # compute appropriate IPW estimator
  psi <- ipw_est(a = a, y = y, gn = gn, est_type = est_type)

  # EIF for the Horvitz-Thompson IPW estimator with _non-stabilized_ weights
  if (est_type == "ht") {
    # optionally, use a different gn in the EIF than the undersmoothed gn
    if (!is.null(gn_eif)) {
      eif <- ((a / gn_eif) * y) - psi
    } else {
      eif <- ((a / gn) * y) - psi
    }

  # EIF for Hajek's IPW estimator with _stabilized_ weights
  } else if (est_type == "hajek") {
    # optionally, use a different gn in the EIF than the undersmoothed gn
    if (!is.null(gn_eif)) {
      eif <- ((a / gn_eif) * (y - psi))
    } else {
      eif <- ((a / gn) * (y - psi))
    }
  }
  return(eif)
}

###############################################################################

#' Compute the Coarsened-at-Random Projection of IPW Estimators
#'
#' @param a Observed values of the treatment.
#' @param y Observed values of the outcome.
#' @param gn Estimates of the propensity score gn used in computing the score
#'  (i.e., (a - gn)) for the treatment mechanism. These same estimates may be
#'  reused in computing the EIF or be replaced instead by \code{gn_dcar}.
#' @param gn_dcar Estimates of the propensity score used in the re-weighting
#'  term (i.e., Qn / gn) of the EIF. The default is \code{NULL}, in which case
#'  \code{gn} is used in the re-weighting term as well as in the score.
#' @param Qn_tsm Counterfactual outcomes, predicted Qn or true Q0, used in the
#'  score equation D_CAR.
#' @param est_type Type of IPW estimator, Hajek's stabilized IPW estimator or
#'  the classical Horvitz-Thompson IPW estimator. Hajek's stabilized estimator
#'  has some advantages in cases where the propensity score may suffer from
#'  near-violations of the assumption of positivity. The default is to use the
#'  classical Horvitz-Thompson IPW estimator.
#'
#' @keywords internal
ipw_dcar <- function(a, y, gn, gn_dcar = NULL, Qn_tsm,
                     est_type = c("ht", "hajek")) {
  # D_CAR for the Horvitz-Thompson IPW estimator with _non-stabilized_ weights
  if (est_type == "ht") {
    if (!is.null(gn_dcar)) {
      dcar <- ((Qn_tsm / gn_dcar) * (a - gn))
    } else {
      dcar <- ((Qn_tsm / gn) * (a - gn))
    }

  # D_CAR for Hajek's IPW estimator with _stabilized_ weights
  } else if (est_type == "hajek") {
    # compute Hajek's IPW estimator
    psi <- ipw_est(a = a, y = y, gn = gn, est_type = est_type)

    # compute D_CAR, which now requires the IPW estimator's realization
    if (!is.null(gn_dcar)) {
      dcar <- ((Qn_tsm - psi) / gn_dcar) * (a - gn)
    } else {
      dcar <- ((Qn_tsm - psi) / gn) * (a - gn)
    }
  }
  return(dcar)
}
