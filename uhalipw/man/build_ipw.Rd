% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npipw.R
\name{build_ipw}
\alias{build_ipw}
\title{Construct undersmoothed IPW estimator for the treatment-specific mean}
\usage{
build_ipw(
  gn_uhal,
  a,
  y,
  Qn_tsm,
  est_type = c("ht", "hajek"),
  bound_gn = NULL,
  gn_dcar = NULL,
  gn_eif = NULL
)
}
\arguments{
\item{gn_uhal}{Undersmoothed estimates of the propensity score gn.}

\item{a}{Observed values of the treatment.}

\item{y}{Observed values of the outcome.}

\item{Qn_tsm}{Counterfactual outcomes, predicted Qn or true Q0, used in the
score equation D_CAR.}

\item{est_type}{Type of IPW estimator, Hajek's stabilized IPW estimator or
the classical Horvitz-Thompson IPW estimator. Hajek's stabilized estimator
has some advantages in cases where the propensity score may suffer from
near-violations of the assumption of positivity. The default is to use the
classical Horvitz-Thompson IPW estimator.}

\item{bound_gn}{A possible lower bound for propensity score truncation.}

\item{gn_dcar}{Estimates of the propensity score used in the re-weighting
term (i.e., Qn / gn) of the EIF. The default is \code{NULL}, in which case
\code{gn} is used in the re-weighting term as well as in the score.}

\item{gn_eif}{Estimates of the propensity score used in the re-weighting
term (i.e., A / gn) of the EIF. The default is \code{NULL}, in which case
\code{gn} is used in the re-weighting term as well as in the IPW estimator.}
}
\description{
While arguments are obvious, note that, by default, the IPW estimator
is computed for the treatment-specific mean E[E[Y | A = 1, W]] under
treatment do(A = 1), while the complementary "control"-specific mean
parameter E[E[Y | A = 0, W]] can be easily computed simply by setting
gn_uhal := 1 - gn_uhal, a := 1 - a, and Qn_tsm := Q(A = 0, W).
}
