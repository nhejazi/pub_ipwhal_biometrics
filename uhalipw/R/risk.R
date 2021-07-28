#' Risk under Cross-Entropy Loss
#'
#' @param obs A \code{numeric} vector of observed binary labels.
#' @param probs A \code{numeric} vector of predicted probabilities of class
#'  assignment, of the form P(obs = 1).
#'
#' @keywords internal
nll <- function(obs, probs) {
  nll_loss <- obs * (-log(probs)) + (1 - obs) * (-log(1 - probs))
  risk <- mean(nll_loss)
  return(risk)
}

#' Mean-squared error
#'
#' @param obs A \code{numeric} vector of observations.
#' @param estim A \code{numeric} vector of predictions.
#'
#' @keywords internal
mse <- function(obs, estim) {
  l2_loss <- (estim - obs)^2
  risk <- mean(l2_loss)
  return(risk)
}
