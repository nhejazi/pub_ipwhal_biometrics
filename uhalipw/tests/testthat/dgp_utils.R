# baseline covariate, treatment mechanism, outcome mechanism
make_simple_data <- function(n_samp = 1000,
                             type = c("arbitrary", "discrete")) {
  # set default to arbitrary covariate types
  type <- match.arg(type)

  # define treatment mechanism
  g0 <- function(W1, W2) {
    plogis(2 * W1 - W2 - W1 * W2)
  }

  # define outcome mechanism
  Q0 <- function(A, W1, W2) {
    A * (W1 + W2 + W1 * W2) + (1 - A) * W1
  }

  # simulate data
  if (type == "arbitrary") {
    W1 <- runif(n_samp, 0.2, 0.8)
  } else if (type == "discrete") {
    W1 <- rbinom(n_samp, 1, 0.3)
  }
  W2 <- rbinom(n_samp, 1, 0.6)
  A <- rbinom(n_samp, 1, g0(W1, W2))
  Y <- Q0(A, W1, W2) + rnorm(n_samp, 0, 0.1)
  data_obs <- as.data.table(list(W1 = W1, W2 = W2, A = A, Y = Y))
  dgp_funs <- list(g0 = g0, Q0 = Q0)
  return(list(data_obs = data_obs, dgp_funs = dgp_funs))
}

# truth for current data-generating mechanism
get_truth <- function(n_samp = 1e7, tsm_contrast = 1,
                      type = c("arbitrary", "discrete")) {

  # set default to arbitrary covariate types
  type <- match.arg(type)

  # compute large data set from data-generating mechanism
  dgp <- make_simple_data(n_samp, type)
  much_data <- dgp$data_obs

  # extract helper functions
  g0 <- dgp$dgp_funs$g0
  Q0 <- dgp$dgp_funs$Q0

  # compute likelihood factors and truth
  g0_est <- with(much_data, g0(W1, W2))
  Q0_est <- Q0(A = tsm_contrast, W1 = much_data$W1, W2 = much_data$W2)
  true_psi <- mean(Q0_est)

  # use influence function to compute efficiency bound
  A_obs <- much_data$A
  ic_ipw <- ((A_obs / g0_est) * Q0_est) - mean((A_obs / g0_est) * Q0_est)
  effic_bound <- var(ic_ipw)

  # output
  return(list(
    true_psi = true_psi,
    effic_bound = effic_bound,
    g0 = g0_est,
    Q0 = Q0_est
  ))
}
