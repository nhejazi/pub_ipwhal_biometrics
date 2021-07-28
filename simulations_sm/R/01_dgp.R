make_simple_data <- function(n_samp = 1000, dgp_type) {
  # Simulation helper for generating data under a variety of data-generating
  # processes, summarized in 3 scenarios. In each case, the treatment mechanism
  # is defined as a function of only two baseline covariates {W1, W2} while the
  # outcome mechanism is defined as a function of {A, W1, W2}. Some scenarios
  # contain extra covariates that impact neither of these mechanisms.

  library(data.table)
  #############################################################################
  # SCENARIO 1: simple DGPs with no positivity issues whatsoever
  #############################################################################

  # 1a: relatively simple DGP with no positivity issues
  if (dgp_type == "1a") {
    # define treatment mechanism
    g0 <- function(W1, W2) {
      plogis(2 * W1 - W2 - W1 * W2)
    }

    # define outcome mechanism
    Q0 <- function(A, W1, W2) {
      (A * (W1 + W2 + W1 * W2)) + ((1 - A) * W1)
    }

    # simulate data
    W1 <- runif(n_samp, 0.2, 0.8)
    W2 <- rbinom(n_samp, 1, 0.6)
    A <- rbinom(n_samp, 1, g0(W1, W2))
    Y <- Q0(A, W1, W2) + rnorm(n_samp, 0, 0.1)
    data_obs <- as.data.table(list(W1 = W1, W2 = W2, A = A, Y = Y))
  }

  # 1b: simple DGP with randomization of treatment
  if (dgp_type == "1b") {
    # define treatment mechanism
    g0 <- function(W1, W2) {
      # randomized, so just ignores covariate information
      rep(0.5, length(W1))
    }

    # define outcome mechanism
    Q0 <- function(A, W1, W2) {
      (A * (W1 + W2 + W1 * W2)) + ((1 - A) * W1)
    }

    # simulate data
    W1 <- runif(n_samp, 0.2, 0.8)
    W2 <- rbinom(n_samp, 1, 0.3)
    A <- rbinom(n_samp, 1, g0(W1, W2))
    Y <- Q0(A, W1, W2) + rnorm(n_samp, 0, 0.1)
    data_obs <- as.data.table(list(W1 = W1, W2 = W2, A = A, Y = Y))
  }

  # 1c: relatively simple DGP with outcome unaffected by exposure
  if (dgp_type == "1c") {
    # define treatment mechanism
    g0 <- function(W1, W2) {
      plogis(2 * W1 - W2 - W1 * W2)
    }

    # define outcome mechanism
    Q0 <- function(A, W1, W2) {
      W1 + W2 + W1 * W2
    }

    # simulate data
    W1 <- runif(n_samp, 0.1, 1.2)
    W2 <- rbinom(n_samp, 1, 0.6)
    A <- rbinom(n_samp, 1, g0(W1, W2))
    Y <- Q0(A, W1, W2) + rnorm(n_samp, 0, 0.1)
    data_obs <- as.data.table(list(W1 = W1, W2 = W2, A = A, Y = Y))
  }

  #############################################################################
  # SCENARIO 2: simple DGPs with moderate to severe positivity issues
  #############################################################################

  # 2a: DGP with severe/significant positivity issues (NOTE: min(g0) ~= 0.017)
  if (dgp_type == "2a") {
    # define treatment mechanism
    g0 <- function(W1, W2) {
      plogis((2 * W1) - (4 * W2) + (W1 * W2))
    }

    # define outcome mechanism
    Q0 <- function(A, W1, W2) {
      (A * (W1 + W2)) + ((1 - A) * W1)
    }

    # simulate data
    W1 <- runif(n_samp, 0, 0.6)
    W2 <- rbinom(n_samp, 1, 0.05)
    A <- rbinom(n_samp, 1, g0(W1, W2))
    Y <- Q0(A, W1, W2) + rnorm(n_samp, 0, 0.1)
    data_obs <- as.data.table(list(W1 = W1, W2 = W2, A = A, Y = Y))
  }

  # 2b: DGP with moderate positivity issues (NOTE: min(g0) ~= 0.0474)
  if (dgp_type == "2b") {
    # define treatment mechanism
    g0 <- function(W1, W2) {
      plogis((4 * W1) - (3 * W2) - (W1 * W2) - (5 * W1^2))
    }

    # define outcome mechanism
    Q0 <- function(A, W1, W2) {
      (A * (W1 - W2)) + ((1 - A) * W1)
    }

    # simulate data
    W1 <- runif(n_samp, 0, 0.6)
    W2 <- rbinom(n_samp, 1, 0.4)
    A <- rbinom(n_samp, 1, g0(W1, W2))
    Y <- Q0(A, W1, W2) + rnorm(n_samp, 0, 0.1)
    data_obs <- as.data.table(list(W1 = W1, W2 = W2, A = A, Y = Y))
  }

  #############################################################################
  # SCENARIO 3: miscellaneous DGPs (e.g., contributed by co-authors)
  #############################################################################

  # DGP from Ashkan's "scenario 1"
  if (dgp_type == "3a") {
    # define treatment mechanism
    g0 <- function(W1, W2) {
      plogis((W1 + W2) / 2)
    }

    # define outcome mechanism
    Q0 <- function(A, W1, W2) {
      W1 + W2
    }

    # simulate data
    W1 <- rnorm(n_samp, 0, 0.5)
    W2 <- rbinom(n_samp, 1, 0.5)
    A <- rbinom(n_samp, 1, g0(W1, W2))
    Y <- Q0(A, W1, W2) + rnorm(n_samp, 0, 0.1)
    data_obs <- as.data.table(list(W1 = W1, W2 = W2, A = A, Y = Y))
  }

  # DGP from Ashkan's "scenario 2"
  if (dgp_type == "3b") {
    # define treatment mechanism
    g0 <- function(W1, W2) {
      plogis((W2^2 - exp(W1 / 2)) / 2)
    }

    # define outcome mechanism
    Q0 <- function(A, W1, W2) {
      -(2 * W2^2) + (2 * W1) + W2 + (W1 * W2) + 0.5
    }

    # simulate data
    W1 <- rnorm(n_samp, 0, 0.5)
    W2 <- rbinom(n_samp, 1, 0.5)
    A <- rbinom(n_samp, 1, g0(W1, W2))
    Y <- Q0(A, W1, W2) + rnorm(n_samp, 0, 0.1)
    data_obs <- as.data.table(list(W1 = W1, W2 = W2, A = A, Y = Y))
  }

  # output
  dgp_funs <- list(g0 = g0, Q0 = Q0)
  return(list(data_obs = data_obs, dgp_funs = dgp_funs))
}

# truth for current data-generating mechanism
# NOTE: g0 and Q0 are functions of {A, W1, W2} only, regardless of scenario
get_truth <- function(n_samp = 1e7, tsm_contrast = 1, dgp_type) {
  # compute large data set from data-generating mechanism
  dgp <- make_simple_data(n_samp = n_samp, dgp_type = dgp_type)
  much_data <- dgp$data_obs

  # extract helper functions
  g0 <- dgp$dgp_funs$g0
  Q0 <- dgp$dgp_funs$Q0

  # extract data components
  A_obs <- much_data$A
  Y_obs <- much_data$Y

  # compute likelihood factors
  g0_est <- with(much_data, g0(W1, W2))
  Qnat_est <- Q0(A = A_obs, W1 = much_data$W1, W2 = much_data$W2)
  Qtsm_est <- Q0(A = tsm_contrast, W1 = much_data$W1, W2 = much_data$W2)

  # compute truth and use influence function to compute efficiency bound
  true_psi <- mean(Qtsm_est)
  eif_ipw <- ((A_obs / g0_est) * (Y_obs - Qnat_est)) + Qtsm_est - true_psi
  effic_bound <- var(eif_ipw)

  # output
  return(list(true_psi = true_psi,
              effic_bound = effic_bound,
              g0 = g0_est,
              Q0 = Qtsm_est))
}
