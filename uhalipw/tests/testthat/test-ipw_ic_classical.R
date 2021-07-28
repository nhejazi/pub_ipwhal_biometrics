context("Check classical IPW estimator using truth and NP-MLE")

# housekeeping and parameters
library(data.table)
library(stringr)
source("dgp_utils.R")
set.seed(72439)
tsm_contrast <- 1
n_obs <- 1000

# generate data for this simulation
sim_truth <- get_truth(n_samp = 1e7, tsm_contrast = tsm_contrast)
dgp <- make_simple_data(n_obs)
data_o <- dgp$data_obs
names_w <- str_subset(colnames(data_o), "W")
names_noty <- colnames(data_o)[!str_detect(colnames(data_o), "Y")]

# also, compute true g_0 and Q_0 from structural equations
g0 <- with(data_o, dgp$dgp_funs$g0(W1 = W1, W2 = W2))
Q0 <- with(data_o, dgp$dgp_funs$Q0(A = tsm_contrast, W1 = W1, W2 = W2))

# for reference, also compute IPW and IC for gn selected by global CV
ipw_g0 <- build_ipw(
  gn_uhal = as.matrix(g0),
  a = data_o$A,
  y = data_o$Y,
  Qn_tsm = Q0,
  est_type = "ht"
)
est_diff <- abs(sim_truth$true_psi - ipw_g0$psi)
true_ic <- abs(mean(((data_o$A / g0) * data_o$Y) - sim_truth$true_psi))
test_that("Classical IPW is AL with known IC for true propensity score", {
  expect_equal(est_diff, true_ic)
})

# also, compute NP-MLE of the propensity score (discrete {W1, W2} only)

# generate data for this simulation
dgp <- make_simple_data(n_obs, type = "discrete")
data_o <- dgp$data_obs
names_w <- str_subset(colnames(data_o), "W")
names_noty <- colnames(data_o)[!str_detect(colnames(data_o), "Y")]
sim_truth <- get_truth(
  n_samp = 1e7, tsm_contrast = tsm_contrast,
  type = "discrete"
)

# also, compute true g_0 and Q_0 from structural equations
g0 <- with(data_o, dgp$dgp_funs$g0(W1 = W1, W2 = W2))
Q0 <- with(data_o, dgp$dgp_funs$Q0(A = tsm_contrast, W1 = W1, W2 = W2))

# compute NP-MLE of gn
g_npmle_fit <- glm(A ~ .^2, data = data_o[, -4], family = "binomial")
g_npmle <- unname(predict(g_npmle_fit, type = "response"))

# for reference, also compute IPW and IC for gn selected by global CV
ipw_gnpmle <- build_ipw(
  gn_uhal = as.matrix(g_npmle),
  a = data_o$A,
  y = data_o$Y,
  Qn_tsm = Q0,
  est_type = "ht"
)
est_diff <- abs(sim_truth$true_psi - ipw_gnpmle$psi)
true_ic <- ((data_o$A / g_npmle) * data_o$Y) - sim_truth$true_psi
dcar_ic <- (Q0 / g_npmle) * (data_o$A - g_npmle)
est_ic <- abs(mean(true_ic - dcar_ic))
test_that("Classical IPW is AL with (IC-D_CAR) for NP-MLE propensity score", {
  expect_equal(est_diff, est_ic)
})
