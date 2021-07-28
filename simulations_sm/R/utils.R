###############################################################################
## HELPERS FOR SUMMARIZING SIMULATIONS
###############################################################################
clean_sim_results <- function(sim_results, type) {
  psi_ipw <- sim_results[paste0("psi_", type)]
  var_ic <- sim_results[paste0("var_", type)]
  sim_cleaned <- list(param_est = psi_ipw, param_var = var_ic) %>%
    as_tibble() %>%
    unnest(cols = c(param_est, param_var))
  return(sim_cleaned)
}

prepare_sim_summary <- function(results, n_samp,
                                type = c("usm_trunc", "usm_mintrunc", "gcv")) {
  # iterate over sample sizes and clean each set of results
  cleaned_results <- lapply(seq_along(results), function(iter) {
      out <- clean_sim_results(results[[iter]], type)
      return(out)
    }) %>%
    set_names(paste0("n", n_samp)) %>%
    bind_rows(.id = "n_obs") %>%
    mutate(
      n_obs = as.numeric(str_remove(n_obs, "n"))
    )
  return(cleaned_results)
}

summarize_sim_results <- function(sim_results, psi_true, ci_level = 0.95) {
  ci_norm_bounds <- c(-1, 1) * abs(stats::qnorm(p = (1 - ci_level) / 2))
  sim_summary <- sim_results %>%
    dplyr::filter(!is.na(param_est)) %>%
    mutate(
      ci_lwr = param_est + (ci_norm_bounds[1] * sqrt(param_var)),
      ci_upr = param_est + (ci_norm_bounds[2] * sqrt(param_var)),
      covers = data.table::between(psi_true, ci_lwr, ci_upr)
    ) %>%
    group_by(n_obs) %>%
    summarise(
      bias = abs(mean(param_est - psi_true)),
      n_bias = bias * sqrt(unique(n_obs)),
      mc_var = var(param_est),
      n_mc_var = mc_var * sqrt(unique(n_obs)),
      avg_var = mean(param_var),
      mse = bias^2 + mc_var,
      nmse = mse * unique(n_obs),
      coverage = mean(covers),
      n_sim = n()
    )
  return(sim_summary)
}

###############################################################################
## summarize performance of selector in terms of size of L1 norm, etc.
###############################################################################
summarize_mn <- function(result, type = c("usm_trunc", "usm_mintrunc")) {
  type <- match.arg(type)
  psi_ipw <- result[[paste0("psi_", type)]]
  var_ic <- result[[paste0("var_", type)]]
  Mn_usm <- result[[paste0("Mn_", type)]]
  Mn_gcv <- result[["Mn_gcv"]]
  gtrunc <- result[[paste0("gtrunc_", type)]]
  result_cleaned <- list(param_est = psi_ipw, param_var = var_ic,
                         Mn_usm = Mn_usm, Mn_gcv = Mn_gcv, gtrunc = gtrunc) %>%
    as_tibble() %>%
    unnest(cols = c("param_est", "param_var", "Mn_usm", "Mn_gcv", "gtrunc"))
  return(result_cleaned)
}

summarize_usm_results <- function(results, n_samp,
                                  type = c("usm_trunc", "usm_mintrunc")) {
  uhal_perf <- lapply(seq_along(results), function(iter) {
      sim_summary <- summarize_mn(results[[iter]], type = type)
      return(sim_summary)
    }) %>%
    set_names(paste0("n", n_samp)) %>%
    bind_rows(.id = "n_samp") %>%
    dplyr::filter(!is.na(Mn_usm)) %>%
    mutate(
      n_samp = as.numeric(str_remove(n_samp, "n")),
      Mn_ratio = Mn_usm / Mn_gcv
    ) %>%
    group_by(n_samp) %>%
    mutate(
      mean_Mn = mean(Mn_usm),
      mean_Mn_ratio = mean(Mn_ratio),
      mean_gtrunc = mean(gtrunc),
      n_sim = n()
    )
  return(uhal_perf)
}
