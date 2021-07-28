# housekeeping
library(here)
library(knitr)
library(tidyverse)
library(data.table)
library(ggsci)
library(patchwork)
library(latex2exp)

# utility functions and settings
source(here("R", "utils.R"))
source(here("R", "visualize.R"))
source(here("R", "01_dgp.R"))
pd <- position_dodge(0.3)

# summarization helper
make_ipw_summary <- function(results_file, tsm_contrast = 1, selector) {
  # get dgp ID from results file
  dgp_type <- unlist(str_split(results_file, "_"))[2]

  # parse other simulation information for file names
  modQn_type <- unlist(str_split(results_file, "_"))[6]
  trunc_type <- unlist(str_split(results_file, "_"))[8]
  l1norm_max <- unlist(str_split(results_file, "_"))[9]

  # get truth and load results data
  sim_truth <- get_truth(n_samp = 1e7, tsm_contrast = tsm_contrast,
                         dgp_type = dgp_type)
  psi_true <- sim_truth$true_psi
  eff_bound <- sim_truth$effic_bound

  # read in data
  results <- readRDS(here("data", selector, results_file))
  n_samp <- as.numeric(str_remove(names(results), "n"))

  # summarize data for undersmoothed and truncated estimator
  usm_results <- prepare_sim_summary(results, n_samp, "usm_trunc")
  usm_summary <- summarize_sim_results(usm_results, psi_true)

  # summarize data for undersmoothed estimator without truncation
  mintrunc_usm_results <- prepare_sim_summary(results, n_samp, "usm_mintrunc")
  mintrunc_usm_summary <- summarize_sim_results(mintrunc_usm_results, psi_true)

  # summarize data for estimator relying on global CV
  gcv_results <- prepare_sim_summary(results, n_samp, "gcv")
  gcv_summary <- summarize_sim_results(gcv_results, psi_true)

  # manually extract results for the unadjusted estimator
  unadj_results <- lapply(seq_along(results), function(iter) {
    psi <- results[[iter]]$psi_unadj
    var_psi <- results[[iter]]$var_unadj
    sim_res <- list(param_est = psi, param_var = var_psi) %>%
      as_tibble() %>%
      unnest(cols = c(param_est, param_var))
  }) %>%
  set_names(paste0("n", n_samp)) %>%
  bind_rows(.id = "n_obs") %>%
  mutate(
    n_obs = as.numeric(str_remove(n_obs, "n"))
  )
  unadj_summary <- summarize_sim_results(unadj_results, psi_true)

  # combine (separately) results and summary information for visualization
  sim_results <- list(usm = usm_results,
                      mintrunc = mintrunc_usm_results,
                      gcv = gcv_results,
                      unadj = unadj_results) %>%
    bind_rows(.id = "est_type") %>%
    mutate(
      label = case_when(est_type == "usm" ~ "undersmooth + truncate",
                        est_type == "mintrunc" ~ "undersmooth",
                        est_type == "gcv" ~ "global CV",
                        est_type == "unadj" ~ "unadjusted"),
      est_category = case_when(est_type == "unadj" ~ "naive",
                               est_type %in% c("usm", "mintrunc") ~
                                 "undersmoothed")
    )

  sim_summary <- list(usm = usm_summary,
                      mintrunc = mintrunc_usm_summary,
                      gcv = gcv_summary,
                      unadj = unadj_summary) %>%
    bind_rows(.id = "est_type") %>%
    mutate(
      label = case_when(est_type == "usm" ~ "undersmooth + truncate",
                        est_type == "mintrunc" ~ "undersmooth",
                        est_type == "gcv" ~ "global CV",
                        est_type == "unadj" ~ "unadjusted"),
      est_category = case_when(est_type %in% c("unadj", "gcv") ~ "naive",
                               est_type %in% c("usm", "mintrunc") ~
                                 "undersmoothed")
    )

  ## VISUALIZE SIMULATION METRICS: SAMPLING DISTRIBUTION, BIAS, VARIANCE, MSE
  sim_results_plotted <- make_sim_plots(sim_results, sim_summary, eff_bound)
  p_panel <- (sim_results_plotted$bias | sim_results_plotted$nbias) /
    (sim_results_plotted$nmse | sim_results_plotted$cov)
  p_panel
  ggsave(filename = here("graphs", paste0("dgp", dgp_type),
                         paste(selector, "est", "trunc", trunc_type, "Qn",
                               modQn_type, l1norm_max, "perf_summary.pdf",
                               sep = "_")),
         plot = p_panel, width = 20, height = 16)

  # examine metrics like L1 norm selected by undersmoothing
  usm_perf_summary <- summarize_usm_results(results, n_samp, "usm_trunc")
  mintrunc_usm_perf_summary <- summarize_usm_results(results, n_samp,
                                                     "usm_mintrunc")
  usm_perf_vis <- make_usm_plots(usm_perf_summary)
  p_usm_perf <- usm_perf_vis$mn + usm_perf_vis$mn_ratio + usm_perf_vis$gtrunc
  mintrunc_usm_perf_vis <- make_usm_plots(mintrunc_usm_perf_summary)
  p_mintrunc_usm_perf <- mintrunc_usm_perf_vis$mn +
    mintrunc_usm_perf_vis$mn_ratio

  # save plots of metrics
  p_usm_perf
  ggsave(filename = here("graphs", paste0("dgp", dgp_type),
                         paste(selector, "metrics_trunc_usm", trunc_type,
                               "Qn", modQn_type, l1norm_max, "perf.pdf",
                               sep = "_")),
         plot = p_usm_perf, width = 20, height = 16)
  p_mintrunc_usm_perf
  ggsave(filename = here("graphs", paste0("dgp", dgp_type),
                         paste(selector, "metrics_mintrunc_usm",
                               trunc_type, "Qn", modQn_type, l1norm_max,
                               "perf.pdf", sep = "_")),
         plot = p_mintrunc_usm_perf, width = 20, height = 16)
}

## D_CAR-based selector
dcar_results_files <-
  c("dgp_1a_ipw_ht_Qn_hal_trunc_joint_140Mncv_2020-02-13_08:31:48.rds",
    "dgp_1a_ipw_ht_Qn_hal_trunc_profile_80Mncv_2020-02-13_03:36:43.rds",
    "dgp_1b_ipw_ht_Qn_hal_trunc_joint_350Mncv_2020-02-13_10:34:52.rds",
    "dgp_1b_ipw_ht_Qn_hal_trunc_profile_350Mncv_2020-02-13_04:39:17.rds",
    "dgp_1c_ipw_ht_Qn_hal_trunc_joint_150Mncv_2020-04-02_11:04:28.rds",
    "dgp_1c_ipw_ht_Qn_hal_trunc_profile_100Mncv_2020-04-02_06:02:34.rds",
    "dgp_2a_ipw_ht_Qn_hal_trunc_joint_80Mncv_2020-03-08_03:35:32.rds",
    "dgp_2a_ipw_ht_Qn_hal_trunc_profile_80Mncv_2020-03-07_21:49:10.rds",
    "dgp_2b_ipw_ht_Qn_hal_trunc_joint_120Mncv_2020-03-08_00:43:11.rds",
    "dgp_2b_ipw_ht_Qn_hal_trunc_profile_120Mncv_2020-03-07_20:52:48.rds",
    "dgp_3a_ipw_ht_Qn_hal_trunc_joint_80Mncv_2020-03-08_05:48:57.rds",
    "dgp_3a_ipw_ht_Qn_hal_trunc_profile_80Mncv_2020-03-07_22:57:10.rds",
    "dgp_3b_ipw_ht_Qn_hal_trunc_joint_400Mncv_2020-03-08_06:31:24.rds",
    "dgp_3b_ipw_ht_Qn_hal_trunc_profile_400Mncv_2020-03-07_23:16:43.rds")

# run summarization and save plots
lapply(dcar_results_files, function(f_results) {
  make_ipw_summary(f_results, selector = "dcar")
})

## score-based selector
score_results_files <-
  c("dgp_1a_ipw_ht_Qn_hal_trunc_profile_200Mncv_2020-03-08_01:37:02.rds",
    "dgp_1b_ipw_ht_Qn_hal_trunc_profile_500Mncv_2020-03-08_05:07:30.rds",
    "dgp_1c_ipw_ht_Qn_hal_trunc_profile_200Mncv_2020-04-02_19:44:02.rds",
    "dgp_2a_ipw_ht_Qn_hal_trunc_profile_200Mncv_2020-03-08_11:05:16.rds",
    "dgp_2b_ipw_ht_Qn_hal_trunc_profile_150Mncv_2020-03-08_23:24:53.rds",
    "dgp_3a_ipw_ht_Qn_hal_trunc_profile_200Mncv_2020-03-08_14:55:09.rds",
    "dgp_3b_ipw_ht_Qn_hal_trunc_profile_400Mncv_2020-03-08_16:05:22.rds")

# run summarization and save plots
lapply(score_results_files, function(f_results) {
  make_ipw_summary(f_results, selector = "score")
})

## plateau selector
plateau_results_files <-
  c("dgp_1a_ipw_ht_Qn_hal_trunc_profile_100Mncv_2020-02-15_22:18:42.rds",
    "dgp_1b_ipw_ht_Qn_hal_trunc_profile_350Mncv_2020-02-15_23:19:43.rds",
    "dgp_2a_ipw_ht_Qn_hal_trunc_profile_100Mncv_2020-03-08_15:49:59.rds",
    "dgp_2b_ipw_ht_Qn_hal_trunc_profile_100Mncv_2020-03-08_06:26:09.rds",
    "dgp_3a_ipw_ht_Qn_hal_trunc_profile_100Mncv_2020-03-08_20:47:31.rds",
    "dgp_3b_ipw_ht_Qn_hal_trunc_profile_200Mncv_2020-03-08_22:18:31.rds")

# run summarization and save plots
lapply(plateau_results_files, function(f_results) {
  make_ipw_summary(f_results, selector = "plateau")
})
