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
pd <- position_dodge(0.4)

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

  # load results form Ashkan's simulations
  ae_results_raw <- readRDS(here("data", "ashkan", "supp_results.rds")) %>%
    as.data.table()
  ae_results_nontrunc <-
    readRDS(here("data", "ashkan", "supp_results_nontrunc.rds")) %>%
    as.data.table() %>%
    dplyr::select(str_subset(colnames(.), "Score")) %>%
    setnames(colnames(.), paste(colnames(.), "nontrunc", sep = "_"))
  ae_results_all <- cbind(ae_results_raw, ae_results_nontrunc)

  # read in data
  results <- readRDS(here("data", selector, results_file))
  n_samp <- as.numeric(str_remove(names(results), "n"))

  # summarize data for undersmoothed and truncated estimator
  usm_results <- prepare_sim_summary(results, n_samp, "usm_trunc")
  usm_summary <- summarize_sim_results(usm_results, psi_true)

  # summarize data for undersmoothed estimator without truncation
  mintrunc_usm_results <- prepare_sim_summary(results, n_samp,
                                              "usm_mintrunc")
  mintrunc_usm_summary <- summarize_sim_results(mintrunc_usm_results,
                                                psi_true)

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

  # manually create results for Ashkan's score-based estimator
  ae_summary_clean <- ae_results_all %>%
    mutate(sim = as.character(sim),
           sim = case_when(sim == "1" ~ "1b",
                           sim == "2" ~ "1a",
                           sim == "3" ~ "2a")
          ) %>%
    setDT() %>%
    melt(
     measure.vars = list(c(1, 2, 12), c(3, 4, 13), c(5, 6, 14), c(10, 11, 15)),
     value.name = c("est", "rn_mse", "coverage", "rnbias"),
     variable.name = "est_type",
     variable.factor = FALSE
    ) %>%
    mutate(
      est_type = case_when(est_type == "1" ~ "dcar",
                           est_type == "2" ~ "score_trunc",
                           est_type == "3" ~ "score_nontrunc")
    ) %>%
    rename(n_obs = n, sim_type = sim, n_bias = rnbias) %>%
    mutate(
      bias = n_bias / sqrt(n_obs),
      mse = rn_mse / sqrt(n_obs),
      nmse = n_obs * (rn_mse / sqrt(n_obs))
    ) %>%
    dplyr::filter(sim_type == dgp_type) %>%
    select(-sim_type, -true, -est, -rn_mse) %>%
    as_tibble()

  # combine summary information for visualization
  sim_summary <- list(usm = usm_summary,
                      mintrunc = mintrunc_usm_summary,
                      gcv = gcv_summary,
                      unadj = unadj_summary) %>%
    bind_rows(.id = "est_type") %>%
    mutate(
      label = case_when(est_type == "usm" ~ "undersmoothed + truncation",
                        est_type == "mintrunc" ~ "undersmoothed only",
                        est_type == "gcv" ~ "global cross-validation",
                        est_type == "unadj" ~ "unadjusted (E[Y|A=1])"),
      est_category = case_when(est_type %in% c("unadj", "gcv") ~ "naive",
                               est_type %in% c("usm", "mintrunc") ~
                                 "undersmoothed")
    )

  # join simulation summary with Ashkan's summary for visualizations
  supp_sim_summary_joined <- sim_summary %>%
    select(-n_sim, -label, -est_category, -mc_var, -n_mc_var, -avg_var) %>%
    full_join(ae_summary_clean) %>%
    mutate(
      label = case_when(est_type == "usm" ~ "DCAR (truncated)",
                        est_type == "mintrunc" ~ "DCAR",
                        est_type == "gcv" ~ "CV",
                        est_type == "unadj" ~ "Unadjusted",
                        est_type == "dcar" ~ "DCAR2",
                        est_type == "score_nontrunc" ~ "Score",
                        est_type == "score_trunc" ~ "Score (truncated)")
    ) %>%
    dplyr::filter(label != "DCAR2")

  ## set axis limits depending on DGP type
  ylim_nbias <- list("1a" = c(0, 5),
                     "1b" = c(0, 0.5),
                     "2a" = c(0, 2))
  ylim_nmse <- list("1a" = c(0, 10),
                    "1b" = c(0, 2),
                    "2a" = c(0, 3.5))

  ## VISUALIZE SIMULATION METRICS: BIAS, SCALED BIAS, SCALED MSE, COVERAGE
  # make plot of absolute bias
  p_bias <- supp_sim_summary_joined %>%
    ggplot(aes(x = as.factor(n_obs), y = bias, group = label,
               shape = label, fill = label)) +
    geom_line(linetype = "dashed", size = 0.5, position = pd) +
    geom_point(alpha = 0.7, size = 8, position = pd) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
    labs(
      x = "",
      y = TeX("|$\\psi - \\hat{\\psi}$|"),
      title = "Raw bias"
    ) +
    theme_bw() +
    theme(legend.background = element_rect(fill = "gray90", size = 0.25,
                                           linetype = "dotted"),
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 28),
          text = element_text(size = 32),
          axis.text.x = element_text(angle = 20, colour = "black",
                                     size = 30, hjust = 1),
          axis.text.y = element_text(colour = "black", size = 30)
         ) +
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    scale_shape_manual(values = c(21, 24, 25, 23, 22, 4)) +
    scale_fill_lancet()

  # make plot of n-scaled variance
  p_nbias <- supp_sim_summary_joined %>%
    ggplot(aes(x = as.factor(n_obs), y = n_bias, group = label,
               shape = label, fill = label)) +
    geom_line(linetype = "dashed", size = 0.5, position = pd) +
    geom_point(alpha = 0.7, size = 8, position = pd) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
    labs(
      x = "Sample size",
      y = TeX("$\\sqrt{n} \\times$ |$\\psi - \\hat{\\psi}$|"),
      title = "Scaled bias"
    ) +
    theme_bw() +
    theme(legend.background = element_rect(fill = "gray90", size = 0.25,
                                           linetype = "dotted"),
          legend.position = "none",
          legend.title = element_blank(),
          text = element_text(size = 32),
          axis.text.x = element_text(angle = 20, colour = "black",
                                     size = 30, hjust = 1),
          axis.text.y = element_text(colour = "black", size = 30)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    scale_shape_manual(values = c(21, 24, 25, 23, 22, 4)) +
    scale_fill_lancet() +
    coord_cartesian(ylim = ylim_nbias[[dgp_type]])

  # make plot of scaled MSE
  p_nmse <- supp_sim_summary_joined %>%
    ggplot(aes(x = as.factor(n_obs), y = nmse, group = label,
               shape = label, fill = label)) +
    geom_line(linetype = "dashed", size = 0.5, position = pd) +
    geom_point(alpha = 0.7, size = 8, position = pd) +
    geom_hline(yintercept = sim_truth$effic_bound, linetype = "dashed",
               colour = "black") +
    labs(
      x = "",
      y = TeX("$n \\times$ ($(\\psi - \\hat{\\psi})^2 + \\hat{\\sigma^2})$"),
      title = "Scaled mean squared error"
    ) +
    theme_bw() +
    theme(legend.background = element_rect(fill = "gray90", size = 0.25,
                                           linetype = "dotted"),
          legend.position = "none",
          legend.title = element_blank(),
          text = element_text(size = 32),
          axis.text.x = element_text(angle = 20, colour = "black",
                                     size = 30, hjust = 1),
          axis.text.y = element_text(colour = "black", size = 30)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    scale_shape_manual(values = c(21, 24, 25, 23, 22, 4)) +
    scale_fill_lancet() +
    coord_cartesian(ylim = ylim_nmse[[dgp_type]])

  # make plot of CI coverage
  p_cov <- supp_sim_summary_joined %>%
    ggplot(aes(x = as.factor(n_obs), y = coverage, group = label,
               shape = label, fill = label)) +
    geom_line(linetype = "dashed", size = 0.5, position = pd) +
    geom_point(alpha = 0.7, size = 8, position = pd) +
    geom_hline(yintercept = 0.95, linetype = "dashed",
               colour = "black") +
    labs(
      x = "",
      y = "Coverage",
      title = "Confidence interval coverage"
    ) +
    theme_bw() +
    theme(legend.background = element_rect(fill = "gray90", size = 0.25,
                                             linetype = "dotted"),
          legend.position = "none",
          legend.title = element_blank(),
          text = element_text(size = 32),
          axis.text.x = element_text(angle = 20, colour = "black",
                                     size = 30, hjust = 1),
          axis.text.y = element_text(colour = "black", size = 30)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    scale_shape_manual(values = c(21, 24, 25, 23, 22, 4)) +
    scale_fill_lancet() +
    ylim(0, 1)

  # create panel plot
  p_panel <- (p_bias + p_nbias) / (p_nmse + p_cov)
  ggsave(filename = here("graphs", "manuscript",
                         paste("supp",
                               paste0("dgp", dgp_type),
                               "comparison.pdf",
                               sep = "_")),
         plot = p_panel, width = 24, height = 16)
}

## results across all selectors
results_files <-
  c(#"dgp_1a_ipw_ht_Qn_hal_trunc_joint_140Mncv_2020-02-13_08:31:48.rds",
    "dgp_1a_ipw_ht_Qn_hal_trunc_profile_80Mncv_2020-02-13_03:36:43.rds",
    "dgp_1b_ipw_ht_Qn_hal_trunc_joint_350Mncv_2020-02-13_10:34:52.rds",
    #"dgp_1b_ipw_ht_Qn_hal_trunc_profile_350Mncv_2020-02-13_04:39:17.rds",
    #"dgp_2a_ipw_ht_Qn_hal_trunc_joint_80Mncv_2020-03-08_03:35:32.rds",
    "dgp_2a_ipw_ht_Qn_hal_trunc_profile_80Mncv_2020-03-07_21:49:10.rds")

# run summarization and save plots
lapply(results_files, function(f_results) {
  make_ipw_summary(f_results, selector = "dcar")
})
