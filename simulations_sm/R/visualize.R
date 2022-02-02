make_sim_plots <- function(sim_results, sim_summary, eff_var) {
  # plot sampling distribution of estimators
  p_samp_dist <- sim_results %>%
    group_by(est_type, n_obs) %>%
    mutate(
      mean_est = mean(param_est, na.rm = TRUE)
    ) %>%
    ungroup(n_obs) %>%
    ggplot(aes(x = param_est, fill = est_type)) +
      geom_histogram(alpha = 0.75, binwidth = 0.01) +
      geom_vline(aes(xintercept = mean_est), linetype = "dashed",
                 colour = "black") +
      geom_vline(aes(xintercept = psi_true), linetype = "dotted",
                 colour = "blue") +
      facet_grid(cols = vars(n_obs), rows = vars(label),
                 scales = "free") +
      labs(
        x = "",
        y = "",
        title = "Sampling distrib. of CV-IPW"
      ) +
      theme_bw() +
      theme(legend.position = "none",
            text = element_text(size = 25),
            axis.text.x = element_text(colour = "black", size = 22),
            axis.text.y = element_text(colour = "black", size = 22))

  # make plot of absolute bias
  p_bias <- sim_summary %>%
    ggplot(aes(x = as.factor(n_obs), y = bias, group = label, shape = label,
               fill = label)) +
    geom_point(alpha = 0.7, size = 8, position = pd) +
    geom_line(linetype = "dashed", size = 0.5, position = pd) +
    geom_hline(aes(yintercept = 0), linetype = "dashed", colour = "black") +
    labs(
      x = "",
      y = TeX("|$\\psi - \\hat{\\psi}$|"),
      title = "Bias"
    ) +
    #facet_grid(vars(est_category), scales = "free_y") +
    theme_bw() +
    theme(legend.background = element_rect(fill = "gray90", size = 0.25,
                                           linetype = "dotted"),
          legend.position = "bottom",
          legend.title = element_blank(),
          text = element_text(size = 30),
          axis.text.x = element_text(angle = 20, colour = "black", size = 26),
          axis.text.y = element_text(colour = "black", size = 26)
         ) +
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    scale_shape_manual(values = c(22, 21, 24, 25)) +
    scale_fill_nejm()

  # make plot of n-scaled variance
  p_nbias <- sim_summary %>%
    ggplot(aes(x = as.factor(n_obs), y = n_bias, group = label, shape = label,
               fill = label)) +
    geom_point(alpha = 0.7, size = 8, position = pd) +
    geom_line(linetype = "dashed", size = 0.5, position = pd) +
    geom_hline(aes(yintercept = 0), linetype = "dashed", colour = "black") +
    labs(
      x = "",
      y = TeX("$\\sqrt{n} \\times$ |$\\psi - \\hat{\\psi}$|"),
      title = "Scaled bias"
    ) +
    #facet_grid(vars(est_category), scales = "free_y") +
    theme_bw() +
    theme(legend.background = element_rect(fill = "gray90", size = 0.25,
                                             linetype = "dotted"),
          legend.position = "none",
          legend.title = element_blank(),
          text = element_text(size = 30),
          axis.text.x = element_text(angle = 20, colour = "black", size = 26),
          axis.text.y = element_text(colour = "black", size = 26)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    scale_shape_manual(values = c(22, 21, 24, 25)) +
    scale_fill_nejm() +
    coord_cartesian(ylim = c(0, 10))

  # make plot of scaled MSE
  p_nmse <- sim_summary %>%
    ggplot(aes(x = as.factor(n_obs), y = nmse, group = label, shape = label,
               fill = label)) +
    geom_point(alpha = 0.7, size = 8, position = pd) +
    geom_line(linetype = "dashed", size = 0.5, position = pd) +
    geom_hline(aes(yintercept = eff_var), linetype = "dashed",
               colour = "black") +
    labs(
      x = "Sample size",
      y = TeX("$n \\times$ ($(\\psi - \\hat{\\psi})^2 + \\hat{\\sigma^2})$"),
      title = "Scaled mean squared error"
    ) +
    #facet_grid(vars(est_category), scales = "free_y") +
    theme_bw() +
    theme(legend.background = element_rect(fill = "gray90", size = 0.25,
                                             linetype = "dotted"),
          legend.position = "none",
          legend.title = element_blank(),
          text = element_text(size = 30),
          axis.text.x = element_text(angle = 20, colour = "black", size = 26),
          axis.text.y = element_text(colour = "black", size = 26)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    scale_shape_manual(values = c(22, 21, 24, 25)) +
    scale_fill_nejm() +
    coord_cartesian(ylim = c(0, 10))

  p_cov <- sim_summary %>%
    ggplot(aes(x = as.factor(n_obs), y = coverage, group = label,
               shape = label, fill = label)) +
    geom_point(alpha = 0.7, size = 8, position = pd) +
    geom_line(linetype = "dashed", size = 0.5, position = pd) +
    geom_hline(yintercept = 0.95, linetype = "dashed", colour = "black") +
    labs(
      x = "Sample size",
      y = "Coverage",
      title = "95% Wald-style confidence intervals"
    ) +
    #facet_grid(vars(est_category), scales = "free_y") +
    theme_bw() +
    theme(legend.background = element_rect(fill = "gray90", size = 0.25,
                                             linetype = "dotted"),
          legend.position = "none",
          legend.title = element_blank(),
          text = element_text(size = 30),
          axis.text.x = element_text(angle = 20, colour = "black", size = 26),
          axis.text.y = element_text(colour = "black", size = 26)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    scale_shape_manual(values = c(22, 21, 24, 25)) +
    scale_fill_nejm() +
    ylim(0, 1)

  # output as list of plots
  return(list(samp_dist = p_samp_dist,
              bias = p_bias,
              nbias = p_nbias,
              nmse = p_nmse,
              cov = p_cov))
}

make_usm_plots <- function(perf_summary) {
  # check undersmoothing-selected Mn
  p_mn <- perf_summary %>%
    ggplot(aes(x = Mn_usm, group = as.factor(n_samp))) +
    geom_histogram(alpha = 0.75) +
    facet_grid(vars(n_samp)) +
    geom_vline(aes(xintercept = mean_Mn), linetype = "dashed") +
    theme_bw() +
    theme(legend.background = element_rect(fill = "gray90", size = 0.25,
                                             linetype = "dotted"),
          legend.position = "none",
          legend.title = element_blank(),
          text = element_text(size = 14),
          axis.text.x = element_text(angle = 20, hjust = 1)) +
    labs(
      title = "Selected L1 norm",
      x = "L1 norm"
    )

  # check ratio of undersmoothing-selected Mn to that chosen by global CV
  p_mn_ratio <- perf_summary %>%
    ggplot(aes(x = Mn_ratio, group = as.factor(n_samp))) +
    geom_histogram(alpha = 0.75) +
    facet_grid(vars(n_samp)) +
    geom_vline(aes(xintercept = mean_Mn_ratio), linetype = "dashed") +
    theme_bw() +
    theme(legend.background = element_rect(fill = "gray90", size = 0.25,
                                             linetype = "dotted"),
          legend.position = "none",
          legend.title = element_blank(),
          text = element_text(size = 14),
          axis.text.x = element_text(angle = 20, hjust = 1)) +
    labs(
      title = "Comparison of selected L1 norms",
      x = "Ratio of L1 norms"
    )

  # check truncation across sample sizes
  p_gtrunc <- perf_summary %>%
    ggplot(aes(x = gtrunc, group = as.factor(n_samp))) +
    geom_histogram(alpha = 0.75) +
    facet_grid(vars(n_samp)) +
    geom_vline(aes(xintercept = mean_gtrunc), linetype = "dashed") +
    theme_bw() +
    theme(legend.background = element_rect(fill = "gray90", size = 0.25,
                                             linetype = "dotted"),
          legend.position = "none",
          legend.title = element_blank(),
          text = element_text(size = 14),
          axis.text.x = element_text(angle = 20, hjust = 1)) +
    labs(
      title = "Selected truncation",
      subtitle = TeX("$\\delta$-truncation of $g(1|w)$"),
      x = "Truncation level"
    )

  # results as list
  return(list(mn = p_mn, mn_ratio = p_mn_ratio, gtrunc = p_gtrunc))
}
