renv::activate(project = here::here())
library(here)
library(cobalt)
library(data.table)
library(tidyverse)
library(patchwork)
library(latex2exp)
library(ggsci)
library(knitr)
library(kableExtra)
library(conflicted)
conflict_prefer("filter", "dplyr")

# load and clean results
ipw_nhefs_out <- readRDS(here("data", "nhefs_ipw_results.rds"))
ipw_nhefs_results <- ipw_nhefs_out$ipw_results %>%
  filter(
    qsmk == "ate",
    gn_estim %in% c("hal_trunc", "hal_mintrunc", "hal_gcv",
                    "glm_unadj", "glm_main", "glm_quad")
  ) %>%
  mutate(
    label = case_when(
      gn_estim == "hal_trunc" ~ "HAL (undersm., trunc.)",
      gn_estim == "hal_mintrunc" ~ "HAL (undersm., min. trunc.)",
      gn_estim == "hal_gcv" ~ "HAL (global CV)",
      gn_estim == "glm_unadj" ~ "GLM (unadjusted)",
      gn_estim == "glm_main" ~ "GLM (main terms)",
      gn_estim == "glm_quad" ~ "GLM (w/ quad. terms)"
    ),
  )

# visualize results for ATE
p_ate <- ipw_nhefs_results %>%
  ggplot(aes(x = label, y = psi_ipw, group = gn_estim, shape = gn_estim,
             fill = gn_estim)) +
  geom_point(size = 9, alpha = 0.75) +
  geom_errorbar(aes(ymin = ci_lwr, ymax = ci_upr), width = 0.2,
                linetype = "dotted") +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 28),
        axis.text.x = element_text(colour = "black", size = 25),
        axis.text.y = element_text(colour = "black", size = 25)) +
  guides(color = guide_legend(title = NULL)) +
  labs(x = TeX("Estimation strategy for propensity score $\\hat{g}(1|W)$"),
       y = TeX("Average treatment effect $\\hat{\\psi}$"),
       title = "Comparison of parametric versus nonparametric IPW estimators",
       subtitle = paste("Assessing impact of smoking cessation on weight gain",
                        "in the NHEFS cohort")
      ) +
  scale_shape_manual(values = c(10, 12, 13, 22, 24, 25)) +
  scale_fill_nejm()
p_ate
ggsave(filename = here("graphs", "ate_nhefs.pdf"), plot = p_ate,
       width = 20, height = 14)

# make ATE results table
sink(here("tables", "ate_nhefs.tex"))
ipw_nhefs_results %>%
  transmute(
    estimator = label,
    ci_lwr = ci_lwr,
    psi = psi_ipw,
    ci_upr = ci_upr
  ) %>%
  kable(col.names = c("Estimator", "Lower 95% CI", "Estimate", "Upper 95% CI"),
        format = "latex",
        booktabs = TRUE,
        digits = 2,
  ) %>%
  kable_styling() %>%
  print()
sink()

# create balance table
sink(here("tables", "ipw_baltab_nhefs.tex"))
ipw_nhefs_baltab <- ipw_nhefs_out$balance$Balance %>%
  #ipw_nhefs_balance$Balance %>%
  select(-Diff.Un) %>%
  kable(col.names = c("Covariate Type", "GLM (unadjusted)", "GLM (main terms)",
                      "GLM (w/ quad. terms)", "HAL (global CV)",
                      "HAL (undersm., min. trunc.)"),
        format = "latex",
        booktabs = TRUE,
        digits = 4,
  ) %>%
  kable_styling() %>%
  print()
sink()

# create plot of propensity score overlap for HAL vs. GLM-based IPW estimators
ipw_nhefs_pscore <- ipw_nhefs_out$nhefs_results %>%
  select(qsmk, starts_with("ipw")) %>%
  select(-ipw_hal_trunc, -ipw_unadj) %>%
  pivot_longer(
    cols = starts_with("ipw"),
    names_to = "pscore_est",
    values_to = "ip_weights"
  ) %>%
  mutate(
    pscore_est = str_remove(pscore_est, "ipw_"),
    pscore_est = case_when(
      pscore_est == "unadj" ~ "GLM (unadjusted)",
      pscore_est == "mainterms" ~ "GLM (main terms)",
      pscore_est == "robins" ~ "GLM (w/ quad. terms)",
      pscore_est == "hal_gcv" ~ "HAL (global CV)",
      pscore_est == "hal_mintrunc" ~ "HAL (undersm., min. trunc.)",
      pscore_est == "hal_trunc" ~ "HAL (undersm., trunc.)"
    ),
    qsmk = case_when(
      qsmk == 0 ~ "Did not quit smoking",
      qsmk == 1 ~ "Quit smoking"
    )
  )
p_ipw_overlap <- ipw_nhefs_pscore %>%
  ggplot(aes(x = ip_weights)) +
  geom_histogram(fill = "blue", color = "white", alpha = 0.5) +
  xlim(1, 10) +
  theme_bw() +
  theme(
    text = element_text(size = 28),
    axis.text.x = element_text(colour = "black", size = 20),
    axis.text.y = element_text(colour = "black", size = 20)
  ) +
  labs(
    x = "Estimated IP Weight",
    y = "",
    title = "Overlap of Empirical Distributions of Estimated IP Weights",
    subtitle = "(across propensity score estimator and treatment condition)"
  ) +
  facet_grid(qsmk ~ pscore_est, scales = "free_y")
ggsave(filename = here("graphs", "ipw_overlap_nhefs.pdf"),
       plot = p_ipw_overlap, width = 22, height = 18)


# create plots of HAL-IPW solution path in lambda
l1_cv <- ipw_nhefs_out$path_plot$l1_norm[1]
dcar_min_idx <- which.min(ipw_nhefs_out$path_plot$dcar_pn)
l1_dcar_min <- ipw_nhefs_out$path_plot$l1_norm[dcar_min_idx]
psi_dcar_min <- ipw_nhefs_out$path_plot$psi[dcar_min_idx]

## 1) IPW solution path vs. lambda
p_ipw_psi_path <- ipw_nhefs_out$path_plot %>%
  filter(l1_norm <= 5 * l1_cv) %>%
  ggplot(aes(x = l1_norm, y = psi)) +
  geom_line(linetype = "dotted") +
  geom_point(size = 8, alpha = 0.75) +
  geom_vline(xintercept = l1_dcar_min, linetype = "dashed") +
  geom_hline(yintercept = psi_dcar_min, linetype = "dashed") +
  labs(
    x = "",
    y = TeX("IPW point estimate $\\hat{\\psi}$"),
    title = "Trajectory of IPW point estimates of average treatment effect",
    subtitle = ""
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 28),
    text = element_text(size = 32),
    axis.text.x = element_blank(),
    axis.text.y = element_text(colour = "black", size = 30)
  )

## 2) D_CAR solution path vs. lambda
p_ipw_dcar_path <- ipw_nhefs_out$path_plot %>%
  filter(l1_norm <= 5 * l1_cv) %>%
  ggplot(aes(x = l1_norm, y = dcar_pn)) +
  geom_line(linetype = "dotted") +
  geom_point(size = 8, alpha = 0.75) +
  geom_vline(xintercept = l1_dcar, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  labs(
    x = TeX("Magnitude of $L_1$ norm of HAL basis coefficients"),
    y = TeX("$P_n |D_{CAR}|$"),
    title = TeX("Solution trajectory of empirical mean of $D_{CAR}$"),
    subtitle = TeX("across increasing degree of undersmoothing of $g_n$")
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 28),
    text = element_text(size = 32),
    axis.text.x = element_text(angle = 20, colour = "black",
                               size = 30, hjust = 1),
    axis.text.y = element_text(colour = "black", size = 30)
  )

## 3) Combine plots and save as single image
p_ipw_path <- p_ipw_psi_path / p_ipw_dcar_path
p_ipw_path
ggsave(filename = here("graphs", "ipw_path_nhefs.pdf"),
       plot = p_ipw_path, width = 22, height = 18)
