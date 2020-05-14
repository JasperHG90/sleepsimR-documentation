# Exploratory analysis of simulation study results

rm(list=ls())
library(sleepsimReval)
library(sleepsimRdata)
library(ggplot2)
library(scales)
library(gridExtra)
library(latex2exp)
library(dplyr)
library(tidyr)
library(purrr)

# Source these scripts
source("helpers/nlp.R")
source("helpers/facet_equal_wrap.R")

# Component distribution fixed effects -----

# Load simulation results

# Load datasets
# (included in 'sleepsimRdata' library)
data("simulation_data_emiss_means")
data("simulation_data_emiss_varmu")
data("simulation_data_gamma_prob")

# Explore Emission means ----

# Get data, parameter results
d <- simulation_data_emiss_means$data_preprocessed
r <- simulation_data_emiss_means$summary_by_scenario
rm(simulation_data_emiss_means)

# Make labels for n/zeta combinations
lbl <- expand.grid(unique(d$n), unique(d$zeta)) %>%
  group_by(Var2) %>%
  arrange(Var2, desc(Var1)) %>%
  mutate(lbl = paste0(Var1, " (z=", Var2, ")"),
         lvl = paste0(Var1, "_", Var2)) %>%
  ungroup() %>%
  select(lbl, lvl)

## Looking for multimodality ----

# Get scenarios for which p-value of multimodal test <= 0.1
mmo <- r %>%
  filter(multimodal <= 0.1) %>%
  mutate(selecvar = paste0(scenario_id, "_", emiss_var_short, "_", state))
mm <- mmo %>%
  select(selecvar) %>%
  pull()

# View param estimates
d %>%
  mutate(selecvar = paste0(scenario_id, "_", emiss_var_short, "_", state)) %>%
  filter(selecvar %in% mm) %>%
  ggplot(aes(x=median)) +
  geom_density() +
  geom_text(data=mmo, aes(x=2, y=0.75, label=multimodal, group = selecvar)) +
  facet_wrap(". ~ selecvar", scales="free_y") 

## Plotting estimates and empirical SE ----

emiss_var_short_nams <- c("EEG mean beta", "EOG median theta", "EOG min beta")
names(emiss_var_short_nams) <- c("EEG_mean_beta", "EOG_median_theta", "EOG_min_beta")

# Plot parameter estimates
# (Across Q)
d %>%
  #filter(n_t == state_lengths[sl]) %>%
  mutate(nsubj = factor(paste0(n, "_", zeta),
                        levels = lbl$lvl,
                        labels = lbl$lbl),
         median = (median-true_val),
         # Need this to label
         n_t_char = factor(paste0("occasions = ", as.character(n_t)),
                           levels = c("occasions = 400", "occasions = 800",
                                      "occasions = 1600"))) %>%
  #select(-Q) %>%
  #filter(Q == 0.4) %>%
  ggplot(aes(x=nsubj, y=median)) +
  geom_hline(yintercept = 0.25, color = "lightgrey") +
  geom_hline(yintercept = 0.5, color = "lightgrey") +
  geom_hline(yintercept = 1, color = "lightgrey") +
  geom_hline(yintercept = 1.5, color = "lightgrey") +
  geom_hline(yintercept = -0.25, color = "lightgrey") +
  geom_hline(yintercept = -0.5, color = "lightgrey") +
  geom_hline(yintercept = -1, color = "lightgrey") +
  geom_hline(yintercept = -1.5, color = "lightgrey") +
  geom_jitter(aes(group = state, color = state,
                  shape = as.factor(zeta)), alpha=0.1) +
  geom_boxplot(aes(group = interaction(nsubj, state), fill=state),
               width = 0.5, position = position_dodge(1),
               outlier.shape = NA) +
  geom_hline(yintercept = 0, color="yellow", size=1,
             linetype = "dashed") +
  coord_flip() +
  theme_thesis(text_size = 6, center_title = TRUE) +
  scale_y_continuous(limits = c(-1.3, 1.3)) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(hjust = 0),
    legend.position = "none"
  ) +
  facet_grid("n_t_char ~ emiss_var_short",
             labeller = labeller(emiss_var_short = emiss_var_short_nams)) +
  ggtitle(TeX(paste0("\\textbf{Mean estimates} ($\\hat{\\theta}_i$)")))

## Collapsed across occasion sizes
d %>%
  #filter(n_t == state_lengths[sl]) %>%
  mutate(nsubj = factor(paste0(n, "_", zeta),
                        levels = lbl$lvl,
                        labels = lbl$lbl),
         median = (median-true_val)) %>%
  #select(-Q) %>%
  #filter(Q == 0.4) %>%
  ggplot(aes(x=nsubj, y=median)) +
  geom_hline(yintercept = 0.25, color = "lightgrey") +
  geom_hline(yintercept = 0.5, color = "lightgrey") +
  geom_hline(yintercept = 1, color = "lightgrey") +
  geom_hline(yintercept = 1.5, color = "lightgrey") +
  geom_hline(yintercept = -0.25, color = "lightgrey") +
  geom_hline(yintercept = -0.5, color = "lightgrey") +
  geom_hline(yintercept = -1, color = "lightgrey") +
  geom_hline(yintercept = -1.5, color = "lightgrey") +
  geom_jitter(aes(color = state),
              alpha=0.2) + #shape = as.factor(zeta)),
  geom_boxplot(aes(group = interaction(nsubj, state), fill=state),
               width = 0.5, position = position_dodge(1),
               outlier.shape = NA) +
  geom_hline(yintercept = 0, color="yellow", size=1,
             linetype = "dashed") +
  coord_flip() +
  theme_thesis(text_size = 9, center_title = TRUE) +
  scale_y_continuous(limits = c(-1.3, 1.3)) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(hjust = 0),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  facet_wrap(". ~ emiss_var_short") +
  ggtitle(TeX(paste0("\\textbf{Mean estimates} ($\\mathbf{\\hat{\\theta}_i$)} \\textbf{. Collapsed accross occasion size and between-subject TPM}")))

# Plot standard error estimates
# Across Q
d %>%
  group_by(var) %>%
  mutate(nsubj = factor(paste0(n, "_", zeta),
                        levels = lbl$lvl,
                        labels = lbl$lbl),
         SE,
         n_t_char = factor(paste0("occasions = ", as.character(n_t)),
                           levels = c("occasions = 400", "occasions = 800",
                                      "occasions = 1600"))) %>%
  #filter(Q == 0.4) %>%
  ggplot(aes(x=nsubj, y=SE)) +
  geom_hline(yintercept = 0.1, color = "lightgrey") +
  geom_hline(yintercept = 0.2, color = "lightgrey") +
  geom_hline(yintercept = 0.3, color = "lightgrey") +
  geom_hline(yintercept = 0.6, color = "lightgrey") +
  geom_jitter(aes(group = state, color = state,
                  shape = as.factor(zeta)), alpha=0.1) +
  geom_boxplot(aes(group = interaction(nsubj, state), fill=state),
               width = 0.5, position = position_dodge(1),
               outlier.shape = NA) +
  coord_flip() +
  theme_thesis(text_size = 6, center_title = TRUE) +
  scale_y_continuous(limits = c(0, 0.8)) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(hjust = 0),
    legend.position = "none"
  ) +
  facet_grid("n_t_char ~ emiss_var_short",
             labeller = labeller(emiss_var_short = emiss_var_short_nams)) +
  ggtitle(TeX(paste0("\\textbf{SE estimates} ($SE(\\hat{\\theta}_i)$)")))

# Plot standard error estimates
# Across Q
# And across n_t
d %>%
  group_by(var) %>%
  mutate(nsubj = factor(paste0(n, "_", zeta),
                        levels = lbl$lvl,
                        labels = lbl$lbl),
         SE) %>%
  ungroup() %>%
  #select(-Q) %>%
  ggplot(aes(x=nsubj, y=SE)) +
  geom_hline(yintercept = 0.1, color = "lightgrey") +
  geom_hline(yintercept = 0.2, color = "lightgrey") +
  geom_hline(yintercept = 0.3, color = "lightgrey") +
  geom_hline(yintercept = 0.6, color = "lightgrey") +
  geom_jitter(aes(group = state, color = state,
                  shape = as.factor(zeta)), alpha=0.1) +
  geom_boxplot(aes(group = interaction(nsubj, state), fill=state),
               width = 0.5, position = position_dodge(1),
               outlier.shape = NA) +
  coord_flip() +
  theme_thesis(text_size = 9, center_title = TRUE) +
  scale_y_continuous(limits = c(0, 0.8)) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(hjust = 0),
    legend.position = "none"
  ) +
  facet_wrap(". ~ emiss_var_short") +
  ggtitle(TeX(paste0("\\textbf{SE estimates} ($SE(\\hat{\\theta}_i)$)")))

## Zipper plot of CCI -----

# Create fractional centile of |z|
d %>%
  # Subset for a variable
  filter(emiss_var_short == "EOG_median_theta") %>%
  # Get modSE for each scenario
  left_join(r %>%
              select(scenario_id, modSE),
            by="scenario_id") %>%
  # Compute the z-score (fc) for each iteration
  # Center the lower/upper CCI
  mutate(fc = abs((median - true_val) / modSE),
         lower = lower - true_val,
         upper = upper - true_val,
         median = median - true_val,
         state_n = paste0(emiss_var_short, "_", n)) %>%
  # Select variables
  select(scenario_id, iteration_id, emiss_var_short,
         state, fc, lower, upper, true_val, median,
         n, n_t, zeta, Q, state_n) %>%
  # Group by emission var and state for each scenario
  group_by(scenario_id, emiss_var_short, state) %>%
  arrange(fc) %>%
  ungroup() %>%
  # Filter for sim settings
  filter(n_t == 800) %>%
  # Compute centiles
  mutate(perc = quantile(fc, seq(0, 1, length.out = n())) %>%
           names() %>%
           stringr::str_replace_all("\\%", "") %>%
           as.numeric(),
         col = ifelse((0 >= lower & 0 <= upper), "blue", "purple")) %>%
  # Filter for quantile
  #filter(perc >= 70) %>%
  # Plot
  ggplot(.) +
  geom_segment(aes(y = perc, yend = perc, x = (lower),
                   xend = (upper), color = col,
                   group = state_n),
               alpha=0.2) +
  scale_x_continuous(name = "95% Central Credible Interval",
                     limits = c(-3, 3)) +
  scale_y_continuous(breaks = c(5, 50, 75, 85, 95),
                     labels = c("5%", "50%", "75%", "85%", "95%")) +
  #geom_hline(yintercept = 5, color = "lightgrey") +
  #geom_hline(yintercept = 50, color = "lightgrey") +
  geom_hline(yintercept = 75, color = "lightgrey") +
  geom_hline(yintercept = 85, color = "lightgrey") +
  geom_hline(yintercept = 95, color = "lightgrey") +
  geom_vline(xintercept=0, color = "yellow",
             linetype = "dashed", size = 0.9) +
  facet_grid("zeta ~ state_n") +
  scale_color_manual(values = c("#2C7085", "#6F447C"),
                     labels = c("Coverers", "Non-coverers")) +
  theme_thesis(text_size=10) +
  theme(legend.title = element_blank(),
        legend.position = "bottom")

## Lollipop plot ----

r1 <- r %>%
  filter(Q == 0.4,
         emiss_var_short == "EEG_mean_beta") %>%
  #select(-modSE, -bias_corr_coverage_mcmc_se) %>%
  mutate(subj_zeta = factor(paste0(n, "_", zeta),
                            levels = lbl$lvl,
                            labels = lbl$lbl))

## Get estimates of parameters of interest
r_est <- r1 %>%
  select(-ends_with("mcmc_se"),
         -multimodal,
         -n, -Q,
         -zeta,
         -true_val) %>%
  gather(estimand, estimate, -scenario_id, -emiss_var_short,
         -state, -subj_zeta, -n_t)

## Get MCMC errors
r_se <- r1 %>%
  select(scenario_id, emiss_var_short, state,
         ends_with("mcmc_se")) %>%
  gather(estimand, mcmc_se, -scenario_id,
         -state, -emiss_var_short) %>%
  mutate(estimand = stringr::str_replace(estimand,
                                         "_mcmc_se",
                                         ""))

## Join datasets
r_out <- r_est %>%
  left_join(r_se,
            by = c("scenario_id",
                   "emiss_var_short",
                   "state",
                   "estimand")) %>%
  mutate(se_lower = estimate - (1.96 * mcmc_se),
         se_upper = estimate + (1.96 * mcmc_se),
         vline_pos = ifelse(estimand == "bias", 0,
                            ifelse(estimand %in% c("coverage", "bias_corr_coverage"), 0.95,
                                   NA)))

# Lollipop plot
r_out %>%
  ggplot(aes(x=subj_zeta, y=estimate, color = state)) +
  #scale_y_continuous(limits = c(-0.5, 0.5)) +
  geom_point(position = position_jitterdodge(0)) +
  geom_linerange(aes(x = subj_zeta,
                     ymin = se_lower,
                     ymax = se_upper),
                 position = position_jitterdodge(0)) +
  coord_flip() +
  geom_hline(aes(yintercept = vline_pos),
             color = "black",
             linetype = "dashed", size=1) +
  theme_thesis(text_size = 5) +
  facet_wrap("estimand~n_t", ncol=3, scales = "free")

# Results are quite similar across n_t, so choose one particular
# (Except coverage, which is a lot worse at n_t == 400)
r_out %>%
  filter(n_t == 800) %>%
  ggplot(aes(x=subj_zeta, y=estimate, color = state)) +
  #scale_y_continuous(limits = c(-0.5, 0.5)) +
  geom_point(position = position_jitterdodge(0)) +
  geom_linerange(aes(x = subj_zeta,
                     ymin = se_lower,
                     ymax = se_upper),
                 position = position_jitterdodge(0)) +
  coord_flip() +
  geom_hline(aes(yintercept = vline_pos),
             color = "black",
             linetype = "dashed", size=1) +
  theme_thesis(text_size = 6) +
  facet_wrap("estimand~state", ncol=3, scales = "free")

# Component distribution random effects -----

# Get data, parameter results
d <-simulation_data_emiss_varmu$data_preprocessed
r <- simulation_data_emiss_varmu$summary_by_scenario
rm(simulation_data_emiss_varmu)

## Emission bet. subject means. Basic checks. ----

## Looking for multimodality ----

# Get scenarios for which p-value of multimodal test <= 0.1
mmo <- r %>%
  filter(multimodal <= 0.1) %>%
  mutate(selecvar = paste0(scenario_id, "_", emiss_var_short, "_", state))
mm <- mmo %>%
  select(selecvar) %>%
  pull()

# View param estimates
d %>%
  mutate(selecvar = paste0(scenario_id, "_", emiss_var_short, "_", state)) %>%
  filter(selecvar %in% mm) %>%
  ggplot(aes(x=median)) +
    geom_density() +
    geom_text(data=mmo, aes(x=2, y=0.75, label=multimodal, group = selecvar)) +
    facet_wrap(". ~ selecvar", scales="free_y") +
    geom_vline(aes(xintercept = zeta), color = "red", linetype="dashed")

## Plotting estimates and empirical SE ----

emiss_var_short_nams <- c("EEG mean beta", "EOG median theta", "EOG min beta")
names(emiss_var_short_nams) <- c("EEG_mean_beta", "EOG_median_theta", "EOG_min_beta")

# Plot parameter estimates
# (Across Q)
d %>%
  #filter(n_t == state_lengths[sl]) %>%
  mutate(nsubj = factor(paste0(n, "_", zeta),
                        levels = lbl$lvl,
                        labels = lbl$lbl),
         median = (median-zeta),
         # Need this to label
         n_t_char = factor(paste0("occasions = ", as.character(n_t)),
                           levels = c("occasions = 400", "occasions = 800",
                                      "occasions = 1600"))) %>%
  #select(-Q) %>%
  #filter(Q == 0.4) %>%
  ggplot(aes(x=nsubj, y=median)) +
  geom_hline(yintercept = 0.25, color = "lightgrey") +
  geom_hline(yintercept = 0.5, color = "lightgrey") +
  geom_hline(yintercept = 1, color = "lightgrey") +
  geom_hline(yintercept = 1.5, color = "lightgrey") +
  geom_hline(yintercept = -0.25, color = "lightgrey") +
  geom_hline(yintercept = -0.5, color = "lightgrey") +
  geom_hline(yintercept = -1, color = "lightgrey") +
  geom_hline(yintercept = -1.5, color = "lightgrey") +
  geom_jitter(aes(group = state, color = state,
                  shape = as.factor(zeta)), alpha=0.1) +
  geom_boxplot(aes(group = interaction(nsubj, state), fill=state),
               width = 0.5, position = position_dodge(1),
               outlier.shape = NA) +
  geom_hline(yintercept = 0, color="yellow", size=1,
             linetype = "dashed") +
  coord_flip() +
  theme_thesis(text_size = 6, center_title = TRUE) +
  scale_y_continuous(limits = c(-1.3, 1.3)) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(hjust = 0),
    legend.position = "none"
  ) +
  facet_grid("n_t_char ~ emiss_var_short",
             labeller = labeller(emiss_var_short = emiss_var_short_nams)) +
  ggtitle(TeX(paste0("\\textbf{Mean estimates} ($\\hat{\\theta}_i$)")))

## Collapsed across occasion sizes
d %>%
  #filter(n_t == state_lengths[sl]) %>%
  mutate(nsubj = factor(paste0(n, "_", zeta),
                        levels = lbl$lvl,
                        labels = lbl$lbl),
         median = (median-zeta)) %>%
  #select(-Q) %>%
  #filter(Q == 0.4) %>%
  ggplot(aes(x=nsubj, y=median)) +
  geom_hline(yintercept = 0.25, color = "lightgrey") +
  geom_hline(yintercept = 0.5, color = "lightgrey") +
  geom_hline(yintercept = 1, color = "lightgrey") +
  geom_hline(yintercept = 1.5, color = "lightgrey") +
  geom_hline(yintercept = -0.25, color = "lightgrey") +
  geom_hline(yintercept = -0.5, color = "lightgrey") +
  geom_hline(yintercept = -1, color = "lightgrey") +
  geom_hline(yintercept = -1.5, color = "lightgrey") +
  geom_jitter(aes(color = state),
              alpha=0.2) + #shape = as.factor(zeta)),
  geom_boxplot(aes(group = interaction(nsubj, state), fill=state),
               width = 0.5, position = position_dodge(1),
               outlier.shape = NA) +
  geom_hline(yintercept = 0, color="yellow", size=1,
             linetype = "dashed") +
  coord_flip() +
  theme_thesis(text_size = 9, center_title = TRUE) +
  scale_y_continuous(limits = c(-2, 2)) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(hjust = 0),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  facet_wrap(". ~ emiss_var_short") +
  ggtitle(TeX(paste0("\\textbf{Mean estimates} ($\\mathbf{\\hat{\\theta}_i$)} \\textbf{. Collapsed accross occasion size and between-subject TPM}")))

# Plot standard error estimates
# Across Q
d %>%
  group_by(var) %>%
  mutate(nsubj = factor(paste0(n, "_", zeta),
                        levels = lbl$lvl,
                        labels = lbl$lbl),
         SE,
         n_t_char = factor(paste0("occasions = ", as.character(n_t)),
                           levels = c("occasions = 400", "occasions = 800",
                                      "occasions = 1600"))) %>%
  #filter(Q == 0.4) %>%
  ggplot(aes(x=nsubj, y=SE)) +
  geom_hline(yintercept = 0.1, color = "lightgrey") +
  geom_hline(yintercept = 0.2, color = "lightgrey") +
  geom_hline(yintercept = 0.3, color = "lightgrey") +
  geom_hline(yintercept = 0.6, color = "lightgrey") +
  geom_jitter(aes(group = state, color = state,
                  shape = as.factor(zeta)), alpha=0.1) +
  geom_boxplot(aes(group = interaction(nsubj, state), fill=state),
               width = 0.5, position = position_dodge(1),
               outlier.shape = NA) +
  coord_flip() +
  theme_thesis(text_size = 6, center_title = TRUE) +
  scale_y_continuous(limits = c(0, 0.8)) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(hjust = 0),
    legend.position = "none"
  ) +
  facet_grid("n_t_char ~ emiss_var_short",
             labeller = labeller(emiss_var_short = emiss_var_short_nams)) +
  ggtitle(TeX(paste0("\\textbf{SE estimates} ($SE(\\hat{\\theta}_i)$)")))

# Plot standard error estimates
# Across Q
# And across n_t
d %>%
  group_by(var) %>%
  mutate(nsubj = factor(paste0(n, "_", zeta),
                        levels = lbl$lvl,
                        labels = lbl$lbl),
         SE) %>%
  ungroup() %>%
  #select(-Q) %>%
  ggplot(aes(x=nsubj, y=SE)) +
  geom_hline(yintercept = 0.1, color = "lightgrey") +
  geom_hline(yintercept = 0.2, color = "lightgrey") +
  geom_hline(yintercept = 0.3, color = "lightgrey") +
  geom_hline(yintercept = 0.6, color = "lightgrey") +
  geom_jitter(aes(group = state, color = state,
                  shape = as.factor(zeta)), alpha=0.1) +
  geom_boxplot(aes(group = interaction(nsubj, state), fill=state),
               width = 0.5, position = position_dodge(1),
               outlier.shape = NA) +
  coord_flip() +
  theme_thesis(text_size = 9, center_title = TRUE) +
  scale_y_continuous(limits = c(0, 0.8)) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(hjust = 0),
    legend.position = "none"
  ) +
  facet_wrap(". ~ emiss_var_short") +
  ggtitle(TeX(paste0("\\textbf{SE estimates} ($SE(\\hat{\\theta}_i)$)")))

## Zipf plot of CCI -----

# Create fractional centile of |z|
d %>%
  # Subset for a variable
  filter(emiss_var_short == "EOG_median_theta") %>%
  # Get modSE for each scenario
  left_join(r %>%
              select(scenario_id, modSE),
            by="scenario_id") %>%
  # Compute the z-score (fc) for each iteration
  # Center the lower/upper CCI
  mutate(fc = abs((median - true_val) / modSE),
         lower = lower - true_val,
         upper = upper - true_val,
         median = median - true_val,
         state_n = paste0(emiss_var_short, "_", n)) %>%
  # Select variables
  select(scenario_id, iteration_id, emiss_var_short,
         state, fc, lower, upper, true_val, median,
         n, n_t, zeta, Q, state_n) %>%
  # Group by emission var and state for each scenario
  group_by(scenario_id, emiss_var_short, state) %>%
  arrange(fc) %>%
  ungroup() %>%
  # Filter for sim settings
  filter(n_t == 800) %>%
  # Compute centiles
  mutate(perc = quantile(fc, seq(0, 1, length.out = n())) %>%
           names() %>%
           stringr::str_replace_all("\\%", "") %>%
           as.numeric(),
         col = ifelse((0 >= lower & 0 <= upper), "blue", "purple")) %>%
  # Filter for quantile
  #filter(perc >= 70) %>%
  # Plot
  ggplot(.) +
  geom_segment(aes(y = perc, yend = perc, x = (lower),
                   xend = (upper), color = col,
                   group = state_n),
               alpha=0.2) +
  scale_x_continuous(name = "95% Central Credible Interval",
                     limits = c(-3, 3)) +
  scale_y_continuous(breaks = c(5, 50, 75, 85, 95),
                     labels = c("5%", "50%", "75%", "85%", "95%")) +
  #geom_hline(yintercept = 5, color = "lightgrey") +
  #geom_hline(yintercept = 50, color = "lightgrey") +
  geom_hline(yintercept = 75, color = "lightgrey") +
  geom_hline(yintercept = 85, color = "lightgrey") +
  geom_hline(yintercept = 95, color = "lightgrey") +
  geom_vline(xintercept=0, color = "yellow",
             linetype = "dashed", size = 0.9) +
  facet_grid("zeta ~ state_n") +
  scale_color_manual(values = c("#2C7085", "#6F447C"),
                     labels = c("Coverers", "Non-coverers")) +
  theme_thesis(text_size=10) +
  theme(legend.title = element_blank(),
        legend.position = "bottom")

## Lollipop plot ----

r1 <- r %>%
  filter(Q == 0.4,
         emiss_var_short == "EEG_mean_beta") %>%
  #select(-modSE, -bias_corr_coverage_mcmc_se) %>%
  mutate(subj_zeta = factor(paste0(n, "_", zeta),
                            levels = lbl$lvl,
                            labels = lbl$lbl))

## Get estimates of parameters of interest
r_est <- r1 %>%
  select(-ends_with("mcmc_se"),
         -multimodal,
         -n, -Q,
         -zeta,
         -true_val) %>%
  gather(estimand, estimate, -scenario_id, -emiss_var_short,
         -state, -subj_zeta, -n_t)

## Get MCMC errors
r_se <- r1 %>%
  select(scenario_id, emiss_var_short, state,
         ends_with("mcmc_se")) %>%
  gather(estimand, mcmc_se, -scenario_id,
         -state, -emiss_var_short) %>%
  mutate(estimand = stringr::str_replace(estimand,
                                         "_mcmc_se",
                                         ""))

## Join datasets
r_out <- r_est %>%
  left_join(r_se,
            by = c("scenario_id",
                   "emiss_var_short",
                   "state",
                   "estimand")) %>%
  mutate(se_lower = estimate - (1.96 * mcmc_se),
         se_upper = estimate + (1.96 * mcmc_se),
         vline_pos = ifelse(estimand == "bias", 0,
                            ifelse(estimand %in% c("coverage", "bias_corr_coverage"), 0.95,
                                   NA)))

# Lollipop plot
r_out %>%
  ggplot(aes(x=subj_zeta, y=estimate, color = state)) +
  #scale_y_continuous(limits = c(-0.5, 0.5)) +
  geom_point(position = position_jitterdodge(0)) +
  geom_linerange(aes(x = subj_zeta,
                     ymin = se_lower,
                     ymax = se_upper),
                 position = position_jitterdodge(0)) +
  coord_flip() +
  geom_hline(aes(yintercept = vline_pos),
             color = "black",
             linetype = "dashed", size=1) +
  theme_thesis(text_size = 5) +
  facet_wrap("estimand~n_t", ncol=3, scales = "free")

# Results are quite similar across n_t, so choose one particular
# (Except coverage, which is a lot worse at n_t == 400)
r_out %>%
  filter(n_t == 800) %>%
  ggplot(aes(x=subj_zeta, y=estimate, color = state)) +
  #scale_y_continuous(limits = c(-0.5, 0.5)) +
  geom_point(position = position_jitterdodge(0)) +
  geom_linerange(aes(x = subj_zeta,
                     ymin = se_lower,
                     ymax = se_upper),
                 position = position_jitterdodge(0)) +
  coord_flip() +
  geom_hline(aes(yintercept = vline_pos),
             color = "black",
             linetype = "dashed", size=1) +
  theme_thesis(text_size = 6) +
  facet_wrap("estimand~state", ncol=3, scales = "free")

# Transition probs ----

## Looking for multimodality ----

d <- simulation_data_gamma_prob$data_preprocessed
r <- simulation_data_gamma_prob$summary_by_scenario

# Get scenarios for which p-value of multimodal test <= 0.1
mmo <- r %>%
  filter(multimodal <= 0.1) %>%
  mutate(selecvar = paste0(scenario_id, "_", transition))
mm <- mmo %>%
  select(selecvar) %>%
  pull()

# View param estimates
d %>%
  mutate(selecvar = paste0(scenario_id, "_", transition)) %>%
  filter(selecvar %in% mm) %>%
  ggplot(aes(x=median)) +
  geom_density() +
  geom_text(data=mmo, aes(x=2, y=0.75, label=multimodal, group = selecvar)) +
  facet_wrap(". ~ selecvar", scales="free_y") 

## Plotting estimates and empirical SE ----

emiss_var_short_nams <- c("EEG mean beta", "EOG median theta", "EOG min beta")
names(emiss_var_short_nams) <- c("EEG_mean_beta", "EOG_median_theta", "EOG_min_beta")

# Plot parameter estimates
# (Across Q)
d %>%
  #filter(n_t == state_lengths[sl]) %>%
  mutate(nsubj = factor(paste0(n, "_", Q)),
         median = (median-true_val),
         # Need this to label
         n_t_char = factor(paste0("occasions = ", as.character(n_t)),
                           levels = c("occasions = 400", "occasions = 800",
                                      "occasions = 1600"))) %>%
  #select(-Q) %>%
  #filter(Q == 0.4) %>%
  ggplot(aes(x=nsubj, y=median)) +
  geom_hline(yintercept = 0.25, color = "lightgrey") +
  geom_hline(yintercept = 0.5, color = "lightgrey") +
  geom_hline(yintercept = 1, color = "lightgrey") +
  geom_hline(yintercept = 1.5, color = "lightgrey") +
  geom_hline(yintercept = -0.25, color = "lightgrey") +
  geom_hline(yintercept = -0.5, color = "lightgrey") +
  geom_hline(yintercept = -1, color = "lightgrey") +
  geom_hline(yintercept = -1.5, color = "lightgrey") +
  geom_jitter(aes(), alpha=0.1) +
  geom_hline(yintercept = 0, color="yellow", size=1,
             linetype = "dashed") +
  coord_flip() +
  theme_thesis(text_size = 6, center_title = TRUE) +
  scale_y_continuous(limits = c(-.1, .1)) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(hjust = 0),
    legend.position = "none"
  ) +
  facet_grid("n_t_char ~ transition") 

# Plot standard error estimates
# Across zeta
d %>%
  #filter(n_t == state_lengths[sl]) %>%
  mutate(nsubj = factor(paste0(n, "_", Q)),
         median = (median-true_val),
         # Need this to label
         n_t_char = factor(paste0("occasions = ", as.character(n_t)),
                           levels = c("occasions = 400", "occasions = 800",
                                      "occasions = 1600"))) %>%
  #filter(Q == 0.4) %>%
  ggplot(aes(x=nsubj, y=SE)) +
  geom_hline(yintercept = 0.1, color = "lightgrey") +
  geom_hline(yintercept = 0.2, color = "lightgrey") +
  geom_hline(yintercept = 0.3, color = "lightgrey") +
  geom_hline(yintercept = 0.6, color = "lightgrey") +
  geom_jitter(alpha=0.1) +
  coord_flip() +
  theme_thesis(text_size = 6, center_title = TRUE) +
  scale_y_continuous(limits = c(0, 0.1)) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(hjust = 0),
    legend.position = "none"
  ) +
  facet_grid("n_t_char ~ transition") 
