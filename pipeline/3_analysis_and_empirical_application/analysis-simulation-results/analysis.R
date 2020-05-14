# Analyze the results of the simulation study

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

#%% Density plots for the emission distributions -----

d <- readRDS("/home/jasper/GitHubProjects/sleepsimR-sleepdata-analysis/app/sleep_data_subset.rds")
vn <- c("EEG mean beta", "EOG median theta", "EOG min beta")
names(vn) <- colnames(d %>% select(-id, -EEG_Fpz_Cz_mean_theta, -sleep_state))
d %>%
  select(-id, -EEG_Fpz_Cz_mean_theta) %>%
  gather(var, val, -sleep_state) %>%
  ggplot(., aes(x=val, fill = sleep_state)) +
    geom_density(alpha = 0.5) +
    facet_wrap(". ~ var",
               labeller = labeller(var = vn)) +
    scale_y_continuous(name = "Density") +
    scale_x_continuous(name = "Logit-transformed EEG/EOG channel value",
                       limits = c(-3.5, 3)) +
    theme_thesis() +
    labs(fill = "Sleep state") +
    theme(legend.position = "bottom",
          legend.title = element_text(face = "bold"))

vn <- c("EEG mean theta", "EOG median theta", "EOG min beta")
names(vn) <- colnames(d %>% select(-id, -EEG_Fpz_Cz_mean_beta, -sleep_state))
d %>%
  select(-id, -EEG_Fpz_Cz_mean_beta) %>%
  gather(var, val, -sleep_state) %>%
  ggplot(., aes(x=val, fill = sleep_state)) +
  geom_density(alpha = 0.5) +
  facet_wrap(". ~ var",
             labeller = labeller(var = vn)) +
  scale_y_continuous(name = "Density") +
  scale_x_continuous(name = "Logit-transformed EEG/EOG channel value",
                     limits = c(-3.5, 3)) +
  theme_thesis() +
  labs(fill = "Sleep state") +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold"))

# Table with summary statistics
#d %>%
  

# https://cran.r-project.org/web/packages/rsimsum/vignettes/D-nlp.html
# https://bmcmedresmethodol.biomedcentral.com/track/pdf/10.1186/1471-2288-14-129

# Load datasets
# (included in 'sleepsimRdata' library)
data("simulation_data_emiss_means")
data("simulation_data_emiss_varmu")
data("simulation_data_gamma_prob")

# Emission means ----

# Get data, parameter results
d <- simulation_data_emiss_means$data_preprocessed
r <- simulation_data_emiss_means$summary_by_scenario
rm(simulation_data_emiss_means)

#%% Conv hist ----

simdata_combined <- simulation_data_emiss_means$data_preprocessed %>%
  mutate(type = "mean") %>%
  bind_rows(simulation_data_emiss_varmu$data_preprocessed %>%
              mutate(var = paste0(emiss_var_short, "_varmu_", state)) %>%
              mutate(type = "var")) %>%
  bind_rows(simulation_data_gamma_prob$data_preprocessed %>%
              rename(var = transition) %>%
              mutate(type = "tp")) %>%
  filter(!(var %in% c("S1toS1", "S2toS1", "S3toS1")))

# Examine convergence history
h <- readRDS("/home/jasper/GitHubProjects/sleepsimR-convergence/history.rds") %>%
  bind_rows()
h[,7:30] <- apply(h[,7:30], 2, function(x) ifelse(x == "Yes", 1, 0))
# Make a mapping
mapping <- c(
  # Means
  "emiss_1_mu_1" = "EEG_mean_beta_median_state1",
  "emiss_1_mu_2" = "EEG_mean_beta_median_state2",
  "emiss_1_mu_3" = "EEG_mean_beta_median_state3",
  "emiss_2_mu_1" = "EOG_median_theta_median_state1",
  "emiss_2_mu_2" = "EOG_median_theta_median_state2",
  "emiss_2_mu_3" = "EOG_median_theta_median_state3",
  "emiss_3_mu_1" = "EOG_min_beta_median_state1",
  "emiss_3_mu_2" = "EOG_min_beta_median_state2",
  "emiss_3_mu_3" = "EOG_min_beta_median_state3",
  # Variances
  "emiss_1_varmu_1" = "EEG_mean_beta_varmu_state1",
  "emiss_1_varmu_2" = "EEG_mean_beta_varmu_state2",
  "emiss_1_varmu_3" = "EEG_mean_beta_varmu_state3",
  "emiss_2_varmu_1" = "EOG_median_theta_varmu_state1",
  "emiss_2_varmu_2" = "EOG_median_theta_varmu_state2",
  "emiss_2_varmu_3" = "EOG_median_theta_varmu_state3",
  "emiss_3_varmu_1" = "EOG_min_beta_varmu_state1",
  "emiss_3_varmu_2" = "EOG_min_beta_varmu_state2",
  "emiss_3_varmu_3" = "EOG_min_beta_varmu_state3",
  # Probs
  "gamma_int_bar_S1toS2" = "S1toS2",
  "gamma_int_bar_S1toS3" = "S1toS3",
  "gamma_int_bar_S2toS2" = "S2toS2",
  "gamma_int_bar_S2toS3" = "S2toS3",
  "gamma_int_bar_S3toS2" = "S3toS2",
  "gamma_int_bar_S3toS3" = "S3toS3"
)
# Gather
hg <- h %>%
  select(-n, -n_t, -zeta, -Q) %>%
  gather(var, param_coverged, -iteration_id, -converged) %>%
  mutate(var = mapping[var]) %>%
  left_join(simdata_combined, by = c("iteration_id", "var"))

# Functions used to set the scales to equal values.
# See https://fishandwhistle.net/post/2018/modifying-facet-scales-in-ggplot2/
FacetEqualWrap <- ggproto(
  "FacetEqualWrap", FacetWrap,
  
  train_scales = function(self, x_scales, y_scales, layout, data, params) {
    
    # doesn't make sense if there is not an x *and* y scale
    if (is.null(x_scales) || is.null(x_scales)) {
      stop("X and Y scales required for facet_equal_wrap")
    }
    
    # regular training of scales
    ggproto_parent(FacetWrap, self)$train_scales(x_scales, y_scales, layout, data, params)
    
    # switched training of scales (x and y and y on x)
    for (layer_data in data) {
      match_id <- match(layer_data$PANEL, layout$PANEL)
      
      x_vars <- intersect(x_scales[[1]]$aesthetics, names(layer_data))
      y_vars <- intersect(y_scales[[1]]$aesthetics, names(layer_data))
      
      SCALE_X <- layout$SCALE_X[match_id]
      ggplot2:::scale_apply(layer_data, y_vars, "train", SCALE_X, x_scales)
      
      SCALE_Y <- layout$SCALE_Y[match_id]
      ggplot2:::scale_apply(layer_data, x_vars, "train", SCALE_Y, y_scales)
    }
    
  }
)
facet_wrap_equal <- function(...) {
  # take advantage of the sanitizing that happens in facet_wrap
  facet_super <- facet_wrap(...)
  
  ggproto(NULL, FacetEqualWrap,
          shrink = facet_super$shrink,
          params = facet_super$params
  )
}
# Look at the difference of the bias estimates
# Mapping for the facet titles
facet_map <- c("mean" = "Comp. dist. fixed effect",
               "tp" = "Bet. subj. \ntransition probability",
               "var" = "Comp. dist. random effect")
pp <- hg %>%
  group_by(var, param_coverged, type, zeta) %>%
  summarize(avgval = mean(median),
            varse = var(median) / sqrt(n())) %>%
  ungroup() %>%
  mutate(param_coverged = ifelse(param_coverged == 1, "val_converged", "val_not_converged")) %>%
  pivot_wider(names_from = param_coverged, values_from = c(avgval, varse))

p <- ggplot(data = pp, aes(x=avgval_val_converged, y=avgval_val_not_converged, color = as.factor(zeta))) +
    geom_point(size = 3, alpha = 0) +
    facet_wrap_equal(~type, scales = "free", labeller = labeller(type = facet_map)) +
    #facet_wrap(". ~ type", scales = "free") +
    geom_abline(intercept =0 , slope = 1, linetype = "dashed") +
    theme_thesis() +
    theme(legend.position = "bottom",
          axis.title = element_text(size=rel(1.2)),
          strip.text = element_text(size=rel(1.1)),
          axis.text = element_text(size=rel(0.9))) +
    scale_x_continuous(name = "Avg. parameter estimate\n(converged)") +
    scale_y_continuous(name = "Avg. parameter estimate\n(not converged)") +
    labs(color = latex2exp::TeX("$\\zeta$"))

p + 
  geom_point(data = pp %>% filter(type == "var"), size = 2.5,
               aes(x=avgval_val_converged, y = avgval_val_not_converged, color = as.factor(zeta)), alpha = 0.6) +
  geom_point(data = pp %>% filter(type == "mean"), size = 2.5,
             aes(x=val_converged, y = val_not_converged), color = "black", alpha = 0.6) +
  geom_point(data = pp %>% filter(type == "tp"), size = 2.5,
             aes(x=val_converged, y = val_not_converged), color = "black", alpha = 0.6) 
  

# Make labels for n/zeta combinations
# Label bottom-to-top
#lbl <- expand.grid(unique(d$n), unique(d$zeta)) %>%
#  group_by(Var2) %>%
#  arrange(Var2, Var1) %>%
#  mutate(lbl = paste0(Var1, " (z=", Var2, ")"),
#         lvl = paste0(Var1, "_", Var2)) %>%
#  ungroup() %>%
#  select(lbl, lvl)

# Make labels for n/zeta combinations
lbl <- expand.grid(unique(d$n), unique(d$zeta)) %>%
  group_by(Var2) %>%
  arrange(Var2, desc(Var1)) %>%
  mutate(lbl = paste0(Var1, " (z=", Var2, ")"),
         lvl = paste0(Var1, "_", Var2)) %>%
  ungroup() %>%
  select(lbl, lvl)

## Emission bet. subject means. Basic checks. ----

# Get data, parameter results
d <- simulation_data_emiss_varmu$data_preprocessed
r <- simulation_data_emiss_varmu$summary_by_scenario
rm(simulation_data_emiss_varmu)

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
    geom_vline(aes(xintercept = true_val), color = "red", linetype="dashed")

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

# Nested loop plot ----

library(rsimsum)
source("nlp.R")
# Colors
cols <- c("state1" = "Awake", "state2" = "NREM", "state3" = "REM")
states <- c("state1", "state2", "state3")

# Subset data
rsub <- r %>%
  # Add percentage bias
  mutate(bias = (bias / true_val)) %>%
  filter(Q == 0.2) %>%
  #subset(emiss_var_short == "EOG_min_beta") %>%
  select(scenario_id, emiss_var_short, state,
         bias, bias_mcmc_se, empirical_se, modSE, n, n_t, zeta, Q) %>%
  arrange(n, n_t, zeta) %>%
  group_by(emiss_var_short, state) %>%
  mutate(dgm = 1:n()) %>%
  ungroup() %>%
  select(-scenario_id) %>%
  mutate(mstate = paste0(emiss_var_short, "_", state))

### NLP for eeg mean beta
tmp <- rsub %>%
  filter(emiss_var_short == "EEG_mean_beta") %>%
  mutate(state = cols[state]) %>%
  rsimsum::simsum(
    data=., estvarname = "bias", true = 0,
    methodvar = "state",
    se = "bias_mcmc_se",
    by = c("n", "n_t", "zeta")
  )
# Plot
df <- rsimsum:::get_data(tmp, stats="bias")
#colnames(df)[5:7] <- c("N", "N_T", "zeta")
nlp(data=df, methodvar = tmp$methodvar, by = tmp$by, stats = "bias", target = 0, top = TRUE) +
  geom_vline(xintercept = seq(13, 13 * 3, 12), color = "grey", linetype = "dotdash") +
  geom_vline(xintercept = seq(1, (48-3), 4), color = "lightgrey", linetype = "dotted") +
  geom_hline(yintercept = 0.05, color = "black", linetype = "dashed", size=0.8) +
  geom_hline(yintercept = -0.05, color = "black", linetype = "dashed", size=0.8) +
  theme_thesis() +
  labs(color = "State") +
  scale_y_continuous(labels = scales::percent,
                     name = "Percent bias") +
  scale_color_discrete(labels = c("Awake\n(-.36)", "NREM\n(-.6)", "REM\n(.7)")) +
  theme(
    panel.grid.major.x = element_blank() ,
    # explicitly set the horizontal lines (or they will disappear too)
    #panel.grid.major.y = element_line( size=.1, color="black" ),
    legend.title = element_text(face = "bold"),
    legend.position = "bottom"
  ) 

### NLP for eog median theta
tmp <- rsub %>%
  filter(emiss_var_short == "EOG_median_theta") %>%
  mutate(state = cols[state]) %>%
  rsimsum::simsum(
    data=., estvarname = "bias", true = 0,
    methodvar = "state",
    se = "bias_mcmc_se",
    by = c("n", "n_t", "zeta")
  )
# Plot
df <- rsimsum:::get_data(tmp, stats="bias")
#colnames(df)[5:7] <- c("N", "N_T", "zeta")
nlp(data=df, methodvar = tmp$methodvar, by = tmp$by, stats = "bias", target = 0, top = TRUE) +
  geom_vline(xintercept = seq(13, 13 * 3, 12), color = "grey", linetype = "dotdash") +
  geom_vline(xintercept = seq(1, (48-3), 4), color = "lightgrey", linetype = "dotted") +
  geom_hline(yintercept = 0.05, color = "black", linetype = "dashed", size=0.8) +
  geom_hline(yintercept = -0.05, color = "black", linetype = "dashed", size=0.8) +
  theme_thesis() +
  labs(color = "State") +
  scale_y_continuous(labels = scales::percent,
                     name = "Percent bias") +
  scale_color_discrete(labels = c("Awake\n(1.01)", "NREM\n(-1.31)", "REM\n(-.24)")) +
  theme(
    panel.grid.major.x = element_blank() ,
    # explicitly set the horizontal lines (or they will disappear too)
    #panel.grid.major.y = element_line( size=.1, color="black" ),
    legend.title = element_text(face = "bold"),
    legend.position = "bottom"
  ) 

# For each var//state combination
# Colors
cols <- c("state1" = "red", "state2" = "green", "state3" = "#619CFF")
states <- c("state1", "state2", "state3")
vars <- c("EOG_min_beta")
var_plot_out <- vector("list", 1)
names(var_plot_out) <- vars
for(var in vars) {
  state_plot_out <- vector("list", 3)
  names(state_plot_out) <- states
  for(state_name in states) {
    tmp <- rsub %>%
      filter(emiss_var_short == var,
             state == state_name) %>%
      rsimsum::simsum(
        data=., estvarname = "bias", true = 0,
        se = "bias_mcmc_se", methodvar = NULL,
        by = c("n", "n_t", "zeta")
      )
    df <- rsimsum:::get_data(tmp, stats="bias")
    outplot <- nlp(data=df, methodvar = tmp$methodvar, by = tmp$by, stats = "bias", target = 0, top = TRUE, col=cols[state_name]) +
      geom_vline(xintercept = seq(13, 13 * 3, 12), color = "grey", linetype = "dotdash") +
      geom_vline(xintercept = seq(1, (48-3), 4), color = "lightgrey", linetype = "dotted") +
      geom_hline(yintercept = 0.05, color = "black", linetype = "dashed", size=0.8) +
      geom_hline(yintercept = -0.05, color = "black", linetype = "dashed", size=0.8) +
      theme_thesis() +
      labs(color = "State") +
      scale_y_continuous(labels = scales::percent,
                         name = "Percent bias") +
      theme(
        panel.grid.major.x = element_blank() ,
        panel.grid.major.y = element_blank(),
        # explicitly set the horizontal lines (or they will disappear too)
        #panel.grid.major.y = element_line( size=.1, color="black" ),
        legend.title = element_text(face = "bold"),
        legend.position = "bottom"
      )
    state_plot_out[[state_name]] <- ggplotGrob(outplot)
  }
  var_plot_out[[var]] <- state_plot_out
}

# Plot each
ggplot() + 
  theme_minimal()
grid::grid.draw(var_plot_out$EOG_min_beta$state1)
ggplot() + 
  theme_minimal()
grid::grid.draw(var_plot_out$EOG_min_beta$state2)
ggplot() + 
  theme_minimal()
grid::grid.draw(var_plot_out$EOG_min_beta$state3)

# Empirical SE
ggplot(r %>% mutate(MSE = sqrt(MSE), MSE_mcmc_se = sqrt(MSE_mcmc_se)) %>%
         mutate(state = state_map[state],
                emiss_var_short = var_map[emiss_var_short]) %>%
         filter(Q==0.2, n_t==800), aes(x=as.factor(n), y=modSE, color = as.factor(zeta)), alpha=0.6) +
  geom_linerange(aes(ymin = modSE - modSE_mcmc_se, 
                     ymax = modSE + modSE_mcmc_se), size=1.1) +
  geom_point(size=3, alpha=0.6) +
  scale_y_continuous(name = "Empirical SE") +
                     #limits = c(0, 1.30)) +
  facet_grid("state ~ emiss_var_short", scales = "free_y") +
  theme_thesis() +
  theme(legend.position = "bottom")

# Empirical SE
state_map <- c("state1" = "Awake", "state2" = "NREM", "state3" = "REM")
var_map <- c("EEG_mean_beta" = "EEG mean beta", "EOG_median_theta" = "EOG median theta",
             "EOG_min_beta" = "EOG min beta")
ggplot(r %>% mutate(MSE = sqrt(MSE), MSE_mcmc_se = sqrt(MSE_mcmc_se)) %>%
         mutate(state = state_map[state],
                emiss_var_short = var_map[emiss_var_short]) %>%
         filter(Q==0.2, n_t==800), aes(x=as.factor(n), y=MSE, color = as.factor(zeta)), alpha=0.6) +
  geom_linerange(aes(ymin = MSE - MSE_mcmc_se, ymax = MSE + MSE_mcmc_se), size=1.1) +
  geom_point(size=3, alpha=0.6) +
  scale_y_continuous(name = "Root MSE",
                     limits = c(0, 1.50)) +
  scale_x_discrete(name = "Number of subjects (N)") +
  facet_grid("state ~ emiss_var_short", scales = "free_y") +
  labs(color = "Bet. subj. variance") +
  theme_thesis() +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold"))

r %>%
  mutate(state = state_map[state],
         emiss_var_short = var_map[emiss_var_short]) %>%
  filter(Q==0.2) %>%
  mutate(empselarger = empirical_se > modSE,
         modselarger = modSE^2 > empirical_se^2) %>%
  group_by(n, zeta) %>%
  summarize(avgcov = mean(coverage))

r %>%
  mutate(state = state_map[state],
         emiss_var_short = var_map[emiss_var_short]) %>%
  filter(Q==0.2, n_t == 400) %>%
  mutate(empselarger = empirical_se > modSE,
         modselarger = modSE^2 > empirical_se^2) %>%
  group_by(n, zeta) %>%
  summarize(propmodselarger = sum(modselarger) / n(),
            propempselarger = sum(empselarger) / n())

ggplot(r %>% mutate(MSE = sqrt(MSE), MSE_mcmc_se = sqrt(MSE_mcmc_se)) %>%
         mutate(state = state_map[state],
                emiss_var_short = var_map[emiss_var_short]) %>%
         filter(Q==0.2, n_t==800), aes(x=as.factor(n), y=coverage, color = as.factor(zeta)), alpha=0.6) +
  geom_linerange(aes(ymin = coverage - coverage_mcmc_se, ymax = coverage + coverage_mcmc_se), size=1.1) +
  geom_point(size=3, alpha=0.6) +
  scale_y_continuous(name = "Root MSE",
                     limits = c(0, 1)) +
  scale_x_discrete(name = "Number of subjects (N)") +
  facet_grid("state ~ emiss_var_short", scales = "free_y") +
  labs(color = "Bet. subj. variance") +
  theme_thesis() +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold"))

# Same for MSE ----

# Subset data
rsub <- r %>%
  #subset(emiss_var_short == "EOG_min_beta") %>%
  select(scenario_id, emiss_var_short, state,
         MSE, MSE_mcmc_se, n, n_t, zeta, Q) %>%
  arrange(n, n_t, zeta) %>%
  group_by(emiss_var_short, state) %>%
  mutate(dgm = 1:n()) %>%
  ungroup() %>%
  select(-scenario_id) %>%
  mutate(mstate = paste0(emiss_var_short, "_", state))

# For each var//state combination
# Colors
cols <- c("state1" = "red", "state2" = "green", "state3" = "blue")
states <- c("state1", "state2", "state3")
vars <- c("EEG_mean_beta", "EOG_median_theta", "EOG_min_beta")
var_plot_out <- vector("list", 3)
names(var_plot_out) <- vars
for(var in vars) {
  state_plot_out <- vector("list", 3)
  names(state_plot_out) <- states
  for(state_name in states) {
    tmp <- rsub %>%
      filter(emiss_var_short == var,
             state == state_name) %>%
      rsimsum::simsum(
        data=., estvarname = "MSE", true = 0,
        se = "MSE_mcmc_se", methodvar = NULL,
        by = c("n", "n_t", "zeta")
      ) %>%
      autoplot(., type="nlp", stats = "mse") +
      geom_step(aes(x=.scenario,
                    y=est), color = cols[state_name]) +
      geom_vline(xintercept = seq(13, 13 * 3, 12), color = "grey", linetype = "dotdash") +
      geom_vline(xintercept = seq(1, (48-3), 4), color = "lightgrey", linetype = "dotted") +
      theme_thesis() +
      theme(
        panel.grid.major.x = element_blank() ,
        # explicitly set the horizontal lines (or they will disappear too)
        panel.grid.major.y = element_line( size=.1, color="black" ),
        legend.title = element_blank(),
        legend.position = "bottom"
      )
    if(state_name != "state1" & var != "EEG_mean_beta") {
      tmp <- tmp +
        theme(
          axis.title.y = element_blank()
        )
    }
    if(var == "EOG_min_beta" & state_name == "state2") {
      tmp <- tmp
    } else {
      tmp <- tmp +
        theme(
          axis.title.x = element_blank()
        )
    }
    state_plot_out[[state_name]] <- ggplotGrob(tmp)
  }
  var_plot_out[[var]] <- state_plot_out
}

# Plot each
ggplot() +
  theme_minimal()
grid::grid.draw(rbind(
  cbind(var_plot_out$EEG_mean_beta$state1,
        var_plot_out$EEG_mean_beta$state2,
        var_plot_out$EEG_mean_beta$state3,
        size = "last"),
  cbind(var_plot_out$EOG_median_theta$state1,
        var_plot_out$EOG_median_theta$state2,
        var_plot_out$EOG_median_theta$state3,
        size = "last"),
  cbind(var_plot_out$EOG_min_beta$state1,
        var_plot_out$EOG_min_beta$state2,
        var_plot_out$EOG_min_beta$state3,
        size = "last"),
  size = "last"
))

# ModSE checks

rsub <- r %>%
  filter(modSE < empirical_se)
rsub2 <- r %>%
  filter(modSE > empirical_se)
# Dealing with:
#  (1) |bias| > 0
#  (2) In the cases where modSE < empSE : zeta is always 1 or 2.
#  (3) Var(theta_i) is too variable in some cases (esp. EOG_min_beta).

# Emission variances -----

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

# Nested loop plot ----

library(rsimsum)

cols <- c("state1" = "Awake", "state2" = "NREM", "state3" = "REM")
states <- c("state1", "state2", "state3")
# Subset data
rsub <- r %>%
  # Add percentage bias
  mutate(bias = (bias / zeta)) %>%
  #subset(emiss_var_short == "EOG_min_beta") %>%
  select(scenario_id, emiss_var_short, state,
         bias, bias_mcmc_se, n, n_t, zeta, Q) %>%
  arrange(n, n_t, zeta) %>%
  group_by(emiss_var_short, state) %>%
  mutate(dgm = 1:n()) %>%
  ungroup() %>%
  select(-scenario_id) %>%
  mutate(mstate = paste0(emiss_var_short, "_", state))

### NLP for eeg mean beta
tmp <- rsub %>%
  filter(emiss_var_short == "EOG_min_beta") %>%
  mutate(state = cols[state]) %>%
  rsimsum::simsum(
    data=., estvarname = "bias", true = 0,
    methodvar = "state",
    se = "bias_mcmc_se",
    by = c("n", "n_t", "zeta")
  )
# Plot
df <- rsimsum:::get_data(tmp, stats="bias")
#colnames(df)[5:7] <- c("N", "N_T", "zeta")
nlp(data=df, methodvar = tmp$methodvar, by = tmp$by, stats = "bias", target = 0, top = TRUE) +
  geom_vline(xintercept = seq(13, 13 * 3, 12), color = "grey", linetype = "dotdash") +
  geom_vline(xintercept = seq(1, (48-3), 4), color = "lightgrey", linetype = "dotted") +
  geom_hline(yintercept = 0.05, color = "black", linetype = "dashed", size=0.8) +
  geom_hline(yintercept = -0.05, color = "black", linetype = "dashed", size=0.8) +
  theme_thesis() +
  labs(color = "State") +
  scale_y_continuous(labels = scales::percent,
                     name = "Percent bias") +
  scale_color_discrete(labels = c("Awake", "NREM", "REM")) +
  theme(
    panel.grid.major.x = element_blank() ,
    # explicitly set the horizontal lines (or they will disappear too)
    #panel.grid.major.y = element_line( size=.1, color="black" ),
    legend.title = element_text(face = "bold"),
    legend.position = "bottom"
  ) 

# For each var//state combination
# Colors
cols <- c("state1" = "red", "state2" = "green", "state3" = "blue")
states <- c("state1", "state2", "state3")
vars <- c("EEG_mean_beta", "EOG_median_theta", "EOG_min_beta")
var_plot_out <- vector("list", 3)
names(var_plot_out) <- vars
for(var in vars) {
  state_plot_out <- vector("list", 3)
  names(state_plot_out) <- states
  for(state_name in states) {
    tmp <- rsub %>%
      filter(emiss_var_short == var,
             state == state_name) %>%
      rsimsum::simsum(
        data=., estvarname = "bias", true = 0,
        se = "bias_mcmc_se", methodvar = NULL,
        by = c("n", "n_t", "zeta")
      ) %>%
      autoplot(., type="nlp", stats = "bias") +
      geom_step(aes(x=.scenario,
                    y=est), color = cols[state_name]) +
      geom_vline(xintercept = seq(13, 13 * 3, 12), color = "grey", linetype = "dotdash") +
      geom_vline(xintercept = seq(1, (48-3), 4), color = "lightgrey", linetype = "dotted") +
      theme_thesis() +
      theme(
        panel.grid.major.x = element_blank() ,
        # explicitly set the horizontal lines (or they will disappear too)
        panel.grid.major.y = element_line( size=.1, color="black" ),
        legend.title = element_blank(),
        legend.position = "bottom"
      )
    if(state_name != "state1" & var != "EEG_mean_beta") {
      tmp <- tmp +
        theme(
          axis.title.y = element_blank()
        )
    }
    if(var == "EOG_min_beta" & state_name == "state2") {
      tmp <- tmp
    } else {
      tmp <- tmp +
        theme(
          axis.title.x = element_blank()
        )
    }
    state_plot_out[[state_name]] <- ggplotGrob(tmp)
  }
  var_plot_out[[var]] <- state_plot_out
}

# Plot each
ggplot() +
  theme_minimal()
grid::grid.draw(rbind(
  cbind(var_plot_out$EEG_mean_beta$state1,
        var_plot_out$EEG_mean_beta$state2,
        var_plot_out$EEG_mean_beta$state3,
        size = "last"),
  cbind(var_plot_out$EOG_median_theta$state1,
        var_plot_out$EOG_median_theta$state2,
        var_plot_out$EOG_median_theta$state3,
        size = "last"),
  cbind(var_plot_out$EOG_min_beta$state1,
        var_plot_out$EOG_min_beta$state2,
        var_plot_out$EOG_min_beta$state3,
        size = "last"),
  size = "last"
))

# Same for MSE ----

# Subset data
rsub <- r %>%
  #subset(emiss_var_short == "EOG_min_beta") %>%
  select(scenario_id, emiss_var_short, state,
         MSE, MSE_mcmc_se, n, n_t, zeta, Q) %>%
  arrange(n, n_t, zeta) %>%
  group_by(emiss_var_short, state) %>%
  mutate(dgm = 1:n()) %>%
  ungroup() %>%
  select(-scenario_id) %>%
  mutate(mstate = paste0(emiss_var_short, "_", state))

# For each var//state combination
# Colors
cols <- c("state1" = "red", "state2" = "green", "state3" = "blue")
states <- c("state1", "state2", "state3")
vars <- c("EEG_mean_beta", "EOG_median_theta", "EOG_min_beta")
var_plot_out <- vector("list", 3)
names(var_plot_out) <- vars
for(var in vars) {
  state_plot_out <- vector("list", 3)
  names(state_plot_out) <- states
  for(state_name in states) {
    tmp <- rsub %>%
      filter(emiss_var_short == var,
             state == state_name) %>%
      rsimsum::simsum(
        data=., estvarname = "MSE", true = 0,
        se = "MSE_mcmc_se", methodvar = NULL,
        by = c("n", "n_t", "zeta")
      ) %>%
      autoplot(., type="nlp", stats = "mse") +
      geom_step(aes(x=.scenario,
                    y=est), color = cols[state_name]) +
      geom_vline(xintercept = seq(13, 13 * 3, 12), color = "grey", linetype = "dotdash") +
      geom_vline(xintercept = seq(1, (48-3), 4), color = "lightgrey", linetype = "dotted") +
      theme_thesis() +
      theme(
        panel.grid.major.x = element_blank() ,
        # explicitly set the horizontal lines (or they will disappear too)
        panel.grid.major.y = element_line( size=.1, color="black" ),
        legend.title = element_blank(),
        legend.position = "bottom"
      )
    if(state_name != "state1" & var != "EEG_mean_beta") {
      tmp <- tmp +
        theme(
          axis.title.y = element_blank()
        )
    }
    if(var == "EOG_min_beta" & state_name == "state2") {
      tmp <- tmp
    } else {
      tmp <- tmp +
        theme(
          axis.title.x = element_blank()
        )
    }
    state_plot_out[[state_name]] <- ggplotGrob(tmp)
  }
  var_plot_out[[var]] <- state_plot_out
}

# Plot each
ggplot() +
  theme_minimal()
grid::grid.draw(rbind(
  cbind(var_plot_out$EEG_mean_beta$state1,
        var_plot_out$EEG_mean_beta$state2,
        var_plot_out$EEG_mean_beta$state3,
        size = "last"),
  cbind(var_plot_out$EOG_median_theta$state1,
        var_plot_out$EOG_median_theta$state2,
        var_plot_out$EOG_median_theta$state3,
        size = "last"),
  cbind(var_plot_out$EOG_min_beta$state1,
        var_plot_out$EOG_min_beta$state2,
        var_plot_out$EOG_min_beta$state3,
        size = "last"),
  size = "last"
))


# Transition probabilities ----

d_out <- simulation_data_gamma_prob$summary_by_scenario
opts <- options("sleepsimR_simulate")
true_values <- data.frame(
  transition = c("S1toS1", "S1toS2", "S1toS3",
                 "S2toS1", "S2toS2", "S2toS3",
                 "S3toS1", "S3toS2", "S3toS3"),
  true_val = as.vector(opts$sleepsimR_simulate$gamma_bar),
  stringsAsFactors = FALSE
)

# Plotting values
plotval <- d_out %>%
  left_join(true_values) %>%
  mutate(bias_mcmc_se = abs((bias_mcmc_se / bias)),
         bias = (bias / true_val),
         mcmc_se_perc = (bias_mcmc_se/bias)) %>%
  select(transition, bias, bias_mcmc_se, mcmc_se_perc, 
         empirical_se, modSE, coverage, coverage_mcmc_se,
         MSE, MSE_mcmc_se,
         bias_corr_coverage, n_t, n, Q, zeta) %>%
  arrange(n, n_t, transition, mcmc_se_perc)

# For labels
state_names <- c("Awake to Awake (.984)", "Awake to NREM (.003)", "Awake to REM (.013)",
                 "NREM to Awake (.007)", "NREM to NREM (.959)", "NREM to REM (.034)",
                 "REM to Awake (.012)", "REM to NREM (.021)", "REM to REM (.967)")
names(state_names) <- c("S1toS1", "S1toS2", "S1toS3",
                        "S2toS1", "S2toS2", "S2toS3",
                        "S3toS1", "S3toS2", "S3toS3")

library(ggplot2)
library(scales)
ggplot(plotval %>% filter(Q == 0.2,
                          zeta == 0.5),
       aes(x=n_t, y=bias, linetype=as.factor(n))) +
  geom_point(size=2.5) +
  geom_line(size=0.7) +
  # Plot MC SE
  geom_ribbon(aes(ymin = bias - bias_mcmc_se,
                  ymax = bias + bias_mcmc_se), alpha = 0.1) +
  geom_hline(yintercept=0, color = "brown") +
  facet_wrap(". ~ transition", scales="free_y",
             labeller = labeller(transition = state_names)) +
  scale_x_continuous(name = "Occasion size",
                     breaks = c(400, 800, 1200, 1600)) +
  scale_y_continuous(name = "Percent bias",
                     breaks = pretty_breaks(n=4),
                     labels = scales::percent) +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash", "dotted")) + 
  labs(linetype = "Subjects (N)") +
  theme_thesis(text_size=10) +
  scale_color_discrete(guide=guide_legend()) +
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        legend.title.align = 0,
        legend.text.align = 0,
        legend.title = element_text(face="bold"),
        axis.text = element_text(size=rel(1.15))) 

# Check modse > empirical se
plotval_sum <- plotval %>%
  mutate(emplarger = empirical_se > modSE,
         empsmaller = empirical_se^2 < modSE^2)

# For labels
state_names <- c("Awake to Awake", "Awake to NREM", "Awake to REM",
                 "NREM to Awake", "NREM to NREM", "NREM to REM",
                 "REM to Awake", "REM to NREM", "REM to REM")
names(state_names) <- c("S1toS1", "S1toS2", "S1toS3",
                        "S2toS1", "S2toS2", "S2toS3",
                        "S3toS1", "S3toS2", "S3toS3")

# Plot coverage
ggplot(plotval %>% filter(Q == 0.2,
                          zeta == 0.25), 
       aes(x=n_t, y=coverage, linetype=as.factor(n))) +
  geom_point(size=2.5) +
  geom_line(size=0.7) +
  # Plot MC SE
  geom_ribbon(aes(ymin = coverage - coverage_mcmc_se,
                  ymax = coverage + coverage_mcmc_se), alpha = 0.1) +
  geom_hline(yintercept = 0.95, color = "brown") +
  facet_wrap(". ~ transition",
             labeller = labeller(transition = state_names)) +
  scale_x_continuous(name = "Occasion size",
                     breaks = c(400, 800, 1200, 1600),
                     limits = c(400,1650)) +
  scale_y_continuous(name = "Percent coverage",
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     limits = c(0, 1),
                     labels = scales::percent) +
  labs(linetype = "Subjects (N)") +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash", "dotted")) + 
  theme_thesis(text_size=10) +
  scale_color_discrete(guide=guide_legend()) +
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        legend.title.align = 0,
        legend.text.align = 0,
        legend.title = element_text(face="bold"),
        axis.text = element_text(size=rel(1.15))) 

# Plot empirical / model SE
ggplot(plotval %>% filter(Q == 0.1), 
       aes(x=n_t, y=empirical_se, linetype=as.factor(n))) +
  geom_point(size=1.7) +
  geom_line() +
  facet_wrap(". ~ transition", scales="free_y",
             labeller = labeller(transition = state_names)) 

# Plot MSE
ggplot(plotval %>% filter(Q == 0.2), 
       aes(x=n_t, y=sqrt(MSE), linetype=as.factor(n))) +
  geom_point(size=2.5) +
  geom_line(size=0.7) +
  geom_hline(yintercept = 0, color = "brown") +
  # Plot MC SE
  geom_ribbon(aes(ymin = sqrt(MSE) - sqrt(MSE_mcmc_se),
                  ymax = sqrt(MSE) + sqrt(MSE_mcmc_se)), alpha = 0.1) +
  facet_wrap(". ~ transition", scales="free_y",
             labeller = labeller(transition = state_names)) +
  scale_x_continuous(name = "Occasion size",
                     breaks = c(400, 800, 1200, 1600)) +
  scale_y_continuous(name = "Root MSE") +
  labs(linetype = "Subjects (N)") +
  theme_thesis(text_size=10) +
  scale_color_discrete(guide=guide_legend()) +
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        legend.title.align = 0,
        legend.text.align = 0,
        legend.title = element_text(face="bold")) 

# Plot coverage
ggplot(plotval %>% filter(Q == 0.2), 
       aes(x=n_t, y=coverage, linetype=as.factor(n))) +
  geom_point(size=1.7) +
  geom_line() +
  geom_hline(yintercept = 0.95, color = "brown") +
  facet_wrap(". ~ transition", scales="free_y",
             labeller = labeller(transition = state_names)) +
  scale_x_continuous(name = "Occasion size",
                     breaks = c(400, 800, 1200, 1600)) +
  scale_y_continuous(name = "Percent coverage",
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     limits = c(0, 1)) +
  labs(linetype = "Subjects (N)") +
  theme_thesis(text_size=10) +
  scale_color_discrete(guide=guide_legend()) +
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        legend.title.align = 0,
        legend.text.align = 0,
        legend.title = element_text(face="bold")) 

# Bias

op2 <- simulation_data_emiss_means$summary_by_scenario %>%
  mutate(bias_mcmc_se = abs((bias_mcmc_se / bias)),
         bias = (bias / true_val),
         mcmc_se_perc = (bias_mcmc_se/bias))

state_map <- c("state1" = "Awake", "state2" = "NREM", "state3" = "REM")
var_map <- c("EEG_mean_beta" = "EEG mean beta", "EOG_median_theta" = "EOG median theta", "EOG_min_beta" = "EOG min beta")
library(ggplot2)
library(latex2exp)
ggplot(op2 %>%
         mutate(emiss_var = paste0(emiss_var_short, "_", state),
                emiss_var_short = var_map[emiss_var_short],
                state = state_map[state]) %>%
         filter(Q == 0.2, n_t ==1600), aes(x=n, y=coverage, linetype=as.factor(zeta))) +
  geom_point(size=2.5) +
  geom_line(size=0.7) +
  # Plot MC SE
  geom_ribbon(aes(ymin = coverage - coverage_mcmc_se,
                  ymax = coverage + coverage_mcmc_se), alpha = 0.1) +
  geom_hline(yintercept = 0.95, color = "brown") +
  facet_grid("state ~ emiss_var_short", scales="free_y") +
  scale_x_continuous(name = "Subject size",
                     breaks = c(10, 20, 40, 60, 80)) +
  scale_y_continuous(name = "Percent coverage",
                     breaks = c(0.70, 0.8, 0.9, 1),
                     limits = c(0.70, 1),
                     labels = scales::percent) +
  scale_linetype_manual(values = c("solid", "dotted", "dashed", "dotdash")) + 
  labs(linetype = TeX("$\\zeta$ (random effect)")) +
  theme_thesis(text_size=10) +
  scale_color_discrete(guide=guide_legend()) +
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        legend.title.align = 0,
        legend.text.align = 0,
        legend.title = element_text(face="bold")) 

op3 <- op2 %>% filter(Q == 0.2, n_t == 1600) %>%
  mutate(modslarger = modSE^2 > empirical_se^2) %>%
  group_by(n, zeta) %>%
  summarize(proplarger = sum(modslarger) / n())

ggplot(op2 %>%
         mutate(emiss_var = paste0(emiss_var_short, "_", state),
                emiss_var_short = var_map[emiss_var_short],
                state = state_map[state]) %>%
         filter(Q == 0.2, n_t ==1600), aes(x=n, y=MSE, linetype=as.factor(zeta))) +
  geom_point(size=2.5) +
  geom_line(size=0.7) +
  geom_hline(yintercept = 0.95, color = "brown") +
  facet_grid("state ~ emiss_var_short", scales="free_y") 
