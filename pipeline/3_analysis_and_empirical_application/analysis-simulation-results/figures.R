# Final plots used in Analysis

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

# 1. FIGURES AND ANALYSES USED IN THESIS -----

#%% Density plots for the emission distributions -----

# Plot figures 15 and 16 in thesis
data("sleepdata")
# Rename the column variables
vn <- c("EEG mean beta", "EOG median theta", "EOG min beta")
names(vn) <- colnames(sleepdata %>% select(-id, -EEG_Fpz_Cz_mean_theta, -sleep_state))
# Plot (figure 14)
sleepdata %>%
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

# Plot (figure 15)
vn <- c("EEG mean theta", "EOG median theta", "EOG min beta")
names(vn) <- colnames(sleepdata %>% select(-id, -EEG_Fpz_Cz_mean_beta, -sleep_state))
sleepdata %>%
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

#%% Convergence history ----

data("simulation_data_emiss_means")
data("simulation_data_emiss_varmu")
data("simulation_data_gamma_prob")

# Make figure 8 (about model convergence)
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
h <- readRDS("../shiny-app-model-convergence/history.rdsg") %>%
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
# Gather (long format)
hg <- h %>%
  select(-n, -n_t, -zeta, -Q) %>%
  gather(var, param_coverged, -iteration_id, -converged) %>%
  mutate(var = mapping[var]) %>%
  left_join(simdata_combined, by = c("iteration_id", "var"))

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

# Plot the scatter plot, but put alpha to 0. (this is an empty plot basically). 
# This way, we can add the plots one by one and color them as in figure 8.
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

# Add the individual plots
p + 
  geom_point(data = pp %>% filter(type == "var"), size = 2.5,
             aes(x=avgval_val_converged, y = avgval_val_not_converged, color = as.factor(zeta)), alpha = 0.6) +
  geom_point(data = pp %>% filter(type == "mean"), size = 2.5,
             aes(x=avgval_val_converged, y = avgval_val_not_converged), color = "black", alpha = 0.6) +
  geom_point(data = pp %>% filter(type == "tp"), size = 2.5,
             aes(x=avgval_val_converged, y = avgval_val_not_converged), color = "black", alpha = 0.6) 


#%% Component distribution fixed effects ----

library(rsimsum)
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
# Figure 7a
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
# Figure 7b
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

# Figure 7c
ggplot() + 
  theme_minimal()
grid::grid.draw(var_plot_out$EOG_min_beta$state3)

# Empirical SE
state_map <- c("state1" = "Awake", "state2" = "NREM", "state3" = "REM")
var_map <- c("EEG_mean_beta" = "EEG mean beta", "EOG_median_theta" = "EOG median theta",
             "EOG_min_beta" = "EOG min beta")
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

# MSE
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

# Look at the evolution of coverage across conditions
r %>%
  mutate(state = state_map[state],
         emiss_var_short = var_map[emiss_var_short]) %>%
  filter(Q==0.2) %>%
  mutate(empselarger = empirical_se > modSE,
         modselarger = modSE^2 > empirical_se^2) %>%
  group_by(n, zeta) %>%
  summarize(avgcov = mean(coverage))

# Proportion where modSE > emp se OR emp se > modSE
r %>%
  mutate(state = state_map[state],
         emiss_var_short = var_map[emiss_var_short]) %>%
  filter(Q==0.2, n_t == 400) %>%
  mutate(empselarger = empirical_se > modSE,
         modselarger = modSE^2 > empirical_se^2) %>%
  group_by(n, zeta) %>%
  summarize(propmodselarger = sum(modselarger) / n(),
            propempselarger = sum(empselarger) / n())

# Plot coverage of the fixed effects
# Figure 7d

op2 <- simulation_data_emiss_means$summary_by_scenario %>%
  mutate(bias_mcmc_se = abs((bias_mcmc_se / bias)),
         bias = (bias / true_val),
         mcmc_se_perc = (bias_mcmc_se/bias))

# Make a mapping
state_map <- c("state1" = "Awake", "state2" = "NREM", "state3" = "REM")
var_map <- c("EEG_mean_beta" = "EEG mean beta", "EOG_median_theta" = "EOG median theta", "EOG_min_beta" = "EOG min beta")
library(ggplot2)
library(latex2exp)
# Make a trellist plot of coverage
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

# Make a trellis plot for MSE
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

#%% Component distribution random effects -----

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
# Plot NLP
# Figure 9a
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

# Same for MSE 

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

#%% Transition probabilities ----

# Get data
d_out <- simulation_data_gamma_prob$summary_by_scenario
# Make mapping with true values
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

# Mapping for labels (with true values)
state_names <- c("Awake to Awake (.984)", "Awake to NREM (.003)", "Awake to REM (.013)",
                 "NREM to Awake (.007)", "NREM to NREM (.959)", "NREM to REM (.034)",
                 "REM to Awake (.012)", "REM to NREM (.021)", "REM to REM (.967)")
names(state_names) <- c("S1toS1", "S1toS2", "S1toS3",
                        "S2toS1", "S2toS2", "S2toS3",
                        "S3toS1", "S3toS2", "S3toS3")

library(ggplot2)
library(scales)
# Percent bias 
# Figure 10a
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
# Figure 10b
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
ggplot(plotval %>% filter(Q == 0.1, zeta==0.5), 
       aes(x=n_t, y=empirical_se, linetype=as.factor(n))) +
  geom_point(size=1.7) +
  geom_line() +
  facet_wrap(". ~ transition", scales="free_y",
             labeller = labeller(transition = state_names)) 

# Plot MSE
ggplot(plotval %>% filter(Q == 0.2, zeta == 0.5), 
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
# Figure 10b
ggplot(plotval %>% filter(Q == 0.2, zeta == 0.25), 
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



