# In this file I analyse two models that have been created using docker application
#
#    sleepsimR-sleepdata-analysis <https://github.com/JasperHG90/sleepsimR-sleepdata-analysis>
#
# It contains the following directories:
#
#  - Models (contains mHMM_cont models used to analyze the sleep data)
#  - PPP (contains data related to posterior predictive checks)
#  - data (contains summary statistics & other data related to the sleep data set)
#
# Written by: Jasper Ginn <j.h.ginn@uu.nl>
# Date: 06/04/2020

rm(list=ls())
# Load libraries
library(sleepsimR)
library(sleepsimReval)
library(tidyverse)
library(sleepsimR)
library(mHMMbayes)
library(sleepsimR)
library(dplyr)
library(ggplot2)
#source("utils_PPP.R")
#source("utils_convergence.R")

# Load original dataset & compute summary statistics etc. for later -----

# Select these channels
selected_channels <- c("EEG_Fpz_Cz_mean_theta","EOG_min_beta", "EOG_median_theta") # "EEG_Fpz_Cz_mean_beta"
# Load data
ait <- readRDS("../simulation_parameters/EEG_data_final.rds") %>%
  filter(age >= 20 & age < 50) %>%
  # Remove these variables
  select(-patient, -gender, -age) %>%
  # From wide to long --> var == channels, val == value
  gather(var, val, -identifier, -sleep_state, -epoch) %>%
  # Group by person, variable and state
  group_by(identifier, var, sleep_state) %>%
  # Perform logit transformation
  mutate(val = log(val / (1-val))) %>%
  ungroup() %>%
  # Group by person and variable
  group_by(identifier, var) %>%
  # Grand-mean center variables
  mutate(val = ((val - mean(val)) / sd(val))) %>%
  ungroup()
# Filter the data
io <- ait %>%
  filter(var %in% selected_channels) %>%
  # Pivot to wide
  pivot_wider(id_cols = c(identifier, epoch, sleep_state), names_from = var, values_from = val) %>%
  # Remove epochs
  select(-epoch)

# Get unique identifiers
# Map to numeric 1-len(identifer))
uid <- unique(as.character(io$identifier))
ids <- 1:length(uid)
names(ids) <- uid
io$id <- unname(ids[as.character(io$identifier)])
io$identifier <- NULL
io <- io[,c(5,1, 2, 3, 4)]
# Plot distributions
io %>%
  gather(var, val, -sleep_state, -id) %>%
  ggplot(aes(x=val, fill=sleep_state)) +
  geom_density(alpha = 0.4) +
  facet_wrap(". ~ var", ncol=3)
# Plot average values
io_avg <- io %>%
  gather(var, val, -sleep_state, -id) %>%
  group_by(id, var, sleep_state) %>%
  summarize(avgch = mean(val))
io_avg %>%
  ggplot(aes(x=avgch, fill=sleep_state)) +
  geom_density(alpha=0.4) +
  facet_wrap("var ~ .")
# Get observed sample between-subject variance
io_avg_betvar <- io_avg %>%
  group_by(var, sleep_state) %>%
  summarize(betvar = var(avgch))
# Between-subject variance is very small on means. (<.1)
# Know that the model is biased upwards on these kinds of values.

# Get summary statistics for each dep var
ss <- io %>%
  gather(variable, value, -id, -sleep_state) %>%
  group_by(id, variable, sleep_state) %>%
  summarize(mvar = mean(value)) %>%
  ungroup() %>%
  group_by(variable, sleep_state) %>%
  summarize(mmvar = mean(mvar),
            vvvar = var(mvar))
# Total variance
tvar <- io %>%
  gather(variable, value, -id, -sleep_state) %>%
  group_by(variable, sleep_state) %>%
  summarize(tvar = var(value))
# Remove sleep states from data
io_orig <- io
io$sleep_state <- NULL

# Save
saveRDS(io_orig, "data/final_dataset_modeling.rds")
saveRDS(ss, "data/summary_stats.rds")
saveRDS(tvar, "data/total_variance.rds")

# Load models and check convergence -----

# Load models
m1 <- readRDS("models/model_5194447c-12ae-4c31-968d-b26b6533b0c7.rds")
m2 <- readRDS("models/model_23862780-9d25-41db-850a-e9d38e8743b9.rds")

# To list
mod_res <- list(
  m1$model,
  m2$model
)

# MAP
m1 <- MAP(mod_res[[1]])
m2 <- MAP(mod_res[[2]])

# Look at trace plots ----

library(gridBase)
library(gridExtra)
library(grid)

# Emiss. dist 1
tpp(mod_res[[1]], mod_res[[2]], "emiss_mu_bar", var=1)
# second base plot
# Density plot
dens_plot(mod_res[[1]], mod_res[[2]], "emiss_mu_bar", var=1)
autocorr_plot(mod_res[[1]], "emiss_mu_bar", var=1, thin = 5)
# Variances
tpp(mod_res[[1]], mod_res[[2]], "emiss_varmu_bar", var=1)
autocorr_plot(mod_res[[1]], "emiss_varmu_bar", var=1)
dens_plot(mod_res[[1]], mod_res[[2]], "emiss_varmu_bar", var=1)

# Emiss. dist 2
tpp(mod_res[[1]], mod_res[[2]], "emiss_mu_bar", var=2)
# second base plot
# Density plot
dens_plot(mod_res[[1]], mod_res[[2]], "emiss_mu_bar", var=2)
autocorr_plot(mod_res[[1]], "emiss_mu_bar", var=2, thin = 5)
# Variances
tpp(mod_res[[1]], mod_res[[2]], "emiss_varmu_bar", var=2)
autocorr_plot(mod_res[[1]], "emiss_varmu_bar", var=2)
dens_plot(mod_res[[1]], mod_res[[2]], "emiss_varmu_bar", var=2)

# Emiss. dist 3
tpp(mod_res[[1]], mod_res[[2]], "emiss_mu_bar", var=3)
# second base plot
# Density plot
dens_plot(mod_res[[1]], mod_res[[2]], "emiss_mu_bar", var=3)
autocorr_plot(mod_res[[1]], "emiss_mu_bar", var=3, thin = 5)
# Variances
tpp(mod_res[[1]], mod_res[[2]], "emiss_varmu_bar", var=3)
autocorr_plot(mod_res[[1]], "emiss_varmu_bar", var=3)
dens_plot(mod_res[[1]], mod_res[[2]], "emiss_varmu_bar", var=3)

# MLR intercepts
tpp(mod_res[[1]], mod_res[[2]], "gamma_int_bar", var=3)
# second base plot
# Density plot
dens_plot(mod_res[[1]], mod_res[[2]], "gamma_int_bar", var=3)
autocorr_plot(mod_res[[1]], "gamma_int_bar", var=3, thin = 5)

# Emission variances

tpp(mod_res[[1]], mod_res[[2]], "emiss_varmu_bar", var=2)
tpp(mod_res[[1]], mod_res[[2]], "emiss_varmu_bar", var=3)
# Gamma intercepts
tpp(mod_res[[1]], mod_res[[2]], "gamma_int_bar", var=1)

# Compute gelman-rubin statistic -----

# Emission means
compute_grs(mod_res[[1]], mod_res[[2]], "emiss_mu_bar", var=1, thin = 5)
compute_grs(mod_res[[1]], mod_res[[2]], "emiss_mu_bar", var=2, thin =5)
compute_grs(mod_res[[1]], mod_res[[2]], "emiss_mu_bar", var=3, thin =5)
# Emission variances
compute_grs(mod_res[[1]], mod_res[[2]], "emiss_varmu_bar", var=1, thin = 5)
compute_grs(mod_res[[1]], mod_res[[2]], "emiss_varmu_bar", var=2, thin = 5)
compute_grs(mod_res[[1]], mod_res[[2]], "emiss_varmu_bar", var=3, thin = 5)
# Gamma intercepts
compute_grs(mod_res[[1]], mod_res[[2]], "gamma_int_bar", var=1, thin = 5)

# Autocorrelation plots -----


autocorr_plot(mod_res[[2]], "emiss_mu_bar", var=1, thin = 5)
autocorr_plot(mod_res[[1]], "emiss_mu_bar", var=2)
autocorr_plot(mod_res[[2]], "emiss_mu_bar", var=2)
autocorr_plot(mod_res[[1]], "emiss_mu_bar", var=3)
autocorr_plot(mod_res[[2]], "emiss_mu_bar", var=3)


autocorr_plot(mod_res[[2]], "emiss_varmu_bar", var=1)
autocorr_plot(mod_res[[1]], "emiss_varmu_bar", var=2)
autocorr_plot(mod_res[[2]], "emiss_varmu_bar", var=2)
autocorr_plot(mod_res[[1]], "emiss_varmu_bar", var=3)
autocorr_plot(mod_res[[2]], "emiss_varmu_bar", var=3)

autocorr_plot(mod_res[[1]], "gamma_int_bar")
autocorr_plot(mod_res[[2]], "gamma_int_bar")

# Density plots -----


dens_plot(mod_res[[1]], mod_res[[2]], "emiss_mu_bar", var=2)
dens_plot(mod_res[[1]], mod_res[[2]], "emiss_mu_bar", var=3)
# Emission variances

dens_plot(mod_res[[1]], mod_res[[2]], "emiss_varmu_bar", var=2)
dens_plot(mod_res[[1]], mod_res[[2]], "emiss_varmu_bar", var=3)
# Gamma intercepts
dens_plot(mod_res[[1]], mod_res[[2]], "gamma_int_bar", var=1)

# MAP Estimates ----

m1$emiss_mu_bar$EEG_Fpz_Cz_mean_theta$mean
m2$emiss_mu_bar$EEG_Fpz_Cz_mean_theta$mean
m1$emiss_mu_bar$EOG_median_theta$mean
m2$emiss_mu_bar$EOG_median_theta$mean
m1$emiss_mu_bar$EOG_min_beta$mean
m2$emiss_mu_bar$EOG_min_beta$mean

m1$emiss_varmu_bar$EEG_Fpz_Cz_mean_theta$mean
m2$emiss_varmu_bar$EEG_Fpz_Cz_mean_theta$mean
m1$emiss_varmu_bar$EOG_median_theta$mean
m2$emiss_varmu_bar$EOG_median_theta$mean
m1$emiss_varmu_bar$EOG_min_beta$mean
m2$emiss_varmu_bar$EOG_min_beta$mean

# Combine posterior distributions and compute summary statistics ----

posterior_combined_out <- combine_posterior_chains(mod_res)
# Save posterior
#saveRDS(posterior_combined_out,"data/posterior_combined_out.rds")
posterior_combined_out <- readRDS("data/posterior_combined_out.rds")

# Get summary statistics
# Emiss mu bar (subject-level means)
t(as.matrix(round(apply(posterior_combined_out$emiss_mu_bar, 2, function(x) {
  c(median(x), sd(x), quantile(x, c(0.025, 0.975)))
}), 3)))

# Emiss varmu bar (between-subject variance)
t(as.matrix(round(apply(posterior_combined_out$emiss_varmu_bar, 2, function(x) {
  c(median(x), sd(x), quantile(x, c(0.025, 0.975)))
}), 3)))

# gamma int bar (subject-level linear predictor on TPMs)
t(as.matrix(round(apply(posterior_combined_out$gamma_int_bar, 2, function(x) {
  c(median(x), sd(x), quantile(x, c(0.025, 0.975)))
}), 3)))

# gamma prob bar (subject-level transition probabilities)
as.matrix(round(apply(posterior_combined_out$gamma_prob_bar, 2, function(x) {
  median(x)
}), 3))

# Plot between-subject values for TPM and emission distributions -----

library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)
library(purrr)
# Between-subject TPM values
o <- posterior_combined_out$gamma_prob_subj %>%
  gather(var, val, -subj_idx) %>%
  group_by(subj_idx, var) %>%
  summarize(mapmed = median(val))
varmap <- c(
  "S1toS1" ="Awake to Awake",
  "S1toS2" ="Awake to REM",
  "S1toS3" ="Awake to NREM",
  "S2toS1" ="NREM to Awake",
  "S2toS2" ="NREM to NREM",
  "S2toS3" ="NREM to REM",
  "S3toS1" ="REM to Awake",
  "S3toS2" ="REM to NREM",
  "S3toS3" ="REM to REM"
)
o %>%
  mutate(var = varmap[var]) %>%
  ggplot(., aes(x=mapmed)) +
  geom_histogram(color = "black",  bins=60,fill = "grey", alpha=0.2) +
  scale_x_continuous(name = "Transition Probability",
                     limits = c(0, 1),
                     breaks = seq(0, 1, .2),
                     labels = seq(0, 1, .2)) +
  scale_y_continuous(name = "Density") +
  facet_wrap(". ~ var", scales = "free_y") +
  theme_thesis()

# Emission means for subjects
k <- posterior_combined_out$PD_subj %>%
  gather(var, val, -subj_idx) %>%
  group_by(subj_idx, var) %>%
  summarize(mapmed = median(val))
k %>%
  ungroup() %>%
  mutate(depvar = map_chr(str_split(k$var, "_"), function(x) x[1]),
         state = map_chr(str_split(k$var, "_"), function(x) paste0(x[2], "_", x[3]))) %>%
  filter(stringr::str_ends(var, "mu_S3")) %>%
  mutate(var=ifelse(
    var == "dep1_mu_S3",
    "EEG mean theta",
    ifelse(
      var == "dep2_mu_S3",
      "EOG median theta",
      "EOG min beta"
    )
  )) %>%
  ggplot(., aes(x=mapmed, group=var)) +
  geom_histogram(alpha=0.2, fill="grey", color="black", bins=20) +
  scale_x_continuous(limits = c(-2, 2),
                     breaks = c(-2, 0, 2),
                     name = "Subject-specific means") +
  scale_y_continuous(name = "Density") +
  facet_wrap(". ~ var", scales= "free_y") +
  theme_thesis()

# Posterior predictive checks -----

# For 2000 iterations, make new data from posterior samples
#PPP_out <- vector("list", 2000)
#for(idx in seq_along(PPP_out)) {
  # Generate data
#  out_data <- simulate_new_data(posterior_combined_out)
#  PPP_out[[idx]] <- out_data
#}
#saveRDS(PPP_out, "PPP/PPP_samples.rds")

# Load sample data
#sample_data <- readRDS("PPP/PPP_samples.rds")
# Compute PPP (these are functions in sleepsimReval library)
#PPP_mean_out <- vector("list", 2000)
#PPP_var_out <- vector("list", 2000)
#PPP_tpm_out <- vector("list", 2000)
#for(idx in seq_along(PPP_mean_out)) {
#  md <- PPP_mean_var(sample_data[[idx]]$data)
#  md2 <- data.frame(matrix(md[[1]], ncol = length(md[[1]])))
#  md3 <- data.frame(matrix(md[[2]], ncol = length(md[[2]])))
#  colnames(md2) <- colnames(md3) <- c(
#    "emiss_distribution_1_state_1", "emiss_distribution_1_state_2", "emiss_distribution_1_state_3",
#    "emiss_distribution_2_state_1", "emiss_distribution_2_state_2", "emiss_distribution_2_state_3",
#    "emiss_distribution_3_state_1", "emiss_distribution_3_state_2", "emiss_distribution_3_state_3"
#  )
#  PPP_mean_out[[idx]] <- md2
#  PPP_var_out[[idx]] <- md3
#  PPP_tpm_out[[idx]] <- PPP_tpm(sample_data[[idx]]$data)
#}
# Save
#saveRDS(list("PPP_mean" = PPP_mean_out,
#             "PPP_var" = PPP_var_out,
#             "PPP_tpm" = PPP_tpm_out),
#        "PPP/PPP_test_data.rds")

# Read
PPPts <- readRDS("PPP/PPP_test_data.rds")
ss <- readRDS("data/summary_stats.rds")
tvar <- readRDS("data/total_variance.rds")

# Plot PPP for mean
ola <- do.call(rbind.data.frame, PPPts$PPP_mean)
# Fceting variable
facvar <- data.frame(
  "val" = ss$mmvar,
  "var" = colnames(ola)
) %>%
  filter(str_ends(var, "_state_3")) %>%
  mutate(var2 = c(
    "EEG mean theta",
    "EOG median theta",
    "EOG min beta"
  ))
# Bind
### Get plot for state three only
library(stringr)
ola %>%
  gather(var, val) %>%
  filter(str_ends(var, "_state_3")) %>%
  mutate(var2 = ifelse(
    str_detect(var, "distribution_1"),
    "EEG mean theta",
    ifelse(
      str_detect(var, "distribution_2"),
      "EOG median theta",
      "EOG min beta"
    )
  )) %>%
  ggplot(aes(x=val, group=var2)) +
  geom_histogram(alpha=0.2, fill = "grey", color = "black") +
  geom_vline(data=facvar, aes(xintercept = val, group=var2), size=1.2, color = "red", linetype = "dashed") +
  facet_wrap(". ~ var2", scales = "free") +
  theme_thesis() +
  scale_x_continuous(name = "Subject-level means") +
  scale_y_continuous(name = "Density") +
  theme(axis.text = element_text(size=rel(1.1)))

# p-values
ola %>%
  gather(var, val) %>%
  left_join(facvar, by="var") %>%
  mutate(eq = val.x >= val.y) %>%
  group_by(var) %>%
  summarize(pval = round(sum(eq) / n(), 2))

# Plot PPP for variance
ola2 <- do.call(rbind.data.frame, PPPts$PPP_var)
# Fceting variable
facvar <- data.frame(
  "val" = tvar$tvar,
  "var" = colnames(ola)
)
# Bind
ola2 %>%
  gather(var, val) %>%
  ggplot(aes(x=val, group=var)) +
  geom_histogram() +
  geom_vline(data=facvar, aes(xintercept = val, group=var),
             size=1.2, color = "red", linetype = "dashed") +
  facet_wrap(". ~ var", scales = "free")

# Transition matrix counts
tpm_counts <- lapply(PPPts$PPP_tpm, function(x) {
  subj_counts <- lapply(x$tpms, function(y) {
    period_count <- lapply(y, function(z) {
      as.data.frame(matrix(z, ncol = length(z)))
    }) %>%
      do.call(rbind.data.frame, .)
  }) %>%
    do.call(rbind.data.frame, .)
  idxsubj <- rep(1:41, 3)
  subj_counts$idxsubj <- idxsubj[order(idxsubj)]
  subj_counts$period <- rep(1:3, 41)
  subj_counts
}) %>%
  do.call(rbind.data.frame, .)
# Column names
colnames(tpm_counts)[1:9] <- c(
  "S1toS1", "S1toS2", "S1toS3",
  "S2toS1", "S2toS2", "S2toS3",
  "S3toS1", "S3toS2", "S3toS3"
)
# Add dataset identifier
tpm_counts$dataset <- map(1:2000, function(x) rep(x, 123)) %>%
  unlist()

# Load dataset used for modeling
io_orig <- readRDS("data/final_dataset_modeling.rds")

# Compute average across all datasets and periods
cp <- tpm_counts %>%
  gather(var, val, -idxsubj, -period, -dataset) %>%
  mutate(from_state = stringr::str_replace(var, "toS[1-3]", ""),
         val = val + runif(n(), 0.0001, 0.0002)) %>%
  group_by(from_state, idxsubj, period, dataset) %>%
  arrange(period, idxsubj, dataset, var) %>%
  mutate(val = val / sum(val)) %>%
  ungroup() %>%
  group_by(var, dataset, period) %>%
  summarize(avgbyvar = mean(val)) %>%
  mutate(avgbyvar = round(avgbyvar, 3)) %>%
  ungroup()

# Real sleep state data
sleep_states <- data.frame(
  "subj" = io_orig$id,
  "states" = as.numeric(io_orig$sleep_state)
)
# TPM values
sleep_states_val <- PPP_tpm(sleep_states)

# Make into df
sleep_states_cc <- sleep_states_val$tpms %>%
  lapply(., function(y) {
    lapply(y, function(z) {
      as.data.frame(matrix(z, ncol = length(z)))
    }) %>%
      do.call(rbind.data.frame, .)
  }) %>%
  do.call(rbind.data.frame, .)
idxsubj <- rep(1:41, 3)
sleep_states_cc$idxsubj <- idxsubj[order(idxsubj)]
sleep_states_cc$period <- rep(1:3, 41)
# Column names
colnames(sleep_states_cc)[1:9] <- c(
  "S1toS1", "S1toS2", "S1toS3",
  "S2toS1", "S2toS2", "S2toS3",
  "S3toS1", "S3toS2", "S3toS3"
)

# Take proportions
cc <- sleep_states_cc %>%
  mutate(dataset = rep(1, n())) %>%
  gather(var, val, -idxsubj, -period, -dataset) %>%
  mutate(from_state = stringr::str_replace(var, "toS[1-3]", ""),
         val = val + runif(n(), 0.0001, 0.0002)) %>%
  group_by(from_state, idxsubj, period, dataset) %>%
  arrange(period, idxsubj, dataset, var) %>%
  mutate(val = val / sum(val)) %>%
  ungroup() %>%
  group_by(var, dataset, period) %>%
  summarize(avgbyvar = mean(val)) %>%
  mutate(avgbyvar = round(avgbyvar, 3)) %>%
  ungroup() %>%
  mutate(Period = as.factor(period))

# Mapping of names
nammap <- c(
  "S1toS1" ="Awake to Awake",
  "S1toS2" ="Awake to REM",
  "S1toS3" ="Awake to NREM",
  "S2toS1" ="NREM to Awake",
  "S2toS2" ="NREM to NREM",
  "S2toS3" ="NREM to REM",
  "S3toS1" ="REM to Awake",
  "S3toS2" ="REM to NREM",
  "S3toS3" ="REM to REM"
)

# Plot
cp %>%
  mutate(var = nammap[var],
         Period = as.factor(period)) %>%
  ggplot(aes(x=avgbyvar, fill = Period)) +
    geom_density(alpha=0.2) +
    facet_wrap(". ~ var", scales = "free_y") +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    geom_vline(data=cc %>%
                 mutate(var = nammap[var]),
               aes(xintercept=avgbyvar, color = Period),
               linetype="dashed", size=1.2) +
    scale_x_continuous(name = "Transition probability") +
    scale_y_continuous(name = "Density") +
    theme_thesis() +
    theme(legend.position = "bottom") +
    theme(legend.title = element_text(face = "bold"),
          axis.text = element_text(size=rel(1.2)))

# P-values
p_values_tpm <- cp %>%
  left_join(cc %>% ungroup() %>% select(-dataset, -period), by="var") %>%
  mutate(is_larger = avgbyvar.x >= avgbyvar.y) %>%
  group_by(var, period) %>%
  summarize(p_value = round(sum(is_larger) / n(), 2))

# State durations
sdsim <- lapply(PPPts$PPP_tpm, function(x) {
  x <- x$lengths
  # For each dataset, go through each subject
  subjmax <- lapply(x, function(y) {
    out_d <- vector("list", 3)
    for(per_idx in seq_along(y)) {
      o <- y[[per_idx]]
      per <- vector("list", 3)
      for(p in seq_along(per)) {
        per[[p]] <- ifelse(length(o[[p]]) == 0, 0, max(o[[p]]))
      }
      out_d[[per_idx]] <- do.call(cbind, per)
    }
    out_em <- as.data.frame(do.call(cbind, out_d))
    colnames(out_em) <- paste0("period_", c(rep(1, 3), rep(2, 3), rep(3,3)),
                               "_state_", rep(1:3, sqrt(ncol(out_em))))
    return(out_em)
  })
  subj_out <- do.call(rbind.data.frame, subjmax)
  # Take average across subjects
  t(as.data.frame(round(apply(subj_out, 2, mean))))
})

# Bind
sdsim <- do.call(rbind.data.frame, sdsim)
row.names(sdsim) <- 1:nrow(sdsim)
sdsim$dataset <- 1:2000

# True values
ssl <- list(sleep_states_val$lengths)
sdsim_cc <- lapply(ssl, function(x) {
  # For each dataset, go through each subject
  subjmax <- lapply(x, function(y) {
    out_d <- vector("list", 3)
    for(per_idx in seq_along(y)) {
      o <- y[[per_idx]]
      per <- vector("list", 3)
      for(p in seq_along(per)) {
        per[[p]] <- ifelse(length(o[[p]]) == 0, 0, max(o[[p]]))
      }
      out_d[[per_idx]] <- do.call(cbind, per)
    }
    out_em <- as.data.frame(do.call(cbind, out_d))
    colnames(out_em) <- paste0("period_", c(rep(1, 3), rep(2, 3), rep(3,3)),
                               "_state_", rep(1:3, sqrt(ncol(out_em))))
    return(out_em)
  })
  subj_out <- do.call(rbind.data.frame, subjmax)
  # Take average across subjects
  tc <- as.data.frame(round(apply(subj_out, 2, mean)))
  tc$var <- row.names(tc)
  colnames(tc) <- c("val", "var")
  row.names(tc) <- 1:nrow(tc)
  tc
}) %>%
  .[[1]]

# Plot
sdsim %>%
  gather(var, val, -dataset) %>%
  ggplot(aes(x=val)) +
  geom_density(fill = "grey", alpha=0.15) +
  geom_vline(data=sdsim_cc, aes(xintercept = val, color=var), linetype="dashed", size=0.9) +
  facet_wrap(". ~ var") +
  theme_bw() +
  theme(legend.position = "none")

# P-values
sdsim %>%
  gather(var, val, -dataset) %>%
  left_join(sdsim_cc, by="var") %>%
  mutate(est_over = val.x >= val.y) %>%
  group_by(var) %>%
  summarize(p_val = sum(est_over) / n())

tpm_counts <- lapply(PPPts$PPP_tpm, function(x) {
  subj_counts <- lapply(x$tpms, function(y) {
    period_count <- lapply(y, function(z) {
      as.data.frame(matrix(z, ncol = length(z)))
    }) %>%
      do.call(rbind.data.frame, .)
  }) %>%
    do.call(rbind.data.frame, .)
  idxsubj <- rep(1:41, 3)
  subj_counts$idxsubj <- idxsubj[order(idxsubj)]
  subj_counts$period <- rep(1:3, 41)
  subj_counts
}) %>%
  do.call(rbind.data.frame, .)
