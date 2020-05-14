library(sleepsimReval)
rm(list=ls())

# Read model
library(sleepsimR)
library(dplyr)
mod <- readRDS("e62f7ffb9ab6541f4c4acb741352d5d7.rds")
mod <- readRDS("6ec4fdffcda16f774359e02da158958f.rds")
trace_plot(mod, "emiss_mu_bar", var=1)

# Read results
set_future_plan("multisession")
res <- parse_sleepsimR_results("results")

# Save as RDS file
saveRDS(res, "res_parsed.rds")

# Read allocations file
#alloc <- read_allocations_file("allocations.json")
# Get uids
#uids <- names(alloc)
#plan("future::multisession")
uids <- lapply(res, function(x) x$iteration_uid) %>%
  unlist()

# Read scenarios
scen <- read.csv("../make_scenarios/scenarios_subs.csv.gz",
                 stringsAsFactors = FALSE)
scen2 <- read.csv("../make_scenarios/scenarios.csv.gz",
                 stringsAsFactors = FALSE)
# Bind rows
scen <- bind_rows(scen, scen2)

# Get means
library(magrittr)
library(dplyr)
library(tidyr)
vmean <- res %>%
  get("emiss_mu_bar") %>%
  bind_rows() %>%
  select_at(vars(contains("mean_state")))
vmean$iteration_id <- res %>%
  get("iteration_uid") %>%
  unname() %>%
  unlist()
vmean_orig <- vmean
# Merge
vmean <- vmean %>%
  left_join(scen %>%
              select(iteration_id,
                     scenario_id))
vmean <- vmean %>%
  # Gather
  gather(var, val, -iteration_id, -scenario_id)
# Real values (bias)
rv <- data.frame(
  "var" = c("EEG_mean_beta_mean_state1", "EEG_mean_beta_mean_state2",
            "EEG_mean_beta_mean_state3", "EOG_median_theta_mean_state1",
            "EOG_median_theta_mean_state2", "EOG_median_theta_mean_state3",
            "EOG_min_beta_mean_state1", "EOG_min_beta_mean_state2",
            "EOG_min_beta_mean_state3"),
  "val" = c(-0.36, -0.6, 0.7, 1.01, -1.31, -.24, .75, -1.31, 0.005),
  stringsAsFactors = FALSE
)
# Summarize by group
library(purrr)
op <- vmean %>%
  split(list(vmean$scenario_id, vmean$var)) %>%
  map(function(x) {
    vari <- unique(x$var)
    gt <- rv %>% filter(var == vari) %>%
      select(val) %>% pull()
    val <- x$val
    bb <- unname(bias(gt, val))
    ese <- unname(emperical_SE(val))
    MSE <- MSE(val, gt)
    return(
      data.frame(
        "bias" = bb[1],
        "emp_se" = ese[1],
        "MSE" = MSE,
        #"MCMC_ERR" = bb[2],
        "scenario_id" = x$scenario_id[1],
        "var" = vari,
        stringsAsFactors = FALSE)
    )
  }) %>%
  do.call(rbind.data.frame, .)
# Get scenarios
op2 <- op %>%
  left_join(., scen %>%
              select(scenario_id, n, n_t, zeta, Q) %>%
              distinct())
library(ggplot2)
ggplot(op2 %>%
         filter(Q == 0.1), aes(x=n, y=bias, linetype=as.factor(n_t), color=as.factor(zeta))) +
  geom_point(size=1.7) +
  geom_line() +
  geom_hline(yintercept=0, color = "brown") +
  facet_wrap(". ~ var", scales="free_y")

# Plot distribution
op3 <- vmean %>%
  left_join(scen %>%
              select(scenario_id, n, n_t, zeta, Q) %>%
              distinct()) %>%
  left_join(rv, by ="var") %>%
  mutate(pbias = ((val.x - val.y)/val.y) * 100)
ggplot(op3 %>%
         filter(Q==0.1, zeta==1, n_t == 1600), aes(x=as.factor(n), y=pbias)) +
    geom_boxplot() +
    facet_wrap(". ~ var", scales = "free_y")

# Transition probabilities
vem <- res %>%
  get("gamma_prob_bar") %>%
  bind_rows()
# Add iteration id
vem$iteration_id <- res %>%
  get("iteration_uid") %>%
  unname() %>%
  unlist()
# Add scenario details
vem2 <- vem %>%
  left_join(scen %>%
              select(iteration_id,
                     scenario_id,
                     n, n_t, zeta, Q)) %>%
  # Gather
  gather(var, val, -iteration_id, -scenario_id,
         -n, -n_t, -Q, -zeta)
# Get ground-truth transition probabilities
library(sleepsimR)
opts <- getOption("sleepsimR_simulate")
tpm <- data.frame(
  "var" = unique(vem2$var),
  "val" = as.vector(opts$gamma_bar)
)
# Compute bias
b <- vem2 %>%
  split(list(vem2$scenario_id, vem2$var)) %>%
  map(function(x) {
    vari <- unique(x$var)
    gt <- tpm %>% filter(var == vari) %>%
      select(val) %>% pull()
    val <- x$val
    bb <- unname(bias(gt, val))
    ese <- unname(emperical_SE(val))
    MSE <- MSE(val, gt)
    return(
      data.frame(
        "bias" = bb[1],
        "emp_se" = ese[1],
        "MSE" = MSE,
        #"MCMC_ERR" = bb[2],
        "scenario_id" = x$scenario_id[1],
        "var" = vari,
        stringsAsFactors = FALSE)
    )
  }) %>%
  do.call(rbind.data.frame, .)
# Merge
op2 <- b %>%
  left_join(scen %>%
              select(scenario_id, n, n_t, Q, zeta) %>%
              distinct(), by="scenario_id")
library(ggplot2)
ggplot(op2 %>%
         filter(zeta == 0.5), aes(x=n, y=bias, linetype=as.factor(n_t), color=as.factor(Q))) +
  geom_point(size=1.7) +
  geom_line() +
  geom_hline(yintercept=0, color = "brown") +
  facet_wrap(". ~ var", scales="free_y")

# Coverage
cover <- res %>%
  get("credible_intervals")
cover <- cover %>%
  bind_rows()

# Emission distribution means
em_dist_means <- lapply(res, function(x) {
  ci <- x$credible_intervals$emiss_mu_bar
  nams <- names(ci)
  for(idx in seq_along(nams)) {
    names(ci[[idx]]) <- paste0("mean_",
                               c("state1_lower", "state1_upper",
                                 "state2_lower", "state2_upper",
                                 "state3_lower", "state3_upper"))
  }
  y <- do.call(cbind.data.frame, ci)
  colnames(y) <- gsub("\\.", "_", colnames(y))
  y
})
# Bind
em_dist_means <- do.call(rbind.data.frame, em_dist_means)
em_dist_means$iteration_id <- vmean_orig$iteration_id
em_dist_means <- em_dist_means %>%
  left_join(scen %>%
              select(scenario_id, iteration_id),
            by = "iteration_id")
# Get scenarios
vmean <- vmean_orig %>%
  left_join(scen %>%
              select(scenario_id, iteration_id),
            by = "iteration_id") %>%
  gather(var, val, -iteration_id, -scenario_id)
# Split by scenario
vmean_split <- split(vmean, list(vmean$scenario_id, vmean$var))
# Compute coverage
out_coverage <- vmean_split %>%
  map(., function(x) {
    ss <- em_dist_means %>%
      filter(scenario_id == x$scenario_id[1]) %>%
      select(starts_with(x$var[1]))
    ci <- lapply(1:nrow(ss), function(z) unname(unlist(ss[z,])))
    gt <- rv %>%
      filter(var == x$var[1]) %>%
      select(val) %>%
      pull()
    covv <- unname(coverage(ci, gt))
    data.frame(
      iteration_id = x$iteration_id[1],
      scenario_id = x$scenario_id[1],
      var = x$var[1],
      coverage = covv[1]
    )
  })
out_coverage <- do.call(rbind.data.frame, out_coverage)
# Combine with scenarios
out_coverage <- out_coverage %>%
  left_join(scen %>%
              select(scenario_id, n, n_t, Q, zeta) %>%
              distinct(),
            by = "scenario_id")
ggplot(out_coverage %>%
         filter(Q==0.4), aes(x=n, y=coverage, linetype=as.factor(n_t), color=as.factor(zeta))) +
  geom_point(size=1.7) +
  geom_line() +
  facet_wrap(". ~ var")

# Transition probabilities
tpm_cov <- lapply(res, function(x) {
  ci <- x$credible_intervals$gamma_prob_bar
  io <- c("int_S1toS1_lower", "int_S1toS1_upper",
           "int_S1toS2_lower", "int_S1toS2_upper",
           "int_S1toS3_lower", "int_S1toS3_upper",
           "int_S2toS1_lower", "int_S2toS1_upper",
           "int_S2toS2_lower", "int_S2toS2_upper",
           "int_S2toS3_lower", "int_S2toS3_upper",
           "int_S3toS1_lower", "int_S3toS1_upper",
           "int_S3toS2_lower", "int_S3toS2_upper",
           "int_S3toS3_lower", "int_S3toS3_upper")
  y <- data.frame(ci)
  colnames(y) <- io
  y
})
# Bind
tpm_cov <- do.call(rbind.data.frame, tpm_cov)
tpm_cov$iteration_id <- vmean_orig$iteration_id
tpm_cov <- tpm_cov %>%
  left_join(scen %>%
              select(scenario_id, iteration_id),
            by = "iteration_id")
# Get scenarios
vem2 <- vem %>%
  left_join(scen %>%
              select(scenario_id, iteration_id),
            by = "iteration_id") %>%
  gather(var, val, -iteration_id, -scenario_id)
# Split by scenario
vem2_split <- split(vem2, list(vem2$scenario_id, vem2$var))
# Compute coverage
out_coverage <- vem2_split %>%
  map(., function(x) {
    ss <- tpm_cov %>%
      filter(scenario_id == x$scenario_id[1]) %>%
      select(starts_with(x$var[1]))
    ci <- lapply(1:nrow(ss), function(z) unname(unlist(ss[z,])))
    gt <- tpm %>%
      filter(var == x$var[1]) %>%
      select(val) %>%
      pull()
    #gt <- mean(x$val)
    covv <- unname(coverage(ci, gt))
    data.frame(
      iteration_id = x$iteration_id[1],
      scenario_id = x$scenario_id[1],
      var = x$var[1],
      coverage = covv[1]
    )
  })
out_coverage <- do.call(rbind.data.frame, out_coverage)
# Combine with scenarios
out_coverage <- out_coverage %>%
  left_join(scen %>%
              select(scenario_id, n, n_t, Q, zeta) %>%
              distinct(),
            by = "scenario_id")
ggplot(out_coverage %>%
         filter(zeta==0.25), aes(x=n, y=coverage, linetype=as.factor(n_t), color=as.factor(Q))) +
  geom_point(size=2) +
  geom_line() +
  facet_wrap(". ~ var")

# Variance
var2 <- res %>%
  get("emiss_varmu_bar") %>%
  bind_rows() %>%
  select_at(vars(contains("mean_state")))
var2$iteration_id <- res %>%
  get("iteration_uid") %>%
  unname() %>%
  unlist() 
varm_orig <- var2
# Merge
var2 <- var2 %>%
  left_join(scen %>%
              select(iteration_id,
                     scenario_id))
var2 <- var2 %>%
  # Gather
  gather(var, val, -iteration_id, -scenario_id)
# Summarize by group
library(purrr)
op <- var2 %>%
  split(list(var2$scenario_id, var2$var))
op <- op %>%
  map(function(x) {
    vari <- unique(x$var)
    gt <- scen %>% filter(iteration_id == x$iteration_id[1]) %>%
      select(zeta) %>% pull()
    val <- x$val
    bb <- unname(bias(gt, val))
    ese <- unname(emperical_SE(val))
    MSE <- MSE(val, gt)
    return(
      data.frame(
        "bias" = bb[1],
        "emp_se" = ese[1],
        "MSE" = MSE,
        #"MCMC_ERR" = bb[2],
        "scenario_id" = x$scenario_id[1],
        "var" = vari,
        stringsAsFactors = FALSE)
    )
  }) %>%
  do.call(rbind.data.frame, .)
# Get scenarios
op2 <- op %>%
  left_join(., scen %>%
              select(scenario_id, n, n_t, zeta, Q) %>%
              distinct())
library(ggplot2)
ggplot(op2 %>%
         filter(Q == 0.4), aes(x=n, y=bias, linetype=as.factor(n_t), color=as.factor(zeta))) +
  geom_point(size=1.7) +
  geom_line() +
  facet_wrap(". ~ var")

# Variance coverage
var_cov <- lapply(res, function(x) {
  ci <- x$credible_intervals$emiss_varmu_bar
  nams <- names(ci)
  for(idx in seq_along(nams)) {
    names(ci[[idx]]) <- paste0("varmu_bar_",
                               c("state1_lower", "state1_upper",
                                 "state2_lower", "state2_upper",
                                 "state3_lower", "state3_upper"))
  }
  y <- do.call(cbind.data.frame, ci)
  colnames(y) <- gsub("\\.", "_", colnames(y))
  y
})
# Bind
var_cov <- do.call(rbind.data.frame, var_cov)
var_cov$iteration_id <- vmean_orig$iteration_id
var_cov2 <- var_cov %>%
  left_join(scen %>%
              select(scenario_id, iteration_id),
            by = "iteration_id")
# Get scenarios
varm <- varm_orig %>%
  left_join(scen %>%
              select(scenario_id, iteration_id),
            by = "iteration_id") %>%
  gather(var, val, -iteration_id, -scenario_id)
# Split by scenario
varm_split <- split(varm, list(varm$scenario_id,varm$var))
# Compute coverage
out_coverage <- varm_split %>%
  map(., function(x) {
    ss <- var_cov2 %>%
      filter(scenario_id == x$scenario_id[1]) %>%
      select(starts_with(x$var[1]))
    ci <- lapply(1:nrow(ss), function(z) unname(unlist(ss[z,])))
    gt <- tpm %>%
      filter(var == x$var[1]) %>%
      select(val) %>%
      pull()
    #gt <- mean(x$val)
    covv <- unname(coverage(ci, gt))
    data.frame(
      iteration_id = x$iteration_id[1],
      scenario_id = x$scenario_id[1],
      var = x$var[1],
      coverage = covv[1]
    )
  })
out_coverage <- do.call(rbind.data.frame, out_coverage)
# Combine with scenarios
out_coverage <- out_coverage %>%
  left_join(scen %>%
              select(scenario_id, n, n_t, Q, zeta) %>%
              distinct(),
            by = "scenario_id")
ggplot(out_coverage %>%
         filter(zeta==0.25), aes(x=n, y=coverage, linetype=as.factor(n_t), color=as.factor(Q))) +
  geom_point(size=2) +
  geom_line() +
  facet_wrap(". ~ var")
