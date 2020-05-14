## Compute sums of squares for each parameter
library(sleepsimReval)
library(sleepsimR)
library(dplyr)
library(tidyr)

# Analysis of the model history
rm(list=ls())
hist <- readRDS("history.rds") %>%
  do.call(rbind.data.frame, .)
hist[,7:30] <- apply(hist[,7:30], 2, function(x) ifelse(x == "Yes", TRUE, FALSE))

# TOtal number of models converged
alconv <- rowSums(hist[,7:30]) == 24
sum(alconv) / length(alconv) # ~ 32% have converged

# By parameter
conv_by_par <- hist[,c(1, 7:30)] %>%
  gather(var,val, -iteration_id) %>%
  group_by(var) %>%
  summarize(pconverged = sum(val) / n()) %>%
  arrange(var)

conv_by_par <- hist[,c(1, 7:30)] %>%
  gather(var,val, -iteration_id) %>%
  group_by(iteration_id) %>%
  summarize(num_params_not_converged = 24 - sum(val)) %>%
  ungroup() %>%
  group_by(num_params_not_converged) %>%
  tally()

hist(conv_by_par$num_params_not_converged)

# Across zeta
conv_by_par <- hist %>%
  select(-converged) %>%
  gather(var,val, -iteration_id, -n, -n_t, -Q, -zeta) %>%
  group_by(var, zeta) %>%
  summarize(pconverged = sum(val) / n()) %>%
  arrange(var, zeta) %>%
  arrange(pconverged)

# Get file names
f <- list.files("data/rerun")
files_sim <- file.path("data/simulations", f)
files_rerun <- file.path("data/rerun", f)

# Compute sums of squares
sse <- function(x, y) {
  sqrt(mean((x - y)^2))
}

# For results
m1 <- m2 <- m3 <- matrix(0, nrow = length(f),
            ncol = (6 + 9 + 9))


# Load
for(idx in 1:length(f)) {
  # Load models
  f1 <- readRDS(files_sim[idx])
  f2 <- readRDS(files_rerun[idx])
  # Burn 
  f1b <- burn(f1)
  f2b <- burn(f2)
  # Compute RMSE
  for(col_idx in 1:6) {
    m1[idx, col_idx] <- sse(f1b$gamma_int_bar[,col_idx], 
                           f2b$gamma_int_bar[,col_idx])
    m2[idx, col_idx] <- sse(f1b$gamma_int_bar[1:1000,col_idx], 
                           f2b$gamma_int_bar[1:1000,col_idx])
    m3[idx, col_idx] <- sse(f1b$gamma_int_bar[1000:2000,col_idx], 
                         f2b$gamma_int_bar[1000:2000,col_idx])
  }
  for(var in c("EEG_mean_beta", "EOG_median_theta", "EOG_min_beta")) {
    for(col_idx in 1:3) {
      m[idx, col_idx] <- sse(f1b$emiss_mu_bar[[var]][,col_idx], 
                             f2b$emiss_mu_bar[[var]][,col_idx])
      m[idx, col_idx] <- sse(f1b$emiss_mu_bar[[var]][1:1000,col_idx], 
                                           f2b$emiss_mu_bar[[var]][1:1000,col_idx])
      m[idx, col_idx] <- sse(f1b$emiss_mu_bar[[var]][1000:2000,col_idx], 
                                             f2b$emiss_mu_bar[[var]][1000:2000,col_idx])
      m[idx, col_idx] <- sse(f1b$emiss_mu_bar[[var]][,col_idx], 
                                 f2b$emiss_mu_bar[[var]][,col_idx])
      m[idx, col_idx] <- sse(f1b$emiss_mu_bar[[var]][1:1000,col_idx], 
                                               f2b$emiss_mu_bar[[var]][1:1000,col_idx])
      m[idx, col_idx] <- sse(f1b$emiss_mu_bar[[var]][1000:2000,col_idx], 
                                                 f2b$emiss_mu_bar[[var]][1000:2000,col_idx])    
    }
  }
}

sleepsimReval::tpp(f1, f2, "emiss_varmu_bar", 2)



RMSE_out_1 <- c()
RMSE_half_1 <- c()
RMSE_half_2 <- c()

for(col_idx in 1:6) {
  RMSE_out_1 <- c(RMSE_out_1,
                  sse(f1b$gamma_int_bar[,col_idx], 
                      f2b$gamma_int_bar[,col_idx]))
  RMSE_half_1 <- c(RMSE_half_1,
                  sse(f1b$gamma_int_bar[1:1000,col_idx], 
                      f2b$gamma_int_bar[1:1000,col_idx]))
  RMSE_half_2 <- c(RMSE_half_2,
                  sse(f1b$gamma_int_bar[1000:2000,col_idx], 
                      f2b$gamma_int_bar[1000:2000,col_idx]))
}

RMSE_out_1 <- c()
RMSE_half_1 <- c()
RMSE_half_2 <- c()
for(col_idx in 1:3) {
  RMSE_out_1 <- c(RMSE_out_1,
                  sse(f1b$emiss_varmu_bar$EOG_median_theta[,col_idx], 
                      f2b$emiss_varmu_bar$EOG_median_theta[,col_idx]))
  RMSE_half_1 <- c(RMSE_half_1,
                   sse(f1b$emiss_varmu_bar$EOG_median_theta[1:1000,col_idx], 
                       f2b$emiss_varmu_bar$EOG_median_theta[1:1000,col_idx]))
  RMSE_half_2 <- c(RMSE_half_2,
                   sse(f1b$emiss_varmu_bar$EOG_median_theta[1000:2000,col_idx], 
                       f2b$emiss_varmu_bar$EOG_median_theta[1000:2000,col_idx]))
}
RMSE_out_1
RMSE_half_1
RMSE_half_2


