library(dplyr)
# Helper used to load models from disk
read_models <- function(file_name) {
  file_sim <- file.path("data/simulations", file_name)
  file_rer <- file.path("data/rerun", file_name)
  # Read and return
  list(
    readRDS(file_sim),
    readRDS(file_rer)
  )
}
# Load scenarios
scen <- read.csv("scenarios.csv.gz")
# Get number of files
files_orig <- list.files("data/simulations")
files_rerun <- list.files("data/rerun")
files_orig <- files_rerun <- intersect(files_orig, files_rerun)
# Number of files
n_o <- length(files_orig)
n_r <- length(files_rerun)
n <- n_o
# Read from history if available
if("history.rds" %in% list.files(".")) {
  history <- readRDS("history.rds")
} else {
  history <- list()
}
# Get names from the list and determine which files still need to be checked
nams_processed <- names(history)
# Setdiff
needs_processing <- setdiff(files_orig, nams_processed)
# Pick a file
f <- needs_processing[1]
# Read
mods <- read_models(f)
# Get model details
scen_cur <- scen %>%
  filter(iteration_id == gsub("\\.rds", "", f))
# Parameter settings
n <- scen_cur$n
n_t <- scen_cur$n_t
zeta <- scen_cur$zeta
Q <- scen_cur$Q
