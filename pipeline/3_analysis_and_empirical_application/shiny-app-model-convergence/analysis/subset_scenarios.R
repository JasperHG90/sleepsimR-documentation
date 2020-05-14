## Subset scenarios
# Global options go here
dd <- read.csv("scenarios.csv.gz")
library(dplyr)
itid <- dd %>% group_by(scenario_id) %>%
  slice(1:3) %>%
  select(iteration_id) %>%
  ungroup() %>%
  pull()

# Download reruns
cmd1 <- "gsutil -m cp -r gs://sleepsimr-baseline/sleepsimrmodels/*.rds data/rerun"

# Download originals
cmd2 <- paste0(
  "gsutil -m cp -r gs://sleepsimr-models/models/simulations/{",
  paste0(
    itid, collapse=","
  ),
  "}.rds data/simulations"
)
