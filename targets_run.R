##########
##### Executing _targets.R
##########
library(targets)
library(parallel)
library(doFuture)
library(tictoc)

# Targets set up
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("Matrix", "TMB", "FishStatsUtils", "VAST", "tidyverse", "lubridate", "sf", "raster", "here", "tools"))

# Clean everything?
clean_start<- TRUE
if(clean_start){
  tar_destroy()
}

# Trying to leverage more compute cores
# cores_avail<- detectCores()
# registerDoFuture()
# plan(multisession, workers = cores_avail-2)

# Checking calls
tar_manifest(fields = "command")

# Graph
# tar_glimpse()
# tar_visnetwork(label = "time", reporter = "forecast", targets_only = TRUE)

# Run it
tic()
tar_make()
toc()


