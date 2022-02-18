##########
##### Executing _targets.R
##########
library(devtools)
# library(gmRi)
library(Matrix)
library(VAST)
library(FishStatsUtils)
library(tidyverse)
library(lubridate)
library(sf)
library(raster)
library(gstat)
library(here)
library(targets)
library(tarchetypes)
library(patchwork)
library(gifski)
library(akima)
library(splines)
library(parallel)
library(doFuture)
library(tools)
library(tictoc)

# Targets set up
options(tidyverse.quiet = TRUE)
tar_option_reset()
tar_option_set(packages = c("Matrix", "TMB", "FishStatsUtils", "VAST", "tidyverse", "lubridate", "sf", "raster", "here", "tools"))

# Clean everything?
clean_start <- TRUE
if (clean_start) {
  tar_destroy()
}

# Trying to leverage more compute cores
cores_avail <- detectCores()
registerDoFuture()
plan(multisession, workers = cores_avail - 2)

# Trying to leverage more compute cores
# Checking calls
tar_manifest(fields = "command", callr_function = NULL)

# Graph
# tar_glimpse()
# tar_visnetwork(label = "time", reporter = "forecast", targets_only = TRUE)

# Run it
tic()
tar_make()
toc()
