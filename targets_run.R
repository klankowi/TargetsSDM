##########
##### Executing _targets.R
##########
library(targets)
library(parallel)
library(doFuture)
library(tictoc)
# 
cores_avail<- detectCores()
registerDoFuture()
plan(multisession, workers = cores_avail-2)

# Clean everything?
clean_start<- TRUE
if(clean_start){
  tar_destroy()
}

# Checking calls
tar_manifest(fields = "command")

# Graph
# tar_glimpse()
# tar_visnetwork(label = "time", reporter = "forecast", targets_only = TRUE)

# Run it
tic()
tar_make()
toc()

# Sync with remote...
system(paste("rclone copy ", "'/home/aallyn/results'", " 'TargetsSDM:Work/TargetsSDM/results'", sep = ""))
system(paste("rclone copy ", "'/home/aallyn/data'", " 'TargetsSDM:Work/TargetsSDM/data'", sep = ""))


# Check on warnings
warning_ind<- which(!is.na(tar_meta(fields = warnings)$warnings))
tar_meta(fields = warnings)[warning_ind, ]

tar_load(vast_fit)
TMBhelper::check_estimability(vast_fit$tmb_list$Obj)


