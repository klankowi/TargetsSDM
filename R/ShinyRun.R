#####
## Shiny script
#####
library(tidyverse)
library(patchwork)
library(sf)
sf_use_s2(FALSE)

nice_times <- as.Date(c(paste0(seq(from = 1985, to = 2100, by = 1), "-03-16"), paste0(seq(from = 1985, to = 2100, by = 1), "-07-16"), paste0(seq(from = 1985, to = 2100, by = 1), "-10-16")))
nice_times <- nice_times[order(nice_times)]


land_sf<- st_read("/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/DFOSDM_app/shp/land.shp")
vast_proj_rds <- readRDS(paste0("/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/targets_output/", "prediction_df/American_lobster_GD2_BetaRW_EpsIIDRED_Temp_SSP5_85_mean_projections.rds"))
vast_fit <- readRDS(paste0("/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/targets_output/", "mod_fits/American_lobster_GD2_BetaRW_EpsIIDRED_Temp_fitted_vast.rds"))

# "Filters"
nice_species_name = "American_lobster"
nice_climate_name = "SSP85_mean"
  
# Summarizing baseline data
base_dat<- make_baseline_shiny(vast_proj_rds = vast_proj_rds, vast_fit = vast_fit, nice_times = nice_times, base_years = c(2015, 2019), nice_species_name = nice_species_name, nice_climate_name = nice_climate_name, out_dir = "/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/DFOSDM_app/baseline_data/")

# Summarizing projections
proj_dat<- make_projection_shiny(vast_proj_rds = vast_proj_rds, vast_fit = vast_fit, nice_times = nice_times, nice_species_name = nice_species_name, nice_climate_name = nice_climate_name, out_dir = "/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/DFOSDM_app/projection_data/")

# Calculate differences
proj_diffs<- make_vast_proj_diff(baseline_df = base_dat, proj_df = proj_dat, out_dir = "/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/DFOSDM_app/projection_data/")

# Plot baseline, future, and difference
base_fall<- plot_base_map(plot_dat_df = base_dat, plot_variable = "Log_Biomass", plot_season = "Fall", land_sf = land, xlim = c(-78.5, -54.5), ylim = c(34.5, 48.25), lab_lat = 36, lab_lon = -67.5, out_dir, land_color = "#d9d9d9", ...)
proj_fall<- plot_base_map(plot_dat_df = proj_dat, plot_variable = "Log_Biomass", plot_season = "Fall", land_sf = land, xlim = c(-78.5, -54.5), ylim = c(34.5, 48.25), lab_lat = 36, lab_lon = -67.5, out_dir, land_color = "#d9d9d9", ...)
prof_diff_fall<- plot_diff_map(plot_dat_df = proj_diffs, plot_variable = "Log_Biomass_Diff", plot_season = "Fall", plot_decade = 2050, land_sf = land, xlim = c(-78.5, -54.5), ylim = c(34.5, 48.25), lab_lat = 36, lab_lon = -67.5, out_dir, land_color = "#d9d9d9", ...)
base_fall + proj_fall + proj_diff_fall