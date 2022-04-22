#####
## Shiny script
#####
library(tidyverse)
library(patchwork)
library(sf)
library(lubridate)
sf_use_s2(FALSE)

nice_times <- as.Date(c(paste0(seq(from = 1985, to = 2100, by = 1), "-03-16"), paste0(seq(from = 1985, to = 2100, by = 1), "-07-16"), paste0(seq(from = 1985, to = 2100, by = 1), "-10-16")))
nice_times <- nice_times[order(nice_times)]

# "Filters"
nice_species_name = "Atlantic_halibut_ST"
nice_climate_name = "SSP5_85_mean"

proj_file<- paste0(nice_species_name, "_", nice_climate_name, "_projections.rds")
fit_file<- paste0(nice_species_name, "_fitted_vast.rds")

land<- st_read("/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/DFOSDM_app/shp/land.shp")
vast_proj_rds <- readRDS(paste0("/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/targets_output/", "prediction_df/", proj_file))
vast_fit <- readRDS(paste0("/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/targets_output/", "mod_fits/", fit_file))

spatial_areas<- st_read("/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/DFOSDM_app/shp/spatial_areas.shp")



# Summarizing baseline data
base_dat<- make_baseline_shiny(vast_proj_rds = vast_proj_rds, vast_fit = vast_fit, nice_times = nice_times, base_years = c(2015, 2019), nice_species_name = nice_species_name, nice_climate_name = nice_climate_name, in_crs = 4326, out_crs = 32619, out_dir = "/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/DFOSDM_app/baseline_data/")
#base_dat<- make_baseline_shiny(vast_proj_rds = vast_proj_rds, vast_fit = vast_fit, nice_times = nice_times, base_years = c(2015, 2019), nice_species_name = nice_species_name, nice_climate_name = nice_climate_name, in_crs = 4326, out_crs = 32619, out_dir = "/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/DFOSDM_app/baseline_no_ST/")

# Summarizing projections
proj_dat<- make_projection_shiny(vast_proj_rds = vast_proj_rds, vast_fit = vast_fit, nice_times = nice_times, nice_species_name = nice_species_name, nice_climate_name = nice_climate_name, in_crs = 4326, out_crs = 32619, out_dir = "/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/DFOSDM_app/projected_data/")
#proj_dat<- make_projection_shiny(vast_proj_rds = vast_proj_rds, vast_fit = vast_fit, nice_times = nice_times, nice_species_name = nice_species_name, nice_climate_name = nice_climate_name, in_crs = 4326, out_crs = 32619, out_dir = "/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/DFOSDM_app/projected_no_ST/")

# Calculate differences
proj_diff<- make_proj_diff(baseline_df = base_dat, proj_df = proj_dat, in_crs = 32619, out_crs = 32619, out_dir = "/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/DFOSDM_app/biomass_diff_data/")
#proj_diff<- make_proj_diff(baseline_df = base_dat, proj_df = proj_dat, in_crs = 32619, out_crs = 32619, out_dir = "/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/DFOSDM_app/biomass_diff_data_no_ST/")

# # Translate to sf?
# base_utm<- df_to_sf(df = base_dat, in_crs = 4326, out_crs = 32619, out_filename = "~/Desktop/temp_baseline_utm.geojson")
# proj_utm<- df_to_sf(df = proj_dat, in_crs = 4326, out_crs = 32619, out_filename = "~/Desktop/temp_proj_utm.geojson")
# diff_utm<- df_to_sf(df = proj_diff, in_crs = 4326, out_crs = 32619, out_filename = "~/Desktop/temp_proj_diff_utm.geojson")
# 
# # Plot baseline, future, and difference
# base_plot<- plot_base_map(plot_dat = base_dat, plot_variable = "Log_Biomass", plot_season = "Fall", land_sf = land, pred_res = 25000, xlim = c(-78.5, -54.5), ylim = c(34.5, 48.25), lab_lat = 36, lab_lon = -67.5, out_filename = "~/Desktop/temp_base.jpg", land_color = "#d9d9d9")
# proj_plot<- plot_proj_map(plot_dat = proj_dat, plot_variable = "Log_Biomass", plot_season = "Fall", plot_decade = 2050, land_sf = land, pred_res = 25000, xlim = c(-78.5, -54.5), ylim = c(34.5, 48.25), lab_lat = 36, lab_lon = -67.5, out_filename = "~/Desktop/temp_2050.jpg", land_color = "#d9d9d9")
# proj_diff_plot<- plot_diff_map(plot_dat = proj_diff, plot_variable = "Log_Biomass_Diff", plot_season = "Fall", plot_decade = 2050, land_sf = land, pred_res = 25000, xlim = c(-78.5, -54.5), ylim = c(34.5, 48.25), lab_lat = 36, lab_lon = -67.5, out_filename = "~/Desktop/temp_diff.jpg", land_color = "#d9d9d9")
# base_plot + proj_plot + proj_diff_plot

# COG and plot
cog_dat<- make_cog_shiny(vast_proj_rds = vast_proj_rds, vast_fit = vast_fit, nice_times = nice_times, nice_species_name = nice_species_name, nice_climate_name = nice_climate_name, out_dir = "/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/DFOSDM_app/center_biomass/")
#cog_dat<- make_cog_shiny(vast_proj_rds = vast_proj_rds, vast_fit = vast_fit, nice_times = nice_times, nice_species_name = nice_species_name, nice_climate_name = nice_climate_name, out_dir = "/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/DFOSDM_app/center_biomass_no_ST/")
# cog_plot<- plot_cog(cog_dat = cog_dat, plot_season = "All", land_sf = land, xlim = c(-78.5, -54.5), ylim = c(34.5, 48.25), out_filename = "/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/DFOSDM_app/center_biomass/American_lobster_SSP85_mean_All_cog.jpg", land_color = "#d9d9d9")

# Index summary and plot
biomass_dat<- make_spatial_index(vast_proj_rds = vast_proj_rds, vast_fit = vast_fit, nice_times = nice_times, nice_species_name = nice_species_name, nice_climate_name = nice_climate_name, out_dir = "/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/DFOSDM_app/spatial_summary/")
#biomass_dat<- make_spatial_index(vast_proj_rds = vast_proj_rds, vast_fit = vast_fit, nice_times = nice_times, nice_species_name = nice_species_name, nice_climate_name = nice_climate_name, out_dir = "/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/DFOSDM_app/spatial_summary_no_ST/")
# bio_plot<- plot_biomass_index(biomass_dat = biomass_dat, plot_season = "All", plot_variable = "Log_Biomass", plot_regions = c("DFO", "NMFS", "GoM", "NMFS_and_DFO", "SNE_and_MAB"), nice_xlab = "Year", nice_ylab = "Total log biomass (kg)", color_pal = c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00'), out_filename = "~/Desktop/temp_log_biomass.jpg")

