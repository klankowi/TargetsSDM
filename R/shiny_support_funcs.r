######
## Collection of functions to support shiny app presentation
######
library(tidyverse)
library(VAST)
library(akima)
library(sf)
sf_use_s2(FALSE)

## Translate to sf function
df_to_sf<- function(df, in_crs, out_crs, out_filename){
  if(FALSE){
    df<- read.csv("/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/DFOSDM_app/baseline_data/American_lobster_SSP85_mean_baseline.csv")
    in_crs<- 4326
    out_crs<- 4326
    out_filename<- "/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/DFOSDM_app/baseline_data/American_lobster_SSP85_mean_baseline_wgs84.shp"
  }
  
  sf_temp<- df %>%
    st_as_sf(., coords = c("Lon", "Lat"), crs = in_crs, remove = FALSE)
  
  if(out_crs != in_crs){
    sf_out<- sf_temp %>%
      st_transform(., crs = out_crs)
  } else {
    sf_out<- sf_temp
  }
  
  # Write and return it
  st_write(sf_out, out_filename, append = FALSE)
  return(sf_out)
}

## Make baseline data function -- for this, we want to have columns for x, y, species, season (spring, summer, fall, all), climate scenario, biomass and log biomass
make_baseline_shiny <- function(vast_proj_rds, vast_fit, nice_times, base_years, nice_species_name, nice_climate_name, in_crs, out_crs, out_dir) {
    if (FALSE) {
        vast_proj_rds <- readRDS(paste0("/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/targets_output/", "prediction_df/American_lobster_GD2_BetaRW_EpsIIDRED_Temp_SSP5_85_mean_projections.rds"))
        vast_fit <- readRDS(paste0("/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/targets_output/", "mod_fits/American_lobster_GD2_BetaRW_EpsIIDRED_Temp_fitted_vast.rds"))
        base_years <- c(2015, 2019)
        nice_species_name <- "American_lobster"
        nice_climate_name <- "SSP85_mean"
        in_crs = 4326
        out_crs = 32619
        out_dir <- "/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/DFOSDM_app/baseline_data"
    }

    # Some general stuff
    # Time series steps
    time_ind <- seq(from = 1, to = vast_proj_rds[[1]]$n_t)
    time_labels <- nice_times

    # Categories
    categories_ind <- seq(from = 1, to = vast_proj_rds[[1]]$n_c)

    # Grid locations
    grid_ind <- seq(from = 1, to = vast_proj_rds[[1]]$n_g)

    for (i in seq_along(vast_proj_rds)) {
        temp_array <- array(vast_proj_rds[[i]][which(names(vast_proj_rds[[i]]) == "D_gct")][[1]], dim = c(unlist(vast_proj_rds[[i]][c("n_g", "n_c", "n_t")])), dimnames = list(grid_ind, categories_ind, time_labels))

        temp_df <- data.frame(aperm(temp_array, c(1, 3, 2)))
        colnames(temp_df) <- nice_times
        temp_df$Lat <- vast_fit$spatial_list$latlon_g[, "Lat"]
        temp_df$Lon <- vast_fit$spatial_list$latlon_g[, "Lon"]

        temp_df <- temp_df %>%
            distinct() %>%
            pivot_longer(., !c(Lat, Lon), names_to = "Time", values_to = "D_gct") %>%
            arrange(Time, Lat, Lon)

        temp_df$Sim_Scenario <- paste0("Sim_", i)

        if (i == 1) {
            temp_out <- temp_df
        } else {
            temp_out <- bind_rows(temp_out, temp_df)
        }
    }

    # Gather up years to filter the data
    min_year <- base_years[1]
    max_year <- base_years[2]

    res_out <- temp_out %>%
        mutate(., Date = as.Date(Time)) %>%
        filter(., format(Date, "%Y") >= min_year & format(Date, "%Y") <= max_year) %>%
        mutate(
            "Year" = format(Date, "%Y"),
            "Month" = format(Date, "%m"),
            "Season" = ifelse(Month == "03", "Spring", ifelse(Month == "10", "Fall", "Summer"))
        )
    
    res_out_seas<- res_out %>%
      group_by(., Lon, Lat, Season) %>%
      summarize_at(., vars(D_gct), list(Biomass = mean), na.rm = FALSE) %>%
      mutate(., "Log_Biomass" = log(Biomass + 1)) %>%
      arrange(., Season, Lon, Lat)
    
    res_out_all<- res_out %>%
      group_by(., Lon, Lat) %>%
      summarize_at(., vars(D_gct), list(Biomass = mean), na.rm = FALSE) %>%
      mutate(., "Log_Biomass" = log(Biomass + 1),
             "Season" = "All") %>%
      arrange(., Season, Lon, Lat) 
    
    res_out_final<- bind_rows(res_out_seas, res_out_all)
    res_out_final$Season<- factor(res_out_final$Season, levels = c("Spring", "Summer", "Fall", "All"), labels = c("Spring", "Summer", "Fall", "All"))
    res_out_final$Species<- nice_species_name
    res_out_final$Climate_Scenario<- nice_climate_name
    
    if(in_crs != out_crs){
      res_out_final<- st_as_sf(res_out_final, coords = c("Lon", "Lat"), crs = in_crs) %>%
        st_transform(., crs = out_crs) %>%
        mutate(., "Lon" = st_coordinates(.)[,1],
               "Lat" = st_coordinates(.)[,2]) %>%
        st_drop_geometry() %>%
        dplyr::select(., Lon, Lat, Season, Climate_Scenario, Species, Biomass, Log_Biomass)
    }
    
    # Save it out
    write.csv(res_out_final, file = paste0(out_dir, paste(nice_species_name, nice_climate_name, "baseline", sep = "_"), ".csv"))
    return(res_out_final)
}

## Make center of gravity data function -- for this, we want to have columns for x, y, season, year
make_cog_shiny<- function(vast_proj_rds, vast_fit, nice_times, nice_species_name, nice_climate_name, out_dir){
  if (FALSE) {
    vast_proj_rds <- readRDS(paste0("/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/targets_output/", "prediction_df/American_lobster_GD2_BetaRW_EpsIIDRED_Temp_SSP5_85_mean_projections.rds"))
    vast_fit <- readRDS(paste0("/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/targets_output/", "mod_fits/American_lobster_GD2_BetaRW_EpsIIDRED_Temp_fitted_vast.rds"))
    nice_species_name <- "American_lobster"
    nice_climate_name <- "SSP85_mean"
    out_dir <- "/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/DFOSDM_app/center_biomass"
  }
  
  # Some general stuff
  # Time series steps
  time_ind <- seq(from = 1, to = vast_proj_rds[[1]]$n_t)
  time_labels <- nice_times
  
  # Categories
  categories_ind <- seq(from = 1, to = vast_proj_rds[[1]]$n_c)
  
  # Grid locations
  grid_ind <- seq(from = 1, to = vast_proj_rds[[1]]$n_g)
  
  for (i in seq_along(vast_proj_rds)) {
    temp_array <- array(vast_proj_rds[[i]][which(names(vast_proj_rds[[i]]) == "D_gct")][[1]], dim = c(unlist(vast_proj_rds[[i]][c("n_g", "n_c", "n_t")])), dimnames = list(grid_ind, categories_ind, time_labels))
    
    temp_df <- data.frame(aperm(temp_array, c(1, 3, 2)))
    colnames(temp_df) <- nice_times
    temp_df$Lat <- vast_fit$spatial_list$latlon_g[, "Lat"]
    temp_df$Lon <- vast_fit$spatial_list$latlon_g[, "Lon"]
    
    temp_df <- temp_df %>%
      distinct() %>%
      pivot_longer(., !c(Lat, Lon), names_to = "Time", values_to = "D_gct") %>%
      arrange(Time, Lat, Lon)
    
    temp_df$Sim_Scenario <- paste0("Sim_", i)
    
    if (i == 1) {
      temp_out <- temp_df
    } else {
      temp_out <- bind_rows(temp_out, temp_df)
    }
  }
  
  res_out <- temp_out %>%
    mutate(., Date = as.Date(Time)) %>%
    mutate(
      "Year" = format(Date, "%Y"),
      "Month" = format(Date, "%m"),
      "Season" = ifelse(Month == "03", "Spring", ifelse(Month == "10", "Fall", "Summer"))
    )
  
  res_out_seas<- res_out %>%
    group_by(., Season, Year) %>%
    nest()
  
  get_cog_df<- function(df) {
    if(FALSE){
      df<- res_out_seas$data[[1]]
    }
    
    cog_lat<- weighted.mean(df[,"Lat"], w = df[, "D_gct"])
    cog_lon<- weighted.mean(df[,"Lon"], w = df[, "D_gct"])
    
    cog_out<- data.frame("Lon" = cog_lon, "Lat" = cog_lat)
    return(cog_out)
  }
  
  res_out_seas<- res_out_seas %>% 
    mutate(., "COG" = map(data, get_cog_df))
  
  res_out_all<- res_out %>%
    group_by(., Year) %>%
    nest() %>%
    mutate(., "COG" = map(data, get_cog_df),
           "Season" = "All")
  
  res_out_final<- bind_rows(res_out_seas, res_out_all)
  res_out_final$Season<- factor(res_out_final$Season, levels = c("Spring", "Summer", "Fall", "All"), labels = c("Spring", "Summer", "Fall", "All"))
  res_out_final$Species<- nice_species_name
  res_out_final$Climate_Scenario<- nice_climate_name
  
  res_out_final<- res_out_final %>%
    dplyr::select(., Year, Season, Species, Climate_Scenario, COG) %>%
    unnest(., cols = c(COG))
  
  res_out_final$Year<- as.numeric(res_out_final$Year)
  
  # Save it out
  write.csv(res_out_final, file = paste0(out_dir, paste(nice_species_name, nice_climate_name, "cog", sep = "_"), ".csv"))
  return(res_out_final)
}

## Make projection data -- for this, we want to have x, y, decade date, species, climate scenario, biomass and log biomass
make_projection_shiny<- function(vast_proj_rds, vast_fit, nice_times, nice_species_name, nice_climate_name, in_crs, out_crs, out_dir){
  if (FALSE) {
  vast_proj_rds <- readRDS(paste0("/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/targets_output/", "prediction_df/American_lobster_GD2_BetaRW_EpsIIDRED_Temp_SSP5_85_mean_projections.rds"))
  vast_fit <- readRDS(paste0("/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/targets_output/", "mod_fits/American_lobster_GD2_BetaRW_EpsIIDRED_Temp_fitted_vast.rds"))
  nice_species_name <- "American_lobster"
  nice_climate_name <- "SSP85_mean"
  in_crs = 4326
  out_crs = 32619
  out_dir <- "/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/DFOSDM_app/projected_data"
  }
  
  # Some general stuff
  # Time series steps
  time_ind <- seq(from = 1, to = vast_proj_rds[[1]]$n_t)
  time_labels <- nice_times
  
  # Categories
  categories_ind <- seq(from = 1, to = vast_proj_rds[[1]]$n_c)
  
  # Grid locations
  grid_ind <- seq(from = 1, to = vast_proj_rds[[1]]$n_g)
  
  for (i in seq_along(vast_proj_rds)) {
    temp_array <- array(vast_proj_rds[[i]][which(names(vast_proj_rds[[i]]) == "D_gct")][[1]], dim = c(unlist(vast_proj_rds[[i]][c("n_g", "n_c", "n_t")])), dimnames = list(grid_ind, categories_ind, time_labels))
    
    temp_df <- data.frame(aperm(temp_array, c(1, 3, 2)))
    colnames(temp_df) <- nice_times
    temp_df$Lat <- vast_fit$spatial_list$latlon_g[, "Lat"]
    temp_df$Lon <- vast_fit$spatial_list$latlon_g[, "Lon"]
    
    temp_df <- temp_df %>%
      distinct() %>%
      pivot_longer(., !c(Lat, Lon), names_to = "Time", values_to = "D_gct") %>%
      arrange(Time, Lat, Lon)
    
    temp_df$Sim_Scenario <- paste0("Sim_", i)
    
    if (i == 1) {
      temp_out <- temp_df
    } else {
      temp_out <- bind_rows(temp_out, temp_df)
    }
  }
  
  # Getting date and season
  res_out <- temp_out %>%
    mutate(., Date = as.Date(Time)) %>%
    mutate(
      "Year" = as.numeric(format(Date, "%Y")),
      "Month" = format(Date, "%m"),
      "Season" = ifelse(Month == "03", "Spring", ifelse(Month == "10", "Fall", "Summer"))
    )
  
  res_out$Decade<- as.numeric(format(floor_date(res_out$Date, years(10)), "%Y"))
  
  res_out<- res_out %>%
    filter(., Year >= 2020 & Year <= 2099)
  
  res_out_seas<- res_out %>%
    group_by(., Lon, Lat, Season, Decade) %>%
    summarize_at(., vars(D_gct), list(Biomass = mean), na.rm = FALSE) %>%
    mutate(., "Log_Biomass" = log(Biomass + 1)) %>%
    arrange(., Season, Decade, Lon, Lat)
  
  res_out_all<- res_out %>%
    group_by(., Lon, Lat, Decade) %>%
    summarize_at(., vars(D_gct), list(Biomass = mean), na.rm = FALSE) %>%
    mutate(., "Log_Biomass" = log(Biomass + 1),
           "Season" = "All") %>%
    arrange(., Season, Decade, Lon, Lat) 
  
  res_out_final<- bind_rows(res_out_seas, res_out_all)
  res_out_final$Season<- factor(res_out_final$Season, levels = c("Spring", "Summer", "Fall", "All"), labels = c("Spring", "Summer", "Fall", "All"))
  res_out_final$Species<- nice_species_name
  res_out_final$Climate_Scenario<- nice_climate_name
  
  if(in_crs != out_crs){
    res_out_final<- st_as_sf(res_out_final, coords = c("Lon", "Lat"), crs = in_crs) %>%
      st_transform(., crs = out_crs) %>%
      ungroup() %>%
      mutate(., "Lon" = st_coordinates(.)[,1],
             "Lat" = st_coordinates(.)[,2]) %>%
      st_drop_geometry() %>%
      dplyr::select(., Lon, Lat, Season, Decade, Climate_Scenario, Species, Biomass, Log_Biomass)
  }
  
  # Save it out
  write.csv(res_out_final, file = paste0(out_dir, paste(nice_species_name, nice_climate_name, "projected", sep = "_"), ".csv"))
  return(res_out_final)
}

## Difference projections from baseline
make_proj_diff<- function(baseline_df, proj_df, in_crs, out_crs, out_dir){
  if(FALSE){
    baseline_df<- base_dat
    proj_df<- proj_dat
    in_crs = 32619
    out_crs = 32619
  }
  
  # Join
  res_out_final<- proj_df %>%
    left_join(., baseline_df, by = c("Lon" = "Lon", "Lat" = "Lat", "Season" = "Season", "Species" = "Species", "Climate_Scenario" = "Climate_Scenario")) %>%
    mutate(., "Biomass_Diff" = Biomass.x - Biomass.y,
           "Log_Biomass_Diff" = Log_Biomass.x - Log_Biomass.y) %>%
    dplyr::select(., Lon, Lat, Season, Decade, Species, Climate_Scenario, Biomass_Diff, Log_Biomass_Diff)
  
  if(in_crs != out_crs){
    res_out_final<- st_as_sf(res_out_final, coords = c("Lon", "Lat"), crs = in_crs) %>%
      st_transform(., crs = out_crs) %>%
      ungroup() %>%
      mutate(., "Lon" = st_coordinates(.)[,1],
             "Lat" = st_coordinates(.)[,2]) %>%
      st_drop_geometry() %>%
      dplyr::select(., Lon, Lat, Season, Decade, Climate_Scenario, Species, Biomass_Diff, Log_Biomass_Diff)
  }
  
  write.csv(res_out_final, file = paste0(out_dir, paste(nice_species_name, nice_climate_name, "projected_diff", sep = "_"), ".csv"))
  return(res_out_final)
  
}

## Make spatial summaries -- for this, we want to have species, season, year, scenario, mean, sd, date, climate scenario
make_spatial_summ<- function(vast_proj_rds, vast_fit, spatial_areas, nice_times, nice_species_name, nice_climate_name, out_dir){
  if(FALSE){
    vast_proj_rds <- readRDS(paste0("/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/targets_output/", "prediction_df/American_lobster_GD2_BetaRW_EpsIIDRED_Temp_SSP5_85_mean_projections.rds"))
    vast_fit <- readRDS(paste0("/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/targets_output/", "mod_fits/American_lobster_GD2_BetaRW_EpsIIDRED_Temp_fitted_vast.rds"))
    spatial_areas<- st_read("/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/DFOSDM_app/shp/spatial_areas.shp")
    nice_species_name <- "American_lobster"
    nice_climate_name <- "SSP85_mean"
    out_dir <- "/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/DFOSDM_app/spatial_summary"
  }

  # Some general stuff
  # Time series steps
  time_ind <- seq(from = 1, to = vast_proj_rds[[1]]$n_t)
  time_labels <- nice_times
  
  # Categories
  categories_ind <- seq(from = 1, to = vast_proj_rds[[1]]$n_c)
  
  # Grid locations
  grid_ind <- seq(from = 1, to = vast_proj_rds[[1]]$n_g)
  
  for (i in seq_along(vast_proj_rds)) {
    temp_array <- array(vast_proj_rds[[i]][which(names(vast_proj_rds[[i]]) == "D_gct")][[1]], dim = c(unlist(vast_proj_rds[[i]][c("n_g", "n_c", "n_t")])), dimnames = list(grid_ind, categories_ind, time_labels))
    
    temp_df <- data.frame(aperm(temp_array, c(1, 3, 2)))
    colnames(temp_df) <- nice_times
    temp_df$Lat <- vast_fit$spatial_list$latlon_g[, "Lat"]
    temp_df$Lon <- vast_fit$spatial_list$latlon_g[, "Lon"]
    
    temp_df <- temp_df %>%
      distinct() %>%
      pivot_longer(., !c(Lat, Lon), names_to = "Time", values_to = "D_gct") %>%
      arrange(Time, Lat, Lon)
    
    temp_df$Sim_Scenario <- paste0("Sim_", i)
    
    temp_df$Biomass<- temp_df$D_gct * 625
    
    if (i == 1) {
      temp_out <- temp_df
    } else {
      temp_out <- bind_rows(temp_out, temp_df)
    }
  }
  
  # Convert prediction df to spatial dataframe...
  pred_sp<- st_as_sf(temp_out, coords = c("Lon", "Lat"), crs = 4326, remove = FALSE)
  
  # Join up with the spatial areas after checking crs
  if(st_crs(pred_sp) != st_crs(spatial_areas)){
    spatial_areas<- st_transform(spatial_areas, st_crs(pred_sp))
  }
  pred_sp<- pred_sp %>%
    st_join(., spatial_areas) %>%
    st_drop_geometry()
  
  res_out <- pred_sp %>%
    mutate(., Date = as.Date(Time),
           "Year" = format(Date, "%Y"),
           "Month" = format(Date, "%m"),
           "Season" = ifelse(Month == "03", "Spring", ifelse(Month == "10", "Fall", "Summer"))
    )
  
  res_out_seas<- res_out %>%
    group_by(., Region, Year, Season) %>%
    summarize_at(., vars(Biomass), list(Biomass = sum), na.rm = FALSE) %>%
    mutate(., "Log_Biomass" = log(Biomass + 1)) %>%
    arrange(., Region, Year, Season)
  
  res_out_all<- res_out %>%
    group_by(., Region, Year) %>%
    summarize_at(., vars(Biomass), list(Biomass = sum), na.rm = FALSE) %>%
    mutate(., "Log_Biomass" = log(Biomass + 1),
           "Season" = "All") %>%
    arrange(., Region, Year, Season) 
  
  res_out_final<- bind_rows(res_out_seas, res_out_all)
  res_out_final$Season<- factor(res_out_final$Season, levels = c("Spring", "Summer", "Fall", "All"), labels = c("Spring", "Summer", "Fall", "All"))
  res_out_final$Species<- nice_species_name
  res_out_final$Climate_Scenario<- nice_climate_name
  
  write.csv(res_out_final, file = paste0(out_dir, paste(nice_species_name, nice_climate_name, "Biomass_Index", sep = "_"), ".csv"))
  return(res_out_final)
  
}

make_spatial_index<- function(vast_proj_rds, vast_fit, spatial_aras, nice_times, nice_species_name, nice_climate_name, out_dir){
  if(FALSE){
    vast_proj_rds <- readRDS(paste0("/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/targets_output/", "prediction_df/American_lobster_ST_SSP5_85_mean_projections.rds"))
    #vast_proj_rds<- read.csv("/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/DFOSDM_app/projected_no_ST/American_lobster_ST_SSP585_mean_projected.csv")
    vast_fit <- readRDS(paste0("/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/targets_output/", "mod_fits/American_lobster_ST_fitted_vast.rds"))
    spatial_areas<- st_read("/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/DFOSDM_app/shp/spatial_areas.shp")
    nice_species_name <- "American_lobster"
    nice_climate_name <- "SSP85_mean"
    out_dir <- "/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/DFOSDM_app/spatial_summary"
  }
  
  # Some general stuff
  # Time series steps
  time_ind <- seq(from = 1, to = vast_proj_rds[[1]]$n_t)
  time_labels <- nice_times
  
  # Categories
  categories_ind <- seq(from = 1, to = vast_proj_rds[[1]]$n_c)
  
  # Index regions
  index_regions_ind <- 1:vast_proj_rds[[1]]$n_l
  index_regions <- vast_fit$settings$strata.limits$STRATA[index_regions_ind]
  
  # Grid index
  grid_ind <- seq(from = 1, to = vast_proj_rds[[1]]$n_g)
  
  
  for (i in seq_along(vast_proj_rds)) {
    temp_array <- array(vast_proj_rds[[i]][which(names(vast_proj_rds[[i]]) == "Index_gctl")][[1]], dim = c(unlist(vast_proj_rds[[i]][c("n_g", "n_c", "n_t", "n_l")], 2)), dimnames = list(grid_ind, categories_ind, time_labels, index_regions))
    
    for(j in 1:dim(temp_array)[4]){
      temp_df<- data.frame("Date" = nice_times, "Biomass" = as.numeric(colSums(temp_array[,,,j])), "Region" = vast_fit$settings$strata.limits$STRATA[index_regions_ind][j])
      
      if(j == 1){
        index_df<- temp_df
      } else {
        index_df<- bind_rows(index_df, temp_df)
      }
    }
    
    index_df$Sim_Scenario <- paste0("Sim_", i)
    
    
    if (i == 1) {
      index_out <- index_df
    } else {
      index_out <- bind_rows(index_out, index_df)
    }
  }
  
  # Add in season
  index_out$Season<- ifelse(format(index_out$Date, "%m") == "03", "Spring", ifelse(format(index_out$Date, "%m") == "10", "Fall", "Summer"))
  
  # Log biomass
  index_out$Log_Biomass<- log(index_out$Biomass)
  
  # Year
  index_out$Year<- format(index_out$Date, "%Y")
  
  write.csv(index_out, file = paste0(out_dir, paste(nice_species_name, nice_climate_name, "Biomass_Index", sep = "_"), ".csv"))
  return(index_out)
}

## Plotting functions
plot_diff_map<- function(plot_dat, plot_variable, plot_season, plot_decade, land_sf, pred_res, xlim = c(-78.5, -54.5), ylim = c(34.5, 48.25), lab_lat = 36, lab_lon = -67.5, out_filename, land_color = "#d9d9d9", ...){
  # For debugging
  if(FALSE){
    plot_dat = base_dat
    plot_variable = "Log_Biomass_Diff"
    plot_season = "All"
    plot_decade = 2050
    and_sf = land
    pred_res = 25000
    xlim = c(-78.5, -54.5)
    ylim = c(34.5, 48.25)
    lab_lat = 36
    lab_lon = -67.5
    out_dir
    land_color = "#d9d9d9"
  }
  
  # Plotting bits
  land_use<- land_sf %>%
    st_transform(., crs = 32619)
  lims_use<- st_as_sf(data.frame("Lon" = xlim, "Lat" = ylim), coords = c("Lon", "Lat"), crs = 4326) %>%
    st_transform(., crs = 32619) %>%
    data.frame(st_coordinates(.))
  lab_use<- st_as_sf(data.frame("Lon" = lab_lon, "Lat" = lab_lat), coords = c("Lon", "Lat"), crs = 4326) %>%
    st_transform(., crs = 32619) %>%
    data.frame(st_coordinates(.))
  
  legend_label<- ifelse(plot_variable == "Log_Biomass_Diff", "Log biomass difference (kg/km2)", ifelse(plot_variable == "Biomass_Diff", "Biomass difference (kg/km2)", ifelse(plot_variable == "Log_Biomass", "Log biomass (kg/km2)", "Biomass (kg/km2)")))
  plot_label<- paste(plot_season, plot_decade, sep = " ")
  # plot_lims<- c(0, max(plot_dat[,plot_variable], na.rm = TRUE))
  
  # Filter to time of interest
  plot_points<- plot_dat %>%
    filter(., Season == plot_season & Decade == plot_decade) %>%
    st_as_sf(., coords = c("Lon", "Lat"), crs = 32619, remove = FALSE) 
  
  bSquare <- function(x, a) {
    a <- sqrt(a)/2
    return(sf::st_buffer(x, dist = a, nQuadSegs=1, endCapStyle = "SQUARE"))
  }
  
  plot_polys<- bSquare(plot_points, pred_res*pred_res)
  
  # Saving bits
  nice_category_names<- unique(plot_points$Species)
  climate_scenario<- unique(plot_points$Climate_Scenario)
  
  if(plot_variable == "Log_Biomass_Diff"){
    plot_out<- ggplot() +
      geom_sf(data = plot_polys, aes(fill = Log_Biomass_Diff, color = Log_Biomass_Diff)) +
      scale_fill_distiller(name = legend_label, type = "div", palette = 'RdBu', na.value = "transparent") + 
      scale_color_distiller(name = legend_label, type = "div", palette = 'RdBu', na.value = "transparent") + 
      annotate("text", x = lab_lon, y = lab_lat, label = plot_label) +
      geom_sf(data = land_sf, fill = land_color, lwd = 0.2, na.rm = TRUE) +
      coord_sf(xlim = xlim, ylim = ylim, expand = FALSE, crs = 4326) +
      theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title = element_blank(), plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "pt"))
  }
  
  if(plot_variable == "Biomass_Diff") {
    plot_out<- ggplot() +
      geom_sf(data = plot_polys, aes(fill = Biomass_Diff, color = Biomass_Diff)) +
      scale_fill_distiller(name = legend_label, type = "div", palette = 'RdBu', na.value = "transparent") + 
      scale_color_distiller(name = legend_label, type = "div", palette = 'RdBu', na.value = "transparent") + 
      annotate("text", x = lab_lon, y = lab_lat, label = plot_label) +
      geom_sf(data = land_sf, fill = land_color, lwd = 0.2, na.rm = TRUE) +
      coord_sf(xlim = xlim, ylim = ylim, expand = FALSE, crs = 4326) +
      theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title = element_blank(), plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "pt"))
  }
  
  # Save it and return it
  ggsave(filename = out_filename, plot_out, width = 8, height = 11, units = "in")
  return(plot_out)
}

plot_base_map<- function(plot_dat, plot_variable, plot_season, land_sf, pred_res = 25000, xlim = c(-78.5, -54.5), ylim = c(34.5, 48.25), lab_lat = 36, lab_lon = -67.5, out_filename, land_color = "#d9d9d9", ...){
  # For debugging
  if(FALSE){
    plot_dat_df<-  base_utm
    plot_variable<- "Log_Biomass"
    plot_season = "Spring"
    land_sf<- land
    pred_res<- 25000
    xlim = c(-78.5, -54.5)
    ylim = c(34.5, 48.25)
    lab_lat = 36
    lab_lon = -67.5
    out_dir = "/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/targets_output/shiny/"
    land_color = "#d9d9d9"
  }
  
  # Plotting bits
  land_use<- land_sf %>%
    st_transform(., crs = 32619)
  lims_use<- st_as_sf(data.frame("Lon" = xlim, "Lat" = ylim), coords = c("Lon", "Lat"), crs = 4326) %>%
    st_transform(., crs = 32619) %>%
    data.frame(st_coordinates(.))
  lab_use<- st_as_sf(data.frame("Lon" = lab_lon, "Lat" = lab_lat), coords = c("Lon", "Lat"), crs = 4326) %>%
    st_transform(., crs = 32619) %>%
    data.frame(st_coordinates(.))
  
  legend_label<- ifelse(plot_variable == "Log_Biomass_Diff", "Log biomass difference (kg/km2)", ifelse(plot_variable == "Biomass_Diff", "Biomass difference (kg/km2)", ifelse(plot_variable == "Log_Biomass", "Log biomass (kg/km2)", "Biomass (kg/km2)")))
  plot_label<- paste(plot_season, "baseline", sep = " ")
  #plot_lims<- c(0, max(plot_dat_df[,plot_variable], na.rm = TRUE))
  
  # Filter to time of interest
  plot_points<- plot_dat %>%
    filter(., Season == plot_season) %>%
    st_as_sf(., coords = c("Lon", "Lat"), crs = 32619, remove = FALSE) 
  
  bSquare <- function(x, a) {
    a <- sqrt(a)/2
    return(sf::st_buffer(x, dist = a, nQuadSegs=1, endCapStyle = "SQUARE"))
  }
  
  plot_polys<- bSquare(plot_points, pred_res*pred_res)
  
  # Saving bits
  nice_category_names<- unique(plot_points$Species)
  climate_scenario<- unique(plot_points$Climate_Scenario)
  
  if(plot_variable == "Log_Biomass"){
    plot_out<- ggplot() +
      geom_sf(data = plot_polys, aes(fill = Log_Biomass, color = Log_Biomass)) +
      scale_fill_viridis_c(name = legend_label, option = "viridis", na.value = "transparent") +
      scale_color_viridis_c(name = legend_label, option = "viridis", na.value = "transparent") +
      annotate("text", x = lab_lon, y = lab_lat, label = plot_label) +
      geom_sf(data = land_use, fill = land_color, lwd = 0.2, na.rm = TRUE) +
      coord_sf(xlim = xlim, ylim = ylim, crs = 4326) +
      theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title = element_blank(), plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "pt"))
  }
  
  if(plot_variable == "Biomass"){
    plot_out<- ggplot() +
      geom_sf(data = plot_polys, aes(fill = Biomass, color = Biomass)) +
      scale_fill_viridis_c(name = legend_label, option = "viridis", na.value = "transparent") +
      scale_color_viridis_c(name = legend_label, option = "viridis", na.value = "transparent") +
      annotate("text", x = lab_lon, y = lab_lat, label = plot_label) +
      geom_sf(data = land_sf, fill = land_color, lwd = 0.2, na.rm = TRUE) +
      coord_sf(xlim = xlim, ylim = ylim, expand = FALSE, crs = 4326) +
      theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title = element_blank(), plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "pt"))
  }
  
  # Save it and return it
  ggsave(filename = out_filename, plot_out, width = 8, height = 11, units = "in")
  return(plot_out)
}

plot_proj_map<- function(plot_dat, plot_variable, plot_season, plot_decade, land_sf, pred_res = 25000, xlim = c(-78.5, -54.5), ylim = c(34.5, 48.25), lab_lat = 36, lab_lon = -67.5, out_filename, land_color = "#d9d9d9", ...){
  # For debugging
  if(FALSE){
    plot_dat_df<-  base_dat
    plot_variable<- "Log_Biomass"
    plot_season = "All"
    land_sf<- land
    xlim = c(-78.5, -54.5)
    ylim = c(34.5, 48.25)
    lab_lat = 36
    lab_lon = -67.5
    out_dir = "/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/targets_output/shiny/"
    land_color = "#d9d9d9"
  }
  
  # Plotting bits
  land_use<- land_sf %>%
    st_transform(., crs = 32619)
  lims_use<- st_as_sf(data.frame("Lon" = xlim, "Lat" = ylim), coords = c("Lon", "Lat"), crs = 4326) %>%
    st_transform(., crs = 32619) %>%
    data.frame(st_coordinates(.))
  lab_use<- st_as_sf(data.frame("Lon" = lab_lon, "Lat" = lab_lat), coords = c("Lon", "Lat"), crs = 4326) %>%
    st_transform(., crs = 32619) %>%
    data.frame(st_coordinates(.))
  
  legend_label<- ifelse(plot_variable == "Log_Biomass_Diff", "Log biomass difference (kg/km2)", ifelse(plot_variable == "Biomass_Diff", "Biomass difference (kg/km2)", ifelse(plot_variable == "Log_Biomass", "Log biomass (kg/km2)", "Biomass (kg/km2)")))
  plot_label<- paste(plot_season, plot_decade, sep = " ")
  # plot_lims<- c(0, max(plot_dat_df[,plot_variable], na.rm = TRUE))
  
  # Filter to time of interest
  plot_points<- plot_dat %>%
    filter(., Season == plot_season & Decade == plot_decade) %>%
    st_as_sf(., coords = c("Lon", "Lat"), crs = 32619, remove = FALSE) 
  
  bSquare <- function(x, a) {
    a <- sqrt(a)/2
    return(sf::st_buffer(x, dist = a, nQuadSegs=1, endCapStyle = "SQUARE"))
  }
  
  plot_polys<- bSquare(plot_points, pred_res*pred_res)
  
  # Saving bits
  nice_category_names<- unique(plot_points$Species)
  climate_scenario<- unique(plot_points$Climate_Scenario)
  
  if(plot_variable == "Log_Biomass"){
    plot_out<- ggplot() +
      geom_sf(data = plot_polys, aes(fill = Log_Biomass, color = Log_Biomass)) +
      scale_fill_viridis_c(name = legend_label, option = "viridis", na.value = "transparent") +
      scale_color_viridis_c(name = legend_label, option = "viridis", na.value = "transparent") +
      annotate("text", x = lab_lon, y = lab_lat, label = plot_label) +
      geom_sf(data = land_sf, fill = land_color, lwd = 0.2, na.rm = TRUE) +
      coord_sf(xlim = xlim, ylim = ylim, expand = FALSE, crs = 4326) +
      theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title = element_blank(), plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "pt"))
  }
  
  if(plot_variable == "Biomass"){
    plot_out<- ggplot() +
      geom_sf(data = plot_polys, aes(fill = Biomass, color = Biomass)) +
      scale_fill_viridis_c(name = legend_label, option = "viridis", na.value = "transparent") +
      scale_color_viridis_c(name = legend_label, option = "viridis", na.value = "transparent") +
      annotate("text", x = lab_lon, y = lab_lat, label = plot_label) +
      geom_sf(data = land_sf, fill = land_color, lwd = 0.2, na.rm = TRUE) +
      coord_sf(xlim = xlim, ylim = ylim, expand = FALSE, crs = 4326) +
      theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title = element_blank(), plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "pt"))
  }
  
  # Save it and return it
  ggsave(filename = out_filename, plot_out, width = 8, height = 11, units = "in")
  return(plot_out)
}

plot_cog<- function(cog_dat, plot_season, land_sf, xlim = c(-78.5, -54.5), ylim = c(34.5, 48.25), out_filename, land_color = "#d9d9d9"){
  
  if(FALSE){
    cog_dat = cog_dat
    plot_season = "All"
    land_sf = land
    xlim = c(-78.5, -54.5)
    ylim = c(34.5, 48.25)
    out_filename = "/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/DFOSDM_app/center_biomass/cog_plot.jpg"
    land_color = "#d9d9d9"
  }
  
  cog_sf<- cog_dat %>%
    filter(., Season == plot_season) %>%
    st_as_sf(., coords = c("Lon", "Lat"), remove = FALSE, crs = 4326)
  
  cog_plot <- ggplot() +
    geom_sf(data = cog_sf, aes(fill = Year), size = 2, shape = 21) +
    scale_fill_viridis_c(name = "Year", limits = c(min(cog_sf$Year), max(cog_sf$Year))) +
    geom_sf(data = land_sf, fill = land_color, lwd = 0.2, na.rm = TRUE) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "pt"))
  
  # Now, the lon/lat time series
  lon_lat_df <- cog_sf %>%
    data.frame(st_coordinates(.))
  lon_lat_df$Date <- as.Date(paste0(lon_lat_df$Year, "-06-15"))
  
  lon_ts <- ggplot() +
    geom_line(data = lon_lat_df, aes(x = Date, y = Lon), color = "#66c2a5", lwd = 2) +
    geom_point(data = lon_lat_df, aes(x = Date, y = Lon), pch = 21, fill = "#66c2a5") +
    scale_x_date(date_breaks = "5 year", date_labels = "%Y", expand = c(0, 0), limits = as.Date(c("1985-01-01", "2100-01-01"))) +
    ylab("Center of longitude") +
    xlab("Date") +
    theme_bw() +
    theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  lat_ts <- ggplot() +
    geom_line(data = lon_lat_df, aes(x = Date, y = Lat), color = "#66c2a5", lwd = 2) +
    geom_point(data = lon_lat_df, aes(x = Date, y = Lat), pch = 21, fill = "#66c2a5") +
    scale_x_date(date_breaks = "5 year", date_labels = "%Y", expand = c(0, 0), limits = as.Date(c("1985-01-01", "2100-01-01"))) +
    ylab("Center of latitude") +
    xlab("Date") +
    theme_bw() +
    theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  plot_out <- (cog_plot) / (lon_ts + lat_ts) + plot_layout(ncol = 1, nrow = 2, widths = c(0.75, 1), heights = c(0.75, 1))
  
  # Save and return it
  ggsave(filename = out_filename, plot = plot_out, width = 8, height = 11, units = "in")
  return(plot_out)
}

plot_biomass_index<- function(biomass_dat, plot_season, plot_variable, plot_regions = c("DFO", "NMFS", "GoM", "SNE_and_MAB"), nice_xlab = "Year", nice_ylab = "Total log biomass (kg)", color_pal = c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00'), out_filename){
  if(FALSE){
    biomass_dat<- index_out
    plot_variable<- "Log_Biomass"
    plot_regions = c("DFO", "NMFS", "GoM", "SNE_and_MAB")
    nice_xlab = "Year"
    nice_ylab = "Total log biomass (kg)"
    color_pal = c('#e41a1c','#377eb8','#4daf4a','#984ea3')
    out_filename 
  }
  
  # Color selection
  if(is.null(plot_regions)){
    plot_regions<- unique(biomass_dat$Region)
  }
  if(!is.null(color_pal)){
    colors_use<- color_pal
  } else {
    color_pal<- c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02')
    colors_use<- color_pal[1:length(plot_regions)]
  }
  
  # Date stuffs
  biomass_dat$Date<- as.Date(paste0(biomass_dat$Year, "-01-01"))
  biomass_dat<- biomass_dat %>%
    filter(., Season == plot_season & Date >= "1985-01-01" & Date <= "2100-01-01" & Region %in% plot_regions)
  date_breaks<- seq.Date(from = as.Date(min(biomass_dat$Date)), to = as.Date(max(biomass_dat$Date)), by = "5 years")

  plot_out<- ggplot() +
    geom_point(data = biomass_dat, aes_string(x = "Date", y = plot_variable, color = "Region")) +
    geom_line(data = biomass_dat, aes_string(x = "Date", y = plot_variable, color = "Region")) +
    scale_color_manual(values = colors_use) +
    scale_x_date(breaks = date_breaks, date_labels = "%Y") +
    xlab({{nice_xlab}}) +
    ylab({{nice_ylab}}) +
    ggtitle(paste0(unique(biomass_dat$Species), " ", unique(biomass_dat$Climate_Scenario))) + 
    theme_bw() +
    theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
  
  # Save and return the plot
  ggsave(filename = out_filename, width = 8, height = 11, units = "in")
  return(plot_out)
}