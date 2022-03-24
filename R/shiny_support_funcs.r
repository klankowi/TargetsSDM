######
## Collection of functions to support shiny app presentation
######

dat <- readRDS(paste0("/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/targets_output/", "prediction_df/American_lobster_GD2_BetaRW_EpsIIDRED_Temp_SSP5_85_mean_projections.rds"))
summary(dat)

# Some general stuff
# Time series steps
time_ind <- seq(from = 1, to = sim_obj[[1]]$n_t)
time_labels <- nice_times

# Index regions
index_regions_ind <- seq(from = 1, to = sim_obj[[1]]$n_l)
index_regions <- vast_fit$settings$strata.limits$STRATA[index_regions_ind]

# Categories
categories_ind <- seq(from = 1, to = sim_obj[[1]]$n_c)

# Grid locations
grid_ind <- seq(from = 1, to = sim_obj[[1]]$n_g)

# Loop through sims and get results
for (i in seq_along(sim_obj)) {

    # Index values?
    if (what == "Index_ctl") {
        temp_array <- array(sim_obj[[i]][which(names(sim_obj[[i]]) == "Index_ctl")][[1]], dim = c(unlist(sim_obj[[i]][c("n_c", "n_t", "n_l")])), dimnames = list(categories_ind, time_labels, index_regions))

        temp_df <- data.frame(aperm(temp_array, c(2, 3, 1)))
        colnames(temp_df) <- gsub(".1", "", colnames(temp_df))

        temp_df$Sim_Scenario <- rep(paste0("Sim_", i), nrow(temp_df))

        temp_df$Time <- nice_times

        if (i == 1) {
            res_out <- temp_df
        } else {
            res_out <- bind_rows(res_out, temp_df)
        }
    }

    # Predicted grid density?
    if (what == "D_gct") {
        temp_array <- array(sim_obj[[i]][which(names(sim_obj[[i]]) == "D_gct")][[1]], dim = c(unlist(sim_obj[[i]][c("n_g", "n_c", "n_t")])), dimnames = list(grid_ind, categories_ind, time_labels))

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
            res_out <- temp_df
        } else {
            res_out <- bind_rows(res_out, temp_df)
        }
    }
}

# Wide to long for "Index_ctl"
if (what == "Index_ctl") {
    # Annual average?
    if (!is.null(out_t_scale)) {
        res_out <- res_out %>%
            pivot_longer(., !c(Sim_Scenario, Time), names_to = "Region", values_to = "Index") %>%
            mutate(., Time = ymd(as.numeric(format(Time, "%Y")), truncated = 2L)) %>%
            group_by(., Time, Region) %>%
            summarise(
                Prob_0.5 = quantile(Index, probs = 0.50),
                Prob_0.1 = quantile(Index, probs = 0.1),
                Prob_0.9 = quantile(Index, probs = 0.9)
            )
    } else {
        res_out <- res_out %>%
            pivot_longer(., !c(Sim_Scenario, Time), names_to = "Region", values_to = "Index") %>%
            group_by(., Time, Region) %>%
            summarise(
                Prob_0.5 = quantile(Index, probs = 0.50),
                Prob_0.1 = quantile(Index, probs = 0.1),
                Prob_0.9 = quantile(Index, probs = 0.9)
            )
    }
}

if (what == "D_gct") {
    # Annual average?
    if (!is.null(out_t_scale)) {
        res_out <- res_out %>%
            mutate(., Time = ymd(as.numeric(format(Time, "%Y")), truncated = 2L)) %>%
            group_by(., Time, Region) %>%
            summarise(
                Prob_0.5 = quantile(Index, probs = 0.50)
            )
    } else {
        res_out <- res_out %>%
            group_by(., Lat, Lon, Time) %>%
            summarise(
                Prob_0.5 = quantile(D_gct, probs = 0.50),
                Prob_0.1 = quantile(D_gct, probs = 0.1),
                Prob_0.9 = quantile(D_gct, probs = 0.9)
            )
    }
}

## Make baseline data function -- for this, we want to have columns for x, y, species, season (spring, summer, fall, all), climate scenario, biomass and log biomass
nice_times <- as.Date(c(paste0(seq(from = 1985, to = 2100, by = 1), "-03-16"), paste0(seq(from = 1985, to = 2100, by = 1), "-07-16"), paste0(seq(from = 1985, to = 2100, by = 1), "-10-16")))
nice_times <- nice_times[order(nice_times)]

make_baseline_shiny <- function(vast_proj_rds, nice_times, base_years, nice_species_name, nice_climate_name, out_crs, out_dir) {
    if (FALSE) {
        vast_proj_rds <- readRDS(paste0("/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/targets_output/", "prediction_df/American_lobster_GD2_BetaRW_EpsIIDRED_Temp_SSP5_85_mean_projections.rds"))
        vast_fit <- readRDS(paste0("/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/targets_output/", "mod_fits/American_lobster_GD2_BetaRW_EpsIIDRED_Temp_fitted_vast.rds"))
        base_years <- c(2015, 2019)
        nice_species_name <- "American_lobster"
        nice_climate_name <- "SSP85_mean"
        out_dir <- "/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/DFOSDM_app/baseline_data"
        out_crs <- 4326
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
}
## Make center of biomass data function -- for this, we want to have columns for x, y, season, year

## Make projection data -- for this, we want to have x, y, decade date, species, climate scenario, biomass and log biomass

## Make spatial summaries -- for this, we want to have species, season, year, scenario, mean, sd, date, climate scenario