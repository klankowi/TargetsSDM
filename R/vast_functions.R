####  Common Resources  ####

pred_template_load <- function(pred_template_dir) {
  if (FALSE) {
    tar_load(pred_template_dir)
  }
  
  # Load the raster template gird
  pred_template_rast <- raster(paste(pred_template_dir, "mod_pred_template.grd", sep = "/"))

  # Convert it to a data frame
  pred_template_df <- as.data.frame(pred_template_rast, xy = TRUE) %>%
    drop_na() %>%
    dplyr::select(., x, y)
  names(pred_template_df) <- c("longitude", "latitude")

  # Return it
  return(pred_template_df)
}


high_res_load <- function(high_res_dir) {
  high_res <- raster(paste(high_res_dir, "HighResTemplate.grd", sep = "/"))
  return(high_res)
}

####  Functions  ####
####

#' @title Make VAST prediction dataframe
#'
#' @description This function creates a dataframe of prediction covariates to combine with the other VAST data
#'
#' @param predict_covariates_stack_agg = The directory holding processed covariate raster stacks
#' @param mask = Shapefile mask
#' @param summarize = Currently, either "annual" or "seasonal" to indicate whether the each dynamic raster stack should be summarized to an annual or seasonal time scale
#' @param ensemble_stat = Either the climate model ensemble statistic to use when working with climate model projections, or NULL. This is only used in naming the output file
#' @param fit_year_min
#' @param fit_year_max
#' @param pred_years
#' @param out_dir = Directory to save the prediction dataframe
#'
#' @return A dataframe with prediction information. This file is also saved in out_dir.
#'
#' @export

make_vast_predict_df <- function(predict_covariates_stack_agg, extra_covariates_stack, covs_rescale = c("Depth", "BS_seasonal", "BT_seasonal", "SS_seasonal", "SST_seasonal"), rescale_params, depth_cut, mask, summarize, ensemble_stat, fit_seasons, fit_year_min, fit_year_max, test_year_max, out_dir) {

  # For debugging
  if (FALSE) {
    tar_load(predict_covariates_stack_agg_out)
    predict_covariates_stack_agg <- predict_covariates_stack_agg_out
    tar_load(static_covariates_stack)
    extra_covariates_stack <- static_covariates_stack
    tar_load(rescale_params)
    tar_load(region_shapefile)
    mask <- region_shapefile
    summarize <- "seasonal"
    ensemble_stat <- "mean"
    fit_year_min <- fit_year_min
    fit_year_max <- fit_year_max
    pred_years <- pred_years
    out_dir <- here::here("scratch/aja/TargetsSDM/data/predict")
    covs_rescale <- c("Depth", "BS_seasonal", "BT_seasonal", "SS_seasonal", "SST_seasonal")
  }

  ####
  ## Need to figure out what to do about depth here!!!

  # Get raster stack covariate files
  rast_files_load <- list.files(predict_covariates_stack_agg, pattern = paste0(summarize, "_", ensemble_stat, ".grd$"), full.names = TRUE)

  # Get variable names
  cov_names_full <- list.files(predict_covariates_stack_agg, pattern = paste0(summarize, "_", ensemble_stat, ".grd$"), full.names = FALSE)
  predict_covs_names <- gsub(paste("_", ensemble_stat, ".grd$", sep = ""), "", gsub("predict_stack_", "", cov_names_full))

  # Looping through prediction stack time steps
  for (i in 1:nlayers(raster::stack(rast_files_load[1]))) {
    # Get the time index
    time_ind <- i

    if ((length(covs_rescale) - 1) == 4) {
      # Load corresponding raster layers matching the time index
      pred_covs_stack_temp <- rotate(raster::stack(raster::stack(rast_files_load[1])[[time_ind]], raster::stack(rast_files_load[2])[[time_ind]], raster::stack(rast_files_load[3])[[time_ind]], raster::stack(rast_files_load[4])[[time_ind]]))
    } else if ((length(covs_rescale) - 1) == 3) {
      # Load corresponding raster layers matching the time index
      pred_covs_stack_temp <- rotate(raster::stack(raster::stack(rast_files_load[1])[[time_ind]], raster::stack(rast_files_load[2])[[time_ind]], raster::stack(rast_files_load[3])[[time_ind]]))
    } else if ((length(covs_rescale) - 1) == 2) {
      # Load corresponding raster layers matching the time index
      pred_covs_stack_temp <- rotate(raster::stack(raster::stack(rast_files_load[1])[[time_ind]], raster::stack(rast_files_load[2])[[time_ind]]))
    } else if ((length(covs_rescale) - 1) == 1) {
      # Load corresponding raster layers matching the time index
      pred_covs_stack_temp <- rotate(raster::stack(raster::stack(rast_files_load[1])[[time_ind]]))
    }
  
    # Mask out values outside area of interest
    pred_covs_stack_temp <- raster::mask(pred_covs_stack_temp, mask = mask)

    # Some processing to keep observations within our area of interest and get things in a "tidy-er" prediction dataframe
    time_name <- sub(".[^.]*$", "", names(pred_covs_stack_temp))
    names(pred_covs_stack_temp) <- paste(time_name, predict_covs_names, sep = "_")
    pred_covs_df_temp <- as.data.frame(pred_covs_stack_temp, xy = TRUE) %>%
      drop_na()
    colnames(pred_covs_df_temp)[2:ncol(pred_covs_df_temp)] <- gsub("X", "", gsub("[.]", "_", colnames(pred_covs_df_temp)[2:ncol(pred_covs_df_temp)]))
    colnames(pred_covs_df_temp)[1:2] <- c("DECDEG_BEGLON", "DECDEG_BEGLAT")

    pred_covs_df_out_temp <- pred_covs_df_temp %>%
      pivot_longer(., -c(DECDEG_BEGLON, DECDEG_BEGLAT), names_to = c("variable"), values_to = "value") %>%
      separate(., variable, into = c("EST_YEAR", "SEASON", "variable"), sep = "_", extra = "merge") %>%
      pivot_wider(., names_from = variable, values_from = value)

    # Adding in some other columns we will want to match up easily with 'vast_data_out'
    pred_covs_df_out_temp <- pred_covs_df_out_temp %>%
      mutate(.,
        EST_YEAR = as.numeric(EST_YEAR),
        DATE = paste(EST_YEAR, case_when(
          SEASON == "Winter" ~ "12-16",
          SEASON == "Spring" ~ "03-16",
          SEASON == "Summer" ~ "07-16",
          SEASON == "Fall" ~ "09-16"
        ), sep = "-"),
        SURVEY = "NMFS",
        SVVESSEL = "DUMMY",
        NMFS_SVSPP = "DUMMY",
        DFO_SPEC = "DUMMY",
        PRESENCE = 1,
        BIOMASS = 1,
        ABUNDANCE = 1,
        ID = paste("DUMMY", DATE, sep = ""),
        PredTF = TRUE
      )

    if (i == 1) {
      pred_covs_out <- pred_covs_df_out_temp
    } else {
      pred_covs_out <- bind_rows(pred_covs_out, pred_covs_df_out_temp)
    }
  }

  # Only going to keep information from fit_year_max through pred_years...
  pred_covs_out_final <- pred_covs_out %>%
    dplyr::filter(., EST_YEAR > fit_year_max & EST_YEAR <= test_year_max)

  # New implementation...
  pred_covs_out_final <- pred_covs_out_final %>%
    mutate(., # VAST_YEAR_COV = EST_YEAR,
      # VAST_YEAR_COV = ifelse(EST_YEAR > fit_year_max, fit_year_max, EST_YEAR),
      VAST_YEAR_COV = ifelse(EST_YEAR > fit_year_max, test_year_max, EST_YEAR),
      VAST_SEASON = case_when(
        SEASON == "Spring" ~ "SPRING",
        SEASON == "Summer" ~ "SUMMER",
        SEASON == "Fall" ~ "FALL"
      ),
      "VAST_YEAR_SEASON" = paste(EST_YEAR, VAST_SEASON, sep = "_")
    )

  # Subset to only seasons of interest...
  pred_covs_out_final <- pred_covs_out_final %>%
    filter(., VAST_SEASON %in% fit_seasons)

  # Need to account for new levels in year season...
  all_years <- seq(from = fit_year_min, to = test_year_max, by = 1)
  all_seasons <- fit_seasons
  year_season_set <- expand.grid("SEASON" = all_seasons, "EST_YEAR" = all_years)
  all_year_season_levels <- apply(year_season_set[, 2:1], MARGIN = 1, FUN = paste, collapse = "_")

  pred_covs_out_final <- pred_covs_out_final %>%
    mutate(.,
      "VAST_YEAR_SEASON" = factor(VAST_YEAR_SEASON, levels = all_year_season_levels),
      "VAST_SEASON" = factor(VAST_SEASON, levels = all_seasons)
    )

  # Name rearrangement!
  # Keep only what we need..
  cov_names <- names(pred_covs_out_final)[-which(names(pred_covs_out_final) %in% c("ID", "DATE", "EST_YEAR", "SEASON", "SURVEY", "SVVESSEL", "DECDEG_BEGLAT", "DECDEG_BEGLON", "NMFS_SVSPP", "DFO_SPEC", "PRESENCE", "BIOMASS", "ABUNDANCE", "PredTF", "VAST_YEAR_COV", "VAST_SEASON", "VAST_YEAR_SEASON"))]
  pred_covs_out_final <- pred_covs_out_final %>%
    dplyr::select(., "ID", "DATE", "EST_YEAR", "SEASON", "SURVEY", "SVVESSEL", "DECDEG_BEGLAT", "DECDEG_BEGLON", "NMFS_SVSPP", "DFO_SPEC", "PRESENCE", "BIOMASS", "ABUNDANCE", "PredTF", "VAST_YEAR_COV", "VAST_SEASON", "VAST_YEAR_SEASON", {{ cov_names }})

  # Any extra covariates will likely be static...
  if (!is.null(extra_covariates_stack)) {
    pred_covs_sf <- points_to_sf(pred_covs_out_final)

    pred_covs_out_final <- static_extract_wrapper(static_covariates_list = extra_covariates_stack, sf_points = pred_covs_sf, date_col_name = "DATE", df_sf = "df", out_dir = NULL)
  }

  # Apply depth cut and drop NAs
  pred_covs_out_final <- pred_covs_out_final %>%
    mutate(.,
      "Depth" = ifelse(Depth > depth_cut, NA, Depth),
      "Summarized" = summarize,
      "Ensemble_Stat" = ensemble_stat
    ) %>%
    drop_na()

  # Rescale
  if (!is.null(rescale_params)) {
    for (i in seq_along(covs_rescale)) {
      match_mean <- rescale_params[which(names(rescale_params) == paste(covs_rescale[i], "Mean", sep = "_"))]
      match_sd <- rescale_params[which(names(rescale_params) == paste(covs_rescale[i], "SD", sep = "_"))]
      pred_covs_out_final <- pred_covs_out_final %>%
        mutate_at(., {{ covs_rescale[i] }}, .funs = covariate_rescale_func, type = "AJA", center = match_mean, scale = match_sd)
    }
  }

  saveRDS(pred_covs_out_final, file = paste(out_dir, "/VAST_pred_df_", summarize, "_", ensemble_stat, ".rds", sep = ""))
  return(pred_covs_out_final)
}

#' @title Make VAST seasonal dataset
#'
#' @description This function reads in a tidy model dataset and does some cleaning and processing to generate a new dataset to accommodate fitting a VAST seasonal (or other intra annual) model. These cleaning and processing steps boil down to creating an ordered, continuous, season-year vector, such that the model can then estimate density even in season-years not surveyed.
#'
#' @param tidy_mod_data = A tidy model datafame with all the information (tows, habitat covariates, species occurrences) needed to fit a species distribution model.
#' @param nmfs_species_code = Numeric NMFS species code
#' @param fit_year_min = Minimum year to keep
#' @param fit_year_max = Maximum year to keep
#' @param pred_df = Either NULL or a dataframe with prediction information as created by `make_vast_predict_df`
#' @param out_dir = Directory to save the tidy model dataframe as an .rds file
#'
#' @return  A VAST seasonal dataset, ready to be split into a `sample data` dataframe and a `covariate data` dataframe. This file is also saved in out_dir.
#'
#' @export

make_vast_seasonal_data <- function(tidy_mod_data, fit_seasons, nmfs_species_code, fit_year_min, fit_year_max, test_year_max, pred_df, out_dir) {

  # For debugging
  if (FALSE) {
    tar_load(tidy_mod_data)
    nmfs_species_code <- nmfs_species_code
    fit_year_min <- fit_year_min
    fit_year_max <- fit_year_max
    fit_seasons <- fit_seasons
    test_year_max<- test_year_max
    tar_load(vast_predict_df)
    pred_df <- NULL
    out_dir <- here::here("scratch/aja/targets_flow/data/combined/")
    fit_seasons
  }

  # Some work on the time span and seasons
  # Previous implementation before trying to include both surveys within a given season
  # data_temp<- tidy_mod_data %>%
  #   filter(., NMFS_SVSPP == nmfs_species_code) %>%
  #   filter(., EST_YEAR >= fit_year_min & EST_YEAR <= fit_year_max) %>%
  #   mutate(., "VAST_SEASON" = case_when(
  #     SURVEY == "DFO" & SEASON == "SPRING" ~ "DFO",
  #     SURVEY == "NMFS" & SEASON == "SPRING" ~ "SPRING",
  #     SURVEY == "DFO" & SEASON == "SUMMER" ~ "SUMMER",
  #     SURVEY == "NMFS" & SEASON == "FALL" ~ "FALL")) %>%
  #   drop_na(VAST_SEASON)

  # New implementatiom...
  data_temp <- tidy_mod_data %>%
    filter(., NMFS_SVSPP == nmfs_species_code) %>%
    # filter(., EST_YEAR >= fit_year_min & EST_YEAR <= fit_year_max) %>%
    filter(., EST_YEAR >= fit_year_min & EST_YEAR <= test_year_max) %>%
    mutate(., "VAST_SEASON" = case_when(
      SURVEY == "DFO" & SEASON == "SPRING" ~ "SPRING",
      SURVEY == "NMFS" & SEASON == "SPRING" ~ "SPRING",
      SURVEY == "DFO" & SEASON == "SUMMER" ~ "SUMMER",
      SURVEY == "NMFS" & SEASON == "FALL" ~ "FALL",
      SURVEY == "DFO" & SEASON == "FALL" ~ as.character("NA")
    )) %>%
    drop_na(VAST_SEASON)

  data_temp <- data_temp %>%
    filter(., VAST_SEASON %in% fit_seasons)

  # Set of years and seasons. The DFO spring survey usually occurs before the NOAA NEFSC spring survey, so ordering accordingly. Pred year max or fit year max??
  all_years<- seq(from = fit_year_min, to = test_year_max, by = 1)
  #all_years <- seq(from = fit_year_min, to = test_year_max, by = 1)
  all_seasons <- fit_seasons
  yearseason_set <- expand.grid("SEASON" = all_seasons, "EST_YEAR" = all_years)
  all_yearseason_levels <- apply(yearseason_set[, 2:1], MARGIN = 1, FUN = paste, collapse = "_")

  # year_set<- sort(unique(data_temp$EST_YEAR))
  # season_set<- c("DFO", "SPRING", "FALL")
  #
  # # Create a grid with all unique combinations of seasons and years and then combine these into one "year_season" variable
  # yearseason_grid<- expand.grid("SEASON" = season_set, "EST_YEAR" = year_set)
  # yearseason_levels<- apply(yearseason_grid[, 2:1], MARGIN = 1, FUN = paste, collapse = "_")
  # yearseason_labels<- round(yearseason_grid$EST_YEAR + (as.numeric(factor(yearseason_grid$VAST_SEASON, levels = season_set))-1)/length(season_set), digits = 1)
  #
  # Similar process, but for the observations
  yearseason_i <- apply(data_temp[, c("EST_YEAR", "VAST_SEASON")], MARGIN = 1, FUN = paste, collapse = "_")
  yearseason_i <- factor(yearseason_i, levels = all_yearseason_levels)

  # Add the year_season factor column to our sampling_data data set
  data_temp$VAST_YEAR_SEASON <- yearseason_i
  data_temp$VAST_SEASON <- factor(data_temp$VAST_SEASON, levels = all_seasons)

  # VAST year
  data_temp$PredTF <- ifelse(data_temp$EST_YEAR <= fit_year_max, FALSE, TRUE)
  # data_temp$VAST_YEAR_COV <- ifelse(data_temp$EST_YEAR > fit_year_max, fit_year_max, data_temp$EST_YEAR)
  # data_temp$VAST_YEAR_COV<- ifelse(data_temp$EST_YEAR > fit_year_max, fit_year_max, data_temp$EST_YEAR)
  data_temp$VAST_YEAR_COV<- data_temp$EST_YEAR

  # Ordering...
  cov_names <- names(data_temp)[-which(names(data_temp) %in% c("ID", "DATE", "EST_YEAR", "SEASON", "SURVEY", "SVVESSEL", "DECDEG_BEGLAT", "DECDEG_BEGLON", "NMFS_SVSPP", "DFO_SPEC", "PRESENCE", "BIOMASS", "ABUNDANCE", "PredTF", "VAST_YEAR_COV", "VAST_SEASON", "VAST_YEAR_SEASON"))]
  # cov_names <- cov_names[-which(cov_names == "Season_Match")]
  data_temp <- data_temp %>%
    dplyr::select("ID", "DATE", "EST_YEAR", "SEASON", "SURVEY", "SVVESSEL", "DECDEG_BEGLAT", "DECDEG_BEGLON", "NMFS_SVSPP", "DFO_SPEC", "PRESENCE", "BIOMASS", "ABUNDANCE", "PredTF", "VAST_YEAR_COV", "VAST_SEASON", "VAST_YEAR_SEASON", {{ cov_names }})

  # Make dummy data for all year_seasons to estimate gaps in sampling if needed
  dummy_data <- data.frame("ID" = sample(data_temp$ID, size = 1), "DATE" = mean(data_temp$DATE, na.rm = TRUE), "EST_YEAR" = yearseason_set[, "EST_YEAR"], "SEASON" = yearseason_set[, "SEASON"], "SURVEY" = "NMFS", "SVVESSEL" = "DUMMY", "DECDEG_BEGLAT" = mean(data_temp$DECDEG_BEGLAT, na.rm = TRUE), "DECDEG_BEGLON" = mean(data_temp$DECDEG_BEGLON, na.rm = TRUE), "NMFS_SVSPP" = "NMFS", "DFO_SPEC" = "DUMMY", "PRESENCE" = 0, "BIOMASS" = 0, "ABUNDANCE" = 0, "PredTF" = TRUE, "VAST_YEAR_COV" = yearseason_set[, "EST_YEAR"], "VAST_SEASON" = yearseason_set[, "SEASON"], "VAST_YEAR_SEASON" = all_yearseason_levels)

  # Add in "covariates"
  col_ind <- ncol(dummy_data)
  for (i in seq_along(cov_names)) {
    col_ind <- col_ind + 1
    cov_vec <- unlist(data_temp[, {{ cov_names }}[i]])
    dummy_data[, col_ind] <- mean(cov_vec, na.rm = TRUE)
    names(dummy_data)[col_ind] <- {{ cov_names }}[i]
  }

  # Combine with original dataset
  vast_data_out <- rbind(data_temp, dummy_data)
  # vast_data_out$VAST_YEAR_COV<- factor(vast_data_out$VAST_YEAR_COV, levels = seq(from = fit_year_min, to = fit_year_max, by = 1))
  vast_data_out$VAST_YEAR_COV <- factor(vast_data_out$VAST_YEAR_COV, levels = seq(from = fit_year_min, to = test_year_max, by = 1))

  # If we have additional years that we want to predict to and NOT Fit too, we aren't quite done just yet...
  if (!is.null(pred_df)) {
    # Name work...
    pred_df <- pred_df %>%
      dplyr::select(., -Summarized, -Ensemble_Stat)

    # Add those -- check names first
    check_names <- all(colnames(pred_df) %in% colnames(vast_data_out)) & all(colnames(vast_data_out) %in% colnames(pred_df))
    if (!check_names) {
      print("Check data and prediction column names, they don't match")
      stop()
    } else {
      pred_df_bind <- pred_df %>%
        dplyr::select(., colnames(vast_data_out))
      # # We only need one observation for each of the times...
      pred_df_bind <- pred_df %>%
        dplyr::select(., colnames(vast_data_out)) %>%
        distinct(., ID, .keep_all = TRUE)
      vast_data_out <- rbind(vast_data_out, pred_df_bind)
    }
  }

  # Save and return it
  saveRDS(vast_data_out, file = paste(out_dir, "vast_data.rds", sep = "/"))
  return(vast_data_out)
}

#' @title Make VAST sample dataset
#'
#' @description This function creates a VAST sample dataset to pass into calls to `VAST::fit_model`.
#'
#' @param vast_seasonal_data = Description
#' @param out_dir = Description
#'
#' @return A sample dataframe that includes all of the "sample" or species occurrence information. This file is also saved in out_dir.
#'
#' @export

make_vast_sample_data <- function(vast_seasonal_data, fit_seasons, out_dir) {

  # For debugging
  if (FALSE) {
    tar_load(vast_seasonal_data)
    out_dir <- here::here("scratch/aja/targets_flow/data/dfo/combined")
  }

  # Select columns we want from the "full" vast_seasonal_data dataset. Area swept Marine fish diversity on the Scotian Shelf, Canada
  vast_samp_dat <- data.frame(
    "Year" = as.numeric(vast_seasonal_data$VAST_YEAR_SEASON) - 1,
    "Lat" = vast_seasonal_data$DECDEG_BEGLAT,
    "Lon" = vast_seasonal_data$DECDEG_BEGLON,
    "Biomass" = vast_seasonal_data$BIOMASS,
    "Swept" = ifelse(vast_seasonal_data$SURVEY == "NMFS", 0.0384, 0.0404),
    "Pred_TF" = vast_seasonal_data$PredTF
  )

  # Save and return it
  saveRDS(vast_samp_dat, file = paste(out_dir, "vast_sample_data.rds", sep = "/"))
  return(vast_samp_dat)
}

#' @title Make VAST covariate dataset
#'
#' @description This function creates a VAST covariate dataset to pass into calls to `VAST::fit_model`.
#'
#' @param vast_seasonal_data = Description
#' @param rescale = Logical indicating whether or not the covariates should be rescaled.
#' @param out_dir = Description
#'
#' @return A sample dataframe that includes all of the covariate information at each unique sample. This file is also saved in out_dir.
#'
#' @export

make_vast_covariate_data <- function(vast_seasonal_data, out_dir) {

  # For debugging
  if (FALSE) {
    tar_load(vast_seasonal_data)
    rescale <-
      out_dir <- here::here("scratch/aja/targets_flow/data/dfo/combined")
  }

  # Some work to make sure that we don't allow covariates for the "DUMMY" observations to be used at the knots...
  vast_seasonal_data_temp <- vast_seasonal_data

  # Select columns we want from the "full" vast_seasonal_data dataset
  vast_cov_dat <- data.frame(
    "Year" = as.numeric(vast_seasonal_data_temp$VAST_YEAR_SEASON) - 1,
    "Year_Cov" = vast_seasonal_data_temp$VAST_YEAR_COV,
    "Season" = vast_seasonal_data_temp$VAST_SEASON,
    "Depth" = vast_seasonal_data_temp$Depth,
    "SST_seasonal" = vast_seasonal_data_temp$SST_seasonal,
    "BT_seasonal" = vast_seasonal_data_temp$BT_seasonal,
    # "BS_seasonal" = vast_seasonal_data_temp$BS_seasonal,
    # "SS_seasonal" = vast_seasonal_data_temp$SS_seasonal,
    "Lat" = vast_seasonal_data_temp$DECDEG_BEGLAT,
    "Lon" = vast_seasonal_data_temp$DECDEG_BEGLON
  )

  # Save and return
  saveRDS(vast_cov_dat, file = paste(out_dir, "vast_covariate_data.rds", sep = "/"))
  return(vast_cov_dat)
}

#' @title Make VAST catachability
#'
#' @description This function creates a VAST catachability dataset to pass into calls to `VAST::fit_model`.
#'
#' @param vast_seasonal_data = Description
#' @param out_dir = Description
#'
#' @return A sample dataframe that includes all of the covariate information at each unique sample. This file is also saved in out_dir.
#'
#' @export

make_vast_catchability_data <- function(vast_seasonal_data, out_dir) {

  # For debugging
  if (FALSE) {
    vast_seasonal_data
    out_dir <- here::here("scratch/aja/targets_flow/data/dfo/combined")
  }

  # Select columns we want from the "full" vast_seasonal_data dataset
  vast_catch_dat <- data.frame(
    "Year" = as.numeric(vast_seasonal_data$VAST_YEAR_SEASON) - 1,
    "Year_Cov" = vast_seasonal_data$VAST_YEAR_COV,
    "Season" = vast_seasonal_data$VAST_SEASON,
    "Lat" = vast_seasonal_data$DECDEG_BEGLAT,
    "Lon" = vast_seasonal_data$DECDEG_BEGLON,
    "Survey" = factor(vast_seasonal_data$SURVEY, levels = c("NMFS", "DFO"))
  )

  # Save and return it
  saveRDS(vast_catch_dat, file = paste(out_dir, "vast_catchability_data.rds", sep = "/"))
  return(vast_catch_dat)
}

#' @title Read in shapefile
#'
#' @description A short function to read in a shapefile given a file path
#'
#' @param polyshape_path = File path to geospatial vector polygon or multipolygon file with .shp extension, specifying the location and shape of the area of interest.
#'
#' @return SF poylgon
#'
#' @export

read_polyshape <- function(polyshape_path) {

  # For debugging
  if (FALSE) {
    polyshape_path <- "~/Box/RES_Data/Shapefiles/NELME_regions/NELME_sf.shp"
  }

  # Read in polygon shapefile from file_path
  shapefile <- st_read(polyshape_path)
  # shapefile$geometry <- shapefile$geometry %>%
  #   s2::s2_rebuild() %>%
  #   sf::st_as_sfc()
  
  # Return it
  return(shapefile)
}

####
#' @title Make VAST extrapolation grid settings from a shapefile
#'
#' @description Create a list of with information defining the extrapolation grid and used by subsequent VAST functions, leveraging code here: https://github.com/James-Thorson-NOAA/VAST/wiki/Creating-an-extrapolation-grid.
#'
#' @param region_shapefile = A geospatial vector sf polygon file, specifying the location and shape of the area of of spatial domain.
#' @param index_shapes = Either NULL or a polygon or multipolygon sf object specifying sub regions of interest. If specified, the function overlays the sf poly/multipolygon with grid locations and assigns each of the grid locations to one (or more) of the polygon areas. The names of these regions should match the names used as `strata.limits`.
#' @param cell_size = The size of grid in meters (since working in UTM). This will control the resolution of the extrapolation grid.
#'
#' @return Tagged list containing extrapolation grid settings needed to fit a VAST model of species occurrence.
#'
#' @export

vast_make_extrap_grid <- function(region_shapefile, index_shapes, cell_size) {

  # For debugging
  if (FALSE) {
    region_shapefile<- st_read(here::here("data/supporting/region_shapefile/full_survey_region.shp"))  
    tar_load(index_shapefiles)
    index_shapes <- index_shapefiles
    cell_size <- NULL
    grid_dim_km = c(10, 10)
    maximum_distance_from_sample = NULL

    region_shapefile <- menh_chull
    index_shapes <- menh_regions_out
    cell_size <- 5000
  }

  # Transform crs of shapefile to common WGS84 lon/lat format.
  region_wgs84 <- region_shapefile

  # Get UTM zone
  lon <- sum(st_bbox(region_wgs84)[c(1, 3)]) / 2
  utm_zone <- floor((lon + 180) / 6) + 1

  # Transform to the UTM zone
  crs_utm <- st_crs(paste0("+proj=utm +zone=", utm_zone, " +ellps=WGS84 +datum=WGS84 +units=m +no_defs "))
  region_utm <- st_transform(region_wgs84, crs = crs_utm)

  # Make extrapolation grid with sf
    region_grid <- st_as_sf(st_make_grid(region_utm, cellsize = cell_size, what = "centers"), crs = crs_utm)

    # Now get only the points that fall within the shape polygon
    points_keep <- data.frame("pt_row" = seq(from = 1, to = nrow(region_grid), by = 1), "in_out" = st_intersects(region_grid, region_utm, sparse = FALSE))
    region_grid <- region_grid %>%
      mutate(., "in_poly" = st_intersects(region_grid, region_utm, sparse = FALSE)) %>%
      filter(., in_poly == TRUE)

    # Convert back to WGS84 lon/lat, as that is what VAST expects.
    extrap_grid <- region_grid %>%
      st_transform(., crs = 4326)

    # Adding in the a strata/region component for stratified abundance. This will depend on index_shapes input.
    if (!is.null(index_shapes)) {
      extrap_grid <- extrap_grid %>%
        st_join(., index_shapes, join = st_within) %>%
        mutate(.,
          "Lon" = as.numeric(st_coordinates(.)[, 1]),
          "Lat" = as.numeric(st_coordinates(.)[, 2])
        ) %>%
        st_drop_geometry() %>%
        dplyr::select(., Lon, Lat, Region) %>%
        mutate(.,
          Area_km2 = ((cell_size / 1000)^2),
          STRATA = factor(Region, levels = index_shapes$Region, labels = index_shapes$Region)
        )
    } else {
      extrap_grid <- extrap_grid %>%
        mutate(.,
          "Lon" = as.numeric(st_coordinates(.)[, 1]),
          "Lat" = as.numeric(st_coordinates(.)[, 2])
        ) %>%
        st_drop_geometry() %>%
        dplyr::select(., Lon, Lat) %>%
        mutate(., Area_km2 = ((cell_size / 1000)^2))
    }

  # Return it
  return(extrap_grid)
}

####
#' @title Make VAST model settings
#'
#' @description Create a list of model settings needed to fit a VAST model for species occurrence, largely copied from VAST::make_settings
#'
#' @param extrap_grid = User created extrapolation grid from vast_make_extrap_grid.
#' @param n_knots = Number of knots to use when geenrating INLA mesh.
#' @param FieldConfig = A vector defining the number of spatial (Omega) and spatio-temporal (Epsilon) factors to include in the model for each of the linear predictors. For each factor, possible values range from 0 (which effectively turns off a given factor), to the number of categories being modeled. If FieldConfig < number of categories, VAST estimates common factors and then loading matrices.
#' @param RhoConfig = A vector defining the temporal structure of intercepts (Beta) and spatio-temporal (Epsilon) variation for each of the linear predictors. See `VAST::make_data` for options.
#' @param bias.correct = Logical boolean determining if Epsilon bias-correction should be done.
#' @param knot_method = The method to use in placing the knots, either "grid" or "samples". Grid places knots on a regular grid, while samples places knots corresponding to the intensity of observations.
#' @param inla_method - The method to pass to INLA to create the mesh, allowing use of the "barrier" mesh method.
#' @param Options = Tagged vector to turn on or off specific options (e.g., SD_site_logdensity, Effective area, etc)
#' @param strata.limits = Dataframe specifying the strata to use. This should align with the `index_shapefile` and ultimately is used by VAST for calculating stratified biomass indices.
#'
#' @return Tagged list containing settings needed to fit a VAST model of species occurrence.
#'
#' @export

vast_make_settings <- function(extrap_grid, n_knots, FieldConfig, RhoConfig, ObsModel, OverdispersionConfig, bias.correct, knot_method, inla_method, Options, strata.limits, version = FishStatsUtils::get_latest_version()) {

  # For debugging
  if (FALSE) {
    tar_load(vast_extrap_grid)
    extrap_grid <- vast_extrap_grid
    FieldConfig <- c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 1)
    RhoConfig <- c("Beta1" = 3, "Beta2" = 3, "Epsilon1" = 2, "Epsilon2" = 2)
    OverdispersionConfig <- c(0, 0)
    bias.correct <- FALSE
    Options <- c("Calculate_Range" = TRUE)
    strata.limits <- strata_use
    n_knots <- 400
    knot_method <- "samples"
    inla_method <- "Barrier"
  }

  # Run FishStatsUtils::make_settings

  settings_out <- make_settings(n_x = n_knots, Region = "User", purpose = "index2", FieldConfig = FieldConfig, RhoConfig = RhoConfig, ObsModel = ObsModel, OverdispersionConfig = OverdispersionConfig, bias.correct = bias.correct, knot_method = knot_method, treat_nonencounter_as_zero = FALSE, strata.limits = strata.limits, Version = version)
  settings_out$Method <- inla_method

  # Adjust options?
  options_new <- settings_out$Options
  if (!is.null(Options)) {
    for (i in seq_along(Options)) {
      options_adjust_i <- Options[i]
      options_new[[which(names(options_new) == names(options_adjust_i))]] <- options_adjust_i
    }
    settings_out <- make_settings(n_x = n_knots, Region = "User", purpose = "index2", FieldConfig = FieldConfig, RhoConfig = RhoConfig, ObsModel = ObsModel, OverdispersionConfig = OverdispersionConfig, bias.correct = bias.correct, knot_method = knot_method, treat_nonencounter_as_zero = FALSE, strata.limits = strata.limits, Options = options_new, Version = version)
    settings_out$Method <- inla_method
  }

  # Return it
  return(settings_out)
}

####
#' @title Make VAST spatial info
#'
#' @description Create a tagged list with VAST spatial information needed
#'
#' @param extrap_grid = User created extrapolation grid from vast_make_extrap_grid.
#' @param vast_settings = A tagged list of VAST settings. Could be either returned by `vast_make_settings` or `make_settings`.
#' @param tidy_mod_data = A tidy model datafame with all the information (tows, habitat covariates, species occurrences) needed to fit a species distribution model.
#' @param out_dir = Directory for VAST to look for certain spatial objects (e.g., knots and extrapolation grid .rds files). This is helpful so that new knots aren't generated every run.
#'
#' @return Returns a tagged list with extrapolation and spatial info in different slots.
#'
#' @export

vast_make_spatial_lists <- function(extrap_grid, vast_settings, tidy_mod_data, out_dir) {

  # For debugging
  if (FALSE) {
    tar_load(vast_extrap_grid)
    extrap_grid <- vast_extrap_grid
    tar_load(vast_settings)
    tar_load(tidy_mod_data)
    inla_method <- "Barrier"
    out_dir <- here::here()
  }


  # Run FishStatsUtiles::make_extrapolation_info
  vast_extrap_info <- make_extrapolation_info(Region = vast_settings$Region, strata.limits = vast_settings$strata.limits, input_grid = extrap_grid, DirPath = out_dir)

  # Run FishStatsUtils::make_spatial_info
  vast_spatial_info <- make_spatial_info(n_x = vast_settings$n_x, Lon_i = tidy_mod_data$DECDEG_BEGLON, Lat_i = tidy_mod_data$DECDEG_BEGLAT, Extrapolation_List = vast_extrap_info, knot_method = vast_settings$knot_method, Method = vast_settings$Method, grid_size_km = vast_settings$grid_size_km, fine_scale = vast_settings$fine_scale, DirPath = out_dir, Save_Results = TRUE)

  # Combine into one list of lists
  spatial_lists_out <- list(vast_extrap_info, vast_spatial_info)
  names(spatial_lists_out) <- c("Extrapolation_List", "Spatial_List")
  return(spatial_lists_out)
}


####
#' @title Reduce VAST prediction dataframe from regular grid to knot locations
#'
#' @description Reduce VAST prediction dataframe to knot locations, where covariates at grid locations are 'snapped' to nearest knot location. As of 2/7/2022 not needed in the workflow.
#'
#' @param vast_predict_df = A dataframe with gridded prediction locations and covariate data information. This is generated by `make_vast_predict_df` in this workflow.
#' @param vast_spatial_lists = A list with VAST extrapolation info and spatial info as generated by `vast_make_spatial_lists`.
#' @param out_dir = Output directory to save reduced prediction dataframe.
#'
#' @return Returns a prediction dataframe with knot locations (instead of original locations) and covariates.
#'
#' @export

reduce_vast_predict_df <- function(vast_predict_df = vast_predict_df, vast_spatial_lists = vast_spatial_lists, out_dir = here::here("data/predict")) {

  # For debugging
  if (FALSE) {
    tar_load(vast_predict_df)
    tar_load(vast_spatial_lists)
  }

  # Knots_sf
  knots_info <- vast_spatial_lists$Spatial_List
  knots_sf <- st_as_sf(data.frame(knots_info$loc_x), coords = c("E_km", "N_km"), crs = attributes(knots_info$loc_i)$projCRS)

  # Get unique prediction locations and assign each prediction location to its nearest knot?
  pred_df_temp <- vast_predict_df %>%
    distinct(., DECDEG_BEGLON, DECDEG_BEGLAT)
  pred_sf <- points_to_sf(pred_df_temp) %>%
    st_transform(., crs = st_crs(knots_sf))

  pred_nearest_knot <- pred_sf %>%
    mutate(., "Nearest_knot" = st_nearest_feature(x = ., y = knots_sf)) %>%
    st_drop_geometry()

  # Merge this with full prediction dataset
  pred_df_out <- vast_predict_df %>%
    left_join(., pred_nearest_knot)

  # Average covariate values based on nearest knot location and output reduced dataframe
  pred_df_out <- pred_df_out %>%
    distinct(., ID, DATE, Nearest_knot, .keep_all = TRUE) %>%
    dplyr::select(-Nearest_knot)

  return(pred_df_out)
}

####
#' @title Make VAST covariate effect objects
#'
#' @description Create covariate effects for both linear predictors. These are passed to the `vast_build_sdm` function.
#'
#' @param X1_coveff_vec = A vector specifying the habitat covariate effects for first linear predictor.
#' @param X2_coveff_vec = A vector specifying the habitat covariate effects for second linear predictor.
#' @param Q1_coveff_vec = A vector specifying the catchability covariate effects for first linear predictor.
#' @param Q2_coveff_vec = A vector specifying the catchability covariate effects for second linear predictor.
#'
#' @return A list with covariate effects for the habitat covariates and first linear predictor (first list slot), habitat covariates and second linear predictor (second list slot), catchability covariates and first linear predictor (third slot) and catchability covariates and second linear predictor (fourth slot).
#'
#' @export

vast_make_coveff <- function(X1_coveff_vec, X2_coveff_vec, Q1_coveff_vec, Q2_coveff_vec) {

  # For debugging
  if (FALSE) {
    X1_coveff_vec <- c(2, 3, 3, 2, rep(3, 32))
    X2_coveff_vec <- c(2, 3, 3, 2, rep(3, 32))
    Q1_coveff_vec <- NULL
    Q2_coveff_vec <- NULL
  }

  # Combine into a list and name it
  if (is.null(Q1_coveff_vec) | is.null(Q2_coveff_vec)) {
    coveff_out <- list("X1config_cp" = matrix(X1_coveff_vec, nrow = 1), "X2config_cp" = matrix(X2_coveff_vec, nrow = 1), "Q1config_k" = NULL, "Q2config_k" = NULL)
  } else {
    coveff_out <- list("X1config_cp" = matrix(X1_coveff_vec, nrow = 1), "X2config_cp" = matrix(X2_coveff_vec, nrow = 1), "Q1config_k" = matrix(Q1_coveff_vec, nrow = 1), "Q2config_k" = matrix(Q2_coveff_vec, nrow = 1))
  }

  # Return it
  return(coveff_out)
}

####
#' @title Build VAST SDM
#'
#' @description Build VAST species distribution model, without running it. This can be helpful to check settings before running `vast_fit_sdm`. Additionally, it can be helpful for making subsequent modifications, particularly to mapping.
#'
#' @param settings = A tagged list with the settings for the model, created with `vast_make_settings`.
#' @param extrap_grid = An extrapolation grid, created with `vast_make_extrap_grid`.
#' @param Method = A character string specifying which Method to use when making the mesh.
#' @param sample_data = A data frame with the biomass sample data for each species at each tow.
#' @param covariate_data = A data frame with the habitat covariate data for each tow.
#' @param X1_formula = A formula for the habitat covariates and first linear predictor.
#' @param X2_formula = A formula for the habitat covariates and second linear predictor.
#' @param X_contrasts = A tagged list specifying the contrasts to use for factor covariates in the model.
#' @param Xconfig_list = A tagged list specifying the habitat and catchability covariate effects for first and second linear predictors.
#' @param catchability_data = A data frame with the catchability data for every sample
#' @param Q1_formula = A formula for the catchability covariates and first linear predictor.
#' @param Q2_formula = A formula for the catchability covariates and second linear predictor.
#' @param index_shapes = An sf polygon or multipolygon object with rows for each of the regions of interest.
#' @param spatial_info_dir = A directory holding the spatial information, including knots and extrapolation locations.
#'
#' @return A VAST `fit_model` object, with the inputs and built TMB object components.
#'
#' @export

vast_build_sdm <- function(settings, spatial_list = NULL, extrapolation_list = NULL, extrap_grid = NULL, sample_data, covariate_data, X1_formula, X2_formula, X_contrasts, Xconfig_list, catchability_data, Q1_formula, Q2_formula, index_shapes = NULL, spatial_info_dir = NULL) {

  # For debugging
  if (FALSE) {
    library(VAST)
    library(tidyverse)
    library(stringr)

    # Seasonal
    tar_load(vast_settings)
    settings <- vast_settings
    tar_load(vast_extrap_grid)
    extrap_grid <- vast_extrap_grid
    tar_load(vast_sample_data)
    sample_data <- vast_sample_data
    tar_load(vast_covariate_data)
    covariate_data <- vast_covariate_data
    X1_formula <- hab_formula
    X2_formula <- hab_formula
    hab_env_coeffs_n <- hab_env_coeffs_n
    tar_load(vast_catchability_data)
    catchability_data <- vast_catchability_data
    catch_formula <- ~Survey
    Q1_formula <- catch_formula
    Q2_formula <- catch_formula
    X_contrasts <- list(Season = contrasts(vast_covariate_data$Season, contrasts = FALSE), Year_Cov = contrasts(vast_covariate_data$Year_Cov, contrasts = FALSE))
    # X_contrasts = list(Year_Cov = contrasts(vast_covariate_data$Year_Cov, contrasts = FALSE))
    tar_load(vast_coveff)
    Xconfig_list <- vast_coveff
    tar_load(index_shapefiles)
    index_shapes <- index_shapefiles
    spatial_info_dir <- here::here("")

    # Annual
    tar_load(vast_settings)
    settings <- vast_settings
    tar_load(vast_extrap_grid)
    extrap_grid <- vast_extrap_grid
    tar_load(vast_sample_data)
    sample_data <- vast_sample_data
    tar_load(vast_covariate_data)
    covariate_data <- vast_covariate_data
    X1_formula <- hab_formula
    X2_formula <- hab_formula
    hab_env_coeffs_n <- hab_env_coeffs_n
    tar_load(vast_catchability_data)
    catchability_data <- vast_catchability_data
    catch_formula <- ~0
    Q1_formula <- catch_formula
    Q2_formula <- catch_formula
    X_contrasts <- list(Year_Cov = contrasts(vast_covariate_data$Year_Cov, contrasts = FALSE))
    tar_load(vast_coveff)
    Xconfig_list <- vast_coveff
    tar_load(index_shapefiles)
    index_shapes <- index_shapefiles
  }

  # Check names
  samp_dat_names <- c("Lat", "Lon", "Year", "Biomass", "Swept", "Pred_TF")
  if (!(all(samp_dat_names %in% names(sample_data)))) {
    stop(paste("Check names in sample data. Must include:", paste0(samp_dat_names, collapse = ","), sep = " "))
  }

  # Covariate data frame names
  if (!is.null(covariate_data)) {
    cov_dat_names1 <- unlist(str_extract_all(format(X1_formula), boundary("word"))[[2]])

    # Remove some stuff associated with the splines...
    spline_words <- c("bs", "degree", "df", "knots", "TRUE", "intercept", unique(as.numeric(unlist(str_extract_all(format(X1_formula), pattern = "[0-9]+", simplify = TRUE)))), "FALSE")
    cov_dat_names1 <- cov_dat_names1[-which(cov_dat_names1 %in% spline_words)]
    cov_dat_names2 <- unlist(str_extract_all(format(X2_formula), boundary("word"))[[2]])
    cov_dat_names2 <- cov_dat_names2[-which(cov_dat_names2 %in% spline_words)]
    cov_dat_names_all <- unique(c(cov_dat_names1, cov_dat_names2))
    if (!(all(cov_dat_names_all %in% names(covariate_data)))) {
      print(names(covariate_data))
      print(names(cov_dat_names_all))
      stop(paste("Check names in covariate data. Must include", paste0(cov_dat_names_all, collapse = ","), sep = " "))
    }
  }

  if (!(all(c("X1config_cp", "X2config_cp", "Q1config_k", "Q2config_k") %in% names(Xconfig_list)))) {
    stop(paste("Check names of Xconfig_list. Must be", paste0(c("X1config_cp", "X2config_cp", "Q1config_k", "Q2config_k"), collapse = ","), sep = ""))
  }

  # Run VAST::fit_model with correct info and settings
  # vast_build_out <- fit_model_aja("settings" = settings, "Method" = settings$Method, "input_grid" = extrap_grid, "Lat_i" = sample_data[, "Lat"], "Lon_i" = sample_data[, "Lon"], "t_i" = sample_data[, "Year"], "c_i" = rep(0, nrow(sample_data)), "b_i" = sample_data[, "Biomass"], "a_i" = sample_data[, "Swept"], "PredTF_i" = sample_data[, "Pred_TF"], "X1config_cp" = Xconfig_list[["X1config_cp"]], "X2config_cp" = Xconfig_list[["X2config_cp"]], "covariate_data" = covariate_data, "X1_formula" = X1_formula, "X2_formula" = X2_formula, "X_contrasts" = X_contrasts, "catchability_data" = catchability_data, "Q1_formula" = Q1_formula, "Q2_formula" = Q2_formula, "Q1config_k" = Xconfig_list[["Q1config_k"]], "Q2config_k" = Xconfig_list[["Q2config_k"]], "newtonsteps" = 1, "getsd" = TRUE, "getReportCovariance" = TRUE, "run_model" = FALSE, "test_fit" = FALSE, "Use_REML" = FALSE, "getJointPrecision" = TRUE, "index_shapes" = index_shapes, "DirPath" = spatial_info_dir)

  vast_build_out <- fit_model("settings" = settings, "Method" = settings$Method, "spatial_list" = spatial_list, "extrapolation_list" = extrapolation_list, "input_grid" = extrap_grid, "Lat_i" = sample_data[, "Lat"], "Lon_i" = sample_data[, "Lon"], "t_i" = sample_data[, "Year"], "c_i" = rep(0, nrow(sample_data)), "b_i" = sample_data[, "Biomass"], "a_i" = sample_data[, "Swept"], "PredTF_i" = sample_data[, "Pred_TF"], "X1config_cp" = Xconfig_list[["X1config_cp"]], "X2config_cp" = Xconfig_list[["X2config_cp"]], "covariate_data" = covariate_data, "X1_formula" = X1_formula, "X2_formula" = X2_formula, "X_contrasts" = X_contrasts, "catchability_data" = catchability_data, "Q1_formula" = Q1_formula, "Q2_formula" = Q2_formula, "Q1config_k" = Xconfig_list[["Q1config_k"]], "Q2config_k" = Xconfig_list[["Q2config_k"]], "newtonsteps" = 1, "getsd" = TRUE, "getReportCovariance" = TRUE, "run_model" = FALSE, "test_fit" = FALSE, "Use_REML" = FALSE, "getJointPrecision" = TRUE)

  # Return it
  return(vast_build_out)
}

####
#' @title Adjust VAST SDM
#'
#' @description Make adjustments to VAST SDM and the model returned in `vast_build_sdm`. This can either be the exact same as the one built using `vast_build_sdm`, or it can update that model with adjustments provided in a tagged list.
#'
#' @param vast_build = A VAST `fit_model` object.
#' @param index_shapes = An sf polygon or multipolygon object with rows for each of the regions of interest.
#' @param spatial_info_dir = A directory holding the spatial information, including knots and extrapolation locations.
#' @param adjustments = Either NULL (default) or a tagged list identifying adjustments that should be made to the vast_build `fit_model` object. If NULL, the identical model defined by the `vast_build` is run and fitted.
#'
#' @return A VAST fit_model object, with the inputs and built TMB object components.
#'
#' @export

vast_make_adjustments <- function(vast_build, index_shapes = NULL, spatial_info_dir = NULL, adjustments = NULL) {

  # For debugging
  if (FALSE) {
    tar_load(vast_build0)
    vast_build <- vast_build0
    tar_load(vast_covariate_data)
    adjustments <- list("log_sigmaXi1_cp" = factor(c(rep(1, length(unique(fit_seasons))), rep(4, nlevels(vast_covariate_data$Year_Cov)), rep(NA, gam_degree * hab_env_coeffs_n))), "log_sigmaXi2_cp" = factor(c(rep(1, length(unique(fit_seasons))), rep(4, nlevels(vast_covariate_data$Year_Cov)), rep(NA, gam_degree * hab_env_coeffs_n))), "lambda1_k" = factor(c(1, NA)), "lambda2_k" = factor(c(1, NA)))
    tar_load(index_shapefiles)
    index_shapes <- index_shapefiles
  }

  # If no adjustments are needed, just need to pull information from vast_build and then set "run_model" to TRUE
  if (is.null(adjustments)) {
    # vast_build_adjust_out <- fit_model_aja("settings" = vast_build$settings, "input_grid" = vast_build$input_args$data_args_input$input_grid, "Method" = vast_build$settings$Method, "Lat_i" = vast_build$data_frame[, "Lat_i"], "Lon_i" = vast_build$data_frame[, "Lon_i"], "t_i" = vast_build$data_frame[, "t_i"], "c_iz" = vast_build$data_frame[, "c_iz"], "b_i" = vast_build$data_frame[, "b_i"], "a_i" = vast_build$data_frame[, "a_i"], "PredTF_i" = vast_build$data_list[["PredTF_i"]], "X1config_cp" = vast_build$input_args$data_args_input[["X1config_cp"]], "X2config_cp" = vast_build$input_args$data_args_input[["X2config_cp"]], "covariate_data" = vast_build$input_args$data_args_input$covariate_data, "X1_formula" = vast_build$input_args$data_args_input$X1_formula, "X2_formula" = vast_build$input_args$data_args_input$X2_formula, "X_contrasts" = vast_build$input_args$data_args_input$X_contrasts, "catchability_data" = vast_build$input_args$data_args_input$catchability_data, "Q1_formula" = vast_build$input_args$data_args_input$Q1_formula, "Q2_formula" = vast_build$input_args$data_args_input$Q2_formula, "Q1config_k" = vast_build$input_args$data_args_input[["Q1config_cp"]], "Q2config_k" = vast_build$input_args$data_args_input[["Q2config_k"]], "newtonsteps" = 1, "getsd" = TRUE, "getReportCovariance" = TRUE, "run_model" = FALSE, "test_fit" = FALSE, "Use_REML" = FALSE, "getJointPrecision" = vast_build$input_args$extra_args$getJointPrecision, "index_shapes" = index_shapes, "DirPath" = spatial_info_dir)

     vast_build_adjust_out <- fit_model("settings" = vast_build$settings, "Method" = vast_build$settings$Method, "spatial_list" = vast_build$spatial_list, "extrapolation_list" = vast_build$extrapolation_list, "input_grid" = vast_build$input_args$data_args_input$input_grid, "Lat_i" = vast_build$data_frame[, "Lat_i"], "Lon_i" = vast_build$data_frame[, "Lon_i"], "t_i" = vast_build$data_frame[, "t_i"], "c_iz" = vast_build$data_frame[, "c_iz"], "b_i" = vast_build$data_frame[, "b_i"], "a_i" = vast_build$data_frame[, "a_i"], "PredTF_i" = vast_build$data_list[["PredTF_i"]], "X1config_cp" = vast_build$input_args$data_args_input[["X1config_cp"]], "X2config_cp" = vast_build$input_args$data_args_input[["X2config_cp"]], "covariate_data" = vast_build$input_args$data_args_input$covariate_data, "X1_formula" = vast_build$input_args$data_args_input$X1_formula, "X2_formula" = vast_build$input_args$data_args_input$X2_formula, "X_contrasts" = vast_build$input_args$data_args_input$X_contrasts, "catchability_data" = vast_build$input_args$data_args_input$catchability_data, "Q1_formula" = vast_build$input_args$data_args_input$Q1_formula, "Q2_formula" = vast_build$input_args$data_args_input$Q2_formula, "Q1config_k" = vast_build$input_args$data_args_input[["Q1config_cp"]], "Q2config_k" = vast_build$input_args$data_args_input[["Q2config_k"]], "newtonsteps" = 1, "getsd" = TRUE, "getReportCovariance" = TRUE, "run_model" = FALSE, "test_fit" = FALSE, "Use_REML" = FALSE, "getJointPrecision" = vast_build$input_args$extra_args$getJointPrecision)

  }

  # If there are adjustments, need to make those and then re run model.
  if (!is.null(adjustments)) {
    # Check names -- trying to think of what the possible adjustment flags would be in the named list
    adjust_names <- c("FieldConfig", "RhoConfig", "X1_formula", "X2_formula", "X1config_cp", "X2config_cp", "X_contrasts", "log_sigmaXi1_cp", "log_sigmaXi2_cp", "lambda1_k", "lambda2_k", "Q1_formula", "Q2_formula", "Q1config_k", "Q2config_k")

    if (!(all(names(adjustments) %in% adjust_names))) {
      stop(paste("Check names in adjustment list. Must be one of", paste0(adjust_names, collapse = ","), sep = " "))
    }

    # First options are going to be in the settings bit..
    if (any(names(adjustments) %in% c("FieldConfig", "RhoConfig"))) {
      # Get just the settings adjustments
      settings_adjusts <- names(adjustments)[which(names(adjustments) %in% names(vast_build$settings))]

      for (i in seq_along(settings_adjusts)) {
        setting_adjust_i <- settings_adjusts[i]
        vast_build$settings[[{{ setting_adjust_i }}]] <- adjustments[[{{ setting_adjust_i }}]]
      }
    }

    # A lot of stuff is going to be in the `vast_build$input_args$data_args_input` object
    if (any(names(adjustments) %in% names(vast_build$input_args$data_args_input))) {

      # Get just the data args adjustments
      data_adjusts <- names(adjustments)[which(names(adjustments) %in% names(vast_build$input_args$data_args_input))]

      for (i in seq_along(data_adjusts)) {
        data_adjust_i <- data_adjusts[i]
        vast_build$input_args$data_args_input[[{{ data_adjust_i }}]] <- adjustments[[{{ data_adjust_i }}]]
      }
    }

    # Only other adjustment (for now) is Map.
    if (any(names(adjustments) %in% c("log_sigmaXi1_cp", "log_sigmaXi2_cp", "lambda1_k", "lambda2_k"))) {
      # Get the original, which we can then edit...
      map_adjust_out <- vast_build$tmb_list$Map

      # Get just the map adjustment names
      map_adjusts <- names(adjustments)[which(names(adjustments) %in% names(vast_build$tmb_list$Map))]

      # Loop over them
      for (i in seq_along(map_adjusts)) {
        map_adjust_i <- map_adjusts[i]
        map_adjust_out[[{{ map_adjust_i }}]] <- adjustments[[{{ map_adjust_i }}]]
      }
    }

    # Now, re-build and fit model. This is slightly different if we have changed map or not...
    if (any(names(adjustments) %in% c("log_sigmaXi1_cp", "log_sigmaXi2_cp", "lambda1_k", "lambda2_k"))) {
      # Adding Map argument
      # vast_build_adjust_out <- fit_model_aja("settings" = vast_build$settings, "input_grid" = vast_build$input_args$data_args_input$input_grid, "Method" = vast_build$settings$Method, "Lat_i" = vast_build$data_frame[, "Lat_i"], "Lon_i" = vast_build$data_frame[, "Lon_i"], "t_i" = vast_build$data_frame[, "t_i"], "c_iz" = vast_build$data_frame[, "c_iz"], "b_i" = vast_build$data_frame[, "b_i"], "a_i" = vast_build$data_frame[, "a_i"], "PredTF_i" = vast_build$data_list[["PredTF_i"]], "X1config_cp" = vast_build$input_args$data_args_input[["X1config_cp"]], "X2config_cp" = vast_build$input_args$data_args_input[["X2config_cp"]], "covariate_data" = vast_build$input_args$data_args_input$covariate_data, "X1_formula" = vast_build$input_args$data_args_input$X1_formula, "X2_formula" = vast_build$input_args$data_args_input$X2_formula, "X_contrasts" = vast_build$input_args$data_args_input$X_contrasts, "catchability_data" = vast_build$input_args$data_args_input$catchability_data, "Q1_formula" = vast_build$input_args$data_args_input$Q1_formula, "Q2_formula" = vast_build$input_args$data_args_input$Q2_formula, "Q1config_k" = vast_build$input_args$data_args_input[["Q1config_k"]], "Q2config_k" = vast_build$input_args$data_args_input[["Q2config_k"]], "Map" = map_adjust_out, "newtonsteps" = 1, "getsd" = TRUE, "getReportCovariance" = TRUE, "run_model" = FALSE, "test_fit" = FALSE, "Use_REML" = FALSE, "getJointPrecision" = FALSE, "index_shapes" = index_shapes, "DirPath" = spatial_info_dir)
      vast_build_adjust_out <- fit_model("settings" = vast_build$settings, "spatial_list" = vast_build$spatial_list, "extrapolation_list" = vast_build$extrapolation_list, "input_grid" = vast_build$input_args$data_args_input$input_grid, "Method" = vast_build$settings$Method, "Lat_i" = vast_build$data_frame[, "Lat_i"], "Lon_i" = vast_build$data_frame[, "Lon_i"], "t_i" = vast_build$data_frame[, "t_i"], "c_iz" = vast_build$data_frame[, "c_iz"], "b_i" = vast_build$data_frame[, "b_i"], "a_i" = vast_build$data_frame[, "a_i"], "PredTF_i" = vast_build$data_list[["PredTF_i"]], "X1config_cp" = vast_build$input_args$data_args_input[["X1config_cp"]], "X2config_cp" = vast_build$input_args$data_args_input[["X2config_cp"]], "covariate_data" = vast_build$input_args$data_args_input$covariate_data, "X1_formula" = vast_build$input_args$data_args_input$X1_formula, "X2_formula" = vast_build$input_args$data_args_input$X2_formula, "X_contrasts" = vast_build$input_args$data_args_input$X_contrasts, "catchability_data" = vast_build$input_args$data_args_input$catchability_data, "Q1_formula" = vast_build$input_args$data_args_input$Q1_formula, "Q2_formula" = vast_build$input_args$data_args_input$Q2_formula, "Q1config_k" = vast_build$input_args$data_args_input[["Q1config_k"]], "Q2config_k" = vast_build$input_args$data_args_input[["Q2config_k"]], "Map" = map_adjust_out, "newtonsteps" = 1, "getsd" = TRUE, "getReportCovariance" = TRUE, "run_model" = FALSE, "test_fit" = FALSE, "Use_REML" = FALSE, "getJointPrecision" = FALSE, "DirPath" = spatial_info_dir)
    } else {
      # No need for Map argument, just build and fit
      # vast_build_adjust_out <- fit_model_aja("settings" = vast_build$settings, "input_grid" = vast_build$input_args$data_args_input$input_grid, "Method" = vast_build$settings$Method, "Lat_i" = vast_build$data_frame[, "Lat_i"], "Lon_i" = vast_build$data_frame[, "Lon_i"], "t_i" = vast_build$data_frame[, "t_i"], "c_iz" = vast_build$data_frame[, "c_iz"], "b_i" = vast_build$data_frame[, "b_i"], "a_i" = vast_build$data_frame[, "a_i"], "PredTF_i" = vast_build$data_list[["PredTF_i"]], "X1config_cp" = vast_build$input_args$data_args_input[["X1config_cp"]], "X2config_cp" = vast_build$input_args$data_args_input[["X2config_cp"]], "covariate_data" = vast_build$input_args$data_args_input$covariate_data, "X1_formula" = vast_build$input_args$data_args_input$X1_formula, "X2_formula" = vast_build$input_args$data_args_input$X2_formula, "X_contrasts" = vast_build$input_args$data_args_input$X_contrasts, "catchability_data" = vast_build$input_args$data_args_input$catchability_data, "Q1_formula" = vast_build$input_args$data_args_input$Q1_formula, "Q2_formula" = vast_build$input_args$data_args_input$Q2_formula, "Q1config_cp" = vast_build$input_args$data_args_input[["Q1config_cp"]], "Q2config_cp" = vast_build$input_args$data_args_input[["Q2config_cp"]], "newtonsteps" = 1, "getsd" = TRUE, "getReportCovariance" = TRUE, "run_model" = FALSE, "test_fit" = FALSE, "Use_REML" = FALSE, "getJointPrecision" = TRUE, "index_shapes" = index_shapes, "DirPath" = spatial_info_dir)
      vast_build_adjust_out <- fit_model("settings" = vast_build$settings, "spatial_list" = vast_build$spatial_list, "extrapolation_list" = vast_build$extrapolation_list, "input_grid" = vast_build$input_args$data_args_input$input_grid, "Method" = vast_build$settings$Method, "Lat_i" = vast_build$data_frame[, "Lat_i"], "Lon_i" = vast_build$data_frame[, "Lon_i"], "t_i" = vast_build$data_frame[, "t_i"], "c_iz" = vast_build$data_frame[, "c_iz"], "b_i" = vast_build$data_frame[, "b_i"], "a_i" = vast_build$data_frame[, "a_i"], "PredTF_i" = vast_build$data_list[["PredTF_i"]], "X1config_cp" = vast_build$input_args$data_args_input[["X1config_cp"]], "X2config_cp" = vast_build$input_args$data_args_input[["X2config_cp"]], "covariate_data" = vast_build$input_args$data_args_input$covariate_data, "X1_formula" = vast_build$input_args$data_args_input$X1_formula, "X2_formula" = vast_build$input_args$data_args_input$X2_formula, "X_contrasts" = vast_build$input_args$data_args_input$X_contrasts, "catchability_data" = vast_build$input_args$data_args_input$catchability_data, "Q1_formula" = vast_build$input_args$data_args_input$Q1_formula, "Q2_formula" = vast_build$input_args$data_args_input$Q2_formula, "Q1config_cp" = vast_build$input_args$data_args_input[["Q1config_cp"]], "Q2config_cp" = vast_build$input_args$data_args_input[["Q2config_cp"]], "newtonsteps" = 1, "getsd" = TRUE, "getReportCovariance" = TRUE, "run_model" = FALSE, "test_fit" = FALSE, "Use_REML" = FALSE, "getJointPrecision" = TRUE, "DirPath" = spatial_info_dir)
    }
  }
  # Return it
  return(vast_build_adjust_out)
}

#' @title Fit VAST SDM
#'
#' @description Fit VAST species distribution model
#'
#' @param vast_build_adjust = A VAST `fit_model` object.
#' @param nice_category_names = A character string to define species/model run and used in naming the fitted files.
#' @param index_shapes = An sf polygon or multipolygon object with rows for each of the regions of interest.
#' @param spatial_info_dir = A directory holding the spatial information, including knots and extrapolation locations.
#' @param out_dir = Output directory to save fitted model object.
#'
#' @return A VAST fit_model object, with the inputs and and outputs, including parameter estimates, extrapolation gid info, spatial list info, data info, and TMB info.
#'
#' @export

vast_fit_sdm <- function(vast_build_adjust, nice_category_names, index_shapes = NULL, spatial_info_dir = NULL, run_final_model = run_final_model, out_dir) {

  # For debugging
  if (FALSE) {
    tar_load(vast_adjust)
    vast_build_adjust <- vast_adjust
    nice_category_names <- nice_category_names
    out_dir <- here::here("results/mod_fits")
    tar_load(index_shapefiles)
    index_shapes <- index_shapefiles
    spatial_info_dir <- here::here("")

    vast_build_adjust <- vast_adjust
    nice_category_names <- nice_category_names
    out_dir <- paste0(res_root, "mod_fits")
    index_shapes <- index_shapefiles
    spatial_info_dir <- here::here("")
  }

  # Build and fit model
  # vast_fit_out <- fit_model_aja("settings" = vast_build_adjust$settings, "input_grid" = vast_build_adjust$input_args$data_args_input$input_grid, "Method" = vast_build_adjust$settings$Method, "Lat_i" = vast_build_adjust$data_frame[, "Lat_i"], "Lon_i" = vast_build_adjust$data_frame[, "Lon_i"], "t_i" = vast_build_adjust$data_frame[, "t_i"], "c_iz" = vast_build_adjust$data_frame[, "c_iz"], "b_i" = vast_build_adjust$data_frame[, "b_i"], "a_i" = vast_build_adjust$data_frame[, "a_i"], "PredTF_i" = vast_build_adjust$data_list[["PredTF_i"]], "X1config_cp" = vast_build_adjust$input_args$data_args_input[["X1config_cp"]], "X2config_cp" = vast_build_adjust$input_args$data_args_input[["X2config_cp"]], "covariate_data" = vast_build_adjust$input_args$data_args_input$covariate_data, "X1_formula" = vast_build_adjust$input_args$data_args_input$X1_formula, "X2_formula" = vast_build_adjust$input_args$data_args_input$X2_formula, "X_contrasts" = vast_build_adjust$input_args$data_args_input$X_contrasts, "catchability_data" = vast_build_adjust$input_args$data_args_input$catchability_data, "Q1_formula" = vast_build_adjust$input_args$data_args_input$Q1_formula, "Q2_formula" = vast_build_adjust$input_args$data_args_input$Q2_formula, "Q1config_cp" = vast_build_adjust$input_args$data_args_input[["Q1config_cp"]], "Q2config_cp" = vast_build_adjust$input_args$data_args_input[["Q2config_cp"]], "Map" = vast_build_adjust$tmb_list$Map, "newtonsteps" = 1, "getsd" = TRUE, "getReportCovariance" = TRUE, "run_model" = run_final_model, "test_fit" = FALSE, "Use_REML" = FALSE, "getJointPrecision" = vast_build_adjust$input_args$extra_args$getJointPrecision, "index_shapes" = index_shapes, "DirPath" = spatial_info_dir)
  vast_fit_out <- fit_model("settings" = vast_build_adjust$settings, "spatial_list" = vast_build_adjust$spatial_list, "extrapolation_list" = vast_build_adjust$extrapolation_list,  "input_grid" = vast_build_adjust$input_args$data_args_input$input_grid, "Method" = vast_build_adjust$settings$Method, "Lat_i" = vast_build_adjust$data_frame[, "Lat_i"], "Lon_i" = vast_build_adjust$data_frame[, "Lon_i"], "t_i" = vast_build_adjust$data_frame[, "t_i"], "c_iz" = vast_build_adjust$data_frame[, "c_iz"], "b_i" = vast_build_adjust$data_frame[, "b_i"], "a_i" = vast_build_adjust$data_frame[, "a_i"], "PredTF_i" = vast_build_adjust$data_list[["PredTF_i"]], "X1config_cp" = vast_build_adjust$input_args$data_args_input[["X1config_cp"]], "X2config_cp" = vast_build_adjust$input_args$data_args_input[["X2config_cp"]], "covariate_data" = vast_build_adjust$input_args$data_args_input$covariate_data, "X1_formula" = vast_build_adjust$input_args$data_args_input$X1_formula, "X2_formula" = vast_build_adjust$input_args$data_args_input$X2_formula, "X_contrasts" = vast_build_adjust$input_args$data_args_input$X_contrasts, "catchability_data" = vast_build_adjust$input_args$data_args_input$catchability_data, "Q1_formula" = vast_build_adjust$input_args$data_args_input$Q1_formula, "Q2_formula" = vast_build_adjust$input_args$data_args_input$Q2_formula, "Q1config_cp" = vast_build_adjust$input_args$data_args_input[["Q1config_cp"]], "Q2config_cp" = vast_build_adjust$input_args$data_args_input[["Q2config_cp"]], "Map" = vast_build_adjust$tmb_list$Map, "newtonsteps" = 1, "getsd" = TRUE, "getReportCovariance" = TRUE, "run_model" = run_final_model, "test_fit" = FALSE, "Use_REML" = FALSE, "getJointPrecision" = vast_build_adjust$input_args$extra_args$getJointPrecision, "DirPath" = spatial_info_dir)

  # Save and return it
  saveRDS(vast_fit_out, file = paste(out_dir, "/", nice_category_names, "_", "fitted_vast.rds", sep = ""))
  return(vast_fit_out)
}

#' @title Predict fitted VAST model
#'
#' @description This function makes predictions from a fitted VAST SDM to new locations using VAST::predict.fit_model. Importantly, to use this feature for new times, at least one location for each time of interest needs to be included during the model fitting process. This dummy observation should have a PredTF value of 1 so that the observation is only used in the predicted probability and NOT estimating the likelihood.
#'
#' @param vast_fitted_sdm = A fitted VAST SDM object, as returned with `vast_fit_sdm`
#' @param nice_category_names = A character string to define species/model run and used in naming the output prediction file.
#' @param predict_variable = Which variable should be predicted, default is density (D_i).
#' @param predict_category = Which category (species/age/size) should be predicted, default is 0.
#' @param predict_vessel = Which sampling category should be predicted, default is 0.
#' @param predict_covariates_df_all = A long data frame with all of the prediction covariates and arranged with lat/long/time columns.
#' @param cov_names = A character string vector with covariate names used in the model fitting.
#' @param coord_names = A character string vector with the names of the longitude and latitude location column of the `predict_covariates_df_all` dataframe, default is c("Lon", "Lat").
#' @param time_col = A character string with the name of the date of the `predict_covariates_df_all` dataframe, default is "DATE".
#' @param out_dir = Output directory to save .rds dataframe with VAST model predictions for each of the lat/long/times in `predict_covariates_df_all`.
#'
#' @return
#'
#' @export

predict_vast <- function(vast_fitted_sdm, nice_category_names, predict_variable = "D_i", predict_category = 0, predict_vessel = 0, predict_covariates_df_all, cov_names, coord_names = c("Lon", "Lat"), time_col = "DATE", out_dir) {

  # For debugging
  if (FALSE) {
    # Targets
    tar_load(vast_fit)
    vast_fit <- readRDS("~/Box/Mills Lab/Projects/sdm_workflow/targets_output/mod_fits/American_lobster_STnoRW_fitted_vast.rds")
    vast_fitted_sdm <- vast_fit
    nice_category_names <- nice_category_names
    nice_category_names <- "American_lobster"
    predict_variable <- "D_i"
    predict_category <- 0
    predict_vessel <- 0
    tar_load(vast_pred_df_post_fit)
    predict_covariates_df_all <- vast_pred_df_post_fit
    predict_covariates_df_all <- readRDS(paste(out_dir, "/VAST_post_fit_pred_df_", summarize, "_", ensemble_stat, ".rds", sep = ""))
    predict_covariates_df_all <- pred_covs_out_final
    coord_names <- c("DECDEG_BEGLON", "DECDEG_BEGLAT")
    cov_names <- c("Depth", "BS_seasonal", "BT_seasonal", "SS_seasonal", "SST_seasonal")

    # Basic example...
    vast_fitted_sdm <- readRDS(here::here("", "results/mod_fits/1011_fitted_vast.rds"))
    nmfs_species_code <- 101
    predict_variable <- "Index_gctl"
    predict_category <- 0
    predict_vessel <- 0
    predict_covariates_df_all <- pred_df
    time_col <- "Year"
    cov_names <- c("Depth", "SST_seasonal", "BT_seasonal")
  }

  # Rename coordinates
  predict_covariates_df_all <- predict_covariates_df_all %>%
    rename_at(vars(all_of(coord_names)), ~ c("Lon", "Lat"))

  # #### Not the biggest fan of this, but for now, building in a work around to resolve some of the memory issues that we were running into by supplying a 0.25 degree grid and trying to predict/project for each season-year from 1980-2100. To overcome this issue, going to try to just make the projections to knots and do the smoothing later.
  # # First, need to get the knot locations
  # knot_locs<- data.frame(vast_fitted_sdm$spatial_list$latlon_g) %>%
  #   st_as_sf(., coords = c("Lon", "Lat"), remove = FALSE) %>%
  #   mutate(., "Pt_Id" = 1:nrow(.))
  #
  # # Nearest knot to each point?
  # pred_sf<- predict_covariates_df_all %>%
  #   st_as_sf(., coords = c("Lon", "Lat"), remove = FALSE)
  #
  # pred_sf<- pred_sf %>%
  #   mutate(., "Nearest_Knot" = st_nearest_feature(., knot_locs))
  #
  # # Average the points...
  # pred_df_knots<- pred_sf %>%
  #   st_drop_geometry()

  if (time_col == "DATE") {
    # pred_df_knots$DATE<- as.Date(pred_df_knots$DATE)
    # group_by_vec<- c({{time_col}}, "Nearest_Knot")
    # pred_df_knots_summ<- pred_df_knots %>%
    #   group_by_at(.vars = group_by_vec) %>%
    #   summarize_at(all_of(c(cov_names, "VAST_YEAR_COV")), mean, na.rm = TRUE) %>%
    #   left_join(., st_drop_geometry(knot_locs), by = c("Nearest_Knot" = "Pt_Id")) %>%
    #   ungroup()
    pred_df_knots_summ <- predict_covariates_df_all
    pred_df_knots_summ$DATE <- as.Date(pred_df_knots_summ$DATE)
    pred_df_knots_summ$Season <- factor(case_when(
      format(pred_df_knots_summ$DATE, "%m-%d") == "12-16" ~ "WINTER",
      format(pred_df_knots_summ$DATE, "%m-%d") == "03-16" ~ "SPRING",
      format(pred_df_knots_summ$DATE, "%m-%d") == "07-16" ~ "SUMMER",
      format(pred_df_knots_summ$DATE, "%m-%d") == "09-16" ~ "FALL"
    ),
    levels = levels(vast_fit$covariate_data$Season)
    )
    pred_df_knots_summ <- pred_df_knots_summ %>%
      rename(., "Year_Cov" = "VAST_YEAR_COV")
    pred_df_knots_summ$Year_Cov <- factor(pred_df_knots_summ$Year_Cov)
    pred_df_knots_summ$Year <- as.numeric(as.character(pred_df_knots_summ$Year_Cov))
  } else {
    # group_by_vec<- c({{time_col}}, "Nearest_Knot")
    # pred_df_knots_summ<- pred_df_knots %>%
    #   group_by_at(.vars = group_by_vec) %>%
    #   summarize_at(all_of(cov_names), mean, na.rm = TRUE) %>%
    #   left_join(., st_drop_geometry(knot_locs), by = c("Nearest_Knot" = "Pt_Id")) %>%
    #   ungroup()
    #
    pred_df_knots_summ <- predict_covariates_df_all
  }

  # Collecting necessary bits from the prediction covariates -- lat, lon, time
  pred_lats <- pred_df_knots_summ$Lat
  pred_lons <- pred_df_knots_summ$Lon
  if (time_col == "DATE") {
    # Need the "VAST" time..
    # pred_times<- pred_df_knots$VAST_YEAR_SEASON[match(pred_df_knots_summ$Nearest_Knot, pred_df_knots$Nearest_Knot)]
    pred_times <- pred_df_knots_summ$VAST_YEAR_SEASON
  } else {
    # pred_times<- as.numeric(unlist(pred_df_knots_summ[{{time_col}}]))
    pred_times <- as.numeric(unlist(pred_df_knots_summ[{{ time_col }}]))
  }

  # Catch stuff...
  pred_sampled_areas <- rep(1, length(pred_lats))
  pred_category <- rep(predict_category, length(pred_lats))
  pred_vessel <- rep(predict_vessel, length(pred_lats))

  # Renaming predict_covariates_df_all to match vast_fit_covariate_data
  cov_dat_name_order <- match(names(vast_fitted_sdm$covariate_data), names(pred_df_knots_summ))
  cov_dat_name_order <- cov_dat_name_order[!is.na(cov_dat_name_order)]
  pred_cov_dat_use <- pred_df_knots_summ[, cov_dat_name_order]

  # Catchability data?
  if (!is.null(vast_fitted_sdm$catchability_data)) {
    pred_catch_dat_use <- pred_cov_dat_use %>%
      dplyr::select(., c(Year, Year_Cov, Season, Lat, Lon))
    pred_catch_dat_use$Survey <- rep("NMFS", nrow(pred_catch_dat_use))
    pred_catch_dat_use$Survey <- factor(pred_catch_dat_use$Survey, levels = c("NMFS", "DFO"))
  } else {
    pred_catch_dat_use <- NULL
  }

  # Make the predictions
  preds_out <- predict.fit_model_aja(x = vast_fitted_sdm, what = predict_variable, Lat_i = pred_lats, Lon_i = pred_lons, t_i = pred_times, a_i = pred_sampled_areas, c_iz = pred_category, NULL, new_covariate_data = pred_cov_dat_use, new_catchability_data = pred_catch_dat_use, do_checks = FALSE)

  # Get everything as a dataframe to make plotting easier...
  if (time_col == "DATE") {
    pred_df_out <- data.frame("Lat" = pred_lats, "Lon" = pred_lons, "Time" = as.Date(predict_covariates_df_all[, {{ time_col }}]), "Pred" = preds_out)
  } else {
    pred_df_out <- data.frame("Lat" = pred_lats, "Lon" = pred_lons, "Time" = predict_covariates_df_all[, {{ time_col }}], "Pred" = preds_out)
  }


  # Save and return it
  saveRDS(pred_df_out, file = paste(out_dir, "/pred_", predict_variable, "_", nice_category_names, ".rds", sep = ""))
  return(pred_df_out)
}

#' @title Plot prediction spatial summaries
#'
#' @description This function plots the total biomass from SDM predictions within spatial area of interest across time.
#'
#' @param pred_df = A dataframe with Lat, Lon, Time and Pred columns.
#' @param spatial_areas = An sf polygon or multipolygon object with rows for each of the spatial areas to calculate average biomass within.
#' @param plot_regions = A character string vector with the names of spatial areas to plot. Needs to be a matching set of areas in the `spatial_areas` object.
#' @param index_scale = A character string of either "raw" or "log" defining the biomass value scale to calculate/plot. In general, log scale is going to be easier to see patterns across different areas.
#' @param year_stop = A year or NULL to plot all years.
#' @param nice_category_names = A character string to define species/model run and used in naming the output prediction file.
#' @param pred_label = A character string to add to the output file plot identifying which prediction data were used to make predictions and calculate total biomass values.
#' @param nice_xlab = A character string for x axis of the plot, default is "Year-Season".
#' @param nice_ylab = A character string for y axis of the plot, default is "Biomass index (metrix tons)"
#' @param paneling = A place holder for eventually creating "panel" plots, default is "none" and has all of the time series plotted on on plot.
#' @param color_pal = A color pallete character string to use for each of the plot_region time series.
#' @param out_dir = Output directory to save biomass time series plot.
#'
#' @return A ggplot object with biomass time series for each of the plot_regions.
#'
#' @export

pred_plot_spatial_summary <- function(pred_df, spatial_areas, plot_regions = c("NMFS", "DFO", "NMFS_and_DFO"), index_scale = "raw", year_stop = NULL, nice_category_names = nice_category_names, pred_label, nice_xlab = "Year-Season", nice_ylab = "Biomass index (metric tons)", paneling = "none", color_pal = NULL, out_dir = paste0(res_root, "plots_maps")) {
  if (FALSE) {
    nice_category_names <- "Atlantic_cod"
    pred_label <- "SSP5_85_mean"
    pred_df <- read.csv(paste(out_dir, "/", nice_category_names, "_", pred_label, "_projections.csv", sep = ""))
    tar_load(index_shapefiles)
    spatial_areas <- index_shapefiles
    plot_regions <- c("NMFS", "DFO", "NMFS_and_DFO")
    index_scale <- "log"
    year_stop <- NULL
    nice_category_names <- nice_category_names
    nice_xlab <- "Year-Season"
    nice_ylab <- "Log biomass index (metric tons)"
    paneling <- "none"
    color_pal <- NULL
    out_dir <- paste0(res_root, "plots_maps")
  }

  # Convert prediction df to spatial dataframe...
  pred_sp <- st_as_sf(pred_df, coords = c("Lon", "Lat"), crs = 4326, remove = FALSE)

  # Join up with the spatial areas
  pred_sp <- pred_sp %>%
    st_join(., spatial_areas) %>%
    st_drop_geometry()

  # Summaries
  index_res_df <- pred_sp %>%
    mutate(., "Year" = format(as.Date(Real_Date), "%Y")) %>%
    group_by(., Year, Region) %>%
    summarize_at(., vars(Dens), c(sum), na.rm = TRUE) %>%
    mutate(., "Date" = as.Date(paste(Year, "06-15", sep = "-")))

  index_res_df_plot <- index_res_df %>%
    filter(., Region %in% plot_regions)

  if (paneling == "none") {
    if (!is.null(color_pal)) {
      colors_use <- color_pal
    } else {
      color_pal <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854")
      colors_use <- color_pal[1:length(plot_regions)]
    }

    # Filter based on years to plot
    if (!is.null(year_stop)) {
      index_res_df_plot <- index_res_df_plot %>%
        filter(., Year < year_stop)
    }

    # Date axis
    date_breaks <- seq.Date(from = as.Date(paste(min(index_res_df_plot$Year), "06-15", sep = "-")), to = as.Date(paste(max(index_res_df_plot$Year), "06-15", sep = "-")), by = "5 years")
    plot_out <- ggplot() +
      # geom_errorbar(data = index_res_df, aes(x = Date, ymin = (Index_Estimate - Index_SD), ymax = (Index_Estimate + Index_SD), color = Index_Region, group = Index_Region), alpha = 0.65) +
      geom_point(data = index_res_df_plot, aes(x = Date, y = Dens, color = Region)) +
      geom_line(data = index_res_df_plot, aes(x = Date, y = Dens, color = Region)) +
      scale_color_manual(values = colors_use) +
      scale_x_date(breaks = date_breaks, date_labels = "%Y") +
      xlab({{ nice_xlab }}) +
      ylab({{ nice_ylab }}) +
      ggtitle({{ nice_category_names }}) +
      theme_bw() +
      theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  }

  # Save and return the plot
  ggsave(plot_out, file = paste(out_dir, "/Post_Fit_Biomass_Index_", index_scale, "_", nice_category_names, "_", pred_label, ".jpg", sep = ""))
  write.csv(index_res_df, file = paste(gsub("plots_maps", "tables", out_dir), "/Post_Fit_Biomass_Index_", index_scale, "_", nice_category_names, "_", pred_label, ".csv", sep = ""))
  return(plot_out)
}

#' @title Plot VAST model predicted density surfaces
#'
#' @description Creates either a panel plot or a gif of VAST model predicted density surfaces
#'
#' @param vast_fit = A VAST `fit_model` object.
#' @param nice_category_names = A character string to define species/model run and used in naming the output prediction file.
#' @param all_times = A vector of all of the unique time steps available from the VAST fitted model
#' @param plot_times = Either NULL to make a plot for each time in `all_times` or a vector of all of the times to plot, which must be a subset of `all_times`
#' @param land_sf = Land sf object
#' @param xlim = A two element vector with the min and max longitudes
#' @param ylim = A two element vector with the min and max latitudes
#' @param panel_or_gif = A character string of either "panel" or "gif" indicating how the multiple plots across time steps should be displayed
#' @param out_dir = Output directory to save the panel plot or gif
#'
#' @return A VAST fit_model object, with the inputs and and outputs, including parameter estimates, extrapolation gid info, spatial list info, data info, and TMB info.
#'
#' @export

vast_fit_plot_density <- function(vast_fit, nice_category_names, mask, all_times = all_times, plot_times = NULL, land_sf, xlim, ylim, panel_or_gif = "gif", out_dir, land_color = "#d9d9d9", panel_cols = NULL, panel_rows = NULL, ...) {
  if (FALSE) {
    tar_load(vast_fit)
    template <- raster("~/GitHub/sdm_workflow/scratch/aja/TargetsSDM/data/supporting/HighResTemplate.grd")
    tar_load(vast_seasonal_data)
    all_times <- as.character(levels(vast_seasonal_data$VAST_YEAR_SEASON))
    plot_times <- NULL
    tar_load(land_sf)
    tar_load(region_shapefile)
    mask <- region_shapefile
    land_color <- "#d9d9d9"
    res_data_path <- "~/Box/RES_Data/"
    xlim <- c(-85, -55)
    ylim <- c(30, 50)
    panel_or_gif <- "gif"
    panel_cols <- NULL
    panel_rows <- NULL
  }

  # Plotting at spatial knots...
  # Getting prediction array
  pred_array <- log(vast_fit$Report$D_gct + 1)

  # Getting time info
  if (!is.null(plot_times)) {
    plot_times <- all_times[which(all_times) %in% plot_times]
  } else {
    plot_times <- all_times
  }

  # Getting spatial information
  spat_data <- vast_fit$extrapolation_list
  loc_g <- spat_data$Data_Extrap[which(spat_data$Data_Extrap[, "Include"] > 0), c("Lon", "Lat")]
  CRS_orig <- sp::CRS("+proj=longlat")
  CRS_proj <- sp::CRS(spat_data$projargs)
  land_sf <- st_crop(land_sf, xmin = xlim[1], ymin = ylim[1], xmax = xlim[2], ymax = ylim[2])

  # Looping through...
  rasts_out <- vector("list", dim(pred_array)[3])
  rasts_range <- pred_array
  rast_lims <- c(0, round(max(rasts_range) + 0.0000001, 2))

  if (dim(pred_array)[3] == 1) {
    data_df <- data.frame(loc_g, z = pred_array[, 1, ])

    # Interpolation
    pred_df <- na.omit(data.frame("x" = data_df$Lon, "y" = data_df$Lat, "layer" = data_df$z))
    pred_df_interp <- interp(pred_df[, 1], pred_df[, 2], pred_df[, 3],
      duplicate = "mean", extrap = TRUE,
      xo = seq(-87.99457, -57.4307, length = 115),
      yo = seq(22.27352, 48.11657, length = 133)
    )
    pred_df_interp_final <- data.frame(expand.grid(x = pred_df_interp$x, y = pred_df_interp$y), z = c(round(pred_df_interp$z, 2)))
    pred_sp <- st_as_sf(pred_df_interp_final, coords = c("x", "y"), crs = CRS_orig)

    pred_df_temp <- pred_sp[which(st_intersects(pred_sp, mask, sparse = FALSE) == TRUE), ]
    coords_keep <- as.data.frame(st_coordinates(pred_df_temp))
    row.names(coords_keep) <- NULL
    pred_df_use <- data.frame(cbind(coords_keep, "z" = as.numeric(pred_df_temp$z)))
    names(pred_df_use) <- c("x", "y", "z")

    # raster_proj<- raster::rasterize(as_Spatial(points_ll), template, field = "z", fun = mean)
    # raster_proj<- as.data.frame(raster_proj, xy = TRUE)
    #
    time_plot_use <- plot_times

    plot_out <- ggplot() +
      geom_tile(data = pred_df_use, aes(x = x, y = y, fill = z)) +
      scale_fill_viridis_c(name = "Log (density+1)", option = "viridis", na.value = "transparent", limits = rast_lims) +
      annotate("text", x = -65, y = 37.5, label = time_plot_use) +
      geom_sf(data = land_sf, fill = land_color, lwd = 0.2, na.rm = TRUE) +
      coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
      theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "pt"))

    ggsave(filename = paste(out_dir, nice_category_names, ".png", sep = "/"), plot_out, width = 11, height = 8, units = "in")
  } else {
    for (tI in 1:dim(pred_array)[3]) {
      data_df <- data.frame(loc_g, z = pred_array[, 1, tI])

      # Interpolation
      pred_df <- na.omit(data.frame("x" = data_df$Lon, "y" = data_df$Lat, "layer" = data_df$z))
      pred_df_interp <- interp(pred_df[, 1], pred_df[, 2], pred_df[, 3],
        duplicate = "mean", extrap = TRUE,
        xo = seq(-87.99457, -57.4307, length = 115),
        yo = seq(22.27352, 48.11657, length = 133)
      )
      pred_df_interp_final <- data.frame(expand.grid(x = pred_df_interp$x, y = pred_df_interp$y), z = c(round(pred_df_interp$z, 2)))
      pred_sp <- st_as_sf(pred_df_interp_final, coords = c("x", "y"), crs = CRS_orig)

      pred_df_temp <- pred_sp[which(st_intersects(pred_sp, mask, sparse = FALSE) == TRUE), ]
      coords_keep <- as.data.frame(st_coordinates(pred_df_temp))
      row.names(coords_keep) <- NULL
      pred_df_use <- data.frame(cbind(coords_keep, "z" = as.numeric(pred_df_temp$z)))
      names(pred_df_use) <- c("x", "y", "z")

      # raster_proj<- raster::rasterize(as_Spatial(points_ll), template, field = "z", fun = mean)
      # raster_proj<- as.data.frame(raster_proj, xy = TRUE)
      #
      time_plot_use <- plot_times[tI]

      rasts_out[[tI]] <- ggplot() +
        geom_tile(data = pred_df_use, aes(x = x, y = y, fill = z)) +
        scale_fill_viridis_c(name = "Log (density+1)", option = "viridis", na.value = "transparent", limits = rast_lims) +
        annotate("text", x = -65, y = 37.5, label = time_plot_use) +
        geom_sf(data = land_sf, fill = land_color, lwd = 0.2, na.rm = TRUE) +
        coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
        theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "pt"))
    }
    if (panel_or_gif == "panel") {
      # Panel plot
      all_plot <- wrap_plots(rasts_out, ncol = panel_cols, nrow = panel_rows, guides = "collect", theme(plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "pt")))
      ggsave(filename = paste0(out_dir, "/", nice_category_names, "_LogDensity.png"), all_plot, width = 11, height = 8, units = "in")
      return(all_plot)
    } else {
      # Make a gif
      plot_loop_func <- function(plot_list) {
        for (i in seq_along(plot_list)) {
          plot_use <- plot_list[[i]]
          print(plot_use)
        }
      }
      invisible(save_gif(plot_loop_func(rasts_out), paste0(out_dir, "/", nice_category_names, "_LogDensity.gif"), delay = 0.75, progress = FALSE))
    }
  }
}

#' @title Plot predicted density surfaces from data frame
#'
#' @description Creates either a panel plot or a gif of predicted density surfaces from a data frame that has location and time information
#'
#' @param pred_df = A dataframe with Lat, Lon, Time and Pred columns
#' @param nice_category_names = A character string to define species/model run and used in naming the output prediction file.
#' @param mask = Land mask
#' @param plot_times = Either NULL to make a plot for each time in `pred_df$Time` or a vector of all of the times to plot, which must be a subset of `pred_df$Time`
#' @param land_sf = Land sf object
#' @param xlim = A two element vector with the min and max longitudes
#' @param ylim = A two element vector with the min and max latitudes
#' @param panel_or_gif = A character string of either "panel" or "gif" indicating how the multiple plots across time steps should be displayed
#' @param out_dir = Output directory to save the panel plot or gif
#'
#' @return NULL. Panel or gif plot is saved in out_dir.
#'
#' @export

vast_df_plot_density <- function(pred_df, nice_category_names, mask, all_times = all_times, plot_times = NULL, land_sf, xlim, ylim, panel_or_gif = "gif", out_dir, land_color = "#d9d9d9", panel_cols = NULL, panel_rows = NULL, ...) {
  if (FALSE) {
    tar_load(vast_predictions)
    pred_df <- vast_predictions
    plot_times <- NULL
    tar_load(land_sf)
    tar_load(region_shapefile)
    mask <- region_shapefile
    land_color <- "#d9d9d9"
    res_data_path <- "~/Box/RES_Data/"
    xlim <- c(-80, -55)
    ylim <- c(35, 50)
    panel_or_gif <- "gif"
    panel_cols <- NULL
    panel_rows <- NULL
  }

  # Time ID column for filtering
  pred_df <- pred_df %>%
    mutate(., "Time_Filter" = as.numeric(Time))

  # Log transform pred_df$Pred
  pred_df$Pred <- log(pred_df$Pred + 1)

  # Getting all unique times
  all_times <- unique(pred_df$Time)

  # Getting time info
  if (!is.null(plot_times)) {
    plot_times <- all_times[which(all_times) %in% plot_times]
  } else {
    plot_times <- all_times
  }

  # Getting spatial information
  land_sf <- st_crop(land_sf, xmin = xlim[1], ymin = ylim[1], xmax = xlim[2], ymax = ylim[2])

  # Looping through...
  rasts_out <- vector("list", length(plot_times))
  rasts_range <- pred_df$Pred
  rast_lims <- c(0, round(max(rasts_range) + 0.0000001, 2))

  for (tI in 1:length(plot_times)) {
    pred_df_temp <- pred_df %>%
      dplyr::filter(., Time_Filter == tI)

    # Interpolation
    pred_df_temp <- na.omit(data.frame("x" = pred_df_temp$Lon, "y" = pred_df_temp$Lat, "layer" = pred_df_temp$Pred))
    pred_df_interp <- interp(pred_df_temp[, 1], pred_df_temp[, 2], pred_df_temp[, 3],
      duplicate = "mean", extrap = TRUE,
      xo = seq(-87.99457, -57.4307, length = 115),
      yo = seq(22.27352, 48.11657, length = 133)
    )
    pred_df_interp_final <- data.frame(expand.grid(x = pred_df_interp$x, y = pred_df_interp$y), z = c(round(pred_df_interp$z, 2)))
    pred_sp <- st_as_sf(pred_df_interp_final, coords = c("x", "y"), crs = 4326)

    pred_df_temp2 <- pred_sp[which(st_intersects(pred_sp, mask, sparse = FALSE) == TRUE), ]
    coords_keep <- as.data.frame(st_coordinates(pred_df_temp2))
    row.names(coords_keep) <- NULL
    pred_df_use <- data.frame(cbind(coords_keep, "z" = as.numeric(pred_df_temp2$z)))
    names(pred_df_use) <- c("x", "y", "z")

    # raster_proj<- raster::rasterize(as_Spatial(points_ll), template, field = "z", fun = mean)
    # raster_proj<- as.data.frame(raster_proj, xy = TRUE)
    #
    time_plot_use <- plot_times[tI]

    rasts_out[[tI]] <- ggplot() +
      geom_tile(data = pred_df_use, aes(x = x, y = y, fill = z)) +
      scale_fill_viridis_c(name = "Log (density+1)", option = "viridis", na.value = "transparent", limits = rast_lims) +
      annotate("text", x = -65, y = 37.5, label = time_plot_use) +
      geom_sf(data = land_sf, fill = land_color, lwd = 0.2, na.rm = TRUE) +
      coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
      theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "pt"))
  }

  if (panel_or_gif == "panel") {
    # Panel plot
    all_plot <- wrap_plots(rasts_out, ncol = panel_cols, nrow = panel_rows, guides = "collect", theme(plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "pt")))
    ggsave(filename = paste0(out_dir, "/", nice_category_names, "_LogDensity.png", sep = ""), all.plot, width = 11, height = 8, units = "in")
  } else {
    # Make a gif
    plot_loop_func <- function(plot_list) {
      for (i in seq_along(plot_list)) {
        plot_use <- plot_list[[i]]
        print(plot_use)
      }
    }
    invisible(save_gif(plot_loop_func(rasts_out), paste0(out_dir, "/", nice_category_names, "_LogDensity.gif"), delay = 0.75, progress = FALSE))
  }
}

vast_post_fit_pred_df <- function(predict_covariates_stack_agg_dir, extra_covariates_stack, covs_rescale = c("Depth", "BS_seasonal", "BT_seasonal", "SS_seasonal", "SST_seasonal"), rescale_params, depth_cut, mask, summarize, ensemble_stat, vast_fit, fit_seasons, pred_year_min = fit_year_min, fit_year_max = pred_years, pred_year_max, out_dir) {

  # For debugging
  if (FALSE) {
    tar_load(vast_fit)
    tar_load(predict_covariates_stack_agg_out)
    predict_covariates_stack_agg_dir <- predict_covariates_stack_agg_out
    tar_load(static_covariates_stack)
    extra_covariates_stack <- static_covariates_stack
    tar_load(rescale_params)
    tar_load(region_shapefile)
    mask <- region_shapefile
    summarize <- "seasonal"
    ensemble_stat <- "mean"
    pred_year_min <- fit_year_min
    fit_year_max <- pred_years
    pred_year_max <- 2100
    fit_seasons <- fit_seasons
    out_dir <- here::here("data/predict")
    covs_rescale <- c("Depth", "BS_seasonal", "BT_seasonal", "SS_seasonal", "SST_seasonal")
  }

  # Get raster stack covariate files
  rast_files_load <- list.files(predict_covariates_stack_agg_dir, pattern = paste0(summarize, "_", ensemble_stat, ".grd$"), full.names = TRUE)

  # Get variable names
  cov_names_full <- list.files(predict_covariates_stack_agg_dir, pattern = paste0(summarize, "_", ensemble_stat, ".grd$"), full.names = FALSE)
  predict_covs_names <- gsub(paste("_", ensemble_stat, ".grd$", sep = ""), "", gsub("predict_stack_", "", cov_names_full))

  # Get the grid locations...
  #####
  ## Start here -- extrap grid in the vast_fit is going to have grid locations duplicated for the different regions. We really just want to get a single estimate for each location and each time. Somehow, need to match the spatial subset up WITH vast_fit$spatial_list$A_gs[i] and [j], as that is going to have the information we need for omega and epsilon....would an ID column on the Data_Extrap match up with the [i] index of the grids??
  #####
  dat_extrap <- vast_fit$extrapolation_list$Data_Extrap %>%
    mutate(., "Extrap_ID" = paste0("GRID_", seq(from = 1, to = nrow(.))))
  extrap_unique_ids <- dat_extrap[which(dat_extrap$STRATA == "All"), ]

  # Coordinates for these unique IDs
  pred_coords <- data.frame("Lon" = extrap_unique_ids$Lon, "Lat" = extrap_unique_ids$Lat, "Extrap_ID" = extrap_unique_ids$Extrap_ID)
  preds_sf <- pred_coords %>%
    st_as_sf(., coords = c("Lon", "Lat"), crs = 4326, remove = FALSE)

  # Looping through prediction stack time steps
  for (i in 1:nlayers(raster::stack(rast_files_load[1]))) {
    # Get the time index
    time_ind <- i

    # Load corresponding raster layers matching the time index
    pred_covs_stack_temp <- rotate(raster::stack(raster::stack(rast_files_load[1])[[time_ind]], raster::stack(rast_files_load[2])[[time_ind]], raster::stack(rast_files_load[3])[[time_ind]], raster::stack(rast_files_load[4])[[time_ind]]))

    # Mask out values outside area of interest
    pred_covs_stack_temp <- raster::mask(pred_covs_stack_temp, mask = mask)

    # Extract values at the knot locations
    pred_covs_df_temp <- data.frame("DECDEG_BEGLON" = pred_coords$Lon, "DECDEG_BEGLAT" = pred_coords$Lat, "Extrap_ID" = pred_coords$Extrap_ID, raster::extract(pred_covs_stack_temp, preds_sf, method = "bilinear")) %>%
      drop_na()

    # Some processing to keep observations within our area of interest and get things in a "tidy-er" prediction dataframe
    time_name <- sub(".[^.]*$", "", names(pred_covs_stack_temp))
    colnames(pred_covs_df_temp)[4:ncol(pred_covs_df_temp)] <- paste(time_name, predict_covs_names, sep = "_")
    colnames(pred_covs_df_temp)[4:ncol(pred_covs_df_temp)] <- gsub("X", "", gsub("[.]", "_", colnames(pred_covs_df_temp)[4:ncol(pred_covs_df_temp)]))

    pred_covs_df_out_temp <- pred_covs_df_temp %>%
      pivot_longer(., -c(DECDEG_BEGLON, DECDEG_BEGLAT, Extrap_ID), names_to = c("variable"), values_to = "value") %>%
      separate(., variable, into = c("EST_YEAR", "SEASON", "variable"), sep = "_", extra = "merge") %>%
      pivot_wider(., names_from = variable, values_from = value)

    # Adding in some other columns we will want to match up easily with 'vast_data_out'
    pred_covs_df_out_temp <- pred_covs_df_out_temp %>%
      mutate(.,
        EST_YEAR = as.numeric(EST_YEAR),
        DATE = paste(EST_YEAR, case_when(
          SEASON == "Winter" ~ "12-16",
          SEASON == "Spring" ~ "03-16",
          SEASON == "Summer" ~ "07-16",
          SEASON == "Fall" ~ "09-16"
        ), sep = "-"),
        SURVEY = "NMFS",
        SVVESSEL = "DUMMY",
        NMFS_SVSPP = "DUMMY",
        DFO_SPEC = "DUMMY",
        PRESENCE = 1,
        BIOMASS = 1,
        ABUNDANCE = 1,
        ID = paste("DUMMY", DATE, sep = ""),
        PredTF = TRUE
      )

    if (i == 1) {
      pred_covs_out <- pred_covs_df_out_temp
    } else {
      pred_covs_out <- bind_rows(pred_covs_out, pred_covs_df_out_temp)
    }
  }

  # Only going to keep information from fit_year_max through pred_years...
  pred_covs_out_final <- pred_covs_out %>%
    dplyr::filter(., EST_YEAR >= pred_year_min & EST_YEAR <= pred_year_max)

  # New implementation...
  pred_covs_out_final <- pred_covs_out_final %>%
    mutate(., # VAST_YEAR_COV = EST_YEAR,
      # VAST_YEAR_COV = ifelse(EST_YEAR > fit_year_max, fit_year_max, EST_YEAR),
      # VAST_YEAR_COV = ifelse(EST_YEAR > pred_year_max, pred_year_max, EST_YEAR),
      VAST_YEAR_COV = EST_YEAR,
      VAST_SEASON = case_when(
        SEASON == "Spring" ~ "SPRING",
        SEASON == "Summer" ~ "SUMMER",
        SEASON == "Fall" ~ "FALL"
      ),
      "VAST_YEAR_SEASON" = paste(VAST_YEAR_COV, VAST_SEASON, sep = "_")
    )

  # Subset to only seasons of interest...
  pred_covs_out_final <- pred_covs_out_final %>%
    filter(., VAST_SEASON %in% fit_seasons)

  # Need to account for new levels in year season...
  # all_years<- seq(from = pred_year_min, to = max(fit_year_max), by = 1)
  all_years <- seq(from = pred_year_min, to = max(pred_year_max), by = 1)
  all_seasons <- fit_seasons
  year_season_set <- expand.grid("SEASON" = all_seasons, "EST_YEAR" = all_years)
  all_year_season_levels <- apply(year_season_set[, 2:1], MARGIN = 1, FUN = paste, collapse = "_")

  pred_covs_out_final <- pred_covs_out_final %>%
    mutate(.,
      "VAST_YEAR_SEASON" = factor(VAST_YEAR_SEASON, levels = all_year_season_levels),
      "VAST_SEASON" = factor(VAST_SEASON, levels = all_seasons)
    )

  # Name rearrangement!
  # Keep only what we need..
  cov_names <- names(pred_covs_out_final)[-which(names(pred_covs_out_final) %in% c("ID", "DATE", "EST_YEAR", "SEASON", "SURVEY", "SVVESSEL", "DECDEG_BEGLAT", "DECDEG_BEGLON", "Extrap_ID", "NMFS_SVSPP", "DFO_SPEC", "PRESENCE", "BIOMASS", "ABUNDANCE", "PredTF", "VAST_YEAR_COV", "VAST_SEASON", "VAST_YEAR_SEASON"))]
  pred_covs_out_final <- pred_covs_out_final %>%
    dplyr::select(., "ID", "DATE", "EST_YEAR", "SEASON", "SURVEY", "SVVESSEL", "DECDEG_BEGLAT", "DECDEG_BEGLON", "Extrap_ID", "NMFS_SVSPP", "DFO_SPEC", "PRESENCE", "BIOMASS", "ABUNDANCE", "PredTF", "VAST_YEAR_COV", "VAST_SEASON", "VAST_YEAR_SEASON", {{ cov_names }})

  # Any extra covariates will likely be static...
  if (!is.null(extra_covariates_stack)) {
    pred_covs_sf <- points_to_sf(pred_covs_out_final)

    pred_covs_out_final <- static_extract_wrapper(static_covariates_list = extra_covariates_stack, sf_points = pred_covs_sf, date_col_name = "DATE", df_sf = "df", out_dir = NULL)
  }

  # Apply depth cut and drop NAs
  pred_covs_out_final <- pred_covs_out_final %>%
    mutate(.,
      "Depth" = ifelse(Depth > depth_cut, NA, Depth),
      "Summarized" = summarize,
      "Ensemble_Stat" = ensemble_stat
    ) %>%
    drop_na()

  # Rescale
  if (!is.null(rescale_params)) {
    for (i in seq_along(covs_rescale)) {
      match_mean <- rescale_params[which(names(rescale_params) == paste(covs_rescale[i], "Mean", sep = "_"))]
      match_sd <- rescale_params[which(names(rescale_params) == paste(covs_rescale[i], "SD", sep = "_"))]
      pred_covs_out_final <- pred_covs_out_final %>%
        mutate_at(., {{ covs_rescale[i] }}, .funs = covariate_rescale_func, type = "AJA", center = match_mean, scale = match_sd)
    }
  }

  # Some formatting...
  pred_covs_out_final <- data.frame(
    "Year" = as.numeric(pred_covs_out_final$VAST_YEAR_SEASON) - 1,
    "Year_Cov" = factor(pred_covs_out_final$VAST_YEAR_COV),
    "Season" = pred_covs_out_final$VAST_SEASON,
    "Depth" = pred_covs_out_final$Depth,
    "SST_seasonal" = pred_covs_out_final$SST_seasonal,
    "BT_seasonal" = pred_covs_out_final$BT_seasonal,
    "BS_seasonal" = pred_covs_out_final$BS_seasonal,
    "SS_seasonal" = pred_covs_out_final$SS_seasonal,
    "Lat" = pred_covs_out_final$DECDEG_BEGLAT,
    "Lon" = pred_covs_out_final$DECDEG_BEGLON,
    "Extrap_ID" = pred_covs_out_final$Extrap_ID,
    "Date" = pred_covs_out_final$DATE
  )

  # Save and return it
  saveRDS(pred_covs_out_final, file = paste(out_dir, "/VAST_post_fit_pred_df_", summarize, "_", ensemble_stat, ".rds", sep = ""))
  return(pred_covs_out_final)
}

make_vast_proj_objects <- function(vast_fit = vast_fit, time_covs, pred_covs_out_final = vast_pred_df_post_fit, gam_degree = gam_degree, hab_env_coeffs_n = hab_env_coeffs_n) {
  if (FALSE) {
    tar_load(vast_fit)
    vast_fit <- vast_fit
    time_covs <- "Season"
    tar_load(vast_pred_df_post_fit)
    pred_covs_out_final <- vast_pred_df_post_fit
    gam_degree <- gam_degree
    hab_env_coeffs_n <- hab_env_coeffs_n
  }
  proj_obs_out <- vector("list", length = 4)
  names(proj_obs_out) <- c("proj_X_contrasts", "proj_X1_config", "proj_X2_config", "proj_map")

  if (length(time_covs) == 2) {
    proj_obs_out[["proj_X_contrasts"]] <- list(Season = contrasts(pred_covs_out_final$Season, contrasts = FALSE), Year_Cov = contrasts(pred_covs_out_final$Year_Cov, contrasts = FALSE))

    ## Move to a new function
    proj_obs_out[["proj_X1_config"]] <- matrix(c(2, rep(3, length(unique(pred_covs_out_final$Season)) - 1), 2, rep(3, nlevels(pred_covs_out_final$Year_Cov) - 1), rep(rep(1, gam_degree), hab_env_coeffs_n)), nrow = 1)
    proj_obs_out[["proj_X2_config"]] <- matrix(c(2, rep(3, length(unique(pred_covs_out_final$Season)) - 1), 2, rep(3, nlevels(pred_covs_out_final$Year_Cov) - 1), rep(rep(1, gam_degree), hab_env_coeffs_n)), nrow = 1)

    proj_map <- vast_fit$tmb_list$Map
    proj_map$gamma1_cp <- factor(c(proj_map$gamma1_cp, seq(from = (1 + max(as.numeric(as.character(levels(proj_map$gamma1_cp))))), to = max(as.numeric(as.character(levels(proj_map$gamma1_cp)))) + length(levels(pred_covs_out_final$Year_Cov)[!levels(pred_covs_out_final$Year_Cov) %in% levels(vast_fit$covariate_data$Year_Cov)]))))
    proj_map$gamma2_cp <- factor(c(proj_map$gamma2_cp, seq(from = (1 + max(as.numeric(as.character(levels(proj_map$gamma2_cp))))), to = max(as.numeric(as.character(levels(proj_map$gamma2_cp)))) + length(levels(pred_covs_out_final$Year_Cov)[!levels(pred_covs_out_final$Year_Cov) %in% levels(vast_fit$covariate_data$Year_Cov)]))))

    proj_map$log_sigmaXi1_cp <- factor(c(rep(1, length(unique(pred_covs_out_final$Season))), rep(4, nlevels(pred_covs_out_final$Year_Cov)), rep(rep(NA, gam_degree), hab_env_coeffs_n)))
    proj_map$log_sigmaXi2_cp <- factor(c(rep(1, length(unique(pred_covs_out_final$Season))), rep(4, nlevels(pred_covs_out_final$Year_Cov)), rep(rep(NA, gam_degree), hab_env_coeffs_n)))
  } else {
    proj_obs_out[["proj_X_contrasts"]] <- list(Season = contrasts(pred_covs_out_final$Season, contrasts = FALSE))

    ## Move to a new function
    proj_obs_out[["proj_X1_config"]] <- vast_fit$X1config_cp
    proj_obs_out[["proj_X2_config"]] <- vast_fit$X2config_cp

    proj_map <- vast_fit$tmb_list$Map
  }

  proj_map$Xiinput1_scp <- factor(seq(from = 1, to = vast_fit$data_list$n_s * dim(proj_obs_out[["proj_X1_config"]])[2]))
  proj_map$Xiinput2_scp <- factor(seq(from = 1, to = vast_fit$data_list$n_s * dim(proj_obs_out[["proj_X2_config"]])[2]))

  proj_map$beta1_ft <- factor(rep(1, length(unique(pred_covs_out_final$Year))))
  proj_map$beta2_ft <- factor(rep(1, length(unique(pred_covs_out_final$Year))))

  proj_map$Beta_mean1_t <- factor(rep(NA, length(unique(pred_covs_out_final$Year))))
  proj_map$Beta_mean2_t <- factor(rep(NA, length(unique(pred_covs_out_final$Year))))

  proj_map$lagrange_tc <- factor(rep(NA, length(unique(pred_covs_out_final$Year))))

  proj_obs_out[["proj_map"]] <- proj_map

  return(proj_obs_out)
}

# project.fit_model_wrapper <- function(n_sims = n_sims, sim_type = 1, vast_fit = vast_fit, time_cov = "Year_Cov", time_cov_method = time_cov_method, index_shapes = index_shapefiles, historical_uncertainty = "none", n_samples = 1, n_proj = n_proj, new_covariate_data = vast_pred_df_post_fit, new_catchability_data = NULL, proj_objects = proj_objects, seed = 123456, nice_category_names = nice_category_names, out_dir = paste0(res_root, "prediction_df")) {
#   if (FALSE) {
#     n_sims <- 1
#     sim_type <- 1
#     tar_load(vast_fit)
#     time_cov <- NULL
#     time_cov_method <- NULL
#     tar_load(index_shapefiles)
#     historical_uncertainty <- "none"
#     n_samples <- 1
#     n_proj <- 3
#     tar_load(vast_pred_df_post_fit)
#     new_covariate_data <- vast_pred_df_post_fit
#     new_catchability_data <- NULL
#     tar_load(proj_objects)
#     proj_objects <- proj_objects
#     seed <- 123456
#     nice_category_names <- nice_category_names
#     out_dir <- paste0(res_root, "prediction_df")
#   }
#   sim_res <- vector("list", length = n_sims)
#   for (rI in 1:n_sims) {
#     sim_res[[rI]] <- project.fit_model(
#       sim_type = sim_type,
#       x = vast_fit,
#       time_cov = time_cov,
#       time_cov_method = time_cov_method,
#       index_shapes = index_shapes,
#       historical_uncertainty = historical_uncertainty,
#       n_samples = 1,
#       n_proj = n_proj,
#       new_covariate_data = new_covariate_data,
#       new_catchability_data = NULL,
#       proj_X_contrasts = proj_objects[["proj_X_contrasts"]],
#       proj_X1_config = proj_objects[["proj_X1_config"]],
#       proj_X2_config = proj_objects[["proj_X2_config"]],
#       proj_map = proj_objects[["proj_map"]],
#       seed = rI
#     )
#     gc()
#   }
#   saveRDS(sim_res, file = paste0(out_dir, "/", nice_category_names, "_", climate_scenario, "_projections.rds"))
#   return(sim_res)
# }

make_new_cov_data<- function(vast_fit = vast_fit, climate_scenario = climate_scenario){
        # Read in data
        cmip6_dat <- read_csv(paste0(here::here("data/predict/"), climate_scenario, "_dat.csv"))
        # new_cov_dat <- cmip6_dat %>%
        #     dplyr::select(., Year, Year_Cov, Season, Depth, BT_seasonal, Lat, Lon)
        new_cov_dat <- cmip6_dat 
        new_cov_dat$Season <- factor(new_cov_dat$Season, levels = levels(vast_fit$covariate_data$Season))
        
        # Make sure we have the right "year" index
        # match_year_cov <- vast_fit$covariate_data %>%
        #   dplyr::select(., Year, Year_Cov, Season) %>%
        #   distinct() %>%
        #   mutate(., "Year_Cov" = as.numeric(as.character(Year_Cov)))
        match_year_cov <- vast_fit$covariate_data %>%
          dplyr::select(., Year, Year_Cov, Season) %>%
          distinct() %>%
          mutate(., "Year_Cov" = as.numeric(as.character(Year_Cov)))

        new_cov_dat <- new_cov_dat %>%
            left_join(., match_year_cov)

        return(new_cov_dat)
}

make_new_catch_data<- function(vast_fit = vast_fit, climate_scenario = climate_scenario){

  if(FALSE){
    climate_scenario = "CMIP6_SSP5_85_mean"
    tar_load(vast_fit_Offshore_hake_full_69)
    vast_fit = vast_fit_Offshore_hake_full_69
  }

   # Read in data
   cmip6_dat <- read.csv(paste0(here::here("data/predict/"), climate_scenario, "_dat.csv"))
   new_cov_dat <- cmip6_dat %>%
     dplyr::select(., Year, Year_Cov, Season, Depth, BT_seasonal, Lat, Lon)
   new_cov_dat$Season <- factor(new_cov_dat$Season, levels = levels(vast_fit$covariate_data$Season))
    
    # Make sure we have the right "year" index
    # Make sure we have the right "year" index
    match_year_cov <- vast_fit$covariate_data %>%
      dplyr::select(., Year, Year_Cov, Season) %>%
      distinct() %>%
      mutate(., "Year_Cov" = as.numeric(as.character(Year_Cov)))

    
    new_cov_dat <- new_cov_dat %>%
      left_join(., match_year_cov)

     # Explicit catchability data
     new_catch_dat <- new_cov_dat %>%
       dplyr::select(., Year, Year_Cov, Season, Lat, Lon) %>%
       mutate(.,
         Survey = factor("NMFS", levels = c("NMFS", "DFO")))
      new_catch_dat$Year_Cov<- factor(new_catch_dat$Year_Cov, levels = unique(new_cov_dat$Year_Cov))
     
     return(new_catch_dat)

}

make_new_spatial_info<- function(vast_fit = vast_fit, new_cov_data = new_cov_data){
  spatial_args_new <- list("n_x" = vast_fit$spatial_list$n_x, "anisotropic_mesh" = vast_fit$spatial_list$MeshList$anisotropic_mesh, "Kmeans" = vast_fit$spatial_list$Kmeans, "Lon_i" = c(vast_fit$data_frame$Lon_i, new_cov_data$Lon), "Lat_i" = c(vast_fit$data_frame$Lat_i, new_cov_data$Lat), "Extrapolation_List" = vast_fit$extrapolation_list, "grid_size_km" = vast_fit$settings$grid_size_km, "fine_scale" = vast_fit$spatial_list$fine_scale)
  spatial_args_input <- combine_lists(input = spatial_args_new, default = vast_fit$input_args$spatial_args_input)
  proj_spatial <- do.call(what = make_spatial_info, args = spatial_args_input)
  return(proj_spatial)
}

project_model_aja<- function (x, what, n_proj = NULL, n_samples, uncert_res = NULL, new_covariate_data = NULL, 
    new_catchability_data = NULL, extrapolation_list = NULL, 
    input_grid = NULL, spatial_list = NULL, historical_uncertainty = "both", 
    seed = 123456, working_dir = paste0(getwd(), "/"), out_dir = NULL, nice_category_names = NULL, climate_scenario = NULL, ...) 
{
    if (FALSE) {
        x = fit
        what = c("Epsilon1_gct", "Epsilon2_gct", "eta1_gct", 
            "eta2_gct", "P1_gct", "D_gct", "Index_ctl", "effective_area_ctl")
        n_proj = 243 - 12
        n_samples = 500
        new_covariate_data = new_cov_dat
        new_catchability_data = new_catch_dat
        historical_uncertainty = "none"
        seed = rI
        working_dir = date_dir
        extrapolation_list = fit$extrapolation_list
        input_grid = fit$extrapolation_list$Data_Extrap[, c("Lon", 
            "Lat", "Region", "Area_km2", "STRATA")]
        spatial_list = proj_spatial
    }
    reload_model(x)
    Obj = x$tmb_list$Obj
    Sdreport = x$parameter_estimates$SD
    if (is.null(new_covariate_data)) {
        new_covariate_data = x$covariate_data
    } else {
        if (!all(colnames(x$covariate_data) %in% colnames(new_covariate_data))) {
            stop("Please ensure that all columns of `x$covariate_data` are present in `new_covariate_data`")
        }
        n_proj_obs = nrow(new_covariate_data)
        t_proj = new_covariate_data$Year
        lat_proj = new_covariate_data$Lat
        lon_proj = new_covariate_data$Lon
        new_covariate_data = new_covariate_data[, match(colnames(x$covariate_data), 
            colnames(new_covariate_data))]
        NN = RANN::nn2(query = x$covariate_data[, c("Lat", "Lon", 
            "Year")], data = new_covariate_data[, c("Lat", "Lon", 
            "Year")], k = 1)
        if (any(NN$nn.dist == 0)) {
            x$covariate_data = x$covariate_data[-which(NN$nn.dist == 
                0), , drop = FALSE]
        }
        if ("Season" %in% names(new_covariate_data)) {
            new_covariate_data$Season <- factor(new_covariate_data$Season, 
                levels = levels(x$covariate_data$Season))
        }
    }
    if (is.null(new_catchability_data)) {
        new_catchability_data = x$catchability_data
    } else {
        if (!all(colnames(x$catchability_data) %in% colnames(new_catchability_data))) {
            stop("Please ensure that all columns of `x$catchability_data` are present in `new_catchability_data`")
        }
        new_catchability_data = new_catchability_data[, match(
            colnames(x$catchability_data),
            colnames(new_catchability_data)
        )]
        new_catchability_data = rbind(x$catchability_data, new_catchability_data)
    }
    if (any(x$data_list$RhoConfig %in% c(0, 3))) {
    }
    if (any(x$data_list$RhoConfig[c("Beta1", "Beta2")] %in% c(0))) {
        stop("`project_model` is currently designed to work only with temporally varying or constant beta terms")
    }
    rmvnorm_prec <- function(mu, prec, n.sims, seed) {
        set.seed(seed)
        z <- matrix(rnorm(length(mu) * n.sims), ncol = n.sims)
        L <- Matrix::Cholesky(prec, super = TRUE)
        z <- Matrix::solve(L, z, system = "Lt")
        z <- Matrix::solve(L, z, system = "Pt")
        z <- as.matrix(z)
        return(mu + z)
    }
    if (historical_uncertainty == "both") {
        u_zr = rmvnorm_prec(mu = Obj$env$last.par.best, prec = Sdreport$jointPrecision, 
            n.sims = n_samples, seed = seed)
    }
    else if (historical_uncertainty == "random") {
        Obj$retape()
        Obj$fn(x$parameter_estimates$par)
        u_zr = Obj$env$last.par.best %o% rep(1, n_samples)
        set.seed(seed)
        MC = Obj$env$MC(keep = TRUE, n = n_samples, antithetic = FALSE)
        u_zr[Obj$env$random, ] = attr(MC, "samples")
    }
    else if (historical_uncertainty == "none") {
        u_zr = Obj$env$last.par.best %o% rep(1, n_samples)
    }
    else {
        stop("Check `historical_uncertainty` argument")
    }
    # t_i = c(x$data_frame$t_i, max(x$data_frame$t_i) + rep(seq_len(n_proj), 
    #     each = 2))
    # b_i = c(x$data_list$b_i, as_units(rep(c(0, mean(x$data_frame$b_i)), 
    #     n_proj), units(x$data_list$b_i)))
    # v_i = c(x$data_frame$v_i, rep(0, 2 * n_proj))
    # Lon_i = c(x$data_frame$Lon_i, rep(mean(x$data_frame$Lon_i), 
    #     2 * n_proj))
    # Lat_i = c(x$data_frame$Lat_i, rep(mean(x$data_frame$Lat_i), 
    #     2 * n_proj))
    # a_i = c(x$data_list$a_i, as_units(rep(mean(x$data_frame$a_i), 
    #     2 * n_proj), units(x$data_list$a_i)))
    # PredTF_i = c(x$data_list$PredTF_i, rep(1, 2 * n_proj))
    # c_iz = rbind(x$data_list$c_iz, x$data_list$c_iz[rep(seq_len(n_proj), 
    #     each = 2), , drop = FALSE])
    # new_catchability_data = rbind(x$catchability_data, new_catchability_data[1:486, 
    #     ])
    # str(b_i)
    # str(new_catchability_data)
    t_i = c(x$data_frame$t_i, t_proj)
    Lon_i = c(x$data_frame$Lon_i, lon_proj)
    Lat_i = c(x$data_frame$Lat_i, lat_proj)
    b_i = c(x$data_list$b_i, as_units(sample(c(0, 1), size = n_proj_obs, 
        replace = TRUE), units(x$data_list$b_i)))
    v_i = c(x$data_frame$v_i, rep(0, n_proj_obs))
    a_i = c(x$data_list$a_i, as_units(rep(mean(x$data_frame$a_i), 
        n_proj_obs), units(x$data_list$a_i)))
    PredTF_i = c(x$data_list$PredTF_i, rep(1, n_proj_obs))
    c_iz = rbind(matrix(x$data_list$c_iz), matrix(0, nrow = n_proj_obs))
    proj_t = x$data_list$n_t + seq_len(n_proj)
    if ("Season" %in% names(new_covariate_data)) {
        x1 = fit_model(settings = x$settings, Lat_i = Lat_i, 
            Lon_i = Lon_i, t_i = t_i, b_i = b_i, a_i = a_i, v_i = v_i, 
            c_iz = c_iz, PredTF_i = PredTF_i, covariate_data = new_covariate_data, 
            X1_formula = x$X1_formula, X2_formula = x$X2_formula, 
            X_contrasts = list(Season = contrasts(new_covariate_data$Season, 
                contrasts = FALSE)), X1config_cp = x$X1config_cp, 
            X2config_cp = x$X2config_cp, catchability_data = new_catchability_data, 
            Q1config_k = x$Q1config_k, Q2config_k = x$Q2config_k, 
            Q1_formula = x$Q1_formula, Q2_formula = x$Q2_formula, 
            build_model = FALSE, working_dir = working_dir, input_grid = input_grid, extrapolation_list = extrapolation_list, spatial_list = spatial_list)
    }
    else {
        x1 = fit_model(settings = x$settings, Lat_i = Lat_i, 
            Lon_i = Lon_i, t_i = t_i, b_i = b_i, a_i = a_i, v_i = v_i, 
            c_iz = c_iz, PredTF_i = PredTF_i, covariate_data = new_covariate_data, 
            X1_formula = x$X1_formula, X2_formula = x$X2_formula, 
            X1config_cp = x$X1config_cp, X2config_cp = x$X2config_cp, 
            catchability_data = new_catchability_data, Q1config_k = x$Q1config_k, 
            Q2config_k = x$Q2config_k, Q1_formula = x$Q1_formula, 
            Q2_formula = x$Q2_formula, build_model = FALSE, working_dir = working_dir)
    }
    out = vector("list", length = n_samples)
    for (sampleI in seq_len(n_samples)) {
        ParList1 = x1$tmb_list$Parameters
        ParList = Obj$env$parList(par = u_zr[, sampleI])
        for (i in seq_along(ParList)) {
            dim = function(x) {
                if (is.vector(x)) {
                  return(length(x))
                }
                else {
                  return(base::dim(x))
                }
            }
            dim_match = (dim(ParList[[i]]) == dim(ParList1[[i]]))
            if (sum(dim_match == FALSE) == 0) {
                ParList1[[i]] = ParList[[i]]
            }
            else if (sum(dim_match == FALSE) == 1) {
                dim_list = lapply(dim(ParList[[i]]), FUN = function(x) {
                  seq_len(x)
                })
                ParList1[[i]][as.matrix(expand.grid(dim_list))] = ParList[[i]][as.matrix(expand.grid(dim_list))]
            }
            else if (sum(dim_match == FALSE) >= 2) {
                stop("Check matching")
            }
        }
        if (x$data_list$RhoConfig["Beta1"] == 3) {
            tmp = ParList1$beta1_ft
            tmp[, proj_t] = NA
            ParList1$beta1_ft = ifelse(is.na(tmp), rowMeans(tmp, 
                na.rm = TRUE) %o% rep(1, ncol(tmp)), ParList1$beta1_ft)
        }
        if (x$data_list$RhoConfig["Beta2"] == 3) {
            tmp = ParList1$beta2_ft
            tmp[, proj_t] = NA
            ParList1$beta2_ft = ifelse(is.na(tmp), rowMeans(tmp, 
                na.rm = TRUE) %o% rep(1, ncol(tmp)), ParList1$beta2_ft)
        }
        if ("Season" %in% names(new_covariate_data)) {
            x2 = fit_model(settings = x$settings, Lat_i = Lat_i, 
                Lon_i = Lon_i, t_i = t_i, b_i = b_i, a_i = a_i, 
                v_i = v_i, c_iz = c_iz, PredTF_i = PredTF_i, 
                covariate_data = new_covariate_data, X1_formula = x$X1_formula, 
                X2_formula = x$X2_formula, X1config_cp = x$X1config_cp, 
                X2config_cp = x$X2config_cp, X_contrasts = list(Season = contrasts(new_covariate_data$Season, 
                  contrasts = FALSE)), catchability_data = new_catchability_data, 
                Q1config_k = x$Q1config_k, Q2config_k = x$Q2config_k, 
                Q1_formula = x$Q1_formula, Q2_formula = x$Q2_formula, 
                run_model = FALSE, Parameters = ParList1, working_dir = working_dir, 
                input_grid = input_grid, extrapolation_list = extrapolation_list, spatial_list = spatial_list)
        }
        else {
            x2 = fit_model(settings = x$settings, Lat_i = Lat_i, 
                Lon_i = Lon_i, t_i = t_i, b_i = b_i, a_i = a_i, 
                v_i = v_i, c_iz = c_iz, PredTF_i = PredTF_i, 
                covariate_data = new_covariate_data, X1_formula = x$X1_formula, 
                X2_formula = x$X2_formula, X1config_cp = x$X1config_cp, 
                X2config_cp = x$X2config_cp, catchability_data = new_catchability_data, 
                Q1config_k = x$Q1config_k, Q2config_k = x$Q2config_k, 
                Q1_formula = x$Q1_formula, Q2_formula = x$Q2_formula, 
                run_model = FALSE, Parameters = ParList1, working_dir = working_dir)
        }
        x2$tmb_list$Obj$env$data$Options_list$simulate_t[] = c(rep(0, 
            x$data_list$n_t), rep(1, n_proj))
        out[[sampleI]] <- simulate_data(fit = x2, type = 1, random_seed = NULL)
        x2$Report = out[[sampleI]]
        out[[sampleI]] = amend_output(x2)
        if (!is.null(what)) {
            out[[sampleI]] <- out[[sampleI]][which(names(out[[sampleI]]) %in% 
                what)]
        }
        gc()
    }
    if (n_samples == 1) {
      out = out[[1]]
    }
    saveRDS(out, file = paste0(out_dir, "/", nice_category_names, "_", climate_scenario, "_projections.rds"))
    return(out)
}


# predict.fit_model_aja <- function(x, what = "D_i", Lat_i, Lon_i, t_i, a_i, c_iz = rep(0, length(t_i)), v_i = rep(0, length(t_i)), new_covariate_data = NULL, new_catchability_data = NULL, do_checks = TRUE, working_dir = paste0(getwd(), "/")) {
#   if (FALSE) {
#     tar_load(vast_fit)
#     x <- vast_fit
#     what <- "D_i"
#     Lat_i <- x$data_frame$Lat_i
#     # Lat_i = pred_cov_dat_use$Lat
#     Lon_i <- x$data_frame$Lon_i
#     # Lon_i = pred_cov_dat_use$Lon
#     t_i <- x$data_frame$t_i
#     # t_i = pred_cov_dat_use$Year
#     a_i <- x$data_frame$a_i
#     # a_i<- rep(unique(pred_sampled_areas), length(Lat_i))
#     c_iz <- rep(0, length(t_i))
#     # c_iz<- rep(unique(predict_category), length(Lat_i))
#     v_i <- rep(0, length(t_i))
#     # v_i<- rep(unique(predict_vessel), length(t_i))
#     new_covariate_data <- NULL
#     # new_covariate_data = pred_cov_dat_use
#     new_catchability_data <- NULL
#     # new_catchability_data = pred_catch_dat_use
#     do_checks <- FALSE

#     x <- vast_fit
#     what <- "Index_gctl"
#     Lat_i <- predict_covariates_df_all[, "DECDEG_BEGLAT"]
#     Lon_i <- predict_covariates_df_all[, "DECDEG_BEGLON"]
#     t_i <- predict_covariates_df_all[, "t_i"]
#     a_i <- predict_covariates_df_all[, "a_i"]
#     c_iz <- predict_covariates_df_all[, "c_iz"]
#     v_i <- predict_covariates_df_all[, "v_i"]
#     new_covariate_data <- pred_cov_dat_use
#     new_catchability_data <- pred_catch_dat_use
#     do_checks <- FALSE
#     working_dir <- paste0(getwd(), "/")

#     # object = vast_fit
#     # x = object
#     # Lat_i = object$data_frame$Lat_i
#     # Lon_i = object$data_frame$Lon_i
#     # t_i = object$data_frame$t_i
#     # a_i = object$data_frame$a_i
#     # c_iz = rep(0,length(t_i))
#     # v_i = rep(0,length(t_i))
#     # what = "P1_iz"
#     # new_covariate_data = object$covariate_data
#     # new_catchability_data = object$catchability_data
#     # do_checks = FALSE

#     x <- vast_fitted_sdm
#     what <- predict_variable
#     Lat_i <- pred_lats
#     Lon_i <- pred_lons
#     t_i <- pred_times
#     a_i <- pred_sampled_areas
#     c_iz <- pred_category
#     v_i <- rep(0, length(t_i))
#     new_covariate_data <- pred_cov_dat_use
#     new_covariate_data <- pred_covs_out_final
#     new_catchability_data <- pred_catch_dat_use
#     do_checks <- FALSE
#     working_dir <- paste0(getwd(), "/")
#   }

#   message("`predict.fit_model(.)` is in beta-testing, and please explore results carefully prior to using")

#   # Check issues
#   if (!(what %in% names(x$Report)) || (length(x$Report[[what]]) != x$data_list$n_i)) {
#     stop("`what` can only take a few options")
#   }
#   if (!is.null(new_covariate_data)) {
#     # Confirm all columns are available
#     if (!all(colnames(x$covariate_data) %in% colnames(new_covariate_data))) {
#       stop("Please ensure that all columns of `x$covariate_data` are present in `new_covariate_data`")
#     }
#     # Eliminate unnecessary columns
#     new_covariate_data <- new_covariate_data[, match(colnames(x$covariate_data), colnames(new_covariate_data))]
#     # Eliminate old-covariates that are also present in new_covariate_data
#     NN <- RANN::nn2(query = x$covariate_data[, c("Lat", "Lon", "Year")], data = new_covariate_data[, c("Lat", "Lon", "Year")], k = 1)
#     if (any(NN$nn.dist == 0)) {
#       x$covariate_data <- x$covariate_data[-which(NN$nn.dist == 0), , drop = FALSE]
#     }
#   }
#   if (!is.null(new_catchability_data)) {
#     # Confirm all columns are available
#     if (!all(colnames(x$catchability_data) %in% colnames(new_catchability_data))) {
#       stop("Please ensure that all columns of `x$catchability_data` are present in `new_covariate_data`")
#     }
#     # Eliminate unnecessary columns
#     new_catchability_data <- new_catchability_data[, match(colnames(x$catchability_data), colnames(new_catchability_data))]
#     # Eliminate old-covariates that are also present in new_covariate_data
#     NN <- RANN::nn2(query = x$catchability_data[, c("Lat", "Lon", "Year")], data = new_catchability_data[, c("Lat", "Lon", "Year")], k = 1)
#     if (any(NN$nn.dist == 0)) {
#       x$catchability_data <- x$catchability_data[-which(NN$nn.dist == 0), , drop = FALSE]
#     }
#   }

#   # Process covariates
#   covariate_data <- rbind(x$covariate_data, new_covariate_data)
#   catchability_data <- rbind(x$catchability_data, new_catchability_data)

#   # Process inputs
#   PredTF_i <- c(x$data_list$PredTF_i, rep(1, length(t_i)))
#   b_i <- c(x$data_frame[, "b_i"], sample(c(0, 1), size = length(t_i), replace = TRUE))
#   c_iz <- rbind(matrix(x$data_frame[, grep("c_iz", names(x$data_frame))]), matrix(c_iz))
#   Lat_i <- c(x$data_frame[, "Lat_i"], Lat_i)
#   Lon_i <- c(x$data_frame[, "Lon_i"], Lon_i)
#   a_i <- c(x$data_frame[, "a_i"], a_i)
#   v_i <- c(x$data_frame[, "v_i"], v_i)
#   t_i <- c(x$data_frame[, "t_i"], t_i)
#   # assign("b_i", b_i, envir=.GlobalEnv)

#   # Build information regarding spatial location and correlation
#   message("\n### Re-making spatial information")
#   spatial_args_new <- list("anisotropic_mesh" = x$spatial_list$MeshList$anisotropic_mesh, "Kmeans" = x$spatial_list$Kmeans, "Lon_i" = Lon_i, "Lat_i" = Lat_i)
#   spatial_args_input <- combine_lists(input = spatial_args_new, default = x$input_args$spatial_args_input)
#   spatial_list <- do.call(what = make_spatial_info, args = spatial_args_input)

#   # Check spatial_list
#   if (!all.equal(spatial_list$MeshList, x$spatial_list$MeshList)) {
#     stop("`MeshList` generated during `predict.fit_model` doesn't match that of original fit; please email package author to report issue")
#   }

#   # Build data
#   # Do *not* restrict inputs to formalArgs(make_data) because other potential inputs are still parsed by make_data for backwards compatibility
#   message("\n### Re-making data object")
#   data_args_new <- list(
#     "c_iz" = c_iz, "b_i" = b_i, "a_i" = a_i, "v_i" = v_i, "PredTF_i" = PredTF_i,
#     "t_i" = t_i, "spatial_list" = spatial_list,
#     "covariate_data" = covariate_data, "catchability_data" = catchability_data
#   )
#   data_args_input <- combine_lists(input = data_args_new, default = x$input_args$data_args_input) # Do *not* use args_to_use
#   data_list <- do.call(what = make_data, args = data_args_input)
#   data_list$n_g <- 0

#   # Build object
#   message("\n### Re-making TMB object")
#   model_args_default <- list("TmbData" = data_list, "RunDir" = working_dir, "Version" = x$settings$Version, "RhoConfig" = x$settings$RhoConfig, "loc_x" = spatial_list$loc_x, "Method" = spatial_list$Method, "Map" = x$tmb_list$Map)
#   model_args_input <- combine_lists(
#     input = list("Parameters" = x$ParHat),
#     default = model_args_default, args_to_use = formalArgs(make_model)
#   )
#   tmb_list <- do.call(what = make_model, args = model_args_input)

#   # Extract output
#   Report <- tmb_list$Obj$report()
#   Y_i <- Report[[what]][(1 + nrow(x$data_frame)):length(Report$D_i)]

#   # sanity check
#   # if( all.equal(covariate_data,x$covariate_data) & Report$jnll!=x$Report$jnll){
#   if (do_checks == TRUE && (Report$jnll != x$Report$jnll)) {
#     message("Problem detected in `predict.fit_model`; returning outputs for diagnostic purposes")
#     Return <- list("Report" = Report, "data_list" = data_list)
#     return(Return)
#   }

#   # return prediction
#   return(Y_i)
# }

# Designing a projection function
#' @title Summarize prediction intervals from simulation results.
#'
#' @description Calculates and returns prediction intervals for a given variable at each projection time step.
#'
#' @param sim_obj = A simulation results list returned by \code{project.fit_model}.
#' @param what = Variable to summarize and plot
#' @param probs = A numeric vector of probabilities to calculate across simulation values at each projection time step.
#'
#' @return A dataframe with location and time information for each sample in `new_projection_data` and the projected \code{what} output variable prob results.
#'
#' @export
summary.sim_results<- function (vast_fit, sim_obj, resp_scale, nice_times = NULL, out_t_scale = NULL, nice_category_names = nice_category_names, climate_scenario = climate_scenario, 
    out_dir) {
    if (FALSE) {
        # tar_load(vast_fit)
        # tar_load(vast_projections)
        # sim_obj <- vast_projections
        # what <- "Index_ctl"
        # nice_times <- nice_times
        # out_t_scale <- "annual"
        # probs <- c(0.1, 0.5, 0.9)
        # mean_instead <- FALSE
        # nice_category_names <- nice_category_names
        # climate_scenario <- climate_scenario
        # out_dir <- paste0(res_root, "prediction_df")
        
        # Capelin
        date_dir<- here::here("2023-02-17/Capelin_BC/")
        vast_fit = fit_full
        sim_obj = uncert_res_full
        resp_scale = "raw"
        nice_times <- nice_times
        out_t_scale = NULL
        probs = c(0.1, 0.5, 0.9)
        mean_instead = FALSE
        nice_category_names = "Capelin_Random"
        climate_scenario = climate_scenario = paste0("gfdl", "_full")
        out_dir = date_dir

        # ## Cod -- COG is off...
        # date_dir <- "~/GitHub/mixedmod_projections/2022-10-25/Cod_BC/"
        # vast_fit = readRDS(paste0(date_dir, "SpST_mod_fit.rds"))
        # sim_obj = readRDS(file = paste0(date_dir, "SpST", "_random_ProjectionsList.rds"))
        # resp_scale = "raw"
        # nice_times <- as.Date(c(paste0(seq(from = 1985, to = 2100, by = 1), "-03-16"), paste0(seq(from = 1985, to = 2100, by = 1), "-07-16"), paste0(seq(from = 1985, to = 2100, by = 1), "-10-16")))
        # nice_times <- nice_times[order(nice_times)]
        # out_t_scale = NULL
        # nice_category_names = paste0("Cod_", "Base_None")
        # climate_scenario = paste0("EnvOnly_Base_5thpercentile", "_SSP5_85")
        # out_dir = date_dir
    }
    
    time_ind <- seq(from = 1, to = length(nice_times))
    time_labels <- nice_times
    index_regions_ind <- seq(from = 1, to = vast_fit$data_list$n_l)
    index_regions <- vast_fit$settings$strata.limits$STRATA[index_regions_ind]
    categories_ind <- seq(from = 1, to = vast_fit$data_list$n_c)
    grid_ind <- seq(from = 1, to = vast_fit$data_list$n_g)

    for (i in seq_along(sim_obj)) {

        if(FALSE){
            # Checking sim_obj 
            summary(sim_obj[[100]][["Index_ctl"]])
        }

        # Going to want the Index...
        ind_array <- array(unlist(sim_obj[[i]][which(names(sim_obj[[i]]) ==
            "Index_ctl")]),
            dim = c(unlist(vast_fit$data_list[c("n_c")]), n_t = length(nice_times), unlist(vast_fit$data_list[c("n_l")])),
            dimnames = list(
                categories_ind, time_labels,
                index_regions
            )  
        )
        ind_df <- data.frame(aperm(ind_array, c(2, 3, 1)))
        colnames(ind_df) <- gsub(".1", "", colnames(ind_df))
        ind_df$Sim_Scenario <- rep(paste0("Sim_", i), nrow(ind_df))
        ind_df$Time <- nice_times
        # ind_df[, 1] <- ifelse(ind_df[, 1] == 0, NA, ind_df[, 1])
        # ind_df <- na.omit(ind_df)
        ind_df <- ind_df %>%
            pivot_longer(., !c(
                Sim_Scenario,
                Time
            ), names_to = "Region", values_to = "Index")

        # Check
        if(FALSE){
            true <- vast_fit$Report$Index_ctl
            str(true)
            str(ind_df)
        }
        
        # Density
        dens_array <- array(unlist(sim_obj[[i]][which(names(sim_obj[[i]]) ==
            "D_gct")]),
            dim = c(unlist(vast_fit$data_list[c("n_g")]), (vast_fit$data_list[c("n_c")]), n_t = length(nice_times)),
            dimnames = list(grid_ind, categories_ind, time_labels)
        )
        dens_df <- data.frame(aperm(dens_array, c(1, 3, 2)))
        colnames(dens_df) <- nice_times
        dens_df$Lat <- vast_fit$spatial_list$latlon_g[, "Lat"]
        dens_df$Lon <- vast_fit$spatial_list$latlon_g[, "Lon"]
        dens_df <- dens_df %>%
            pivot_longer(.,
                !c(Lat, Lon),
                names_to = "Time", values_to = "D_gct"
            ) %>%
            arrange(Time, Lat, Lon)
        dens_df$Time <- as.Date(dens_df$Time)
        dens_df$Sim_Scenario <- paste0("Sim_", i)

        # Check
        if(FALSE){
            true <- vast_fit$Report$D_gct
            str(true)
            str(dens_df)
        }

        # Center of gravity
        cog_array <- array(unlist(sim_obj[[i]][which(names(sim_obj[[i]]) ==
            "mean_Z_ctm")]),
        dim = c(unlist(vast_fit$data_list[c("n_c")]), n_t = length(nice_times), 2),
        dimnames = list(
            categories_ind, time_labels,
            c("Lon", "Lat")
        )
        )
        cog_true_df <- data.frame(aperm(cog_array, c(2, 3, 1)))
        names(cog_true_df)<- c("Eastings", "Northings")
        cog_true_df$Time<- nice_times
        cog_true_df$Time <- as.Date(cog_true_df$Time)
        cog_true_df$Sim_Scenario <- paste0("Sim_", i)

        cog_df <- cog_from_dens(vast_fit = vast_fit, proj_dens = dens_df, proj_ind = ind_df)

        # Check
        if(FALSE){
            true_lon <- vast_fit$Report$mean_Z_ctm[,,1]
            str(true_lon)
            true_lat<- vast_fit$Report$mean_Z_ctm[,,2]
            str(true_lat)
            str(cog_df)
        }

        # Effective area
        eff_area_array <- array(unlist(sim_obj[[i]][which(names(sim_obj[[i]]) ==
            "effective_area_ctl")]),
            dim = c(unlist(vast_fit$data_list[c("n_c")]), n_t = length(nice_times), unlist(vast_fit$data_list[c("n_l")])),
            dimnames = list(
                categories_ind, time_labels,
                index_regions
            )  
        )
        eff_area_true_df <- data.frame(aperm(eff_area_array, c(2, 3, 1)))
        colnames(eff_area_true_df) <- gsub(".1", "", colnames(eff_area_true_df))
        eff_area_true_df$Sim_Scenario <- rep(paste0("Sim_", i), nrow(eff_area_true_df))
        eff_area_true_df$Time <- nice_times
        # ind_df[, 1] <- ifelse(ind_df[, 1] == 0, NA, ind_df[, 1])
        # ind_df <- na.omit(ind_df)
        eff_area_true_df <- eff_area_true_df %>%
            pivot_longer(., !c(
                Sim_Scenario,
                Time
            ), names_to = "Region", values_to = "Eff_Area")

        
        eff_area_df <- eff_area_from_dens(vast_fit = vast_fit, proj_dens = dens_df, proj_ind = ind_df)

        # Check
        if (FALSE) {
            vast_fit$Report$effective_area_ctl
            str(true_eff_area)
            str(eff_area_df)
        }
        
        if(resp_scale == "log"){
            ind_df$Index <- log(ind_df$Index)
            dens_df$D_gct <- log(dens_df$D_gct)
        }
        
        if (i == 1) {
            res_out_ind <- ind_df
            res_out_dens <- dens_df
            res_out_cog <- cog_df
            res_out_cog_true<- cog_true_df
            res_out_eff_area <- eff_area_df
            res_out_eff_area_true<- eff_area_true_df
        } else {
            res_out_ind <- bind_rows(res_out_ind, ind_df)
            res_out_dens <- bind_rows(res_out_dens, dens_df)
            res_out_cog <- bind_rows(res_out_cog, cog_df)
            res_out_cog_true <- bind_rows(res_out_cog_true, cog_true_df)
            res_out_eff_area <- bind_rows(res_out_eff_area, eff_area_df)
            res_out_eff_area_true<- bind_rows(res_out_eff_area_true, eff_area_true_df)
        }
    }

    # Calculate summaries across all runs
    res_out_ind <- res_out_ind %>%
        group_by(., Time, Region) %>%
        summarise(
            Prob_0.5 = quantile(Index, probs = 0.5, na.rm = TRUE),
            Prob_0.1 = quantile(Index, probs = 0.1, na.rm = TRUE),
            Prob_0.9 = quantile(Index, probs = 0.9, na.rm = TRUE)
        )
    
    res_out_dens <- res_out_dens %>%
        group_by(., Lat, Lon, Time) %>%
        summarise(
            Prob_0.5 = quantile(D_gct, probs = 0.5, na.rm = TRUE),
            Prob_0.1 = quantile(D_gct, probs = 0.1, na.rm = TRUE),
            Prob_0.9 = quantile(D_gct, probs = 0.9, na.rm = TRUE)
        )
    
    res_out_cog <- res_out_cog %>%
        group_by(., Time) %>%
        summarise(
            Lon_Prob_0.5 = quantile(Mean_Lon, probs = 0.5, na.rm = TRUE),
            Lon_Prob_0.1 = quantile(Mean_Lon, probs = 0.1, na.rm = TRUE),
            Lon_Prob_0.9 = quantile(Mean_Lon, probs = 0.9, na.rm = TRUE),
            Lat_Prob_0.5 = quantile(Mean_Lat, probs = 0.5, na.rm = TRUE),
            Lat_Prob_0.1 = quantile(Mean_Lat, probs = 0.1, na.rm = TRUE),
            Lat_Prob_0.9 = quantile(Mean_Lat, probs = 0.9, na.rm = TRUE)
        )
    
    res_out_cog_true<- res_out_cog_true %>%
        group_by(., Time) %>%
        summarise(
            Lon_Prob_0.5 = quantile(Eastings, probs = 0.5, na.rm = TRUE),
            Lon_Prob_0.1 = quantile(Eastings, probs = 0.1, na.rm = TRUE),
            Lon_Prob_0.9 = quantile(Eastings, probs = 0.9, na.rm = TRUE),
            Lat_Prob_0.5 = quantile(Northings, probs = 0.5, na.rm = TRUE),
            Lat_Prob_0.1 = quantile(Northings, probs = 0.1, na.rm = TRUE),
            Lat_Prob_0.9 = quantile(Northings, probs = 0.9, na.rm = TRUE)
        )
    
    res_out_eff_area <- res_out_eff_area %>%
        group_by(., Time, Region) %>%
        summarise(
            Prob_0.5 = quantile(Eff_Area, probs = 0.5, na.rm = TRUE),
            Prob_0.1 = quantile(Eff_Area, probs = 0.1, na.rm = TRUE),
            Prob_0.9 = quantile(Eff_Area, probs = 0.9, na.rm = TRUE)
        )
    
    res_out_eff_area_true <- res_out_eff_area_true %>%
        group_by(., Time, Region) %>%
        summarise(
            Prob_0.5 = quantile(Eff_Area, probs = 0.5, na.rm = TRUE),
            Prob_0.1 = quantile(Eff_Area, probs = 0.1, na.rm = TRUE),
            Prob_0.9 = quantile(Eff_Area, probs = 0.9, na.rm = TRUE)
        )
    
    res_out<- list("Index" = res_out_ind, "Dens" = res_out_dens, "COG" = res_out_cog, "COG_True" = res_out_cog_true, "EffArea" = res_out_eff_area, "EffArea_True" = res_out_eff_area_true)
 
    saveRDS(res_out, file = paste0(out_dir, "/", nice_category_names, 
        "_", climate_scenario, ".rds"))
    return(res_out)
}
    

match_strata_fn_aja <- function(points, strata_dataframe, index_shapes) {
  if (FALSE) {
    points <- Tmp
    l <- 1
    strata_dataframe <- strata.limits[l, , drop = FALSE]
    index_shapes <- index_shapes
  }
  if (is.null(index_shapes)) {
    # Default all strata
    match_latitude_TF <- match_longitude_TF <- match_depth_TF <- rep(TRUE, nrow(strata_dataframe))
    if (all(c("south_border", "north_border") %in% names(strata_dataframe))) {
      match_latitude_TF <- as.numeric(x["BEST_LAT_DD"]) > strata_dataframe[, "south_border"] & as.numeric(x["BEST_LAT_DD"]) <= strata_dataframe[, "north_border"]
    }
    if (all(c("west_border", "east_border") %in% names(strata_dataframe))) {
      match_longitude_TF <- as.numeric(x["BEST_LON_DD"]) > strata_dataframe[, "west_border"] & as.numeric(x["BEST_LON_DD"]) <= strata_dataframe[, "east_border"]
    }
    if (all(c("shallow_border", "deep_border") %in% names(strata_dataframe))) {
      match_depth_TF <- as.numeric(x["BEST_DEPTH_M"]) > strata_dataframe[, "shallow_border"] & as.numeric(x["BEST_DEPTH_M"]) <= strata_dataframe[, "deep_border"]
    }
    # Return stuff
    Char <- as.character(strata_dataframe[match_latitude_TF & match_longitude_TF & match_depth_TF, "STRATA"])
    return(ifelse(length(Char) == 0, NA, Char))
  }

  # Andrew edit...
  if (!is.null(index_shapes)) {
    Tmp_sf <- data.frame(points) %>%
      st_as_sf(., coords = c("BEST_LON_DD", "BEST_LAT_DD"), crs = st_crs(index_shapes), remove = FALSE)
    match_shape <- Tmp_sf %>%
      st_join(., index_shapes, join = st_within) %>%
      mutate(., "Row_ID" = seq(from = 1, to = nrow(.))) %>%
      st_drop_geometry() %>%
      dplyr::select(., Region) %>%
      as.vector()
    return(match_shape)
  }
}

Prepare_User_Extrapolation_Data_Fn_aja <- function(input_grid, strata.limits = NULL, projargs = NA, zone = NA, flip_around_dateline = TRUE, index_shapes, ...) {
  if (FALSE) {
    # Run make_extrapolation_info_aja first...
    strata.limits <- strata.limits
    input_grid <- input_grid
    projargs <- projargs
    zone <- zone
    flip_around_dateline <- flip_around_dateline
    index_shapes <- index_shapes
  }
  if (is.null(strata.limits)) {
    strata.limits <- data.frame(STRATA = "All_areas")
  }
  message("Using strata ", strata.limits)
  Data_Extrap <- input_grid
  Area_km2_x <- Data_Extrap[, "Area_km2"]
  Tmp <- cbind(BEST_LAT_DD = Data_Extrap[, "Lat"], BEST_LON_DD = Data_Extrap[, "Lon"])
  if ("Depth" %in% colnames(Data_Extrap)) {
    Tmp <- cbind(Tmp, BEST_DEPTH_M = Data_Extrap[, "Depth"])
  }
  a_el <- as.data.frame(matrix(NA, nrow = nrow(Data_Extrap), ncol = nrow(strata.limits), dimnames = list(NULL, strata.limits[, "STRATA"])))
  for (l in 1:ncol(a_el)) {
    a_el[, l] <- match_strata_fn_aja(points = Tmp, strata_dataframe = strata.limits[l, , drop = FALSE], index_shapes = index_shapes[index_shapes$Region == as.character(strata.limits[l, , drop = FALSE]), ])
    a_el[, l] <- ifelse(is.na(a_el[, l]), 0, Area_km2_x)
  }
  tmpUTM <- project_coordinates(X = Data_Extrap[, "Lon"], Y = Data_Extrap[, "Lat"], projargs = projargs, zone = zone, flip_around_dateline = flip_around_dateline)
  Data_Extrap <- cbind(Data_Extrap, Include = 1)
  if (all(c("E_km", "N_km") %in% colnames(Data_Extrap))) {
    Data_Extrap[, c("E_km", "N_km")] <- tmpUTM[, c("X", "Y")]
  } else {
    Data_Extrap <- cbind(Data_Extrap, E_km = tmpUTM[, "X"], N_km = tmpUTM[, "Y"])
  }
  Return <- list(a_el = a_el, Data_Extrap = Data_Extrap, zone = attr(tmpUTM, "zone"), projargs = attr(tmpUTM, "projargs"), flip_around_dateline = flip_around_dateline, Area_km2_x = Area_km2_x)
  return(Return)
}

make_extrapolation_info_aja <- function(Region, projargs = NA, zone = NA, strata.limits = data.frame(STRATA = "All_areas"), create_strata_per_region = FALSE, max_cells = NULL, input_grid = NULL, observations_LL = NULL, grid_dim_km = c(2, 2), maximum_distance_from_sample = NULL, grid_in_UTM = TRUE, grid_dim_LL = c(0.1, 0.1), region = c("south_coast", "west_coast"), strata_to_use = c("SOG", "WCVI", "QCS", "HS", "WCHG"), epu_to_use = c("All", "Georges_Bank", "Mid_Atlantic_Bight", "Scotian_Shelf", "Gulf_of_Maine", "Other")[1], survey = "Chatham_rise", surveyname = "propInWCGBTS", flip_around_dateline, nstart = 100, area_tolerance = 0.05, backwards_compatible_kmeans = FALSE, DirPath = paste0(getwd(), "/"), index_shapes, ...) {
  if (FALSE) {
    # First run fit_model_aja...
    Region <- settings$Region
    projargs <- NA
    zone <- settings$zone
    strata.limits <- settings$strata.limits
    create_strata_per_region <- FALSE
    max_cells <- settings$max_cells
    input_grid <- input_grid
    observations_LL <- NULL
    grid_dim_km <- settings$grid_size_km
    maximum_distance_from_sample <- NULL
    index_shapes <- index_shapes
  }
  if (is.null(max_cells)) {
    max_cells <- Inf
  }
  for (rI in seq_along(Region)) {
    Extrapolation_List <- NULL
    if (tolower(Region[rI]) == "user") {
      if (is.null(input_grid)) {
        stop("Because you're using a user-supplied region, please provide 'input_grid' input")
      }
      if (!(all(c("Lat", "Lon", "Area_km2") %in% colnames(input_grid)))) {
        stop("'input_grid' must contain columns named 'Lat', 'Lon', and 'Area_km2'")
      }
      if (missing(flip_around_dateline)) {
        flip_around_dateline <- FALSE
      }
      Extrapolation_List <- Prepare_User_Extrapolation_Data_Fn_aja(strata.limits = strata.limits, input_grid = input_grid, projargs = projargs, zone = zone, flip_around_dateline = flip_around_dateline, index_shapes = index_shapes, ...)
    }
    if (is.null(Extrapolation_List)) {
      if (is.null(observations_LL)) {
        stop("Because you're using a new Region[rI], please provide 'observations_LL' input with columns named `Lat` and `Lon`")
      }
      if (missing(flip_around_dateline)) {
        flip_around_dateline <- FALSE
      }
      Extrapolation_List <- Prepare_Other_Extrapolation_Data_Fn(strata.limits = strata.limits, observations_LL = observations_LL, grid_dim_km = grid_dim_km, maximum_distance_from_sample = maximum_distance_from_sample, grid_in_UTM = grid_in_UTM, grid_dim_LL = grid_dim_LL, projargs = projargs, zone = zone, flip_around_dateline = flip_around_dateline, ...)
    }
    if (rI == 1) {
      Return <- Extrapolation_List
    } else {
      Return <- combine_extrapolation_info(Return, Extrapolation_List, create_strata_per_region = create_strata_per_region)
    }
  }
  if (max_cells < nrow(Return$Data_Extrap)) {
    message("# Reducing extrapolation-grid from ", nrow(Return$Data_Extrap), " to ", max_cells, " cells for Region(s): ", paste(Region, collapse = ", "))
    loc_orig <- Return$Data_Extrap[, c("E_km", "N_km")]
    loc_orig <- loc_orig[which(Return$Area_km2_x > 0), ]
    Kmeans <- make_kmeans(n_x = max_cells, loc_orig = loc_orig, nstart = nstart, randomseed = 1, iter.max = 1000, DirPath = DirPath, Save_Results = TRUE, kmeans_purpose = "extrapolation", backwards_compatible_kmeans = backwards_compatible_kmeans)
    Kmeans[["cluster"]] <- RANN::nn2(data = Kmeans[["centers"]], query = Return$Data_Extrap[, c("E_km", "N_km")], k = 1)$nn.idx[, 1]
    aggregate_vector <- function(values_x, index_x, max_index, FUN = sum) {
      tapply(values_x, INDEX = factor(index_x, levels = 1:max_index), FUN = FUN)
    }
    a_el <- matrix(NA, nrow = max_cells, ncol = ncol(Return$a_el))
    for (lI in 1:ncol(Return$a_el)) {
      a_el[, lI] <- aggregate_vector(values_x = Return$a_el[, lI], index_x = Kmeans$cluster, max_index = max_cells)
    }
    Area_km2_x <- aggregate_vector(values_x = Return$Area_km2_x, index_x = Kmeans$cluster, max_index = max_cells)
    Include <- aggregate_vector(values_x = Return$Data_Extrap[, "Include"], index_x = Kmeans$cluster, max_index = max_cells, FUN = function(vec) {
      any(vec > 0)
    })
    lonlat_g <- project_coordinates(X = Kmeans$centers[, "E_km"], Y = Kmeans$centers[, "N_km"], projargs = "+proj=longlat +ellps=WGS84", origargs = Return$projargs)
    Data_Extrap <- cbind(Lon = lonlat_g[, 1], Lat = lonlat_g[, 2], Include = Include, Kmeans$centers)
    Return <- list(a_el = a_el, Data_Extrap = Data_Extrap, zone = Return$zone, projargs = Return$projargs, flip_around_dateline = Return$flip_around_dateline, Area_km2_x = Area_km2_x)
  }
  if (length(Region) > 1 & create_strata_per_region == TRUE) {
    Return$a_el <- cbind(Total = rowSums(Return$a_el), Return$a_el)
  }
  class(Return) <- "make_extrapolation_info"
  return(Return)
}

fit_model_aja <- function(settings, Method, Lat_i, Lon_i, t_i, b_i, a_i, c_iz = rep(0, length(b_i)), v_i = rep(0, length(b_i)), category_names = NULL, working_dir = paste0(getwd(), "/"), X1config_cp = NULL, X2config_cp = NULL, covariate_data, X1_formula = ~0, X2_formula = ~0, Q1config_k = NULL, Q2config_k = NULL, catchability_data, Q1_formula = ~0, Q2_formula = ~0, newtonsteps = 1, silent = TRUE, build_model = TRUE, run_model = TRUE, test_fit = TRUE, framework = "TMBad", use_new_epsilon = TRUE, ...) {
  if (FALSE) {
    # Run vast_fit_sdm first...

    "settings" <- settings
    "input_grid" <- extrap_grid
    "Lat_i" <- sample_data[, "Lat"]
    "Lon_i" <- sample_data[, "Lon"]
    "t_i" <- sample_data[, "Year"]
    "c_i" <- rep(0, nrow(sample_data))
    "b_i" <- sample_data[, "Biomass"]
    "v_i" <- rep(0, length(b_i))
    "a_i" <- sample_data[, "Swept"]
    "PredTF_i" <- sample_data[, "Pred_TF"]
    "X1config_cp" <- Xconfig_list[["X1config_cp"]]
    "X2config_cp" <- Xconfig_list[["X2config_cp"]]
    "covariate_data" <- covariate_data
    "X1_formula" <- X1_formula
    "X2_formula" <- X2_formula
    "X_contrasts" <- X_contrasts
    "catchability_data" <- catchability_data
    "Q1_formula" <- Q1_formula
    "Q2_formula" <- Q2_formula
    "Q1config_k" <- Xconfig_list[["Q1config_k"]]
    "Q2config_k" <- Xconfig_list[["Q2config_k"]]
    "newtonsteps" <- 1
    "getsd" <- TRUE
    "getReportCovariance" <- TRUE
    "run_model" <- FALSE
    "test_fit" <- FALSE
    "Use_REML" <- FALSE
    "getJointPrecision" <- FALSE
    "index_shapes" <- index_shapes

    # Now, go into make_extrapolation_info_aja
  }
  extra_args <- list(...)
  extra_args <- c(extra_args, extra_args$extrapolation_args, extra_args$spatial_args, extra_args$optimize_args, extra_args$model_args)
  start_time = Sys.time()
  data_frame <- data.frame(Lat_i = Lat_i, Lon_i = Lon_i, a_i = a_i, v_i = v_i, b_i = b_i, t_i = t_i, c_iz = c_iz)
  year_labels <- seq(min(t_i), max(t_i))
  years_to_plot <- which(year_labels %in% t_i)
  if(is.null(category_names)) category_names = paste0( 1:(max(c_iz,na.rm=TRUE)+1) )
  message("\n### Writing output from `fit_model` in directory: ", working_dir)
  dir.create(working_dir, showWarnings = FALSE, recursive = TRUE)
  capture.output(settings, file = file.path(working_dir, "settings.txt"))
  if (is.null(extra_args$extrapolation_list)) {
    message("\n### Making extrapolation-grid")
    extrapolation_args_default <- list(Region = settings$Region, strata.limits = settings$strata.limits, zone = settings$zone, max_cells = settings$max_cells, DirPath = working_dir)
    extrapolation_args_input <- combine_lists(input = extra_args, default = extrapolation_args_default, args_to_use = formalArgs(make_extrapolation_info_aja))
    extrapolation_list <- do.call(what = make_extrapolation_info_aja, args = extrapolation_args_input)
  } else {
    extrapolation_args_input = NULL
    extrapolation_list = extra_args$extrapolation_list
  }

  if (is.null(extra_args$spatial_list)) {
    message("\n### Making spatial information")
    spatial_args_default <- list(grid_size_km = settings$grid_size_km, n_x = settings$n_x, Method = Method, Lon_i = Lon_i, Lat_i = Lat_i, Extrapolation_List = extrapolation_list, DirPath = working_dir, Save_Results = TRUE, fine_scale = settings$fine_scale, knot_method = settings$knot_method)
    spatial_args_input <- combine_lists(input = extra_args, default = spatial_args_default, args_to_use = c(formalArgs(make_spatial_info), formalArgs(INLA::inla.mesh.create)))
    spatial_list <- do.call(what = make_spatial_info, args = spatial_args_input)
  } else {
    spatial_args_input = NULL
    spatial_list = extra_args$spatial_list
  }

  if(is.null(extra_args$data_list)){
    message("\n### Making data object")
    if (missing(covariate_data)) {
      covariate_data <- NULL
    }
    if (missing(catchability_data)) {
      catchability_data <- NULL
    }
    data_args_default <- list(Version = settings$Version, FieldConfig = settings$FieldConfig, OverdispersionConfig = settings$OverdispersionConfig, RhoConfig = settings$RhoConfig, VamConfig = settings$VamConfig, ObsModel = settings$ObsModel, c_iz = c_iz, b_i = b_i, a_i = a_i, v_i = v_i, s_i = spatial_list$knot_i - 1, t_i = t_i, spatial_list = spatial_list, Options = settings$Options, Aniso = settings$use_anisotropy, X1config_cp = X1config_cp, X2config_cp = X2config_cp, covariate_data = covariate_data, X1_formula = X1_formula, X2_formula = X2_formula, Q1config_k = Q1config_k, Q2config_k = Q2config_k, catchability_data = catchability_data, Q1_formula = Q1_formula, Q2_formula = Q2_formula)
    data_args_input <- combine_lists(input = extra_args, default = data_args_default)
    data_list <- do.call(what = make_data, args = data_args_input)
  } else {
    data_args_input = NULL
    data_list = extra_args$data_list
  }
  
  message("\n### Making TMB object")
  model_args_default <- list(TmbData = data_list, RunDir = working_dir, Version = settings$Version, RhoConfig = settings$RhoConfig, loc_x = spatial_list$loc_x, Method = spatial_list$Method, build_model = build_model, framework = framework)
  model_args_input <- combine_lists(input = extra_args, default = model_args_default, args_to_use = formalArgs(make_model))
  tmb_list <- do.call(what = make_model, args = model_args_input)
  if (run_model == FALSE | build_model == FALSE) {
    input_args <- list(extra_args = extra_args, extrapolation_args_input = extrapolation_args_input, model_args_input = model_args_input, spatial_args_input = spatial_args_input, data_args_input = data_args_input)
    Return <- list(data_frame = data_frame, extrapolation_list = extrapolation_list, spatial_list = spatial_list, data_list = data_list, tmb_list = tmb_list, year_labels = year_labels, years_to_plot = years_to_plot, settings = settings, input_args = input_args)
    class(Return) <- "fit_model"
    return(Return)
  }
  if (silent == TRUE) {
    tmb_list$Obj$env$beSilent()
  }
  if (test_fit == TRUE) {
    message("\n### Testing model at initial values")
    LogLike0 <- tmb_list$Obj$fn(tmb_list$Obj$par)
    Gradient0 <- tmb_list$Obj$gr(tmb_list$Obj$par)
    if (any(Gradient0 == 0)) {
      message("\n")
      stop("Please check model structure; some parameter has a gradient of zero at starting values\n",
        call. = FALSE
      )
    } else {
      message("Looks good: All fixed effects have a nonzero gradient")
    }
  }
  message("\n### Estimating parameters")
  optimize_args_default1 <- list(lower = tmb_list$Lower, upper = tmb_list$Upper, loopnum = 1)
  optimize_args_default1 <- combine_lists(default = optimize_args_default1, input = extra_args, args_to_use = formalArgs(TMBhelper::fit_tmb))
  optimize_args_input1 <- list(obj = tmb_list$Obj, savedir = NULL, newtonsteps = 0, bias.correct = FALSE, control = list(eval.max = 50000, iter.max = 50000, trace = 1), quiet = TRUE, getsd = FALSE)
  optimize_args_input1 <- combine_lists(default = optimize_args_default1, input = optimize_args_input1, args_to_use = formalArgs(TMBhelper::fit_tmb))
  parameter_estimates1 <- do.call(what = TMBhelper::fit_tmb, args = optimize_args_input1)
  if (exists("check_fit") & test_fit == TRUE) {
    problem_found <- VAST::check_fit(parameter_estimates1)
    if (problem_found == TRUE) {
      message("\n")
      stop("Please change model structure to avoid problems with parameter estimates and then re-try; see details in `?check_fit`\n", call. = FALSE)
    }
  }
  if ((use_new_epsilon == TRUE) & (settings$bias.correct == TRUE) & (framework == "TMBad") & ("Index_ctl" %in% settings$vars_to_correct)) {
    settings$vars_to_correct = setdiff(settings$vars_to_correct, c("Index_ctl", "Index_cyl"))
    if (length(settings$vars_to_correct) == 0) {
      settings$bias.correct = FALSE
    }
    settings$vars_to_correct = c(settings$vars_to_correct, "eps_Index_ctl")
  }
  
  optimize_args_default2 <- list(obj = tmb_list$Obj, lower = tmb_list$Lower, upper = tmb_list$Upper, savedir = working_dir, bias.correct = settings$bias.correct, newtonsteps = newtonsteps, bias.correct.control = list(sd = FALSE, split = NULL, nsplit = 1, vars_to_correct = settings$vars_to_correct), control = list(eval.max = 10000, iter.max = 10000, trace = 1), loopnum = 1, getJointPrecision = TRUE, start_time_elapsed = parameter_estimates1$time_for_run)
  optimize_args_input2 <- combine_lists(input = extra_args, default = optimize_args_default2, args_to_use = formalArgs(TMBhelper::fit_tmb))
  optimize_args_input2 <- combine_lists(input = list(startpar = parameter_estimates1$par), default = optimize_args_input2)
  parameter_estimates2 <- do.call(what = TMBhelper::fit_tmb, args = optimize_args_input2)

  if ((use_new_epsilon == TRUE) & (framework == "TMBad") & ("eps_Index_ctl" %in% settings$vars_to_correct) & !is.null(parameter_estimates2$SD)) {
    message("\n### Applying faster epsilon bias-correction estimator")
    fit = list(parameter_estimates = parameter_estimates2, tmb_list = tmb_list, input_args = list(model_args_input = model_args_input))
    parameter_estimates2$SD = apply_epsilon(fit)
  }

  if ("par" %in% names(parameter_estimates2)) {
    if (!is.null(tmb_list$Obj$env$intern) && tmb_list$Obj$env$intern == TRUE) {
      Report = as.list(tmb_list$Obj$env$reportenv)
    } else {
      Report = tmb_list$Obj$report()
    }
    ParHat = tmb_list$Obj$env$parList(parameter_estimates2$par)
    Report = amend_output(Report = Report, TmbData = data_list, Map = tmb_list$Map, Sdreport = parameter_estimates2$SD, year_labels = year_labels, category_names = category_names, extrapolation_list = extrapolation_list)
  } else {
    Report = ParHat = "Model is not converged"
  }
  

  input_args <- list(extra_args = extra_args, extrapolation_args_input = extrapolation_args_input, model_args_input = model_args_input, spatial_args_input = spatial_args_input, optimize_args_input1 = optimize_args_input1, optimize_args_input2 = optimize_args_input2, data_args_input = data_args_input)
  Return <- list(data_frame = data_frame, extrapolation_list = extrapolation_list, spatial_list = spatial_list, data_list = data_list, tmb_list = tmb_list, parameter_estimates = parameter_estimates2, Report = Report, ParHat = ParHat, year_labels = year_labels, years_to_plot = years_to_plot, settings = settings, input_args = input_args, X1config_cp = X1config_cp, X2config_cp = X2config_cp, covariate_data = covariate_data, X1_formula = X1_formula, X2_formula = X2_formula, Q1config_k = Q1config_k, Q2config_k = Q1config_k, catchability_data = catchability_data, Q1_formula = Q1_formula, Q2_formula = Q2_formula, total_time = Sys.time() - start_time)

  Return$effects <- list()
  if (!is.null(catchability_data)) {
    catchability_data_full <- data.frame(catchability_data, linear_predictor = 0)
    Q1_formula_full <- update.formula(Q1_formula, linear_predictor ~ . + 0)
    call_Q1 <- lm(Q1_formula_full, data = catchability_data_full)$call
    Q2_formula_full <- update.formula(Q2_formula, linear_predictor ~ . + 0)
    call_Q2 <- lm(Q2_formula_full, data = catchability_data_full)$call
    Return$effects <- c(Return$effects, list(call_Q1 = call_Q1, call_Q2 = call_Q2, catchability_data_full = catchability_data_full))
  }
  if (!is.null(covariate_data)) {
    covariate_data_full <- data.frame(covariate_data, linear_predictor = 0)
    X1_formula_full <- update.formula(X1_formula, linear_predictor ~ . + 0)
    call_X1 <- lm(X1_formula_full, data = covariate_data_full)$call
    X2_formula_full <- update.formula(X2_formula, linear_predictor ~ . + 0)
    call_X2 <- lm(X2_formula_full, data = covariate_data_full)$call
    Return$effects <- c(Return$effects, list(call_X1 = call_X1, call_X2 = call_X2, covariate_data_full = covariate_data_full))
  }

  # Add stuff for marginaleffects package
  Return$last.par.best<- tmb_list$Obj$env$last.par.best

  # Class and return
  class(Return) <- "fit_model"
  return(Return)
}

vast_read_region_shape <- function(region_shapefile_dir) {
  region_file <- list.files(region_shapefile_dir, pattern = ".shp$", full.names = TRUE)
  region_sf <- st_read(region_file)
  # region_sf$geometry <- region_sf$geometry %>%
  #   s2::s2_rebuild() %>%
  #   sf::st_as_sfc()
  return(region_sf)
}

vast_read_index_shapes <- function(index_shapefiles_dir) {
  if (FALSE) {
    index_shapefiles_dir <- "~/GitHub/sdm_workflow/scratch/aja/TargetsSDM/data/supporting/index_shapefiles/"
    index_shapefiles_dir <- "~/data/supporting/index_shapefiles/"
  }

  index_files <- list.files(index_shapefiles_dir, pattern = ".shp$", full.names = TRUE)

  for (i in seq_along(index_files)) {
    index_shapes_temp <- st_read(index_files[i])
    # index_shapes_temp$geometry <- index_shapes_temp$geometry %>%
    #   s2::s2_rebuild() %>%
    #   sf::st_as_sfc()
    
    if (i == 1) {
      index_shapes_out <- index_shapes_temp
    } else {
      index_shapes_out <- bind_rows(index_shapes_out, index_shapes_temp)
    }
  }
  return(index_shapes_out)
}

######
## Getting and plotting biomass index time series
######
get_vast_index_timeseries <- function(vast_fit, all_times, nice_category_names, index_scale = c("raw", "log"), out_dir) {
  if (FALSE) {
    tar_load(vast_fit)
    tar_load(vast_seasonal_data)
    all_times <- levels(vast_seasonal_data$VAST_YEAR_SEASON)
    nice_category_names <- nice_category_names
    index_scale <- "raw"
    out_dir <- paste0(res_root, "tables")

    tar_load(vast_fit)
    vast_fit <- vast_fitted
    nice_category_names <- "Atlantic halibut"
    index_scale <- "raw"
    out_dir <- here::here("scratch/aja/TargetsSDM/results/tables")

    vast_fit <- vast_fit
    all_times <- unique(vast_samp_dat$Year)
    nice_category_names <- "Atlantic_halibut_habcovs"
    index_scale <- c("raw")
    out_dir <- here::here("", "results/tables")

    vast_fit = fit
    all_times = hist_times
    nice_category_names = "Atlantic cod"
    index_scale = "log"
    out_dir = date_dir
  }

  TmbData <- vast_fit$data_list
  Sdreport <- vast_fit$parameter_estimates$SD

  # Time series steps
  time_ind <- 1:TmbData$n_t
  time_labels <- sort(unique(vast_fit$data_frame$t_i)[time_ind])

  # Missing year issue?
  if(length(seq(from = min(time_labels), to = max(time_labels))) != length(time_labels)){
    time_labels<- seq(from = min(unique(vast_fit$data_frame$t_i)[time_ind], na.rm = TRUE), to = max(unique(vast_fit$data_frame$t_i)[time_ind], na.rm = TRUE))
  }

  # Index regions
  index_regions_ind <- 1:TmbData$n_l
  index_regions <- vast_fit$settings$strata.limits$STRATA[index_regions_ind]

  # Categories
  categories_ind <- 1:TmbData$n_c

  # Get the index information
  SD <- TMB::summary.sdreport(Sdreport)
  # SD_stderr<- TMB:::as.list.sdreport(Sdreport, what = "Std. Error", report = TRUE)
  # SD_estimate<- TMB:::as.list.sdreport(Sdreport, what = "Estimate", report = TRUE)
  # if(vast_fit$settings$bias.correct == TRUE && "unbiased" %in% names(Sdreport)){
  #   SD_estimate_biascorrect<- TMB:::as.list.sdreport(Sdreport, what = "Std. (bias.correct)", report = TRUE)
  # }
  #
  # Now, populate array with values
  Index_ctl <- log_Index_ctl <- array(NA, dim = c(unlist(TmbData[c("n_c", "n_t", "n_l")]), 2), dimnames = list(categories_ind, time_labels, index_regions, c("Estimate", "Std. Error")))

  if (index_scale == "raw") {
    if ("unbiased" %in% names(Sdreport)) {
      Index_ctl[] <- SD[which(rownames(SD) == "Index_ctl"), c("Est. (bias.correct)", "Std. Error")]
    } else {
      Index_ctl[] <- SD[which(rownames(SD) == "Index_ctl"), c("Estimate", "Std. Error")]
    }

    index_res_array <- Index_ctl
  } else {
    if ("unbiased" %in% names(Sdreport)) {
      log_Index_ctl[] <- log(SD[which(rownames(SD) == "Index_ctl"), c("Est. (bias.correct)", "Std. Error")])
    } else {
      log_Index_ctl[] <- SD[which(rownames(SD) == "ln_Index_ctl"), c("Estimate", "Std. Error")]
    }
    index_res_array <- log_Index_ctl
  }
  

  # Data manipulation to get out out the array and to something more "plottable"
  for (i in seq_along(categories_ind)) {
    index_array_temp <- index_res_array[i, , , ]

    if (length(dim(index_array_temp)) == 2) {
      index_res_temp_est <- data.frame("Time" = as.numeric(rownames(index_array_temp)), "Category" = categories_ind[i], index_array_temp[, 1])
      colnames(index_res_temp_est)[3] <- "Index_Estimate"
      index_res_temp_est$Index_Region <- index_regions

      index_res_temp_sd <- data.frame("Time" = as.numeric(rownames(index_array_temp)), "Category" = categories_ind[i], index_array_temp[, 2])
      colnames(index_res_temp_sd)[3] <- "Index_SD"
      index_res_temp_sd$Index_Region <- index_regions

      index_res_temp_out <- index_res_temp_est %>%
        left_join(., index_res_temp_sd)
    } else {
      index_res_temp_est <- data.frame("Time" = as.numeric(rownames(index_array_temp[, , 1])), "Category" = categories_ind[i], index_array_temp[, , 1]) %>%
        pivot_longer(cols = -c(Time, Category), names_to = "Index_Region", values_to = "Index_Estimate")
      index_res_temp_sd <- data.frame("Time" = as.numeric(rownames(index_array_temp[, , 1])), "Category" = categories_ind[i], index_array_temp[, , 2]) %>%
        pivot_longer(cols = -c(Time, Category), names_to = "Index_Region", values_to = "Index_SD")
      index_res_temp_out <- index_res_temp_est %>%
        left_join(., index_res_temp_sd)
    }

    if (i == 1) {
      index_res_out <- index_res_temp_out
    } else {
      index_res_out <- bind_rows(index_res_out, index_res_temp_out)
    }


    # if(dim(index_array_temp)[2] == 3){
    #   index_res_temp_est<- data.frame("Time" = as.numeric(rownames(index_array_temp[,,1])), "Category" = categories_ind[i], index_array_temp[,,1]) %>%
    #     pivot_longer(cols = -c(Time, Category), names_to = "Index_Region", values_to = "Index_Estimate")
    #   index_res_temp_sd<- data.frame("Time" = as.numeric(rownames(index_array_temp[,,1])), "Category" = categories_ind[i], index_array_temp[,,2]) %>%
    #     pivot_longer(cols = -c(Time, Category), names_to = "Index_Region", values_to = "Index_SD")
    #   index_res_temp_out<- index_res_temp_est %>%
    #     left_join(., index_res_temp_sd)
    #
    #   if(i == 1){
    #     index_res_out<- index_res_temp_out
    #   } else {
    #     index_res_out<- bind_rows(index_res_out, index_res_temp_out)
    #   }
    # } else if(as.numeric(dim(index_array_temp)[2]) == 2){
    #   index_res_temp_est<- data.frame("Time" = as.numeric(rownames(index_array_temp[,,1])), "Category" = categories_ind[i], index_array_temp[,,1]) %>%
    #     pivot_longer(cols = -c(Time, Category), names_to = "Index_Region", values_to = "Index_Estimate")
    #   index_res_temp_sd<- data.frame("Time" = as.numeric(rownames(index_array_temp[,,1])), "Category" = categories_ind[i], index_array_temp[,,2]) %>%
    #     pivot_longer(cols = -c(Time, Category), names_to = "Index_Region", values_to = "Index_SD")
    #   index_res_temp_out<- index_res_temp_est %>%
    #     left_join(., index_res_temp_sd)
    #
    #   if(i == 1){
    #     index_res_out<- index_res_temp_out
    #   } else {
    #     index_res_out<- bind_rows(index_res_out, index_res_temp_out)
    #   }
    # }
  }

  # if(!is.null(vast_fit$covariate_data)){
  #   year_start<- min(as.numeric(as.character(vast_fit$covariate_data$Year_Cov)))
  #
  #   if(any(grepl("Season", vast_fit$X1_formula))){
  #     seasons<- nlevels(unique(vast_fit$covariate_data$Season))
  #     if(seasons == 3 & max(time_labels) == 347){
  #       time_labels_use<- paste(rep(seq(from = year_start, to = 2100), each = 3), rep(c("SPRING", "SUMMER", "FALL")), sep = "-")
  #     }
  #   } else {
  #     time_labels_use<- paste(rep(seq(from = year_start, to = 2100), each = 1), rep(c("FALL")), sep = "-")
  #   }
  #
  #   index_res_out$Date<- factor(rep(time_labels_use, length(index_regions)), levels = time_labels_use)
  #
  # } else {
  #   # Just basic years...
  #   time_labels_use<- seq(from = min(vast_fit$year_labels), to = max(vast_fit$year_labels))
  #   index_res_out$Date<- factor(rep(time_labels_use, each = length(index_regions)), levels = time_labels_use)
  # }
  #

  index_res_out$Date <- rep(factor(as.character(all_times), levels = as.character(all_times)), each = length(unique(index_res_out$Index_Region)))

  # Date info
  index_res_out <- index_res_out %>%
    mutate(., Year = as.numeric(gsub("([0-9]+).*$", "\\1", Date)))

  if (any(str_detect(as.character(index_res_out$Date), "[:upper:]"))) {
    index_res_out$Date <- as.Date(paste(index_res_out$Year, ifelse(grepl("SPRING", index_res_out$Date), "-04-15",
      ifelse(grepl("SUMMER", index_res_out$Date), "-07-15", "-10-15")
    ), sep = ""))
  } else {
    # Staggered data sequence...
    month_day_seqs <- seq.Date(from = as.Date(paste0(min(index_res_out$Year), "-05-01")), to = as.Date(paste0(min(index_res_out$Year), "-07-30")), by = "day")
    index_res_out$Date <- as.Date(paste(index_res_out$Year, format(sample(month_day_seqs, size = nrow(index_res_out), replace = TRUE), "%m-%d"), sep = "-"))
  }

  # Save and return it
  write.csv(index_res_out, file = paste(out_dir, "/Biomass_Index_", index_scale, "_", nice_category_names, ".csv", sep = ""))
  return(index_res_out)
}

plot_vast_index_timeseries <- function(index_res_df, year_stop = NULL, index_scale, nice_category_names, nice_xlab, nice_ylab, paneling = c("category", "index_region", "none"), color_pal = c("#66c2a5", "#fc8d62", "#8da0cb"), out_dir) {
  if (FALSE) {
    tar_load(biomass_indices)
    index_res_df <- index_res_out
    index_res_df <- biomass_indices
    nice_category_names <- "American lobster"
    nice_xlab <- "Year-Season"
    nice_ylab <- "Biomass index (metric tons)"
    color_pal <- NULL
    paneling <- "none"
    date_breaks <- "5 year"
    out_dir <- paste0(res_root, "plots_maps")

    index_res_df <- biomass_indices
    index_scale <- "raw"
    nice_category_names <- nice_category_names
    nice_xlab <- "Year"
    nice_ylab <- "Biomass index (metric tons)"
    paneling <- "none"
    color_pal <- NULL
    out_dir <- here::here("", "Objective 3/Temp_Results")
    
    index_res_df = bio_index_df
    index_scale = "log"
    nice_category_names = "Atlantic cod"
    nice_xlab = "Year-Season"
    nice_ylab = "Biomass index (metric tons)"
    paneling = "none"
    color_pal = NULL
    out_dir = date_dir
  }

  if (paneling == "none") {
    if (!is.null(color_pal)) {
      colors_use <- color_pal
    } else {
      color_pal <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854")
      colors_use <- color_pal[1:length(unique(index_res_df$Index_Region))]
    }

    # Filter based on years to plot
    if (!is.null(year_stop)) {
      index_res_df <- index_res_df %>%
        filter(., Year < year_stop)
    }

    # Date axis
    date_breaks <- seq.Date(from = as.Date(paste(min(index_res_df$Year), "06-15", sep = "-")), to = as.Date(paste(max(index_res_df$Year), "06-15", sep = "-")), by = "year")
    plot_out <- ggplot() +
      geom_errorbar(data = index_res_df, aes(x = Date, ymin = (Index_Estimate - Index_SD), ymax = (Index_Estimate + Index_SD), color = Index_Region, group = Index_Region), alpha = 0.65) +
      geom_point(data = index_res_df, aes(x = Date, y = Index_Estimate, color = Index_Region)) +
      geom_line(data = index_res_df, aes(x = Date, y = Index_Estimate, color = Index_Region)) +
      scale_color_manual(values = colors_use) +
      scale_x_date(breaks = date_breaks, date_labels = "%Y") +
      xlab({{ nice_xlab }}) +
      ylab({{ nice_ylab }}) +
      ggtitle({{ nice_category_names }}) +
      theme_bw() +
      theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  }

  # Save and return the plot
  ggsave(plot_out, file = paste(out_dir, "/Biomass_Index_", index_scale, "_", nice_category_names, ".jpg", sep = ""))
  return(plot_out)
}

plot_vast_projected_index <- function(vast_projections, year_stop = NULL, index_scale, nice_category_names = nice_category_names, climate_scenario = climate_scenario, nice_times = nice_times, region_keep = c("DFO", "NMFS", "GoM", "SNE_and_MAB"), nice_xlab, nice_ylab, paneling = c("category", "index_region", "none"), color_pal = c("#66c2a5", "#fc8d62", "#8da0cb"), out_dir) {
  if (FALSE) {
    tar_load(vast_projection_summ_index)
    vast_projections <- vast_projection_summ_index
    vast_projections <- res_out
    year_stop <- NULL
    index_scale <- "log"
    nice_category_names <- nice_category_names
    climate_scenario <- climate_scenario
    nice_times <- nice_times
    nice_xlab <- "Year-Season"
    nice_ylab <- "Log Biomass index (metric tons)"
    color_pal <- NULL
    paneling <- "none"
    date_breaks <- "5 year"
    out_dir <- paste0(res_root, "plots_maps")

    vast_projections <- res_out
    region_keep<- NULL
    year_stop <- NULL
    index_scale <- "log"
    nice_category_names <- "Atlantic cod"
    climate_scenario <- "SSP5_85_Mean"
    nice_times <- nice_times
    nice_xlab <- "time"
    nice_ylab<- "Log Biomass index (metric tons)"
    color_pal = c("#66c2a5", "#fc8d62", "#8da0cb")
    paneling <- "none"
    date_breaks <- "10 year"
    out_dir = date_dir
  }

  if (paneling == "none") {
    if (!is.null(color_pal)) {
      colors_use <- color_pal
    } else {
      color_pal <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854")
      colors_use <- color_pal[1:length(unique(vast_projections$Region))]
    }

    # Filter based on years to plot
    if (!is.null(year_stop)) {
      index_res_df <- vast_projections %>%
        filter(., Time < year_stop)
    } else {
      index_res_df <- vast_projections
    }

    if (!is.null(region_keep)) {
      index_res_df <- index_res_df %>%
        dplyr::filter(., Region %in% region_keep)
      index_res_df$Region <- factor(index_res_df$Region, levels = region_keep)
    }

    if (index_scale == "log") {
      index_res_df <- index_res_df %>%
        mutate_at(., c("Prob_0.1", "Prob_0.5", "Prob_0.9"), log)
    }

    # Date axis
    date_breaks <- seq.Date(from = as.Date(paste(min(index_res_df$Time), "06-15", sep = "-")), to = as.Date(paste(max(index_res_df$Time), "06-15", sep = "-")), by = date_breaks)
    
    plot_out <- ggplot(data = index_res_df, aes(x = Time, ymin = Prob_0.1, ymax = Prob_0.9, fill = Region)) +
      geom_ribbon(alpha = 0.75) +
      geom_line(data = index_res_df, aes(x = Time, y = Prob_0.5, color = Region)) +
      scale_color_manual(values = colors_use) +
      scale_fill_manual(values = colors_use) +
      # scale_x_date(breaks = date_breaks, date_labels = "%Y") +
      xlab({{ nice_xlab }}) +
      ylab({{ nice_ylab }}) +
      ggtitle({{ nice_category_names }}) +
      theme_bw() +
      theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      facet_wrap(~Region)
  }

  # Save and return the plot
  ggsave(plot_out, file = paste(out_dir, "/Biomass_Index_", index_scale, "_", nice_category_names, "_", climate_scenario, ".jpg", sep = ""))
  return(plot_out)
}

######
## Plot parameter effects...
######
#' @title Adapts package \code{effects}
#'
#' @inheritParams effects::Effect
#' @param which_formula which formula to use e.g., \code{"X1"}
#'
#' @rawNamespace S3method(effects::Effect, fit_model)
#' @export
Effect.fit_model_aja <- function(focal.predictors, mod, which_formula = "X1", pad_values = c(), 
                                 category_to_use, ncat,
                                 ...) {
  # if (FALSE) {
  #   tar_load(vast_fit)
  #   focal.predictors <- c("Depth", "SST_seasonal", "BT_seasonal")
  #   mod <- fit_base
  #   which_formula <- "X1"
  #   xlevels <- 100
  #   pad_values <- c(1)
  # 
  #   covariate_data_full <- mod$effects$covariate_data_full
  #   catchability_data_full <- mod$effects$catchability_data_full
  # }

  # Error checks
  #if (mod$data_list$n_c > 1 & which_formula %in% c("X1", "X2")) {
  #  stop("`Effect.fit_model` is not currently designed for multivariate models using density covariates")
  #}
  if (!all(c("covariate_data_full", "catchability_data_full") %in% ls(.GlobalEnv))) {
    stop("Please load `covariate_data_full` and `catchability_data_full` into global memory")
  }
  if (!requireNamespace("effects")) {
    stop("please install the effects package")
  }
  if (!("effects" %in% names(mod))) {
    stop("`effects` slot not detected in input to `Effects.fit_model`. Please update model using later package version.")
  }

  # Identify formula-specific stuff
  if (which_formula == "X1") {
    formula_orig <- mod$X1_formula
    parname <- "gamma1_cp"
    mod$call <- mod$effects$call_X1
  } else if (which_formula == "X2") {
    formula_orig <- mod$X2_formula
    parname <- "gamma2_cp"
    mod$call <- mod$effects$call_X2
  } else if (which_formula == "Q1") {
    formula_orig <- mod$Q1_formula
    parname <- "lambda1_k"
    mod$call <- mod$effects$call_Q1
  } else if (which_formula == "Q2") {
    formula_orig <- mod$Q2_formula
    parname <- "lambda2_k"
    mod$call <- mod$effects$call_Q2
  } else {
    stop("Check `which_formula` input")
  }
  
  # Adjust for category to use
  # Extract parameters / covariance
  whichnum <- which(names(mod$parameter_estimates$par) == parname)
  whichnum <- whichnum[seq(category_to_use, length(whichnum), ncat)]
  mod$parhat <- mod$parameter_estimates$par[whichnum]
  #mod$parhat <- mod$parhat[seq(category_to_use, length(mod$parhat), ncat)]
  
  if (is.null(mod$parameter_estimates$SD$cov.fixed)) {
    mod$covhat <- array(0, dim = rep(length(mod$parhat), 2))
  } else {
    mod$covhat <- mod$parameter_estimates$SD$cov.fixed[whichnum, whichnum, drop = FALSE]
  }

  # # Fill in values that are mapped off
  # if(parname %in% names(mod$tmb_list$Obj$env$map)){
  #   mod$parhat = mod$parhat[mod$tmb_list$Obj$env$map[[parname]]]
  #   mod$covhat = mod$covhat[mod$tmb_list$Obj$env$map[[parname]], mod$tmb_list$Obj$env$map[[parname]], drop = FALSE]
  #   mod$parhat = ifelse(is.na(mod$parhat), 0, mod$parhat)
  #   mod$covhat = ifelse(is.na(mod$covhat), 0, mod$covhat)
  # }

  # add names
  names(mod$parhat)[] <- parname
  if (length(pad_values) != 0) {
    parhat <- rep(NA, length(mod$parhat) + length(pad_values))
    parhat[setdiff(1:length(parhat), pad_values)] <- mod$parhat
    covhat <- array(NA, dim = dim(mod$covhat) + rep(length(pad_values), 2))
    covhat[setdiff(1:length(parhat), pad_values), setdiff(1:length(parhat), pad_values)] <- mod$covhat
    mod$parhat <- ifelse(is.na(parhat), 0, parhat)
    mod$covhat <- ifelse(is.na(covhat), 0, covhat)
    # parname = c("padded_intercept", parname)
  }
  # rownames(mod$covhat) = colnames(mod$covhat) = names(mod$parhat)

  # Augment stuff
  formula_full <- stats::update.formula(formula_orig, linear_predictor ~ . + 0)
  mod$coefficients <- mod$parhat
  mod$vcov <- mod$covhat
  mod$formula <- formula_full
  mod$family <- stats::gaussian(link = "identity")

  if (FALSE) {
    Tmp <- model.matrix(formula_full, data = fit$effects$catchability_data)
  }

  # Functions for package
  family.fit_model <- function(x, ...) x$family
  vcov.fit_model <- function(x, ...) x$vcov

  # dummy functions to make Effect.default work
  dummyfuns <- list(variance = function(mu) mu, initialize = expression(mustart = y + 0.1), dev.resids = function(...) stats::poisson()$dev.res(...))

  # Replace family (for reasons I don't really understand)
  fam <- mod$family
  for (i in names(dummyfuns)) {
    if (is.null(fam[[i]])) fam[[i]] <- dummyfuns[[i]]
  }

  # allow calculation of effects ...
  if (length(formals(fam$variance)) > 1) {
    warning("overriding variance function for effects: computed variances may be incorrect")
    fam$variance <- dummyfuns$variance
  }

  # Bundle arguments
  args <- list(call = mod$call, coefficients = mod$coefficients, vcov = mod$vcov, family = fam, formula = formula_full)

  # Do call
  effects::Effect.default(focal.predictors, mod, ..., sources = args)
}

get_vast_covariate_effects <- function(vast_fit, params_plot, params_plot_levels, effects_pad_values, nice_category_names, out_dir,
                                       category_to_use, ncat, 
                                       ...) {
  # if (FALSE) {
  #   tar_load(vast_fit)
  #   params_plot <- c(
  #     "Depth", "SST_seasonal", "BT_seasonal",
  #     "SS_seasonal", "BS_seasonal"
  #   )
  #   params_plot_levels <- 100
  #   effects_pad_values <- c(1)
  #   nice_category_names <- nice_category_names
  #   out_dir <- paste0(res_root, "tables")
  # 
  #   vast_fit = mod_comp_res$Fitted_Mod[[1]]
  #   params_plot = c("index")
  #   params_plot_levels = 100
  #   effects_pad_values = c(1)
  #   nice_category_names = "Capelin"
  # }
  assign("covariate_data_full", vast_fit$effects$covariate_data_full[seq(category_to_use, nrow(vast_fit$effects$covariate_data_full), ncat),], 
        envir = .GlobalEnv)
  assign("catchability_data_full", vast_fit$effects$catchability_data_full, 
        envir = .GlobalEnv)
  x1_rescale <- function(x) plogis(x)
  x2_rescale <- function(x) exp(x)
  
  for (i in seq_along(params_plot)) {
    if (any(grepl(params_plot[i], labels(terms(vast_fit$X1_formula))))) {
      pred_dat_temp_X1 <- data.frame(Effect.fit_model_aja(focal.predictors = params_plot[i], 
                                                          mod = vast_fit, 
                                                          which_formula = "X1", 
                                                          xlevels = params_plot_levels, 
                                                          pad_values = effects_pad_values,
                                                          category_to_use = category_to_use,
                                                          ncat = ncat
                                                          )) %>%
        mutate(., Lin_pred = "X1")
    }
    
    if(any(grepl(params_plot[i], labels(terms(vast_fit$X2_formula))))){
      pred_dat_temp_X2 <- data.frame(Effect.fit_model_aja(focal.predictors = params_plot[i], 
                                                          mod = vast_fit, which_formula = "X2", 
                                                          xlevels = params_plot_levels, 
                                                          pad_values = effects_pad_values,
                                                          category_to_use = category_to_use,
                                                          ncat = ncat
                                                          )) %>% 
      mutate(., Lin_pred = "X2")
    }
    
    if(exists("pred_dat_temp_X1") & !exists("pred_dat_temp_X2")){
      pred_dat_out_temp <- pred_dat_temp_X1
      rm(pred_dat_temp_X1)
    }
    if(!exists("pred_dat_temp_X1") & exists("pred_dat_temp_X2")){
      pred_dat_out_temp <- pred_dat_temp_X2
      rm(pred_dat_temp_X2)
    }
    if(exists("pred_dat_temp_X1") & exists("pred_dat_temp_X2")) {
      pred_dat_out_temp <- bind_rows(pred_dat_temp_X1, pred_dat_temp_X2)
      rm(pred_dat_temp_X1)
      rm(pred_dat_temp_X2)
    }
    
    if (i == 1) {
      pred_dat_out <- pred_dat_out_temp
    } else {
      pred_dat_out <- bind_rows(pred_dat_out, pred_dat_out_temp)
    }
  }
  pred_dat_out <- pred_dat_out %>%
    pivot_longer(., !c(fit, se, lower, upper, Lin_pred), names_to = "Covariate", values_to = "Value") %>%
    drop_na()
  saveRDS(pred_dat_out, file = paste(out_dir, "/", nice_category_names, "_covariate_effects.rds", sep = ""))
  return(pred_dat_out)
}


plot_vast_covariate_effects <- function (vast_covariate_effects, vast_fit, nice_category_names, 
    out_dir, ...) 
{
    if (FALSE) {
        tar_load(vast_covariate_effects_Offshore_hake_full_69)
        vast_covariate_effects<- vast_covariate_effects_Offshore_hake_full_69
        tar_load(vast_fit_Offshore_hake_full_69)
        vast_fit<- vast_fit_Offshore_hake_full_69
        nice_category_names <- "American lobster"
        plot_rows <- 2
        res_root <- "/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/targets_output/"
        out_dir <- paste0(res_root, "tables")
    }
    names_stay <- c("fit", "se", "lower", "upper", "Lin_pred")
    # vast_cov_eff_l <- vast_covariate_effects %>% pivot_longer(.,
    #     names_to = "Variable", values_to = "Covariate_Value",
    #     -{
    #         {
    #             names_stay
    #         }
    #     }) %>% drop_na(Covariate_Value)
    vast_cov_eff_l <- vast_covariate_effects %>%
        drop_na(Value)
    # ylim_dat <- vast_cov_eff_l %>% group_by(., Lin_pred, Variable) %>%
    #     summarize(., Min = min(lower, na.rm = TRUE), Max = max(upper,
    #         na.rm = TRUE))
    ylim_dat <- vast_cov_eff_l %>%
        group_by(., Lin_pred, Covariate) %>%
        summarize(., Min = min(lower, na.rm = TRUE), Max = max(upper, na.rm = TRUE))
    
    plot_out <- ggplot() +
        geom_ribbon(data = vast_cov_eff_l, aes(x = Value, ymin = lower, ymax = upper), fill = "#bdbdbd") +
        geom_line(data = vast_cov_eff_l, aes(x = Value, y = fit)) +
        xlab("Scaled covariate value") +
        ylab("Linear predictor fitted value") +
        facet_grid(Lin_pred ~ Covariate, scales = "free") +
        theme_bw() +
        theme(strip.background = element_blank())
    
    names_keep <- unique(vast_cov_eff_l$Covariate)
    samp_dat <- vast_fit$covariate_data %>% dplyr::select({
        {
            names_keep
        }
    }) %>% gather(., "Covariate", "Value")
    plot_out2 <- plot_out +
        geom_rug(data = samp_dat, aes(x = Value))
    ggsave(plot_out2, filename = paste(out_dir, "/", nice_category_names, 
        "_covariate_effects.jpg", sep = ""))
    return(plot_out2)
}
######
## Plot samples, knots and mesh
######
vast_plot_design <- function(vast_fit, land, spat_grid, xlim = c(-80, -55), ylim = c(35, 50), land_color = "#f0f0f0", out_dir) {
  if (FALSE) {
    vast_fit <- readRDS(paste(res_root, "mod_fits/American_lobster_STnoRW_fitted_vast.rds", sep = ""))
    tar_load(land_sf)
    spat_grid <- "~/GitHub/sdm_workflow/scratch/aja/TargetsSDM/data/predict/predict_stack_SST_seasonal_mean.grd"
    land <- land_sf
    xlim <- c(-80, -55)
    ylim <- c(35, 50)
    land_color <- "#f0f0f0"

    vast_fit <- vast_fitted
    land <- land_use
    spat_grid <- spat_grid
    xlim <- xlim_use
    ylim <- ylim_use
    land_color <- "#f0f0f0"
    out_dir <- main_dir
  }

  # Read in raster
  spat_grid <- rotate(raster::stack(spat_grid)[[1]])

  # Intensity surface of sample locations and then a plot of the knot locations/mesh over the top?
  samp_dat <- vast_fit$data_frame %>%
    distinct(., Lon_i, Lat_i, .keep_all = TRUE) %>%
    st_as_sf(., coords = c("Lon_i", "Lat_i"), remove = FALSE, crs = st_crs(land))

  cell_samps <- table(cellFromXY(spat_grid, data.frame("x" = samp_dat$Lon_i, "y" = samp_dat$Lat_i)))

  # Put back into raster...
  spat_grid[] <- 0
  spat_grid[as.numeric(names(cell_samps))] <- cell_samps
  spat_grid_plot <- as.data.frame(spat_grid, xy = TRUE)
  names(spat_grid_plot)[3] <- "Samples"
  spat_grid_plot$Samples <- ifelse(spat_grid_plot$Samples == 0, NA, spat_grid_plot$Samples)

  tow_samps <- ggplot() +
    geom_tile(data = spat_grid_plot, aes(x = x, y = y, fill = Samples)) +
    scale_fill_gradient2(name = "Tow samples", low = "#bdbdbd", high = "#525252", na.value = "white") +
    geom_sf(data = land, fill = land_color, lwd = 0.2, na.rm = TRUE) +
    coord_sf(xlim = xlim, ylim = ylim, expand = 0) +
    theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "pt")) +
    ggtitle("Tow samples")

  # Knots and mesh...
  # Getting spatial information
  spat_data <- vast_fit$extrapolation_list
  extrap_grid <- data.frame("Lon" = as.numeric(spat_data$Data_Extrap$Lon), "Lat" = as.numeric(spat_data$Data_Extrap$Lat)) %>%
    distinct(., Lon, Lat)

  tow_samps_grid <- tow_samps +
    geom_point(data = extrap_grid, aes(x = Lon, y = Lat), fill = "#41ab5d", pch = 21, size = 0.75) +
    ggtitle("VAST spatial extrapolation grid")

  # Get mesh as sf
  mesh_sf <- vast_mesh_to_sf(vast_fit, crs_transform = "+proj=longlat +datum=WGS84 +no_defs")$triangles
  tow_samps_mesh <- tow_samps +
    geom_sf(data = land, fill = land_color, lwd = 0.2, na.rm = TRUE) +
    geom_sf(data = mesh_sf, fill = NA, color = "#41ab5d") +
    coord_sf(xlim = xlim, ylim = ylim, expand = 0) +
    theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "pt")) +
    ggtitle("INLA Mesh")

  # Plot em together
  plot_out <- tow_samps + tow_samps_grid + tow_samps_mesh

  # Save it
  ggsave(plot_out, file = paste(out_dir, "/", "samples_grid_knots_plot.jpg", sep = ""), height = 8, width = 11)
  return(plot_out)
}

#####
## Plot covariate values
#####

plot_spattemp_cov_ts <- function(predict_covariates_stack_agg, summarize = "seasonal", ensemble_stat = "mean", all_tows_with_all_covs, regions, land, out_dir) {
  if (FALSE) {
    tar_load(predict_covariates_stack_agg_out)
    predict_covariates_stack_agg <- predict_covariates_stack_agg_out
    summarize <- "seasonal"
    ensemble_stat <- "mean"
    tar_load(all_tows_with_all_covs)
    tar_load(land_sf)
    land <- land_sf
    tar_load(index_shapefiles)
    out_dir <- paste0(res_root, "plots_maps")
  }

  # Get raster stack covariate files
  rast_files_load <- list.files(predict_covariates_stack_agg, pattern = paste0(summarize, "_", ensemble_stat, ".grd$"), full.names = TRUE)

  # Get variable names
  cov_names_full <- list.files(predict_covariates_stack_agg, pattern = paste0(summarize, "_", ensemble_stat, ".grd$"), full.names = FALSE)
  predict_covs_names <- gsub(paste("_", ensemble_stat, ".grd", sep = ""), "", gsub("predict_stack_", "", cov_names_full))

  # Loop through
  for (i in seq_along(rast_files_load)) {
    # Get variable names
    cov_names_full <- list.files(predict_covariates_stack_agg, pattern = paste0(summarize, "_", ensemble_stat, ".grd$"), full.names = FALSE)[i]
    predict_covs_names <- gsub(paste("_", ensemble_stat, ".grd", sep = ""), "", gsub("predict_stack_", "", cov_names_full))

    # Prediction values
    spattemp_summs <- data.frame(raster::extract(raster::rotate(raster::stack(rast_files_load[i])), index_shapefiles, fun = mean))
    spattemp_summs$Region <- factor(unique(as.character(index_shapefiles$Region)), levels = c("DFO", "Georges_Bank", "GoM", "Mid_Atlantic_Bight", "NMFS_and_DFO", "NMFS", "SNE_and_MAB", "Southern_New_England", "Scotian_Shelf"))
    spattemp_summs <- spattemp_summs %>%
      drop_na(., Region) %>%
      filter(., Region %in% c("NMFS_and_DFO", "DFO", "GoM", "SNE_and_MAB"))

    # Gather
    spattemp_summs_df <- spattemp_summs %>%
      pivot_longer(., names_to = "Time", values_to = "Value", -Region)

    # Formatting Time
    spattemp_summs_df <- spattemp_summs_df %>%
      mutate(., Date = gsub("X", "", gsub("[.]", "-", Time)))
    spattemp_summs_df$Date <- as.Date(paste(as.numeric(gsub("([0-9]+).*$", "\\1", spattemp_summs_df$Date)), ifelse(grepl("Spring", spattemp_summs_df$Date), "-04-15", ifelse(grepl("Summer", spattemp_summs_df$Date), "-07-15", ifelse(grepl("Winter", spattemp_summs_df$Date), "-12-15", "-10-15"))), sep = ""))

    # Data values
    cov_dat <- all_tows_with_all_covs %>%
      dplyr::select(., Season_Match, DECDEG_BEGLON, DECDEG_BEGLAT, {{ predict_covs_names }})
    cov_dat$Date <- as.Date(paste(as.numeric(gsub("([0-9]+).*$", "\\1", cov_dat$Season_Match)), ifelse(grepl("Spring", cov_dat$Season_Match), "-04-15", ifelse(grepl("Summer", cov_dat$Season_Match), "-07-15", ifelse(grepl("Winter", cov_dat$Season_Match), "-12-15", "-10-15"))), sep = ""))

    # Get summary by region...
    cov_dat <- cov_dat %>%
      st_as_sf(., coords = c("DECDEG_BEGLON", "DECDEG_BEGLAT"), crs = st_crs(index_shapefiles), remove = FALSE) %>%
      st_join(., index_shapefiles, join = st_within) %>%
      st_drop_geometry()

    cov_dat_plot <- cov_dat %>%
      group_by(., Date, Region) %>%
      summarize_at(., .vars = {{ predict_covs_names }}, .funs = mean, na.rm = TRUE)

    cov_dat_plot$Region <- factor(cov_dat_plot$Region, levels = c("DFO", "Georges_Bank", "GoM", "Mid_Atlantic_Bight", "NMFS_and_DFO", "NMFS", "SNE_and_MAB", "Southern_New_England", "Scotian_Shelf"))


    cov_dat_plot <- cov_dat_plot %>%
      drop_na(., c({{ predict_covs_names }}, Region)) %>%
      filter(., Region %in% c("NMFS_and_DFO", "DFO", "GoM", "SNE_and_MAB"))

    # Plot
    if (predict_covs_names == "Depth") {
      plot_out <- ggplot() +
        geom_line(data = spattemp_summs_df, aes(x = Date, y = Value, color = Region)) +
        scale_color_manual(name = "Region", values = c("#1b9e77", "#d95f02", "#7570b3", "#666666")) +
        geom_point(data = cov_dat_plot, aes(x = Date, y = Depth), fill = "black", pch = 21, alpha = 0.2) +
        facet_wrap(~Region, nrow = 2) +
        theme_bw()
    }

    if (predict_covs_names == "BS_seasonal") {
      plot_out <- ggplot() +
        geom_line(data = spattemp_summs_df, aes(x = Date, y = Value, color = Region)) +
        scale_color_manual(name = "Region", values = c("#1b9e77", "#d95f02", "#7570b3", "#666666")) +
        geom_point(data = cov_dat_plot, aes(x = Date, y = BS_seasonal), fill = "black", pch = 21, alpha = 0.2) +
        facet_wrap(~Region, nrow = 2) +
        theme_bw()
    }

    if (predict_covs_names == "SS_seasonal") {
      plot_out <- ggplot() +
        geom_line(data = spattemp_summs_df, aes(x = Date, y = Value, color = Region)) +
        scale_color_manual(name = "Region", values = c("#1b9e77", "#d95f02", "#7570b3", "#666666")) +
        geom_point(data = cov_dat_plot, aes(x = Date, y = SS_seasonal), fill = "black", pch = 21, alpha = 0.2) +
        facet_wrap(~Region, nrow = 2) +
        theme_bw()
    }

    if (predict_covs_names == "BT_seasonal") {
      plot_out <- ggplot() +
        geom_line(data = spattemp_summs_df, aes(x = Date, y = Value, color = Region)) +
        scale_color_manual(name = "Region", values = c("#1b9e77", "#d95f02", "#7570b3", "#666666")) +
        geom_point(data = cov_dat_plot, aes(x = Date, y = BT_seasonal), fill = "black", pch = 21, alpha = 0.2) +
        facet_wrap(~Region, nrow = 2) +
        theme_bw()
    }

    if (predict_covs_names == "SST_seasonal") {
      plot_out <- ggplot() +
        geom_line(data = spattemp_summs_df, aes(x = Date, y = Value, color = Region)) +
        scale_color_manual(name = "Region", values = c("#1b9e77", "#d95f02", "#7570b3", "#666666")) +
        geom_point(data = cov_dat_plot, aes(x = Date, y = SST_seasonal), fill = "black", pch = 21, alpha = 0.2) +
        facet_wrap(~Region, nrow = 2) +
        theme_bw()
    }

    ggsave(paste(out_dir, "/", predict_covs_names, "_covariate_plot.jpg", sep = ""), plot_out)
  }
}

plot_scaled_cov_ts <- function(vast_pred_df_post_fit, all_tows_with_all_covs, covs = c("Depth", "SST_seasonal", "BT_seasonal", "BS_seasonal", "SS_seasonal"), regions, land, out_dir) {
  if (FALSE) {
    tar_load(vast_pred_df_post_fit)
    tar_load(all_tows_with_all_covs_rescale)
    covs <- c("SST_seasonal", "BT_seasonal", "BS_seasonal", "SS_seasonal")
    tar_load(land_sf)
    land <- land_sf
    tar_load(index_shapefiles)
    out_dir <- paste0(res_root, "plots_maps")
  }

  # Spatial to get summaries by regions
  pred_sf <- st_as_sf(vast_pred_df_post_fit, coords = c("Lon", "Lat"), crs = 4326)
  fit_sf <- st_as_sf(all_tows_with_all_covs_rescale, coords = c("DECDEG_BEGLON", "DECDEG_BEGLAT"), crs = 4326)

  # Join with region shapefiles
  pred_df <- pred_sf %>%
    st_join(., index_shapefiles, join = st_within) %>%
    st_drop_geometry() %>%
    filter(., Region %in% c("NMFS_and_DFO", "DFO", "GoM", "SNE_and_MAB")) %>%
    group_by(., Date, Region) %>%
    summarize_at(., .vars = {{ covs }}, .funs = mean, na.rm = TRUE)
  pred_df$Region <- factor(pred_df$Region, levels = c("NMFS_and_DFO", "DFO", "GoM", "SNE_and_MAB"))

  fit_df <- fit_sf %>%
    st_join(., index_shapefiles, join = st_within) %>%
    st_drop_geometry() %>%
    filter(., Region %in% c("NMFS_and_DFO", "DFO", "GoM", "SNE_and_MAB")) %>%
    filter(., DATE >= "1985-01-01")
  fit_df$DATE <- as.Date(paste(as.numeric(gsub("([0-9]+).*$", "\\1", fit_df$Season_Match)), ifelse(grepl("Spring", fit_df$Season_Match), "-04-15", ifelse(grepl("Summer", fit_df$Season_Match), "-07-15", ifelse(grepl("Winter", fit_df$Season_Match), "-12-15", "-10-15"))), sep = ""))

  fit_df <- fit_df %>%
    group_by(., DATE, Region) %>%
    summarize_at(., .vars = {{ covs }}, .funs = mean, na.rm = TRUE)
  fit_df$Region <- factor(fit_df$Region, levels = c("NMFS_and_DFO", "DFO", "GoM", "SNE_and_MAB"))


  # Loop through
  for (i in seq_along(covs)) {
    cov_plot <- covs[i]

    pred_plot <- pred_df %>%
      dplyr::select(., Date, {{ cov_plot }}, Region) %>%
      mutate(.,
        "Region" = factor(Region, levels(pred_df$Region)),
        "Date" = as.Date(Date)
      )


    fit_plot <- fit_df %>%
      dplyr::select(., DATE, {{ cov_plot }}, Region) %>%
      rename_at(vars(DATE), str_to_sentence) %>%
      mutate(.,
        "Region" = factor(Region, levels(pred_df$Region)),
        "Date" = as.Date(Date)
      )


    # Plot
    if (cov_plot == "BS_seasonal") {
      plot_out <- ggplot() +
        geom_line(data = pred_plot, aes(x = Date, y = BS_seasonal, color = Region)) +
        scale_color_manual(name = "Region", values = c("#1b9e77", "#d95f02", "#7570b3", "#666666")) +
        geom_point(data = fit_plot, aes(x = Date, y = BS_seasonal), fill = "black", pch = 21, alpha = 0.2) +
        facet_wrap(~Region, nrow = 2) +
        theme_bw()
    }

    if (cov_plot == "SS_seasonal") {
      plot_out <- ggplot() +
        geom_line(data = pred_plot, aes(x = Date, y = SS_seasonal, color = Region)) +
        scale_color_manual(name = "Region", values = c("#1b9e77", "#d95f02", "#7570b3", "#666666")) +
        geom_point(data = fit_plot, aes(x = Date, y = SS_seasonal), fill = "black", pch = 21, alpha = 0.2) +
        facet_wrap(~Region, nrow = 2) +
        theme_bw()
    }

    if (cov_plot == "BT_seasonal") {
      plot_out <- ggplot() +
        geom_line(data = pred_plot, aes(x = Date, y = BT_seasonal, color = Region)) +
        scale_color_manual(name = "Region", values = c("#1b9e77", "#d95f02", "#7570b3", "#666666")) +
        geom_point(data = fit_plot, aes(x = Date, y = BT_seasonal), fill = "black", pch = 21, alpha = 0.2) +
        facet_wrap(~Region, nrow = 2) +
        theme_bw()
    }

    if (cov_plot == "SST_seasonal") {
      plot_out <- ggplot() +
        geom_line(data = pred_plot, aes(x = Date, y = SST_seasonal, color = Region)) +
        scale_color_manual(name = "Region", values = c("#1b9e77", "#d95f02", "#7570b3", "#666666")) +
        geom_point(data = fit_plot, aes(x = Date, y = SST_seasonal), fill = "black", pch = 21, alpha = 0.2) +
        facet_wrap(~Region, nrow = 2) +
        theme_bw()
    }

    ggsave(paste(out_dir, "/", cov_plot, "_scaled_covariate_plot.jpg", sep = ""), plot_out)
  }
}


#####
## VAST inla mesh to sf object
#####
#' @title Convert VAST INLA mesh to sf object
#'
#' @description Convert inla.mesh to sp objects, totally taken from David Keith here https://github.com/Dave-Keith/Paper_2_SDMs/blob/master/mesh_build_example/convert_inla_mesh_to_sf.R and Finn Lindgren here
# # https://groups.google.com/forum/#!topic/r-inla-discussion-group/z1n1exlZrKM
#'
#' @param vast_fit A fitted VAST model
#' @param crs_transform Optional crs to transform mesh into
#' @return A list with \code{sp} objects for triangles and vertices:
# \describe{
#   \item{triangles}{\code{SpatialPolygonsDataFrame} object with the triangles in
#   the same order as in the original mesh, but each triangle looping through
#   the vertices in clockwise order (\code{sp} standard) instead of
#   counterclockwise order (\code{inla.mesh} standard). The \code{data.frame}
#   contains the vertex indices for each triangle, which is needed to link to
#   functions defined on the vertices of the triangulation.
#   \item{vertices}{\code{SpatialPoints} object with the vertex coordinates,
#   in the same order as in the original mesh.}
# }
#' @export
#
vast_mesh_to_sf <- function(vast_fit, crs_transform = "+proj=longlat +datum=WGS84 +no_defs") {
  if (FALSE) {
    tar_load(vast_fit)
    crs_transform <- "+proj=longlat +datum=WGS84 +no_defs"
  }
  require(sp) || stop("Install sp, else thine code shan't work for thee")
  require(sf) || stop("Install sf or this code will be a mess")
  require(INLA) || stop("You need the R-INLA package for this, note that it's not crantastic...
                        install.packages('INLA', repos=c(getOption('repos'), INLA='https://inla.r-inla-download.org/R/stable'), dep=TRUE)")

  # Get the extrapolation mesh information from the vast_fitted object
  mesh <- vast_fit$spatial_list$MeshList$anisotropic_mesh
  mesh["crs"] <- vast_fit$extrapolation_list$projargs

  # Grab the CRS if it exists, NA is fine (NULL spits a warning, but is also fine)
  crs <- sp::CRS(mesh$crs)

  # Make sure the CRS isn't a geocentric one, which is won't be if yo look up geocentric..
  # isgeocentric <- identical(inla.as.list.CRS(crs)[["proj"]], "geocent")
  isgeocentric <- inla.crs_is_geocent(mesh$crs)
  # Look up geo-centric coordinate systems, nothing we'll need to worry about, but stop if so
  if (isgeocentric || (mesh$manifold == "S2")) {
    stop(paste0(
      "'sp and sf' don't support storing polygons in geocentric coordinates.\n",
      "Convert to a map projection with inla.spTransform() before calling inla.mesh2sf()."
    ))
  }
  # This pulls out from the mesh the triangles as polygons, this was the piece I couldn't figure out.
  triangles <- SpatialPolygonsDataFrame(
    Sr = SpatialPolygons(
      lapply(
        1:nrow(mesh$graph$tv),
        function(x) {
          tv <- mesh$graph$tv[x, , drop = TRUE]
          Polygons(list(Polygon(mesh$loc[tv[c(1, 3, 2, 1)], 1:2, drop = FALSE])), ID = x)
        }
      ),
      proj4string = crs
    ),
    data = as.data.frame(mesh$graph$tv[, c(1, 3, 2), drop = FALSE]),
    match.ID = FALSE
  )

  # This one is easy, just grab the vertices (points)
  vertices <- SpatialPoints(mesh$loc[, 1:2, drop = FALSE], proj4string = crs)

  # Make these sf objects
  triangles <- st_as_sf(triangles)
  vertices <- st_as_sf(vertices)

  # Transform?
  if (!is.null(crs_transform)) {
    triangles <- st_transform(triangles, crs = crs_transform)
    vertices <- st_transform(vertices, crs = crs_transform)
  }

  # Add your output list.
  return_sf <- list(triangles = triangles, vertices = vertices)
  return(return_sf)
}

#' @title Plot VAST model spatial and spatio-temporal surfaces
#'
#' @description Creates either a panel plot or a gif of VAST model spatial or spatio-temporal parameter surfaces or derived quantities
#'
#' @param vast_fit = A VAST `fit_model` object
#' @param manaul_pred_df = Data frame with location, time and predicted value information OR NULL if using vast_fit
#' @param spatial_var = An estimated spatial coefficient or predicted value or NULL. Currently works for `D_gct`, `R1_gct`, `R2_gct`, `P1_gct`, `P2_gct`, `Omega1_gc`, `Omega2_gc`, `Epsilon1_gct`, `Epsilon2_gct`.
#' @param nice_category_names = A character string to define species/model run and used in naming the output prediction file.
#' @param pred_label = Explain
#' @param all_times = A vector of all of the unique time steps available from the VAST fitted model
#' @param plot_times = Either NULL to make a plot for each time in `all_times` or a vector of all of the times to plot, which must be a subset of `all_times`
#' @param land_sf = Land sf object
#' @param xlim = A two element vector with the min and max longitudes
#' @param ylim = A two element vector with the min and max latitudes
#' @param panel_or_gif = A character string of either "panel" or "gif" indicating how the multiple plots across time steps should be displayed
#' @param out_dir = Output directory to save the panel plot or gif
#'
#' @return A VAST fit_model object, with the inputs and and outputs, including parameter estimates, extrapolation gid info, spatial list info, data info, and TMB info.
#'
#' @export

vast_fit_plot_spatial <- function(vast_fit, manual_pred_df, pred_grid, spatial_var, nice_category_names, pred_label, mask, all_times = all_times, plot_times = NULL, land_sf, xlim, ylim, lab_lat = 33.75, lab_lon = -67.5, panel_or_gif = "gif", out_dir, land_color = "#d9d9d9", panel_cols = NULL, panel_rows = NULL, ...) {
  if (FALSE) {
    tar_load(vast_fit)
    template <- raster("~/GitHub/sdm_workflow/scratch/aja/TargetsSDM/data/supporting/HighResTemplate.grd")
    tar_load(vast_seasonal_data)
    all_times <- as.character(levels(vast_seasonal_data$VAST_YEAR_SEASON))
    plot_times <- NULL
    tar_load(land_sf)
    tar_load(region_shapefile)
    mask <- region_shapefile
    land_color <- "#d9d9d9"
    res_data_path <- "~/Box/RES_Data/"
    xlim <- c(-78.5, -56)
    ylim <- c(35, 48)
    panel_or_gif <- "gif"
    panel_cols <- NULL
    panel_rows <- NULL
    grid_space_use <- 10
    lab_lat <- 36
    lab_lon <- -60
    spatial_var <- "D_gct"

    # Manual
    vast_fit <- readRDS("/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/targets_output/mod_fits/American_lobster_STnoRW_fitted_vast.rds")
    manual_pred_df <- read.csv(paste0(res_root, "prediction_df/American_lobster_SSP5_85_mean_projections.csv"))
    tar_load(region_shapefile)
    mask <- region_shapefile
    nice_category_names <- "American_lobster"
    pred_label <- "SSP5_85_mean"
    all_times <- unique(manual_pred_df$Real_Date)
    plot_times <- NULL
    tar_load(land_sf)
    xlim <- c(-78.5, -55)
    ylim <- c(34.5, 48.25)
    panel_or_gif <- "gif"
    panel_cols <- NULL
    panel_rows <- NULL
    grid_space_utm <- 25000
    lab_lat <- 36
    lab_lon <- -67.5
    out_dir <- paste0(res_root, "plots_maps")
    land_color <- "#d9d9d9"

    vast_fit = fit
    manual_pred_df = NULL
    spatial_var = "D_gct"
    nice_category_names = "Atlantic cod"
    pred_label = "Obs"
    mask = mask
    all_times = all_times
    plot_times = NULL
    land_sf = land_sf
    xlim = c(-78.5, -56)
    ylim = c(35, 48)
    lab_lat = 36
    lab_lon = -60
    panel_or_gif = "gif"
    out_dir = date_dir
    land_color = "#d9d9d9"


  }

  if (is.null(manual_pred_df)) {
    # Plotting at spatial knots...
    # First check the spatial_var, only a certain subset are being used...
    if (!spatial_var %in% c("D_gct", "R1_gct", "R2_gct", "P1_gct", "P2_gct", "Omega1_gc", "Omega2_gc", "Epsilon1_gct", "Epsilon2_gct")) {
      stop(print("Check `spatial_var` input. Currently must be one of `D_gct`, `R1_gct`, `R2_gct`, `P1_gct`, `P2_gct`, `Omega1_gc`, `Omega2_gc`, `Epsilon1_gct`, `Epsilon2_gct`."))
    }

    # Getting prediction array
    pred_array <- vast_fit$Report[[{{ spatial_var }}]]
    if (spatial_var == "D_gct") {
      pred_array <- log(pred_array)
      pred_array<- drop_units(pred_array)
    }

    # Getting time info
    if (!is.null(plot_times)) {
      plot_times <- all_times[which(all_times) %in% plot_times]
    } else {
      plot_times <- all_times
    }

    # Getting spatial information
    if (vast_fit$spatial_list$fine_scale == TRUE) {
      spat_data <- vast_fit$extrapolation_list
      locs <- data.frame(spat_data$Data_Extrap[which(spat_data$Data_Extrap[, "Include"] > 0), c("Lon", "Lat")]) %>%
        distinct()
    } else {
      spat_data <- vast_fit$spatial_list
      locs <- spat_data$latlon_s[1:spat_data$n_x, ]
    }
    
    
    CRS_orig <- sp::CRS("+proj=longlat")
    CRS_proj <- sp::CRS(spat_data$projargs)
    land_sf <- st_crop(land_sf, xmin = xlim[1], ymin = ylim[1], xmax = xlim[2], ymax = ylim[2])

    # Looping through...
    rasts_out <- vector("list", dim(pred_array)[length(dim(pred_array))])
    rasts_range <- pred_array
    # rast_lims_min <- ifelse(spatial_var %in% c("D_gct", "R1_gct", "R2_gct", "P1_gct", "P2_gct"), 0, min(rasts_range))
    # rast_lims_max <- ifelse(spatial_var %in% c("D_gct", "R1_gct", "R2_gct", "P1_gct", "P2_gct"), round(max(rasts_range) + 0.0000001, 2), max(rasts_range))
    rast_lims <- c(min(rasts_range), max(rasts_range))

    if (length(dim(pred_array)) == 2) {
      data_df <- data.frame(locs, z = pred_array)

      # Interpolation
      pred_df <- na.omit(data.frame("x" = data_df$Lon, "y" = data_df$Lat, "layer" = data_df$z))
      pred_df_interp <- interp(pred_df[, 1], pred_df[, 2], pred_df[, 3],
        duplicate = "mean", extrap = TRUE,
        xo = seq(-87.99457, -57.4307, length = x_dim_length),
        yo = seq(22.27352, 48.11657, length = y_dim_length)
      )
      pred_df_interp_final <- data.frame(expand.grid(x = pred_df_interp$x, y = pred_df_interp$y), z = c(round(pred_df_interp$z, 2)))
      pred_sp <- st_as_sf(pred_df_interp_final, coords = c("x", "y"), crs = CRS_orig)

      pred_df_temp <- pred_sp[which(st_intersects(pred_sp, mask, sparse = FALSE) == TRUE), ]
      coords_keep <- as.data.frame(st_coordinates(pred_df_temp))
      row.names(coords_keep) <- NULL
      pred_df_use <- data.frame(cbind(coords_keep, "z" = as.numeric(pred_df_temp$z)))
      names(pred_df_use) <- c("x", "y", "z")

      plot_out <- ggplot() +
        geom_tile(data = pred_df_use, aes(x = x, y = y, fill = z)) +
        scale_fill_viridis_c(name = spatial_var, option = "viridis", na.value = "transparent", limits = rast_lims) +
        annotate("text", x = lab_lon, y = lab_lat, label = spatial_var) +
        geom_sf(data = land_sf, fill = land_color, lwd = 0.2, na.rm = TRUE) +
        coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
        theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "pt"))

      ggsave(filename = paste(out_dir, "/", nice_category_names, "_", climate_scenario, "_", spatial_var, ".png", sep = ""), plot_out, width = 11, height = 8, units = "in")
      return(plot_out)
    } else {
      for (tI in 1:dim(pred_array)[3]) {
        data_df <- data.frame(locs, z = pred_array[, 1, tI])

        # Interpolation
        # pred_df <- na.omit(data.frame("x" = data_df$Lon, "y" = data_df$Lat, "layer" = data_df$z))
        # pred_df_interp <- interp(pred_df[, 1], pred_df[, 2], pred_df[, 3],
        #   duplicate = "mean", extrap = TRUE,
        #   xo = seq(-87.99457, -57.4307, length = 113),
        #   yo = seq(22.27352, 48.11657, length = 133)
        # )
        # pred_df_interp_final <- data.frame(expand.grid(x = pred_df_interp$x, y = pred_df_interp$y), z = c(round(pred_df_interp$z, 2)))
        # pred_sp <- st_as_sf(pred_df_interp_final, coords = c("x", "y"), crs = CRS_orig)

        # pred_df_temp <- pred_sp[which(st_intersects(pred_sp, mask, sparse = FALSE) == TRUE), ]
        # coords_keep <- as.data.frame(st_coordinates(pred_df_temp))
        # row.names(coords_keep) <- NULL
        # pred_df_use <- data.frame(cbind(coords_keep, "z" = as.numeric(pred_df_temp$z)))
        # names(pred_df_use) <- c("x", "y", "z")

        krig_mod <- Krig(data.frame("x" = data_df$Lon, "y" = data_df$Lat), data_df$z)
        pred_df_interp <- as.data.frame(interpolate(pred_grid, krig_mod), xy = TRUE)
        pred_sp <- st_as_sf(pred_df_interp, coords = c("x", "y"), crs = CRS_orig)

        pred_df_temp <- pred_sp[which(st_intersects(pred_sp, mask, sparse = FALSE) == TRUE), ]
        pred_df_use <- data.frame(st_coordinates(pred_df_temp), "z" = as.numeric(pred_df_temp$layer))
        names(pred_df_use) <- c("x", "y", "z")

        time_plot_use <- plot_times[tI]

        rasts_out[[tI]] <- ggplot() +
          geom_tile(data = pred_df_use, aes(x = x, y = y, fill = z)) +
          scale_fill_viridis_c(name = spatial_var, option = "viridis", na.value = "transparent", limits = rast_lims) +
          annotate("text", x = lab_lon, y = lab_lat, label = time_plot_use) +
          geom_sf(data = land_sf, fill = land_color, lwd = 0.2, na.rm = TRUE) +
          coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
          theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "pt"))
      }
      if (panel_or_gif == "panel") {
        # Panel plot
        all_plot <- wrap_plots(rasts_out, ncol = panel_cols, nrow = panel_rows, guides = "collect", theme(plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "pt")))
        ggsave(filename = paste0(out_dir, "/", nice_category_names, "_", pred_label, "_", spatial_var, ".png"), all_plot, width = 11, height = 8, units = "in")
        return(all_plot)
      } else {
        # Make a gif
        plot_loop_func <- function(plot_list) {
          for (i in seq_along(plot_list)) {
            plot_use <- plot_list[[i]]
            print(plot_use)
          }
        }
        invisible(save_gif(plot_loop_func(rasts_out), paste0(out_dir, "/", nice_category_names, "_", pred_label, "_", spatial_var, ".gif"), delay = 0.75, progress = FALSE))
      }
    }
  } else {
    # Using manual post fit predictions data frame
    # Getting time info
    if (!is.null(plot_times)) {
      plot_times <- all_times[which(all_times) %in% plot_times]
    } else {
      plot_times <- all_times
    }

    # Getting spatial information
    land_sf <- st_transform(land_sf, crs = 32619)

    # Converting limits...
    plot_coords <- st_sfc(st_point(c(xlim[1], ylim[1])), st_point(c(xlim[2], ylim[2])), crs = 4326) %>%
      st_transform(., crs = 32619) %>%
      st_coordinates(.)

    lab_coords <- st_sfc(st_point(c(lab_lon[1], lab_lat[1])), crs = 4326) %>%
      st_transform(., crs = 32619) %>%
      st_coordinates(.)

    # We want a prediction grid...
    mask_utm <- vast_mesh_to_sf(vast_fit)$triangles %>%
      st_union() %>%
      st_transform(., crs = 32619)
    pred_grid <- mask_utm %>%
      st_make_grid(., cellsize = grid_space_utm, what = "polygons", square = TRUE)

    pred_grid2 <- pred_grid %>%
      st_intersection(., mask_utm)

    # Raster storage and limits
    rasts_out <- vector("list", length(plot_times))
    rast_lims <- c(0, max(log((625 * manual_pred_df$Dens) + 1)))

    for (tI in seq_along(plot_times)) {
      time_plot_use <- plot_times[tI]

      plot_label <- paste(format(as.Date(time_plot_use), "%Y"), ifelse(grepl("03-16", time_plot_use), "SPRING", ifelse(grepl("07-16", time_plot_use), "SUMMER", "FALL")), sep = " ")

      data_df <- manual_pred_df %>%
        dplyr::filter(., Real_Date == time_plot_use) %>%
        st_as_sf(., coords = c("Lon", "Lat"), crs = 4326, remove = FALSE) %>%
        st_transform(., crs = 32619) %>%
        dplyr::select(., Lon, Lat, Dens, geometry)

      pred_idw <- idw(Dens ~ 1, data_df, pred_grid2) %>%
        mutate(.,
          "Biomass" = var1.pred * 625,
          "Log_Biomass" = log(Biomass + 1)
        )

      plot_out <- ggplot() +
        geom_sf(data = pred_idw, aes(fill = Log_Biomass, color = Log_Biomass)) +
        scale_fill_viridis_c(name = "Log biomass (kg)", option = "viridis", na.value = "transparent", limits = rast_lims) +
        scale_color_viridis_c(name = "Log biomass (kg)", option = "viridis", na.value = "transparent", limits = rast_lims) +
        annotate("text", x = lab_lon, y = lab_lat, label = plot_label) +
        geom_sf(data = land_sf, fill = land_color, lwd = 0.2, na.rm = TRUE) +
        coord_sf(xlim = xlim, ylim = ylim, expand = FALSE, crs = 4326) +
        theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "pt"))
    }
    if (panel_or_gif == "panel") {
      # Panel plot
      all_plot <- wrap_plots(rasts_out, ncol = panel_cols, nrow = panel_rows, guides = "collect", theme(plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "pt")))
      ggsave(filename = paste0(out_dir, "/", nice_category_names, "_", pred_label, "_", "density.png"), all_plot, width = 11, height = 8, units = "in")
      return(all_plot)
    } else {
      # Make a gif
      plot_loop_func <- function(plot_list) {
        for (i in seq_along(plot_list)) {
          plot_use <- plot_list[[i]]
          print(plot_use)
        }
      }
      invisible(save_gif(plot_loop_func(rasts_out), paste0(out_dir, "/", nice_category_names, "_", pred_label, "_", "density.gif"), delay = 0.75, progress = FALSE))
    }
  }
}


#' @title Get VAST point predictions
#'
#' @description Generates a dataframe with observed and VAST model predictions at sample locations
#'
#' @param vast_fit = A VAST `fit_model` object.
#' @param use_PredTF_only = Logical TRUE/FALSE. If TRUE, then only the locations specified as PredTF == 1 will be extracted. Otherwise, all points will be included.
#' @param nice_category_names = A character string to define species/model run and used in naming the output prediction file.
#' @param out_dir = Output directory to save the dataset
#'
#' @return A dataframe with lat, lon, observations and model predictions
#'
#' @export

vast_get_point_preds<- function(vast_fit, use_PredTF_only, nice_category_names, out_dir) {
  if (FALSE) {
    vast_fit <- mod_comp_res$Fitted_Mod[[1]]
    use_PredTF_only <- TRUE
    nice_category_names <- "Atlantic_cod_test"
    out_dir <- here::here("", "results/tables")
  }

  # Collecting the sample data
  samp_dat <- vast_fit$data_frame %>%
      dplyr::select(., Lat_i, Lon_i, b_i, t_i) %>%
      drop_units()
  names(samp_dat) <- c("Lat", "Lon", "Biomass", "Year")
  samp_dat$Presence <- ifelse(samp_dat$Biomass > 0, 1, 0)

  # Now, getting the model predictions
  pred_dat <- vast_fit$Report

  # Combine em
  samp_pred_out <- data.frame(samp_dat, "Predicted_ProbPresence" = pred_dat$R1_i, "Predicted_Biomass" = pred_dat$D_i)

  # Add PredTF column -- this is 1 if the sample is only going to be used in predicted probability and NOT in estimating the likelihood
  samp_pred_out$PredTF_i <- vast_fit$data_list$PredTF_i

  # Subset if use_PredTF_only is TRUE
  if (use_PredTF_only) {
    samp_pred_out <- samp_pred_out %>%
      dplyr::filter(., PredTF_i == 1)
  }

  # Save and return it
  saveRDS(samp_pred_out, paste0(out_dir, "/", nice_category_names, "_obs_pred.rds"))
  return(samp_pred_out)
}

#' @title Get VAST knot predictions for spatial or spatio-temporal parameters/derived quantities
#'
#' @description Generates a dataframe with VAST model spatial or spatio-temporal parameters/derived quantities at each knot location
#'
#' @param vast_fit = A VAST `fit_model` object.
#' @param spatial_var = An estimated spatial coefficient or predicted value. Currently works for `D_gct`, `R1_gct`, `R2_gct`, `P1_gct`, `P2_gct`, `Omega1_gc`, `Omega2_gc`, `Epsilon1_gct`, `Epsilon2_gct`.
#' @param nice_category_names = A character string to define species/model run and used in naming the output prediction file.
#' @param out_dir = Output directory to save the dataframe
#'
#' @return A dataframe with lat, lon, observations and model predictions
#'
#' @export

vast_get_extrap_spatial <- function(vast_fit, spatial_var, nice_category_names, out_dir) {
  if (FALSE) {
    vast_fit <- vast_fitted
    spatial_var <- "D_gct"
    nice_category_names <- "Atlantic_halibut"
    out_dir <- here::here("", "results/tables")
  }

  # First check the spatial_var, only a certain subset are being used...
  if (!spatial_var %in% c("D_gct", "R1_gct", "R2_gct", "P1_gct", "P2_gct", "Omega1_gc", "Omega2_gc", "Epsilon1_gct", "Epsilon2_gct")) {
    stop(print("Check `spatial_var` input. Currently must be one of `D_gct`, `R1_gct`, `R2_gct`, `P1_gct`, `P2_gct`, `Omega1_gc`, `Omega2_gc`, `Epsilon1_gct`, `Epsilon2_gct`."))
  }

  # Getting prediction array
  pred_array <- vast_fit$Report[[{{ spatial_var }}]]
  if (spatial_var == "D_gct") {
    pred_array <- log(pred_array + 1)
  }

  # Getting time info
  times <- as.character(vast_fit$year_labels)

  # Getting extrapolation grid locations
  spat_data <- vast_fit$extrapolation_list
  loc_g <- spat_data$Data_Extrap[which(spat_data$Data_Extrap[, "Include"] > 0), c("Lon", "Lat")]

  # Creating the dataframe to save...
  df_out_temp <- as.data.frame(pred_array)
  colnames(df_out_temp) <- paste0("Time_", times)
  df_out_temp <- cbind(loc_g, df_out_temp)

  df_out <- df_out_temp %>%
    pivot_longer(., cols = !c("Lon", "Lat"), names_to = "Time", values_to = {{ spatial_var }}) %>%
    arrange(., Time, Lon, Lat)

  # Save and return it
  saveRDS(df_out, paste0(out_dir, "/", nice_category_names, "_", spatial_var, "_df.rds"))
  return(df_out)
}

#' @title Get VAST center of gravity
#'
#' @description Blah
#'
#' @param vast_fit = A VAST `fit_model` object.
#' @param land_sf = Land sf object
#' @param xlim = A two element vector with the min and max longitudes
#' @param ylim = A two element vector with the min and max latitudes
#' @param nice_category_names = Species name
#' @param out_dir = Output directory to save the dataset
#'
#' @return Blah
#'
#' @export

vast_get_cog <- function(vast_fit, all_times, nice_category_names, out_dir) {
  if (FALSE) {
    tar_load(vast_fit)
    tar_load(vast_seasonal_data)
    all_times <- levels(vast_seasonal_data$VAST_YEAR_SEASON)
    summarize <- TRUE
    tar_load(land_sf)
    xlim <- c(-80, -55)
    ylim <- c(35, 50)
    nice_category_names <- nice_category_names
    land_color <- "#d9d9d9"
    out_dir <- paste0(res_root, "plots_maps")
  }

  TmbData <- vast_fit$data_list
  Sdreport <- vast_fit$parameter_estimates$SD

  # Time series steps
  time_ind <- seq(from = 1, to = TmbData$n_t)
  time_labels <- sort(unique(vast_fit$data_frame$t_i)[time_ind])

  # Categories
  categories_ind <- seq(from = 1, to = TmbData$n_c)

  # Get the index information
  SD <- TMB::summary.sdreport(Sdreport)

  # Now, populate array with values
  mean_Z_ctm <- array(NA, dim = c(unlist(TmbData[c("n_c", "n_t")]), 2, 2), dimnames = list(categories_ind, time_labels, c("Lon", "Lat"), c("Estimate", "Std. Error")))
  mean_Z_ctm[] <- SD[which(rownames(SD) == "mean_Z_ctm"), c("Estimate", "Std. Error")]

  index_res_array <- mean_Z_ctm

  # Data manipulation to get out out the array and to something more "plottable"
  for (i in seq_along(categories_ind)) {
    index_array_temp <- index_res_array[i, , , ]

    index_res_temp_est <- data.frame("Time" = as.numeric(rownames(index_array_temp[, , 1])), "Category" = categories_ind[i], index_array_temp[, , 1])
    index_res_temp_sd <- data.frame("Time" = as.numeric(rownames(index_array_temp[, , 1])), "Category" = categories_ind[i], index_array_temp[, , 2])
    names(index_res_temp_sd)[3:4] <- c("Lon_SD", "Lat_SD")
    index_res_temp_out <- index_res_temp_est %>%
      left_join(., index_res_temp_sd) %>%
      mutate(.,
        "Lon_Min" = Lon - Lon_SD,
        "Lon_Max" = Lon + Lon_SD,
        "Lat_Min" = Lat - Lat_SD,
        "Lat_Max" = Lat + Lat_SD
      )

    if (i == 1) {
      index_res_out <- index_res_temp_out
    } else {
      index_res_out <- bind_rows(index_res_out, index_res_temp_out)
    }
  }

  index_res_out$Date <- factor(all_times, levels = all_times)

  # Date info
  index_res_out <- index_res_out %>%
    mutate(., Year = as.numeric(gsub("([0-9]+).*$", "\\1", Date)))

  if (any(str_detect(as.character(index_res_out$Date), LETTERS))) {
    index_res_out$Date <- as.Date(paste(index_res_out$Year, ifelse(grepl("SPRING", index_res_out$Date), "-04-15",
      ifelse(grepl("SUMMER", index_res_out$Date), "-07-15", "-10-15")
    ), sep = ""))
  } else {
    index_res_out$Date <- as.Date(paste(index_res_out$Year, "-06-15", sep = ""))
  }

  # Save and return the full table
  write.csv(index_res_out, paste0(out_dir, "/COG_", nice_category_names, ".csv"))
  return(index_res_out)
}

plot_cog <- function(cog_df, df_crs = "+proj=utm +zone=19 +datum=WGS84 +units=km +no_defs", plot_crs = 4326, summarize = TRUE, land_sf, xlim, ylim, nice_category_names, land_color = "#d9d9d9", color_pal = NULL, out_dir) {
  if (FALSE) {
    tar_load(vast_cog)
    cog_df <- vast_cog
    df_crs <- "+proj=utm +zone=19 +datum=WGS84 +units=km +no_defs"
    plot_crs <- 4326
    summarize <- TRUE
    tar_load(land_sf)
    xlim <- c(-80, -55)
    ylim <- c(35, 50)
    nice_category_names <- nice_category_names
    land_color <- "#d9d9d9"
    out_dir <- paste0(res_root, "plots_maps")
  }

  # Plot is either going to be annual or each individual season
  if (summarize) {
    cog_df <- cog_df %>%
      group_by(., Year, Category, .drop = FALSE) %>%
      summarize_at(., vars(c("Lon", "Lat", "Lon_Min", "Lon_Max", "Lat_Min", "Lat_Max")), mean, na.rm = TRUE)

    # First, the map.
    cog_sf <- st_as_sf(cog_df, coords = c("Lon", "Lat"), crs = df_crs)

    # Transform to be in WGS84
    cog_sf_wgs84 <- st_transform(cog_sf, crs = plot_crs)

    # Base map
    cog_plot <- ggplot() +
      geom_sf(data = cog_sf_wgs84, aes(fill = Year), size = 2, shape = 21) +
      scale_fill_viridis_c(name = "Year", limits = c(min(cog_sf_wgs84$Year), max(cog_sf_wgs84$Year))) +
      geom_sf(data = land_sf, fill = land_color, lwd = 0.2, na.rm = TRUE) +
      coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
      theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "pt"))

    # Now, the lon/lat time series
    lon_lat_df <- cog_sf_wgs84 %>%
      data.frame(st_coordinates(.))
    lon_lat_min <- st_as_sf(cog_df, coords = c("Lon_Min", "Lat_Min"), crs = df_crs) %>%
      st_transform(., plot_crs) %>%
      data.frame(st_coordinates(.)) %>%
      dplyr::select(c("X", "Y"))

    names(lon_lat_min) <- c("Lon_Min_WGS", "Lat_Min_WGS")
    lon_lat_max <- st_as_sf(cog_df, coords = c("Lon_Max", "Lat_Max"), crs = df_crs) %>%
      st_transform(., plot_crs) %>%
      data.frame(st_coordinates(.)) %>%
      dplyr::select(c("X", "Y"))
    names(lon_lat_max) <- c("Lon_Max_WGS", "Lat_Max_WGS")

    lon_lat_df <- cbind(lon_lat_df, lon_lat_min, lon_lat_max)
    names(lon_lat_df)[8:9] <- c("Lon", "Lat")
    lon_lat_df$Date <- as.Date(paste0(lon_lat_df$Year, "-06-15"))

    if (!is.null(color_pal)) {
      colors_use <- color_pal
    } else {
      color_pal <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854")
      colors_use <- color_pal[1:length(unique(lon_lat_df$Category))]
    }

    lon_ts <- ggplot() +
      geom_ribbon(data = lon_lat_df, aes(x = Date, ymin = Lon_Min_WGS, ymax = Lon_Max_WGS), fill = "#66c2a5", alpha = 0.3) +
      geom_line(data = lon_lat_df, aes(x = Date, y = Lon), color = "#66c2a5", lwd = 2) +
      # scale_fill_manual(name = "Category", values = '#66c2a5') +
      scale_x_date(date_breaks = "5 year", date_labels = "%Y") +
      ylab("Center of longitude") +
      xlab("Date") +
      theme_bw() +
      theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

    lat_ts <- ggplot() +
      geom_ribbon(data = lon_lat_df, aes(x = Date, ymin = Lat_Min_WGS, ymax = Lat_Max_WGS), fill = "#66c2a5", alpha = 0.3) +
      geom_line(data = lon_lat_df, aes(x = Date, y = Lat), color = "#66c2a5", lwd = 2) +
      # scale_fill_manual(name = "Category", values = '#66c2a5') +
      scale_x_date(date_breaks = "5 year", date_labels = "%Y") +
      ylab("Center of latitude") +
      xlab("Date") +
      theme_bw() +
      theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

    plot_out <- (cog_plot) / (lon_ts + lat_ts) + plot_layout(ncol = 1, nrow = 2, widths = c(0.75, 1), heights = c(0.75, 1))
  }

  # Save and return it
  ggsave(plot_out, file = paste(out_dir, "/COG_", "_", nice_category_names, ".jpg", sep = ""))
  return(plot_out)
}

plot_proj_dens_gif <- function(vast_projections, plot_scale = "log", land_sf, mask = region_shapefile, xlim = c(-78.5, -54.5), ylim = c(34.5, 48.25), lab_lat = 36, lab_lon = -67.5, out_dir, land_color = "#d9d9d9", ...) {
  # For debugging
  if (FALSE) {
    tar_load(vast_projection_summ_dens)
    vast_projections <- vast_projection_summ_dens
    plot_scale <- "log"
    tar_load(land_sf)
    tar_load(region_shapefile)
    mask <- region_shapefile
    xlim <- c(-78.5, -54.5)
    ylim <- c(34.5, 48.25)
    lab_lat <- 36
    lab_lon <- -67.5
    land_color <- "#d9d9d9"
    out_dir <- paste0(res_root, "plots_maps")
  }

  if (plot_scale == "log") {
    rasts_range <- log(vast_projections[, 4:6] + 1)
    rast_lims_min <- min(c(0, min(rasts_range)))
    rast_lims_max <- max(c(round(max(rasts_range) + 0.0000001, 2), max(rasts_range)))
    rast_lims <- c(rast_lims_min, rast_lims_max)
  } else {
    rasts_range <- vast_projections[, 4:6]
    rast_lims_min <- min(c(0, min(rasts_range)))
    rast_lims_max <- max(c(round(max(rasts_range) + 0.0000001, 2), max(rasts_range)))
    rast_lims <- c(rast_lims_min, rast_lims_max)
  }

  rasts_out <- vector("list", length(unique(vast_projections$Time)))

  for (tI in 1:length(unique(vast_projections$Time))) {
    data_df <- vast_projections[vast_projections$Time == unique(vast_projections$Time)[tI], ]

    # Interpolation
    pred_df <- na.omit(data.frame("x" = data_df$Lon, "y" = data_df$Lat, "layer" = data_df$Prob_0.5))
    pred_df_interp <- interp(pred_df[, 1], pred_df[, 2], pred_df[, 3],
      duplicate = "mean", extrap = TRUE,
      xo = seq(-87.99457, -57.4307, length = 115),
      yo = seq(22.27352, 48.11657, length = 133)
    )

    pred_df_interp_final <- data.frame(expand.grid(x = pred_df_interp$x, y = pred_df_interp$y), z = c(round(pred_df_interp$z, 2)))
    pred_sp <- st_as_sf(pred_df_interp_final, coords = c("x", "y"), crs = 4326)

    pred_df_temp <- pred_sp[which(st_intersects(pred_sp, mask, sparse = FALSE) == TRUE), ]
    coords_keep <- as.data.frame(st_coordinates(pred_df_temp))
    row.names(coords_keep) <- NULL
    pred_df_use <- data.frame(cbind(coords_keep, "z" = as.numeric(pred_df_temp$z)))
    names(pred_df_use) <- c("x", "y", "z")

    time_plot_use <- unique(vast_projections$Time)[tI]

    if (plot_scale == "log") {
      pred_df_use$z <- log(pred_df_use$z + 1)
    }

    rasts_out[[tI]] <- ggplot() +
      geom_tile(data = pred_df_use, aes(x = x, y = y, fill = z, color = z)) +
      scale_fill_viridis_c(name = "Density kg/km2", option = "viridis", na.value = "transparent", limits = rast_lims) +
      scale_color_viridis_c(name = "Density kg/km2", option = "viridis", na.value = "transparent", limits = rast_lims) +
      annotate("text", x = lab_lon, y = lab_lat, label = time_plot_use) +
      geom_sf(data = land_sf, fill = land_color, lwd = 0.2, na.rm = TRUE) +
      coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
      theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "pt"))
  }

  # Make a gif
  plot_loop_func <- function(plot_list) {
    for (i in seq_along(plot_list)) {
      plot_use <- plot_list[[i]]
      print(plot_use)
    }
  }

  invisible(save_gif(plot_loop_func(rasts_out), paste0(out_dir, "/", nice_category_names, "_", climate_scenario, "_", "Density.gif"), delay = 0.75, progress = FALSE))
}




#' @title VAST manual predictions
#'
#' @description Make predictions from a fitted VAST model given new covariate data, with a few key assumptions. In particular, right now, just the "environment_only" option is implemented. This will assume that the spatio-temporal surfaces are constant in the "future" and all projected occurrence changes arise mainly from changes in environmental conditions.
#'
#' @param vast_fit = A VAST `fit_model` object.
#' @param pred_type = Character string signaling the prediction type to use. Currently only implemented as "environment_only".
#' @param obs_or_gird = Either "obs" or "grid" signaling if the predictions should be done at an observation level OR at a grid level. For climate projections, we are using the grid option as the "new_covariate_data" uses extrapolation grid location and then bilinear interpolation to calculate covariate values at these locations.
#' @param pred_label = Character string to characterize the predictions, included in the output file name. For example, this could be based on the climate scenario used.
#' @param set_re_zero = Logical flag. IF true, this will set all random effects to 0.
#' @param new_covariate_data = A dataframe with new covariate data to use to make the predictions. Has to include all of the necessary columns to match those used during the model fitting process.
#' @param new_catchability_data = A dataframe with new catchability data to use to make the predictions. Has to include all of the necessary columns to match those used during the model fitting process.
#' @param nice_category_names = A character string to define species/model run and used in naming the output prediction file.
#' @param out_dir = Output directory to save the dataset
#'
#' @return A dataframe with columns for time, space, all of the data used to make the predictions, and the final predicted density at each location.
#'
#' @export

vast_predict_manual <- function(vast_fit, pred_type = "environment_only", obs_or_grid = "grid", pred_label, set_re_zero = TRUE, new_covariate_data, new_catchability_data, nice_category_names = nice_category_names, out_dir = paste0(res_root, "prediction_df")) {
  if (FALSE) {
    # For debugging
    vast_fit <- readRDS("/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/targets_output/mod_fits/Atlantic_cod_ST1noRW_fitted_vast.rds")
    nice_category_names <- "Atlantic_cod"
    pred_type <- "environment_only"
    obs_or_grid <- "grid"
    pred_label <- "SSP5_85_mean"
    set_re_zero <- TRUE
    new_covariate_data <- readRDS("~/GitHub/TargetsSDM/data/predict/VAST_post_fit_pred_df_seasonal_mean.rds")
    new_catchability_data <- NULL
    out_dir <- paste0(res_root, "prediction_df")

    # From VAST code...
    # P1_gct(g,c,t) = Omega1_gc(g,c) + beta1_tc(t,c) + Epsilon1_gct(g,c,t) + eta1_gct(g,c,t) + iota_ct(c,t)
    # P2_gct(g,c,t) = Omega2_gc(g,c) + beta2_tc(t,c) + Epsilon2_gct(g,c,t) + eta2_gct(g,c,t)
    # R1_gct(g,c,t) = Type(1.0) - exp( -exp(P1_gct(g,c,t)) )
    # R2_gct(g,c,t) = exp(P1_gct(g,c,t)) / R1_gct(g,c,t) * exp( P2_gct(g,c,t) )
    # D_gct(g,c,t) = exp( P1_gct(g,c,t) + P2_gct(g,c,t) ) # Use this line to prevent numerical over/underflow
  }

  preds_out_full <- data.frame("Model_Date" = factor(new_covariate_data$Year, levels = levels(new_covariate_data$Year)), "Real_Date" = new_covariate_data$Date, "Lat" = new_covariate_data$Lat, "Lon" = new_covariate_data$Lon, "Extrap_ID" = new_covariate_data$Extrap_ID)

  # Something going on here with the predictions...we aren't getting a year covariate estimate for 2015-2019. Need to ask Jim about this...
  new_covariate_data$Year_Cov <- factor(ifelse(as.numeric(new_covariate_data$Year_Cov) > 31, "2015", as.character(new_covariate_data$Year_Cov)), levels = levels(vast_fit$covariate_data$Year_Cov))

  if (pred_type == "environment_only") {
    # Here, we are fixing everything else in the model and only looking at the influence of environmental changes on species occurrence.
    # Extracting what we need from the fitted model object.
    preds_out_sp <- data.frame("Lat" = new_covariate_data$Lat, "Lon" = new_covariate_data$Lon, "Extrap_ID" = new_covariate_data$Extrap_ID) %>%
      distinct(., Extrap_ID, .keep_all = TRUE)

    ## Random effects
    # Omega first, this will be a surface and is constant through time. These values are going to change depending on if the fine_scale option is used. When they are used, values at knot locations (or grid locations) are estimated using a bilinear interpolation. When they are not, just nearest neighbor are used. To account for this, we just need to make sure to update value at locations as we go. If we are interested in observations, we will want to use the A_is and related information. If we are interested in extrapolation grid locations, we will want to use the A_gs matrix. With the A_is matrix, slot i corresponds to the observation ID, slot j to the knot ID, and slot x to the value.
    if (obs_or_grid == "grid") {
      # Keeping track of the grid locations we want to retain
      extrap_ids_check <- readr::parse_number(preds_out_sp$Extrap_ID)

      if (!set_re_zero) {
        preds_out_sp$omega1 <- vast_fit$Report[[{{ "Omega1_gc" }}]][extrap_ids_check, 1]
        preds_out_sp$omega2 <- vast_fit$Report[[{{ "Omega2_gc" }}]][extrap_ids_check, 1]
      } else {
        preds_out_sp$omega1 <- vast_fit$Report[[{{ "Omega1_gc" }}]][extrap_ids_check, 1]
        preds_out_sp$omega2 <- vast_fit$Report[[{{ "Omega2_gc" }}]][extrap_ids_check, 1]
        # preds_out_sp$omega1<- 0
        # preds_out_sp$omega2<- 0
      }

      # Bring that over to preds_out_full -- join based on location...
      preds_out_full <- preds_out_full %>%
        left_join(., preds_out_sp, by = c("Lon" = "Lon", "Lat" = "Lat", "Extrap_ID" = "Extrap_ID"))

      # Average intercepts, these are constant across knots (sometimes, time too)
      beta1 <- data.frame("Model_Date" = factor(unique(new_covariate_data$Year), levels = levels(new_covariate_data$Year)), "beta1" = vast_fit$Report[[{{ "beta1_tc" }}]])
      beta2 <- data.frame("Model_Date" = factor(unique(new_covariate_data$Year), levels = levels(new_covariate_data$Year)), "beta2" = vast_fit$Report[[{{ "beta2_tc" }}]])
      betas <- beta1 %>%
        left_join(., beta2)

      preds_out_full <- preds_out_full %>%
        left_join(., betas, by = c("Model_Date" = "Model_Date"))

      # Spatio_temporal -- similar to omega, but now with a temporal dimension...
      if (!set_re_zero) {
        eps1_proj <- vast_fit$Report[[{{ "Epsilon1_gct" }}]][extrap_ids_check, 1, ]
        eps2_proj <- vast_fit$Report[[{{ "Epsilon2_gct" }}]][extrap_ids_check, 1, ]
      } else {
        eps1_proj <- matrix(0, nrow = nrow(preds_out_sp), ncol = vast_fit$data_list$n_t)
        eps2_proj <- matrix(0, nrow = nrow(preds_out_sp), ncol = vast_fit$data_list$n_t)
      }
    } else {
      # Add knot to observation implementation with A_is here
    }

    # Some work on column names...
    eps1_wide <- data.frame("Extrap_ID" = preds_out_sp$Extrap_ID, eps1_proj)
    colnames(eps1_wide) <- c("Extrap_ID", as.character(unique(new_covariate_data$Year)))

    # Come back to this and make it cleaner!!!
    if (!set_re_zero) {
      for (i in 95:106) {
        col_match <- paste0("2015_", str_extract(names(eps1_wide)[i], "[A-Z]+"))
        eps1_wide[, i] <- eps1_wide[, {{ col_match }}]
      }
    }
    eps1_long <- eps1_wide %>%
      pivot_longer(-Extrap_ID, names_to = "Model_Date", values_to = "epsilon1") %>%
      mutate(., "Model_Date" = factor(Model_Date, levels = levels(new_covariate_data$Year))) %>%
      arrange(Model_Date)

    eps2_wide <- data.frame("Extrap_ID" = preds_out_sp$Extrap_ID, eps2_proj)
    colnames(eps2_wide) <- c("Extrap_ID", as.character(unique(new_covariate_data$Year)))

    if (!set_re_zero) {
      # Come back to this and make it cleaner!!!
      for (i in 95:106) {
        col_match <- paste0("2015_", str_extract(names(eps2_wide)[i], "[A-Z]+"))
        eps2_wide[, i] <- eps2_wide[, {{ col_match }}]
      }
    }

    eps2_long <- eps2_wide %>%
      pivot_longer(-Extrap_ID, names_to = "Model_Date", values_to = "epsilon2") %>%
      mutate(., "Model_Date" = factor(Model_Date, levels = levels(new_covariate_data$Year))) %>%
      arrange(Model_Date)

    # Combine them
    eps_out <- left_join(eps1_long, eps2_long)

    # Add to preds...
    preds_out_full <- preds_out_full %>%
      left_join(., eps_out, by = c("Extrap_ID" = "Extrap_ID", "Model_Date" = "Model_Date"))

    ## Habitat covariates using gridded environmental data
    # For this, I think it is going to be easier to do the predictions and then supply those single eta1 and eta2 values to a simpler function
    x1_covs <- vast_fit$ParHat[[{{ "gamma1_cp" }}]]
    colnames(x1_covs) <- attributes(vast_fit$data_list$X1_gctp)$dimnames[[4]]
    x2_covs <- vast_fit$ParHat[[{{ "gamma2_cp" }}]]
    colnames(x2_covs) <- attributes(vast_fit$data_list$X2_gctp)$dimnames[[4]]

    # Create a model matrix
    mod_mat <- model.matrix(object = vast_fit$X1_formula, data = new_covariate_data, contrasts.arg = list(Season = contrasts(new_covariate_data$Season, contrasts = FALSE), Year_Cov = contrasts(new_covariate_data$Year_Cov, contrasts = FALSE)))

    # Not necessary in this simple example, but matters when we have new years in the new_covariate_data
    cols_keep <- colnames(mod_mat) %in% colnames(x1_covs)
    mod_mat <- mod_mat[, cols_keep]

    # Generate predictions from covariate values and model matrix
    preds_out_full$eta1 <- as.numeric(t(as.vector(x1_covs) %*% t(mod_mat)))
    preds_out_full$eta2 <- as.numeric(t(as.vector(x2_covs) %*% t(mod_mat)))

    ## Make the calculations...
    preds_out_full <- preds_out_full %>%
      mutate(.,
        "P1" = omega1 + beta1 + epsilon1 + eta1,
        "P2" = omega2 + beta2 + epsilon2 + eta2,
        "R1" = 1 - exp(-exp(P1)),
        "R2" = (exp(P1) / R1) * exp(P2),
        "Dens" = R1 * R2
      )

    # Save and return it
    write.csv(preds_out_full, file = paste(out_dir, "/", nice_category_names, "_", pred_label, "_projections.csv", sep = ""))
    return(preds_out_full)
  }
}

vast_predict_idw <- function(vast_fit, manual_pred_df, grid_space_utm, nice_category_names, climate_scenario, all_times = all_times, out_dir, ...) {
  if (FALSE) {
    vast_fit <- readRDS("/Users/aallyn/Library/CloudStorage/Box-Box/Mills Lab/Projects/sdm_workflow/targets_output/mod_fits/American_lobster_STnoRW_fitted_vast.rds")
    manual_pred_df <- read.csv(paste0(res_root, "prediction_df/American_lobster_SSP5_85_mean_projections.csv"))
    nice_category_names <- "American_lobster"
    climate_scenario <- "SSP5_85_mean"
    all_times <- unique(manual_pred_df$Real_Date)
    plot_times <- NULL
    grid_space_utm <- 25000
  }

  # We want a prediction grid...
  mask_utm <- vast_mesh_to_sf(vast_fit)$triangles %>%
    st_union() %>%
    st_transform(., crs = 32619)
  pred_grid <- mask_utm %>%
    st_make_grid(., cellsize = grid_space_utm, what = "polygons", square = TRUE)

  pred_grid2 <- pred_grid %>%
    st_intersection(., mask_utm)

  pred_idw_out_all <- vector("list", length = length(all_times))

  for (tI in seq_along(all_times)) {
    time_use <- all_times[tI]

    data_df <- manual_pred_df %>%
      dplyr::filter(., Real_Date == time_use) %>%
      st_as_sf(., coords = c("Lon", "Lat"), crs = 4326, remove = FALSE) %>%
      st_transform(., crs = 32619) %>%
      dplyr::select(., Lon, Lat, Dens, geometry)

    pred_idw_sf <- idw(Dens ~ 1, data_df, pred_grid) %>%
      mutate(.,
        "Biomass" = var1.pred * 625,
        "Log_Biomass" = log(Biomass + 1),
        "Date" = time_use,
        "Species" = nice_category_names,
        "Climate_Scenario" = pred_label
      )

    pred_idw_temp <- pred_idw_sf %>%
      dplyr::select(., x, y, Date, Species, Climate_Scenario, Biomass, Log_Biomass)

    na_ids <- seq(from = 1, to = nrow(pred_idw_sf))
    intersect_list <- st_intersects(pred_idw_sf, mask_utm)
    na_ids <- na_ids[-which(na_ids %in% unlist(t(intersect_list)))]

    pred_idw_temp$Biomass[na_ids] <- NA
    pred_idw_temp$Log_Biomass[na_ids] <- NA

    pred_idw_temp <- pred_idw_temp %>%
      st_drop_geometry(.)

    pred_idw_out_all[[tI]] <- pred_idw_temp
  }

  # Save and return that
  pred_idw_out <- do.call("rbind", pred_idw_out_all)
  write.csv(pred_idw_out, file = paste0(out_dir, nice_category_names, "_", climate_scenario, "_idw_predictions.csv"))
  return(pred_idw_out)
}

pred_idw_to_ncdf <- function(pred_idw_list, nice_category_names, climate_scenario, out_dir) {
  # Into a dataframe...
  pred_idw_wide <- do.call("rbind", pred_idw_out_all) %>%
    dplyr::select(., -Log_Biomass) %>%
    pivot_wider(., names_from = Date, values_from = Biomass)

  # NetCDF bits
  lon <- pred_idw_wide$x
  lat <- pred_idw_wide$y
  time <- seq(from = 1, to = length(all_times)) - 1
  tunits <- "Year-season since spring 1985"

  nlon <- length(lon)
  nlat <- length(lat)
  ntime <- length(time)

  bio_mat <- as.matrix(pred_idw_wide[5:(5 + ntime - 1)])
  bio_array <- array(bio_mat, dim = c(nlon, nlat, ntime))

  ncpath <- out_dir
  ncname <- paste(nice_category_names, "_", climate_scenario, "_ncdf4", sep = "")
  ncfname <- paste(ncpath, ncname, ".nc", sep = "")
  dname <- "biomass"

  # define dimensions
  londim <- ncdim_def("lon", "degrees_east_crs_32619", as.double(lon))
  latdim <- ncdim_def("lat", "degrees_north_crs_32619", as.double(lat))
  timedim <- ncdim_def("time", tunits, as.double(time))

  # define variables
  fillvalue <- 1e32
  dlname <- "Predicted biomass"
  bio_def <- ncvar_def("biomass", "kg", list(londim, latdim, timedim), fillvalue, dlname, prec = "single")

  # Create netCDF and put in arrays
  ncout <- nc_create(ncfname, list(bio_def), force_v4 = TRUE)

  # put variables
  ncvar_put(ncout, bio_def, bio_array)

  # put additional attributes into dimension and data variables
  ncatt_put(ncout, "lon", "axis", "X") # ,verbose=FALSE) #,definemode=FALSE)
  ncatt_put(ncout, "lat", "axis", "Y")
  ncatt_put(ncout, "time", "axis", "T")

  # add global attributes
  ncatt_put(ncout, 0, "title", "SDM Workflow VAST model predicted biomass")
  ncatt_put(ncout, 0, "institution", "Gulf of Maine Research Institute")
  # ncatt_put(ncout,0,"source",datasource$value)
  # ncatt_put(ncout,0,"references",references$value)
  history <- paste("Andrew Allyn", date(), sep = ", ")
  ncatt_put(ncout, 0, "history", history)
  # ncatt_put(ncout,0, "Conventions", Conventions$value)

  # Get a summary of the created file:
  nc_close(ncout)
}

get_range_edges<- function(vast_fit, category_names = nice_category_names, strata_names = NULL, all_times, n_samples = 100, quantiles = c(0.05, 0.5, 0.95), calculate_relative_to_average = FALSE){

if(FALSE){
   vast_fit = fit
   all_times = all_times
   category_names = nice_category_names
   strata_names = NULL
   n_samples = 100
   quantiles = c(0.05, 0.5, 0.95)
   calculate_relative_to_average = FALSE
}
  
  # Unpack
  Sdreport = vast_fit$parameter_estimates$SD 
  Report = vast_fit$Report 
  TmbData = vast_fit$data_list 
  Obj = vast_fit$tmb_list$Obj
  year_labels = all_times

  # Default inputs
  if( is.null(strata_names) ) strata_names = 1:TmbData$n_l
  if( is.null(category_names) ) category_names = 1:TmbData$n_c
  if( is.null(colnames(TmbData$Z_gm)) ){
    m_labels = paste0("axis", 1:ncol(TmbData$Z_gm))
  }else{
    m_labels = colnames(TmbData$Z_gm)
  }
  
  
  ##### Local function
  D_gcyr = sample_variable( Sdreport=Sdreport, Obj=Obj, variable_name="D_gct", n_samples=n_samples )
  
  # Calculate quantiles from observed and sampled densities D_gcy
  E_zctm = array(NA, dim=c(length(quantiles),dim(Report$D_gct)[2:3],ncol(TmbData$Z_gm)) )
  E_zctmr = array(NA, dim=c(length(quantiles),dim(Report$D_gct)[2:3],ncol(TmbData$Z_gm),n_samples) )
  Mean_cmr = array(NA, dim=c(dim(Report$D_gct)[2],ncol(TmbData$Z_gm),n_samples) )
  prop_zctm = array(NA, dim=c(dim(Report$D_gct)[1:3],ncol(TmbData$Z_gm)) )
  prop_zctmr = array(NA, dim=c(dim(Report$D_gct)[1:3],ncol(TmbData$Z_gm),n_samples) )
  for( rI in 0:n_samples ){
    for( mI in 1:ncol(TmbData$Z_gm) ){
      order_g = order(TmbData$Z_gm[,mI], decreasing=FALSE)
      if(rI==0) prop_zctm[,,,mI] = apply( Report$D_gct, MARGIN=2:3, FUN=function(vec){cumsum(vec[order_g])/sum(vec)} )
      if(rI>=0) prop_zctmr[,,,mI,rI] = apply( D_gcyr[,,,rI,drop=FALSE], MARGIN=2:3, FUN=function(vec){cumsum(vec[order_g])/sum(vec)} )
      
      # Calculate edge
      for( cI in 1:dim(E_zctm)[2] ){
        if(rI>=1){
          if( calculate_relative_to_average==TRUE ){
            Mean_cmr[cI,mI,rI] = weighted.mean( as.vector(TmbData$Z_gm[,mI]%o%rep(1,dim(Report$D_gct)[3])), w=as.vector(D_gcyr[,cI,,rI]) )
          }else{
            Mean_cmr[cI,mI,rI] = 0
          }
        }
        for( zI in 1:dim(E_zctm)[1] ){
          for( tI in 1:dim(E_zctm)[3] ){
            if(rI==0){
              index_tmp = which.min( (prop_zctm[,cI,tI,mI]-quantiles[zI])^2 )
              E_zctm[zI,cI,tI,mI] = TmbData$Z_gm[order_g[index_tmp],mI]
            }
            if(rI>=1){
              index_tmp = which.min( (prop_zctmr[,cI,tI,mI,rI]-quantiles[zI])^2 )
              E_zctmr[zI,cI,tI,mI,rI] = TmbData$Z_gm[order_g[index_tmp],mI] - Mean_cmr[cI,mI,rI]
            }
          }}
      }
    }}
  SE_zctm = apply( E_zctmr, MARGIN=1:4, FUN=sd )
  Edge_zctm = abind::abind( "Estimate"=E_zctm, "Std_Error"=SE_zctm, along=5 )
  dimnames(Edge_zctm)[[1]] = paste0("quantile_",quantiles)
  
  # transform matrix into a dataframe
  Edge_df <- Edge_zctm %>%
      as_tibble(., rownames = "Quantile") %>%
      pivot_longer(., cols = -Quantile) %>%
      separate(., name, into = c("Category", "Year", "Axis", "Quantity"), sep = "[.]")
 
  Edge_df$Year <- year_labels[as.numeric(Edge_df$Year)] # DANGER--would prefer to carry through real year values
  
  # Rename axes
  Edge_df$Axis <- ifelse(Edge_df$Axis == "1", "Longitude", "Latitude")

  # Slight formatting to get things to be a bit easier to plot
  Edge_df_mu <- Edge_df %>%
      filter(., Quantity == "Estimate") %>%
      pivot_wider(., names_from = Axis, values_from = value) %>%
      dplyr::select(., -Quantity)
  colnames(Edge_df_mu)[4:5]<- c("Longitude_Mean", "Latitude_Mean")

  Edge_df_sd <- Edge_df %>%
      filter(., Quantity == "Std_Error") %>%
      pivot_wider(., names_from = Axis, values_from = value) %>%
      dplyr::select(., -Quantity)
  colnames(Edge_df_sd)[4:5]<- c("Longitude_SD", "Latitude_SD")

  Edge_df_out <- Edge_df_mu %>%
      left_join(., Edge_df_sd)
  
  # Return
  return(Edge_df_out)
}


plot_cog <- function(cog_df, df_crs = "+proj=utm +zone=19 +datum=WGS84 +units=km +no_defs", plot_crs = 4326, summarize = TRUE, land_sf, xlim, ylim, nice_category_names, land_color = "#d9d9d9", color_pal = NULL, out_dir) {
  if (FALSE) {
    tar_load(vast_cog)
    cog_df <- vast_cog
    df_crs <- "+proj=utm +zone=19 +datum=WGS84 +units=km +no_defs"
    plot_crs <- 4326
    summarize <- TRUE
    tar_load(land_sf)
    xlim <- c(-80, -55)
    ylim <- c(35, 50)
    nice_category_names <- nice_category_names
    land_color <- "#d9d9d9"
    out_dir <- paste0(res_root, "plots_maps")

    cog_df = cog_res
    df_crs = "+proj=utm +zone=19 +datum=WGS84 +units=km +no_defs"
    plot_crs = 4326
    summarize = TRUE
    land_sf = land_use
    xlim = xlim_use
    ylim = ylim_use
    nice_category_names = nice_category_names
    land_color = "#d9d9d9"
    color_pal = NULL
    out_dir = here("results/vast")
  }

  # Plot is either going to be annual or each individual season
  if (summarize) {
    cog_df <- cog_df %>%
      group_by(., Year, Category, .drop = FALSE) %>%
      summarize_at(., vars(c("Lon", "Lat", "Lon_Min", "Lon_Max", "Lat_Min", "Lat_Max")), mean, na.rm = TRUE)

    # First, the map.
    cog_sf <- st_as_sf(cog_df, coords = c("Lon", "Lat"), crs = df_crs)

    # Transform to be in WGS84
    cog_sf_wgs84 <- st_transform(cog_sf, crs = plot_crs)

    # Base map
    cog_plot <- ggplot() +
      geom_sf(data = cog_sf_wgs84, aes(fill = Year), size = 2, shape = 21) +
      scale_fill_viridis_c(name = "Year", limits = c(min(cog_sf_wgs84$Year), max(cog_sf_wgs84$Year))) +
      geom_sf(data = land_sf, fill = land_color, lwd = 0.2, na.rm = TRUE) +
      coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
      theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "pt"))

    # Now, the lon/lat time series
    lon_lat_df <- cog_sf_wgs84 %>%
      data.frame(st_coordinates(.))
    lon_lat_min <- st_as_sf(cog_df, coords = c("Lon_Min", "Lat_Min"), crs = df_crs) %>%
      st_transform(., plot_crs) %>%
      data.frame(st_coordinates(.)) %>%
      dplyr::select(c("X", "Y"))

    names(lon_lat_min) <- c("Lon_Min_WGS", "Lat_Min_WGS")
    lon_lat_max <- st_as_sf(cog_df, coords = c("Lon_Max", "Lat_Max"), crs = df_crs) %>%
      st_transform(., plot_crs) %>%
      data.frame(st_coordinates(.)) %>%
      dplyr::select(c("X", "Y"))
    names(lon_lat_max) <- c("Lon_Max_WGS", "Lat_Max_WGS")

    lon_lat_df <- cbind(lon_lat_df, lon_lat_min, lon_lat_max)
    names(lon_lat_df)[8:9] <- c("Lon", "Lat")
    lon_lat_df$Date <- as.Date(paste0(lon_lat_df$Year, "-06-15"))

    if (!is.null(color_pal)) {
      colors_use <- color_pal
    } else {
      color_pal <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854")
      colors_use <- color_pal[1:length(unique(lon_lat_df$Category))]
    }

    lon_ts <- ggplot() +
      geom_ribbon(data = lon_lat_df, aes(x = Date, ymin = Lon_Min_WGS, ymax = Lon_Max_WGS), fill = "#66c2a5", alpha = 0.3) +
      geom_line(data = lon_lat_df, aes(x = Date, y = Lon), color = "#66c2a5", lwd = 2) +
      # scale_fill_manual(name = "Category", values = '#66c2a5') +
      scale_x_date(date_breaks = "5 year", date_labels = "%Y") +
      ylab("Center of longitude") +
      xlab("Date") +
      theme_bw() +
      theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

    lat_ts <- ggplot() +
      geom_ribbon(data = lon_lat_df, aes(x = Date, ymin = Lat_Min_WGS, ymax = Lat_Max_WGS), fill = "#66c2a5", alpha = 0.3) +
      geom_line(data = lon_lat_df, aes(x = Date, y = Lat), color = "#66c2a5", lwd = 2) +
      # scale_fill_manual(name = "Category", values = '#66c2a5') +
      scale_x_date(date_breaks = "5 year", date_labels = "%Y") +
      ylab("Center of latitude") +
      xlab("Date") +
      theme_bw() +
      theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

    plot_out <- (cog_plot) / (lon_ts + lat_ts) + plot_layout(ncol = 1, nrow = 2, widths = c(0.75, 1), heights = c(0.75, 1))
  }

  # Save and return it
  ggsave(plot_out, file = paste(out_dir, "/COG_", "_", nice_category_names, ".jpg", sep = ""))
  return(plot_out)
}


check_estimability_aja<- function (obj, h, version) {
  if(FALSE){
    obj = vast_fit$tmbl_list$obj
    h =  substitute()
    version = vast_fit$settings$version
  }
    ParHat = TMBhelper:::extract_fixed(obj)
    Gr = obj$gr(ParHat)
    if (any(Gr > 0.01)) 
        stop("Some gradients are high, please improve optimization and only then use `Check_Identifiable`")
    List = NULL
    if (missing(h)) {
        List[["Hess"]] = optimHess(par = ParHat, fn = obj$fn, 
            gr = obj$gr)
    }
    else {
        List[["Hess"]] = h
    }
    List[["Eigen"]] = eigen(List[["Hess"]])
    List[["WhichBad"]] = which(List[["Eigen"]]$values < sqrt(.Machine$double.eps))
    if (length(List[["WhichBad"]]) == 0) {
        message("All parameters are estimable")
    }
    else {
        RowMax = apply(List[["Eigen"]]$vectors[, List[["WhichBad"]], 
            drop = FALSE], MARGIN = 1, FUN = function(vec) {
            max(abs(vec))
        })
        List[["BadParams"]] = data.frame(Param = names(obj$par), 
            MLE = ParHat, Param_check = ifelse(RowMax > 0.1, 
                "Bad", "OK"))
        print(List[["BadParams"]])
    }
    return(invisible(List))
}

###############################################################################################################################
##
##  DHARMa utility functions -- from A Gruss
##
###############################################################################################################################

######## Define the "getQuantile" function 
getQuantile <- function ( simulations, observed, integerResponse, method = c( "PIT", "traditional" ) ) {

	method = match.arg( method )
	n = length( observed )
  	if ( nrow( simulations ) != n ) stop( "Wrong dimension of simulations" )
  	nSim = ncol( simulations )
	if ( method == "traditional" ) {
		if ( integerResponse == F ) {
			if( any( duplicated( observed ) ) ) 
				message( paste0( "Model family was recognized or set as continuous, ", 
					"but duplicate values were detected in the response. ", 
					"Consider whether you are fitting an appropriate model." ) )
			values = as.vector( simulations )[duplicated( as.vector( simulations ) )]
      		if ( length( values ) > 0 ) {
        			if ( all( values %% 1 == 0 ) ) {
          				integerResponse = T
          				message( paste0( "Model family was recognized or set as continuous, ", 
						"but duplicate values were detected in the simulation - ", 
						"changing to integer residuals (see ?simulateResiduals for details)" ) )
        			} else {
          				message( paste0( "Duplicate non-integer values found in the simulation. ", 
						"If this is because you are fitting a non-inter valued discrete response model, ", 
						"note that DHARMa does not perform appropriate randomization for such cases." ) )
        			}
			}
		}
		scaledResiduals = rep( NA, n )
    		for ( i in 1 : n ) {
 			if ( integerResponse == T ) {
				scaledResiduals[i] <- DHARMa.ecdf( simulations[i,] + 
					runif( nSim, -0.5, 0.5 ) ) ( observed[i] + runif( 1, -0.5, 0.5 ) )
			} else {
				scaledResiduals[i] <- DHARMa.ecdf( simulations[i,] )( observed[i] )
			}
    		}
	} else {
		scaledResiduals = rep( NA, n ) 
    		for ( i in 1 : n ) {
      		minSim <- mean( simulations[i,] < observed[i] ) 
      		maxSim <- mean( simulations[i,] <= observed[i] ) 
      		if ( minSim == maxSim ) scaledResiduals[i] = minSim
      		else scaledResiduals[i] = runif( 1, minSim, maxSim )
    		}
  	}
  	return( scaledResiduals )

}

######## Define the "getRandomState" function 
getRandomState <- function ( seed = NULL ) {
  
	current = mget( ".Random.seed", envir = .GlobalEnv, ifnotfound = list( NULL ) )[[1]]
	if ( !is.null( seed ) && is.logical( seed ) && seed == F ) {
    			restoreCurrent <- function() { }    
  	} else {
    		restoreCurrent <- function() {
      		if ( is.null( current ) ) 
				rm( ".Random.seed", envir = .GlobalEnv ) 
			else 
				assign( ".Random.seed", current , envir = .GlobalEnv )
    		}    
  	}

  	#### Set seed
  	if ( is.numeric( seed ) ) set.seed( seed )

  	#### Ensure that RNG has been initialized
  	if ( is.null( current ) ) runif( 1 ) 
  
  	randomState = list( seed, state = get( ".Random.seed", globalenv() ), 
		kind = RNGkind(), restoreCurrent = restoreCurrent )  
  	return( randomState )

}

######## Define the "create_DHARMa" function 
create_DHARMa <- function ( simulatedResponse, observedResponse, fittedPredictedResponse = NULL, 
	integerResponse = F, seed = 123, method = c( "PIT", "traditional" ) ) {

		randomState <- getRandomState( seed )
  		on.exit( { randomState$restoreCurrent() } )
  		match.arg( method )
		out = list()
  		out$simulatedResponse = simulatedResponse
  		out$refit = F
  		out$integerResponse = integerResponse
  		out$observedResponse = observedResponse
		if ( !is.matrix( simulatedResponse ) & !is.null( observedResponse ) ) 
			stop( "Either scaled residuals or simulations and observations have to be provided" )
  		if ( ncol( simulatedResponse ) < 2 ) 
			stop( "simulatedResponse with less than 2 simulations provided - cannot calculate residuals on that." )
		if ( ncol( simulatedResponse ) < 10 ) 
			warning( "simulatedResponse with less than 10 simulations provided. This rarely makes sense" )
		out$nObs = length( observedResponse )
		if ( out$nObs < 3 ) 
			stop( "Warning - number of observations < 3 ... this rarely makes sense" )
		if ( ! ( out$nObs == nrow( simulatedResponse ) ) ) 
			stop( "Dimensions of observedResponse and simulatedResponse do not match" )
		out$nSim = ncol( simulatedResponse )
		out$scaledResiduals = getQuantile( simulations = simulatedResponse, observed = observedResponse, 
			integerResponse = integerResponse, method = method )
		if ( is.null( fittedPredictedResponse ) ) {
    			message( "No fitted predicted response provided, using the mean of the simulations" )
    			fittedPredictedResponse = apply( simulatedResponse, 1, mean )
  		}
  		out$fittedPredictedResponse = fittedPredictedResponse
  		out$randomState = randomState
  		class( out ) = "DHARMa"
  		return( out )

}

plot_DHARMa_res<- function(n_samples, fit, response_units =  NULL, out_dir){

    # Extracting objects from fitted model
    reload_model(fit)
    Obj = fit$tmb_list$Obj
    n_g_orig = Obj$env$data$n_g
    Obj$env$data$n_g = 0
    
    # Simulating response
    b_iz = matrix(NA, nrow = length(fit$data_frame$b_i), ncol = n_samples)
    for (zI in 1:n_samples) {
        if (zI %% max(1, floor(n_samples / 10)) == 0) {
            message("  Finished sample ", zI, " of ", n_samples)
        }
        b_iz[, zI] = simulate_data(fit = list(tmb_list = list(Obj = Obj)), type = 1)$b_i
    }
    
    if (any(is.na(b_iz))) {
        stop("Check simulated residuals for NA values")
    }
    
    if (is.null(response_units)) {
        b_iz <- as_units(b_iz, units(fit$data_frame$b_i))
        dharmaRes = create_DHARMa(simulatedResponse = b_iz, observedResponse = fit$data_frame$b_i, fittedPredictedResponse = fit$Report, integer = FALSE)
    } else {
        b_iz <- as_units(b_iz, response_units)
        dharmaRes = create_DHARMa(simulatedResponse = b_iz, observedResponse = as_units(fit$data_frame$b_i, response_units), fittedPredictedResponse = fit$Report, integer = FALSE)
    }
    
    prop_lessthan_i = apply( as.numeric( b_iz ) < outer( as.numeric( fit$data_frame$b_i ), rep( 1, n_samples ) ), MARGIN = 1, FUN = mean )
    prop_lessthanorequalto_i = apply(as.numeric(b_iz) <= outer(
        as.numeric(fit$data_frame$b_i),
        rep(1, n_samples)
    ), MARGIN = 1, FUN = mean)
    PIT_i = runif(min = prop_lessthan_i, max = prop_lessthanorequalto_i, n = length(prop_lessthan_i))
    dharmaRes$scaledResiduals = PIT_i
    
    #### Save the DHARMa residuals
    save(dharmaRes, file = paste0(out_dir, "dharmaRes.RData"))
    
    #### Produce and save an histogram of DHARMa residuals
    val = dharmaRes$scaledResiduals
    val[val == 0] = -0.01
    val[val == 1] = 1.01
    # val = data.frame(dharmaRes$scaledResiduals)
    # val[val == 0] = -0.01
    # val[val == 1] = 1.01

    # out_hist <- ggplot() +
    #     geom_histogram(data = val, aes(x = dharmaRes.scaledResiduals), breaks = seq(-0.02, 1.02, len = 53))
    jpeg(paste0(out_dir, "Histogram_of_DHARMa_residuals.png"), width = 6, height = 7, units = "in", res = 600)
    hist(val,
        breaks = seq(-0.02, 1.02, len = 53), col = c("red", rep("lightgrey", 50), "red"),
        main = "", xlab = "Residuals", cex.main = 2.5
    )
    dev.off()
    
    #### Produce and save a QQ-plot of DHARMa residuals
    jpeg(paste0(out_dir, "QQplot_of_DHARMa_residuals.png"), width = 6, height = 7, units = "in", res = 600)
    gap::qqunif(dharmaRes$scaledResiduals, pch = 2, bty = "n", logscale = F, col = "black", cex = 0.6, main = "", cex.main = 2.5)
    dev.off()
}

Prepare_User_Extrapolation_Data_Fn<- function (input_grid, strata.limits = NULL, projargs = NA, zone = NA,  flip_around_dateline = TRUE, ...) 
{

    if (FALSE) {
        input_grid = extrap_df
        strata.limits = data.frame("STRATA" = c("All", "NMFS", "DFO"))
        projargs = NA
        zone = NA
        flip_around_dateline = TRUE
    }
    
    if (is.null(strata.limits)) {
        strata.limits = data.frame(STRATA = "All_areas")
    }
    message("Using strata ", strata.limits)
    Data_Extrap <- input_grid
    Area_km2_x = Data_Extrap[, "Area_km2"]
    Tmp = cbind(BEST_LAT_DD = Data_Extrap[, "Lat"], BEST_LON_DD = Data_Extrap[,  "Lon"])
    if ("Depth" %in% colnames(Data_Extrap)) {
        Tmp = cbind(Tmp, BEST_DEPTH_M = Data_Extrap[, "Depth"])
    }
    if("STRATA" %in% colnames(Data_Extrap)){
        Tmp = cbind(Tmp, BEST_STRATA = as.character(Data_Extrap[, "STRATA"]))
    }
    a_el = as.data.frame(matrix(NA, nrow = nrow(Data_Extrap), ncol = nrow(strata.limits), dimnames = list(NULL, strata.limits[, "STRATA"])))
    for (l in 1:ncol(a_el)) {
        a_el[, l] = apply(Tmp, MARGIN = 1, FUN = match_strata_fn, 
            strata_dataframe = strata.limits[l, , drop = FALSE])
        a_el[, l] = ifelse(is.na(a_el[, l]), 0, Area_km2_x)
    }
    tmpUTM = project_coordinates(X = Data_Extrap[, "Lon"], Y = Data_Extrap[, 
        "Lat"], projargs = projargs, zone = zone, flip_around_dateline = flip_around_dateline)
    Data_Extrap = cbind(Data_Extrap, Include = 1)
    if (all(c("E_km", "N_km") %in% colnames(Data_Extrap))) {
        Data_Extrap[, c("E_km", "N_km")] = tmpUTM[, c("X", "Y")]
    } else {
        Data_Extrap = cbind(Data_Extrap, E_km = tmpUTM[, "X"], N_km = tmpUTM[, "Y"])
    }
    Return = list(a_el = a_el, Data_Extrap = Data_Extrap, zone = attr(tmpUTM, 
        "zone"), projargs = attr(tmpUTM, "projargs"), flip_around_dateline = flip_around_dateline, 
        Area_km2_x = Area_km2_x)
    return(Return)
}

match_strata_fn<-  function (x, strata_dataframe)  {
    match_latitude_TF = match_longitude_TF = match_depth_TF = match_strata_TF = rep(TRUE, nrow(strata_dataframe))
    if (all(c("south_border", "north_border") %in% names(strata_dataframe))) {
        match_latitude_TF = as.numeric(x["BEST_LAT_DD"]) > strata_dataframe[, "south_border"] & as.numeric(x["BEST_LAT_DD"]) <=  strata_dataframe[, "north_border"]
    }
    if (all(c("west_border", "east_border") %in% names(strata_dataframe))) {
        match_longitude_TF = as.numeric(x["BEST_LON_DD"]) > strata_dataframe[, "west_border"] & as.numeric(x["BEST_LON_DD"]) <=  strata_dataframe[, "east_border"]
    }
    if (all(c("shallow_border", "deep_border") %in% names(strata_dataframe))) {
        match_depth_TF = as.numeric(x["BEST_DEPTH_M"]) > strata_dataframe[, "shallow_border"] & as.numeric(x["BEST_DEPTH_M"]) <= strata_dataframe[, "deep_border"]
    }
    if(names(strata_dataframe) == "STRATA"){
        match_strata_TF = as.character(x["BEST_STRATA"]) == strata_dataframe[, "STRATA"]
    }
    Char = as.character(strata_dataframe[match_latitude_TF & match_longitude_TF &  match_depth_TF & match_strata_TF, "STRATA"])
    return(ifelse(length(Char) == 0, NA, Char))
}

cog_from_dens<- function(vast_fit, proj_dens, proj_ind){

    if (FALSE) {
        # vast_fit = vast_fit
        # proj_dens = dens_df
        # proj_ind = ind_df
    }
    
    # Check units
    units(proj_dens$D_gct) <- units(vast_fit$Report$D_gct)
    units(proj_ind$Index)<- units(vast_fit$Report$Index_ctl)
    
    ## Calculate Index_gctl (D_gct * a_gl)
    # Get extrapolation info
    if (vast_fit$data_list$n_l > 1) {
        extrap_dat <- data.frame(vast_fit$extrapolation_list$Data_Extrap, vast_fit$extrapolation_list$a_el) %>%
            filter(., Region == "All") %>%
            rename(., "a_gl" = "All")
        proj_ind <- proj_ind %>%
            filter(., Region == "All")
    } else {
        extrap_dat <- data.frame(vast_fit$extrapolation_list$Data_Extrap, "a_gl" = vast_fit$extrapolation_list$a_el)
    }
    
    # Join to area and multiply --- this is different!
    if(FALSE){
        check <- data.frame(vast_fit$Report$Index_gctl[, , 1, ]) 
        index_gctl <- proj_dens %>%
            left_join(., extrap_dat) %>%
            drop_na() %>%
            mutate(.,
                "Index_gctl" = D_gct * a_gl
            )
            mutate("RowID" = seq(from = 1, to = nrow(.))) %>%
            pivot_longer(., !RowID, names_to = "Strata", values_to = "Index_gctl") %>%
            filter(., Strata == "Stratum_1")
        str(check)
    }

    index_gctl <- proj_dens %>%
        left_join(., extrap_dat) %>%
        drop_na() %>%
        mutate(.,
            "Index_gctl" = D_gct * a_gl
        )
    
    # Join to overall index ctl
    index_gctl <- index_gctl %>%
        left_join(., proj_ind) %>%
        mutate(., "Index_Prop" = Index_gctl / Index) %>%
        distinct(., Lat, Lon, Time, Sim_Scenario, D_gct, E_km, N_km, Index_gctl, Index, Index_Prop) %>%
        mutate(.,
            "Lon_wt" = E_km * Index_Prop,
            "Lat_wt" = N_km * Index_Prop
        )
    
    # Now get COG measures by multiplying Index_Prob by lon/lat
    cog_out <- index_gctl %>%
        group_by(Time, Sim_Scenario) %>%
        summarize(.,
            "Mean_Lon" = sum(Lon_wt),
            "Mean_Lat" = sum(Lat_wt)
        )
    return(cog_out)
}


eff_area_from_dens<- function(vast_fit, proj_dens, proj_ind){

    if (FALSE) {
        # vast_fit = vast_fit
        # proj_dens = dens_df
        # proj_ind = ind_df
    }
    
    # Check units
    proj_dens <- proj_dens %>%
        drop_units()
    proj_ind <- proj_ind %>%
        drop_units
    
    # One region?
    if(vast_fit$data_list$n_l == 1){
        # Getting mean density
        # Extrapolation info
        extrap_dat <- data.frame(vast_fit$extrapolation_list$Data_Extrap, "a_gl" = vast_fit$extrapolation_list$a_el)
        
        # First, need index_gctl - density multiplied by area of knot
        index_gctl <- proj_dens %>%
            left_join(., extrap_dat, by = c("Lon" = "Lon", "Lat" = "Lat")) %>%
            drop_na() %>%
            mutate(.,
                "Index_gctl" = D_gct * a_gl
            )
        
        # Next, need information on the total index within the area
        index_gctl <- index_gctl %>%
            left_join(., proj_ind)
        
        # Finally, get the mean density and then use that to get the effective area occupied Summarize across knots and get effective area
        eff_area <- index_gctl %>%
            mutate(., "mean_d_ctl_temp" = D_gct * Index_gctl / Index) %>%
            group_by(., Time, Sim_Scenario, Region) %>%
            summarize(.,
                "mean_d_ctl" = sum(mean_d_ctl_temp, na.rm = TRUE)
            ) %>%
            left_join(., proj_ind) %>%
                mutate(., "Eff_Area" = Index / mean_d_ctl) %>%
                dplyr::select(., Time, Sim_Scenario, Region, Eff_Area)
    } else {
        # Multiple regions...
        regs_vec <- as.vector(unlist(vast_fit$settings$strata.limits))
        
        for(i in seq_along(regs_vec)){
            # Extrapolation info
            extrap_dat <- data.frame(vast_fit$extrapolation_list$Data_Extrap, vast_fit$extrapolation_list$a_el) %>%
                filter(., Region == regs_vec[i]) %>%
                rename(., "a_gl" = regs_vec[i])
            proj_ind_use2 <- proj_ind %>%
                filter(., Region == regs_vec[i])
            
            # First, need index_gctl - density multiplied by area of knot
            index_gctl <- proj_dens %>%
                left_join(., extrap_dat, by = c("Lon" = "Lon", "Lat" = "Lat")) %>%
                drop_na() %>%
                mutate(.,
                    "Index_gctl" = D_gct * a_gl
                )
            
            # Next, need information on the total index within the area
            index_gctl <- index_gctl %>%
                left_join(., proj_ind_use2)
                   
            # Finally, get the mean density and then use that to get the effective area occupied Summarize across knots and get effective area
            eff_area_temp <- index_gctl %>%
                mutate(., "mean_d_ctl_temp" = D_gct * Index_gctl / Index) %>%
                group_by(., Time, Sim_Scenario, Region) %>%
                summarize(.,
                    "mean_d_ctl" = sum(mean_d_ctl_temp, na.rm = TRUE)
                ) %>%
                left_join(., proj_ind_use2) %>%
                    mutate(., "Eff_Area" = Index / mean_d_ctl) %>%
                    dplyr::select(., Time, Sim_Scenario, Region, Eff_Area)
                
            if(i == 1){
                eff_area<- eff_area_temp
            } else {
                eff_area<- bind_rows(eff_area, eff_area_temp)
            }
        }
    }
    return(eff_area)
}

# summary.sim_results<- function (vast_fit, sim_obj, resp_scale, nice_times = NULL, out_t_scale = NULL, nice_category_names = nice_category_names, climate_scenario = climate_scenario, 
#     out_dir) {
#     if (FALSE) {
#         # tar_load(vast_fit)
#         # tar_load(vast_projections)
#         # sim_obj <- vast_projections
#         # what <- "Index_ctl"
#         # nice_times <- nice_times
#         # out_t_scale <- "annual"
#         # probs <- c(0.1, 0.5, 0.9)
#         # mean_instead <- FALSE
#         # nice_category_names <- nice_category_names
#         # climate_scenario <- climate_scenario
#         # out_dir <- paste0(res_root, "prediction_df")
        
#         # Capelin
#         date_dir<- here::here("2023-02-17/Capelin_BC/")
#         vast_fit = fit_full
#         sim_obj = uncert_res_full
#         resp_scale = "raw"
#         nice_times <- nice_times
#         out_t_scale = NULL
#         probs = c(0.1, 0.5, 0.9)
#         mean_instead = FALSE
#         nice_category_names = "Capelin_Random"
#         climate_scenario = climate_scenario = paste0("gfdl", "_full")
#         out_dir = date_dir

#         # ## Cod -- COG is off...
#         # date_dir <- "~/GitHub/mixedmod_projections/2022-10-25/Cod_BC/"
#         # vast_fit = readRDS(paste0(date_dir, "SpST_mod_fit.rds"))
#         # sim_obj = readRDS(file = paste0(date_dir, "SpST", "_random_ProjectionsList.rds"))
#         # resp_scale = "raw"
#         # nice_times <- as.Date(c(paste0(seq(from = 1985, to = 2100, by = 1), "-03-16"), paste0(seq(from = 1985, to = 2100, by = 1), "-07-16"), paste0(seq(from = 1985, to = 2100, by = 1), "-10-16")))
#         # nice_times <- nice_times[order(nice_times)]
#         # out_t_scale = NULL
#         # nice_category_names = paste0("Cod_", "Base_None")
#         # climate_scenario = paste0("EnvOnly_Base_5thpercentile", "_SSP5_85")
#         # out_dir = date_dir
#     }
    
#     time_ind <- seq(from = 1, to = length(nice_times))
#     time_labels <- nice_times
#     index_regions_ind <- seq(from = 1, to = vast_fit$data_list$n_l)
#     index_regions <- vast_fit$settings$strata.limits$STRATA[index_regions_ind]
#     categories_ind <- seq(from = 1, to = vast_fit$data_list$n_c)
#     grid_ind <- seq(from = 1, to = vast_fit$data_list$n_g)

#     for (i in seq_along(sim_obj)) {

#         if(FALSE){
#             # Checking sim_obj 
#             summary(sim_obj[[100]][["Index_ctl"]])
#         }

#         # Going to want the Index...
#         ind_array <- array(unlist(sim_obj[[i]][which(names(sim_obj[[i]]) ==
#             "Index_ctl")]),
#             dim = c(unlist(vast_fit$data_list[c("n_c")]), n_t = length(nice_times), unlist(vast_fit$data_list[c("n_l")])),
#             dimnames = list(
#                 categories_ind, time_labels,
#                 index_regions
#             )  
#         )
#         ind_df <- data.frame(aperm(ind_array, c(2, 3, 1)))
#         colnames(ind_df) <- gsub(".1", "", colnames(ind_df))
#         ind_df$Sim_Scenario <- rep(paste0("Sim_", i), nrow(ind_df))
#         ind_df$Time <- nice_times
#         # ind_df[, 1] <- ifelse(ind_df[, 1] == 0, NA, ind_df[, 1])
#         # ind_df <- na.omit(ind_df)
#         ind_df <- ind_df %>%
#             pivot_longer(., !c(
#                 Sim_Scenario,
#                 Time
#             ), names_to = "Region", values_to = "Index")

#         # Check
#         if(FALSE){
#             true <- vast_fit$Report$Index_ctl
#             str(true)
#             str(ind_df)
#         }
        
#         # Density
#         dens_array <- array(unlist(sim_obj[[i]][which(names(sim_obj[[i]]) ==
#             "D_gct")]),
#             dim = c(unlist(vast_fit$data_list[c("n_g")]), (vast_fit$data_list[c("n_c")]), n_t = length(nice_times)),
#             dimnames = list(grid_ind, categories_ind, time_labels)
#         )
#         dens_df <- data.frame(aperm(dens_array, c(1, 3, 2)))
#         colnames(dens_df) <- nice_times
#         dens_df$Lat <- vast_fit$spatial_list$latlon_g[, "Lat"]
#         dens_df$Lon <- vast_fit$spatial_list$latlon_g[, "Lon"]
#         dens_df <- dens_df %>%
#             pivot_longer(.,
#                 !c(Lat, Lon),
#                 names_to = "Time", values_to = "D_gct"
#             ) %>%
#             arrange(Time, Lat, Lon)
#         dens_df$Time <- as.Date(dens_df$Time)
#         dens_df$Sim_Scenario <- paste0("Sim_", i)

#         # Check
#         if(FALSE){
#             true <- vast_fit$Report$D_gct
#             str(true)
#             str(dens_df)
#         }

#         # Center of gravity
#         cog_array <- array(unlist(sim_obj[[i]][which(names(sim_obj[[i]]) ==
#             "mean_Z_ctm")]),
#         dim = c(unlist(vast_fit$data_list[c("n_c")]), n_t = length(nice_times), 2),
#         dimnames = list(
#             categories_ind, time_labels,
#             c("Lon", "Lat")
#         )
#         )
#         cog_true_df <- data.frame(aperm(cog_array, c(2, 3, 1)))
#         names(cog_true_df)<- c("Eastings", "Northings")
#         cog_true_df$Time<- nice_times
#         cog_true_df$Time <- as.Date(cog_true_df$Time)
#         cog_true_df$Sim_Scenario <- paste0("Sim_", i)

#         cog_df <- cog_from_dens(vast_fit = vast_fit, proj_dens = dens_df, proj_ind = ind_df)

#         # Check
#         if(FALSE){
#             true_lon <- vast_fit$Report$mean_Z_ctm[,,1]
#             str(true_lon)
#             true_lat<- vast_fit$Report$mean_Z_ctm[,,2]
#             str(true_lat)
#             str(cog_df)
#         }

#         # Effective area
#         eff_area_array <- array(unlist(sim_obj[[i]][which(names(sim_obj[[i]]) ==
#             "effective_area_ctl")]),
#             dim = c(unlist(vast_fit$data_list[c("n_c")]), n_t = length(nice_times), unlist(vast_fit$data_list[c("n_l")])),
#             dimnames = list(
#                 categories_ind, time_labels,
#                 index_regions
#             )  
#         )
#         eff_area_true_df <- data.frame(aperm(eff_area_array, c(2, 3, 1)))
#         colnames(eff_area_true_df) <- gsub(".1", "", colnames(eff_area_true_df))
#         eff_area_true_df$Sim_Scenario <- rep(paste0("Sim_", i), nrow(eff_area_true_df))
#         eff_area_true_df$Time <- nice_times
#         # ind_df[, 1] <- ifelse(ind_df[, 1] == 0, NA, ind_df[, 1])
#         # ind_df <- na.omit(ind_df)
#         eff_area_true_df <- eff_area_true_df %>%
#             pivot_longer(., !c(
#                 Sim_Scenario,
#                 Time
#             ), names_to = "Region", values_to = "Eff_Area")

        
#         eff_area_df <- eff_area_from_dens(vast_fit = vast_fit, proj_dens = dens_df, proj_ind = ind_df)

#         # Check
#         if (FALSE) {
#             vast_fit$Report$effective_area_ctl
#             str(true_eff_area)
#             str(eff_area_df)
#         }
        
#         if(resp_scale == "log"){
#             ind_df$Index <- log(ind_df$Index)
#             dens_df$D_gct <- log(dens_df$D_gct)
#         }
        
#         if (i == 1) {
#             res_out_ind <- ind_df
#             res_out_dens <- dens_df
#             res_out_cog <- cog_df
#             res_out_cog_true<- cog_true_df
#             res_out_eff_area <- eff_area_df
#             res_out_eff_area_true<- eff_area_true_df
#         } else {
#             res_out_ind <- bind_rows(res_out_ind, ind_df)
#             res_out_dens <- bind_rows(res_out_dens, dens_df)
#             res_out_cog <- bind_rows(res_out_cog, cog_df)
#             res_out_cog_true <- bind_rows(res_out_cog_true, cog_true_df)
#             res_out_eff_area <- bind_rows(res_out_eff_area, eff_area_df)
#             res_out_eff_area_true<- bind_rows(res_out_eff_area_true, eff_area_true_df)
#         }
#     }

#     # Calculate summaries across all runs
#     res_out_ind <- res_out_ind %>%
#         group_by(., Time, Region) %>%
#         summarise(
#             Prob_0.5 = quantile(Index, probs = 0.5, na.rm = TRUE),
#             Prob_0.1 = quantile(Index, probs = 0.1, na.rm = TRUE),
#             Prob_0.9 = quantile(Index, probs = 0.9, na.rm = TRUE)
#         )
    
#     res_out_dens <- res_out_dens %>%
#         group_by(., Lat, Lon, Time) %>%
#         summarise(
#             Prob_0.5 = quantile(D_gct, probs = 0.5, na.rm = TRUE),
#             Prob_0.1 = quantile(D_gct, probs = 0.1, na.rm = TRUE),
#             Prob_0.9 = quantile(D_gct, probs = 0.9, na.rm = TRUE)
#         )
    
#     res_out_cog <- res_out_cog %>%
#         group_by(., Time) %>%
#         summarise(
#             Lon_Prob_0.5 = quantile(Mean_Lon, probs = 0.5, na.rm = TRUE),
#             Lon_Prob_0.1 = quantile(Mean_Lon, probs = 0.1, na.rm = TRUE),
#             Lon_Prob_0.9 = quantile(Mean_Lon, probs = 0.9, na.rm = TRUE),
#             Lat_Prob_0.5 = quantile(Mean_Lat, probs = 0.5, na.rm = TRUE),
#             Lat_Prob_0.1 = quantile(Mean_Lat, probs = 0.1, na.rm = TRUE),
#             Lat_Prob_0.9 = quantile(Mean_Lat, probs = 0.9, na.rm = TRUE)
#         )
    
#     res_out_cog_true<- res_out_cog_true %>%
#         group_by(., Time) %>%
#         summarise(
#             Lon_Prob_0.5 = quantile(Eastings, probs = 0.5, na.rm = TRUE),
#             Lon_Prob_0.1 = quantile(Eastings, probs = 0.1, na.rm = TRUE),
#             Lon_Prob_0.9 = quantile(Eastings, probs = 0.9, na.rm = TRUE),
#             Lat_Prob_0.5 = quantile(Northings, probs = 0.5, na.rm = TRUE),
#             Lat_Prob_0.1 = quantile(Northings, probs = 0.1, na.rm = TRUE),
#             Lat_Prob_0.9 = quantile(Northings, probs = 0.9, na.rm = TRUE)
#         )
    
#     res_out_eff_area <- res_out_eff_area %>%
#         group_by(., Time, Region) %>%
#         summarise(
#             Prob_0.5 = quantile(Eff_Area, probs = 0.5, na.rm = TRUE),
#             Prob_0.1 = quantile(Eff_Area, probs = 0.1, na.rm = TRUE),
#             Prob_0.9 = quantile(Eff_Area, probs = 0.9, na.rm = TRUE)
#         )
    
#     res_out_eff_area_true <- res_out_eff_area_true %>%
#         group_by(., Time, Region) %>%
#         summarise(
#             Prob_0.5 = quantile(Eff_Area, probs = 0.5, na.rm = TRUE),
#             Prob_0.1 = quantile(Eff_Area, probs = 0.1, na.rm = TRUE),
#             Prob_0.9 = quantile(Eff_Area, probs = 0.9, na.rm = TRUE)
#         )
    
#     res_out<- list("Index" = res_out_ind, "Dens" = res_out_dens, "COG" = res_out_cog, "COG_True" = res_out_cog_true, "EffArea" = res_out_eff_area, "EffArea_True" = res_out_eff_area_true)
 
#     saveRDS(res_out, file = paste0(out_dir, "/", nice_category_names, 
#         "_", climate_scenario, ".rds"))
#     return(res_out)
# }
    