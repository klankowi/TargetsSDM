#####
## Making projections from fitted VAST model with new data
#####

# Overview ----------------------------------------------------------------
# The goal of this short example is to show options for making projections with a fitted VAST model and particular focus on handling estimated random effects (beta, epsilon). There are a few different reasons for doing this. First, and most generally, as far as I am aware there isn't any great "guidelines" out there for what species distribution modelers should be doing with these random effects when making forecasts or climate scale projections. This seems like an inevitable question as mixed models become the new standard and there is a priority on producing distribution and abundance forecasts and projections to increase our understanding of potential ecosystem changes on species occurrence and advise conservation and management decisions. Second, and most practically, we are facing this challenge on a specific project where the time scales and extent of interest present computational challenges such that we can't just include all of the forecast/projection years of interest and have the model estimate these effects during the model fitting process leveraging some autoregressive process (e.g., AR1 or RW). Without this computational challenge, things are practically more straight forward. We modify `rho_config` to meet our needs and then include dummy observations for each time step of interest. Internally, then, VAST (and other mixed modeling software) will estimate the random effects for those time steps. Although practically straightforward, there's still a general question there on if that is what we *should* be doing for long-term, climate scale projections. The computational issue arises as it is asking a lot for the machine to estimate N_knots * N_time steps random effects with an autoregressive structure. For instance in our work, this would amount to 400 knots * 3 seasons * 120 years.

# Preliminaries -----------------------------------------------------------
library(tidyverse)

library(VAST)

# Run for the first time? This is the seasonal wiki example
first_run <- TRUE

if (first_run) {
  # Load data and quick exploration of structure
  # Set of years and seasons. The DFO spring survey usually occurs before the NOAA NEFSC spring survey, so ordering accordingly.
  example <- load_example(data_set = "NWA_yellowtail_seasons")
  example$sampling_data <- example$sampling_data %>%
    filter(., year >= 2000)
  year_set <- sort(unique(example$sampling_data[, "year"]))
  season_set <- c("DFO", "SPRING", "FALL")

  # Create a grid with all unique combinations of seasons and years and then combine these into one "year_season" variable
  yearseason_grid <- expand.grid("season" = season_set, "year" = year_set)
  yearseason_levels <- apply(yearseason_grid[, 2:1], MARGIN = 1, FUN = paste, collapse = "_")
  yearseason_labels <- round(yearseason_grid[, "year"] + (as.numeric(factor(yearseason_grid[, "season"], levels = season_set)) - 1) / length(season_set), digits = 1)

  # Similar process, but for the observations
  yearseason_i <- apply(example$sampling_data[, c("year", "season")], MARGIN = 1, FUN = paste, collapse = "_")
  yearseason_i <- factor(yearseason_i, levels = yearseason_levels)

  # Add the year_season factor column to our sampling_data data set
  example$sampling_data$year_season <- yearseason_i
  example$sampling_data$season <- factor(example$sampling_data$season, levels = season_set)

  # Some last processing steps
  example$sampling_data <- example$sampling_data[, c("year", "season", "year_season", "latitude", "longitude", "swept", "weight")]

  # Make dummy observation for each season-year combination
  dummy_data <- data.frame(
    year = yearseason_grid[, "year"],
    season = yearseason_grid[, "season"],
    year_season = yearseason_levels,
    latitude = mean(example$sampling_data[, "latitude"]),
    longitude = mean(example$sampling_data[, "longitude"]),
    swept = mean(example$sampling_data[, "swept"]),
    weight = 0,
    dummy = TRUE
  )

  # Combine with sampling data
  full_data <- rbind(cbind(example$sampling_data, dummy = FALSE), dummy_data)

  # Create sample data
  samp_dat <- data.frame(
    "year_season" = as.numeric(full_data$year_season) - 1,
    "Lat" = full_data$latitude,
    "Lon" = full_data$longitude,
    "weight" = full_data$weight,
    "Swept" = full_data$swept,
    "Dummy" = full_data$dummy
  )

  # Covariate data. Note here, case sensitive!
  cov_dat <- data.frame(
    "Year" = as.numeric(full_data$year_season) - 1,
    "Year_Cov" = factor(full_data$year, levels = year_set),
    "Season" = full_data$season,
    "Lat" = full_data$latitude,
    "Lon" = full_data$longitude
  )

  # Inspect
  table("year_season" = cov_dat$Year, "Actual_year" = cov_dat$Year_Cov)
  table("year_season" = cov_dat$Year, "Actual_season" = cov_dat$Season)

  #####
  ## Model settings
  #####
  # Make settings
  settings <- make_settings(
    n_x = 250,
    Region = example$Region,
    strata.limits = example$strata.limits,
    purpose = "index2",
    FieldConfig = c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 1),
    RhoConfig = c("Beta1" = 3, "Beta2" = 3, "Epsilon1" = 4, "Epsilon2" = 4),
    ObsModel = c(1, 1),
    bias.correct = FALSE,
    Options = c("treat_nonencounter_as_zero" = TRUE)
  )

  # Creating model formula
  formula_use <- ~ Season + Year_Cov

  # Implement corner constraint for linear effect but not spatially varying effect:
  # * one level for each term is 2 (just spatially varying)
  # * all other levels for each term is 3 (spatialy varying plus linear effect)
  X1config_cp_use <- matrix(c(2, rep(3, nlevels(cov_dat$Season) - 1), 2, rep(3, nlevels(cov_dat$Year_Cov) - 1)), nrow = 1)
  X2config_cp_use <- matrix(c(2, rep(3, nlevels(cov_dat$Season) - 1), 2, rep(3, nlevels(cov_dat$Year_Cov) - 1)), nrow = 1)

  #####
  ## Model fit -- make sure to use new functions
  #####

  fit_orig <- fit_model(
    "settings" = settings,
    "Lat_i" = samp_dat[, "Lat"],
    "Lon_i" = samp_dat[, "Lon"],
    "t_i" = samp_dat[, "year_season"],
    "b_i" = samp_dat[, "weight"],
    "a_i" = samp_dat[, "Swept"],
    "X1config_cp" = X1config_cp_use,
    "X2config_cp" = X2config_cp_use,
    "covariate_data" = cov_dat,
    "X1_formula" = formula_use,
    "X2_formula" = formula_use,
    "X_contrasts" = list(Season = contrasts(cov_dat$Season, contrasts = FALSE), Year_Cov = contrasts(cov_dat$Year_Cov, contrasts = FALSE)),
    "run_model" = FALSE,
    "PredTF_i" = samp_dat[, "Dummy"]
  )

  # Adjust mapping for log_sigmaXi and fitting final model -- pool variance for all seasons and then set year's to NA
  Map_adjust <- fit_orig$tmb_list$Map

  # Pool variances for each term to a single value
  Map_adjust$log_sigmaXi1_cp <- factor(c(
    rep(as.numeric(Map_adjust$log_sigmaXi1_cp[1]), nlevels(cov_dat$Season)),
    rep(
      as.numeric(Map_adjust$log_sigmaXi1_cp[nlevels(cov_dat$Season) + 1]),
      nlevels(cov_dat$Year_Cov)
    )
  ))

  Map_adjust$log_sigmaXi2_cp <- factor(c(
    rep(as.numeric(Map_adjust$log_sigmaXi2_cp[1]), nlevels(cov_dat$Season)),
    rep(
      as.numeric(Map_adjust$log_sigmaXi2_cp[nlevels(cov_dat$Season) + 1]),
      nlevels(cov_dat$Year_Cov)
    )
  ))

  settings_new <- settings
  settings_new$Options["report_additional_variables"] <- TRUE
  fit <- fit_model(
    "settings" = settings_new,
    "Lat_i" = samp_dat[, "Lat"],
    "Lon_i" = samp_dat[, "Lon"],
    "t_i" = samp_dat[, "year_season"],
    "b_i" = samp_dat[, "weight"],
    "a_i" = samp_dat[, "Swept"],
    "X1config_cp" = X1config_cp_use,
    "X2config_cp" = X2config_cp_use,
    "covariate_data" = cov_dat,
    "X1_formula" = formula_use,
    "X2_formula" = formula_use,
    "X_contrasts" = list(Season = contrasts(cov_dat$Season, contrasts = FALSE), Year_Cov = contrasts(cov_dat$Year_Cov, contrasts = FALSE)),
    "newtonsteps" = 1,
    "PredTF_i" = samp_dat[, "Dummy"],
    "Map" = Map_adjust,
    "getJointPrecision" = TRUE,
    "run_model" = TRUE,
    "build_model" = TRUE
  )
}

# Addressing the challenge ----------------
# Key steps:
# 1. Fit the model with observational data, turning on options slot 16 (double check this in c++ code) and `getJointPrecision = TRUE` so that we get out the precision matrices for Q1 and Q2. Q1 and Q2 are the inverses of the covariance matrices (for first and second linear predictors?), which describe the spatial correlation among knot locations.
if (!first_run) {
  # Read in the model fit
  fit <- readRDS("~/Desktop/vast_seasonal_example.rds")
}

# Helper function for multivariate normal based on a precision matrix rather than variance:covariance matrix.
rmvnorm_prec <- function(mu, prec, n.sims, seed) {
  set.seed(seed)
  z <- matrix(rnorm(length(mu) * n.sims), ncol = n.sims)
  L <- Matrix::Cholesky(prec, super = TRUE)
  z <- Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
  z <- Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
  z <- as.matrix(z)
  return(mu + z)
}

# Padding list time elements with additional 0s for each additional time step. Assumes time dimension is the LAST dimention of the list element
pad_list_notXi_element <- function(x, pad_dim, pad_value = 0, pad_times = n_proj_steps) {
  if (FALSE) {
    x <- pred_draw$sigmaXi1_cp
    pad_dim <- 2
    pad_value <- pred_draw$sigmaXi1_cp[, ncol(pred_draw$sigmaXi1_cp)]
    pad_times <- 2
  }

  # Array?
  if (length(dim(x)) >= 3) {
    out_list_element <- unname(abind::abind(x, array(pad_value, replace(dim(x), pad_dim, pad_times)), along = pad_dim))
  }

  # Matrix?
  if (length(dim(x)) == 2) {
    if (nrow(x) < ncol(x)) {
      out_list_element <- matrix(c(x, rep(pad_value, pad_times)), nrow = 1)
    }
    if (nrow(x) == ncol(x)) {
      out_list_element <- diag(x = 1, nrow = nrow(x) + pad_times, ncol = ncol(x) + pad_times)
    }
    if (nrow(x) > ncol(x)) {
      out_list_element <- matrix(c(x, rep(pad_value, pad_times)), nrow = nrow(x) + pad_times)
    }
  }

  # Vector
  if (length(dim(x)) == 0) {
    out_list_element <- c(x, rep(pad_value, pad_times))
  }

  return(out_list_element)
}

pad_list_Xi_element <- function(x, pad_dim, pad_value, pad_start = time_cov_start, pad_times = time_cov_add) {

  # Arrays
  if (length(dim(x)) == 3) {
    # Just appending to the "end"
    if (pad_start + pad_times > dim(x)[3]) {
      out_list_element <- unname(abind::abind(array(x[, , 1:(pad_start - 1)], dim = c(dim(x)[1], dim(x)[2], (pad_start - 1))), array(pad_value, replace(dim(x), pad_dim, pad_times)), along = 3))
    } else {
      # "Inserting" work in progress...
    }
  }

  # Matrices
  if (length(dim(x)) == 2) {
    if (nrow(x) < ncol(x)) {
      # Just appending to the "end"
      if (pad_start + pad_times > ncol(x)) {
        out_list_element <- matrix(c(x[1:(pad_start - 1)], rep(pad_value, pad_times)), nrow = 1)
      } else {
        # "Inserting" work in progress...
        out_list_element <- matrix(c(x[, 1:(pad_start - 1)], rep(pad_value, pad_times), x[, (pad_value + pad_times + 1):(ncol(x) + pad_value + pad_times)]), nrow = 1)
      }
    }
  }
  return(out_list_element)
}

# Also going to need some fake "projection" data -- should have created this before fitting the seasonal model as it would have been faster.
temp_df <- fit$covariate_data[fit$covariate_data["Year_Cov"] == "2016" | fit$covariate_data["Year_Cov"] == "2017", ]
temp_df$Year_Cov <- ifelse(temp_df$Year_Cov == "2016", "2018", "2019")

# Bind it to original data to get correct factor levels
temp_df2 <- rbind.data.frame(fit$covariate_data, temp_df)

# Season and year sets
year_set <- sort(unique(as.numeric(as.character(temp_df2[, "Year_Cov"]))))
season_set <- c("DFO", "SPRING", "FALL")

# Reworking factor levels...
yearseason_grid <- expand.grid("season" = season_set, "year" = year_set)
yearseason_levels <- apply(yearseason_grid[, 2:1], MARGIN = 1, FUN = paste, collapse = "_")
yearseason_labels <- round(yearseason_grid[, "year"] + (as.numeric(factor(yearseason_grid[, "season"], levels = season_set)) - 1) / length(season_set), digits = 1)

# Similar process, but for the observations
yearseason_i <- apply(temp_df2[, c("Year_Cov", "Season")], MARGIN = 1, FUN = paste, collapse = "_")
yearseason_i <- factor(yearseason_i, levels = yearseason_levels)

# This should be the new time index...
temp_df2$Year <- as.numeric(yearseason_i) - 1

# Extract last two years. Not sure yet still on how best to provide this. For now, mimicing predict function implementation
new_projection_data <- temp_df2[temp_df2$Year > max(fit$covariate_data$Year), ]

# Designing a projection function
#' @title Project density from fitted VAST model
#'
#' @description Calculates density projections from a VAST `fit_model` object for `new_projection_data`.
#'
#' @param x = A VAST \code{fit_model} object.
#' @param time_cov = A character string or NULL defining a time covariate fitted as a spatially varying term as in the seasonal model. The current implementation can ONLY handle one of these.
#' @param category = Integer for which category to project.
#' @param what = Which output from \code{fit_model$Report} should be extracted; by default projected density.
#' @param historical_uncertainty = Character string indicating how to sample uncertainty during the historical period, where "both" samples fixed and random effects, "random" samples only random effects, and "none" assumes no uncertainty during historical period and uses MLE estimates for fixed effects and random effects conditioned on fixed effects (?)
#' @param n_samples = Integer with the number of times to sample
#' @param new_covariate_data = A dataframe with new covariate values for projection data.
#' @param new_catchability_data =
#' @param seed = 
#'
#' @return A dataframe with location and time information for each sample in `new_projection_data` and the projected \code{what} output variable.
#'
#' @export
project.fit_model <- function(x, time_cov, category = 1, what = "D_gct", historical_uncertainty = "both", n_samples, new_covariate_data, new_catchability_data = NULL, seed = 123456) {

  # For debugging
  if (FALSE) {
    x = fit
    time_cov = c("Year_Cov")
    category = 1
    what = "D_gct"
    historical_uncertainty = "both"
    n_samples= 2
    new_covariate_data <- new_projection_data
    new_catchability_data <- NULL
    seed <- 123456
  }

  #####
  ## Sanity checks
  #####

  # Check columns of new_covariate_data -- from predict function
  if (!is.null(new_covariate_data)) {
    # Confirm all columns are available
    if (!all(colnames(x$covariate_data) %in% colnames(new_covariate_data))) {
      stop("Please ensure that all columns of `x$covariate_data` are present in `new_covariate_data`")
    }
    # Eliminate unnecessary columns
    new_covariate_data <- new_covariate_data[, match(colnames(x$covariate_data), colnames(new_covariate_data))]
    # Eliminate old-covariates that are also present in new_covariate_data
    NN <- RANN::nn2(query = x$covariate_data[, c("Lat", "Lon", "Year")], data = new_covariate_data[, c("Lat", "Lon", "Year")], k = 1)
    if (any(NN$nn.dist == 0)) {
      x$covariate_data <- x$covariate_data[-which(NN$nn.dist == 0), , drop = FALSE]
    }
  }

  # Check factor levels. This is a likely stumbling block if people just supply the `new_covariate_data` without taking into account where the new "years" or "seasons" fit relative to those used in the model fitting.
  projection_levels <- lapply(new_covariate_data[sapply(new_covariate_data, is.factor)], levels)
  fit_levels <- lapply(x$covariate_data[sapply(x$covariate_data, is.factor)], levels)

  fac_level_check <- sapply(names(fit_levels), function(x) {
    all(unlist(fit_levels[x]) %in% unlist(projection_levels[x]))
  })

  if (!any(fac_level_check)) {
    stop("Please check that all factor levels in `x$covariate_data` are present in `new_covariate_data`.")
  }

  # Check if time_cov is only one string
  if (length(time_cov) != 1) {
    stop("Projection function is implemented only for one time covariate")
  }
  # Check if time_cov is not null that X1_formula and X2_formula match
  if (!is.null(time_cov) & x$X1_formula != x$X2_formula) {
    stop("Projection function for time covariates is currently implemented only for models with matching X1_formula and X2_formula")
  }

  #####
  ## Generating predictions
  #####

  # Figuring out how many projection time steps there are by comparing new_covariate_data to fit$covariate_data. Not very robust...
  n_proj <- max(new_covariate_data$Year) - max(x$covariate_data$Year)

  # Get dimensions of each of the variables, which we will need later on
  var_dims <- unlist(lapply(x$Report, FUN = function(x) length(dim(x))))

  # Get the index for time_cov from model matrix column names now if there is a time_cov, which we will need later on. Accounting for intercept would mean subtracting 1, but we actually want to start padding after the last value, so no adjustment.
  if (!is.null(time_cov)) {
    mod_mat_cols <- colnames(model.matrix(x$X1_formula, data = x$covariate_data, contrasts.arg = x$input_args$extra_args$X_contrasts))
    time_cov_start <- max(which(grepl(time_cov, mod_mat_cols)))
    time_cov_add <- length(projection_levels[[time_cov]]) - length(fit_levels[[time_cov]])
  }

  # Unpack fitted model
  Obj = x$tmb_list$Obj
  Sdreport = x$parameter_estimates$SD

  if (n_samples > 1) {

    # Going to loop over samples from the jointprecision of fitted model
    for (i in 1:n_samples) {
      if (historical_uncertainty == "both") {
        u_z = rmvnorm_prec(mu = Obj$env$last.par.best, prec = Sdreport$jointPrecision, n.sims = 1, seed = seed)[, 1]
      } else if (historical_uncertainty == "random") {
        set.seed(i)
        u_z = Obj$env$last.par.best
        MC = Obj$env$MC(keep = TRUE, n = 1, antithetic = FALSE)
        u_z[Obj$env$random] = attr(MC, "samples")[, 1]
      } else if (historical_uncertainty == "none") {
        u_z = Obj$env$last.par.best
      } else {
        stop("Check `historical_uncertainty` argument")
      }
      
      # Get ParList
      ParList = Obj$env$parList(par = u_z)

      #####
      ## Padding
      #####
      # Pad beta1_ft/beta2_ft -- pad with value from ParList? Or 0?
      pad_vars <- c("beta1_ft")
      pad_value <- ParList$beta1_ft[1, ncol(ParList$beta1_ft)]
      ParList[pad_vars] <- lapply(ParList[pad_vars], FUN = pad_list_notXi_element, pad_dim = 2, pad_value = pad_value, pad_times = n_proj)

      pad_vars <- c("beta2_ft")
      pad_value <- ParList$beta2_ft[1, ncol(ParList$beta2_ft)]
      ParList[pad_vars] <- lapply(ParList[pad_vars], FUN = pad_list_notXi_element, pad_dim = 2, pad_value = pad_value, pad_times = n_proj)

      # Pad Beta_mean1_t/mean2_t, both with 0
      pad_vars <- c("Beta_mean1_t", "Beta_mean2_t")
      ParList[pad_vars] <- lapply(ParList[pad_vars], FUN = pad_list_notXi_element, pad_dim = 1, pad_value = 0, pad_times = n_proj)

      # Pad Epsiloninput1/2_sff
      pad_vars <- c("Epsiloninput1_sff", "Epsiloninput2_sff")
      ParList[pad_vars] <- lapply(ParList[pad_vars], FUN = pad_list_notXi_element, pad_dim = 3, pad_value = 0, pad_times = n_proj)

      # Pad "time" Xivariables. We will need to account for added factor levels NOT added years as n_proj that are fitted as spatially varying effects (in this example, season and year) in:
      # gamma1/2_cp (pad with 0) [1, 1:21] 2nd dimension
      # log_sigmaXi1/2_cp (pad with estimated variance instead of 0) [1, 1:21] 2nd dimension
      # Xiinput1/2_scp (pad with 0) [1:304, 1, 1:21] 3rd dimension
      # All of them have "time" as the last dimension, which simplifies things a bit...maybe
      if (!is.null(time_covs)) {
        # gamma1/2_cp (pad with 0) [1, 1:21] 2nd dimension
        pad_vars <- c("gamma1_cp", "gamma2_cp")
        ParList[pad_vars] <- lapply(ParList[pad_vars], FUN = pad_list_Xi_element, pad_dim = 2, pad_value = 0, pad_start = time_cov_start, pad_times = time_cov_add)

        # log_sigmaXi1/2_cp (pad with parList last value, assume fixed in the future) [1, 1:21] 2nd dimension
        pad_vars <- c("log_sigmaXi1_cp")
        ParList[pad_vars] <- lapply(ParList[pad_vars], FUN = pad_list_Xi_element, pad_dim = 2, pad_value = ParList$log_sigmaXi1_cp[time_cov_start - 1], pad_start = time_cov_start, pad_times = time_cov_add)

        pad_vars <- c("log_sigmaXi2_cp")
        ParList[pad_vars] <- lapply(ParList[pad_vars], FUN = pad_list_Xi_element, pad_dim = 2, pad_value = ParList$log_sigmaXi2_cp[time_cov_start - 1], pad_start = time_cov_start, pad_times = time_cov_add)

        # Xiinput1/2 scp (pad with 0), time is the third dimension
        pad_vars <- c("Xiinput1_scp", "Xiinput2_scp")
        ParList[pad_vars] <- lapply(ParList[pad_vars], FUN = pad_list_Xi_element, pad_dim = 3, pad_value = 0, pad_start = time_cov_start, pad_times = time_cov_add)
      }

      # Done padding, could add some sort of sanity check here on dimensions.

      #####
      ## Rebuilding model
      #####
      # Processing inputs 
      t_i<- c(x$data_frame$t_i, new_covariate_data$Year)
      b_i <- c(x$data_frame$b_i, sample(c(0, mean(x$data_frame$b_i)), size = nrow(new_covariate_data), replace = TRUE))
      v_i<- c(x$data_frame$v_i, rep(0, nrow(new_covariate_data)))
      Lon_i <- c(x$data_frame$Lon_i, new_covariate_data$Lon)
      Lat_i <- c(x$data_frame$Lat_i, new_covariate_data$Lat)
      a_i<- c(x$data_frame$a_i, rep(mean(x$data_frame$a_i), nrow(new_covariate_data)))
      PredTF_i<- c(x$data_list$PredTF_i, rep(1, nrow(new_covariate_data)))
      c_iz <- rbind(x$data_list$c_iz, x$data_list$c_iz[sample(1:nrow(x$data_list$c_iz), size = nrow(new_covariate_data), replace = TRUE), , drop = FALSE])
      all_covariate_data <- rbind(x$covariate_data, new_covariate_data)
      if (is.null(new_catchability_data)) {
        all_catchability_data <- rbind(x$catchability_data, x$catchability_data[sample(1:nrow(x$data_list$c_iz), size = nrow(new_covariate_data), replace = TRUE), , drop = FALSE])
      } else {
        all_catchability_data <- rbind(x$catchability_data, new_catchability_data)
      }

      # Calling `fit_model` prompted a crash and eigen complaint.
      crash <- FALSE
      if (crash) {
        # Rebuild model
        proj_mod <- fit_model(
          settings = x$settings,
          Lat_i = Lat_i,
          Lon_i = Lon_i,
          t_i = t_i,
          b_i = b_i,
          a_i = a_i,
          v_i = v_i,
          c_iz = c_iz,
          PredTF_i = PredTF_i,
          covariate_data = new_covariate_data,
          X1_formula = x$X1_formula,
          X2_formula = x$X2_formula,
          X1config_cp = x$X1config_cp,
          X2config_cp = x$X2config_cp,
          catchability_data = new_catchability_data,
          Q1config_k = x$Q1config_k,
          Q2config_k = x$Q2config_k,
          Q1_formula = x$Q1_formula,
          Q2_formula = x$Q2_formula,
          run_model = FALSE,
          Parameters = ParList
        )
      }
      
      # Working up to using MakeADFun directly instead? Seems like we are going to need to adjust some other things -- X1/X2 config and Map??
      Xi_check <- c("Xiinput1_scp", "Xiinput2_scp")
      if (any(names(ParList) %in% Xi_check)) {
        if (Xi_check[1] %in% names(ParList)) {
          # Update X1 config...
          X1config_cp_proj <- matrix(c(2, rep(3, length(unique(all_covariate_data$Season)) - 1), 2, rep(3, nlevels(all_covariate_data$Year_Cov) - 1)), nrow = 1)
        }

        if (Xi_check[2] %in% names(ParList)) {
          # Update X2 config...
          X2config_cp_proj <- matrix(c(2, rep(3, length(unique(all_covariate_data$Season)) - 1), 2, rep(3, nlevels(all_covariate_data$Year_Cov) - 1)), nrow = 1)
        }
      } else {
        # Do nothing, use original fit Xiconfigs
        X1config_cp_proj <- x$X1config_cp
        X2config_cp_proj <- x$X2config_cp
      }
       
       # Adjusting`Map`? gamma1/gamma2 cp, log_sigmaXi1/2_cp, Xiinput1/2_scp, beta1/2_ft. Going to be different if we have additional covariates. Also still only dealing with years...
       #!!NOT ROBUST!!
       proj_map <- x$tmb_list$Map

       proj_map$gamma1_cp <- c(proj_map$gamma1_cp, seq(from = max(as.numeric(as.character(levels(proj_map$gamma1_cp)))), to = max(as.numeric(as.character(levels(proj_map$gamma1_cp)))) + length(levels(all_covariate_data$Year_Cov)[!levels(all_covariate_data$Year_Cov) %in% levels(x$covariate_data$Year_Cov)])))
       proj_map$gamma2_cp <- c(proj_map$gamma2_cp, seq(from = max(as.numeric(as.character(levels(proj_map$gamma2_cp)))), to = max(as.numeric(as.character(levels(proj_map$gamma2_cp)))) + length(levels(all_covariate_data$Year_Cov)[!levels(all_covariate_data$Year_Cov) %in% levels(x$covariate_data$Year_Cov)])))

       proj_map$log_sigmaXi1_cp <- factor(c(rep(1, length(unique(all_covariate_data$Season))), rep(4, nlevels(all_covariate_data$Year_Cov))))
       proj_map$log_sigmaXi2_cp <- factor(c(rep(1, length(unique(all_covariate_data$Season))), rep(4, nlevels(all_covariate_data$Year_Cov))))

       proj_map$Xiinput1_scp <- factor(seq(from = 1, to = dim(ParList$Xiinput1_scp)[1] * dim(ParList$Xiinput1_scp)[3]))
       proj_map$Xiinput2_scp <- factor(seq(from = 1, to = dim(ParList$Xiinput2_scp)[1] * dim(ParList$Xiinput2_scp)[3]))

       proj_map$beta1_ft <- factor(rep(1, length(unique(all_covariate_data$Year))))
       proj_map$beta2_ft <- factor(rep(1, length(unique(all_covariate_data$Year))))

       # Remaking spatial information
       message("\n### Re-making spatial information")
       spatial_args_new <- list("anisotropic_mesh" = x$spatial_list$MeshList$anisotropic_mesh, "Kmeans" = x$spatial_list$Kmeans, "Lon_i" = Lon_i, "Lat_i" = Lat_i)
       spatial_args_input <- combine_lists(input = spatial_args_new, default = x$input_args$spatial_args_input)
       spatial_list <- do.call(what = make_spatial_info, args = spatial_args_input)

      # Remaking -- from fit_model
      data_args_default = list(
        Version = x$settings$Version, FieldConfig = x$settings$FieldConfig, OverdispersionConfig = x$settings$OverdispersionConfig, RhoConfig = x$settings$RhoConfig, VamConfig = x$settings$VamConfig, ObsModel = x$settings$ObsModel, c_iz = c_iz, b_i = b_i, a_i = a_i, v_i = v_i, s_i = spatial_list$knot_i - 1, t_i = t_i, spatial_list = spatial_list, Options = x$settings$Options, Aniso = x$settings$use_anisotropy, X1config_cp = X1config_cp_proj, X2config_cp = X2config_cp_proj, covariate_data = all_covariate_data, X1_formula = x$X1_formula, X2_formula = x$X2_formula, Q1config_k = x$Q1config_k, Q2config_k = x$Q2config_k, catchability_data = all_catchability_data, Q1_formula = x$Q1_formula, Q2_formula = x$Q2_formula
      )
      extra_args <- list("X_contrasts" = list(Season = contrasts(all_covariate_data$Season, contrasts = FALSE), Year_Cov = contrasts(all_covariate_data$Year_Cov, contrasts = FALSE)))
      data_args_input = combine_lists(input = extra_args, default = data_args_default)
      data_list = do.call(what = make_data, args = data_args_input)

      message("\n### Re-making TMB object")
      model_args_default <- list("TmbData" = data_list, "RunDir" = getwd(), "Version" = x$settings$Version, "RhoConfig" = x$settings$RhoConfig, "loc_x" = spatial_list$loc_x, "Method" = spatial_list$Method, "Map" = x$tmb_list$Map)
      model_args_input <- combine_lists(
        input = list("Parameters" = par_list_new, "Map" = proj_map, "Random" = NULL, "run_model" = FALSE, "build_model" = TRUE, "RunDir" = x$input_args$model_args_input$RunDir),
        default = model_args_default, args_to_use = formalArgs(make_model)
      )
      tmb_list <- do.call(what = make_model, args = model_args_input)

      #####
      ## Simulating with the new model
      #####
      # May need to change above this as simulate_data wants Obj from `fit$tmb_list$Obj`. Provided we figure that out, then we call simulate. Which type would we use? type = 3 as this will simulate the random effects values for the time steps we padded??
      sim_temp <- simulate_data(
        fit = tmb_list$Obj,
        type = 2,
        random_seed = seed
      )

      # From that, we are going to want to...save everything? Save just our derived values of interest (e.g., D_gct or D_i?)
      if (i == 1) {
        res <- res_to_save_from_sim_temp
      } else {
        # Bind them?
        res <- bind_somehow(res, res_to_save_from_sim_temp)
      }
      
      # Print update
      message('Sample ', i, ' out of ' n_samples, ' is done')
    }









    




















      # Draw fram pedictive distribution
      pred_draw <- sample_variable(
        Sdreport = x$parameter_estimates$SD,
        Obj = x$tmb_list$Obj,
        rmvnorm_mu = Obj$env$last.par.best,
        variable_name = NULL,
        n_samples = 1,
        seed = i,
        sample_fixed <- TRUE
      )

      # Pad Epsilon1/2 g/s_t array variables that are indexed by "t" in third dimension
      pad_vars <- c("Epsilon1_gct", "Epsilon1_sct", "Epsiloninput1_sft", "eta1_gct", "P1_gct", "R1_gct", "Epsilon2_gct", "Epsilon2_sct", "Epsiloninput2_sft", "eta2_gct", "P2_gct", "R2_gct")
      pred_draw[pad_vars] <- lapply(pred_draw[pad_vars], FUN = pad_list_notXi_element, pad_dim = 3, pad_value = 0, pad_times = n_proj_steps)

      # Pad matrix variables that are index by "t" Ltime_epsilon1/2_tf [n_t,n_t]
      pad_vars <- c("Ltime_epsilon1_tf", "Ltime_epsilon2_tf")
      pred_draw[pad_vars] <- lapply(pred_draw[pad_vars], FUN = pad_list_notXi_element, pad_dim = 2, pad_value = NULL, pad_times = n_proj_steps)

      # Pad index_ctl
      pad_vars <- c("Index_ctl")
      pred_draw[pad_vars] <- lapply(pred_draw[pad_vars], FUN = pad_list_notXi_element, pad_dim = 2, pad_value = 0, pad_times = n_proj_steps)

      # Pad index_gctl
      pad_vars <- c("Index_gctl")
      pred_draw[pad_vars] <- lapply(pred_draw[pad_vars], FUN = pad_list_notXi_element, pad_dim = 3, pad_value = 0, pad_times = n_proj_steps)

      # Pad Beta_mean1_t, Beta_mean2_t
      pad_vars <- c("Beta_mean1_t", "Beta_mean2_t")
      pred_draw[pad_vars] <- lapply(pred_draw[pad_vars], FUN = pad_list_notXi_element, pad_dim = 1, pad_value = 0, pad_times = n_proj_steps)

      # Pad beta1_mean_tf/beta2_mean_tf
      pad_vars <- c("beta1_mean_tf", "beta2_mean_tf")
      pred_draw[pad_vars] <- lapply(pred_draw[pad_vars], FUN = pad_list_notXi_element, pad_dim = 2, pad_value = 0, pad_times = n_proj_steps)

      # Pad beta1_tc
      pad_vars <- c("beta1_tc")
      pred_draw[pad_vars] <- lapply(pred_draw[pad_vars], FUN = pad_list_notXi_element, pad_dim = 2, pad_value = 2.48, pad_times = n_proj_steps)

      # Pad beta2_tc
      pad_vars <- c("beta2_tc")
      pred_draw[pad_vars] <- lapply(pred_draw[pad_vars], FUN = pad_list_notXi_element, pad_dim = 2, pad_value = 1.42, pad_times = n_proj_steps)

      # Pad iota_ct
      pad_vars <- c("iota_ct")
      pred_draw[pad_vars] <- lapply(pred_draw[pad_vars], FUN = pad_list_notXi_element, pad_dim = 2, pad_value = 0, pad_times = n_proj_steps)

      # Pad "time" Xivariables. We will need to account for added factor levels NOT added years as n_proj_steps that are fitted as spatially varying effects (in this example, season and year) in:
      # sigmaXi1/2_cp (pad with estimated variance instead of 0) [1, 1:21] 2nd dimension
      # Xi1/2_gcp (pad with 0) [1:2500, 1, 1:21] 3rd dimension
      # Xi1/2_scp (pad with 0) [1:304, 1, 1:21] 3rd dimension
      # All of them have "time" as the last dimension, which simplifies things a bit...maybe
      if (!is.null(time_covs)) {
        # sigmaXi1/2_cp (pad with estimated variance instead, assume fixed in the future) [1, 1:21] 2nd dimension
        pad_vars <- c("sigmaXi1_cp")
        pred_draw[pad_vars] <- lapply(pred_draw[pad_vars], FUN = pad_list_Xi_element, pad_dim = 2, pad_value = x$Report$sigmaXi1_cp[time_cov_start - 1], pad_start = time_cov_start, pad_times = time_cov_add)

        pad_vars <- c("sigmaXi2_cp")
        pred_draw[pad_vars] <- lapply(pred_draw[pad_vars], FUN = pad_list_Xi_element, pad_dim = 2, pad_value = x$Report$sigmaXi2_cp[time_cov_start - 1], pad_start = time_cov_start, pad_times = time_cov_add)

        # Xi1/2 gcp/scp (pad with 0), time is the third dimension
        pred_vars <- c("Xi1_scp", "Xi1_gcp", "Xi2_scp", "Xi2_gcp")
        pred_draw[pad_vars] <- lapply(pred_draw[pad_vars], FUN = pad_list_Xi_element, pad_dim = 3, pad_value = 0, pad_start = time_cov_start, pad_times = time_cov_add)
      }

      # Done padding, could add some sort of sanity check here?

      # I'm a bit confused here between working with report and working with parList. This was mentioned briefly, though for some it sounds like there is a 1:1 match directly. For others, this is not the case.
    }
  }

  ## If n_samples = 1 OR n_samples = NULL, then we do something else
  if (n_samples == 1) {
    # Just run sample_variable once, pad time variables, run makeADFun, and simulate
  }

  if (is.null(n_samples)) {
    # This option would use the MLE and never run sample_variable, pad time variables, run makeADFun and simulate
  }
}






  # Generate a new `parList` object that will have additional storage to account for the projection time steps.
  par_list_new <- x$tmb_list$Obj$env$parList()

  #####
  ## New implementation
  #####
  # Figuring out how many projection time steps there are by comparing new_covariate_data to fit$covariate_data
  n_proj_steps <- max(new_covariate_data$Year) - max(x$covariate_data$Year)

  # Sampling gamma1/2, which are not in `x$Report()`? These are going to be main effects for each season and year, though could also include habitat covariates. No additional seasons in this case, so just need to draw main year effects for "new" years -- what about something like `  gamma1_samps <- rnorm(n_gammas_add, mean = unique(x$Report[["beta1_tc"]]), sd = x$Report$sigmaXi1_cp[4])`?

  # Determining which variables we need to sample -- must be a creative way of doing this -- check slots if they are empty or not (e.g., sapply(par_list_new, is_empty))? Names difference from par_list_new and what is expected by sample_variable, which uses names from `Report`. Also do we want to be sampling "s" or "g" -- depend on settings$fine_scale -- and "c" or "f" -- depend on n_c? Finally, beta pieces? Hardwiring this for now
  variables_df <- data.frame("report_name" = c("Beta_mean1_t", "Beta_mean2_t", "sigmaXi1_cp", "sigmaXi2_cp", "Epsiloninput1_sft", "Epsiloninput2_sft", "Xi1_scp", "Xi2_scp"), "par_list_new_name" = c("Beta_mean1_t", "Beta_mean2_t", "log_sigmaXi1_cp", "log_sigmaXi2_cp", "Epsiloninput1_sff", "Epsiloninput2_sff", "Xi1_scp", "Xi2_scp"))

  # Starting sample -- does this account for the AR structure internally?
  pred_samps <- sample_variables(
    Sdreport <- x$parameter_estimates$SD,
    Obj <- x$tmb_list$Obj,
    rmvnorm_mu = Obj$env$last.par.best,
    variable_names <- variables_df$report_name,
    n_samples <- 1,
    seed <- seed,
    sample_fixed <- TRUE
  )

  # For each projection time step, run new `sample_variables` function. How would this change if we have AR structure?? Do we want to bootstrap this?
  for (i in 1:n_proj_steps) {
    pred_samps <- sample_variables(
      Sdreport <- x$parameter_estimates$SD,
      Obj <- x$tmb_list$Obj,
      rmvnorm_mu = last_sample,
      variable_names <- variables_df$report_name,
      n_samples <- 1,
      seed <- seed,
      sample_fixed <- TRUE
    )

    # Append values to par_list_new
    # Add that here
  }

  # Bring together fitted covariate data and new_covariate_data
  all_covariate_data <- rbind(x$covariate_data, new_covariate_data)

  # Along with updating `parList`, we also are going to need to change some of the other components to the model structure !! COME BACK TO THIS TO MAKE IT MORE ROBUST !!
  # X1/X2config parts
  Xi_check <- c("Xiinput1_scp", "Xiinput2_scp")
  if (any(res_fit$par_list_name %in% Xi_check)) {
    if (Xi_check[1] %in% res_fit$par_list_name) {
      # Update X1 config...
      X1config_cp_proj <- matrix(c(2, rep(3, length(unique(all_covariate_data$Season)) - 1), 2, rep(3, nlevels(all_covariate_data$Year_Cov) - 1)), nrow = 1)
    }

    if (Xi_check[2] %in% res_fit$par_list_name) {
      # Update X2 config...
      X2config_cp_proj <- matrix(c(2, rep(3, length(unique(all_covariate_data$Season)) - 1), 2, rep(3, nlevels(all_covariate_data$Year_Cov) - 1)), nrow = 1)
    }
  } else {
    # Do nothing, use original fit Xiconfigs
    X1config_cp_proj <- x$X1config_cp
    X2config_cp_proj <- x$X2config_cp
  }

  # Adjusting `Map`? slopts
  proj_map <- x$tmb_list$Map

  proj_map$gamma1_cp <- c(proj_map$gamma1_cp, seq(from = max(as.numeric(as.character(levels(proj_map$gamma1_cp)))), to = max(as.numeric(as.character(levels(proj_map$gamma1_cp)))) + length(levels(new_projection_data$Year_Cov)[!levels(new_projection_data$Year_Cov) %in% levels(x$covariate_data$Year_Cov)])))
  proj_map$gamma2_cp <- c(proj_map$gamma2_cp, seq(from = max(as.numeric(as.character(levels(proj_map$gamma2_cp)))), to = max(as.numeric(as.character(levels(proj_map$gamma2_cp)))) + length(levels(new_projection_data$Year_Cov)[!levels(new_projection_data$Year_Cov) %in% levels(x$covariate_data$Year_Cov)])))

  proj_map$log_sigmaXi1_cp <- factor(c(rep(1, length(unique(new_projection_data$Season))), rep(4, nlevels(new_projection_data$Year_Cov))))
  proj_map$log_sigmaXi2_cp <- factor(c(rep(1, length(unique(new_projection_data$Season))), rep(4, nlevels(new_projection_data$Year_Cov))))

  proj_map$Xiinput1_scp <- factor(seq(from = 1, to = dim(par_list_new$Xiinput1_scp)[1] * dim(par_list_new$Xiinput1_scp)[3]))
  proj_map$Xiinput2_scp <- factor(seq(from = 1, to = dim(par_list_new$Xiinput2_scp)[1] * dim(par_list_new$Xiinput2_scp)[3]))

  proj_map$beta1_ft <- factor(rep(1, length(unique(new_projection_data$Year))))
  proj_map$beta2_ft <- factor(rep(1, length(unique(new_projection_data$Year))))

  # With the added values in `parList`, adjustments to config and map objects, rebuild the model.
  # Process inputs
  PredTF_i <- c(x$data_list$PredTF_i, rep(1, length(t_proj_i)))
  b_i <- c(x$data_frame[, "b_i"], sample(c(0, 1), size = length(t_proj_i), replace = TRUE))
  c_iz <- rbind(matrix(x$data_frame[, grep("c_iz", names(x$data_frame))]), matrix(c_proj_iz))
  Lat_i <- c(x$data_frame[, "Lat_i"], Lat_proj_i)
  Lon_i <- c(x$data_frame[, "Lon_i"], Lon_proj_i)
  a_i <- c(x$data_frame[, "a_i"], a_proj_i)
  v_i <- c(x$data_frame[, "v_i"], v_proj_i)
  t_i <- c(x$data_frame[, "t_i"], t_proj_i)

  # # Build information regarding spatial location and correlation
  message("\n### Re-making spatial information")
  spatial_args_new <- list("anisotropic_mesh" = x$spatial_list$MeshList$anisotropic_mesh, "Kmeans" = x$spatial_list$Kmeans, "Lon_i" = Lon_i, "Lat_i" = Lat_i)
  spatial_args_input <- combine_lists(input = spatial_args_new, default = x$input_args$spatial_args_input)
  spatial_list <- do.call(what = make_spatial_info, args = spatial_args_input)

  # # Check spatial_list
  if (!all.equal(spatial_list$MeshList, x$spatial_list$MeshList)) {
    stop("`MeshList` generated during `predict.fit_model` doesn't match that of original fit; please email package author to report issue")
  }

  # # Build data
  # # Do *not* restrict inputs to formalArgs(make_data) because other potential inputs are still parsed by make_data for backwards compatibility
  # message("\n### Re-making data object")
  data_args_new <- list(
    "c_iz" = c_iz, "b_i" = b_i, "a_i" = a_i, "v_i" = v_i, "PredTF_i" = PredTF_i, "t_i" = t_i, "spatial_list" = spatial_list, "covariate_data" = all_covariate_data, "catchability_data" = catchability_data, "X_contrasts" = list(Season = contrasts(all_covariate_data$Season, contrasts = FALSE), Year_Cov = contrasts(all_covariate_data$Year_Cov, contrasts = FALSE)), "X1config_cp" = X1config_cp_proj, "X2config_cp" = X2config_cp_proj
  )
  data_args_input <- combine_lists(input = data_args_new, default = x$input_args$data_args_input) # Do *not* use args_to_use
  data_list <- do.call(what = make_data, args = data_args_input)
  data_list$n_g <- 0

  # Doing some checking here for potential indexing issues...
  n_is <- data_list$n_i
  n_ss <- data_list$n_s
  n_gs <- data_list$n_g
  n_ts <- data_list$n_t
  n_bis <- length(data_list$b_i)
  n_tis <- length(data_list$t_i)
  PredTF_is <- length(data_list$PredTF_i)

  ncol((data_list$X1config_cp)) == length(levels(all_covariate_data$Season)) + length(levels(all_covariate_data$Year_Cov))
  str(data_list$a_gl) # Not sure about this one
  nrow(data_list$X1_ip) == n_is & ncol(data_list$X1_ip) == length(levels(all_covariate_data$Season)) + length(levels(all_covariate_data$Year_Cov))
  dim(data_list$X1_gctp)[3] == n_ts & dim(data_list$X1_gctp)[4] == length(levels(all_covariate_data$Season)) + length(levels(all_covariate_data$Year_Cov))

  str(data_list$Ais_ij)

  # Build object -- ERROR, something about TMB has received an error from eigen and then it crashes...
  message("\n### Re-making TMB object")
  model_args_default <- list("TmbData" = data_list, "RunDir" = getwd(), "Version" = x$settings$Version, "RhoConfig" = x$settings$RhoConfig, "loc_x" = spatial_list$loc_x, "Method" = spatial_list$Method, "Map" = x$tmb_list$Map)
  model_args_input <- combine_lists(
    input = list("Parameters" = par_list_new, "Map" = proj_map, "Random" = NULL, "run_model" = FALSE, "build_model" = TRUE, "RunDir" = x$input_args$model_args_input$RunDir),
    default = model_args_default, args_to_use = formalArgs(make_model)
  )
  tmb_list <- do.call(what = make_model, args = model_args_input)










  # Start with the relatively simple gamma1/gamma2 bits?
  # QUESTION: SHOULD I INSTEAD LOOK AT MODIFYING `SAMPLE_VARIABLE` TO SAMPLE ALL OF THEM?
  # How would we address when gamma includes other covariates? This will influence the log_sigmaXi indexing and where the gamma values get added in.
  # How should we generate these? Normal draw using mean beta1_tc (expected intercepts among years?) and sigma values as SD? Should this be from Report or par_list_new?
  # !!!! COME  BACK TO THIS AS ITS HARDWIRED NOW !!!!
  n_gammas_add <- length(unlist(projection_levels["Year_Cov"])[-which(unlist(projection_levels["Year_Cov"]) %in% unlist(fit_levels["Year_Cov"]))])

  gamma1_samps <- rnorm(n_gammas_add, mean = unique(x$Report[["beta1_tc"]]), sd = x$Report$sigmaXi1_cp[4])
  par_list_new[[{{ "gamma1_cp" }}]] <- matrix(c(par_list_new[[{{ "gamma1_cp" }}]], gamma1_samps), nrow = 1)

  gamma2_samps <- rnorm(n_gammas_add, mean = unique(x$Report[["beta2_tc"]]), sd = x$Report$sigmaXi2_cp[4])
  par_list_new[[{{ "gamma1_cp" }}]] <- matrix(c(par_list_new[[{{ "gamma2_cp" }}]], gamma2_samps), nrow = 1)

  # log_sigmaXi parts
  par_list_new[[{{ "log_sigmaXi1_cp" }}]] <- matrix(c(par_list_new[[{{ "log_sigmaXi1_cp" }}]], rep(par_list_new[[{{ "log_sigmaXi1_cp" }}]][4], n_gammas_add)), nrow = 1)
  par_list_new[[{{ "log_sigmaXi2_cp" }}]] <- matrix(c(par_list_new[[{{ "log_sigmaXi2_cp" }}]], rep(par_list_new[[{{ "log_sigmaXi2_cp" }}]][4], n_gammas_add)), nrow = 1)

  # Now, moving onto the more complicated random effects. There's likely a clever way of figuring out which of these we need to deal with and it will depend on the model swttings -- check slots if they are empty or not (e.g., sapply(par_list_new, is_empty))? Names difference from par_list_new and what is expected by sample_variable, which uses names from `Report`. Also switch depending on if extrapolation grid vs. knots?
  # !!!! COME  BACK TO THIS AS ITS HARDWIRED NOW !!!!

  # Figuring out how many projection time steps there are by comparing new_covariate_data to fit$covariate_data
  n_proj_steps <- max(new_covariate_data$Year) - max(x$covariate_data$Year)

  re_match_df <- data.frame("par_list_name" = c("Epsiloninput1_sff", "Epsiloninput2_sff", "Xiinput1_scp", "Xiinput2_scp"), "sample_variable_name" = c("Epsiloninput1_sft", "Epsiloninput2_sft", "Xi1_scp", "Xi2_scp"))
  res_fit <- re_match_df[1:4, ]

  # Next, get precision matrices...
  Q1 <- x$Report$Q1
  Q2 <- x$Report$Q2

  # Likely a much better way of doing this next part, but my brain thinks first in loops. So, for any of these random effects in the original fit, we are going to want to append values for each of the projection time steps.
  for (i in 1:nrow(res_fit)) {

    # Get random effect name info
    re_use <- res_fit[i, ]

    # Add storage to that re effect slot in par_list_new.
    par_list_new[[{{ re_use$par_list_name }}]] <- unname(abind::abind(par_list_new[[{{ re_use$par_list_name }}]], array(0, replace(dim(par_list_new[[{{ re_use$par_list_name }}]]), 3, n_proj_steps)), along = 3))

    # Run the `sample_variable` function to get one sample from the predictive distribution for the `re_use` random effect.
    # This is a bit odd, as samples across res are going to be independent?
    post_samps <- sample_variable(
      Sdreport = x$parameter_estimates$SD,
      Obj = x$tmb_list$Obj,
      variable_name = re_use$sample_variable_name,
      n_samples = 1,
      seed = seed,
      sample_fixed = TRUE
    )

    # Going to keep the last one and use this as the mean for the multivariate normal distribution?
    # Do we even need to do this? Could we just use the last fitted value instead (e.g., fit$Report$Epsiloninput1_sft[,category, dim(fit$Report$Epsiloninput1_sft)[3]])?
    mu_use <- post_samps[, category, dim(post_samps)[3]]

    # Inner loop over n_proj_steps
    for (j in 1:n_proj_steps) {

      # Multivariate normal draw. Precision matrix depends on first or second linear predictor...
      if (grepl("1", re_use$par_list_name)) {
        mv_samp <- rmvnorm_prec(mu = mu_use, prec = Q1, n.sims = 1, seed = seed)
      } else {
        mv_samp <- rmvnorm_prec(mu = mu_use, prec = Q2, n.sims = 1, seed = seed)
      }

      # I'm a bit lost here when thinking about adding in the autoregressive component if needed. Looking at the cpp code, it looks like this comes in for Epsilon1 (as an example) as `Epsilon_rho2_f(f) * Epsiloninput2_sft(s,f,t-1)`. It seems like that structure would then remove the need to do the MV draws -- the time series of epsilon values would just move forward from the last fitted surface and Epsilon_rho_1_f? Skipping this for now.

      # Appending to par_list_new object. This is going to be different for Xi1 and Xi2 as these have a third dimension of parameter rather than time. Only adding new years -- would be sensitive to if model formula had year first or second? Come back to
      if (grepl("sff", re_use$par_list_name) | grepl("sft", re_use$par_list_name)) {
        par_list_new[[re_use$par_list_name]][, category, x$data_list$n_t + j] <- mv_samp
      } else {
        if (grepl("1", re_use$par_list_name)) {
          par_list_new[[re_use$par_list_name]][, category, dim(x$X1config_cp)[2] + j] <- mv_samp
        } else {
          par_list_new[[re_use$par_list_name]][, category, dim(x$X2config_cp)[2] + j] <- mv_samp
        }
      }

      # End of inner loop over n_proj_steps
    }

    # End of outer loop over random effects
  }

  # Some checks here...
  length(par_list_new$gamma1_cp) - length(x$tmb_list$Obj$env$parList()$gamma1_cp) == n_gammas_add
  dim(par_list_new$Epsiloninput1_sff)[3] - dim(x$tmb_list$Obj$env$parList()$Epsiloninput1_sff)[3] == n_proj_steps

  # Bring together fitted covariate data and new_covariate_data
  all_covariate_data <- rbind(x$covariate_data, new_covariate_data)

  # Along with updating `parList`, we also might need to change the X1 or X2 configs. !! COME BACK TO THIS TO MAKE IT MORE ROBUST !!
  Xi_check <- c("Xiinput1_scp", "Xiinput2_scp")
  if (any(res_fit$par_list_name %in% Xi_check)) {
    if (Xi_check[1] %in% res_fit$par_list_name) {
      # Update X1 config...
      X1config_cp_proj <- matrix(c(2, rep(3, length(unique(all_covariate_data$Season)) - 1), 2, rep(3, nlevels(all_covariate_data$Year_Cov) - 1)), nrow = 1)
    }

    if (Xi_check[2] %in% res_fit$par_list_name) {
      # Update X2 config...
      X2config_cp_proj <- matrix(c(2, rep(3, length(unique(all_covariate_data$Season)) - 1), 2, rep(3, nlevels(all_covariate_data$Year_Cov) - 1)), nrow = 1)
    }
  } else {
    # Do nothing, use original fit Xiconfigs
    X1config_cp_proj <- x$X1config_cp
    X2config_cp_proj <- x$X2config_cp
  }

  # Finally, adjust `Map`? gamma1/gamma2 cp, log_sigmaXi1/2_cp, Xiinput1/2_scp, beta1/2_ft. Going to be different if we have additional covariates. Also still only dealing with years...
  proj_map <- x$tmb_list$Map

  proj_map$gamma1_cp <- c(proj_map$gamma1_cp, seq(from = max(as.numeric(as.character(levels(proj_map$gamma1_cp)))), to = max(as.numeric(as.character(levels(proj_map$gamma1_cp)))) + length(levels(new_projection_data$Year_Cov)[!levels(new_projection_data$Year_Cov) %in% levels(x$covariate_data$Year_Cov)])))
  proj_map$gamma2_cp <- c(proj_map$gamma2_cp, seq(from = max(as.numeric(as.character(levels(proj_map$gamma2_cp)))), to = max(as.numeric(as.character(levels(proj_map$gamma2_cp)))) + length(levels(new_projection_data$Year_Cov)[!levels(new_projection_data$Year_Cov) %in% levels(x$covariate_data$Year_Cov)])))

  proj_map$log_sigmaXi1_cp <- factor(c(rep(1, length(unique(new_projection_data$Season))), rep(4, nlevels(new_projection_data$Year_Cov))))
  proj_map$log_sigmaXi2_cp <- factor(c(rep(1, length(unique(new_projection_data$Season))), rep(4, nlevels(new_projection_data$Year_Cov))))

  proj_map$Xiinput1_scp <- factor(seq(from = 1, to = dim(par_list_new$Xiinput1_scp)[1] * dim(par_list_new$Xiinput1_scp)[3]))
  proj_map$Xiinput2_scp <- factor(seq(from = 1, to = dim(par_list_new$Xiinput2_scp)[1] * dim(par_list_new$Xiinput2_scp)[3]))

  proj_map$beta1_ft <- factor(rep(1, length(unique(new_projection_data$Year))))
  proj_map$beta2_ft <- factor(rep(1, length(unique(new_projection_data$Year))))

  # With the added values in `parList`, adjustments to config and map objects, rebuild the model.
  # Process inputs
  PredTF_i <- c(x$data_list$PredTF_i, rep(1, length(t_proj_i)))
  b_i <- c(x$data_frame[, "b_i"], sample(c(0, 1), size = length(t_proj_i), replace = TRUE))
  c_iz <- rbind(matrix(x$data_frame[, grep("c_iz", names(x$data_frame))]), matrix(c_proj_iz))
  Lat_i <- c(x$data_frame[, "Lat_i"], Lat_proj_i)
  Lon_i <- c(x$data_frame[, "Lon_i"], Lon_proj_i)
  a_i <- c(x$data_frame[, "a_i"], a_proj_i)
  v_i <- c(x$data_frame[, "v_i"], v_proj_i)
  t_i <- c(x$data_frame[, "t_i"], t_proj_i)

  # # Build information regarding spatial location and correlation
  message("\n### Re-making spatial information")
  spatial_args_new <- list("anisotropic_mesh" = x$spatial_list$MeshList$anisotropic_mesh, "Kmeans" = x$spatial_list$Kmeans, "Lon_i" = Lon_i, "Lat_i" = Lat_i)
  spatial_args_input <- combine_lists(input = spatial_args_new, default = x$input_args$spatial_args_input)
  spatial_list <- do.call(what = make_spatial_info, args = spatial_args_input)

  # # Check spatial_list
  if (!all.equal(spatial_list$MeshList, x$spatial_list$MeshList)) {
    stop("`MeshList` generated during `predict.fit_model` doesn't match that of original fit; please email package author to report issue")
  }

  # # Build data
  # # Do *not* restrict inputs to formalArgs(make_data) because other potential inputs are still parsed by make_data for backwards compatibility
  # message("\n### Re-making data object")
  data_args_new <- list(
    "c_iz" = c_iz, "b_i" = b_i, "a_i" = a_i, "v_i" = v_i, "PredTF_i" = PredTF_i, "t_i" = t_i, "spatial_list" = spatial_list, "covariate_data" = all_covariate_data, "catchability_data" = catchability_data, "X_contrasts" = list(Season = contrasts(all_covariate_data$Season, contrasts = FALSE), Year_Cov = contrasts(all_covariate_data$Year_Cov, contrasts = FALSE)), "X1config_cp" = X1config_cp_proj, "X2config_cp" = X2config_cp_proj
  )
  data_args_input <- combine_lists(input = data_args_new, default = x$input_args$data_args_input) # Do *not* use args_to_use
  data_list <- do.call(what = make_data, args = data_args_input)
  data_list$n_g <- 0

  # Doing some checking here for potential indexing issues...
  n_is <- data_list$n_i
  n_ss <- data_list$n_s
  n_gs <- data_list$n_g
  n_ts <- data_list$n_t
  n_bis <- length(data_list$b_i)
  n_tis <- length(data_list$t_i)
  PredTF_is <- length(data_list$PredTF_i)

  ncol((data_list$X1config_cp)) == length(levels(all_covariate_data$Season)) + length(levels(all_covariate_data$Year_Cov))
  str(data_list$a_gl) # Not sure about this one
  nrow(data_list$X1_ip) == n_is & ncol(data_list$X1_ip) == length(levels(all_covariate_data$Season)) + length(levels(all_covariate_data$Year_Cov))
  dim(data_list$X1_gctp)[3] == n_ts & dim(data_list$X1_gctp)[4] == length(levels(all_covariate_data$Season)) + length(levels(all_covariate_data$Year_Cov))

  str(data_list$Ais_ij)

  # Build object -- ERROR, something about TMB has received an error from eigen and then it crashes...
  message("\n### Re-making TMB object")
  model_args_default <- list("TmbData" = data_list, "RunDir" = getwd(), "Version" = x$settings$Version, "RhoConfig" = x$settings$RhoConfig, "loc_x" = spatial_list$loc_x, "Method" = spatial_list$Method, "Map" = x$tmb_list$Map)
  model_args_input <- combine_lists(
    input = list("Parameters" = par_list_new, "Map" = proj_map, "Random" = NULL, "run_model" = FALSE, "build_model" = TRUE, "RunDir" = x$input_args$model_args_input$RunDir),
    default = model_args_default, args_to_use = formalArgs(make_model)
  )
  tmb_list <- do.call(what = make_model, args = model_args_input)

  # Extract category and variable from report
}