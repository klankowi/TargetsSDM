#####
## Making projections from fitted VAST model with new data
#####

# Overview ----------------------------------------------------------------
# The goal of this short example is to show options for making projections with a fitted VAST model and particular focus on handling estimated random effects (beta, epsilon). There are a few different reasons for doing this. First, and most generally, as far as I am aware there isn't any great "guidelines" out there for what species distribution modelers should be doing with these random effects when making forecasts or climate scale projections. This seems like an inevitable question as mixed models become the new standard and there is a priority on producing distribution and abundance forecasts and projections to increase our understanding of potential ecosystem changes on species occurrence and advise conservation and management decisions. Second, and most practically, we are facing this challenge on a specific project where the time scales and extent of interest present computational challenges such that we can't just include all of the forecast/projection years of interest and have the model estimate these effects during the model fitting process leveraging some autoregressive process (e.g., AR1 or RW). Without this computational challenge, things are practically more straight forward. We modify `rho_config` to meet our needs and then include dummy observations for each time step of interest. Internally, then, VAST (and other mixed modeling software) will estimate the random effects for those time steps. Although practically straightforward, there's still a general question there on if that is what we *should* be doing for long-term, climate scale projections. The computational issue arises as it is asking a lot for the machine to estimate N_knots * N_time steps random effects with an autoregressive structure. For instance in our work, this would amount to 400 knots * 3 seasons * 120 years.

# Preliminaries -----------------------------------------------------------
library(VAST)

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
  rep(as.numeric(Map_adjust$log_sigmaXi1_cp[nlevels(cov_dat$Season) + 1]), nlevels(cov_dat$Year_Cov))
))
Map_adjust$log_sigmaXi2_cp <- factor(c(
  rep(as.numeric(Map_adjust$log_sigmaXi2_cp[1]), nlevels(cov_dat$Season)),
  rep(as.numeric(Map_adjust$log_sigmaXi2_cp[nlevels(cov_dat$Season) + 1]), nlevels(cov_dat$Year_Cov))
))

# Addressing the challenge ----------------
# Key steps:
# 1. Fit the model with observational data, turning on options slot 16 (double check this in c++ code) and `getJointPrecision = TRUE` so that we get out the precision matrices for Q1 and Q2. Q1 and Q2 are the inverses of the covariance matrices (for first and second linear predictors?), which describe the spatial correlation among knot locations.
# Fit final model with new mapping
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

# 2. Set up storage `ParList` object, which should have same structure as the #1, with added storage in epsilon, Xi slots where added storage depends on the projection time steps (e.g., if only projecting one year, epsilon slot would have an additional knots * categories * time (1) of storage).
n_proj_years <- 2 # User input

# Internal bits of a function...
n_knots <- fit$data_list$n_s
n_cats <- fit$data_list$n_c

par_list_new <- fit$tmb_list$Parameters
str(par_list_new)

# Add place holders for simulated/approximated epsilon values. Come back to this, could we make it more itneractive and have it check settings and then modify things as needed?
par_list_new[["Epsiloninput1_sff"]] <- abind(par_list_new[["Epsiloninput1_sff"]], array(0, replace(dim(par_list_new[["Epsiloninput1_sff"]]), 3, n_proj_years)), along = 3)
par_list_new[["Epsiloninput2_sff"]] <- abind(par_list_new[["Epsiloninput2_sff"]], array(0, replace(dim(par_list_new[["Epsiloninput2_sff"]]), 3, n_proj_years)), along = 3)

# Add place holders for Xi's and account for spatially-varying year effect. Come back to this, could we make ti more interactive and have it check Xconfig for "3"'s.'
par_list_new[["Xiinput1_scp"]] <- abind(par_list_new[["Xiinput1_scp"]], array(0, replace(dim(par_list_new[["Xiinput1_scp"]]), 3, n_proj_steps)), along = 3)
par_list_new[["Xiinput2_scp"]] <- abind(par_list_new[["Xiinput2_scp"]], array(0, replace(dim(par_list_new[["Xiinput2_scp"]]), 3, n_proj_steps)), along = 3)

# Do we need beta??

# Then if doing a single time steps projection:
# 3. Run the `sample_variable` function to get one sample from the predictive distribution for fixed and random effects,
post_samps <- sample_variable(
  Sdreport = fit$parameter_estimates$SD,
  Obj = fit$tmb_list$Obj
)

# 4. Grab the epsilon and Xi bits from #3,
# 5. Use `rmvnorm_prec` function to calculate spatially correlated epsilon and Xi values, supplying values from #3 as “mu” and #1 as “prec” arguments. These are surfaces, though save as a vector indexed by sct (or gct),
# Sample from GMRF using sparse precision
rmvnorm_prec <- function(mu, prec, n.sims, seed) {
  set.seed(seed)
  z <- matrix(rnorm(length(mu) * n.sims), ncol = n.sims)
  L <- Matrix::Cholesky(prec, super = TRUE)
  z <- Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
  z <- Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
  z <- as.matrix(z)
  return(mu + z)
}

tmb_sd <- fit$parameter_estimates$SD
tmb_obj <- fit$tmb_list$Obj
n_sims <- 1
samps <- rmvnorm_prec(mu = tmb_obj$env$last.par.best, prec = tmb_sd, n_sims)

# 6. For epsilon and if we have temporal structure for epsilon, adjust values from #4 given auto-regressive structure (see INLA docs, as I think there is an example simulating spatial-correlation and then adding temporal dependence to those values),
# 7. Save values in ParList,

# If doing multi-year projection
# Do 5-6 as many times sequentially as needed to get the full padding of times t.