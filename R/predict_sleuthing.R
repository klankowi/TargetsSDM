# Predict function crash sleuthing #

# On restart...
library(devtools)
library(VAST)
library(tidyverse)
library(splines)  # Used to include basis-splines
library(effects)  # Used to visualize covariate effects
library(sf)
library(raster)
library(patchwork)
library(akima)
library(lubridate)
library(profvis)
library(googledrive)
library(PresenceAbsence)
library(MLmetrics)
library(gmRi)
library(transformr)
library(gganimate)

# Source a bunch of "vast" helper function. These are all located here: https://github.com/aallyn/TargetsSDM. Working on getting this into a nicer package style. I think this could also change if you fork/clone the repo, then you can source from the local copies...
source(here::here("", "R/nmfs_functions.R"))
source(here::here("", "R/dfo_functions.R"))
source(here::here("", "R/combo_functions.R"))
source(here::here("", "R/covariate_functions.R"))
source(here::here("", "R/vast_functions.R"))
source(here::here("", "R/SDM_PredValidation_Functions.R"))

# Run?
eval_use<- FALSE

# Read in tidy model data and create a fake prediction dataframe
tidy_mod_data<- readRDS(here::here("", "data/combined/tidy_mod_data.rds"))

# Filter tidy model data to keep only one species of interest -- Halibut in this case
tidy_mod_data_run<- tidy_mod_data %>%
  dplyr::filter(., DFO_SPEC == 30)

# Prep the sample dataframe, we also include a "Pred_TF" column -- when Pred_TF == 1, the observation is only included in the predictions and NOT in the model fitting. This is a helpful way to do some quick model validation.
vast_sample_data<- data.frame("Year" = tidy_mod_data_run$EST_YEAR, "Lat" = tidy_mod_data_run$DECDEG_BEGLAT, "Lon" = tidy_mod_data_run$DECDEG_BEGLON, "Biomass" = tidy_mod_data_run$BIOMASS, "Swept" = ifelse(tidy_mod_data_run$SURVEY == "NMFS", 0.0384, 0.0404), "Pred_TF" = rep(0, nrow(tidy_mod_data_run)))

# We will need a covariate data frame
vast_covariate_data<- data.frame("Year" = tidy_mod_data_run$EST_YEAR, "Depth" = tidy_mod_data_run$Depth, "SST_seasonal" = tidy_mod_data_run$SST_seasonal, "BT_seasonal" = tidy_mod_data_run$BT_seasonal, "Lat" = tidy_mod_data_run$DECDEG_BEGLAT, "Lon" = tidy_mod_data_run$DECDEG_BEGLON)

dummy_dat_cov<- vast_covariate_data[nrow(vast_covariate_data),]
dummy_dat_cov$Year<- 2020
vast_covariate_data_aug<- rbind(vast_covariate_data, dummy_dat_cov)

dummy_dat<- vast_sample_data[nrow(vast_sample_data),]
dummy_dat$Year<- 2020
dummy_dat$Pred_TF<- 1

vast_sample_data_aug<- rbind(vast_sample_data, dummy_dat)

# Read in fitted model
vast_fitted_hab_covs<- readRDS(here::here("", "results/mod_fits/1011_fitted_vast.rds"))

# Grab a sample of data from the last three years of observations to use as "fake" temperature forecasts
fitted_dat<- vast_fitted_hab_covs$covariate_data
pred_df<- fitted_dat %>%
  filter(., Year >= 2015)
pred_df<- pred_df[sample(1:nrow(pred_df), size = 250, replace = FALSE), ]

# Manually fudge the year to pretend like it comes from "the future"
pred_df$Year<- 2020

summary(pred_df)

# Now, break into...predict_vast function to get our "prediction" dataframe formatted correctly
nmfs_species_code = 101
predict_variable = "D_i"
predict_category = 0
predict_vessel = 0
predict_covariates_df_all<- pred_df
time_col = "Year"
cov_names = c("Depth", "SST_seasonal", "BT_seasonal")

#### Not the biggest fan of this, but for now, building in a work around to resolve some of the memory issues that we were running into by supplying a 0.25 degree grid and trying to predict/project for each season-year from 1980-2100. To overcome this issue, going to try to just make the projections to knots and do the smoothing later.
# First, need to get the knot locations
knot_locs<- data.frame(vast_fitted_sdm$spatial_list$latlon_g) %>%
  st_as_sf(., coords = c("Lon", "Lat"), remove = FALSE) %>%
  mutate(., "Pt_Id" = 1:nrow(.))

# Nearest knot to each point?
pred_sf<- predict_covariates_df_all %>%
  st_as_sf(., coords = c("Lon", "Lat"), remove = FALSE)

pred_sf<- pred_sf %>%
  mutate(., "Nearest_Knot" = st_nearest_feature(., knot_locs))

# Average the points...
pred_df_knots<- pred_sf %>%
  st_drop_geometry()
group_by_vec<- c({{time_col}}, "Nearest_Knot")
pred_df_knots<- pred_df_knots %>%
  group_by_at(.vars = group_by_vec) %>%
  summarize_at(all_of(cov_names), mean, na.rm = TRUE) %>%
  left_join(., st_drop_geometry(knot_locs), by = c("Nearest_Knot" = "Pt_Id")) %>%
  ungroup()

# Collecting necessary bits from the prediction covariates -- lat, lon, time
pred_lats<- pred_df_knots$Lat
pred_lons<- pred_df_knots$Lon
pred_times<- as.numeric(unlist(pred_df_knots[{{time_col}}]))

# Catch stuff...
pred_sampled_areas<- rep(1, length(pred_lats))
pred_category<- rep(predict_category, length(pred_lats))
pred_vessel<- rep(predict_vessel, length(pred_lats))

# Renaming predict_covariates_df_all to match vast_fit_covariate_data
pred_cov_dat_name_order<- which(names(pred_df_knots) %in% names(vast_fitted_sdm$covariate_data))
pred_cov_dat_use<- pred_df_knots[,pred_cov_dat_name_order]

# Catchability data?
if(!is.null(vast_fitted_sdm$catchability_data)){
  pred_catch_dat_use<- pred_cov_dat_use %>%
    dplyr::select(., c(Year, Year_Cov, Season, Lat, Lon, Survey)
    )
  pred_catch_dat_use$Survey<- rep("NMFS", nrow(pred_catch_dat_use))
  pred_catch_dat_use$Survey<- factor(pred_catch_dat_use$Survey, levels = c("NMFS", "DFO", "DUMMY"))
} else {
  pred_catch_dat_use<- NULL
}

# Now break into predict.fit_model_aja
x = vast_fitted_hab_covs
what = predict_variable
Lat_i = pred_lats
Lon_i = pred_lons
t_i = pred_times
a_i = pred_sampled_areas
c_iz = pred_category
v_i = rep(0,length(t_i)) 
new_covariate_data = pred_cov_dat_use
new_catchability_data = pred_catch_dat_use
do_checks = FALSE
working_dir = paste0(getwd(), "/")

message("`predict.fit_model(.)` is in beta-testing, and please explore results carefully prior to using")

# Check issues
if( !(what%in%names(x$Report)) || (length(x$Report[[what]])!=x$data_list$n_i) ){
  stop("`what` can only take a few options")
}
if( !is.null(new_covariate_data) ){
  # Confirm all columns are available
  if( !all(colnames(x$covariate_data) %in% colnames(new_covariate_data)) ){
    stop("Please ensure that all columns of `x$covariate_data` are present in `new_covariate_data`")
  }
  # Eliminate unnecessary columns
  new_covariate_data = new_covariate_data[,match(colnames(x$covariate_data),colnames(new_covariate_data))]
  # Eliminate old-covariates that are also present in new_covariate_data
  NN = RANN::nn2( query=x$covariate_data[,c('Lat','Lon','Year')], data=new_covariate_data[,c('Lat','Lon','Year')], k=1 )
  if( any(NN$nn.dist==0) ){
    x$covariate_data = x$covariate_data[-which(NN$nn.dist==0),,drop=FALSE]
  }
}
if( !is.null(new_catchability_data) ){
  # Confirm all columns are available
  if( !all(colnames(x$catchability_data) %in% colnames(new_catchability_data)) ){
    stop("Please ensure that all columns of `x$catchability_data` are present in `new_covariate_data`")
  }
  # Eliminate unnecessary columns
  new_catchability_data = new_catchability_data[,match(colnames(x$catchability_data),colnames(new_catchability_data))]
  # Eliminate old-covariates that are also present in new_covariate_data
  NN = RANN::nn2( query=x$catchability_data[,c('Lat','Lon','Year')], data=new_catchability_data[,c('Lat','Lon','Year')], k=1 )
  if( any(NN$nn.dist==0) ){
    x$catchability_data = x$catchability_data[-which(NN$nn.dist==0),,drop=FALSE]
  }
}

# Process covariates
covariate_data = rbind( x$covariate_data, new_covariate_data )
catchability_data = rbind( x$catchability_data, new_catchability_data )

# Process inputs
PredTF_i = c( x$data_list$PredTF_i, rep(1,length(t_i)) )
b_i = c( x$data_frame[,"b_i"], sample(c(0, 1), size = length(t_i), replace = TRUE))
c_iz = rbind( matrix(x$data_frame[,grep("c_iz",names(x$data_frame))]), matrix(c_iz) )
Lat_i = c( x$data_frame[,"Lat_i"], Lat_i )
Lon_i = c( x$data_frame[,"Lon_i"], Lon_i )
a_i = c( x$data_frame[,"a_i"], a_i )
v_i = c( x$data_frame[,"v_i"], v_i )
t_i = c( x$data_frame[,"t_i"], t_i )
#assign("b_i", b_i, envir=.GlobalEnv)

# Build information regarding spatial location and correlation
message("\n### Re-making spatial information")
spatial_args_new = list("anisotropic_mesh"=x$spatial_list$MeshList$anisotropic_mesh, "Kmeans"=x$spatial_list$Kmeans, "Lon_i"=Lon_i, "Lat_i"=Lat_i )
spatial_args_input = combine_lists( input=spatial_args_new, default=x$input_args$spatial_args_input )
spatial_list = do.call( what=make_spatial_info, args=spatial_args_input )

# Check spatial_list
if( !all.equal(spatial_list$MeshList,x$spatial_list$MeshList) ){
  stop("`MeshList` generated during `predict.fit_model` doesn't match that of original fit; please email package author to report issue")
}

# Build data
# Do *not* restrict inputs to formalArgs(make_data) because other potential inputs are still parsed by make_data for backwards compatibility
message("\n### Re-making data object")
data_args_new = list( "c_iz"=c_iz, "b_i"=b_i, "a_i"=a_i, "v_i"=v_i, "PredTF_i"=PredTF_i,
                      "t_i"=t_i, "spatial_list"=spatial_list,
                      "covariate_data"=covariate_data, "catchability_data"=catchability_data )
data_args_input = combine_lists( input=data_args_new, default=x$input_args$data_args_input )  # Do *not* use args_to_use
data_list = do.call( what=make_data, args=data_args_input )
data_list$n_g = 0

# Build object
message("\n### Re-making TMB object")
model_args_default = list("TmbData"=data_list, "RunDir"=working_dir, "Version"=x$settings$Version, "RhoConfig"=x$settings$RhoConfig, "loc_x"=spatial_list$loc_x, "Method"=spatial_list$Method, "Map" = x$tmb_list$Map)
model_args_input = combine_lists( input=list("Parameters"=x$ParHat),
                                  default=model_args_default, args_to_use=formalArgs(make_model) )
tmb_list = do.call( what=make_model, args=model_args_input )

# Extract output
TmbData<- tmb_list$data_list
Sdreport<- tmb_list$parameter_estimates$SD


Report = tmb_list$Obj$report()
Y_i = Report[[what]][(1+nrow(x$data_frame)):length(Report$D_i)]

# sanity check
#if( all.equal(covariate_data,x$covariate_data) & Report$jnll!=x$Report$jnll){
if( do_checks==TRUE && (Report$jnll!=x$Report$jnll) ){
  message("Problem detected in `predict.fit_model`; returning outputs for diagnostic purposes")
  Return = list("Report"=Report, "data_list"=data_list)
  return(Return)
}

# return prediction
return(Y_i)






#####
## What about sdmTMB's predict since we can't get population level effects...
#####
newdata = pred_cov_dat_use
object = vast_fitted_hab_covs

library(ggplot2)
pcod_spde <- sdmTMB::make_mesh(d, c("X", "Y"), cutoff = 15)
unique(d$year)
m <- sdmTMB(data = d, formula = density ~ 1, ar1_fields = TRUE,  extra_time = 2019L, include_spatial = FALSE, time = "year", spde = pcod_spde, family = tweedie(link = "log"))

# Add a year to our grid:
grid2019 <- qcs_grid[qcs_grid$year == max(qcs_grid$year), ]
grid2019$year <- 2019L # `L` because `year` is an integer in the data
qcsgrid_forecast <- rbind(qcs_grid, grid2019)

predictions <- predict(m, newdata = qcsgrid_forecast, re_form = NA)



re_form<- NA # Used for population level predictions, which Index_gctl would be
xy_cols<- attributes(object$spatial_list$latlon_g)$dimnames[[2]]

test <- suppressWarnings(tryCatch(object$tmb_list$Obj$report(), error = function(e) NA))
if (all(is.na(test))) object <- update_model(object)

# from glmmTMB:
pop_pred<- (!is.null(re_form) && ((re_form == ~0) || identical(re_form, NA)))

tmb_data <- object$data_list
tmb_data$do_predict <- 1L

if (!is.null(newdata)) {
  if (any(!xy_cols %in% names(newdata)) && isFALSE(pop_pred))
    stop("`xy_cols` (the column names for the x and y coordinates) ",
         "are not in `newdata`. Did you miss specifying the argument ",
         "`xy_cols` to match your data? The newer `make_mesh()` ",
         "(vs. `make_spde()`) takes care of this for you.", call. = FALSE
    )
  
  original_time <- sort(unique(object$data_frame$t_i))
  new_data_time <- sort(unique(newdata$Year))
  
  if (!all(new_data_time %in% original_time))
    stop("Some new time elements were found in `newdata`. ",
         "For now, make sure only time elements from the original dataset ",
         "are present. If you would like to predict on new time elements, see ",
         "the `extra_time` argument in `?sdmTMB:::predict.sdmTMB`.",
         call. = FALSE
    )

  # If making population predictions (with standard errors), we don't need to worry about space, so fill in dummy values if the user hasn't made any:
  fake_spatial_added <- FALSE
  if (pop_pred) {
    for (i in c(1, 2)) {
      if (!xy_cols[[i]] %in% names(newdata)) {
        newdata[[xy_cols[[i]]]] <- mean(object$data[[xy_cols[[i]]]], na.rm = TRUE)
        fake_spatial_added <- TRUE
      }
    }
  }
  
  if (sum(is.na(new_data_time)) > 0)
    stop("There is at least one NA value in the time column. ",
         "Please remove it.", call. = FALSE)
  
  newdata$sdm_orig_id <- seq(1, nrow(newdata))
  fake_newdata <- unique(newdata[,xy_cols])
  fake_newdata[["sdm_spatial_id"]] <- seq(1, nrow(fake_newdata)) - 1L
  
  newdata <- base::merge(newdata, fake_newdata, by = xy_cols, all.x = TRUE, all.y = FALSE)
  newdata <- newdata[order(newdata$sdm_orig_id),, drop=FALSE]
  
  proj_mesh <- INLA::inla.spde.make.A(object$spatial_list$MeshList$anisotropic_mesh, loc = as.matrix(fake_newdata[,xy_cols, drop = FALSE]))
  
  # this formula has breakpt() etc. in it:
  #thresh <- check_and_parse_thresh_params(object$formula, newdata)
  formula <- object$X1_formula # this one does not
  
  nd <- newdata
  response <- "Biomass"
  sdmTMB_fake_response <- FALSE
  if (!response %in% names(nd)) {
    nd[[response]] <- 0 # fake for model.matrix
    sdmTMB_fake_response <- TRUE
  }
  
  # deal with prediction IID random intercepts:
  RE_names <- barnames(object$split_formula$reTrmFormulas)
  ## not checking so that not all factors need to be in prediction:
  # fct_check <- vapply(RE_names, function(x) check_valid_factor_levels(data[[x]], .name = x), TRUE)
  proj_RE_indexes <- vapply(RE_names, function(x) as.integer(nd[[x]]) - 1L, rep(1L, nrow(nd)))
  nobs_RE <- unname(apply(proj_RE_indexes, 2L, max)) + 1L
  if (length(nobs_RE) == 0L) nobs_RE <- 0L
  ln_tau_G_index <- unlist(lapply(seq_along(nobs_RE), function(i) rep(i, each = nobs_RE[i]))) - 1L
  
  #if (!"mgcv" %in% names(object)) object[["mgcv"]] <- FALSE
  proj_X_ij <- matrix(999)
  if (!object$mgcv) {
    proj_X_ij <- tryCatch({model.matrix(object$formula, data = nd)}, error = function(e) NA)
  }
  if (object$mgcv || identical(proj_X_ij, NA)) {
    proj_X_ij <- mgcv::predict.gam(object$mgcv_mod, type = "lpmatrix", newdata = nd)
  }
  if (!is.null(object$time_varying))
    proj_X_rw_ik <- model.matrix(object$time_varying, data = nd)
  else
    proj_X_rw_ik <- matrix(0, ncol = 1, nrow = 1) # dummy
  
  if (length(area) != nrow(proj_X_ij) && length(area) != 1L) {
    stop("`area` should be of the same length as `nrow(newdata)` or of length 1.", call. = FALSE)
  }
  
  tmb_data$proj_X_threshold <- thresh$X_threshold
  tmb_data$area_i <- if (length(area) == 1L) rep(area, nrow(proj_X_ij)) else area
  tmb_data$proj_mesh <- proj_mesh
  tmb_data$proj_X_ij <- proj_X_ij
  tmb_data$ln_tau_G_index <- ln_tau_G_index
  tmb_data$proj_X_rw_ik <- proj_X_rw_ik
  tmb_data$proj_RE_indexes <- proj_RE_indexes
  tmb_data$proj_year <- make_year_i(nd[[object$time]])
  tmb_data$proj_lon <- newdata[[xy_cols[[1]]]]
  tmb_data$proj_lat <- newdata[[xy_cols[[2]]]]
  tmb_data$calc_se <- as.integer(se_fit)
  tmb_data$pop_pred <- as.integer(pop_pred)
  tmb_data$calc_time_totals <- as.integer(!se_fit)
  tmb_data$proj_spatial_index <- newdata$sdm_spatial_id
  tmb_data$proj_t_i <- as.numeric(newdata[[object$time]])
  tmb_data$proj_t_i <- tmb_data$proj_t_i - mean(unique(tmb_data$proj_t_i)) # center on mean
  
  epsilon_covariate <- rep(0, length(unique(newdata[[object$time]])))
  if (tmb_data$est_epsilon_model) {
    # covariate vector dimensioned by number of time steps
    time_steps <- unique(newdata[[object$time]])
    for (i in seq_along(time_steps)) {
      epsilon_covariate[i] <- newdata[newdata[[object$time]] == time_steps[i],
                                      object$epsilon_predictor, drop = TRUE][[1]]
    }
  }
  tmb_data$epsilon_predictor <- epsilon_covariate
  
  new_tmb_obj <- TMB::MakeADFun(
    data = tmb_data,
    parameters = object$tmb_obj$env$parList(),
    map = object$tmb_map,
    random = object$tmb_random,
    DLL = "sdmTMB",
    silent = TRUE
  )
  
  old_par <- object$model$par
  # need to initialize the new TMB object once:
  new_tmb_obj$fn(old_par)
  lp <- new_tmb_obj$env$last.par.best
  
  r <- new_tmb_obj$report(lp)
  
  if (isFALSE(pop_pred)) {
    nd$est <- r$proj_eta
    nd$est_non_rf <- r$proj_fe
    nd$est_rf <- r$proj_rf
    nd$omega_s <- r$proj_re_sp_st
    nd$zeta_s <- r$proj_re_sp_slopes
    nd$epsilon_st <- r$proj_re_st_vector
  }
  
  nd$sdm_spatial_id <- NULL
  nd$sdm_orig_id <- NULL
  
  obj <- new_tmb_obj
  
  if (se_fit) {
    sr <- TMB::sdreport(new_tmb_obj, bias.correct = FALSE)
    ssr <- summary(sr, "report")
    if (pop_pred) {
      proj_eta <- ssr[row.names(ssr) == "proj_fe", , drop = FALSE]
    } else {
      proj_eta <- ssr[row.names(ssr) == "proj_eta", , drop = FALSE]
    }
    row.names(proj_eta) <- NULL
    d <- as.data.frame(proj_eta)
    names(d) <- c("est", "se")
    nd$est <- d$est
    nd$est_se <- d$se
  }
  
  if ("sdmTMB_fake_year" %in% names(nd)) {
    nd <- nd[!nd$sdmTMB_fake_year,,drop=FALSE]
    nd$sdmTMB_fake_year <- NULL
  }
  if (fake_spatial_added) {
    for (i in 1:2) nd[[xy_cols[[i]]]] <- NULL
  }
  if (sdmTMB_fake_response) {
    nd[[response]] <- NULL
  }