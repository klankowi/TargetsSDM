# Helper function for multivariate normal based on a precision matrix rather than variance:covariance matrix.
rmvnorm_prec <- function(mu, prec, n.sims, seed) {
  if (FALSE) {
    mu <- Obj$env$last.par.best
    prec <- Sdreport$jointPrecision
    n.sims <- 1
    seed <- seed
  }
  set.seed(seed)
  z <- matrix(rnorm(length(mu) * n.sims), ncol = n.sims)
  L <- Matrix::Cholesky(prec, super = TRUE)
  z <- Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
  z <- Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
  z <- as.matrix(z)
  return(mu + z)
}

# Padding list time elements with additional 0s for each additional time step. Assumes time dimension is the LAST dimension of the list element
pad_list_notXi_element <- function(x, pad_dim, pad_value = 0, pad_times = n_proj_steps) {
  if (FALSE) {
    x <- ParList$Beta_mean1_t
    pad_dim <- 1
    pad_value <- 0
    pad_times <- n_proj
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
  if (length(dim(x)) == 1 | length(dim(x)) == 0) {
    out_list_element <- c(x, rep(pad_value, pad_times))
  }

  return(out_list_element)
}

pad_list_Xi_element <- function(x, pad_dim, pad_value, pad_start = time_cov_start, pad_times = time_cov_add) {
  if (FALSE) {
    pad_vars <- c("Xiinput1_scp")
    x <- ParList[[pad_vars]]
    pad_dim <- 3
    pad_value <- mean(ParList$gamma1_cp[gamma_mean_start:gamma_mean_end])
    pad_value <- 0
    pad_start <- time_cov_start
    pad_times <- time_cov_add
  }

  # Arrays
  if (length(dim(x)) == 3) {
    # Just appending to the "end"
    if (pad_start >= dim(x)[3]) {
      out_list_element <- unname(abind::abind(array(x[, , 1:(pad_start - 1)], dim = c(dim(x)[1], dim(x)[2], (pad_start - 1))), array(pad_value, replace(dim(x), pad_dim, pad_times)), along = 3))
    } else {
      out_list_element_temp <- unname(abind::abind(array(x[, , 1:(pad_start - 1)], dim = c(dim(x)[1], dim(x)[2], (pad_start - 1))), array(pad_value, replace(dim(x), pad_dim, pad_times)), along = 3))
      out_list_element <- unname(abind::abind(out_list_element_temp, array(x[, , pad_start:dim(x)[3]], dim = c(dim(x)[1], dim(x)[2], length(pad_start:dim(x)[3]))), along = 3))
    }
  }

  # Matrices
  if (length(dim(x)) == 2) {
    if (nrow(x) < ncol(x)) {
      # Just appending to the "end"
      if (pad_start >= dim(x)[2]) {
        out_list_element <- matrix(c(x[1:(pad_start - 1)], rep(pad_value, pad_times)), nrow = 1)
      } else {
        # "Inserting" work in progress...
        out_list_element <- matrix(c(x[, 1:(pad_start - 1)], rep(pad_value, pad_times), x[, pad_start:ncol(x)]), nrow = 1)
      }
    }
  }
  return(out_list_element)
}

# Designing a projection function
#' @title Project density from fitted VAST model.
#'
#' @description Calculates density projections from a VAST `fit_model` object for `new_covariate_data`.
#'
#' @param x = A VAST \code{fit_model} object.
#' @param time_cov = A character string or NULL defining a time covariate fitted as a spatially varying term as in the seasonal model. The current implementation can ONLY handle one of these.
#' @param historical_uncertainty = Character string indicating how to sample uncertainty during the historical period, where "both" samples fixed and random effects, "random" samples only random effects, and "none" assumes no uncertainty during historical period and uses MLE estimates for fixed effects and random effects conditioned on fixed effects (?)
#' @param n_samples = Integer with the number of times to sample uncertainty from the historical period.
#' @param new_covariate_data = A dataframe with new covariate values for projection data.
#' @param new_catchability_data = A dataframe with new catchability values for projection data.
#' @param seed = Random seed number to use for the projections, default = 123456.
#'
#' @return A list with ...
#'
#' @export

project.fit_model <- function(sim_type = 1, x, time_cov, index_shapes, historical_uncertainty = "both", n_samples, n_proj, new_covariate_data = NULL, new_catchability_data = NULL, proj_X_contrasts = NULL, proj_X1_config = NULL, proj_X2_config = NULL, proj_map = NULL, seed = 123456) {

  # For debugging
  if (FALSE) {
    #   library(VAST)
    # sim_type <- 1
    tar_load(vast_fit)
    x <- vast_fit
    time_cov <- "Year_Cov"
    tar_load(index_shapefiles)
    index_shapes <- index_shapefiles
    historical_uncertainty <- "random"
    n_samples <- 1
    n_proj <- 243
    # new_covariate_data <- new_projection_data
    pred_covs_out_final <- readRDS("~/GitHub/TargetsSDM/data/predict/VAST_post_fit_pred_df_seasonal_mean.rds")
    new_covariate_data <- pred_covs_out_final
    new_catchability_data <- NULL

    # sim_type <- 1
    # x <- fit
    # time_cov <- NULL
    # historical_uncertainty <- "both"
    # n_samples <- 1
    # n_proj <- 6
    # new_covariate_data <- NULL
    # new_catchability_data <- NULL
    # proj_X_contrasts <- NULL
    # proj_X1_config <- NULL
    # proj_X2_config <- NULL
    # proj_map <- proj_map
    # seed <- rI
  }

  #####
  ## Should be able to delete this...
  #####
  # Load package
  library(VAST)

  # load data set
  # see `?load_example` for list of stocks with example data
  # that are installed automatically with `FishStatsUtils`.
  example <- load_example(data_set = "EBS_pollock")

  # Make settings (turning off bias.correct to save time for example)
  settings <- make_settings(
    n_x = 100,
    Region = example$Region,
    purpose = "index2",
    bias.correct = FALSE
  )

  # Run model
  fit <- fit_model(
    settings = settings,
    Lat_i = example$sampling_data[, "Lat"],
    Lon_i = example$sampling_data[, "Lon"],
    t_i = example$sampling_data[, "Year"],
    b_i = example$sampling_data[, "Catch_KG"],
    a_i = example$sampling_data[, "AreaSwept_km2"], run_model = FALSE
  )

  #####
  ## Sanity checks
  #####

  # Check columns of new_covariate_data -- from predict function
  if (!is.null(new_covariate_data) & !is.null(x$covariate_data)) {
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

    # Now check factor levels
    projection_levels <- lapply(new_covariate_data[sapply(new_covariate_data, is.factor)], levels)
    fit_levels <- lapply(x$covariate_data[sapply(x$covariate_data, is.factor)], levels)

    fac_level_check <- sapply(names(fit_levels), function(x) {
      all(unlist(fit_levels[x]) %in% unlist(projection_levels[x]))
    })

    if (!any(fac_level_check)) {
      stop("Please check that all factor levels in `x$covariate_data` are present in `new_covariate_data`.")
    }
  }

  # Check if time_cov is only one string
  if (!is.null(time_cov) & length(time_cov) != 1) {
    stop("Projection function is implemented only for one time covariate")
  }

  # Check if time_cov is not null that X1_formula and X2_formula match
  if (!is.null(time_cov) & x$X1_formula != x$X2_formula) {
    stop("Projection function for time covariates is currently implemented only for models with matching X1_formula and X2_formula")
  }

  #####
  ## Generating predictions
  #####

  # Get dimensions of each of the variables, which we will need later on
  var_dims <- unlist(lapply(x$Report, FUN = function(x) length(dim(x))))

  # Get the index for time_cov from model matrix column names now if there is a time_cov, which we will need later on. Accounting for intercept would mean subtracting 1, but we actually want to start padding after the last value, so no adjustment.
  if (!is.null(time_cov)) {
    mod_mat_cols <- colnames(model.matrix(x$X1_formula, data = x$covariate_data, contrasts.arg = x$input_args$extra_args$X_contrasts))
    time_cov_start <- max(which(grepl(time_cov, mod_mat_cols)))
    time_cov_add <- length(projection_levels[[time_cov]]) - length(fit_levels[[time_cov]])
  }

  # Unpack fitted model
  Obj <- x$tmb_list$Obj
  Sdreport <- x$parameter_estimates$SD

  if (historical_uncertainty == "both") {
    # # Sample random effects
    # Obj$retape()
    # Obj$fn(x$parameter_estimates$par)
    # u_z <- Obj$env$last.par.best

    # # Simulate random effects
    # set.seed(seed)
    # MC <- Obj$env$MC(keep = TRUE, n = 1, antithetic = FALSE)
    # u_z[Obj$env$random] <- attr(MC, "samples")[, 1]

    # # Now get the fixed effects
    u_z <- rmvnorm_prec(mu = Obj$env$last.par.best, prec = Sdreport$jointPrecision, n.sims = 1, seed = seed)[, 1]
    # u_z[setdiff(1:length(u_z_fixed), Obj$env$random)] <- u_z_fixed[setdiff(1:length(u_z_fixed), Obj$env$random)]
  } else if (historical_uncertainty == "random") {
    # Retape and call once to get last.par.best to work
    Obj$retape()
    Obj$fn(x$parameter_estimates$par)
    u_z <- Obj$env$last.par.best
    # Simulate random effects
    set.seed(seed)
    MC <- Obj$env$MC(keep = TRUE, n = 1, antithetic = FALSE)
    u_z[Obj$env$random] <- attr(MC, "samples")[, 1]
  } else if (historical_uncertainty == "none") {
    u_z <- Obj$env$last.par.best
  } else {
    stop("Check `historical_uncertainty` argument")
  }

  # Get ParList
  ParList <- Obj$env$parList(par = u_z)

  #####
  ## Padding
  #####
  # Pad beta1_ft/beta2_ft -- pad with value from ParList? Or 0?
  pad_vars <- c("beta1_ft")
  # pad_value <- ParList$beta1_ft[1, ncol(ParList$beta1_ft)]
  pad_value <- 0
  ParList[pad_vars] <- lapply(ParList[pad_vars], FUN = pad_list_notXi_element, pad_dim = 2, pad_value = pad_value, pad_times = n_proj)

  pad_vars <- c("beta2_ft")
  # pad_value <- ParList$beta2_ft[1, ncol(ParList$beta2_ft)]
  pad_value <- 0
  ParList[pad_vars] <- lapply(ParList[pad_vars], FUN = pad_list_notXi_element, pad_dim = 2, pad_value = pad_value, pad_times = n_proj)

  # Pad Beta_mean1_t/mean2_t, both with 0
  pad_vars <- c("Beta_mean1_t", "Beta_mean2_t")
  ParList[pad_vars] <- lapply(ParList[pad_vars], FUN = pad_list_notXi_element, pad_dim = 1, pad_value = 0, pad_times = n_proj)

  # Pad Epsiloninput1/2_sff
  pad_vars <- c("Epsiloninput1_sff", "Epsiloninput2_sff")
  ParList[pad_vars] <- lapply(ParList[pad_vars], FUN = pad_list_notXi_element, pad_dim = 3, pad_value = 0, pad_times = n_proj)

  # Pad lagrange_tc -- what is this?
  pad_vars <- c("lagrange_tc")
  pad_value <- 0
  ParList[pad_vars] <- lapply(ParList[pad_vars], FUN = pad_list_notXi_element, pad_dim = 2, pad_value = pad_value, pad_times = n_proj)

  # Pad "time" Xivariables. We will need to account for added factor levels NOT added years as n_proj that are fitted as spatially varying effects (in this example, season and year) in:
  # gamma1/2_cp (pad with 0) [1, 1:21] 2nd dimension
  # log_sigmaXi1/2_cp (pad with estimated variance instead of 0) [1, 1:21] 2nd dimension
  # Xiinput1/2_scp (pad with 0) [1:304, 1, 1:21] 3rd dimension
  # All of them have "time" as the last dimension, which simplifies things a bit...maybe
  if (!is.null(time_cov)) {
    gamma_mean_end <- length(unique(x$covariate_data$Season)) + length(unique(x$covariate_data$Year_Cov))
    gamma_mean_start <- length(unique(x$covariate_data$Season)) + length(unique(x$covariate_data$Year_Cov))

    # gamma1/2_cp (pad with 0) [1, 1:21] 2nd dimension
    pad_vars <- c("gamma1_cp")
    ParList[pad_vars] <- lapply(ParList[pad_vars], FUN = pad_list_Xi_element, pad_dim = 2, pad_value = mean(ParList$gamma1_cp[gamma_mean_start:gamma_mean_end]), pad_start = time_cov_start, pad_times = time_cov_add)

    pad_vars <- c("gamma2_cp")
    ParList[pad_vars] <- lapply(ParList[pad_vars], FUN = pad_list_Xi_element, pad_dim = 2, pad_value = mean(ParList$gamma1_cp[gamma_mean_start:gamma_mean_end]), pad_start = time_cov_start, pad_times = time_cov_add)

    # log_sigmaXi1/2_cp (pad with parList last value, assume fixed in the future) [1, 1:21] 2nd dimension
    pad_vars <- c("log_sigmaXi1_cp")
    ParList[pad_vars] <- lapply(ParList[pad_vars], FUN = pad_list_Xi_element, pad_dim = 2, pad_value = ParList$log_sigmaXi1_cp[gamma_mean_end], pad_start = time_cov_start, pad_times = time_cov_add)

    pad_vars <- c("log_sigmaXi2_cp")
    ParList[pad_vars] <- lapply(ParList[pad_vars], FUN = pad_list_Xi_element, pad_dim = 2, pad_value = ParList$log_sigmaXi2_cp[gamma_mean_end], pad_start = time_cov_start, pad_times = time_cov_add)

    # Xiinput1/2 scp (pad with 0), time is the third dimension
    pad_vars <- c("Xiinput1_scp", "Xiinput2_scp")
    ParList[pad_vars] <- lapply(ParList[pad_vars], FUN = pad_list_Xi_element, pad_dim = 3, pad_value = 0, pad_start = time_cov_start, pad_times = time_cov_add)
  }

  # Done padding, could add some sort of sanity check here on dimensions.

  #####
  ## Rebuilding model
  #####
  # Processing inputs
  if (!is.null(new_covariate_data)) {
    t_i <- c(x$data_frame$t_i, new_covariate_data$Year)
    b_i <- c(x$data_frame$b_i, sample(c(0, mean(x$data_frame$b_i)), size = nrow(new_covariate_data), replace = TRUE))
    v_i <- c(x$data_frame$v_i, rep(0, nrow(new_covariate_data)))
    Lon_i <- c(x$data_frame$Lon_i, new_covariate_data$Lon)
    Lat_i <- c(x$data_frame$Lat_i, new_covariate_data$Lat)
    a_i <- c(x$data_frame$a_i, rep(mean(x$data_frame$a_i), nrow(new_covariate_data)))
    PredTF_i <- c(x$data_list$PredTF_i, rep(1, nrow(new_covariate_data)))
    c_iz <- rbind(x$data_list$c_iz, x$data_list$c_iz[sample(1:nrow(x$data_list$c_iz), size = nrow(new_covariate_data), replace = TRUE), , drop = FALSE])
    all_covariate_data <- rbind(x$covariate_data, new_covariate_data)
    all_catchability_data <- NULL
  } else {
    t_i <- c(x$data_frame$t_i, max(x$data_frame$t_i) + rep(1:n_proj, each = 2))
    b_i <- c(x$data_frame$b_i, rep(c(0, mean(x$data_frame$b_i)), n_proj))
    v_i <- c(x$data_frame$v_i, rep(0, 2 * n_proj))
    Lon_i <- c(x$data_frame$Lon_i, rep(mean(x$data_frame$Lon_i), 2 * n_proj))
    Lat_i <- c(x$data_frame$Lat_i, rep(mean(x$data_frame$Lat_i), 2 * n_proj))
    a_i <- c(x$data_frame$a_i, rep(mean(x$data_frame$a_i), 2 * n_proj))
    PredTF_i <- c(x$data_list$PredTF_i, rep(1, 2 * n_proj))
    c_iz <- rbind(x$data_list$c_iz, x$data_list$c_iz[rep(1:n_proj, each = 2), , drop = FALSE])
    all_covariate_data <- NULL
    all_catchability_data <- NULL
  }

  # Rebuild
  proj_mod <- fit_model_aja(
    settings = x$settings,
    input_grid = x$input_args$data_args_input$input_grid,
    index_shapes = index_shapes,
    Method = x$settings$Method,
    Lat_i = Lat_i,
    Lon_i = Lon_i,
    t_i = t_i,
    b_i = b_i,
    a_i = a_i,
    v_i = v_i,
    c_iz = c_iz,
    PredTF_i = PredTF_i,
    covariate_data = all_covariate_data,
    X1_formula = x$X1_formula,
    X2_formula = x$X2_formula,
    X1config_cp = proj_X1_config,
    X2config_cp = proj_X2_config,
    X_contrasts = proj_X_contrasts,
    catchability_data = all_catchability_data,
    Q1config_k = x$Q1config_k,
    Q2config_k = x$Q2config_k,
    Q1_formula = x$Q1_formula,
    Q2_formula = x$Q2_formula,
    run_model = FALSE,
    Parameters = ParList,
    Map = proj_map
    # Random = NULL
    # "getJointPrecision" = TRUE,
    # "report_additional_variables" = TRUE
  )

  #####
  ## Simulating with the new model
  #####
  # Which type would we use? I would have thought type = 3 as this will simulate the random effects values for the time steps we padded, but that doesn't work? Complains about getJointprecision, which is set to true in original fit?
  # Simulate Epsiloninput / Betainput for projection years
  proj_mod$tmb_list$Obj$env$data$Options_list$simulate_t <- c(rep(0, x$data_list$n_t), rep(1, n_proj))

  # Simulate type=1 so Omegas and other random effects are held fixed
  Sim <- simulate_data.aja(
    fit = proj_mod,
    type = 1,
    joint_prec = x$parameter_estimates$SD$jointPrecision,
    par_use = ParList,
    random_seed = NULL
  )

  proj_mod$Report <- Sim
  Sim <- amend_output(proj_mod)

  return(Sim)

  # End function
}

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
summary.sim_results <- function(vast_fit, sim_obj, what, nice_times = NULL, probs = c(0.1, 0.5, 0.9), mean_instead = FALSE, nice_category_names = nice_category_names, climate_scenario = climate_scenario, out_dir) {
  if (FALSE) {
    tar_load(vast_fit)
    tar_load(vast_projections)
    sim_obj <- vast_projections
    what <- "D_gct"
    probs <- c(0.1, 0.5, 0.9)
    nice_times <- seq(as.Date("1980-06-15"), as.Date("2092-06-15"), by = "year")
    nice_times <- seq(from = 1, to = sim_obj[[1]]$n_t)
    mean_instead <- FALSE
  }

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
        res_out<- bind_rows(res_out, temp_df)
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
        res_out<- bind_rows(res_out, temp_df)
      }
    }
  }

  # Wide to long for "Index_ctl"
  if (what == "Index_ctl") {
    res_out <- res_out %>%
      pivot_longer(., !c(Sim_Scenario, Time), names_to = "Region", values_to = "Index") %>%
      group_by(., Time, Region) %>%
      summarise(
        Prob_0.5 = quantile(Index, probs = 0.50),
        Prob_0.1 = quantile(Index, probs = 0.1),
        Prob_0.9 = quantile(Index, probs = 0.9)
      )
  }
  
  if (what == "D_gct") {
    res_out <- res_out %>%
      group_by(., Lat, Lon, Time) %>%
      summarise(
        Prob_0.5 = quantile(D_gct, probs = 0.50),
        Prob_0.1 = quantile(D_gct, probs = 0.1),
        Prob_0.9 = quantile(D_gct, probs = 0.9)
      )
  }
  
  
  # Save and return it
  saveRDS(res_out, file = paste0(out_dir, "/", nice_category_names, "_", climate_scenario, "_", what, ".rds"))
  return(res_out)
}

simulate_data.aja <- function(fit, type = 1, par_use, joint_prec = NULL, random_seed = NULL) {
  rmvnorm_prec <- function(mu, prec, n.sims, random_seed) {
    set.seed(random_seed)
    z <- matrix(rnorm(length(mu) * n.sims), ncol = n.sims)
    L <- Matrix::Cholesky(prec, super = TRUE)
    z <- Matrix::solve(L, z, system = "Lt")
    z <- Matrix::solve(L, z, system = "Pt")
    z <- as.matrix(z)
    return(mu + z)
  }
  dlls <- getLoadedDLLs()
  isTMBdll <- function(dll) !is(try(getNativeSymbolInfo("MakeADFunObject", dll), TRUE), "try-error")
  TMBdll <- sapply(dlls, isTMBdll)
  if (sum(TMBdll) == 0) {
    stop("VAST is not linked as a DLL, so `simulate_data` will not work.\n    Please re-run model (potentially from informative starting values to save time) to use `simulate_data`")
  } else if (sum(TMBdll) >= 2) {
    warning("VAST is linked to multiple DLLs. Please consider using dyn.unload() to unload\n    earlier VAST runs to avoid potentially ambiguous behavior when running `simulate_data`")
  }
  Obj <- fit$tmb_list$Obj
  simulate_random_effects_orig <- Obj$env$data$Options_list$Options["simulate_random_effects"]
  revert_settings <- function(simulate_random_effects) {
    Obj$env$data$Options_list$Options["simulate_random_effects"] <- simulate_random_effects
  }
  on.exit(revert_settings(simulate_random_effects_orig))
  if (type == 1) {
    Obj$env$data$Options_list$Options["simulate_random_effects"] <- FALSE
    set.seed(random_seed)
    Return <- Obj$simulate(complete = TRUE)
  }
  if (type == 2) {
    Obj$env$data$Options_list$Options["simulate_random_effects"] <- TRUE
    set.seed(random_seed)
    Return <- Obj$simulate(complete = TRUE)
  }
  if (type == 3) {
    # Informative error messages
    if (!("jointPrecision" %in% names(fit$parameter_estimates$SD))) {
      stop("jointPrecision not present in fit$parameter_estimates$SD; please re-run with `getJointPrecision=TRUE`")
    }

    # Sample from joint distribution
    newpar <- rmvnorm_prec(mu = Obj$env$last.par.best, prec = fit$parameter_estimates$SD$jointPrecision, n.sims = 1, random_seed = random_seed)[, 1]

    # Simulate
    Obj$env$data$Options_list$Options["simulate_random_effects"] <- FALSE
    Return <- Obj$simulate(par = newpar, complete = TRUE)
  }
  if (type == 4) {
    set.seed(random_seed)
    warning("Type-4 residuals are still under development, please use with care and note that they may change at any point.")
    newpar <- Obj$env$last.par.best
    MC <- Obj$env$MC(keep = TRUE, n = 1, antithetic = FALSE)
    newpar[Obj$env$random] <- attr(MC, "samples")
    Obj$env$data$Options_list$Options["simulate_random_effects"] <- FALSE
    Return <- Obj$simulate(par = newpar, complete = TRUE)
  }
  return(Return)
}

