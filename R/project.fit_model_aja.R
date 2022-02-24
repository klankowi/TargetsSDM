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

# Designing a projection function
#' @title Project density from fitted VAST model.
#'
#' @description Calculates density projections from a VAST `fit_model` object for `new_covariate_data`.
#'
#' @param x = A VAST \code{fit_model} object.
#' @param time_cov = A character string or NULL defining a time covariate fitted as a spatially varying term as in the seasonal model. The current implementation can ONLY handle one of these.
#' @param category = Integer for which category to project.
#' @param what = Which output from \code{fit_model$Report} should be extracted; by default projected density.
#' @param historical_uncertainty = Character string indicating how to sample uncertainty during the historical period, where "both" samples fixed and random effects, "random" samples only random effects, and "none" assumes no uncertainty during historical period and uses MLE estimates for fixed effects and random effects conditioned on fixed effects (?)
#' @param n_samples = Integer with the number of times to sample uncertainty from the historical period.
#' @param new_covariate_data = A dataframe with new covariate values for projection data.
#' @param new_catchability_data = A dataframe with new catchability values for projection data.
#' @param seed = Random seed number to use for the projections, default = 123456.
#'
#' @return A dataframe with location and time information for each sample in `new_projection_data` and the projected \code{what} output variable.
#'
#' @export
project.fit_model <- function(x, time_cov, category = 1, what = "D_gct", historical_uncertainty = "both", n_samples, new_covariate_data, new_catchability_data = NULL, seed = 123456) {

  # For debugging
  if (FALSE) {
      library(VAST)
      x = readRDS("~/Desktop/YTseasonalmod.rds")
      time_cov = c("Year_Cov")
      category = 1
      what = "D_gct"
      historical_uncertainty = "both"
      n_samples = 2
      new_covariate_data <- readRDS("~/Desktop/MockProjectionData.rds")
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
      
      # I'm not entirely sure why the crash happens. It does seems like we are going to need to adjust some other things -- X1/X2 config and Map -- as the dimensions of these will not match from original model fit to what we want with the re-built model.
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
      
      proj_map$gamma1_cp <- factor(c(proj_map$gamma1_cp, seq(from = max(as.numeric(as.character(levels(proj_map$gamma1_cp))))+1, to = max(as.numeric(as.character(levels(proj_map$gamma1_cp)))) + length(levels(all_covariate_data$Year_Cov)[!levels(all_covariate_data$Year_Cov) %in% levels(x$covariate_data$Year_Cov)]))))
      proj_map$gamma2_cp <- factor(c(proj_map$gamma2_cp, seq(from = max(as.numeric(as.character(levels(proj_map$gamma2_cp))))+1, to = max(as.numeric(as.character(levels(proj_map$gamma2_cp)))) + length(levels(all_covariate_data$Year_Cov)[!levels(all_covariate_data$Year_Cov) %in% levels(x$covariate_data$Year_Cov)]))))
      
      proj_map$log_sigmaXi1_cp <- factor(c(rep(1, length(unique(all_covariate_data$Season))), rep(4, nlevels(all_covariate_data$Year_Cov))))
      proj_map$log_sigmaXi2_cp <- factor(c(rep(1, length(unique(all_covariate_data$Season))), rep(4, nlevels(all_covariate_data$Year_Cov))))
      
      proj_map$Xiinput1_scp <- factor(seq(from = 1, to = dim(ParList$Xiinput1_scp)[1] * dim(ParList$Xiinput1_scp)[3]))
      proj_map$Xiinput2_scp <- factor(seq(from = 1, to = dim(ParList$Xiinput2_scp)[1] * dim(ParList$Xiinput2_scp)[3]))
      
      proj_map$beta1_ft <- factor(rep(1, length(unique(all_covariate_data$Year))))
      proj_map$beta2_ft <- factor(rep(1, length(unique(all_covariate_data$Year))))
      
      proj_map$Beta_mean1_t <- factor(rep(NA, length(unique(all_covariate_data$Year))))
      proj_map$Beta_mean2_t <- factor(rep(NA, length(unique(all_covariate_data$Year))))

      # New call to fit_model...
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
          covariate_data = all_covariate_data,
          X1_formula = x$X1_formula,
          X2_formula = x$X2_formula,
          X1config_cp = X1config_cp_proj,
          X2config_cp = X2config_cp_proj,
          X_contrasts = list(Season = contrasts(all_covariate_data$Season, contrasts = FALSE), Year_Cov = contrasts(all_covariate_data$Year_Cov, contrasts = FALSE)),
          catchability_data = all_catchability_data,
          Q1config_k = x$Q1config_k,
          Q2config_k = x$Q2config_k,
          Q1_formula = x$Q1_formula,
          Q2_formula = x$Q2_formula,
          run_model = FALSE,
          Parameters = ParList,
          Map = proj_map, 
          Random = NULL,
          "getJointPrecision" = TRUE
      )
      
      #####
      ## Simulating with the new model
      #####
      # Which type would we use? I would have thought type = 3 as this will simulate the random effects values for the time steps we padded, but that doesn't work? Complains about getJointprecision, which is set to true in original fit?
      sim_temp <- simulate_data(
        fit = proj_mod,
        type = 1,
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
  }
  # Not sure if other options make sense for n_samples?

  # End function
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
