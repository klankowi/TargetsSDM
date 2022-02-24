sample_variable <- function(Sdreport, Obj, rmvnorm_mu = Obj$env$last.par.best, variable_name, n_samples = 100, sample_fixed = TRUE, seed = 123456) {
    if (FALSE) {
        # For debugging
        Sdreport <- x$parameter_estimates$SD
        Obj <- x$tmb_list$Obj
        variable_name <- re_match_df$sample_variable_name[1:2]
        n_samples <- 2
        seed <- seed
        sample_fixed <- TRUE
    }

    if (!("jointPrecision" %in% names(Sdreport))) {
        stop("jointPrecision not present in Sdreport; please re-run with `getJointPrecision=TRUE`")
    }

    if (!is.null(variable_name) & any(!(variable_name %in% names(Obj$report())))) {
        stop(
            "At least one of ", variable_name, " not found in report(.); please choose check your requested variable name from available list: ",
            paste(names(Obj$report()), collapse = ", ")
        )
    }

    if (sample_fixed == TRUE) {
        u_zr <- rmvnorm_prec(mu = rmvnorm_mu, prec = Sdreport$jointPrecision, n.sims = n_samples, seed = seed)
    } else {
        u_zr <- rmvnorm_mu %o% rep(1, n_samples)
        MC <- Obj$env$MC(keep = TRUE, n = n_samples, antithetic = FALSE)
        u_zr[Obj$env$random, ] <- attr(MC, "samples")
    }
    if (!is.null(variable_name)) {
        message(
            "# Obtaining samples from predictive distribution for variable ",
            variable_name
        )
    } else {
        message(
            "# Obtaining samples from predictive distribution for all variables"
        )
    }

    for (rI in 1:n_samples) {
        if (rI %% max(1, floor(n_samples / 10)) == 0) {
            message("  Finished sample ", rI, " of ", n_samples)
        }
        if (!is.null(variable_name)) {
            Var <- Obj$report(par = u_zr[, rI])[variable_name]
        } else {
            Var <- Obj$report(par = u_zr[, rI])
        }
        if (rI == 1) {
            Var_zr <- Var
        }
        if (rI >= 2) {
            Var_zr <- sapply(names(Var_zr), FUN = function(x) {
                abind::abind(Var_zr[[x]], Var[[x]], along = length(dim(Var[[x]])) + 1)
            })
        }
    }
    return(Var_zr)
}