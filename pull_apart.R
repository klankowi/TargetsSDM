
rm(list=ls())

library(VAST)
library(here)
library(tidyverse)

source(here('R/VAST_functions.R'))

load("C:/Users/klankowicz/Documents/GitHub/Atlantic-Cod-Habitat-VAST/VAST_runs/refine_effort/refine_effort.RData")

ncat <- fit$data_list$n_c
cat.to.use <- 1

smalls <- fit
smalls$data_frame <- subset(smalls$data_frame, c_iz == '0')
smalls$covariate_data
smalls$extrapolation_list
smalls$spatial_list
smalls$data_list$n_i <- nrow(survs[survs$Age == '0',])
smalls$data_list$n_c <- 1
smalls$data_list$n_e <- 1
smalls$data_list$Options_list$Expansion_cz <- smalls$data_list$Options_list$Expansion_cz[cat.to.use,]
smalls$data_list$Options_list$metadata_ctz <- smalls$data_list$Options_list$metadata_ctz[cat.to.use,,]
smalls$data_list$ObsModel_ez <- smalls$data_list$ObsModel_ez[cat.to.use,]
smalls$data_list$X1config_cp <- smalls$data_list$X1config_cp[cat.to.use,]
smalls$data_list$X2config_cp <- smalls$data_list$X2config_cp[cat.to.use,]
smalls$data_list$b_i <- smalls$data_list$b_i[seq(1, length(smalls$data_list$b_i), ncat)]
smalls$data_list$a_i <- smalls$data_list$a_i[seq(1, length(smalls$data_list$a_i), ncat)]
smalls$data_list$c_iz <- matrix(ncol=1, nrow=length(smalls$data_list$a_i),
                                data=rep(0, length(smalls$data_list$a_i)))
smalls$data_list$e_i <- rep(0, length(smalls$data_list$a_i))
smalls$data_list$t_i <- smalls$data_list$t_i[seq(1, length(smalls$data_list$t_i), ncat)]
smalls$data_list$v_i <- rep(0, length(smalls$data_list$a_i))
smalls$data_list$PredTF_i <- rep(0, length(smalls$data_list$a_i))
smalls$data_list$X1_ip <- smalls$data_list$X1_ip[seq(1, nrow(smalls$data_list$X1_ip), ncat), ]
smalls$data_list$X1_gctp <- smalls$data_list$X1_gctp[,1,,]
smalls$data_list$X2_ip <- smalls$data_list$X2_ip[seq(1, nrow(smalls$data_list$X2_ip), ncat), ]
smalls$data_list$X2_gctp <- smalls$data_list$X2_gctp[,1,,]
smalls$data_list$Q1_ik <- smalls$data_list$Q1_ik[seq(1, nrow(smalls$data_list$Q1_ik), ncat),]
smalls$data_list$Q2_ik <- smalls$data_list$Q2_ik[seq(1, nrow(smalls$data_list$Q2_ik), ncat),]

tmb.par <- smalls$tmb_list$Obj$par
tmb.par.beta1_ft <- tmb.par[names(tmb.par) == 'beta1_ft']
tmb.par.beta1_ft <- tmb.par.beta1_ft[seq(1, length(tmb.par.beta1_ft), ncat)]
tmb.par.beta2_ft <- tmb.par[names(tmb.par) == 'beta2_ft']
tmb.par.beta2_ft <- tmb.par.beta2_ft[seq(1, length(tmb.par.beta2_ft), ncat)]

tmb.par.gamma1_cp <- tmb.par[names(tmb.par) == 'gamma1_cp']
tmb.par.gamma1_cp <- tmb.par.gamma1_cp[seq(1, length(tmb.par.gamma1_cp), ncat)]

tmb.par.gamma2_cp <- tmb.par[names(tmb.par) == 'gamma2_cp']
tmb.par.gamma2_cp <- tmb.par.gamma2_cp[seq(1, length(tmb.par.gamma2_cp), ncat)]

tmb.par.L_epsilon1_z <- tmb.par[names(tmb.par) == 'L_epsilon1_z']
tmb.par.L_epsilon1_z <- tmb.par.L_epsilon1_z[seq(1, length(tmb.par.L_epsilon1_z), ncat)]

tmb.par.L_epsilon2_z <- tmb.par[names(tmb.par) == 'L_epsilon2_z']
tmb.par.L_epsilon2_z <- tmb.par.L_epsilon2_z[seq(1, length(tmb.par.L_epsilon2_z), ncat)]

tmb.par.L_omega1_z <- tmb.par[names(tmb.par) == 'L_omega1_z']
tmb.par.L_omega1_z <- tmb.par.L_omega1_z[seq(1, length(tmb.par.L_omega1_z), ncat)]

tmb.par.L_omega2_z <- tmb.par[names(tmb.par) == 'L_omega2_z']
tmb.par.L_omega2_z <- tmb.par.L_omega2_z[seq(1, length(tmb.par.L_omega2_z), ncat)]

tmb.par.ln_H_input <- tmb.par[names(tmb.par) == 'ln_H_input']

tmb.par.logkappa1 <- tmb.par[names(tmb.par) == 'logkappa1']
tmb.par.logkappa2 <- tmb.par[names(tmb.par) == 'logkappa2']

tmb.par.logSigmaM <- tmb.par[names(tmb.par) == 'logSigmaM']
tmb.par.logSigmaM <- tmb.par.logSigmaM[seq(1, length(tmb.par.logSigmaM), ncat)]

tmb.par <- c(tmb.par.ln_H_input, 
             tmb.par.beta1_ft, tmb.par.gamma1_cp, tmb.par.L_omega1_z, tmb.par.L_epsilon1_z, tmb.par.logkappa1,
             tmb.par.beta2_ft, tmb.par.gamma2_cp, tmb.par.L_omega2_z, tmb.par.L_epsilon2_z, tmb.par.logkappa2,
             tmb.par.logSigmaM)
smalls$tmb_list$Obj$par <- tmb.par

tmb.par <- smalls$tmb_list$Upper
tmb.par.beta1_ft <- tmb.par[names(tmb.par) == 'beta1_ft']
tmb.par.beta1_ft <- tmb.par.beta1_ft[seq(1, length(tmb.par.beta1_ft), ncat)]
tmb.par.beta2_ft <- tmb.par[names(tmb.par) == 'beta2_ft']
tmb.par.beta2_ft <- tmb.par.beta2_ft[seq(1, length(tmb.par.beta2_ft), ncat)]

tmb.par.gamma1_cp <- tmb.par[names(tmb.par) == 'gamma1_cp']
tmb.par.gamma1_cp <- tmb.par.gamma1_cp[seq(1, length(tmb.par.gamma1_cp), ncat)]

tmb.par.gamma2_cp <- tmb.par[names(tmb.par) == 'gamma2_cp']
tmb.par.gamma2_cp <- tmb.par.gamma2_cp[seq(1, length(tmb.par.gamma2_cp), ncat)]

tmb.par.L_epsilon1_z <- tmb.par[names(tmb.par) == 'L_epsilon1_z']
tmb.par.L_epsilon1_z <- tmb.par.L_epsilon1_z[seq(1, length(tmb.par.L_epsilon1_z), ncat)]

tmb.par.L_epsilon2_z <- tmb.par[names(tmb.par) == 'L_epsilon2_z']
tmb.par.L_epsilon2_z <- tmb.par.L_epsilon2_z[seq(1, length(tmb.par.L_epsilon2_z), ncat)]

tmb.par.L_omega1_z <- tmb.par[names(tmb.par) == 'L_omega1_z']
tmb.par.L_omega1_z <- tmb.par.L_omega1_z[seq(1, length(tmb.par.L_omega1_z), ncat)]

tmb.par.L_omega2_z <- tmb.par[names(tmb.par) == 'L_omega2_z']
tmb.par.L_omega2_z <- tmb.par.L_omega2_z[seq(1, length(tmb.par.L_omega2_z), ncat)]

tmb.par.ln_H_input <- tmb.par[names(tmb.par) == 'ln_H_input']

tmb.par.logkappa1 <- tmb.par[names(tmb.par) == 'logkappa1']
tmb.par.logkappa2 <- tmb.par[names(tmb.par) == 'logkappa2']

tmb.par.logSigmaM <- tmb.par[names(tmb.par) == 'logSigmaM']
tmb.par.logSigmaM <- tmb.par.logSigmaM[seq(1, length(tmb.par.logSigmaM), ncat)]

tmb.par <- c(tmb.par.ln_H_input, 
             tmb.par.beta1_ft, tmb.par.gamma1_cp, tmb.par.L_omega1_z, tmb.par.L_epsilon1_z, tmb.par.logkappa1,
             tmb.par.beta2_ft, tmb.par.gamma2_cp, tmb.par.L_omega2_z, tmb.par.L_epsilon2_z, tmb.par.logkappa2,
             tmb.par.logSigmaM)
smalls$tmb_list$Upper <- tmb.par

tmb.par <- smalls$tmb_list$Lower
tmb.par.beta1_ft <- tmb.par[names(tmb.par) == 'beta1_ft']
tmb.par.beta1_ft <- tmb.par.beta1_ft[seq(1, length(tmb.par.beta1_ft), ncat)]
tmb.par.beta2_ft <- tmb.par[names(tmb.par) == 'beta2_ft']
tmb.par.beta2_ft <- tmb.par.beta2_ft[seq(1, length(tmb.par.beta2_ft), ncat)]

tmb.par.gamma1_cp <- tmb.par[names(tmb.par) == 'gamma1_cp']
tmb.par.gamma1_cp <- tmb.par.gamma1_cp[seq(1, length(tmb.par.gamma1_cp), ncat)]

tmb.par.gamma2_cp <- tmb.par[names(tmb.par) == 'gamma2_cp']
tmb.par.gamma2_cp <- tmb.par.gamma2_cp[seq(1, length(tmb.par.gamma2_cp), ncat)]

tmb.par.L_epsilon1_z <- tmb.par[names(tmb.par) == 'L_epsilon1_z']
tmb.par.L_epsilon1_z <- tmb.par.L_epsilon1_z[seq(1, length(tmb.par.L_epsilon1_z), ncat)]

tmb.par.L_epsilon2_z <- tmb.par[names(tmb.par) == 'L_epsilon2_z']
tmb.par.L_epsilon2_z <- tmb.par.L_epsilon2_z[seq(1, length(tmb.par.L_epsilon2_z), ncat)]

tmb.par.L_omega1_z <- tmb.par[names(tmb.par) == 'L_omega1_z']
tmb.par.L_omega1_z <- tmb.par.L_omega1_z[seq(1, length(tmb.par.L_omega1_z), ncat)]

tmb.par.L_omega2_z <- tmb.par[names(tmb.par) == 'L_omega2_z']
tmb.par.L_omega2_z <- tmb.par.L_omega2_z[seq(1, length(tmb.par.L_omega2_z), ncat)]

tmb.par.ln_H_input <- tmb.par[names(tmb.par) == 'ln_H_input']

tmb.par.logkappa1 <- tmb.par[names(tmb.par) == 'logkappa1']
tmb.par.logkappa2 <- tmb.par[names(tmb.par) == 'logkappa2']

tmb.par.logSigmaM <- tmb.par[names(tmb.par) == 'logSigmaM']
tmb.par.logSigmaM <- tmb.par.logSigmaM[seq(1, length(tmb.par.logSigmaM), ncat)]

tmb.par <- c(tmb.par.ln_H_input, 
             tmb.par.beta1_ft, tmb.par.gamma1_cp, tmb.par.L_omega1_z, tmb.par.L_epsilon1_z, tmb.par.logkappa1,
             tmb.par.beta2_ft, tmb.par.gamma2_cp, tmb.par.L_omega2_z, tmb.par.L_epsilon2_z, tmb.par.logkappa2,
             tmb.par.logSigmaM)
smalls$tmb_list$Lower <- tmb.par

smalls$tmb_list$Parameters$beta1_ft <- t(as.matrix(smalls$tmb_list$Parameters$beta1_ft[1,]))
smalls$tmb_list$Parameters$gamma1_cp<- t(as.matrix(smalls$tmb_list$Parameters$gamma1_cp[1,]))
smalls$tmb_list$Parameters$L_omega1_z <- smalls$tmb_list$Parameters$L_omega1_z[1]
smalls$tmb_list$Parameters$L_epsilon1_z <- smalls$tmb_list$Parameters$L_epsilon1_z[1]
smalls$tmb_list$Parameters$L_beta1_z <- smalls$tmb_list$Parameters$L_beta1_z[1]
smalls$tmb_list$Parameters$Beta_mean1_c <- smalls$tmb_list$Parameters$Beta_mean1_c[1]
smalls$tmb_list$Parameters$Beta_rho1_f <- smalls$tmb_list$Parameters$Beta_rho1_f[1]
smalls$tmb_list$Parameters$Epsilon_rho1_f <- smalls$tmb_list$Parameters$Epsilon_rho1_f[1]
smalls$tmb_list$Parameters$log_sigmaXi1_cp <- t(as.matrix(smalls$tmb_list$Parameters$log_sigmaXi1_cp[1,]))
smalls$tmb_list$Parameters$Xiinput1_scp <- smalls$tmb_list$Parameters$Xiinput1_scp[,1,]
smalls$tmb_list$Parameters$Omegainput1_sf <- as.matrix(smalls$tmb_list$Parameters$Omegainput1_sf[,1])
smalls$tmb_list$Parameters$Epsiloninput1_sff <- smalls$tmb_list$Parameters$Epsiloninput1_sff[,1,]

smalls$tmb_list$Parameters$beta2_ft <- t(as.matrix(smalls$tmb_list$Parameters$beta2_ft[1,]))
smalls$tmb_list$Parameters$gamma2_cp<- t(as.matrix(smalls$tmb_list$Parameters$gamma2_cp[1,]))
smalls$tmb_list$Parameters$L_omega2_z <- smalls$tmb_list$Parameters$L_omega2_z[1]
smalls$tmb_list$Parameters$L_epsilon2_z <- smalls$tmb_list$Parameters$L_epsilon2_z[1]
smalls$tmb_list$Parameters$L_beta2_z <- smalls$tmb_list$Parameters$L_beta2_z[1]
smalls$tmb_list$Parameters$Beta_mean2_c <- smalls$tmb_list$Parameters$Beta_mean2_c[1]
smalls$tmb_list$Parameters$Beta_rho2_f <- smalls$tmb_list$Parameters$Beta_rho2_f[1]
smalls$tmb_list$Parameters$Epsilon_rho2_f <- smalls$tmb_list$Parameters$Epsilon_rho2_f[1]
smalls$tmb_list$Parameters$log_sigmaXi2_cp <- t(as.matrix(smalls$tmb_list$Parameters$log_sigmaXi2_cp[1,]))
smalls$tmb_list$Parameters$Xiinput2_scp <- smalls$tmb_list$Parameters$Xiinput2_scp[,1,]
smalls$tmb_list$Parameters$Omegainput2_sf <- as.matrix(smalls$tmb_list$Parameters$Omegainput2_sf[,1])
smalls$tmb_list$Parameters$Epsiloninput2_sff <- smalls$tmb_list$Parameters$Epsiloninput2_sff[,1,]

smalls$tmb_list$Parameters$logSigmaM <- t(as.matrix(smalls$tmb_list$Parameters$logSigmaM[1,]))
smalls$tmb_list$Parameters$lagrange_tc <- as.matrix(smalls$tmb_list$Parameters$lagrange_tc[,1])

smalls$tmb_list$Map$logSigmaM <- factor(c(1, NA, NA))
smalls$tmb_list$Map$Epsilon_rho1_f <- factor(NA)
smalls$tmb_list$Map$Epsilon_rho2_f <- factor(NA)
smalls$tmb_list$Map$gamma1_cp <- factor(c(1,2,3,4,5,6,7))
smalls$tmb_list$Map$gamma2_cp <- factor(c(1,2,3,4,5,6,7))
smalls$tmb_list$Map$log_sigmaXi1_cp <- factor(c(NA, NA, NA, NA, NA, NA, NA))
smalls$tmb_list$Map$log_sigmaXi2_cp <- factor(c(NA, NA, NA, NA, NA, NA, NA))
smalls$tmb_list$Map$lagrange_tc <- factor(rep(NA, 80))
smalls$tmb_list$Map$Beta_mean1_c <- factor(NA)
smalls$tmb_list$Map$Beta_rho1_f <- factor(NA)
smalls$tmb_list$Map$L_beta1_z <- factor(NA)
smalls$tmb_list$Map$Beta_mean2_c <- factor(NA)
smalls$tmb_list$Map$Beta_rho2_f <- factor(NA)
smalls$tmb_list$Map$L_beta2_z <- factor(NA)

tmb.par <- smalls$parameter_estimates$diagnostics
tmb.par.beta1_ft <- subset(tmb.par, Param == 'beta1_ft')
tmb.par.beta1_ft <- tmb.par.beta1_ft[seq(1, nrow(tmb.par.beta1_ft), ncat),]

tmb.par.beta2_ft <- subset(tmb.par, Param == 'beta2_ft')
tmb.par.beta2_ft <- tmb.par.beta2_ft[seq(1, nrow(tmb.par.beta2_ft), ncat),]

tmb.par.gamma1_cp <- subset(tmb.par, Param == 'gamma1_cp')
tmb.par.gamma1_cp <- tmb.par.gamma1_cp[seq(1, nrow(tmb.par.gamma1_cp), ncat),]

tmb.par.gamma2_cp <- subset(tmb.par, Param == 'gamma2_cp')
tmb.par.gamma2_cp <- tmb.par.gamma2_cp[seq(1, nrow(tmb.par.gamma2_cp), ncat),]

tmb.par.L_epsilon1_z <- subset(tmb.par, Param == 'L_epsilon1_z')
tmb.par.L_epsilon1_z <- tmb.par.L_epsilon1_z[seq(1, nrow(tmb.par.L_epsilon1_z), ncat),]

tmb.par.L_epsilon2_z <- subset(tmb.par, Param == 'L_epsilon2_z')
tmb.par.L_epsilon2_z <- tmb.par.L_epsilon2_z[seq(1, nrow(tmb.par.L_epsilon2_z), ncat),]

tmb.par.L_omega1_z <- subset(tmb.par, Param == 'L_omega1_z')
tmb.par.L_omega1_z <- tmb.par.L_omega1_z[seq(1, nrow(tmb.par.L_omega1_z), ncat),]

tmb.par.L_omega2_z <- subset(tmb.par, Param == 'L_omega2_z')
tmb.par.L_omega2_z <- tmb.par.L_omega2_z[seq(1, nrow(tmb.par.L_omega2_z), ncat),]

tmb.par.ln_H_input <- subset(tmb.par, Param == 'ln_H_input')

tmb.par.logkappa1 <- subset(tmb.par, Param == 'logkappa1')
tmb.par.logkappa2 <- subset(tmb.par, Param == 'logkappa2')

tmb.par.logSigmaM <- subset(tmb.par, Param == 'logSigmaM')
tmb.par.logSigmaM <- tmb.par.logSigmaM[seq(1, nrow(tmb.par.logSigmaM), ncat),]

tmb.par <- rbind(tmb.par.ln_H_input, 
             tmb.par.beta1_ft, tmb.par.gamma1_cp, tmb.par.L_omega1_z, tmb.par.L_epsilon1_z, tmb.par.logkappa1,
             tmb.par.beta2_ft, tmb.par.gamma2_cp, tmb.par.L_omega2_z, tmb.par.L_epsilon2_z, tmb.par.logkappa2,
             tmb.par.logSigmaM)
smalls$parameter_estimates$diagnostics <- tmb.par

tmb.par <- smalls$parameter_estimates$par
tmb.par.beta1_ft <- tmb.par[names(tmb.par) == 'beta1_ft']
tmb.par.beta1_ft <- tmb.par.beta1_ft[seq(1, length(tmb.par.beta1_ft), ncat)]
tmb.par.beta2_ft <- tmb.par[names(tmb.par) == 'beta2_ft']
tmb.par.beta2_ft <- tmb.par.beta2_ft[seq(1, length(tmb.par.beta2_ft), ncat)]

tmb.par.gamma1_cp <- tmb.par[names(tmb.par) == 'gamma1_cp']
tmb.par.gamma1_cp <- tmb.par.gamma1_cp[seq(1, length(tmb.par.gamma1_cp), ncat)]

tmb.par.gamma2_cp <- tmb.par[names(tmb.par) == 'gamma2_cp']
tmb.par.gamma2_cp <- tmb.par.gamma2_cp[seq(1, length(tmb.par.gamma2_cp), ncat)]

tmb.par.L_epsilon1_z <- tmb.par[names(tmb.par) == 'L_epsilon1_z']
tmb.par.L_epsilon1_z <- tmb.par.L_epsilon1_z[seq(1, length(tmb.par.L_epsilon1_z), ncat)]

tmb.par.L_epsilon2_z <- tmb.par[names(tmb.par) == 'L_epsilon2_z']
tmb.par.L_epsilon2_z <- tmb.par.L_epsilon2_z[seq(1, length(tmb.par.L_epsilon2_z), ncat)]

tmb.par.L_omega1_z <- tmb.par[names(tmb.par) == 'L_omega1_z']
tmb.par.L_omega1_z <- tmb.par.L_omega1_z[seq(1, length(tmb.par.L_omega1_z), ncat)]

tmb.par.L_omega2_z <- tmb.par[names(tmb.par) == 'L_omega2_z']
tmb.par.L_omega2_z <- tmb.par.L_omega2_z[seq(1, length(tmb.par.L_omega2_z), ncat)]

tmb.par.ln_H_input <- tmb.par[names(tmb.par) == 'ln_H_input']

tmb.par.logkappa1 <- tmb.par[names(tmb.par) == 'logkappa1']
tmb.par.logkappa2 <- tmb.par[names(tmb.par) == 'logkappa2']

tmb.par.logSigmaM <- tmb.par[names(tmb.par) == 'logSigmaM']
tmb.par.logSigmaM <- tmb.par.logSigmaM[seq(1, length(tmb.par.logSigmaM), ncat)]

tmb.par <- c(tmb.par.ln_H_input, 
             tmb.par.beta1_ft, tmb.par.gamma1_cp, tmb.par.L_omega1_z, tmb.par.L_epsilon1_z, tmb.par.logkappa1,
             tmb.par.beta2_ft, tmb.par.gamma2_cp, tmb.par.L_omega2_z, tmb.par.L_epsilon2_z, tmb.par.logkappa2,
             tmb.par.logSigmaM)
smalls$parameter_estimates$par <- tmb.par

smalls$parameter_estimates$SD$value

vast_habcovs_effs<- get_vast_covariate_effects(vast_fit = fit, 
                                               params_plot = c("BATHY.DEPTH"), 
                                               params_plot_levels = 100, effects_pad_values = c(),
                                               nice_category_names = "Atlantic halibut", 
                                               out_dir = here::here("", "results/tables"))

# Warnings are...interesting....
vast_habcovs_plot<- plot_vast_covariate_effects(vast_covariate_effects = vast_habcovs_effs, 
                                                vast_fit = vast_fitted_hab_covs, 
                                                nice_category_names = "Atlantic halibut", 
                                                out_dir = here::here("", "results/plots_maps"))
vast_habcovs_plot