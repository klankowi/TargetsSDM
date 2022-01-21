#####
## Making predictions from fitted VAST model with new data
# The goal of this short example is to show how we can extract what we need from a fitted VAST model to make predictions using new data. In particular, we are looking to generate a dataframe that has a row for each knot/time of interest. Then, we want to have columns with the components needed to calculate `P1` and `P2` as in the VAST.cpp files. For example, with the first linear predictor we are going to want to have Omega1_gc(g,c) + beta1_tc(t,c) + Epsilon1_gct(g,c,t) + eta1_gct(g,c,t).
#####

# Load packages
library(VAST)
library(splines)
library(tidyverse)
library(sf)

# load data set
# see `?load_example` for list of stocks with example data
# that are installed automatically with `FishStatsUtils`.
example = load_example( data_set="covariate_example" )

# Test settings
field_config_use<- c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 1) # Modify here to narrow down where issue is occurring...

# Log and rescale BOT_DEPTH now...
example$covariate_data[,'BOT_DEPTH']<- log(example$covariate_data[,'BOT_DEPTH'])
# Rescale covariates being used to have an SD >0.1 and <10 (for numerical stability)
example$covariate_data[,'BOT_DEPTH'] = example$covariate_data[,'BOT_DEPTH'] / 100

# Make settings (turning off bias.correct to save time for example). **NOTE: When spatial/spatio-temporal components are turned off, manual predictions match those from `predict_fit.model`
settings = make_settings( n_x=200,
                          ObsModel = c(2,1),
                          Region=example$Region,
                          purpose="index2",
                          use_anisotropy=FALSE,
                          bias.correct=FALSE,
                          fine_scale=TRUE)
settings$FieldConfig<- field_config_use
settings$RhoConfig<- c("Beta1" = 1, "Beta2" = 1, "Epsilon1" = 2, "Epsilon2" = 2)

# Define formula -- a simple one for now
X1_formula = ~ bs( BOT_DEPTH, degree=2, intercept=FALSE)
X2_formula = X1_formula

# In covariate example, all years are set to NA as depth doesn't vary temporally. This creates an issue with predict_fit.model as it wants to have lat, lon, year columns. If we just proceed, `fit_model` complains because years aren't sequential. So, just adjusting these to be sequential for this example. I'm sure there's a faster way with a look up table, brain was working slowly...
example$covariate_data[,'Year']<- ifelse(example$covariate_data[,'Year'] == 2001, 2001,
                                             ifelse(example$covariate_data[,'Year'] == 2003, 2002,
                                                    ifelse(example$covariate_data[,'Year'] == 2005, 2003,
                                                           ifelse(example$covariate_data[,'Year'] == 2007, 2004,
                                                                  ifelse(example$covariate_data[,'Year'] == 2009, 2005, 
                                                                         ifelse(example$covariate_data[,'Year'] == 2011, 2006,
                                                                                ifelse(example$covariate_data[,'Year'] == 2013, 2007,
                                                                                       ifelse(example$covariate_data[,'Year'] == 2015, 2008,
                                                                                              ifelse(example$covariate_data[,'Year'] == 2017, 2009, NA)))))))))
example$covariate_data[,'Year']<- as.numeric(factor(example$covariate_data[,'Year']))-1

example$sampling_data[,'Year']<- ifelse(example$sampling_data[,'Year'] == 2001, 2001,
                                             ifelse(example$sampling_data[,'Year'] == 2003, 2002,
                                                    ifelse(example$sampling_data[,'Year'] == 2005, 2003,
                                                           ifelse(example$sampling_data[,'Year'] == 2007, 2004,
                                                                  ifelse(example$sampling_data[,'Year'] == 2009, 2005, 
                                                                         ifelse(example$sampling_data[,'Year'] == 2011, 2006,
                                                                                ifelse(example$sampling_data[,'Year'] == 2013, 2007,
                                                                                       ifelse(example$sampling_data[,'Year'] == 2015, 2008,
                                                                                              ifelse(example$sampling_data[,'Year'] == 2017, 2009, NA)))))))))
example$sampling_data[,'Year']<- as.numeric(factor(example$sampling_data[,'Year']))-1


# Run model
fit = fit_model( "settings" = settings,
                 Lat_i = example$sampling_data[,'Lat'],
                 Lon_i = example$sampling_data[,'Lon'],
                 t_i = example$sampling_data[,'Year'],
                 b_i = example$sampling_data[,'Catch_KG'],
                 a_i = example$sampling_data[,'AreaSwept_km2'],
                 X1_formula = X1_formula,
                 X2_formula = X2_formula,
                 covariate_data = example$covariate_data,
                 run_model = FALSE)

# Predictions using `predict.fit_model`
object = fit
p1_true<- predict( x=object,
         Lat_i = object$data_frame$Lat_i,
         Lon_i = object$data_frame$Lon_i,
         t_i = object$data_frame$t_i,
         a_i = object$data_frame$a_i,
         what = "P1_iz",
         new_covariate_data = example$covariate_data,
         do_checks = FALSE )
p2_true<- predict( x=object,
                   Lat_i = object$data_frame$Lat_i,
                   Lon_i = object$data_frame$Lon_i,
                   t_i = object$data_frame$t_i,
                   a_i = object$data_frame$a_i,
                   what = "P2_iz",
                   new_covariate_data = example$covariate_data,
                   do_checks = FALSE )
d_true<- predict( x=object,
                  Lat_i = object$data_frame$Lat_i,
                  Lon_i = object$data_frame$Lon_i,
                  t_i = object$data_frame$t_i,
                  a_i = object$data_frame$a_i,
                  what = "D_i",
                  new_covariate_data = example$covariate_data,
                  do_checks = FALSE )


# Manual predictions...long working
new_covariate_data = example$covariate_data
all_times = unique(example$covariate_data[,'Year'])

# Here, we are fixing everything else in the model and only looking at the influence of environmental changes on species occurrence.
# Extracting what we need from the fitted model object. This will depend on the "use_years" -- if more than one year are supplied, then the average for the intercept/spatio-temporal components are taken. If just one year, than that year's estimate is used.
preds_out<- data.frame("Date" = factor(new_covariate_data$Year), "Obs" = seq(from = 1, to = fit$data_list$n_i, by = 1), "Lat" = new_covariate_data$Lat, "Lon" = new_covariate_data$Lon)
                       
## Random effects 
# Omega first, this will be a surface and is constant through time. These values are going to change depending on if the fine_scale option is used. When they are used, values at knot locations (or grid locations) are estimated using a bilinear interpolation. When they are not, just nearest neighbor are used. To account for this, we just need to make sure to update value at locations as we go. With the A_is matrix, slot i corresponds to the observation ID, slot j to the knot ID, and slot x to the value. 
omega1_proj<- matrix(data = 0, nrow = fit$data_list$n_i, ncol = 1)
for(i in seq_along(fit$spatial_list$A_is@i)){
  obs_ind<- fit$spatial_list$A_is@i[i]+1
  knot_ind<- fit$spatial_list$A_is@j[i]+1
  omega1_proj[obs_ind,1]<- omega1_proj[obs_ind,1] + (fit$spatial_list$A_is@x[i] * fit$Report[[{{"Omega1_sc"}}]][knot_ind,1])
}

preds_out$omega1<- omega1_proj[,1]

omega2_proj<- matrix(data = 0, nrow = fit$data_list$n_i, ncol = 1)
for(i in seq_along(fit$spatial_list$A_is@i)){
  obs_ind<- fit$spatial_list$A_is@i[i]+1
  knot_ind<- fit$spatial_list$A_is@j[i]+1
  omega2_proj[obs_ind,1]<- omega2_proj[obs_ind,1] + (fit$spatial_list$A_is@x[i] * fit$Report[[{{"Omega2_sc"}}]][knot_ind,1])
}

preds_out$omega2<- omega2_proj[,1]

# Average intercepts, these are constant across knots (sometimes, time too)
beta1<- data.frame("Date" = factor(all_times), "beta1" = fit$Report[[{{"beta1_tc"}}]])
beta2<- data.frame("Date" = factor(all_times), "beta2" = fit$Report[[{{"beta2_tc"}}]])
betas<- beta1 %>%
  left_join(., beta2)

preds_out<- preds_out %>%
  left_join(., betas, by = c("Date" = "Date"))

# Spatio_temporal -- similar to omega, but now with a temporal dimension...
eps1_proj<- matrix(data = 0, nrow = fit$data_list$n_i, ncol = fit$data_list$n_t)
for(i in seq_along(fit$spatial_list$A_is@i)){
  for(t in 1:fit$data_list$n_t){
    obs_ind<- fit$spatial_list$A_is@i[i]+1
    knot_ind<- fit$spatial_list$A_is@j[i]+1
    eps1_proj[obs_ind,t]<- eps1_proj[obs_ind,t] + (fit$spatial_list$A_is@x[i] * fit$Report[[{{"Epsilon1_sct"}}]][knot_ind,1,t])
  }
}

eps2_proj<- matrix(data = 0, nrow = fit$data_list$n_i, ncol = fit$data_list$n_t)
for(i in seq_along(fit$spatial_list$A_is@i)){
  for(t in 1:fit$data_list$n_t){
    obs_ind<- fit$spatial_list$A_is@i[i]+1
    knot_ind<- fit$spatial_list$A_is@j[i]+1
    eps2_proj[obs_ind,t]<- eps2_proj[obs_ind,t] + (fit$spatial_list$A_is@x[i] * fit$Report[[{{"Epsilon2_sct"}}]][knot_ind,1,t])
  }
}

# Some work on column names...
eps1_wide<- data.frame("Obs" = seq(from = 1, to = fit$data_list$n_i, by = 1), eps1_proj)
colnames(eps1_wide)<- c("Obs", all_times)
eps1_long<- eps1_wide %>%
  pivot_longer(-Obs, names_to = "Date", values_to = "epsilon1") %>%
  mutate(., "Date" = factor(Date, levels = all_times)) %>%
  arrange(Date)

eps2_wide<- data.frame("Obs" = seq(from = 1, to = fit$data_list$n_i, by = 1), eps2_proj)
colnames(eps2_wide)<- c("Obs", all_times)
eps2_long<- eps2_wide %>%
  pivot_longer(-Obs, names_to = "Date", values_to = "epsilon2") %>%
  mutate(., "Date" = factor(Date, levels = all_times)) %>%
  arrange(Date)

# Combine them
eps_out<- left_join(eps1_long, eps2_long)

# Add to preds...
preds_out<- preds_out %>%
  left_join(., eps_out, by = c("Obs" = "Obs", "Date" = "Date"))

## Habitat covariates using gridded environmental data
# For this, I think it is going to be easier to do the predictions and then supply those single eta1 and eta2 values to a simpler function
x1_covs<- fit$ParHat[[{{"gamma1_cp"}}]]
colnames(x1_covs)<- attributes(fit$data_list$X1_gctp)$dimnames[[4]]
x2_covs<- fit$ParHat[[{{"gamma2_cp"}}]]
colnames(x2_covs)<- attributes(fit$data_list$X2_gctp)$dimnames[[4]]

# Create a model matrix
mod_mat<- model.matrix(object = fit$X1_formula, data = new_covariate_data)

# Not necessary in this simple example, but matters when we have new years in the new_covariate_data
cols_keep<- colnames(mod_mat) %in% colnames(x1_covs)
mod_mat<- mod_mat[,cols_keep]

# Generate predictions from covariate values and model matrix
preds_out$eta1<- as.numeric(t(as.vector(x1_covs) %*% t(mod_mat)))
preds_out$eta2<- as.numeric(t(as.vector(x2_covs) %*% t(mod_mat)))

## Make the calculations...
preds_out<- preds_out %>% 
  mutate(., "P1" = omega1 + beta1 + epsilon1 + eta1,
         "P2" = omega2 + beta2 + epsilon2 + eta2,
         # "R1" = exp(P1)/(1+exp(P1)), #DLogNorm implementation
         # "R2" = exp(P2), #DLogNorm implementation
         "R1" = 1- exp(-exp(P1)), # Poisson-link implementation
         "R2" = (exp(P1)/R1) * exp(P2), # Poisson-link implementation
         "Dens" = R1*R2)

# Check
all.equal(p1_true, preds_out$P1)
all.equal(d_true, preds_out$Dens) # Without spatial/spatio-temporal variation, density values match but off in units 
head(d_true)
head(preds_out$Dens)
