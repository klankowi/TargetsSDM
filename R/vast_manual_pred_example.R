#####
## Making predictions from fitted VAST model with new data
# The goal of this short example is to show how we can extract what we need from a fitted VAST model to make predictions using new data. In particular, we are looking to generate a dataframe that has a row for each knot/time of interest. Then, we want to have columns with the components needed to calculate `P1` and `P2` as in the VAST.cpp files. For example, with the first linear predictor we are going to want to have Omega1_gc(g,c) + beta1_tc(t,c) + Epsilon1_gct(g,c,t) + eta1_gct(g,c,t).
#####

# Load packages
library(VAST)
library(splines)

# load data set
# see `?load_example` for list of stocks with example data
# that are installed automatically with `FishStatsUtils`.
example = load_example( data_set="covariate_example" )

# Test settings
field_config_use<- c("Omega1" = 0, "Epsilon1" = 0, "Omega2" = 0, "Epsilon2" = 0) # Modify here to narrow down where issue is occurring...

# Log BOT_DEPTH now...
example$covariate_data[,'BOT_DEPTH']<- log(example$covariate_data[,'BOT_DEPTH'])
# Rescale covariates being used to have an SD >0.1 and <10 (for numerical stability)
example$covariate_data[,'BOT_DEPTH'] = example$covariate_data[,'BOT_DEPTH'] / 100

# Make settings (turning off bias.correct to save time for example). **NOTE: When spatial/spatio-temporal components are turned off, manual predictions match those from `predict_fit.model`
settings = make_settings( n_x=200,
                          ObsModel = c(1,0),
                          Region=example$Region,
                          purpose="index2",
                          use_anisotropy=FALSE,
                          bias.correct=FALSE,
                          fine_scale=TRUE )
settings$FieldConfig<- field_config_use

# Define formula -- a simple one for now
X1_formula = ~ bs( BOT_DEPTH, degree=2, intercept=FALSE)
X2_formula = X1_formula

# NOW necessary to duplicate across years
year_set = min(example$covariate_data[,'Year']):max(example$covariate_data[,'Year'])
tmp = example$covariate_data
for( year in year_set ){
  if(year == year_set[1]){
    example$covariate_data[,'Year'] = year
  }else{
    tmp[,'Year'] = year
    example$covariate_data = rbind(example$covariate_data, tmp)
  }
}

# Make new factor
example$covariate_data = cbind( example$covariate_data,
                                "Year_factor"=factor(example$covariate_data[,'Year'],levels=sort(unique(example$covariate_data[,'Year']))) )


# Run model
fit = fit_model( "settings" = settings,
                 Lat_i = example$sampling_data[,'Lat'],
                 Lon_i = example$sampling_data[,'Lon'],
                 t_i = example$sampling_data[,'Year'],
                 b_i = example$sampling_data[,'Catch_KG'],
                 a_i = example$sampling_data[,'AreaSwept_km2'],
                 X1_formula = X1_formula,
                 X2_formula = X2_formula,
                 covariate_data = example$covariate_data )

# Predictions using `predict.fit_model`
object = fit
p1_true<- predict( x=object,
         Lat_i = object$data_frame$Lat_i,
         Lon_i = object$data_frame$Lon_i,
         t_i = object$data_frame$t_i,
         a_i = object$data_frame$a_i,
         what = "P1_i",
         new_covariate_data = example$covariate_data,
         do_checks = FALSE )
dens_true<- predict( x=object,
                  Lat_i = object$data_frame$Lat_i,
                  Lon_i = object$data_frame$Lon_i,
                  t_i = object$data_frame$t_i,
                  a_i = object$data_frame$a_i,
                  what = "D_gct",
                  new_covariate_data = example$covariate_data,
                  do_checks = FALSE )


# Manual predictions...
new_covariate_data = example$covariate_data
all_times = unique(example$covariate_data[,'Year'])

# Here, we are fixing everything else in the model and only looking at the influence of environmental changes on species occurrence.
# Extracting what we need from the fitted model object. This will depend on the "use_years" -- if more than one year are supplied, then the average for the intercept/spatio-temporal components are taken. If just one year, than that year's estimate is used.
## Estimates at knot locations (s)
knot_preds_out<- data.frame("Knot" = seq(from = 1, to = fit$data_list$n_s, by = 1))

# Omega first, this will be a surface and is constant through time..
knot_preds_out$omega1<- fit$Report[[{{"Omega1_sc"}}]][,1]
knot_preds_out$omega2<- fit$Report[[{{"Omega2_sc"}}]][,1]

# Average intercepts, these are constant across knots (sometimes, time too)
beta1<- data.frame("Date" = factor(all_times), "beta1" = fit$Report[[{{"beta1_tc"}}]])
beta2<- data.frame("Date" = factor(all_times), "beta2" = fit$Report[[{{"beta2_tc"}}]])
betas<- merge(beta1, beta2)

# Spatio_temporal array.
eps1_wide<- data.frame("Knot" = seq(from = 1, to = fit$data_list$n_s, by = 1), fit$Report[[{{"Epsilon1_sct"}}]])
colnames(eps1_wide)<- c("Knot", all_times)
eps2_wide<- data.frame("Knot" = seq(from = 1, to = fit$data_list$n_s, by = 1), fit$Report[[{{"Epsilon2_sct"}}]])
colnames(eps2_wide)<- c("Knot", all_times)

# Long dataframe...
eps1_long<- data.frame(reshape(eps1_wide, direction = "long", varying = names(eps1_wide)[-1], v.names = "epsilon1", idvar = c("Knot"), timevar = "Date", times = names(eps1_wide)[-1]))
eps1_long$Date<- factor(eps1_long$Date, levels = all_times)
rownames(eps1_long)<- NULL

eps2_long<- data.frame(reshape(eps2_wide, direction = "long", varying = names(eps2_wide)[-1], v.names = "epsilon2", idvar = c("Knot"), timevar = "Date", times = names(eps2_wide)[-1]))
eps2_long$Date<- factor(eps2_long$Date, levels = all_times)
rownames(eps2_long)<- NULL

# Combine them
eps_out<- merge(eps1_long, eps2_long)

# Bring in other components
knot_preds_out_temp<- merge(eps_out, knot_preds_out, by = "Knot")
knot_preds_out<- merge(knot_preds_out_temp, betas, by = "Date")

# Arrange
knot_preds_out<- knot_preds_out[order(knot_preds_out$Date, knot_preds_out$Knot), ]

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
eta1<- t(as.vector(x1_covs) %*% t(mod_mat))
eta2<- t(as.vector(x2_covs) %*% t(mod_mat))

# Let's get this into a bit nicer format to keep track of time/space, too...
preds_out<- data.frame("Date" = factor(new_covariate_data$Year), "Lat" = new_covariate_data$Lat, "Lon" = new_covariate_data$Lon, "eta1" = eta1, "eta2" = eta2)

## Bringing everything together. For the random effects measured at knot locations, we need to snap these values to the prediction lat/long grid cells. 
# Find nearest knot to each observation
predict_nn = RANN::nn2(query = preds_out[,c('Lat','Lon')], data = fit$spatial_list$latlon_x[,c('Lat','Lon')], k=1)$nn.idx

# Bring over knot estimate...
preds_out<- cbind(preds_out, knot_preds_out[match(predict_nn, knot_preds_out$Knot), ])

## Make the calculations...
p1_test<- preds_out$omega1 + preds_out$beta1 + preds_out$epsilon1 + preds_out$eta1
p2_test<- preds_out$omega2 + preds_out$beta2 + preds_out$epsilon2 + preds_out$eta2
r1_test<- exp(P1)/(1+exp(P1))
r2_test<- exp(P2)
dens_test<- r1_test * r2_test

# Check
all.equal(p1_true, p1_test, tolerance = 0.001)
all.equal(dens_true, dens_test, tolerance = 0.001)
