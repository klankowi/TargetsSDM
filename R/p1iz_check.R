#####
## Changing "_iz" predicted values depending on fine_scale TRUE/FALSE
#####

# Load packages
library(VAST)
library(splines)

# load data set
# see `?load_example` for list of stocks with example data
# that are installed automatically with `FishStatsUtils`.
example = load_example( data_set="covariate_example" )

# Test settings
field_config_use<- c("Omega1" = 1, "Epsilon1" = 0, "Omega2" = 1, "Epsilon2" = 0) # Modify here to narrow down where issue is occurring...

# Log and rescale BOT_DEPTH now...
example$covariate_data[,'BOT_DEPTH']<- log(example$covariate_data[,'BOT_DEPTH'])
# Rescale covariates being used to have an SD >0.1 and <10 (for numerical stability)
example$covariate_data[,'BOT_DEPTH'] = example$covariate_data[,'BOT_DEPTH'] / 100

year_set = min(example$covariate_data[,'Year']):max(example$covariate_data[,'Year'])
example$covariate_data[,'Year'] = NA

# NOW necessary to duplicate across years
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


# Make settings (turning off bias.correct to save time for example). 
settings_true = settings_false = make_settings( n_x=200,
                          ObsModel = c(1,0),
                          Region=example$Region,
                          purpose="index2",
                          use_anisotropy=FALSE,
                          bias.correct=FALSE,
                          fine_scale=FALSE,
                          FieldConfig = field_config_use)
# Adjust settings_true fine_scale
settings_true['fine_scale']<- TRUE

# Define formula -- a simple one for now
X1_formula = ~ bs( BOT_DEPTH, degree=2, intercept=FALSE)
X2_formula = X1_formula

# Fit the models
# Run model
fit_true = fit_model( "settings" = settings_true,
                 Lat_i = example$sampling_data[,'Lat'],
                 Lon_i = example$sampling_data[,'Lon'],
                 t_i = example$sampling_data[,'Year'],
                 b_i = example$sampling_data[,'Catch_KG'],
                 a_i = example$sampling_data[,'AreaSwept_km2'],
                 X1_formula = X1_formula,
                 X2_formula = X2_formula,
                 covariate_data = example$covariate_data )
fit_false = fit_model( "settings" = settings_false,
                       Lat_i = example$sampling_data[,'Lat'],
                       Lon_i = example$sampling_data[,'Lon'],
                       t_i = example$sampling_data[,'Year'],
                       b_i = example$sampling_data[,'Catch_KG'],
                       a_i = example$sampling_data[,'AreaSwept_km2'],
                       X1_formula = X1_formula,
                       X2_formula = X2_formula,
                       covariate_data = example$covariate_data )

# Predictions using `predict.fit_model`
p1_true = predict( x=fit_true,
                   Lat_i = fit_true$data_frame$Lat_i,
                   Lon_i = fit_true$data_frame$Lon_i,
                   t_i = fit_true$data_frame$t_i,
                   a_i = fit_true$data_frame$a_i,
                   what = "P1_iz",
                   new_covariate_data = example$covariate_data,
                   do_checks = FALSE )
p1_false = predict( x=fit_false,
                    Lat_i = fit_false$data_frame$Lat_i,
                    Lon_i = fit_false$data_frame$Lon_i,
                    t_i = fit_false$data_frame$t_i,
                    a_i = fit_false$data_frame$a_i,
                    what = "P1_iz",
                    new_covariate_data = example$covariate_data,
                    do_checks = FALSE )

# Test
all.equal(p1_true, p1_false, tolerance = 0.01)

#####
## What is going on with Ais and Ags?
#####

str(fit_true$spatial_list$A_is) # Weights relating knots to the observations
str(fit_false$spatial_list$A_is) # No weights, just a match of nearest neighbors


summary(fit_true$Report$Omega1_sc)
summary(fit_false$Report$Omega1_sc)
