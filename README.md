# TargetsSDM
## Overview
This a repository that supports fitting the Vector Autoregressive Spatio-Temporal Model (Thorson and Barnett 2017)[https://doi.org/10.1093/icesjms/fsw193] (Thorson 2019)[https://doi.org/10.1016/j.fishres.2018.10.013] as a species distribution model to accomplish two objectives. The first objective is to create a vignette that walks through the modeling approach and slowly build model complexity with well documented code using example data. The second objective is to deploy a species distribution modeling workflow using the (Targets package)[https://github.com/ropensci/targets] in R. 

## The Vignette
The (Vignette)[[https://github.com/aallyn/TargetsSDM/blob/main/Vignette.Rmd] is a learning space. Within the Vignette, we start off with basic single species annual timestep model that estimates yearly intercepts as fixed effects and persistent spatial variability and ephemeral spatio-temporal variability. After fitting this basic model, we then demonstrate some functions that we have written to make inferences from the model (e.g., get parameter estimates, check model fit to the data, validate model to testing data) and visualize derived quantities (e.g., plotting statified biomass indices, mapping predicted density across space/time). From this basic model, we then slowly add complexity to the model with an aim towards ending with a single species, seasonal time step model that includes habitat covariates, catchability covariates, persistent spatial-variability, ephemeral spatio-temporal variability and temporal correlation in model intercepts as well as ephemeral spatio-temporal variability. 

Although we hope this provides useful "how to" information to potential users, the Vignette does not dive deeply into the statistical nitty-gritty that is going on behind the scenes. For that level of information, we'd encourage folks to read the many papers linked at the GitHub VAST page. In addition, we also provide our own (VAST for Beginners)[] document to accompany these other materials. This document comes with the disclaimer that none of us are statisticians. Although, we have spent considerable time trying to better understand the VAST modeling approach and more general spatio-temporal modeling concepts, we are sure that there are things we have not understood entirely. If you do come across something that stands out, please let us know as we are always trying to learn!

## The Targets workflow
![Workflow](SDMWorkflow.drawio.png “Targets Workflow”)
### Background
The Targets workflow of the repository more product focused than the Vignette. While we have tried to write functions and implement a workflow that is relatively flexibile, it is really tailored to a project that is designed to produce marine species distribution and abundance projections within the northwest Altantic ocean. To reach this goal, this Targets workflow

1.  integrates fisheries-independent bottom trawl data from the NOAA Northeast Fisheries Science Center spring/fall and Department of Fisheries and Oceans Canada spring/summer surveys
2.  fits a seasonal model VAST model that accounts for habitat covariates, catchability covariates, persistent spatial variability, ephemeral spatio-temporal variability, and potential temporal correlations in species occurrence across the study domain as well as between ephemeral spatio-temporal variability in successive seasons
3.  makes inferences from the fitted model, including vizualizing species-environment habitat relationships, evaluating model fit to the data and validating the predictive skill of the model using holdout, testing data
3.  uses the fitted model to project species distribution and abundance from 1982 to 2100 with potential environmental conditions estimated from an ensemble of statistically downscaled, biased corrected, CMIP6 global climate models

### Executing the workflow on your local machine
**WARNING: THIS TAKES A LONG TIME (DAYS) TO RUN AS IT IS CURRENTLY SET UP -- DAYS**
Although there are a number of functions within this repository that the workflow relies on, ultimately, executing the workflow involves interacting with two scripts. The first is the (_targets.R)[https://github.com/aallyn/TargetsSDM/blob/main/_targets.R] script. This is the playbook. At the top of the script, it has all of the library dependencies and then a few `source` calls to load required functions. Then, there is a code chunk ("Model fitting stuffs"). This sets up what play should be run -- the species, the years to use in model fitting, the years to project out to, which habitat model forumla to use, which catchability formula to use, etc. Following this chunk, there is the actual workflow that gets written as a list of connected "targets", where each target has a specific name and executes one command or function. The second script is the (targets_run.R)[https://github.com/aallyn/TargetsSDM/blob/main/targets_run.R] script. This script is really just a wrapper that calls the `_targets` file and executes the workflow by running the `tar_make()` function. 

To run everything as is on your local machine, you first want to fork and then clone this repo using your preferred Git/GitHub tool. When doing so, you would need to select where to clone the repo on your local machine. For this example, let's say you decided to clone the repo to your desktop, resulting in a new directory at `Desktop/TargetsSDM`. Next, you will want to make sure that you have installed DVC and they provide some (nice instructions)[https://dvc.org/doc/install] for how to do this depending on your OS. After installing DVC, you will then want to open up a new terminal and then change directories so that you go from the root directory to the directory of your cloned copy of this repo. On a Mac, this might look something like `cd ~/Desktop/TargetsSDM`. Once there, just to confirm you are in the right spot, type `ls` and you should see the list of files associated with this repo. Finally, once you have confirmed you are where you need to be, you can type `dvc pull` in the terminal window. This might take a moment as it is pulling down data from a remote google drive directory to your local machine. To make sure this step worked, type `ls ./data/covariates/static/` in your terminal window. You should see four different files: `Depth.grd`, `Depth.grd.dvc`, `Depth.gri`, `Depth.gri.dvc`. If you do, you are all set to continue. If you don't, feel free to [reach out to me](mailto:aallyn@gmri.org). 

Now that you have all of the code and data, you are ready to run the workflow. If you wanted to run everything "as is", you could open up the `targets_run.R` code, and run everything through line 29 (the `toc()` command following `tar_make()`). That said, it's unlikely that you will want to do that. Instead, you probably want to modify some component to the "Model fitting stuffs", which are in the `_targets.R` file. The most likely situation is you want to run a different species. To do that, you would just need to change the `nmfs_species_code` and `nice_category_names` to match. After making these edits, you would save the file and then return to the `targets_run.R` code and run everything through line 29. Other relatively easy changes include:

- Changing `fit_year_min` or `fit_year_max` to adjust the training data period
- Changing `pred_years` to adjust the projection time horizon
- Changing `hab_formula` to alter which covariates are included, with the caveat that if you remove `Season` or `Year_Cov` OR if you change `gam_degree` there's more editing to do
- Changing `field_config` or `rho_config` to adjust some of the VAST settings 
- Changing `strata_use` strata to use to get statified biomass indices 

Editing any of these objects is allowed and facilitated in the same way as the change to run a new species -- make the edits in `_targets.R`, save the file, and then run through line 29 of `targets_run.R`. There is some more flexibility, too (e.g., running the workflow and using an annual model rather than seasonal). That said, accounting for all these different options was making the code a bit messy. If you are interested in some modification beyond the changes listed above, please [email me](mailto:aallyn@gmri.org) and I'll do what I can to help. 

### Executing the workflow remotely
Coming soon...


