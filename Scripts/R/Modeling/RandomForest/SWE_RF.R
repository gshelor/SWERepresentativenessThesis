##### SWE Random Forest Model for CONUS #####
### script by Griffin Shelor

##### Loading packages #####
library(pacman)
p_load(here, tidyverse, randomForest, ranger, ModelMetrics, gt, gtExtras, sf, terra, lime, sfext, rsample, caret, parallel, grateful, ModelMetrics, Metrics)
options(mc.cores = parallel::detectCores())

##### reading in data #####
### CONUS AOI
CONUS_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "CONUS", "CONUS_AOI.gpkg"))
### SNOTEL Climatologies used to fit the model
Snotel_Clim_sf <- read_sf(here("Data", "SNOTEL", "Combined", "GIS", "SnotelCombined_ClimatologyCovars.gpkg"))
### SNOTEL Annual values for testing
Snotel_sf <- read_sf(here("Data", "SNOTEL", "Combined", "GIS", "SnotelCombined_AnnualCovars.gpkg"))

### converting to dfs for model fitting process
Snotel_Clim_df <- sf_to_df(Snotel_Clim_sf) |>
  drop_na()
Snotel_df <- sf_to_df(Snotel_sf) |>
  drop_na()
### making sure landcover is type integer
Snotel_Clim_df$landcover <- as.integer(Snotel_Clim_df$landcover)
Snotel_Clim_df$landcover_triclass <- as.integer(Snotel_Clim_df$landcover_triclass)
Snotel_df$landcover <- as.integer(Snotel_df$landcover)
Snotel_df$landcover_triclass <- as.integer(Snotel_df$landcover_triclass)

##### splitting data into training and testing datasets #####
# set.seed(802)
# Snotel_split <- initial_split(Snotel_Clim_df, prop = 0.7)
# Snotel_train <- training(Snotel_split)
# Snotel_test <- testing(Snotel_split)
### keeping as dataframes for caret
### OctApr aggregates only
train_OctApr_x <- Snotel_Clim_df[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "OctApr_tminMean", "OctApr_tmaxMean", "OctApr_sradMean", "landcover_triclass")]
train_y <- Snotel_Clim_df$peak_swe
test_OctApr_x <- Snotel_df[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "OctApr_tminMean", "OctApr_tmaxMean", "OctApr_sradMean", "landcover_triclass")]
test_y <- Snotel_df$peak_swe
### OctApr aggregates only but with aspect
train_OctAprAspect_x <- Snotel_Clim_df[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tminMean", "OctApr_tmaxMean", "OctApr_sradMean", "landcover_triclass")]
test_OctAprAspect_x <- Snotel_df[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tminMean", "OctApr_tmaxMean", "OctApr_sradMean", "landcover_triclass")]
### OctApr aggregates only but with aspect, and srad removed
train_OctAprAspect_Nosrad_x <- Snotel_Clim_df[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tminMean", "OctApr_tmaxMean", "landcover_triclass")]
test_OctAprAspect_Nosrad_x <- Snotel_df[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tminMean", "OctApr_tmaxMean", "landcover_triclass")]
### OctApr aggregates only but with aspect, and srad, tmin removed
train_OctAprAspect_Nosradtmin_x <- Snotel_Clim_df[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tmaxMean", "landcover_triclass")]
test_OctAprAspect_Nosradtmin_x <- Snotel_df[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tmaxMean", "landcover_triclass")]
### OctApr aggregates only but with aspect, and srad, tmin, tmax removed
train_OctAprAspect_Nosradtmintmax_x <- Snotel_Clim_df[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "landcover_triclass")]
test_OctAprAspect_Nosradtmintmax_x <- Snotel_df[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "landcover_triclass")]
### OctApr aggregates only but with aspect, and no tmin
train_OctAprAspect_Notmin_x <- Snotel_Clim_df[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tmaxMean", "OctApr_sradMean", "landcover_triclass")]
test_OctAprAspect_Notmin_x <- Snotel_df[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tmaxMean", "OctApr_sradMean", "landcover_triclass")]
### OctApr aggregates only but with aspect, and no tmin, tmax
train_OctAprAspect_Notmintmax_x <- Snotel_Clim_df[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_sradMean", "landcover_triclass")]
test_OctAprAspect_Notmintmax_x <- Snotel_df[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_sradMean", "landcover_triclass")]
### OctApr aggregates only but with aspect, and no tmax
train_OctAprAspect_Notmax_x <- Snotel_Clim_df[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tminMean", "OctApr_sradMean", "landcover_triclass")]
test_OctAprAspect_Notmax_x <- Snotel_df[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tminMean", "OctApr_sradMean", "landcover_triclass")]
### OctApr aggregates and DecFeb aggregates for prcpSum and tmean
train_OctAprDecFeb_x <- Snotel_Clim_df[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "OctApr_tminMean", "OctApr_tmaxMean", "OctApr_sradMean", "landcover_triclass")]
test_OctAprDecFeb_x <- Snotel_df[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "OctApr_tminMean", "OctApr_tmaxMean", "OctApr_sradMean", "landcover_triclass")]
### Testing OctMay aggregations in model instead of OctApr
train_OctMay_x <- Snotel_Clim_df[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "elevation", "slope", "OctMay_tminMean", "OctMay_tmaxMean", "OctMay_sradMean", "landcover_triclass")]
test_OctMay_x <- Snotel_df[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "elevation", "slope", "OctMay_tminMean", "OctMay_tmaxMean", "OctMay_sradMean", "landcover_triclass")]
### OctMay aggregates only but with aspect, and srad, tmin removed
train_OctMayAspect_Nosradtmin_x <- Snotel_Clim_df[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "elevation", "slope", "aspect", "OctMay_tmaxMean", "landcover_triclass")]
test_OctMayAspect_Nosradtmin_x <- Snotel_df[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "elevation", "slope", "aspect", "OctMay_tmaxMean", "landcover_triclass")]
### OctMay aggregates only but with aspect, and srad, tmin, tmax removed
train_OctMayAspect_Nosradtmintmax_x <- Snotel_Clim_df[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "elevation", "slope", "aspect", "landcover_triclass")]
test_OctMayAspect_Nosradtmintmax_x <- Snotel_df[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "elevation", "slope", "aspect", "landcover_triclass")]
### OctMay aggregates only but with aspect, and no tmin, tmax
train_OctMayAspect_Notmintmax_x <- Snotel_Clim_df[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "elevation", "slope", "aspect", "OctMay_sradMean", "landcover_triclass")]
test_OctMayAspect_Notmintmax_x <- Snotel_df[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "elevation", "slope", "aspect", "OctMay_sradMean", "landcover_triclass")]
### OctMay aggregates and DecFeb aggregates for prcpSum and tmean
train_OctMayDecFeb_x <- Snotel_Clim_df[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "OctMay_tminMean", "OctMay_tmaxMean", "OctMay_sradMean", "landcover_triclass")]
test_OctMayDecFeb_x <- Snotel_df[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "OctMay_tminMean", "OctMay_tmaxMean", "OctMay_sradMean", "landcover_triclass")]
### Testing SepMay aggregations in model instead of OctApr
train_SepMay_x <- Snotel_Clim_df[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "elevation", "slope", "SepMay_tminMean", "SepMay_tmaxMean", "SepMay_sradMean", "landcover_triclass")]
test_SepMay_x <- Snotel_df[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "elevation", "slope", "SepMay_tminMean", "SepMay_tmaxMean", "SepMay_sradMean", "landcover_triclass")]
### SepMay aggregates and DecFeb aggregates for prcpSum and tmean
train_SepMayDecFeb_x <- Snotel_Clim_df[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "SepMay_tminMean", "SepMay_tmaxMean", "SepMay_sradMean", "landcover_triclass")]
test_SepMayDecFeb_x <- Snotel_df[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "SepMay_tminMean", "SepMay_tmaxMean", "SepMay_sradMean", "landcover_triclass")]
### SepMay aggregates only but with aspect, and srad, tmin removed
train_SepMayAspect_Nosradtmin_x <- Snotel_Clim_df[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "elevation", "slope", "aspect", "SepMay_tmaxMean", "landcover_triclass")]
test_SepMayAspect_Nosradtmin_x <- Snotel_df[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "elevation", "slope", "aspect", "SepMay_tmaxMean", "landcover_triclass")]
### SepMay aggregates only but with aspect, and srad, tmin, tmax removed
train_SepMayAspect_Nosradtmintmax_x <- Snotel_Clim_df[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "elevation", "slope", "aspect", "landcover_triclass")]
test_SepMayAspect_Nosradtmintmax_x <- Snotel_df[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "elevation", "slope", "aspect", "landcover_triclass")]
### SepMay aggregates only but with aspect, and no tmin, tmax
train_SepMayAspect_Notmintmax_x <- Snotel_Clim_df[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "elevation", "slope", "aspect", "SepMay_sradMean", "landcover_triclass")]
test_SepMayAspect_Notmintmax_x <- Snotel_df[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "elevation", "slope", "aspect", "SepMay_sradMean", "landcover_triclass")]

##### reading in rasters for spatial predictions #####
prcpSum_OctAprClimSum_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpSum_OctAprSum_ClimatologyCDM.tif"))
prcpMean_OctAprClimMean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpMean_OctAprMean_ClimatologyCDM.tif"))
tmean_OctAprClimSum_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmean_OctAprSum_ClimatologyCDM.tif"))
tmin_OctAprClimMean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmin_OctAprMean_ClimatologyCDM.tif"))
tmax_OctAprClimMean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmax_OctAprMean_ClimatologyCDM.tif"))
srad_OctAprClimMean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "srad_OctAprMean_ClimatologyCDM.tif"))
### OctMay Rasters
prcpSum_OctMayClimSum_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpSum_OctMaySum_ClimatologyCDM.tif"))
prcpMean_OctMayClimMean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpMean_OctMayMean_ClimatologyCDM.tif"))
tmean_OctMayClimSum_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmean_OctMaySum_ClimatologyCDM.tif"))
tmin_OctMayClimMean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmin_OctMayMean_ClimatologyCDM.tif"))
tmax_OctMayClimMean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmax_OctMayMean_ClimatologyCDM.tif"))
srad_OctMayClimMean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "srad_OctMayMean_ClimatologyCDM.tif"))
### SepMay Rasters
prcpSum_SepMayClimSum_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpSum_SepMaySum_ClimatologyCDM.tif"))
prcpMean_SepMayClimMean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpMean_SepMayMean_ClimatologyCDM.tif"))
tmean_SepMayClimSum_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmean_SepMaySum_ClimatologyCDM.tif"))
tmin_SepMayClimMean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmin_SepMayMean_ClimatologyCDM.tif"))
tmax_SepMayClimMean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmax_SepMayMean_ClimatologyCDM.tif"))
srad_SepMayClimMean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "srad_SepMayMean_ClimatologyCDM.tif"))

### CONUS topography
CONUS_DEM <- rast(here("Data", "DEM", "CONUS", "CONUSDEMMosaic.tif"))
CONUS_Slope <- rast(here("Data", "DEM", "CONUS", "CONUSSlope.tif"))
CONUS_Aspect <- rast(here("Data", "DEM", "CONUS", "CONUSAspect.tif"))
### CONUS Landcover
CONUS_LCClim_rast <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Climatology", "LC_CONUS_Climatology.tif"))


##### Fitting RF Model with Caret #####
### Hyperparameter Tuning Grid
tune_grid_OctApr_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(train_OctApr_x) - 1, by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 25, by = 1)
)

tune_grid_OctAprAspect_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(train_OctAprAspect_x) - 1, by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 25, by = 1)
)

tune_grid_OctAprAspect_Nosrad_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(train_OctAprAspect_Nosrad_x) - 1, by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 25, by = 1)
)

tune_grid_OctAprAspect_Nosradtmin_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(train_OctAprAspect_Nosradtmin_x) - 1, by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 25, by = 1)
)

##### Tune grid final model #####
tune_grid_OctAprAspect_Nosradtmintmax_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(train_OctAprAspect_Nosradtmintmax_x) - 1, by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 25, by = 1)
)

tune_grid_OctAprAspect_Notmin_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(train_OctAprAspect_Notmin_x) - 1, by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 25, by = 1)
)

tune_grid_OctAprAspect_Notmintmax_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(train_OctAprAspect_Notmintmax_x) - 1, by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 25, by = 1)
)

tune_grid_OctAprAspect_Notmax_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(train_OctAprAspect_Notmax_x) - 1, by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 25, by = 1)
)

tune_grid_OctMay_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(train_OctMay_x) - 1, by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 25, by = 1)
)

tune_grid_OctMayAspect_Nosradtmin_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(train_OctMayAspect_Nosradtmin_x) - 1, by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 25, by = 1)
)

tune_grid_OctMayAspect_Nosradtmintmax_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(train_OctMayAspect_Nosradtmintmax_x) - 1, by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 25, by = 1)
)

tune_grid_OctMayAspect_Notmintmax_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(train_OctMayAspect_Notmintmax_x) - 1, by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 25, by = 1)
)

tune_grid_SepMay_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(train_SepMay_x) - 1, by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 25, by = 1)
)

tune_grid_OctAprDecFeb_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(train_OctAprDecFeb_x) - 1, by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 25, by = 1)
)

tune_grid_OctMayDecFeb_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(train_OctMayDecFeb_x) - 1, by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 25, by = 1)
)

tune_grid_SepMayDecFeb_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(train_SepMayDecFeb_x) - 1, by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 25, by = 1)
)

tune_grid_SepMayAspect_Nosradtmin_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(train_SepMayAspect_Nosradtmin_x) - 1, by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 25, by = 1)
)

tune_grid_SepMayAspect_Nosradtmintmax_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(train_SepMayAspect_Nosradtmintmax_x) - 1, by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 25, by = 1)
)

tune_grid_SepMayAspect_Notmintmax_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(train_SepMayAspect_Notmintmax_x) - 1, by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 25, by = 1)
)

##### Training Control with Cross-Validation #####
tune_control_rf <- trainControl(
  method = "oob", # out of bag
  number = 20, # 10-fold cross-validation
  # repeats = 3, # Repeat cross-validation 3 times (more robust)
  verboseIter = TRUE, # Print progress during training
  returnData = FALSE, # Don't save the training data (saves memory)
  savePredictions = "none", # Save predictions from the best model
  returnResamp = "final", # Save resampling results from the best model
  allowParallel = TRUE # Enable parallel processing (if available)
)

tune_control_loocv_rf <- trainControl(
  method = "LOOCV", # leave one out CV
  number = 20, # 10-fold cross-validation
  # repeats = 3, # Repeat cross-validation 3 times (more robust)
  verboseIter = TRUE, # Print progress during training
  returnData = FALSE, # Don't save the training data (saves memory)
  savePredictions = "none", # Save predictions from the best model
  returnResamp = "final", # Save resampling results from the best model
  allowParallel = TRUE, # Enable parallel processing (if available)
  predictionBounds = c(0, NA)
)

tune_control_boot_rf <- trainControl(
  method = "boot_all", # out of bag
  number = 20, # 10-fold cross-validation
  # repeats = 3, # Repeat cross-validation 3 times (more robust)
  verboseIter = TRUE, # Print progress during training
  returnData = FALSE, # Don't save the training data (saves memory)
  savePredictions = "none", # Save predictions from the best model
  returnResamp = "final", # Save resampling results from the best model
  allowParallel = TRUE, # Enable parallel processing (if available)
  predictionBounds = c(0, NA)
)

tune_control_rf <- trainControl(
  method = "cv", # out of bag
  number = 20, # 10-fold cross-validation
  # repeats = 3, # Repeat cross-validation 3 times (more robust)
  verboseIter = TRUE, # Print progress during training
  returnData = FALSE, # Don't save the training data (saves memory)
  savePredictions = "none", # Save predictions from the best model
  returnResamp = "final", # Save resampling results from the best model
  allowParallel = TRUE, # Enable parallel processing (if available)
  predictionBounds = c(0, NA)
)

### Train the Random Forest Model
set.seed(802)
RF_model1 <- train(
  x = train_OctApr_x,
  y = train_y,
  method = "rf",
  trControl = tune_control_rf,
  tuneGrid = tune_grid_OctApr_rf,
  importance = TRUE,
  n.tree = 2500
)
set.seed(802)
RF_model2 <- train(
  x = train_OctMay_x,
  y = train_y,
  method = "rf",
  trControl = tune_control_rf,
  tuneGrid = tune_grid_OctMay_rf,
  importance = TRUE,
  n.tree = 2500
)
set.seed(802)
RF_model3 <- train(
  x = train_SepMay_x,
  y = train_y,
  method = "rf",
  trControl = tune_control_rf,
  tuneGrid = tune_grid_SepMay_rf,
  importance = TRUE,
  n.tree = 2500
)
### Model 4, OctApr aggregates and DecFeb aggregates for prcpSum and tmean
set.seed(802)
RF_model4 <- train(
  x = train_OctAprDecFeb_x,
  y = train_y,
  method = "rf",
  trControl = tune_control_rf,
  tuneGrid = tune_grid_OctAprDecFeb_rf,
  importance = TRUE,
  n.tree = 2500
)
### Model 5, OctMay aggregates and DecFeb aggregates for prcpSum and tmean
set.seed(802)
RF_model5 <- train(
  x = train_OctMayDecFeb_x,
  y = train_y,
  method = "rf",
  trControl = tune_control_rf,
  tuneGrid = tune_grid_OctMayDecFeb_rf,
  importance = TRUE,
  n.tree = 2500
)
### Model 6, SepMay aggregates and DecFeb aggregates for prcpSum and tmean
set.seed(802)
RF_model6 <- train(
  x = train_SepMayDecFeb_x,
  y = train_y,
  method = "rf",
  trControl = tune_control_rf,
  tuneGrid = tune_grid_SepMayDecFeb_rf,
  importance = TRUE,
  n.tree = 2500
)
### Model 7, OctApr aggregates and aspect
set.seed(802)
RF_model7 <- train(
  x = train_OctAprAspect_x,
  y = train_y,
  method = "rf",
  trControl = tune_control_rf,
  tuneGrid = tune_grid_OctAprAspect_rf,
  importance = TRUE,
  n.tree = 2500
)
### Model 8, OctApr aggregates and aspect, no srad
set.seed(802)
RF_model8 <- train(
  x = train_OctAprAspect_Nosrad_x,
  y = train_y,
  method = "rf",
  trControl = tune_control_rf,
  tuneGrid = tune_grid_OctAprAspect_Nosrad_rf,
  importance = TRUE,
  n.tree = 2500
)
### Model 9, OctApr aggregates and aspect, no srad or tmin
set.seed(802)
RF_model9 <- train(
  x = train_OctAprAspect_Nosradtmin_x,
  y = train_y,
  method = "rf",
  trControl = tune_control_rf,
  tuneGrid = tune_grid_OctAprAspect_Nosradtmin_rf,
  importance = TRUE,
  n.tree = 2500
)
##### Model 10, OctApr aggregates and aspect, no srad or tmin or tmax #####
set.seed(802)
RF_model10 <- train(
  x = train_OctAprAspect_Nosradtmintmax_x,
  y = train_y,
  method = "rf",
  trControl = tune_control_rf,
  tuneGrid = tune_grid_OctAprAspect_Nosradtmintmax_rf,
  importance = TRUE,
  n.tree = 2500
)
### Model 11, OctApr aggregates and aspect, no tmin
set.seed(802)
RF_model11 <- train(
  x = train_OctAprAspect_Notmin_x,
  y = train_y,
  method = "rf",
  trControl = tune_control_rf,
  tuneGrid = tune_grid_OctAprAspect_Notmin_rf,
  importance = TRUE,
  n.tree = 2500
)
### Model 12, OctApr aggregates and aspect, no tmin or tmax
set.seed(802)
RF_model12 <- train(
  x = train_OctAprAspect_Notmintmax_x,
  y = train_y,
  method = "rf",
  trControl = tune_control_rf,
  tuneGrid = tune_grid_OctAprAspect_Notmintmax_rf,
  importance = TRUE,
  n.tree = 2500
)
### Model 13, OctApr aggregates and aspect, no tmax
set.seed(802)
RF_model13 <- train(
  x = train_OctAprAspect_Notmax_x,
  y = train_y,
  method = "rf",
  trControl = tune_control_rf,
  tuneGrid = tune_grid_OctAprAspect_Notmax_rf,
  importance = TRUE,
  n.tree = 2500
)
### Model 14, OctMay aggregates and aspect, no srad, tmin
set.seed(802)
RF_model14 <- train(
  x = train_OctMayAspect_Nosradtmin_x,
  y = train_y,
  method = "rf",
  trControl = tune_control_rf,
  tuneGrid = tune_grid_OctMayAspect_Nosradtmin_rf,
  importance = TRUE,
  n.tree = 2500
)
### Model 15, OctMay aggregates and aspect, no srad, tmin, tmax
set.seed(802)
RF_model15 <- train(
  x = train_OctMayAspect_Nosradtmintmax_x,
  y = train_y,
  method = "rf",
  trControl = tune_control_rf,
  tuneGrid = tune_grid_OctMayAspect_Nosradtmintmax_rf,
  importance = TRUE,
  n.tree = 2500
)
### Model 16, OctMay aggregates and aspect, no tmin or tmax
set.seed(802)
RF_model16 <- train(
  x = train_OctMayAspect_Notmintmax_x,
  y = train_y,
  method = "rf",
  trControl = tune_control_rf,
  tuneGrid = tune_grid_OctMayAspect_Notmintmax_rf,
  importance = TRUE,
  n.tree = 2500
)
### Model 17, SepMay aggregates and aspect, no srad, tmin
set.seed(802)
RF_model17 <- train(
  x = train_SepMayAspect_Nosradtmin_x,
  y = train_y,
  method = "rf",
  trControl = tune_control_rf,
  tuneGrid = tune_grid_SepMayAspect_Nosradtmin_rf,
  importance = TRUE,
  n.tree = 2500
)
### Model 18, SepMay aggregates and aspect, no srad, tmin, tmax
set.seed(802)
RF_model18 <- train(
  x = train_SepMayAspect_Nosradtmintmax_x,
  y = train_y,
  method = "rf",
  trControl = tune_control_rf,
  tuneGrid = tune_grid_SepMayAspect_Nosradtmintmax_rf,
  importance = TRUE,
  n.tree = 2500
)
### Model 19, SepMay aggregates and aspect, no tmin or tmax
set.seed(802)
RF_model19 <- train(
  x = train_SepMayAspect_Notmintmax_x,
  y = train_y,
  method = "rf",
  trControl = tune_control_rf,
  tuneGrid = tune_grid_SepMayAspect_Notmintmax_rf,
  importance = TRUE,
  n.tree = 2500
)
##### saving model to avoid fitting it again #####
saveRDS(RF_model1, file = here("Data", "FittedModels", "Caret", "RF", "OctApr_RFModel_1km.rds"))
saveRDS(RF_model2, file = here("Data", "FittedModels", "Caret", "RF", "OctMay_RFModel_1km.rds"))
saveRDS(RF_model3, file = here("Data", "FittedModels", "Caret", "RF", "SepMay_RFModel_1km.rds"))
saveRDS(RF_model4, file = here("Data", "FittedModels", "Caret", "RF", "OctAprDecFeb_RFModel_1km.rds"))
saveRDS(RF_model5, file = here("Data", "FittedModels", "Caret", "RF", "OctMayDecFeb_RFModel_1km.rds"))
saveRDS(RF_model6, file = here("Data", "FittedModels", "Caret", "RF", "SepMayDecFeb_RFModel_1km.rds"))
saveRDS(RF_model7, file = here("Data", "FittedModels", "Caret", "RF", "OctAprAspect_RFModel_1km.rds"))
saveRDS(RF_model8, file = here("Data", "FittedModels", "Caret", "RF", "OctAprAspect_Nosrad_RFModel_1km.rds"))
saveRDS(RF_model9, file = here("Data", "FittedModels", "Caret", "RF", "OctAprAspect_Nosradtmin_RFModel_1km.rds"))
saveRDS(RF_model10, file = here("Data", "FittedModels", "Caret", "RF", "OctAprAspect_Nosradtmintmax_RFModel_1km.rds"))
saveRDS(RF_model11, file = here("Data", "FittedModels", "Caret", "RF", "OctAprAspect_Notmin_RFModel_1km.rds"))
saveRDS(RF_model12, file = here("Data", "FittedModels", "Caret", "RF", "OctAprAspect_Notmintmax_RFModel_1km.rds"))
saveRDS(RF_model13, file = here("Data", "FittedModels", "Caret", "RF", "OctAprAspect_Notmax_RFModel_1km.rds"))
saveRDS(RF_model14, file = here("Data", "FittedModels", "Caret", "RF", "OctMayAspect_Nosradtmin_RFModel_1km.rds"))
saveRDS(RF_model15, file = here("Data", "FittedModels", "Caret", "RF", "OctMayAspect_Nosradtmintmax_RFModel_1km.rds"))
saveRDS(RF_model16, file = here("Data", "FittedModels", "Caret", "RF", "OctMayAspect_Notmintmax_RFModel_1km.rds"))
saveRDS(RF_model17, file = here("Data", "FittedModels", "Caret", "RF", "SepMayAspect_Nosradtmin_RFModel_1km.rds"))
saveRDS(RF_model18, file = here("Data", "FittedModels", "Caret", "RF", "SepMayAspect_Nosradtmintmax_RFModel_1km.rds"))
saveRDS(RF_model19, file = here("Data", "FittedModels", "Caret", "RF", "SepMayAspect_Notmintmax_RFModel_1km.rds"))

##### Examine the Best Parameters and Results #####
### Model 1
print(RF_model1$bestTune)
print(RF_model1$results)
set.seed(802)
RF_model1_preds <- predict(RF_model1, test_OctApr_x)
rmse(test_y, RF_model1_preds)
### Model 2
print(RF_model2$bestTune)
print(RF_model2$results)
set.seed(802)
RF_model2_preds <- predict(RF_model2, test_OctMay_x)
rmse(test_y, RF_model2_preds)
### Model 3
print(RF_model3$bestTune)
print(RF_model3$results)
set.seed(802)
RF_model3_preds <- predict(RF_model3, test_SepMay_x)
rmse(test_y, RF_model3_preds)
### Model 4
print(RF_model4$bestTune)
print(RF_model4$results)
set.seed(802)
RF_model4_preds <- predict(RF_model4, test_OctAprDecFeb_x)
rmse(test_y, RF_model4_preds)
### Model 5
print(RF_model5$bestTune)
print(RF_model5$results)
set.seed(802)
RF_model5_preds <- predict(RF_model5, test_OctMayDecFeb_x)
rmse(test_y, RF_model5_preds)
### Model 6
print(RF_model6$bestTune)
print(RF_model6$results)
set.seed(802)
RF_model6_preds <- predict(RF_model6, test_SepMayDecFeb_x)
rmse(test_y, RF_model6_preds)
### Model 7
print(RF_model7$bestTune)
print(RF_model7$results)
set.seed(802)
RF_model7_preds <- predict(RF_model7, test_OctAprAspect_x)
rmse(test_y, RF_model7_preds)
### Model 8
print(RF_model8$bestTune)
print(RF_model8$results)
set.seed(802)
RF_model8_preds <- predict(RF_model8, test_OctAprAspect_Nosrad_x)
rmse(test_y, RF_model8_preds)
### Model 9
print(RF_model9$bestTune)
print(RF_model9$results)
set.seed(802)
RF_model9_preds <- predict(RF_model9, test_OctAprAspect_Nosradtmin_x)
rmse(test_y, RF_model9_preds)
### Model 10
##### Final Model #####
RF_model10 <- read_rds(here("Data", "FittedModels", "Caret", "RF", "OctAprAspect_Nosradtmintmax_RFModel_1km.rds"))
print(RF_model10$bestTune)
print(RF_model10$results)
set.seed(802)
RF_model10_preds <- predict(RF_model10, test_OctAprAspect_Nosradtmintmax_x)
set.seed(802)
Snotel_SWEpreds <- Snotel_df |>
  mutate(SWE_preds = predict(RF_model10, Snotel_df),
         SWE_ae = ae(peak_swe, SWE_preds)) |>
  select(site_id, site_name, state, WaterYear, peak_swe, SWE_preds, SWE_ae)
Snotel_SWEpreds_WYErrorMeans <- Snotel_SWEpreds |>
  dplyr::group_by(WaterYear) |>
  summarise(mean_SWE_ae = mean(SWE_ae),
            median_SWE_ae = median(SWE_ae))
plot(Snotel_SWEpreds$peak_swe, Snotel_SWEpreds$SWE_ae, main = "SNOTEL Peak SWE Value on X-axis, Absolute Error of Predictions on Y")
hist(Snotel_SWEpreds$SWE_ae, main = "Histogram of Absolute Error at All Annual SNOTEL Sites")
plot(Snotel_SWEpreds$WaterYear, Snotel_SWEpreds$SWE_ae, main = "Water Year on X-axis, Absolute Error of Predictions on Y")
plot(Snotel_SWEpreds_WYErrorMeans$WaterYear, Snotel_SWEpreds_WYErrorMeans$mean_SWE_ae, main = "Water Year on X-axis, Mean Absolute Error of Predictions on Y")
plot(Snotel_SWEpreds_WYErrorMeans$WaterYear, Snotel_SWEpreds_WYErrorMeans$median_SWE_ae, main = "Water Year on X-axis, Median Absolute Error of Predictions on Y")
rmse(test_y, RF_model10_preds)
### Model 11
print(RF_model11$bestTune)
print(RF_model11$results)
set.seed(802)
RF_model11_preds <- predict(RF_model11, test_OctAprAspect_Notmin_x)
rmse(test_y, RF_model11_preds)
### Model 12
print(RF_model12$bestTune)
print(RF_model12$results)
set.seed(802)
RF_model12_preds <- predict(RF_model12, test_OctAprAspect_Notmintmax_x)
rmse(test_y, RF_model12_preds)
### Model 13
print(RF_model13$bestTune)
print(RF_model13$results)
set.seed(802)
RF_model13_preds <- predict(RF_model13, test_OctAprAspect_Notmax_x)
rmse(test_y, RF_model13_preds)
### Model 14
print(RF_model14$bestTune)
print(RF_model14$results)
set.seed(802)
RF_model14_preds <- predict(RF_model14, test_OctMayAspect_Nosradtmin_x)
rmse(test_y, RF_model14_preds)
### Model 15
print(RF_model15$bestTune)
print(RF_model15$results)
set.seed(802)
RF_model15_preds <- predict(RF_model15, test_OctMayAspect_Nosradtmintmax_x)
rmse(test_y, RF_model15_preds)
### Model 16
print(RF_model16$bestTune)
print(RF_model16$results)
set.seed(802)
RF_model16_preds <- predict(RF_model16, test_OctMayAspect_Notmintmax_x)
rmse(test_y, RF_model16_preds)
### Model 17
print(RF_model17$bestTune)
print(RF_model17$results)
set.seed(802)
RF_model17_preds <- predict(RF_model17, test_SepMayAspect_Nosradtmin_x)
rmse(test_y, RF_model17_preds)
### Model 18
print(RF_model18$bestTune)
print(RF_model18$results)
set.seed(802)
RF_model18_preds <- predict(RF_model18, test_SepMayAspect_Nosradtmintmax_x)
rmse(test_y, RF_model18_preds)
### Model 19
print(RF_model19$bestTune)
print(RF_model19$results)
set.seed(802)
RF_model19_preds <- predict(RF_model19, test_SepMayAspect_Notmintmax_x)
rmse(test_y, RF_model19_preds)


