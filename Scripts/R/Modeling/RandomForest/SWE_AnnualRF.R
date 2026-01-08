##### SWE Random Forest Model for CONUS #####
### script by Griffin Shelor

##### Loading packages #####
library(pacman)
p_load(here, tidyverse, randomForest, sf, terra, sfext, rsample, caret, parallel, grateful, ModelMetrics, Metrics, mcprogress)
# options(mc.cores = parallel::detectCores())

##### reading in data #####
### CONUS AOI
CONUS_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "CONUS", "CONUS_AOI.gpkg"))
### Alaska AOI
# AK_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "Alaska", "AK_AOI.gpkg"))
### Combined SNOTEL Annual values for model fitting and evaluation
# Snotel_sf <- read_sf(here("Data", "SNOTEL", "Combined", "GIS", "SnotelCombined_AnnualCovars_AnnualZeroPts.gpkg"))

### Reading in CONUS
SnotelCONUS_sf <- read_sf(here("Data", "SNOTEL", "CONUS", "GIS", "GPKG", "Annual", "SnotelCONUSAnnualZeroPts_PeakSWE_Covars.gpkg")) |>
  drop_na(OctApr_prcpSumCDMSum, OctApr_tmeanCDMSum, elevation, slope, aspect, OctApr_tminMean, OctApr_tmaxMean, OctApr_sradMean, landcover_triclass, OctMay_prcpSumCDMSum, OctMay_tmeanCDMSum, OctMay_tminMean, OctMay_tmaxMean, OctMay_sradMean, SepMay_prcpSumCDMSum, SepMay_tmeanCDMSum, SepMay_tminMean, SepMay_tmaxMean, SepMay_sradMean, DecFeb_prcpSumCDMSum, DecFeb_tmeanCDMSum, DecFeb_tminMean, DecFeb_tmaxMean, DecFeb_sradMean)

### Extracting just Alaska data
# SnotelAK_sf <- read_sf(here("Data", "SNOTEL", "Alaska", "GIS", "GPKG", "Annual", "SnotelAK_PeakSWE_Covars.gpkg")) |>
#   drop_na(OctApr_prcpSumCDMSum, OctApr_tmeanCDMSum, elevation, slope, aspect, OctApr_tminMean, OctApr_tmaxMean, OctApr_sradMean, landcover_triclass, OctMay_prcpSumCDMSum, OctMay_tmeanCDMSum, OctMay_tminMean, OctMay_tmaxMean, OctMay_sradMean, SepMay_prcpSumCDMSum, SepMay_tmeanCDMSum, SepMay_tminMean, SepMay_tmaxMean, SepMay_sradMean, DecFeb_prcpSumCDMSum, DecFeb_tmeanCDMSum, DecFeb_tminMean, DecFeb_tmaxMean, DecFeb_sradMean)

### converting to df for model fitting process
# Snotel_df <- sf_to_df(Snotel_sf) |>
#   drop_na(OctApr_prcpSumCDMSum, OctApr_tmeanCDMSum, elevation, slope, aspect, OctApr_tminMean, OctApr_tmaxMean, OctApr_sradMean, landcover_triclass, OctMay_prcpSumCDMSum, OctMay_tmeanCDMSum, OctMay_tminMean, OctMay_tmaxMean, OctMay_sradMean, SepMay_prcpSumCDMSum, SepMay_tmeanCDMSum, SepMay_tminMean, SepMay_tmaxMean, SepMay_sradMean)
SnotelCONUS_df <- sf_to_df(SnotelCONUS_sf)
# SnotelAK_df <- sf_to_df(SnotelAK_sf)
### making sure landcover is type integer
SnotelCONUS_df$landcover <- as.integer(SnotelCONUS_df$landcover)
SnotelCONUS_df$landcover_triclass <- as.integer(SnotelCONUS_df$landcover_triclass)
# SnotelAK_df$landcover <- as.integer(SnotelAK_df$landcover)
# SnotelAK_df$landcover_triclass <- as.integer(SnotelAK_df$landcover_triclass)
# 
# ##### splitting data into training and testing datasets #####
# ### splitting Alaska
# set.seed(802)
# SnotelAK_split <- initial_split(SnotelAK_df, prop = 0.75, strata = peak_swe)
# SnotelAK_train <- training(SnotelAK_split)
# SnotelAK_test <- testing(SnotelAK_split)
### splitting CONUS
set.seed(802)
SnotelCONUS_split <- initial_split(SnotelCONUS_df, prop = 0.75, strata = peak_swe)
SnotelCONUS_train <- training(SnotelCONUS_split)
SnotelCONUS_test <- testing(SnotelCONUS_split)
# ##### Alaska training and testing DFs #####
# ### keeping as dataframes for caret
# ### OctApr aggregates only
# AKTrain_OctApr_x <- SnotelAK_train[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "OctApr_tminMean", "OctApr_tmaxMean", "OctApr_sradMean", "landcover_triclass")]
# AKTrain_y <- SnotelAK_train$peak_swe
# AKTest_OctApr_x <- SnotelAK_test[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "OctApr_tminMean", "OctApr_tmaxMean", "OctApr_sradMean", "landcover_triclass")]
# AKTest_y <- SnotelAK_test$peak_swe
# ### OctApr aggregates only but with aspect
# AKTrain_OctAprAspect_x <- SnotelAK_train[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tminMean", "OctApr_tmaxMean", "OctApr_sradMean", "landcover_triclass")]
# AKTest_OctAprAspect_x <- SnotelAK_test[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tminMean", "OctApr_tmaxMean", "OctApr_sradMean", "landcover_triclass")]
# ### OctApr aggregates only but with aspect, and srad removed
# AKTrain_OctAprAspect_Nosrad_x <- SnotelAK_train[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tminMean", "OctApr_tmaxMean", "landcover_triclass")]
# AKTest_OctAprAspect_Nosrad_x <- SnotelAK_test[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tminMean", "OctApr_tmaxMean", "landcover_triclass")]
# ### OctApr aggregates only but with aspect, and srad, tmin removed
# AKTrain_OctAprAspect_Nosradtmin_x <- SnotelAK_train[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tmaxMean", "landcover_triclass")]
# AKTest_OctAprAspect_Nosradtmin_x <- SnotelAK_test[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tmaxMean", "landcover_triclass")]
# ### OctApr aggregates only but with aspect, and srad, tmin, tmax removed
# AKTrain_OctAprAspect_Nosradtmintmax_x <- SnotelAK_train[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "landcover_triclass")]
# AKTest_OctAprAspect_Nosradtmintmax_x <- SnotelAK_test[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "landcover_triclass")]
# ### OctApr aggregates only but with aspect, and no tmin
# AKTrain_OctAprAspect_Notmin_x <- SnotelAK_train[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tmaxMean", "OctApr_sradMean", "landcover_triclass")]
# AKTest_OctAprAspect_Notmin_x <- SnotelAK_test[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tmaxMean", "OctApr_sradMean", "landcover_triclass")]
# ### OctApr aggregates only but with aspect, and no tmin, tmax
# AKTrain_OctAprAspect_Notmintmax_x <- SnotelAK_train[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_sradMean", "landcover_triclass")]
# AKTest_OctAprAspect_Notmintmax_x <- SnotelAK_test[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_sradMean", "landcover_triclass")]
# ### OctApr aggregates only but with aspect, and no tmax
# AKTrain_OctAprAspect_Notmax_x <- SnotelAK_train[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tminMean", "OctApr_sradMean", "landcover_triclass")]
# AKTest_OctAprAspect_Notmax_x <- SnotelAK_test[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tminMean", "OctApr_sradMean", "landcover_triclass")]
# ### OctApr aggregates and DecFeb aggregates for prcpSum and tmean
# AKTrain_OctAprDecFeb_x <- SnotelAK_train[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "OctApr_tminMean", "OctApr_tmaxMean", "OctApr_sradMean", "landcover_triclass")]
# AKTest_OctAprDecFeb_x <- SnotelAK_test[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "OctApr_tminMean", "OctApr_tmaxMean", "OctApr_sradMean", "landcover_triclass")]
# ### Testing OctMay aggregations in model instead of OctApr
# AKTrain_OctMay_x <- SnotelAK_train[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "elevation", "slope", "OctMay_tminMean", "OctMay_tmaxMean", "OctMay_sradMean", "landcover_triclass")]
# AKTest_OctMay_x <- SnotelAK_test[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "elevation", "slope", "OctMay_tminMean", "OctMay_tmaxMean", "OctMay_sradMean", "landcover_triclass")]
# ### OctMay aggregates only but with aspect, and srad, tmin removed
# AKTrain_OctMayAspect_Nosradtmin_x <- SnotelAK_train[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "elevation", "slope", "aspect", "OctMay_tmaxMean", "landcover_triclass")]
# AKTest_OctMayAspect_Nosradtmin_x <- SnotelAK_test[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "elevation", "slope", "aspect", "OctMay_tmaxMean", "landcover_triclass")]
# ### OctMay aggregates only but with aspect, and srad, tmin, tmax removed
# AKTrain_OctMayAspect_Nosradtmintmax_x <- SnotelAK_train[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "elevation", "slope", "aspect", "landcover_triclass")]
# AKTest_OctMayAspect_Nosradtmintmax_x <- SnotelAK_test[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "elevation", "slope", "aspect", "landcover_triclass")]
# ### OctMay aggregates only but with aspect, and no tmin, tmax
# AKTrain_OctMayAspect_Notmintmax_x <- SnotelAK_train[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "elevation", "slope", "aspect", "OctMay_sradMean", "landcover_triclass")]
# AKTest_OctMayAspect_Notmintmax_x <- SnotelAK_test[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "elevation", "slope", "aspect", "OctMay_sradMean", "landcover_triclass")]
# ### OctMay aggregates and DecFeb aggregates for prcpSum and tmean
# AKTrain_OctMayDecFeb_x <- SnotelAK_train[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "OctMay_tminMean", "OctMay_tmaxMean", "OctMay_sradMean", "landcover_triclass")]
# AKTest_OctMayDecFeb_x <- SnotelAK_test[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "OctMay_tminMean", "OctMay_tmaxMean", "OctMay_sradMean", "landcover_triclass")]
# ### Testing SepMay aggregations in model instead of OctApr
# AKTrain_SepMay_x <- SnotelAK_train[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "elevation", "slope", "SepMay_tminMean", "SepMay_tmaxMean", "SepMay_sradMean", "landcover_triclass")]
# AKTest_SepMay_x <- SnotelAK_test[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "elevation", "slope", "SepMay_tminMean", "SepMay_tmaxMean", "SepMay_sradMean", "landcover_triclass")]
# ### SepMay aggregates and DecFeb aggregates for prcpSum and tmean
# AKTrain_SepMayDecFeb_x <- SnotelAK_train[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "SepMay_tminMean", "SepMay_tmaxMean", "SepMay_sradMean", "landcover_triclass")]
# AKTest_SepMayDecFeb_x <- SnotelAK_test[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "SepMay_tminMean", "SepMay_tmaxMean", "SepMay_sradMean", "landcover_triclass")]
# ### SepMay aggregates only but with aspect, and srad, tmin removed
# AKTrain_SepMayAspect_Nosradtmin_x <- SnotelAK_train[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "elevation", "slope", "aspect", "SepMay_tmaxMean", "landcover_triclass")]
# AKTest_SepMayAspect_Nosradtmin_x <- SnotelAK_test[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "elevation", "slope", "aspect", "SepMay_tmaxMean", "landcover_triclass")]
# ### SepMay aggregates only but with aspect, and srad, tmin, tmax removed
# AKTrain_SepMayAspect_Nosradtmintmax_x <- SnotelAK_train[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "elevation", "slope", "aspect", "landcover_triclass")]
# AKTest_SepMayAspect_Nosradtmintmax_x <- SnotelAK_test[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "elevation", "slope", "aspect", "landcover_triclass")]
# ### SepMay aggregates only but with aspect, and no tmin, tmax
# AKTrain_SepMayAspect_Notmintmax_x <- SnotelAK_train[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "elevation", "slope", "aspect", "SepMay_sradMean", "landcover_triclass")]
# AKTest_SepMayAspect_Notmintmax_x <- SnotelAK_test[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "elevation", "slope", "aspect", "SepMay_sradMean", "landcover_triclass")]
# ### OctMay aggregates and DecFeb aggregates for prcpSum and tmean
# AKTrain_OctMayDecFebAspect_x <- SnotelAK_train[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "aspect", "OctMay_tminMean", "OctMay_tmaxMean", "OctMay_sradMean", "landcover_triclass")]
# AKTest_OctMayDecFebAspect_x <- SnotelAK_test[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "aspect", "OctMay_tminMean", "OctMay_tmaxMean", "OctMay_sradMean", "landcover_triclass")]


##### CONUS training and testing DFs #####
### keeping as dataframes for caret
### OctApr aggregates only
CONUSTrain_OctApr_x <- SnotelCONUS_train[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "OctApr_tminMean", "OctApr_tmaxMean", "OctApr_sradMean", "landcover_triclass")]
CONUSTrain_y <- SnotelCONUS_train$peak_swe
CONUSTest_OctApr_x <- SnotelCONUS_test[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "OctApr_tminMean", "OctApr_tmaxMean", "OctApr_sradMean", "landcover_triclass")]
CONUSTest_y <- SnotelCONUS_test$peak_swe
### OctApr aggregates only but with aspect
CONUSTrain_OctAprAspect_x <- SnotelCONUS_train[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tminMean", "OctApr_tmaxMean", "OctApr_sradMean", "landcover_triclass")]
CONUSTest_OctAprAspect_x <- SnotelCONUS_test[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tminMean", "OctApr_tmaxMean", "OctApr_sradMean", "landcover_triclass")]
### OctApr aggregates only but with aspect, and srad removed
CONUSTrain_OctAprAspect_Nosrad_x <- SnotelCONUS_train[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tminMean", "OctApr_tmaxMean", "landcover_triclass")]
CONUSTest_OctAprAspect_Nosrad_x <- SnotelCONUS_test[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tminMean", "OctApr_tmaxMean", "landcover_triclass")]
### OctApr aggregates only but with aspect, and srad, tmin removed
CONUSTrain_OctAprAspect_Nosradtmin_x <- SnotelCONUS_train[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tmaxMean", "landcover_triclass")]
CONUSTest_OctAprAspect_Nosradtmin_x <- SnotelCONUS_test[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tmaxMean", "landcover_triclass")]
### OctApr aggregates only but with aspect, and srad, tmin, tmax removed
CONUSTrain_OctAprAspect_Nosradtmintmax_x <- SnotelCONUS_train[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "landcover_triclass")]
CONUSTest_OctAprAspect_Nosradtmintmax_x <- SnotelCONUS_test[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "landcover_triclass")]
### OctApr aggregates only but with aspect, and no tmin
CONUSTrain_OctAprAspect_Notmin_x <- SnotelCONUS_train[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tmaxMean", "OctApr_sradMean", "landcover_triclass")]
CONUSTest_OctAprAspect_Notmin_x <- SnotelCONUS_test[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tmaxMean", "OctApr_sradMean", "landcover_triclass")]
### OctApr aggregates only but with aspect, and no tmin, tmax
CONUSTrain_OctAprAspect_Notmintmax_x <- SnotelCONUS_train[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_sradMean", "landcover_triclass")]
CONUSTest_OctAprAspect_Notmintmax_x <- SnotelCONUS_test[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_sradMean", "landcover_triclass")]
### OctApr aggregates only but with aspect, and no tmax
CONUSTrain_OctAprAspect_Notmax_x <- SnotelCONUS_train[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tminMean", "OctApr_sradMean", "landcover_triclass")]
CONUSTest_OctAprAspect_Notmax_x <- SnotelCONUS_test[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tminMean", "OctApr_sradMean", "landcover_triclass")]
### OctApr aggregates and DecFeb aggregates for prcpSum and tmean
CONUSTrain_OctAprDecFeb_x <- SnotelCONUS_train[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "OctApr_tminMean", "OctApr_tmaxMean", "OctApr_sradMean", "landcover_triclass")]
CONUSTest_OctAprDecFeb_x <- SnotelCONUS_test[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "OctApr_tminMean", "OctApr_tmaxMean", "OctApr_sradMean", "landcover_triclass")]
### Testing OctMay aggregations in model instead of OctApr
CONUSTrain_OctMay_x <- SnotelCONUS_train[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "elevation", "slope", "OctMay_tminMean", "OctMay_tmaxMean", "OctMay_sradMean", "landcover_triclass")]
CONUSTest_OctMay_x <- SnotelCONUS_test[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "elevation", "slope", "OctMay_tminMean", "OctMay_tmaxMean", "OctMay_sradMean", "landcover_triclass")]
### OctMay aggregates only but with aspect, and srad, tmin removed
CONUSTrain_OctMayAspect_Nosradtmin_x <- SnotelCONUS_train[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "elevation", "slope", "aspect", "OctMay_tmaxMean", "landcover_triclass")]
CONUSTest_OctMayAspect_Nosradtmin_x <- SnotelCONUS_test[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "elevation", "slope", "aspect", "OctMay_tmaxMean", "landcover_triclass")]
### OctMay aggregates only but with aspect, and srad, tmin, tmax removed
CONUSTrain_OctMayAspect_Nosradtmintmax_x <- SnotelCONUS_train[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "elevation", "slope", "aspect", "landcover_triclass")]
CONUSTest_OctMayAspect_Nosradtmintmax_x <- SnotelCONUS_test[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "elevation", "slope", "aspect", "landcover_triclass")]
### OctMay aggregates only but with aspect, and no tmin, tmax
CONUSTrain_OctMayAspect_Notmintmax_x <- SnotelCONUS_train[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "elevation", "slope", "aspect", "OctMay_sradMean", "landcover_triclass")]
CONUSTest_OctMayAspect_Notmintmax_x <- SnotelCONUS_test[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "elevation", "slope", "aspect", "OctMay_sradMean", "landcover_triclass")]
### OctMay aggregates and DecFeb aggregates for prcpSum and tmean
CONUSTrain_OctMayDecFeb_x <- SnotelCONUS_train[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "OctMay_tminMean", "OctMay_tmaxMean", "OctMay_sradMean", "landcover_triclass")]
CONUSTest_OctMayDecFeb_x <- SnotelCONUS_test[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "OctMay_tminMean", "OctMay_tmaxMean", "OctMay_sradMean", "landcover_triclass")]
### Testing SepMay aggregations in model instead of OctApr
CONUSTrain_SepMay_x <- SnotelCONUS_train[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "elevation", "slope", "SepMay_tminMean", "SepMay_tmaxMean", "SepMay_sradMean", "landcover_triclass")]
CONUSTest_SepMay_x <- SnotelCONUS_test[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "elevation", "slope", "SepMay_tminMean", "SepMay_tmaxMean", "SepMay_sradMean", "landcover_triclass")]
### SepMay aggregates and DecFeb aggregates for prcpSum and tmean
CONUSTrain_SepMayDecFeb_x <- SnotelCONUS_train[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "SepMay_tminMean", "SepMay_tmaxMean", "SepMay_sradMean", "landcover_triclass")]
CONUSTest_SepMayDecFeb_x <- SnotelCONUS_test[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "SepMay_tminMean", "SepMay_tmaxMean", "SepMay_sradMean", "landcover_triclass")]
### SepMay aggregates only but with aspect, and srad, tmin removed
CONUSTrain_SepMayAspect_Nosradtmin_x <- SnotelCONUS_train[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "elevation", "slope", "aspect", "SepMay_tmaxMean", "landcover_triclass")]
CONUSTest_SepMayAspect_Nosradtmin_x <- SnotelCONUS_test[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "elevation", "slope", "aspect", "SepMay_tmaxMean", "landcover_triclass")]
### SepMay aggregates only but with aspect, and srad, tmin, tmax removed
CONUSTrain_SepMayAspect_Nosradtmintmax_x <- SnotelCONUS_train[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "elevation", "slope", "aspect", "landcover_triclass")]
CONUSTest_SepMayAspect_Nosradtmintmax_x <- SnotelCONUS_test[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "elevation", "slope", "aspect", "landcover_triclass")]
### SepMay aggregates only but with aspect, and no tmin, tmax
CONUSTrain_SepMayAspect_Notmintmax_x <- SnotelCONUS_train[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "elevation", "slope", "aspect", "SepMay_sradMean", "landcover_triclass")]
CONUSTest_SepMayAspect_Notmintmax_x <- SnotelCONUS_test[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "elevation", "slope", "aspect", "SepMay_sradMean", "landcover_triclass")]
### OctMay aggregates and DecFeb aggregates for prcpSum and tmean
CONUSTrain_OctMayDecFebAspect_x <- SnotelCONUS_train[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "aspect", "OctMay_tminMean", "OctMay_tmaxMean", "OctMay_sradMean", "landcover_triclass")]
CONUSTest_OctMayDecFebAspect_x <- SnotelCONUS_test[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "aspect", "OctMay_tminMean", "OctMay_tmaxMean", "OctMay_sradMean", "landcover_triclass")]
### OctMay aggregates and DecFeb aggregates for prcpSum and tmean, srad removed
CONUSTrain_OctMayDecFebAspect_Nosrad_x <- SnotelCONUS_train[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "aspect", "OctMay_tminMean", "OctMay_tmaxMean", "landcover_triclass")]
CONUSTest_OctMayDecFebAspect_Nosrad_x <- SnotelCONUS_test[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "aspect", "OctMay_tminMean", "OctMay_tmaxMean", "landcover_triclass")]



##### Hyperparameter Tuning Grid #####
tune_grid_OctApr_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(CONUSTrain_OctApr_x), by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 15, by = 1)
)

tune_grid_OctAprAspect_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(CONUSTrain_OctAprAspect_x), by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 15, by = 1)
)

tune_grid_OctAprAspect_Nosrad_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(CONUSTrain_OctAprAspect_Nosrad_x), by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 15, by = 1)
)

tune_grid_OctAprAspect_Nosradtmin_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(CONUSTrain_OctAprAspect_Nosradtmin_x), by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 15, by = 1)
)

##### Tune grid final model #####
tune_grid_OctAprAspect_Nosradtmintmax_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(CONUSTrain_OctAprAspect_Nosradtmintmax_x), by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 15, by = 1)
)

tune_grid_OctAprAspect_Notmin_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(CONUSTrain_OctAprAspect_Notmin_x), by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 15, by = 1)
)

tune_grid_OctAprAspect_Notmintmax_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(CONUSTrain_OctAprAspect_Notmintmax_x), by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 15, by = 1)
)

tune_grid_OctAprAspect_Notmax_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(CONUSTrain_OctAprAspect_Notmax_x), by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 15, by = 1)
)

tune_grid_OctMay_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(CONUSTrain_OctMay_x), by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 15, by = 1)
)

tune_grid_OctMayAspect_Nosradtmin_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(CONUSTrain_OctMayAspect_Nosradtmin_x), by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 15, by = 1)
)

tune_grid_OctMayAspect_Nosradtmintmax_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(CONUSTrain_OctMayAspect_Nosradtmintmax_x), by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 15, by = 1)
)

tune_grid_OctMayAspect_Notmintmax_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(CONUSTrain_OctMayAspect_Notmintmax_x), by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 15, by = 1)
)

tune_grid_SepMay_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(CONUSTrain_SepMay_x), by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 15, by = 1)
)

tune_grid_OctAprDecFeb_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(CONUSTrain_OctAprDecFeb_x), by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 15, by = 1)
)

tune_grid_OctMayDecFeb_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(CONUSTrain_OctMayDecFeb_x), by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 15, by = 1)
)

tune_grid_SepMayDecFeb_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(CONUSTrain_SepMayDecFeb_x), by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 15, by = 1)
)

tune_grid_SepMayAspect_Nosradtmin_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(CONUSTrain_SepMayAspect_Nosradtmin_x), by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 15, by = 1)
)

tune_grid_SepMayAspect_Nosradtmintmax_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(CONUSTrain_SepMayAspect_Nosradtmintmax_x), by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 15, by = 1)
)

tune_grid_SepMayAspect_Notmintmax_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(CONUSTrain_SepMayAspect_Notmintmax_x), by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 15, by = 1)
)

tune_grid_OctMayDecFebAspect_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(CONUSTrain_OctMayDecFebAspect_x), by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 15, by = 1)
)

tune_grid_OctMayDecFebAspect_Nosrad_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(CONUSTrain_OctMayDecFebAspect_Nosrad_x), by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 15, by = 1)
)

##### specifying train control objects #####
tune_control_oob_rf <- trainControl(
  method = "oob", # out of bag
  number = 10, # 10-fold cross-validation
  # repeats = 3, # Repeat cross-validation 3 times (more robust)
  verboseIter = TRUE, # Print progress during training
  returnData = FALSE, # Don't save the training data (saves memory)
  savePredictions = "none", # Save predictions from the best model
  returnResamp = "final", # Save resampling results from the best model
  allowParallel = TRUE, # Enable parallel processing (if available)
  predictionBounds = c(0, NA)
)

# tune_control_loocv_rf <- trainControl(
#   method = "LOOCV", # leave one out CV
#   # number = 10, # 10-fold cross-validation
#   # repeats = 3, # Repeat cross-validation 3 times (more robust)
#   verboseIter = TRUE, # Print progress during training
#   returnData = FALSE, # Don't save the training data (saves memory)
#   savePredictions = "none", # Save predictions from the best model
#   returnResamp = "final", # Save resampling results from the best model
#   allowParallel = TRUE, # Enable parallel processing (if available)
#   predictionBounds = c(0, NA)
# )

# tune_control_boot_rf <- trainControl(
#   method = "boot", # bootstrapping
#   number = 10,
#   # repeats = 3, # Repeat cross-validation 3 times (more robust)
#   verboseIter = TRUE, # Print progress during training
#   returnData = FALSE, # Don't save the training data (saves memory)
#   savePredictions = "none", # Save predictions from the best model
#   returnResamp = "final", # Save resampling results from the best model
#   allowParallel = TRUE, # Enable parallel processing (if available)
#   predictionBounds = c(0, NA)
# )

tune_control_cv_rf <- trainControl(
  method = "cv", # out of bag
  number = 10, # 10-fold cross-validation
  # repeats = 3, # Repeat cross-validation 3 times (more robust)
  verboseIter = TRUE, # Print progress during training
  returnData = FALSE, # Don't save the training data (saves memory)
  savePredictions = "none", # Save predictions from the best model
  returnResamp = "final", # Save resampling results from the best model
  allowParallel = TRUE, # Enable parallel processing (if available)
  predictionBounds = c(0, NA)
)

##### Train the Random Forest Models #####
### Function to parallelize the model-fitting process
fit_RF_models <- function(i, train_data, train_y, tune_grid, tune_control_grid_rf){
  ### training model
  set.seed(802)
  train(
    x = train_data[[i]],
    y = train_y,
    method = "rf",
    trControl = tune_control_grid_rf,
    tuneGrid = tune_grid[[i]],
    # importance = "impurity",
    importance = TRUE,
    ntree = 2500
    # num.trees = 2500
  )
}

### making lists of training dfs and tune_grids to input to function
# AKTrain_df_list <- list(AKTrain_OctApr_x, AKTrain_OctMay_x, AKTrain_SepMay_x, AKTrain_OctAprDecFeb_x, AKTrain_OctMayDecFeb_x, AKTrain_SepMayDecFeb_x, AKTrain_OctAprAspect_x, AKTrain_OctAprAspect_Nosrad_x, AKTrain_OctAprAspect_Nosradtmin_x, AKTrain_OctAprAspect_Nosradtmintmax_x, AKTrain_OctAprAspect_Notmin_x, AKTrain_OctAprAspect_Notmintmax_x, AKTrain_OctAprAspect_Notmax_x, AKTrain_OctMayAspect_Nosradtmin_x, AKTrain_OctMayAspect_Nosradtmintmax_x, AKTrain_OctMayAspect_Notmintmax_x, AKTrain_SepMayAspect_Nosradtmin_x, AKTrain_SepMayAspect_Nosradtmintmax_x, AKTrain_SepMayAspect_Notmintmax_x, AKTrain_OctMayDecFebAspect_x)
CONUSTrain_df_list <- list(CONUSTrain_OctApr_x, CONUSTrain_OctMay_x, CONUSTrain_SepMay_x, CONUSTrain_OctAprDecFeb_x, CONUSTrain_OctMayDecFeb_x, CONUSTrain_SepMayDecFeb_x, CONUSTrain_OctAprAspect_x, CONUSTrain_OctAprAspect_Nosrad_x, CONUSTrain_OctAprAspect_Nosradtmin_x, CONUSTrain_OctAprAspect_Nosradtmintmax_x, CONUSTrain_OctAprAspect_Notmin_x, CONUSTrain_OctAprAspect_Notmintmax_x, CONUSTrain_OctAprAspect_Notmax_x, CONUSTrain_OctMayAspect_Nosradtmin_x, CONUSTrain_OctMayAspect_Nosradtmintmax_x, CONUSTrain_OctMayAspect_Notmintmax_x, CONUSTrain_SepMayAspect_Nosradtmin_x, CONUSTrain_SepMayAspect_Nosradtmintmax_x, CONUSTrain_SepMayAspect_Notmintmax_x, CONUSTrain_OctMayDecFebAspect_x, CONUSTrain_OctMayDecFebAspect_Nosrad_x)
### list of tune grids
tune_grid_list <- list(tune_grid_OctApr_rf, tune_grid_OctMay_rf, tune_grid_SepMay_rf, tune_grid_OctAprDecFeb_rf, tune_grid_OctMayDecFeb_rf, tune_grid_SepMayDecFeb_rf, tune_grid_OctAprAspect_rf, tune_grid_OctAprAspect_Nosrad_rf, tune_grid_OctAprAspect_Nosradtmin_rf, tune_grid_OctAprAspect_Nosradtmintmax_rf, tune_grid_OctAprAspect_Notmin_rf, tune_grid_OctAprAspect_Notmintmax_rf, tune_grid_OctAprAspect_Notmax_rf, tune_grid_OctMayAspect_Nosradtmin_rf, tune_grid_OctMayAspect_Nosradtmintmax_rf, tune_grid_OctMayAspect_Notmintmax_rf, tune_grid_SepMayAspect_Nosradtmin_rf, tune_grid_SepMayAspect_Nosradtmintmax_rf, tune_grid_SepMayAspect_Notmintmax_rf, tune_grid_OctMayDecFebAspect_rf, tune_grid_OctMayDecFebAspect_Nosrad_rf)
### list of test dfs to evaluate model
# AKTest_df_list <- list(AKTest_OctApr_x, AKTest_OctMay_x, AKTest_SepMay_x, AKTest_OctAprDecFeb_x, AKTest_OctMayDecFeb_x, AKTest_SepMayDecFeb_x, AKTest_OctAprAspect_x, AKTest_OctAprAspect_Nosrad_x, AKTest_OctAprAspect_Nosradtmin_x, AKTest_OctAprAspect_Nosradtmintmax_x, AKTest_OctAprAspect_Notmin_x, AKTest_OctAprAspect_Notmintmax_x, AKTest_OctAprAspect_Notmax_x, AKTest_OctMayAspect_Nosradtmin_x, AKTest_OctMayAspect_Nosradtmintmax_x, AKTest_OctMayAspect_Notmintmax_x, AKTest_SepMayAspect_Nosradtmin_x, AKTest_SepMayAspect_Nosradtmintmax_x, AKTest_SepMayAspect_Notmintmax_x, AKTest_OctMayDecFebAspect_x)
CONUSTest_df_list <- list(CONUSTest_OctApr_x, CONUSTest_OctMay_x, CONUSTest_SepMay_x, CONUSTest_OctAprDecFeb_x, CONUSTest_OctMayDecFeb_x, CONUSTest_SepMayDecFeb_x, CONUSTest_OctAprAspect_x, CONUSTest_OctAprAspect_Nosrad_x, CONUSTest_OctAprAspect_Nosradtmin_x, CONUSTest_OctAprAspect_Nosradtmintmax_x, CONUSTest_OctAprAspect_Notmin_x, CONUSTest_OctAprAspect_Notmintmax_x, CONUSTest_OctAprAspect_Notmax_x, CONUSTest_OctMayAspect_Nosradtmin_x, CONUSTest_OctMayAspect_Nosradtmintmax_x, CONUSTest_OctMayAspect_Notmintmax_x, CONUSTest_SepMayAspect_Nosradtmin_x, CONUSTest_SepMayAspect_Nosradtmintmax_x, CONUSTest_SepMayAspect_Notmintmax_x, CONUSTest_OctMayDecFebAspect_x, CONUSTest_OctMayDecFebAspect_Nosrad_x)


##### training Alaska models resampled using OOB error #####
# AKOOBModels <- pmclapply(X = seq(1:length(AKTrain_df_list)), FUN = fit_RF_models, train_data = AKTrain_df_list, train_y = AKTrain_y, tune_grid = tune_grid_list, tune_control_grid_rf = tune_control_oob_rf, mc.cores = ifelse(length(tune_grid_list) > detectCores(), ceiling(detectCores() / 1.5), length(tune_grid_list)), mc.silent = FALSE)
# ### evaluating test metrics
# print("printing model metrics for Alaska OOB models")
# for (x in 1:length(AKOOBModels)){
#   print(x)
#   temp_model <- AKOOBModels[[x]]
#   print(temp_model$bestTune)
#   temp_results <- temp_model$results
#   temp_best_model <- temp_results |>
#     filter(mtry == temp_model$bestTune$mtry) # & min.node.size == temp_model$bestTune$min.node.size)
#   print(temp_best_model)
#   set.seed(802)
#   RF_model_preds <- predict(temp_model, AKTest_df_list[[x]])
#   print(rmse(AKTest_y, RF_model_preds))
#   print(mae(AKTest_y, RF_model_preds))
# }
##### training models resampled using bootstrapping #####
# BootModels <- pmclapply(X = 1:19, FUN = fit_RF_models, train_data = train_df_list, tune_grid = tune_grid_list, tune_control_grid_rf = tune_control_boot_rf, mc.cores = length(tune_grid_list), mc.silent = FALSE)
# ### evaluating test metrics
# for (x in 1:length(BootModels)){
#   print(x)
#   temp_model <- BootModels[[x]]
#   print(temp_model$bestTune)
#   print(temp_model$results)
#   set.seed(802)
#   RF_model_preds <- predict(temp_model, test_df_list[[x]])
#   print(rmse(test_y, RF_model_preds))
#   print(mae(AKTest_y, RF_model_preds))
# }
##### training Alaska models resampled using cross-validation #####
# AKCVModels <- pmclapply(X = seq(1:length(AKTrain_df_list)), FUN = fit_RF_models, train_data = AKTrain_df_list, train_y = AKTrain_y, tune_grid = tune_grid_list, tune_control_grid_rf = tune_control_cv_rf, mc.cores = ifelse(length(tune_grid_list) > detectCores(), ceiling(detectCores() / 1.5), length(tune_grid_list)), mc.silent = FALSE)
# ### evaluating test metrics
# print("printing model metrics for Alaska Cross Validation models")
# for (x in 1:length(AKCVModels)){
#   print(x)
#   temp_model <- AKCVModels[[x]]
#   print(temp_model$bestTune)
#   temp_results <- temp_model$results
#   temp_best_model <- temp_results |>
#     filter(mtry == temp_model$bestTune$mtry) # & min.node.size == temp_model$bestTune$min.node.size)
#   print(temp_best_model)
#   set.seed(802)
#   RF_model_preds <- predict(temp_model, AKTest_df_list[[x]])
#   print(rmse(AKTest_y, RF_model_preds))
#   print(mae(AKTest_y, RF_model_preds))
# }

##### training Alaska models resampled using Leave One Out cross-validation (LOOCV) #####
# AKLOOCVModels <- pmclapply(X = seq(1:length(AKTrain_df_list)), FUN = fit_RF_models, train_data = AKTrain_df_list, train_y = AKTrain_y, tune_grid = tune_grid_list, tune_control_grid_rf = tune_control_loocv_rf, mc.cores = ifelse(length(tune_grid_list) > detectCores(), ceiling(detectCores() / 1.5), length(tune_grid_list)), mc.silent = FALSE)
# ### evaluating test metrics
# print("printing model metrics for Cross Validation models")
# for (x in 1:length(AKLOOCVModels)){
#   print(x)
#   temp_model <- AKLOOCVModels[[x]]
#   print(temp_model$bestTune)
#   temp_results <- temp_model$results
#   temp_best_model <- temp_results |>
#     filter(mtry == temp_model$bestTune$mtry) # & min.node.size == temp_model$bestTune$min.node.size)
#   print(temp_best_model)
#   set.seed(802)
#   RF_model_preds <- predict(temp_model, AKTest_df_list[[x]])
#   print(rmse(AKTest_y, RF_model_preds))
# }


##### training CONUS models resampled using OOB error #####
# print("training CONUS models")
# CONUSOOBModels <- pmclapply(X = seq(1:length(CONUSTrain_df_list)), FUN = fit_RF_models, train_data = CONUSTrain_df_list, train_y = CONUSTrain_y, tune_grid = tune_grid_list, tune_control_grid_rf = tune_control_oob_rf, mc.cores = ifelse(length(tune_grid_list) > detectCores(), ceiling(detectCores() / 1.5), length(tune_grid_list)), mc.silent = FALSE)
# ### evaluating test metrics
# print("printing model metrics for CONUS OOB models")
# for (x in 1:length(CONUSOOBModels)){
#   print(x)
#   temp_model <- CONUSOOBModels[[x]]
#   print(temp_model$bestTune)
#   temp_results <- temp_model$results
#   temp_best_model <- temp_results |>
#     filter(mtry == temp_model$bestTune$mtry) # & min.node.size == temp_model$bestTune$min.node.size)
#   print(temp_best_model)
#   set.seed(802)
#   RF_model_preds <- predict(temp_model, CONUSTest_df_list[[x]])
#   print(rmse(CONUSTest_y, RF_model_preds))
#   print(mae(CONUSTest_y, RF_model_preds))
# }

##### training CONUS models resampled using cross-validation #####
# CONUSCVModels <- pmclapply(X = seq(1:length(CONUSTrain_df_list)), FUN = fit_RF_models, train_data = CONUSTrain_df_list, train_y = CONUSTrain_y, tune_grid = tune_grid_list, tune_control_grid_rf = tune_control_cv_rf, mc.cores = ifelse(length(tune_grid_list) > detectCores(), ceiling(detectCores() / 1.5), length(tune_grid_list)), mc.silent = FALSE)
CONUSCVModels <- pmclapply(X = seq(1:length(CONUSTrain_df_list)), FUN = fit_RF_models, train_data = CONUSTrain_df_list, train_y = CONUSTrain_y, tune_grid = tune_grid_list, tune_control_grid_rf = tune_control_cv_rf, mc.cores = 5, mc.silent = FALSE)
### evaluating test metrics
print("printing model metrics for Cross Validation models")
for (x in 1:length(CONUSCVModels)){
  print(x)
  temp_model <- CONUSCVModels[[x]]
  print(temp_model$bestTune)
  temp_results <- temp_model$results
  temp_best_model <- temp_results |>
    filter(mtry == temp_model$bestTune$mtry) # & min.node.size == temp_model$bestTune$min.node.size)
  print(temp_best_model)
  set.seed(802)
  RF_model_preds <- predict(temp_model, CONUSTest_df_list[[x]])
  print(rmse(CONUSTest_y, RF_model_preds))
  print(mae(CONUSTest_y, RF_model_preds))
}

##### training CONUS models resampled using Leave One Out cross-validation (LOOCV) #####
# CONUSLOOCVModels <- pmclapply(X = seq(1:length(CONUSTrain_df_list)), FUN = fit_RF_models, train_data = CONUSTrain_df_list, train_y = CONUSTrain_y, tune_grid = tune_grid_list, tune_control_grid_rf = tune_control_loocv_rf, mc.cores = ifelse(length(tune_grid_list) > detectCores(), ceiling(detectCores() / 1.5), length(tune_grid_list)), mc.silent = FALSE)
# ### evaluating test metrics
# print("printing model metrics for Cross Validation models")
# for (x in 1:length(CONUSLOOCVModels)){
#   print(x)
#   temp_model <- CONUSLOOCVModels[[x]]
#   print(temp_model$bestTune)
#   print(temp_model$results)
#   set.seed(802)
#   RF_model_preds <- predict(temp_model, CONUSTest_df_list[[x]])
#   print(rmse(test_y, RF_model_preds))
# }

# write_rds(AKOOBModels, here("Data", "FittedModels", "Caret", "RF", "AK_OOBModelList.rds"), compress = "gz")
# write_rds(BootModels, here("Data", "FittedModels", "Caret", "RF", "BootModelList.rds"), compress = "gz")
# write_rds(AKCVModels, here("Data", "FittedModels", "Caret", "RF", "AK_CVModelList.rds"), compress = "gz")
# write_rds(AKLOOCVModels, here("Data", "FittedModels", "Caret", "RF", "AK_LOOCVModelList.rds"), compress = "gz")
# write_rds(CONUSOOBModels, here("Data", "FittedModels", "Caret", "RF", "CONUS_OOBModelList.rds"), compress = "gz")
# write_rds(CONUSCVModels, here("Data", "FittedModels", "Caret", "RF", "CONUS_CVModelList.rds"), compress = "gz")
# write_rds(CONUSLOOCVModels, here("Data", "FittedModels", "Caret", "RF", "CONUS_LOOCVModelList.rds"), compress = "gz")
# saveRDS(OOBModels, here("Data", "FittedModels", "Caret", "RF", "OOBModelList.rds"))
# saveRDS(BootModels, here("Data", "FittedModels", "Caret", "RF", "BootModelList.rds"))
# saveRDS(CVModels, here("Data", "FittedModels", "Caret", "RF", "CVModelList.rds"))



##### saving model to avoid fitting it again #####
# saveRDS(RF_model1, file = here("Data", "FittedModels", "Caret", "RF", "OctApr_RFModel_1km.rds"))
# saveRDS(RF_model2, file = here("Data", "FittedModels", "Caret", "RF", "OctMay_RFModel_1km.rds"))
# saveRDS(RF_model3, file = here("Data", "FittedModels", "Caret", "RF", "SepMay_RFModel_1km.rds"))
# saveRDS(RF_model4, file = here("Data", "FittedModels", "Caret", "RF", "OctAprDecFeb_RFModel_1km.rds"))
# saveRDS(RF_model5, file = here("Data", "FittedModels", "Caret", "RF", "OctMayDecFeb_RFModel_1km.rds"))
# saveRDS(RF_model6, file = here("Data", "FittedModels", "Caret", "RF", "SepMayDecFeb_RFModel_1km.rds"))
# saveRDS(RF_model7, file = here("Data", "FittedModels", "Caret", "RF", "OctAprAspect_RFModel_1km.rds"))
# saveRDS(RF_model8, file = here("Data", "FittedModels", "Caret", "RF", "OctAprAspect_Nosrad_RFModel_1km.rds"))
# saveRDS(RF_model9, file = here("Data", "FittedModels", "Caret", "RF", "OctAprAspect_Nosradtmin_RFModel_1km.rds"))
# saveRDS(RF_model10, file = here("Data", "FittedModels", "Caret", "RF", "OctAprAspect_Nosradtmintmax_RFModel_1km.rds"))
# saveRDS(RF_model11, file = here("Data", "FittedModels", "Caret", "RF", "OctAprAspect_Notmin_RFModel_1km.rds"))
# saveRDS(RF_model12, file = here("Data", "FittedModels", "Caret", "RF", "OctAprAspect_Notmintmax_RFModel_1km.rds"))
# saveRDS(RF_model13, file = here("Data", "FittedModels", "Caret", "RF", "OctAprAspect_Notmax_RFModel_1km.rds"))
# saveRDS(RF_model14, file = here("Data", "FittedModels", "Caret", "RF", "OctMayAspect_Nosradtmin_RFModel_1km.rds"))
# saveRDS(RF_model15, file = here("Data", "FittedModels", "Caret", "RF", "OctMayAspect_Nosradtmintmax_RFModel_1km.rds"))
# saveRDS(RF_model16, file = here("Data", "FittedModels", "Caret", "RF", "OctMayAspect_Notmintmax_RFModel_1km.rds"))
# saveRDS(RF_model17, file = here("Data", "FittedModels", "Caret", "RF", "SepMayAspect_Nosradtmin_RFModel_1km.rds"))
# saveRDS(RF_model18, file = here("Data", "FittedModels", "Caret", "RF", "SepMayAspect_Nosradtmintmax_RFModel_1km.rds"))
# saveRDS(RF_model19, file = here("Data", "FittedModels", "Caret", "RF", "SepMayAspect_Notmintmax_RFModel_1km.rds"))
# 
# ##### Examine the Best Parameters and Results #####
# ### Model 1
# print(RF_model1$bestTune)
# print(RF_model1$results)
# set.seed(802)
# RF_model1_preds <- predict(RF_model1, test_OctApr_x)
# rmse(test_y, RF_model1_preds)
# ### Model 2
# print(RF_model2$bestTune)
# print(RF_model2$results)
# set.seed(802)
# RF_model2_preds <- predict(RF_model2, test_OctMay_x)
# rmse(test_y, RF_model2_preds)
# ### Model 3
# print(RF_model3$bestTune)
# print(RF_model3$results)
# set.seed(802)
# RF_model3_preds <- predict(RF_model3, test_SepMay_x)
# rmse(test_y, RF_model3_preds)
# ### Model 4
# print(RF_model4$bestTune)
# print(RF_model4$results)
# set.seed(802)
# RF_model4_preds <- predict(RF_model4, test_OctAprDecFeb_x)
# rmse(test_y, RF_model4_preds)
# ### Model 5
# print(RF_model5$bestTune)
# print(RF_model5$results)
# set.seed(802)
# RF_model5_preds <- predict(RF_model5, test_OctMayDecFeb_x)
# rmse(test_y, RF_model5_preds)
# ### Model 6
# print(RF_model6$bestTune)
# print(RF_model6$results)
# set.seed(802)
# RF_model6_preds <- predict(RF_model6, test_SepMayDecFeb_x)
# rmse(test_y, RF_model6_preds)
# ### Model 7
# print(RF_model7$bestTune)
# print(RF_model7$results)
# set.seed(802)
# RF_model7_preds <- predict(RF_model7, test_OctAprAspect_x)
# rmse(test_y, RF_model7_preds)
# ### Model 8
# print(RF_model8$bestTune)
# print(RF_model8$results)
# set.seed(802)
# RF_model8_preds <- predict(RF_model8, test_OctAprAspect_Nosrad_x)
# rmse(test_y, RF_model8_preds)
# ### Model 9
# print(RF_model9$bestTune)
# print(RF_model9$results)
# set.seed(802)
# RF_model9_preds <- predict(RF_model9, test_OctAprAspect_Nosradtmin_x)
# rmse(test_y, RF_model9_preds)
# ### Model 10
# # RF_model10 <- read_rds(here("Data", "FittedModels", "Caret", "RF", "OctAprAspect_Nosradtmintmax_RFModel_1km.rds"))
# print(RF_model10$bestTune)
# print(RF_model10$results)
# set.seed(802)
# RF_model10_preds <- predict(RF_model10, test_OctAprAspect_Nosradtmintmax_x)
# set.seed(802)
# Snotel_SWEpreds <- Snotel_test |>
#   mutate(SWE_preds = predict(RF_model10, Snotel_test),
#          SWE_ae = ae(peak_swe, SWE_preds)) |>
#   select(site_id, site_name, state, WaterYear, peak_swe, SWE_preds, SWE_ae)
# Snotel_SWEpreds_WYErrorMeans <- Snotel_SWEpreds |>
#   dplyr::group_by(WaterYear) |>
#   summarise(mean_SWE_ae = mean(SWE_ae),
#             median_SWE_ae = median(SWE_ae))
# plot(Snotel_SWEpreds$peak_swe, Snotel_SWEpreds$SWE_ae, main = "SNOTEL Peak SWE Value on X-axis, Absolute Error of Predictions on Y")
# hist(Snotel_SWEpreds$SWE_ae, main = "Histogram of Absolute Error at All Annual SNOTEL Sites")
# plot(Snotel_SWEpreds$WaterYear, Snotel_SWEpreds$SWE_ae, main = "Water Year on X-axis, Absolute Error of Predictions on Y")
# plot(Snotel_SWEpreds_WYErrorMeans$WaterYear, Snotel_SWEpreds_WYErrorMeans$mean_SWE_ae, main = "Water Year on X-axis, Mean Absolute Error of Predictions on Y")
# plot(Snotel_SWEpreds_WYErrorMeans$WaterYear, Snotel_SWEpreds_WYErrorMeans$median_SWE_ae, main = "Water Year on X-axis, Median Absolute Error of Predictions on Y")
# rmse(test_y, RF_model10_preds)
# ### Model 11
# print(RF_model11$bestTune)
# print(RF_model11$results)
# set.seed(802)
# RF_model11_preds <- predict(RF_model11, test_OctAprAspect_Notmin_x)
# rmse(test_y, RF_model11_preds)
# ### Model 12
# print(RF_model12$bestTune)
# print(RF_model12$results)
# set.seed(802)
# RF_model12_preds <- predict(RF_model12, test_OctAprAspect_Notmintmax_x)
# rmse(test_y, RF_model12_preds)
# ### Model 13
# print(RF_model13$bestTune)
# print(RF_model13$results)
# set.seed(802)
# RF_model13_preds <- predict(RF_model13, test_OctAprAspect_Notmax_x)
# rmse(test_y, RF_model13_preds)
# ### Model 14
# print(RF_model14$bestTune)
# print(RF_model14$results)
# set.seed(802)
# RF_model14_preds <- predict(RF_model14, test_OctMayAspect_Nosradtmin_x)
# rmse(test_y, RF_model14_preds)
# ### Model 15
# print(RF_model15$bestTune)
# print(RF_model15$results)
# set.seed(802)
# RF_model15_preds <- predict(RF_model15, test_OctMayAspect_Nosradtmintmax_x)
# rmse(test_y, RF_model15_preds)
# ### Model 16
# print(RF_model16$bestTune)
# print(RF_model16$results)
# set.seed(802)
# RF_model16_preds <- predict(RF_model16, test_OctMayAspect_Notmintmax_x)
# rmse(test_y, RF_model16_preds)
# ### Model 17
# print(RF_model17$bestTune)
# print(RF_model17$results)
# set.seed(802)
# RF_model17_preds <- predict(RF_model17, test_SepMayAspect_Nosradtmin_x)
# rmse(test_y, RF_model17_preds)
# ### Model 18
# print(RF_model18$bestTune)
# print(RF_model18$results)
# set.seed(802)
# RF_model18_preds <- predict(RF_model18, test_SepMayAspect_Nosradtmintmax_x)
# rmse(test_y, RF_model18_preds)
# ### Model 19
# print(RF_model19$bestTune)
# print(RF_model19$results)
# set.seed(802)
# RF_model19_preds <- predict(RF_model19, test_SepMayAspect_Notmintmax_x)
# rmse(test_y, RF_model19_preds)
