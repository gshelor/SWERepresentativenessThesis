##### SWE Random Forest Error Analysis #####
##### Loading packages #####
library(pacman)
p_load(here, tidyverse, randomForest, gt, gtExtras, sf, terra, sfext, rsample, caret, parallel, grateful, ModelMetrics, Metrics)
# options(mc.cores = parallel::detectCores())

##### reading in data #####
### CONUS AOI
CONUS_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "CONUS", "CONUS_AOI.gpkg"))
### SNOTEL Climatologies used to fit the model
# Snotel_Clim_sf <- read_sf(here("Data", "SNOTEL", "Combined", "GIS", "SnotelCombined_ClimatologyCovars.gpkg"))

### SNOTEL Annual values
SnotelCONUS_sf <- read_sf(here("Data", "SNOTEL", "CONUS", "GIS", "GPKG", "Annual", "SnotelCONUSAnnualZeroPts_PeakSWE_Covars.gpkg")) |>
  drop_na(OctApr_prcpSumCDMSum, OctApr_tmeanCDMSum, elevation, slope, aspect, landcover_triclass)
### SNOTEL Annual values
SnotelAK_sf <- read_sf(here("Data", "SNOTEL", "Alaska", "GIS", "GPKG", "Annual", "SnotelAK_PeakSWE_Covars.gpkg")) |>
  drop_na(OctMay_prcpSumCDMSum, OctMay_tmeanCDMSum, DecFeb_prcpSumCDMSum, DecFeb_tmeanCDMSum, elevation, slope, aspect, OctMay_tminMean, OctMay_tmaxMean, landcover_triclass)

### converting to df for model fitting process
SnotelCONUS_df <- sf_to_df(SnotelCONUS_sf) |>
  drop_na(OctApr_prcpSumCDMSum, OctApr_tmeanCDMSum, elevation, slope, aspect, landcover_triclass)
### making sure landcover is type integer
SnotelCONUS_df$landcover <- as.integer(SnotelCONUS_df$landcover)
SnotelCONUS_df$landcover_triclass <- as.integer(SnotelCONUS_df$landcover_triclass)

##### splitting data into training and testing datasets #####
set.seed(802)
SnotelCONUS_split <- initial_split(SnotelCONUS_df, prop = 0.75, strata = peak_swe)
SnotelCONUS_train <- training(SnotelCONUS_split)
SnotelCONUS_test <- testing(SnotelCONUS_split)

### OctApr aggregates only but with aspect, and srad, tmin, tmax removed
train_OctAprAspect_Nosradtmintmax_x <- SnotelCONUS_train[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "landcover_triclass")]
test_OctAprAspect_Nosradtmintmax_x <- SnotelCONUS_test[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "landcover_triclass")]
train_y <- SnotelCONUS_train$peak_swe
test_y <- SnotelCONUS_test$peak_swe



##### Final Model #####
FinalModel <- read_rds(here("Data", "FittedModels", "Caret", "RF", "AnnualCONUSRFModel_oob.rds"))
print(FinalModel$bestTune)
print(FinalModel$results)
set.seed(802)
FinalModel_preds <- predict(FinalModel, test_OctAprAspect_Nosradtmintmax_x)
set.seed(802)
Snotel_SWEpreds <- SnotelCONUS_test |>
  mutate(SWE_preds = predict(FinalModel, SnotelCONUS_test),
         SWE_error = peak_swe - SWE_preds,
         SWE_ae = ae(peak_swe, SWE_preds),
         SWE_se = se(peak_swe, SWE_preds),
         SWE_relative_error = SWE_ae / peak_swe,
         SWE_pct_error = SWE_relative_error * 100) |>
  select(site_id, site_name, state, WaterYear, peak_swe, SWE_preds, SWE_error, SWE_ae, SWE_se, SWE_relative_error, OctApr_prcpSumCDMSum, OctApr_tmeanCDMSum, elevation, slope, aspect, landcover_triclass)
Snotel_SWEpreds_WYErrorMeans <- Snotel_SWEpreds |>
  dplyr::group_by(WaterYear) |>
  summarise(mean_SWE_ae = mean(SWE_ae),
            median_SWE_ae = median(SWE_ae),
            mean_SWE_se = mean(SWE_se),
            median_SWE_se = median(SWE_se))

### grouping error metrics by site
Snotel_SWEpreds_SiteMeans <- Snotel_SWEpreds |>
  dplyr::group_by(site_name) |>
  summarise(site_id = mean(site_id),
            mean_SWE_ae = mean(SWE_ae),
            median_SWE_ae = median(SWE_ae),
            mean_SWE_se = mean(SWE_se),
            median_SWE_se = median(SWE_se),
            SWE_sd = sd(peak_swe),
            SWE_mean = mean(peak_swe),
            SWE_COV = SWE_sd / SWE_mean,
            elevation = mean(elevation),
            mean_tmeanCDM = mean(OctApr_tmeanCDMSum))

plot(Snotel_SWEpreds_SiteMeans$SWE_mean, Snotel_SWEpreds_SiteMeans$SWE_COV, main = "Mean SWE_max for Each SNOTEL Site vs SWE Coeffiecient of Variation")

plot(Snotel_SWEpreds_SiteMeans$elevation, Snotel_SWEpreds_SiteMeans$SWE_COV, main = "Elevation vs SWE Coeffiecient of Variation")
plot(Snotel_SWEpreds_SiteMeans$mean_tmeanCDM, Snotel_SWEpreds_SiteMeans$SWE_COV, main = "Mean tmean CDM Values vs SWE Coeffiecient of Variation")

plot(Snotel_SWEpreds_SiteMeans$site_id, Snotel_SWEpreds_SiteMeans$mean_SWE_ae)

### error metrics
median(Snotel_SWEpreds$SWE_ae)
mean(Snotel_SWEpreds$SWE_ae)
median(Snotel_SWEpreds$SWE_error)
mean(Snotel_SWEpreds$SWE_error)

hist(Snotel_SWEpreds$SWE_ae, main = "Histogram of Absolute Error at Test Annual SNOTEL Sites")
hist(Snotel_SWEpreds$SWE_error, main = "Histogram of Error at Test Annual SNOTEL Sites")
plot(Snotel_SWEpreds$peak_swe, Snotel_SWEpreds$SWE_ae, main = "SNOTEL Peak SWE Value on X-axis, Absolute Error of Test Predictions on Y")
plot(Snotel_SWEpreds$peak_swe, Snotel_SWEpreds$SWE_relative_error, main = "SNOTEL Peak SWE Value on X-axis, Relative Error of Test Predictions on Y")

plot(Snotel_SWEpreds$WaterYear, Snotel_SWEpreds$SWE_ae, main = "Water Year on X-axis, Absolute Error of Test Predictions on Y")
plot(Snotel_SWEpreds$elevation, Snotel_SWEpreds$SWE_ae, main = "Elevation vs Test Absolute Error")
plot(Snotel_SWEpreds$slope, Snotel_SWEpreds$SWE_ae, main = "Slope vs Test Absolute Error")
plot(Snotel_SWEpreds$OctApr_tmeanCDMSum, Snotel_SWEpreds$SWE_ae, main = "tmean CDMs vs Test Absolute Error")
plot(Snotel_SWEpreds$OctApr_prcpSumCDMSum, Snotel_SWEpreds$SWE_ae, main = "prcpSum CDMs vs Test Absolute Error")
plot(Snotel_SWEpreds$landcover_triclass, Snotel_SWEpreds$SWE_ae, main = "landcover vs Test Absolute Error")

cor(Snotel_SWEpreds$peak_swe, Snotel_SWEpreds$SWE_ae)
# cor(Snotel_SWEpreds$WaterYear, Snotel_SWEpreds$SWE_ae)
cor(Snotel_SWEpreds$elevation, Snotel_SWEpreds$SWE_ae)
cor(Snotel_SWEpreds$slope, Snotel_SWEpreds$SWE_ae)
cor(Snotel_SWEpreds$OctApr_tmeanCDMSum, Snotel_SWEpreds$SWE_ae)
cor(Snotel_SWEpreds$OctApr_prcpSumCDMSum, Snotel_SWEpreds$SWE_ae)
cor(Snotel_SWEpreds$landcover_triclass, Snotel_SWEpreds$SWE_ae)

plot(Snotel_SWEpreds_WYErrorMeans$WaterYear, Snotel_SWEpreds_WYErrorMeans$mean_SWE_ae, main = "Water Year on X-axis, Mean Absolute Error of Predictions on Y")
plot(Snotel_SWEpreds_WYErrorMeans$WaterYear, Snotel_SWEpreds_WYErrorMeans$median_SWE_ae, main = "Water Year on X-axis, Median Absolute Error of Predictions on Y")

##### Overall Model RMSE, zero pts included #####
rmse(test_y, FinalModel_preds)



### looking at highest error sites
HighErrorSites <- Snotel_SWEpreds_SiteMeans |>
  filter(mean_SWE_ae > quantile(mean_SWE_ae, 0.95))
# StateHighErrorCounts <- as.data.frame(table(HighErrorSites$state))
plot(HighErrorSites$SWE_mean, HighErrorSites$mean_SWE_ae)
# plot(HighErrorSites$elevation, HighErrorSites$mean_SWE_ae)
# plot(HighErrorSites$SWE_mean, HighErrorSites$SWE_COV)

### sum of error
AE_Sum <- sum(Snotel_SWEpreds$SWE_ae)
Snotel_TopSWE <- Snotel_SWEpreds |>
  arrange(desc(peak_swe))
Snotel_BottomSWE <- Snotel_SWEpreds |>
  arrange(peak_swe)
Snotel_Top5pctSWE <- Snotel_TopSWE |>
  filter(peak_swe >= quantile(peak_swe, 0.95))
Snotel_Top10pctSWE <- Snotel_TopSWE |>
  filter(peak_swe >= quantile(peak_swe, 0.9))
SWE_Top5pctErrorSum <- sum(Snotel_Top5pctSWE$SWE_ae)
SWE_Top10pctErrorSum <- sum(Snotel_Top10pctSWE$SWE_ae)

SWE_Top5pctErrorSum / AE_Sum
SWE_Top10pctErrorSum / AE_Sum

### looking at the proportion of error made up by sites in the bottom 10% of CONUS sites in peak SWE
## pretty sure it's just zero points
Snotel_Bottom5pctSWE <- Snotel_BottomSWE |>
  filter(peak_swe <= quantile(peak_swe, 0.05))
Snotel_Bottom10pctSWE <- Snotel_BottomSWE |>
  filter(peak_swe <= quantile(peak_swe, 0.1))
SWE_Bottom10pctErrorSum <- sum(Snotel_Bottom10pctSWE$SWE_ae)
SWE_Bottom5pctErrorSum <- sum(Snotel_Bottom5pctSWE$SWE_ae)
SWE_Bottom10pctErrorSum / AE_Sum
SWE_Bottom5pctErrorSum / AE_Sum

##### Evaluating Error excluding zero pts #####
Snotel_SWEpreds_SnotelOnly <- Snotel_SWEpreds |>
  filter(site_id != -999)


### error metrics
median(Snotel_SWEpreds_SnotelOnly$SWE_ae)
mean(Snotel_SWEpreds_SnotelOnly$SWE_ae)
median(Snotel_SWEpreds_SnotelOnly$SWE_error)
mean(Snotel_SWEpreds_SnotelOnly$SWE_error)
rmse(Snotel_SWEpreds_SnotelOnly$peak_swe, Snotel_SWEpreds_SnotelOnly$SWE_preds)

### sum of error
SnotelOnlyAE_Sum <- sum(Snotel_SWEpreds_SnotelOnly$SWE_ae)
SnotelOnly_TopSWE <- Snotel_SWEpreds_SnotelOnly |>
  arrange(desc(peak_swe))
SnotelOnly_BottomSWE <- Snotel_SWEpreds_SnotelOnly |>
  arrange(peak_swe)
SnotelOnly_Top5pctSWE <- SnotelOnly_TopSWE |>
  filter(peak_swe >= quantile(peak_swe, 0.95))
SnotelOnly_Top10pctSWE <- SnotelOnly_TopSWE |>
  filter(peak_swe >= quantile(peak_swe, 0.9))
SnotelOnlySWE_Top5pctErrorSum <- sum(SnotelOnly_Top5pctSWE$SWE_ae)
SnotelOnlySWE_Top10pctErrorSum <- sum(SnotelOnly_Top10pctSWE$SWE_ae)

SnotelOnlySWE_Top5pctErrorSum / SnotelOnlyAE_Sum
SnotelOnlySWE_Top10pctErrorSum / SnotelOnlyAE_Sum

SnotelOnly_Bottom5pctSWE <- SnotelOnly_BottomSWE |>
  filter(peak_swe <= quantile(peak_swe, 0.05))
SnotelOnly_Bottom10pctSWE <- SnotelOnly_BottomSWE |>
  filter(peak_swe <= quantile(peak_swe, 0.1))
SnotelOnlySWE_Bottom10pctErrorSum <- sum(SnotelOnly_Bottom10pctSWE$SWE_ae)
SnotelOnlySWE_Bottom5pctErrorSum <- sum(SnotelOnly_Bottom5pctSWE$SWE_ae)
SnotelOnlySWE_Bottom10pctErrorSum / SnotelOnlyAE_Sum
SnotelOnlySWE_Bottom5pctErrorSum / SnotelOnlyAE_Sum
