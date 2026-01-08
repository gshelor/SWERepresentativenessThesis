##### SWE Random Forest Model for CONUS #####
### script by Griffin Shelor

##### Loading packages #####
library(pacman)
p_load(here, tidyverse, randomForest, sf, terra, sfext, rsample, caret, parallel, grateful, ModelMetrics, Metrics, mcprogress)
# options(mc.cores = parallel::detectCores())

##### reading in data #####
### Alaska AOI
AK_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "Alaska", "AK_AOI.gpkg"))
### Combined SNOTEL Annual values for model fitting and evaluation
# Snotel_sf <- read_sf(here("Data", "SNOTEL", "Combined", "GIS", "SnotelCombined_AnnualCovars_AnnualZeroPts.gpkg"))


### Extracting just Alaska data
SnotelAK_sf <- read_sf(here("Data", "SNOTEL", "Alaska", "GIS", "GPKG", "Annual", "SnotelAK_PeakSWE_Covars.gpkg")) |>
  drop_na(OctApr_prcpSumCDMSum, OctApr_tmeanCDMSum, elevation, slope, aspect, OctApr_tminMean, OctApr_tmaxMean, OctApr_sradMean, landcover_triclass, OctMay_prcpSumCDMSum, OctMay_tmeanCDMSum, OctMay_tminMean, OctMay_tmaxMean, OctMay_sradMean, SepMay_prcpSumCDMSum, SepMay_tmeanCDMSum, SepMay_tminMean, SepMay_tmaxMean, SepMay_sradMean)

### converting to df for model fitting process
# Snotel_df <- sf_to_df(Snotel_sf) |>
#   drop_na(OctApr_prcpSumCDMSum, OctApr_tmeanCDMSum, elevation, slope, aspect, OctApr_tminMean, OctApr_tmaxMean, OctApr_sradMean, landcover_triclass, OctMay_prcpSumCDMSum, OctMay_tmeanCDMSum, OctMay_tminMean, OctMay_tmaxMean, OctMay_sradMean, SepMay_prcpSumCDMSum, SepMay_tmeanCDMSum, SepMay_tminMean, SepMay_tmaxMean, SepMay_sradMean)
SnotelAK_df <- sf_to_df(SnotelAK_sf)
### making sure landcover is type integer
SnotelAK_df$landcover <- as.integer(SnotelAK_df$landcover)
SnotelAK_df$landcover_triclass <- as.integer(SnotelAK_df$landcover_triclass)

##### splitting data into training and testing datasets #####
### splitting Alaska
set.seed(802)
SnotelAK_split <- initial_split(SnotelAK_df, prop = 0.75, strata = peak_swe)
SnotelAK_train <- training(SnotelAK_split)
SnotelAK_test <- testing(SnotelAK_split)
##### Alaska training and testing DFs #####
### keeping as dataframes for caret
### OctMay aggregates and DecFeb aggregates for prcpSum and tmean
### ran a second test, took srad out
AKTrain_OctMayDecFebAspect_x <- SnotelAK_train[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "aspect", "OctMay_tminMean", "OctMay_tmaxMean", "landcover_triclass")]
AKTrain_y <- SnotelAK_train$peak_swe
AKTest_OctMayDecFebAspect_x <- SnotelAK_test[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "aspect", "OctMay_tminMean", "OctMay_tmaxMean", "landcover_triclass")]
AKTest_y <- SnotelAK_test$peak_swe

##### Hyperparameter Tuning Grid #####
tune_grid_OctMayDecFebAspect_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(AKTrain_OctMayDecFebAspect_x), by = 1) #,
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
AKTrain_df_list <- list(AKTrain_OctMayDecFebAspect_x)
### list of tune grids
tune_grid_list <- list(tune_grid_OctMayDecFebAspect_rf)
### list of test dfs to evaluate model
AKTest_df_list <- list(AKTest_OctMayDecFebAspect_x)


##### training Alaska models resampled using OOB error #####
AKOOBModels <- pmclapply(X = seq(1:length(AKTrain_df_list)), FUN = fit_RF_models, train_data = AKTrain_df_list, train_y = AKTrain_y, tune_grid = tune_grid_list, tune_control_grid_rf = tune_control_oob_rf, mc.cores = ifelse(length(tune_grid_list) > detectCores(), ceiling(detectCores() / 1.5), length(tune_grid_list)), mc.silent = FALSE)
### evaluating test metrics
print("printing model metrics for Alaska OOB models")
for (x in 1:length(AKOOBModels)){
  print(x)
  temp_model <- AKOOBModels[[x]]
  print(temp_model$bestTune)
  temp_results <- temp_model$results
  temp_best_model <- temp_results |>
    filter(mtry == temp_model$bestTune$mtry) # & min.node.size == temp_model$bestTune$min.node.size)
  print(temp_best_model)
  set.seed(802)
  RF_model_preds <- predict(temp_model, AKTest_df_list[[x]])
  print(rmse(AKTest_y, RF_model_preds))
  print(mae(AKTest_y, RF_model_preds))
}
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
AKCVModels <- pmclapply(X = seq(1:length(AKTrain_df_list)), FUN = fit_RF_models, train_data = AKTrain_df_list, train_y = AKTrain_y, tune_grid = tune_grid_list, tune_control_grid_rf = tune_control_cv_rf, mc.cores = ifelse(length(tune_grid_list) > detectCores(), ceiling(detectCores() / 1.5), length(tune_grid_list)), mc.silent = FALSE)
### evaluating test metrics
print("printing model metrics for Alaska Cross Validation models")
for (x in 1:length(AKCVModels)){
  print(x)
  temp_model <- AKCVModels[[x]]
  print(temp_model$bestTune)
  temp_results <- temp_model$results
  temp_best_model <- temp_results |>
    filter(mtry == temp_model$bestTune$mtry) # & min.node.size == temp_model$bestTune$min.node.size)
  print(temp_best_model)
  set.seed(802)
  RF_model_preds <- predict(temp_model, AKTest_df_list[[x]])
  print(rmse(AKTest_y, RF_model_preds))
  print(mae(AKTest_y, RF_model_preds))
}

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
