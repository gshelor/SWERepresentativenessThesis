##### Trying XGBoost for modeling maximum SWE (mSWE) #####

##### Loading packages, reading in data #####
library(pacman)
p_load(tidyverse, xgboost, caret, sf, terra, here, sfext, rsample, parallel, mcprogress, fastDummies, ModelMetrics, Metrics)
# options(mc.cores = parallel::detectCores())

### AOI
CONUS_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "CONUS", "CONUS_AOI.gpkg"))
### SNOTEL Annual values for model fitting and evaluation
Snotel_sf <- read_sf(here("Data", "SNOTEL", "Combined", "GIS", "SnotelCombined_AnnualCovars_AnnualZeroPts.gpkg"))

### converting to df for model fitting process
Snotel_df <- sf_to_df(Snotel_sf) |>
  drop_na(OctApr_prcpSumCDMSum, OctApr_tmeanCDMSum, elevation, slope, aspect, OctApr_tminMean, OctApr_tmaxMean, OctApr_sradMean, landcover_triclass, OctMay_prcpSumCDMSum, OctMay_tmeanCDMSum, OctMay_tminMean, OctMay_tmaxMean, OctMay_sradMean, SepMay_prcpSumCDMSum, SepMay_tmeanCDMSum, SepMay_tminMean, SepMay_tmaxMean, SepMay_sradMean)
### making sure landcover is type integer
Snotel_df$landcover <- as.integer(Snotel_df$landcover)
Snotel_df$landcover_triclass <- as.integer(Snotel_df$landcover_triclass)

##### splitting data into training and testing datasets #####
set.seed(802)
Snotel_split <- initial_split(Snotel_df, prop = 0.75)
Snotel_train <- training(Snotel_split)
Snotel_test <- testing(Snotel_split)

##### Separating variables for modeling and model eval #####
### keeping as dataframes for caret
### OctApr aggregates only
train_OctApr_x <- Snotel_train[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "OctApr_tminMean", "OctApr_tmaxMean", "OctApr_sradMean", "landcover_triclass")]
train_y <- Snotel_train$peak_swe
test_OctApr_x <- Snotel_test[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "OctApr_tminMean", "OctApr_tmaxMean", "OctApr_sradMean", "landcover_triclass")]
test_y <- Snotel_test$peak_swe
### OctApr aggregates only but with aspect
train_OctAprAspect_x <- Snotel_train[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tminMean", "OctApr_tmaxMean", "OctApr_sradMean", "landcover_triclass")]
test_OctAprAspect_x <- Snotel_test[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tminMean", "OctApr_tmaxMean", "OctApr_sradMean", "landcover_triclass")]
### OctApr aggregates only but with aspect, and srad removed
train_OctAprAspect_Nosrad_x <- Snotel_train[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tminMean", "OctApr_tmaxMean", "landcover_triclass")]
test_OctAprAspect_Nosrad_x <- Snotel_test[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tminMean", "OctApr_tmaxMean", "landcover_triclass")]
### OctApr aggregates only but with aspect, and srad, tmin removed
train_OctAprAspect_Nosradtmin_x <- Snotel_train[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tmaxMean", "landcover_triclass")]
test_OctAprAspect_Nosradtmin_x <- Snotel_test[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tmaxMean", "landcover_triclass")]
### OctApr aggregates only but with aspect, and srad, tmin, tmax removed
train_OctAprAspect_Nosradtmintmax_x <- Snotel_train[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "landcover_triclass")]
test_OctAprAspect_Nosradtmintmax_x <- Snotel_test[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "landcover_triclass")]
### OctApr aggregates only but with aspect, and no tmin
train_OctAprAspect_Notmin_x <- Snotel_train[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tmaxMean", "OctApr_sradMean", "landcover_triclass")]
test_OctAprAspect_Notmin_x <- Snotel_test[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tmaxMean", "OctApr_sradMean", "landcover_triclass")]
### OctApr aggregates only but with aspect, and no tmin, tmax
train_OctAprAspect_Notmintmax_x <- Snotel_train[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_sradMean", "landcover_triclass")]
test_OctAprAspect_Notmintmax_x <- Snotel_test[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_sradMean", "landcover_triclass")]
### OctApr aggregates only but with aspect, and no tmax
train_OctAprAspect_Notmax_x <- Snotel_train[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tminMean", "OctApr_sradMean", "landcover_triclass")]
test_OctAprAspect_Notmax_x <- Snotel_test[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "OctApr_tminMean", "OctApr_sradMean", "landcover_triclass")]
### OctApr aggregates and DecFeb aggregates for prcpSum and tmean
train_OctAprDecFeb_x <- Snotel_train[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "OctApr_tminMean", "OctApr_tmaxMean", "OctApr_sradMean", "landcover_triclass")]
test_OctAprDecFeb_x <- Snotel_test[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "OctApr_tminMean", "OctApr_tmaxMean", "OctApr_sradMean", "landcover_triclass")]
### Testing OctMay aggregations in model instead of OctApr
train_OctMay_x <- Snotel_train[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "elevation", "slope", "OctMay_tminMean", "OctMay_tmaxMean", "OctMay_sradMean", "landcover_triclass")]
test_OctMay_x <- Snotel_test[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "elevation", "slope", "OctMay_tminMean", "OctMay_tmaxMean", "OctMay_sradMean", "landcover_triclass")]
### OctMay aggregates only but with aspect, and srad, tmin removed
train_OctMayAspect_Nosradtmin_x <- Snotel_train[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "elevation", "slope", "aspect", "OctMay_tmaxMean", "landcover_triclass")]
test_OctMayAspect_Nosradtmin_x <- Snotel_test[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "elevation", "slope", "aspect", "OctMay_tmaxMean", "landcover_triclass")]
### OctMay aggregates only but with aspect, and srad, tmin, tmax removed
train_OctMayAspect_Nosradtmintmax_x <- Snotel_train[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "elevation", "slope", "aspect", "landcover_triclass")]
test_OctMayAspect_Nosradtmintmax_x <- Snotel_test[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "elevation", "slope", "aspect", "landcover_triclass")]
### OctMay aggregates only but with aspect, and no tmin, tmax
train_OctMayAspect_Notmintmax_x <- Snotel_train[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "elevation", "slope", "aspect", "OctMay_sradMean", "landcover_triclass")]
test_OctMayAspect_Notmintmax_x <- Snotel_test[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "elevation", "slope", "aspect", "OctMay_sradMean", "landcover_triclass")]
### OctMay aggregates and DecFeb aggregates for prcpSum and tmean
train_OctMayDecFeb_x <- Snotel_train[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "OctMay_tminMean", "OctMay_tmaxMean", "OctMay_sradMean", "landcover_triclass")]
test_OctMayDecFeb_x <- Snotel_test[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "OctMay_tminMean", "OctMay_tmaxMean", "OctMay_sradMean", "landcover_triclass")]
### Testing SepMay aggregations in model instead of OctApr
train_SepMay_x <- Snotel_train[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "elevation", "slope", "SepMay_tminMean", "SepMay_tmaxMean", "SepMay_sradMean", "landcover_triclass")]
test_SepMay_x <- Snotel_test[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "elevation", "slope", "SepMay_tminMean", "SepMay_tmaxMean", "SepMay_sradMean", "landcover_triclass")]
### SepMay aggregates and DecFeb aggregates for prcpSum and tmean
train_SepMayDecFeb_x <- Snotel_train[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "SepMay_tminMean", "SepMay_tmaxMean", "SepMay_sradMean", "landcover_triclass")]
test_SepMayDecFeb_x <- Snotel_test[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "SepMay_tminMean", "SepMay_tmaxMean", "SepMay_sradMean", "landcover_triclass")]
### SepMay aggregates only but with aspect, and srad, tmin removed
train_SepMayAspect_Nosradtmin_x <- Snotel_train[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "elevation", "slope", "aspect", "SepMay_tmaxMean", "landcover_triclass")]
test_SepMayAspect_Nosradtmin_x <- Snotel_test[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "elevation", "slope", "aspect", "SepMay_tmaxMean", "landcover_triclass")]
### SepMay aggregates only but with aspect, and srad, tmin, tmax removed
train_SepMayAspect_Nosradtmintmax_x <- Snotel_train[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "elevation", "slope", "aspect", "landcover_triclass")]
test_SepMayAspect_Nosradtmintmax_x <- Snotel_test[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "elevation", "slope", "aspect", "landcover_triclass")]
### SepMay aggregates only but with aspect, and no tmin, tmax
train_SepMayAspect_Notmintmax_x <- Snotel_train[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "elevation", "slope", "aspect", "SepMay_sradMean", "landcover_triclass")]
test_SepMayAspect_Notmintmax_x <- Snotel_test[, c("SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "elevation", "slope", "aspect", "SepMay_sradMean", "landcover_triclass")]



##### Training xgBoost Model #####
### converting matrices to be compatible with xgBoost package
# dtrain <- xgb.DMatrix(data = train_x, label = train_y)
# dtest <- xgb.DMatrix(data = test_x, label = test_y)

# params <- list(
#   objective = "reg:squarederror", # Regression task
#   eta = 0.1, # Learning rate
#   max_depth = 6, # Maximum depth of trees
#   min_child_weight = 1, # Minimum sum of instance weight needed in a child
#   subsample = 0.8, # Subsample ratio of the training instances
#   colsample_bytree = 0.8 # Subsample ratio of columns when constructing each tree
# )
# 
# set.seed(802)
# xgb_model <- xgb.train(
#   params = params,
#   data = dtrain,
#   nrounds = 500, # Number of boosting rounds
#   watchlist = list(train = dtrain, eval = dtest),
#   early_stopping_rounds = 10 # Stop if no improvement after 10 rounds
# )


tune_grid_xgb <- expand.grid(
  nrounds = seq(50, 1000, by = 50),
  max_depth = c(2:8),
  eta = c(0.01, 0.1, 0.3),
  gamma = 0,
  colsample_bytree = seq(0.5, 1, by = 0.1),
  min_child_weight = 1,
  subsample = seq(0, 1, by = 0.25)
)

tune_control_xgb <- trainControl(
  method = "cv",
  number = 10, # 10-fold cross-validation
  verboseIter = TRUE,
  predictionBounds = c(0, NA)
)

##### Train the Random Forest Models #####
### Function to parallelize the model-fitting process
fit_xgb_models <- function(i, train_data, tune_grid, tune_control){
  ### training model
  set.seed(802)
  train(
    x = train_data[[i]],
    y = train_y,
    method = "xgbTree",
    trControl = tune_control,
    tuneGrid = tune_grid
  )
}

### making lists of training dfs and tune_grids to input to function
train_df_list <- list(train_OctApr_x, train_OctMay_x, train_SepMay_x, train_OctAprDecFeb_x, train_OctMayDecFeb_x, train_SepMayDecFeb_x, train_OctAprAspect_x, train_OctAprAspect_Nosrad_x, train_OctAprAspect_Nosradtmin_x, train_OctAprAspect_Nosradtmintmax_x, train_OctAprAspect_Notmin_x, train_OctAprAspect_Notmintmax_x, train_OctAprAspect_Notmax_x, train_OctMayAspect_Nosradtmin_x, train_OctMayAspect_Nosradtmintmax_x, train_OctMayAspect_Notmintmax_x, train_SepMayAspect_Nosradtmin_x, train_SepMayAspect_Nosradtmintmax_x, train_SepMayAspect_Notmintmax_x)
### list of test dfs to evaluate model
test_df_list <- list(test_OctApr_x, test_OctMay_x, test_SepMay_x, test_OctAprDecFeb_x, test_OctMayDecFeb_x, test_SepMayDecFeb_x, test_OctAprAspect_x, test_OctAprAspect_Nosrad_x, test_OctAprAspect_Nosradtmin_x, test_OctAprAspect_Nosradtmintmax_x, test_OctAprAspect_Notmin_x, test_OctAprAspect_Notmintmax_x, test_OctAprAspect_Notmax_x, test_OctMayAspect_Nosradtmin_x, test_OctMayAspect_Nosradtmintmax_x, test_OctMayAspect_Notmintmax_x, test_SepMayAspect_Nosradtmin_x, test_SepMayAspect_Nosradtmintmax_x, test_SepMayAspect_Notmintmax_x)

### Calling Function
xgbModels <- pmclapply(X = 1:19, FUN = fit_xgb_models, train_data = train_df_list, tune_grid = tune_grid_xgb, tune_control = tune_control_xgb, mc.cores = ceiling(detectCores() / 1.5), mc.silent = FALSE)

write_rds(xgbModels, here("Data", "FittedModels", "Caret", "xgBoost", "Annual", "xgbModels.rds"), compress = "gz")
### evaluating test metrics
for (x in 1:length(xgbModels)){
  print(paste("xgB Model:", x))
  temp_model <- xgbModels[[x]]
  print(temp_model$bestTune)
  print(temp_model$results)
  set.seed(802)
  xgb_model_preds <- predict(temp_model, test_df_list[[x]])
  print(rmse(test_y, xgb_model_preds))
  print(mae(test_y, xgb_model_preds))
}

# ### train the model with cross validation and different possible hyperparameters
# trainstarttime <- Sys.time()
# set.seed(802)
# xgb_model1 <- train(
#   x = train_OctApr_x,
#   y = train_y,
#   method = "xgbTree",
#   trControl = tune_control_xgb,
#   tuneGrid = tune_grid_xgb
# )
# set.seed(802)
# xgb_model2 <- train(
#   x = train_OctMay_x,
#   y = train_y,
#   method = "xgbTree",
#   trControl = tune_control_xgb,
#   tuneGrid = tune_grid_xgb
# )
# set.seed(802)
# xgb_model3 <- train(
#   x = train_SepMay_x,
#   y = train_y,
#   method = "xgbTree",
#   trControl = tune_control_xgb,
#   tuneGrid = tune_grid_xgb
# )
# set.seed(802)
# xgb_model4 <- train(
#   x = train_OctAprDecFeb_x,
#   y = train_y,
#   method = "xgbTree",
#   trControl = tune_control_xgb,
#   tuneGrid = tune_grid_xgb
# )
# ### xgb model 5
# set.seed(802)
# xgb_model5 <- train(
#   x = train_OctMayDecFeb_x,
#   y = train_y,
#   method = "xgbTree",
#   trControl = tune_control_xgb,
#   tuneGrid = tune_grid_xgb
# )
# ### xgb model 6
# set.seed(802)
# xgb_model6 <- train(
#   x = train_SepMayDecFeb_x,
#   y = train_y,
#   method = "xgbTree",
#   trControl = tune_control_xgb,
#   tuneGrid = tune_grid_xgb
# )
# trainendtime <- Sys.time()
# trainendtime - trainstarttime
# ### saving model to avoid fitting it again
# saveRDS(xgb_model1, file = here("Data", "FittedModels", "Caret", "xgBoost", "OctApr_xgBModel.rds"))
# saveRDS(xgb_model2, file = here("Data", "FittedModels", "Caret", "xgBoost", "OctMay_xgBModel.rds"))
# saveRDS(xgb_model3, file = here("Data", "FittedModels", "Caret", "xgBoost", "SepMay_xgBModel.rds"))
# saveRDS(xgb_model4, file = here("Data", "FittedModels", "Caret", "xgBoost", "OctAprDecFeb_xgBModel.rds"))
# saveRDS(xgb_model5, file = here("Data", "FittedModels", "Caret", "xgBoost", "OctMayDecFeb_xgBModel.rds"))
# saveRDS(xgb_model6, file = here("Data", "FittedModels", "Caret", "xgBoost", "SepMayDecFeb_xgBModel.rds"))
# xgb_model1 <- readRDS(here("Data", "FittedModels", "Caret", "xgBoost", "OctApr_xgBModel.rds"))
# xgb_model2 <- readRDS(here("Data", "FittedModels", "Caret", "xgBoost", "OctMay_xgBModel.rds"))
# xgb_model3 <- readRDS(here("Data", "FittedModels", "Caret", "xgBoost", "SepMay_xgBModel.rds"))
# xgb_model4 <- readRDS(here("Data", "FittedModels", "Caret", "xgBoost", "OctAprDecFeb_xgBModel.rds"))
# xgb_model5 <- readRDS(here("Data", "FittedModels", "Caret", "xgBoost", "OctMayDecFeb_xgBModel.rds"))
# xgb_model6 <- readRDS(here("Data", "FittedModels", "Caret", "xgBoost", "SepMayDecFeb_xgBModel.rds"))
# 
# ### Display the best parameters and results for best models
# ### Model 1
# print(xgb_model1$bestTune)
# xgb_model1_BestTune_results <- xgb_model1$results |>
#   filter(RMSE == min(RMSE, na.rm = TRUE))
# xgb_model1_BestTune_results$RMSE
# xgb_model1_BestTune_results$Rsquared
# set.seed(802)
# xgb_model1_preds <- predict(xgb_model1, test_OctApr_x)
# rmse(test_y, xgb_model1_preds)
# ### Model 2
# print(xgb_model2$bestTune)
# xgb_model2_BestTune_results <- xgb_model2$results |>
#   filter(RMSE == min(RMSE, na.rm = TRUE))
# xgb_model2_BestTune_results$RMSE
# xgb_model2_BestTune_results$Rsquared
# set.seed(802)
# xgb_model2_preds <- predict(xgb_model2, test_OctMay_x)
# rmse(test_y, xgb_model2_preds)
# ### Model 3
# print(xgb_model3$bestTune)
# xgb_model3_BestTune_results <- xgb_model3$results |>
#   filter(RMSE == min(RMSE, na.rm = TRUE))
# xgb_model3_BestTune_results$RMSE
# xgb_model3_BestTune_results$Rsquared
# set.seed(802)
# xgb_model3_preds <- predict(xgb_model3, test_SepMay_x)
# rmse(test_y, xgb_model3_preds)
# ### Model 4
# print(xgb_model4$bestTune)
# xgb_model4_BestTune_results <- xgb_model4$results |>
#   filter(RMSE == min(RMSE, na.rm = TRUE))
# xgb_model4_BestTune_results$RMSE
# xgb_model4_BestTune_results$Rsquared
# set.seed(802)
# xgb_model4_preds <- predict(xgb_model4, test_OctAprDecFeb_x)
# rmse(test_y, xgb_model4_preds)
# ### Model 5
# print(xgb_model5$bestTune)
# xgb_model5_BestTune_results <- xgb_model5$results |>
#   filter(RMSE == min(RMSE, na.rm = TRUE))
# xgb_model5_BestTune_results$RMSE
# xgb_model5_BestTune_results$Rsquared
# set.seed(802)
# xgb_model5_preds <- predict(xgb_model5, test_OctMayDecFeb_x)
# rmse(test_y, xgb_model5_preds)
# ### Model 6
# print(xgb_model6$bestTune)
# xgb_model6_BestTune_results <- xgb_model6$results |>
#   filter(RMSE == min(RMSE, na.rm = TRUE))
# xgb_model6_BestTune_results$RMSE
# xgb_model6_BestTune_results$Rsquared
# set.seed(802)
# xgb_model6_preds <- predict(xgb_model6, test_SepMayDecFeb_x)
# rmse(test_y, xgb_model6_preds)




