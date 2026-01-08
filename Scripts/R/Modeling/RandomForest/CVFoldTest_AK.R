##### testing a number of folds for my RF model in caret #####
##### Loading packages #####
library(pacman)
p_load(here, tidyverse, randomForest, sf, terra, sfext, rsample, caret, parallel, grateful, ModelMetrics, Metrics, mcprogress)
# options(mc.cores = parallel::detectCores())

##### reading in data #####
### Alaska AOI
AK_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "Alaska", "AK_AOI.gpkg"))
### SNOTEL Annual values
SnotelAK_sf <- read_sf(here("Data", "SNOTEL", "Alaska", "GIS", "GPKG", "Annual", "SnotelAK_PeakSWE_Covars.gpkg")) |>
  drop_na(OctMay_prcpSumCDMSum, OctMay_tmeanCDMSum, DecFeb_prcpSumCDMSum, DecFeb_tmeanCDMSum, elevation, slope, aspect, OctMay_tminMean, OctMay_tmaxMean, landcover_triclass)

### converting to df for model fitting process
SnotelAK_df <- sf_to_df(SnotelAK_sf) |>
  drop_na(OctMay_prcpSumCDMSum, OctMay_tmeanCDMSum, DecFeb_prcpSumCDMSum, DecFeb_tmeanCDMSum, elevation, slope, aspect, OctMay_tminMean, OctMay_tmaxMean, landcover_triclass)
### making sure landcover is type integer
SnotelAK_df$landcover <- as.integer(SnotelAK_df$landcover)
SnotelAK_df$landcover_triclass <- as.integer(SnotelAK_df$landcover_triclass)

##### splitting data into training and testing datasets #####
set.seed(802)
Snotel_split <- initial_split(SnotelAK_df, prop = 0.75, strata = peak_swe)
SnotelAK_train <- training(Snotel_split)
SnotelAK_test <- testing(Snotel_split)
### keeping as dataframes for caret
### OctApr aggregates only
### OctApr aggregates only but with aspect, and srad, tmin, tmax removed
# train_OctAprAspect_Nosradtmintmax_x <- SnotelAK_train[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "landcover_triclass")]
# test_OctAprAspect_Nosradtmintmax_x <- SnotelAK_test[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "landcover_triclass")]
### OctMay aggregates and DecFeb aggregates for prcpSum and tmean, plus tmin, tmax, srad
AKTrain_OctMayDecFebAspect_Nosrad_x <- SnotelAK_train[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "aspect", "OctMay_tminMean", "OctMay_tmaxMean", "landcover_triclass")]
AKTest_OctMayDecFebAspect_Nosrad_x <- SnotelAK_test[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "aspect", "OctMay_tminMean", "OctMay_tmaxMean", "landcover_triclass")]
train_y <- SnotelAK_train$peak_swe
test_y <- SnotelAK_test$peak_swe



##### Tune grid final model #####
tune_grid_OctMayDecFebAspect_Nosrad_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(AKTrain_OctMayDecFebAspect_Nosrad_x), by = 1) #,
  # splitrule = "variance",
  # min.node.size = seq(5, 25, by = 1)
)


cv_func <- function(fold){
  tune_control_rf <- trainControl(
    method = "cv", # out of bag
    number = fold, # k-fold cross-validation
    # repeats = 3, # Repeat cross-validation 3 times (more robust)
    verboseIter = TRUE, # Print progress during training
    returnData = FALSE, # Don't save the training data (saves memory)
    savePredictions = "none", # Save predictions from the best model
    returnResamp = "final", # Save resampling results from the best model
    allowParallel = TRUE, # Enable parallel processing (if available)
    predictionBounds = c(0, NA)
  )
  
  ##### Model 10, OctApr aggregates and aspect, no srad or tmin or tmax #####
  set.seed(802)
  RF_model <- train(
    x = AKTrain_OctMayDecFebAspect_Nosrad_x,
    y = train_y,
    method = "rf",
    trControl = tune_control_rf,
    tuneGrid = tune_grid_OctMayDecFebAspect_Nosrad_rf,
    importance = TRUE,
    # importance = "impurity",
    # num.trees = 2500
    ntree = 2500
  )
  return(RF_model)
}




folds <- 2:50
FoldsTest <- pmclapply(folds, cv_func, mc.cores = ifelse(length(folds) > detectCores(), ceiling(detectCores() / 1.5), length(folds)), mc.silent = FALSE)

### storing key metrics from each round of CV
rsquared_vals <- c()
rmse_vals <- c()
mae_vals <- c()
test_rmse_vals <- c()
test_mae_vals <- c()
for (x in 1:length(FoldsTest)){
  temp_model <- FoldsTest[[x]]
  temp_results <- temp_model$results
  rsquared_vals <- c(rsquared_vals, max(temp_results$Rsquared))
  rmse_vals <- c(rmse_vals, min(temp_results$RMSE))
  mae_vals <- c(mae_vals, min(temp_results$MAE))
  set.seed(802)
  temp_preds <- predict(temp_model, AKTest_OctMayDecFebAspect_Nosrad_x)
  test_rmse_vals <- c(test_rmse_vals, rmse(test_y, temp_preds))
  test_mae_vals <- c(test_mae_vals, mae(test_y, temp_preds))
}

print(which.max(rsquared_vals))
# which.min(rmse_vals)
# which.min(mae_vals)
print("number of folds for best test error metrics (RMSE, MAE)")
print(which.min(test_rmse_vals))
print(which.min(test_mae_vals))

### max R^2
print("max R^2")
rsquared_vals[which.max(rsquared_vals)]
# rmse_vals[which.min(rmse_vals)]
# mae_vals[which.min(mae_vals)]
# rmse_vals[which.min(test_rmse_vals)]
### model metrics for model with best test RMSE
print("model metrics for model with best test RMSE")
print(rsquared_vals[which.min(test_rmse_vals)])
print(test_rmse_vals[which.min(test_rmse_vals)])
print(test_mae_vals[which.min(test_rmse_vals)])
### model metrics for model with best test MAE
print("model metrics for model with best test MAE")
print(rsquared_vals[which.min(test_mae_vals)])
print(test_rmse_vals[which.min(test_mae_vals)])
print(test_mae_vals[which.min(test_mae_vals)])

hist(rsquared_vals)
# hist(rmse_vals)
hist(test_mae_vals)
hist(test_rmse_vals)

plot(folds, rsquared_vals)
plot(folds, test_rmse_vals)
plot(folds, test_mae_vals)


##### Final Model Extraction #####
print("FINAL MODEL")
FinalModel <- FoldsTest[[which.min(test_rmse_vals)]]
print(FinalModel)
print(FinalModel$bestTune)
print(FinalModel$results)
### model metrics for model with best test RMSE
print("model metrics for model with best test RMSE")
print(rsquared_vals[which.min(test_rmse_vals)])
print(test_rmse_vals[which.min(test_rmse_vals)])
print(test_mae_vals[which.min(test_rmse_vals)])

write_rds(FinalModel, here("Data", "FittedModels", "Caret", "RF", "AnnualAKRFModel.rds"), compress = "gz")
