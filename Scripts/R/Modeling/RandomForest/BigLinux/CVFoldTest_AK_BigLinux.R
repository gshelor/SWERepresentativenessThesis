##### testing a number of folds for my RF model in caret #####
##### Loading packages #####
library(pacman)
p_load(tidyverse, randomForest, sf, terra, sfext, rsample, caret, parallel, grateful, ModelMetrics, Metrics, mcprogress)
# options(mc.cores = parallel::detectCores())
setwd("/media/Research/Morafkai/GriffinS")

##### reading in data #####
### Alaska AOI
# AK_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "Alaska", "AK_AOI.gpkg"))
### SNOTEL Annual values
SnotelAK_sf <- read_sf("Combined/GIS/SnotelAK_PeakSWE_Covars.gpkg") |>
  drop_na(OctApr_prcpSumCDMSum, OctApr_tmeanCDMSum, elevation, slope, aspect, landcover_triclass)

### converting to df for model fitting process
SnotelAK_df <- sf_to_df(SnotelAK_sf) |>
  drop_na(OctApr_prcpSumCDMSum, OctApr_tmeanCDMSum, elevation, slope, aspect, landcover_triclass)
### making sure landcover is type integer
SnotelAK_df$landcover <- as.integer(SnotelAK_df$landcover)
SnotelAK_df$landcover_triclass <- as.integer(SnotelAK_df$landcover_triclass)

##### splitting data into training and testing datasets #####
set.seed(802)
Snotel_split <- initial_split(SnotelAK_df, prop = 0.75)
Snotel_train <- training(Snotel_split)
Snotel_test <- testing(Snotel_split)
### keeping as dataframes for caret
### OctApr aggregates only
### OctApr aggregates only but with aspect, and srad, tmin, tmax removed
train_OctAprAspect_Nosradtmintmax_x <- Snotel_train[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "landcover_triclass")]
test_OctAprAspect_Nosradtmintmax_x <- Snotel_test[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "landcover_triclass")]
train_y <- Snotel_train$peak_swe
test_y <- Snotel_test$peak_swe



##### Tune grid final model #####
tune_grid_OctAprAspect_Nosradtmintmax_rf <- expand.grid(
  ### Tune mtry (number of variables randomly sampled)
  mtry = seq(2, ncol(train_OctAprAspect_Nosradtmintmax_x), by = 1),
  splitrule = "variance",
  min.node.size = seq(5, 25, by = 1)
)


cv_func <- function(fold){
  tune_control_rf <- trainControl(
    method = "cv", # out of bag
    number = fold, # 10-fold cross-validation
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
    x = train_OctAprAspect_Nosradtmintmax_x,
    y = train_y,
    method = "ranger",
    trControl = tune_control_rf,
    tuneGrid = tune_grid_OctAprAspect_Nosradtmintmax_rf,
    # importance = TRUE,
    importance = "impurity",
    num.trees = 2500
  )
  
  return(RF_model)
}




folds <- 2:50
FoldsTest <- pmclapply(folds, cv_func, mc.cores = length(folds), mc.silent = FALSE)

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
  temp_preds <- predict(temp_model, test_OctAprAspect_Nosradtmintmax_x)
  test_rmse_vals <- c(test_rmse_vals, rmse(test_y, temp_preds))
  test_mae_vals <- c(test_mae_vals, mae(test_y, temp_preds))
}

which.max(rsquared_vals)
which.min(rmse_vals)
which.min(mae_vals)
which.min(test_rmse_vals)
which.min(test_mae_vals)

rsquared_vals[which.max(rsquared_vals)]
rmse_vals[which.min(rmse_vals)]
mae_vals[which.min(mae_vals)]
rmse_vals[which.min(test_rmse_vals)]
### model metrics for model with best test RMSE
rsquared_vals[which.min(test_rmse_vals)]
rmse_vals[which.min(test_rmse_vals)]
test_rmse_vals[which.min(test_rmse_vals)]
test_mae_vals[which.min(test_rmse_vals)]
### model metrics for model with best test MAE
rsquared_vals[which.min(test_mae_vals)]
test_rmse_vals[which.min(test_mae_vals)]
test_mae_vals[which.min(test_mae_vals)]

hist(rsquared_vals)
# hist(rmse_vals)
hist(test_mae_vals)
hist(test_rmse_vals)

plot(folds, rsquared_vals)
plot(folds, test_rmse_vals)
plot(folds, test_mae_vals)


##### Final Model Extraction #####
FinalModel <- FoldsTest[[which.min(test_rmse_vals)]]
FinalModel
FinalModel$bestTune
FinalModel$results

write_rds(FinalModel, "FittedModels/AnnualAKRFModel.rds", compress = "gz")
