##### SHAP #####
totalstarttime <- Sys.time()
##### Loading packages #####
library(pacman)
p_load(tidyverse, randomForest, sf, terra, sfext, caret, parallel, data.table, rsample, mcprogress, fastshap, shapviz)
# options(mc.cores = parallel::detectCores())
setwd("/media/Research/Morafkai/GriffinS")

##### reading in data #####
AKRFModel <- read_rds("FittedModels/AnnualAKRFModel.rds")
### CONUS AOI
# CONUS_AOI <- read_sf("SHAPAnalysis/Data/L3_Ecoregions_USB/CONUS/CONUS_AOI.gpkg")
### Alaska AOI
AK_AOI <- read_sf("SHAPAnalysis/Data/L3_Ecoregions_USB/Alaska/AK_AOI.gpkg")
### SNOTEL Data used to fit the model
SnotelAK_sf <- read_sf("Combined/GIS/SnotelAK_PeakSWE_Covars.gpkg") |>
  drop_na(OctMay_prcpSumCDMSum, OctMay_tmeanCDMSum, DecFeb_prcpSumCDMSum, DecFeb_tmeanCDMSum, elevation, slope, aspect, OctMay_tminMean, OctMay_tmaxMean, landcover_triclass)


### converting to dfs for model fitting process
SnotelAK_df <- sf_to_df(SnotelAK_sf) |>
  drop_na(OctMay_prcpSumCDMSum, OctMay_tmeanCDMSum, DecFeb_prcpSumCDMSum, DecFeb_tmeanCDMSum, elevation, slope, aspect, OctMay_tminMean, OctMay_tmaxMean, landcover_triclass)
### making sure landcover is type integer
# SnotelAK_df$landcover <- as.integer(SnotelAK_df$landcover)
SnotelAK_df$landcover_triclass <- as.integer(SnotelAK_df$landcover_triclass)

##### splitting data into training and testing datasets #####
set.seed(802)
SnotelAK_split <- initial_split(SnotelAK_df, prop = 0.75, strata = peak_swe)
SnotelAK_train <- training(SnotelAK_split)
SnotelAK_test <- testing(SnotelAK_split)

### OctMay aggregates and DecFeb aggregates for prcpSum and tmean, plus tmin, tmax
AKTrain_OctMayDecFebAspect_Nosrad_x <- SnotelAK_train[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "aspect", "OctMay_tminMean", "OctMay_tmaxMean", "landcover_triclass")]
AKTest_OctMayDecFebAspect_Nosrad_x <- SnotelAK_test[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "aspect", "OctMay_tminMean", "OctMay_tmaxMean", "landcover_triclass")]
trainAK_y <- SnotelAK_train$peak_swe
testAK_y <- SnotelAK_test$peak_swe

### saving variable names to be listed in ranking order later
AKvars <- c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "aspect", "OctMay_tminMean", "OctMay_tmaxMean", "landcover_triclass")
### creating folders to save feature dependence plots if they don't exist
for (x in AKvars){
  if (dir.exists(paste0("SHAPAnalysis/Outputs/Plots/Fastshap/FeatureDependence/Alaska/", x)) == FALSE){
    dir.create(paste0("SHAPAnalysis/Outputs/Plots/Fastshap/FeatureDependence/Alaska/", x))
  }
}


##### SHAP stuff #####
### in order to run treeshap, the model must be "unified"
# ?unify
# UnifiedRF <- unify(AKRFModel$finalModel, data = train_OctAprAspect_Nosradtmintmax_x)

### years of study period
years <- 1993:2020
### function to specify prediction method for fastshap explain()
pfun <- function(object, newdata) {  # needs to return a numeric vector
  predict(object, newdata)  
}
##### Alaska Function to call fastshap and create SHAP-based classes #####
AKSHAPFunc <- function(year){
  # temp_starttime <- Sys.time()
  # print(temp_starttime)
  
  ##### Alaska #####
  # print("Alaska")
  ### reading in csv versions of data frames so Fastshap can understand them
  ### SWE_max predictions
  # print("reading in prediction csv")
  AKPredictions <- fread(file = paste0("SHAPAnalysis/Data/SHAP/PredictionCSVs/Alaska/", paste0("AK_SWEmaxPredictions", year, ".csv")), data.table = FALSE)
  ### removing xy
  # AKPredictions_NoXY <- AKPredictions[,3]
  colnames(AKPredictions) <- c("x", "y", "peak_swe")
  
  # print("reading in predictor csv")
  AKPredictors <- fread(file = paste0("SHAPAnalysis/Data/SHAP/PredictorCSVs/Alaska/", paste0("AK_SWEmaxPredictors", year, ".csv")), data.table = FALSE) |>
    select(OctMay_prcpSumCDMSum, OctMay_tmeanCDMSum, DecFeb_prcpSumCDMSum, DecFeb_tmeanCDMSum, elevation, slope, aspect, OctMay_tminMean, OctMay_tmaxMean, landcover_triclass) #|>
  # filter(paste0(x,y) %in% paste0(AKPredictions$x, AKPredictions$y))
  ### removing xy data for Fastshap
  # AKPredictors_NoXY <- AKPredictors |>
  #   select(OctApr_prcpSumCDMSum, OctApr_tmeanCDMSum, elevation, slope, aspect, landcover_triclass)
  ### storing xy separately
  # AKPredictors_XY <- AKPredictors |>
  #   select(x,y)
  
  ### binding predictors with associated predictions
  # print("binding predictors with associated predictions")
  # AKVals <- cbind(AKPredictors_NoXY, AKPredictions_NoXY)
  
  ### Running Fastshap
  # print("running Fastshap")
  set.seed(802)
  AK_Fastshap <- explain(AKRFModel$finalModel, X = AKTrain_OctMayDecFebAspect_Nosrad_x, nsim = 100, pred_wrapper = pfun, newdata = as.data.frame(AKPredictors), adjust = TRUE)
  AK_Fastshap_df <- as.data.frame(AK_Fastshap)
  write_rds(AK_Fastshap, paste0("SHAPAnalysis/Data/SHAP/Fastshap/Alaska/", paste0("AKSWEmax_SHAPs", year, ".rds")), compress = "gz")
  
  ### making plots of feature dependence
  ## only making them and saving them, not displaying because it takes forever
  sv <- shapviz(object = AK_Fastshap, X = AKPredictors)
  plot1 <- sv_dependence(sv, v = AKvars[1], color_var = AKvars[1], jitter_width = 0) +
    ggtitle("Oct-May Precipitation CDM Sum Feature Dependence", subtitle =
              paste("Alaska", as.character(year))) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
    xlab(label = "Oct-May Total Precipitation CDM Sum (mm)") +
    ylab(label = "SHAP Value")
  plot2 <- sv_dependence(sv, v = AKvars[2], color_var = AKvars[2], jitter_width = 0) +
    ggtitle("Oct-May Mean Temperature CDM Sum Feature Dependence", subtitle = paste("Alaska", as.character(year))) + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
    xlab(label = expression("Oct-May Mean Temp CDM Sum ("*degree*C*")")) +
    ylab(label = "SHAP Value")
  plot3 <- sv_dependence(sv, v = AKvars[3], color_var = AKvars[3], jitter_width = 0) +
    ggtitle("Dec-Feb Precipitation CDM Sum Feature Dependence", subtitle =
              paste("Alaska", as.character(year))) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
    xlab(label = "Dec-Feb Total Precipitation CDM Sum (mm)") +
    ylab(label = "SHAP Value")
  plot4 <- sv_dependence(sv, v = AKvars[4], color_var = AKvars[4], jitter_width = 0) +
    ggtitle("Dec-Feb Mean Temperature CDM Sum Feature Dependence", subtitle = paste("Alaska", as.character(year))) + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
    xlab(label = expression("Dec-Feb Mean Temp CDM Sum ("*degree*C*")")) +
    ylab(label = "SHAP Value")
  plot5 <- sv_dependence(sv, v = AKvars[5], color_var = AKvars[5], jitter_width = 0) +
    ggtitle("Elevation Feature Dependence", subtitle = paste("Alaska", as.character(year))) + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
    xlab(label = "Elevation (m)") +
    ylab(label = "SHAP Value")
  plot6 <- sv_dependence(sv, v = AKvars[6], color_var = AKvars[6], jitter_width = 0) +
    ggtitle("Slope Feature Dependence", subtitle = paste("Alaska", as.character(year))) + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
    xlab(label = "Slope") +
    ylab(label = "SHAP Value")
  plot7 <- sv_dependence(sv, v = AKvars[7], color_var = AKvars[7], jitter_width = 0) +
    ggtitle("Aspect Feature Dependence", subtitle = paste("Alaska", as.character(year))) + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
    xlab(label = "Aspect") +
    ylab(label = "SHAP Value")
  plot8 <- sv_dependence(sv, v = AKvars[8], color_var = AKvars[8], jitter_width = 0) +
    ggtitle("Oct-May Mean Minimum Temperature Feature Dependence", subtitle = paste("Alaska", as.character(year))) + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
    xlab(label = expression("Oct-May Mean Min Temp ("*degree*C*")")) +
    ylab(label = "SHAP Value")
  plot9 <- sv_dependence(sv, v = AKvars[9], color_var = AKvars[9], jitter_width = 0) +
    ggtitle("Oct-May Mean Maximum Temperature Feature Dependence", subtitle = paste("Alaska", as.character(year))) + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
    xlab(label = expression("Oct-May Mean Max Temp ("*degree*C*")")) +
    ylab(label = "SHAP Value")
  plot10 <- sv_dependence(sv, v = AKvars[10], color_var = AKvars[10], jitter_width = 0) +
    ggtitle("Landcover Feature Dependence", subtitle = paste("Alaska", as.character(year))) + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
    xlab(label = "Landcover") +
    ylab(label = "SHAP Value")
  
  ### Saving Feature Dependence plots
  ggsave(paste0("SHAPAnalysis/Outputs/Plots/Fastshap/FeatureDependence/Alaska/", AKvars[1], "/", paste0("Alaska", year, "_", AKvars[1], "_FeatureDependence.png")), plot = plot1)
  ggsave(paste0("SHAPAnalysis/Outputs/Plots/Fastshap/FeatureDependence/Alaska/", AKvars[2], "/", paste0("Alaska", year, "_", AKvars[2], "_FeatureDependence.png")), plot = plot2)
  ggsave(paste0("SHAPAnalysis/Outputs/Plots/Fastshap/FeatureDependence/Alaska/", AKvars[3], "/", paste0("Alaska", year, "_", AKvars[3], "_FeatureDependence.png")), plot = plot3)
  ggsave(paste0("SHAPAnalysis/Outputs/Plots/Fastshap/FeatureDependence/Alaska/", AKvars[4], "/", paste0("Alaska", year, "_", AKvars[4], "_FeatureDependence.png")), plot = plot4)
  ggsave(paste0("SHAPAnalysis/Outputs/Plots/Fastshap/FeatureDependence/Alaska/", AKvars[5], "/", paste0("Alaska", year, "_", AKvars[5], "_FeatureDependence.png")), plot = plot5)
  ggsave(paste0("SHAPAnalysis/Outputs/Plots/Fastshap/FeatureDependence/Alaska/", AKvars[6], "/", paste0("Alaska", year, "_", AKvars[6], "_FeatureDependence.png")), plot = plot6)
  ggsave(paste0("SHAPAnalysis/Outputs/Plots/Fastshap/FeatureDependence/Alaska/", AKvars[7], "/", paste0("Alaska", year, "_", AKvars[7], "_FeatureDependence.png")), plot = plot7)
  ggsave(paste0("SHAPAnalysis/Outputs/Plots/Fastshap/FeatureDependence/Alaska/", AKvars[8], "/", paste0("Alaska", year, "_", AKvars[8], "_FeatureDependence.png")), plot = plot8)
  ggsave(paste0("SHAPAnalysis/Outputs/Plots/Fastshap/FeatureDependence/Alaska/", AKvars[9], "/", paste0("Alaska", year, "_", AKvars[9], "_FeatureDependence.png")), plot = plot9)
  ggsave(paste0("SHAPAnalysis/Outputs/Plots/Fastshap/FeatureDependence/Alaska/", AKvars[10], "/", paste0("Alaska", year, "_", AKvars[10], "_FeatureDependence.png")), plot = plot10)
  
  
  ### Extracting SHAP values, adding dummy rank columns which will be filled in later
  AK_SHAPs <- AK_Fastshap_df |>
    mutate(Rank_OctMayprcpSumCDM_varimp = -999,
           Rank_OctMaytmeanCDM_varimp = -999,
           Rank_DecFebprcpSumCDM_varimp = -999,
           Rank_DecFebtmeanCDM_varimp = -999,
           Rank_elevation_varimp = -999,
           Rank_slope_varimp = -999,
           Rank_aspect_varimp = -999,
           Rank_OctMaytminMean_varimp = -999,
           Rank_OctMaytmaxMean_varimp = -999,
           Rank_landcover_varimp = -999,
           class = "dummy string",
           class_top1 = "dummy string",
           class_top2 = "dummy string",
           class_top3 = "dummy string",
           class_top4 = "dummy string",
           class_top5 = "dummy string")
  
  
  ### extracting only columns containing shap values
  shap_values <- AK_SHAPs[, c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "aspect", "OctMay_tminMean", "OctMay_tmaxMean", "landcover_triclass")]
  # AKvars <- names(shap_values)
  
  ### ranking SHAP values by absolute value
  ranked_data <- t(apply(shap_values, 1, function(row) {
    abs_vals <- abs(row)
    ranks <- row_number(desc(abs_vals)) # Rank in descending order of absolute value
    ranked_vars <- vars[order(ranks)]
    class_all <- paste(ranked_vars, collapse = ", ")
    class_top1 <- paste(ranked_vars[1], collapse = ", ")
    class_top2 <- paste(ranked_vars[1:2], collapse = ", ")
    class_top3 <- paste(ranked_vars[1:3], collapse = ", ")
    class_top4 <- paste(ranked_vars[1:4], collapse = ", ")
    class_top5 <- paste(ranked_vars[1:5], collapse = ", ")
    
    return(c(ranks, class_all, class_top1, class_top2, class_top3, class_top4, class_top5))
  }))
  
  ### Assign the ranks back to the original data frame
  AK_SHAPs$Rank_OctMayprcpSumCDM_varimp <- as.numeric(ranked_data[, 1])
  AK_SHAPs$Rank_OctMaytmeanCDM_varimp <- as.numeric(ranked_data[, 2])
  AK_SHAPs$Rank_DecFebprcpSumCDM_varimp <- as.numeric(ranked_data[, 3])
  AK_SHAPs$Rank_DecFebtmeanCDM_varimp <- as.numeric(ranked_data[, 4])
  AK_SHAPs$Rank_elevation_varimp <- as.numeric(ranked_data[, 5])
  AK_SHAPs$Rank_slope_varimp <- as.numeric(ranked_data[, 6])
  AK_SHAPs$Rank_aspect_varimp <- as.numeric(ranked_data[, 7])
  AK_SHAPs$Rank_OctMaytminMean_varimp <- as.numeric(ranked_data[, 8])
  AK_SHAPs$Rank_OctMaytmaxMean_varimp <- as.numeric(ranked_data[, 9])
  AK_SHAPs$Rank_landcover_varimp <- as.numeric(ranked_data[, 10])
  AK_SHAPs$class <- ranked_data[, 11]
  AK_SHAPs$class_top1 <- ranked_data[, 12]
  AK_SHAPs$class_top2 <- ranked_data[, 13]
  AK_SHAPs$class_top3 <- ranked_data[, 14]
  AK_SHAPs$class_top4 <- ranked_data[, 15]
  AK_SHAPs$class_top5 <- ranked_data[, 16]
  
  ### binding Shapley values and initial classes to predicted peak SWE df to be rasterized by class below
  AK_SHAPs <- cbind(AKPredictions, AK_SHAPs)
  ### saving SHAPs
  fwrite(AK_SHAPs, file = paste0("SHAPAnalysis/Data/SHAP/Fastshap/SHAPs/Alaska/", paste0("Alaska", year, "SHAPs.csv")))
  
  ### reading in prediction raster to use as basis for rasterization
  PredictionRast <- rast(paste0("SHAPAnalysis/Outputs/MapTIFFs/Alaska/", paste0("AK_SWEmax", year, ".tif")))
  AK_SHAPClasses <- terra::rasterize(vect(st_as_sf(AK_SHAPs, coords = c("x", "y"), crs = crs(AK_AOI))), y = PredictionRast, field = "class")
  ### top 1
  AK_SHAPClasses_top1 <- terra::rasterize(vect(st_as_sf(AK_SHAPs, coords = c("x", "y"), crs = crs(AK_AOI))), y = PredictionRast, field = "class_top1")
  ### top 2
  AK_SHAPClasses_top2 <- terra::rasterize(vect(st_as_sf(AK_SHAPs, coords = c("x", "y"), crs = crs(AK_AOI))), y = PredictionRast, field = "class_top2")
  ### top 3
  AK_SHAPClasses_top3 <- terra::rasterize(vect(st_as_sf(AK_SHAPs, coords = c("x", "y"), crs = crs(AK_AOI))), y = PredictionRast, field = "class_top3")
  ### top 4
  AK_SHAPClasses_top4 <- terra::rasterize(vect(st_as_sf(AK_SHAPs, coords = c("x", "y"), crs = crs(AK_AOI))), y = PredictionRast, field = "class_top4")
  ### top 5
  AK_SHAPClasses_top5 <- terra::rasterize(vect(st_as_sf(AK_SHAPs, coords = c("x", "y"), crs = crs(AK_AOI))), y = PredictionRast, field = "class_top5")
  
  ### Saving Rasters
  writeRaster(AK_SHAPClasses, file = paste0(getwd(), "/SHAPAnalysis/Data/SHAP/Classes/Alaska/", paste0("Alaska", year, "SHAPClasses.tif")) , overwrite = TRUE)
  writeRaster(AK_SHAPClasses_top1, file = paste0(getwd(), "/SHAPAnalysis/Data/SHAP/Classes/Alaska/", paste0("Alaska", year, "SHAPClasses_top1.tif")) , overwrite = TRUE)
  writeRaster(AK_SHAPClasses_top2, file = paste0(getwd(), "/SHAPAnalysis/Data/SHAP/Classes/Alaska/", paste0("Alaska", year, "SHAPClasses_top2.tif")) , overwrite = TRUE)
  writeRaster(AK_SHAPClasses_top3, file = paste0(getwd(), "/SHAPAnalysis/Data/SHAP/Classes/Alaska/", paste0("Alaska", year, "SHAPClasses_top3.tif")) , overwrite = TRUE)
  writeRaster(AK_SHAPClasses_top4, file = paste0(getwd(), "/SHAPAnalysis/Data/SHAP/Classes/Alaska/", paste0("Alaska", year, "SHAPClasses_top4.tif")) , overwrite = TRUE)
  writeRaster(AK_SHAPClasses_top5, file = paste0(getwd(), "/SHAPAnalysis/Data/SHAP/Classes/Alaska/", paste0("Alaska", year, "SHAPClasses_top5.tif")) , overwrite = TRUE)
  
  # temp_endtime <- Sys.time()
  # print(temp_endtime)
  # temp_endtime - temp_starttime
}

### calling function to run treeshap on all years
## change number of cores if running on personal computer
# AKSWE_SHAPs <- pmclapply(X = years, FUN = AKSHAPFunc, mc.cores = ifelse(length(years) > detectCores(), ceiling(detectCores() / 1.5), length(years)), mc.silent = FALSE, mc.set.seed = FALSE)
AKSWE_SHAPs <- pmclapply(X = years, FUN = AKSHAPFunc, mc.cores = 14, mc.silent = FALSE, mc.set.seed = FALSE)

# SWE_SHAPs[[1]]
# str(SWE_SHAPs[[1]])



totalendtime <- Sys.time()
totaltime <- totalendtime - totalstarttime
print(totaltime)
