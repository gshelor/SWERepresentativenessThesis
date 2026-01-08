##### Making Spatial Predictions with Best Random Forest Model #####
starttime <- Sys.time()
##### Loading packages, Reading in model/data #####
library(pacman)
p_load(here, caret, randomForest, tidyverse, terra, sf, tidyterra, RColorBrewer, viridis, ModelMetrics, ggthemes, ggspatial, mcprogress, data.table, parallel)

### CONUS AOI
CONUS_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "CONUS", "CONUSEcoregionsNoStates.gpkg"))
### Alaska AOI
# AK_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "Alaska", "AK_StateBoundary.gpkg"))

### reading in topography rasters since they are considered constant throughout the study period
## CONUS
CONUSDEM_rast <- rast(here("Data", "DEM", "CONUS", "CONUSDEMMosaic.tif"))
CONUSSlope_rast <- rast(here("Data", "DEM", "CONUS", "CONUSSlope.tif"))
CONUSAspect_rast <- rast(here("Data", "DEM", "CONUS", "CONUSAspect.tif"))
## Alaska
# AKDEM_rast <- rast(here("Data", "DEM", "Alaska", "AKDEMMosaic.tif"))
# AKSlope_rast <- rast(here("Data", "DEM", "Alaska", "AKSlope.tif"))
# AKAspect_rast <- rast(here("Data", "DEM", "Alaska", "AKAspect.tif"))

### loading model (pre-saved)
CONUSRFModel <- readRDS(here("Data", "FittedModels", "Caret", "RF", "AnnualCONUSRFModel_oob.rds"))
# AKRFModel <- readRDS(here("Data", "FittedModels", "Caret", "RF", "AnnualAKRFModel.rds"))

##### Function which will make predictions for each year, for both CONUS and Alaska #####
### storing vector of years
years <- 1993:2020
### making predictions
RFPredict <- function(year){
  # print(year)
  ### reading in rasters dependent on year
  prcpSum_OctAprCDMSum_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpSum_OctApr", year, "CDMSum.tif")))
  prcpMean_OctAprMean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpMean_OctApr", year, "CDMMean.tif")))
  tmean_OctAprCDMSum_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmean_OctApr", year, "CDMSum.tif")))
  tmin_OctAprMean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmin_OctApr", year, "Mean.tif")))
  tmax_OctAprMean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmax_OctApr", year, "Mean.tif")))
  srad_OctAprMean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("srad_OctApr", year, "Mean.tif")))
  Landcover_rast <- as.int(rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", paste0("LC_CONUS_", year, ".tif"))))
  
  ### Alaska
  ### Alaska variables
  ## c(OctMay_prcpSumCDMSum, OctMay_tmeanCDMSum, DecFeb_prcpSumCDMSum, DecFeb_tmeanCDMSum, elevation, slope, aspect, OctMay_tminMean, OctMay_tmaxMean, OctMay_sradMean, landcover_triclass)
  # AKprcpSum_OctMayCDMSum_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpSum_OctMay", year, "CDMSum.tif")))
  # AKprcpSum_DecFebCDMSum_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpSum_DecFeb", year, "CDMSum.tif")))
  # # AKprcpMean_OctMayMean_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpMean_OctMay", year, "CDMMean.tif")))
  # AKtmean_OctMayCDMSum_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmean_OctMay", year, "CDMSum.tif")))
  # AKtmean_DecFebCDMSum_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmean_DecFeb", year, "CDMSum.tif")))
  # AKtmin_OctMayMean_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmin_OctMay", year, "Mean.tif")))
  # AKtmax_OctMayMean_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmax_OctMay", year, "Mean.tif")))
  # # AKsrad_OctMayMean_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("srad_OctMay", year, "Mean.tif")))
  # AKLandcover_rast <- as.int(rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", paste0("LC_AK_", year, ".tif"))))
  
  
  ### Stacking rasters for predictions
  ## CONUS
  # print("stacking CONUS rasters")
  CONUSPredictors_rasts <- c(prcpSum_OctAprCDMSum_rast, tmean_OctAprCDMSum_rast, CONUSDEM_rast, CONUSSlope_rast, CONUSAspect_rast)
  ### filling NAs from numeric rasters with dummy value so the prediction runs
  ## they're all getting cropped/masked out at the end anyway so it shouldn't change predictions
  CONUSPredictors_rasts <- subst(CONUSPredictors_rasts, NA, 0)
  ### filling NAs with 0 for landcover because the substituted values need to be one of the factors (0,1,2)
  Landcover_rast <- as.int(subst(Landcover_rast, NA, 0))
  CONUSPredictors_rasts <- c(CONUSPredictors_rasts, Landcover_rast)
  names(CONUSPredictors_rasts) <- c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "landcover_triclass")
  
  ## Alaska
  # print("stacking Alaska rasters")
  # AKPredictors_rasts <- c(AKprcpSum_OctMayCDMSum_rast, AKtmean_OctMayCDMSum_rast,  AKprcpSum_DecFebCDMSum_rast,  AKtmean_DecFebCDMSum_rast, AKDEM_rast, AKSlope_rast, AKAspect_rast, AKtmin_OctMayMean_rast, AKtmax_OctMayMean_rast)
  # ### filling NAs from numeric rasters with dummy value so the prediction runs
  # ## they're all getting cropped/masked out at the end anyway so it shouldn't change predictions
  # AKPredictors_rasts <- subst(AKPredictors_rasts, NA, 0)
  # ### filling NAs with 0 for landcover because the substituted values need to be one of the factors (0,1,2)
  # AKLandcover_rast <- as.int(subst(AKLandcover_rast, NA, 0))
  # AKPredictors_rasts <- c(AKPredictors_rasts, AKLandcover_rast)
  # names(AKPredictors_rasts) <- c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "aspect", "OctMay_tminMean", "OctMay_tmaxMean", "landcover_triclass")
  
  
  ### Making predictions
  # print("Making Predictions for CONUS")
  set.seed(802)
  CONUSRF_PredictionMap <- terra::predict(CONUSPredictors_rasts, CONUSRFModel)
  ### cropping model because it predicts outside the land parts
  CONUSRF_PredictionMap_cropmask <- crop(CONUSRF_PredictionMap, CONUS_AOI, mask = TRUE)
  ### ensuring no values below 0
  CONUSRF_PredictionMap_cropmask <- app(CONUSRF_PredictionMap_cropmask, fun = function(x){x[x < 0] <- 0; return(x)})
  ### writing raster out as tif
  writeRaster(CONUSRF_PredictionMap_cropmask, here("Outputs", "MapTIFFs", "CONUS", paste0("CONUS_SWEmax", year, ".tif")), overwrite = TRUE)
  ### converting rasters to dataframe
  CONUSRF_PredictionMap_df <- as.data.frame(CONUSRF_PredictionMap_cropmask, xy = TRUE)
  fwrite(CONUSRF_PredictionMap_df, here("Data", "SHAP", "PredictionCSVs", "CONUS", paste0("CONUS_SWEmaxPredictions", year, ".csv")))
  CONUSRF_Predictors_df <- as.data.frame(CONUSPredictors_rasts, xy = TRUE) |>
    filter(paste0(x, y) %in% paste0(CONUSRF_PredictionMap_df$x, CONUSRF_PredictionMap_df$y))
  fwrite(CONUSRF_Predictors_df, file = here("Data", "SHAP", "PredictorCSVs", "CONUS", paste0("CONUS_SWEmaxPredictors", year, ".csv")))

  ##### Making ggplots of Predicted Peak SWE #####
  # print("making ggplot for CONUS")
  CONUSSWE_plot <- ggplot() +
    theme_bw() +
    geom_spatraster(data = CONUSRF_PredictionMap_cropmask) +
    scale_fill_whitebox_c(palette = "deep", direction = 1, name = "Peak SWE (mm)") +
    geom_sf(data = CONUS_AOI, fill = NA, color = "black") +
    ggtitle(label = paste0("Predicted Peak SWE for Water Year ", year)) +
    ggspatial::annotation_scale(location = "br") +
    ggspatial::annotation_north_arrow(location = "bl", height = unit(1.25, "cm"), width = unit(1.25, "cm")) +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
  CONUSSWE_plot
  ### saving ggplot
  ggsave(here("Outputs", "Plots", "SWEMaxMaps", "CONUS", "1Km", paste0("CONUSmSWE", year, ".png")))

  ### removing title for plots going into figures
  CONUSSWE_plot2 <- ggplot() +
    theme_bw() +
    geom_spatraster(data = CONUSRF_PredictionMap_cropmask) +
    scale_fill_whitebox_c(palette = "deep", direction = 1, name = "Peak SWE (mm)") +
    geom_sf(data = CONUS_AOI, fill = NA, color = "black") +
    # ggtitle(label = paste0("Predicted Peak SWE for Water Year ", year)) +
    ggspatial::annotation_scale(location = "br") +
    ggspatial::annotation_north_arrow(location = "bl", height = unit(1.25, "cm"), width = unit(1.25, "cm")) +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
  CONUSSWE_plot2
  ### saving ggplot
  ggsave(here("Outputs", "Plots", "SWEMaxMaps", "PaperFigures", "CONUS", paste0("CONUSmSWE", year, ".png")))
  
  ### Making Alaska predictions
  # set.seed(802)
  # # print("Making Predictions for Alaska")
  # AKRF_PredictionMap <- terra::predict(AKPredictors_rasts, AKRFModel)
  # ### cropping model because it predicts outside the land parts
  # AKRF_PredictionMap_cropmask <- crop(AKRF_PredictionMap, AK_AOI, mask = TRUE)
  # ### ensuring no values below 0
  # AKRF_PredictionMap_cropmask <- app(AKRF_PredictionMap_cropmask, fun = function(x){x[x < 0] <- 0; return(x)})
  # ### writing raster out as tif
  # writeRaster(AKRF_PredictionMap_cropmask, here("Outputs", "MapTIFFs", "Alaska", paste0("AK_SWEmax", year, ".tif")), overwrite = TRUE)
  # ### converting rasters to dataframe
  # AKRF_Predictors_df <- as.data.frame(AKPredictors_rasts, xy = TRUE)
  # fwrite(AKRF_Predictors_df, here("Data", "SHAP", "PredictorCSVs", "Alaska", paste0("AK_SWEmaxPredictors", year, ".csv")))
  # AKRF_PredictionMap_df <- as.data.frame(AKRF_PredictionMap_cropmask, xy = TRUE)
  # fwrite(AKRF_PredictionMap_df, here("Data", "SHAP", "PredictionCSVs", "Alaska", paste0("AK_SWEmaxPredictions", year, ".csv")))
  # AKRF_Predictors_df <- as.data.frame(AKPredictors_rasts, xy = TRUE) |>
  #   filter(paste0(x, y) %in% paste0(AKRF_PredictionMap_df$x, AKRF_PredictionMap_df$y))
  # fwrite(AKRF_Predictors_df, here("Data", "SHAP", "PredictorCSVs", "Alaska", paste0("AK_SWEmaxPredictors", year, ".csv")))
  # 
  # ### Alaska
  # # print("making ggplot for Alaska")
  # AKSWE_plot <- ggplot() +
  #   theme_bw() +
  #   geom_spatraster(data = AKRF_PredictionMap_cropmask) +
  #   scale_fill_whitebox_c(palette = "deep", direction = 1, name = "Peak SWE (mm)") +
  #   geom_sf(data = AK_AOI, fill = NA, color = "black") +
  #   ggtitle(label = paste0("Predicted Peak SWE for Water Year ", year)) +
  #   ggspatial::annotation_scale(location = "br") +
  #   ggspatial::annotation_north_arrow(location = "tr", height = unit(1.25, "cm"), width = unit(1.25, "cm")) +
  #   theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
  # AKSWE_plot
  # ### saving ggplot
  # ggsave(here("Outputs", "Plots", "SWEMaxMaps", "Alaska", "1Km", paste0("AKmSWE", year, "_1km.png")))
  # 
  # ### making plot without title so it will be better suited for publication (captions used in publication to describe/title plot)
  # AKSWE_plot2 <- ggplot() +
  #   theme_bw() +
  #   geom_spatraster(data = AKRF_PredictionMap_cropmask) +
  #   scale_fill_whitebox_c(palette = "deep", direction = 1, name = "Peak SWE (mm)") +
  #   geom_sf(data = AK_AOI, fill = NA, color = "black") +
  #   # ggtitle(label = paste0("Predicted Peak SWE for Water Year ", year)) +
  #   ggspatial::annotation_scale(location = "br") +
  #   ggspatial::annotation_north_arrow(location = "tr", height = unit(1.25, "cm"), width = unit(1.25, "cm")) +
  #   theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
  # AKSWE_plot2
  # ### saving ggplot
  # ggsave(here("Outputs", "Plots", "SWEMaxMaps", "PaperFigures", "Alaska", paste0("AKmSWE", year, "_1km.png")))
}

### making predictions on all years using parallel package
## change value of mc.cores if running on personal computer
# Predictions <- pmclapply(years, RFPredict, mc.cores = ifelse(length(years) > detectCores(), ceiling(detectCores() / 1.5), length(years)), mc.silent = FALSE)
Predictions <- pmclapply(years, RFPredict, mc.cores = 10, mc.silent = FALSE)

endtime <- Sys.time()
totaltime <- endtime - starttime
totaltime
