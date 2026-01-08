##### SHAP Analysis #####
##### SCRIPT NO LONGER NECESSARY, FUNCTION INCLUDED IN SWE_SHAP.R #####
### loading packages, reading in data
library(pacman)
p_load(here, tidyverse, sf, terra, treeshap, parallel, mcprogress, data.table)


### setting numberyears I'm going to run my class creation function for
years <- 1993:2020

### Reading in AOIs to use to set CRS of SHAP-based classes
### CONUS AOI
CONUS_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "CONUS", "CONUS_AOI.gpkg"))
### Alaska AOI
AK_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "Alaska", "AK_StateBoundary.gpkg"))

##### Function to Define SWE-based Snow Classes #####
SWEclass_creation <- function(year){
  ### reading in CONUS treeshap object
  CONUS_Treeshap <- read_rds(here("Data", "SHAP", "Treeshap", "CONUS", paste0("CONUSSWEMax_SHAPs", year, ".rds")))
  
  ### making plots of feature dependence
  ## only making them and saving them, not displaying because it takes forever
  plot1 <- plot_feature_dependence(CONUS_Treeshap, variable = vars[1], title = "Oct-Apr Precipitation CDM Sum Feature Dependence") + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab(label = "Oct-Apr Total Precipitation CDM Sum (mm)") +
    ylab(label = "SHAP Value")
  plot2 <- plot_feature_dependence(CONUS_Treeshap, variable = vars[2], title = "Oct-Apr Mean Temperature CDM Sum Feature Dependence") + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab(label = expression("Oct-Apr Mean Temp CDM Sum ("*degree*C*")")) +
    ylab(label = "SHAP Value")
  plot3 <- plot_feature_dependence(CONUS_Treeshap, variable = vars[3], title = "Elevation Feature Dependence") + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab(label = "Elevation (m)") +
    ylab(label = "SHAP Value")
  plot4 <- plot_feature_dependence(CONUS_Treeshap, variable = vars[4], title = "Slope Feature Dependence") + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab(label = "Slope") +
    ylab(label = "SHAP Value")
  plot5 <- plot_feature_dependence(CONUS_Treeshap, variable = vars[5], title = "Aspect Feature Dependence") + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab(label = "Aspect") +
    ylab(label = "SHAP Value")
  plot6 <- plot_feature_dependence(CONUS_Treeshap, variable = vars[6], title = "Landcover Feature Dependence") + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab(label = "Landcover") +
    ylab(label = "SHAP Value")
  
  ### Saving Feature Dependence plots
  ggsave(here("Outputs", "Plots", "Treeshap", "FeatureDependence", "CONUS", vars[1], paste0("CONUS", year, "_", vars[1], "_FeatureDependence.png")), plot = plot1)
  ggsave(here("Outputs", "Plots", "Treeshap", "FeatureDependence", "CONUS", vars[2], paste0("CONUS", year, "_", vars[2], "_FeatureDependence.png")), plot = plot2)
  ggsave(here("Outputs", "Plots", "Treeshap", "FeatureDependence", "CONUS", vars[3], paste0("CONUS", year, "_", vars[3], "_FeatureDependence.png")), plot = plot3)
  ggsave(here("Outputs", "Plots", "Treeshap", "FeatureDependence", "CONUS", vars[4], paste0("CONUS", year, "_", vars[4], "_FeatureDependence.png")), plot = plot4)
  ggsave(here("Outputs", "Plots", "Treeshap", "FeatureDependence", "CONUS", vars[5], paste0("CONUS", year, "_", vars[5], "_FeatureDependence.png")), plot = plot5)
  ggsave(here("Outputs", "Plots", "Treeshap", "FeatureDependence", "CONUS", vars[6], paste0("CONUS", year, "_", vars[6], "_FeatureDependence.png")), plot = plot6)
  
  ### binding Classes to XY values from predictions, converting them to points, rasterizing them
  CONUSPredictions <- fread(file = here("Data", "SHAP", "PredictionCSVs", "CONUS", paste0("CONUS_SWEmaxPredictions", year, ".csv")), data.table = FALSE)
  
  # print("reading in predictor csv")
  CONUSPredictors <- fread(file = here("Data", "SHAP", "PredictorCSVs", "CONUS", paste0("CONUS_SWEmaxPredictors", year, ".csv")), data.table = FALSE) |>
    filter(paste0(x,y) %in% paste0(CONUSPredictions$x, CONUSPredictions$y))
  
  ### Extracting SHAP values, adding dummy rank columns which will be filled in later
  CONUS_SHAPs <- CONUS_Treeshap$shaps |>
    mutate(Rank_prcpSumCDM_varimp = -999,
           Rank_tmeanCDM_varimp = -999,
           Rank_elevation_varimp = -999,
           Rank_slope_varimp = -999,
           Rank_aspect_varimp = -999,
           Rank_landcover_varimp = -999,
           class = "dummy string")
  
  ### extracting only columns containing shap values
  shap_values <- CONUS_SHAPs[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "landcover_triclass")]
  # vars <- names(shap_values)
  
  ### ranking SHAP values by absolute value
  ranked_data <- t(apply(shap_values, 1, function(row) {
    abs_vals <- abs(row)
    ranks <- row_number(desc(abs_vals)) # Rank in descending order of absolute value
    ranked_vars <- vars[order(ranks)]
    return(c(ranks, paste(ranked_vars, collapse = ", ")))
  }))
  
  ### Assign the ranks back to the original data frame
  CONUS_SHAPs$Rank_prcpSumCDM_varimp <- as.numeric(ranked_data[, 1])
  CONUS_SHAPs$Rank_tmeanCDM_varimp <- as.numeric(ranked_data[, 2])
  CONUS_SHAPs$Rank_elevation_varimp <- as.numeric(ranked_data[, 3])
  CONUS_SHAPs$Rank_slope_varimp <- as.numeric(ranked_data[, 4])
  CONUS_SHAPs$Rank_aspect_varimp <- as.numeric(ranked_data[, 5])
  CONUS_SHAPs$Rank_landcover_varimp <- as.numeric(ranked_data[, 6])
  CONUS_SHAPs$class <- ranked_data[, 7]
  
  # for (i in 1:nrow(CONUS_SHAPs)){
  #   ### Ranking absolute values of SHAP values (distance from 0 == greater importance to prediction)
  #   ranks <- row_number(desc(c(abs(CONUS_SHAPs$OctApr_prcpSumCDMSum[i]), abs(CONUS_SHAPs$OctApr_tmeanCDMSum[i]), abs(CONUS_SHAPs$elevation[i]), abs(CONUS_SHAPs$slope[i]), abs(CONUS_SHAPs$aspect[i]), abs(CONUS_SHAPs$landcover_triclass[i]))))
  #   
  #   ### Assigning true variable importance ranks (according to SHAP)
  #   CONUS_SHAPs$Rank_prcpSumCDM_varimp[i] = ranks[1]
  #   CONUS_SHAPs$Rank_tmeanCDM_varimp[i] = ranks[2]
  #   CONUS_SHAPs$Rank_elevation_varimp[i] = ranks[3]
  #   CONUS_SHAPs$Rank_slope_varimp[i] = ranks[4]
  #   CONUS_SHAPs$Rank_aspect_varimp[i] = ranks[5]
  #   CONUS_SHAPs$Rank_landcover_varimp[i] = ranks[6]
  #   
  #   ### Assigning Class based on order of variables in ranking
  #   CONUS_SHAPs$class[i] <- paste(vars[which(ranks == 1)], vars[which(ranks == 2)], vars[which(ranks == 3)], vars[which(ranks == 4)], vars[which(ranks == 5)], vars[which(ranks == 6)], sep = ", ")
  # }
  
  ### saving SHAPs
  fwrite(CONUS_SHAPs, file = here("Data", "SHAP", "Treeshap", "SHAPs", "CONUS", paste0("CONUS", year, "SHAPs.csv")))
  
  ### reading in prediction raster to use as basis for rasterization
  PredictionRast <- rast(here("Outputs", "MapTIFFs", "CONUS", paste0("CONUS_SWEmax", year, ".tif")))
  CONUS_SHAPClasses <- terra::rasterize(vect(st_as_sf(cbind(CONUSPredictors, CONUS_SHAPs), coords = c("x", "y"), crs = crs(CONUS_AOI))), y = PredictionRast, field = "class")
  
  ### Saving Raster
  writeRaster(CONUS_SHAPClasses, file = here("Data", "SHAP", "Classes", "CONUS", paste0("CONUS", year, "SHAPClasses.tif")))
  
  
  
  ### Alaska
  ### reading in AK treeshap object
  AK_Treeshap <- read_rds(here("Data", "SHAP", "Treeshap", "Alaska", paste0("AKSWEMax_SHAPs", year, ".rds")))
  
  ### making plots of feature dependence
  ## only making them and saving them, not displaying because it takes forever
  plot1 <- plot_feature_dependence(AK_Treeshap, variable = vars[1]) + theme_bw()
  plot2 <- plot_feature_dependence(AK_Treeshap, variable = vars[2]) + theme_bw()
  plot3 <- plot_feature_dependence(AK_Treeshap, variable = vars[3]) + theme_bw()
  plot4 <- plot_feature_dependence(AK_Treeshap, variable = vars[4]) + theme_bw()
  plot5 <- plot_feature_dependence(AK_Treeshap, variable = vars[5]) + theme_bw()
  plot6 <- plot_feature_dependence(AK_Treeshap, variable = vars[6]) + theme_bw()
  
  ### Saving Feature Dependence plots
  ggsave(here("Outputs", "Plots", "Treeshap", "FeatureDependence", "Alaska", vars[1], paste0("Alaska", year, "_", vars[1], "_FeatureDependence.png")), plot = plot1)
  ggsave(here("Outputs", "Plots", "Treeshap", "FeatureDependence", "Alaska", vars[2], paste0("Alaska", year, "_", vars[2], "_FeatureDependence.png")), plot = plot2)
  ggsave(here("Outputs", "Plots", "Treeshap", "FeatureDependence", "Alaska", vars[3], paste0("Alaska", year, "_", vars[3], "_FeatureDependence.png")), plot = plot3)
  ggsave(here("Outputs", "Plots", "Treeshap", "FeatureDependence", "Alaska", vars[4], paste0("Alaska", year, "_", vars[4], "_FeatureDependence.png")), plot = plot4)
  ggsave(here("Outputs", "Plots", "Treeshap", "FeatureDependence", "Alaska", vars[5], paste0("Alaska", year, "_", vars[5], "_FeatureDependence.png")), plot = plot5)
  ggsave(here("Outputs", "Plots", "Treeshap", "FeatureDependence", "Alaska", vars[6], paste0("Alaska", year, "_", vars[6], "_FeatureDependence.png")), plot = plot6)
  
  ### Extracting SHAP values, adding dummy rank columns which will be filled in later
  AK_SHAPs <- AK_Treeshap$shaps |>
    mutate(Rank_prcpSumCDM_varimp = -999,
           Rank_tmeanCDM_varimp = -999,
           Rank_elevation_varimp = -999,
           Rank_slope_varimp = -999,
           Rank_aspect_varimp = -999,
           Rank_landcover_varimp = -999,
           class = "dummy string")
  
  
  ### extracting only columns containing shap values
  shap_values <- AK_SHAPs[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "landcover_triclass")]
  # vars <- names(shap_values)

  ### ranking SHAP values by absolute value
  ranked_data <- t(apply(shap_values, 1, function(row) {
    abs_vals <- abs(row)
    ranks <- row_number(desc(abs_vals)) # Rank in descending order of absolute value
    ranked_vars <- vars[order(ranks)]
    return(c(ranks, paste(ranked_vars, collapse = ", ")))
  }))

  ### Assign the ranks back to the original data frame
  AK_SHAPs$Rank_prcpSumCDM_varimp <- as.numeric(ranked_data[, 1])
  AK_SHAPs$Rank_tmeanCDM_varimp <- as.numeric(ranked_data[, 2])
  AK_SHAPs$Rank_elevation_varimp <- as.numeric(ranked_data[, 3])
  AK_SHAPs$Rank_slope_varimp <- as.numeric(ranked_data[, 4])
  AK_SHAPs$Rank_aspect_varimp <- as.numeric(ranked_data[, 5])
  AK_SHAPs$Rank_landcover_varimp <- as.numeric(ranked_data[, 6])
  AK_SHAPs$class <- ranked_data[, 7]
  
  # for (i in 1:nrow(AK_SHAPs)){
  #   ### Ranking absolute values of SHAP values (distance from 0 == greater importance to prediction)
  #   ranks <- row_number(desc(c(abs(AK_SHAPs$OctApr_prcpSumCDMSum[i]), abs(AK_SHAPs$OctApr_tmeanCDMSum[i]), abs(AK_SHAPs$elevation[i]), abs(AK_SHAPs$slope[i]), abs(AK_SHAPs$aspect[i]), abs(AK_SHAPs$landcover_triclass[i]))))
  #   
  #   ### Assigning true variable importance ranks (according to SHAP)
  #   AK_SHAPs$Rank_prcpSumCDM_varimp[i] = ranks[1]
  #   AK_SHAPs$Rank_tmeanCDM_varimp[i] = ranks[2]
  #   AK_SHAPs$Rank_elevation_varimp[i] = ranks[3]
  #   AK_SHAPs$Rank_slope_varimp[i] = ranks[4]
  #   AK_SHAPs$Rank_aspect_varimp[i] = ranks[5]
  #   AK_SHAPs$Rank_landcover_varimp[i] = ranks[6]
  #   
  #   ### Assigning Class based on order of variables in ranking
  #   AK_SHAPs$class[i] <- paste(vars[which(ranks == 1)], vars[which(ranks == 2)], vars[which(ranks == 3)], vars[which(ranks == 4)], vars[which(ranks == 5)], vars[which(ranks == 6)], sep = ", ")
  # }
  
  ### saving SHAPs
  write_rds(AK_SHAPs, file = here("Data", "SHAP", "Treeshap", "SHAPs", "Alaska", paste0("Alaska", year, "SHAPs.rds")), compress = "gz")
  
  ### binding Classes to XY values from predictions, converting them to points, rasterizing them
  AKPredictions <- fread(file = here("Data", "SHAP", "PredictionCSVs", "Alaska", paste0("AK_SWEmaxPredictions", year, ".csv")), data.table = FALSE)
  
  # print("reading in predictor csv")
  AKPredictors <- fread(file = here("Data", "SHAP", "PredictorCSVs", "Alaska", paste0("AK_SWEmaxPredictors", year, ".csv")), data.table = FALSE) |>
    filter(paste0(x,y) %in% paste0(AKPredictions$x, AKPredictions$y))
  
  ### reading in prediction raster to use as basis for rasterization
  PredictionRast <- rast(here("Outputs", "MapTIFFs", "Alaska", paste0("AK_SWEmax", year, ".tif")))
  AK_SHAPClasses <- terra::rasterize(vect(st_as_sf(cbind(AKPredictors, AK_SHAPs), coords = c("x", "y"), crs = crs(AK_AOI))), y = PredictionRast, field = "class")
  
  ### Saving Raster
  writeRaster(AK_SHAPClasses, file = here("Data", "SHAP", "Classes", "Alaska", paste0("Alaska", year, "SHAPClasses.tif")))
  
}


### Creating Classes
CreatedClasses <- pmclapply(years, SWEclass_creation, mc.cores = detectCores() - 2, mc.silent = FALSE)


