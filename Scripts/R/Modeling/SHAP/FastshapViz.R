##### redoing feature importance plots from fastshap #####

##### loading packages #####
library(pacman)
p_load(here, tidyverse, data.table, shapviz, fastshap, ggthemes)



##### Function to read in Fastshap Object, create shapviz object, make FD plots #####
years <- 1993:2020
### saving variable names to be listed in ranking order later
vars <- c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "landcover_triclass")
AKvars <- c("OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "DecFeb_prcpSumCDMSum", "DecFeb_tmeanCDMSum", "elevation", "slope", "aspect", "OctMay_tminMean", "OctMay_tmaxMean", "landcover_triclass")

### creating folders to save feature dependence plots if they don't exist
# for (x in vars){
#   if (dir.exists(here("Outputs", "Plots", "Fastshap", "FeatureDependence", "CONUS", paste0(x, "/"))) == FALSE){
#     dir.create(here("Outputs", "Plots", "Fastshap", "FeatureDependence", "CONUS",  paste0(x, "/")))
#   }
#   if (dir.exists(here("Outputs", "Plots", "Fastshap", "FeatureDependence", "PaperFigures", "CONUS",  paste0(x, "/"))) == FALSE){
#     dir.create(here("Outputs", "Plots", "Fastshap", "FeatureDependence", "PaperFigures", "CONUS",  paste0(x, "/")))
#   }
# }
# for (x in AKvars){
#   if (dir.exists(here("Outputs", "Plots", "Fastshap", "FeatureDependence", "Alaska",  paste0(x, "/"))) == FALSE){
#     dir.create(here("Outputs", "Plots", "Fastshap", "FeatureDependence", "Alaska",  paste0(x, "/")))
#   }
#   if (dir.exists(here("Outputs", "Plots", "Fastshap", "FeatureDependence", "PaperFigures", "Alaska",  paste0(x, "/"))) == FALSE){
#     dir.create(here("Outputs", "Plots", "Fastshap", "FeatureDependence", "PaperFigures", "Alaska",  paste0(x, "/")))
#   }
# }

ShapViz_func <- function(year){
  ### reading in CONUS predictors because ShapViz needs that for some fucking reason
  # CONUSPredictors <- fread(file = here("Data", "SHAP", "PredictorCSVs", "CONUS", paste0("CONUS_SWEmaxPredictors", year, ".csv")), data.table = FALSE) |>
  #   select(OctApr_prcpSumCDMSum, OctApr_tmeanCDMSum, elevation, slope, aspect, landcover_triclass)
  # ### reading in fastshap object
  # CONUSFastshap <- read_rds(here("Data", "SHAP", "Fastshap", "CONUS", paste0("CONUSSWEmax_SHAPs", year, ".rds")))
  # 
  # ### making plots of feature dependence
  # ## only making them and saving them, not displaying because it takes forever
  # ### making shapviz object to create feature dependence plot
  # sv <- shapviz(object = CONUSFastshap, X = CONUSPredictors)
  # plot1 <- sv_dependence(sv, v = vars[1], color_var = vars[1], jitter_width = 0) +
  #   ggtitle("Oct-Apr Precipitation CDM Sum Feature Dependence", subtitle = paste("Western US", as.character(year))) + 
  #   theme_bw() +
  #   theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  #   xlab(label = "Oct-Apr Total Precipitation CDM Sum (mm)") +
  #   ylab(label = "SHAP Value")
  # plot2 <- sv_dependence(sv, v = vars[2], color_var = vars[2], jitter_width = 0) +
  #   ggtitle("Oct-Apr Mean Temperature CDM Sum Feature Dependence", subtitle = paste("Western US", as.character(year))) + 
  #   theme_bw() +
  #   theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  #   xlab(label = expression("Oct-Apr Mean Temp CDM Sum ("*degree*C*")")) +
  #   ylab(label = "SHAP Value")
  # plot3 <- sv_dependence(sv, v = vars[3], color_var = vars[3], jitter_width = 0) +
  #   ggtitle("Elevation Feature Dependence", subtitle = paste("Western US", as.character(year))) + 
  #   theme_bw() +
  #   theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  #   xlab(label = "Elevation (m)") +
  #   ylab(label = "SHAP Value")
  # plot4 <- sv_dependence(sv, v = vars[4], color_var = vars[4], jitter_width = 0) +
  #   ggtitle("Slope Feature Dependence", subtitle = paste("Western US", as.character(year))) + 
  #   theme_bw() +
  #   theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  #   xlab(label = "Slope") +
  #   ylab(label = "SHAP Value")
  # plot5 <- sv_dependence(sv, v = vars[5], color_var = vars[5], jitter_width = 0) +
  #   ggtitle("Aspect Feature Dependence", subtitle = paste("Western US", as.character(year))) + 
  #   theme_bw() +
  #   theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  #   xlab(label = "Aspect") +
  #   ylab(label = "SHAP Value")
  # plot6 <- sv_dependence(sv, v = vars[6], color_var = vars[6], jitter_width = 0) +
  #   ggtitle("Landcover Feature Dependence", subtitle = paste("Western US", as.character(year))) + 
  #   theme_bw() +
  #   theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  #   xlab(label = "Landcover") +
  #   ylab(label = "SHAP Value")
  # 
  # ### Saving Feature Dependence plots
  # ggsave(here("Outputs", "Plots", "Fastshap", "FeatureDependence", "CONUS", vars[1], paste0("CONUS", year, "_", vars[1],  "_FeatureDependence.png")), plot = plot1)
  # ggsave(here("Outputs", "Plots", "Fastshap", "FeatureDependence", "CONUS", vars[2], paste0("CONUS", year, "_", vars[2], "_FeatureDependence.png")), plot = plot2)
  # ggsave(here("Outputs", "Plots", "Fastshap", "FeatureDependence", "CONUS", vars[3], paste0("CONUS", year, "_", vars[3], "_FeatureDependence.png")), plot = plot3)
  # ggsave(here("Outputs", "Plots", "Fastshap", "FeatureDependence", "CONUS", vars[4], paste0("CONUS", year, "_", vars[4], "_FeatureDependence.png")), plot = plot4)
  # ggsave(here("Outputs", "Plots", "Fastshap", "FeatureDependence", "CONUS", vars[5], paste0("CONUS", year, "_", vars[5], "_FeatureDependence.png")), plot = plot5)
  # ggsave(here("Outputs", "Plots", "Fastshap", "FeatureDependence", "CONUS", vars[6], paste0("CONUS", year, "_", vars[6], "_FeatureDependence.png")), plot = plot6)
  # 
  # 
  # 
  # ##### now FD plots but without the title (for the paper) #####
  # plot1 <- sv_dependence(sv, v = vars[1], color_var = vars[1], jitter_width = 0) +
  #   ggtitle("", subtitle = paste("Western US", as.character(year))) + 
  #   theme_bw() +
  #   theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  #   xlab(label = "Oct-Apr Total Precipitation CDM Sum (mm)") +
  #   ylab(label = "SHAP Value")
  # plot2 <- sv_dependence(sv, v = vars[2], color_var = vars[2], jitter_width = 0) +
  #   ggtitle("", subtitle = paste("Western US", as.character(year))) + 
  #   theme_bw() +
  #   theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  #   xlab(label = expression("Oct-Apr Mean Temp CDM Sum ("*degree*C*")")) +
  #   ylab(label = "SHAP Value")
  # plot3 <- sv_dependence(sv, v = vars[3], color_var = vars[3], jitter_width = 0) +
  #   ggtitle("", subtitle = paste("Western US", as.character(year))) + 
  #   theme_bw() +
  #   theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  #   xlab(label = "Elevation (m)") +
  #   ylab(label = "SHAP Value")
  # plot4 <- sv_dependence(sv, v = vars[4], color_var = vars[4], jitter_width = 0) +
  #   ggtitle("", subtitle = paste("Western US", as.character(year))) + 
  #   theme_bw() +
  #   theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  #   xlab(label = "Slope") +
  #   ylab(label = "SHAP Value")
  # plot5 <- sv_dependence(sv, v = vars[5], color_var = vars[5], jitter_width = 0) +
  #   ggtitle("", subtitle = paste("Western US", as.character(year))) + 
  #   theme_bw() +
  #   theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  #   xlab(label = "Aspect") +
  #   ylab(label = "SHAP Value")
  # plot6 <- sv_dependence(sv, v = vars[6], color_var = vars[6], jitter_width = 0) +
  #   ggtitle("", subtitle = paste("Western US", as.character(year))) + 
  #   theme_bw() +
  #   theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  #   xlab(label = "Landcover") +
  #   ylab(label = "SHAP Value")
  # 
  # ### Saving Feature Dependence plots
  # ggsave(here("Outputs", "Plots", "Fastshap", "FeatureDependence", "PaperFigures", "CONUS", vars[1], paste0("CONUS", year, "_", vars[1],  "_FeatureDependence.png")), plot = plot1)
  # ggsave(here("Outputs", "Plots", "Fastshap", "FeatureDependence", "PaperFigures", "CONUS", vars[2], paste0("CONUS", year, "_", vars[2], "_FeatureDependence.png")), plot = plot2)
  # ggsave(here("Outputs", "Plots", "Fastshap", "FeatureDependence", "PaperFigures", "CONUS", vars[3], paste0("CONUS", year, "_", vars[3], "_FeatureDependence.png")), plot = plot3)
  # ggsave(here("Outputs", "Plots", "Fastshap", "FeatureDependence", "PaperFigures", "CONUS", vars[4], paste0("CONUS", year, "_", vars[4], "_FeatureDependence.png")), plot = plot4)
  # ggsave(here("Outputs", "Plots", "Fastshap", "FeatureDependence", "PaperFigures", "CONUS", vars[5], paste0("CONUS", year, "_", vars[5], "_FeatureDependence.png")), plot = plot5)
  # ggsave(here("Outputs", "Plots", "Fastshap", "FeatureDependence", "PaperFigures", "CONUS", vars[6], paste0("CONUS", year, "_", vars[6], "_FeatureDependence.png")), plot = plot6)
  
  
  ##### Alaska shapviz #####
  ### reading in predictors
  AKPredictors <- fread(file = here("Data", "SHAP", "PredictorCSVs", "Alaska", paste0("AK_SWEmaxPredictors", year, ".csv")), data.table = FALSE) |>
    select(OctMay_prcpSumCDMSum, OctMay_tmeanCDMSum, DecFeb_prcpSumCDMSum, DecFeb_tmeanCDMSum, elevation, slope, aspect, OctMay_tminMean, OctMay_tmaxMean, landcover_triclass)
  
  ### reading in AK fastshap object
  AK_Fastshap <- read_rds(here("Data", "SHAP", "Fastshap", "Alaska", paste0("AKSWEmax_SHAPs", year, ".rds")))
  
  ### making plots of feature dependence
  ## only making them and saving them, not displaying because it takes forever
  sv <- shapviz(object = AK_Fastshap, X = AKPredictors)
  
  ### now without the title
  plot1 <- sv_dependence(sv, v = AKvars[1], color_var = AKvars[1], jitter_width = 0) +
    ggtitle("", subtitle = paste("Alaska", as.character(year))) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
    xlab(label = "Oct-May Total Precipitation CDM Sum (mm)") +
    ylab(label = "SHAP Value")
  plot2 <- sv_dependence(sv, v = AKvars[2], color_var = AKvars[2], jitter_width = 0) +
    ggtitle("", subtitle = paste("Alaska", as.character(year))) + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
    xlab(label = expression("Oct-May Mean Temp CDM Sum ("*degree*C*")")) +
    ylab(label = "SHAP Value")
  plot3 <- sv_dependence(sv, v = AKvars[3], color_var = AKvars[3], jitter_width = 0) +
    ggtitle("", subtitle = paste("Alaska", as.character(year))) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
    xlab(label = "Dec-Feb Total Precipitation CDM Sum (mm)") +
    ylab(label = "SHAP Value")
  plot4 <- sv_dependence(sv, v = AKvars[4], color_var = AKvars[4], jitter_width = 0) +
    ggtitle("", subtitle = paste("Alaska", as.character(year))) + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
    xlab(label = expression("Dec-Feb Mean Temp CDM Sum ("*degree*C*")")) +
    ylab(label = "SHAP Value")
  plot5 <- sv_dependence(sv, v = AKvars[5], color_var = AKvars[5], jitter_width = 0) +
    ggtitle("", subtitle = paste("Alaska", as.character(year))) + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
    xlab(label = "Elevation (m)") +
    ylab(label = "SHAP Value")
  plot6 <- sv_dependence(sv, v = AKvars[6], color_var = AKvars[6], jitter_width = 0) +
    ggtitle("", subtitle = paste("Alaska", as.character(year))) + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
    xlab(label = "Slope") +
    ylab(label = "SHAP Value")
  plot7 <- sv_dependence(sv, v = AKvars[7], color_var = AKvars[7], jitter_width = 0) +
    ggtitle("", subtitle = paste("Alaska", as.character(year))) + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
    xlab(label = "Aspect") +
    ylab(label = "SHAP Value")
  plot8 <- sv_dependence(sv, v = AKvars[8], color_var = AKvars[8], jitter_width = 0) +
    ggtitle("", subtitle = paste("Alaska", as.character(year))) + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
    xlab(label = expression("Oct-May Mean Min Temp ("*degree*C*")")) +
    ylab(label = "SHAP Value")
  plot9 <- sv_dependence(sv, v = AKvars[9], color_var = AKvars[9], jitter_width = 0) +
    ggtitle("", subtitle = paste("Alaska", as.character(year))) + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
    xlab(label = expression("Oct-May Mean Max Temp ("*degree*C*")")) +
    ylab(label = "SHAP Value")
  plot10 <- sv_dependence(sv, v = AKvars[10], color_var = AKvars[10], jitter_width = 0) +
    ggtitle("", subtitle = paste("Alaska", as.character(year))) + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
    xlab(label = "Landcover") +
    ylab(label = "SHAP Value")
  
  ### Saving Feature Dependence plots
  ggsave(here("Outputs", "Plots", "Fastshap", "FeatureDependence", "PaperFigures", "Alaska", AKvars[1], paste0("Alaska", year, "_", AKvars[1], "_FeatureDependence.png")), plot = plot1)
  ggsave(here("Outputs", "Plots", "Fastshap", "FeatureDependence", "PaperFigures", "Alaska", AKvars[2], paste0("Alaska", year, "_", AKvars[2], "_FeatureDependence.png")), plot = plot2)
  ggsave(here("Outputs", "Plots", "Fastshap", "FeatureDependence", "PaperFigures", "Alaska", AKvars[3], paste0("Alaska", year, "_", AKvars[3], "_FeatureDependence.png")), plot = plot3)
  ggsave(here("Outputs", "Plots", "Fastshap", "FeatureDependence", "PaperFigures", "Alaska", AKvars[4], paste0("Alaska", year, "_", AKvars[4], "_FeatureDependence.png")), plot = plot4)
  ggsave(here("Outputs", "Plots", "Fastshap", "FeatureDependence", "PaperFigures", "Alaska", AKvars[5], paste0("Alaska", year, "_", AKvars[5], "_FeatureDependence.png")), plot = plot5)
  ggsave(here("Outputs", "Plots", "Fastshap", "FeatureDependence", "PaperFigures", "Alaska", AKvars[6], paste0("Alaska", year, "_", AKvars[6], "_FeatureDependence.png")), plot = plot6)
  ggsave(here("Outputs", "Plots", "Fastshap", "FeatureDependence", "PaperFigures", "Alaska", AKvars[7], paste0("Alaska", year, "_", AKvars[7], "_FeatureDependence.png")), plot = plot7)
  ggsave(here("Outputs", "Plots", "Fastshap", "FeatureDependence", "PaperFigures", "Alaska", AKvars[8], paste0("Alaska", year, "_", AKvars[8], "_FeatureDependence.png")), plot = plot8)
  ggsave(here("Outputs", "Plots", "Fastshap", "FeatureDependence", "PaperFigures", "Alaska", AKvars[9], paste0("Alaska", year, "_", AKvars[9], "_FeatureDependence.png")), plot = plot9)
  ggsave(here("Outputs", "Plots", "Fastshap", "FeatureDependence", "PaperFigures", "Alaska", AKvars[10], paste0("Alaska", year, "_", AKvars[10], "_FeatureDependence.png")), plot = plot10)
}


### ok now actually creating and saving the figures
pmclapply(X = years, FUN = ShapViz_func, mc.cores = length(years), mc.silent = FALSE, mc.set.seed = FALSE)



