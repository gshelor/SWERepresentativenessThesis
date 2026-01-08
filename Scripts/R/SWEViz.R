##### Fixing plots from SWE_max maps so the scale is the same #####

library(pacman)
p_load(here, tidyverse, terra, sf, tidyterra, RColorBrewer, viridis, ggthemes, ggspatial, basemaps)
years <- 1993:2020

### CONUS AOI
CONUS_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "CONUS", "CONUSEcoregionsNoStates.gpkg"))
### Alaska AOI
AK_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "Alaska", "AK_AOI.gpkg"))
### basemaps
CONUS_bmap <- basemap_terra(ext = st_bbox(CONUS_AOI), map_service = "carto", map_type = "light")
AK_bmap <- basemap_terra(ext = st_bbox(AK_AOI), map_service = "carto", map_type = "light")

##### calculating max SWE value to use as the max value in plot's color scales #####
### reading in rasters, calculating max value, checking other rasters to see if their max is bigger
## setting max to 0 as the first value to beat
AKswemax_vals <- c()
AKswemean_vals <- c()
for (i in 1:length(years)){
  temp_rast <- rast(here("Outputs", "MapTIFFs", "Alaska", paste0("AK_SWEmax", years[i], ".tif")))
  temp_max <- max(values(temp_rast), na.rm = TRUE)
  temp_mean <- mean(values(temp_rast), na.rm = TRUE)
  AKswemax_vals <- c(AKswemax_vals, temp_max)
  AKswemean_vals <- c(AKswemean_vals, temp_mean)
}
### rounding highest SWE value to nearest number divisible by 100
AKswemax <- ceiling(mean(AKswemax_vals, na.rm = TRUE) / 100) * 100

## setting max to 0 as the first value to beat
CONUSswemax_vals <- c()
CONUSswemean_vals <- c()
for (i in 1:length(years)){
  temp_rast <- rast(here("Outputs", "MapTIFFs", "CONUS", paste0("CONUS_SWEmax", years[i], ".tif")))
  temp_max <- max(values(temp_rast), na.rm = TRUE)
  temp_mean <- mean(values(temp_rast), na.rm = TRUE)
  CONUSswemax_vals <- c(CONUSswemax_vals, temp_max)
  CONUSswemean_vals <- c(CONUSswemean_vals, temp_mean)
}
### rounding highest SWE value to nearest number divisible by 100
CONUSswemax <- ceiling(mean(CONUSswemax_vals, na.rm = TRUE) / 100) * 100

### making plots
for (i in 1:length(years)){
  ##### Alaska plots #####
  temp_rast <- rast(here("Outputs", "MapTIFFs", "Alaska", paste0("AK_SWEmax", years[i], ".tif")))
  
  AKSWE_plot <- ggplot() +
    theme_bw() +
    # geom_spatraster_rgb(data = AK_bmap) +
    geom_spatraster(data = temp_rast) +
    scale_fill_whitebox_c(palette = "deep", direction = 1, name = "Peak SWE (mm)", limits = c(0, AKswemax)) +
    geom_sf(data = AK_AOI, fill = NA, color = "black") +
    ggtitle(label = paste0("Predicted Peak SWE for Water Year ", years[i])) +
    ggspatial::annotation_scale(location = "br") +
    ggspatial::annotation_north_arrow(location = "tr", height = unit(1.25, "cm"), width = unit(1.25, "cm")) +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
  AKSWE_plot
  ### saving ggplot
  ggsave(here("Outputs", "Plots", "SWEMaxMaps", "Alaska", "1Km", paste0("AKmSWE", years[i], "_1km.png")))

  ### making plot without title so it will be better suited for publication (captions used in publication to describe/title plot)
  AKSWE_plot2 <- ggplot() +
    theme_bw() +
    # geom_spatraster_rgb(data = AK_bmap) +
    geom_spatraster(data = temp_rast) +
    scale_fill_whitebox_c(palette = "deep", direction = 1, name = "Peak SWE (mm)", limits = c(0, AKswemax)) +
    geom_sf(data = AK_AOI, fill = NA, color = "black") +
    # ggtitle(label = paste0("Predicted Peak SWE for Water Year ", year)) +
    ggspatial::annotation_scale(location = "br") +
    ggspatial::annotation_north_arrow(location = "tr", height = unit(1.25, "cm"), width = unit(1.25, "cm")) +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
  AKSWE_plot2
  ### saving ggplot
  ggsave(here("Outputs", "Plots", "SWEMaxMaps", "PaperFigures", "Alaska", paste0("AKmSWE", years[i], "_1km.png")))
  
  
  ##### Making CONUS plots #####
  temp_rast <- rast(here("Outputs", "MapTIFFs", "CONUS", paste0("CONUS_SWEmax", years[i], ".tif")))
  
  CONUSSWE_plot <- ggplot() +
    theme_bw() +
    # geom_spatraster_rgb(data = CONUS_bmap) +
    geom_spatraster(data = temp_rast) +
    scale_fill_whitebox_c(palette = "deep", direction = 1, name = "Peak SWE (mm)", limits = c(0, CONUSswemax)) +
    geom_sf(data = CONUS_AOI, fill = NA, color = "black") +
    ggtitle(label = paste0("Predicted Peak SWE for Water Year ", years[i])) +
    ggspatial::annotation_scale(location = "br") +
    ggspatial::annotation_north_arrow(location = "bl", height = unit(1.25, "cm"), width = unit(1.25, "cm")) +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
  CONUSSWE_plot
  ### saving ggplot
  ggsave(here("Outputs", "Plots", "SWEMaxMaps", "CONUS", "1Km", paste0("CONUSmSWE", years[i], ".png")))
  
  ### removing title for plots going into figures
  CONUSSWE_plot2 <- ggplot() +
    theme_bw() +
    # geom_spatraster_rgb(data = CONUS_bmap) +
    geom_spatraster(data = temp_rast) +
    scale_fill_whitebox_c(palette = "deep", direction = 1, name = "Peak SWE (mm)", limits = c(0, CONUSswemax)) +
    geom_sf(data = CONUS_AOI, fill = NA, color = "black") +
    # ggtitle(label = paste0("Predicted Peak SWE for Water Year ", year)) +
    ggspatial::annotation_scale(location = "br") +
    ggspatial::annotation_north_arrow(location = "bl", height = unit(1.25, "cm"), width = unit(1.25, "cm")) +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
  CONUSSWE_plot2
  ### saving ggplot
  ggsave(here("Outputs", "Plots", "SWEMaxMaps", "PaperFigures", "CONUS", paste0("CONUSmSWE", years[i], ".png")))
}



##### Calculating Average max SWE prediction for each cell #####
### Alaska
for (i in 1:length(years)){
  temp_rast <- rast(here("Outputs", "MapTIFFs", "Alaska", paste0("AK_SWEmax", years[i], ".tif")))
  if (i == 1){
    AKSWE_rasts <- temp_rast
  } else{
    AKSWE_rasts <- c(AKSWE_rasts, temp_rast)
  }
}

### calculating mean
AKSWEMean_rast <- app(AKSWE_rasts, fun = "mean")
plot(AKSWEMean_rast)
writeRaster(AKSWEMean_rast, here("Outputs", "MapTIFFs", "Alaska", "Summary", "AKSWEmaxMean.tif"))

### CONUS
for (i in 1:length(years)){
  temp_rast <- rast(here("Outputs", "MapTIFFs", "CONUS", paste0("CONUS_SWEmax", years[i], ".tif")))
  if (i == 1){
    CONUSSWE_rasts <- temp_rast
  } else{
    CONUSSWE_rasts <- c(CONUSSWE_rasts, temp_rast)
  }
}

### calculating mean
CONUSSWEMean_rast <- app(CONUSSWE_rasts, fun = "mean")
plot(CONUSSWEMean_rast)
max(CONUSSWEMean_rast, na.rm = TRUE)
writeRaster(CONUSSWEMean_rast, here("Outputs", "MapTIFFs", "CONUS", "Summary", "CONUSSWEmaxMean.tif"))
