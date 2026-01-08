##### Daymet Climatologies #####
### only using 1992-2020 since to calculate these since those are the years of my study period
### Total start time
StartTime <- Sys.time()
### loading packages
library(pacman)
p_load(here, tidyverse, sf, terra, ncdf4, parallel, tidyterra)


AlaskaStart <- Sys.time()
##### creating Alaska climatologies using mosaics of monthly variable aggregates #####
print("starting climatology calculation for Alaska")
for (variable in c("prcp", "tmin", "tmean", "tmax", "srad")){
  print(variable)
  ### creating vector to store file paths
  meanfilepaths <- c()
  medianfilepaths <- c()
  for (daymetyear in 1992:2023){
    currentyear <- as.character(daymetyear)
    variablemeanfiles_temp <- list.files(path = here("Data", "Daymet", "Alaska", "Mosaics"), pattern = paste0("^", variable, "_", currentyear, "MonthlyMean", "*"), all.files= T, full.names= T)
    # variablemedianfiles_temp <- list.files(path = here("Data", "Daymet", "Alaska", "Mosaics"), pattern = paste0("^", variable, "_", currentyear, "MonthlyMedian", "*"), all.files= T, full.names= T)
    meanfilepaths <- c(meanfilepaths, variablemeanfiles_temp)
    # medianfilepaths <- c(medianfilepaths, variablemedianfiles_temp)
  }
  rast1 <- rast(meanfilepaths[1])
  rast2 <- rast(meanfilepaths[2])
  rast3 <- rast(meanfilepaths[3])
  rast4 <- rast(meanfilepaths[4])
  rast5 <- rast(meanfilepaths[5])
  rast6 <- rast(meanfilepaths[6])
  rast7 <- rast(meanfilepaths[7])
  rast8 <- rast(meanfilepaths[8])
  rast9 <- rast(meanfilepaths[9])
  rast10 <- rast(meanfilepaths[10])
  rast11 <- rast(meanfilepaths[11])
  rast12 <- rast(meanfilepaths[12])
  rast13 <- rast(meanfilepaths[13])
  rast14 <- rast(meanfilepaths[14])
  rast15 <- rast(meanfilepaths[15])
  rast16 <- rast(meanfilepaths[16])
  rast17 <- rast(meanfilepaths[17])
  rast18 <- rast(meanfilepaths[18])
  rast19 <- rast(meanfilepaths[19])
  rast20 <- rast(meanfilepaths[20])
  rast21 <- rast(meanfilepaths[21])
  rast22 <- rast(meanfilepaths[22])
  rast23 <- rast(meanfilepaths[23])
  rast24 <- rast(meanfilepaths[24])
  rast25 <- rast(meanfilepaths[25])
  rast26 <- rast(meanfilepaths[26])
  rast27 <- rast(meanfilepaths[27])
  rast28 <- rast(meanfilepaths[28])
  rast29 <- rast(meanfilepaths[29])
  rast30 <- rast(meanfilepaths[30])
  rast31 <- rast(meanfilepaths[31])
  rast32 <- rast(meanfilepaths[32])
  
  ### concatenating specific months of each raster to create climatological average over study period
  for (lyr in 1:nlyr(rast1)){
    if (lyr < 10){
      climaterast <- c(rast2[[lyr]], rast3[[lyr]], rast4[[lyr]], rast5[[lyr]], rast6[[lyr]], rast7[[lyr]], rast8[[lyr]], rast9[[lyr]], rast10[[lyr]], rast11[[lyr]], rast12[[lyr]], rast13[[lyr]], rast14[[lyr]], rast15[[lyr]], rast16[[lyr]], rast17[[lyr]], rast18[[lyr]], rast19[[lyr]], rast20[[lyr]], rast21[[lyr]], rast22[[lyr]], rast23[[lyr]], rast24[[lyr]], rast25[[lyr]], rast26[[lyr]], rast27[[lyr]], rast28[[lyr]], rast29[[lyr]], rast30[[lyr]])
      
      ### taking climatological average
      climatologyrast <- app(climaterast, fun = mean)
      # climatologyrast <- app(climaterast, fun = median)
      
      ### saving climatology raster as tif
      writeRaster(climatologyrast, here("Data", "Daymet", "Alaska", "Climatology", "FinalTIFFs", paste0(variable, "_Month", as.character(lyr), "MeanClimatology.tif")), overwrite = TRUE)
      # writeRaster(climatologyrast, here("Data", "Daymet", "Alaska", "Climatology", "FinalTIFFs", paste0(variable, "_Month", as.character(lyr), "MedianClimatology.tif")))
    } else {
      climaterast <- c(rast1[[lyr]], rast2[[lyr]], rast3[[lyr]], rast4[[lyr]], rast5[[lyr]], rast6[[lyr]], rast7[[lyr]], rast8[[lyr]], rast9[[lyr]], rast10[[lyr]], rast11[[lyr]], rast12[[lyr]], rast13[[lyr]], rast14[[lyr]], rast15[[lyr]], rast16[[lyr]], rast17[[lyr]], rast18[[lyr]], rast19[[lyr]], rast20[[lyr]], rast21[[lyr]], rast22[[lyr]], rast23[[lyr]], rast24[[lyr]], rast25[[lyr]], rast26[[lyr]], rast27[[lyr]], rast28[[lyr]], rast29[[lyr]])
      
      ### taking climatological average
      climatologyrast <- app(climaterast, fun = mean)
      # climatologyrast <- app(climaterast, fun = median)
      
      ### saving climatology raster as tif
      writeRaster(climatologyrast, here("Data", "Daymet", "Alaska", "Climatology", "FinalTIFFs", paste0(variable, "_Month", as.character(lyr), "MeanClimatology.tif")), overwrite = TRUE)
      # writeRaster(climatologyrast, here("Data", "Daymet", "Alaska", "Climatology", "FinalTIFFs", paste0(variable, "_Month", as.character(lyr), "MedianClimatology.tif")))
    }
  }
}


print("starting PrcpSum climatology calculation for Alaska")
for (variable in c("prcp")){
  print(variable)
  ### creating vector to store file paths
  meanfilepaths <- c()
  medianfilepaths <- c()
  for (daymetyear in 1992:2023){
    currentyear <- as.character(daymetyear)
    variablemeanfiles_temp <- list.files(path = here("Data", "Daymet", "Alaska", "Mosaics"), pattern = paste0("^", variable, "_", currentyear, "MonthlySum", "*"), all.files= T, full.names= T)
    # variablemedianfiles_temp <- list.files(path = here("Data", "Daymet", "Alaska", "Mosaics"), pattern = paste0("^", variable, "_", currentyear, "MonthlyMedian", "*"), all.files= T, full.names= T)
    meanfilepaths <- c(meanfilepaths, variablemeanfiles_temp)
    # medianfilepaths <- c(medianfilepaths, variablemedianfiles_temp)
  }
  rast1 <- rast(meanfilepaths[1])
  rast2 <- rast(meanfilepaths[2])
  rast3 <- rast(meanfilepaths[3])
  rast4 <- rast(meanfilepaths[4])
  rast5 <- rast(meanfilepaths[5])
  rast6 <- rast(meanfilepaths[6])
  rast7 <- rast(meanfilepaths[7])
  rast8 <- rast(meanfilepaths[8])
  rast9 <- rast(meanfilepaths[9])
  rast10 <- rast(meanfilepaths[10])
  rast11 <- rast(meanfilepaths[11])
  rast12 <- rast(meanfilepaths[12])
  rast13 <- rast(meanfilepaths[13])
  rast14 <- rast(meanfilepaths[14])
  rast15 <- rast(meanfilepaths[15])
  rast16 <- rast(meanfilepaths[16])
  rast17 <- rast(meanfilepaths[17])
  rast18 <- rast(meanfilepaths[18])
  rast19 <- rast(meanfilepaths[19])
  rast20 <- rast(meanfilepaths[20])
  rast21 <- rast(meanfilepaths[21])
  rast22 <- rast(meanfilepaths[22])
  rast23 <- rast(meanfilepaths[23])
  rast24 <- rast(meanfilepaths[24])
  rast25 <- rast(meanfilepaths[25])
  rast26 <- rast(meanfilepaths[26])
  rast27 <- rast(meanfilepaths[27])
  rast28 <- rast(meanfilepaths[28])
  rast29 <- rast(meanfilepaths[29])
  rast30 <- rast(meanfilepaths[30])
  rast31 <- rast(meanfilepaths[31])
  rast32 <- rast(meanfilepaths[32])
  # rast33 <- rast(meanfilepaths[33])
  
  ### concatenating specific months of each raster to create climatological average over study period
  for (lyr in 1:nlyr(rast1)){
    if (lyr < 10){
      climaterast <- c(rast2[[lyr]], rast3[[lyr]], rast4[[lyr]], rast5[[lyr]], rast6[[lyr]], rast7[[lyr]], rast8[[lyr]], rast9[[lyr]], rast10[[lyr]], rast11[[lyr]], rast12[[lyr]], rast13[[lyr]], rast14[[lyr]], rast15[[lyr]], rast16[[lyr]], rast17[[lyr]], rast18[[lyr]], rast19[[lyr]], rast20[[lyr]], rast21[[lyr]], rast22[[lyr]], rast23[[lyr]], rast24[[lyr]], rast25[[lyr]], rast26[[lyr]], rast27[[lyr]], rast28[[lyr]], rast29[[lyr]], rast30[[lyr]])
      
      ### taking climatological average
      climatologyrast <- app(climaterast, fun = mean)
      # climatologyrast <- app(climaterast, fun = median)
      
      ### saving climatology raster as tif
      writeRaster(climatologyrast, here("Data", "Daymet", "Alaska", "Climatology", "FinalTIFFs", paste0(variable, "Sum_Month", as.character(lyr), "MeanClimatology.tif")), overwrite = TRUE)
      # writeRaster(climatologyrast, here("Data", "Daymet", "Alaska", "Climatology", "FinalTIFFs", paste0(variable, "Sum_Month", as.character(lyr), "MedianClimatology.tif")))
    } else {
      climaterast <- c(rast1[[lyr]], rast2[[lyr]], rast3[[lyr]], rast4[[lyr]], rast5[[lyr]], rast6[[lyr]], rast7[[lyr]], rast8[[lyr]], rast9[[lyr]], rast10[[lyr]], rast11[[lyr]], rast12[[lyr]], rast13[[lyr]], rast14[[lyr]], rast15[[lyr]], rast16[[lyr]], rast17[[lyr]], rast18[[lyr]], rast19[[lyr]], rast20[[lyr]], rast21[[lyr]], rast22[[lyr]], rast23[[lyr]], rast24[[lyr]], rast25[[lyr]], rast26[[lyr]], rast27[[lyr]], rast28[[lyr]], rast29[[lyr]])
      
      ### taking climatological average
      climatologyrast <- app(climaterast, fun = mean)
      # climatologyrast <- app(climaterast, fun = median)
      
      ### saving climatology raster as tif
      writeRaster(climatologyrast, here("Data", "Daymet", "Alaska", "Climatology", "FinalTIFFs", paste0(variable, "Sum_Month", as.character(lyr), "MeanClimatology.tif")), overwrite = TRUE)
      # writeRaster(climatologyrast, here("Data", "Daymet", "Alaska", "Climatology", "FinalTIFFs", paste0(variable, "_Month", as.character(lyr), "MedianClimatology.tif")))
    }
  }
}
AlaskaEnd <- Sys.time()







ConusStart <- Sys.time()
##### creating CONUS climatologies using mosaics of variables #####
print("starting climatology calculation for CONUS")
for (variable in c("prcp", "tmin", "tmean", "tmax", "srad")){
  print(variable)
  ### creating vector to store file paths
  meanfilepaths <- c()
  medianfilepaths <- c()
  for (daymetyear in 1992:2023){
    currentyear <- as.character(daymetyear)
    variablemeanfiles_temp <- list.files(path = here("Data", "Daymet", "CONUS", "Mosaics", "NetCDFs"), pattern = paste0("^", variable, "_", currentyear, "MonthlyMean", "*"), all.files= T, full.names= T)
    # variablemedianfiles_temp <- list.files(path = here("Data", "Daymet", "CONUS", "Mosaics"), pattern = paste0("^", variable, "_", currentyear, "MonthlyMedian", "*"), all.files= T, full.names= T)
    meanfilepaths <- c(meanfilepaths, variablemeanfiles_temp)
    # medianfilepaths <- c(medianfilepaths, variablemedianfiles_temp)
  }
  rast1 <- rast(meanfilepaths[1])
  rast2 <- rast(meanfilepaths[2])
  rast3 <- rast(meanfilepaths[3])
  rast4 <- rast(meanfilepaths[4])
  rast5 <- rast(meanfilepaths[5])
  rast6 <- rast(meanfilepaths[6])
  rast7 <- rast(meanfilepaths[7])
  rast8 <- rast(meanfilepaths[8])
  rast9 <- rast(meanfilepaths[9])
  rast10 <- rast(meanfilepaths[10])
  rast11 <- rast(meanfilepaths[11])
  rast12 <- rast(meanfilepaths[12])
  rast13 <- rast(meanfilepaths[13])
  rast14 <- rast(meanfilepaths[14])
  rast15 <- rast(meanfilepaths[15])
  rast16 <- rast(meanfilepaths[16])
  rast17 <- rast(meanfilepaths[17])
  rast18 <- rast(meanfilepaths[18])
  rast19 <- rast(meanfilepaths[19])
  rast20 <- rast(meanfilepaths[20])
  rast21 <- rast(meanfilepaths[21])
  rast22 <- rast(meanfilepaths[22])
  rast23 <- rast(meanfilepaths[23])
  rast24 <- rast(meanfilepaths[24])
  rast25 <- rast(meanfilepaths[25])
  rast26 <- rast(meanfilepaths[26])
  rast27 <- rast(meanfilepaths[27])
  rast28 <- rast(meanfilepaths[28])
  rast29 <- rast(meanfilepaths[29])
  rast30 <- rast(meanfilepaths[30])
  rast31 <- rast(meanfilepaths[31])
  rast32 <- rast(meanfilepaths[32])
  # rast33 <- rast(meanfilepaths[33])
  
  ### concatenating specific months of each raster to create climatological average over study period
  for (lyr in 1:nlyr(rast1)){
    if (lyr < 10){
      climaterast <- c(rast2[[lyr]], rast3[[lyr]], rast4[[lyr]], rast5[[lyr]], rast6[[lyr]], rast7[[lyr]], rast8[[lyr]], rast9[[lyr]], rast10[[lyr]], rast11[[lyr]], rast12[[lyr]], rast13[[lyr]], rast14[[lyr]], rast15[[lyr]], rast16[[lyr]], rast17[[lyr]], rast18[[lyr]], rast19[[lyr]], rast20[[lyr]], rast21[[lyr]], rast22[[lyr]], rast23[[lyr]], rast24[[lyr]], rast25[[lyr]], rast26[[lyr]], rast27[[lyr]], rast28[[lyr]], rast29[[lyr]], rast30[[lyr]])
      
      ### taking climatological average
      climatologyrast <- app(climaterast, fun = mean)
      # climatologyrast <- app(climaterast, fun = median)
      
      ### saving climatology raster as tif
      writeRaster(climatologyrast, here("Data", "Daymet", "CONUS", "Climatology", "FinalTIFFs", paste0(variable, "_Month", as.character(lyr), "MeanClimatology.tif")))
      # writeRaster(climatologyrast, here("Data", "Daymet", "CONUS", "Climatology", "FinalTIFFs", paste0(variable, "_Month", as.character(lyr), "MedianClimatology.tif")))
    } else {
      climaterast <- c(rast1[[lyr]], rast2[[lyr]], rast3[[lyr]], rast4[[lyr]], rast5[[lyr]], rast6[[lyr]], rast7[[lyr]], rast8[[lyr]], rast9[[lyr]], rast10[[lyr]], rast11[[lyr]], rast12[[lyr]], rast13[[lyr]], rast14[[lyr]], rast15[[lyr]], rast16[[lyr]], rast17[[lyr]], rast18[[lyr]], rast19[[lyr]], rast20[[lyr]], rast21[[lyr]], rast22[[lyr]], rast23[[lyr]], rast24[[lyr]], rast25[[lyr]], rast26[[lyr]], rast27[[lyr]], rast28[[lyr]], rast29[[lyr]])
      
      ### taking climatological average
      climatologyrast <- app(climaterast, fun = mean)
      # climatologyrast <- app(climaterast, fun = median)
      
      ### saving climatology raster as tif
      writeRaster(climatologyrast, here("Data", "Daymet", "CONUS", "Climatology", "FinalTIFFs", paste0(variable, "_Month", as.character(lyr), "MeanClimatology.tif")), overwrite = TRUE)
      # writeRaster(climatologyrast, here("Data", "Daymet", "CONUS", "Climatology", "FinalTIFFs", paste0(variable, "_Month", as.character(lyr), "MedianClimatology.tif")))
    }
  }
}


print("starting PrcpSum climatology calculation for CONUS")
for (variable in c("prcp")){
  print(variable)
  ### creating vector to store file paths
  meanfilepaths <- c()
  medianfilepaths <- c()
  for (daymetyear in 1992:2023){
    currentyear <- as.character(daymetyear)
    variablemeanfiles_temp <- list.files(path = here("Data", "Daymet", "CONUS", "Mosaics", "NetCDFs"), pattern = paste0("^", variable, "_", currentyear, "MonthlySum", "*"), all.files= T, full.names= T)
    # variablemedianfiles_temp <- list.files(path = here("Data", "Daymet", "CONUS", "Mosaics"), pattern = paste0("^", variable, "_", currentyear, "MonthlyMedian", "*"), all.files= T, full.names= T)
    meanfilepaths <- c(meanfilepaths, variablemeanfiles_temp)
    # medianfilepaths <- c(medianfilepaths, variablemedianfiles_temp)
  }
  rast1 <- rast(meanfilepaths[1])
  rast2 <- rast(meanfilepaths[2])
  rast3 <- rast(meanfilepaths[3])
  rast4 <- rast(meanfilepaths[4])
  rast5 <- rast(meanfilepaths[5])
  rast6 <- rast(meanfilepaths[6])
  rast7 <- rast(meanfilepaths[7])
  rast8 <- rast(meanfilepaths[8])
  rast9 <- rast(meanfilepaths[9])
  rast10 <- rast(meanfilepaths[10])
  rast11 <- rast(meanfilepaths[11])
  rast12 <- rast(meanfilepaths[12])
  rast13 <- rast(meanfilepaths[13])
  rast14 <- rast(meanfilepaths[14])
  rast15 <- rast(meanfilepaths[15])
  rast16 <- rast(meanfilepaths[16])
  rast17 <- rast(meanfilepaths[17])
  rast18 <- rast(meanfilepaths[18])
  rast19 <- rast(meanfilepaths[19])
  rast20 <- rast(meanfilepaths[20])
  rast21 <- rast(meanfilepaths[21])
  rast22 <- rast(meanfilepaths[22])
  rast23 <- rast(meanfilepaths[23])
  rast24 <- rast(meanfilepaths[24])
  rast25 <- rast(meanfilepaths[25])
  rast26 <- rast(meanfilepaths[26])
  rast27 <- rast(meanfilepaths[27])
  rast28 <- rast(meanfilepaths[28])
  rast29 <- rast(meanfilepaths[29])
  rast30 <- rast(meanfilepaths[30])
  rast31 <- rast(meanfilepaths[31])
  rast32 <- rast(meanfilepaths[32])
  # rast33 <- rast(meanfilepaths[33])
  
  ### concatenating specific months of each raster to create climatological average over study period
  for (lyr in 1:nlyr(rast1)){
    if (lyr < 10){
      climaterast <- c(rast2[[lyr]], rast3[[lyr]], rast4[[lyr]], rast5[[lyr]], rast6[[lyr]], rast7[[lyr]], rast8[[lyr]], rast9[[lyr]], rast10[[lyr]], rast11[[lyr]], rast12[[lyr]], rast13[[lyr]], rast14[[lyr]], rast15[[lyr]], rast16[[lyr]], rast17[[lyr]], rast18[[lyr]], rast19[[lyr]], rast20[[lyr]], rast21[[lyr]], rast22[[lyr]], rast23[[lyr]], rast24[[lyr]], rast25[[lyr]], rast26[[lyr]], rast27[[lyr]], rast28[[lyr]], rast29[[lyr]], rast30[[lyr]])
      
      ### taking climatological average
      climatologyrast <- app(climaterast, fun = mean)
      # climatologyrast <- app(climaterast, fun = median)
      
      ### saving climatology raster as tif
      writeRaster(climatologyrast, here("Data", "Daymet", "CONUS", "Climatology", "FinalTIFFs", paste0(variable, "Sum_Month", as.character(lyr), "MeanClimatology.tif")))
      # writeRaster(climatologyrast, here("Data", "Daymet", "CONUS", "Climatology", "FinalTIFFs", paste0(variable, "Sum_Month", as.character(lyr), "MedianClimatology.tif")))
    } else {
      climaterast <- c(rast1[[lyr]], rast2[[lyr]], rast3[[lyr]], rast4[[lyr]], rast5[[lyr]], rast6[[lyr]], rast7[[lyr]], rast8[[lyr]], rast9[[lyr]], rast10[[lyr]], rast11[[lyr]], rast12[[lyr]], rast13[[lyr]], rast14[[lyr]], rast15[[lyr]], rast16[[lyr]], rast17[[lyr]], rast18[[lyr]], rast19[[lyr]], rast20[[lyr]], rast21[[lyr]], rast22[[lyr]], rast23[[lyr]], rast24[[lyr]], rast25[[lyr]], rast26[[lyr]], rast27[[lyr]], rast28[[lyr]], rast29[[lyr]])
      
      ### taking climatological average
      climatologyrast <- app(climaterast, fun = mean)
      # climatologyrast <- app(climaterast, fun = median)
      
      ### saving climatology raster as tif
      writeRaster(climatologyrast, here("Data", "Daymet", "CONUS", "Climatology", "FinalTIFFs", paste0(variable, "Sum_Month", as.character(lyr), "MeanClimatology.tif")))
      # writeRaster(climatologyrast, here("Data", "Daymet", "CONUS", "Climatology", "FinalTIFFs", paste0(variable, "Sum_Month", as.character(lyr), "MedianClimatology.tif")))
    }
  }
}


ConusEnd <- Sys.time()



##### making some plots #####
Ecoregions_sf <- read_sf(here("Data", "L3_Ecoregions_USB", "CONUS", "CONUS_AOI.gpkg"))
Ecoregions_bbox <- st_bbox(Ecoregions_sf)

unique(Ecoregions_sf$STATE_NAME)

AK_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "Alaska", "AK_AOI.gpkg"))

SnotelPoints <- read_sf(here("Data", "SNOTEL", "CONUS", "GIS", "SnotelData_CONUS_PeakSWE.gpkg")) |>
  group_by(site_id) |>
  summarise(swe = max(peak_swe))
SnotelPoints_AK <- read_sf(here("Data", "SNOTEL", "Alaska", "GIS", "SnotelData_AK_PeakSWE.gpkg")) |>
  group_by(site_id) |>
  summarise(swe = max(peak_swe))

SnotelPoints_AK_bbox <- st_bbox(SnotelPoints_AK)
SnotelPoints_AK_bbox[[1]]

CONUS_DEM <- rast(here("Data", "DEM", "CONUS", "CONUSDEMMosaic.tif"))
AK_DEM <- rast(here("Data", "DEM", "Alaska", "AKDEMMosaic.tif"))

AK_tminDecClimatology <- rast(here("Data", "Daymet", "Alaska", "Climatology", "FinalTIFFs", "tmin_Month12MeanClimatology.tif"))

CONUS_tminDecClimatology <- rast(here("Data", "Daymet", "CONUS", "Climatology", "FinalTIFFs", "tmin_Month12MeanClimatology.tif"))

### plotting to see what happens
testplot <- ggplot() +
  theme_bw() +
  geom_spatraster(data = CONUS_DEM) +
  scale_fill_wiki_c(na.value = NA) +
  geom_sf(data = Ecoregions_sf, fill = NA, color = "black") +
  geom_sf(data = SnotelPoints, color = 'red')
testplot

testplot2 <- ggplot() +
  theme_bw() +
  geom_spatraster(data = CONUS_tminDecClimatology) +
  scale_fill_wiki_c(na.value = NA) +
  geom_sf(data = Ecoregions_sf, fill = NA, color = "black") +
  geom_sf(data = SnotelPoints, color = 'red')
testplot2



### plotting to see what happens
testplot_AK <- ggplot() +
  theme_bw() +
  geom_spatraster(data = AK_DEM) +
  scale_fill_wiki_c(na.value = NA) +
  geom_sf(data = AK_AOI, fill = NA, color = "black") +
  geom_sf(data = SnotelPoints_AK, color = 'red')
testplot_AK

### plotting to see what happens
testplot_AK2 <- ggplot() +
  theme_bw() +
  geom_spatraster(data = AK_tminDecClimatology) +
  scale_fill_wiki_c(na.value = NA) +
  geom_sf(data = AK_AOI, fill = NA, color = "black") +
  geom_sf(data = SnotelPoints_AK, color = 'red')
testplot_AK2



AlaskaTime <- AlaskaEnd - AlaskaStart
ConusTime <- ConusEnd - ConusStart
EndTime <- Sys.time()
TotalTime <- EndTime - StartTime
AlaskaTime
ConusTime
TotalTime