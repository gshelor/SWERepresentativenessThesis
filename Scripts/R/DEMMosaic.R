##### Script for Mosaicing DEM Datasets and Deriving Slope and Aspect #####
### Script by Griffin Shelor
### start time
TotalStartTime <- Sys.time()

##### loading packages, reading in AOIs and landcover #####
library(pacman)
p_load(here, tidyverse, sf, terra, parallel, tidyterra)
# options(mc.cores = parallel::detectCores())

### reading in AOIs
CONUS_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "CONUS", "CONUS_AOI.gpkg"))
AK_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "Alaska", "AK_AOI.gpkg"))
st_crs(AK_AOI)
# AK_AOI <- st_transform(AK_AOI, 6393)

### reading in land cover to use for resampling
AK_2020LC <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2020.tif"))
plot(AK_2020LC)
plot(st_geometry(AK_AOI), add = T)
CONUS_2020LC <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "LC_CONUS_2020.tif"))


##### Alaska DEM Mosaicing #####
### making list of all the different tif tiles
# AKrastlist <- list.files(path = here("Data", "DEM", "Alaska", "FABDEM"), pattern='^N.*\\.tif$', all.files= T, full.names= T)
# 
# ### converting all the files to spatraster
# AK_rasts <- lapply(AKrastlist, rast)
# 
# ### making a spatraster collection
# AKdem_sprc <- sprc(AK_rasts)
# 
# ### mosaicing dem_sprc instead of using for loop to merge
# AKdem_rast <- mosaic(AKdem_sprc, fun = "max")
# plot(AKdem_rast, main = "mosaiced alaska dem")
# 
# ### projecting data into Alaska Albers so it can be cropped to AOI
# AKdem_rast_project <- project(AKdem_rast, "epsg:6393")
# 
# plot(AKdem_rast_project, main = "reprojected alaska dem")
# ### for some reason the crop doesn't work? So I'm writing out the reprojected raster and cropping it in QGIS
# writeRaster(AKdem_rast_project, here("Data", "DEM", "Alaska", "AKReprojectedDEM.tif"), overwrite = TRUE)
# AKdem_rast_project <- rast(here("Data", "DEM", "Alaska", "AKReprojectedDEM.tif"))
# AKdem_rast_crop <- terra::crop(AKdem_rast_project, AK_AOI, mask = TRUE)
# 
# plot(AKdem_rast_crop, main = "reprojected and cropped alaska dem")

### this raster cropped in QGIS
AKdem_rast_crop <- rast(here("Data", "DEM", "Alaska", "AKClippedDEM.tif"))
AKdem_rast_resample <- terra::resample(AKdem_rast_crop, AK_2020LC, method = "cubicspline", threads = TRUE)

plot(AKdem_rast_resample, main = "final resampled alaska dem")

### creating slope raster
AKSlope_rast <- terrain(AKdem_rast_resample, v = "slope", neighbors = 8, unit = "degrees")

### creating aspect raster
AKAspect_rast <- terrain(AKdem_rast_resample, v = "aspect", neighbors = 8, unit = "degrees")

### writing dem, slope, and aspect rasters to a tif file
writeRaster(AKdem_rast_resample, here("Data", "DEM", "Alaska", "AKDEMMosaic.tif"), overwrite = TRUE)
writeRaster(AKSlope_rast, here("Data", "DEM", "Alaska", "AKSlope.tif"), overwrite = TRUE)
writeRaster(AKAspect_rast, here("Data", "DEM", "Alaska", "AKAspect.tif"), overwrite = TRUE)



##### CONUS DEM Mosaicing #####
### making lists of all the different tif tiles for each state
CONUSrastlist <- list.files(path = here("Data", "DEM", "CONUS", "Washington"), pattern='.tif$', all.files= T, full.names= T)
CONUSrastlist2 <- list.files(path = here("Data", "DEM", "CONUS", "Oregon"), pattern='.tif$', all.files= T, full.names= T)
CONUSrastlist3 <- list.files(path = here("Data", "DEM", "CONUS", "California"), pattern='.tif$', all.files= T, full.names= T)
CONUSrastlist4 <- list.files(path = here("Data", "DEM", "CONUS", "NewMexico"), pattern='.tif$', all.files= T, full.names= T)
CONUSrastlist5 <- list.files(path = here("Data", "DEM", "CONUS", "Arizona"), pattern='.tif$', all.files= T, full.names= T)
CONUSrastlist6 <- list.files(path = here("Data", "DEM", "CONUS", "Colorado"), pattern='.tif$', all.files= T, full.names= T)
CONUSrastlist7 <- list.files(path = here("Data", "DEM", "CONUS", "Wyoming"), pattern='.tif$', all.files= T, full.names= T)
CONUSrastlist8 <- list.files(path = here("Data", "DEM", "CONUS", "Nevada"), pattern='.tif$', all.files= T, full.names= T)
CONUSrastlist9 <- list.files(path = here("Data", "DEM", "CONUS", "Utah"), pattern='.tif$', all.files= T, full.names= T)
CONUSrastlist10 <- list.files(path = here("Data", "DEM", "CONUS", "Idaho"), pattern='.tif$', all.files= T, full.names= T)
CONUSrastlist11 <- list.files(path = here("Data", "DEM", "CONUS", "Montana"), pattern='.tif$', all.files= T, full.names= T)
CONUSrastlist12 <- list.files(path = here("Data", "DEM", "CONUS", "SouthDakota"), pattern='.tif$', all.files= T, full.names= T)
CONUSrastlist13 <- list.files(path = here("Data", "DEM", "CONUS", "NorthDakota"), pattern='.tif$', all.files= T, full.names= T)
# CONUSrastlist14 <- list.files(path = here("Data", "DEM", "CONUS", "Nebraska"), pattern='.tif$', all.files= T, full.names= T)

### converting all the files to spatraster
CONUS_rasts <- lapply(c(CONUSrastlist, CONUSrastlist2, CONUSrastlist3, CONUSrastlist4, CONUSrastlist5, CONUSrastlist6, CONUSrastlist7, CONUSrastlist8, CONUSrastlist9, CONUSrastlist10, CONUSrastlist11, CONUSrastlist12, CONUSrastlist13), rast)

### making a spatraster collection
CONUSdem_sprc <- sprc(CONUS_rasts)
### mosaicing spatraster collection
CONUSdem_rast <- mosaic(CONUSdem_sprc, fun = "max")

### projecting data into USGS contiguous albers equal area conis (USGS version, it's what the ecoregions used to create the AOI were in) so it can be cropped to AOI
CONUSdem_rast <- project(CONUSdem_rast, crs(CONUS_AOI))

CONUSdem_rast <- terra::crop(CONUSdem_rast, CONUS_AOI, mask = TRUE)
writeRaster(CONUSdem_rast, here("Data", "DEM", "CONUS", "CONUSDEM_Cropped.tif"), overwrite = TRUE)

CONUSdem_rast_resample <- terra::resample(CONUSdem_rast, CONUS_2020LC, method = "cubicspline", threads = TRUE)

### creating slope raster
CONUSSlope_rast <- terrain(CONUSdem_rast_resample, v = "slope", neighbors = 8, unit = "degrees")

### creating aspect raster
CONUSAspect_rast <- terrain(CONUSdem_rast_resample, v = "aspect", neighbors = 8, unit = "degrees")

### writing dem, slope, and aspect rasters to a tif file
writeRaster(CONUSdem_rast_resample, here("Data", "DEM", "CONUS", "CONUSDEMMosaic.tif"), overwrite = TRUE)
writeRaster(CONUSSlope_rast, here("Data", "DEM", "CONUS", "CONUSSlope.tif"), overwrite = TRUE)
writeRaster(CONUSAspect_rast, here("Data", "DEM", "CONUS", "CONUSAspect.tif"), overwrite = TRUE)


### reading in AOIs and a landcover tif to test plots
CONUS_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "CONUS", "CONUS_AOI.gpkg"))
AK_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "Alaska", "AK_AOI.gpkg"))
CONUS_2020LC <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "LC_CONUS_2020.tif"))
AK_2020LC <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2020.tif"))
CONUSdem_rast_resample <- rast(here("Data", "DEM", "CONUS", "CONUSDEMMosaic.tif"))
AKdem_rast_resample <- rast(here("Data", "DEM", "Alaska", "AKDEMMosaic.tif"))
CONUSSlope_rast <- rast(here("Data", "DEM", "CONUS", "CONUSSlope.tif"))
CONUSAspect_rast <- rast(here("Data", "DEM", "CONUS", "CONUSAspect.tif"))
AKSlope_rast <- rast(here("Data", "DEM", "Alaska", "AKSlope.tif"))
AKAspect_rast <- rast(here("Data", "DEM", "Alaska", "AKAspect.tif"))

### plotting to see what happens
dem_plot <- ggplot() +
  theme_bw() +
  geom_spatraster(data = CONUSdem_rast_resample) +
  scale_fill_wiki_c(na.value = NA) +
  geom_sf(data = CONUS_AOI, fill = NA, color = "black")

lc_plot <- ggplot() +
  theme_bw() +
  geom_spatraster(data = as.int(CONUS_2020LC)) +
  scale_fill_wiki_d(na.value = NA) +
  geom_sf(data = CONUS_AOI, fill = NA, color = "black")

slope_plot <- ggplot() +
  theme_bw() +
  geom_spatraster(data = CONUSSlope_rast) +
  scale_fill_wiki_c(na.value = NA) +
  geom_sf(data = CONUS_AOI, fill = NA, color = "black")

aspect_plot <- ggplot() +
  theme_bw() +
  geom_spatraster(data = CONUSAspect_rast) +
  scale_fill_wiki_c(na.value = NA) +
  geom_sf(data = CONUS_AOI, fill = NA, color = "black")


### alaska plots
dem_plot_ak <- ggplot() +
  theme_bw() +
  geom_spatraster(data = AKdem_rast_resample) +
  scale_fill_wiki_c(na.value = NA) +
  geom_sf(data = AK_AOI, fill = NA, color = "black")

lc_plot_ak <- ggplot() +
  theme_bw() +
  geom_spatraster(data = as.int(AK_2020LC)) +
  scale_fill_wiki_d(na.value = NA) +
  geom_sf(data = AK_AOI, fill = NA, color = "black")

slope_plot_ak <- ggplot() +
  theme_bw() +
  geom_spatraster(data = AKSlope_rast) +
  scale_fill_wiki_c(na.value = NA) +
  geom_sf(data = AK_AOI, fill = NA, color = "black")

aspect_plot_ak <- ggplot() +
  theme_bw() +
  geom_spatraster(data = AKAspect_rast) +
  scale_fill_wiki_c(na.value = NA) +
  geom_sf(data = AK_AOI, fill = NA, color = "black")


### calling plots
dem_plot
lc_plot
slope_plot
aspect_plot
dem_plot_ak
lc_plot_ak
slope_plot_ak
aspect_plot_ak

### DEM mosaicing end time
TotalEndTime <- Sys.time()
TotalTime <- TotalEndTime - TotalStartTime
TotalTime
