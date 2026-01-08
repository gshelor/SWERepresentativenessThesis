##### Script for Daymet resampling #####

### Daymet mosaics need to be resampled to align with land cover

##### loading packages, reading in data to borrow CRS from #####
library(pacman)
p_load(here, tidyverse, terra, sf, daymetr)

### reading in AOIs
AK_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "Alaska", "AK_AOI.gpkg"))
CONUS_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "CONUS", "CONUS_AOI.gpkg"))

### reading in landcovers to use as basis for resampling and cropping
AK_2020LC <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2020.tif"))
CONUS_2020LC <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2020.tif"))

##### Converting CONUS netCDF files to geotiff so I can reproject and crop and resample them in QGIS or ArcPro or python maybe because the terra crop and mask functions didn't work for me for some reason #####
### please god let this work I am sad and have been in my office for many hours begging and pleading for this to get easier instead of harder

### CONUS first now, so I can start modeling it sooner
# CONUS_Mosaic_list <- list.files(path = here("Data", "Daymet", "CONUS", "Mosaics", "NetCDFs"), pattern = paste0("MonthlyMean_Mosaic.nc$"), all.files= T, full.names= T)
# CONUS_Mosaic_list_filenames <- list.files(path = here("Data", "Daymet", "CONUS", "Mosaics", "NetCDFs"), pattern = paste0("MonthlyMean_Mosaic.nc$"), all.files= T)
# CONUS_PrcpSum_Mosaic_list <- list.files(path = here("Data", "Daymet", "CONUS", "Mosaics", "NetCDFs"), pattern = paste0("^.*\\MonthlySum_Mosaic.nc$"), all.files= T, full.names= T)
# CONUS_PrcpSum_Mosaic_list_filenames <- list.files(path = here("Data", "Daymet", "CONUS", "Mosaics", "NetCDFs"), pattern = paste0("^.*\\MonthlySum_Mosaic.nc$"), all.files= T)
# 
# 
# for(file in 1:length(CONUS_Mosaic_list)){
#   temp_rast <- rast(CONUS_Mosaic_list[file])
#   temp_filename <- CONUS_Mosaic_list_filenames[file]
#   ### grabbing the variable name and year for when I save the tif
#   temp_var <- str_split_1(temp_filename, pattern = "_")[[1]]
#   temp_yearMonthlyMean <- str_split_1(temp_filename, pattern = "_")[[2]]
#   temp_year <- str_split_1(temp_yearMonthlyMean, pattern = "MonthlyMean")[[1]]
#   for (x in 1:12){
#     temp_lyr_rast <- temp_rast[[x]]
#     if (x == 1){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "JanMean.tif")), overwrite = TRUE)
#     } else if (x == 2){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_FebMean.tif")), overwrite = TRUE)
#     } else if (x == 3){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_MarMean.tif")), overwrite = TRUE)
#     } else if (x == 4){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_AprMean.tif")), overwrite = TRUE)
#     } else if (x == 5){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_MayMean.tif")), overwrite = TRUE)
#     } else if (x == 6){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_JuneMean.tif")), overwrite = TRUE)
#     } else if (x == 7){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_JulyMean.tif")), overwrite = TRUE)
#     } else if (x == 8){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_AugMean.tif")), overwrite = TRUE)
#     } else if (x == 9){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_SepMean.tif")), overwrite = TRUE)
#     } else if (x == 10){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_OctMean.tif")), overwrite = TRUE)
#     } else if (x == 11){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_NovMean.tif")), overwrite = TRUE)
#     } else{
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_DecMean.tif")), overwrite = TRUE)
#     }
#   }
# }
# 
# ### writing out TIFF file of CONUS Prcp sum data
# for(file in 1:length(CONUS_PrcpSum_Mosaic_list)){
#   temp_rast <- rast(CONUS_PrcpSum_Mosaic_list[file])
#   temp_filename <- CONUS_PrcpSum_Mosaic_list_filenames[file]
#   ### grabbing the variable name and year for when I save the tif
#   temp_var <- str_split_1(temp_filename, pattern = "_")[[1]]
#   temp_yearMonthlyMean <- str_split_1(temp_filename, pattern = "_")[[2]]
#   temp_year <- str_split_1(temp_yearMonthlyMean, pattern = "MonthlySum")[[1]]
#   for (x in 1:12){
#     temp_lyr_rast <- temp_rast[[x]]
#     if (x == 1){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_JanSum.tif")), overwrite = TRUE)
#     } else if (x == 2){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_FebSum.tif")), overwrite = TRUE)
#     } else if (x == 3){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_MarSum.tif")), overwrite = TRUE)
#     } else if (x == 4){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_AprSum.tif")), overwrite = TRUE)
#     } else if (x == 5){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_MaySum.tif")), overwrite = TRUE)
#     } else if (x == 6){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_JuneSum.tif")), overwrite = TRUE)
#     } else if (x == 7){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_JulySum.tif")), overwrite = TRUE)
#     } else if (x == 8){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_AugSum.tif")), overwrite = TRUE)
#     } else if (x == 9){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_SepSum.tif")), overwrite = TRUE)
#     } else if (x == 10){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_OctSum.tif")), overwrite = TRUE)
#     } else if (x == 11){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_NovSum.tif")), overwrite = TRUE)
#     } else{
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_DecSum.tif")), overwrite = TRUE)
#     }
#   }
# }
# 
# 
# 
# 
# 
### making list of files to be reprojected, resampled, and cropped
### also printing lists of just the file names without the full file pathways to use when saving the resampled, cropped rasters
### CONUS
CONUS_climatology_temp_rast_list <- list.files(path = here("Data", "Daymet", "CONUS", "Climatology"), pattern='^t.*\\.tif$', all.files= T, full.names= T)
CONUS_climatology_srad_rast_list <- list.files(path = here("Data", "Daymet", "CONUS", "Climatology"), pattern='^s.*\\.tif$', all.files= T, full.names= T)
CONUS_climatology_prcp_rast_list <- list.files(path = here("Data", "Daymet", "CONUS", "Climatology"), pattern='^p.*\\.tif$', all.files= T, full.names= T)
CONUS_climatology_temp_rast_list_filenames <- list.files(path = here("Data", "Daymet", "CONUS", "Climatology"), pattern='^t.*\\.tif$', all.files= T)
CONUS_climatology_srad_rast_list_filenames <- list.files(path = here("Data", "Daymet", "CONUS", "Climatology"), pattern='^s.*\\.tif$', all.files= T)
CONUS_climatology_prcp_rast_list_filenames <- list.files(path = here("Data", "Daymet", "CONUS", "Climatology"), pattern='^p.*\\.tif$', all.files= T)
### mosaic file lists
CONUS_mosaic_temp_rast_list <- list.files(path = here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs"), pattern='^t.*\\.tif$', all.files= T, full.names= T)
CONUS_mosaic_srad_rast_list <- list.files(path = here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs"), pattern='^s.*\\.tif$', all.files= T, full.names= T)
CONUS_mosaic_prcp_rast_list <- list.files(path = here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs"), pattern='^p.*\\.tif$', all.files= T, full.names= T)
CONUS_mosaic_temp_rast_list_filenames <- list.files(path = here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs"), pattern='^t.*\\.tif$', all.files= T)
CONUS_mosaic_srad_rast_list_filenames <- list.files(path = here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs"), pattern='^s.*\\.tif$', all.files= T)
CONUS_mosaic_prcp_rast_list_filenames <- list.files(path = here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs"), pattern='^p.*\\.tif$', all.files= T)



##### reprojecting, resampling, cropping, and writing out new files for CONUS #####
print("reprojecting and resampling and cropping CONUS temperature climatology tifs")
for (file in 1:length(CONUS_climatology_temp_rast_list)){
  print(file/length(CONUS_climatology_temp_rast_list))
  print(CONUS_climatology_temp_rast_list_filenames[file])
  temp_rast <- rast(CONUS_climatology_temp_rast_list[file])
  temp_rast <- project(temp_rast, crs(CONUS_AOI))
  temp_rast_crop <- crop(temp_rast, CONUS_AOI, mask = TRUE)
  temp_rast_resampled <- resample(temp_rast_crop, CONUS_2020LC, method = "cubicspline", threads = TRUE)
  writeRaster(temp_rast_resampled, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", CONUS_climatology_temp_rast_list_filenames[file]), overwrite = TRUE)
}

print("reprojecting and resampling and cropping CONUS srad climatology tifs")
for (file in 1:length(CONUS_climatology_srad_rast_list)){
  print(file/length(CONUS_climatology_srad_rast_list))
  print(CONUS_climatology_srad_rast_list_filenames[file])
  temp_rast <- rast(CONUS_climatology_srad_rast_list[file])
  temp_rast <- project(temp_rast, crs(CONUS_AOI))
  temp_rast_crop <- crop(temp_rast, CONUS_AOI, mask = TRUE)
  temp_rast_resampled <- resample(temp_rast_crop, CONUS_2020LC, method = "cubicspline", threads = TRUE)
  writeRaster(temp_rast_resampled, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", CONUS_climatology_srad_rast_list_filenames[file]), overwrite = TRUE)
}

print("reprojecting and resampling and cropping CONUS prcp climatology tifs")
for (file in 1:length(CONUS_climatology_prcp_rast_list)){
  print(file/length(CONUS_climatology_prcp_rast_list))
  print(CONUS_climatology_prcp_rast_list_filenames[file])
  temp_rast <- rast(CONUS_climatology_prcp_rast_list[file])
  temp_rast <- project(temp_rast, crs(CONUS_AOI))
  temp_rast_crop <- crop(temp_rast, CONUS_AOI, mask = TRUE)
  temp_rast_resampled <- resample(temp_rast_crop, CONUS_2020LC, method = "cubicspline", threads = TRUE)
  writeRaster(temp_rast_resampled, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", CONUS_climatology_prcp_rast_list_filenames[file]), overwrite = TRUE)
}

print("reprojecting and resampling and cropping CONUS temperature (tmin, tmean, tmax) mosaiced tiffs (Separated by var, year, month)")
for (file in 1:length(CONUS_mosaic_temp_rast_list)){
  print(file/length(CONUS_mosaic_temp_rast_list))
  print(CONUS_mosaic_temp_rast_list_filenames[file])
  temp_rast <- rast(CONUS_mosaic_temp_rast_list[file])
  temp_rast <- project(temp_rast, crs(CONUS_AOI))
  temp_rast_crop <- crop(temp_rast, CONUS_AOI, mask = TRUE)
  temp_rast_resampled <- resample(temp_rast_crop, CONUS_2020LC, method = "cubicspline", threads = TRUE)
  writeRaster(temp_rast_resampled, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", CONUS_mosaic_temp_rast_list_filenames[file]), overwrite = TRUE)
}

print("reprojecting and resampling and cropping srad CONUS mosaiced tiffs (Separated by year, month)")
for (file in 1:length(CONUS_mosaic_srad_rast_list)){
  print(file/length(CONUS_mosaic_srad_rast_list))
  print(CONUS_mosaic_srad_rast_list_filenames[file])
  temp_rast <- rast(CONUS_mosaic_srad_rast_list[file])
  temp_rast <- project(temp_rast, crs(CONUS_AOI))
  temp_rast_crop <- crop(temp_rast, CONUS_AOI, mask = TRUE)
  temp_rast_resampled <- resample(temp_rast_crop, CONUS_2020LC, method = "cubicspline", threads = TRUE)
  writeRaster(temp_rast_resampled, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", CONUS_mosaic_srad_rast_list_filenames[file]), overwrite = TRUE)
}

print("reprojecting and resampling and cropping prcp CONUS mosaiced tiffs (Separated by year, month)")
for (file in 1:length(CONUS_mosaic_prcp_rast_list)){
  print(file/length(CONUS_mosaic_prcp_rast_list))
  print(CONUS_mosaic_prcp_rast_list_filenames[file])
  temp_rast <- rast(CONUS_mosaic_prcp_rast_list[file])
  temp_rast <- project(temp_rast, crs(CONUS_AOI))
  temp_rast_crop <- crop(temp_rast, CONUS_AOI, mask = TRUE)
  temp_rast_resampled <- resample(temp_rast_crop, CONUS_2020LC, method = "cubicspline", threads = TRUE)
  writeRaster(temp_rast_resampled, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", CONUS_mosaic_prcp_rast_list_filenames[file]), overwrite = TRUE)
}


##### Converting Alaska netCDF files to geotiff so I can reproject and crop and resample them in QGIS or ArcPro or python maybe because the terra crop and mask functions didn't work for me for some reason #####
### Alaska
# AK_Mosaic_list <- list.files(path = here("Data", "Daymet", "Alaska", "Mosaics", "StackedTIFFs"), pattern = paste0("^[A-Za-z].*\\MonthlyMean_Mosaic.tif$"), all.files= T, full.names= T)
# AK_Mosaic_list_filenames <- list.files(path = here("Data", "Daymet", "Alaska", "Mosaics", "StackedTIFFs"), pattern = paste0("^[A-Za-z].*\\MonthlyMean_Mosaic.tif$"), all.files= T)
# AK_PrcpSum_Mosaic_list <- list.files(path = here("Data", "Daymet", "Alaska", "Mosaics", "StackedTIFFs"), pattern = paste0("^[A-Za-z].*\\MonthlySum_Mosaic.tif$"), all.files= T, full.names= T)
# AK_PrcpSum_Mosaic_list_filenames <- list.files(path = here("Data", "Daymet", "Alaska", "Mosaics", "StackedTIFFs"), pattern = paste0("^[A-Za-z].*\\MonthlySum_Mosaic.tif$"), all.files= T)

### Monthly mean mosaics
# for(file in 1:length(AK_Mosaic_list)){
#   temp_rast <- rast(AK_Mosaic_list[file])
#   temp_filename <- AK_Mosaic_list_filenames[file]
#   ### grabbing the variable name and year for when I save the tif
#   temp_var <- str_split_1(temp_filename, pattern = "_")[[1]]
#   temp_yearMonthlyMean <- str_split_1(temp_filename, pattern = "_")[[2]]
#   temp_year <- str_split_1(temp_yearMonthlyMean, pattern = "MonthlyMean")[[1]]
#   for (x in 1:12){
#     temp_lyr_rast <- temp_rast[[x]]
#     if (x == 1){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "JanMean.tif")), overwrite = TRUE)
#     } else if (x == 2){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_FebMean.tif")), overwrite = TRUE)
#     } else if (x == 3){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_MarMean.tif")), overwrite = TRUE)
#     } else if (x == 4){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_AprMean.tif")), overwrite = TRUE)
#     } else if (x == 5){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_MayMean.tif")), overwrite = TRUE)
#     } else if (x == 6){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_JuneMean.tif")), overwrite = TRUE)
#     } else if (x == 7){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_JulyMean.tif")), overwrite = TRUE)
#     } else if (x == 8){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_AugMean.tif")), overwrite = TRUE)
#     } else if (x == 9){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_SepMean.tif")), overwrite = TRUE)
#     } else if (x == 10){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_OctMean.tif")), overwrite = TRUE)
#     } else if (x == 11){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_NovMean.tif")), overwrite = TRUE)
#     } else{
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_DecMean.tif")), overwrite = TRUE)
#     }
#   }
# }
# 
# ### converting Prcp Monthly Sum to tiff files
# for(file in 1:length(AK_PrcpSum_Mosaic_list)){
#   temp_rast <- rast(AK_PrcpSum_Mosaic_list[file])
#   temp_filename <- AK_PrcpSum_Mosaic_list_filenames[file]
#   ### grabbing the variable name and year for when I save the tif
#   temp_var <- str_split_1(temp_filename, pattern = "_")[[1]]
#   temp_yearMonthlyMean <- str_split_1(temp_filename, pattern = "_")[[2]]
#   temp_year <- str_split_1(temp_yearMonthlyMean, pattern = "MonthlyMean")[[1]]
#   for (x in 1:12){
#     temp_lyr_rast <- temp_rast[[x]]
#     if (x == 1){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "JanSum.tif")), overwrite = TRUE)
#     } else if (x == 2){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_FebSum.tif")), overwrite = TRUE)
#     } else if (x == 3){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_MarSum.tif")), overwrite = TRUE)
#     } else if (x == 4){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_AprSum.tif")), overwrite = TRUE)
#     } else if (x == 5){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_MaySum.tif")), overwrite = TRUE)
#     } else if (x == 6){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_JuneSum.tif")), overwrite = TRUE)
#     } else if (x == 7){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_JulySum.tif")), overwrite = TRUE)
#     } else if (x == 8){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_AugSum.tif")), overwrite = TRUE)
#     } else if (x == 9){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_SepSum.tif")), overwrite = TRUE)
#     } else if (x == 10){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_OctSum.tif")), overwrite = TRUE)
#     } else if (x == 11){
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_NovSum.tif")), overwrite = TRUE)
#     } else{
#       writeRaster(temp_lyr_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", paste0(temp_var, temp_year, "_DecSum.tif")), overwrite = TRUE)
#     }
#   }
# }


### making list of files to be reprojected, resampled, and cropped
### also printing lists of just the file names without the full file pathways to use when saving the resampled, cropped rasters
### Alaska
AK_climatology_temp_rast_list <- list.files(path = here("Data", "Daymet", "Alaska", "Climatology", "TIFFs"), pattern='^t.*\\.tif$', all.files= T, full.names= T)
AK_climatology_srad_rast_list <- list.files(path = here("Data", "Daymet", "Alaska", "Climatology", "TIFFs"), pattern='^s.*\\.tif$', all.files= T, full.names= T)
AK_climatology_prcp_rast_list <- list.files(path = here("Data", "Daymet", "Alaska", "Climatology", "TIFFs"), pattern='^p.*\\.tif$', all.files= T, full.names= T)
AK_climatology_temp_rast_list_filenames <- list.files(path = here("Data", "Daymet", "Alaska", "Climatology", "TIFFs"), pattern='^t.*\\.tif$', all.files= T)
AK_climatology_srad_rast_list_filenames <- list.files(path = here("Data", "Daymet", "Alaska", "Climatology", "TIFFs"), pattern='^s.*\\.tif$', all.files= T)
AK_climatology_prcp_rast_list_filenames <- list.files(path = here("Data", "Daymet", "Alaska", "Climatology", "TIFFs"), pattern='^p.*\\.tif$', all.files= T)
### mosaic file lists
AK_mosaic_temp_rast_list <- list.files(path = here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs"), pattern='^t.*\\.tif$', all.files= T, full.names= T)
AK_mosaic_srad_rast_list <- list.files(path = here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs"), pattern='^s.*\\.tif$', all.files= T, full.names= T)
AK_mosaic_prcp_rast_list <- list.files(path = here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs"), pattern='^p.*\\.tif$', all.files= T, full.names= T)
AK_mosaic_temp_rast_list_filenames <- list.files(path = here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs"), pattern='^t.*\\.tif$', all.files= T)
AK_mosaic_srad_rast_list_filenames <- list.files(path = here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs"), pattern='^s.*\\.tif$', all.files= T)
AK_mosaic_prcp_rast_list_filenames <- list.files(path = here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs"), pattern='^p.*\\.tif$', all.files= T)



##### reprojecting, resampling, cropping, and writing out new files for Alaska #####
print("reprojecting and resampling and cropping Alaska temperature climatology tifs")
for (file in 1:length(AK_climatology_temp_rast_list)){
  print(file/length(AK_climatology_temp_rast_list))
  print(AK_climatology_temp_rast_list_filenames[file])
  temp_rast <- rast(AK_climatology_temp_rast_list[file])
  temp_rast <- project(temp_rast, crs(AK_AOI))
  temp_rast_crop <- crop(temp_rast, AK_AOI, mask = TRUE)
  temp_rast_resampled <- resample(temp_rast_crop, AK_2020LC, method = "cubicspline", threads = TRUE)
  writeRaster(temp_rast_resampled, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", AK_climatology_temp_rast_list_filenames[file]), overwrite = TRUE)
}

print("reprojecting and resampling and cropping Alaska srad climatology tifs")
for (file in 1:length(AK_climatology_srad_rast_list)){
  print(file/length(AK_climatology_srad_rast_list))
  print(AK_climatology_srad_rast_list_filenames[file])
  temp_rast <- rast(AK_climatology_srad_rast_list[file])
  temp_rast <- project(temp_rast, crs(AK_AOI))
  temp_rast_crop <- crop(temp_rast, AK_AOI, mask = TRUE)
  temp_rast_resampled <- resample(temp_rast_crop, AK_2020LC, method = "cubicspline", threads = TRUE)
  writeRaster(temp_rast_resampled, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", AK_climatology_srad_rast_list_filenames[file]), overwrite = TRUE)
}

print("reprojecting and resampling and cropping Alaska prcp climatology tifs")
for (file in 1:length(AK_climatology_prcp_rast_list)){
  print(file/length(AK_climatology_prcp_rast_list))
  print(AK_climatology_prcp_rast_list_filenames[file])
  temp_rast <- rast(AK_climatology_prcp_rast_list[file])
  temp_rast <- project(temp_rast, crs(AK_AOI))
  temp_rast_crop <- crop(temp_rast, AK_AOI, mask = TRUE)
  temp_rast_resampled <- resample(temp_rast_crop, AK_2020LC, method = "cubicspline", threads = TRUE)
  writeRaster(temp_rast_resampled, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", AK_climatology_prcp_rast_list_filenames[file]), overwrite = TRUE)
}

print("reprojecting and resampling and cropping Alaska temperature (tmin, tmean, tmax) mosaiced tiffs (Separated by var, year, month)")
for (file in 1:length(AK_mosaic_temp_rast_list)){
  print(file/length(AK_mosaic_temp_rast_list))
  print(AK_mosaic_temp_rast_list_filenames[file])
  temp_rast <- rast(AK_mosaic_temp_rast_list[file])
  temp_rast <- project(temp_rast, crs(AK_AOI))
  temp_rast_crop <- crop(temp_rast, AK_AOI, mask = TRUE)
  temp_rast_resampled <- resample(temp_rast_crop, AK_2020LC, method = "cubicspline", threads = TRUE)
  writeRaster(temp_rast_resampled, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", AK_mosaic_temp_rast_list_filenames[file]), overwrite = TRUE)
}

print("reprojecting and resampling and cropping srad Alaska mosaiced tiffs (Separated by year, month)")
for (file in 1:length(AK_mosaic_srad_rast_list)){
  print(file/length(AK_mosaic_srad_rast_list))
  print(AK_mosaic_srad_rast_list_filenames[file])
  temp_rast <- rast(AK_mosaic_srad_rast_list[file])
  temp_rast <- project(temp_rast, crs(AK_AOI))
  temp_rast_crop <- crop(temp_rast, AK_AOI, mask = TRUE)
  temp_rast_resampled <- resample(temp_rast_crop, AK_2020LC, method = "cubicspline", threads = TRUE)
  writeRaster(temp_rast_resampled, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", AK_mosaic_srad_rast_list_filenames[file]), overwrite = TRUE)
}

print("reprojecting and resampling and cropping prcp Alaska mosaiced tiffs (Separated by year, month)")
for (file in 1:length(AK_mosaic_prcp_rast_list)){
  print(file/length(AK_mosaic_prcp_rast_list))
  print(AK_mosaic_prcp_rast_list_filenames[file])
  temp_rast <- rast(AK_mosaic_prcp_rast_list[file])
  temp_rast <- project(temp_rast, crs(AK_AOI))
  temp_rast_crop <- crop(temp_rast, AK_AOI, mask = TRUE)
  temp_rast_resampled <- resample(temp_rast_crop, AK_2020LC, method = "cubicspline", threads = TRUE)
  writeRaster(temp_rast_resampled, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", AK_mosaic_prcp_rast_list_filenames[file]), overwrite = TRUE)
}


