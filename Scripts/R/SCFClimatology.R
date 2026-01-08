##### Snow Cover Frequency Mean Calculation #####
starttime <- Sys.time()
##### loading packages, reading in data #####
library(pacman)
p_load(here, tidyverse, sf, terra, parallel)
# options(mc.cores = parallel::detectCores())

CONUS_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "CONUS", "CONUS_AOI.gpkg"))
### reading in LC raster to serve as basis for reprojections later
LC_CONUS_2020 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2020.tif"))
Snotel_CONUS_Annual <- read_sf(here("Data", "SNOTEL", "CONUS", "GIS", "GPKG", "Annual", "SnotelData_CONUS_PeakSWE.gpkg"))

### listing SCF rasters
# SCF_files <- list.files(path = here("Data", "SCF"), pattern = "*.tif")
# 
# ##### Reprojecting, cropping, and resampling SCF rasters before averaging #####
# for (file in SCF_files){
#   print(file)
#   ### reading in raster
#   temp_scf_rast <- rast(here("Data", "SCF", file))
#   
#   ### reprojecting raster
#   temp_scf_project_rast <- project(temp_scf_rast, LC_CONUS_2020, method = "cubicspline")
#   
#   temp_scf_subst_rast <- subst(temp_scf_project_rast, from = NA, to = -1)
#   
#   temp_scf_crop_rast <- crop(temp_scf_subst_rast, CONUS_AOI, mask = TRUE)
#   
#   writeRaster(temp_scf_crop_rast, here("Data", "SCF", "Resampled", file), overwrite = TRUE)
# }

### 2009, 2018 and also all years before 2003 not available due to either failure of the SnowCloudMetrics GEE code or due to MODIS data not being available because MODIS wasn't launched until 1999
SCF_2003 <- rast(here("Data", "SCF", "Resampled", "SCF2003US.tif"))
plot(SCF_2003, main = "Snow Cover Frequency 2003")
SCF_2004 <- rast(here("Data", "SCF", "Resampled", "SCF2004US.tif"))
plot(SCF_2004, main = "Snow Cover Frequency 2004")
SCF_2005 <- rast(here("Data", "SCF", "Resampled", "SCF2005US.tif"))
plot(SCF_2005, main = "Snow Cover Frequency 2005")
SCF_2006 <- rast(here("Data", "SCF", "Resampled", "SCF2006US.tif"))
plot(SCF_2006, main = "Snow Cover Frequency 2006")
SCF_2007 <- rast(here("Data", "SCF", "Resampled", "SCF2007US.tif"))
plot(SCF_2007, main = "Snow Cover Frequency 2007")
SCF_2008 <- rast(here("Data", "SCF", "Resampled", "SCF2008US.tif"))
plot(SCF_2008, main = "Snow Cover Frequency 2008")
# SCF_2009 <- rast(here("Data", "SCF", "Resampled", "SCF2009US.tif"))
# plot(SCF_2009, main = "Snow Cover Frequency 2009")
SCF_2010 <- rast(here("Data", "SCF", "Resampled", "SCF2010US.tif"))
plot(SCF_2010, main = "Snow Cover Frequency 2010")
SCF_2011 <- rast(here("Data", "SCF", "Resampled", "SCF2011US.tif"))
plot(SCF_2011, main = "Snow Cover Frequency 2011")
SCF_2012 <- rast(here("Data", "SCF", "Resampled", "SCF2012US.tif"))
plot(SCF_2012, main = "Snow Cover Frequency 2012")
SCF_2013 <- rast(here("Data", "SCF", "Resampled", "SCF2013US.tif"))
plot(SCF_2013, main = "Snow Cover Frequency 2013")
SCF_2014 <- rast(here("Data", "SCF", "Resampled", "SCF2014US.tif"))
plot(SCF_2014, main = "Snow Cover Frequency 2014")
SCF_2015 <- rast(here("Data", "SCF", "Resampled", "SCF2015US.tif"))
plot(SCF_2015, main = "Snow Cover Frequency 2015")
SCF_2016 <- rast(here("Data", "SCF", "Resampled", "SCF2016US.tif"))
plot(SCF_2016, main = "Snow Cover Frequency 2016")
SCF_2017 <- rast(here("Data", "SCF", "Resampled", "SCF2017US.tif"))
plot(SCF_2017, main = "Snow Cover Frequency 2017")
SCF_2019 <- rast(here("Data", "SCF", "Resampled", "SCF2019US.tif"))
plot(SCF_2019, main = "Snow Cover Frequency 2019")
SCF_2020 <- rast(here("Data", "SCF", "Resampled", "SCF2020US.tif"))
plot(SCF_2020, main = "Snow Cover Frequency 2020")

SCF_rasts_allyears <- c(SCF_2003, SCF_2004, SCF_2005, SCF_2006, SCF_2007, SCF_2008, SCF_2010, SCF_2011, SCF_2012, SCF_2013, SCF_2014, SCF_2015, SCF_2016, SCF_2017, SCF_2019, SCF_2020)

SCF_Climatology <- app(SCF_rasts_allyears, fun = "max")

plot(SCF_Climatology, main = "Snow Cover Frequency Max Between 2003-2020, excluding 2009 and 2018")

writeRaster(SCF_Climatology, here("Data", "SCF", "Resampled", "SCF2003_2020.tif"), overwrite = TRUE)

SCF_Climatology_NoSnow <- app(SCF_Climatology, fun = function(x){x[x > 0] <- NA; return(x)})
plot(SCF_Climatology_NoSnow, main = "Areas of No Snow Cover in MODIS Between 2003-2020, excluding 2009 and 2018")

writeRaster(SCF_Climatology_NoSnow, here("Data", "SCF", "Resampled", "SCF2003_2020NoSnow.tif"), overwrite = TRUE)


endtime <- Sys.time()
totaltime <- endtime - starttime
totaltime
