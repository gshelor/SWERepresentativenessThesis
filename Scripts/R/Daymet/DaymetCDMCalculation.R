##### Creating Cooling Degree Month Rasters #####
##### loading packages #####
starttime <- Sys.time()
library(pacman)
p_load(here, tidyverse, sf, terra, parallel)
# options(mc.cores = parallel::detectCores())

##### Calculating CDMs from Climatologies #####
for (month_num in 1:12){
  if (month_num == 1){
    month_char <- "Jan"
  } else if (month_num == 2){
    month_char <- "Feb"
  } else if (month_num == 3){
    month_char <- "Mar"
  } else if (month_num == 4){
    month_char <- "Apr"
  } else if (month_num == 5){
    month_char <- "May"
  } else if (month_num == 6){
    month_char <- "June"
  } else if (month_num == 7){
    month_char <- "July"
  } else if (month_num == 8){
    month_char <- "Aug"
  } else if (month_num == 9){
    month_char <- "Sep"
  } else if (month_num == 10){
    month_char <- "Oct"
  } else if (month_num == 11){
    month_char <- "Nov"
  } else{
    month_char <- "Dec"
  }
  print(paste("Creating Climatology CDMs for CONUS, month:", month_char))
  ### CONUS
  temp_prcpMean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", paste0("prcp_", "Month", month_num, "MeanClimatology.tif")))
  temp_prcpSum_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", paste0("prcpSum_", "Month", month_num, "MeanClimatology.tif")))
  temp_tmean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", paste0("tmean_", "Month", month_num, "MeanClimatology.tif")))

  ### compiling prcp and tmean rasters to one spatraster for CDM creation function below
  temp_prcpMean_tmean_rast <- c(temp_prcpMean_rast, temp_tmean_rast)
  temp_prcpSum_tmean_rast <- c(temp_prcpSum_rast, temp_tmean_rast)

  ### creating CDM
  temp_prcpMean_CDM_rast <- lapp(temp_prcpMean_tmean_rast, fun=function(x, y){x[y >= 10] <- 0; return(x)})
  temp_prcpSum_CDM_rast <- lapp(temp_prcpSum_tmean_rast, fun=function(x, y){x[y >= 10] <- 0; return(x)})
  ## tmean CDM
  ## subtracting 10 first, because otherwise areas that should be 0s result with values of 10
  temp_tmean_rast_sub10 <- 10 - temp_tmean_rast
  ## setting negative values to be 0 (these areas have tmean values above 10)
  temp_tmean_CDM_rast <- app(temp_tmean_rast_sub10, fun = function(x){x[x <= 0] <- 0; return(x)})

  ### writing rasters
  writeRaster(temp_prcpMean_CDM_rast, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", paste0("prcpMean_", month_char, "ClimatologyCDM.tif")), overwrite = TRUE)
  writeRaster(temp_prcpSum_CDM_rast, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", paste0("prcpSum_", month_char, "ClimatologyCDM.tif")), overwrite = TRUE)
  writeRaster(temp_tmean_CDM_rast, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", paste0("tmean_", month_char, "ClimatologyCDM.tif")), overwrite = TRUE)
  
  ### Alaska
  print(paste("Creating Climatology CDMs for Alaska, month:", month_char))
  temp_prcpMean_AK_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", paste0("prcp_", "Month", month_num, "MeanClimatology.tif")))
  temp_prcpSum_AK_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", paste0("prcpSum_", "Month", month_num, "MeanClimatology.tif")))
  temp_tmean_AK_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", paste0("tmean_", "Month", month_num, "MeanClimatology.tif")))
  
  
  ### compiling prcp and tmean rasters to one spatraster for CDM creation functions below
  temp_prcpMean_tmean_AK_rast <- c(temp_prcpMean_AK_rast, temp_tmean_AK_rast)
  temp_prcpSum_tmean_AK_rast <- c(temp_prcpSum_AK_rast, temp_tmean_AK_rast)
  
  ### creating CDM
  temp_prcpMean_AK_CDM_rast <- lapp(temp_prcpMean_tmean_AK_rast, fun=function(x, y){x[y >= 10] <- 0; return(x)})
  temp_prcpSum_AK_CDM_rast <- lapp(temp_prcpSum_tmean_AK_rast, fun=function(x, y){x[y >= 10] <- 0; return(x)})
  ## tmean CDM
  ## subtracting 10 first, because otherwise areas that should be 0s result with values of 10
  temp_tmean_AK_rast_sub10 <- 10 - temp_tmean_AK_rast
  ## setting negative values to be 0 (these areas have tmean values above 10)
  temp_tmean_AK_CDM_rast <- app(temp_tmean_AK_rast_sub10, fun = function(x){x[x <= 0] <- 0; return(x)})
  
  ### writing rasters
  writeRaster(temp_prcpMean_AK_CDM_rast, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", paste0("prcpMean_", month_char, "ClimatologyCDM.tif")), overwrite = TRUE)
  writeRaster(temp_prcpSum_AK_CDM_rast, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", paste0("prcpSum_", month_char, "ClimatologyCDM.tif")), overwrite = TRUE)
  writeRaster(temp_tmean_AK_CDM_rast, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", paste0("tmean_", month_char, "ClimatologyCDM.tif")), overwrite = TRUE)
}



##### Calculating CDMs for individual years #####
for (year_num in 1992:2020){
    for (month_num in 1:12){
      ### setting month as a string for when I save the CDM tif
      if (month_num == 1){
        month_char <- "Jan"
      } else if (month_num == 2){
        month_char <- "Feb"
      } else if (month_num == 3){
        month_char <- "Mar"
      } else if (month_num == 4){
        month_char <- "Apr"
      } else if (month_num == 5){
        month_char <- "May"
      } else if (month_num == 6){
        month_char <- "June"
      } else if (month_num == 7){
        month_char <- "July"
      } else if (month_num == 8){
        month_char <- "Aug"
      } else if (month_num == 9){
        month_char <- "Sep"
      } else if (month_num == 10){
        month_char <- "Oct"
      } else if (month_num == 11){
        month_char <- "Nov"
      } else{
        month_char <- "Dec"
      }
      print(paste("Creating CDMs for CONUS, year:", year_num, "month:", month_char))
      ### reading in CONUS rasters
      ## somehow I fucked up during the resampling of the mean rasters and left out an underscore in the filename. The data is fine, but the file name was intended to follow the same pattern as the other Daymet variables and I messed that one up. oops
      if (month_num == 1){
        temp_prcpMean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", year_num, month_char, "Mean.tif")))
        temp_tmean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmean", year_num, month_char, "Mean.tif")))
      } else{
        temp_prcpMean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", year_num, "_", month_char, "Mean.tif")))
        temp_tmean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmean", year_num, "_", month_char, "Mean.tif")))
      }
      ### reading in prcpSum raster since it always has an underscore
      temp_prcpSum_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", year_num, "_", month_char, "Sum.tif")))


      ### compiling prcp and tmean rasters to one spatraster for CDM creation function below
      temp_prcpMean_tmean_rast <- c(temp_prcpMean_rast, temp_tmean_rast)
      temp_prcpSum_tmean_rast <- c(temp_prcpSum_rast, temp_tmean_rast)

      ### creating CDM
      temp_prcpMean_CDM_rast <- lapp(temp_prcpMean_tmean_rast, fun=function(x, y){x[y >= 10] <- 0; return(x)})
      temp_prcpSum_CDM_rast <- lapp(temp_prcpSum_tmean_rast, fun=function(x, y){x[y >= 10] <- 0; return(x)})
      ## tmean CDM
      ## subtracting 10 first, because otherwise areas that should be 0s result with values of 10
      temp_tmean_rast_sub10 <- 10 - temp_tmean_rast
      ## setting negative values to be 0 (these areas have tmean values above 10)
      temp_tmean_CDM_rast <- app(temp_tmean_rast_sub10, fun = function(x){x[x <= 0] <- 0; return(x)})

      ### writing rasters
      writeRaster(temp_prcpMean_CDM_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", paste0("prcpMean_", month_char, year_num, "CDM.tif")), overwrite = TRUE)
      writeRaster(temp_prcpSum_CDM_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", paste0("prcpSum_", month_char, year_num, "CDM.tif")), overwrite = TRUE)
      writeRaster(temp_tmean_CDM_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", paste0("tmean_", month_char, year_num, "CDM.tif")), overwrite = TRUE)
      
      ### Alaska
      print(paste("Creating CDMs for Alaska, year:", year_num, "month:", month_char))
      ### reading in Alaska rasters
      if (month_num == 1){
        temp_prcpMean_AK_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", year_num, month_char, "Mean.tif")))
        temp_tmean_AK_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmean", year_num, month_char, "Mean.tif")))
        temp_prcpSum_AK_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", year_num, "MonthlySum", month_char, "Sum.tif")))
      } else{
        temp_prcpMean_AK_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", year_num, "_", month_char, "Mean.tif")))
        temp_tmean_AK_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmean", year_num, "_", month_char, "Mean.tif")))
        temp_prcpSum_AK_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", year_num, "MonthlySum_", month_char, "Sum.tif")))
      }
      
      
      
      ### compiling prcp and tmean rasters to one spatraster for CDM creation function below
      temp_prcpMean_tmean_AK_rast <- c(temp_prcpMean_AK_rast, temp_tmean_AK_rast)
      temp_prcpSum_tmean_AK_rast <- c(temp_prcpSum_AK_rast, temp_tmean_AK_rast)
      
      ### creating CDM
      temp_prcpMean_AK_CDM_rast <- lapp(temp_prcpMean_tmean_AK_rast, fun=function(x, y){x[y >= 10] <- 0; return(x)})
      temp_prcpSum_AK_CDM_rast <- lapp(temp_prcpSum_tmean_AK_rast, fun=function(x, y){x[y >= 10] <- 0; return(x)})
      ## tmean CDM
      ## subtracting 10 first, because otherwise areas that should be 0s result with values of 10
      temp_tmean_AK_rast_sub10 <- 10 - temp_tmean_AK_rast
      ## setting negative values to be 0 (these areas have tmean values above 10)
      temp_tmean_AK_CDM_rast <- app(temp_tmean_AK_rast_sub10, fun = function(x){x[x <= 0] <- 0; return(x)})
      
      ### writing rasters
      writeRaster(temp_prcpMean_AK_CDM_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", paste0("prcpMean_", month_char, year_num, "CDM.tif")), overwrite = TRUE)
      writeRaster(temp_prcpSum_AK_CDM_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", paste0("prcpSum_", month_char, year_num, "CDM.tif")), overwrite = TRUE)
      writeRaster(temp_tmean_AK_CDM_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", paste0("tmean_", month_char, year_num, "CDM.tif")), overwrite = TRUE)
  }
}


endtime <- Sys.time()
totaltime <- endtime - starttime
totaltime


