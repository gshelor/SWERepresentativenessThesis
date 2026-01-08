##### Script for Mosaicing Daymet data #####
### Daymet mosaicing start time
DaymetStartTime <- Sys.time()
library(pacman)
p_load(here, tidyverse, sf, terra, ncdf4, parallel, daymetr, tidyterra)


### Tmean Start time
# TmeanStartTime <- Sys.time()
##### Calculating tmean #####
##### Alaska tmean calculation #####
# print("calculating tmean for Alaska")
# for (daymetyear in 1992:2020){
#   currentyear <- as.character(daymetyear)
#   print(currentyear)
#   ### identifying tiles that are the correct variable and in the same year
#   tminfiles_temp <- list.files(path = here("Data", "Daymet", "Alaska", "Tiles", "tmin_tmax"), pattern = paste0("^", "tmin", "_", currentyear, "_", "*"), all.files= T, full.names= T)
#   tmaxfiles_temp <- list.files(path = here("Data", "Daymet", "Alaska", "Tiles", "tmin_tmax"), pattern = paste0("^", "tmax", "_", currentyear, "_", "*"), all.files= T, full.names= T)
#   for (i in 1:length(tminfiles_temp)){
#     ### grabbing the tile id
#     ## only need to do this once, it's for when I save the tmean netCDF
#     tmin_filename_temp <- tminfiles_temp[i]
#     tmin_fileextension <- str_split_1(tmin_filename_temp, pattern = paste0("tmin_", currentyear, "_"))[[2]]
#     tileid_char <- str_split_1(tmin_fileextension, pattern = ".nc")[[1]]
#     tileid <- as.numeric(tileid_char)
#     ### creating tmean layer and writing netcdf file
#     tmean_temp <- daymet_grid_tmean(path = here("Data", "Daymet", "Alaska", "Tiles", "tmin_tmax"), product = tileid, year = daymetyear, internal = TRUE)
#     writeRaster(tmean_temp, here("Data", "Daymet", "Alaska", "Tiles", "tmean", paste0("tmean_", currentyear, "_", tileid_char, ".tif")), overwrite = TRUE)
#   }
# }
# 
# 
# ##### CONUS tmean calculation #####
# # print("calculating tmean for CONUS")
# # for (daymetyear in 1992:2020){
# #   currentyear <- as.character(daymetyear)
# #   print(currentyear)
# #   ### identifying tiles that are the correct variable and in the same year
# #   tminfiles_temp <- list.files(path = here("Data", "Daymet", "CONUS", "tmin_tmax"), pattern = paste0("^", "tmin", "_", currentyear, "_", "*"), all.files= T, full.names= T)
# #   tmaxfiles_temp <- list.files(path = here("Data", "Daymet", "CONUS", "tmin_tmax"), pattern = paste0("^", "tmax", "_", currentyear, "_", "*"), all.files= T, full.names= T)
# #   for (i in 1:length(tminfiles_temp)){
# #     ### grabbing the tile id
# #     ## only need to do this once, it's for when I save the tmean netCDF
# #     tmin_filename_temp <- tminfiles_temp[i]
# #     tmin_fileextension <- str_split_1(tmin_filename_temp, pattern = paste0("tmin_", currentyear, "_"))[[2]]
# #     tileid_char <- str_split_1(tmin_fileextension, pattern = ".nc")[[1]]
# #     tileid <- as.numeric(tileid_char)
# #     ### creating tmean layer and writing netCDF file
# #     tmean_temp <- daymet_grid_tmean(path = here("Data", "Daymet", "CONUS", "tmin_tmax"), product = tileid, year = daymetyear, internal = TRUE)
# #     writeRaster(tmean_temp, here("Data", "Daymet", "CONUS", "tmean", paste0("tmean_", currentyear, "_", tileid_char, ".tif")), overwrite = TRUE)
# #   }
# # }
# ### tmean end time
# TmeanEndTime <- Sys.time()


### Aggregation Start time
MonthlyAggStartTime <- Sys.time()
##### Aggregating Daymet to get monthly means & sums #####
##### Alaska Daymet monthly aggregate calculation #####
# print("calculating monthly aggregates for Alaska")
# for (daymetyear in 1992:2020){
#   currentyear <- as.character(daymetyear)
#   print(currentyear)
#   ### identifying tiles that are the correct variable and in the same year
#   tminfiles_temp <- list.files(path = here("Data", "Daymet", "Alaska", "Tiles", "tmin_tmax"), pattern = paste0("^", "tmin", "_", currentyear, "_", "*"), all.files= T, full.names= T)
#   tmaxfiles_temp <- list.files(path = here("Data", "Daymet", "Alaska", "Tiles", "tmin_tmax"), pattern = paste0("^", "tmax", "_", currentyear, "_", "*"), all.files= T, full.names= T)
#   prcpfiles_temp <- list.files(path = here("Data", "Daymet", "Alaska", "Tiles", "prcp"), pattern = paste0("^", "prcp", "_", currentyear, "_", "*"), all.files= T, full.names= T)
#   sradfiles_temp <- list.files(path = here("Data", "Daymet", "Alaska", "Tiles", "srad"), pattern = paste0("^", "srad", "_", currentyear, "_", "*"), all.files= T, full.names= T)
#   for (i in 1:length(tminfiles_temp)){
#     ### grabbing the tile id
#     ## only need to do this once, it's for when I save the tmean netCDF
#     tmin_filename_temp <- tminfiles_temp[i]
#     tmin_fileextension <- str_split_1(tmin_filename_temp, pattern = paste0("tmin_", currentyear, "_"))[[2]]
#     tileid_char <- str_split_1(tmin_fileextension, pattern = ".nc")[[1]]
#     ### tmin mean and medians
#     tmin_mean_temp <- daymet_grid_agg(tminfiles_temp[i], int = "monthly", fun = "mean", internal = TRUE)
#     # tmin_median_temp <- daymet_grid_agg(tminfiles_temp[i], int = "monthly", fun = "median", internal = TRUE)
#     ### writing netCDF files of monthly aggregates
#     writeRaster(tmin_mean_temp, here("Data", "Daymet", "Alaska", "Aggregates", paste0("tmin_", currentyear, "_MonthlyMean_", tileid_char, ".tif")), overwrite = TRUE)
#     # writeRaster(tmin_median_temp, here("Data", "Daymet", "Alaska", "Aggregates", paste0("tmin_", currentyear, "_MonthlyMedian_", tileid_char, ".tif")), overwrite = TRUE)
# 
#     ### tmax mean and medians
#     tmax_mean_temp <- daymet_grid_agg(tmaxfiles_temp[i], int = "monthly", fun = "mean", internal = TRUE)
#     # tmax_median_temp <- daymet_grid_agg(tmaxfiles_temp[i], int = "monthly", fun = "median", internal = TRUE)
#     ### writing netCDF files of monthly aggregates
#     writeRaster(tmax_mean_temp, here("Data", "Daymet", "Alaska", "Aggregates", paste0("tmax_", currentyear, "_MonthlyMean_", tileid_char, ".tif")), overwrite = TRUE)
#     # writeRaster(tmax_median_temp, here("Data", "Daymet", "Alaska", "Aggregates", paste0("tmax_", currentyear, "_MonthlyMedian_", tileid_char, ".tif")), overwrite = TRUE)
# 
#     ### prcp mean and medians
#     prcp_mean_temp <- daymet_grid_agg(prcpfiles_temp[i], int = "monthly", fun = "mean", internal = TRUE)
#     # prcp_median_temp <- daymet_grid_agg(prcpfiles_temp[i], int = "monthly", fun = "median", internal = TRUE)
#     prcp_sum_temp <- daymet_grid_agg(prcpfiles_temp[i], int = "monthly", fun = "sum", internal = TRUE)
#     ### writing netCDF files of monthly aggregates
#     writeRaster(prcp_mean_temp, here("Data", "Daymet", "Alaska", "Aggregates", paste0("prcp_", currentyear, "_MonthlyMean_", tileid_char, ".tif")), overwrite = TRUE)
#     # writeRaster(prcp_median_temp, here("Data", "Daymet", "Alaska", "Aggregates", paste0("prcp_", currentyear, "_MonthlyMedian_", tileid_char, ".tif")))
#     writeRaster(prcp_sum_temp, here("Data", "Daymet", "Alaska", "Aggregates", paste0("prcp_", currentyear, "_MonthlySum_", tileid_char, ".tif")), overwrite = TRUE)
# 
#     ### srad mean and medians
#     srad_mean_temp <- daymet_grid_agg(sradfiles_temp[i], int = "monthly", fun = "mean", internal = TRUE)
#     # srad_median_temp <- daymet_grid_agg(sradfiles_temp[i], int = "monthly", fun = "median", internal = TRUE)
#     ### writing netCDF files of monthly aggregates
#     writeRaster(srad_mean_temp, here("Data", "Daymet", "Alaska", "Aggregates", paste0("srad_", currentyear, "_MonthlyMean_", tileid_char, ".tif")), overwrite = TRUE)
#     # writeRaster(srad_median_temp, here("Data", "Daymet", "Alaska", "Aggregates", paste0("srad_", currentyear, "_MonthlyMedian_", tileid_char, ".tif")))
#   }
# }

##### CONUS Daymet monthly aggregation #####
# print("calculating monthly aggregates for CONUS")
# for (daymetyear in 1992:2020){
#   currentyear <- as.character(daymetyear)
#   print(currentyear)
#   ### identifying tiles that are the correct variable and in the same year
#   tminfiles_temp <- list.files(path = here("Data", "Daymet", "CONUS", "tmin_tmax"), pattern = paste0("^", "tmin", "_", currentyear, "_", "*"), all.files= T, full.names= T)
#   tmaxfiles_temp <- list.files(path = here("Data", "Daymet", "CONUS", "tmin_tmax"), pattern = paste0("^", "tmax", "_", currentyear, "_", "*"), all.files= T, full.names= T)
#   prcpfiles_temp <- list.files(path = here("Data", "Daymet", "CONUS", "prcp"), pattern = paste0("^", "prcp", "_", currentyear, "_", "*"), all.files= T, full.names= T)
#   sradfiles_temp <- list.files(path = here("Data", "Daymet", "CONUS", "srad"), pattern = paste0("^", "srad", "_", currentyear, "_", "*"), all.files= T, full.names= T)
#   for (i in 1:length(tminfiles_temp)){
#     ### grabbing the tile id
#     ## only need to do this once, it's for when I save the aggregated netCDF
#     tmin_filename_temp <- tminfiles_temp[i]
#     tmin_fileextension <- str_split_1(tmin_filename_temp, pattern = paste0("tmin_", currentyear, "_"))[[2]]
#     tileid_char <- str_split_1(tmin_fileextension, pattern = ".nc")[[1]]
#     ### tmin mean and medians
#     tmin_mean_temp <- daymet_grid_agg(tminfiles_temp[i], int = "monthly", fun = "mean", internal = TRUE)
#     tmin_median_temp <- daymet_grid_agg(tminfiles_temp[i], int = "monthly", fun = "median", internal = TRUE)
#     ### writing netCDF files of monthly aggregates
#     writeRaster(tmin_mean_temp, here("Data", "Daymet", "CONUS", "Aggregates", paste0("tmin_", currentyear, "_MonthlyMean_", tileid_char, ".tif")), overwrite = TRUE)
#     writeRaster(tmin_median_temp, here("Data", "Daymet", "CONUS", "Aggregates", paste0("tmin_", currentyear, "_MonthlyMedian_", tileid_char, ".tif")), overwrite = TRUE)
#     
#     ### tmax mean and medians
#     tmax_mean_temp <- daymet_grid_agg(tmaxfiles_temp[i], int = "monthly", fun = "mean", internal = TRUE)
#     tmax_median_temp <- daymet_grid_agg(tmaxfiles_temp[i], int = "monthly", fun = "median", internal = TRUE)
#     ### writing netCDF files of monthly aggregates
#     writeRaster(tmax_mean_temp, here("Data", "Daymet", "CONUS", "Aggregates", paste0("tmax_", currentyear, "_MonthlyMean_", tileid_char, ".tif")), overwrite = TRUE)
#     writeRaster(tmax_median_temp, here("Data", "Daymet", "CONUS", "Aggregates", paste0("tmax_", currentyear, "_MonthlyMedian_", tileid_char, ".tif")), overwrite = TRUE)
#     
#     ### prcp mean and medians
#     prcp_mean_temp <- daymet_grid_agg(prcpfiles_temp[i], int = "monthly", fun = "mean", internal = TRUE)
#     prcp_median_temp <- daymet_grid_agg(prcpfiles_temp[i], int = "monthly", fun = "median", internal = TRUE)
#     prcp_sum_temp <- daymet_grid_agg(prcpfiles_temp[i], int = "monthly", fun = "sum", internal = TRUE)
#     ### writing netCDF files of monthly aggregates
#     writeRaster(prcp_mean_temp, here("Data", "Daymet", "CONUS", "Aggregates", paste0("prcp_", currentyear, "_MonthlyMean_", tileid_char, ".tif")), overwrite = TRUE)
#     writeRaster(prcp_median_temp, here("Data", "Daymet", "CONUS", "Aggregates", paste0("prcp_", currentyear, "_MonthlyMedian_", tileid_char, ".tif")), overwrite = TRUE)
#     writeRaster(prcp_sum_temp, here("Data", "Daymet", "CONUS", "Aggregates", paste0("prcp_", currentyear, "_MonthlySum_", tileid_char, ".tif")), overwrite = TRUE)
#     
#     ### srad mean and medians
#     srad_mean_temp <- daymet_grid_agg(sradfiles_temp[i], int = "monthly", fun = "mean", internal = TRUE)
#     srad_median_temp <- daymet_grid_agg(sradfiles_temp[i], int = "monthly", fun = "median", internal = TRUE)
#     ### writing netCDF files of monthly aggregates
# writeRaster(srad_mean_temp, here("Data", "Daymet", "CONUS", "Aggregates", paste0("srad_", currentyear, "_MonthlyMean_", tileid_char, ".tif")), overwrite = TRUE)
# writeRaster(srad_median_temp, here("Data", "Daymet", "CONUS", "Aggregates", paste0("srad_", currentyear, "_MonthlyMedian_", tileid_char, ".tif")), overwrite = TRUE)
#   }
# }


##### Aggregating tmean separately because the daymetr function produced an error #####
##### Alaska tmean monthly aggregation #####
# print("calculating monthly aggregate of tmean for Alaska")
# for (daymetyear in 1992:2020){
#   currentyear <- as.character(daymetyear)
#   print(currentyear)
#   ### identifying tiles that are the correct variable and in the same year
#   tmeanfiles_temp <- list.files(path = here("Data", "Daymet", "Alaska", "Tiles", "tmean"), pattern = paste0("^", "tmean", "_", currentyear, "_", "*"), all.files= T, full.names= T)
# 
#   for (i in 1:length(tmeanfiles_temp)){
#     ### getting tile id for when I save
#     tmean_filename_temp <- tmeanfiles_temp[i]
#     tmean_fileextension <- str_split_1(tmean_filename_temp, pattern = paste0("tmean_", currentyear, "_"))[[2]]
#     tileid_char <- str_split_1(tmean_fileextension, pattern = ".tif")[[1]]
#     ### reading in files one by one as spatrasters
#     tmean_temp <- rast(tmeanfiles_temp[i])
#     ### separating months out individually to calculate a monthly mean
#     if (daymetyear %in% c(1992, 1996, 2000, 2004, 2008, 2012, 2016, 2020)){
#       ### need a special formula for leap years to separate out months
#       ## daymet excludes Dec 31 in leap years
#       tmean_temp_jan <- tmean_temp[[1:31]]
#       tmean_temp_feb <- tmean_temp[[(1+31):(31+29)]]
#       tmean_temp_mar <- tmean_temp[[(1+31+29):(31+29+31)]]
#       tmean_temp_apr <- tmean_temp[[(1+31+29+31):(31+29+31+30)]]
#       tmean_temp_may <- tmean_temp[[(1+31+29+31+30):(31+29+31+30+31)]]
#       tmean_temp_june <- tmean_temp[[(1+31+29+31+30+31):(31+29+31+30+31+30)]]
#       tmean_temp_july <- tmean_temp[[(1+31+29+31+30+31+30):(31+29+31+30+31+30+31)]]
#       tmean_temp_aug <- tmean_temp[[(1+31+29+31+30+31+30+31):(31+29+31+30+31+30+31+31)]]
#       tmean_temp_sep <- tmean_temp[[(1+31+29+31+30+31+30+31+31):(31+29+31+30+31+30+31+31+30)]]
#       tmean_temp_oct <- tmean_temp[[(1+31+29+31+30+31+30+31+31+30):(31+29+31+30+31+30+31+31+30+31)]]
#       tmean_temp_nov <- tmean_temp[[(1+31+29+31+30+31+30+31+31+30+31):(31+29+31+30+31+30+31+31+30+31+30)]]
#       tmean_temp_dec <- tmean_temp[[(1+31+29+31+30+31+30+31+31+30+31+30):365]]
#     } else{
#       tmean_temp_jan <- tmean_temp[[1:31]]
#       tmean_temp_feb <- tmean_temp[[(1+31):(31+28)]]
#       tmean_temp_mar <- tmean_temp[[(1+31+28):(31+28+31)]]
#       tmean_temp_apr <- tmean_temp[[(1+31+28+31):(31+28+31+30)]]
#       tmean_temp_may <- tmean_temp[[(1+31+28+31+30):(31+28+31+30+31)]]
#       tmean_temp_june <- tmean_temp[[(1+31+28+31+30+31):(31+28+31+30+31+30)]]
#       tmean_temp_july <- tmean_temp[[(1+31+28+31+30+31+30):(31+28+31+30+31+30+31)]]
#       tmean_temp_aug <- tmean_temp[[(1+31+28+31+30+31+30+31):(31+28+31+30+31+30+31+31)]]
#       tmean_temp_sep <- tmean_temp[[(1+31+28+31+30+31+30+31+31):(31+28+31+30+31+30+31+31+30)]]
#       tmean_temp_oct <- tmean_temp[[(1+31+28+31+30+31+30+31+31+30):(31+28+31+30+31+30+31+31+30+31)]]
#       tmean_temp_nov <- tmean_temp[[(1+31+28+31+30+31+30+31+31+30+31):(31+28+31+30+31+30+31+31+30+31+30)]]
#       tmean_temp_dec <- tmean_temp[[(1+31+28+31+30+31+30+31+31+30+31+30):365]]
#     }
#     ### performing aggregation
#     tmean_temp_jan_agg_mean <- app(tmean_temp_jan, fun = mean)
#     # tmean_temp_jan_agg_med <- app(tmean_temp_jan, fun = median)
#     tmean_temp_feb_agg_mean <- app(tmean_temp_feb, fun = mean)
#     # tmean_temp_feb_agg_med <- app(tmean_temp_feb, fun = median)
#     tmean_temp_mar_agg_mean <- app(tmean_temp_mar, fun = mean)
#     # tmean_temp_mar_agg_med <- app(tmean_temp_mar, fun = median)
#     tmean_temp_apr_agg_mean <- app(tmean_temp_apr, fun = mean)
#     # tmean_temp_apr_agg_med <- app(tmean_temp_apr, fun = median)
#     tmean_temp_may_agg_mean <- app(tmean_temp_may, fun = mean)
#     # tmean_temp_may_agg_med <- app(tmean_temp_may, fun = median)
#     tmean_temp_june_agg_mean <- app(tmean_temp_june, fun = mean)
#     # tmean_temp_june_agg_med <- app(tmean_temp_june, fun = median)
#     tmean_temp_july_agg_mean <- app(tmean_temp_july, fun = mean)
#     # tmean_temp_july_agg_med <- app(tmean_temp_july, fun = median)
#     tmean_temp_aug_agg_mean <- app(tmean_temp_aug, fun = mean)
#     # tmean_temp_aug_agg_med <- app(tmean_temp_aug, fun = median)
#     tmean_temp_sep_agg_mean <- app(tmean_temp_sep, fun = mean)
#     # tmean_temp_sep_agg_med <- app(tmean_temp_sep, fun = median)
#     tmean_temp_oct_agg_mean <- app(tmean_temp_oct, fun = mean)
#     # tmean_temp_oct_agg_med <- app(tmean_temp_oct, fun = median)
#     tmean_temp_nov_agg_mean <- app(tmean_temp_nov, fun = mean)
#     # tmean_temp_nov_agg_med <- app(tmean_temp_nov, fun = median)
#     tmean_temp_dec_agg_mean <- app(tmean_temp_dec, fun = mean)
#     # tmean_temp_dec_agg_med <- app(tmean_temp_dec, fun = median)
# 
#     ### concatenating aggregated months back together
#     tmean_temp_mean_agg <- c(tmean_temp_jan_agg_mean, tmean_temp_feb_agg_mean, tmean_temp_mar_agg_mean, tmean_temp_apr_agg_mean, tmean_temp_may_agg_mean, tmean_temp_june_agg_mean, tmean_temp_july_agg_mean, tmean_temp_aug_agg_mean, tmean_temp_sep_agg_mean, tmean_temp_oct_agg_mean, tmean_temp_nov_agg_mean, tmean_temp_dec_agg_mean)
#     names(tmean_temp_mean_agg) <- c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")
# 
#     # tmean_temp_med_agg <- c(tmean_temp_jan_agg_med, tmean_temp_feb_agg_med, tmean_temp_mar_agg_med, tmean_temp_apr_agg_med, tmean_temp_may_agg_med, tmean_temp_june_agg_med, tmean_temp_july_agg_med, tmean_temp_aug_agg_med, tmean_temp_sep_agg_med, tmean_temp_oct_agg_med, tmean_temp_nov_agg_med, tmean_temp_dec_agg_med)
#     # names(tmean_temp_med_agg) <- c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")
# 
#     ### writing netCDF files of monthly aggregates
#     writeRaster(tmean_temp_mean_agg, here("Data", "Daymet", "Alaska", "Aggregates", paste0("tmean_", currentyear, "_MonthlyMean_", tileid_char, ".tif")), overwrite = TRUE)
#     # writeRaster(tmean_temp_med_agg, here("Data", "Daymet", "Alaska", "Aggregates", paste0("tmean_", currentyear, "_MonthlyMedian_", tileid_char, ".tif")), overwrite = TRUE)
#   }
# }


##### CONUS tmean monthly aggregation #####
# print("calculating monthly aggregate of tmean for CONUS")
# for (daymetyear in 1992:2020){
#   currentyear <- as.character(daymetyear)
#   print(currentyear)
#   ### identifying tiles that are the correct variable and in the same year
#   tmeanfiles_temp <- list.files(path = here("Data", "Daymet", "CONUS", "tmean"), pattern = paste0("^", "tmean", "_", currentyear, "_", "*"), all.files= T, full.names= T)
#   
#   for (i in 1:length(tmeanfiles_temp)){
#     ### getting tile id for when I save
#     tmean_filename_temp <- tmeanfiles_temp[i]
#     tmean_fileextension <- str_split_1(tmean_filename_temp, pattern = paste0("tmean_", currentyear, "_"))[[2]]
#     tileid_char <- str_split_1(tmean_fileextension, pattern = ".nc")[[1]]
#     ### reading in files one by one as spatrasters
#     tmean_temp <- rast(tmeanfiles_temp[i])
#     ### separating months out individually to calculate a monthly mean
#     if (daymetyear %in% c(1992, 1996, 2000, 2004, 2008, 2012, 2016, 2020)){
#       ### need a special formula for leap years to separate out months
#       ## daymet excludes Dec 31 in leap years
#       tmean_temp_jan <- tmean_temp[[1:31]]
#       tmean_temp_feb <- tmean_temp[[(1+31):(31+29)]]
#       tmean_temp_mar <- tmean_temp[[(1+31+29):(31+29+31)]]
#       tmean_temp_apr <- tmean_temp[[(1+31+29+31):(31+29+31+30)]]
#       tmean_temp_may <- tmean_temp[[(1+31+29+31+30):(31+29+31+30+31)]]
#       tmean_temp_june <- tmean_temp[[(1+31+29+31+30+31):(31+29+31+30+31+30)]]
#       tmean_temp_july <- tmean_temp[[(1+31+29+31+30+31+30):(31+29+31+30+31+30+31)]]
#       tmean_temp_aug <- tmean_temp[[(1+31+29+31+30+31+30+31):(31+29+31+30+31+30+31+31)]]
#       tmean_temp_sep <- tmean_temp[[(1+31+29+31+30+31+30+31+31):(31+29+31+30+31+30+31+31+30)]]
#       tmean_temp_oct <- tmean_temp[[(1+31+29+31+30+31+30+31+31+30):(31+29+31+30+31+30+31+31+30+31)]]
#       tmean_temp_nov <- tmean_temp[[(1+31+29+31+30+31+30+31+31+30+31):(31+29+31+30+31+30+31+31+30+31+30)]]
#       tmean_temp_dec <- tmean_temp[[(1+31+29+31+30+31+30+31+31+30+31+30):365]]
#     } else{
#       tmean_temp_jan <- tmean_temp[[1:31]]
#       tmean_temp_feb <- tmean_temp[[(1+31):(31+28)]]
#       tmean_temp_mar <- tmean_temp[[(1+31+28):(31+28+31)]]
#       tmean_temp_apr <- tmean_temp[[(1+31+28+31):(31+28+31+30)]]
#       tmean_temp_may <- tmean_temp[[(1+31+28+31+30):(31+28+31+30+31)]]
#       tmean_temp_june <- tmean_temp[[(1+31+28+31+30+31):(31+28+31+30+31+30)]]
#       tmean_temp_july <- tmean_temp[[(1+31+28+31+30+31+30):(31+28+31+30+31+30+31)]]
#       tmean_temp_aug <- tmean_temp[[(1+31+28+31+30+31+30+31):(31+28+31+30+31+30+31+31)]]
#       tmean_temp_sep <- tmean_temp[[(1+31+28+31+30+31+30+31+31):(31+28+31+30+31+30+31+31+30)]]
#       tmean_temp_oct <- tmean_temp[[(1+31+28+31+30+31+30+31+31+30):(31+28+31+30+31+30+31+31+30+31)]]
#       tmean_temp_nov <- tmean_temp[[(1+31+28+31+30+31+30+31+31+30+31):(31+28+31+30+31+30+31+31+30+31+30)]]
#       tmean_temp_dec <- tmean_temp[[(1+31+28+31+30+31+30+31+31+30+31+30):365]]
#     }
#     ### performing aggregation
#     tmean_temp_jan_agg_mean <- app(tmean_temp_jan, fun = mean)
#     tmean_temp_jan_agg_med <- app(tmean_temp_jan, fun = median)
#     tmean_temp_feb_agg_mean <- app(tmean_temp_feb, fun = mean)
#     tmean_temp_feb_agg_med <- app(tmean_temp_feb, fun = median)
#     tmean_temp_mar_agg_mean <- app(tmean_temp_mar, fun = mean)
#     tmean_temp_mar_agg_med <- app(tmean_temp_mar, fun = median)
#     tmean_temp_apr_agg_mean <- app(tmean_temp_apr, fun = mean)
#     tmean_temp_apr_agg_med <- app(tmean_temp_apr, fun = median)
#     tmean_temp_may_agg_mean <- app(tmean_temp_may, fun = mean)
#     tmean_temp_may_agg_med <- app(tmean_temp_may, fun = median)
#     tmean_temp_june_agg_mean <- app(tmean_temp_june, fun = mean)
#     tmean_temp_june_agg_med <- app(tmean_temp_june, fun = median)
#     tmean_temp_july_agg_mean <- app(tmean_temp_july, fun = mean)
#     tmean_temp_july_agg_med <- app(tmean_temp_july, fun = median)
#     tmean_temp_aug_agg_mean <- app(tmean_temp_aug, fun = mean)
#     tmean_temp_aug_agg_med <- app(tmean_temp_aug, fun = median)
#     tmean_temp_sep_agg_mean <- app(tmean_temp_sep, fun = mean)
#     tmean_temp_sep_agg_med <- app(tmean_temp_sep, fun = median)
#     tmean_temp_oct_agg_mean <- app(tmean_temp_oct, fun = mean)
#     tmean_temp_oct_agg_med <- app(tmean_temp_oct, fun = median)
#     tmean_temp_nov_agg_mean <- app(tmean_temp_nov, fun = mean)
#     tmean_temp_nov_agg_med <- app(tmean_temp_nov, fun = median)
#     tmean_temp_dec_agg_mean <- app(tmean_temp_dec, fun = mean)
#     tmean_temp_dec_agg_med <- app(tmean_temp_dec, fun = median)
#     
#     ### concatenating aggregated months back together
#     tmean_temp_mean_agg <- c(tmean_temp_jan_agg_mean, tmean_temp_feb_agg_mean, tmean_temp_mar_agg_mean, tmean_temp_apr_agg_mean, tmean_temp_may_agg_mean, tmean_temp_june_agg_mean, tmean_temp_july_agg_mean, tmean_temp_aug_agg_mean, tmean_temp_sep_agg_mean, tmean_temp_oct_agg_mean, tmean_temp_nov_agg_mean, tmean_temp_dec_agg_mean)
#     names(tmean_temp_mean_agg) <- c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")
#     
#     tmean_temp_med_agg <- c(tmean_temp_jan_agg_med, tmean_temp_feb_agg_med, tmean_temp_mar_agg_med, tmean_temp_apr_agg_med, tmean_temp_may_agg_med, tmean_temp_june_agg_med, tmean_temp_july_agg_med, tmean_temp_aug_agg_med, tmean_temp_sep_agg_med, tmean_temp_oct_agg_med, tmean_temp_nov_agg_med, tmean_temp_dec_agg_med)
#     names(tmean_temp_med_agg) <- c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")
#     
#     ### writing netCDF files of monthly aggregates
#     writeRaster(tmean_temp_mean_agg, here("Data", "Daymet", "CONUS", "Aggregates", paste0("tmean_", currentyear, "_MonthlyMean_", tileid_char, ".tif")), overwrite = TRUE)
#     writeRaster(tmean_temp_med_agg, here("Data", "Daymet", "CONUS", "Aggregates", paste0("tmean_", currentyear, "_MonthlyMedian_", tileid_char, ".tif")), overwrite = TRUE)
#   }
# }

### Aggregation end time
MonthlyAggEndTime <- Sys.time()


##### Daymet Mosaicing #####
MosaicStartTime <- Sys.time()
# print("mosaicing monthly aggregates for Alaska")
# ##### Mosaicing monthly aggregated daymet tiles by variable for Alaska #####
# for (daymetyear in 1992:2020){
#   currentyear <- as.character(daymetyear)
#   print(currentyear)
#   for (variable in c("prcp", "tmin", "tmean", "tmax", "srad")){
#     ### identifying tiles that are the correct variable and in the same year
#     if (variable == "prcp"){
#       prcpsumfiles_temp <- list.files(path = here("Data", "Daymet", "Alaska", "Aggregates"), pattern = paste0("^", variable, "_", currentyear, "_MonthlySum", "*"), all.files= T, full.names= T)
#       meanfiles_temp <- list.files(path = here("Data", "Daymet", "Alaska", "Aggregates"), pattern = paste0("^", variable, "_", currentyear, "_MonthlyMean", "*"), all.files= T, full.names= T)
#       # medianfiles_temp <- list.files(path = here("Data", "Daymet", "Alaska", "Aggregates"), pattern = paste0("^", variable, "_", currentyear, "_MonthlyMedian", "*"), all.files= T, full.names= T)
#       ### reading in tiles as spatrasters, making into spatraster collection, writing mosaiced netCDF file to new folder for easier access
#       prcpsum_rasttemp <- lapply(prcpsumfiles_temp, rast)
#       mean_rasttemp <- lapply(meanfiles_temp, rast)
#       # median_rasttemp <- lapply(medianfiles_temp, rast)
# 
#       prcpsum_sprctemp <- sprc(prcpsum_rasttemp)
#       mean_sprctemp <- sprc(mean_rasttemp)
#       # median_sprctemp <- sprc(median_rasttemp)
# 
#       prcpsum_mosaictemp <- mosaic(prcpsum_sprctemp, overwrite = TRUE)
#       mean_mosaictemp <- mosaic(mean_sprctemp, overwrite = TRUE)
#       # median_mosaictemp <- mosaic(median_sprctemp, overwrite = TRUE)
# 
#       writeRaster(prcpsum_mosaictemp, here("Data", "Daymet", "Alaska", "Mosaics", "StackedTIFFs", paste0(variable, "_", currentyear, "MonthlySum_Mosaic.tif")), overwrite = TRUE)
#       writeRaster(mean_mosaictemp, here("Data", "Daymet", "Alaska", "Mosaics", "StackedTIFFs", paste0(variable, "_", currentyear, "MonthlyMean_Mosaic.tif")), overwrite = TRUE)
#       # writeRaster(median_mosaictemp, here("Data", "Daymet", "Alaska", "Mosaics", "StackedTIFFs", paste0(variable, "_", currentyear, "MonthlyMedian_Mosaic.tif")))
#     } else {
#       meanfiles_temp <- list.files(path = here("Data", "Daymet", "Alaska", "Aggregates"), pattern = paste0("^", variable, "_", currentyear, "_MonthlyMean", "*"), all.files= T, full.names= T)
#       # medianfiles_temp <- list.files(path = here("Data", "Daymet", "Alaska", "Aggregates"), pattern = paste0("^", variable, "_", currentyear, "_MonthlyMedian", "*"), all.files= T, full.names= T)
#       ### reading in tiles as spatrasters, making into spatraster collection, writing mosaiced netCDF file to new folder for easier access
#       mean_rasttemp <- lapply(meanfiles_temp, rast)
#       # median_rasttemp <- lapply(medianfiles_temp, rast)
# 
#       mean_sprctemp <- sprc(mean_rasttemp)
#       # median_sprctemp <- sprc(median_rasttemp)
# 
#       mean_mosaictemp <- mosaic(mean_sprctemp, overwrite = TRUE)
#       # median_mosaictemp <- mosaic(median_sprctemp, overwrite = TRUE)
# 
#       writeRaster(mean_mosaictemp, here("Data", "Daymet", "Alaska", "Mosaics", "StackedTIFFs", paste0(variable, "_", currentyear, "MonthlyMean_Mosaic.tif")), overwrite = TRUE)
#       # writeRaster(median_mosaictemp, here("Data", "Daymet", "Alaska", "Mosaics", "StackedTIFFs", paste0(variable, "_", currentyear, "MonthlyMedian_Mosaic.tif")), overwrite = TRUE)
#     }
#   }
# }


##### Mosaicing monthly aggregated daymet tiles by variable for CONUS #####
# print("mosaicing monthly aggregates for CONUS")
# for (daymetyear in 1992:2020){
# currentyear <- as.character(daymetyear)
# print(currentyear)
# for (variable in c("tmean")){
#     ### identifying tiles that are the correct variable and in the same year
#     if (variable == "prcp"){
#       prcpsumfiles_temp <- list.files(path = here("Data", "Daymet", "CONUS", "Aggregates"), pattern = paste0("^", variable, "_", currentyear, "_MonthlySum", "*"), all.files= T, full.names= T)
#       meanfiles_temp <- list.files(path = here("Data", "Daymet", "CONUS", "Aggregates"), pattern = paste0("^", variable, "_", currentyear, "_MonthlyMean", "*"), all.files= T, full.names= T)
#       medianfiles_temp <- list.files(path = here("Data", "Daymet", "CONUS", "Aggregates"), pattern = paste0("^", variable, "_", currentyear, "_MonthlyMedian", "*"), all.files= T, full.names= T)
#       ### reading in tiles as spatrasters, making into spatraster collection, writing mosaiced netCDF file to new folder for easier access
#       prcpsum_rasttemp <- lapply(prcpsumfiles_temp, rast)
#       mean_rasttemp <- lapply(meanfiles_temp, rast)
#       median_rasttemp <- lapply(medianfiles_temp, rast)
# 
#       prcpsum_sprctemp <- sprc(prcpsum_rasttemp)
#       mean_sprctemp <- sprc(mean_rasttemp)
#       median_sprctemp <- sprc(median_rasttemp)
# 
#       prcpsum_mosaictemp <- mosaic(prcpsum_sprctemp, overwrite = TRUE)
#       mean_mosaictemp <- mosaic(mean_sprctemp, overwrite = TRUE)
#       median_mosaictemp <- mosaic(median_sprctemp, overwrite = TRUE)
# 
#       writeRaster(prcpsum_mosaictemp, here("Data", "Daymet", "CONUS", "Mosaics", paste0(variable, "_", currentyear, "MonthlySum_Mosaic.tif")), overwrite = TRUE)
#       writeRaster(mean_mosaictemp, here("Data", "Daymet", "CONUS", "Mosaics", paste0(variable, "_", currentyear, "MonthlyMean_Mosaic.tif")), overwrite = TRUE)
#       writeRaster(median_mosaictemp, here("Data", "Daymet", "CONUS", "Mosaics", paste0(variable, "_", currentyear, "MonthlyMedian_Mosaic.tif")), overwrite = TRUE)
#     } else {
#       meanfiles_temp <- list.files(path = here("Data", "Daymet", "CONUS", "Aggregates"), pattern = paste0("^", variable, "_", currentyear, "_MonthlyMean", "*"), all.files= T, full.names= T)
#       # medianfiles_temp <- list.files(path = here("Data", "Daymet", "CONUS", "Aggregates"), pattern = paste0("^", variable, "_", currentyear, "_MonthlyMedian", "*"), all.files= T, full.names= T)
#       ### reading in tiles as spatrasters, making into spatraster collection, writing mosaiced netCDF file to new folder for easier access
#       mean_rasttemp <- lapply(meanfiles_temp, rast)
#       # median_rasttemp <- lapply(medianfiles_temp, rast)
# 
#       mean_sprctemp <- sprc(mean_rasttemp)
#       # median_sprctemp <- sprc(median_rasttemp)
# 
#       mean_mosaictemp <- mosaic(mean_sprctemp, overwrite = TRUE)
#       # median_mosaictemp <- mosaic(median_sprctemp, overwrite = TRUE)
# 
#       writeRaster(mean_mosaictemp, here("Data", "Daymet", "CONUS", "Mosaics", paste0(variable, "_", currentyear, "MonthlyMean_Mosaic.tif")), overwrite = TRUE)
#       # writeRaster(median_mosaictemp, here("Data", "Daymet", "CONUS", "Mosaics", paste0(variable, "_", currentyear, "MonthlyMedian_Mosaic.tif")), overwrite = TRUE)
#     }
#   }
# }

### Mosaicing End Time
MosaicEndTime <- Sys.time()

##### Creating Climatologies across study period #####
AlaskaClimateStart <- Sys.time()
##### creating Alaska climatologies using mosaics of monthly variable aggregates #####
print("starting climatology calculation for Alaska")
for (variable in c("prcp", "tmin", "tmean", "tmax", "srad")){
  print(variable)
  ### creating vector to store file paths
  meanfilepaths <- c()
  # medianfilepaths <- c()
  for (daymetyear in 1992:2020){
    currentyear <- as.character(daymetyear)
    print(currentyear)
    variablemeanfiles_temp <- list.files(path = here("Data", "Daymet", "Alaska", "Mosaics", "StackedTIFFs"), pattern = paste0("^", variable, "_", currentyear, "MonthlyMean", "*"), all.files= T, full.names= T)
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
  # rast30 <- rast(meanfilepaths[30])
  # rast31 <- rast(meanfilepaths[31])
  # rast32 <- rast(meanfilepaths[32])

  ### concatenating specific months of each raster to create climatological average over study period
  for (lyr in 1:nlyr(rast1)){
    if (lyr < 10){
      climaterast <- c(rast2[[lyr]], rast3[[lyr]], rast4[[lyr]], rast5[[lyr]], rast6[[lyr]], rast7[[lyr]], rast8[[lyr]], rast9[[lyr]], rast10[[lyr]], rast11[[lyr]], rast12[[lyr]], rast13[[lyr]], rast14[[lyr]], rast15[[lyr]], rast16[[lyr]], rast17[[lyr]], rast18[[lyr]], rast19[[lyr]], rast20[[lyr]], rast21[[lyr]], rast22[[lyr]], rast23[[lyr]], rast24[[lyr]], rast25[[lyr]], rast26[[lyr]], rast27[[lyr]], rast28[[lyr]], rast29[[lyr]])

      ### taking climatological average
      climatologyrast <- app(climaterast, fun = mean)
      # climatologyrast <- app(climaterast, fun = median)

      ### saving climatology raster as tif
      writeRaster(climatologyrast, here("Data", "Daymet", "Alaska", "Climatology", "TIFFs", paste0(variable, "_Month", as.character(lyr), "MeanClimatology.tif")), overwrite = TRUE)
      # writeRaster(climatologyrast, here("Data", "Daymet", "Alaska", "Climatology", paste0(variable, "_Month", as.character(lyr), "MedianClimatology.tif")))
    } else {
      climaterast <- c(rast1[[lyr]], rast2[[lyr]], rast3[[lyr]], rast4[[lyr]], rast5[[lyr]], rast6[[lyr]], rast7[[lyr]], rast8[[lyr]], rast9[[lyr]], rast10[[lyr]], rast11[[lyr]], rast12[[lyr]], rast13[[lyr]], rast14[[lyr]], rast15[[lyr]], rast16[[lyr]], rast17[[lyr]], rast18[[lyr]], rast19[[lyr]], rast20[[lyr]], rast21[[lyr]], rast22[[lyr]], rast23[[lyr]], rast24[[lyr]], rast25[[lyr]], rast26[[lyr]], rast27[[lyr]], rast28[[lyr]])

      ### taking climatological average
      climatologyrast <- app(climaterast, fun = mean)
      # climatologyrast <- app(climaterast, fun = median)

      ### saving climatology raster as tif
      writeRaster(climatologyrast, here("Data", "Daymet", "Alaska", "Climatology", "TIFFs", paste0(variable, "_Month", as.character(lyr), "MeanClimatology.tif")), overwrite = TRUE)
      # writeRaster(climatologyrast, here("Data", "Daymet", "Alaska", "Climatology", paste0(variable, "_Month", as.character(lyr), "MedianClimatology.tif")))
    }
  }
}
AlaskaClimateEnd <- Sys.time()



##### creating CONUS climatologies using mosaics of variables #####
# CONUSClimateStart <- Sys.time()
# print("starting climatology calculation for CONUS")
# for (variable in c("tmean")){
#   print(variable)
#   ### creating vector to store file paths
#   meanfilepaths <- c()
#   medianfilepaths <- c()
#   for (daymetyear in 1992:2020){
#     currentyear <- as.character(daymetyear)
#     variablemeanfiles_temp <- list.files(path = here("Data", "Daymet", "CONUS", "Mosaics"), pattern = paste0("^", variable, "_", currentyear, "MonthlyMean", "*"), all.files= T, full.names= T)
#     # variablemedianfiles_temp <- list.files(path = here("Data", "Daymet", "CONUS", "Mosaics"), pattern = paste0("^", variable, "_", currentyear, "MonthlyMedian", "*"), all.files= T, full.names= T)
#     meanfilepaths <- c(meanfilepaths, variablemeanfiles_temp)
#     # medianfilepaths <- c(medianfilepaths, variablemedianfiles_temp)
#   }
#   rast1 <- rast(meanfilepaths[1])
#   rast2 <- rast(meanfilepaths[2])
#   rast3 <- rast(meanfilepaths[3])
#   rast4 <- rast(meanfilepaths[4])
#   rast5 <- rast(meanfilepaths[5])
#   rast6 <- rast(meanfilepaths[6])
#   rast7 <- rast(meanfilepaths[7])
#   rast8 <- rast(meanfilepaths[8])
#   rast9 <- rast(meanfilepaths[9])
#   rast10 <- rast(meanfilepaths[10])
#   rast11 <- rast(meanfilepaths[11])
#   rast12 <- rast(meanfilepaths[12])
#   rast13 <- rast(meanfilepaths[13])
#   rast14 <- rast(meanfilepaths[14])
#   rast15 <- rast(meanfilepaths[15])
#   rast16 <- rast(meanfilepaths[16])
#   rast17 <- rast(meanfilepaths[17])
#   rast18 <- rast(meanfilepaths[18])
#   rast19 <- rast(meanfilepaths[19])
#   rast20 <- rast(meanfilepaths[20])
#   rast21 <- rast(meanfilepaths[21])
#   rast22 <- rast(meanfilepaths[22])
#   rast23 <- rast(meanfilepaths[23])
#   rast24 <- rast(meanfilepaths[24])
#   rast25 <- rast(meanfilepaths[25])
#   rast26 <- rast(meanfilepaths[26])
#   rast27 <- rast(meanfilepaths[27])
#   rast28 <- rast(meanfilepaths[28])
#   rast29 <- rast(meanfilepaths[29])
#   rast30 <- rast(meanfilepaths[30])
#   rast31 <- rast(meanfilepaths[31])
#   rast32 <- rast(meanfilepaths[32])
#   
#   ### concatenating specific months of each raster to create climatological average over study period
#   for (lyr in 1:nlyr(rast1)){
#     if (lyr < 10){
#       climaterast <- c(rast2[[lyr]], rast3[[lyr]], rast4[[lyr]], rast5[[lyr]], rast6[[lyr]], rast7[[lyr]], rast8[[lyr]], rast9[[lyr]], rast10[[lyr]], rast11[[lyr]], rast12[[lyr]], rast13[[lyr]], rast14[[lyr]], rast15[[lyr]], rast16[[lyr]], rast17[[lyr]], rast18[[lyr]], rast19[[lyr]], rast20[[lyr]], rast21[[lyr]], rast22[[lyr]], rast23[[lyr]], rast24[[lyr]], rast25[[lyr]], rast26[[lyr]], rast27[[lyr]], rast28[[lyr]], rast29[[lyr]], rast30[[lyr]], rast31[[lyr]], rast32[[lyr]])
#       
#       ### taking climatological average
#       climatologyrast <- app(climaterast, fun = mean)
#       # climatologyrast <- app(climaterast, fun = median)
#       
#       ### saving climatology raster as tif
#       writeRaster(climatologyrast, here("Data", "Daymet", "CONUS", "Climatology", paste0(variable, "_Month", as.character(lyr), "MeanClimatology.tif")), overwrite = TRUE)
#       # writeRaster(climatologyrast, here("Data", "Daymet", "CONUS", "Climatology", paste0(variable, "_Month", as.character(lyr), "MedianClimatology.tif")))
#     } else {
#       climaterast <- c(rast1[[lyr]], rast2[[lyr]], rast3[[lyr]], rast4[[lyr]], rast5[[lyr]], rast6[[lyr]], rast7[[lyr]], rast8[[lyr]], rast9[[lyr]], rast10[[lyr]], rast11[[lyr]], rast12[[lyr]], rast13[[lyr]], rast14[[lyr]], rast15[[lyr]], rast16[[lyr]], rast17[[lyr]], rast18[[lyr]], rast19[[lyr]], rast20[[lyr]], rast21[[lyr]], rast22[[lyr]], rast23[[lyr]], rast24[[lyr]], rast25[[lyr]], rast26[[lyr]], rast27[[lyr]], rast28[[lyr]], rast29[[lyr]], rast30[[lyr]], rast31[[lyr]])
#       
#       ### taking climatological average
#       climatologyrast <- app(climaterast, fun = mean)
#       # climatologyrast <- app(climaterast, fun = median)
#       
#       ### saving climatology raster as tif
#       writeRaster(climatologyrast, here("Data", "Daymet", "CONUS", "Climatology", paste0(variable, "_Month", as.character(lyr), "MeanClimatology.tif")), overwrite = TRUE)
#       # writeRaster(climatologyrast, here("Data", "Daymet", "CONUS", "Climatology", paste0(variable, "_Month", as.character(lyr), "MedianClimatology.tif")))
#     }
#   }
# }
# CONUSClimateEnd <- Sys.time()


##### making some plots #####
# CONUS_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "CONUS", "CONUS_AOI.gpkg"))

AK_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "Alaska", "AK_AOI.gpkg"))

# SnotelPoints <- read_sf(here("Data", "SNOTEL", "CONUS", "GIS", "SnotelData_CONUS_PeakSWE.gpkg")) |>
#   group_by(site_id) |>
#   summarise(swe = max(peak_swe))
SnotelPoints_AK <- read_sf(here("Data", "SNOTEL", "Alaska", "GIS", "GPKG", "Climatology", "SnotelData_AK_PeakSWEClimatology.gpkg"))

# CONUS_DEM <- rast(here("Data", "DEM", "CONUS", "CONUSDEMMosaic.tif"))
# AK_DEM <- rast(here("Data", "DEM", "Alaska", "AKDEMMosaic.tif"))

AK_tmeanDecClimatology <- rast(here("Data", "Daymet", "Alaska", "Climatology", "TIFFs", "tmean_Month12MeanClimatology.tif"))

# CONUS_tminDecClimatology <- rast(here("Data", "Daymet", "CONUS", "Climatology", "tmin_Month12MeanClimatology.tif"))

### plotting to see what happens
# testplot <- ggplot() +
#   theme_bw() +
#   geom_spatraster(data = CONUS_DEM) +
#   scale_fill_wiki_c(na.value = NA) +
#   geom_sf(data = CONUS_AOI, fill = NA, color = "black") +
#   geom_sf(data = SnotelPoints, color = 'red')

# testplot2 <- ggplot() +
#   theme_bw() +
#   geom_spatraster(data = CONUS_tminDecClimatology) +
#   scale_fill_wiki_c(na.value = NA) +
#   geom_sf(data = CONUS_AOI, fill = NA, color = "black") +
#   geom_sf(data = SnotelPoints, color = 'red')


### plotting to see what happens
# testplot_AK <- ggplot() +
#   theme_bw() +
#   geom_spatraster(data = AK_DEM) +
#   scale_fill_wiki_c(na.value = NA) +
#   geom_sf(data = AK_AOI, fill = NA, color = "black") +
#   geom_sf(data = SnotelPoints_AK, color = 'red')


### plotting to see what happens
testplot_AK2 <- ggplot() +
  theme_bw() +
  geom_spatraster(data = AK_tmeanDecClimatology) +
  scale_fill_wiki_c(na.value = NA) +
  geom_sf(data = AK_AOI, fill = NA, color = "black") +
  geom_sf(data = SnotelPoints_AK, color = 'red')


# testplot
# testplot2
# testplot_AK
testplot_AK2


### Daymet mosaicing end time
# TmeanTime <- TmeanEndTime - TmeanStartTime
# MonthlyAggTime <- MonthlyAggEndTime - MonthlyAggStartTime
# MosaicTime <- MosaicEndTime - MosaicStartTime
AlaskaClimateTime <- AlaskaClimateEnd - AlaskaClimateStart
# CONUSClimateTime <- CONUSClimateEnd - CONUSClimateStart
# TmeanTime
# MonthlyAggTime
# MosaicTime
AlaskaClimateTime
# CONUSClimateTime
DaymetEndTime <- Sys.time()
DaymetTime <- DaymetEndTime - DaymetStartTime
DaymetTime

##### End of Script #####

