##### Daymet Cooling Degree Month (CDM) value aggregation script #####
### CDM rasters created in R, aggregated across different time ranges here

starttime <- Sys.time()
##### loading in packages #####
library(pacman)
p_load(here, tidyverse, parallel, sf, terra)
# options(mc.cores = parallel::detectCores())

CONUSstarttime <- Sys.time()
##### reading in rasters for CONUS climatology CDM aggregates #####
### reading in rasters
## prcpMean
temp_sep_prcpMean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "prcpMean_SepClimatologyCDM.tif"))
temp_oct_prcpMean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "prcpMean_OctClimatologyCDM.tif"))
temp_nov_prcpMean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "prcpMean_NovClimatologyCDM.tif"))
temp_dec_prcpMean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "prcpMean_DecClimatologyCDM.tif"))
temp_jan_prcpMean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "prcpMean_JanClimatologyCDM.tif"))
temp_feb_prcpMean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "prcpMean_FebClimatologyCDM.tif"))
temp_mar_prcpMean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "prcpMean_MarClimatologyCDM.tif"))
temp_apr_prcpMean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "prcpMean_AprClimatologyCDM.tif"))
temp_may_prcpMean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "prcpMean_MayClimatologyCDM.tif"))
## prcpSum
temp_sep_prcpSum_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "prcpSum_SepClimatologyCDM.tif"))
temp_oct_prcpSum_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "prcpSum_OctClimatologyCDM.tif"))
temp_nov_prcpSum_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "prcpSum_NovClimatologyCDM.tif"))
temp_dec_prcpSum_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "prcpSum_DecClimatologyCDM.tif"))
temp_jan_prcpSum_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "prcpSum_JanClimatologyCDM.tif"))
temp_feb_prcpSum_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "prcpSum_FebClimatologyCDM.tif"))
temp_mar_prcpSum_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "prcpSum_MarClimatologyCDM.tif"))
temp_apr_prcpSum_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "prcpSum_AprClimatologyCDM.tif"))
temp_may_prcpSum_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "prcpSum_MayClimatologyCDM.tif"))
## tmean
temp_sep_tmean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "tmean_SepClimatologyCDM.tif"))
temp_oct_tmean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "tmean_OctClimatologyCDM.tif"))
temp_nov_tmean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "tmean_NovClimatologyCDM.tif"))
temp_dec_tmean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "tmean_DecClimatologyCDM.tif"))
temp_jan_tmean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "tmean_JanClimatologyCDM.tif"))
temp_feb_tmean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "tmean_FebClimatologyCDM.tif"))
temp_mar_tmean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "tmean_MarClimatologyCDM.tif"))
temp_apr_tmean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "tmean_AprClimatologyCDM.tif"))
temp_may_tmean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "tmean_MayClimatologyCDM.tif"))
## srad
temp_sep_srad_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "srad_Month9MeanClimatology.tif"))
temp_oct_srad_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "srad_Month10MeanClimatology.tif"))
temp_nov_srad_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "srad_Month11MeanClimatology.tif"))
temp_dec_srad_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "srad_Month12MeanClimatology.tif"))
temp_jan_srad_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "srad_Month1MeanClimatology.tif"))
temp_feb_srad_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "srad_Month2MeanClimatology.tif"))
temp_mar_srad_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "srad_Month3MeanClimatology.tif"))
temp_apr_srad_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "srad_Month4MeanClimatology.tif"))
temp_may_srad_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "srad_Month5MeanClimatology.tif"))
## tmin
temp_sep_tmin_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmin_Month9MeanClimatology.tif"))
temp_oct_tmin_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmin_Month10MeanClimatology.tif"))
temp_nov_tmin_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmin_Month11MeanClimatology.tif"))
temp_dec_tmin_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmin_Month12MeanClimatology.tif"))
temp_jan_tmin_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmin_Month1MeanClimatology.tif"))
temp_feb_tmin_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmin_Month2MeanClimatology.tif"))
temp_mar_tmin_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmin_Month3MeanClimatology.tif"))
temp_apr_tmin_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmin_Month4MeanClimatology.tif"))
temp_may_tmin_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmin_Month5MeanClimatology.tif"))
## tmax
temp_sep_tmax_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmax_Month9MeanClimatology.tif"))
temp_oct_tmax_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmax_Month10MeanClimatology.tif"))
temp_nov_tmax_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmax_Month11MeanClimatology.tif"))
temp_dec_tmax_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmax_Month12MeanClimatology.tif"))
temp_jan_tmax_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmax_Month1MeanClimatology.tif"))
temp_feb_tmax_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmax_Month2MeanClimatology.tif"))
temp_mar_tmax_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmax_Month3MeanClimatology.tif"))
temp_apr_tmax_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmax_Month4MeanClimatology.tif"))
temp_may_tmax_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmax_Month5MeanClimatology.tif"))

##### creating climatology CDM aggregates #####
print("Creating climatology CDM aggregates for CONUS")
## full "snow season" aggregate
## prcpMean sum
temp_prcpMean_sep_may_sum_rast <- temp_sep_prcpMean_rast + temp_oct_prcpMean_rast + temp_nov_prcpMean_rast + temp_dec_prcpMean_rast + temp_jan_prcpMean_rast + temp_feb_prcpMean_rast + temp_mar_prcpMean_rast + temp_apr_prcpMean_rast + temp_may_prcpMean_rast
writeRaster(temp_prcpMean_sep_may_sum_rast, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpMean_SepMaySum_ClimatologyCDM.tif"), overwrite = TRUE)
## prcpMean mean
temp_prcpMean_sep_may_mean_rast <- temp_prcpMean_sep_may_sum_rast / 9
writeRaster(temp_prcpMean_sep_may_mean_rast, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpMean_SepMayMean_ClimatologyCDM.tif"), overwrite = TRUE)
## prcpSum sum
temp_prcpSum_sep_may_sum_rast <- temp_sep_prcpSum_rast + temp_oct_prcpSum_rast + temp_nov_prcpSum_rast + temp_dec_prcpSum_rast + temp_jan_prcpSum_rast + temp_feb_prcpSum_rast + temp_mar_prcpSum_rast + temp_apr_prcpSum_rast + temp_may_prcpSum_rast
writeRaster(temp_prcpSum_sep_may_sum_rast, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpSum_SepMaySum_ClimatologyCDM.tif"), overwrite = TRUE)
## prcpSum mean
temp_prcpSum_sep_may_mean_rast <- temp_prcpSum_sep_may_sum_rast / 9
writeRaster(temp_prcpSum_sep_may_mean_rast, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpSum_SepMayMean_ClimatologyCDM.tif"), overwrite = TRUE)
## tmean sum
temp_tmean_sep_may_sum_rast <- temp_sep_tmean_rast + temp_oct_tmean_rast + temp_nov_tmean_rast + temp_dec_tmean_rast + temp_jan_tmean_rast + temp_feb_tmean_rast + temp_mar_tmean_rast + temp_apr_tmean_rast + temp_may_tmean_rast
writeRaster(temp_tmean_sep_may_sum_rast, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmean_SepMaySum_ClimatologyCDM.tif"), overwrite = TRUE)
## tmean mean
temp_tmean_sep_may_mean_rast <- temp_tmean_sep_may_sum_rast / 9
writeRaster(temp_tmean_sep_may_mean_rast, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmean_SepMayMean_ClimatologyCDM.tif"), overwrite = TRUE)
## tmin sum
temp_tmin_sep_may_sum_rast <- temp_sep_tmin_rast + temp_oct_tmin_rast + temp_nov_tmin_rast + temp_dec_tmin_rast + temp_jan_tmin_rast + temp_feb_tmin_rast + temp_mar_tmin_rast + temp_apr_tmin_rast + temp_may_tmin_rast
## tmin mean
temp_tmin_sep_may_mean_rast <- temp_tmin_sep_may_sum_rast / 9
writeRaster(temp_tmin_sep_may_mean_rast, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmin_SepMayMean_ClimatologyCDM.tif"), overwrite = TRUE)
## tmax sum
temp_tmax_sep_may_sum_rast <- temp_sep_tmax_rast + temp_oct_tmax_rast + temp_nov_tmax_rast + temp_dec_tmax_rast + temp_jan_tmax_rast + temp_feb_tmax_rast + temp_mar_tmax_rast + temp_apr_tmax_rast + temp_may_tmax_rast
## tmax mean
temp_tmax_sep_may_mean_rast <- temp_tmax_sep_may_sum_rast / 9
writeRaster(temp_tmax_sep_may_mean_rast, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmax_SepMayMean_ClimatologyCDM.tif"), overwrite = TRUE)
## srad sum
temp_srad_sep_may_sum_rast <- temp_sep_srad_rast + temp_oct_srad_rast + temp_nov_srad_rast + temp_dec_srad_rast + temp_jan_srad_rast + temp_feb_srad_rast + temp_mar_srad_rast + temp_apr_srad_rast + temp_may_srad_rast
## srad mean
temp_srad_sep_may_mean_rast <- temp_srad_sep_may_sum_rast / 9
writeRaster(temp_srad_sep_may_mean_rast, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "srad_SepMayMean_ClimatologyCDM.tif"), overwrite = TRUE)


## full snow season but just from start of water year (October)
## prcpMean sum
temp_prcpMean_oct_may_sum_rast <- temp_oct_prcpMean_rast + temp_nov_prcpMean_rast + temp_dec_prcpMean_rast + temp_jan_prcpMean_rast + temp_feb_prcpMean_rast + temp_mar_prcpMean_rast + temp_apr_prcpMean_rast + temp_may_prcpMean_rast
writeRaster(temp_prcpMean_oct_may_sum_rast, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpMean_OctMaySum_ClimatologyCDM.tif"), overwrite = TRUE)
## prcpMean mean
temp_prcpMean_oct_may_mean_rast <- temp_prcpMean_oct_may_sum_rast / 8
writeRaster(temp_prcpMean_oct_may_mean_rast, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpMean_OctMayMean_ClimatologyCDM.tif"), overwrite = TRUE)
## prcpSum sum
temp_prcpSum_oct_may_sum_rast <- temp_oct_prcpSum_rast + temp_nov_prcpSum_rast + temp_dec_prcpSum_rast + temp_jan_prcpSum_rast + temp_feb_prcpSum_rast + temp_mar_prcpSum_rast + temp_apr_prcpSum_rast + temp_may_prcpSum_rast
writeRaster(temp_prcpSum_oct_may_sum_rast, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpSum_OctMaySum_ClimatologyCDM.tif"), overwrite = TRUE)
## prcpSum mean
temp_prcpSum_oct_may_mean_rast <- temp_prcpSum_oct_may_sum_rast / 8
writeRaster(temp_prcpSum_oct_may_mean_rast, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpSum_OctMayMean_ClimatologyCDM.tif"), overwrite = TRUE)
## tmean sum
temp_tmean_oct_may_sum_rast <- temp_oct_tmean_rast + temp_nov_tmean_rast + temp_dec_tmean_rast + temp_jan_tmean_rast + temp_feb_tmean_rast + temp_mar_tmean_rast + temp_apr_tmean_rast + temp_may_tmean_rast
writeRaster(temp_tmean_oct_may_sum_rast, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmean_OctMaySum_ClimatologyCDM.tif"), overwrite = TRUE)
## tmean mean
temp_tmean_oct_may_mean_rast <- temp_tmean_oct_may_sum_rast / 8
writeRaster(temp_tmean_oct_may_mean_rast, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmean_OctMayMean_ClimatologyCDM.tif"), overwrite = TRUE)
## tmin sum
temp_tmin_oct_may_sum_rast <- temp_oct_tmin_rast + temp_nov_tmin_rast + temp_dec_tmin_rast + temp_jan_tmin_rast + temp_feb_tmin_rast + temp_mar_tmin_rast + temp_apr_tmin_rast + temp_may_tmin_rast
## tmin mean
temp_tmin_oct_may_mean_rast <- temp_tmin_oct_may_sum_rast / 8
writeRaster(temp_tmin_oct_may_mean_rast, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmin_OctMayMean_ClimatologyCDM.tif"), overwrite = TRUE)
## tmax sum
temp_tmax_oct_may_sum_rast <- temp_oct_tmax_rast + temp_nov_tmax_rast + temp_dec_tmax_rast + temp_jan_tmax_rast + temp_feb_tmax_rast + temp_mar_tmax_rast + temp_apr_tmax_rast + temp_may_tmax_rast
## tmax mean
temp_tmax_oct_may_mean_rast <- temp_tmax_oct_may_sum_rast / 8
writeRaster(temp_tmax_oct_may_mean_rast, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmax_OctMayMean_ClimatologyCDM.tif"), overwrite = TRUE)
## srad sum
temp_srad_oct_may_sum_rast <- temp_oct_srad_rast + temp_nov_srad_rast + temp_dec_srad_rast + temp_jan_srad_rast + temp_feb_srad_rast + temp_mar_srad_rast + temp_apr_srad_rast + temp_may_srad_rast
## srad mean
temp_srad_oct_may_mean_rast <- temp_srad_oct_may_sum_rast / 8
writeRaster(temp_srad_oct_may_mean_rast, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "srad_OctMayMean_ClimatologyCDM.tif"), overwrite = TRUE)



## "full" snow season but just from start of water year (October) to April
## prcpMean sum
temp_prcpMean_oct_apr_sum_rast <- temp_oct_prcpMean_rast + temp_nov_prcpMean_rast + temp_dec_prcpMean_rast + temp_jan_prcpMean_rast + temp_feb_prcpMean_rast + temp_mar_prcpMean_rast + temp_apr_prcpMean_rast
writeRaster(temp_prcpMean_oct_apr_sum_rast, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpMean_OctAprSum_ClimatologyCDM.tif"), overwrite = TRUE)
## prcpMean mean
temp_prcpMean_oct_apr_mean_rast <- temp_prcpMean_oct_apr_sum_rast / 7
writeRaster(temp_prcpMean_oct_apr_mean_rast, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpMean_OctAprMean_ClimatologyCDM.tif"), overwrite = TRUE)
## prcpSum sum
temp_prcpSum_oct_apr_sum_rast <- temp_oct_prcpSum_rast + temp_nov_prcpSum_rast + temp_dec_prcpSum_rast + temp_jan_prcpSum_rast + temp_feb_prcpSum_rast + temp_mar_prcpSum_rast + temp_apr_prcpSum_rast
writeRaster(temp_prcpSum_oct_apr_sum_rast, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpSum_OctAprSum_ClimatologyCDM.tif"), overwrite = TRUE)
## prcpSum mean
temp_prcpSum_oct_apr_mean_rast <- temp_prcpSum_oct_apr_sum_rast / 7
writeRaster(temp_prcpSum_oct_apr_mean_rast, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpSum_OctAprMean_ClimatologyCDM.tif"), overwrite = TRUE)
## tmean sum
temp_tmean_oct_apr_sum_rast <- temp_oct_tmean_rast + temp_nov_tmean_rast + temp_dec_tmean_rast + temp_jan_tmean_rast + temp_feb_tmean_rast + temp_mar_tmean_rast + temp_apr_tmean_rast
writeRaster(temp_tmean_oct_apr_sum_rast, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmean_OctAprSum_ClimatologyCDM.tif"), overwrite = TRUE)
## tmean mean
temp_tmean_oct_apr_mean_rast <- temp_tmean_oct_apr_sum_rast / 7
writeRaster(temp_tmean_oct_apr_mean_rast, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmean_OctAprMean_ClimatologyCDM.tif"), overwrite = TRUE)
## tmin sum
temp_tmin_oct_apr_sum_rast <- temp_oct_tmin_rast + temp_nov_tmin_rast + temp_dec_tmin_rast + temp_jan_tmin_rast + temp_feb_tmin_rast + temp_mar_tmin_rast + temp_apr_tmin_rast
## tmin mean
temp_tmin_oct_apr_mean_rast <- temp_tmin_oct_apr_sum_rast / 7
writeRaster(temp_tmin_oct_apr_mean_rast, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmin_OctAprMean_ClimatologyCDM.tif"), overwrite = TRUE)
## tmax sum
temp_tmax_oct_apr_sum_rast <- temp_oct_tmax_rast + temp_nov_tmax_rast + temp_dec_tmax_rast + temp_jan_tmax_rast + temp_feb_tmax_rast + temp_mar_tmax_rast + temp_apr_tmax_rast
## tmax mean
temp_tmax_oct_apr_mean_rast <- temp_tmax_oct_apr_sum_rast / 7
writeRaster(temp_tmax_oct_apr_mean_rast, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmax_OctAprMean_ClimatologyCDM.tif"), overwrite = TRUE)
## srad sum
temp_srad_oct_apr_sum_rast <- temp_oct_srad_rast + temp_nov_srad_rast + temp_dec_srad_rast + temp_jan_srad_rast + temp_feb_srad_rast + temp_mar_srad_rast + temp_apr_srad_rast
## srad mean
temp_srad_oct_apr_mean_rast <- temp_srad_oct_apr_sum_rast / 7
writeRaster(temp_srad_oct_apr_mean_rast, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "srad_OctAprMean_ClimatologyCDM.tif"), overwrite = TRUE)





## Winter month aggregates (Dec-Feb)
## prcpMean sum
temp_prcpMean_dec_feb_sum_rast <- temp_dec_prcpMean_rast + temp_jan_prcpMean_rast + temp_feb_prcpMean_rast
writeRaster(temp_prcpMean_dec_feb_sum_rast, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpMean_DecFebSum_ClimatologyCDM.tif"), overwrite = TRUE)
## prcpMean mean
temp_prcpMean_dec_feb_mean_rast <- temp_prcpMean_dec_feb_sum_rast / 3
writeRaster(temp_prcpMean_dec_feb_mean_rast, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpMean_DecFebMean_ClimatologyCDM.tif"), overwrite = TRUE)
## prcpSum sum
temp_prcpSum_dec_feb_sum_rast <- temp_dec_prcpSum_rast + temp_jan_prcpSum_rast + temp_feb_prcpSum_rast
writeRaster(temp_prcpSum_dec_feb_sum_rast, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpSum_DecFebSum_ClimatologyCDM.tif"), overwrite = TRUE)
## prcpSum mean
temp_prcpSum_dec_feb_mean_rast <- temp_prcpSum_dec_feb_sum_rast / 3
writeRaster(temp_prcpSum_dec_feb_mean_rast, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpSum_DecFebMean_ClimatologyCDM.tif"), overwrite = TRUE)
## tmean sum
temp_tmean_dec_feb_sum_rast <- temp_dec_tmean_rast + temp_jan_tmean_rast + temp_feb_tmean_rast
writeRaster(temp_tmean_dec_feb_sum_rast, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmean_DecFebSum_ClimatologyCDM.tif"), overwrite = TRUE)
## tmean mean
temp_tmean_dec_feb_mean_rast <- temp_tmean_dec_feb_sum_rast / 3
writeRaster(temp_tmean_dec_feb_mean_rast, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmean_DecFebMean_ClimatologyCDM.tif"), overwrite = TRUE)
## tmin sum
temp_tmin_dec_feb_sum_rast <- temp_dec_tmin_rast + temp_jan_tmin_rast + temp_feb_tmin_rast
## tmin mean
temp_tmin_dec_feb_mean_rast <- temp_tmin_dec_feb_sum_rast / 3
writeRaster(temp_tmin_dec_feb_mean_rast, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmin_DecFebMean_ClimatologyCDM.tif"), overwrite = TRUE)
## tmax sum
temp_tmax_dec_feb_sum_rast <- temp_dec_tmax_rast + temp_jan_tmax_rast + temp_feb_tmax_rast
## tmax mean
temp_tmax_dec_feb_mean_rast <- temp_tmax_dec_feb_sum_rast / 3
writeRaster(temp_tmax_dec_feb_mean_rast, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmax_DecFebMean_ClimatologyCDM.tif"), overwrite = TRUE)
## srad sum
temp_srad_dec_feb_sum_rast <- temp_dec_srad_rast + temp_jan_srad_rast + temp_feb_srad_rast
## srad mean
temp_srad_dec_feb_mean_rast <- temp_srad_dec_feb_sum_rast / 3
writeRaster(temp_srad_dec_feb_mean_rast, here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "srad_DecFebMean_ClimatologyCDM.tif"), overwrite = TRUE)

##### annual CDM aggregates #####
for (year_num in 1993:2020){
  print("Creating annual CDM aggregates for CONUS")
  print(year_num)
  ### reading in prcpMean rasters for given year
  temp_sep_prcpMean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("prcpMean_Sep", year_num - 1, "CDM.tif")))
  temp_oct_prcpMean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("prcpMean_Oct", year_num - 1, "CDM.tif")))
  temp_nov_prcpMean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("prcpMean_Nov", year_num - 1, "CDM.tif")))
  temp_dec_prcpMean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("prcpMean_Dec", year_num - 1, "CDM.tif")))
  temp_jan_prcpMean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("prcpMean_Jan", year_num, "CDM.tif")))
  temp_feb_prcpMean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("prcpMean_Feb", year_num, "CDM.tif")))
  temp_mar_prcpMean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("prcpMean_Mar", year_num, "CDM.tif")))
  temp_apr_prcpMean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("prcpMean_Apr", year_num, "CDM.tif")))
  temp_may_prcpMean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("prcpMean_May", year_num, "CDM.tif")))
  ### reading in prcpSum rasters for given year
  temp_sep_prcpSum_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("prcpSum_Sep", year_num - 1, "CDM.tif")))
  temp_oct_prcpSum_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("prcpSum_Oct", year_num - 1, "CDM.tif")))
  temp_nov_prcpSum_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("prcpSum_Nov", year_num - 1, "CDM.tif")))
  temp_dec_prcpSum_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("prcpSum_Dec", year_num - 1, "CDM.tif")))
  temp_jan_prcpSum_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("prcpSum_Jan", year_num, "CDM.tif")))
  temp_feb_prcpSum_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("prcpSum_Feb", year_num, "CDM.tif")))
  temp_mar_prcpSum_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("prcpSum_Mar", year_num, "CDM.tif")))
  temp_apr_prcpSum_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("prcpSum_Apr", year_num, "CDM.tif")))
  temp_may_prcpSum_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("prcpSum_May", year_num, "CDM.tif")))
  ### reading in tmean rasters for given year
  temp_sep_tmean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("tmean_Sep", year_num - 1, "CDM.tif")))
  temp_oct_tmean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("tmean_Oct", year_num - 1, "CDM.tif")))
  temp_nov_tmean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("tmean_Nov", year_num - 1, "CDM.tif")))
  temp_dec_tmean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("tmean_Dec", year_num - 1, "CDM.tif")))
  temp_jan_tmean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("tmean_Jan", year_num, "CDM.tif")))
  temp_feb_tmean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("tmean_Feb", year_num, "CDM.tif")))
  temp_mar_tmean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("tmean_Mar", year_num, "CDM.tif")))
  temp_apr_tmean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("tmean_Apr", year_num, "CDM.tif")))
  temp_may_tmean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("tmean_May", year_num, "CDM.tif")))

  ### reading in tmin rasters for given year
  temp_sep_tmin_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", paste0("tmin", year_num - 1, "_SepMean.tif")))
  temp_oct_tmin_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", paste0("tmin", year_num - 1, "_OctMean.tif")))
  temp_nov_tmin_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", paste0("tmin", year_num - 1, "_NovMean.tif")))
  temp_dec_tmin_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", paste0("tmin", year_num - 1, "_DecMean.tif")))
  temp_jan_tmin_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", paste0("tmin", year_num, "JanMean.tif")))
  temp_feb_tmin_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", paste0("tmin", year_num, "_FebMean.tif")))
  temp_mar_tmin_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", paste0("tmin", year_num, "_MarMean.tif")))
  temp_apr_tmin_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", paste0("tmin", year_num, "_AprMean.tif")))
  temp_may_tmin_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", paste0("tmin", year_num, "_MayMean.tif")))

  ### reading in tmax rasters for given year
  temp_sep_tmax_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", paste0("tmax", year_num - 1, "_SepMean.tif")))
  temp_oct_tmax_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", paste0("tmax", year_num - 1, "_OctMean.tif")))
  temp_nov_tmax_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", paste0("tmax", year_num - 1, "_NovMean.tif")))
  temp_dec_tmax_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", paste0("tmax", year_num - 1, "_DecMean.tif")))
  temp_jan_tmax_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", paste0("tmax", year_num, "JanMean.tif")))
  temp_feb_tmax_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", paste0("tmax", year_num, "_FebMean.tif")))
  temp_mar_tmax_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", paste0("tmax", year_num, "_MarMean.tif")))
  temp_apr_tmax_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", paste0("tmax", year_num, "_AprMean.tif")))
  temp_may_tmax_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", paste0("tmax", year_num, "_MayMean.tif")))
  ### reading in srad rasters for given year
  temp_sep_srad_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", paste0("srad", year_num - 1, "_SepMean.tif")))
  temp_oct_srad_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", paste0("srad", year_num - 1, "_OctMean.tif")))
  temp_nov_srad_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", paste0("srad", year_num - 1, "_NovMean.tif")))
  temp_dec_srad_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", paste0("srad", year_num - 1, "_DecMean.tif")))
  temp_jan_srad_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", paste0("srad", year_num, "JanMean.tif")))
  temp_feb_srad_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", paste0("srad", year_num, "_FebMean.tif")))
  temp_mar_srad_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", paste0("srad", year_num, "_MarMean.tif")))
  temp_apr_srad_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", paste0("srad", year_num, "_AprMean.tif")))
  temp_may_srad_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics",  "TIFFs", "ClippedResampled", paste0("srad", year_num, "_MayMean.tif")))


  ### creating CDM aggregates using WY-specific rasters
  ##### creating annual CDM aggregates #####
  ## full "snow season" aggregate
  ## prcpMean sum
  temp_prcpMean_sep_may_sum_rast <- temp_sep_prcpMean_rast + temp_oct_prcpMean_rast + temp_nov_prcpMean_rast + temp_dec_prcpMean_rast + temp_jan_prcpMean_rast + temp_feb_prcpMean_rast + temp_mar_prcpMean_rast + temp_apr_prcpMean_rast + temp_may_prcpMean_rast
  writeRaster(temp_prcpMean_sep_may_sum_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpMean_SepMay", year_num, "CDMSum.tif")), overwrite = TRUE)
  ## prcpMean mean
  temp_prcpMean_sep_may_mean_rast <- temp_prcpMean_sep_may_sum_rast / 9
  writeRaster(temp_prcpMean_sep_may_mean_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpMean_SepMay", year_num, "CDMMean.tif")), overwrite = TRUE)
  ## prcpSum sum
  temp_prcpSum_sep_may_sum_rast <- temp_sep_prcpSum_rast + temp_oct_prcpSum_rast + temp_nov_prcpSum_rast + temp_dec_prcpSum_rast + temp_jan_prcpSum_rast + temp_feb_prcpSum_rast + temp_mar_prcpSum_rast + temp_apr_prcpSum_rast + temp_may_prcpSum_rast
  writeRaster(temp_prcpSum_sep_may_sum_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpSum_SepMay", year_num, "CDMSum.tif")), overwrite = TRUE)
  ## prcpSum mean
  temp_prcpSum_sep_may_mean_rast <- temp_prcpSum_sep_may_sum_rast / 9
  writeRaster(temp_prcpSum_sep_may_mean_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpSum_SepMay", year_num, "CDMMean.tif")), overwrite = TRUE)
  ## tmean sum
  temp_tmean_sep_may_sum_rast <- temp_sep_tmean_rast + temp_oct_tmean_rast + temp_nov_tmean_rast + temp_dec_tmean_rast + temp_jan_tmean_rast + temp_feb_tmean_rast + temp_mar_tmean_rast + temp_apr_tmean_rast + temp_may_tmean_rast
  writeRaster(temp_tmean_sep_may_sum_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmean_SepMay", year_num, "CDMSum.tif")), overwrite = TRUE)
  ## tmean mean
  temp_tmean_sep_may_mean_rast <- temp_tmean_sep_may_sum_rast / 9
  writeRaster(temp_tmean_sep_may_mean_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmean_SepMay", year_num, "CDMMean.tif")), overwrite = TRUE)
  ## tmin sum
  temp_tmin_sep_may_sum_rast <- temp_sep_tmin_rast + temp_oct_tmin_rast + temp_nov_tmin_rast + temp_dec_tmin_rast + temp_jan_tmin_rast + temp_feb_tmin_rast + temp_mar_tmin_rast + temp_apr_tmin_rast + temp_may_tmin_rast
  ## tmin mean
  temp_tmin_sep_may_mean_rast <- temp_tmin_sep_may_sum_rast / 9
  writeRaster(temp_tmin_sep_may_mean_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmin_SepMay", year_num, "Mean.tif")), overwrite = TRUE)
  ## tmax sum
  temp_tmax_sep_may_sum_rast <- temp_sep_tmax_rast + temp_oct_tmax_rast + temp_nov_tmax_rast + temp_dec_tmax_rast + temp_jan_tmax_rast + temp_feb_tmax_rast + temp_mar_tmax_rast + temp_apr_tmax_rast + temp_may_tmax_rast
  ## tmax mean
  temp_tmax_sep_may_mean_rast <- temp_tmax_sep_may_sum_rast / 9
  writeRaster(temp_tmax_sep_may_mean_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmax_SepMay", year_num, "Mean.tif")), overwrite = TRUE)
  ## srad sum
  temp_srad_sep_may_sum_rast <- temp_sep_srad_rast + temp_oct_srad_rast + temp_nov_srad_rast + temp_dec_srad_rast + temp_jan_srad_rast + temp_feb_srad_rast + temp_mar_srad_rast + temp_apr_srad_rast + temp_may_srad_rast
  ## srad mean
  temp_srad_sep_may_mean_rast <- temp_srad_sep_may_sum_rast / 9
  writeRaster(temp_srad_sep_may_mean_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("srad_SepMay", year_num, "Mean.tif")), overwrite = TRUE)


  ## full snow season but just from start of water year (October)
  ## prcpMean sum
  temp_prcpMean_oct_may_sum_rast <- temp_oct_prcpMean_rast + temp_nov_prcpMean_rast + temp_dec_prcpMean_rast + temp_jan_prcpMean_rast + temp_feb_prcpMean_rast + temp_mar_prcpMean_rast + temp_apr_prcpMean_rast + temp_may_prcpMean_rast
  writeRaster(temp_prcpMean_oct_may_sum_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpMean_OctMay", year_num, "CDMSum.tif")), overwrite = TRUE)
  ## prcpMean mean
  temp_prcpMean_oct_may_mean_rast <- temp_prcpMean_oct_may_sum_rast / 8
  writeRaster(temp_prcpMean_oct_may_mean_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpMean_OctMay", year_num, "CDMMean.tif")), overwrite = TRUE)
  ## prcpSum sum
  temp_prcpSum_oct_may_sum_rast <- temp_oct_prcpSum_rast + temp_nov_prcpSum_rast + temp_dec_prcpSum_rast + temp_jan_prcpSum_rast + temp_feb_prcpSum_rast + temp_mar_prcpSum_rast + temp_apr_prcpSum_rast + temp_may_prcpSum_rast
  writeRaster(temp_prcpSum_oct_may_sum_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpSum_OctMay", year_num, "CDMSum.tif")), overwrite = TRUE)
  ## prcpSum mean
  temp_prcpSum_oct_may_mean_rast <- temp_prcpSum_oct_may_sum_rast / 8
  writeRaster(temp_prcpSum_oct_may_mean_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpSum_OctMay", year_num, "CDMMean.tif")), overwrite = TRUE)
  ## tmean sum
  temp_tmean_oct_may_sum_rast <- temp_oct_tmean_rast + temp_nov_tmean_rast + temp_dec_tmean_rast + temp_jan_tmean_rast + temp_feb_tmean_rast + temp_mar_tmean_rast + temp_apr_tmean_rast + temp_may_tmean_rast
  writeRaster(temp_tmean_oct_may_sum_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmean_OctMay", year_num, "CDMSum.tif")), overwrite = TRUE)
  ## tmean mean
  temp_tmean_oct_may_mean_rast <- temp_tmean_oct_may_sum_rast / 8
  writeRaster(temp_tmean_oct_may_mean_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmean_OctMay", year_num, "CDMMean.tif")), overwrite = TRUE)
  ## tmin sum
  temp_tmin_oct_may_sum_rast <- temp_oct_tmin_rast + temp_nov_tmin_rast + temp_dec_tmin_rast + temp_jan_tmin_rast + temp_feb_tmin_rast + temp_mar_tmin_rast + temp_apr_tmin_rast + temp_may_tmin_rast
  ## tmin mean
  temp_tmin_oct_may_mean_rast <- temp_tmin_oct_may_sum_rast / 8
  writeRaster(temp_tmin_oct_may_mean_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmin_OctMay", year_num, "Mean.tif")), overwrite = TRUE)
  ## tmax sum
  temp_tmax_oct_may_sum_rast <- temp_oct_tmax_rast + temp_nov_tmax_rast + temp_dec_tmax_rast + temp_jan_tmax_rast + temp_feb_tmax_rast + temp_mar_tmax_rast + temp_apr_tmax_rast + temp_may_tmax_rast
  ## tmax mean
  temp_tmax_oct_may_mean_rast <- temp_tmax_oct_may_sum_rast / 8
  writeRaster(temp_tmax_oct_may_mean_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmax_OctMay", year_num, "Mean.tif")), overwrite = TRUE)
  ## srad sum
  temp_srad_oct_may_sum_rast <- temp_oct_srad_rast + temp_nov_srad_rast + temp_dec_srad_rast + temp_jan_srad_rast + temp_feb_srad_rast + temp_mar_srad_rast + temp_apr_srad_rast + temp_may_srad_rast
  ## srad mean
  temp_srad_oct_may_mean_rast <- temp_srad_oct_may_sum_rast / 8
  writeRaster(temp_srad_oct_may_mean_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("srad_OctMay", year_num, "Mean.tif")), overwrite = TRUE)



  ## "full" snow season but just from start of water year (October) to April
  ## prcpMean sum
  temp_prcpMean_oct_apr_sum_rast <- temp_oct_prcpMean_rast + temp_nov_prcpMean_rast + temp_dec_prcpMean_rast + temp_jan_prcpMean_rast + temp_feb_prcpMean_rast + temp_mar_prcpMean_rast + temp_apr_prcpMean_rast
  writeRaster(temp_prcpMean_oct_apr_sum_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpMean_OctApr", year_num, "CDMSum.tif")), overwrite = TRUE)
  ## prcpMean mean
  temp_prcpMean_oct_apr_mean_rast <- temp_prcpMean_oct_apr_sum_rast / 7
  writeRaster(temp_prcpMean_oct_apr_mean_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpMean_OctApr", year_num, "CDMMean.tif")), overwrite = TRUE)
  ## prcpSum sum
  temp_prcpSum_oct_apr_sum_rast <- temp_oct_prcpSum_rast + temp_nov_prcpSum_rast + temp_dec_prcpSum_rast + temp_jan_prcpSum_rast + temp_feb_prcpSum_rast + temp_mar_prcpSum_rast + temp_apr_prcpSum_rast
  writeRaster(temp_prcpSum_oct_apr_sum_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpSum_OctApr", year_num, "CDMSum.tif")), overwrite = TRUE)
  ## prcpSum mean
  temp_prcpSum_oct_apr_mean_rast <- temp_prcpSum_oct_apr_sum_rast / 7
  writeRaster(temp_prcpSum_oct_apr_mean_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpSum_OctApr", year_num, "CDMMean.tif")), overwrite = TRUE)
  ## tmean sum
  temp_tmean_oct_apr_sum_rast <- temp_oct_tmean_rast + temp_nov_tmean_rast + temp_dec_tmean_rast + temp_jan_tmean_rast + temp_feb_tmean_rast + temp_mar_tmean_rast + temp_apr_tmean_rast
  writeRaster(temp_tmean_oct_apr_sum_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmean_OctApr", year_num, "CDMSum.tif")), overwrite = TRUE)
  ## tmean mean
  temp_tmean_oct_apr_mean_rast <- temp_tmean_oct_apr_sum_rast / 7
  writeRaster(temp_tmean_oct_apr_mean_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmean_OctApr", year_num, "CDMMean.tif")), overwrite = TRUE)
  ## tmin sum
  temp_tmin_oct_apr_sum_rast <- temp_oct_tmin_rast + temp_nov_tmin_rast + temp_dec_tmin_rast + temp_jan_tmin_rast + temp_feb_tmin_rast + temp_mar_tmin_rast + temp_apr_tmin_rast
  ## tmin mean
  temp_tmin_oct_apr_mean_rast <- temp_tmin_oct_apr_sum_rast / 7
  writeRaster(temp_tmin_oct_apr_mean_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmin_OctApr", year_num, "Mean.tif")), overwrite = TRUE)
  ## tmax sum
  temp_tmax_oct_apr_sum_rast <- temp_oct_tmax_rast + temp_nov_tmax_rast + temp_dec_tmax_rast + temp_jan_tmax_rast + temp_feb_tmax_rast + temp_mar_tmax_rast + temp_apr_tmax_rast
  ## tmax mean
  temp_tmax_oct_apr_mean_rast <- temp_tmax_oct_apr_sum_rast / 7
  writeRaster(temp_tmax_oct_apr_mean_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmax_OctApr", year_num, "Mean.tif")), overwrite = TRUE)
  ## srad sum
  temp_srad_oct_apr_sum_rast <- temp_oct_srad_rast + temp_nov_srad_rast + temp_dec_srad_rast + temp_jan_srad_rast + temp_feb_srad_rast + temp_mar_srad_rast + temp_apr_srad_rast
  ## srad mean
  temp_srad_oct_apr_mean_rast <- temp_srad_oct_apr_sum_rast / 7
  writeRaster(temp_srad_oct_apr_mean_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("srad_OctApr", year_num, "Mean.tif")), overwrite = TRUE)

  ## Winter month aggregates (Dec-Feb)
  ## prcpMean sum
  temp_prcpMean_dec_feb_sum_rast <- temp_dec_prcpMean_rast + temp_jan_prcpMean_rast + temp_feb_prcpMean_rast
  writeRaster(temp_prcpMean_dec_feb_sum_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpMean_DecFeb", year_num, "CDMSum.tif")), overwrite = TRUE)
  ## prcpMean mean
  temp_prcpMean_dec_feb_mean_rast <- temp_prcpMean_dec_feb_sum_rast / 3
  writeRaster(temp_prcpMean_dec_feb_mean_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpMean_DecFeb", year_num, "CDMMean.tif")), overwrite = TRUE)
  ## prcpSum sum
  temp_prcpSum_dec_feb_sum_rast <- temp_dec_prcpSum_rast + temp_jan_prcpSum_rast + temp_feb_prcpSum_rast
  writeRaster(temp_prcpSum_dec_feb_sum_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpSum_DecFeb", year_num, "CDMSum.tif")), overwrite = TRUE)
  ## prcpSum mean
  temp_prcpSum_dec_feb_mean_rast <- temp_prcpSum_dec_feb_sum_rast / 3
  writeRaster(temp_prcpSum_dec_feb_mean_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpSum_DecFeb", year_num, "CDMMean.tif")), overwrite = TRUE)
  ## tmean sum
  temp_tmean_dec_feb_sum_rast <- temp_dec_tmean_rast + temp_jan_tmean_rast + temp_feb_tmean_rast
  writeRaster(temp_tmean_dec_feb_sum_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmean_DecFeb", year_num, "CDMSum.tif")), overwrite = TRUE)
  ## tmean mean
  temp_tmean_dec_feb_mean_rast <- temp_tmean_dec_feb_sum_rast / 3
  writeRaster(temp_tmean_dec_feb_mean_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmean_DecFeb", year_num, "CDMMean.tif")), overwrite = TRUE)
  ## tmin sum
  temp_tmin_dec_feb_sum_rast <- temp_dec_tmin_rast + temp_jan_tmin_rast + temp_feb_tmin_rast
  ## tmin mean
  temp_tmin_dec_feb_mean_rast <- temp_tmin_dec_feb_sum_rast / 3
  writeRaster(temp_tmin_dec_feb_mean_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmin_DecFeb", year_num, "Mean.tif")), overwrite = TRUE)
  ## tmax sum
  temp_tmax_dec_feb_sum_rast <- temp_dec_tmax_rast + temp_jan_tmax_rast + temp_feb_tmax_rast
  ## tmax mean
  temp_tmax_dec_feb_mean_rast <- temp_tmax_dec_feb_sum_rast / 3
  writeRaster(temp_tmax_dec_feb_mean_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmax_DecFeb", year_num, "Mean.tif")), overwrite = TRUE)
  ## srad sum
  temp_srad_dec_feb_sum_rast <- temp_dec_srad_rast + temp_jan_srad_rast + temp_feb_srad_rast
  ## srad mean
  temp_srad_dec_feb_mean_rast <- temp_srad_dec_feb_sum_rast / 3
  writeRaster(temp_srad_dec_feb_mean_rast, here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("srad_DecFeb", year_num, "Mean.tif")), overwrite = TRUE)
}
CONUSendtime <- Sys.time()


Alaskastarttime <- Sys.time()
##### ALASKA #####
##### reading in rasters for Alaska climatology CDM aggregates #####
### reading in rasters
## prcpMean
temp_sep_prcpMean_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "prcpMean_SepClimatologyCDM.tif"))
temp_oct_prcpMean_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "prcpMean_OctClimatologyCDM.tif"))
temp_nov_prcpMean_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "prcpMean_NovClimatologyCDM.tif"))
temp_dec_prcpMean_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "prcpMean_DecClimatologyCDM.tif"))
temp_jan_prcpMean_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "prcpMean_JanClimatologyCDM.tif"))
temp_feb_prcpMean_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "prcpMean_FebClimatologyCDM.tif"))
temp_mar_prcpMean_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "prcpMean_MarClimatologyCDM.tif"))
temp_apr_prcpMean_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "prcpMean_AprClimatologyCDM.tif"))
temp_may_prcpMean_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "prcpMean_MayClimatologyCDM.tif"))
## prcpSum
temp_sep_prcpSum_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "prcpSum_SepClimatologyCDM.tif"))
temp_oct_prcpSum_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "prcpSum_OctClimatologyCDM.tif"))
temp_nov_prcpSum_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "prcpSum_NovClimatologyCDM.tif"))
temp_dec_prcpSum_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "prcpSum_DecClimatologyCDM.tif"))
temp_jan_prcpSum_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "prcpSum_JanClimatologyCDM.tif"))
temp_feb_prcpSum_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "prcpSum_FebClimatologyCDM.tif"))
temp_mar_prcpSum_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "prcpSum_MarClimatologyCDM.tif"))
temp_apr_prcpSum_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "prcpSum_AprClimatologyCDM.tif"))
temp_may_prcpSum_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "prcpSum_MayClimatologyCDM.tif"))
## tmean
temp_sep_tmean_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "tmean_SepClimatologyCDM.tif"))
temp_oct_tmean_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "tmean_OctClimatologyCDM.tif"))
temp_nov_tmean_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "tmean_NovClimatologyCDM.tif"))
temp_dec_tmean_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "tmean_DecClimatologyCDM.tif"))
temp_jan_tmean_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "tmean_JanClimatologyCDM.tif"))
temp_feb_tmean_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "tmean_FebClimatologyCDM.tif"))
temp_mar_tmean_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "tmean_MarClimatologyCDM.tif"))
temp_apr_tmean_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "tmean_AprClimatologyCDM.tif"))
temp_may_tmean_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "tmean_MayClimatologyCDM.tif"))
## srad
temp_sep_srad_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "srad_Month9MeanClimatology.tif"))
temp_oct_srad_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "srad_Month10MeanClimatology.tif"))
temp_nov_srad_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "srad_Month11MeanClimatology.tif"))
temp_dec_srad_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "srad_Month12MeanClimatology.tif"))
temp_jan_srad_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "srad_Month1MeanClimatology.tif"))
temp_feb_srad_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "srad_Month2MeanClimatology.tif"))
temp_mar_srad_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "srad_Month3MeanClimatology.tif"))
temp_apr_srad_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "srad_Month4MeanClimatology.tif"))
temp_may_srad_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "srad_Month5MeanClimatology.tif"))
## tmin
temp_sep_tmin_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmin_Month9MeanClimatology.tif"))
temp_oct_tmin_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmin_Month10MeanClimatology.tif"))
temp_nov_tmin_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmin_Month11MeanClimatology.tif"))
temp_dec_tmin_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmin_Month12MeanClimatology.tif"))
temp_jan_tmin_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmin_Month1MeanClimatology.tif"))
temp_feb_tmin_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmin_Month2MeanClimatology.tif"))
temp_mar_tmin_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmin_Month3MeanClimatology.tif"))
temp_apr_tmin_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmin_Month4MeanClimatology.tif"))
temp_may_tmin_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmin_Month5MeanClimatology.tif"))
## tmax
temp_sep_tmax_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmax_Month9MeanClimatology.tif"))
temp_oct_tmax_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmax_Month10MeanClimatology.tif"))
temp_nov_tmax_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmax_Month11MeanClimatology.tif"))
temp_dec_tmax_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmax_Month12MeanClimatology.tif"))
temp_jan_tmax_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmax_Month1MeanClimatology.tif"))
temp_feb_tmax_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmax_Month2MeanClimatology.tif"))
temp_mar_tmax_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmax_Month3MeanClimatology.tif"))
temp_apr_tmax_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmax_Month4MeanClimatology.tif"))
temp_may_tmax_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmax_Month5MeanClimatology.tif"))

##### creating Alaska climatology CDM aggregates #####
print("Creating climatology CDM aggregates for Alaska")
## full "snow season" aggregate
## prcpMean sum
temp_prcpMean_sep_may_sum_rast <- temp_sep_prcpMean_rast + temp_oct_prcpMean_rast + temp_nov_prcpMean_rast + temp_dec_prcpMean_rast + temp_jan_prcpMean_rast + temp_feb_prcpMean_rast + temp_mar_prcpMean_rast + temp_apr_prcpMean_rast + temp_may_prcpMean_rast
writeRaster(temp_prcpMean_sep_may_sum_rast, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpMean_SepMaySum_ClimatologyCDM.tif"), overwrite = TRUE)
## prcpMean mean
temp_prcpMean_sep_may_mean_rast <- temp_prcpMean_sep_may_sum_rast / 9
writeRaster(temp_prcpMean_sep_may_mean_rast, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpMean_SepMayMean_ClimatologyCDM.tif"), overwrite = TRUE)
## prcpSum sum
temp_prcpSum_sep_may_sum_rast <- temp_sep_prcpSum_rast + temp_oct_prcpSum_rast + temp_nov_prcpSum_rast + temp_dec_prcpSum_rast + temp_jan_prcpSum_rast + temp_feb_prcpSum_rast + temp_mar_prcpSum_rast + temp_apr_prcpSum_rast + temp_may_prcpSum_rast
writeRaster(temp_prcpSum_sep_may_sum_rast, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpSum_SepMaySum_ClimatologyCDM.tif"), overwrite = TRUE)
## prcpSum mean
temp_prcpSum_sep_may_mean_rast <- temp_prcpSum_sep_may_sum_rast / 9
writeRaster(temp_prcpSum_sep_may_mean_rast, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpSum_SepMayMean_ClimatologyCDM.tif"), overwrite = TRUE)
## tmean sum
temp_tmean_sep_may_sum_rast <- temp_sep_tmean_rast + temp_oct_tmean_rast + temp_nov_tmean_rast + temp_dec_tmean_rast + temp_jan_tmean_rast + temp_feb_tmean_rast + temp_mar_tmean_rast + temp_apr_tmean_rast + temp_may_tmean_rast
writeRaster(temp_tmean_sep_may_sum_rast, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmean_SepMaySum_ClimatologyCDM.tif"), overwrite = TRUE)
## tmean mean
temp_tmean_sep_may_mean_rast <- temp_tmean_sep_may_sum_rast / 9
writeRaster(temp_tmean_sep_may_mean_rast, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmean_SepMayMean_ClimatologyCDM.tif"), overwrite = TRUE)
## tmin sum
temp_tmin_sep_may_sum_rast <- temp_sep_tmin_rast + temp_oct_tmin_rast + temp_nov_tmin_rast + temp_dec_tmin_rast + temp_jan_tmin_rast + temp_feb_tmin_rast + temp_mar_tmin_rast + temp_apr_tmin_rast + temp_may_tmin_rast
## tmin mean
temp_tmin_sep_may_mean_rast <- temp_tmin_sep_may_sum_rast / 9
writeRaster(temp_tmin_sep_may_mean_rast, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmin_SepMayMean_ClimatologyCDM.tif"), overwrite = TRUE)
## tmax sum
temp_tmax_sep_may_sum_rast <- temp_sep_tmax_rast + temp_oct_tmax_rast + temp_nov_tmax_rast + temp_dec_tmax_rast + temp_jan_tmax_rast + temp_feb_tmax_rast + temp_mar_tmax_rast + temp_apr_tmax_rast + temp_may_tmax_rast
## tmax mean
temp_tmax_sep_may_mean_rast <- temp_tmax_sep_may_sum_rast / 9
writeRaster(temp_tmax_sep_may_mean_rast, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmax_SepMayMean_ClimatologyCDM.tif"), overwrite = TRUE)
## srad sum
temp_srad_sep_may_sum_rast <- temp_sep_srad_rast + temp_oct_srad_rast + temp_nov_srad_rast + temp_dec_srad_rast + temp_jan_srad_rast + temp_feb_srad_rast + temp_mar_srad_rast + temp_apr_srad_rast + temp_may_srad_rast
## srad mean
temp_srad_sep_may_mean_rast <- temp_srad_sep_may_sum_rast / 9
writeRaster(temp_srad_sep_may_mean_rast, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "srad_SepMayMean_ClimatologyCDM.tif"), overwrite = TRUE)


## full snow season but just from start of water year (October)
## prcpMean sum
temp_prcpMean_oct_may_sum_rast <- temp_oct_prcpMean_rast + temp_nov_prcpMean_rast + temp_dec_prcpMean_rast + temp_jan_prcpMean_rast + temp_feb_prcpMean_rast + temp_mar_prcpMean_rast + temp_apr_prcpMean_rast + temp_may_prcpMean_rast
writeRaster(temp_prcpMean_oct_may_sum_rast, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpMean_OctMaySum_ClimatologyCDM.tif"), overwrite = TRUE)
## prcpMean mean
temp_prcpMean_oct_may_mean_rast <- temp_prcpMean_oct_may_sum_rast / 8
writeRaster(temp_prcpMean_oct_may_mean_rast, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpMean_OctMayMean_ClimatologyCDM.tif"), overwrite = TRUE)
## prcpSum sum
temp_prcpSum_oct_may_sum_rast <- temp_oct_prcpSum_rast + temp_nov_prcpSum_rast + temp_dec_prcpSum_rast + temp_jan_prcpSum_rast + temp_feb_prcpSum_rast + temp_mar_prcpSum_rast + temp_apr_prcpSum_rast + temp_may_prcpSum_rast
writeRaster(temp_prcpSum_oct_may_sum_rast, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpSum_OctMaySum_ClimatologyCDM.tif"), overwrite = TRUE)
## prcpSum mean
temp_prcpSum_oct_may_mean_rast <- temp_prcpSum_oct_may_sum_rast / 8
writeRaster(temp_prcpSum_oct_may_mean_rast, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpSum_OctMayMean_ClimatologyCDM.tif"), overwrite = TRUE)
## tmean sum
temp_tmean_oct_may_sum_rast <- temp_oct_tmean_rast + temp_nov_tmean_rast + temp_dec_tmean_rast + temp_jan_tmean_rast + temp_feb_tmean_rast + temp_mar_tmean_rast + temp_apr_tmean_rast + temp_may_tmean_rast
writeRaster(temp_tmean_oct_may_sum_rast, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmean_OctMaySum_ClimatologyCDM.tif"), overwrite = TRUE)
## tmean mean
temp_tmean_oct_may_mean_rast <- temp_tmean_oct_may_sum_rast / 8
writeRaster(temp_tmean_oct_may_mean_rast, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmean_OctMayMean_ClimatologyCDM.tif"), overwrite = TRUE)
## tmin sum
temp_tmin_oct_may_sum_rast <- temp_oct_tmin_rast + temp_nov_tmin_rast + temp_dec_tmin_rast + temp_jan_tmin_rast + temp_feb_tmin_rast + temp_mar_tmin_rast + temp_apr_tmin_rast + temp_may_tmin_rast
## tmin mean
temp_tmin_oct_may_mean_rast <- temp_tmin_oct_may_sum_rast / 8
writeRaster(temp_tmin_oct_may_mean_rast, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmin_OctMayMean_ClimatologyCDM.tif"), overwrite = TRUE)
## tmax sum
temp_tmax_oct_may_sum_rast <- temp_oct_tmax_rast + temp_nov_tmax_rast + temp_dec_tmax_rast + temp_jan_tmax_rast + temp_feb_tmax_rast + temp_mar_tmax_rast + temp_apr_tmax_rast + temp_may_tmax_rast
## tmax mean
temp_tmax_oct_may_mean_rast <- temp_tmax_oct_may_sum_rast / 8
writeRaster(temp_tmax_oct_may_mean_rast, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmax_OctMayMean_ClimatologyCDM.tif"), overwrite = TRUE)
## srad sum
temp_srad_oct_may_sum_rast <- temp_oct_srad_rast + temp_nov_srad_rast + temp_dec_srad_rast + temp_jan_srad_rast + temp_feb_srad_rast + temp_mar_srad_rast + temp_apr_srad_rast + temp_may_srad_rast
## srad mean
temp_srad_oct_may_mean_rast <- temp_srad_oct_may_sum_rast / 8
writeRaster(temp_srad_oct_may_mean_rast, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "srad_OctMayMean_ClimatologyCDM.tif"), overwrite = TRUE)



## "full" snow season but just from start of water year (October) to April
## prcpMean sum
temp_prcpMean_oct_apr_sum_rast <- temp_oct_prcpMean_rast + temp_nov_prcpMean_rast + temp_dec_prcpMean_rast + temp_jan_prcpMean_rast + temp_feb_prcpMean_rast + temp_mar_prcpMean_rast + temp_apr_prcpMean_rast
writeRaster(temp_prcpMean_oct_apr_sum_rast, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpMean_OctAprSum_ClimatologyCDM.tif"), overwrite = TRUE)
## prcpMean mean
temp_prcpMean_oct_apr_mean_rast <- temp_prcpMean_oct_apr_sum_rast / 7
writeRaster(temp_prcpMean_oct_apr_mean_rast, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpMean_OctAprMean_ClimatologyCDM.tif"), overwrite = TRUE)
## prcpSum sum
temp_prcpSum_oct_apr_sum_rast <- temp_oct_prcpSum_rast + temp_nov_prcpSum_rast + temp_dec_prcpSum_rast + temp_jan_prcpSum_rast + temp_feb_prcpSum_rast + temp_mar_prcpSum_rast + temp_apr_prcpSum_rast
writeRaster(temp_prcpSum_oct_apr_sum_rast, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpSum_OctAprSum_ClimatologyCDM.tif"), overwrite = TRUE)
## prcpSum mean
temp_prcpSum_oct_apr_mean_rast <- temp_prcpSum_oct_apr_sum_rast / 7
writeRaster(temp_prcpSum_oct_apr_mean_rast, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpSum_OctAprMean_ClimatologyCDM.tif"), overwrite = TRUE)
## tmean sum
temp_tmean_oct_apr_sum_rast <- temp_oct_tmean_rast + temp_nov_tmean_rast + temp_dec_tmean_rast + temp_jan_tmean_rast + temp_feb_tmean_rast + temp_mar_tmean_rast + temp_apr_tmean_rast
writeRaster(temp_tmean_oct_apr_sum_rast, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmean_OctAprSum_ClimatologyCDM.tif"), overwrite = TRUE)
## tmean mean
temp_tmean_oct_apr_mean_rast <- temp_tmean_oct_apr_sum_rast / 7
writeRaster(temp_tmean_oct_apr_mean_rast, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmean_OctAprMean_ClimatologyCDM.tif"), overwrite = TRUE)
## tmin sum
temp_tmin_oct_apr_sum_rast <- temp_oct_tmin_rast + temp_nov_tmin_rast + temp_dec_tmin_rast + temp_jan_tmin_rast + temp_feb_tmin_rast + temp_mar_tmin_rast + temp_apr_tmin_rast
## tmin mean
temp_tmin_oct_apr_mean_rast <- temp_tmin_oct_apr_sum_rast / 7
writeRaster(temp_tmin_oct_apr_mean_rast, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmin_OctAprMean_ClimatologyCDM.tif"), overwrite = TRUE)
## tmax sum
temp_tmax_oct_apr_sum_rast <- temp_oct_tmax_rast + temp_nov_tmax_rast + temp_dec_tmax_rast + temp_jan_tmax_rast + temp_feb_tmax_rast + temp_mar_tmax_rast + temp_apr_tmax_rast
## tmax mean
temp_tmax_oct_apr_mean_rast <- temp_tmax_oct_apr_sum_rast / 7
writeRaster(temp_tmax_oct_apr_mean_rast, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmax_OctAprMean_ClimatologyCDM.tif"), overwrite = TRUE)
## srad sum
temp_srad_oct_apr_sum_rast <- temp_oct_srad_rast + temp_nov_srad_rast + temp_dec_srad_rast + temp_jan_srad_rast + temp_feb_srad_rast + temp_mar_srad_rast + temp_apr_srad_rast
## srad mean
temp_srad_oct_apr_mean_rast <- temp_srad_oct_apr_sum_rast / 7
writeRaster(temp_srad_oct_apr_mean_rast, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "srad_OctAprMean_ClimatologyCDM.tif"), overwrite = TRUE)





## Winter month aggregates (Dec-Feb)
## prcpMean sum
temp_prcpMean_dec_feb_sum_rast <- temp_dec_prcpMean_rast + temp_jan_prcpMean_rast + temp_feb_prcpMean_rast
writeRaster(temp_prcpMean_dec_feb_sum_rast, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpMean_DecFebSum_ClimatologyCDM.tif"), overwrite = TRUE)
## prcpMean mean
temp_prcpMean_dec_feb_mean_rast <- temp_prcpMean_dec_feb_sum_rast / 3
writeRaster(temp_prcpMean_dec_feb_mean_rast, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpMean_DecFebMean_ClimatologyCDM.tif"), overwrite = TRUE)
## prcpSum sum
temp_prcpSum_dec_feb_sum_rast <- temp_dec_prcpSum_rast + temp_jan_prcpSum_rast + temp_feb_prcpSum_rast
writeRaster(temp_prcpSum_dec_feb_sum_rast, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpSum_DecFebSum_ClimatologyCDM.tif"), overwrite = TRUE)
## prcpSum mean
temp_prcpSum_dec_feb_mean_rast <- temp_prcpSum_dec_feb_sum_rast / 3
writeRaster(temp_prcpSum_dec_feb_mean_rast, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpSum_DecFebMean_ClimatologyCDM.tif"), overwrite = TRUE)
## tmean sum
temp_tmean_dec_feb_sum_rast <- temp_dec_tmean_rast + temp_jan_tmean_rast + temp_feb_tmean_rast
writeRaster(temp_tmean_dec_feb_sum_rast, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmean_DecFebSum_ClimatologyCDM.tif"), overwrite = TRUE)
## tmean mean
temp_tmean_dec_feb_mean_rast <- temp_tmean_dec_feb_sum_rast / 3
writeRaster(temp_tmean_dec_feb_mean_rast, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmean_DecFebMean_ClimatologyCDM.tif"), overwrite = TRUE)
## tmin sum
temp_tmin_dec_feb_sum_rast <- temp_dec_tmin_rast + temp_jan_tmin_rast + temp_feb_tmin_rast
## tmin mean
temp_tmin_dec_feb_mean_rast <- temp_tmin_dec_feb_sum_rast / 3
writeRaster(temp_tmin_dec_feb_mean_rast, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmin_DecFebMean_ClimatologyCDM.tif"), overwrite = TRUE)
## tmax sum
temp_tmax_dec_feb_sum_rast <- temp_dec_tmax_rast + temp_jan_tmax_rast + temp_feb_tmax_rast
## tmax mean
temp_tmax_dec_feb_mean_rast <- temp_tmax_dec_feb_sum_rast / 3
writeRaster(temp_tmax_dec_feb_mean_rast, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmax_DecFebMean_ClimatologyCDM.tif"), overwrite = TRUE)
## srad sum
temp_srad_dec_feb_sum_rast <- temp_dec_srad_rast + temp_jan_srad_rast + temp_feb_srad_rast
## srad mean
temp_srad_dec_feb_mean_rast <- temp_srad_dec_feb_sum_rast / 3
writeRaster(temp_srad_dec_feb_mean_rast, here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "srad_DecFebMean_ClimatologyCDM.tif"), overwrite = TRUE)

##### annual Alaska CDM aggregates #####
for (year_num in 1993:2020){
  print("Creating annual CDM aggregates for Alaska")
  print(year_num)
  ### reading in prcpMean rasters for given year
  temp_sep_prcpMean_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("prcpMean_Sep", year_num - 1, "CDM.tif")))
  temp_oct_prcpMean_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("prcpMean_Oct", year_num - 1, "CDM.tif")))
  temp_nov_prcpMean_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("prcpMean_Nov", year_num - 1, "CDM.tif")))
  temp_dec_prcpMean_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("prcpMean_Dec", year_num - 1, "CDM.tif")))
  temp_jan_prcpMean_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("prcpMean_Jan", year_num, "CDM.tif")))
  temp_feb_prcpMean_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("prcpMean_Feb", year_num, "CDM.tif")))
  temp_mar_prcpMean_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("prcpMean_Mar", year_num, "CDM.tif")))
  temp_apr_prcpMean_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("prcpMean_Apr", year_num, "CDM.tif")))
  temp_may_prcpMean_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("prcpMean_May", year_num, "CDM.tif")))
  ### reading in prcpSum rasters for given year
  temp_sep_prcpSum_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("prcpSum_Sep", year_num - 1, "CDM.tif")))
  temp_oct_prcpSum_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("prcpSum_Oct", year_num - 1, "CDM.tif")))
  temp_nov_prcpSum_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("prcpSum_Nov", year_num - 1, "CDM.tif")))
  temp_dec_prcpSum_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("prcpSum_Dec", year_num - 1, "CDM.tif")))
  temp_jan_prcpSum_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("prcpSum_Jan", year_num, "CDM.tif")))
  temp_feb_prcpSum_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("prcpSum_Feb", year_num, "CDM.tif")))
  temp_mar_prcpSum_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("prcpSum_Mar", year_num, "CDM.tif")))
  temp_apr_prcpSum_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("prcpSum_Apr", year_num, "CDM.tif")))
  temp_may_prcpSum_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("prcpSum_May", year_num, "CDM.tif")))
  ### reading in tmean rasters for given year
  temp_sep_tmean_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("tmean_Sep", year_num - 1, "CDM.tif")))
  temp_oct_tmean_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("tmean_Oct", year_num - 1, "CDM.tif")))
  temp_nov_tmean_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("tmean_Nov", year_num - 1, "CDM.tif")))
  temp_dec_tmean_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("tmean_Dec", year_num - 1, "CDM.tif")))
  temp_jan_tmean_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("tmean_Jan", year_num, "CDM.tif")))
  temp_feb_tmean_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("tmean_Feb", year_num, "CDM.tif")))
  temp_mar_tmean_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("tmean_Mar", year_num, "CDM.tif")))
  temp_apr_tmean_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("tmean_Apr", year_num, "CDM.tif")))
  temp_may_tmean_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", "CDM", paste0("tmean_May", year_num, "CDM.tif")))
  
  ### reading in tmin rasters for given year
  temp_sep_tmin_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", paste0("tmin", year_num - 1, "_SepMean.tif")))
  temp_oct_tmin_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", paste0("tmin", year_num - 1, "_OctMean.tif")))
  temp_nov_tmin_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", paste0("tmin", year_num - 1, "_NovMean.tif")))
  temp_dec_tmin_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", paste0("tmin", year_num - 1, "_DecMean.tif")))
  temp_jan_tmin_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", paste0("tmin", year_num, "JanMean.tif")))
  temp_feb_tmin_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", paste0("tmin", year_num, "_FebMean.tif")))
  temp_mar_tmin_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", paste0("tmin", year_num, "_MarMean.tif")))
  temp_apr_tmin_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", paste0("tmin", year_num, "_AprMean.tif")))
  temp_may_tmin_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", paste0("tmin", year_num, "_MayMean.tif")))
  
  ### reading in tmax rasters for given year
  temp_sep_tmax_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", paste0("tmax", year_num - 1, "_SepMean.tif")))
  temp_oct_tmax_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", paste0("tmax", year_num - 1, "_OctMean.tif")))
  temp_nov_tmax_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", paste0("tmax", year_num - 1, "_NovMean.tif")))
  temp_dec_tmax_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", paste0("tmax", year_num - 1, "_DecMean.tif")))
  temp_jan_tmax_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", paste0("tmax", year_num, "JanMean.tif")))
  temp_feb_tmax_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", paste0("tmax", year_num, "_FebMean.tif")))
  temp_mar_tmax_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", paste0("tmax", year_num, "_MarMean.tif")))
  temp_apr_tmax_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", paste0("tmax", year_num, "_AprMean.tif")))
  temp_may_tmax_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", paste0("tmax", year_num, "_MayMean.tif")))
  ### reading in srad rasters for given year
  temp_sep_srad_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", paste0("srad", year_num - 1, "_SepMean.tif")))
  temp_oct_srad_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", paste0("srad", year_num - 1, "_OctMean.tif")))
  temp_nov_srad_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", paste0("srad", year_num - 1, "_NovMean.tif")))
  temp_dec_srad_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", paste0("srad", year_num - 1, "_DecMean.tif")))
  temp_jan_srad_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", paste0("srad", year_num, "JanMean.tif")))
  temp_feb_srad_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", paste0("srad", year_num, "_FebMean.tif")))
  temp_mar_srad_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", paste0("srad", year_num, "_MarMean.tif")))
  temp_apr_srad_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", paste0("srad", year_num, "_AprMean.tif")))
  temp_may_srad_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics",  "TIFFs", "ClippedResampled", paste0("srad", year_num, "_MayMean.tif")))
  
  
  ### creating CDM aggregates using WY-specific rasters
  ##### creating annual Alaska CDM aggregates #####
  ## full "snow season" aggregate
  ## prcpMean sum
  temp_prcpMean_sep_may_sum_rast <- temp_sep_prcpMean_rast + temp_oct_prcpMean_rast + temp_nov_prcpMean_rast + temp_dec_prcpMean_rast + temp_jan_prcpMean_rast + temp_feb_prcpMean_rast + temp_mar_prcpMean_rast + temp_apr_prcpMean_rast + temp_may_prcpMean_rast
  writeRaster(temp_prcpMean_sep_may_sum_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpMean_SepMay", year_num, "CDMSum.tif")), overwrite = TRUE)
  ## prcpMean mean
  temp_prcpMean_sep_may_mean_rast <- temp_prcpMean_sep_may_sum_rast / 9
  writeRaster(temp_prcpMean_sep_may_mean_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpMean_SepMay", year_num, "CDMMean.tif")), overwrite = TRUE)
  ## prcpSum sum
  temp_prcpSum_sep_may_sum_rast <- temp_sep_prcpSum_rast + temp_oct_prcpSum_rast + temp_nov_prcpSum_rast + temp_dec_prcpSum_rast + temp_jan_prcpSum_rast + temp_feb_prcpSum_rast + temp_mar_prcpSum_rast + temp_apr_prcpSum_rast + temp_may_prcpSum_rast
  writeRaster(temp_prcpSum_sep_may_sum_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpSum_SepMay", year_num, "CDMSum.tif")), overwrite = TRUE)
  ## prcpSum mean
  temp_prcpSum_sep_may_mean_rast <- temp_prcpSum_sep_may_sum_rast / 9
  writeRaster(temp_prcpSum_sep_may_mean_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpSum_SepMay", year_num, "CDMMean.tif")), overwrite = TRUE)
  ## tmean sum
  temp_tmean_sep_may_sum_rast <- temp_sep_tmean_rast + temp_oct_tmean_rast + temp_nov_tmean_rast + temp_dec_tmean_rast + temp_jan_tmean_rast + temp_feb_tmean_rast + temp_mar_tmean_rast + temp_apr_tmean_rast + temp_may_tmean_rast
  writeRaster(temp_tmean_sep_may_sum_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmean_SepMay", year_num, "CDMSum.tif")), overwrite = TRUE)
  ## tmean mean
  temp_tmean_sep_may_mean_rast <- temp_tmean_sep_may_sum_rast / 9
  writeRaster(temp_tmean_sep_may_mean_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmean_SepMay", year_num, "CDMMean.tif")), overwrite = TRUE)
  ## tmin sum
  temp_tmin_sep_may_sum_rast <- temp_sep_tmin_rast + temp_oct_tmin_rast + temp_nov_tmin_rast + temp_dec_tmin_rast + temp_jan_tmin_rast + temp_feb_tmin_rast + temp_mar_tmin_rast + temp_apr_tmin_rast + temp_may_tmin_rast
  ## tmin mean
  temp_tmin_sep_may_mean_rast <- temp_tmin_sep_may_sum_rast / 9
  writeRaster(temp_tmin_sep_may_mean_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmin_SepMay", year_num, "Mean.tif")), overwrite = TRUE)
  ## tmax sum
  temp_tmax_sep_may_sum_rast <- temp_sep_tmax_rast + temp_oct_tmax_rast + temp_nov_tmax_rast + temp_dec_tmax_rast + temp_jan_tmax_rast + temp_feb_tmax_rast + temp_mar_tmax_rast + temp_apr_tmax_rast + temp_may_tmax_rast
  ## tmax mean
  temp_tmax_sep_may_mean_rast <- temp_tmax_sep_may_sum_rast / 9
  writeRaster(temp_tmax_sep_may_mean_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmax_SepMay", year_num, "Mean.tif")), overwrite = TRUE)
  ## srad sum
  temp_srad_sep_may_sum_rast <- temp_sep_srad_rast + temp_oct_srad_rast + temp_nov_srad_rast + temp_dec_srad_rast + temp_jan_srad_rast + temp_feb_srad_rast + temp_mar_srad_rast + temp_apr_srad_rast + temp_may_srad_rast
  ## srad mean
  temp_srad_sep_may_mean_rast <- temp_srad_sep_may_sum_rast / 9
  writeRaster(temp_srad_sep_may_mean_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("srad_SepMay", year_num, "Mean.tif")), overwrite = TRUE)
  
  
  ## full snow season but just from start of water year (October)
  ## prcpMean sum
  temp_prcpMean_oct_may_sum_rast <- temp_oct_prcpMean_rast + temp_nov_prcpMean_rast + temp_dec_prcpMean_rast + temp_jan_prcpMean_rast + temp_feb_prcpMean_rast + temp_mar_prcpMean_rast + temp_apr_prcpMean_rast + temp_may_prcpMean_rast
  writeRaster(temp_prcpMean_oct_may_sum_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpMean_OctMay", year_num, "CDMSum.tif")), overwrite = TRUE)
  ## prcpMean mean
  temp_prcpMean_oct_may_mean_rast <- temp_prcpMean_oct_may_sum_rast / 8
  writeRaster(temp_prcpMean_oct_may_mean_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpMean_OctMay", year_num, "CDMMean.tif")), overwrite = TRUE)
  ## prcpSum sum
  temp_prcpSum_oct_may_sum_rast <- temp_oct_prcpSum_rast + temp_nov_prcpSum_rast + temp_dec_prcpSum_rast + temp_jan_prcpSum_rast + temp_feb_prcpSum_rast + temp_mar_prcpSum_rast + temp_apr_prcpSum_rast + temp_may_prcpSum_rast
  writeRaster(temp_prcpSum_oct_may_sum_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpSum_OctMay", year_num, "CDMSum.tif")), overwrite = TRUE)
  ## prcpSum mean
  temp_prcpSum_oct_may_mean_rast <- temp_prcpSum_oct_may_sum_rast / 8
  writeRaster(temp_prcpSum_oct_may_mean_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpSum_OctMay", year_num, "CDMMean.tif")), overwrite = TRUE)
  ## tmean sum
  temp_tmean_oct_may_sum_rast <- temp_oct_tmean_rast + temp_nov_tmean_rast + temp_dec_tmean_rast + temp_jan_tmean_rast + temp_feb_tmean_rast + temp_mar_tmean_rast + temp_apr_tmean_rast + temp_may_tmean_rast
  writeRaster(temp_tmean_oct_may_sum_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmean_OctMay", year_num, "CDMSum.tif")), overwrite = TRUE)
  ## tmean mean
  temp_tmean_oct_may_mean_rast <- temp_tmean_oct_may_sum_rast / 8
  writeRaster(temp_tmean_oct_may_mean_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmean_OctMay", year_num, "CDMMean.tif")), overwrite = TRUE)
  ## tmin sum
  temp_tmin_oct_may_sum_rast <- temp_oct_tmin_rast + temp_nov_tmin_rast + temp_dec_tmin_rast + temp_jan_tmin_rast + temp_feb_tmin_rast + temp_mar_tmin_rast + temp_apr_tmin_rast + temp_may_tmin_rast
  ## tmin mean
  temp_tmin_oct_may_mean_rast <- temp_tmin_oct_may_sum_rast / 8
  writeRaster(temp_tmin_oct_may_mean_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmin_OctMay", year_num, "Mean.tif")), overwrite = TRUE)
  ## tmax sum
  temp_tmax_oct_may_sum_rast <- temp_oct_tmax_rast + temp_nov_tmax_rast + temp_dec_tmax_rast + temp_jan_tmax_rast + temp_feb_tmax_rast + temp_mar_tmax_rast + temp_apr_tmax_rast + temp_may_tmax_rast
  ## tmax mean
  temp_tmax_oct_may_mean_rast <- temp_tmax_oct_may_sum_rast / 8
  writeRaster(temp_tmax_oct_may_mean_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmax_OctMay", year_num, "Mean.tif")), overwrite = TRUE)
  ## srad sum
  temp_srad_oct_may_sum_rast <- temp_oct_srad_rast + temp_nov_srad_rast + temp_dec_srad_rast + temp_jan_srad_rast + temp_feb_srad_rast + temp_mar_srad_rast + temp_apr_srad_rast + temp_may_srad_rast
  ## srad mean
  temp_srad_oct_may_mean_rast <- temp_srad_oct_may_sum_rast / 8
  writeRaster(temp_srad_oct_may_mean_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("srad_OctMay", year_num, "Mean.tif")), overwrite = TRUE)
  
  
  
  ## "full" snow season but just from start of water year (October) to April
  ## prcpMean sum
  temp_prcpMean_oct_apr_sum_rast <- temp_oct_prcpMean_rast + temp_nov_prcpMean_rast + temp_dec_prcpMean_rast + temp_jan_prcpMean_rast + temp_feb_prcpMean_rast + temp_mar_prcpMean_rast + temp_apr_prcpMean_rast
  writeRaster(temp_prcpMean_oct_apr_sum_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpMean_OctApr", year_num, "CDMSum.tif")), overwrite = TRUE)
  ## prcpMean mean
  temp_prcpMean_oct_apr_mean_rast <- temp_prcpMean_oct_apr_sum_rast / 7
  writeRaster(temp_prcpMean_oct_apr_mean_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpMean_OctApr", year_num, "CDMMean.tif")), overwrite = TRUE)
  ## prcpSum sum
  temp_prcpSum_oct_apr_sum_rast <- temp_oct_prcpSum_rast + temp_nov_prcpSum_rast + temp_dec_prcpSum_rast + temp_jan_prcpSum_rast + temp_feb_prcpSum_rast + temp_mar_prcpSum_rast + temp_apr_prcpSum_rast
  writeRaster(temp_prcpSum_oct_apr_sum_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpSum_OctApr", year_num, "CDMSum.tif")), overwrite = TRUE)
  ## prcpSum mean
  temp_prcpSum_oct_apr_mean_rast <- temp_prcpSum_oct_apr_sum_rast / 7
  writeRaster(temp_prcpSum_oct_apr_mean_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpSum_OctApr", year_num, "CDMMean.tif")), overwrite = TRUE)
  ## tmean sum
  temp_tmean_oct_apr_sum_rast <- temp_oct_tmean_rast + temp_nov_tmean_rast + temp_dec_tmean_rast + temp_jan_tmean_rast + temp_feb_tmean_rast + temp_mar_tmean_rast + temp_apr_tmean_rast
  writeRaster(temp_tmean_oct_apr_sum_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmean_OctApr", year_num, "CDMSum.tif")), overwrite = TRUE)
  ## tmean mean
  temp_tmean_oct_apr_mean_rast <- temp_tmean_oct_apr_sum_rast / 7
  writeRaster(temp_tmean_oct_apr_mean_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmean_OctApr", year_num, "CDMMean.tif")), overwrite = TRUE)
  ## tmin sum
  temp_tmin_oct_apr_sum_rast <- temp_oct_tmin_rast + temp_nov_tmin_rast + temp_dec_tmin_rast + temp_jan_tmin_rast + temp_feb_tmin_rast + temp_mar_tmin_rast + temp_apr_tmin_rast
  ## tmin mean
  temp_tmin_oct_apr_mean_rast <- temp_tmin_oct_apr_sum_rast / 7
  writeRaster(temp_tmin_oct_apr_mean_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmin_OctApr", year_num, "Mean.tif")), overwrite = TRUE)
  ## tmax sum
  temp_tmax_oct_apr_sum_rast <- temp_oct_tmax_rast + temp_nov_tmax_rast + temp_dec_tmax_rast + temp_jan_tmax_rast + temp_feb_tmax_rast + temp_mar_tmax_rast + temp_apr_tmax_rast
  ## tmax mean
  temp_tmax_oct_apr_mean_rast <- temp_tmax_oct_apr_sum_rast / 7
  writeRaster(temp_tmax_oct_apr_mean_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmax_OctApr", year_num, "Mean.tif")), overwrite = TRUE)
  ## srad sum
  temp_srad_oct_apr_sum_rast <- temp_oct_srad_rast + temp_nov_srad_rast + temp_dec_srad_rast + temp_jan_srad_rast + temp_feb_srad_rast + temp_mar_srad_rast + temp_apr_srad_rast
  ## srad mean
  temp_srad_oct_apr_mean_rast <- temp_srad_oct_apr_sum_rast / 7
  writeRaster(temp_srad_oct_apr_mean_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("srad_OctApr", year_num, "Mean.tif")), overwrite = TRUE)
  
  ## Winter month aggregates (Dec-Feb)
  ## prcpMean sum
  temp_prcpMean_dec_feb_sum_rast <- temp_dec_prcpMean_rast + temp_jan_prcpMean_rast + temp_feb_prcpMean_rast
  writeRaster(temp_prcpMean_dec_feb_sum_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpMean_DecFeb", year_num, "CDMSum.tif")), overwrite = TRUE)
  ## prcpMean mean
  temp_prcpMean_dec_feb_mean_rast <- temp_prcpMean_dec_feb_sum_rast / 3
  writeRaster(temp_prcpMean_dec_feb_mean_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpMean_DecFeb", year_num, "CDMMean.tif")), overwrite = TRUE)
  ## prcpSum sum
  temp_prcpSum_dec_feb_sum_rast <- temp_dec_prcpSum_rast + temp_jan_prcpSum_rast + temp_feb_prcpSum_rast
  writeRaster(temp_prcpSum_dec_feb_sum_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpSum_DecFeb", year_num, "CDMSum.tif")), overwrite = TRUE)
  ## prcpSum mean
  temp_prcpSum_dec_feb_mean_rast <- temp_prcpSum_dec_feb_sum_rast / 3
  writeRaster(temp_prcpSum_dec_feb_mean_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpSum_DecFeb", year_num, "CDMMean.tif")), overwrite = TRUE)
  ## tmean sum
  temp_tmean_dec_feb_sum_rast <- temp_dec_tmean_rast + temp_jan_tmean_rast + temp_feb_tmean_rast
  writeRaster(temp_tmean_dec_feb_sum_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmean_DecFeb", year_num, "CDMSum.tif")), overwrite = TRUE)
  ## tmean mean
  temp_tmean_dec_feb_mean_rast <- temp_tmean_dec_feb_sum_rast / 3
  writeRaster(temp_tmean_dec_feb_mean_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmean_DecFeb", year_num, "CDMMean.tif")), overwrite = TRUE)
  ## tmin sum
  temp_tmin_dec_feb_sum_rast <- temp_dec_tmin_rast + temp_jan_tmin_rast + temp_feb_tmin_rast
  ## tmin mean
  temp_tmin_dec_feb_mean_rast <- temp_tmin_dec_feb_sum_rast / 3
  writeRaster(temp_tmin_dec_feb_mean_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmin_DecFeb", year_num, "Mean.tif")), overwrite = TRUE)
  ## tmax sum
  temp_tmax_dec_feb_sum_rast <- temp_dec_tmax_rast + temp_jan_tmax_rast + temp_feb_tmax_rast
  ## tmax mean
  temp_tmax_dec_feb_mean_rast <- temp_tmax_dec_feb_sum_rast / 3
  writeRaster(temp_tmax_dec_feb_mean_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmax_DecFeb", year_num, "Mean.tif")), overwrite = TRUE)
  ## srad sum
  temp_srad_dec_feb_sum_rast <- temp_dec_srad_rast + temp_jan_srad_rast + temp_feb_srad_rast
  ## srad mean
  temp_srad_dec_feb_mean_rast <- temp_srad_dec_feb_sum_rast / 3
  writeRaster(temp_srad_dec_feb_mean_rast, here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("srad_DecFeb", year_num, "Mean.tif")), overwrite = TRUE)
}

Alaskaendtime <- Sys.time()
CONUStime <- CONUSendtime - CONUSstarttime
Alaskatime <- Alaskaendtime - Alaskastarttime
Endtime <- Sys.time()
Totaltime <- Endtime - starttime
CONUStime
Alaskatime
Totaltime
##### End of Script #####
