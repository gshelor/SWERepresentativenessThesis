##### Preparation for random forest modeling of Peak SWE in CONUS #####
Totalstarttime <- Sys.time()
##### loading packages, reading in data #####
library(pacman)
p_load(here, tidyverse, sf, terra, parallel, rsample, mcprogress)
# options(mc.cores = parallel::detectCores())

CONUSClimstarttime <- Sys.time()
### all peak swe points for each site for all 28 years
# Snotel_CONUS <- read_sf(here("Data", "SNOTEL", "CONUS", "GIS", "GPKG", "Climatology", "SnotelData_CONUS_PeakSWEClimatology.gpkg"))
# Snotel_AK <- read_sf(here("Data", "SNOTEL", "Alaska", "GIS", "GPKG", "Climatology", "SnotelData_AK_PeakSWEClimatology.gpkg"))
CONUS_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "CONUS", "CONUS_AOI.gpkg"))
AK_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "Alaska", "AK_AOI.gpkg"))


##### reading in CONUS climatologies #####
### prcpSum
# prcpSum_JanClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "prcpSum_Month1MeanClimatology.tif"))
# prcpSum_FebClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "prcpSum_Month2MeanClimatology.tif"))
# prcpSum_MarClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "prcpSum_Month3MeanClimatology.tif"))
# prcpSum_AprClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "prcpSum_Month4MeanClimatology.tif"))
# prcpSum_MayClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "prcpSum_Month5MeanClimatology.tif"))
# prcpSum_JuneClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "prcpSum_Month6MeanClimatology.tif"))
# prcpSum_JulyClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "prcpSum_Month7MeanClimatology.tif"))
# prcpSum_AugClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "prcpSum_Month8MeanClimatology.tif"))
# prcpSum_SepClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "prcpSum_Month9MeanClimatology.tif"))
# prcpSum_OctClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "prcpSum_Month10MeanClimatology.tif"))
# prcpSum_NovClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "prcpSum_Month11MeanClimatology.tif"))
# prcpSum_DecClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "prcpSum_Month12MeanClimatology.tif"))
# ### prcpMean
# prcpMean_JanClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "prcp_Month1MeanClimatology.tif"))
# prcpMean_FebClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "prcp_Month2MeanClimatology.tif"))
# prcpMean_MarClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "prcp_Month3MeanClimatology.tif"))
# prcpMean_AprClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "prcp_Month4MeanClimatology.tif"))
# prcpMean_MayClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "prcp_Month5MeanClimatology.tif"))
# prcpMean_JuneClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "prcp_Month6MeanClimatology.tif"))
# prcpMean_JulyClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "prcp_Month7MeanClimatology.tif"))
# prcpMean_AugClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "prcp_Month8MeanClimatology.tif"))
# prcpMean_SepClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "prcp_Month9MeanClimatology.tif"))
# prcpMean_OctClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "prcp_Month10MeanClimatology.tif"))
# prcpMean_NovClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "prcp_Month11MeanClimatology.tif"))
# prcpMean_DecClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "prcp_Month12MeanClimatology.tif"))
# ## tmin
# tmin_JanClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmin_Month1MeanClimatology.tif"))
# tmin_FebClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmin_Month2MeanClimatology.tif"))
# tmin_MarClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmin_Month3MeanClimatology.tif"))
# tmin_AprClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmin_Month4MeanClimatology.tif"))
# tmin_MayClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmin_Month5MeanClimatology.tif"))
# tmin_JuneClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmin_Month6MeanClimatology.tif"))
# tmin_JulyClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmin_Month7MeanClimatology.tif"))
# tmin_AugClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmin_Month8MeanClimatology.tif"))
# tmin_SepClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmin_Month9MeanClimatology.tif"))
# tmin_OctClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmin_Month10MeanClimatology.tif"))
# tmin_NovClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmin_Month11MeanClimatology.tif"))
# tmin_DecClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmin_Month12MeanClimatology.tif"))
# 
# ## tmean
# tmean_JanClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmean_Month1MeanClimatology.tif"))
# tmean_FebClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmean_Month2MeanClimatology.tif"))
# tmean_MarClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmean_Month3MeanClimatology.tif"))
# tmean_AprClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmean_Month4MeanClimatology.tif"))
# tmean_MayClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmean_Month5MeanClimatology.tif"))
# tmean_JuneClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmean_Month6MeanClimatology.tif"))
# tmean_JulyClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmean_Month7MeanClimatology.tif"))
# tmean_AugClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmean_Month8MeanClimatology.tif"))
# tmean_SepClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmean_Month9MeanClimatology.tif"))
# tmean_OctClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmean_Month10MeanClimatology.tif"))
# tmean_NovClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmean_Month11MeanClimatology.tif"))
# tmean_DecClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmean_Month12MeanClimatology.tif"))
# 
# ## tmax
# tmax_JanClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmax_Month1MeanClimatology.tif"))
# tmax_FebClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmax_Month2MeanClimatology.tif"))
# tmax_MarClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmax_Month3MeanClimatology.tif"))
# tmax_AprClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmax_Month4MeanClimatology.tif"))
# tmax_MayClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmax_Month5MeanClimatology.tif"))
# tmax_JuneClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmax_Month6MeanClimatology.tif"))
# tmax_JulyClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmax_Month7MeanClimatology.tif"))
# tmax_AugClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmax_Month8MeanClimatology.tif"))
# tmax_SepClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmax_Month9MeanClimatology.tif"))
# tmax_OctClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmax_Month10MeanClimatology.tif"))
# tmax_NovClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmax_Month11MeanClimatology.tif"))
# tmax_DecClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "tmax_Month12MeanClimatology.tif"))
# 
# ## srad
# srad_JanClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "srad_Month1MeanClimatology.tif"))
# srad_FebClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "srad_Month2MeanClimatology.tif"))
# srad_MarClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "srad_Month3MeanClimatology.tif"))
# srad_AprClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "srad_Month4MeanClimatology.tif"))
# srad_MayClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "srad_Month5MeanClimatology.tif"))
# srad_JuneClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "srad_Month6MeanClimatology.tif"))
# srad_JulyClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "srad_Month7MeanClimatology.tif"))
# srad_AugClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "srad_Month8MeanClimatology.tif"))
# srad_SepClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "srad_Month9MeanClimatology.tif"))
# srad_OctClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "srad_Month10MeanClimatology.tif"))
# srad_NovClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "srad_Month11MeanClimatology.tif"))
# srad_DecClim <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "srad_Month12MeanClimatology.tif"))
# 
# 
# ##### reading in CONUS climatology CDMs #####
# ### Oct-Apr Sums
# prcpSum_OctAprClimSum_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpSum_OctAprSum_ClimatologyCDM.tif"))
# prcpMean_OctAprClimMean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpMean_OctAprMean_ClimatologyCDM.tif"))
# tmean_OctAprClimSum_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmean_OctAprSum_ClimatologyCDM.tif"))
# tmin_OctAprClimMean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmin_OctAprMean_ClimatologyCDM.tif"))
# tmax_OctAprClimMean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmax_OctAprMean_ClimatologyCDM.tif"))
# srad_OctAprClimMean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "srad_OctAprMean_ClimatologyCDM.tif"))
# ### Oct-May Sums
# prcpSum_OctMayClimSum_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpSum_OctMaySum_ClimatologyCDM.tif"))
# prcpMean_OctMayClimMean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpMean_OctMayMean_ClimatologyCDM.tif"))
# tmean_OctMayClimSum_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmean_OctMaySum_ClimatologyCDM.tif"))
# tmin_OctMayClimMean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmin_OctMayMean_ClimatologyCDM.tif"))
# tmax_OctMayClimMean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmax_OctMayMean_ClimatologyCDM.tif"))
# srad_OctMayClimMean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "srad_OctMayMean_ClimatologyCDM.tif"))
# ### Sep-May Sums
# prcpSum_SepMayClimSum_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpSum_SepMaySum_ClimatologyCDM.tif"))
# prcpMean_SepMayClimMean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpMean_SepMayMean_ClimatologyCDM.tif"))
# tmean_SepMayClimSum_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmean_SepMaySum_ClimatologyCDM.tif"))
# tmin_SepMayClimMean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmin_SepMayMean_ClimatologyCDM.tif"))
# tmax_SepMayClimMean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmax_SepMayMean_ClimatologyCDM.tif"))
# srad_SepMayClimMean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "srad_SepMayMean_ClimatologyCDM.tif"))
# ### Dec-Feb Sums
# prcpSum_DecFebClimSum_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpSum_DecFebSum_ClimatologyCDM.tif"))
# prcpMean_DecFebClimMean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpMean_DecFebMean_ClimatologyCDM.tif"))
# tmean_DecFebClimSum_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmean_DecFebSum_ClimatologyCDM.tif"))
# tmin_DecFebClimMean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmin_DecFebMean_ClimatologyCDM.tif"))
# tmax_DecFebClimMean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmax_DecFebMean_ClimatologyCDM.tif"))
# srad_DecFebClimMean_rast <- rast(here("Data", "Daymet", "CONUS", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "srad_DecFebMean_ClimatologyCDM.tif"))


##### reading in CONUS DEM, slope, aspect #####
### CONUS
CONUS_DEM <- rast(here("Data", "DEM", "CONUS", "CONUSDEMMosaic.tif"))
CONUS_Slope <- rast(here("Data", "DEM", "CONUS", "CONUSSlope.tif"))
CONUS_Aspect <- rast(here("Data", "DEM", "CONUS", "CONUSAspect.tif"))


##### reading in CONUS Landcover #####
### CONUS
LC_CONUS_1992 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "LC_CONUS_1992.tif"))
LC_CONUS_1993 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "LC_CONUS_1993.tif"))
LC_CONUS_1994 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "LC_CONUS_1994.tif"))
LC_CONUS_1995 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "LC_CONUS_1995.tif"))
LC_CONUS_1996 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "LC_CONUS_1996.tif"))
LC_CONUS_1997 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "LC_CONUS_1997.tif"))
LC_CONUS_1998 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "LC_CONUS_1998.tif"))
LC_CONUS_1999 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "LC_CONUS_1999.tif"))
LC_CONUS_2000 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "LC_CONUS_2000.tif"))
LC_CONUS_2001 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "LC_CONUS_2001.tif"))
LC_CONUS_2002 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "LC_CONUS_2002.tif"))
LC_CONUS_2003 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "LC_CONUS_2003.tif"))
LC_CONUS_2004 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "LC_CONUS_2004.tif"))
LC_CONUS_2005 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "LC_CONUS_2005.tif"))
LC_CONUS_2006 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "LC_CONUS_2006.tif"))
LC_CONUS_2007 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "LC_CONUS_2007.tif"))
LC_CONUS_2008 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "LC_CONUS_2008.tif"))
LC_CONUS_2009 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "LC_CONUS_2009.tif"))
LC_CONUS_2010 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "LC_CONUS_2010.tif"))
LC_CONUS_2011 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "LC_CONUS_2011.tif"))
LC_CONUS_2012 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "LC_CONUS_2012.tif"))
LC_CONUS_2013 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "LC_CONUS_2013.tif"))
LC_CONUS_2014 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "LC_CONUS_2014.tif"))
LC_CONUS_2015 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "LC_CONUS_2015.tif"))
LC_CONUS_2016 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "LC_CONUS_2016.tif"))
LC_CONUS_2017 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "LC_CONUS_2017.tif"))
LC_CONUS_2018 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "LC_CONUS_2018.tif"))
LC_CONUS_2019 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "LC_CONUS_2019.tif"))
LC_CONUS_2020 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "LC_CONUS_2020.tif"))


### making a list of landcover rasters for when landcover classes get extracted
lc_CONUS_rast_list <- list(LC_CONUS_1993, LC_CONUS_1994, LC_CONUS_1995, LC_CONUS_1996, LC_CONUS_1997, LC_CONUS_1998, LC_CONUS_1999, LC_CONUS_2000, LC_CONUS_2001, LC_CONUS_2002, LC_CONUS_2003, LC_CONUS_2004, LC_CONUS_2005, LC_CONUS_2006, LC_CONUS_2007, LC_CONUS_2008, LC_CONUS_2009, LC_CONUS_2010, LC_CONUS_2011, LC_CONUS_2012, LC_CONUS_2013, LC_CONUS_2014, LC_CONUS_2015, LC_CONUS_2016, LC_CONUS_2017, LC_CONUS_2018, LC_CONUS_2019, LC_CONUS_2020)
### concatenating landcover rasters for landcover climatology
## CONUS
# lc_CONUS_rasts <- c(LC_CONUS_1992, LC_CONUS_1993, LC_CONUS_1994, LC_CONUS_1995, LC_CONUS_1996, LC_CONUS_1997, LC_CONUS_1998, LC_CONUS_1999, LC_CONUS_2000, LC_CONUS_2001, LC_CONUS_2002, LC_CONUS_2003, LC_CONUS_2004, LC_CONUS_2005, LC_CONUS_2006, LC_CONUS_2007, LC_CONUS_2008, LC_CONUS_2009, LC_CONUS_2010, LC_CONUS_2011, LC_CONUS_2012, LC_CONUS_2013, LC_CONUS_2014, LC_CONUS_2015, LC_CONUS_2016, LC_CONUS_2017, LC_CONUS_2018, LC_CONUS_2019, LC_CONUS_2020)

### creating single raster which represents most common landcover class for each cell across all years in study period
# LC_CONUS_Climatology_rast <- app(lc_CONUS_rasts, fun = "modal")
# writeRaster(LC_CONUS_Climatology_rast, here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Climatology", "LC_CONUS_Climatology.tif"), overwrite = TRUE)

##### reading in RECLASSIFIED CONUS Landcover #####
### CONUS
LC_CONUS_TriClass_1992 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_1992.tif"))
LC_CONUS_TriClass_1993 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_1993.tif"))
LC_CONUS_TriClass_1994 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_1994.tif"))
LC_CONUS_TriClass_1995 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_1995.tif"))
LC_CONUS_TriClass_1996 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_1996.tif"))
LC_CONUS_TriClass_1997 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_1997.tif"))
LC_CONUS_TriClass_1998 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_1998.tif"))
LC_CONUS_TriClass_1999 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_1999.tif"))
LC_CONUS_TriClass_2000 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2000.tif"))
LC_CONUS_TriClass_2001 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2001.tif"))
LC_CONUS_TriClass_2002 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2002.tif"))
LC_CONUS_TriClass_2003 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2003.tif"))
LC_CONUS_TriClass_2004 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2004.tif"))
LC_CONUS_TriClass_2005 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2005.tif"))
LC_CONUS_TriClass_2006 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2006.tif"))
LC_CONUS_TriClass_2007 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2007.tif"))
LC_CONUS_TriClass_2008 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2008.tif"))
LC_CONUS_TriClass_2009 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2009.tif"))
LC_CONUS_TriClass_2010 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2010.tif"))
LC_CONUS_TriClass_2011 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2011.tif"))
LC_CONUS_TriClass_2012 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2012.tif"))
LC_CONUS_TriClass_2013 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2013.tif"))
LC_CONUS_TriClass_2014 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2014.tif"))
LC_CONUS_TriClass_2015 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2015.tif"))
LC_CONUS_TriClass_2016 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2016.tif"))
LC_CONUS_TriClass_2017 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2017.tif"))
LC_CONUS_TriClass_2018 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2018.tif"))
LC_CONUS_TriClass_2019 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2019.tif"))
LC_CONUS_TriClass_2020 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2020.tif"))

### making a list of landcover rasters for when landcover classes get extracted
LC_CONUS_TriClass_rast_list <- list(LC_CONUS_TriClass_1993, LC_CONUS_TriClass_1994, LC_CONUS_TriClass_1995, LC_CONUS_TriClass_1996, LC_CONUS_TriClass_1997, LC_CONUS_TriClass_1998, LC_CONUS_TriClass_1999, LC_CONUS_TriClass_2000, LC_CONUS_TriClass_2001, LC_CONUS_TriClass_2002, LC_CONUS_TriClass_2003, LC_CONUS_TriClass_2004, LC_CONUS_TriClass_2005, LC_CONUS_TriClass_2006, LC_CONUS_TriClass_2007, LC_CONUS_TriClass_2008, LC_CONUS_TriClass_2009, LC_CONUS_TriClass_2010, LC_CONUS_TriClass_2011, LC_CONUS_TriClass_2012, LC_CONUS_TriClass_2013, LC_CONUS_TriClass_2014, LC_CONUS_TriClass_2015, LC_CONUS_TriClass_2016, LC_CONUS_TriClass_2017, LC_CONUS_TriClass_2018, LC_CONUS_TriClass_2019, LC_CONUS_TriClass_2020)
### concatenating landcover rasters for landcover climatology
## CONUS
# LC_CONUS_TriClass_rasts <- c(LC_CONUS_TriClass_1992, LC_CONUS_TriClass_1993, LC_CONUS_TriClass_1994, LC_CONUS_TriClass_1995, LC_CONUS_TriClass_1996, LC_CONUS_TriClass_1997, LC_CONUS_TriClass_1998, LC_CONUS_TriClass_1999, LC_CONUS_TriClass_2000, LC_CONUS_TriClass_2001, LC_CONUS_TriClass_2002, LC_CONUS_TriClass_2003, LC_CONUS_TriClass_2004, LC_CONUS_TriClass_2005, LC_CONUS_TriClass_2006, LC_CONUS_TriClass_2007, LC_CONUS_TriClass_2008, LC_CONUS_TriClass_2009, LC_CONUS_TriClass_2010, LC_CONUS_TriClass_2011, LC_CONUS_TriClass_2012, LC_CONUS_TriClass_2013, LC_CONUS_TriClass_2014, LC_CONUS_TriClass_2015, LC_CONUS_TriClass_2016, LC_CONUS_TriClass_2017, LC_CONUS_TriClass_2018, LC_CONUS_TriClass_2019, LC_CONUS_TriClass_2020)

### creating single raster which represents most common landcover class for each cell across all years in study period
# LC_CONUS_TriClass_Climatology_rast <- app(LC_CONUS_TriClass_rasts, fun = "modal")
# writeRaster(LC_CONUS_TriClass_Climatology_rast, here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "Climatology", "LC_CONUS_TriClass_Climatology.tif"), overwrite = TRUE)
# LC_CONUS_TriClass_Climatology_rast <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "Climatology", "LC_CONUS_Climatology_Reclass.tif"))
# 
# ##### Generating Random "Zero points" to be assigned Max SWE value of 0 #####
# ### useful for reducing or eliminating overestimation of SWE in no/low snow areas
# ### Created Raster of no snow areas using output from SnowCloudMetrics GEE code
# SCF_Climatology <- rast(here("Data", "SCF", "Resampled", "SCF2003_2020.tif"))
# SCF_Climatology_NoSnow <- rast(here("Data", "SCF", "Resampled", "SCF2003_2020NoSnow.tif"))
# 
# ### Randomly sampling points within areas of no snow
# set.seed(802)
# Zero_points <- spatSample(SCF_Climatology_NoSnow, size = nrow(Snotel_CONUS), method = "random", replace = FALSE, na.rm = TRUE, as.points = TRUE, values = FALSE)
# 
# Zero_pts_sf <- st_as_sf(Zero_points)
# st_geometry(Zero_pts_sf) <- "geom"
# Zero_pts_sf <- Zero_pts_sf |>
#   mutate(peak_swe = 0,
#          site_id = -999,
#          site_name = "zero values",
#          description = "randomly generated points to be assigned 0 for peak SWE",
#          state = "randomly generated points to be assigned 0 for peak SWE",
#          start = min(Snotel_CONUS$start),
#          end = max(Snotel_CONUS$end), .before = "geom")
# 
# ### making plot of points on SCF raster as a check
# plot(SCF_Climatology, main = "SCF Climatology, SNOTEL pts in red")
# points(vect(Snotel_CONUS), col = "red")
# points(vect(Zero_pts_sf), col = "black")
# 
# ### writing out zero points
# write_sf(Zero_pts_sf, here("Data", "SNOTEL", "ZeroPts", "CONUSZeroPts.gpkg"), append = FALSE)
# 
# ### binding zero points to Actual SNOTEL data
# Snotel_CONUS <- rbind(Snotel_CONUS, Zero_pts_sf)
# 
# ##### Extracting climatology and topography values to CONUS SNOTEL points #####
# ### extracting topographic variables
# Snotel_CONUS$elevation <- terra::extract(CONUS_DEM, Snotel_CONUS)[,2]
# Snotel_CONUS$slope <- terra::extract(CONUS_Slope, Snotel_CONUS)[,2]
# Snotel_CONUS$aspect <- terra::extract(CONUS_Aspect, Snotel_CONUS)[,2]
# ### extracting climatologies
# ## prcpSum
# Snotel_CONUS$jan_prcpSum <- terra::extract(prcpSum_JanClim, Snotel_CONUS)[,2]
# Snotel_CONUS$feb_prcpSum <- terra::extract(prcpSum_FebClim, Snotel_CONUS)[,2]
# Snotel_CONUS$mar_prcpSum <- terra::extract(prcpSum_MarClim, Snotel_CONUS)[,2]
# Snotel_CONUS$apr_prcpSum <- terra::extract(prcpSum_AprClim, Snotel_CONUS)[,2]
# Snotel_CONUS$may_prcpSum <- terra::extract(prcpSum_MayClim, Snotel_CONUS)[,2]
# Snotel_CONUS$june_prcpSum <- terra::extract(prcpSum_JuneClim, Snotel_CONUS)[,2]
# Snotel_CONUS$july_prcpSum <- terra::extract(prcpSum_JulyClim, Snotel_CONUS)[,2]
# Snotel_CONUS$aug_prcpSum <- terra::extract(prcpSum_AugClim, Snotel_CONUS)[,2]
# Snotel_CONUS$sep_prcpSum <- terra::extract(prcpSum_SepClim, Snotel_CONUS)[,2]
# Snotel_CONUS$oct_prcpSum <- terra::extract(prcpSum_OctClim, Snotel_CONUS)[,2]
# Snotel_CONUS$nov_prcpSum <- terra::extract(prcpSum_NovClim, Snotel_CONUS)[,2]
# Snotel_CONUS$dec_prcpSum <- terra::extract(prcpSum_DecClim, Snotel_CONUS)[,2]
# ## prcpMean
# Snotel_CONUS$jan_prcpMean <- terra::extract(prcpMean_JanClim, Snotel_CONUS)[,2]
# Snotel_CONUS$feb_prcpMean <- terra::extract(prcpMean_FebClim, Snotel_CONUS)[,2]
# Snotel_CONUS$mar_prcpMean <- terra::extract(prcpMean_MarClim, Snotel_CONUS)[,2]
# Snotel_CONUS$apr_prcpMean <- terra::extract(prcpMean_AprClim, Snotel_CONUS)[,2]
# Snotel_CONUS$may_prcpMean <- terra::extract(prcpMean_MayClim, Snotel_CONUS)[,2]
# Snotel_CONUS$june_prcpMean <- terra::extract(prcpMean_JuneClim, Snotel_CONUS)[,2]
# Snotel_CONUS$july_prcpMean <- terra::extract(prcpMean_JulyClim, Snotel_CONUS)[,2]
# Snotel_CONUS$aug_prcpMean <- terra::extract(prcpMean_AugClim, Snotel_CONUS)[,2]
# Snotel_CONUS$sep_prcpMean <- terra::extract(prcpMean_SepClim, Snotel_CONUS)[,2]
# Snotel_CONUS$oct_prcpMean <- terra::extract(prcpMean_OctClim, Snotel_CONUS)[,2]
# Snotel_CONUS$nov_prcpMean <- terra::extract(prcpMean_NovClim, Snotel_CONUS)[,2]
# Snotel_CONUS$dec_prcpMean <- terra::extract(prcpMean_DecClim, Snotel_CONUS)[,2]
# 
# ## tmin
# Snotel_CONUS$jan_tmin <- terra::extract(tmin_JanClim, Snotel_CONUS)[,2]
# Snotel_CONUS$feb_tmin <- terra::extract(tmin_FebClim, Snotel_CONUS)[,2]
# Snotel_CONUS$mar_tmin <- terra::extract(tmin_MarClim, Snotel_CONUS)[,2]
# Snotel_CONUS$apr_tmin <- terra::extract(tmin_AprClim, Snotel_CONUS)[,2]
# Snotel_CONUS$may_tmin <- terra::extract(tmin_MayClim, Snotel_CONUS)[,2]
# Snotel_CONUS$june_tmin <- terra::extract(tmin_JuneClim, Snotel_CONUS)[,2]
# Snotel_CONUS$july_tmin <- terra::extract(tmin_JulyClim, Snotel_CONUS)[,2]
# Snotel_CONUS$aug_tmin <- terra::extract(tmin_AugClim, Snotel_CONUS)[,2]
# Snotel_CONUS$sep_tmin <- terra::extract(tmin_SepClim, Snotel_CONUS)[,2]
# Snotel_CONUS$oct_tmin <- terra::extract(tmin_OctClim, Snotel_CONUS)[,2]
# Snotel_CONUS$nov_tmin <- terra::extract(tmin_NovClim, Snotel_CONUS)[,2]
# Snotel_CONUS$dec_tmin <- terra::extract(tmin_DecClim, Snotel_CONUS)[,2]
# 
# ## tmean
# Snotel_CONUS$jan_tmean <- terra::extract(tmean_JanClim, Snotel_CONUS)[,2]
# Snotel_CONUS$feb_tmean <- terra::extract(tmean_FebClim, Snotel_CONUS)[,2]
# Snotel_CONUS$mar_tmean <- terra::extract(tmean_MarClim, Snotel_CONUS)[,2]
# Snotel_CONUS$apr_tmean <- terra::extract(tmean_AprClim, Snotel_CONUS)[,2]
# Snotel_CONUS$may_tmean <- terra::extract(tmean_MayClim, Snotel_CONUS)[,2]
# Snotel_CONUS$june_tmean <- terra::extract(tmean_JuneClim, Snotel_CONUS)[,2]
# Snotel_CONUS$july_tmean <- terra::extract(tmean_JulyClim, Snotel_CONUS)[,2]
# Snotel_CONUS$aug_tmean <- terra::extract(tmean_AugClim, Snotel_CONUS)[,2]
# Snotel_CONUS$sep_tmean <- terra::extract(tmean_SepClim, Snotel_CONUS)[,2]
# Snotel_CONUS$oct_tmean <- terra::extract(tmean_OctClim, Snotel_CONUS)[,2]
# Snotel_CONUS$nov_tmean <- terra::extract(tmean_NovClim, Snotel_CONUS)[,2]
# Snotel_CONUS$dec_tmean <- terra::extract(tmean_DecClim, Snotel_CONUS)[,2]
# 
# ## tmax
# Snotel_CONUS$jan_tmax <- terra::extract(tmax_JanClim, Snotel_CONUS)[,2]
# Snotel_CONUS$feb_tmax <- terra::extract(tmax_FebClim, Snotel_CONUS)[,2]
# Snotel_CONUS$mar_tmax <- terra::extract(tmax_MarClim, Snotel_CONUS)[,2]
# Snotel_CONUS$apr_tmax <- terra::extract(tmax_AprClim, Snotel_CONUS)[,2]
# Snotel_CONUS$may_tmax <- terra::extract(tmax_MayClim, Snotel_CONUS)[,2]
# Snotel_CONUS$june_tmax <- terra::extract(tmax_JuneClim, Snotel_CONUS)[,2]
# Snotel_CONUS$july_tmax <- terra::extract(tmax_JulyClim, Snotel_CONUS)[,2]
# Snotel_CONUS$aug_tmax <- terra::extract(tmax_AugClim, Snotel_CONUS)[,2]
# Snotel_CONUS$sep_tmax <- terra::extract(tmax_SepClim, Snotel_CONUS)[,2]
# Snotel_CONUS$oct_tmax <- terra::extract(tmax_OctClim, Snotel_CONUS)[,2]
# Snotel_CONUS$nov_tmax <- terra::extract(tmax_NovClim, Snotel_CONUS)[,2]
# Snotel_CONUS$dec_tmax <- terra::extract(tmax_DecClim, Snotel_CONUS)[,2]
# 
# ## srad
# Snotel_CONUS$jan_srad <- terra::extract(srad_JanClim, Snotel_CONUS)[,2]
# Snotel_CONUS$feb_srad <- terra::extract(srad_FebClim, Snotel_CONUS)[,2]
# Snotel_CONUS$mar_srad <- terra::extract(srad_MarClim, Snotel_CONUS)[,2]
# Snotel_CONUS$apr_srad <- terra::extract(srad_AprClim, Snotel_CONUS)[,2]
# Snotel_CONUS$may_srad <- terra::extract(srad_MayClim, Snotel_CONUS)[,2]
# Snotel_CONUS$june_srad <- terra::extract(srad_JuneClim, Snotel_CONUS)[,2]
# Snotel_CONUS$july_srad <- terra::extract(srad_JulyClim, Snotel_CONUS)[,2]
# Snotel_CONUS$aug_srad <- terra::extract(srad_AugClim, Snotel_CONUS)[,2]
# Snotel_CONUS$sep_srad <- terra::extract(srad_SepClim, Snotel_CONUS)[,2]
# Snotel_CONUS$oct_srad <- terra::extract(srad_OctClim, Snotel_CONUS)[,2]
# Snotel_CONUS$nov_srad <- terra::extract(srad_NovClim, Snotel_CONUS)[,2]
# Snotel_CONUS$dec_srad <- terra::extract(srad_DecClim, Snotel_CONUS)[,2]
# 
# 
# ### Extracting climatology/CDM aggregates
# ## Oct-May aggregates
# Snotel_CONUS$OctMay_prcpSumCDMSum <- terra::extract(prcpSum_OctMayClimSum_rast, Snotel_CONUS)[,2]
# Snotel_CONUS$OctMay_tmeanCDMSum <- terra::extract(tmean_OctMayClimSum_rast, Snotel_CONUS)[,2]
# Snotel_CONUS$OctMay_tminMean <- terra::extract(tmin_OctMayClimMean_rast, Snotel_CONUS)[,2]
# Snotel_CONUS$OctMay_tmaxMean <- terra::extract(tmax_OctMayClimMean_rast, Snotel_CONUS)[,2]
# Snotel_CONUS$OctMay_sradMean <- terra::extract(srad_OctMayClimMean_rast, Snotel_CONUS)[,2]
# ## Oct-Apr aggregates
# Snotel_CONUS$OctApr_prcpSumCDMSum <- terra::extract(prcpSum_OctAprClimSum_rast, Snotel_CONUS)[,2]
# Snotel_CONUS$OctApr_tmeanCDMSum <- terra::extract(tmean_OctAprClimSum_rast, Snotel_CONUS)[,2]
# Snotel_CONUS$OctApr_tminMean <- terra::extract(tmin_OctAprClimMean_rast, Snotel_CONUS)[,2]
# Snotel_CONUS$OctApr_tmaxMean <- terra::extract(tmax_OctAprClimMean_rast, Snotel_CONUS)[,2]
# Snotel_CONUS$OctApr_sradMean <- terra::extract(srad_OctAprClimMean_rast, Snotel_CONUS)[,2]
# ## Sep-May aggregates
# Snotel_CONUS$SepMay_prcpSumCDMSum <- terra::extract(prcpSum_SepMayClimSum_rast, Snotel_CONUS)[,2]
# Snotel_CONUS$SepMay_tmeanCDMSum <- terra::extract(tmean_SepMayClimSum_rast, Snotel_CONUS)[,2]
# Snotel_CONUS$SepMay_tminMean <- terra::extract(tmin_SepMayClimMean_rast, Snotel_CONUS)[,2]
# Snotel_CONUS$SepMay_tmaxMean <- terra::extract(tmax_SepMayClimMean_rast, Snotel_CONUS)[,2]
# Snotel_CONUS$SepMay_sradMean <- terra::extract(srad_SepMayClimMean_rast, Snotel_CONUS)[,2]
# ## Dec-Feb aggregates
# Snotel_CONUS$DecFeb_prcpSumCDMSum <- terra::extract(prcpSum_DecFebClimSum_rast, Snotel_CONUS)[,2]
# Snotel_CONUS$DecFeb_tmeanCDMSum <- terra::extract(tmean_DecFebClimSum_rast, Snotel_CONUS)[,2]
# Snotel_CONUS$DecFeb_tminMean <- terra::extract(tmin_DecFebClimMean_rast, Snotel_CONUS)[,2]
# Snotel_CONUS$DecFeb_tmaxMean <- terra::extract(tmax_DecFebClimMean_rast, Snotel_CONUS)[,2]
# Snotel_CONUS$DecFeb_sradMean <- terra::extract(srad_DecFebClimMean_rast, Snotel_CONUS)[,2]
# ### Extracting Climatological (modal) land cover value from land cover raster
# Snotel_CONUS$landcover <- terra::extract(LC_CONUS_Climatology_rast, Snotel_CONUS)[,2]
# Snotel_CONUS$landcover_triclass <- terra::extract(LC_CONUS_TriClass_Climatology_rast, Snotel_CONUS)[,2]
# 
# 
# ##### writing out SNOTEL CONUS SWE Climatology shapefile/gpkg with covariates #####
# # write_sf(Snotel_CONUS, here("Data", "SNOTEL", "CONUS", "GIS", "SHP", "SnotelCONUS_PeakSWEClimatology_Covars.shp"))
# write_sf(Snotel_CONUS, here("Data", "SNOTEL", "CONUS", "GIS", "GPKG", "Climatology", "SnotelCONUS_PeakSWEClimatology_Covars.gpkg"), append = FALSE)
# CONUSClimendtime <- Sys.time()
# CONUSClimtime <- CONUSClimendtime - CONUSClimstarttime




##### Alaska Climatology #####
# AKClimStartTime <- Sys.time()
# ### reading in Alaska data
# Snotel_AK <- read_sf(here("Data", "SNOTEL", "Alaska", "GIS", "GPKG", "Climatology", "SnotelData_AK_PeakSWEClimatology.gpkg"))
AK_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "Alaska", "AK_AOI.gpkg"))

##### Reading in Alaska Climatology rasters #####
### prcpSum
# prcpSum_JanClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "prcpSum_Month1MeanClimatology.tif"))
# prcpSum_FebClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "prcpSum_Month2MeanClimatology.tif"))
# prcpSum_MarClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "prcpSum_Month3MeanClimatology.tif"))
# prcpSum_AprClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "prcpSum_Month4MeanClimatology.tif"))
# prcpSum_MayClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "prcpSum_Month5MeanClimatology.tif"))
# prcpSum_JuneClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "prcpSum_Month6MeanClimatology.tif"))
# prcpSum_JulyClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "prcpSum_Month7MeanClimatology.tif"))
# prcpSum_AugClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "prcpSum_Month8MeanClimatology.tif"))
# prcpSum_SepClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "prcpSum_Month9MeanClimatology.tif"))
# prcpSum_OctClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "prcpSum_Month10MeanClimatology.tif"))
# prcpSum_NovClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "prcpSum_Month11MeanClimatology.tif"))
# prcpSum_DecClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "prcpSum_Month12MeanClimatology.tif"))
# ### prcpMean
# prcpMean_JanClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "prcp_Month1MeanClimatology.tif"))
# prcpMean_FebClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "prcp_Month2MeanClimatology.tif"))
# prcpMean_MarClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "prcp_Month3MeanClimatology.tif"))
# prcpMean_AprClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "prcp_Month4MeanClimatology.tif"))
# prcpMean_MayClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "prcp_Month5MeanClimatology.tif"))
# prcpMean_JuneClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "prcp_Month6MeanClimatology.tif"))
# prcpMean_JulyClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "prcp_Month7MeanClimatology.tif"))
# prcpMean_AugClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "prcp_Month8MeanClimatology.tif"))
# prcpMean_SepClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "prcp_Month9MeanClimatology.tif"))
# prcpMean_OctClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "prcp_Month10MeanClimatology.tif"))
# prcpMean_NovClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "prcp_Month11MeanClimatology.tif"))
# prcpMean_DecClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "prcp_Month12MeanClimatology.tif"))
# ## tmin
# tmin_JanClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmin_Month1MeanClimatology.tif"))
# tmin_FebClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmin_Month2MeanClimatology.tif"))
# tmin_MarClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmin_Month3MeanClimatology.tif"))
# tmin_AprClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmin_Month4MeanClimatology.tif"))
# tmin_MayClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmin_Month5MeanClimatology.tif"))
# tmin_JuneClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmin_Month6MeanClimatology.tif"))
# tmin_JulyClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmin_Month7MeanClimatology.tif"))
# tmin_AugClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmin_Month8MeanClimatology.tif"))
# tmin_SepClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmin_Month9MeanClimatology.tif"))
# tmin_OctClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmin_Month10MeanClimatology.tif"))
# tmin_NovClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmin_Month11MeanClimatology.tif"))
# tmin_DecClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmin_Month12MeanClimatology.tif"))
# 
# ## tmean
# tmean_JanClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmean_Month1MeanClimatology.tif"))
# tmean_FebClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmean_Month2MeanClimatology.tif"))
# tmean_MarClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmean_Month3MeanClimatology.tif"))
# tmean_AprClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmean_Month4MeanClimatology.tif"))
# tmean_MayClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmean_Month5MeanClimatology.tif"))
# tmean_JuneClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmean_Month6MeanClimatology.tif"))
# tmean_JulyClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmean_Month7MeanClimatology.tif"))
# tmean_AugClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmean_Month8MeanClimatology.tif"))
# tmean_SepClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmean_Month9MeanClimatology.tif"))
# tmean_OctClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmean_Month10MeanClimatology.tif"))
# tmean_NovClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmean_Month11MeanClimatology.tif"))
# tmean_DecClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmean_Month12MeanClimatology.tif"))
# 
# ## tmax
# tmax_JanClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmax_Month1MeanClimatology.tif"))
# tmax_FebClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmax_Month2MeanClimatology.tif"))
# tmax_MarClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmax_Month3MeanClimatology.tif"))
# tmax_AprClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmax_Month4MeanClimatology.tif"))
# tmax_MayClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmax_Month5MeanClimatology.tif"))
# tmax_JuneClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmax_Month6MeanClimatology.tif"))
# tmax_JulyClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmax_Month7MeanClimatology.tif"))
# tmax_AugClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmax_Month8MeanClimatology.tif"))
# tmax_SepClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmax_Month9MeanClimatology.tif"))
# tmax_OctClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmax_Month10MeanClimatology.tif"))
# tmax_NovClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmax_Month11MeanClimatology.tif"))
# tmax_DecClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "tmax_Month12MeanClimatology.tif"))
# 
# ## srad
# srad_JanClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "srad_Month1MeanClimatology.tif"))
# srad_FebClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "srad_Month2MeanClimatology.tif"))
# srad_MarClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "srad_Month3MeanClimatology.tif"))
# srad_AprClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "srad_Month4MeanClimatology.tif"))
# srad_MayClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "srad_Month5MeanClimatology.tif"))
# srad_JuneClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "srad_Month6MeanClimatology.tif"))
# srad_JulyClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "srad_Month7MeanClimatology.tif"))
# srad_AugClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "srad_Month8MeanClimatology.tif"))
# srad_SepClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "srad_Month9MeanClimatology.tif"))
# srad_OctClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "srad_Month10MeanClimatology.tif"))
# srad_NovClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "srad_Month11MeanClimatology.tif"))
# srad_DecClim <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "srad_Month12MeanClimatology.tif"))
# 
# 
# ##### reading in Alaska climatology CDMs #####
# ### Oct-Apr Sums
# prcpSum_OctAprClimSum_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpSum_OctAprSum_ClimatologyCDM.tif"))
# prcpMean_OctAprClimMean_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpMean_OctAprMean_ClimatologyCDM.tif"))
# tmean_OctAprClimSum_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmean_OctAprSum_ClimatologyCDM.tif"))
# tmin_OctAprClimMean_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmin_OctAprMean_ClimatologyCDM.tif"))
# tmax_OctAprClimMean_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmax_OctAprMean_ClimatologyCDM.tif"))
# srad_OctAprClimMean_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "srad_OctAprMean_ClimatologyCDM.tif"))
# ### Oct-May Sums
# prcpSum_OctMayClimSum_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpSum_OctMaySum_ClimatologyCDM.tif"))
# prcpMean_OctMayClimMean_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpMean_OctMayMean_ClimatologyCDM.tif"))
# tmean_OctMayClimSum_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmean_OctMaySum_ClimatologyCDM.tif"))
# tmin_OctMayClimMean_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmin_OctMayMean_ClimatologyCDM.tif"))
# tmax_OctMayClimMean_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmax_OctMayMean_ClimatologyCDM.tif"))
# srad_OctMayClimMean_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "srad_OctMayMean_ClimatologyCDM.tif"))
# ### Sep-May Sums
# prcpSum_SepMayClimSum_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpSum_SepMaySum_ClimatologyCDM.tif"))
# prcpMean_SepMayClimMean_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpMean_SepMayMean_ClimatologyCDM.tif"))
# tmean_SepMayClimSum_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmean_SepMaySum_ClimatologyCDM.tif"))
# tmin_SepMayClimMean_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmin_SepMayMean_ClimatologyCDM.tif"))
# tmax_SepMayClimMean_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmax_SepMayMean_ClimatologyCDM.tif"))
# srad_SepMayClimMean_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "srad_SepMayMean_ClimatologyCDM.tif"))
# ### Dec-Feb Sums
# prcpSum_DecFebClimSum_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpSum_DecFebSum_ClimatologyCDM.tif"))
# prcpMean_DecFebClimMean_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "prcpMean_DecFebMean_ClimatologyCDM.tif"))
# tmean_DecFebClimSum_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmean_DecFebSum_ClimatologyCDM.tif"))
# tmin_DecFebClimMean_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmin_DecFebMean_ClimatologyCDM.tif"))
# tmax_DecFebClimMean_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "tmax_DecFebMean_ClimatologyCDM.tif"))
# srad_DecFebClimMean_rast <- rast(here("Data", "Daymet", "Alaska", "Climatology", "Resampled", "CDM", "CDM_Aggregates", "srad_DecFebMean_ClimatologyCDM.tif"))
# 
##### Reading in Alaska topographical variables #####
### Alaska
AK_DEM <- rast(here("Data", "DEM", "Alaska", "AKDEMMosaic.tif"))
AK_Slope <- rast(here("Data", "DEM", "Alaska", "AKSlope.tif"))
AK_Aspect <- rast(here("Data", "DEM", "Alaska", "AKAspect.tif"))


##### Reading in Alaska landcover #####
### Alaska
LC_AK_1992 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "LC_AK_1992.tif"))
LC_AK_1993 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "LC_AK_1993.tif"))
LC_AK_1994 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "LC_AK_1994.tif"))
LC_AK_1995 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "LC_AK_1995.tif"))
LC_AK_1996 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "LC_AK_1996.tif"))
LC_AK_1997 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "LC_AK_1997.tif"))
LC_AK_1998 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "LC_AK_1998.tif"))
LC_AK_1999 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "LC_AK_1999.tif"))
LC_AK_2000 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "LC_AK_2000.tif"))
LC_AK_2001 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "LC_AK_2001.tif"))
LC_AK_2002 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "LC_AK_2002.tif"))
LC_AK_2003 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "LC_AK_2003.tif"))
LC_AK_2004 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "LC_AK_2004.tif"))
LC_AK_2005 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "LC_AK_2005.tif"))
LC_AK_2006 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "LC_AK_2006.tif"))
LC_AK_2007 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "LC_AK_2007.tif"))
LC_AK_2008 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "LC_AK_2008.tif"))
LC_AK_2009 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "LC_AK_2009.tif"))
LC_AK_2010 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "LC_AK_2010.tif"))
LC_AK_2011 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "LC_AK_2011.tif"))
LC_AK_2012 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "LC_AK_2012.tif"))
LC_AK_2013 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "LC_AK_2013.tif"))
LC_AK_2014 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "LC_AK_2014.tif"))
LC_AK_2015 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "LC_AK_2015.tif"))
LC_AK_2016 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "LC_AK_2016.tif"))
LC_AK_2017 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "LC_AK_2017.tif"))
LC_AK_2018 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "LC_AK_2018.tif"))
LC_AK_2019 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "LC_AK_2019.tif"))
LC_AK_2020 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "LC_AK_2020.tif"))

### making a list of landcover rasters for when landcover classes get extracted
lc_AK_rast_list <- list(LC_AK_1993, LC_AK_1994, LC_AK_1995, LC_AK_1996, LC_AK_1997, LC_AK_1998, LC_AK_1999, LC_AK_2000, LC_AK_2001, LC_AK_2002, LC_AK_2003, LC_AK_2004, LC_AK_2005, LC_AK_2006, LC_AK_2007, LC_AK_2008, LC_AK_2009, LC_AK_2010, LC_AK_2011, LC_AK_2012, LC_AK_2013, LC_AK_2014, LC_AK_2015, LC_AK_2016, LC_AK_2017, LC_AK_2018, LC_AK_2019, LC_AK_2020)
### stacking landcover rasters for land cover climatology
# lc_AK_rasts <- c(LC_AK_1992, LC_AK_1993, LC_AK_1994, LC_AK_1995, LC_AK_1996, LC_AK_1997, LC_AK_1998, LC_AK_1999, LC_AK_2000, LC_AK_2001, LC_AK_2002, LC_AK_2003, LC_AK_2004, LC_AK_2005, LC_AK_2006, LC_AK_2007, LC_AK_2008, LC_AK_2009, LC_AK_2010, LC_AK_2011, LC_AK_2012, LC_AK_2013, LC_AK_2014, LC_AK_2015, LC_AK_2016, LC_AK_2017, LC_AK_2018, LC_AK_2019, LC_AK_2020)
### creating landcover raster with each cell representing most common landcover type across the study period
# LC_AK_Climatology_rast <- app(lc_AK_rasts, fun = "modal")
# writeRaster(LC_AK_Climatology_rast, here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Climatology", "LC_AK_Climatology.tif"), overwrite = TRUE)


##### reading in RECLASSIFIED Alaska Landcover #####
### Alaska
LC_AK_TriClass_1992 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_1992.tif"))
LC_AK_TriClass_1993 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_1993.tif"))
LC_AK_TriClass_1994 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_1994.tif"))
LC_AK_TriClass_1995 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_1995.tif"))
LC_AK_TriClass_1996 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_1996.tif"))
LC_AK_TriClass_1997 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_1997.tif"))
LC_AK_TriClass_1998 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_1998.tif"))
LC_AK_TriClass_1999 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_1999.tif"))
LC_AK_TriClass_2000 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2000.tif"))
LC_AK_TriClass_2001 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2001.tif"))
LC_AK_TriClass_2002 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2002.tif"))
LC_AK_TriClass_2003 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2003.tif"))
LC_AK_TriClass_2004 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2004.tif"))
LC_AK_TriClass_2005 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2005.tif"))
LC_AK_TriClass_2006 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2006.tif"))
LC_AK_TriClass_2007 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2007.tif"))
LC_AK_TriClass_2008 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2008.tif"))
LC_AK_TriClass_2009 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2009.tif"))
LC_AK_TriClass_2010 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2010.tif"))
LC_AK_TriClass_2011 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2011.tif"))
LC_AK_TriClass_2012 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2012.tif"))
LC_AK_TriClass_2013 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2013.tif"))
LC_AK_TriClass_2014 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2014.tif"))
LC_AK_TriClass_2015 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2015.tif"))
LC_AK_TriClass_2016 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2016.tif"))
LC_AK_TriClass_2017 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2017.tif"))
LC_AK_TriClass_2018 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2018.tif"))
LC_AK_TriClass_2019 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2019.tif"))
LC_AK_TriClass_2020 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2020.tif"))

### making a list of landcover rasters for when landcover classes get extracted
LC_AK_TriClass_rast_list <- list(LC_AK_TriClass_1993, LC_AK_TriClass_1994, LC_AK_TriClass_1995, LC_AK_TriClass_1996, LC_AK_TriClass_1997, LC_AK_TriClass_1998, LC_AK_TriClass_1999, LC_AK_TriClass_2000, LC_AK_TriClass_2001, LC_AK_TriClass_2002, LC_AK_TriClass_2003, LC_AK_TriClass_2004, LC_AK_TriClass_2005, LC_AK_TriClass_2006, LC_AK_TriClass_2007, LC_AK_TriClass_2008, LC_AK_TriClass_2009, LC_AK_TriClass_2010, LC_AK_TriClass_2011, LC_AK_TriClass_2012, LC_AK_TriClass_2013, LC_AK_TriClass_2014, LC_AK_TriClass_2015, LC_AK_TriClass_2016, LC_AK_TriClass_2017, LC_AK_TriClass_2018, LC_AK_TriClass_2019, LC_AK_TriClass_2020)
### concatenating landcover rasters for landcover climatology
## Alaska
# LC_AK_TriClass_rasts <- c(LC_AK_TriClass_1992, LC_AK_TriClass_1993, LC_AK_TriClass_1994, LC_AK_TriClass_1995, LC_AK_TriClass_1996, LC_AK_TriClass_1997, LC_AK_TriClass_1998, LC_AK_TriClass_1999, LC_AK_TriClass_2000, LC_AK_TriClass_2001, LC_AK_TriClass_2002, LC_AK_TriClass_2003, LC_AK_TriClass_2004, LC_AK_TriClass_2005, LC_AK_TriClass_2006, LC_AK_TriClass_2007, LC_AK_TriClass_2008, LC_AK_TriClass_2009, LC_AK_TriClass_2010, LC_AK_TriClass_2011, LC_AK_TriClass_2012, LC_AK_TriClass_2013, LC_AK_TriClass_2014, LC_AK_TriClass_2015, LC_AK_TriClass_2016, LC_AK_TriClass_2017, LC_AK_TriClass_2018, LC_AK_TriClass_2019, LC_AK_TriClass_2020)

### creating single raster which represents most common landcover class for each cell across all years in study period
# # LC_AK_TriClass_Climatology_rast <- app(LC_AK_TriClass_rasts, fun = "modal")
# # writeRaster(LC_AK_TriClass_Climatology_rast, here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "Climatology", "LC_AK_TriClass_Climatology.tif"), overwrite = TRUE)
# LC_AK_TriClass_Climatology_rast <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "Climatology", "LC_AK_Climatology_Reclass.tif"))
# 
# ##### Extracting climatology and topography values to Alaska SNOTEL points #####
# ### extracting topographic variables
# Snotel_AK$elevation <- terra::extract(AK_DEM, Snotel_AK)[,2]
# Snotel_AK$slope <- terra::extract(AK_Slope, Snotel_AK)[,2]
# Snotel_AK$aspect <- terra::extract(AK_Aspect, Snotel_AK)[,2]
# ### extracting climatologies
# ## prcpSum
# Snotel_AK$jan_prcpSum <- terra::extract(prcpSum_JanClim, Snotel_AK)[,2]
# Snotel_AK$feb_prcpSum <- terra::extract(prcpSum_FebClim, Snotel_AK)[,2]
# Snotel_AK$mar_prcpSum <- terra::extract(prcpSum_MarClim, Snotel_AK)[,2]
# Snotel_AK$apr_prcpSum <- terra::extract(prcpSum_AprClim, Snotel_AK)[,2]
# Snotel_AK$may_prcpSum <- terra::extract(prcpSum_MayClim, Snotel_AK)[,2]
# Snotel_AK$june_prcpSum <- terra::extract(prcpSum_JuneClim, Snotel_AK)[,2]
# Snotel_AK$july_prcpSum <- terra::extract(prcpSum_JulyClim, Snotel_AK)[,2]
# Snotel_AK$aug_prcpSum <- terra::extract(prcpSum_AugClim, Snotel_AK)[,2]
# Snotel_AK$sep_prcpSum <- terra::extract(prcpSum_SepClim, Snotel_AK)[,2]
# Snotel_AK$oct_prcpSum <- terra::extract(prcpSum_OctClim, Snotel_AK)[,2]
# Snotel_AK$nov_prcpSum <- terra::extract(prcpSum_NovClim, Snotel_AK)[,2]
# Snotel_AK$dec_prcpSum <- terra::extract(prcpSum_DecClim, Snotel_AK)[,2]
# ## prcpMean
# Snotel_AK$jan_prcpMean <- terra::extract(prcpMean_JanClim, Snotel_AK)[,2]
# Snotel_AK$feb_prcpMean <- terra::extract(prcpMean_FebClim, Snotel_AK)[,2]
# Snotel_AK$mar_prcpMean <- terra::extract(prcpMean_MarClim, Snotel_AK)[,2]
# Snotel_AK$apr_prcpMean <- terra::extract(prcpMean_AprClim, Snotel_AK)[,2]
# Snotel_AK$may_prcpMean <- terra::extract(prcpMean_MayClim, Snotel_AK)[,2]
# Snotel_AK$june_prcpMean <- terra::extract(prcpMean_JuneClim, Snotel_AK)[,2]
# Snotel_AK$july_prcpMean <- terra::extract(prcpMean_JulyClim, Snotel_AK)[,2]
# Snotel_AK$aug_prcpMean <- terra::extract(prcpMean_AugClim, Snotel_AK)[,2]
# Snotel_AK$sep_prcpMean <- terra::extract(prcpMean_SepClim, Snotel_AK)[,2]
# Snotel_AK$oct_prcpMean <- terra::extract(prcpMean_OctClim, Snotel_AK)[,2]
# Snotel_AK$nov_prcpMean <- terra::extract(prcpMean_NovClim, Snotel_AK)[,2]
# Snotel_AK$dec_prcpMean <- terra::extract(prcpMean_DecClim, Snotel_AK)[,2]
# 
# ## tmin
# Snotel_AK$jan_tmin <- terra::extract(tmin_JanClim, Snotel_AK)[,2]
# Snotel_AK$feb_tmin <- terra::extract(tmin_FebClim, Snotel_AK)[,2]
# Snotel_AK$mar_tmin <- terra::extract(tmin_MarClim, Snotel_AK)[,2]
# Snotel_AK$apr_tmin <- terra::extract(tmin_AprClim, Snotel_AK)[,2]
# Snotel_AK$may_tmin <- terra::extract(tmin_MayClim, Snotel_AK)[,2]
# Snotel_AK$june_tmin <- terra::extract(tmin_JuneClim, Snotel_AK)[,2]
# Snotel_AK$july_tmin <- terra::extract(tmin_JulyClim, Snotel_AK)[,2]
# Snotel_AK$aug_tmin <- terra::extract(tmin_AugClim, Snotel_AK)[,2]
# Snotel_AK$sep_tmin <- terra::extract(tmin_SepClim, Snotel_AK)[,2]
# Snotel_AK$oct_tmin <- terra::extract(tmin_OctClim, Snotel_AK)[,2]
# Snotel_AK$nov_tmin <- terra::extract(tmin_NovClim, Snotel_AK)[,2]
# Snotel_AK$dec_tmin <- terra::extract(tmin_DecClim, Snotel_AK)[,2]
# 
# ## tmean
# Snotel_AK$jan_tmean <- terra::extract(tmean_JanClim, Snotel_AK)[,2]
# Snotel_AK$feb_tmean <- terra::extract(tmean_FebClim, Snotel_AK)[,2]
# Snotel_AK$mar_tmean <- terra::extract(tmean_MarClim, Snotel_AK)[,2]
# Snotel_AK$apr_tmean <- terra::extract(tmean_AprClim, Snotel_AK)[,2]
# Snotel_AK$may_tmean <- terra::extract(tmean_MayClim, Snotel_AK)[,2]
# Snotel_AK$june_tmean <- terra::extract(tmean_JuneClim, Snotel_AK)[,2]
# Snotel_AK$july_tmean <- terra::extract(tmean_JulyClim, Snotel_AK)[,2]
# Snotel_AK$aug_tmean <- terra::extract(tmean_AugClim, Snotel_AK)[,2]
# Snotel_AK$sep_tmean <- terra::extract(tmean_SepClim, Snotel_AK)[,2]
# Snotel_AK$oct_tmean <- terra::extract(tmean_OctClim, Snotel_AK)[,2]
# Snotel_AK$nov_tmean <- terra::extract(tmean_NovClim, Snotel_AK)[,2]
# Snotel_AK$dec_tmean <- terra::extract(tmean_DecClim, Snotel_AK)[,2]
# 
# ## tmax
# Snotel_AK$jan_tmax <- terra::extract(tmax_JanClim, Snotel_AK)[,2]
# Snotel_AK$feb_tmax <- terra::extract(tmax_FebClim, Snotel_AK)[,2]
# Snotel_AK$mar_tmax <- terra::extract(tmax_MarClim, Snotel_AK)[,2]
# Snotel_AK$apr_tmax <- terra::extract(tmax_AprClim, Snotel_AK)[,2]
# Snotel_AK$may_tmax <- terra::extract(tmax_MayClim, Snotel_AK)[,2]
# Snotel_AK$june_tmax <- terra::extract(tmax_JuneClim, Snotel_AK)[,2]
# Snotel_AK$july_tmax <- terra::extract(tmax_JulyClim, Snotel_AK)[,2]
# Snotel_AK$aug_tmax <- terra::extract(tmax_AugClim, Snotel_AK)[,2]
# Snotel_AK$sep_tmax <- terra::extract(tmax_SepClim, Snotel_AK)[,2]
# Snotel_AK$oct_tmax <- terra::extract(tmax_OctClim, Snotel_AK)[,2]
# Snotel_AK$nov_tmax <- terra::extract(tmax_NovClim, Snotel_AK)[,2]
# Snotel_AK$dec_tmax <- terra::extract(tmax_DecClim, Snotel_AK)[,2]
# 
# ## srad
# Snotel_AK$jan_srad <- terra::extract(srad_JanClim, Snotel_AK)[,2]
# Snotel_AK$feb_srad <- terra::extract(srad_FebClim, Snotel_AK)[,2]
# Snotel_AK$mar_srad <- terra::extract(srad_MarClim, Snotel_AK)[,2]
# Snotel_AK$apr_srad <- terra::extract(srad_AprClim, Snotel_AK)[,2]
# Snotel_AK$may_srad <- terra::extract(srad_MayClim, Snotel_AK)[,2]
# Snotel_AK$june_srad <- terra::extract(srad_JuneClim, Snotel_AK)[,2]
# Snotel_AK$july_srad <- terra::extract(srad_JulyClim, Snotel_AK)[,2]
# Snotel_AK$aug_srad <- terra::extract(srad_AugClim, Snotel_AK)[,2]
# Snotel_AK$sep_srad <- terra::extract(srad_SepClim, Snotel_AK)[,2]
# Snotel_AK$oct_srad <- terra::extract(srad_OctClim, Snotel_AK)[,2]
# Snotel_AK$nov_srad <- terra::extract(srad_NovClim, Snotel_AK)[,2]
# Snotel_AK$dec_srad <- terra::extract(srad_DecClim, Snotel_AK)[,2]
# 
# 
# ### Extracting climatology/CDM aggregates
# ## Oct-May aggregates
# Snotel_AK$OctMay_prcpSumCDMSum <- terra::extract(prcpSum_OctMayClimSum_rast, Snotel_AK)[,2]
# Snotel_AK$OctMay_tmeanCDMSum <- terra::extract(tmean_OctMayClimSum_rast, Snotel_AK)[,2]
# Snotel_AK$OctMay_tminMean <- terra::extract(tmin_OctMayClimMean_rast, Snotel_AK)[,2]
# Snotel_AK$OctMay_tmaxMean <- terra::extract(tmax_OctMayClimMean_rast, Snotel_AK)[,2]
# Snotel_AK$OctMay_sradMean <- terra::extract(srad_OctMayClimMean_rast, Snotel_AK)[,2]
# ## Oct-Apr aggregates
# Snotel_AK$OctApr_prcpSumCDMSum <- terra::extract(prcpSum_OctAprClimSum_rast, Snotel_AK)[,2]
# Snotel_AK$OctApr_tmeanCDMSum <- terra::extract(tmean_OctAprClimSum_rast, Snotel_AK)[,2]
# Snotel_AK$OctApr_tminMean <- terra::extract(tmin_OctAprClimMean_rast, Snotel_AK)[,2]
# Snotel_AK$OctApr_tmaxMean <- terra::extract(tmax_OctAprClimMean_rast, Snotel_AK)[,2]
# Snotel_AK$OctApr_sradMean <- terra::extract(srad_OctAprClimMean_rast, Snotel_AK)[,2]
# ## Sep-May aggregates
# Snotel_AK$SepMay_prcpSumCDMSum <- terra::extract(prcpSum_SepMayClimSum_rast, Snotel_AK)[,2]
# Snotel_AK$SepMay_tmeanCDMSum <- terra::extract(tmean_SepMayClimSum_rast, Snotel_AK)[,2]
# Snotel_AK$SepMay_tminMean <- terra::extract(tmin_SepMayClimMean_rast, Snotel_AK)[,2]
# Snotel_AK$SepMay_tmaxMean <- terra::extract(tmax_SepMayClimMean_rast, Snotel_AK)[,2]
# Snotel_AK$SepMay_sradMean <- terra::extract(srad_SepMayClimMean_rast, Snotel_AK)[,2]
# ## Dec-Feb aggregates
# Snotel_AK$DecFeb_prcpSumCDMSum <- terra::extract(prcpSum_DecFebClimSum_rast, Snotel_AK)[,2]
# Snotel_AK$DecFeb_tmeanCDMSum <- terra::extract(tmean_DecFebClimSum_rast, Snotel_AK)[,2]
# Snotel_AK$DecFeb_tminMean <- terra::extract(tmin_DecFebClimMean_rast, Snotel_AK)[,2]
# Snotel_AK$DecFeb_tmaxMean <- terra::extract(tmax_DecFebClimMean_rast, Snotel_AK)[,2]
# Snotel_AK$DecFeb_sradMean <- terra::extract(srad_DecFebClimMean_rast, Snotel_AK)[,2]
# ### Extracting Climatological (modal) land cover value from Alaska land cover raster
# Snotel_AK$landcover <- terra::extract(LC_AK_Climatology_rast, Snotel_AK)[,2]
# Snotel_AK$landcover_triclass <- terra::extract(LC_AK_TriClass_Climatology_rast, Snotel_AK)[,2]
# 
# 
# ##### writing out SNOTEL Alaska SWE shapefile/gpkg with covariates #####
# # write_sf(Snotel_AK, here("Data", "SNOTEL", "Alaska", "GIS", "SHP", "SnotelAK_PeakSWEClimatology_Covars.shp"))
# write_sf(Snotel_AK, here("Data", "SNOTEL", "Alaska", "GIS", "GPKG", "Climatology", "SnotelAK_PeakSWEClimatology_Covars.gpkg"), append = FALSE)
# 
# AKClimendtime <- Sys.time()
# AKClimtime <- AKClimendtime - AKClimStartTime

##### CONUS ANNUAL STUFF #####
CONUSStarttime <- Sys.time()
##### Reading in CONUS variables and extracting them by year #####
Snotel_CONUS <- read_sf(here("Data", "SNOTEL", "CONUS", "GIS", "GPKG", "Annual", "SnotelData_CONUS_PeakSWE.gpkg"))


### creating zero points for individual years
ZeroPts_func <- function(year){
  if(year < 2003 | year == 2009 | year == 2018){
    SCF_NoSnow_rast <- rast(here("Data", "SCF", "Resampled", "SCF2003_2020NoSnow.tif"))
  } else{
    ### reading in SCF raster, filtering 
    SCF_rast <- rast(here("Data", "SCF", "Resampled", paste0("SCF", year, "US.tif")))
    
    SCF_max_rast <- app(SCF_rast, fun = "max")
    
    SCF_NoSnow_rast <- app(SCF_max_rast, fun = function(x){x[x > 0] <- NA; return(x)})
  }
  
  ### filtering Snotel points to just year zero pts are being created for so the correct number of zero pts are created
  temp_snotel <- Snotel_CONUS |>
    filter(WaterYear == as.numeric(year))
  ### Randomly sampling points within areas of no snow
  set.seed(802)
  Zero_points <- spatSample(SCF_NoSnow_rast, size = nrow(temp_snotel), method = "random", replace = FALSE, na.rm = TRUE, as.points = TRUE, values = FALSE)
  
  Zero_pts_sf <- st_as_sf(Zero_points)
  st_geometry(Zero_pts_sf) <- "geom"
  Zero_pts_sf <- Zero_pts_sf |>
    mutate(WaterYear = year,
           peak_swe = 0,
           site_id = -999,
           site_name = "zero values",
           description = "randomly generated points to be assigned 0 for peak SWE",
           state = "randomly generated points to be assigned 0 for peak SWE",
           start = min(temp_snotel$start),
           end = max(temp_snotel$end), .before = "geom")
  
  return(Zero_pts_sf)
}

### running function to create zero points
years <- 1993:2020
AllZeroPts <- pmclapply(years, ZeroPts_func, mc.cores = ifelse(length(years) > detectCores(), ceiling(detectCores() / 1.5), length(years)), mc.silent = FALSE)

### binding zero points to Snotel CONUS
for (x in 1:length(AllZeroPts)){
  if(x == 1){
    Snotel_CONUS <- rbind(Snotel_CONUS, AllZeroPts[[x]])
  }else{
    Snotel_CONUS <- rbind(Snotel_CONUS, AllZeroPts[[x]])
  }
}

### filtering out years so that appropriate values can be extracted
years <- 1993:2020
for (year in 1:length(years)){
  print(paste("CONUS:", years[year]))
  
  ### prcpSum
  prcpSum_Jan <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_JanSum.tif")))
  prcpSum_Feb <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_FebSum.tif")))
  prcpSum_Mar <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_MarSum.tif")))
  prcpSum_Apr <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_AprSum.tif")))
  prcpSum_May <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_MaySum.tif")))
  prcpSum_June <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_JuneSum.tif")))
  prcpSum_July <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_JulySum.tif")))
  prcpSum_Aug <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_AugSum.tif")))
  prcpSum_Sep <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_SepSum.tif")))
  prcpSum_Oct <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_OctSum.tif")))
  prcpSum_Nov <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_NovSum.tif")))
  prcpSum_Dec <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_DecSum.tif")))
  ### prcpMean
  prcpMean_Jan <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "JanMean.tif")))
  prcpMean_Feb <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_FebMean.tif")))
  prcpMean_Mar <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_MarMean.tif")))
  prcpMean_Apr <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_AprMean.tif")))
  prcpMean_May <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_MayMean.tif")))
  prcpMean_June <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_JuneMean.tif")))
  prcpMean_July <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_JulyMean.tif")))
  prcpMean_Aug <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_AugMean.tif")))
  prcpMean_Sep <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_SepMean.tif")))
  prcpMean_Oct <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_OctMean.tif")))
  prcpMean_Nov <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_NovMean.tif")))
  prcpMean_Dec <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_DecMean.tif")))
  ## tmin
  tmin_Jan <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmin", years[year], "JanMean.tif")))
  tmin_Feb <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmin", years[year], "_FebMean.tif")))
  tmin_Mar <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmin", years[year], "_MarMean.tif")))
  tmin_Apr <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmin", years[year], "_AprMean.tif")))
  tmin_May <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmin", years[year], "_MayMean.tif")))
  tmin_June <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmin", years[year], "_JuneMean.tif")))
  tmin_July <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmin", years[year], "_JulyMean.tif")))
  tmin_Aug <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmin", years[year], "_AugMean.tif")))
  tmin_Sep <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmin", years[year], "_SepMean.tif")))
  tmin_Oct <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmin", years[year], "_OctMean.tif")))
  tmin_Nov <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmin", years[year], "_NovMean.tif")))
  tmin_Dec <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmin", years[year], "_DecMean.tif")))
  
  ## tmean
  tmean_Jan <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmean", years[year], "JanMean.tif")))
  tmean_Feb <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmean", years[year], "_FebMean.tif")))
  tmean_Mar <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmean", years[year], "_MarMean.tif")))
  tmean_Apr <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmean", years[year], "_AprMean.tif")))
  tmean_May <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmean", years[year], "_MayMean.tif")))
  tmean_June <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmean", years[year], "_JuneMean.tif")))
  tmean_July <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmean", years[year], "_JulyMean.tif")))
  tmean_Aug <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmean", years[year], "_AugMean.tif")))
  tmean_Sep <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmean", years[year], "_SepMean.tif")))
  tmean_Oct <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmean", years[year], "_OctMean.tif")))
  tmean_Nov <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmean", years[year], "_NovMean.tif")))
  tmean_Dec <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmean", years[year], "_DecMean.tif")))
  
  ## tmax
  tmax_Jan <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmax", years[year], "JanMean.tif")))
  tmax_Feb <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmax", years[year], "_FebMean.tif")))
  tmax_Mar <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmax", years[year], "_MarMean.tif")))
  tmax_Apr <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmax", years[year], "_AprMean.tif")))
  tmax_May <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmax", years[year], "_MayMean.tif")))
  tmax_June <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmax", years[year], "_JuneMean.tif")))
  tmax_July <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmax", years[year], "_JulyMean.tif")))
  tmax_Aug <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmax", years[year], "_AugMean.tif")))
  tmax_Sep <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmax", years[year], "_SepMean.tif")))
  tmax_Oct <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmax", years[year], "_OctMean.tif")))
  tmax_Nov <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmax", years[year], "_NovMean.tif")))
  tmax_Dec <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmax", years[year], "_DecMean.tif")))
  
  ## srad
  srad_Jan <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("srad", years[year], "JanMean.tif")))
  srad_Feb <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("srad", years[year], "_FebMean.tif")))
  srad_Mar <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("srad", years[year], "_MarMean.tif")))
  srad_Apr <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("srad", years[year], "_AprMean.tif")))
  srad_May <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("srad", years[year], "_MayMean.tif")))
  srad_June <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("srad", years[year], "_JuneMean.tif")))
  srad_July <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("srad", years[year], "_JulyMean.tif")))
  srad_Aug <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("srad", years[year], "_AugMean.tif")))
  srad_Sep <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("srad", years[year], "_SepMean.tif")))
  srad_Oct <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("srad", years[year], "_OctMean.tif")))
  srad_Nov <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("srad", years[year], "_NovMean.tif")))
  srad_Dec <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", paste0("srad", years[year], "_DecMean.tif")))
  
  
  ##### reading in CONUS annual CDMs #####
  ### Oct-Apr Sums
  prcpSum_OctAprSum_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpSum_OctApr", years[year], "CDMSum.tif")))
  prcpMean_OctAprMean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpMean_OctApr", years[year], "CDMMean.tif")))
  tmean_OctAprSum_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmean_OctApr", years[year], "CDMSum.tif")))
  tmin_OctAprMean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmin_OctApr", years[year], "Mean.tif")))
  tmax_OctAprMean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmax_OctApr", years[year], "Mean.tif")))
  srad_OctAprMean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("srad_OctApr", years[year], "Mean.tif")))
  ### Oct-May Sums
  prcpSum_OctMaySum_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpSum_OctMay", years[year], "CDMSum.tif")))
  prcpMean_OctMayMean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpMean_OctMay", years[year], "CDMMean.tif")))
  tmean_OctMaySum_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmean_OctMay", years[year], "CDMSum.tif")))
  tmin_OctMayMean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmin_OctMay", years[year], "Mean.tif")))
  tmax_OctMayMean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmax_OctMay", years[year], "Mean.tif")))
  srad_OctMayMean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("srad_OctMay", years[year], "Mean.tif")))
  ### Sep-May Sums
  prcpSum_SepMaySum_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpSum_SepMay", years[year], "CDMSum.tif")))
  prcpMean_SepMayMean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpMean_SepMay", years[year], "CDMMean.tif")))
  tmean_SepMaySum_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmean_SepMay", years[year], "CDMSum.tif")))
  tmin_SepMayMean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmin_SepMay", years[year], "Mean.tif")))
  tmax_SepMayMean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmax_SepMay", years[year], "Mean.tif")))
  srad_SepMayMean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("srad_SepMay", years[year], "Mean.tif")))
  ### Dec-Feb Sums
  prcpSum_DecFebSum_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpSum_DecFeb", years[year], "CDMSum.tif")))
  prcpMean_DecFebMean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpMean_DecFeb", years[year], "CDMMean.tif")))
  tmean_DecFebSum_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmean_DecFeb", years[year], "CDMSum.tif")))
  tmin_DecFebMean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmin_DecFeb", years[year], "Mean.tif")))
  tmax_DecFebMean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmax_DecFeb", years[year], "Mean.tif")))
  srad_DecFebMean_rast <- rast(here("Data", "Daymet", "CONUS", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("srad_DecFeb", years[year], "Mean.tif")))
  
  
  ##### Extracting weather and topography values to Annual CONUS SNOTEL points #####
  ### extracting topographic variables
  Snotel_CONUS$elevation <- terra::extract(CONUS_DEM, Snotel_CONUS)[,2]
  Snotel_CONUS$slope <- terra::extract(CONUS_Slope, Snotel_CONUS)[,2]
  Snotel_CONUS$aspect <- terra::extract(CONUS_Aspect, Snotel_CONUS)[,2]
  
  if (year == 1){
    temp_snotel <- Snotel_CONUS |>
      filter(WaterYear == years[year])
    temp_lc_rast <- lc_CONUS_rast_list[[year]]
    temp_lc_TriClass_rast <- LC_CONUS_TriClass_rast_list[[year]]
    ### extracting climatologies
    ## prcpSum
    temp_snotel$jan_prcpSum <- terra::extract(prcpSum_Jan, temp_snotel)[,2]
    temp_snotel$feb_prcpSum <- terra::extract(prcpSum_Feb, temp_snotel)[,2]
    temp_snotel$mar_prcpSum <- terra::extract(prcpSum_Mar, temp_snotel)[,2]
    temp_snotel$apr_prcpSum <- terra::extract(prcpSum_Apr, temp_snotel)[,2]
    temp_snotel$may_prcpSum <- terra::extract(prcpSum_May, temp_snotel)[,2]
    temp_snotel$june_prcpSum <- terra::extract(prcpSum_June, temp_snotel)[,2]
    temp_snotel$july_prcpSum <- terra::extract(prcpSum_July, temp_snotel)[,2]
    temp_snotel$aug_prcpSum <- terra::extract(prcpSum_Aug, temp_snotel)[,2]
    temp_snotel$sep_prcpSum <- terra::extract(prcpSum_Sep, temp_snotel)[,2]
    temp_snotel$oct_prcpSum <- terra::extract(prcpSum_Oct, temp_snotel)[,2]
    temp_snotel$nov_prcpSum <- terra::extract(prcpSum_Nov, temp_snotel)[,2]
    temp_snotel$dec_prcpSum <- terra::extract(prcpSum_Dec, temp_snotel)[,2]
    ## prcpMean
    temp_snotel$jan_prcpMean <- terra::extract(prcpMean_Jan, temp_snotel)[,2]
    temp_snotel$feb_prcpMean <- terra::extract(prcpMean_Feb, temp_snotel)[,2]
    temp_snotel$mar_prcpMean <- terra::extract(prcpMean_Mar, temp_snotel)[,2]
    temp_snotel$apr_prcpMean <- terra::extract(prcpMean_Apr, temp_snotel)[,2]
    temp_snotel$may_prcpMean <- terra::extract(prcpMean_May, temp_snotel)[,2]
    temp_snotel$june_prcpMean <- terra::extract(prcpMean_June, temp_snotel)[,2]
    temp_snotel$july_prcpMean <- terra::extract(prcpMean_July, temp_snotel)[,2]
    temp_snotel$aug_prcpMean <- terra::extract(prcpMean_Aug, temp_snotel)[,2]
    temp_snotel$sep_prcpMean <- terra::extract(prcpMean_Sep, temp_snotel)[,2]
    temp_snotel$oct_prcpMean <- terra::extract(prcpMean_Oct, temp_snotel)[,2]
    temp_snotel$nov_prcpMean <- terra::extract(prcpMean_Nov, temp_snotel)[,2]
    temp_snotel$dec_prcpMean <- terra::extract(prcpMean_Dec, temp_snotel)[,2]
    
    ## tmin
    temp_snotel$jan_tmin <- terra::extract(tmin_Jan, temp_snotel)[,2]
    temp_snotel$feb_tmin <- terra::extract(tmin_Feb, temp_snotel)[,2]
    temp_snotel$mar_tmin <- terra::extract(tmin_Mar, temp_snotel)[,2]
    temp_snotel$apr_tmin <- terra::extract(tmin_Apr, temp_snotel)[,2]
    temp_snotel$may_tmin <- terra::extract(tmin_May, temp_snotel)[,2]
    temp_snotel$june_tmin <- terra::extract(tmin_June, temp_snotel)[,2]
    temp_snotel$july_tmin <- terra::extract(tmin_July, temp_snotel)[,2]
    temp_snotel$aug_tmin <- terra::extract(tmin_Aug, temp_snotel)[,2]
    temp_snotel$sep_tmin <- terra::extract(tmin_Sep, temp_snotel)[,2]
    temp_snotel$oct_tmin <- terra::extract(tmin_Oct, temp_snotel)[,2]
    temp_snotel$nov_tmin <- terra::extract(tmin_Nov, temp_snotel)[,2]
    temp_snotel$dec_tmin <- terra::extract(tmin_Dec, temp_snotel)[,2]
    
    ## tmean
    temp_snotel$jan_tmean <- terra::extract(tmean_Jan, temp_snotel)[,2]
    temp_snotel$feb_tmean <- terra::extract(tmean_Feb, temp_snotel)[,2]
    temp_snotel$mar_tmean <- terra::extract(tmean_Mar, temp_snotel)[,2]
    temp_snotel$apr_tmean <- terra::extract(tmean_Apr, temp_snotel)[,2]
    temp_snotel$may_tmean <- terra::extract(tmean_May, temp_snotel)[,2]
    temp_snotel$june_tmean <- terra::extract(tmean_June, temp_snotel)[,2]
    temp_snotel$july_tmean <- terra::extract(tmean_July, temp_snotel)[,2]
    temp_snotel$aug_tmean <- terra::extract(tmean_Aug, temp_snotel)[,2]
    temp_snotel$sep_tmean <- terra::extract(tmean_Sep, temp_snotel)[,2]
    temp_snotel$oct_tmean <- terra::extract(tmean_Oct, temp_snotel)[,2]
    temp_snotel$nov_tmean <- terra::extract(tmean_Nov, temp_snotel)[,2]
    temp_snotel$dec_tmean <- terra::extract(tmean_Dec, temp_snotel)[,2]
    
    ## tmax
    temp_snotel$jan_tmax <- terra::extract(tmax_Jan, temp_snotel)[,2]
    temp_snotel$feb_tmax <- terra::extract(tmax_Feb, temp_snotel)[,2]
    temp_snotel$mar_tmax <- terra::extract(tmax_Mar, temp_snotel)[,2]
    temp_snotel$apr_tmax <- terra::extract(tmax_Apr, temp_snotel)[,2]
    temp_snotel$may_tmax <- terra::extract(tmax_May, temp_snotel)[,2]
    temp_snotel$june_tmax <- terra::extract(tmax_June, temp_snotel)[,2]
    temp_snotel$july_tmax <- terra::extract(tmax_July, temp_snotel)[,2]
    temp_snotel$aug_tmax <- terra::extract(tmax_Aug, temp_snotel)[,2]
    temp_snotel$sep_tmax <- terra::extract(tmax_Sep, temp_snotel)[,2]
    temp_snotel$oct_tmax <- terra::extract(tmax_Oct, temp_snotel)[,2]
    temp_snotel$nov_tmax <- terra::extract(tmax_Nov, temp_snotel)[,2]
    temp_snotel$dec_tmax <- terra::extract(tmax_Dec, temp_snotel)[,2]
    
    ## srad
    temp_snotel$jan_srad <- terra::extract(srad_Jan, temp_snotel)[,2]
    temp_snotel$feb_srad <- terra::extract(srad_Feb, temp_snotel)[,2]
    temp_snotel$mar_srad <- terra::extract(srad_Mar, temp_snotel)[,2]
    temp_snotel$apr_srad <- terra::extract(srad_Apr, temp_snotel)[,2]
    temp_snotel$may_srad <- terra::extract(srad_May, temp_snotel)[,2]
    temp_snotel$june_srad <- terra::extract(srad_June, temp_snotel)[,2]
    temp_snotel$july_srad <- terra::extract(srad_July, temp_snotel)[,2]
    temp_snotel$aug_srad <- terra::extract(srad_Aug, temp_snotel)[,2]
    temp_snotel$sep_srad <- terra::extract(srad_Sep, temp_snotel)[,2]
    temp_snotel$oct_srad <- terra::extract(srad_Oct, temp_snotel)[,2]
    temp_snotel$nov_srad <- terra::extract(srad_Nov, temp_snotel)[,2]
    temp_snotel$dec_srad <- terra::extract(srad_Dec, temp_snotel)[,2]
    
    
    ### Extracting climatology/CDM aggregates
    ## Oct-May aggregates
    temp_snotel$OctMay_prcpSumCDMSum <- terra::extract(prcpSum_OctMaySum_rast, temp_snotel)[,2]
    temp_snotel$OctMay_tmeanCDMSum <- terra::extract(tmean_OctMaySum_rast, temp_snotel)[,2]
    temp_snotel$OctMay_tminMean <- terra::extract(tmin_OctMayMean_rast, temp_snotel)[,2]
    temp_snotel$OctMay_tmaxMean <- terra::extract(tmax_OctMayMean_rast, temp_snotel)[,2]
    temp_snotel$OctMay_sradMean <- terra::extract(srad_OctMayMean_rast, temp_snotel)[,2]
    ## Oct-Apr aggregates
    temp_snotel$OctApr_prcpSumCDMSum <- terra::extract(prcpSum_OctAprSum_rast, temp_snotel)[,2]
    temp_snotel$OctApr_tmeanCDMSum <- terra::extract(tmean_OctAprSum_rast, temp_snotel)[,2]
    temp_snotel$OctApr_tminMean <- terra::extract(tmin_OctAprMean_rast, temp_snotel)[,2]
    temp_snotel$OctApr_tmaxMean <- terra::extract(tmax_OctAprMean_rast, temp_snotel)[,2]
    temp_snotel$OctApr_sradMean <- terra::extract(srad_OctAprMean_rast, temp_snotel)[,2]
    ## Sep-May aggregates
    temp_snotel$SepMay_prcpSumCDMSum <- terra::extract(prcpSum_SepMaySum_rast, temp_snotel)[,2]
    temp_snotel$SepMay_tmeanCDMSum <- terra::extract(tmean_SepMaySum_rast, temp_snotel)[,2]
    temp_snotel$SepMay_tminMean <- terra::extract(tmin_SepMayMean_rast, temp_snotel)[,2]
    temp_snotel$SepMay_tmaxMean <- terra::extract(tmax_SepMayMean_rast, temp_snotel)[,2]
    temp_snotel$SepMay_sradMean <- terra::extract(srad_SepMayMean_rast, temp_snotel)[,2]
    ## Dec-Feb aggregates
    temp_snotel$DecFeb_prcpSumCDMSum <- terra::extract(prcpSum_DecFebSum_rast, temp_snotel)[,2]
    temp_snotel$DecFeb_tmeanCDMSum <- terra::extract(tmean_DecFebSum_rast, temp_snotel)[,2]
    temp_snotel$DecFeb_tminMean <- terra::extract(tmin_DecFebMean_rast, temp_snotel)[,2]
    temp_snotel$DecFeb_tmaxMean <- terra::extract(tmax_DecFebMean_rast, temp_snotel)[,2]
    temp_snotel$DecFeb_sradMean <- terra::extract(srad_DecFebMean_rast, temp_snotel)[,2]
    temp_snotel$landcover = terra::extract(temp_lc_rast, temp_snotel)[,2]
    temp_snotel$landcover_triclass = terra::extract(temp_lc_TriClass_rast, temp_snotel)[,2]
  } else{
    temp_snotel2 <- Snotel_CONUS |>
      filter(WaterYear == years[year])
    temp_lc_rast <- lc_CONUS_rast_list[[year]]
    temp_lc_TriClass_rast <- LC_CONUS_TriClass_rast_list[[year]]
    ### extracting climatologies
    ## prcpSum
    temp_snotel2$jan_prcpSum <- terra::extract(prcpSum_Jan, temp_snotel2)[,2]
    temp_snotel2$feb_prcpSum <- terra::extract(prcpSum_Feb, temp_snotel2)[,2]
    temp_snotel2$mar_prcpSum <- terra::extract(prcpSum_Mar, temp_snotel2)[,2]
    temp_snotel2$apr_prcpSum <- terra::extract(prcpSum_Apr, temp_snotel2)[,2]
    temp_snotel2$may_prcpSum <- terra::extract(prcpSum_May, temp_snotel2)[,2]
    temp_snotel2$june_prcpSum <- terra::extract(prcpSum_June, temp_snotel2)[,2]
    temp_snotel2$july_prcpSum <- terra::extract(prcpSum_July, temp_snotel2)[,2]
    temp_snotel2$aug_prcpSum <- terra::extract(prcpSum_Aug, temp_snotel2)[,2]
    temp_snotel2$sep_prcpSum <- terra::extract(prcpSum_Sep, temp_snotel2)[,2]
    temp_snotel2$oct_prcpSum <- terra::extract(prcpSum_Oct, temp_snotel2)[,2]
    temp_snotel2$nov_prcpSum <- terra::extract(prcpSum_Nov, temp_snotel2)[,2]
    temp_snotel2$dec_prcpSum <- terra::extract(prcpSum_Dec, temp_snotel2)[,2]
    ## prcpMean
    temp_snotel2$jan_prcpMean <- terra::extract(prcpMean_Jan, temp_snotel2)[,2]
    temp_snotel2$feb_prcpMean <- terra::extract(prcpMean_Feb, temp_snotel2)[,2]
    temp_snotel2$mar_prcpMean <- terra::extract(prcpMean_Mar, temp_snotel2)[,2]
    temp_snotel2$apr_prcpMean <- terra::extract(prcpMean_Apr, temp_snotel2)[,2]
    temp_snotel2$may_prcpMean <- terra::extract(prcpMean_May, temp_snotel2)[,2]
    temp_snotel2$june_prcpMean <- terra::extract(prcpMean_June, temp_snotel2)[,2]
    temp_snotel2$july_prcpMean <- terra::extract(prcpMean_July, temp_snotel2)[,2]
    temp_snotel2$aug_prcpMean <- terra::extract(prcpMean_Aug, temp_snotel2)[,2]
    temp_snotel2$sep_prcpMean <- terra::extract(prcpMean_Sep, temp_snotel2)[,2]
    temp_snotel2$oct_prcpMean <- terra::extract(prcpMean_Oct, temp_snotel2)[,2]
    temp_snotel2$nov_prcpMean <- terra::extract(prcpMean_Nov, temp_snotel2)[,2]
    temp_snotel2$dec_prcpMean <- terra::extract(prcpMean_Dec, temp_snotel2)[,2]
    
    ## tmin
    temp_snotel2$jan_tmin <- terra::extract(tmin_Jan, temp_snotel2)[,2]
    temp_snotel2$feb_tmin <- terra::extract(tmin_Feb, temp_snotel2)[,2]
    temp_snotel2$mar_tmin <- terra::extract(tmin_Mar, temp_snotel2)[,2]
    temp_snotel2$apr_tmin <- terra::extract(tmin_Apr, temp_snotel2)[,2]
    temp_snotel2$may_tmin <- terra::extract(tmin_May, temp_snotel2)[,2]
    temp_snotel2$june_tmin <- terra::extract(tmin_June, temp_snotel2)[,2]
    temp_snotel2$july_tmin <- terra::extract(tmin_July, temp_snotel2)[,2]
    temp_snotel2$aug_tmin <- terra::extract(tmin_Aug, temp_snotel2)[,2]
    temp_snotel2$sep_tmin <- terra::extract(tmin_Sep, temp_snotel2)[,2]
    temp_snotel2$oct_tmin <- terra::extract(tmin_Oct, temp_snotel2)[,2]
    temp_snotel2$nov_tmin <- terra::extract(tmin_Nov, temp_snotel2)[,2]
    temp_snotel2$dec_tmin <- terra::extract(tmin_Dec, temp_snotel2)[,2]
    
    ## tmean
    temp_snotel2$jan_tmean <- terra::extract(tmean_Jan, temp_snotel2)[,2]
    temp_snotel2$feb_tmean <- terra::extract(tmean_Feb, temp_snotel2)[,2]
    temp_snotel2$mar_tmean <- terra::extract(tmean_Mar, temp_snotel2)[,2]
    temp_snotel2$apr_tmean <- terra::extract(tmean_Apr, temp_snotel2)[,2]
    temp_snotel2$may_tmean <- terra::extract(tmean_May, temp_snotel2)[,2]
    temp_snotel2$june_tmean <- terra::extract(tmean_June, temp_snotel2)[,2]
    temp_snotel2$july_tmean <- terra::extract(tmean_July, temp_snotel2)[,2]
    temp_snotel2$aug_tmean <- terra::extract(tmean_Aug, temp_snotel2)[,2]
    temp_snotel2$sep_tmean <- terra::extract(tmean_Sep, temp_snotel2)[,2]
    temp_snotel2$oct_tmean <- terra::extract(tmean_Oct, temp_snotel2)[,2]
    temp_snotel2$nov_tmean <- terra::extract(tmean_Nov, temp_snotel2)[,2]
    temp_snotel2$dec_tmean <- terra::extract(tmean_Dec, temp_snotel2)[,2]
    
    ## tmax
    temp_snotel2$jan_tmax <- terra::extract(tmax_Jan, temp_snotel2)[,2]
    temp_snotel2$feb_tmax <- terra::extract(tmax_Feb, temp_snotel2)[,2]
    temp_snotel2$mar_tmax <- terra::extract(tmax_Mar, temp_snotel2)[,2]
    temp_snotel2$apr_tmax <- terra::extract(tmax_Apr, temp_snotel2)[,2]
    temp_snotel2$may_tmax <- terra::extract(tmax_May, temp_snotel2)[,2]
    temp_snotel2$june_tmax <- terra::extract(tmax_June, temp_snotel2)[,2]
    temp_snotel2$july_tmax <- terra::extract(tmax_July, temp_snotel2)[,2]
    temp_snotel2$aug_tmax <- terra::extract(tmax_Aug, temp_snotel2)[,2]
    temp_snotel2$sep_tmax <- terra::extract(tmax_Sep, temp_snotel2)[,2]
    temp_snotel2$oct_tmax <- terra::extract(tmax_Oct, temp_snotel2)[,2]
    temp_snotel2$nov_tmax <- terra::extract(tmax_Nov, temp_snotel2)[,2]
    temp_snotel2$dec_tmax <- terra::extract(tmax_Dec, temp_snotel2)[,2]
    
    ## srad
    temp_snotel2$jan_srad <- terra::extract(srad_Jan, temp_snotel2)[,2]
    temp_snotel2$feb_srad <- terra::extract(srad_Feb, temp_snotel2)[,2]
    temp_snotel2$mar_srad <- terra::extract(srad_Mar, temp_snotel2)[,2]
    temp_snotel2$apr_srad <- terra::extract(srad_Apr, temp_snotel2)[,2]
    temp_snotel2$may_srad <- terra::extract(srad_May, temp_snotel2)[,2]
    temp_snotel2$june_srad <- terra::extract(srad_June, temp_snotel2)[,2]
    temp_snotel2$july_srad <- terra::extract(srad_July, temp_snotel2)[,2]
    temp_snotel2$aug_srad <- terra::extract(srad_Aug, temp_snotel2)[,2]
    temp_snotel2$sep_srad <- terra::extract(srad_Sep, temp_snotel2)[,2]
    temp_snotel2$oct_srad <- terra::extract(srad_Oct, temp_snotel2)[,2]
    temp_snotel2$nov_srad <- terra::extract(srad_Nov, temp_snotel2)[,2]
    temp_snotel2$dec_srad <- terra::extract(srad_Dec, temp_snotel2)[,2]
    
    
    ### Extracting climatology/CDM aggregates
    ## Oct-May aggregates
    temp_snotel2$OctMay_prcpSumCDMSum <- terra::extract(prcpSum_OctMaySum_rast, temp_snotel2)[,2]
    temp_snotel2$OctMay_tmeanCDMSum <- terra::extract(tmean_OctMaySum_rast, temp_snotel2)[,2]
    temp_snotel2$OctMay_tminMean <- terra::extract(tmin_OctMayMean_rast, temp_snotel2)[,2]
    temp_snotel2$OctMay_tmaxMean <- terra::extract(tmax_OctMayMean_rast, temp_snotel2)[,2]
    temp_snotel2$OctMay_sradMean <- terra::extract(srad_OctMayMean_rast, temp_snotel2)[,2]
    ## Oct-Apr aggregates
    temp_snotel2$OctApr_prcpSumCDMSum <- terra::extract(prcpSum_OctAprSum_rast, temp_snotel2)[,2]
    temp_snotel2$OctApr_tmeanCDMSum <- terra::extract(tmean_OctAprSum_rast, temp_snotel2)[,2]
    temp_snotel2$OctApr_tminMean <- terra::extract(tmin_OctAprMean_rast, temp_snotel2)[,2]
    temp_snotel2$OctApr_tmaxMean <- terra::extract(tmax_OctAprMean_rast, temp_snotel2)[,2]
    temp_snotel2$OctApr_sradMean <- terra::extract(srad_OctAprMean_rast, temp_snotel2)[,2]
    ## Sep-May aggregates
    temp_snotel2$SepMay_prcpSumCDMSum <- terra::extract(prcpSum_SepMaySum_rast, temp_snotel2)[,2]
    temp_snotel2$SepMay_tmeanCDMSum <- terra::extract(tmean_SepMaySum_rast, temp_snotel2)[,2]
    temp_snotel2$SepMay_tminMean <- terra::extract(tmin_SepMayMean_rast, temp_snotel2)[,2]
    temp_snotel2$SepMay_tmaxMean <- terra::extract(tmax_SepMayMean_rast, temp_snotel2)[,2]
    temp_snotel2$SepMay_sradMean <- terra::extract(srad_SepMayMean_rast, temp_snotel2)[,2]
    ## Dec-Feb aggregates
    temp_snotel2$DecFeb_prcpSumCDMSum <- terra::extract(prcpSum_DecFebSum_rast, temp_snotel2)[,2]
    temp_snotel2$DecFeb_tmeanCDMSum <- terra::extract(tmean_DecFebSum_rast, temp_snotel2)[,2]
    temp_snotel2$DecFeb_tminMean <- terra::extract(tmin_DecFebMean_rast, temp_snotel2)[,2]
    temp_snotel2$DecFeb_tmaxMean <- terra::extract(tmax_DecFebMean_rast, temp_snotel2)[,2]
    temp_snotel2$DecFeb_sradMean <- terra::extract(srad_DecFebMean_rast, temp_snotel2)[,2]
    temp_snotel2$landcover = terra::extract(temp_lc_rast, temp_snotel2)[,2]
    temp_snotel2$landcover_triclass = terra::extract(temp_lc_TriClass_rast, temp_snotel2)[,2]
    temp_snotel <- rbind(temp_snotel, temp_snotel2)
  }
}
### storing combined years with landcover added back in the Snotel_CONUS variable
Snotel_CONUS <- temp_snotel

##### writing out annual CONUS gpkg with covariates #####
# write_sf(Snotel_CONUS, here("Data", "SNOTEL", "CONUS", "GIS", "SHP", "SnotelCONUS_PeakSWE_Covars.shp"))
# write_sf(Snotel_CONUS, here("Data", "SNOTEL", "CONUS", "GIS", "GPKG", "Annual", "SnotelCONUS_PeakSWE_Covars.gpkg"), append = FALSE)
write_sf(Snotel_CONUS, here("Data", "SNOTEL", "CONUS", "GIS", "GPKG", "Annual", "SnotelCONUSAnnualZeroPts_PeakSWE_Covars.gpkg"), append = FALSE)

CONUSEndtime <- Sys.time()
CONUStime <- CONUSEndtime - CONUSStarttime

##### ALASKA Annual stuff #####
AKStarttime <- Sys.time()

##### Reading in Alaska variables and extracting them by year #####
Snotel_AK <- read_sf(here("Data", "SNOTEL", "Alaska", "GIS", "GPKG", "Annual", "SnotelData_AK_PeakSWE.gpkg"))
# ### filtering out years so that appropriate values can be extracted
years <- 1993:2020
for (year in 1:length(years)){
  print(paste("Alaska:", years[year]))

  ### prcpSum
  prcpSum_Jan <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "JanSum.tif")))
  prcpSum_Feb <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_FebSum.tif")))
  prcpSum_Mar <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_MarSum.tif")))
  prcpSum_Apr <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_AprSum.tif")))
  prcpSum_May <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_MaySum.tif")))
  prcpSum_June <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_JuneSum.tif")))
  prcpSum_July <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_JulySum.tif")))
  prcpSum_Aug <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_AugSum.tif")))
  prcpSum_Sep <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_SepSum.tif")))
  prcpSum_Oct <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_OctSum.tif")))
  prcpSum_Nov <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_NovSum.tif")))
  prcpSum_Dec <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_DecSum.tif")))
  ### prcpMean
  prcpMean_Jan <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "JanMean.tif")))
  prcpMean_Feb <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_FebMean.tif")))
  prcpMean_Mar <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_MarMean.tif")))
  prcpMean_Apr <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_AprMean.tif")))
  prcpMean_May <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_MayMean.tif")))
  prcpMean_June <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_JuneMean.tif")))
  prcpMean_July <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_JulyMean.tif")))
  prcpMean_Aug <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_AugMean.tif")))
  prcpMean_Sep <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_SepMean.tif")))
  prcpMean_Oct <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_OctMean.tif")))
  prcpMean_Nov <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_NovMean.tif")))
  prcpMean_Dec <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("prcp", years[year], "_DecMean.tif")))
  ## tmin
  tmin_Jan <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmin", years[year], "JanMean.tif")))
  tmin_Feb <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmin", years[year], "_FebMean.tif")))
  tmin_Mar <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmin", years[year], "_MarMean.tif")))
  tmin_Apr <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmin", years[year], "_AprMean.tif")))
  tmin_May <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmin", years[year], "_MayMean.tif")))
  tmin_June <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmin", years[year], "_JuneMean.tif")))
  tmin_July <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmin", years[year], "_JulyMean.tif")))
  tmin_Aug <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmin", years[year], "_AugMean.tif")))
  tmin_Sep <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmin", years[year], "_SepMean.tif")))
  tmin_Oct <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmin", years[year], "_OctMean.tif")))
  tmin_Nov <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmin", years[year], "_NovMean.tif")))
  tmin_Dec <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmin", years[year], "_DecMean.tif")))

  ## tmean
  tmean_Jan <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmean", years[year], "JanMean.tif")))
  tmean_Feb <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmean", years[year], "_FebMean.tif")))
  tmean_Mar <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmean", years[year], "_MarMean.tif")))
  tmean_Apr <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmean", years[year], "_AprMean.tif")))
  tmean_May <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmean", years[year], "_MayMean.tif")))
  tmean_June <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmean", years[year], "_JuneMean.tif")))
  tmean_July <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmean", years[year], "_JulyMean.tif")))
  tmean_Aug <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmean", years[year], "_AugMean.tif")))
  tmean_Sep <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmean", years[year], "_SepMean.tif")))
  tmean_Oct <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmean", years[year], "_OctMean.tif")))
  tmean_Nov <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmean", years[year], "_NovMean.tif")))
  tmean_Dec <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmean", years[year], "_DecMean.tif")))

  ## tmax
  tmax_Jan <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmax", years[year], "JanMean.tif")))
  tmax_Feb <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmax", years[year], "_FebMean.tif")))
  tmax_Mar <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmax", years[year], "_MarMean.tif")))
  tmax_Apr <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmax", years[year], "_AprMean.tif")))
  tmax_May <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmax", years[year], "_MayMean.tif")))
  tmax_June <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmax", years[year], "_JuneMean.tif")))
  tmax_July <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmax", years[year], "_JulyMean.tif")))
  tmax_Aug <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmax", years[year], "_AugMean.tif")))
  tmax_Sep <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmax", years[year], "_SepMean.tif")))
  tmax_Oct <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmax", years[year], "_OctMean.tif")))
  tmax_Nov <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmax", years[year], "_NovMean.tif")))
  tmax_Dec <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("tmax", years[year], "_DecMean.tif")))

  ## srad
  srad_Jan <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("srad", years[year], "JanMean.tif")))
  srad_Feb <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("srad", years[year], "_FebMean.tif")))
  srad_Mar <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("srad", years[year], "_MarMean.tif")))
  srad_Apr <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("srad", years[year], "_AprMean.tif")))
  srad_May <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("srad", years[year], "_MayMean.tif")))
  srad_June <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("srad", years[year], "_JuneMean.tif")))
  srad_July <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("srad", years[year], "_JulyMean.tif")))
  srad_Aug <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("srad", years[year], "_AugMean.tif")))
  srad_Sep <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("srad", years[year], "_SepMean.tif")))
  srad_Oct <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("srad", years[year], "_OctMean.tif")))
  srad_Nov <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("srad", years[year], "_NovMean.tif")))
  srad_Dec <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", paste0("srad", years[year], "_DecMean.tif")))


  ##### reading in Alaska annual CDMs #####
  ### Oct-Apr Sums
  prcpSum_OctAprSum_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpSum_OctApr", years[year], "CDMSum.tif")))
  prcpMean_OctAprMean_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpMean_OctApr", years[year], "CDMMean.tif")))
  tmean_OctAprSum_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmean_OctApr", years[year], "CDMSum.tif")))
  tmin_OctAprMean_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmin_OctApr", years[year], "Mean.tif")))
  tmax_OctAprMean_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmax_OctApr", years[year], "Mean.tif")))
  srad_OctAprMean_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("srad_OctApr", years[year], "Mean.tif")))
  ### Oct-May Sums
  prcpSum_OctMaySum_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpSum_OctMay", years[year], "CDMSum.tif")))
  prcpMean_OctMayMean_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpMean_OctMay", years[year], "CDMMean.tif")))
  tmean_OctMaySum_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmean_OctMay", years[year], "CDMSum.tif")))
  tmin_OctMayMean_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmin_OctMay", years[year], "Mean.tif")))
  tmax_OctMayMean_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmax_OctMay", years[year], "Mean.tif")))
  srad_OctMayMean_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("srad_OctMay", years[year], "Mean.tif")))
  ### Sep-May Sums
  prcpSum_SepMaySum_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpSum_SepMay", years[year], "CDMSum.tif")))
  prcpMean_SepMayMean_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpMean_SepMay", years[year], "CDMMean.tif")))
  tmean_SepMaySum_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmean_SepMay", years[year], "CDMSum.tif")))
  tmin_SepMayMean_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmin_SepMay", years[year], "Mean.tif")))
  tmax_SepMayMean_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmax_SepMay", years[year], "Mean.tif")))
  srad_SepMayMean_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("srad_SepMay", years[year], "Mean.tif")))
  ### Dec-Feb Sums
  prcpSum_DecFebSum_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpSum_DecFeb", years[year], "CDMSum.tif")))
  prcpMean_DecFebMean_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("prcpMean_DecFeb", years[year], "CDMMean.tif")))
  tmean_DecFebSum_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmean_DecFeb", years[year], "CDMSum.tif")))
  tmin_DecFebMean_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmin_DecFeb", years[year], "Mean.tif")))
  tmax_DecFebMean_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("tmax_DecFeb", years[year], "Mean.tif")))
  srad_DecFebMean_rast <- rast(here("Data", "Daymet", "Alaska", "Mosaics", "TIFFs", "ClippedResampled", "CDM", "CDM_Aggregates", paste0("srad_DecFeb", years[year], "Mean.tif")))


  ##### Extracting weather and topography values to Annual Alaska SNOTEL points #####
  ### extracting topographic variables
  Snotel_AK$elevation <- terra::extract(AK_DEM, Snotel_AK)[,2]
  Snotel_AK$slope <- terra::extract(AK_Slope, Snotel_AK)[,2]
  Snotel_AK$aspect <- terra::extract(AK_Aspect, Snotel_AK)[,2]

  if (year == 1){
    temp_snotel <- Snotel_AK |>
      filter(WaterYear == years[year])
    temp_lc_rast <- lc_AK_rast_list[[year]]
    temp_lc_TriClass_rast <- LC_AK_TriClass_rast_list[[year]]
    ### extracting climatologies
    ## prcpSum
    temp_snotel$jan_prcpSum <- terra::extract(prcpSum_Jan, temp_snotel)[,2]
    temp_snotel$feb_prcpSum <- terra::extract(prcpSum_Feb, temp_snotel)[,2]
    temp_snotel$mar_prcpSum <- terra::extract(prcpSum_Mar, temp_snotel)[,2]
    temp_snotel$apr_prcpSum <- terra::extract(prcpSum_Apr, temp_snotel)[,2]
    temp_snotel$may_prcpSum <- terra::extract(prcpSum_May, temp_snotel)[,2]
    temp_snotel$june_prcpSum <- terra::extract(prcpSum_June, temp_snotel)[,2]
    temp_snotel$july_prcpSum <- terra::extract(prcpSum_July, temp_snotel)[,2]
    temp_snotel$aug_prcpSum <- terra::extract(prcpSum_Aug, temp_snotel)[,2]
    temp_snotel$sep_prcpSum <- terra::extract(prcpSum_Sep, temp_snotel)[,2]
    temp_snotel$oct_prcpSum <- terra::extract(prcpSum_Oct, temp_snotel)[,2]
    temp_snotel$nov_prcpSum <- terra::extract(prcpSum_Nov, temp_snotel)[,2]
    temp_snotel$dec_prcpSum <- terra::extract(prcpSum_Dec, temp_snotel)[,2]
    ## prcpMean
    temp_snotel$jan_prcpMean <- terra::extract(prcpMean_Jan, temp_snotel)[,2]
    temp_snotel$feb_prcpMean <- terra::extract(prcpMean_Feb, temp_snotel)[,2]
    temp_snotel$mar_prcpMean <- terra::extract(prcpMean_Mar, temp_snotel)[,2]
    temp_snotel$apr_prcpMean <- terra::extract(prcpMean_Apr, temp_snotel)[,2]
    temp_snotel$may_prcpMean <- terra::extract(prcpMean_May, temp_snotel)[,2]
    temp_snotel$june_prcpMean <- terra::extract(prcpMean_June, temp_snotel)[,2]
    temp_snotel$july_prcpMean <- terra::extract(prcpMean_July, temp_snotel)[,2]
    temp_snotel$aug_prcpMean <- terra::extract(prcpMean_Aug, temp_snotel)[,2]
    temp_snotel$sep_prcpMean <- terra::extract(prcpMean_Sep, temp_snotel)[,2]
    temp_snotel$oct_prcpMean <- terra::extract(prcpMean_Oct, temp_snotel)[,2]
    temp_snotel$nov_prcpMean <- terra::extract(prcpMean_Nov, temp_snotel)[,2]
    temp_snotel$dec_prcpMean <- terra::extract(prcpMean_Dec, temp_snotel)[,2]

    ## tmin
    temp_snotel$jan_tmin <- terra::extract(tmin_Jan, temp_snotel)[,2]
    temp_snotel$feb_tmin <- terra::extract(tmin_Feb, temp_snotel)[,2]
    temp_snotel$mar_tmin <- terra::extract(tmin_Mar, temp_snotel)[,2]
    temp_snotel$apr_tmin <- terra::extract(tmin_Apr, temp_snotel)[,2]
    temp_snotel$may_tmin <- terra::extract(tmin_May, temp_snotel)[,2]
    temp_snotel$june_tmin <- terra::extract(tmin_June, temp_snotel)[,2]
    temp_snotel$july_tmin <- terra::extract(tmin_July, temp_snotel)[,2]
    temp_snotel$aug_tmin <- terra::extract(tmin_Aug, temp_snotel)[,2]
    temp_snotel$sep_tmin <- terra::extract(tmin_Sep, temp_snotel)[,2]
    temp_snotel$oct_tmin <- terra::extract(tmin_Oct, temp_snotel)[,2]
    temp_snotel$nov_tmin <- terra::extract(tmin_Nov, temp_snotel)[,2]
    temp_snotel$dec_tmin <- terra::extract(tmin_Dec, temp_snotel)[,2]

    ## tmean
    temp_snotel$jan_tmean <- terra::extract(tmean_Jan, temp_snotel)[,2]
    temp_snotel$feb_tmean <- terra::extract(tmean_Feb, temp_snotel)[,2]
    temp_snotel$mar_tmean <- terra::extract(tmean_Mar, temp_snotel)[,2]
    temp_snotel$apr_tmean <- terra::extract(tmean_Apr, temp_snotel)[,2]
    temp_snotel$may_tmean <- terra::extract(tmean_May, temp_snotel)[,2]
    temp_snotel$june_tmean <- terra::extract(tmean_June, temp_snotel)[,2]
    temp_snotel$july_tmean <- terra::extract(tmean_July, temp_snotel)[,2]
    temp_snotel$aug_tmean <- terra::extract(tmean_Aug, temp_snotel)[,2]
    temp_snotel$sep_tmean <- terra::extract(tmean_Sep, temp_snotel)[,2]
    temp_snotel$oct_tmean <- terra::extract(tmean_Oct, temp_snotel)[,2]
    temp_snotel$nov_tmean <- terra::extract(tmean_Nov, temp_snotel)[,2]
    temp_snotel$dec_tmean <- terra::extract(tmean_Dec, temp_snotel)[,2]

    ## tmax
    temp_snotel$jan_tmax <- terra::extract(tmax_Jan, temp_snotel)[,2]
    temp_snotel$feb_tmax <- terra::extract(tmax_Feb, temp_snotel)[,2]
    temp_snotel$mar_tmax <- terra::extract(tmax_Mar, temp_snotel)[,2]
    temp_snotel$apr_tmax <- terra::extract(tmax_Apr, temp_snotel)[,2]
    temp_snotel$may_tmax <- terra::extract(tmax_May, temp_snotel)[,2]
    temp_snotel$june_tmax <- terra::extract(tmax_June, temp_snotel)[,2]
    temp_snotel$july_tmax <- terra::extract(tmax_July, temp_snotel)[,2]
    temp_snotel$aug_tmax <- terra::extract(tmax_Aug, temp_snotel)[,2]
    temp_snotel$sep_tmax <- terra::extract(tmax_Sep, temp_snotel)[,2]
    temp_snotel$oct_tmax <- terra::extract(tmax_Oct, temp_snotel)[,2]
    temp_snotel$nov_tmax <- terra::extract(tmax_Nov, temp_snotel)[,2]
    temp_snotel$dec_tmax <- terra::extract(tmax_Dec, temp_snotel)[,2]

    ## srad
    temp_snotel$jan_srad <- terra::extract(srad_Jan, temp_snotel)[,2]
    temp_snotel$feb_srad <- terra::extract(srad_Feb, temp_snotel)[,2]
    temp_snotel$mar_srad <- terra::extract(srad_Mar, temp_snotel)[,2]
    temp_snotel$apr_srad <- terra::extract(srad_Apr, temp_snotel)[,2]
    temp_snotel$may_srad <- terra::extract(srad_May, temp_snotel)[,2]
    temp_snotel$june_srad <- terra::extract(srad_June, temp_snotel)[,2]
    temp_snotel$july_srad <- terra::extract(srad_July, temp_snotel)[,2]
    temp_snotel$aug_srad <- terra::extract(srad_Aug, temp_snotel)[,2]
    temp_snotel$sep_srad <- terra::extract(srad_Sep, temp_snotel)[,2]
    temp_snotel$oct_srad <- terra::extract(srad_Oct, temp_snotel)[,2]
    temp_snotel$nov_srad <- terra::extract(srad_Nov, temp_snotel)[,2]
    temp_snotel$dec_srad <- terra::extract(srad_Dec, temp_snotel)[,2]


    ### Extracting climatology/CDM aggregates
    ## Oct-May aggregates
    temp_snotel$OctMay_prcpSumCDMSum <- terra::extract(prcpSum_OctMaySum_rast, temp_snotel)[,2]
    temp_snotel$OctMay_tmeanCDMSum <- terra::extract(tmean_OctMaySum_rast, temp_snotel)[,2]
    temp_snotel$OctMay_tminMean <- terra::extract(tmin_OctMayMean_rast, temp_snotel)[,2]
    temp_snotel$OctMay_tmaxMean <- terra::extract(tmax_OctMayMean_rast, temp_snotel)[,2]
    temp_snotel$OctMay_sradMean <- terra::extract(srad_OctMayMean_rast, temp_snotel)[,2]
    ## Oct-Apr aggregates
    temp_snotel$OctApr_prcpSumCDMSum <- terra::extract(prcpSum_OctAprSum_rast, temp_snotel)[,2]
    temp_snotel$OctApr_tmeanCDMSum <- terra::extract(tmean_OctAprSum_rast, temp_snotel)[,2]
    temp_snotel$OctApr_tminMean <- terra::extract(tmin_OctAprMean_rast, temp_snotel)[,2]
    temp_snotel$OctApr_tmaxMean <- terra::extract(tmax_OctAprMean_rast, temp_snotel)[,2]
    temp_snotel$OctApr_sradMean <- terra::extract(srad_OctAprMean_rast, temp_snotel)[,2]
    ## Sep-May aggregates
    temp_snotel$SepMay_prcpSumCDMSum <- terra::extract(prcpSum_SepMaySum_rast, temp_snotel)[,2]
    temp_snotel$SepMay_tmeanCDMSum <- terra::extract(tmean_SepMaySum_rast, temp_snotel)[,2]
    temp_snotel$SepMay_tminMean <- terra::extract(tmin_SepMayMean_rast, temp_snotel)[,2]
    temp_snotel$SepMay_tmaxMean <- terra::extract(tmax_SepMayMean_rast, temp_snotel)[,2]
    temp_snotel$SepMay_sradMean <- terra::extract(srad_SepMayMean_rast, temp_snotel)[,2]
    ## Dec-Feb aggregates
    temp_snotel$DecFeb_prcpSumCDMSum <- terra::extract(prcpSum_DecFebSum_rast, temp_snotel)[,2]
    temp_snotel$DecFeb_tmeanCDMSum <- terra::extract(tmean_DecFebSum_rast, temp_snotel)[,2]
    temp_snotel$DecFeb_tminMean <- terra::extract(tmin_DecFebMean_rast, temp_snotel)[,2]
    temp_snotel$DecFeb_tmaxMean <- terra::extract(tmax_DecFebMean_rast, temp_snotel)[,2]
    temp_snotel$DecFeb_sradMean <- terra::extract(srad_DecFebMean_rast, temp_snotel)[,2]
    temp_snotel$landcover = terra::extract(temp_lc_rast, temp_snotel)[,2]
    temp_snotel$landcover_triclass = terra::extract(temp_lc_TriClass_rast, temp_snotel)[,2]
  } else{
    temp_snotel2 <- Snotel_AK |>
      filter(WaterYear == years[year])
    temp_lc_rast <- lc_AK_rast_list[[year]]
    temp_lc_TriClass_rast <- LC_AK_TriClass_rast_list[[year]]
    ### extracting climatologies
    ## prcpSum
    temp_snotel2$jan_prcpSum <- terra::extract(prcpSum_Jan, temp_snotel2)[,2]
    temp_snotel2$feb_prcpSum <- terra::extract(prcpSum_Feb, temp_snotel2)[,2]
    temp_snotel2$mar_prcpSum <- terra::extract(prcpSum_Mar, temp_snotel2)[,2]
    temp_snotel2$apr_prcpSum <- terra::extract(prcpSum_Apr, temp_snotel2)[,2]
    temp_snotel2$may_prcpSum <- terra::extract(prcpSum_May, temp_snotel2)[,2]
    temp_snotel2$june_prcpSum <- terra::extract(prcpSum_June, temp_snotel2)[,2]
    temp_snotel2$july_prcpSum <- terra::extract(prcpSum_July, temp_snotel2)[,2]
    temp_snotel2$aug_prcpSum <- terra::extract(prcpSum_Aug, temp_snotel2)[,2]
    temp_snotel2$sep_prcpSum <- terra::extract(prcpSum_Sep, temp_snotel2)[,2]
    temp_snotel2$oct_prcpSum <- terra::extract(prcpSum_Oct, temp_snotel2)[,2]
    temp_snotel2$nov_prcpSum <- terra::extract(prcpSum_Nov, temp_snotel2)[,2]
    temp_snotel2$dec_prcpSum <- terra::extract(prcpSum_Dec, temp_snotel2)[,2]
    ## prcpMean
    temp_snotel2$jan_prcpMean <- terra::extract(prcpMean_Jan, temp_snotel2)[,2]
    temp_snotel2$feb_prcpMean <- terra::extract(prcpMean_Feb, temp_snotel2)[,2]
    temp_snotel2$mar_prcpMean <- terra::extract(prcpMean_Mar, temp_snotel2)[,2]
    temp_snotel2$apr_prcpMean <- terra::extract(prcpMean_Apr, temp_snotel2)[,2]
    temp_snotel2$may_prcpMean <- terra::extract(prcpMean_May, temp_snotel2)[,2]
    temp_snotel2$june_prcpMean <- terra::extract(prcpMean_June, temp_snotel2)[,2]
    temp_snotel2$july_prcpMean <- terra::extract(prcpMean_July, temp_snotel2)[,2]
    temp_snotel2$aug_prcpMean <- terra::extract(prcpMean_Aug, temp_snotel2)[,2]
    temp_snotel2$sep_prcpMean <- terra::extract(prcpMean_Sep, temp_snotel2)[,2]
    temp_snotel2$oct_prcpMean <- terra::extract(prcpMean_Oct, temp_snotel2)[,2]
    temp_snotel2$nov_prcpMean <- terra::extract(prcpMean_Nov, temp_snotel2)[,2]
    temp_snotel2$dec_prcpMean <- terra::extract(prcpMean_Dec, temp_snotel2)[,2]

    ## tmin
    temp_snotel2$jan_tmin <- terra::extract(tmin_Jan, temp_snotel2)[,2]
    temp_snotel2$feb_tmin <- terra::extract(tmin_Feb, temp_snotel2)[,2]
    temp_snotel2$mar_tmin <- terra::extract(tmin_Mar, temp_snotel2)[,2]
    temp_snotel2$apr_tmin <- terra::extract(tmin_Apr, temp_snotel2)[,2]
    temp_snotel2$may_tmin <- terra::extract(tmin_May, temp_snotel2)[,2]
    temp_snotel2$june_tmin <- terra::extract(tmin_June, temp_snotel2)[,2]
    temp_snotel2$july_tmin <- terra::extract(tmin_July, temp_snotel2)[,2]
    temp_snotel2$aug_tmin <- terra::extract(tmin_Aug, temp_snotel2)[,2]
    temp_snotel2$sep_tmin <- terra::extract(tmin_Sep, temp_snotel2)[,2]
    temp_snotel2$oct_tmin <- terra::extract(tmin_Oct, temp_snotel2)[,2]
    temp_snotel2$nov_tmin <- terra::extract(tmin_Nov, temp_snotel2)[,2]
    temp_snotel2$dec_tmin <- terra::extract(tmin_Dec, temp_snotel2)[,2]

    ## tmean
    temp_snotel2$jan_tmean <- terra::extract(tmean_Jan, temp_snotel2)[,2]
    temp_snotel2$feb_tmean <- terra::extract(tmean_Feb, temp_snotel2)[,2]
    temp_snotel2$mar_tmean <- terra::extract(tmean_Mar, temp_snotel2)[,2]
    temp_snotel2$apr_tmean <- terra::extract(tmean_Apr, temp_snotel2)[,2]
    temp_snotel2$may_tmean <- terra::extract(tmean_May, temp_snotel2)[,2]
    temp_snotel2$june_tmean <- terra::extract(tmean_June, temp_snotel2)[,2]
    temp_snotel2$july_tmean <- terra::extract(tmean_July, temp_snotel2)[,2]
    temp_snotel2$aug_tmean <- terra::extract(tmean_Aug, temp_snotel2)[,2]
    temp_snotel2$sep_tmean <- terra::extract(tmean_Sep, temp_snotel2)[,2]
    temp_snotel2$oct_tmean <- terra::extract(tmean_Oct, temp_snotel2)[,2]
    temp_snotel2$nov_tmean <- terra::extract(tmean_Nov, temp_snotel2)[,2]
    temp_snotel2$dec_tmean <- terra::extract(tmean_Dec, temp_snotel2)[,2]

    ## tmax
    temp_snotel2$jan_tmax <- terra::extract(tmax_Jan, temp_snotel2)[,2]
    temp_snotel2$feb_tmax <- terra::extract(tmax_Feb, temp_snotel2)[,2]
    temp_snotel2$mar_tmax <- terra::extract(tmax_Mar, temp_snotel2)[,2]
    temp_snotel2$apr_tmax <- terra::extract(tmax_Apr, temp_snotel2)[,2]
    temp_snotel2$may_tmax <- terra::extract(tmax_May, temp_snotel2)[,2]
    temp_snotel2$june_tmax <- terra::extract(tmax_June, temp_snotel2)[,2]
    temp_snotel2$july_tmax <- terra::extract(tmax_July, temp_snotel2)[,2]
    temp_snotel2$aug_tmax <- terra::extract(tmax_Aug, temp_snotel2)[,2]
    temp_snotel2$sep_tmax <- terra::extract(tmax_Sep, temp_snotel2)[,2]
    temp_snotel2$oct_tmax <- terra::extract(tmax_Oct, temp_snotel2)[,2]
    temp_snotel2$nov_tmax <- terra::extract(tmax_Nov, temp_snotel2)[,2]
    temp_snotel2$dec_tmax <- terra::extract(tmax_Dec, temp_snotel2)[,2]

    ## srad
    temp_snotel2$jan_srad <- terra::extract(srad_Jan, temp_snotel2)[,2]
    temp_snotel2$feb_srad <- terra::extract(srad_Feb, temp_snotel2)[,2]
    temp_snotel2$mar_srad <- terra::extract(srad_Mar, temp_snotel2)[,2]
    temp_snotel2$apr_srad <- terra::extract(srad_Apr, temp_snotel2)[,2]
    temp_snotel2$may_srad <- terra::extract(srad_May, temp_snotel2)[,2]
    temp_snotel2$june_srad <- terra::extract(srad_June, temp_snotel2)[,2]
    temp_snotel2$july_srad <- terra::extract(srad_July, temp_snotel2)[,2]
    temp_snotel2$aug_srad <- terra::extract(srad_Aug, temp_snotel2)[,2]
    temp_snotel2$sep_srad <- terra::extract(srad_Sep, temp_snotel2)[,2]
    temp_snotel2$oct_srad <- terra::extract(srad_Oct, temp_snotel2)[,2]
    temp_snotel2$nov_srad <- terra::extract(srad_Nov, temp_snotel2)[,2]
    temp_snotel2$dec_srad <- terra::extract(srad_Dec, temp_snotel2)[,2]


    ### Extracting climatology/CDM aggregates
    ## Oct-May aggregates
    temp_snotel2$OctMay_prcpSumCDMSum <- terra::extract(prcpSum_OctMaySum_rast, temp_snotel2)[,2]
    temp_snotel2$OctMay_tmeanCDMSum <- terra::extract(tmean_OctMaySum_rast, temp_snotel2)[,2]
    temp_snotel2$OctMay_tminMean <- terra::extract(tmin_OctMayMean_rast, temp_snotel2)[,2]
    temp_snotel2$OctMay_tmaxMean <- terra::extract(tmax_OctMayMean_rast, temp_snotel2)[,2]
    temp_snotel2$OctMay_sradMean <- terra::extract(srad_OctMayMean_rast, temp_snotel2)[,2]
    ## Oct-Apr aggregates
    temp_snotel2$OctApr_prcpSumCDMSum <- terra::extract(prcpSum_OctAprSum_rast, temp_snotel2)[,2]
    temp_snotel2$OctApr_tmeanCDMSum <- terra::extract(tmean_OctAprSum_rast, temp_snotel2)[,2]
    temp_snotel2$OctApr_tminMean <- terra::extract(tmin_OctAprMean_rast, temp_snotel2)[,2]
    temp_snotel2$OctApr_tmaxMean <- terra::extract(tmax_OctAprMean_rast, temp_snotel2)[,2]
    temp_snotel2$OctApr_sradMean <- terra::extract(srad_OctAprMean_rast, temp_snotel2)[,2]
    ## Sep-May aggregates
    temp_snotel2$SepMay_prcpSumCDMSum <- terra::extract(prcpSum_SepMaySum_rast, temp_snotel2)[,2]
    temp_snotel2$SepMay_tmeanCDMSum <- terra::extract(tmean_SepMaySum_rast, temp_snotel2)[,2]
    temp_snotel2$SepMay_tminMean <- terra::extract(tmin_SepMayMean_rast, temp_snotel2)[,2]
    temp_snotel2$SepMay_tmaxMean <- terra::extract(tmax_SepMayMean_rast, temp_snotel2)[,2]
    temp_snotel2$SepMay_sradMean <- terra::extract(srad_SepMayMean_rast, temp_snotel2)[,2]
    ## Dec-Feb aggregates
    temp_snotel2$DecFeb_prcpSumCDMSum <- terra::extract(prcpSum_DecFebSum_rast, temp_snotel2)[,2]
    temp_snotel2$DecFeb_tmeanCDMSum <- terra::extract(tmean_DecFebSum_rast, temp_snotel2)[,2]
    temp_snotel2$DecFeb_tminMean <- terra::extract(tmin_DecFebMean_rast, temp_snotel2)[,2]
    temp_snotel2$DecFeb_tmaxMean <- terra::extract(tmax_DecFebMean_rast, temp_snotel2)[,2]
    temp_snotel2$DecFeb_sradMean <- terra::extract(srad_DecFebMean_rast, temp_snotel2)[,2]
    temp_snotel2$landcover = terra::extract(temp_lc_rast, temp_snotel2)[,2]
    temp_snotel2$landcover_triclass = terra::extract(temp_lc_TriClass_rast, temp_snotel2)[,2]
    temp_snotel <- rbind(temp_snotel, temp_snotel2)
  }
}
### storing combined years with landcover added back in the Snotel_AK variable
Snotel_AK <- temp_snotel

##### writing out annual Alaska gpkg with covariates #####
# write_sf(Snotel_AK, here("Data", "SNOTEL", "Alaska", "GIS", "SHP", "SnotelAK_PeakSWE_Covars.shp"), append = FALSE)
write_sf(Snotel_AK, here("Data", "SNOTEL", "Alaska", "GIS", "GPKG", "Annual", "SnotelAK_PeakSWE_Covars.gpkg"), append = FALSE)

AKEndtime <- Sys.time()
AKtime <- AKEndtime - AKStarttime


##### Combining Snotel_CONUS and Snotel_AK for model fitting/evaluation #####
combinestarttime <- Sys.time()
### reading in files
# Snotel_CONUS_clim <- read_sf(here("Data", "SNOTEL", "CONUS", "GIS", "GPKG", "Climatology", "SnotelCONUS_PeakSWEClimatology_Covars.gpkg"))
# Snotel_AK_clim <- read_sf(here("Data", "SNOTEL", "Alaska", "GIS", "GPKG", "Climatology", "SnotelAK_PeakSWEClimatology_Covars.gpkg"))
# Snotel_CONUS <- read_sf(here("Data", "SNOTEL", "CONUS", "GIS", "GPKG", "Annual", "SnotelCONUS_PeakSWE_Covars.gpkg"))
Snotel_CONUS <- read_sf(here("Data", "SNOTEL", "CONUS", "GIS", "GPKG", "Annual", "SnotelCONUSAnnualZeroPts_PeakSWE_Covars.gpkg"))
Snotel_AK <- read_sf(here("Data", "SNOTEL", "Alaska", "GIS", "GPKG", "Annual", "SnotelAK_PeakSWE_Covars.gpkg"))

### binding appropriate files
# Snotel_Combined_Clim <- rbind(Snotel_CONUS_clim, st_transform(Snotel_AK_clim, crs(Snotel_CONUS_clim)))
Snotel_Combined <- rbind(Snotel_CONUS, st_transform(Snotel_AK, crs(Snotel_CONUS)))

# write_sf(Snotel_Combined_Clim, here("Data", "SNOTEL", "Combined", "GIS", "SnotelCombined_ClimatologyCovars.gpkg"), append = FALSE)
# write_sf(Snotel_Combined, here("Data", "SNOTEL", "Combined", "GIS", "SnotelCombined_AnnualCovars.gpkg"), append = FALSE)
write_sf(Snotel_Combined, here("Data", "SNOTEL", "Combined", "GIS", "SnotelCombined_AnnualCovars_AnnualZeroPts.gpkg"), append = FALSE)

combineendtime <- Sys.time()
combinetime <- combineendtime - combinestarttime

# CONUSClimtime
# AKClimtime
CONUStime
AKtime
combinetime
Totalendtime <- Sys.time()
Totaltime <- Totalendtime - Totalstarttime
Totaltime
##### End of Script #####