##### Script to be used for downloading SNOTEL, Daymet, Elevation, and NLCD data for my thesis work #####
### script by Griffin Shelor
TotalStartTime <- Sys.time()

##### installing/loading packages #####
# install.packages(c("pacman", "here", "tidyverse", "daymetr", "snotelr", "rstac", "rsi"))
library(pacman)
p_load(here, tidyverse, daymetr, snotelr, rstac, rsi, ncdf4, parallel, svMisc, sf, terra)
# here()

### putting this function here, gonna use it later
`%nin%` = Negate(`%in%`)

### Snotel Download Start time
SnotelStartTime <- Sys.time()
##### Downloading SNOTEL data #####
### calling snotel_info first to initialize what sites I actually want to download
SnotelInfo <- snotel_info()


### ensuring start and end columns in SnotelInfo are in date format
SnotelInfo$start <- ymd(SnotelInfo$start)
SnotelInfo$end <- ymd(SnotelInfo$end)

### filtering to only include sites that cover from October 1, 1992 to September 30, 2020
SnotelInfo_Clean <- SnotelInfo |>
  # filter(site_id != 201 & site_id != 1315 & site_id != 623 & site_id != 870 & site_id != 1326 & site_id != 453) |>
  # filter(start <= "1992-10-01" & end >= "2020-09-30") |>
  arrange(site_id)

### separating Snotel info by whether it is in Alaska or not
SnotelInfo_AK <- SnotelInfo_Clean |>
  filter(state == "AK") |>
  ### site 201 had issues
  filter(site_id != 201)

SnotelInfo_CONUS <- SnotelInfo_Clean |>
  filter(state != "AK") |>
  ### spratt creek had issues downloading
  filter(site_id != 778)

### downloading snotel data for conus fails at site 657 for some reason, but it works when I download it separately????
# SnotelData_657 <- snotel_download(657, internal = TRUE)
# SnotelData_778 <- snotel_download(778, internal = TRUE)


### now downloading SNOTEL data only for sites that cover my whole study period
SnotelData_AK <- snotel_download(site_id = SnotelInfo_AK$site_id, internal = TRUE)
SnotelData_CONUS <- snotel_download(site_id = SnotelInfo_CONUS$site_id, internal = TRUE)

### binding site 657 and rest of CONUS Snotel data back together
SnotelData_CONUS <- rbind(SnotelData_CONUS, SnotelData_778)

### removing NAs, arranging by both site_id and date
SnotelData_AK <- SnotelData_AK |>
  drop_na(snow_water_equivalent) |>
  separate(date, into = c("Year", "Month", "Day"), sep = "-") |>
  mutate(date = ymd(paste0(Year, "-", Month, "-", Day)), .before = 15)

### Converting individual year, month, and day columns to be numeric
SnotelData_AK$Year <- as.numeric(SnotelData_AK$Year)
SnotelData_AK$Month <- as.numeric(SnotelData_AK$Month)
SnotelData_AK$Day <- as.numeric(SnotelData_AK$Day)

### filtering out just relevant dates, adding new columns to assist with grabbing daymet data
SnotelData_AK <- SnotelData_AK |>
  filter(date > "1992-09-30") |>
  filter(date < "2020-10-01")

### Daymet is annoyingly weird about leap years so I'm separating those out and creating the yearday column differently for those years
### it also drops Dec 31 for some reason but I'm going to keep SNOTEL obs from that date since I'm really after SWE_max, but why is Oak Ridge National Lab doing this to me
SnotelData_AK_LeapYears <- SnotelData_AK |>
  filter(Year %in% c(1992, 1996, 2000, 2004, 2008, 2012, 2016, 2020)) |>
  mutate(yearday = case_when(Month == 1 ~ Day,
                             Month == 2 ~ 31 + Day,
                             Month == 3 ~ 31 + 29 + Day,
                             Month == 4 ~ 31 + 29 + 31 + Day,
                             Month == 5 ~ 31 + 29 + 31 + 30 + Day,
                             Month == 6 ~ 31 + 29 + 31 + 30 + 31 + Day,
                             Month == 7 ~ 31 + 29 + 31 + 30 + 31 + 30 + Day,
                             Month == 8 ~ 31 + 29 + 31 + 30 + 31 + 30 + 31 + Day,
                             Month == 9 ~ 31 + 29 + 31 + 30 + 31 + 30 + 31 + 31 + Day,
                             Month == 10 ~ 31 + 29 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + Day,
                             Month == 11 ~ 31 + 29 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + Day,
                             TRUE ~ 31 + 29 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + 30 + Day), .before = 15)

SnotelData_AK_NonLeapYears <- SnotelData_AK |>
  filter(Year %nin% c(1992, 1996, 2000, 2004, 2008, 2012, 2016, 2020)) |>
  mutate(yearday = case_when(Month == 1 ~ Day,
                             Month == 2 ~ 31 + Day,
                             Month == 3 ~ 31 + 28 + Day,
                             Month == 4 ~ 31 + 28 + 31 + Day,
                             Month == 5 ~ 31 + 28 + 31 + 30 + Day,
                             Month == 6 ~ 31 + 28 + 31 + 30 + 31 + Day,
                             Month == 7 ~ 31 + 28 + 31 + 30 + 31 + 30 + Day,
                             Month == 8 ~ 31 + 28 + 31 + 30 + 31 + 30 + 31 + Day,
                             Month == 9 ~ 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + Day,
                             Month == 10 ~ 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + Day,
                             Month == 11 ~ 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + Day,
                             TRUE ~ 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + 30 + Day), .before = 15)

### binding leap year data to non leap year data
SnotelData_AK <- rbind(SnotelData_AK_LeapYears, SnotelData_AK_NonLeapYears) |>
  mutate(WaterYear = case_when(Month >= 10 ~ Year + 1,
                               TRUE ~ Year),
         study_year = WaterYear - min(WaterYear)) |>
  arrange(date, site_id)

### filtering out observations from each year for each site where the site does not have observations for most of the water year
AKSiteWYObsCounts <- SnotelData_AK |>
  group_by(site_id, WaterYear) |>
  summarise(n_obs = n()) |>
  ungroup() |>
  filter(n_obs > 300)

AKunique_ids <- unique(AKSiteWYObsCounts$site_id)
for (i in 1:length(AKunique_ids)){
  if (i == 1){
    temp_years_df <- AKSiteWYObsCounts |>
      filter(site_id == AKunique_ids[i])
    temp_years <- temp_years_df$WaterYear
    SnotelAK_clean <- SnotelData_AK |>
      filter(site_id == AKunique_ids[i] & WaterYear %in% temp_years)
  } else{
    temp_years_df <- AKSiteWYObsCounts |>
      filter(site_id == AKunique_ids[i])
    temp_years <- temp_years_df$WaterYear
    SnotelAK2 <- SnotelData_AK |>
      filter(site_id == AKunique_ids[i] & WaterYear %in% temp_years)
    SnotelAK_clean <- rbind(SnotelAK_clean, SnotelAK2)
  }
}

SnotelData_AK <- SnotelAK_clean
rm(SnotelAK_clean)

### now repeating process for CONUS
### removing NAs, arranging by both site_id and date
SnotelData_CONUS <- SnotelData_CONUS |>
  drop_na(snow_water_equivalent) |>
  separate(date, into = c("Year", "Month", "Day"), sep = "-") |>
  mutate(date = ymd(paste0(Year, "-", Month, "-", Day)), .before = 15)

### Converting individual year, month, and day columns to be numeric
SnotelData_CONUS$Year <- as.numeric(SnotelData_CONUS$Year)
SnotelData_CONUS$Month <- as.numeric(SnotelData_CONUS$Month)
SnotelData_CONUS$Day <- as.numeric(SnotelData_CONUS$Day)

### filtering out just relevant dates, adding new columns to assist with grabbing daymet data
SnotelData_CONUS <- SnotelData_CONUS |>
  filter(date > "1992-09-30") |>
  filter(date < "2020-10-01")

### Daymet is annoyingly weird about leap years so I'm separating those out and creating the yearday column differently for those years
### it drops Dec 31 for some reason, why is Oak Ridge National Lab doing this to me
SnotelData_CONUS_LeapYears <- SnotelData_CONUS |>
  filter(Year %in% c(1992, 1996, 2000, 2004, 2008, 2012, 2016, 2020)) |>
  mutate(yearday = case_when(Month == 1 ~ Day,
                             Month == 2 ~ 31 + Day,
                             Month == 3 ~ 31 + 29 + Day,
                             Month == 4 ~ 31 + 29 + 31 + Day,
                             Month == 5 ~ 31 + 29 + 31 + 30 + Day,
                             Month == 6 ~ 31 + 29 + 31 + 30 + 31 + Day,
                             Month == 7 ~ 31 + 29 + 31 + 30 + 31 + 30 + Day,
                             Month == 8 ~ 31 + 29 + 31 + 30 + 31 + 30 + 31 + Day,
                             Month == 9 ~ 31 + 29 + 31 + 30 + 31 + 30 + 31 + 31 + Day,
                             Month == 10 ~ 31 + 29 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + Day,
                             Month == 11 ~ 31 + 29 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + Day,
                             TRUE ~ 31 + 29 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + 30 + Day), .before = 15)


SnotelData_CONUS_NonLeapYears <- SnotelData_CONUS |>
  filter(Year %nin% c(1992, 1996, 2000, 2004, 2008, 2012, 2016, 2020)) |>
  mutate(yearday = case_when(Month == 1 ~ Day,
                             Month == 2 ~ 31 + Day,
                             Month == 3 ~ 31 + 28 + Day,
                             Month == 4 ~ 31 + 28 + 31 + Day,
                             Month == 5 ~ 31 + 28 + 31 + 30 + Day,
                             Month == 6 ~ 31 + 28 + 31 + 30 + 31 + Day,
                             Month == 7 ~ 31 + 28 + 31 + 30 + 31 + 30 + Day,
                             Month == 8 ~ 31 + 28 + 31 + 30 + 31 + 30 + 31 + Day,
                             Month == 9 ~ 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + Day,
                             Month == 10 ~ 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + Day,
                             Month == 11 ~ 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + Day,
                             TRUE ~ 31 + 28 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + 30 + Day), .before = 15)
### binding leap year data to non leap year data
SnotelData_CONUS <- rbind(SnotelData_CONUS_LeapYears, SnotelData_CONUS_NonLeapYears) |>
  mutate(WaterYear = case_when(Month >= 10 ~ Year + 1,
                               TRUE ~ Year),
         study_year = WaterYear - min(WaterYear)) |>
  arrange(date, site_id)

### filtering out observations from each year for each site where the site does not have observations for most of the water year
CONUSSiteWYObsCounts <- SnotelData_CONUS |>
  group_by(site_id, WaterYear) |>
  summarise(n_obs = n()) |>
  ungroup() |>
  filter(n_obs > 300)

CONUSunique_ids <- unique(CONUSSiteWYObsCounts$site_id)
for (i in 1:length(CONUSunique_ids)){
  if (i == 1){
    temp_years_df <- CONUSSiteWYObsCounts |>
      filter(site_id == CONUSunique_ids[i])
    temp_years <- temp_years_df$WaterYear
    SnotelCONUS_clean <- SnotelData_CONUS |>
      filter(site_id == CONUSunique_ids[i] & WaterYear %in% temp_years)
  } else{
    temp_years_df <- CONUSSiteWYObsCounts |>
      filter(site_id == CONUSunique_ids[i])
    temp_years <- temp_years_df$WaterYear
    SnotelCONUS2 <- SnotelData_CONUS |>
      filter(site_id == CONUSunique_ids[i] & WaterYear %in% temp_years)
    SnotelCONUS_clean <- rbind(SnotelCONUS_clean, SnotelCONUS2)
  }
}

SnotelData_CONUS <- SnotelCONUS_clean
rm(SnotelCONUS_clean)

### writing SNOTEL data as csvs
write_csv(SnotelData_AK, here("Data", "SNOTEL", "Alaska", "SnotelData_AK.csv"))
write_csv(SnotelData_CONUS, here("Data", "SNOTEL", "CONUS", "SnotelData_CONUS.csv"))

### reading in SNOTEL data from csvs which were saved in a previous run of this script
SnotelData_AK <- read_csv(here("Data", "SNOTEL", "Alaska", "SnotelData_AK.csv"))
SnotelData_CONUS <- read_csv(here("Data", "SNOTEL", "CONUS", "SnotelData_CONUS.csv"))

### converting SNOTEL data to extract peak SWE at each SNOTEL site for each year
### Alaska first
### separating out each site, grouping by Water year, identifying peak SWE
uniqueAKids <- unique(SnotelData_AK$site_id)
for (id in 1:length(uniqueAKids)){
  if (id == 1){
    temp_id <- SnotelData_AK |>
      filter(site_id == uniqueAKids[id])
    relevant_info <- SnotelInfo_AK |>
      filter(site_id == uniqueAKids[id]) |>
      select(site_id, site_name, description, state, start, end, longitude, latitude)
    SnotelData_AK_PeakSWE <- temp_id |>
      group_by(WaterYear) |>
      summarise(peak_swe = max(snow_water_equivalent)) |>
      ### adding columns back
      # site_id, site_name, description, state, start, end, longitude, latitude
      mutate(site_id = relevant_info$site_id,
             site_name = relevant_info$site_name,
             description = relevant_info$description,
             state = relevant_info$state,
             start = relevant_info$start,
             end = relevant_info$end,
             longitude = relevant_info$longitude,
             latitude = relevant_info$latitude)
  } else{
    temp_id <- SnotelData_AK |>
      filter(site_id == uniqueAKids[id])
    relevant_info <- SnotelInfo_AK |>
      filter(site_id == uniqueAKids[id]) |>
      select(site_id, site_name, description, state, start, end, longitude, latitude)
    temp_id_peakSWE  <- temp_id |>
      group_by(WaterYear) |>
      summarise(peak_swe = max(snow_water_equivalent)) |>
      ### adding columns back
      # site_id, site_name, description, state, start, end, longitude, latitude
      mutate(site_id = relevant_info$site_id,
             site_name = relevant_info$site_name,
             description = relevant_info$description,
             state = relevant_info$state,
             start = relevant_info$start,
             end = relevant_info$end,
             longitude = relevant_info$longitude,
             latitude = relevant_info$latitude)
    SnotelData_AK_PeakSWE <- rbind(SnotelData_AK_PeakSWE, temp_id_peakSWE)
  }
}


### Now CONUS
### separating out each site, grouping by Water year, identifying peak SWE
uniqueCONUSids <- unique(SnotelData_CONUS$site_id)
for (id in 1:length(uniqueCONUSids)){
  if (id == 1){
    temp_id <- SnotelData_CONUS |>
      filter(site_id == uniqueCONUSids[id])
    relevant_info <- SnotelInfo_CONUS |>
      filter(site_id == uniqueCONUSids[id]) |>
      select(site_id, site_name, description, state, start, end, longitude, latitude)
    SnotelData_CONUS_PeakSWE <- temp_id |>
      group_by(WaterYear) |>
      summarise(peak_swe = max(snow_water_equivalent)) |>
      ### adding columns back
      # site_id, site_name, description, state, start, end, longitude, latitude
      mutate(site_id = relevant_info$site_id,
             site_name = relevant_info$site_name,
             description = relevant_info$description,
             state = relevant_info$state,
             start = relevant_info$start,
             end = relevant_info$end,
             longitude = relevant_info$longitude,
             latitude = relevant_info$latitude)
  } else{
    temp_id <- SnotelData_CONUS |>
      filter(site_id == uniqueCONUSids[id])
    relevant_info <- SnotelInfo_CONUS |>
      filter(site_id == uniqueCONUSids[id]) |>
      select(site_id, site_name, description, state, start, end, longitude, latitude)
    temp_id_peakSWE  <- temp_id |>
      group_by(WaterYear) |>
      summarise(peak_swe = max(snow_water_equivalent)) |>
      ### adding columns back
      # site_id, site_name, description, state, start, end, longitude, latitude
      mutate(site_id = relevant_info$site_id,
             site_name = relevant_info$site_name,
             description = relevant_info$description,
             state = relevant_info$state,
             start = relevant_info$start,
             end = relevant_info$end,
             longitude = relevant_info$longitude,
             latitude = relevant_info$latitude)
    SnotelData_CONUS_PeakSWE <- rbind(SnotelData_CONUS_PeakSWE, temp_id_peakSWE)
  }
}


### writing SNOTEL data as csvs
write_csv(SnotelData_AK_PeakSWE, here("Data", "SNOTEL", "Alaska", "SnotelData_AK_WY9320_PeakSWE.csv"))
write_csv(SnotelData_CONUS_PeakSWE, here("Data", "SNOTEL", "CONUS", "SnotelData_CONUS_WY9320_PeakSWE.csv"))

### reading in SNOTEL data from csvs which were saved in a previous run of this script (or right above this line)
SnotelData_AK_PeakSWE <- read_csv(here("Data", "SNOTEL", "Alaska", "SnotelData_AK_WY9320_PeakSWE.csv"))
SnotelData_CONUS_PeakSWE <- read_csv(here("Data", "SNOTEL", "CONUS", "SnotelData_CONUS_WY9320_PeakSWE.csv"))

### calculating SWE Climatologies
## Alaska first
Snotel_AK_PeakSWEClim <- SnotelData_AK_PeakSWE |>
  group_by(site_id) |>
  reframe(peak_swe = mean(peak_swe),
          site_name = site_name,
          description = description,
          state = state,
          start = start,
          end = end,
          longitude = longitude,
          latitude = latitude) |>
  ## extracting unique rows (one average peak SWE per site_id, otherwise each site ID is repeated for each year it has observations)
  unique()
## now CONUS
Snotel_CONUS_PeakSWEClim <- SnotelData_CONUS_PeakSWE |>
  group_by(site_id) |>
  reframe(peak_swe = mean(peak_swe),
          site_name = site_name,
          description = description,
          state = state,
          start = start,
          end = end,
          longitude = longitude,
          latitude = latitude) |>
  ## extracting unique rows (one average peak SWE per site_id, otherwise each site ID is repeated for each year it has observations)
  unique()

### writing SNOTEL climatologies to csv
write_csv(Snotel_AK_PeakSWEClim, here("Data", "SNOTEL", "Alaska", "Snotel_AK_WY9320_PeakSWEClimatology.csv"))
write_csv(Snotel_CONUS_PeakSWEClim, here("Data", "SNOTEL", "CONUS", "Snotel_CONUS_WY9320_PeakSWEClimatology.csv"))

### Converting csvs (both climatology and annual data) to sf object, then writing to shapefile and geopackage
### converting annual peak SWE values to sf
## Alaska first
SnotelData_AK_PeakSWE_sf <- st_as_sf(SnotelData_AK_PeakSWE, coords = c("longitude", "latitude"), crs = 4326)
## CONUS
SnotelData_CONUS_PeakSWE_sf <- st_as_sf(SnotelData_CONUS_PeakSWE, coords = c("longitude", "latitude"), crs = 4326)

### converting SWE climatologies to sf
## Alaska first
Snotel_AK_PeakSWEClim_sf <- st_as_sf(Snotel_AK_PeakSWEClim, coords = c("longitude", "latitude"), crs = 4326)
## CONUS
Snotel_CONUS_PeakSWEClim_sf <- st_as_sf(Snotel_CONUS_PeakSWEClim, coords = c("longitude", "latitude"), crs = 4326)

### reprojecting AK and CONUS peak SWE layers to correct CRS (EPSG:6393 (alaska albers) for AK, ESRI:102039 (usa contiguous albers equal area conic usgs version) for CONUS)
### reading in AOIs
AK_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "Alaska", "AK_AOI.gpkg"))
CONUS_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "CONUS", "CONUS_AOI.gpkg"))
## Alaska first
SnotelData_AK_PeakSWE_sf <- st_transform(SnotelData_AK_PeakSWE_sf, crs(AK_AOI))
## CONUS
SnotelData_CONUS_PeakSWE_sf <- st_transform(SnotelData_CONUS_PeakSWE_sf, crs(CONUS_AOI))

### reprojecting SWE climatology layers to correct projections
## Alaska first
Snotel_AK_PeakSWEClim_sf <- st_transform(Snotel_AK_PeakSWEClim_sf, crs(AK_AOI))
## CONUS
Snotel_CONUS_PeakSWEClim_sf <- st_transform(Snotel_CONUS_PeakSWEClim_sf, crs(CONUS_AOI))

### exporting annual peak SWE values as geopackage, and as shapefile
## Alaska
st_write(SnotelData_AK_PeakSWE_sf, here("Data", "SNOTEL", "Alaska", "GIS", "GPKG", "Annual", "SnotelData_AK_PeakSWE.gpkg"), append = FALSE)
## CONUS
st_write(SnotelData_CONUS_PeakSWE_sf, here("Data", "SNOTEL", "CONUS", "GIS", "GPKG", "Annual", "SnotelData_CONUS_PeakSWE.gpkg"), append = FALSE)
## now shapefile
## Alaska
st_write(SnotelData_AK_PeakSWE_sf, here("Data", "SNOTEL", "Alaska", "GIS", "SHP", "Annual", "SnotelData_AK_PeakSWE.shp"), append = FALSE)
## CONUS
st_write(SnotelData_CONUS_PeakSWE_sf, here("Data", "SNOTEL", "CONUS", "GIS", "SHP", "Annual", "SnotelData_CONUS_PeakSWE.shp"), append = FALSE)

### exporting climatologies as geopackage, and as shapefile
## Alaska
# st_write(Snotel_AK_PeakSWEClim_sf, here("Data", "SNOTEL", "Alaska", "GIS", "SnotelData_AK_PeakSWEClimatology.gpkg"), append = FALSE)
## CONUS
# st_write(Snotel_CONUS_PeakSWEClim_sf, here("Data", "SNOTEL", "CONUS", "GIS", "SnotelData_CONUS_PeakSWEClimatology.gpkg"), append = FALSE)
## now shapefile
## Alaska
# st_write(Snotel_AK_PeakSWEClim_sf, here("Data", "SNOTEL", "Alaska", "GIS", "SnotelData_AK_PeakSWEClimatology.shp"), append = FALSE)
## CONUS
# st_write(Snotel_CONUS_PeakSWEClim_sf, here("Data", "SNOTEL", "CONUS", "GIS", "SnotelData_CONUS_PeakSWEClimatology.shp"), append = FALSE)

### Snotel download end time
SnotelEndTime <- Sys.time()
SnotelDownloadTime <- SnotelEndTime - SnotelStartTime


##### Downloading Daymet netcdf files for whole study area #####
### Daymet NetCDF tile download start time
DaymetTileStartTime <- Sys.time()
DaymetTileStartTime
### Alaska first
### ended up splitting daymet downloads between like 8 different computers, but in total I am downloading daymet data from 1992-2023

### using AOI polygons to calculate min/max latitude and longitude for downloading Daymet data
AK_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "Alaska", "AK_AOI.gpkg"))
AK_Ecoregions <- read_sf(here("Data", "L3_Ecoregions_USB", "Alaska", "AKEcoregions.gpkg"))
### transforming to WGS84 for downloading daymet data
AK_AOI_WGS <- st_transform(AK_AOI, 4326)

# CONUS_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "CONUS", "CONUS_AOI.gpkg"))
# CONUS_AOI_WGS <- st_transform(CONUS_AOI, 4326)

### using st_bbox to get min/max of lat/lon
AK_bbox <- st_bbox(AK_AOI_WGS)
AK_minlon <- AK_bbox[[1]]
AK_minlat <- AK_bbox[[2]]
AK_maxlon <- AK_bbox[[3]]
AK_maxlat <- AK_bbox[[4]]
### CONUS
# CONUS_bbox <- st_bbox(CONUS_AOI)
# CONUS_minlon <- CONUS_bbox[[1]]
# CONUS_minlat <- CONUS_bbox[[2]]
# CONUS_maxlon <- CONUS_bbox[[3]]
# CONUS_maxlat <- CONUS_bbox[[4]]

### downloading variables separately
# download_daymet_tiles(
#   location = c(AK_minlat, AK_minlon, AK_maxlat, AK_maxlon),
#   # tiles = c(12601, 12602, 12603, ),
#   start = 1992,
#   end = 2020,
#   path = here("Data", "Daymet", "Alaska", "Tiles", "tmin_tmax"),
#   param = "tmin",
#   silent = FALSE,
#   force = FALSE
# )
# 
# download_daymet_tiles(
#   location = c(AK_minlat, AK_minlon, AK_maxlat, AK_maxlon),
#   # tiles,
#   start = 1992,
#   end = 2020,
#   path = here("Data", "Daymet", "Alaska", "Tiles", "tmin_tmax"),
#   param = "tmax",
#   silent = FALSE,
#   force = FALSE
# )
# 
# download_daymet_tiles(
#   location = c(AK_minlat, AK_minlon, AK_maxlat, AK_maxlon),
#   # tiles,
#   start = 1992,
#   end = 2020,
#   path = here("Data", "Daymet", "Alaska", "Tiles", "srad"),
#   param = "srad",
#   silent = FALSE,
#   force = FALSE
# )
# 
# download_daymet_tiles(
#   location = c(AK_minlat, AK_minlon, AK_maxlat, AK_maxlon),
#   # tiles,
#   start = 1992,
#   end = 2020,
#   path = here("Data", "Daymet", "Alaska", "Tiles", "prcp"),
#   param = "prcp",
#   silent = FALSE,
#   force = FALSE
# )

### CONUS next
# download_daymet_tiles(
#   location = c(CONUS_minlat, CONUS_minlon, CONUS_maxlat, CONUS_maxlon),
#   # tiles,
#   start = 1992,
#   end = 2023,
#   path = here("Data", "Daymet", "CONUS", "tmin_tmax"),
#   param = "tmin",
#   silent = FALSE,
#   force = FALSE
# )
# 
# download_daymet_tiles(
#   location = c(CONUS_minlat, CONUS_minlon, CONUS_maxlat, CONUS_maxlon),
#   # tiles,
#   start = 1992,
#   end = 2023,
#   path = here("Data", "Daymet", "CONUS", "tmin_tmax"),
#   param = "tmax",
#   silent = FALSE,
#   force = FALSE
# )
# 
# download_daymet_tiles(
#   location = c(CONUS_minlat, CONUS_minlon, CONUS_maxlat, CONUS_maxlon),
#   # tiles,
#   start = 1992,
#   end = 2023,
#   path = here("Data", "Daymet", "CONUS", "srad"),
#   param = "srad",
#   silent = FALSE,
#   force = FALSE
# )
# 
# download_daymet_tiles(
#   location = c(CONUS_minlat, CONUS_minlon, CONUS_maxlat, CONUS_maxlon),
#   # tiles,
#   start = 1992,
#   end = 2023,
#   path = here("Data", "Daymet", "CONUS", "prcp"),
#   param = "prcp",
#   silent = FALSE,
#   force = FALSE
# )

### Daymet netcdf tile end time
DaymetTileEndTime <- Sys.time()
DaymetTileDownloadTime <- DaymetTileEndTime - DaymetTileStartTime


### DEM Download Start time
DEMStartTime <- Sys.time()
##### Downloading DEM data #####
### making lists of states so the DEM will be downloaded in chunks so as not to break the download function/crash my R session
### Made a copy of the Alaska AOI and split it into chunks so the DEM download would be easier/less likely to crash
# AK_AOI_Sliced <- read_sf(here("Data", "L3_Ecoregions_USB", "Alaska", "AK_AOI_Sliced.gpkg"))
# AK_AOI_Sliced <- AK_AOI_Sliced |>
#   mutate(slice = 1:nrow(AK_AOI_Sliced))
### projecting AK_AOI_Sliced into Alaska Albers projection because get_dem() can throw an error in WGS 84
# AK_AOI_Sliced <- st_transform(AK_AOI_Sliced, 6393)

### re-reading in CONUS AOI so that DEM download won't be grabbed in WGS 84
# CONUS_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "CONUS", "CONUS_AOI.gpkg"))
# states <- c("Alaska", unique(CONUS_AOI$STATE_NAME))

### Fuck it, I'm apparently back to trying the rsi package
### downloading copernicus with the RSI package, getting USGS 3DEP data with rstac because it fails with RSI (and maybe with rstac too, it's hard to mosaic)
# for (state in states){
#   print(state)
#   if (state == "Alaska"){
#     for (slice_id in 1:nrow(AK_AOI_Sliced)){
#       print(as.character(slice_id))
      ### creating bounding box which gets fed into rstac process
      # bbox_temp <- AK_AOI_Sliced |>
      #   filter(slice == slice_id)
      
      ### creating separate object reprojected to Alaska Albers since rsi package expects projected data
      # bbox_temp_rsi <- st_transform(bbox_temp, 6393)
      
      # temp_dem <- get_dem(bbox_temp_rsi, pixel_x_size = 30, pixel_y_size = 30, output_filename = here("Data", "DEM", "Alaska", "Slices", paste0("AKDEMSlice", slice_id, ".tif")))
      
      ### ok now back to rstac method, creating bounding box in WGS 84 that rstac expects
      # bbox_temp <- st_bbox(bbox_temp)
      
      # temp_dem <- get_stac_data(bbox_temp, start_date = NULL, end_date = NULL, pixel_x_size = 30, pixel_y_size = 30, asset_names = "data", stac_source = "https://planetarycomputer.microsoft.com/api/stac/v1/", collection = "3dep-seamless", output_filename = here("Data", "DEM", "Alaska", "Slices", paste0("AKDEMSlice", slice_id, ".tif")))
      ### connecting to the appropriate stac API
      # stac_source <- rstac::stac("https://planetarycomputer.microsoft.com/api/stac/v1")
      # stac_source
      # 
      # rstac::get_request(stac_source)
      # 
      # collections_query <- stac_source |>
      #   rstac::collections()
      # collections_query
      # 
      # class(stac_source)
      # class(collections_query)
      # 
      # available_collections <- rstac::get_request(collections_query)
      # available_collections
      
      ## setting up STAC search for NASA DEM data
      # rstac::stac_search(
      #   q = stac_source,
      #   collections = "3dep-seamless",
      #   bbox = bbox_temp) |>
      #   rstac::get_request()
      # 
      # 
      # stac_query <- rstac::stac_search(
      #   q = stac_source,
      #   collections = "3dep-seamless",
      #   bbox = bbox_temp)
      # 
      # executed_stac_query <- rstac::get_request(stac_query)
      # executed_stac_query
      
      ### signing in to MPC in order to access STAC with rstac
      # signed_stac_query <- rstac::items_sign(
      #   executed_stac_query,
      #   rstac::sign_planetary_computer())
      # signed_stac_query
      
      ### downloading elevation data
      # rstac::assets_download(signed_stac_query, "data", output_dir = here("Data", "DEM", "Alaska", "USGSTiles"), overwrite = TRUE)
    # }
  # } else if (state == "New Mexico"){
    ### setting bounding box for DEM download
    # bbox_temp <- CONUS_AOI |>
    #   filter(STATE_NAME == state) #|>
      # st_transform(4326) |>
      # st_bbox()
    
    # temp_dem <- get_dem(bbox_temp, pixel_x_size = 30, pixel_y_size = 30, output_filename = here("Data", "DEM", "CONUS", "NewMexico", "NewMexicoDEM.tif"))
    
    # temp_dem <- get_stac_data(bbox_temp, start_date = NULL, end_date = NULL, pixel_x_size = 30, pixel_y_size = 30, asset_names = "elevation", stac_source = "https://planetarycomputer.microsoft.com/api/stac/v1/", collection = "nasadem", output_filename = here("Data", "DEM", "CONUS", "NewMexico", "NewMexicoDEM.tif"))
    
    ### connecting to the appropriate stac API
    # stac_source <- rstac::stac("https://planetarycomputer.microsoft.com/api/stac/v1")
    # stac_source
    # 
    # rstac::get_request(stac_source)
    # 
    # collections_query <- stac_source |>
    #   rstac::collections()
    # collections_query
    # 
    # class(stac_source)
    # class(collections_query)
    # 
    # available_collections <- rstac::get_request(collections_query)
    # available_collections
    
    ## setting up STAC search for Copernicus DEM data
    # rstac::stac_search(
    #   q = stac_source,
    #   collections = "cop-dem-glo-30",
    #   bbox = bbox_temp) |>
    #   rstac::get_request()
    # 
    # 
    # stac_query <- rstac::stac_search(
    #   q = stac_source,
    #   collections = "cop-dem-glo-30",
    #   bbox = bbox_temp)
    # 
    # executed_stac_query <- rstac::get_request(stac_query)
    # executed_stac_query
    
    ### signing in to MPC in order to access STAC with rstac
    # signed_stac_query <- rstac::items_sign(
    #   executed_stac_query,
    #   rstac::sign_planetary_computer())
    # signed_stac_query
    
    ### downloading elevation data
    # rstac::assets_download(signed_stac_query, "data", output_dir = here("Data", "DEM", "CONUS", "NewMexico"), overwrite = TRUE)
  # } else if (state == "South Dakota") {
    ### setting bounding box for DEM download
    # bbox_temp <- CONUS_AOI |>
    #   filter(STATE_NAME == state) #|>
      # st_transform(4326) |>
      # st_bbox()
    
    # temp_dem <- get_dem(bbox_temp, pixel_x_size = 30, pixel_y_size = 30, output_filename = here("Data", "DEM", "CONUS", "SouthDakota", "SouthDakotaDEM.tif"))
    
    # temp_dem <- get_stac_data(bbox_temp, start_date = NULL, end_date = NULL, pixel_x_size = 30, pixel_y_size = 30, asset_names = "elevation", stac_source = "https://planetarycomputer.microsoft.com/api/stac/v1/", collection = "nasadem", output_filename = here("Data", "DEM", "CONUS", "SouthDakota", "SouthDakotaDEM.tif"))
    
    ### connecting to the appropriate stac API
    # stac_source <- rstac::stac("https://planetarycomputer.microsoft.com/api/stac/v1")
    # stac_source
    # 
    # rstac::get_request(stac_source)
    # 
    # collections_query <- stac_source |>
    #   rstac::collections()
    # collections_query
    # 
    # class(stac_source)
    # class(collections_query)
    # 
    # available_collections <- rstac::get_request(collections_query)
    # available_collections
    
    ## setting up STAC search for Copernicus DEM data
    # rstac::stac_search(
    #   q = stac_source,
    #   collections = "cop-dem-glo-30",
    #   bbox = bbox_temp) |>
    #   rstac::get_request()
    # 
    # 
    # stac_query <- rstac::stac_search(
    #   q = stac_source,
    #   collections = "cop-dem-glo-30",
    #   bbox = bbox_temp)
    # 
    # executed_stac_query <- rstac::get_request(stac_query)
    # executed_stac_query
    
    ### signing in to MPC in order to access STAC with rstac
    # signed_stac_query <- rstac::items_sign(
    #   executed_stac_query,
    #   rstac::sign_planetary_computer())
    # signed_stac_query

    ### downloading elevation data
    # rstac::assets_download(signed_stac_query, "data", output_dir = here("Data", "DEM", "CONUS", "SouthDakota"), overwrite = TRUE)
  # } else if (state == "North Dakota"){
    ### setting bounding box for DEM download
    # bbox_temp <- CONUS_AOI |>
    #   filter(STATE_NAME == state) #|>
      # st_transform(4326) |>
      # st_bbox()
    
    # temp_dem <- get_dem(bbox_temp, pixel_x_size = 30, pixel_y_size = 30, output_filename = here("Data", "DEM", "CONUS", "NorthDakota", "NorthDakotaDEM.tif"))
    # temp_dem <- get_stac_data(bbox_temp, start_date = NULL, end_date = NULL, pixel_x_size = 30, pixel_y_size = 30, asset_names = "elevation", stac_source = "https://planetarycomputer.microsoft.com/api/stac/v1/", collection = "nasadem", output_filename = here("Data", "DEM", "CONUS", "NorthDakota", "NorthDakotaDEM.tif"))  
    
      ### connecting to the appropriate stac API
      # stac_source <- rstac::stac("https://planetarycomputer.microsoft.com/api/stac/v1")
      # stac_source
      # 
      # rstac::get_request(stac_source)
      # 
      # collections_query <- stac_source |>
      #   rstac::collections()
      # collections_query
      # 
      # class(stac_source)
      # class(collections_query)
      # 
      # available_collections <- rstac::get_request(collections_query)
      # available_collections
      
      ## setting up STAC search for Copernicus DEM data
      # rstac::stac_search(
      #   q = stac_source,
      #   collections = "cop-dem-glo-30",
      #   bbox = bbox_temp) |>
      #   rstac::get_request()
      # 
      # 
      # stac_query <- rstac::stac_search(
      #   q = stac_source,
      #   collections = "cop-dem-glo-30",
      #   bbox = bbox_temp)
      # 
      # executed_stac_query <- rstac::get_request(stac_query)
      # executed_stac_query
      
      ### signing in to MPC in order to access STAC with rstac
      # signed_stac_query <- rstac::items_sign(
      #   executed_stac_query,
      #   rstac::sign_planetary_computer())
      # signed_stac_query
      
      ### downloading elevation data
      # rstac::assets_download(signed_stac_query, "data", output_dir = here("Data", "DEM", "CONUS", "SouthDakota"), overwrite = TRUE)
  # } else {
  #   bbox_temp <- CONUS_AOI |>
  #     filter(STATE_NAME == state) #|>
      # st_transform(4326) |>
      # st_bbox()
    
    # temp_dem <- get_dem(bbox_temp, pixel_x_size = 30, pixel_y_size = 30, output_filename = here("Data", "DEM", "CONUS", paste0(state), paste0(state, "DEM.tif")))
    
    # temp_dem <- get_stac_data(bbox_temp, start_date = NULL, end_date = NULL, pixel_x_size = 30, pixel_y_size = 30, asset_names = "elevation", stac_source = "https://planetarycomputer.microsoft.com/api/stac/v1/", collection = "nasadem", output_filename = here("Data", "DEM", "CONUS", paste0(state), paste0(state, "DEM.tif")))
    
    ### connecting to the appropriate stac API
    # stac_source <- rstac::stac("https://planetarycomputer.microsoft.com/api/stac/v1")
    # stac_source
    # 
    # rstac::get_request(stac_source)
    # 
    # collections_query <- stac_source |>
    #   rstac::collections()
    # collections_query
    # 
    # class(stac_source)
    # class(collections_query)
    # 
    # available_collections <- rstac::get_request(collections_query)
    # available_collections
    
    ## setting up STAC search for Copernicus DEM data
    # rstac::stac_search(
    #   q = stac_source,
    #   collections = "cop-dem-glo-30",
    #   bbox = bbox_temp) |>
    #   rstac::get_request()
    # 
    # 
    # stac_query <- rstac::stac_search(
    #   q = stac_source,
    #   collections = "cop-dem-glo-30",
    #   bbox = bbox_temp)
    # 
    # executed_stac_query <- rstac::get_request(stac_query)
    # executed_stac_query
    
    ### signing in to MPC in order to access STAC with rstac
    # signed_stac_query <- rstac::items_sign(
    #   executed_stac_query,
    #   rstac::sign_planetary_computer())
    # signed_stac_query
    
    ### downloading elevation data
    # rstac::assets_download(signed_stac_query, "data", output_dir = here("Data", "DEM", "CONUS",  paste0(state)), overwrite = TRUE)
#   }
# }

### DEM Download End Time
# DEMEndTime <- Sys.time()
# DEMDownloadTime <- DEMEndTime - DEMStartTime
# 

##### Land cover download #####
### Land Cover start time
LCStartTime <- Sys.time()
### starting with Alaska
### transforming AOI to WGS84 because rstac seems to require it
# AK_AOI_WGS <- st_transform(AK_AOI, 4326)

### creating bounding box in WGS 84 that rstac expects
# bbox_temp <- st_bbox(AK_AOI_WGS)
# 
# 
# ### connecting to the appropriate stac API
# stac_source <- rstac::stac("https://planetarycomputer.microsoft.com/api/stac/v1")
# stac_source
# 
# rstac::get_request(stac_source)
# 
# collections_query <- stac_source |>
#   rstac::collections()
# collections_query
# 
# class(stac_source)
# class(collections_query)
# 
# available_collections <- rstac::get_request(collections_query)
# available_collections
# 
# ## setting up STAC search for ESA CCI landcover data
# rstac::stac_search(
#   q = stac_source,
#   collections = "esa-cci-lc",
#   bbox = bbox_temp) |>
#   rstac::get_request()
# 
# 
# stac_query <- rstac::stac_search(
#   q = stac_source,
#   collections = "esa-cci-lc",
#   bbox = bbox_temp)
# 
# executed_stac_query <- rstac::get_request(stac_query)
# executed_stac_query
# 
# ### signing in to MPC in order to access STAC with rstac
# signed_stac_query <- rstac::items_sign(
#   executed_stac_query,
#   rstac::sign_planetary_computer())
# signed_stac_query
# 
# ### downloading landcover data
# rstac::assets_download(signed_stac_query, "lccs_class", output_dir = here("Data", "Landcover", "Alaska", "ESACCI"), overwrite = TRUE)


### Now Conus
### transforming AOI to WGS84 because rstac seems to require it
# CONUS_AOI_WGS <- st_transform(CONUS_AOI, 4326)
# 
# ### creating bounding box in WGS 84 that rstac expects
# bbox_temp <- st_bbox(CONUS_AOI_WGS)
# 
# 
# ### connecting to the appropriate stac API
# stac_source <- rstac::stac("https://planetarycomputer.microsoft.com/api/stac/v1")
# stac_source
# 
# rstac::get_request(stac_source)
# 
# collections_query <- stac_source |>
#   rstac::collections()
# collections_query
# 
# class(stac_source)
# class(collections_query)
# 
# available_collections <- rstac::get_request(collections_query)
# available_collections
# 
# ## setting up STAC search for NASA DEM data
# rstac::stac_search(
#   q = stac_source,
#   collections = "esa-cci-lc",
#   bbox = bbox_temp) |>
#   rstac::get_request()
# 
# 
# stac_query <- rstac::stac_search(
#   q = stac_source,
#   collections = "esa-cci-lc",
#   bbox = bbox_temp)
# 
# executed_stac_query <- rstac::get_request(stac_query)
# executed_stac_query
# 
# ### signing in to MPC in order to access STAC with rstac
# signed_stac_query <- rstac::items_sign(
#   executed_stac_query,
#   rstac::sign_planetary_computer())
# signed_stac_query
# 
# ### downloading landcover data
# rstac::assets_download(signed_stac_query, "lccs_class", output_dir = here("Data", "Landcover", "CONUS", "ESACCI"), overwrite = TRUE)


### Land Cover end time
LCEndTime <- Sys.time()
LCDownloadTime <- LCEndTime - LCStartTime

##### End of Script #####
# SnotelDownloadTime
DaymetTileDownloadTime
# DEMDownloadTime
LCDownloadTime
TotalEndTime <- Sys.time()
TotalTime <- TotalEndTime - TotalStartTime
TotalTime
