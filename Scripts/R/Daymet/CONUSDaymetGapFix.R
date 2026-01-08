##### CONUS daymet gap fix #####
TotalStartTime <- Sys.time()
library(pacman)
p_load(here, tidyverse, daymetr, ncdf4, parallel, svMisc, sf, terra)
# future::plan("multisession")


### CONUS next
download_daymet_tiles(
  # location = c(CONUS_minlat, CONUS_minlon, CONUS_maxlat, CONUS_maxlon),
  tiles = c(10834, 10835, 10836),
  start = 1992,
  end = 2023,
  path = here("Data", "Daymet", "CONUS", "tmin_tmax"),
  param = "tmin",
  silent = FALSE,
  force = FALSE
)

download_daymet_tiles(
  # location = c(CONUS_minlat, CONUS_minlon, CONUS_maxlat, CONUS_maxlon),
  tiles = c(10834, 10835, 10836),
  start = 1992,
  end = 2023,
  path = here("Data", "Daymet", "CONUS", "tmin_tmax"),
  param = "tmax",
  silent = FALSE,
  force = FALSE
)

download_daymet_tiles(
  # location = c(CONUS_minlat, CONUS_minlon, CONUS_maxlat, CONUS_maxlon),
  tiles = c(10834, 10835, 10836),
  start = 1992,
  end = 2023,
  path = here("Data", "Daymet", "CONUS", "srad"),
  param = "srad",
  silent = FALSE,
  force = FALSE
)

download_daymet_tiles(
  # location = c(CONUS_minlat, CONUS_minlon, CONUS_maxlat, CONUS_maxlon),
  tiles = c(10834, 10835, 10836),
  start = 1992,
  end = 2023,
  path = here("Data", "Daymet", "CONUS", "prcp"),
  param = "prcp",
  silent = FALSE,
  force = FALSE
)



