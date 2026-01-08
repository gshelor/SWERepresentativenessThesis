##### Reclassifying the ESA Landcover #####

### loading packages, reading in AOIs
library(pacman)
p_load(here, tidyverse, terra, sf)
CONUS_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "CONUS", "CONUS_AOI.gpkg"))
AK_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "Alaska", "AK_AOI.gpkg"))

##### Creating Reclassification Matrix #####
### 0 = "No Vegetation"
### 1 = "Short Vegetation"
### 2 = "Tall Vegetation"
ReclassMatrix <- matrix(data = c(0, 0,
                                 10, 1,
                                 11, 1,
                                 12, 1,
                                 20, 1,
                                 30, 1,
                                 40, 1,
                                 50, 2,
                                 60, 2,
                                 61, 2,
                                 62, 1,
                                 70, 2,
                                 71, 2,
                                 72, 1,
                                 80, 2,
                                 81, 2,
                                 82, 1,
                                 90, 2,
                                 100, 2,
                                 110, 1,
                                 120, 1,
                                 121, 1,
                                 122, 1,
                                 130, 1,
                                 140, 1,
                                 150, 1,
                                 151, 1,
                                 152, 1,
                                 153, 1,
                                 160, 1,
                                 170, 1,
                                 180, 1,
                                 190, 0,
                                 200, 0,
                                 201, 0,
                                 202, 0,
                                 210, 0,
                                 220, 0), ncol = 2, byrow = TRUE)

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
lc_CONUS_rast_list <- list(LC_CONUS_1992, LC_CONUS_1993, LC_CONUS_1994, LC_CONUS_1995, LC_CONUS_1996, LC_CONUS_1997, LC_CONUS_1998, LC_CONUS_1999, LC_CONUS_2000, LC_CONUS_2001, LC_CONUS_2002, LC_CONUS_2003, LC_CONUS_2004, LC_CONUS_2005, LC_CONUS_2006, LC_CONUS_2007, LC_CONUS_2008, LC_CONUS_2009, LC_CONUS_2010, LC_CONUS_2011, LC_CONUS_2012, LC_CONUS_2013, LC_CONUS_2014, LC_CONUS_2015, LC_CONUS_2016, LC_CONUS_2017, LC_CONUS_2018, LC_CONUS_2019, LC_CONUS_2020)

years <- 1992:2020
for (x in 1:length(years)){
  print(paste("CONUS:", years[x]))
  temp_lc <- lc_CONUS_rast_list[[x]]
  if (ext(temp_lc) != ext(LC_CONUS_2020)){
    print(paste(years[x], ext(temp_lc)))
  }
  temp_lc <- project(temp_lc, crs(CONUS_AOI), method = "near")
  temp_reclassified <- classify(temp_lc, rcl = ReclassMatrix)
  writeRaster(temp_reclassified, here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", paste0("LC_CONUS_", years[x], ".tif")), overwrite = TRUE)
}


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

### For some reason, the ArcPro modelbuilder tool I made spits out the 1997 landcover with a different extent than all of the others. I have looked at everything I can think of and I see no reason for this anywhere in the tool
## I figured it out, I left out a zero in the reproject tool in my ArcPro model. This is fixed as of March 29, 2025
## so I resample it to match 2020
# LC_AK_1997 <- resample(LC_AK_1997, LC_AK_2020, method = "near")
# writeRaster(LC_AK_1997, here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "LC_AK_1997.tif"), overwrite = TRUE)

### making a list of landcover rasters for when landcover classes get extracted
lc_AK_rast_list <- list(LC_AK_1992, LC_AK_1993, LC_AK_1994, LC_AK_1995, LC_AK_1996, LC_AK_1997, LC_AK_1998, LC_AK_1999, LC_AK_2000, LC_AK_2001, LC_AK_2002, LC_AK_2003, LC_AK_2004, LC_AK_2005, LC_AK_2006, LC_AK_2007, LC_AK_2008, LC_AK_2009, LC_AK_2010, LC_AK_2011, LC_AK_2012, LC_AK_2013, LC_AK_2014, LC_AK_2015, LC_AK_2016, LC_AK_2017, LC_AK_2018, LC_AK_2019, LC_AK_2020)


for (x in 1:length(years)){
  print(paste("Alaska:", years[x]))
  temp_lc <- lc_AK_rast_list[[x]]
  if (ext(temp_lc) != ext(LC_AK_2020)){
    print(paste(years[x], ext(temp_lc)))
  }
  temp_lc <- project(temp_lc, crs(AK_AOI), method = "near")
  temp_reclassified <- classify(temp_lc, rcl = ReclassMatrix)
  writeRaster(temp_reclassified, here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", paste0("LC_AK_", years[x], ".tif")), overwrite = TRUE)
}


##### Creating Climatologies from Reclassified Landcover rasters #####
### CONUS
LC_CONUS_1992 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_1992.tif"))
LC_CONUS_1993 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_1993.tif"))
LC_CONUS_1994 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_1994.tif"))
LC_CONUS_1995 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_1995.tif"))
LC_CONUS_1996 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_1996.tif"))
LC_CONUS_1997 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_1997.tif"))
LC_CONUS_1998 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_1998.tif"))
LC_CONUS_1999 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_1999.tif"))
LC_CONUS_2000 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2000.tif"))
LC_CONUS_2001 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2001.tif"))
LC_CONUS_2002 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2002.tif"))
LC_CONUS_2003 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2003.tif"))
LC_CONUS_2004 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2004.tif"))
LC_CONUS_2005 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2005.tif"))
LC_CONUS_2006 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2006.tif"))
LC_CONUS_2007 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2007.tif"))
LC_CONUS_2008 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2008.tif"))
LC_CONUS_2009 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2009.tif"))
LC_CONUS_2010 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2010.tif"))
LC_CONUS_2011 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2011.tif"))
LC_CONUS_2012 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2012.tif"))
LC_CONUS_2013 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2013.tif"))
LC_CONUS_2014 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2014.tif"))
LC_CONUS_2015 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2015.tif"))
LC_CONUS_2016 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2016.tif"))
LC_CONUS_2017 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2017.tif"))
LC_CONUS_2018 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2018.tif"))
LC_CONUS_2019 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2019.tif"))
LC_CONUS_2020 <- rast(here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "LC_CONUS_2020.tif"))

### making a list of landcover rasters for when landcover classes get extracted
lc_CONUS_rasts <- c(LC_CONUS_1992, LC_CONUS_1993, LC_CONUS_1994, LC_CONUS_1995, LC_CONUS_1996, LC_CONUS_1997, LC_CONUS_1998, LC_CONUS_1999, LC_CONUS_2000, LC_CONUS_2001, LC_CONUS_2002, LC_CONUS_2003, LC_CONUS_2004, LC_CONUS_2005, LC_CONUS_2006, LC_CONUS_2007, LC_CONUS_2008, LC_CONUS_2009, LC_CONUS_2010, LC_CONUS_2011, LC_CONUS_2012, LC_CONUS_2013, LC_CONUS_2014, LC_CONUS_2015, LC_CONUS_2016, LC_CONUS_2017, LC_CONUS_2018, LC_CONUS_2019, LC_CONUS_2020)
### creating climatology
lc_CONUS_climatology <- app(lc_CONUS_rasts, fun = "modal")
writeRaster(lc_CONUS_climatology, here("Data", "Landcover", "CONUS", "ESACCI", "Clipped", "Reclassified", "Climatology", "LC_CONUS_Climatology_Reclass.tif"))


### Alaska
LC_AK_1992 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_1992.tif"))
LC_AK_1993 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_1993.tif"))
LC_AK_1994 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_1994.tif"))
LC_AK_1995 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_1995.tif"))
LC_AK_1996 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_1996.tif"))
LC_AK_1997 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_1997.tif"))
LC_AK_1998 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_1998.tif"))
LC_AK_1999 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_1999.tif"))
LC_AK_2000 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2000.tif"))
LC_AK_2001 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2001.tif"))
LC_AK_2002 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2002.tif"))
LC_AK_2003 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2003.tif"))
LC_AK_2004 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2004.tif"))
LC_AK_2005 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2005.tif"))
LC_AK_2006 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2006.tif"))
LC_AK_2007 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2007.tif"))
LC_AK_2008 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2008.tif"))
LC_AK_2009 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2009.tif"))
LC_AK_2010 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2010.tif"))
LC_AK_2011 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2011.tif"))
LC_AK_2012 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2012.tif"))
LC_AK_2013 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2013.tif"))
LC_AK_2014 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2014.tif"))
LC_AK_2015 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2015.tif"))
LC_AK_2016 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2016.tif"))
LC_AK_2017 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2017.tif"))
LC_AK_2018 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2018.tif"))
LC_AK_2019 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2019.tif"))
LC_AK_2020 <- rast(here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "LC_AK_2020.tif"))


lc_AK_rasts <- c(LC_AK_1992, LC_AK_1993, LC_AK_1994, LC_AK_1995, LC_AK_1996, LC_AK_1997, LC_AK_1998, LC_AK_1999, LC_AK_2000, LC_AK_2001, LC_AK_2002, LC_AK_2003, LC_AK_2004, LC_AK_2005, LC_AK_2006, LC_AK_2007, LC_AK_2008, LC_AK_2009, LC_AK_2010, LC_AK_2011, LC_AK_2012, LC_AK_2013, LC_AK_2014, LC_AK_2015, LC_AK_2016, LC_AK_2017, LC_AK_2018, LC_AK_2019, LC_AK_2020)

lc_AK_climatology <- app(lc_AK_rasts, fun = "modal")
writeRaster(lc_AK_climatology, here("Data", "Landcover", "Alaska", "ESACCI", "Clipped", "Reclassified", "Climatology", "LC_AK_Climatology_Reclass.tif"), overwrite = TRUE)


plot(lc_CONUS_climatology)
plot(lc_AK_climatology)
