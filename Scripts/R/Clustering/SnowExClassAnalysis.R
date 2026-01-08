##### Analyzing SnowEx Site Presence in SWE-based classes #####

##### loading packages, reading in data #####
library(pacman)
p_load(here, tidyverse, sf, terra, gt, gtExtras)

### reading in AOIs
AK_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "Alaska", "AK_AOI.gpkg"))
CONUS_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "CONUS", "CONUSEcoregionsNoStates.gpkg"))

### reading in SnowEx site csv
SnowEx_df <- read_csv(here("Data", "SnowEx", "SnowExSiteLocations.csv"))

### filtering by study area
SnowExAK_df <- SnowEx_df |>
  filter(state == "AK")
SnowExCONUS_df <- SnowEx_df |>
  filter(state != "AK")

### converting to sf objects
SnowExAK_sf <- st_as_sf(SnowExAK_df, coords = c("lon", "lat"), crs = 4326) |>
  st_transform(crs = crs(AK_AOI))
SnowExCONUS_sf <- st_as_sf(SnowExCONUS_df, coords = c("lon", "lat"), crs = 4326) |>
  st_transform(crs = crs(CONUS_AOI))


### reading in class rasters to calculate most common class and number of unique classes
CONUS1993Classes_rast <- rast(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "CONUS1993_KClasses.tif"))
CONUS1994Classes_rast <- rast(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "CONUS1994_KClasses.tif"))
CONUS1995Classes_rast <- rast(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "CONUS1995_KClasses.tif"))
CONUS1996Classes_rast <- rast(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "CONUS1996_KClasses.tif"))
CONUS1997Classes_rast <- rast(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "CONUS1997_KClasses.tif"))
CONUS1998Classes_rast <- rast(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "CONUS1998_KClasses.tif"))
CONUS1999Classes_rast <- rast(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "CONUS1999_KClasses.tif"))
CONUS2000Classes_rast <- rast(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "CONUS2000_KClasses.tif"))
CONUS2001Classes_rast <- rast(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "CONUS2001_KClasses.tif"))
CONUS2002Classes_rast <- rast(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "CONUS2002_KClasses.tif"))
CONUS2003Classes_rast <- rast(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "CONUS2003_KClasses.tif"))
CONUS2004Classes_rast <- rast(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "CONUS2004_KClasses.tif"))
CONUS2005Classes_rast <- rast(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "CONUS2005_KClasses.tif"))
CONUS2006Classes_rast <- rast(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "CONUS2006_KClasses.tif"))
CONUS2007Classes_rast <- rast(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "CONUS2007_KClasses.tif"))
CONUS2008Classes_rast <- rast(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "CONUS2008_KClasses.tif"))
CONUS2009Classes_rast <- rast(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "CONUS2009_KClasses.tif"))
CONUS2010Classes_rast <- rast(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "CONUS2010_KClasses.tif"))
CONUS2011Classes_rast <- rast(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "CONUS2011_KClasses.tif"))
CONUS2012Classes_rast <- rast(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "CONUS2012_KClasses.tif"))
CONUS2013Classes_rast <- rast(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "CONUS2013_KClasses.tif"))
CONUS2014Classes_rast <- rast(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "CONUS2014_KClasses.tif"))
CONUS2015Classes_rast <- rast(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "CONUS2015_KClasses.tif"))
CONUS2016Classes_rast <- rast(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "CONUS2016_KClasses.tif"))
CONUS2017Classes_rast <- rast(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "CONUS2017_KClasses.tif"))
CONUS2018Classes_rast <- rast(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "CONUS2018_KClasses.tif"))
CONUS2019Classes_rast <- rast(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "CONUS2019_KClasses.tif"))
CONUS2020Classes_rast <- rast(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "CONUS2020_KClasses.tif"))

### reading in class rasters
AK1993Classes_rast <- rast(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "AK1993_KClasses.tif"))
AK1994Classes_rast <- rast(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "AK1994_KClasses.tif"))
AK1995Classes_rast <- rast(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "AK1995_KClasses.tif"))
AK1996Classes_rast <- rast(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "AK1996_KClasses.tif"))
AK1997Classes_rast <- rast(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "AK1997_KClasses.tif"))
AK1998Classes_rast <- rast(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "AK1998_KClasses.tif"))
AK1999Classes_rast <- rast(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "AK1999_KClasses.tif"))
AK2000Classes_rast <- rast(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "AK2000_KClasses.tif"))
AK2001Classes_rast <- rast(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "AK2001_KClasses.tif"))
AK2002Classes_rast <- rast(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "AK2002_KClasses.tif"))
AK2003Classes_rast <- rast(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "AK2003_KClasses.tif"))
AK2004Classes_rast <- rast(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "AK2004_KClasses.tif"))
AK2005Classes_rast <- rast(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "AK2005_KClasses.tif"))
AK2006Classes_rast <- rast(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "AK2006_KClasses.tif"))
AK2007Classes_rast <- rast(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "AK2007_KClasses.tif"))
AK2008Classes_rast <- rast(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "AK2008_KClasses.tif"))
AK2009Classes_rast <- rast(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "AK2009_KClasses.tif"))
AK2010Classes_rast <- rast(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "AK2010_KClasses.tif"))
AK2011Classes_rast <- rast(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "AK2011_KClasses.tif"))
AK2012Classes_rast <- rast(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "AK2012_KClasses.tif"))
AK2013Classes_rast <- rast(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "AK2013_KClasses.tif"))
AK2014Classes_rast <- rast(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "AK2014_KClasses.tif"))
AK2015Classes_rast <- rast(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "AK2015_KClasses.tif"))
AK2016Classes_rast <- rast(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "AK2016_KClasses.tif"))
AK2017Classes_rast <- rast(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "AK2017_KClasses.tif"))
AK2018Classes_rast <- rast(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "AK2018_KClasses.tif"))
AK2019Classes_rast <- rast(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "AK2019_KClasses.tif"))
AK2020Classes_rast <- rast(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "AK2020_KClasses.tif"))


### Reading in Modal class raster and number of unique classes raster
CONUSModalClass_rast <- rast(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "Summary", "CONUSModalClass.tif"))
plot(CONUSModalClass_rast)
plot(st_geometry(SnowExCONUS_sf), add = TRUE)
AKModalClass_rast <- rast(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "Summary", "AKModalClass.tif"))
plot(AKModalClass_rast)
plot(st_geometry(SnowExAK_sf), add = TRUE)
CONUSNumUniqueClass_rast <- rast(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "Summary", "CONUSNumUniqueClass.tif"))
AKNumUniqueClass_rast <- rast(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "Summary", "AKNumUniqueClass.tif"))

### reading in Sturm classes
SturmListonCONUS_rast <- rast(here("Data", "SturmListon2021Classes", "Resampled", "CONUSSturmClasses.tif"))
plot(SturmListonCONUS_rast, main = "Sturm WUS classes")
plot(st_geometry(SnowExCONUS_sf), add = TRUE)
SturmListonAK_rast <- rast(here("Data", "SturmListon2021Classes", "Resampled", "AKSturmClasses.tif"))
plot(SturmListonAK_rast, main = "Sturm AK classes")
plot(st_geometry(SnowExAK_sf), add = TRUE)

##### Extracting Classes to SnowEx Points #####
### Alaska first
SnowExAK_sf$Class1993 <- terra::extract(AK1993Classes_rast, SnowExAK_sf)[,2]
SnowExAK_sf$Class1994 <- terra::extract(AK1994Classes_rast, SnowExAK_sf)[,2]
SnowExAK_sf$Class1995 <- terra::extract(AK1995Classes_rast, SnowExAK_sf)[,2]
SnowExAK_sf$Class1996 <- terra::extract(AK1996Classes_rast, SnowExAK_sf)[,2]
SnowExAK_sf$Class1997 <- terra::extract(AK1997Classes_rast, SnowExAK_sf)[,2]
SnowExAK_sf$Class1998 <- terra::extract(AK1998Classes_rast, SnowExAK_sf)[,2]
SnowExAK_sf$Class1999 <- terra::extract(AK1999Classes_rast, SnowExAK_sf)[,2]
SnowExAK_sf$Class2000 <- terra::extract(AK2000Classes_rast, SnowExAK_sf)[,2]
SnowExAK_sf$Class2001 <- terra::extract(AK2001Classes_rast, SnowExAK_sf)[,2]
SnowExAK_sf$Class2002 <- terra::extract(AK2002Classes_rast, SnowExAK_sf)[,2]
SnowExAK_sf$Class2003 <- terra::extract(AK2003Classes_rast, SnowExAK_sf)[,2]
SnowExAK_sf$Class2004 <- terra::extract(AK2004Classes_rast, SnowExAK_sf)[,2]
SnowExAK_sf$Class2005 <- terra::extract(AK2005Classes_rast, SnowExAK_sf)[,2]
SnowExAK_sf$Class2006 <- terra::extract(AK2006Classes_rast, SnowExAK_sf)[,2]
SnowExAK_sf$Class2007 <- terra::extract(AK2007Classes_rast, SnowExAK_sf)[,2]
SnowExAK_sf$Class2008 <- terra::extract(AK2008Classes_rast, SnowExAK_sf)[,2]
SnowExAK_sf$Class2009 <- terra::extract(AK2009Classes_rast, SnowExAK_sf)[,2]
SnowExAK_sf$Class2010 <- terra::extract(AK2010Classes_rast, SnowExAK_sf)[,2]
SnowExAK_sf$Class2011 <- terra::extract(AK2011Classes_rast, SnowExAK_sf)[,2]
SnowExAK_sf$Class2012 <- terra::extract(AK2012Classes_rast, SnowExAK_sf)[,2]
SnowExAK_sf$Class2013 <- terra::extract(AK2013Classes_rast, SnowExAK_sf)[,2]
SnowExAK_sf$Class2014 <- terra::extract(AK2014Classes_rast, SnowExAK_sf)[,2]
SnowExAK_sf$Class2015 <- terra::extract(AK2015Classes_rast, SnowExAK_sf)[,2]
SnowExAK_sf$Class2016 <- terra::extract(AK2016Classes_rast, SnowExAK_sf)[,2]
SnowExAK_sf$Class2017 <- terra::extract(AK2017Classes_rast, SnowExAK_sf)[,2]
SnowExAK_sf$Class2018 <- terra::extract(AK2018Classes_rast, SnowExAK_sf)[,2]
SnowExAK_sf$Class2019 <- terra::extract(AK2019Classes_rast, SnowExAK_sf)[,2]
SnowExAK_sf$Class2020 <- terra::extract(AK2020Classes_rast, SnowExAK_sf)[,2]
SnowExAK_sf$ModalClass <- terra::extract(AKModalClass_rast, SnowExAK_sf)[,2]
SnowExAK_sf$NumUniqueClasses <- terra::extract(AKNumUniqueClass_rast, SnowExAK_sf)[,2]
SnowExAK_sf$SturmClass <- terra::extract(SturmListonAK_rast, SnowExAK_sf)[,2]

### sturm classes read into R as numbers which is incredibly frustrating, so I figured out which number is which class, yippee
SnowExAK_sf <- SnowExAK_sf |>
  mutate(SturmClass_string = case_when(SturmClass == 1 ~ "Tundra",
                                       SturmClass == 2 ~ "Boreal Forest",
                                       SturmClass == 3 ~ "Maritime",
                                       SturmClass == 4 ~ "Ephemeral (or no snow)",
                                       SturmClass == 5 ~ "Prairie",
                                       SturmClass == 6 ~ "Montane Forest",
                                       SturmClass == 7 ~ "Ice",
                                       SturmClass == 8 ~ "Ocean/Water"))

### CONUS first
SnowExCONUS_sf$Class1993 <- terra::extract(CONUS1993Classes_rast, SnowExCONUS_sf)[,2]
SnowExCONUS_sf$Class1994 <- terra::extract(CONUS1994Classes_rast, SnowExCONUS_sf)[,2]
SnowExCONUS_sf$Class1995 <- terra::extract(CONUS1995Classes_rast, SnowExCONUS_sf)[,2]
SnowExCONUS_sf$Class1996 <- terra::extract(CONUS1996Classes_rast, SnowExCONUS_sf)[,2]
SnowExCONUS_sf$Class1997 <- terra::extract(CONUS1997Classes_rast, SnowExCONUS_sf)[,2]
SnowExCONUS_sf$Class1998 <- terra::extract(CONUS1998Classes_rast, SnowExCONUS_sf)[,2]
SnowExCONUS_sf$Class1999 <- terra::extract(CONUS1999Classes_rast, SnowExCONUS_sf)[,2]
SnowExCONUS_sf$Class2000 <- terra::extract(CONUS2000Classes_rast, SnowExCONUS_sf)[,2]
SnowExCONUS_sf$Class2001 <- terra::extract(CONUS2001Classes_rast, SnowExCONUS_sf)[,2]
SnowExCONUS_sf$Class2002 <- terra::extract(CONUS2002Classes_rast, SnowExCONUS_sf)[,2]
SnowExCONUS_sf$Class2003 <- terra::extract(CONUS2003Classes_rast, SnowExCONUS_sf)[,2]
SnowExCONUS_sf$Class2004 <- terra::extract(CONUS2004Classes_rast, SnowExCONUS_sf)[,2]
SnowExCONUS_sf$Class2005 <- terra::extract(CONUS2005Classes_rast, SnowExCONUS_sf)[,2]
SnowExCONUS_sf$Class2006 <- terra::extract(CONUS2006Classes_rast, SnowExCONUS_sf)[,2]
SnowExCONUS_sf$Class2007 <- terra::extract(CONUS2007Classes_rast, SnowExCONUS_sf)[,2]
SnowExCONUS_sf$Class2008 <- terra::extract(CONUS2008Classes_rast, SnowExCONUS_sf)[,2]
SnowExCONUS_sf$Class2009 <- terra::extract(CONUS2009Classes_rast, SnowExCONUS_sf)[,2]
SnowExCONUS_sf$Class2010 <- terra::extract(CONUS2010Classes_rast, SnowExCONUS_sf)[,2]
SnowExCONUS_sf$Class2011 <- terra::extract(CONUS2011Classes_rast, SnowExCONUS_sf)[,2]
SnowExCONUS_sf$Class2012 <- terra::extract(CONUS2012Classes_rast, SnowExCONUS_sf)[,2]
SnowExCONUS_sf$Class2013 <- terra::extract(CONUS2013Classes_rast, SnowExCONUS_sf)[,2]
SnowExCONUS_sf$Class2014 <- terra::extract(CONUS2014Classes_rast, SnowExCONUS_sf)[,2]
SnowExCONUS_sf$Class2015 <- terra::extract(CONUS2015Classes_rast, SnowExCONUS_sf)[,2]
SnowExCONUS_sf$Class2016 <- terra::extract(CONUS2016Classes_rast, SnowExCONUS_sf)[,2]
SnowExCONUS_sf$Class2017 <- terra::extract(CONUS2017Classes_rast, SnowExCONUS_sf)[,2]
SnowExCONUS_sf$Class2018 <- terra::extract(CONUS2018Classes_rast, SnowExCONUS_sf)[,2]
SnowExCONUS_sf$Class2019 <- terra::extract(CONUS2019Classes_rast, SnowExCONUS_sf)[,2]
SnowExCONUS_sf$Class2020 <- terra::extract(CONUS2020Classes_rast, SnowExCONUS_sf)[,2]
SnowExCONUS_sf$ModalClass <- terra::extract(CONUSModalClass_rast, SnowExCONUS_sf)[,2]
SnowExCONUS_sf$NumUniqueClasses <- terra::extract(CONUSNumUniqueClass_rast, SnowExCONUS_sf)[,2]
SnowExCONUS_sf$SturmClass <- terra::extract(SturmListonCONUS_rast, SnowExCONUS_sf)[,2]

### sturm classes read into R as numbers which is incredibly frustrating, so I figured out which number is which class, yippee
SnowExCONUS_sf <- SnowExCONUS_sf |>
  mutate(SturmClass_string = case_when(SturmClass == 1 ~ "Tundra",
                                       SturmClass == 2 ~ "Boreal Forest",
                                       SturmClass == 3 ~ "Maritime",
                                       SturmClass == 4 ~ "Ephemeral (or no snow)",
                                       SturmClass == 5 ~ "Prairie",
                                       SturmClass == 6 ~ "Montane Forest",
                                       SturmClass == 7 ~ "Ice",
                                       SturmClass == 8 ~ "Ocean/Water"))


##### writing out SnowEx sf objects #####
write_sf(SnowExCONUS_sf, here("Data", "SnowEx", "SnowExCONUS.gpkg"))
write_sf(SnowExAK_sf, here("Data", "SnowEx", "SnowExAK.gpkg"))
write_sf(SnowExCONUS_sf, here("Data", "SnowEx", "SHP", "SnowExCONUS.shp"))
write_sf(SnowExAK_sf, here("Data", "SnowEx", "SHP", "SnowExAK.shp"))

##### Examining how many SnowEx sites are in each class #####
### figuring out how many times different classes are observed in SnowEx sites
SnowExAK_ClassCounts <- SnowExAK_sf |> 
  count(ModalClass)
SnowExCONUS_ClassCounts <- SnowExCONUS_sf |> 
  count(ModalClass)
### reading in class average csvs to use as dataframes to aggregate by
AKClassesAvg <- read_csv(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "CSVs", "Summary", "AKClassSHAPAvgs.csv")) |>
  select(class_cluster, mean_swe_max, median_swe_max) |>
  mutate(SnowExSites = numeric(nrow(AKClassesAvg)))
  # mutate(SnowExSites = case_when(class_cluster %in% SnowExAK_ClassCounts$ModalClass ~ SnowExAK_ClassCounts$n[which(SnowExAK_ClassCounts$ModalClass == class_cluster)],
  #                                TRUE ~ 0))
### the tidyverse is failing me, so
for (i in 1:nrow(AKClassesAvg)){
  temp_cluster_sites <- sum(SnowExAK_ClassCounts$n[SnowExAK_ClassCounts$ModalClass == i])
  AKClassesAvg$SnowExSites[i] = temp_cluster_sites
}
### CONUS
CONUSClassesAvg <- read_csv(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "CSVs", "Summary", "CONUSClassSHAPAvgs.csv")) |>
  select(class_cluster, mean_swe_max, median_swe_max) |>
  mutate(SnowExSites = numeric(nrow(CONUSClassesAvg)))
# mutate(SnowExSites = case_when(class_cluster %in% SnowExCONUS_ClassCounts$ModalClass ~ SnowExCONUS_ClassCounts$n[which(SnowExCONUS_ClassCounts$ModalClass == class_cluster)],
#                                TRUE ~ 0))
### the tidyverse is failing me, so
for (i in 1:nrow(CONUSClassesAvg)){
  temp_cluster_sites <- sum(SnowExCONUS_ClassCounts$n[SnowExCONUS_ClassCounts$ModalClass == i])
  CONUSClassesAvg$SnowExSites[i] = temp_cluster_sites
}


##### Making Tables showing class number, median peak swe, and number of SnowEx sites #####
### Alaska
AKClasses_SnowEx_gt <- AKClassesAvg |>
  gt() |> # use 'gt' to make an awesome table...
  gt_theme_guardian() |>
  # tab_header(
  #   title = paste(), # ...with this title
  #   subtitle = "")  |>  # and this subtitle
  ##tab_style(style = cell_fill("bisque"),
  ##        locations = cells_body()) |>  # add fill color to table
  fmt_number( # A column (numeric data)
    columns = c(median_swe_max), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # A column (numeric data)
    columns = c(SnowExSites), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 0 # With four decimal places
  ) |>
  data_color( # Update cell colors, testing different color palettes
    columns = c(median_swe_max),
    fn = scales::col_numeric( # <- bc it's numeric
      palette = brewer.pal(11, "RdBu"), # A color scheme (gradient)
      domain = c(), # Column scale endpoints
      reverse = FALSE
    )
  ) |>
  data_color( # Update cell colors, testing different color palettes
    columns = c(SnowExSites), # ...for dose column
    # direction = "row",
    fn = scales::col_numeric( # <- bc it's numeric
      palette = brewer.pal(9, "Blues"), # A color scheme (gradient)
      domain = c(), # Column scale endpoints
      reverse = FALSE
    )
  ) |>
  cols_label(class_cluster = "Class", median_swe_max = "Median Max SWE (mm)", SnowExSites = "# of SnowEx Sites") |> # Update labels
  # cols_move_to_end(columns = c("mean_AbsSHAP_DecFeb_prcpSumCDMSum", "mean_AbsSHAP_DecFeb_tmeanCDMSum", "mean_AbsSHAP_elevation", "mean_AbsSHAP_slope", "mean_AbsSHAP_aspect", "mean_AbsSHAP_landcover_triclass")) |>
  cols_hide(c(mean_swe_max)) |>
  opt_vertical_padding(scale = 1)
# tab_footnote(
#   footnote = ""
# )
AKClasses_SnowEx_gt
### saving table
AKClasses_SnowEx_gt |>
  gtsave(
    "AKClusterSnowExSiteCounts.png", expand = 5,
    path = here("Outputs", "ChartsGraphsTables", "ClusterSnowEx")
  )

### now an alaska one with a title
AKClasses_SnowExTitle_gt <- AKClassesAvg |>
  gt() |> # use 'gt' to make an awesome table...
  gt_theme_guardian() |>
  tab_header(
    title = "Number of SnowEx Sites in Each Class", # ...with this title
    subtitle = "Alaska")  |>  # and this subtitle
  opt_align_table_header(align = "center") |>
  ##tab_style(style = cell_fill("bisque"),
  ##        locations = cells_body()) |>  # add fill color to table
  fmt_number( # A column (numeric data)
    columns = c(median_swe_max), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # A column (numeric data)
    columns = c(SnowExSites), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 0 # With four decimal places
  ) |>
  data_color( # Update cell colors, testing different color palettes
    columns = c(median_swe_max),
    fn = scales::col_numeric( # <- bc it's numeric
      palette = brewer.pal(11, "RdBu"), # A color scheme (gradient)
      domain = c(), # Column scale endpoints
      reverse = FALSE
    )
  ) |>
  data_color( # Update cell colors, testing different color palettes
    columns = c(SnowExSites), # ...for dose column
    # direction = "row",
    fn = scales::col_numeric( # <- bc it's numeric
      palette = brewer.pal(9, "Blues"), # A color scheme (gradient)
      domain = c(), # Column scale endpoints
      reverse = FALSE
    )
  ) |>
  cols_label(class_cluster = "Class", median_swe_max = "Median Max SWE (mm)", SnowExSites = "# of SnowEx Sites") |> # Update labels
  # cols_move_to_end(columns = c("mean_AbsSHAP_DecFeb_prcpSumCDMSum", "mean_AbsSHAP_DecFeb_tmeanCDMSum", "mean_AbsSHAP_elevation", "mean_AbsSHAP_slope", "mean_AbsSHAP_aspect", "mean_AbsSHAP_landcover_triclass")) |>
  cols_hide(c(mean_swe_max)) |>
  opt_vertical_padding(scale = 1)
# tab_footnote(
#   footnote = ""
# )
AKClasses_SnowExTitle_gt
### saving table
AKClasses_SnowExTitle_gt |>
  gtsave(
    "AKClusterSnowExSiteCounts_Title.png", expand = 5,
    path = here("Outputs", "ChartsGraphsTables", "ClusterSnowEx")
  )


### CONUS
CONUSClasses_SnowEx_gt <- CONUSClassesAvg |>
  gt() |> # use 'gt' to make an awesome table...
  gt_theme_guardian() |>
  # tab_header(
  #   title = paste(), # ...with this title
  #   subtitle = "")  |>  # and this subtitle
  ##tab_style(style = cell_fill("bisque"),
  ##        locations = cells_body()) |>  # add fill color to table
  fmt_number( # A column (numeric data)
    columns = c(median_swe_max), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # A column (numeric data)
    columns = c(SnowExSites), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 0 # With four decimal places
  ) |>
  data_color( # Update cell colors, testing different color palettes
    columns = c(median_swe_max),
    fn = scales::col_numeric( # <- bc it's numeric
      palette = brewer.pal(11, "RdBu"), # A color scheme (gradient)
      domain = c(), # Column scale endpoints
      reverse = FALSE
    )
  ) |>
  data_color( # Update cell colors, testing different color palettes
    columns = c(SnowExSites), # ...for dose column
    # direction = "row",
    fn = scales::col_numeric( # <- bc it's numeric
      palette = brewer.pal(9, "Blues"), # A color scheme (gradient)
      domain = c(), # Column scale endpoints
      reverse = FALSE
    )
  ) |>
  cols_label(class_cluster = "Class", median_swe_max = "Median Max SWE (mm)", SnowExSites = "# of SnowEx Sites") |> # Update labels
  # cols_move_to_end(columns = c("mean_AbsSHAP_DecFeb_prcpSumCDMSum", "mean_AbsSHAP_DecFeb_tmeanCDMSum", "mean_AbsSHAP_elevation", "mean_AbsSHAP_slope", "mean_AbsSHAP_aspect", "mean_AbsSHAP_landcover_triclass")) |>
  cols_hide(c(mean_swe_max)) |>
  opt_vertical_padding(scale = 1)
# tab_footnote(
#   footnote = ""
# )
CONUSClasses_SnowEx_gt
### saving table
CONUSClasses_SnowEx_gt |>
  gtsave(
    "CONUSClusterSnowExSiteCounts.png", expand = 5,
    path = here("Outputs", "ChartsGraphsTables", "ClusterSnowEx")
  )

### now an alaska one with a title
CONUSClasses_SnowExTitle_gt <- CONUSClassesAvg |>
  gt() |> # use 'gt' to make an awesome table...
  gt_theme_guardian() |>
  tab_header(
    title = "Number of SnowEx Sites in Each Class", # ...with this title
    subtitle = "Western US")  |>  # and this subtitle
  opt_align_table_header(align = "center") |>
  ##tab_style(style = cell_fill("bisque"),
  ##        locations = cells_body()) |>  # add fill color to table
  fmt_number( # A column (numeric data)
    columns = c(median_swe_max), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # A column (numeric data)
    columns = c(SnowExSites), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 0 # With four decimal places
  ) |>
  data_color( # Update cell colors, testing different color palettes
    columns = c(median_swe_max),
    fn = scales::col_numeric( # <- bc it's numeric
      palette = brewer.pal(11, "RdBu"), # A color scheme (gradient)
      domain = c(), # Column scale endpoints
      reverse = FALSE
    )
  ) |>
  data_color( # Update cell colors, testing different color palettes
    columns = c(SnowExSites), # ...for dose column
    # direction = "row",
    fn = scales::col_numeric( # <- bc it's numeric
      palette = brewer.pal(9, "Blues"), # A color scheme (gradient)
      domain = c(), # Column scale endpoints
      reverse = FALSE
    )
  ) |>
  cols_label(class_cluster = "Class", median_swe_max = "Median Max SWE (mm)", SnowExSites = "# of SnowEx Sites") |> # Update labels
  # cols_move_to_end(columns = c("mean_AbsSHAP_DecFeb_prcpSumCDMSum", "mean_AbsSHAP_DecFeb_tmeanCDMSum", "mean_AbsSHAP_elevation", "mean_AbsSHAP_slope", "mean_AbsSHAP_aspect", "mean_AbsSHAP_landcover_triclass")) |>
  cols_hide(c(mean_swe_max)) |>
  opt_vertical_padding(scale = 1)
# tab_footnote(
#   footnote = ""
# )
CONUSClasses_SnowExTitle_gt
### saving table
CONUSClasses_SnowExTitle_gt |>
  gtsave(
    "CONUSClusterSnowExSiteCounts_Title.png", expand = 5,
    path = here("Outputs", "ChartsGraphsTables", "ClusterSnowEx")
  )




##### Making Tables Showing info about SnowEx Sites and class info #####
### Alaska
AKSnowExClass_gt <- SnowExAK_sf |>
  gt() |> # use 'gt' to make an awesome table...
  gt_theme_guardian() |>
  # tab_header(
  #   title = paste(), # ...with this title
  #   subtitle = "")  |>  # and this subtitle
  ##tab_style(style = cell_fill("bisque"),
  ##        locations = cells_body()) |>  # add fill color to table
  fmt_number( # A column (numeric data)
    columns = c(ModalClass), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 0 # With four decimal places
  ) |>
  fmt_number( # A column (numeric data)
    columns = c(NumUniqueClasses), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 0 # With four decimal places
  ) |>
  # data_color( # Update cell colors, testing different color palettes
  #   columns = c(median_swe_max),
  #   fn = scales::col_numeric( # <- bc it's numeric
  #     palette = brewer.pal(11, "RdBu"), # A color scheme (gradient)
  #     domain = c(), # Column scale endpoints
  #     reverse = FALSE
  #   )
  # ) |>
  # data_color( # Update cell colors, testing different color palettes
  #   columns = c(SnowExSites), # ...for dose column
  #   # direction = "row",
  #   fn = scales::col_numeric( # <- bc it's numeric
  #     palette = brewer.pal(9, "Blues"), # A color scheme (gradient)
  #     domain = c(), # Column scale endpoints
  #     reverse = FALSE
  #   )
  # ) |>
  cols_label(site_name = "Site Name", ModalClass = "Modal SWE Class", NumUniqueClasses = "# of Unique SWE Classes", SturmClass_string = "Sturm Class") |> # Update labels
  # cols_move_to_end(columns = c("mean_AbsSHAP_DecFeb_prcpSumCDMSum", "mean_AbsSHAP_DecFeb_tmeanCDMSum", "mean_AbsSHAP_elevation", "mean_AbsSHAP_slope", "mean_AbsSHAP_aspect", "mean_AbsSHAP_landcover_triclass")) |>
  cols_hide(c(state, geometry, Class1993, Class1994, Class1995, Class1996, Class1997, Class1998, Class1999, Class2000, Class2001, Class2002, Class2003, Class2004, Class2005, Class2006, Class2007, Class2008, Class2009, Class2010, Class2011, Class2012, Class2013, Class2014, Class2015, Class2016, Class2017, Class2018, Class2019, Class2020, SturmClass)) |>
  opt_vertical_padding(scale = 2)
# tab_footnote(
#   footnote = ""
# )
AKSnowExClass_gt
### saving table
AKSnowExClass_gt |>
  gtsave(
    "AKSnowExClassInfo.png", expand = 5,
    path = here("Outputs", "ChartsGraphsTables", "ClusterSnowEx")
  )

### now an alaska one with a title
AKSnowExClassTitle_gt <- SnowExAK_sf |>
  gt() |> # use 'gt' to make an awesome table...
  gt_theme_guardian() |>
  tab_header(
    title = "SnowEx SWE and Sturm Classes", # ...with this title
    subtitle = "Alaska")  |>  # and this subtitle
  opt_align_table_header(align = "center") |>
  ##tab_style(style = cell_fill("bisque"),
  ##        locations = cells_body()) |>  # add fill color to table
  fmt_number( # A column (numeric data)
    columns = c(ModalClass), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 0 # With four decimal places
  ) |>
  fmt_number( # A column (numeric data)
    columns = c(NumUniqueClasses), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 0 # With four decimal places
  ) |>
  # data_color( # Update cell colors, testing different color palettes
  #   columns = c(median_swe_max),
  #   fn = scales::col_numeric( # <- bc it's numeric
  #     palette = brewer.pal(11, "RdBu"), # A color scheme (gradient)
  #     domain = c(), # Column scale endpoints
  #     reverse = FALSE
  #   )
  # ) |>
  # data_color( # Update cell colors, testing different color palettes
  #   columns = c(SnowExSites), # ...for dose column
  #   # direction = "row",
  #   fn = scales::col_numeric( # <- bc it's numeric
  #     palette = brewer.pal(9, "Blues"), # A color scheme (gradient)
  #     domain = c(), # Column scale endpoints
  #     reverse = FALSE
  #   )
  # ) |>
  cols_label(site_name = "Site Name", ModalClass = "Modal SWE Class", NumUniqueClasses = "# of Unique SWE Classes", SturmClass_string = "Sturm Class") |> # Update labels
  # cols_move_to_end(columns = c("mean_AbsSHAP_DecFeb_prcpSumCDMSum", "mean_AbsSHAP_DecFeb_tmeanCDMSum", "mean_AbsSHAP_elevation", "mean_AbsSHAP_slope", "mean_AbsSHAP_aspect", "mean_AbsSHAP_landcover_triclass")) |>
  cols_hide(c(state, geometry, Class1993, Class1994, Class1995, Class1996, Class1997, Class1998, Class1999, Class2000, Class2001, Class2002, Class2003, Class2004, Class2005, Class2006, Class2007, Class2008, Class2009, Class2010, Class2011, Class2012, Class2013, Class2014, Class2015, Class2016, Class2017, Class2018, Class2019, Class2020, SturmClass)) |>
  opt_vertical_padding(scale = 1)
# tab_footnote(
#   footnote = ""
# )
AKSnowExClassTitle_gt
### saving table
AKSnowExClassTitle_gt |>
  gtsave(
    "AKSnowExClassInfo_Title.png", expand = 5,
    path = here("Outputs", "ChartsGraphsTables", "ClusterSnowEx")
  )


### CONUS
CONUSSnowExClass_gt <- SnowExCONUS_sf |>
  gt() |> # use 'gt' to make an awesome table...
  gt_theme_guardian() |>
  # tab_header(
  #   title = paste(), # ...with this title
  #   subtitle = "")  |>  # and this subtitle
  ##tab_style(style = cell_fill("bisque"),
  ##        locations = cells_body()) |>  # add fill color to table
  fmt_number( # A column (numeric data)
    columns = c(ModalClass), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 0 # With four decimal places
  ) |>
  fmt_number( # A column (numeric data)
    columns = c(NumUniqueClasses), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 0 # With four decimal places
  ) |>
  # data_color( # Update cell colors, testing different color palettes
  #   columns = c(median_swe_max),
  #   fn = scales::col_numeric( # <- bc it's numeric
  #     palette = brewer.pal(11, "RdBu"), # A color scheme (gradient)
  #     domain = c(), # Column scale endpoints
  #     reverse = FALSE
  #   )
  # ) |>
  # data_color( # Update cell colors, testing different color palettes
  #   columns = c(SnowExSites), # ...for dose column
  #   # direction = "row",
  #   fn = scales::col_numeric( # <- bc it's numeric
  #     palette = brewer.pal(9, "Blues"), # A color scheme (gradient)
  #     domain = c(), # Column scale endpoints
  #     reverse = FALSE
  #   )
  # ) |>
  cols_label(site_name = "Site Name", ModalClass = "Modal SWE Class", NumUniqueClasses = "# of Unique SWE Classes", SturmClass_string = "Sturm Class") |> # Update labels
  # cols_move_to_end(columns = c("mean_AbsSHAP_DecFeb_prcpSumCDMSum", "mean_AbsSHAP_DecFeb_tmeanCDMSum", "mean_AbsSHAP_elevation", "mean_AbsSHAP_slope", "mean_AbsSHAP_aspect", "mean_AbsSHAP_landcover_triclass")) |>
  cols_hide(c(state, geometry, Class1993, Class1994, Class1995, Class1996, Class1997, Class1998, Class1999, Class2000, Class2001, Class2002, Class2003, Class2004, Class2005, Class2006, Class2007, Class2008, Class2009, Class2010, Class2011, Class2012, Class2013, Class2014, Class2015, Class2016, Class2017, Class2018, Class2019, Class2020, SturmClass)) |>
  opt_vertical_padding(scale = 1)
# tab_footnote(
#   footnote = ""
# )
CONUSSnowExClass_gt
### saving table
CONUSSnowExClass_gt |>
  gtsave(
    "CONUSSnowExClassInfo.png", expand = 5,
    path = here("Outputs", "ChartsGraphsTables", "ClusterSnowEx")
  )

### now an alaska one with a title
CONUSSnowExClassTitle_gt <- SnowExCONUS_sf|>
  gt() |> # use 'gt' to mCONUSe an awesome table...
  gt_theme_guardian() |>
  tab_header(
    title = "SnowEx SWE and Sturm Classes", # ...with this title
    subtitle = "Western US")  |>  # and this subtitle
  opt_align_table_header(align = "center") |>
  ##tab_style(style = cell_fill("bisque"),
  ##        locations = cells_body()) |>  # add fill color to table
  fmt_number( # A column (numeric data)
    columns = c(ModalClass), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 0 # With four decimal places
  ) |>
  fmt_number( # A column (numeric data)
    columns = c(NumUniqueClasses), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 0 # With four decimal places
  ) |>
  # data_color( # Update cell colors, testing different color palettes
  #   columns = c(median_swe_max),
  #   fn = scales::col_numeric( # <- bc it's numeric
  #     palette = brewer.pal(11, "RdBu"), # A color scheme (gradient)
  #     domain = c(), # Column scale endpoints
  #     reverse = FALSE
  #   )
  # ) |>
  # data_color( # Update cell colors, testing different color palettes
  #   columns = c(SnowExSites), # ...for dose column
  #   # direction = "row",
  #   fn = scales::col_numeric( # <- bc it's numeric
  #     palette = brewer.pal(9, "Blues"), # A color scheme (gradient)
  #     domain = c(), # Column scale endpoints
  #     reverse = FALSE
  #   )
  # ) |>
  cols_label(site_name = "Site Name", ModalClass = "Modal SWE Class", NumUniqueClasses = "# of Unique SWE Classes", SturmClass_string = "Sturm Class") |> # Update labels
  # cols_move_to_end(columns = c("mean_AbsSHAP_DecFeb_prcpSumCDMSum", "mean_AbsSHAP_DecFeb_tmeanCDMSum", "mean_AbsSHAP_elevation", "mean_AbsSHAP_slope", "mean_AbsSHAP_aspect", "mean_AbsSHAP_landcover_triclass")) |>
  cols_hide(c(state, geometry, Class1993, Class1994, Class1995, Class1996, Class1997, Class1998, Class1999, Class2000, Class2001, Class2002, Class2003, Class2004, Class2005, Class2006, Class2007, Class2008, Class2009, Class2010, Class2011, Class2012, Class2013, Class2014, Class2015, Class2016, Class2017, Class2018, Class2019, Class2020, SturmClass)) |>
  opt_vertical_padding(scale = 1)
# tab_footnote(
#   footnote = ""
# )
CONUSSnowExClassTitle_gt
### saving table
CONUSSnowExClassTitle_gt |>
  gtsave(
    "CONUSSnowExClassInfo_Title.png", expand = 5,
    path = here("Outputs", "ChartsGraphsTables", "ClusterSnowEx")
  )


