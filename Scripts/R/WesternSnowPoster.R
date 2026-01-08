##### Getting Stuff Ready for Western Snow Conference Poster #####

##### loading packages, reading in Data #####
library(pacman)
p_load(tidyverse, here, sf, terra, randomForest, caret, treeshap, parallel, mcprogress, tidyterra, gt, gtExtras, webshot2)

### CONUS AOI
CONUS_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "CONUS", "CONUS_AOI.gpkg"))
CONUS_States <- read_sf(here("Data", "L3_Ecoregions_USB", "CONUS", "CONUSEcoregionsNoStates.gpkg"))
### Alaska AOIs
AK_Ecoregions <- read_sf(here("Data", "L3_Ecoregions_USB", "Alaska", "AK_StateBoundary.gpkg"))
AK_State <- read_sf(here("Data", "L3_Ecoregions_USB", "Alaska", "AK_StateBoundary.gpkg"))
### SNOTEL Annual values for model fitting and evaluation
Snotel_CONUS <- read_sf(here("Data", "SNOTEL", "CONUS", "GIS", "GPKG", "Annual", "SnotelData_CONUS_PeakSWE.gpkg"))
Snotel_AK <- read_sf(here("Data", "SNOTEL", "Alaska", "GIS", "GPKG", "Annual", "SnotelData_AK_PeakSWE.gpkg"))
Snotel_sf <- rbind(Snotel_CONUS, st_transform(Snotel_AK, crs(Snotel_CONUS)))
UniqueCONUSSites <- Snotel_sf |>
  filter(site_id != -999) |>
  filter(state != "AK") |>
  group_by(site_id) |>
  summarise(mean_swe = mean(peak_swe)) |>
  mutate(Snotel_field = "SNOTEL Station")

UniqueAKSites <- Snotel_sf |>
  filter(site_id != -999) |>
  filter(state == "AK") |>
  group_by(site_id) |>
  summarise(mean_swe = mean(peak_swe)) |>
  mutate(Snotel_field = "SNOTEL Station")

### reading in CONUS DEM to make plot
CONUS_DEM <- rast(here("Data", "DEM", "CONUS", "CONUSDEMMosaic.tif"))
### Alaska
AK_DEM <- rast(here("Data", "DEM", "Alaska", "AKDEMMosaic.tif"))

##### Making Study Area Plots #####
# CONUSSite_plot <- ggplot() +
#   theme_bw() +
#   geom_spatraster(data = CONUS_DEM, show.legend = FALSE) +
#   scale_fill_whitebox_c(palette = "high_relief", direction = 1) +
#   geom_sf(data = CONUS_States, fill = NA, color = "black") +
#   geom_sf(data = UniqueCONUSSites, fill = UniqueCONUSSites$Snotel_field) +
#   ggtitle(label = "SNOTEL Sites") +
#   ggspatial::annotation_scale(location = "br") +
#   ggspatial::annotation_north_arrow(location = "bl", height = unit(1.25, "cm"), width = unit(1.25, "cm")) +
#   theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
# CONUSSite_plot
### saving ggplot
# ggsave(here("Outputs", "Plots", "WesternSnow", "CONUS", "1Km", paste0("CONUSmSWE", year, ".png")))



##### Extracting Snow Cover Classes to SnowEx sites #####
### reading in CSV of SnowEx site locations
SnowEx_df <- read_csv(here("Data", "SnowEx", "SnowExSiteLocations.csv"))
SnowExCONUS_df <- SnowEx_df |>
  filter(State != "AK") |>
  ### filtering out this site because I don't have good coordinates for it
  filter(SiteName != "Fraser Experimental Forest, CO" & SiteName != "Snodgrass, CO")
SnowExAK_df <- SnowEx_df |>
  filter(State == "AK")
### converting df to sf object
SnowExAK_sf <- st_transform(st_as_sf(SnowExAK_df, coords = c("Lon", "Lat"), crs = 4326), crs(AK_Ecoregions)) |>
  mutate(Class1993 = "dummystring",
         Class1994 = "dummystring",
         Class1995 = "dummystring",
         Class1996 = "dummystring",
         Class1997 = "dummystring",
         Class1998 = "dummystring",
         Class1999 = "dummystring",
         Class2000 = "dummystring",
         Class2001 = "dummystring",
         Class2002 = "dummystring",
         Class2003 = "dummystring",
         Class2004 = "dummystring",
         Class2005 = "dummystring",
         Class2006 = "dummystring",
         Class2007 = "dummystring",
         Class2008 = "dummystring",
         Class2009 = "dummystring",
         Class2010 = "dummystring",
         Class2011 = "dummystring",
         Class2012 = "dummystring",
         Class2013 = "dummystring",
         Class2014 = "dummystring",
         Class2015 = "dummystring",
         Class2016 = "dummystring",
         Class2017 = "dummystring",
         Class2018 = "dummystring",
         Class2019 = "dummystring",
         Class2020 = "dummystring")
SnowExCONUS_sf <- st_transform(st_as_sf(SnowExCONUS_df, coords = c("Lon", "Lat"), crs = 4326), crs(CONUS_AOI)) |>
  mutate(Class1993 = "dummystring",
         Class1994 = "dummystring",
         Class1995 = "dummystring",
         Class1996 = "dummystring",
         Class1997 = "dummystring",
         Class1998 = "dummystring",
         Class1999 = "dummystring",
         Class2000 = "dummystring",
         Class2001 = "dummystring",
         Class2002 = "dummystring",
         Class2003 = "dummystring",
         Class2004 = "dummystring",
         Class2005 = "dummystring",
         Class2006 = "dummystring",
         Class2007 = "dummystring",
         Class2008 = "dummystring",
         Class2009 = "dummystring",
         Class2010 = "dummystring",
         Class2011 = "dummystring",
         Class2012 = "dummystring",
         Class2013 = "dummystring",
         Class2014 = "dummystring",
         Class2015 = "dummystring",
         Class2016 = "dummystring",
         Class2017 = "dummystring",
         Class2018 = "dummystring",
         Class2019 = "dummystring",
         Class2020 = "dummystring")

##### Reading in Class Rasters and Extracting classes to SnowEx objects #####
years <- 1993:2020
for (i in 1:length(years)){
  CONUSClasses_rast <- rast(here("WesternSnow", "Classes", "CONUS", paste0("CONUS", years[i], "SHAPClasses.tif")))
  AKClasses_rast <- rast(here("WesternSnow", "Classes", "Alaska", paste0("Alaska", years[i], "SHAPClasses.tif")))
  SnowExCONUS_sf[,i + 3] = terra::extract(CONUSClasses_rast, SnowExCONUS_sf)[,2]
  SnowExAK_sf[,i + 3] = terra::extract(AKClasses_rast, SnowExAK_sf)[,2]
}

### Calculating most common class and number of unique classes
SnowExCONUS_sf <- SnowExCONUS_sf %>%
  mutate(
    most_common_class = apply(.[, 4:31], 1, function(row) {
      # Flatten the row into a single vector of strings
      all_strings <- unlist(row)
      
      # Calculate the frequency of each string
      string_counts <- table(all_strings)
      
      # Find the string with the maximum frequency
      if (length(string_counts) > 0) {
        most_common <- names(string_counts)[which.max(string_counts)]
        return(most_common)
      } else {
        return(NA_character_) # Or some other default if the row is all NA/empty
      }
    }),
    unique_classes = apply(.[, 4:31], 1, function(row) {
      all_strings <- unlist(row)
      unique_count <- length(unique(all_strings))
      return(unique_count)
    })
  )

SnowExAK_sf <- SnowExAK_sf %>%
  mutate(
    most_common_class = apply(.[, 4:31], 1, function(row) {
      # Flatten the row into a single vector of strings
      all_strings <- unlist(row)
      
      # Calculate the frequency of each string
      string_counts <- table(all_strings)
      
      # Find the string with the maximum frequency
      if (length(string_counts) > 0) {
        most_common <- names(string_counts)[which.max(string_counts)]
        return(most_common)
      } else {
        return(NA_character_) # Or some other default if the row is all NA/empty
      }
    }),
    unique_classes = apply(.[, 4:31], 1, function(row) {
      all_strings <- unlist(row)
      unique_count <- length(unique(all_strings))
      return(unique_count)
    })
  )

### saving sf objects
write_sf(SnowExAK_sf, here("Data", "SnowEx", "SnowExAK_TempClasses.gpkg"), append = FALSE)
write_sf(SnowExCONUS_sf, here("Data", "SnowEx", "SnowExCONUS_TempClasses.gpkg"), append = FALSE)
### combining sf objects to make table of most common class
SnowExCombined_df <- rbind(st_drop_geometry(SnowExCONUS_sf), st_drop_geometry(SnowExAK_sf)) |>
  filter(SiteName != "Southern Sierra CZO")

SnowExCommonClass_gt <- SnowExCombined_df |>
  arrange(State) |>
  gt() |> # use 'gt' to make an awesome table...
  gt_theme_538() |>
  tab_header(
    title = "Most Common SWE-based Classes", # ...with this title
    subtitle = "Class is defined as variables listed in order of absolute value of SHAP values")  |>  # and this subtitle
  tab_style(style = cell_fill("bisque"),
            locations = cells_body()) |>  # add fill color to table
  cols_label(SiteName = "Site Name", most_common_class = "Most Common Class/ Variable Importance Order", unique_classes = "Number of Unique Classes") |> # Update labels
  # cols_move_to_end(columns = "VoA_Rating") |>
  cols_hide(c(State, Class1993, Class1994, Class1995, Class1996, Class1997, Class1998, Class1999, Class2000, Class2001, Class2002, Class2003, Class2004, Class2005,Class2006, Class2007, Class2008, Class2009, Class2010,Class2011, Class2012, Class2013, Class2014, Class2015,Class2016, Class2017, Class2018, Class2019, Class2020)) |>
  tab_footnote(
    footnote = "CDM = Cooling Degree Month, where prcpSumCDM & tmeanCDM = 0 if mean temp > 10 C, otherwise prcpSum value is kept and tmeanCDM = 10 - tmean. prcpSum and tmean CDM values are then summed from October-April of a given water year. landcover_triclass = either 0, 1, or 2 corresponding with no vegetation, short vegetation, or tall vegetation"
  ) |>
  tab_options(footnotes.font.size = 12, table.font.size = 20)
SnowExCommonClass_gt
### saving table
SnowExCommonClass_gt |>
  gtsave(
    "SnowExMostCommonClassTable.png", expand = 5,
    path = here("WesternSnow", "PlotsMaps")
  )


##### Making Table of Final Model Metrics #####
FinalModelMetrics_df <- data.frame(FinalR2 = 0.8522428,
                                   TestRMSE = 132.4932)

FinalModelMetrics_gt <- FinalModelMetrics_df |>
  gt() |>
  gt_theme_538() |>
  tab_header(
    title = "Final Model Evaluation")  |>
  tab_style(style = cell_fill("bisque"),
            locations = cells_body()) |>
  cols_label(FinalR2 = "Final R^2", TestRMSE = "Test RMSE")
FinalModelMetrics_gt
### saving table
FinalModelMetrics_gt |>
  gtsave(
    "FinalModelMetricsTable.png", expand = 5,
    path = here("WesternSnow", "PlotsMaps")
  )

##### Making Maps of SWE classes for 2020 #####
CONUSClasses2020_rast <- rast(here("WesternSnow", "Classes", "CONUS", "CONUS2020SHAPClasses.tif"))
AKClasses2020_rast <- rast(here("WesternSnow", "Classes", "Alaska", "Alaska2020SHAPClasses.tif"))


# CONUSSWEClass_plot <- ggplot() +
#   theme_bw() +
#   geom_spatraster(data = CONUSClasses2020_rast, show.legend = FALSE, na.rm = TRUE) +
#   # scale_fill_wiki_d(direction = 1) +
#   geom_sf(data = CONUS_States, fill = NA, color = "black") +
#   ggtitle(label = "SWE-based Classes for Water Year 2020") +
#   ggspatial::annotation_scale(location = "br") +
#   ggspatial::annotation_north_arrow(location = "bl", height = unit(1.25, "cm"), width = unit(1.25, "cm")) +
#   theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
# CONUSSWEClass_plot
# ggsave(here("WesternSnow", "PlotsMaps", "CONUS2020ClassMap.png"))
