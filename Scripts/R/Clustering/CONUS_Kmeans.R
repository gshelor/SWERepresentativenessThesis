##### K means classification based on fastshap #####
totalstarttime <- Sys.time()
##### Loading packages #####
library(pacman)
p_load(tidyverse, here, sf, terra, sfext, parallel, data.table, mcprogress, dtplyr, gt, gtExtras, RColorBrewer)
# options(mc.cores = parallel::detectCores())
# setwd("/media/Research/Morafkai/GriffinS")

### reading in csvs
CONUSSHAP_csvs <- list.files(path = here("Data", "SHAP", "Fastshap", "SHAPs", "CONUS"), pattern = ".csv", full.names = TRUE)

years <- 1993:2020

ReadAllYears_func <- function(year){
  temp_shaps <- fread(file = paste0(here("Data", "SHAP", "Fastshap", "SHAPs", "CONUS"), paste0("CONUS", year, "SHAPs.csv")))
}
CONUSSHAPs <- pmclapply(X = years, FUN = ReadAllYears_func, mc.cores = length(years))

ReadAllYearsPreds_func <- function(year){
  temp_preds <- fread(file = here("Data", "SHAP", "PredictorCSVs", "CONUS", paste0("CONUS_SWEmaxPredictors", year, ".csv")))
}
CONUSPredictors <- pmclapply(X = years, FUN = ReadAllYearsPreds_func, mc.cores = length(years))

for (x in 1:length(CONUSSHAPs)){
  temp_dt <- CONUSSHAPs[[x]]
  year <- years[x]
  temp_dt <- temp_dt |>
    mutate(WaterYear = year,
           pastexyWY = paste0(x, y, WaterYear))
  if (x == 1){
    CONUS_SHAPs_All <- temp_dt
  } else{
    CONUS_SHAPs_All <- rbind(CONUS_SHAPs_All, temp_dt)
  }
}

for (x in 1:length(CONUSPredictors)){
  temp_dt <- CONUSPredictors[[x]]
  year <- years[x]
  temp_dt <- temp_dt |>
    mutate(WaterYear = year,
           pastexyWY = paste0(x, y, WaterYear))
  if (x == 1){
    CONUS_Predictors_All <- temp_dt
  } else{
    CONUS_Predictors_All <- rbind(CONUS_Predictors_All, temp_dt)
  }
}


CONUS_SHAPs_All <- CONUS_SHAPs_All |>
  mutate(AbsSHAP_OctApr_prcpSumCDMSum = abs(OctApr_prcpSumCDMSum),
         AbsSHAP_OctApr_tmeanCDMSum = abs(OctApr_tmeanCDMSum),
         AbsSHAP_elevation = abs(elevation),
         AbsSHAP_slope = abs(slope),
         AbsSHAP_aspect = abs(aspect),
         AbsSHAP_landcover_triclass = abs(landcover_triclass))

# fwrite(CONUS_SHAPs_All, file = "SHAPAnalysis/Data/SHAP/Fastshap/SHAPs/CONUS/AllYearsCombined/CONUSAllSHAPs.csv")
# CONUS_SHAPs_All <- fread("SHAPAnalysis/Data/SHAP/Fastshap/SHAPs/CONUS/AllYearsCombined/CONUSAllSHAPs.csv")


##### Kmeans Clustering CONUS SHAP absolute values #####
### initializing df to store clustering error metrics
KmeansErrors_df <- data.frame(k = 1:25,
                              sse = numeric(25))

Kmeans_func <- function(i){
  set.seed(802)
  kmeans_error <- kmeans(CONUS_SHAPs_All[, .(peak_swe,AbsSHAP_OctApr_prcpSumCDMSum, AbsSHAP_OctApr_tmeanCDMSum, AbsSHAP_elevation, AbsSHAP_slope, AbsSHAP_aspect, AbsSHAP_landcover_triclass)], centers = i)$tot.withinss
  return(kmeans_error)
}

Kmeans_errors <- pmclapply(X = KmeansErrors_df$k, FUN = Kmeans_func, mc.cores = 6, mc.silent = TRUE, mc.set.seed = FALSE)

KmeansErrors_df$sse <- unlist(Kmeans_errors)

KmeansErrors_df$pct_error_change <- numeric(25)
for (i in 1:25){
  if (i == 1){
    KmeansErrors_df$pct_error_change[i] = 1
  } else{
    KmeansErrors_df$pct_error_change[i] = (KmeansErrors_df$sse[i-1] - KmeansErrors_df$sse[i]) / KmeansErrors_df$sse[i - 1]
  }
}
fwrite(KmeansErrors_df, here("SHAPAnalysis", "Data", "SHAP", "Classes", "CONUS", "KMeans", "CSVs", "Summary", "KmeansErrors.csv"))
KmeansError_plot <- ggplot(data = KmeansErrors_df, mapping = aes(x = k, y = sse)) +
  theme_bw() +
  geom_point() +
  geom_line() +
  # ggtitle(label = "Total Sum of Squared Error for Different Numbers of K-Means Clusters", subtitle = "WUS Study Area") +
  ggtitle(label = "WUS Study Area") +
  xlab("K Number of Clusters") +
  ylab("Total SSE") +
  theme(plot.title = element_text(hjust = 0.5, size = 12), plot.subtitle = element_text(hjust = 0.5))
KmeansError_plot
ggsave(plot = KmeansError_plot, here("Outputs", "Plots", "KmeansError", "WUSKmeansSSEPlot.png"))

# break
##### Final KMeans with appropriate number of clusters #####
set.seed(802)
FinalKmeans <- kmeans(CONUS_SHAPs_All[, .(peak_swe, AbsSHAP_OctApr_prcpSumCDMSum, AbsSHAP_OctApr_tmeanCDMSum, AbsSHAP_elevation, AbsSHAP_slope, AbsSHAP_aspect, AbsSHAP_landcover_triclass)], centers = 9)

### assigning initial class number
CONUS_SHAPs_All <- CONUS_SHAPs_All |>
  mutate(class_cluster = FinalKmeans$cluster)

### calculating averages of each cluster so I can reset class numbers to be based on max SWE
CONUSCluster_Avg <- CONUS_SHAPs_All |>
  group_by(class_cluster) |>
  summarise(mean_swe_max = mean(peak_swe),
            mean_AbsSHAP_OctApr_prcpSumCDMSum = mean(AbsSHAP_OctApr_prcpSumCDMSum), 
            mean_AbsSHAP_OctApr_tmeanCDMSum = mean(AbsSHAP_OctApr_tmeanCDMSum),
            mean_AbsSHAP_elevation = mean(AbsSHAP_elevation), 
            mean_AbsSHAP_slope = mean(AbsSHAP_slope),
            mean_AbsSHAP_aspect = mean(AbsSHAP_aspect), 
            mean_AbsSHAP_landcover_triclass = mean(AbsSHAP_landcover_triclass)) |>
  mutate(swe_max_rank = dense_rank(desc(mean_swe_max)))

### reassigning class number based on descending order of average max swe for the classes
CONUS_SHAPs_All <- CONUS_SHAPs_All |>
  mutate(class_cluster2 = case_when(class_cluster == 1 ~ CONUSCluster_Avg$swe_max_rank[CONUSCluster_Avg$class_cluster == 1],
                                    class_cluster == 2 ~ CONUSCluster_Avg$swe_max_rank[CONUSCluster_Avg$class_cluster == 2],
                                    class_cluster == 3 ~ CONUSCluster_Avg$swe_max_rank[CONUSCluster_Avg$class_cluster == 3],
                                    class_cluster == 4 ~ CONUSCluster_Avg$swe_max_rank[CONUSCluster_Avg$class_cluster == 4],
                                    class_cluster == 5 ~ CONUSCluster_Avg$swe_max_rank[CONUSCluster_Avg$class_cluster == 5],
                                    class_cluster == 6 ~ CONUSCluster_Avg$swe_max_rank[CONUSCluster_Avg$class_cluster == 6],
                                    class_cluster == 7 ~ CONUSCluster_Avg$swe_max_rank[CONUSCluster_Avg$class_cluster == 7],
                                    class_cluster == 8 ~ CONUSCluster_Avg$swe_max_rank[CONUSCluster_Avg$class_cluster == 8],
                                    class_cluster == 9 ~ CONUSCluster_Avg$swe_max_rank[CONUSCluster_Avg$class_cluster == 9]),
         class_cluster = class_cluster2)

### dropping class_cluster2 since it was only supposed to be temporary
CONUS_SHAPs_All <- CONUS_SHAPs_All |>
  select(-one_of("class_cluster2")) |>
  mutate(pastexyWY = paste0(x, y, WaterYear))

### dropping xy from conus predictors since it's already in the SHAPs df
CONUS_Predictors_All <- CONUS_Predictors_All |>
  select(OctApr_prcpSumCDMSum, OctApr_tmeanCDMSum, elevation, slope, aspect, landcover_triclass, WaterYear, pastexyWY)

### binding SHAPs to predictors so I can make a table later
## changing colnames first since fastshap calls the columns containing the SHAP values for each variable by just the name of the variable, which is kind of inconvenient
colnames(CONUS_SHAPs_All) <- c("x", "y", "peak_swe", "SHAPOctApr_prcpSumCDMSum", "SHAPOctApr_tmeanCDMSum", "SHAPelevation", "SHAPslope", "SHAPaspect", "SHAPlandcover_triclass", "Rank_prcpSumCDM_varimp", "Rank_tmeanCDM_varimp", "Rank_elevation_varimp", "Rank_slope_varimp", "Rank_aspect_varimp", "Rank_landcover_varimp", "class", "class_top1", "class_top2", "class_top3", "class_top4", "class_top5", "WaterYear", "AbsSHAP_OctApr_prcpSumCDMSum", "AbsSHAP_OctApr_tmeanCDMSum", "AbsSHAP_elevation", "AbsSHAP_slope", "AbsSHAP_aspect", "AbsSHAP_landcover_triclass", "class_cluster", "pastexyWY")


# CONUS_SHAPs_All <- full_join(CONUS_SHAPs_All, CONUS_Predictors_All, by = "pastexyWY")
CONUS_Predictors_All <- full_join(CONUS_SHAPs_All[, .(x, y, peak_swe, class_cluster, pastexyWY)], CONUS_Predictors_All, by = "pastexyWY") |>
  select(-one_of("pastexyWY"))

CONUS_SHAPs_All <- CONUS_SHAPs_All |>
  select(x, y, peak_swe, SHAPOctApr_prcpSumCDMSum, SHAPOctApr_tmeanCDMSum, SHAPelevation, SHAPslope, SHAPaspect, SHAPlandcover_triclass, Rank_prcpSumCDM_varimp, Rank_tmeanCDM_varimp, Rank_elevation_varimp, Rank_slope_varimp, Rank_aspect_varimp, Rank_landcover_varimp, WaterYear, AbsSHAP_OctApr_prcpSumCDMSum, AbsSHAP_OctApr_tmeanCDMSum, AbsSHAP_elevation, AbsSHAP_slope, AbsSHAP_aspect, AbsSHAP_landcover_triclass, class_cluster)

### writing out csvs so I never have to do this again even though they're huge files that take up a shitton of space
fwrite(CONUS_SHAPs_All, file = here("Data", "SHAP", "Fastshap", "SHAPs", "CONUS", "AllYearsCombined", "CONUSAllSHAPs_KClasses.csv"))
fwrite(CONUS_Predictors_All, file = here("Data", "SHAP", "Fastshap", "SHAPs", "CONUS", "AllYearsCombined", "CONUSAllPredictors_KClasses.csv"))
# rm()
CONUS_SHAPs_All <- fread(here("Data", "SHAP", "Fastshap", "SHAPs", "CONUS", "AllYearsCombined", "CONUSAllSHAPs_KClasses.csv"))
##### Rasterizing K-means-based classes #####
CONUS_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "CONUS", "CONUSEcoregionsNoStates.gpkg"))
years <- 1993:2020
SaveClasses_func <- function(year){
  PredictionRast <- rast(here("Outputs", "MapTIFFs", "CONUS", paste0("CONUS_SWEmax", year, ".tif")))
  temp_classes <- CONUS_SHAPs_All |>
    filter(WaterYear == year)
  fwrite(temp_classes, file = here("Data", "SHAP", "Classes", "CONUS", "KMeans", "CSVs", paste0("CONUS", year, "_KClasses.csv")))
  ### converting year's classes to raster
  temp_classes_rast <- terra::rasterize(vect(st_as_sf(temp_classes, coords = c("x", "y"), crs = crs(CONUS_AOI))), y = PredictionRast, field = "class_cluster")
  writeRaster(temp_classes_rast, here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", paste0("CONUS", year, "_KClasses.tif")), overwrite = TRUE)
  # return(temp_classes_rast)
}

### rasterizing classes for each year
SWEClassesByYear <- pmclapply(X = years, FUN = SaveClasses_func, mc.cores = 12, mc.silent = TRUE, mc.set.seed = FALSE)

### reading in class rasters to calculate most common class and number of unique classes
CONUS1993Classes_rast <- rast(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "CONUS1993_KClasses.tif"))
CONUS1994Classes_rast <- rast(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "CONUS1994_KClasses.tif"))
CONUS1995Classes_rast <- rast(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "CONUS1995_KClasses.tif"))
# plot(CONUS1994Classes_rast)
CONUS1996Classes_rast <- rast(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "CONUS1996_KClasses.tif"))
# plot(CONUS1994Classes_rast)
CONUS1997Classes_rast <- rast(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "CONUS1997_KClasses.tif"))
# plot(CONUS1994Classes_rast)
CONUS1998Classes_rast <- rast(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "CONUS1998_KClasses.tif"))
# plot(CONUS1994Classes_rast)
CONUS1999Classes_rast <- rast(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "CONUS1999_KClasses.tif"))
# plot(CONUS1994Classes_rast)
CONUS2000Classes_rast <- rast(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "CONUS2000_KClasses.tif"))
# plot(CONUS1994Classes_rast)
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

CONUSAllYearsClasses_rast <- c(CONUS1993Classes_rast, CONUS1994Classes_rast, CONUS1995Classes_rast, CONUS1996Classes_rast, CONUS1997Classes_rast, CONUS1998Classes_rast, CONUS1999Classes_rast, CONUS2000Classes_rast, CONUS2001Classes_rast, CONUS2002Classes_rast, CONUS2003Classes_rast, CONUS2004Classes_rast, CONUS2005Classes_rast, CONUS2006Classes_rast, CONUS2007Classes_rast, CONUS2008Classes_rast, CONUS2009Classes_rast, CONUS2010Classes_rast, CONUS2011Classes_rast, CONUS2012Classes_rast, CONUS2013Classes_rast, CONUS2014Classes_rast, CONUS2015Classes_rast, CONUS2016Classes_rast, CONUS2017Classes_rast, CONUS2018Classes_rast, CONUS2019Classes_rast, CONUS2020Classes_rast)

CONUSModalClass_rast <- as.int(app(CONUSAllYearsClasses_rast, fun = "modal"))
CONUSNumUniqueClass_rast <- app(CONUSAllYearsClasses_rast, fun = function(x){length(unique(x))})
CONUSNumUniqueClass_rast <- crop(CONUSNumUniqueClass_rast, CONUS_AOI, mask = TRUE)
# CONUSUniqueClasses_rast <- app(CONUSAllYearsClasses_rast, fun = function(x){as.factor(unique(x))})
# CONUSUniqueClasses_rast <- crop(CONUSUniqueClasses_rast, CONUS_AOI, mask = TRUE)
plot(as.factor(CONUSModalClass_rast), main = "most common class", reverse = TRUE)
plot(st_geometry(CONUS_AOI), add = TRUE)
plot(CONUSNumUniqueClass_rast, main = "number of unique classes")
plot(st_geometry(CONUS_AOI), add = TRUE)

### saving most common class raster and number of unique classes raster
writeRaster(CONUSModalClass_rast, here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "Summary", "CONUSModalClass.tif"), overwrite = TRUE)
writeRaster(CONUSNumUniqueClass_rast, here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "Summary", "CONUSNumUniqueClass.tif"), overwrite = TRUE)

##### Calculating Averages by Cluster #####
CONUSCluster_Avg <- CONUS_SHAPs_All |>
  group_by(class_cluster) |>
  summarise(mean_swe_max = mean(peak_swe),
            median_swe_max = median(peak_swe),
            mean_SHAPOctApr_prcpSumCDMSum = mean(SHAPOctApr_prcpSumCDMSum), 
            mean_SHAPOctApr_tmeanCDMSum = mean(SHAPOctApr_tmeanCDMSum),
            mean_SHAPelevation = mean(SHAPelevation), 
            mean_SHAPslope = mean(SHAPslope),
            mean_SHAPaspect = mean(SHAPaspect), 
            mean_SHAPlandcover_triclass = mean(SHAPlandcover_triclass),
            mean_AbsSHAP_OctApr_prcpSumCDMSum = mean(AbsSHAP_OctApr_prcpSumCDMSum), 
            mean_AbsSHAP_OctApr_tmeanCDMSum = mean(AbsSHAP_OctApr_tmeanCDMSum),
            mean_AbsSHAP_elevation = mean(AbsSHAP_elevation), 
            mean_AbsSHAP_slope = mean(AbsSHAP_slope),
            mean_AbsSHAP_aspect = mean(AbsSHAP_aspect), 
            mean_AbsSHAP_landcover_triclass = mean(AbsSHAP_landcover_triclass),
            median_SHAPOctApr_prcpSumCDMSum = median(SHAPOctApr_prcpSumCDMSum), 
            median_SHAPOctApr_tmeanCDMSum = median(SHAPOctApr_tmeanCDMSum),
            median_SHAPelevation = median(SHAPelevation), 
            median_SHAPslope = median(SHAPslope),
            median_SHAPaspect = median(SHAPaspect), 
            median_SHAPlandcover_triclass = median(SHAPlandcover_triclass),
            median_AbsSHAP_OctApr_prcpSumCDMSum = median(AbsSHAP_OctApr_prcpSumCDMSum), 
            median_AbsSHAP_OctApr_tmeanCDMSum = median(AbsSHAP_OctApr_tmeanCDMSum),
            median_AbsSHAP_elevation = median(AbsSHAP_elevation), 
            median_AbsSHAP_slope = median(AbsSHAP_slope),
            median_AbsSHAP_aspect = median(AbsSHAP_aspect), 
            median_AbsSHAP_landcover_triclass = median(AbsSHAP_landcover_triclass))
### writing group averages as csv
fwrite(CONUSCluster_Avg, here("Data", "SHAP", "Classes", "CONUS", "KMeans", "CSVs", "Summary", "CONUSClassSHAPAvgs.csv"))



### making table of cluster averages
CONUSCluster_Avg_gt <- CONUSCluster_Avg |>
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
    columns = c(mean_AbsSHAP_OctApr_prcpSumCDMSum), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # A column (numeric data)
    columns = c(mean_AbsSHAP_OctApr_tmeanCDMSum), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # Another column (also numeric data)
    columns = c(mean_AbsSHAP_elevation), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
  ) |>
  fmt_number( # Another column (also numeric data)
    columns = c(mean_AbsSHAP_slope), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
  ) |>
  fmt_number( # Another column (also numeric data)
    columns = c(mean_AbsSHAP_aspect), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
  ) |>
  fmt_number( # Another column (also numeric data)
    columns = c(mean_AbsSHAP_landcover_triclass), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
  ) |>
  data_color( # Update cell colors, testing different color palettes
    columns = c(median_swe_max), # ...for dose column
    fn = scales::col_numeric( # <- bc it's numeric
      palette = brewer.pal(11, "RdBu"), # A color scheme (gradient)
      domain = c(), # Column scale endpoints
      reverse = FALSE
    )
  ) |>
  data_color( # Update cell colors, testing different color palettes
    columns = c(mean_AbsSHAP_OctApr_prcpSumCDMSum, mean_AbsSHAP_OctApr_tmeanCDMSum, mean_AbsSHAP_elevation, mean_AbsSHAP_slope, mean_AbsSHAP_aspect, mean_AbsSHAP_landcover_triclass), # ...for dose column
    direction = "row",
    fn = scales::col_numeric( # <- bc it's numeric
      palette = brewer.pal(11, "PRGn"), # A color scheme (gradient)
      domain = c(), # Column scale endpoints
      reverse = FALSE
    )
  ) |>
  # data_color( # Update cell colors, testing different color palettes
  #   columns = c(mean_AbsSHAP_OctApr_tmeanCDMSum), # ...for dose column
  #   fn = scales::col_numeric( # <- bc it's numeric
  #     palette = brewer.pal(11, "PRGn"), # A color scheme (gradient)
  #     domain = c(), # Column scale endpoints
  #     reverse = FALSE
  #   )
  # ) |>
  # data_color( # Update cell colors, testing different color palettes
  #   columns = c(mean_AbsSHAP_elevation), # ...for dose column
  #   fn = scales::col_numeric( # <- bc it's numeric
  #     palette = brewer.pal(11, "PRGn"), # A color scheme (gradient)
  #     domain = c(), # Column scale endpoints
  #     reverse = FALSE
  #   )
  # ) |>
  # data_color( # Update cell colors, testing different color palettes
  #   columns = c(mean_AbsSHAP_slope), # ...for dose column
  #   fn = scales::col_numeric( # <- bc it's numeric
  #     palette = brewer.pal(11, "PRGn"), # A color scheme (gradient)
  #     domain = c(), # Column scale endpoints
  #     reverse = FALSE
  #   )
  # ) |>
  # data_color( # Update cell colors, testing different color palettes
  #   columns = c(mean_AbsSHAP_aspect), # ...for dose column
  #   fn = scales::col_numeric( # <- bc it's numeric
  #     palette = brewer.pal(11, "PRGn"), # A color scheme (gradient)
  #     domain = c(), # Column scale endpoints
  #     reverse = FALSE
  #   )
  # ) |>
  # data_color( # Update cell colors, testing different color palettes
  #   columns = c(mean_AbsSHAP_landcover_triclass), # ...for dose column
  #   fn = scales::col_numeric( # <- bc it's numeric
  #     palette = brewer.pal(11, "PRGn"), # A color scheme (gradient)
  #     domain = c(), # Column scale endpoints
  #     reverse = FALSE
  #   )
  # ) |>
  cols_label(class_cluster = "Class", median_swe_max = "Median Max SWE (mm)", mean_AbsSHAP_OctApr_prcpSumCDMSum = "Oct-Apr prcp CDM Sum", mean_AbsSHAP_OctApr_tmeanCDMSum = "Oct-Apr tmean CDM Sum", mean_AbsSHAP_elevation = "Elevation", mean_AbsSHAP_slope = "Slope", mean_AbsSHAP_aspect = "Aspect", mean_AbsSHAP_landcover_triclass = "Landcover") |> # Update labels
  cols_move_to_end(columns = c("mean_AbsSHAP_elevation", "mean_AbsSHAP_slope", "mean_AbsSHAP_aspect", "mean_AbsSHAP_landcover_triclass")) |>
  cols_hide(c(mean_SHAPOctApr_prcpSumCDMSum, mean_SHAPOctApr_tmeanCDMSum, mean_SHAPelevation, mean_SHAPslope, mean_SHAPaspect, mean_SHAPlandcover_triclass, mean_swe_max, median_SHAPOctApr_prcpSumCDMSum, median_SHAPOctApr_tmeanCDMSum, median_SHAPelevation, median_SHAPslope, median_SHAPaspect, median_SHAPlandcover_triclass, median_AbsSHAP_OctApr_prcpSumCDMSum, median_AbsSHAP_OctApr_tmeanCDMSum, median_AbsSHAP_elevation, median_AbsSHAP_slope, median_AbsSHAP_aspect, median_AbsSHAP_landcover_triclass)) |>
  opt_vertical_padding(scale = 2)
# tab_footnote(
#   footnote = ""
# )
CONUSCluster_Avg_gt
### saving table
CONUSCluster_Avg_gt |>
  gtsave(
    "CONUSClusterAvg.png", expand = 5,
    path = here("Outputs", "ChartsGraphsTables", "Kmeans")
  )

### making table of cluster averages
CONUSCluster_AvgTitle_gt <- CONUSCluster_Avg |>
  gt() |> # use 'gt' to make an awesome table...
  gt_theme_guardian() |>
  tab_header(
    title = "SWE Class Absolute SHAP Value Averages", # ...with this title
    subtitle = "Western US")  |>  # and this subtitle
  ##tab_style(style = cell_fill("bisque"),
  ##        locations = cells_body()) |>  # add fill color to table
  opt_align_table_header(align = "center") |>
  fmt_number( # A column (numeric data)
    columns = c(median_swe_max), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # A column (numeric data)
    columns = c(mean_AbsSHAP_OctApr_prcpSumCDMSum), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # A column (numeric data)
    columns = c(mean_AbsSHAP_OctApr_tmeanCDMSum), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # Another column (also numeric data)
    columns = c(mean_AbsSHAP_elevation), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
  ) |>
  fmt_number( # Another column (also numeric data)
    columns = c(mean_AbsSHAP_slope), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
  ) |>
  fmt_number( # Another column (also numeric data)
    columns = c(mean_AbsSHAP_aspect), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
  ) |>
  fmt_number( # Another column (also numeric data)
    columns = c(mean_AbsSHAP_landcover_triclass), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
  ) |>
  data_color( # Update cell colors, testing different color palettes
    columns = c(median_swe_max), # ...for dose column
    fn = scales::col_numeric( # <- bc it's numeric
      palette = brewer.pal(11, "RdBu"), # A color scheme (gradient)
      domain = c(), # Column scale endpoints
      reverse = FALSE
    )
  ) |>
  data_color( # Update cell colors, testing different color palettes
    columns = c(mean_AbsSHAP_OctApr_prcpSumCDMSum, mean_AbsSHAP_OctApr_tmeanCDMSum, mean_AbsSHAP_elevation, mean_AbsSHAP_slope, mean_AbsSHAP_aspect, mean_AbsSHAP_landcover_triclass), # ...for dose column
    direction = "row",
    fn = scales::col_numeric( # <- bc it's numeric
      palette = brewer.pal(11, "PRGn"), # A color scheme (gradient)
      domain = c(), # Column scale endpoints
      reverse = FALSE
    )
  ) |>
  # data_color( # Update cell colors, testing different color palettes
  #   columns = c(mean_AbsSHAP_OctApr_tmeanCDMSum), # ...for dose column
  #   fn = scales::col_numeric( # <- bc it's numeric
  #     palette = brewer.pal(11, "PRGn"), # A color scheme (gradient)
  #     domain = c(), # Column scale endpoints
  #     reverse = FALSE
  #   )
  # ) |>
  # data_color( # Update cell colors, testing different color palettes
  #   columns = c(mean_AbsSHAP_elevation), # ...for dose column
  #   fn = scales::col_numeric( # <- bc it's numeric
  #     palette = brewer.pal(11, "PRGn"), # A color scheme (gradient)
  #     domain = c(), # Column scale endpoints
  #     reverse = FALSE
  #   )
  # ) |>
  # data_color( # Update cell colors, testing different color palettes
  #   columns = c(mean_AbsSHAP_slope), # ...for dose column
  #   fn = scales::col_numeric( # <- bc it's numeric
  #     palette = brewer.pal(11, "PRGn"), # A color scheme (gradient)
  #     domain = c(), # Column scale endpoints
  #     reverse = FALSE
  #   )
  # ) |>
  # data_color( # Update cell colors, testing different color palettes
  #   columns = c(mean_AbsSHAP_aspect), # ...for dose column
  #   fn = scales::col_numeric( # <- bc it's numeric
  #     palette = brewer.pal(11, "PRGn"), # A color scheme (gradient)
  #     domain = c(), # Column scale endpoints
  #     reverse = FALSE
  #   )
  # ) |>
  # data_color( # Update cell colors, testing different color palettes
  #   columns = c(mean_AbsSHAP_landcover_triclass), # ...for dose column
  #   fn = scales::col_numeric( # <- bc it's numeric
  #     palette = brewer.pal(11, "PRGn"), # A color scheme (gradient)
  #     domain = c(), # Column scale endpoints
  #     reverse = FALSE
  #   )
  # ) |>
  cols_label(class_cluster = "Class", median_swe_max = "Median Max SWE (mm)", mean_AbsSHAP_OctApr_prcpSumCDMSum = "Oct-Apr prcp CDM Sum", mean_AbsSHAP_OctApr_tmeanCDMSum = "Oct-Apr tmean CDM Sum", mean_AbsSHAP_elevation = "Elevation", mean_AbsSHAP_slope = "Slope", mean_AbsSHAP_aspect = "Aspect", mean_AbsSHAP_landcover_triclass = "Landcover") |> # Update labels
  cols_move_to_end(columns = c("mean_AbsSHAP_elevation", "mean_AbsSHAP_slope", "mean_AbsSHAP_aspect", "mean_AbsSHAP_landcover_triclass")) |>
  cols_hide(c(mean_SHAPOctApr_prcpSumCDMSum, mean_SHAPOctApr_tmeanCDMSum, mean_SHAPelevation, mean_SHAPslope, mean_SHAPaspect, mean_SHAPlandcover_triclass, mean_swe_max, median_SHAPOctApr_prcpSumCDMSum, median_SHAPOctApr_tmeanCDMSum, median_SHAPelevation, median_SHAPslope, median_SHAPaspect, median_SHAPlandcover_triclass, median_AbsSHAP_OctApr_prcpSumCDMSum, median_AbsSHAP_OctApr_tmeanCDMSum, median_AbsSHAP_elevation, median_AbsSHAP_slope, median_AbsSHAP_aspect, median_AbsSHAP_landcover_triclass)) |>
  opt_vertical_padding(scale = 2)
# tab_footnote(
#   footnote = ""
# )
CONUSCluster_AvgTitle_gt
### saving table
CONUSCluster_AvgTitle_gt |>
  gtsave(
    "CONUSClusterAvg_Title.png", expand = 5,
    path = here("Outputs", "ChartsGraphsTables", "Kmeans")
  )


### making table of cluster averages
CONUSCluster_SHAPAvg_gt <- CONUSCluster_Avg |>
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
    columns = c(mean_SHAPOctApr_prcpSumCDMSum), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # A column (numeric data)
    columns = c(mean_SHAPOctApr_tmeanCDMSum), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # Another column (also numeric data)
    columns = c(mean_SHAPelevation), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
  ) |>
  fmt_number( # Another column (also numeric data)
    columns = c(mean_SHAPslope), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
  ) |>
  fmt_number( # Another column (also numeric data)
    columns = c(mean_SHAPaspect), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
  ) |>
  fmt_number( # Another column (also numeric data)
    columns = c(mean_SHAPlandcover_triclass), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
  ) |>
  data_color( # Update cell colors, testing different color palettes
    columns = c(median_swe_max), # ...for dose column
    fn = scales::col_numeric( # <- bc it's numeric
      palette = brewer.pal(11, "RdBu"), # A color scheme (gradient)
      domain = c(), # Column scale endpoints
      reverse = FALSE
    )
  ) |>
  data_color( # Update cell colors, testing different color palettes
    columns = c(mean_SHAPOctApr_prcpSumCDMSum, mean_SHAPOctApr_tmeanCDMSum, mean_SHAPelevation, mean_SHAPslope, mean_SHAPaspect, mean_SHAPlandcover_triclass), # ...for dose column
    direction = "row",
    fn = scales::col_numeric( # <- bc it's numeric
      palette = brewer.pal(11, "PRGn"), # A color scheme (gradient)
      domain = c(), # Column scale endpoints
      reverse = FALSE
    )
  ) |>
  # data_color( # Update cell colors, testing different color palettes
  #   columns = c(mean_AbsSHAP_OctApr_tmeanCDMSum), # ...for dose column
  #   fn = scales::col_numeric( # <- bc it's numeric
  #     palette = brewer.pal(11, "PRGn"), # A color scheme (gradient)
  #     domain = c(), # Column scale endpoints
  #     reverse = FALSE
  #   )
  # ) |>
  # data_color( # Update cell colors, testing different color palettes
  #   columns = c(mean_AbsSHAP_elevation), # ...for dose column
  #   fn = scales::col_numeric( # <- bc it's numeric
  #     palette = brewer.pal(11, "PRGn"), # A color scheme (gradient)
  #     domain = c(), # Column scale endpoints
  #     reverse = FALSE
  #   )
  # ) |>
  # data_color( # Update cell colors, testing different color palettes
  #   columns = c(mean_AbsSHAP_slope), # ...for dose column
  #   fn = scales::col_numeric( # <- bc it's numeric
  #     palette = brewer.pal(11, "PRGn"), # A color scheme (gradient)
  #     domain = c(), # Column scale endpoints
  #     reverse = FALSE
  #   )
  # ) |>
  # data_color( # Update cell colors, testing different color palettes
  #   columns = c(mean_AbsSHAP_aspect), # ...for dose column
  #   fn = scales::col_numeric( # <- bc it's numeric
  #     palette = brewer.pal(11, "PRGn"), # A color scheme (gradient)
  #     domain = c(), # Column scale endpoints
  #     reverse = FALSE
  #   )
  # ) |>
  # data_color( # Update cell colors, testing different color palettes
  #   columns = c(mean_AbsSHAP_landcover_triclass), # ...for dose column
  #   fn = scales::col_numeric( # <- bc it's numeric
  #     palette = brewer.pal(11, "PRGn"), # A color scheme (gradient)
  #     domain = c(), # Column scale endpoints
  #     reverse = FALSE
  #   )
  # ) |>
  cols_label(class_cluster = "Class", median_swe_max = "Median Max SWE (mm)", mean_SHAPOctApr_prcpSumCDMSum = "Oct-Apr prcp CDM Sum", mean_SHAPOctApr_tmeanCDMSum = "Oct-Apr tmean CDM Sum", mean_SHAPelevation = "Elevation", mean_SHAPslope = "Slope", mean_SHAPaspect = "Aspect", mean_SHAPlandcover_triclass = "Landcover") |> # Update labels
  cols_move_to_end(columns = c("mean_SHAPelevation", "mean_SHAPslope", "mean_SHAPaspect", "mean_SHAPlandcover_triclass")) |>
  cols_hide(c(mean_AbsSHAP_OctApr_prcpSumCDMSum, mean_AbsSHAP_OctApr_tmeanCDMSum, mean_AbsSHAP_elevation, mean_AbsSHAP_slope, mean_AbsSHAP_aspect, mean_AbsSHAP_landcover_triclass, mean_swe_max, median_SHAPOctApr_prcpSumCDMSum, median_SHAPOctApr_tmeanCDMSum, median_SHAPelevation, median_SHAPslope, median_SHAPaspect, median_SHAPlandcover_triclass, median_AbsSHAP_OctApr_prcpSumCDMSum, median_AbsSHAP_OctApr_tmeanCDMSum, median_AbsSHAP_elevation, median_AbsSHAP_slope, median_AbsSHAP_aspect, median_AbsSHAP_landcover_triclass)) |>
  opt_vertical_padding(scale = 2)
# tab_footnote(
#   footnote = ""
# )
CONUSCluster_SHAPAvg_gt
### saving table
CONUSCluster_SHAPAvg_gt |>
  gtsave(
    "CONUSClusterSHAPAvg.png", expand = 5,
    path = here("Outputs", "ChartsGraphsTables", "Kmeans")
  )

### making table of cluster averages
CONUSCluster_SHAPAvgTitle_gt <- CONUSCluster_Avg |>
  gt() |> # use 'gt' to make an awesome table...
  gt_theme_guardian() |>
  tab_header(
    title = "SWE Class SHAP Value Averages", # ...with this title
    subtitle = "Western US")  |>  # and this subtitle
  ##tab_style(style = cell_fill("bisque"),
  ##        locations = cells_body()) |>  # add fill color to table
  opt_align_table_header(align = "center") |>
  fmt_number( # A column (numeric data)
    columns = c(median_swe_max), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # A column (numeric data)
    columns = c(mean_SHAPOctApr_prcpSumCDMSum), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # A column (numeric data)
    columns = c(mean_SHAPOctApr_tmeanCDMSum), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # Another column (also numeric data)
    columns = c(mean_SHAPelevation), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
  ) |>
  fmt_number( # Another column (also numeric data)
    columns = c(mean_SHAPslope), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
  ) |>
  fmt_number( # Another column (also numeric data)
    columns = c(mean_SHAPaspect), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
  ) |>
  fmt_number( # Another column (also numeric data)
    columns = c(mean_SHAPlandcover_triclass), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
  ) |>
  data_color( # Update cell colors, testing different color palettes
    columns = c(median_swe_max), # ...for dose column
    fn = scales::col_numeric( # <- bc it's numeric
      palette = brewer.pal(11, "RdBu"), # A color scheme (gradient)
      domain = c(), # Column scale endpoints
      reverse = FALSE
    )
  ) |>
  data_color( # Update cell colors, testing different color palettes
    columns = c(mean_SHAPOctApr_prcpSumCDMSum, mean_SHAPOctApr_tmeanCDMSum, mean_SHAPelevation, mean_SHAPslope, mean_SHAPaspect, mean_SHAPlandcover_triclass),
    direction = "row",
    fn = scales::col_numeric( # <- bc it's numeric
      palette = brewer.pal(11, "PRGn"), # A color scheme (gradient)
      domain = c(), # Column scale endpoints
      reverse = FALSE
    )
  ) |>
  cols_label(class_cluster = "Class", median_swe_max = "Median Max SWE (mm)", mean_SHAPOctApr_prcpSumCDMSum = "Oct-Apr prcp CDM Sum", mean_SHAPOctApr_tmeanCDMSum = "Oct-Apr tmean CDM Sum", mean_SHAPelevation = "Elevation", mean_SHAPslope = "Slope", mean_SHAPaspect = "Aspect", mean_SHAPlandcover_triclass = "Landcover") |> # Update labels
  cols_move_to_end(columns = c("mean_SHAPelevation", "mean_SHAPslope", "mean_SHAPaspect", "mean_SHAPlandcover_triclass")) |>
  cols_hide(c(mean_AbsSHAP_OctApr_prcpSumCDMSum, mean_AbsSHAP_OctApr_tmeanCDMSum, mean_AbsSHAP_elevation, mean_AbsSHAP_slope, mean_AbsSHAP_aspect, mean_AbsSHAP_landcover_triclass, mean_swe_max, median_SHAPOctApr_prcpSumCDMSum, median_SHAPOctApr_tmeanCDMSum, median_SHAPelevation, median_SHAPslope, median_SHAPaspect, median_SHAPlandcover_triclass, median_AbsSHAP_OctApr_prcpSumCDMSum, median_AbsSHAP_OctApr_tmeanCDMSum, median_AbsSHAP_elevation, median_AbsSHAP_slope, median_AbsSHAP_aspect, median_AbsSHAP_landcover_triclass)) |>
  opt_vertical_padding(scale = 2)
# tab_footnote(
#   footnote = ""
# )
CONUSCluster_SHAPAvgTitle_gt
### saving table
CONUSCluster_SHAPAvgTitle_gt |>
  gtsave(
    "CONUSClusterSHAPAvg_Title.png", expand = 5,
    path = here("Outputs", "ChartsGraphsTables", "Kmeans")
  )


#rm()

##### grouping predictors by cluster so I can make a table for them too #####
CONUS_Predictors_All <- fread(here("Data", "SHAP", "Fastshap", "SHAPs", "CONUS", "AllYearsCombined", "CONUSAllPredictors_KClasses.csv"))

### define function to calculate mode
find_mode <- function(x) {
  u <- unique(x)
  tab <- tabulate(match(x, u))
  u[tab == max(tab)]
}

CONUSCluster_PredAvg <- CONUS_Predictors_All |>
  group_by(class_cluster) |>
  summarise(mean_swe_max = mean(peak_swe),
            median_swe_max = median(peak_swe),
            mean_OctApr_prcpSumCDMSum = mean(OctApr_prcpSumCDMSum), 
            mean_OctApr_tmeanCDMSum = mean(OctApr_tmeanCDMSum),
            mean_elevation = mean(elevation), 
            mean_slope = mean(slope),
            mean_aspect = mean(aspect), 
            mean_landcover_triclass = find_mode(landcover_triclass),
            median_OctApr_prcpSumCDMSum = median(OctApr_prcpSumCDMSum), 
            median_OctApr_tmeanCDMSum = median(OctApr_tmeanCDMSum),
            median_elevation = median(elevation), 
            median_slope = median(slope),
            median_aspect = median(aspect), 
            median_landcover_triclass = median(landcover_triclass))

### writing group averages as csv
fwrite(CONUSCluster_PredAvg, here("Data", "SHAP", "Classes", "CONUS", "KMeans", "CSVs", "Summary", "CONUSClassPredAvgs.csv"))

### making table of cluster averages
CONUSCluster_PredAvg_gt <- CONUSCluster_PredAvg |>
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
    columns = c(mean_OctApr_prcpSumCDMSum), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # A column (numeric data)
    columns = c(mean_OctApr_tmeanCDMSum), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # Another column (also numeric data)
    columns = c(mean_elevation), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
  ) |>
  fmt_number( # Another column (also numeric data)
    columns = c(mean_slope), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
  ) |>
  fmt_number( # Another column (also numeric data)
    columns = c(mean_aspect), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
  ) |>
  fmt_number( # Another column (also numeric data)
    columns = c(mean_landcover_triclass), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 0 # I want this column to have zero decimal places
  ) |>
  data_color( # Update cell colors, testing different color palettes
    columns = c(median_swe_max), # ...for dose column
    fn = scales::col_numeric( # <- bc it's numeric
      palette = brewer.pal(11, "RdBu"), # A color scheme (gradient)
      domain = c(), # Column scale endpoints
      reverse = FALSE
    )
  ) |>
  data_color( # Update cell colors, testing different color palettes
    columns = c(mean_OctApr_prcpSumCDMSum, mean_OctApr_tmeanCDMSum, mean_elevation, mean_slope, mean_aspect, mean_landcover_triclass), # ...for dose column
    # direction = "row",
    fn = scales::col_numeric( # <- bc it's numeric
      palette = brewer.pal(11, "PRGn"), # A color scheme (gradient)
      domain = c(), # Column scale endpoints
      reverse = FALSE
    )
  ) |>
  cols_label(class_cluster = "Class", median_swe_max = "Median Max SWE (mm)", mean_OctApr_prcpSumCDMSum = "Oct-Apr prcp CDM Sum", mean_OctApr_tmeanCDMSum = "Oct-Apr tmean CDM Sum", mean_elevation = "Elevation (m)", mean_slope = "Slope", mean_aspect = "Aspect", mean_landcover_triclass = "Landcover") |> # Update labels
  cols_move_to_end(columns = c("mean_elevation", "mean_slope", "mean_aspect", "mean_landcover_triclass")) |>
  cols_hide(c(mean_swe_max)) |>
  opt_vertical_padding(scale = 2) #|>
# tab_footnote(
#   footnote = "Landcover values are modes, all other values are means"
# )
CONUSCluster_PredAvg_gt
### saving table
CONUSCluster_PredAvg_gt |>
  gtsave(
    "CONUSClusterPredAvg.png", expand = 5,
    path = here("Outputs", "ChartsGraphsTables", "Kmeans")
  )

### making table of cluster averages
CONUSCluster_PredAvgTitle_gt <- CONUSCluster_PredAvg |>
  gt() |> # use 'gt' to make an awesome table...
  gt_theme_guardian() |>
  tab_header(
    title = "SWE Class Variable Averages", # ...with this title
    subtitle = "Western US")  |>  # and this subtitle
  ##tab_style(style = cell_fill("bisque"),
  ##        locations = cells_body()) |>  # add fill color to table
  opt_align_table_header(align = "center") |>
  fmt_number( # A column (numeric data)
    columns = c(median_swe_max), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # A column (numeric data)
    columns = c(mean_OctApr_prcpSumCDMSum), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # A column (numeric data)
    columns = c(mean_OctApr_tmeanCDMSum), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # Another column (also numeric data)
    columns = c(mean_elevation), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
  ) |>
  fmt_number( # Another column (also numeric data)
    columns = c(mean_slope), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
  ) |>
  fmt_number( # Another column (also numeric data)
    columns = c(mean_aspect), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
  ) |>
  fmt_number( # Another column (also numeric data)
    columns = c(mean_landcover_triclass), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 0 # I want this column to have zero decimal places
  ) |>
  data_color( # Update cell colors, testing different color palettes
    columns = c(median_swe_max), # ...for dose column
    fn = scales::col_numeric( # <- bc it's numeric
      palette = brewer.pal(11, "RdBu"), # A color scheme (gradient)
      domain = c(), # Column scale endpoints
      reverse = FALSE
    )
  ) |>
  data_color( # Update cell colors, testing different color palettes
    columns = c(mean_OctApr_prcpSumCDMSum, mean_OctApr_tmeanCDMSum, mean_elevation, mean_slope, mean_aspect, mean_landcover_triclass), # ...for dose column
    # direction = "row",
    fn = scales::col_numeric( # <- bc it's numeric
      palette = brewer.pal(11, "PRGn"), # A color scheme (gradient)
      domain = c(), # Column scale endpoints
      reverse = FALSE
    )
  ) |>
  cols_label(class_cluster = "Class", median_swe_max = "Median Max SWE (mm)", mean_OctApr_prcpSumCDMSum = "Oct-Apr prcp CDM Sum", mean_OctApr_tmeanCDMSum = "Oct-Apr tmean CDM Sum", mean_elevation = "Elevation (m)", mean_slope = "Slope", mean_aspect = "Aspect", mean_landcover_triclass = "Landcover") |> # Update labels
  cols_move_to_end(columns = c("mean_elevation", "mean_slope", "mean_aspect", "mean_landcover_triclass")) |>
  cols_hide(c(mean_swe_max)) |>
  opt_vertical_padding(scale = 2) |>
  tab_footnote(
    footnote = "Landcover values are modes, all other values are means except for Max SWE"
    )
CONUSCluster_PredAvgTitle_gt
### saving table
CONUSCluster_PredAvgTitle_gt |>
  gtsave(
    "CONUSClusterPredAvg_Title.png", expand = 5,
    path = here("Outputs", "ChartsGraphsTables", "Kmeans")
  )

##### Making stacked bar chart showing each covariate's percentage of total importance for each class #####
### reading group averages
CONUSCluster_Avg <- read_csv(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "CSVs", "Summary", "CONUSClassSHAPAvgs.csv")) |>
  mutate(total_absshap = mean_AbsSHAP_OctApr_prcpSumCDMSum + mean_AbsSHAP_OctApr_tmeanCDMSum + mean_AbsSHAP_elevation + mean_AbsSHAP_slope + mean_AbsSHAP_aspect + mean_AbsSHAP_landcover_triclass,
         pct_AbsSHAP_OctApr_prcpSumCDMSum = mean_AbsSHAP_OctApr_prcpSumCDMSum / total_absshap * 100, 
         pct_AbsSHAP_OctApr_tmeanCDMSum = mean_AbsSHAP_OctApr_tmeanCDMSum / total_absshap * 100,
         pct_AbsSHAP_elevation = mean_AbsSHAP_elevation / total_absshap * 100,
         pct_AbsSHAP_slope = mean_AbsSHAP_slope / total_absshap * 100,
         pct_AbsSHAP_aspect = mean_AbsSHAP_aspect / total_absshap * 100,
         pct_AbsSHAP_landcover_triclass = mean_AbsSHAP_landcover_triclass / total_absshap * 100) |>
  select(class_cluster, total_absshap, pct_AbsSHAP_OctApr_prcpSumCDMSum, pct_AbsSHAP_OctApr_tmeanCDMSum, pct_AbsSHAP_elevation, pct_AbsSHAP_slope, pct_AbsSHAP_aspect, pct_AbsSHAP_landcover_triclass)

### ggplot time
# c(pct_AbsSHAP_OctApr_prcpSumCDMSum, pct_AbsSHAP_OctApr_tmeanCDMSum, pct_AbsSHAP_DecFeb_prcpSumCDMSum, pct_AbsSHAP_DecFeb_tmeanCDMSum, pct_AbsSHAP_elevation, pct_AbsSHAP_slope, pct_AbsSHAP_aspect, pct_AbsSHAP_OctApr_tminMean, pct_AbsSHAP_OctApr_tmaxMean, pct_AbsSHAP_landcover_triclass)))
CONUSCluster_Avg_Pivot = CONUSCluster_Avg |>
  pivot_longer(
    cols = starts_with("pct_"),
  )

### setting name to be as a factor so I can set what order the bars are in in the stacked bar chart below
CONUSCluster_Avg_Pivot$name = factor(CONUSCluster_Avg_Pivot$name, levels = c("pct_AbsSHAP_landcover_triclass",  "pct_AbsSHAP_aspect", "pct_AbsSHAP_slope", "pct_AbsSHAP_elevation", "pct_AbsSHAP_OctApr_tmeanCDMSum", "pct_AbsSHAP_OctApr_prcpSumCDMSum"))

PctAbsSHAP_plot <- ggplot(data = CONUSCluster_Avg_Pivot, mapping = aes(x = class_cluster, y = value, fill = name)) +
  theme_bw() +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_brewer("", palette = "PRGn", labels = c("Landcover", "Aspect", "Slope", "Elevation", "Oct-Apr Tmean CDM Sum", "Oct-Apr prcp CDM Sum")) +
  ggtitle(label = "", subtitle = "Western US") +
  xlab("Class") +
  ylab("% of Absolute SHAP Value Sum") +
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9)) +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

PctAbsSHAP_plot
ggsave(plot = PctAbsSHAP_plot, here("Outputs", "ChartsGraphsTables", "SHAPPcts", "CONUSPctAbsSHAP.png"))

### now with a title
PctAbsSHAPTitle_plot <- ggplot(data = CONUSCluster_Avg_Pivot, mapping = aes(x = class_cluster, y = value, fill = name)) +
  theme_bw() +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_brewer("", palette = "PRGn", labels = c("Landcover", "Aspect", "Slope", "Elevation", "Oct-Apr Tmean CDM Sum", "Oct-Apr prcp CDM Sum")) +
  ggtitle(label = "Predictor Proportion of Total Absolute SHAP Value by Class", subtitle = "Western US") +
  xlab("Class") +
  ylab("% of Absolute SHAP Value Sum") +
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9)) +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

PctAbsSHAPTitle_plot
ggsave(plot = PctAbsSHAPTitle_plot, here("Outputs", "ChartsGraphsTables", "SHAPPcts", "CONUSPctAbsSHAP_Title.png"))
