##### K means classification based on fastshap #####
totalstarttime <- Sys.time()
##### Loading packages #####
library(pacman)
p_load(here, tidyverse, sf, terra, sfext, parallel, data.table, mcprogress, dtplyr, gt, gtExtras, RColorBrewer)
# options(mc.cores = parallel::detectCores())
# setwd("/media/Research/Morafkai/GriffinS")

### reading in csvs
# AKSHAP_csvs <- list.files(path = here("Data", "SHAP", "Fastshap", "SHAPs", "Alaska"), pattern = ".csv", full.names = TRUE)
# 
years <- 1993:2020
# 
# ReadAllYears_func <- function(year){
#   temp_shaps <- fread(file = paste0(here("Data", "SHAP", "Fastshap", "SHAPs", "Alaska"), paste0("Alaska", year, "SHAPs.csv")))
# }
# AKSHAPs <- pmclapply(X = years, FUN = ReadAllYears_func, mc.cores = length(years))
# 
ReadAllYearsPreds_func <- function(year){
  temp_preds <- fread(file = here("Data", "SHAP", "PredictorCSVs", "Alaska", paste0("AK_SWEmaxPredictors", year, ".csv")))
}
### reading in csvs of predictor rasters
AKPredictors <- pmclapply(X = years, FUN = ReadAllYearsPreds_func, mc.cores = length(years))

### adding water year column so I can bind it with max SWE values and class clusters later
for (x in 1:length(AKPredictors)){
  temp_dt <- AKPredictors[[x]]
  year <- years[x]
  temp_dt <- temp_dt |>
    mutate(WaterYear = year,
           pastexyWY = paste0(x, y, WaterYear))
  if (x == 1){
    AK_Predictors_All <- temp_dt
  } else{
    AK_Predictors_All <- rbind(AK_Predictors_All, temp_dt)
  }
}

 
# for (x in 1:length(AKSHAPs)){
#   temp_dt <- AKSHAPs[[x]]
#   year <- years[x]
#   temp_dt <- temp_dt |>
#     mutate(WaterYear = year)
#   if (x == 1){
#     AK_SHAPs_All <- temp_dt
#   } else{
#     AK_SHAPs_All <- rbind(AK_SHAPs_All, temp_dt)
#   }
# }
# 
# 
# AK_SHAPs_All <- AK_SHAPs_All |>
#   mutate(AbsSHAP_OctMay_prcpSumCDMSum = abs(OctMay_prcpSumCDMSum),
#          AbsSHAP_OctMay_tmeanCDMSum = abs(OctMay_tmeanCDMSum),
#          AbsSHAP_DecFeb_prcpSumCDMSum = abs(DecFeb_prcpSumCDMSum),
#          AbsSHAP_DecFeb_tmeanCDMSum = abs(DecFeb_tmeanCDMSum),
#          AbsSHAP_elevation = abs(elevation),
#          AbsSHAP_slope = abs(slope),
#          AbsSHAP_aspect = abs(aspect),
#          AbsSHAP_OctMay_tminMean = abs(OctMay_tminMean),
#          AbsSHAP_OctMay_tmaxMean = abs(OctMay_tmaxMean),
#          AbsSHAP_landcover_triclass = abs(landcover_triclass))
# 
# fwrite(AK_SHAPs_All, file = here("Data", "SHAP", "Fastshap", "SHAPs", "Alaska", "AllYearsCombined", "AlaskaAllSHAPs.csv"))
# AK_SHAPs_All <- fread(here("Data", "SHAP", "Fastshap", "SHAPs", "Alaska", "AllYearsCombined", "AlaskaAllSHAPs.csv"))


##### Kmeans Clustering Alaska SHAP absolute values #####
### initializing df to store clustering error metrics
AKKmeansErrors_df <- data.frame(k = 1:25,
                              sse = numeric(25))

AKKmeans_func <- function(i){
  set.seed(802)
  kmeans_error <- kmeans(AK_SHAPs_All[, .(peak_swe, AbsSHAP_OctMay_prcpSumCDMSum, AbsSHAP_OctMay_tmeanCDMSum, AbsSHAP_DecFeb_prcpSumCDMSum, AbsSHAP_DecFeb_tmeanCDMSum, AbsSHAP_elevation, AbsSHAP_slope, AbsSHAP_aspect, AbsSHAP_OctMay_tminMean, AbsSHAP_OctMay_tmaxMean, AbsSHAP_landcover_triclass)], centers = i)$tot.withinss
  return(kmeans_error)
}

AKKmeans_errors <- pmclapply(X = AKKmeansErrors_df$k, FUN = AKKmeans_func, mc.cores = 8, mc.silent = TRUE, mc.set.seed = FALSE)

AKKmeansErrors_df$sse <- unlist(AKKmeans_errors)

AKKmeansErrors_df$pct_error_change <- numeric(25)
for (i in 1:25){
  if (i == 1){
    AKKmeansErrors_df$pct_error_change[i] = 1
  } else{
    AKKmeansErrors_df$pct_error_change[i] = (AKKmeansErrors_df$sse[i-1] - AKKmeansErrors_df$sse[i]) / AKKmeansErrors_df$sse[i - 1]
  }
}

AKKmeansError_plot <- ggplot(data = AKKmeansErrors_df, mapping = aes(x = k, y = sse)) +
  theme_bw() +
  geom_point() +
  geom_line() +
  # geom_smooth(method = "gam") +
  # ggtitle(label = "Total Sum of Squared Error for Different Numbers of K-Means Clusters", subtitle = "WUS Study Area") +
  ggtitle(label = "Alaska Study Area") +
  xlab("K Number of Clusters") +
  ylab("Total SSE") +
  theme(plot.title = element_text(hjust = 0.5, size = 12), plot.subtitle = element_text(hjust = 0.5))
AKKmeansError_plot
ggsave(plot = AKKmeansError_plot, here("Outputs", "Plots", "KmeansError", "AKKmeansSSEPlot.png"))

##### Final KMeans with appropriate number of clusters #####
set.seed(802)
FinalKmeans <- kmeans(AK_SHAPs_All[, .(peak_swe, AbsSHAP_OctMay_prcpSumCDMSum, AbsSHAP_OctMay_tmeanCDMSum, AbsSHAP_DecFeb_prcpSumCDMSum, AbsSHAP_DecFeb_tmeanCDMSum, AbsSHAP_elevation, AbsSHAP_slope, AbsSHAP_aspect, AbsSHAP_OctMay_tminMean, AbsSHAP_OctMay_tmaxMean, AbsSHAP_landcover_triclass)], centers = 5)

AK_SHAPs_All <- AK_SHAPs_All |>
  mutate(class_cluster = FinalKmeans$cluster)

### Calculating Averages by Cluster so I can re-assign cluster numbers (not recluster) by descending order of average max swe
AKCluster_Avg <- AK_SHAPs_All |>
  group_by(class_cluster) |>
  summarise(mean_swe_max = mean(peak_swe)) |>
  mutate(swe_max_rank = dense_rank(desc(mean_swe_max)))

### reassigning class number based on descending order of average max swe for the classes
AK_SHAPs_All <- AK_SHAPs_All |>
  mutate(class_cluster2 = case_when(class_cluster == 1 ~ AKCluster_Avg$swe_max_rank[AKCluster_Avg$class_cluster == 1],
                                    class_cluster == 2 ~ AKCluster_Avg$swe_max_rank[AKCluster_Avg$class_cluster == 2],
                                    class_cluster == 3 ~ AKCluster_Avg$swe_max_rank[AKCluster_Avg$class_cluster == 3],
                                    class_cluster == 4 ~ AKCluster_Avg$swe_max_rank[AKCluster_Avg$class_cluster == 4],
                                    class_cluster == 5 ~ AKCluster_Avg$swe_max_rank[AKCluster_Avg$class_cluster == 5]),
         class_cluster = class_cluster2,
         pastexyWY = paste0(x, y, WaterYear)) |>
  ### dropping class_cluster2 since it was only supposed to be temporary
  select(-one_of("class_cluster2"))

### dropping xy from AK predictors since it's already in the SHAPs df
AK_Predictors_All <- AK_Predictors_All |>
  select(OctMay_prcpSumCDMSum, OctMay_tmeanCDMSum, DecFeb_prcpSumCDMSum, DecFeb_tmeanCDMSum, elevation, slope, aspect, OctMay_tminMean, OctMay_tmaxMean, landcover_triclass, WaterYear, pastexyWY)

### binding SHAPs to predictors so I can make a table later
## changing colnames first since fastshap calls the columns containing the SHAP values for each variable by just the name of the variable, which is kind of inconvenient
colnames(AK_SHAPs_All) <- c("x", "y", "peak_swe", "SHAPOctMay_prcpSumCDMSum", "SHAPOctMay_tmeanCDMSum", "SHAPDecFeb_prcpSumCDMSum", "SHAPDecFeb_tmeanCDMSum", "SHAPelevation", "SHAPslope", "SHAPaspect", "SHAPOctMay_tminMean", "SHAPOctMay_tmaxMean", "SHAPlandcover_triclass", "Rank_OctMayprcpSumCDM_varimp", "Rank_OctMaytmeanCDM_varimp", "Rank_DecFebprcpSumCDM_varimp", "Rank_DecFebtmeanCDM_varimp", "Rank_elevation_varimp", "Rank_slope_varimp", "Rank_aspect_varimp", "Rank_OctMaytminMean_varimp", "Rank_OctMaytmaxMean_varimp", "Rank_landcover_varimp", "class", "class_top1", "class_top2", "class_top3", "class_top4", "class_top5", "WaterYear", "AbsSHAP_OctMay_prcpSumCDMSum", "AbsSHAP_OctMay_tmeanCDMSum", "AbsSHAP_DecFeb_prcpSumCDMSum", "AbsSHAP_DecFeb_tmeanCDMSum", "AbsSHAP_elevation", "AbsSHAP_slope", "AbsSHAP_aspect", "AbsSHAP_OctMay_tminMean", "AbsSHAP_OctMay_tmaxMean", "AbsSHAP_landcover_triclass", "class_cluster", "pastexyWY")


# AK_SHAPs_All <- full_join(AK_SHAPs_All, AK_Predictors_All, by = "pastexyWY")
AK_Predictors_All <- full_join(AK_SHAPs_All[, .(x, y, peak_swe, class_cluster, pastexyWY)], AK_Predictors_All, by = "pastexyWY") |>
  select(-one_of("pastexyWY"))

### selecting only columns that I'm actually going to use later
AK_SHAPs_All <- AK_SHAPs_All |>
  select(x, y, peak_swe, SHAPOctMay_prcpSumCDMSum, SHAPOctMay_tmeanCDMSum, SHAPDecFeb_prcpSumCDMSum, SHAPDecFeb_tmeanCDMSum, SHAPelevation, SHAPslope, SHAPaspect, SHAPOctMay_tminMean, SHAPOctMay_tmaxMean, SHAPlandcover_triclass, Rank_OctMayprcpSumCDM_varimp, Rank_OctMaytmeanCDM_varimp, Rank_DecFebprcpSumCDM_varimp, Rank_DecFebtmeanCDM_varimp, Rank_elevation_varimp, Rank_slope_varimp, Rank_aspect_varimp, Rank_OctMaytminMean_varimp, Rank_OctMaytmaxMean_varimp, Rank_landcover_varimp, WaterYear, AbsSHAP_OctMay_prcpSumCDMSum, AbsSHAP_OctMay_tmeanCDMSum, AbsSHAP_DecFeb_prcpSumCDMSum, AbsSHAP_DecFeb_tmeanCDMSum, AbsSHAP_elevation, AbsSHAP_slope, AbsSHAP_aspect, AbsSHAP_OctMay_tminMean, AbsSHAP_OctMay_tmaxMean, AbsSHAP_landcover_triclass, class_cluster)


### writing csvs of this because I don't want to have to do this again
fwrite(AK_SHAPs_All, file = here("Data", "SHAP", "Fastshap", "SHAPs", "Alaska", "AllYearsCombined", "AlaskaAllSHAPs_KClasses.csv"))
fwrite(AK_Predictors_All, file = here("Data", "SHAP", "Fastshap", "SHAPs", "Alaska", "AllYearsCombined", "AKAllPredictors_KClasses.csv"))
rm()

##### Rasterizing K-means-based classes #####
AK_SHAPs_All <- fread(here("Data", "SHAP", "Fastshap", "SHAPs", "Alaska", "AllYearsCombined", "AlaskaAllSHAPs_KClasses.csv"))
AK_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "Alaska", "AK_AOI.gpkg"))
years <- 1993:2020
SaveClasses_func <- function(year){
  PredictionRast <- rast(here("Outputs", "MapTIFFs", "Alaska", paste0("AK_SWEmax", year, ".tif")))
  temp_classes <- AK_SHAPs_All |>
    filter(WaterYear == year)
  fwrite(temp_classes, file = here("Data", "SHAP", "Classes", "Alaska", "KMeans", "CSVs", paste0("AK", year, "_KClasses.csv")))
  ### converting year's classes to raster
  temp_classes_rast <- terra::rasterize(vect(st_as_sf(temp_classes, coords = c("x", "y"), crs = crs(AK_AOI))), y = PredictionRast, field = "class_cluster")
  writeRaster(temp_classes_rast, here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", paste0("AK", year, "_KClasses.tif")), overwrite = TRUE)
  # return(temp_classes_rast)
}

SWEClassesByYear <- pmclapply(X = years, FUN = SaveClasses_func, mc.cores = 10, mc.silent = TRUE, mc.set.seed = FALSE)

### reading in class rasters
AK1993Classes_rast <- rast(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "AK1993_KClasses.tif"))
AK1994Classes_rast <- rast(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "AK1994_KClasses.tif"))
AK1995Classes_rast <- rast(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "AK1995_KClasses.tif"))
AK1996Classes_rast <- rast(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "AK1996_KClasses.tif"))
# plot(AK1996Classes_rast)
AK1997Classes_rast <- rast(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "AK1997_KClasses.tif"))
# plot(AK1997Classes_rast)
AK1998Classes_rast <- rast(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "AK1998_KClasses.tif"))
# plot(AK1998Classes_rast)
AK1999Classes_rast <- rast(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "AK1999_KClasses.tif"))
# plot(AK1999Classes_rast)
AK2000Classes_rast <- rast(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "AK2000_KClasses.tif"))
# plot(AK2000Classes_rast)
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

AKAllYearsClasses_rast <- c(AK1993Classes_rast, AK1994Classes_rast, AK1995Classes_rast, AK1996Classes_rast, AK1997Classes_rast, AK1998Classes_rast, AK1999Classes_rast, AK2000Classes_rast, AK2001Classes_rast, AK2002Classes_rast, AK2003Classes_rast, AK2004Classes_rast, AK2005Classes_rast, AK2006Classes_rast, AK2007Classes_rast, AK2008Classes_rast, AK2009Classes_rast, AK2010Classes_rast, AK2011Classes_rast, AK2012Classes_rast, AK2013Classes_rast, AK2014Classes_rast, AK2015Classes_rast, AK2016Classes_rast, AK2017Classes_rast, AK2018Classes_rast, AK2019Classes_rast, AK2020Classes_rast)

AKModalClass_rast <- app(AKAllYearsClasses_rast, fun = "modal")
AKNumUniqueClass_rast <- app(AKAllYearsClasses_rast, fun = function(x){length(unique(x))})
AKNumUniqueClass_rast <- crop(AKNumUniqueClass_rast, AK_AOI, mask = TRUE)
# AKUniqueClasses_rast <- app(AKAllYearsClasses_rast, fun = function(x){as.factor(unique(x))})
# AKUniqueClasses_rast <- crop(AKUniqueClasses_rast, AK_AOI, mask = TRUE)
plot(AKModalClass_rast, main = "most common class")
plot(st_geometry(AK_AOI), add = TRUE)
plot(AKNumUniqueClass_rast, main = "number of unique classes")
plot(st_geometry(AK_AOI), add = TRUE)

### saving most common class raster and number of unique classes raster
writeRaster(AKModalClass_rast, here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "Summary", "AKModalClass.tif"), overwrite = TRUE)
writeRaster(AKNumUniqueClass_rast, here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "Summary", "AKNumUniqueClass.tif"), overwrite = TRUE)

##### Calculating Averages by Cluster #####
AKCluster_Avg <- AK_SHAPs_All |>
  group_by(class_cluster) |>
  summarise(mean_swe_max = mean(peak_swe),
            median_swe_max = median(peak_swe),
            mean_SHAPOctMay_prcpSumCDMSum = mean(SHAPOctMay_prcpSumCDMSum), 
            mean_SHAPOctMay_tmeanCDMSum = mean(SHAPOctMay_tmeanCDMSum), 
            mean_SHAPDecFeb_prcpSumCDMSum = mean(SHAPDecFeb_prcpSumCDMSum), 
            mean_SHAPDecFeb_tmeanCDMSum = mean(SHAPDecFeb_tmeanCDMSum), 
            mean_SHAPelevation = mean(SHAPelevation), 
            mean_SHAPslope = mean(SHAPslope), 
            mean_SHAPaspect = mean(SHAPaspect), 
            mean_SHAPOctMay_tminMean = mean(SHAPOctMay_tminMean), 
            mean_SHAPOctMay_tmaxMean = mean(SHAPOctMay_tmaxMean), 
            mean_SHAPlandcover_triclass = mean(SHAPlandcover_triclass),
            mean_AbsSHAP_OctMay_prcpSumCDMSum = mean(AbsSHAP_OctMay_prcpSumCDMSum), 
            mean_AbsSHAP_OctMay_tmeanCDMSum = mean(AbsSHAP_OctMay_tmeanCDMSum), 
            mean_AbsSHAP_DecFeb_prcpSumCDMSum = mean(AbsSHAP_DecFeb_prcpSumCDMSum), 
            mean_AbsSHAP_DecFeb_tmeanCDMSum = mean(AbsSHAP_DecFeb_tmeanCDMSum), 
            mean_AbsSHAP_elevation = mean(AbsSHAP_elevation), 
            mean_AbsSHAP_slope = mean(AbsSHAP_slope), 
            mean_AbsSHAP_aspect = mean(AbsSHAP_aspect), 
            mean_AbsSHAP_OctMay_tminMean = mean(AbsSHAP_OctMay_tminMean), 
            mean_AbsSHAP_OctMay_tmaxMean = mean(AbsSHAP_OctMay_tmaxMean), 
            mean_AbsSHAP_landcover_triclass = mean(AbsSHAP_landcover_triclass),
            median_SHAPOctMay_prcpSumCDMSum = median(SHAPOctMay_prcpSumCDMSum), 
            median_SHAPOctMay_tmeanCDMSum = median(SHAPOctMay_tmeanCDMSum), 
            median_SHAPDecFeb_prcpSumCDMSum = median(SHAPDecFeb_prcpSumCDMSum), 
            median_SHAPDecFeb_tmeanCDMSum = median(SHAPDecFeb_tmeanCDMSum), 
            median_SHAPelevation = median(SHAPelevation), 
            median_SHAPslope = median(SHAPslope), 
            median_SHAPaspect = median(SHAPaspect), 
            median_SHAPOctMay_tminMean = median(SHAPOctMay_tminMean), 
            median_SHAPOctMay_tmaxMean = median(SHAPOctMay_tmaxMean), 
            median_SHAPlandcover_triclass = median(SHAPlandcover_triclass),
            median_AbsSHAP_OctMay_prcpSumCDMSum = median(AbsSHAP_OctMay_prcpSumCDMSum), 
            median_AbsSHAP_OctMay_tmeanCDMSum = median(AbsSHAP_OctMay_tmeanCDMSum), 
            median_AbsSHAP_DecFeb_prcpSumCDMSum = median(AbsSHAP_DecFeb_prcpSumCDMSum), 
            median_AbsSHAP_DecFeb_tmeanCDMSum = median(AbsSHAP_DecFeb_tmeanCDMSum), 
            median_AbsSHAP_elevation = median(AbsSHAP_elevation), 
            median_AbsSHAP_slope = median(AbsSHAP_slope), 
            median_AbsSHAP_aspect = median(AbsSHAP_aspect), 
            median_AbsSHAP_OctMay_tminMean = median(AbsSHAP_OctMay_tminMean), 
            median_AbsSHAP_OctMay_tmaxMean = median(AbsSHAP_OctMay_tmaxMean), 
            median_AbsSHAP_landcover_triclass = median(AbsSHAP_landcover_triclass))
### writing group averages as csv
fwrite(AKCluster_Avg, here("Data", "SHAP", "Classes", "Alaska", "KMeans", "CSVs", "Summary", "AKClassSHAPAvgs.csv"))

##### making table of cluster absolute SHAP value averages #####
AKCluster_Avg_gt <- AKCluster_Avg |>
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
    columns = c(mean_AbsSHAP_OctMay_prcpSumCDMSum), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # A column (numeric data)
    columns = c(mean_AbsSHAP_OctMay_tmeanCDMSum), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # A column (numeric data)
    columns = c(mean_AbsSHAP_OctMay_tminMean), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # Another column (also numeric data)
    columns = c(mean_AbsSHAP_OctMay_tmaxMean), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
  ) |> 
  fmt_number( # Another column (also numeric data)
    columns = c(mean_AbsSHAP_DecFeb_prcpSumCDMSum), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
  ) |>
  fmt_number( # Another column (also numeric data)
    columns = c(mean_AbsSHAP_DecFeb_tmeanCDMSum), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
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
    columns = c(median_swe_max),
    fn = scales::col_numeric( # <- bc it's numeric
      palette = brewer.pal(11, "RdBu"), # A color scheme (gradient)
      domain = c(), # Column scale endpoints
      reverse = FALSE
    )
  ) |>
  data_color( # Update cell colors, testing different color palettes
    columns = c(mean_AbsSHAP_OctMay_prcpSumCDMSum, mean_AbsSHAP_OctMay_tmeanCDMSum, mean_AbsSHAP_OctMay_tminMean, mean_AbsSHAP_OctMay_tmaxMean, mean_AbsSHAP_DecFeb_prcpSumCDMSum, mean_AbsSHAP_DecFeb_tmeanCDMSum, mean_AbsSHAP_elevation, mean_AbsSHAP_slope, mean_AbsSHAP_aspect, mean_AbsSHAP_landcover_triclass), # ...for dose column
    direction = "row",
    fn = scales::col_numeric( # <- bc it's numeric
      palette = brewer.pal(11, "PRGn"), # A color scheme (gradient)
      domain = c(), # Column scale endpoints
      reverse = FALSE
    )
  ) |>
  cols_label(class_cluster = "Class", median_swe_max = "Median Max SWE (mm)", mean_AbsSHAP_OctMay_prcpSumCDMSum = "Oct-May prcp CDM Sum", mean_AbsSHAP_OctMay_tmeanCDMSum = "Oct-May tmean CDM Sum", mean_AbsSHAP_OctMay_tminMean = "Oct-May tmin Mean", mean_AbsSHAP_OctMay_tmaxMean = "Oct-May tmax Mean", mean_AbsSHAP_DecFeb_prcpSumCDMSum = "Dec-Feb prcp CDM Sum", mean_AbsSHAP_DecFeb_tmeanCDMSum = "Dec-Feb tmean CDM Sum", mean_AbsSHAP_elevation = "Elevation", mean_AbsSHAP_slope = "Slope", mean_AbsSHAP_aspect = "Aspect", mean_AbsSHAP_landcover_triclass = "Landcover") |> # Update labels
  cols_move_to_end(columns = c("mean_AbsSHAP_DecFeb_prcpSumCDMSum", "mean_AbsSHAP_DecFeb_tmeanCDMSum", "mean_AbsSHAP_elevation", "mean_AbsSHAP_slope", "mean_AbsSHAP_aspect", "mean_AbsSHAP_landcover_triclass")) |>
  cols_hide(c(mean_SHAPOctMay_prcpSumCDMSum, mean_SHAPOctMay_tmeanCDMSum, mean_SHAPOctMay_tminMean, mean_SHAPOctMay_tmaxMean, mean_SHAPDecFeb_prcpSumCDMSum, mean_SHAPDecFeb_tmeanCDMSum, mean_SHAPelevation, mean_SHAPslope, mean_SHAPaspect, mean_SHAPlandcover_triclass, mean_swe_max, median_SHAPOctMay_prcpSumCDMSum, median_SHAPOctMay_tmeanCDMSum, median_SHAPDecFeb_prcpSumCDMSum, median_SHAPDecFeb_tmeanCDMSum, median_SHAPelevation, median_SHAPslope, median_SHAPaspect, median_SHAPOctMay_tminMean, median_SHAPOctMay_tmaxMean, median_SHAPlandcover_triclass, median_AbsSHAP_OctMay_prcpSumCDMSum, median_AbsSHAP_OctMay_tmeanCDMSum, median_AbsSHAP_DecFeb_prcpSumCDMSum, median_AbsSHAP_DecFeb_tmeanCDMSum, median_AbsSHAP_elevation, median_AbsSHAP_slope, median_AbsSHAP_aspect, median_AbsSHAP_OctMay_tminMean, median_AbsSHAP_OctMay_tmaxMean, median_AbsSHAP_landcover_triclass)) |>
  opt_vertical_padding(scale = 3)
  # tab_footnote(
  #   footnote = ""
  # )
AKCluster_Avg_gt
### saving table
AKCluster_Avg_gt |>
  gtsave(
    "AKClusterAvg.png", expand = 5,
    path = here("Outputs", "ChartsGraphsTables", "Kmeans")
  )

### same table but with a title for the presentation
AKCluster_AvgTitle_gt <- AKCluster_Avg |>
  gt() |> # use 'gt' to make an awesome table...
  gt_theme_guardian() |>
  tab_header(
    title = "SWE Class Absolute SHAP Value Averages", # ...with this title
    subtitle = "Alaska")  |>  # and this subtitle
  opt_align_table_header(align = "center") |>
  ##tab_style(style = cell_fill("bisque"),
  ##        locations = cells_body()) |>  # add fill color to table
  fmt_number( # A column (numeric data)
    columns = c(median_swe_max), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # A column (numeric data)
    columns = c(mean_AbsSHAP_OctMay_prcpSumCDMSum), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # A column (numeric data)
    columns = c(mean_AbsSHAP_OctMay_tmeanCDMSum), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # A column (numeric data)
    columns = c(mean_AbsSHAP_OctMay_tminMean), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # Another column (also numeric data)
    columns = c(mean_AbsSHAP_OctMay_tmaxMean), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
  ) |> 
  fmt_number( # Another column (also numeric data)
    columns = c(mean_AbsSHAP_DecFeb_prcpSumCDMSum), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
  ) |>
  fmt_number( # Another column (also numeric data)
    columns = c(mean_AbsSHAP_DecFeb_tmeanCDMSum), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
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
    columns = c(mean_AbsSHAP_OctMay_prcpSumCDMSum, mean_AbsSHAP_OctMay_tmeanCDMSum, mean_AbsSHAP_OctMay_tminMean, mean_AbsSHAP_OctMay_tmaxMean, mean_AbsSHAP_DecFeb_prcpSumCDMSum, mean_AbsSHAP_DecFeb_tmeanCDMSum, mean_AbsSHAP_elevation, mean_AbsSHAP_slope, mean_AbsSHAP_aspect, mean_AbsSHAP_landcover_triclass), # ...for dose column
    direction = "row",
    fn = scales::col_numeric( # <- bc it's numeric
      palette = brewer.pal(11, "PRGn"), # A color scheme (gradient)
      domain = c(), # Column scale endpoints
      reverse = FALSE
    )
  ) |>
  cols_label(class_cluster = "Class", median_swe_max = "Median Max SWE (mm)", mean_AbsSHAP_OctMay_prcpSumCDMSum = "Oct-May prcp CDM Sum", mean_AbsSHAP_OctMay_tmeanCDMSum = "Oct-May tmean CDM Sum", mean_AbsSHAP_OctMay_tminMean = "Oct-May tmin Mean", mean_AbsSHAP_OctMay_tmaxMean = "Oct-May tmax Mean", mean_AbsSHAP_DecFeb_prcpSumCDMSum = "Dec-Feb prcp CDM Sum", mean_AbsSHAP_DecFeb_tmeanCDMSum = "Dec-Feb tmean CDM Sum", mean_AbsSHAP_elevation = "Elevation", mean_AbsSHAP_slope = "Slope", mean_AbsSHAP_aspect = "Aspect", mean_AbsSHAP_landcover_triclass = "Landcover") |> # Update labels
  cols_move_to_end(columns = c("mean_AbsSHAP_DecFeb_prcpSumCDMSum", "mean_AbsSHAP_DecFeb_tmeanCDMSum", "mean_AbsSHAP_elevation", "mean_AbsSHAP_slope", "mean_AbsSHAP_aspect", "mean_AbsSHAP_landcover_triclass")) |>
  cols_hide(c(mean_SHAPOctMay_prcpSumCDMSum, mean_SHAPOctMay_tmeanCDMSum, mean_SHAPOctMay_tminMean, mean_SHAPOctMay_tmaxMean, mean_SHAPDecFeb_prcpSumCDMSum, mean_SHAPDecFeb_tmeanCDMSum, mean_SHAPelevation, mean_SHAPslope, mean_SHAPaspect, mean_SHAPlandcover_triclass, mean_swe_max, median_SHAPOctMay_prcpSumCDMSum, median_SHAPOctMay_tmeanCDMSum, median_SHAPDecFeb_prcpSumCDMSum, median_SHAPDecFeb_tmeanCDMSum, median_SHAPelevation, median_SHAPslope, median_SHAPaspect, median_SHAPOctMay_tminMean, median_SHAPOctMay_tmaxMean, median_SHAPlandcover_triclass, median_AbsSHAP_OctMay_prcpSumCDMSum, median_AbsSHAP_OctMay_tmeanCDMSum, median_AbsSHAP_DecFeb_prcpSumCDMSum, median_AbsSHAP_DecFeb_tmeanCDMSum, median_AbsSHAP_elevation, median_AbsSHAP_slope, median_AbsSHAP_aspect, median_AbsSHAP_OctMay_tminMean, median_AbsSHAP_OctMay_tmaxMean, median_AbsSHAP_landcover_triclass)) |>
  opt_vertical_padding(scale = 3)
# tab_footnote(
#   footnote = ""
# )
AKCluster_AvgTitle_gt
### saving table
AKCluster_AvgTitle_gt |>
  gtsave(
    "AKClusterAvgTitle.png", expand = 5,
    path = here("Outputs", "ChartsGraphsTables", "Kmeans")
  )

##### making table of cluster raw SHAP value averages #####
AKCluster_SHAPAvg_gt <- AKCluster_Avg |>
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
    columns = c(mean_SHAPOctMay_prcpSumCDMSum), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # A column (numeric data)
    columns = c(mean_SHAPOctMay_tmeanCDMSum), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # A column (numeric data)
    columns = c(mean_SHAPOctMay_tminMean), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # Another column (also numeric data)
    columns = c(mean_SHAPOctMay_tmaxMean), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
  ) |> 
  fmt_number( # Another column (also numeric data)
    columns = c(mean_SHAPDecFeb_prcpSumCDMSum), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
  ) |>
  fmt_number( # Another column (also numeric data)
    columns = c(mean_SHAPDecFeb_tmeanCDMSum), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
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
    columns = c(median_swe_max),
    fn = scales::col_numeric( # <- bc it's numeric
      palette = brewer.pal(11, "RdBu"), # A color scheme (gradient)
      domain = c(), # Column scale endpoints
      reverse = FALSE
    )
  ) |>
  data_color( # Update cell colors, testing different color palettes
    columns = c(mean_SHAPOctMay_prcpSumCDMSum, mean_SHAPOctMay_tmeanCDMSum, mean_SHAPOctMay_tminMean, mean_SHAPOctMay_tmaxMean, mean_SHAPDecFeb_prcpSumCDMSum, mean_SHAPDecFeb_tmeanCDMSum, mean_SHAPelevation, mean_SHAPslope, mean_SHAPaspect, mean_SHAPlandcover_triclass), # ...for dose column
    direction = "row",
    fn = scales::col_numeric( # <- bc it's numeric
      palette = brewer.pal(11, "PRGn"), # A color scheme (gradient)
      domain = c(), # Column scale endpoints
      reverse = FALSE
    )
  ) |>
  cols_label(class_cluster = "Class", median_swe_max = "Median Max SWE (mm)", mean_SHAPOctMay_prcpSumCDMSum = "Oct-May prcp CDM Sum", mean_SHAPOctMay_tmeanCDMSum = "Oct-May tmean CDM Sum", mean_SHAPOctMay_tminMean = "Oct-May tmin Mean", mean_SHAPOctMay_tmaxMean = "Oct-May tmax Mean", mean_SHAPDecFeb_prcpSumCDMSum = "Dec-Feb prcp CDM Sum", mean_SHAPDecFeb_tmeanCDMSum = "Dec-Feb tmean CDM Sum", mean_SHAPelevation = "Elevation", mean_SHAPslope = "Slope", mean_SHAPaspect = "Aspect", mean_SHAPlandcover_triclass = "Landcover") |> # Update labels
  cols_move_to_end(columns = c("mean_SHAPDecFeb_prcpSumCDMSum", "mean_SHAPDecFeb_tmeanCDMSum", "mean_SHAPelevation", "mean_SHAPslope", "mean_SHAPaspect", "mean_SHAPlandcover_triclass")) |>
  cols_hide(c(mean_AbsSHAP_OctMay_prcpSumCDMSum, mean_AbsSHAP_OctMay_tmeanCDMSum, mean_AbsSHAP_OctMay_tminMean, mean_AbsSHAP_OctMay_tmaxMean, mean_AbsSHAP_DecFeb_prcpSumCDMSum, mean_AbsSHAP_DecFeb_tmeanCDMSum, mean_AbsSHAP_elevation, mean_AbsSHAP_slope, mean_AbsSHAP_aspect, mean_AbsSHAP_landcover_triclass, mean_swe_max, median_SHAPOctMay_prcpSumCDMSum, median_SHAPOctMay_tmeanCDMSum, median_SHAPDecFeb_prcpSumCDMSum, median_SHAPDecFeb_tmeanCDMSum, median_SHAPelevation, median_SHAPslope, median_SHAPaspect, median_SHAPOctMay_tminMean, median_SHAPOctMay_tmaxMean, median_SHAPlandcover_triclass, median_AbsSHAP_OctMay_prcpSumCDMSum, median_AbsSHAP_OctMay_tmeanCDMSum, median_AbsSHAP_DecFeb_prcpSumCDMSum, median_AbsSHAP_DecFeb_tmeanCDMSum, median_AbsSHAP_elevation, median_AbsSHAP_slope, median_AbsSHAP_aspect, median_AbsSHAP_OctMay_tminMean, median_AbsSHAP_OctMay_tmaxMean, median_AbsSHAP_landcover_triclass)) |>
  opt_vertical_padding(scale = 3)
# tab_footnote(
#   footnote = ""
# )
AKCluster_SHAPAvg_gt
### saving table
AKCluster_SHAPAvg_gt |>
  gtsave(
    "AKClusterSHAPAvg.png", expand = 5,
    path = here("Outputs", "ChartsGraphsTables", "Kmeans")
  )

### same table but with a title for the presentation
AKCluster_SHAPAvgTitle_gt <- AKCluster_Avg |>
  gt() |> # use 'gt' to make an awesome table...
  gt_theme_guardian() |>
  tab_header(
    title = "SWE Class SHAP Value Averages", # ...with this title
    subtitle = "Alaska")  |>  # and this subtitle
  opt_align_table_header(align = "center") |>
  ##tab_style(style = cell_fill("bisque"),
  ##        locations = cells_body()) |>  # add fill color to table
  fmt_number( # A column (numeric data)
    columns = c(median_swe_max), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # A column (numeric data)
    columns = c(mean_SHAPOctMay_prcpSumCDMSum), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # A column (numeric data)
    columns = c(mean_SHAPOctMay_tmeanCDMSum), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # A column (numeric data)
    columns = c(mean_SHAPOctMay_tminMean), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # Another column (also numeric data)
    columns = c(mean_SHAPOctMay_tmaxMean), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
  ) |> 
  fmt_number( # Another column (also numeric data)
    columns = c(mean_SHAPDecFeb_prcpSumCDMSum), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
  ) |>
  fmt_number( # Another column (also numeric data)
    columns = c(mean_SHAPDecFeb_tmeanCDMSum), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
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
    columns = c(mean_SHAPOctMay_prcpSumCDMSum, mean_SHAPOctMay_tmeanCDMSum, mean_SHAPOctMay_tminMean, mean_SHAPOctMay_tmaxMean, mean_SHAPDecFeb_prcpSumCDMSum, mean_SHAPDecFeb_tmeanCDMSum, mean_SHAPelevation, mean_SHAPslope, mean_SHAPaspect, mean_SHAPlandcover_triclass), # ...for dose column
    direction = "row",
    fn = scales::col_numeric( # <- bc it's numeric
      palette = brewer.pal(11, "PRGn"), # A color scheme (gradient)
      domain = c(), # Column scale endpoints
      reverse = FALSE
    )
  ) |>
  cols_label(class_cluster = "Class", median_swe_max = "Median Max SWE (mm)", mean_SHAPOctMay_prcpSumCDMSum = "Oct-May prcp CDM Sum", mean_SHAPOctMay_tmeanCDMSum = "Oct-May tmean CDM Sum", mean_SHAPOctMay_tminMean = "Oct-May tmin Mean", mean_SHAPOctMay_tmaxMean = "Oct-May tmax Mean", mean_SHAPDecFeb_prcpSumCDMSum = "Dec-Feb prcp CDM Sum", mean_SHAPDecFeb_tmeanCDMSum = "Dec-Feb tmean CDM Sum", mean_SHAPelevation = "Elevation", mean_SHAPslope = "Slope", mean_SHAPaspect = "Aspect", mean_SHAPlandcover_triclass = "Landcover") |> # Update labels
  cols_move_to_end(columns = c("mean_SHAPDecFeb_prcpSumCDMSum", "mean_SHAPDecFeb_tmeanCDMSum", "mean_SHAPelevation", "mean_SHAPslope", "mean_SHAPaspect", "mean_SHAPlandcover_triclass")) |>
  cols_hide(c(mean_AbsSHAP_OctMay_prcpSumCDMSum, mean_AbsSHAP_OctMay_tmeanCDMSum, mean_AbsSHAP_OctMay_tminMean, mean_AbsSHAP_OctMay_tmaxMean, mean_AbsSHAP_DecFeb_prcpSumCDMSum, mean_AbsSHAP_DecFeb_tmeanCDMSum, mean_AbsSHAP_elevation, mean_AbsSHAP_slope, mean_AbsSHAP_aspect, mean_AbsSHAP_landcover_triclass, mean_swe_max, median_SHAPOctMay_prcpSumCDMSum, median_SHAPOctMay_tmeanCDMSum, median_SHAPDecFeb_prcpSumCDMSum, median_SHAPDecFeb_tmeanCDMSum, median_SHAPelevation, median_SHAPslope, median_SHAPaspect, median_SHAPOctMay_tminMean, median_SHAPOctMay_tmaxMean, median_SHAPlandcover_triclass, median_AbsSHAP_OctMay_prcpSumCDMSum, median_AbsSHAP_OctMay_tmeanCDMSum, median_AbsSHAP_DecFeb_prcpSumCDMSum, median_AbsSHAP_DecFeb_tmeanCDMSum, median_AbsSHAP_elevation, median_AbsSHAP_slope, median_AbsSHAP_aspect, median_AbsSHAP_OctMay_tminMean, median_AbsSHAP_OctMay_tmaxMean, median_AbsSHAP_landcover_triclass)) |>
  opt_vertical_padding(scale = 3)
# tab_footnote(
#   footnote = ""
# )
AKCluster_SHAPAvgTitle_gt
### saving table
AKCluster_SHAPAvgTitle_gt |>
  gtsave(
    "AKClusterSHAPAvgTitle.png", expand = 5,
    path = here("Outputs", "ChartsGraphsTables", "Kmeans")
  )


##### grouping predictors by cluster so I can make a table for them too #####
AK_Predictors_All <- fread(here("Data", "SHAP", "Fastshap", "SHAPs", "Alaska", "AllYearsCombined", "AKAllPredictors_KClasses.csv"))

### define function to calculate mode
find_mode <- function(x) {
  u <- unique(x)
  tab <- tabulate(match(x, u))
  u[tab == max(tab)]
}



AKCluster_PredAvg <- AK_Predictors_All |>
  group_by(class_cluster) |>
  summarise(mean_swe_max = mean(peak_swe),
            median_swe_max = median(peak_swe),
            mean_OctMay_prcpSumCDMSum = mean(OctMay_prcpSumCDMSum), 
            mean_OctMay_tmeanCDMSum = mean(OctMay_tmeanCDMSum), 
            mean_DecFeb_prcpSumCDMSum = mean(DecFeb_prcpSumCDMSum), 
            mean_DecFeb_tmeanCDMSum = mean(DecFeb_tmeanCDMSum), 
            mean_elevation = mean(elevation), 
            mean_slope = mean(slope), 
            mean_aspect = mean(aspect), 
            mean_OctMay_tminMean = mean(OctMay_tminMean), 
            mean_OctMay_tmaxMean = mean(OctMay_tmaxMean),
            median_OctMay_prcpSumCDMSum = median(OctMay_prcpSumCDMSum), 
            median_OctMay_tmeanCDMSum = median(OctMay_tmeanCDMSum), 
            median_DecFeb_prcpSumCDMSum = median(DecFeb_prcpSumCDMSum), 
            median_DecFeb_tmeanCDMSum = median(DecFeb_tmeanCDMSum), 
            median_elevation = median(elevation), 
            median_slope = median(slope), 
            median_aspect = median(aspect), 
            median_OctMay_tminMean = median(OctMay_tminMean), 
            median_OctMay_tmaxMean = median(OctMay_tmaxMean),
            mean_landcover_triclass = find_mode(landcover_triclass))

### writing group averages as csv
fwrite(AKCluster_PredAvg, here("Data", "SHAP", "Classes", "Alaska", "KMeans", "CSVs", "Summary", "AKClassPredAvgs.csv"))


### making table of cluster averages
AKCluster_PredAvg_gt <- AKCluster_PredAvg |>
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
    columns = c(mean_OctMay_prcpSumCDMSum), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # A column (numeric data)
    columns = c(mean_OctMay_tmeanCDMSum), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # A column (numeric data)
    columns = c(mean_OctMay_tminMean), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # Another column (also numeric data)
    columns = c(mean_OctMay_tmaxMean), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
  ) |> 
  fmt_number( # Another column (also numeric data)
    columns = c(mean_DecFeb_prcpSumCDMSum), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
  ) |>
  fmt_number( # Another column (also numeric data)
    columns = c(mean_DecFeb_tmeanCDMSum), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
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
    columns = c(median_swe_max),
    fn = scales::col_numeric( # <- bc it's numeric
      palette = brewer.pal(11, "RdBu"), # A color scheme (gradient)
      domain = c(), # Column scale endpoints
      reverse = FALSE
    )
  ) |>
  data_color( # Update cell colors, testing different color palettes
    columns = c(mean_OctMay_prcpSumCDMSum, mean_OctMay_tmeanCDMSum, mean_OctMay_tminMean, mean_OctMay_tmaxMean, mean_DecFeb_prcpSumCDMSum, mean_DecFeb_tmeanCDMSum, mean_elevation, mean_slope, mean_aspect, mean_landcover_triclass), # ...for dose column
    # direction = "row",
    fn = scales::col_numeric( # <- bc it's numeric
      palette = brewer.pal(11, "PRGn"), # A color scheme (gradient)
      domain = c(), # Column scale endpoints
      reverse = FALSE
    )
  ) |>
  cols_label(class_cluster = "Class", median_swe_max = "Median Max SWE (mm)", mean_OctMay_prcpSumCDMSum = "Oct-May prcp CDM Sum", mean_OctMay_tmeanCDMSum = "Oct-May tmean CDM Sum", mean_OctMay_tminMean = "Oct-May tmin Mean", mean_OctMay_tmaxMean = "Oct-May tmax Mean", mean_DecFeb_prcpSumCDMSum = "Dec-Feb prcp CDM Sum", mean_DecFeb_tmeanCDMSum = "Dec-Feb tmean CDM Sum", mean_elevation = "Elevation", mean_slope = "Slope", mean_aspect = "Aspect", mean_landcover_triclass = "Landcover") |> # Update labels
  cols_move_to_end(columns = c("mean_DecFeb_prcpSumCDMSum", "mean_DecFeb_tmeanCDMSum", "mean_elevation", "mean_slope", "mean_aspect", "mean_landcover_triclass")) |>
  cols_hide(c(median_OctMay_prcpSumCDMSum, median_OctMay_tmeanCDMSum, median_DecFeb_prcpSumCDMSum, median_DecFeb_tmeanCDMSum, median_elevation, median_slope, median_aspect, median_OctMay_tminMean, median_OctMay_tmaxMean, mean_swe_max)) |>
  opt_vertical_padding(scale = 3)
# tab_footnote(
#   footnote = ""
# )
AKCluster_PredAvg_gt
### saving table
AKCluster_PredAvg_gt |>
  gtsave(
    "AKClusterPredAvg.png", expand = 5,
    path = here("Outputs", "ChartsGraphsTables", "Kmeans")
  )

### same table but with a title for the presentation
AKCluster_PredAvgTitle_gt <- AKCluster_PredAvg |>
  gt() |> # use 'gt' to make an awesome table...
  gt_theme_guardian() |>
  tab_header(
    title = "SWE Class Variable Averages", # ...with this title
    subtitle = "Alaska")  |>  # and this subtitle
  opt_align_table_header(align = "center") |>
  ##tab_style(style = cell_fill("bisque"),
  ##        locations = cells_body()) |>  # add fill color to table
  fmt_number( # A column (numeric data)
    columns = c(median_swe_max), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # A column (numeric data)
    columns = c(mean_OctMay_prcpSumCDMSum), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # A column (numeric data)
    columns = c(mean_OctMay_tmeanCDMSum), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # A column (numeric data)
    columns = c(mean_OctMay_tminMean), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  fmt_number( # Another column (also numeric data)
    columns = c(mean_OctMay_tmaxMean), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
  ) |> 
  fmt_number( # Another column (also numeric data)
    columns = c(mean_DecFeb_prcpSumCDMSum), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
  ) |>
  fmt_number( # Another column (also numeric data)
    columns = c(mean_DecFeb_tmeanCDMSum), # What column variable? FinalVoATop25$VoA_Ranking
    decimals = 2 # I want this column to have zero decimal places
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
    columns = c(mean_OctMay_prcpSumCDMSum, mean_OctMay_tmeanCDMSum, mean_OctMay_tminMean, mean_OctMay_tmaxMean, mean_DecFeb_prcpSumCDMSum, mean_DecFeb_tmeanCDMSum, mean_elevation, mean_slope, mean_aspect, mean_landcover_triclass), # ...for dose column
    # direction = "row",
    fn = scales::col_numeric( # <- bc it's numeric
      palette = brewer.pal(11, "PRGn"), # A color scheme (gradient)
      domain = c(), # Column scale endpoints
      reverse = FALSE
    )
  ) |>
  cols_label(class_cluster = "Class", median_swe_max = "Median Max SWE (mm)", mean_OctMay_prcpSumCDMSum = "Oct-May prcp CDM Sum", mean_OctMay_tmeanCDMSum = "Oct-May tmean CDM Sum", mean_OctMay_tminMean = "Oct-May tmin Mean", mean_OctMay_tmaxMean = "Oct-May tmax Mean", mean_DecFeb_prcpSumCDMSum = "Dec-Feb prcp CDM Sum", mean_DecFeb_tmeanCDMSum = "Dec-Feb tmean CDM Sum", mean_elevation = "Elevation", mean_slope = "Slope", mean_aspect = "Aspect", mean_landcover_triclass = "Landcover") |> # Update labels
  cols_move_to_end(columns = c("mean_DecFeb_prcpSumCDMSum", "mean_DecFeb_tmeanCDMSum", "mean_elevation", "mean_slope", "mean_aspect", "mean_landcover_triclass")) |>
  cols_hide(c(median_OctMay_prcpSumCDMSum, median_OctMay_tmeanCDMSum, median_DecFeb_prcpSumCDMSum, median_DecFeb_tmeanCDMSum, median_elevation, median_slope, median_aspect, median_OctMay_tminMean, median_OctMay_tmaxMean, mean_swe_max)) |>
  opt_vertical_padding(scale = 3) |>
  tab_footnote(
    footnote = "Landcover values are modes, all other values are means except max SWE"
    )
AKCluster_PredAvgTitle_gt
### saving table
AKCluster_PredAvgTitle_gt |>
  gtsave(
    "AKClusterPredAvgTitle.png", expand = 5,
    path = here("Outputs", "ChartsGraphsTables", "Kmeans")
  )


##### Making stacked bar chart showing each covariate's percentage of total importance for each class #####
### reading group averages
AKCluster_Avg <- read_csv(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "CSVs", "Summary", "AKClassSHAPAvgs.csv")) |>
  mutate(total_absshap = mean_AbsSHAP_OctMay_prcpSumCDMSum + mean_AbsSHAP_OctMay_tmeanCDMSum + mean_AbsSHAP_DecFeb_prcpSumCDMSum + mean_AbsSHAP_DecFeb_tmeanCDMSum + mean_AbsSHAP_elevation + mean_AbsSHAP_slope + mean_AbsSHAP_aspect + mean_AbsSHAP_OctMay_tminMean + mean_AbsSHAP_OctMay_tmaxMean + mean_AbsSHAP_landcover_triclass,
         pct_AbsSHAP_OctMay_prcpSumCDMSum = mean_AbsSHAP_OctMay_prcpSumCDMSum / total_absshap * 100, 
         pct_AbsSHAP_OctMay_tmeanCDMSum = mean_AbsSHAP_OctMay_tmeanCDMSum / total_absshap * 100,
         pct_AbsSHAP_DecFeb_prcpSumCDMSum = mean_AbsSHAP_DecFeb_prcpSumCDMSum / total_absshap * 100,
         pct_AbsSHAP_DecFeb_tmeanCDMSum = mean_AbsSHAP_DecFeb_tmeanCDMSum / total_absshap * 100,
         pct_AbsSHAP_elevation = mean_AbsSHAP_elevation / total_absshap * 100,
         pct_AbsSHAP_slope = mean_AbsSHAP_slope / total_absshap * 100,
         pct_AbsSHAP_aspect = mean_AbsSHAP_aspect / total_absshap * 100,
         pct_AbsSHAP_OctMay_tminMean = mean_AbsSHAP_OctMay_tminMean / total_absshap * 100,
         pct_AbsSHAP_OctMay_tmaxMean = mean_AbsSHAP_OctMay_tmaxMean / total_absshap * 100,
         pct_AbsSHAP_landcover_triclass = mean_AbsSHAP_landcover_triclass / total_absshap * 100) |>
  select(class_cluster, total_absshap, pct_AbsSHAP_OctMay_prcpSumCDMSum, pct_AbsSHAP_OctMay_tmeanCDMSum, pct_AbsSHAP_DecFeb_prcpSumCDMSum, pct_AbsSHAP_DecFeb_tmeanCDMSum, pct_AbsSHAP_elevation, pct_AbsSHAP_slope, pct_AbsSHAP_aspect, pct_AbsSHAP_OctMay_tminMean, pct_AbsSHAP_OctMay_tmaxMean, pct_AbsSHAP_landcover_triclass)

### ggplot time
# c(pct_AbsSHAP_OctMay_prcpSumCDMSum, pct_AbsSHAP_OctMay_tmeanCDMSum, pct_AbsSHAP_DecFeb_prcpSumCDMSum, pct_AbsSHAP_DecFeb_tmeanCDMSum, pct_AbsSHAP_elevation, pct_AbsSHAP_slope, pct_AbsSHAP_aspect, pct_AbsSHAP_OctMay_tminMean, pct_AbsSHAP_OctMay_tmaxMean, pct_AbsSHAP_landcover_triclass)))
AKCluster_Avg_Pivot = AKCluster_Avg |>
  pivot_longer(
    cols = starts_with("pct_"),
  )

### setting name to be as a factor so I can set what order the bars are in in the stacked bar chart below
AKCluster_Avg_Pivot$name = factor(AKCluster_Avg_Pivot$name, levels = c("pct_AbsSHAP_landcover_triclass", "pct_AbsSHAP_elevation", "pct_AbsSHAP_aspect", "pct_AbsSHAP_slope", "pct_AbsSHAP_OctMay_tminMean", "pct_AbsSHAP_OctMay_tmaxMean", "pct_AbsSHAP_DecFeb_tmeanCDMSum", "pct_AbsSHAP_OctMay_tmeanCDMSum", "pct_AbsSHAP_DecFeb_prcpSumCDMSum", "pct_AbsSHAP_OctMay_prcpSumCDMSum"))

PctAbsSHAP_plot <- ggplot(data = AKCluster_Avg_Pivot, mapping = aes(x = class_cluster, y = value, fill = name)) +
  theme_bw() +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_brewer("", palette = "PRGn", labels = c("Landcover", "Elevation", "Aspect", "Slope", "Oct-May Tmin Mean", "Oct-May Tmax Mean", "Dec-Feb Tmean CDM Sum", "Oct-May Tmean CDM Sum", "Dec-Feb prcp CDM Sum", "Oct-May prcp CDM Sum")) +
  ggtitle(label = "", subtitle = "Alaska") +
  xlab("Class") +
  ylab("% of Absolute SHAP Value Sum") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
  
PctAbsSHAP_plot
ggsave(plot = PctAbsSHAP_plot, here("Outputs", "ChartsGraphsTables", "SHAPPcts", "AKPctAbsSHAP.png"))

### now with a title
PctAbsSHAPTitle_plot <- ggplot(data = AKCluster_Avg_Pivot, mapping = aes(x = class_cluster, y = value, fill = name)) +
  theme_bw() +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_brewer("", palette = "PRGn", labels = c("Landcover", "Elevation", "Aspect", "Slope", "Oct-May Tmin Mean", "Oct-May Tmax Mean", "Dec-Feb Tmean CDM Sum", "Oct-May Tmean CDM Sum", "Dec-Feb prcp CDM Sum", "Oct-May prcp CDM Sum")) +
  ggtitle(label = "Predictor Proportion of Total Absolute SHAP Value by Class", subtitle = "Alaska") +
  xlab("Class") +
  ylab("% of Absolute SHAP Value Sum") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

PctAbsSHAPTitle_plot
ggsave(plot = PctAbsSHAPTitle_plot, here("Outputs", "ChartsGraphsTables", "SHAPPcts", "AKPctAbsSHAP_Title.png"))
