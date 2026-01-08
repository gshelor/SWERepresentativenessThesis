##### Comparing My SWE-based classes to Sturm and Liston #####

##### Loading packages, reading in data #####
library(pacman)
p_load(tidyverse, here, sf, terra, tidyterra, sabre, data.table, bespatial, parallel, mcprogress)

### reading in AOIs
AK_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "Alaska", "AK_AOI.gpkg"))
CONUS_AOI <- read_sf(here("Data", "L3_Ecoregions_USB", "CONUS", "CONUSEcoregionsNoStates.gpkg"))


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


### Reading in Modal class raster #####
CONUSModalClass_rast <- rast(here("Data", "SHAP", "Classes", "CONUS", "KMeans", "Rasters", "Summary", "CONUSModalClass.tif"))
plot(as.factor(CONUSModalClass_rast), main = "Modal Class")
AKModalClass_rast <- rast(here("Data", "SHAP", "Classes", "Alaska", "KMeans", "Rasters", "Summary", "AKModalClass.tif"))
plot(AKModalClass_rast, main = "AK Modal Class")

### reading in Sturm classes
SturmListonCONUS_rast <- rast(here("Data", "SturmListon2021Classes", "Resampled", "CONUSSturmClasses.tif"))
plot(SturmListonCONUS_rast, main = "Sturm WUS classes")
SturmListonAK_rast <- rast(here("Data", "SturmListon2021Classes", "Resampled", "AKSturmClasses.tif"))
plot(SturmListonAK_rast, main = "Sturm AK classes")

##### Mapcurves Map Comparison #####
AKGOF_modal <- sabre::mapcurves_calc(as.factor(SturmListonAK_rast), as.factor(AKModalClass_rast))
AKGOF_modal$gof

CONUSGOF_modal <- sabre::mapcurves_calc(as.factor(SturmListonCONUS_rast), as.factor(CONUSModalClass_rast))
CONUSGOF_modal$gof
# plot(CONUSGOF_modal$map2)



##### V-measure Map Comparison #####
CONUSVMes_modal <- sabre::vmeasure_calc(as.factor(SturmListonCONUS_rast), as.factor(CONUSModalClass_rast))
CONUSVMes_modal

AKVmes_modal <- vmeasure_calc(as.factor(SturmListonAK_rast), as.factor(AKModalClass_rast))
AKVmes_modal

### making lists of rasters for each region
CONUSAllYearsClasses_list <- list(CONUS1993Classes_rast, CONUS1994Classes_rast, CONUS1995Classes_rast, CONUS1996Classes_rast, CONUS1997Classes_rast, CONUS1998Classes_rast, CONUS1999Classes_rast, CONUS2000Classes_rast, CONUS2001Classes_rast, CONUS2002Classes_rast, CONUS2003Classes_rast, CONUS2004Classes_rast, CONUS2005Classes_rast, CONUS2006Classes_rast, CONUS2007Classes_rast, CONUS2008Classes_rast, CONUS2009Classes_rast, CONUS2010Classes_rast, CONUS2011Classes_rast, CONUS2012Classes_rast, CONUS2013Classes_rast, CONUS2014Classes_rast, CONUS2015Classes_rast, CONUS2016Classes_rast, CONUS2017Classes_rast, CONUS2018Classes_rast, CONUS2019Classes_rast, CONUS2020Classes_rast)
AKAllYearsClasses_list <- list(AK1993Classes_rast, AK1994Classes_rast, AK1995Classes_rast, AK1996Classes_rast, AK1997Classes_rast, AK1998Classes_rast, AK1999Classes_rast, AK2000Classes_rast, AK2001Classes_rast, AK2002Classes_rast, AK2003Classes_rast, AK2004Classes_rast, AK2005Classes_rast, AK2006Classes_rast, AK2007Classes_rast, AK2008Classes_rast, AK2009Classes_rast, AK2010Classes_rast, AK2011Classes_rast, AK2012Classes_rast, AK2013Classes_rast, AK2014Classes_rast, AK2015Classes_rast, AK2016Classes_rast, AK2017Classes_rast, AK2018Classes_rast, AK2019Classes_rast, AK2020Classes_rast)

### functions to calculate v-measure for each year
CONUSVmes_func <- function(classraster){
  CONUSVmes <- vmeasure_calc(as.factor(SturmListonCONUS_rast), as.factor(classraster))
}
AKVmes_func <- function(AKclassraster){
  AKVmes <- vmeasure_calc(as.factor(SturmListonAK_rast), as.factor(AKclassraster))
}

### calculating Vmeasure for CONUS classes
CONUS_Vmeasures <- pmclapply(X = CONUSAllYearsClasses_list, FUN = CONUSVmes_func, mc.cores = length(CONUSAllYearsClasses_list), mc.set.seed = FALSE, spinner = TRUE)
### calculating v-measure for Alaska classes
AK_Vmeasures <- pmclapply(X = AKAllYearsClasses_list, FUN = AKVmes_func, mc.cores = length(AKAllYearsClasses_list), mc.set.seed = FALSE, spinner = TRUE)



##### Analyzing V-measures #####
### WUS/CONUS
CONUS_VmeasureScores <- c()
CONUS_HomogeneityScores <- c()
CONUS_CompletenessScores <- c()
for (i in 1:28){
  temp_vmes <- CONUS_Vmeasures[[i]]
  CONUS_VmeasureScores <- c(CONUS_VmeasureScores, temp_vmes$v_measure)
  CONUS_HomogeneityScores <- c(CONUS_HomogeneityScores, temp_vmes$homogeneity)
  CONUS_CompletenessScores <- c(CONUS_CompletenessScores, temp_vmes$completeness)
}

### putting Vmeasure metrics into df
CONUS_Vmeasure_df <- data.frame(years = years,
                                v_measures = CONUS_VmeasureScores,
                                homogeneity = CONUS_HomogeneityScores,
                                completeness = CONUS_CompletenessScores)
### saving df as csv
write_csv(CONUS_Vmeasure_df, here("Outputs", "CSVs", "Vmeasure", "CONUSVmeasureScores.csv"))
### plotting scores by year
plot(1993:2020, CONUS_VmeasureScores, main = "CONUS V-Measure Values", ylim = c(min(CONUS_CompletenessScores),max(CONUS_HomogeneityScores)))
points(1993:2020, CONUS_HomogeneityScores, col = "blue")
points(1993:2020, CONUS_CompletenessScores, col = "red")

### need to pivot for the plot
CONUS_Vmeasure_pivot_df <- CONUS_Vmeasure_df |>
  pivot_longer(
    cols = c(v_measures, homogeneity, completeness),
    names_to = "metric",
    values_to = "scores"
  )

CONUSVmeasure_plot <- ggplot(data = CONUS_Vmeasure_pivot_df, mapping = aes(x = years, y = scores)) +
  theme_bw() +
  geom_point(aes(color = metric)) +
  geom_line(aes(color = metric)) +
  # geom_point(mapping = aes(y = homogeneity, col = homogeneity)) +
  # geom_line(mapping = aes(y = homogeneity, col = homogeneity)) +
  # geom_point(mapping = aes(y = completeness, col = completeness)) +
  # geom_line(mapping = aes(y = completeness, col = completeness)) +
  scale_color_manual("", labels = c("Completeness", "Homogeneity", "V-measure"), values = c("red", "blue", "black")) +
  # ggtitle(label = "Total Sum of Squared Error for Different Numbers of K-Means Clusters", subtitle = "WUS Study Area") +
  ggtitle(label = "WUS Study Area") +
  xlab("Year") +
  ylab("V-Measure Scores") +
  theme(plot.title = element_text(hjust = 0.5, size = 12), plot.subtitle = element_text(hjust = 0.5))
CONUSVmeasure_plot
ggsave(plot = CONUSVmeasure_plot, here("Outputs", "Plots", "Vmeasure", "CONUSVMeasureScores.png"))


### using linear model to evaluate R^2 between years and v-measure scores
years <- 1993:2020
set.seed(802)
CONUSVMesModel <- lm(CONUS_VmeasureScores ~ years)
summary(CONUSVMesModel)

### Alaska
AK_VmeasureScores <- c()
AK_HomogeneityScores <- c()
AK_CompletenessScores <-c()
for (i in 1:28){
  temp_vmes <- AK_Vmeasures[[i]]
  AK_VmeasureScores <- c(AK_VmeasureScores, temp_vmes$v_measure)
  AK_HomogeneityScores <- c(AK_HomogeneityScores, temp_vmes$homogeneity)
  AK_CompletenessScores <- c(AK_CompletenessScores, temp_vmes$completeness)
}

### putting Vmeasure metrics into df
AK_Vmeasure_df <- data.frame(years = years,
                             v_measures = AK_VmeasureScores,
                             homogeneity = AK_HomogeneityScores,
                             completeness = AK_CompletenessScores)
### saving df as csv
write_csv(AK_Vmeasure_df, here("Outputs", "CSVs", "Vmeasure", "AKVmeasureScores.csv"))

plot(1993:2020, AK_VmeasureScores, main = "AK V-Measure Values", ylim = c(min(c(AK_CompletenessScores, AK_HomogeneityScores)),max(c(AK_HomogeneityScores, AK_CompletenessScores))))
points(1993:2020, AK_HomogeneityScores, col = "blue")
points(1993:2020, AK_CompletenessScores, col = "red")
### need to pivot for the plot
AK_Vmeasure_pivot_df <- AK_Vmeasure_df |>
  pivot_longer(
    cols = c(v_measures, homogeneity, completeness),
    names_to = "metric",
    values_to = "scores"
  )

AKVmeasure_plot <- ggplot(data = AK_Vmeasure_pivot_df, mapping = aes(x = years, y = scores)) +
  theme_bw() +
  geom_point(aes(color = metric)) +
  geom_line(aes(color = metric)) +
  # geom_point(mapping = aes(y = homogeneity, col = homogeneity)) +
  # geom_line(mapping = aes(y = homogeneity, col = homogeneity)) +
  # geom_point(mapping = aes(y = completeness, col = completeness)) +
  # geom_line(mapping = aes(y = completeness, col = completeness)) +
  scale_color_manual("", labels = c("Completeness", "Homogeneity", "V-measure"), values = c("red", "blue", "black")) +
  # ggtitle(label = "Total Sum of Squared Error for Different Numbers of K-Means Clusters", subtitle = "WUS Study Area") +
  ggtitle(label = "Alaska Study Area") +
  xlab("Year") +
  ylab("V-Measure Scores") +
  theme(plot.title = element_text(hjust = 0.5, size = 12), plot.subtitle = element_text(hjust = 0.5))
AKVmeasure_plot
ggsave(plot = AKVmeasure_plot, here("Outputs", "Plots", "Vmeasure", "AKVMeasureScores.png"))


### using linear model to evaluate R^2 between years and v-measure scores
years <- 1993:2020
set.seed(802)
AKVMesModel <- lm(AK_VmeasureScores ~ years)
summary(AKVMesModel)
