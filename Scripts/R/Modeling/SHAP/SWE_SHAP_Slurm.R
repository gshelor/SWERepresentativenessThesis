##### SHAP #####
### as opposed to running this script locally, intending to run it on UNR's HPC Pronghorn using Slurm
totalstarttime <- Sys.time()
##### Loading packages #####
# install.packages("pacman")
library(tidyverse)
# library(future)
install.packages("treeshap", lib = "/data/gpfs/home/gshelor")
library(treeshap, pos = "/data/gpfs/home/gshelor/treeshap")
library(sf)
library(sfext)
# p_load(caret, parallel, treeshap, future)
# options(mc.cores = parallel::detectCores())
# future::plan("multisession", workers = 16)

##### reading in data #####
RFModel <- readRDS("OctAprAspect_Nosradtmintmax_RFModel.rds")
### SNOTEL Climatologies used to fit the model
Snotel_Clim_sf <- read_sf(here("Data", "SNOTEL", "Combined", "GIS", "SnotelCombined_ClimatologyCovars.gpkg"))
### converting to dfs for model fitting process
Snotel_Clim_df <- sf_to_df(Snotel_Clim_sf) |>
  drop_na()
### making sure landcover is type factor
Snotel_Clim_df$landcover <- as.integer(Snotel_Clim_df$landcover)
Snotel_Clim_df$landcover_triclass <- as.integer(Snotel_Clim_df$landcover_triclass)
### OctApr aggregates only but with aspect, and srad, tmin, tmax removed
train_OctAprAspect_Nosradtmintmax_x <- Snotel_Clim_df[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "landcover_triclass")]


##### SHAP stuff #####
### in order to run treeshap, the model must be "unified"
UnifiedRF <- unify(RFModel$finalModel, data = train_OctAprAspect_Nosradtmintmax_x)

### years of study period
# years <- 1993:2020
# for (i in 1:length(years)){
#   temp_starttime <- Sys.time()
#   print(temp_starttime)
#   print(years[i])
#   print("CONUS")
#   ### reading in csv versions of data frames so treeshap can understand them
#   ### SWE_max predictions
#   print("reading in prediction csv")
#   CONUSPredictions <- read_csv(here("Data", "SHAP", "PredictionCSVs", "CONUS", paste0("CONUS_SWEmaxPredictions", years[i], ".csv")))
#   ### removing xy
#   CONUSPredictions_NoXY <- CONUSPredictions[,3]
#   colnames(CONUSPredictions_NoXY) <- "peak_SWE"
#   
#   print("reading in predictor csv")
#   CONUSPredictors <- read_csv(here("Data", "SHAP", "PredictorCSVs", "CONUS", paste0("CONUS_SWEmaxPredictors", years[i], ".csv"))) |>
#     filter(paste0(x,y) %in% paste0(CONUSPredictions$x, CONUSPredictions$y))
#   ### removing xy data for treeshap
#   CONUSPredictors_NoXY <- CONUSPredictors |>
#     select(OctApr_prcpSumCDMSum, OctApr_tmeanCDMSum, elevation, slope, aspect, landcover_triclass)
#   ### storing xy separately
#   CONUSPredictors_XY <- CONUSPredictors |>
#     select(x,y)
#   
#   ### Running Treeshap
#   print("running treeshap")
#   set.seed(802)
#   CONUSTreeshap <- treeshap(unified_model = UnifiedRF, x = CONUSPredictors[1:500,])
#   write_rds(CONUSTreeshap, here("Data", "SHAP", "Treeshap", "CONUS", paste0("CONUSSWEmax_SHAPs", years[i], ".rds")))
#   
#   
#   ### Alaska
#   print("Alaska")
#   ### reading in csv versions of data frames so treeshap can understand them
#   ### SWE_max predictions
#   print("reading in prediction csv")
#   AKPredictions <- read_csv(here("Data", "SHAP", "PredictionCSVs", "Alaska", paste0("AK_SWEmaxPredictions", years[i], ".csv")))
#   ### removing xy
#   AKPredictions_NoXY <- AKPredictions[,3]
#   colnames(AKPredictions_NoXY) <- "peak_SWE"
#   
#   print("reading in predictor csv")
#   AKPredictors <- read_csv(here("Data", "SHAP", "PredictorCSVs", "Alaska", paste0("AK_SWEmaxPredictors", years[i], ".csv"))) |>
#     filter(paste0(x,y) %in% paste0(AKPredictions$x, AKPredictions$y))
#   ### removing xy data for treeshap
#   AKPredictors_NoXY <- AKPredictors |>
#     select(OctApr_prcpSumCDMSum, OctApr_tmeanCDMSum, elevation, slope, aspect, landcover_triclass)
#   
#   ### Running Treeshap
#   print("running treeshap")
#   AKTreeshap <- treeshap(unified_model = UnifiedRF, x = AKPredictors)
#   write_rds(AKTreeshap, here("Data", "SHAP", "Treeshap", "Alaska", paste0("AKSWEmax_SHAPs", years[i], ".rds")))
#   
#   temp_endtime <- Sys.time()
#   print(temp_endtime)
#   temp_endtime - temp_starttime
# }

### Alaska
print("Alaska")
### reading in csv versions of data frames so treeshap can understand them
### SWE_max predictions
print("reading in prediction csv")
AKPredictions <- read_csv("AK_SWEmaxPredictions1993.csv")
### removing xy
# AKPredictions_NoXY <- AKPredictions[,3]
# colnames(AKPredictions_NoXY) <- "peak_SWE"

print("reading in predictor csv")
AKPredictors <- read_csv("AK_SWEmaxPredictors1993.csv") |>
  filter(paste0(x,y) %in% paste0(AKPredictions$x, AKPredictions$y))
### removing xy data for treeshap
# AKPredictors_NoXY <- AKPredictors |>
#   select(OctApr_prcpSumCDMSum, OctApr_tmeanCDMSum, elevation, slope, aspect, landcover_triclass)

### Running Treeshap
print("running treeshap")
AKTreeshap <- treeshap(unified_model = UnifiedRF, x = AKPredictors)
write_rds(AKTreeshap, "AKSWEmax_SHAPs1993.rds")

totalendtime <- Sys.time()
totaltime <- totalendtime - totalstarttime
totaltime
