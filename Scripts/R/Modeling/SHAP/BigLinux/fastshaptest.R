##### Loading packages #####
library(pacman)
p_load(tidyverse, randomForest, sf, terra, sfext, caret, parallel, treeshap, data.table, rsample, mcprogress, fastshap, doParallel, shapviz)
# options(mc.cores = parallel::detectCores())
# future::plan("multisession")
setwd("/media/Research/Morafkai/GriffinS")

##### reading in data #####
RFModel <- read_rds("FittedModels/AnnualRFModel10.rds")
### CONUS AOI
# CONUS_AOI <- read_sf("SHAPAnalysis/Data/L3_Ecoregions_USB/CONUS/CONUS_AOI.gpkg")
### Alaska AOI
# AK_AOI <- read_sf("SHAPAnalysis/Data/L3_Ecoregions_USB/Alaska/AK_AOI.gpkg")
### SNOTEL Data used to fit the model
Snotel_sf <- read_sf("Combined/GIS/SnotelCombined_AnnualCovars_AnnualZeroPts.gpkg")


### converting to dfs for model fitting process
Snotel_df <- sf_to_df(Snotel_sf) |>
  drop_na(OctApr_prcpSumCDMSum, OctApr_tmeanCDMSum, elevation, slope, aspect, landcover_triclass)
### making sure landcover is type factor
Snotel_df$landcover <- as.integer(Snotel_df$landcover)
Snotel_df$landcover_triclass <- as.integer(Snotel_df$landcover_triclass)

##### splitting data into training and testing datasets #####
set.seed(802)
Snotel_split <- initial_split(Snotel_df, prop = 0.75)
Snotel_train <- training(Snotel_split)
Snotel_test <- testing(Snotel_split)

### OctApr aggregates only but with aspect, and srad, tmin, tmax removed
train_OctAprAspect_Nosradtmintmax_x <- Snotel_train[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "landcover_triclass")]
test_OctAprAspect_Nosradtmintmax_x <- Snotel_test[, c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "landcover_triclass")]
train_y <- Snotel_train$peak_swe
test_y <- Snotel_test$peak_swe

# UnifiedRF <- unify(RFModel$finalModel, data = train_OctAprAspect_Nosradtmintmax_x)


AKPredictors <- fread(file = paste0("SHAPAnalysis/Data/SHAP/PredictorCSVs/Alaska/", paste0("AK_SWEmaxPredictors", 2020, ".csv")), data.table = FALSE)
AKPredictors_fastshap <- AKPredictors |>
  select(OctApr_prcpSumCDMSum, OctApr_tmeanCDMSum, elevation, slope, aspect, landcover_triclass)
CONUSPredictors <- fread(file = paste0("SHAPAnalysis/Data/SHAP/PredictorCSVs/CONUS/", paste0("CONUS_SWEmaxPredictors", 2020, ".csv")), data.table = FALSE)
CONUSPredictors_fastshap <- CONUSPredictors |>
  select(OctApr_prcpSumCDMSum, OctApr_tmeanCDMSum, elevation, slope, aspect, landcover_triclass)

set.seed(802)
random_row_nums <- sort(sample(1:nrow(AKPredictors_fastshap), size = 250, replace = FALSE))
set.seed(802)
CONUSrandom_row_nums <- sort(sample(1:nrow(CONUSPredictors_fastshap), size = 250, replace = FALSE))

RandomRows_df <- AKPredictors_fastshap[random_row_nums,]
CONUSRandomRows_df <- CONUSPredictors_fastshap[CONUSrandom_row_nums,]

##### Treeshap #####
# set.seed(802)
# treeshapstarttime <- Sys.time()
# AK_Treeshap <- treeshap(unified_model = UnifiedRF, x = as.data.frame(head(AKPredictors, 200)))
# treeshapendtime <- Sys.time()
# treeshap_time <- as.numeric(treeshapendtime - treeshapstarttime) * 60
### treeshap typically takes ~ 7 minutes for 200 rows

##### Fastshap #####
### function to specify prediction method for fastshap explain()
pfun <- function(object, newdata) {  # needs to return a numeric vector
  predict(object, newdata)  
}

fastshapstarttime <- Sys.time()
set.seed(802)
AK_Fastshap <- explain(RFModel$finalModel, X = train_OctAprAspect_Nosradtmintmax_x, nsim = 2, pred_wrapper = pfun, newdata = RandomRows_df, adjust = TRUE, parallel = FALSE)
set.seed(802)
AK_Fastshap2 <- explain(RFModel$finalModel, X = train_OctAprAspect_Nosradtmintmax_x, nsim = 5, pred_wrapper = pfun, newdata = RandomRows_df, adjust = TRUE, parallel = FALSE)
Fastshap_diff <- AK_Fastshap2 - AK_Fastshap





nsims <- c(50, 100, 250, 300, 500, 1000, 2500, 5000, 10000)
sim_secondstep <- c(100, 250, 300, 500, 1000, 2500, 5000, 10000)
SHAPChangeTest <- function(sim_num){
  # registerDoParallel(cores = 2)
  set.seed(802)
  AK_Fastshap <- explain(RFModel$finalModel, X = train_OctAprAspect_Nosradtmintmax_x, nsim = sim_num, pred_wrapper = pfun, newdata = RandomRows_df, adjust = TRUE, parallel = FALSE)
}

SHAPs <- pmclapply(nsims, SHAPChangeTest, progress = TRUE, spinner = TRUE, mc.set.seed = FALSE, mc.cores = length(nsims))

PctChanges <- list()
for (x in 2:length(nsims)){
  Fastshap_diff_pct <- abs((SHAPs[[x]] - SHAPs[[x-1]]) / SHAPs[[x-1]]) * 100
  PctChanges[[x-1]] <- Fastshap_diff_pct
}

for (i in 1:length(PctChanges)){
  temp_matrix <- PctChanges[[i]]
  temp_shaps <- SHAPs[[i + 1]]
  for (s in 1:nrow(temp_matrix)){
    for (x in 1:ncol(temp_matrix)){
      if (is.infinite(temp_matrix[s,x])){
        temp_matrix[s,x] = temp_shaps[s,x]
      } else if (is.nan(temp_matrix[s,x]) | is.na(temp_matrix[s,x])){
        temp_matrix[s,x] = 0
      } else{
        temp_matrix[s,x] = temp_matrix[s,x]
      }
    }
  }
  ### assigning corrected matrix back to PctChanges list
  PctChanges[[i]] <- temp_matrix
}


tmean_pctchange <- c()
prcp_pctchange <- c()
elev_pctchange <- c()
slope_pctchange <- c()
aspect_pctchange <- c()
lc_pctchange <- c()
for (i in 1:length(PctChanges)){
  # print(i)
  temp_matrix <- PctChanges[[i]]
  # print(colnames(temp_matrix)[x])
  # print(mean(temp_matrix[,x]))
  tmean_pctchange <- c(tmean_pctchange, median(temp_matrix[,1]))
  prcp_pctchange <- c(prcp_pctchange, median(temp_matrix[,2]))
  elev_pctchange <- c(elev_pctchange, median(temp_matrix[,3]))
  slope_pctchange <- c(slope_pctchange, median(temp_matrix[,4]))
  aspect_pctchange <- c(aspect_pctchange, median(temp_matrix[,5]))
  lc_pctchange <- c(lc_pctchange, median(temp_matrix[,6]))
}

tmean_pctchange <- c(tmean_pctchange[1], tmean_pctchange)
prcp_pctchange <- c(prcp_pctchange[1], prcp_pctchange)
elev_pctchange <- c(elev_pctchange[1], elev_pctchange)
slope_pctchange <- c(slope_pctchange[1], slope_pctchange)
aspect_pctchange <- c(aspect_pctchange[1], aspect_pctchange)
lc_pctchange <- c(lc_pctchange[1], lc_pctchange)

plot(nsims, tmean_pctchange, type = "l", xlab = "nsims", ylab = "pct change in shapley values", ylim = c(0, max(c(tmean_pctchange, prcp_pctchange, elev_pctchange, slope_pctchange, aspect_pctchange, lc_pctchange))), main = "Alaska")
lines(nsims, prcp_pctchange, col = "blue")
lines(nsims, elev_pctchange, col = "red")
lines(nsims, slope_pctchange, col = "orange")
lines(nsims, aspect_pctchange, col = "green")
lines(nsims, lc_pctchange, col = "purple")



### testing changes in Shapley values for CONUS predictions
CONUSSHAPChangeTest <- function(sim_num){
  # registerDoParallel(cores = 2)
  set.seed(802)
  Fastshap <- explain(RFModel$finalModel, X = train_OctAprAspect_Nosradtmintmax_x, nsim = sim_num, pred_wrapper = pfun, newdata = CONUSRandomRows_df, adjust = TRUE, parallel = FALSE)
}
CONUSSHAPs <- pmclapply(nsims, CONUSSHAPChangeTest, progress = TRUE, spinner = TRUE, mc.set.seed = FALSE, mc.cores = length(nsims) * 2)

CONUSPctChanges <- list()
for (x in 2:length(nsims)){
  Fastshap_diff_pct <- abs((CONUSSHAPs[[x]] - CONUSSHAPs[[x-1]]) / CONUSSHAPs[[x-1]]) * 100
  CONUSPctChanges[[x-1]] <- Fastshap_diff_pct
}

for (i in 1:length(CONUSPctChanges)){
  temp_matrix <- CONUSPctChanges[[i]]
  temp_shaps <- CONUSSHAPs[[i + 1]]
  for (s in 1:nrow(temp_matrix)){
    for (x in 1:ncol(temp_matrix)){
      if (is.infinite(temp_matrix[s,x])){
        temp_matrix[s,x] = temp_shaps[s,x]
      } else if (is.nan(temp_matrix[s,x]) | is.na(temp_matrix[s,x])){
        temp_matrix[s,x] = 0
      } else{
        temp_matrix[s,x] = temp_matrix[s,x]
      }
    }
  }
  ### assigning corrected matrix back to PctChanges list
  CONUSPctChanges[[i]] <- temp_matrix
}



CONUStmean_pctchange <- c()
CONUSprcp_pctchange <- c()
CONUSelev_pctchange <- c()
CONUSslope_pctchange <- c()
CONUSaspect_pctchange <- c()
CONUSlc_pctchange <- c()
for (i in 1:length(CONUSPctChanges)){
  # print(i)
  temp_matrix <- CONUSPctChanges[[i]]
  # print(colnames(temp_matrix)[x])
  # print(mean(temp_matrix[,x]))
  CONUStmean_pctchange <- c(CONUStmean_pctchange, median(temp_matrix[,1]))
  CONUSprcp_pctchange <- c(CONUSprcp_pctchange, median(temp_matrix[,2]))
  CONUSelev_pctchange <- c(CONUSelev_pctchange, median(temp_matrix[,3]))
  CONUSslope_pctchange <- c(CONUSslope_pctchange, median(temp_matrix[,4]))
  CONUSaspect_pctchange <- c(CONUSaspect_pctchange, median(temp_matrix[,5]))
  CONUSlc_pctchange <- c(CONUSlc_pctchange, median(temp_matrix[,6]))
}

CONUStmean_pctchange <- c(CONUStmean_pctchange[1], CONUStmean_pctchange)
CONUSprcp_pctchange <- c(CONUSprcp_pctchange[1], CONUSprcp_pctchange)
CONUSelev_pctchange <- c(CONUSelev_pctchange[1], CONUSelev_pctchange)
CONUSslope_pctchange <- c(CONUSslope_pctchange[1], CONUSslope_pctchange)
CONUSaspect_pctchange <- c(CONUSaspect_pctchange[1], CONUSaspect_pctchange)
CONUSlc_pctchange <- c(CONUSlc_pctchange[1], CONUSlc_pctchange)

plot(nsims, CONUStmean_pctchange, type = "l", xlab = "nsims", ylab = "pct change in shapley values", ylim = c(0, max(c(CONUStmean_pctchange, CONUSprcp_pctchange, CONUSelev_pctchange, CONUSslope_pctchange, CONUSaspect_pctchange, CONUSlc_pctchange))), main = "CONUS")
lines(nsims, CONUSprcp_pctchange, col = "blue")
lines(nsims, CONUSelev_pctchange, col = "red")
lines(nsims, CONUSslope_pctchange, col = "orange")
lines(nsims, CONUSaspect_pctchange, col = "green")
lines(nsims, CONUSlc_pctchange, col = "purple")



set.seed(802)
model <- lm(time_lengths ~ nsims)
summary(model)

est_days_func <- function(s, aoi_rows){
  est_days <- s / 60 / 60 * aoi_rows / 200 / 24
}

ak_est_days <- est_days_func(time_lengths, 1535592)
conus_est_days <- est_days_func(time_lengths, 3311878)
model2 <- lm(conus_est_days ~ nsims)
summary(model2)
model3 <- lm(ak_est_days ~ nsims)
summary(model3)
conus_est_treeshap_days <- est_days_func(treeshap_time, 3311878)

plot(nsims, conus_est_days, col = "blue", main = "number of days to run 1 year vs number of simulations used to estimate Shapley value")
points(nsims, ak_est_days,, col = "red")
abline(h = est_days_func(baseline_length * 60, 1535592), col = "red")
abline(h = est_days_func(baseline_length * 60, 3311878), col = "blue")
abline(h = conus_est_treeshap_days, col = "green")


# set.seed(802)
# fastshapstarttime2 <- Sys.time()
# AK_Fastshap2 <- explain(RFModel$finalModel, X = train_OctAprAspect_Nosradtmintmax_x, nsim = 50, pred_wrapper = pfun, newdata = as.data.frame(head(AKPredictors_fastshap, 200)), adjust = FALSE)
# fastshapendtime2 <- Sys.time()



# print("treeshap time")
# print(treeshapendtime - treeshapstarttime)
# print("fastshap time")
# print(fastshapendtime - fastshapstarttime)
# print("fastshap time with no adjustment")
# print(fastshapendtime2 - fastshapstarttime2)

# write_csv(AK_Treeshap$shaps, "treeshaptest.csv")
# write_csv(as.data.frame(AK_Fastshap), "fastshaptest.csv")
# write_csv(as.data.frame(AK_Fastshap2), "fastshaptest2.csv")


##### Comparing Treeshap with Fastshap ####
Fastshap50sims <- read_csv("SHAPTests/NsimTests/SHAPTest_nsim50.csv")
Fastshap25sims <- read_csv("SHAPTests/NsimTests/SHAPTest_nsim25.csv")
Fastshap100sims <- read_csv("SHAPTests/NsimTests/SHAPTest_nsim100.csv")
Fastshap1000sims <- read_csv("SHAPTests/NsimTests/SHAPTest_nsim1000.csv")
Treeshap_test <- read_csv("treeshaptest.csv")

set.seed(802)
t.test(Fastshap1000sims$landcover_triclass, Treeshap_test$landcover_triclass)
set.seed(802)
t.test(Fastshap100sims$landcover_triclass, Treeshap_test$landcover_triclass)
set.seed(802)
t.test(Fastshap50sims$landcover_triclass, Treeshap_test$landcover_triclass)
set.seed(802)
t.test(Fastshap25sims$landcover_triclass, Treeshap_test$landcover_triclass)

plot(Fastshap1000sims$landcover_triclass, Treeshap_test$landcover_triclass)
plot(Fastshap100sims$landcover_triclass, Treeshap_test$landcover_triclass)
plot(Fastshap50sims$landcover_triclass, Treeshap_test$landcover_triclass)
plot(Fastshap25sims$landcover_triclass, Treeshap_test$landcover_triclass)

plot(Fastshap1000sims$aspect, Treeshap_test$aspect)
plot(Fastshap100sims$aspect, Treeshap_test$aspect)
plot(Fastshap50sims$aspect, Treeshap_test$aspect)
plot(Fastshap25sims$aspect, Treeshap_test$aspect)

plot(Fastshap1000sims$landcover_triclass, Treeshap_test$landcover_triclass)
plot(Fastshap100sims$landcover_triclass, Treeshap_test$landcover_triclass)
plot(Fastshap50sims$landcover_triclass, Treeshap_test$landcover_triclass)
plot(Fastshap25sims$landcover_triclass, Treeshap_test$landcover_triclass)
plot(Fastshap1000sims$landcover_triclass, Treeshap_test$landcover_triclass)
plot(Fastshap100sims$landcover_triclass, Treeshap_test$landcover_triclass)
plot(Fastshap50sims$landcover_triclass, Treeshap_test$landcover_triclass)
plot(Fastshap25sims$landcover_triclass, Treeshap_test$landcover_triclass)
plot(Fastshap1000sims$landcover_triclass, Treeshap_test$landcover_triclass)
plot(Fastshap100sims$landcover_triclass, Treeshap_test$landcover_triclass)
plot(Fastshap50sims$landcover_triclass, Treeshap_test$landcover_triclass)
plot(Fastshap25sims$landcover_triclass, Treeshap_test$landcover_triclass)
plot(Fastshap1000sims$landcover_triclass, Treeshap_test$landcover_triclass)
plot(Fastshap100sims$landcover_triclass, Treeshap_test$landcover_triclass)
plot(Fastshap50sims$landcover_triclass, Treeshap_test$landcover_triclass)
plot(Fastshap25sims$landcover_triclass, Treeshap_test$landcover_triclass)



##### Testing Shapviz for feature dependence plots #####
### saving variable names to be listed in ranking order later
vars <- c("OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect", "landcover_triclass")
### dummy value for year so plot code works
year <- 2020
sv <- shapviz(object = AK_Fastshap, X = head(AKPredictors_fastshap, 200))
plot1 <- sv_dependence(sv, v = vars[1], color_var = vars[1]) +
  ggtitle("Oct-Apr Precipitation CDM Sum Feature Dependence", subtitle =
            paste("Alaska", as.character(year))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  xlab(label = "Oct-Apr Total Precipitation CDM Sum (mm)") +
  ylab(label = "SHAP Value")
plot1

plot2 <- sv_dependence(sv, v = vars[2], color_var = vars[2]) +
  ggtitle("Oct-Apr Precipitation CDM Sum Feature Dependence", subtitle =
            paste("Alaska", as.character(year))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  xlab(label = "Oct-Apr Total Precipitation CDM Sum (mm)") +
  ylab(label = "SHAP Value")
plot2

plot3 <- sv_dependence(sv, v = vars[3], color_var = vars[3]) +
  ggtitle("Oct-Apr Precipitation CDM Sum Feature Dependence", subtitle =
            paste("Alaska", as.character(year))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  xlab(label = "Oct-Apr Total Precipitation CDM Sum (mm)") +
  ylab(label = "SHAP Value")
plot3

plot4 <- sv_dependence(sv, v = vars[4], color_var = vars[4]) +
  ggtitle("Oct-Apr Precipitation CDM Sum Feature Dependence", subtitle =
            paste("Alaska", as.character(year))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  xlab(label = "Oct-Apr Total Precipitation CDM Sum (mm)") +
  ylab(label = "SHAP Value")
plot4

plot5 <- sv_dependence(sv, v = vars[5], color_var = vars[5]) +
  ggtitle("Oct-Apr Precipitation CDM Sum Feature Dependence", subtitle =
            paste("Alaska", as.character(year))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  xlab(label = "Oct-Apr Total Precipitation CDM Sum (mm)") +
  ylab(label = "SHAP Value")
plot5

plot6 <- sv_dependence(sv, v = vars[6], color_var = vars[6], jitter_width = 0) +
  ggtitle("Oct-Apr Precipitation CDM Sum Feature Dependence", subtitle =
            paste("Alaska", as.character(year))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  xlab(label = "Oct-Apr Total Precipitation CDM Sum (mm)") +
  ylab(label = "SHAP Value")
plot6
