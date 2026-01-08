##### MRoS #####

library(pacman)
p_load(here, tidyverse, rsample)

MRoS <- read_csv(here("Data", "MRoS", "mros_output_20240607_nodupe.csv")) |>
  select(name, latitude, longitude, temp_air_idw_lapse_const, temp_air_idw_lapse_var,  temp_air_nearest_site_const, temp_air_nearest_site_var, temp_air_avg_obs, temp_air_nearest, plp, elevation, eco_level3, state, temp_air_flag, rh_flag, dist_temp_air_flag, dist_temp_dew_flag, closest_temp_air_flag, closest_temp_dew_flag, nstation_temp_air_flag, nstation_temp_dew_flag, pval_temp_air_flag, pval_temp_dew_flag, phase_flag) |>
  filter(state %in% c("Colorado", "Montana", "Nevada", "Oregon", "California", "Washington", "Utah", "Idaho", "New Mexico", "Wyoming", "Arizona", "South Dakota")) |>
  filter((name == "Rain" | name == "Snow") & plp <= 100) |>
  mutate(name_factor = as.factor(name)) |>
  drop_na(name, temp_air_idw_lapse_const, temp_air_idw_lapse_var, plp, elevation)
  
  
  
plot(MRoS$longitude, MRoS$latitude)

MRoS_flag_pass <- MRoS |>
  filter(temp_air_flag == "Pass" & rh_flag == "Pass" & dist_temp_air_flag == "Pass" & dist_temp_dew_flag == "Pass" & closest_temp_air_flag == "Pass" & closest_temp_dew_flag == "Pass" & nstation_temp_air_flag == "Pass" & nstation_temp_dew_flag == "Pass" & pval_temp_air_flag == "Pass" & pval_temp_dew_flag == "Pass" & phase_flag == "Pass")

### splitting data for logistic regression model
set.seed(802)
data_split <- initial_split(MRoS, prop = 0.75)
TrainingData_alldata <- training(data_split)
TestingData_alldata <- testing(data_split)

set.seed(802)
data_split2 <- initial_split(MRoS_flag_pass, prop = 0.75)
TrainingData_flagpass <- training(data_split2)
TestingData_flagpass <- testing(data_split2)


#fit logistic regression model
set.seed(802)
model_alldata<- glm(name_factor ~ temp_air_idw_lapse_var, family="binomial", data=TrainingData_alldata)
summary(model_alldata)
# plot(model)

#fit logistic regression model
set.seed(802)
model_flagpass <- glm(name_factor ~ temp_air_idw_lapse_var, family="binomial", data = TrainingData_flagpass)
summary(model_flagpass)

### predicting model
predictions_alldata<- predict(model_alldata, TestingData_alldata, type = "response")
plot(TestingData_alldata$temp_air_idw_lapse_var, predictions_alldata)
# hist(predictions_alldata)


predictions_flagpass <- predict(model_flagpass, TestingData_flagpass, type = "response")
plot(TestingData_flagpass$temp_air_idw_lapse_var, predictions_flagpass)
# hist(predictions_flagpass)
summary(predictions_flagpass)

TestingData_flagpass <- TestingData_flagpass |>
  mutate(prob_of_snow_mod1 = predictions_flagpass)

TestingData_flagpass_snowprob_sub50 <- TestingData_flagpass |>
  filter(prob_of_snow_mod1 < 0.5)

### fit logistic regression model
set.seed(802)
model_flagpass_const <- glm(name_factor ~ temp_air_idw_lapse_const, family="binomial", data = TrainingData_flagpass)
summary(model_flagpass_const)


predictions_flagpass_const <- predict(model_flagpass_const, TestingData_flagpass, type = "response")
plot(TestingData_flagpass$temp_air_idw_lapse_const, predictions_flagpass_const)


TestingData_flagpass <- TestingData_flagpass |>
  mutate(prob_of_snow_mod2 = predictions_flagpass_const)

TestingData_flagpass_const_snowprob_sub50 <- TestingData_flagpass |>
  filter(prob_of_snow_mod2 < 0.5)


### converting prediction probability to actual precip phase
TestingData_flagpass <- TestingData_flagpass |>
  mutate(phase_pred_mod1 = as.factor(ifelse(prob_of_snow_mod1 >= 50, "Snow", "Rain")),
         phase_pred_mod2 = as.factor(ifelse(prob_of_snow_mod2 >= 50, "Snow", "Rain")))



### fit logistic regression model
set.seed(802)
model_flagpass_nearest <- glm(name_factor ~ temp_air_nearest, family="binomial", data = TrainingData_flagpass)
summary(model_flagpass_nearest)


predictions_flagpass_nearest <- predict(model_flagpass_nearest, TestingData_flagpass, type = "response")
plot(TestingData_flagpass$temp_air_nearest, predictions_flagpass_nearest)


TestingData_flagpass <- TestingData_flagpass |>
  mutate(prob_of_snow_mod3 = predictions_flagpass_nearest)

TestingData_flagpass_nearest_snowprob_sub50 <- TestingData_flagpass |>
  filter(prob_of_snow_mod3 < 0.5)


### converting prediction probability to actual precip phase
TestingData_flagpass <- TestingData_flagpass |>
  mutate(phase_pred_mod3 = as.factor(ifelse(prob_of_snow_mod3 >= 50, "Snow", "Rain")))
