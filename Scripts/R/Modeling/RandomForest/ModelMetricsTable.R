##### Making gt tables of model metrics to match gt tables used in rest of thesis #####

### loading packages, reading in data
library(pacman)
p_load(here, tidyverse, gt, gtExtras)


### Alaska Model Metrics
AKModelMetrics <- read_csv(here("Outputs", "CSVs", "ModelMetrics", "AKModelMetrics.csv")) |>
  filter(resampling_method != "LOOCV") |>
  # filter(is_min_test_RMSE != 1) |>
  drop_na(best_model_R2, testing_RMSE, testing_MAE) |>
  arrange(desc(resampling_method)) |>
  select(AOI, resampling_method, num_folds, vars_included, best_model_R2, testing_RMSE, testing_MAE, is_min_test_RMSE) |>
  mutate(resampling_method = case_when(resampling_method == "oob" ~ "OOB",
                                       TRUE ~ resampling_method))

### gonna make separate tables for each resampling method because it looks weird in latex all stretched out
AKOOBModelMetrics <- AKModelMetrics |>
  filter(resampling_method == "OOB")
AKCVModelMetrics <- AKModelMetrics |>
  filter(resampling_method == "CV")

### CONUS Model Metrics
CONUSModelMetrics <- read_csv(here("Outputs", "CSVs", "ModelMetrics", "CONUSModelMetrics.csv")) |>
  filter(resampling_method != "LOOCV") |>
  filter(resampling_method != "CV") |>
  # filter(is_min_test_RMSE != 1) |>
  drop_na(best_model_R2, testing_RMSE, testing_MAE) |>
  arrange(desc(resampling_method)) |>
  select(AOI, resampling_method, num_folds, vars_included, best_model_R2, testing_RMSE, testing_MAE, is_min_test_RMSE) |>
  mutate(resampling_method = case_when(resampling_method == "oob" ~ "OOB",
                                       TRUE ~ resampling_method),
         AOI = "WUS")

### final model for each region
AKFinalModel <- AKModelMetrics |>
  filter(is_min_test_RMSE == 1)

CONUSFinalModel <- CONUSModelMetrics |>
  filter(is_min_test_RMSE == 1)

FinalModels <- rbind(CONUSFinalModel, AKFinalModel)


##### Making gt tables #####
### Alaska OOB
AKOOBModelMetrics_gt <- AKOOBModelMetrics |>
  gt() |> # use 'gt' to make an awesome table...
  gt_theme_guardian() |>
  # tab_header(
  #   title = paste(), # ...with this title
  #   subtitle = "")  |>  # and this subtitle
  ##tab_style(style = cell_fill("bisque"),
  ##        locations = cells_body()) |>  # add fill color to table
  fmt_number( # A column (numeric data)
    columns = c(best_model_R2, testing_RMSE, testing_MAE), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  # fmt_number( # A column (numeric data)
  #   columns = c(), # What column variable? FinalVoATop25$VoA_Rating
  #   decimals = 0 # With four decimal places
  # ) |>
  data_color( # Update cell colors, testing different color palettes
    columns = c(best_model_R2),
    fn = scales::col_numeric( # <- bc it's numeric
      palette = brewer.pal(11, "RdYlGn"), # A color scheme (gradient)
      domain = c(), # Column scale endpoints
      reverse = FALSE
    )
  ) |>
  data_color( # Update cell colors, testing different color palettes
    columns = c(testing_RMSE, testing_MAE), # ...for dose column
    # direction = "row",
    fn = scales::col_numeric( # <- bc it's numeric
      palette = brewer.pal(11, "RdYlGn"), # A color scheme (gradient)
      domain = c(), # Column scale endpoints
      reverse = TRUE
    )
  ) |>
  cols_label(AOI = "Region", vars_included = "Covariates", resampling_method = "Resampling", num_folds = "# Folds", best_model_R2 = html("R<sup>2</sup>"), testing_RMSE = "RMSE", testing_MAE = "MAE") |> # Update labels
  # cols_move_to_end(columns = c("mean_AbsSHAP_DecFeb_prcpSumCDMSum", "mean_AbsSHAP_DecFeb_tmeanCDMSum", "mean_AbsSHAP_elevation", "mean_AbsSHAP_slope", "mean_AbsSHAP_aspect", "mean_AbsSHAP_landcover_triclass")) |>
  cols_hide(c(is_min_test_RMSE, num_folds, AOI, resampling_method)) |>
  opt_horizontal_padding(scale = 3) |>
  opt_vertical_padding(scale = 0.5)

AKOOBModelMetrics_gt
### saving table
AKOOBModelMetrics_gt |>
  gtsave(
    "AKOOBModelMetrics.png", expand = 5,
    path = here("Outputs", "ChartsGraphsTables", "ModelMetrics")
  )

### Alaska CV models
AKCVModelMetrics_gt <- AKCVModelMetrics |>
  gt() |> # use 'gt' to make an awesome table...
  gt_theme_guardian() |>
  # tab_header(
  #   title = paste(), # ...with this title
  #   subtitle = "")  |>  # and this subtitle
  ##tab_style(style = cell_fill("bisque"),
  ##        locations = cells_body()) |>  # add fill color to table
  fmt_number( # A column (numeric data)
    columns = c(best_model_R2, testing_RMSE, testing_MAE), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  # fmt_number( # A column (numeric data)
  #   columns = c(), # What column variable? FinalVoATop25$VoA_Rating
  #   decimals = 0 # With four decimal places
  # ) |>
  data_color( # Update cell colors, testing different color palettes
    columns = c(best_model_R2),
    fn = scales::col_numeric( # <- bc it's numeric
      palette = brewer.pal(11, "RdYlGn"), # A color scheme (gradient)
      domain = c(), # Column scale endpoints
      reverse = FALSE
    )
  ) |>
  data_color( # Update cell colors, testing different color palettes
    columns = c(testing_RMSE, testing_MAE), # ...for dose column
    # direction = "row",
    fn = scales::col_numeric( # <- bc it's numeric
      palette = brewer.pal(11, "RdYlGn"), # A color scheme (gradient)
      domain = c(), # Column scale endpoints
      reverse = TRUE
    )
  ) |>
  cols_label(AOI = "Region", vars_included = "Covariates", resampling_method = "Resampling", num_folds = "# Folds", best_model_R2 = html("R<sup>2</sup>"), testing_RMSE = "RMSE", testing_MAE = "MAE") |> # Update labels
  # cols_move_to_end(columns = c("mean_AbsSHAP_DecFeb_prcpSumCDMSum", "mean_AbsSHAP_DecFeb_tmeanCDMSum", "mean_AbsSHAP_elevation", "mean_AbsSHAP_slope", "mean_AbsSHAP_aspect", "mean_AbsSHAP_landcover_triclass")) |>
  cols_hide(c(is_min_test_RMSE, AOI, resampling_method)) |>
  opt_horizontal_padding(scale = 3) |>
  opt_vertical_padding(scale = 0.5)

AKCVModelMetrics_gt
### saving table
AKCVModelMetrics_gt |>
  gtsave(
    "AKCVModelMetrics.png", expand = 5,
    path = here("Outputs", "ChartsGraphsTables", "ModelMetrics")
  )


### CONUS
CONUSModelMetrics_gt <- CONUSModelMetrics |>
  gt() |> # use 'gt' to make an awesome table...
  gt_theme_guardian() |>
  # tab_header(
  #   title = paste(), # ...with this title
  #   subtitle = "")  |>  # and this subtitle
  ##tab_style(style = cell_fill("bisque"),
  ##        locations = cells_body()) |>  # add fill color to table
  fmt_number( # A column (numeric data)
    columns = c(best_model_R2, testing_RMSE, testing_MAE), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  # fmt_number( # A column (numeric data)
  #   columns = c(), # What column variable? FinalVoATop25$VoA_Rating
  #   decimals = 0 # With four decimal places
  # ) |>
  data_color( # Update cell colors, testing different color palettes
    columns = c(best_model_R2),
    fn = scales::col_numeric( # <- bc it's numeric
      palette = brewer.pal(11, "RdYlGn"), # A color scheme (gradient)
      domain = c(), # Column scale endpoints
      reverse = FALSE
    )
  ) |>
  data_color( # Update cell colors, testing different color palettes
    columns = c(testing_RMSE, testing_MAE), # ...for dose column
    # direction = "row",
    fn = scales::col_numeric( # <- bc it's numeric
      palette = brewer.pal(11, "RdYlGn"), # A color scheme (gradient)
      domain = c(), # Column scale endpoints
      reverse = TRUE
    )
  ) |>
  cols_label(AOI = "Region", vars_included = "Covariates", resampling_method = "Resampling", num_folds = "# Folds", best_model_R2 = html("R<sup>2</sup>"), testing_RMSE = "RMSE", testing_MAE = "MAE") |> # Update labels
  # cols_move_to_end(columns = c("mean_AbsSHAP_DecFeb_prcpSumCDMSum", "mean_AbsSHAP_DecFeb_tmeanCDMSum", "mean_AbsSHAP_elevation", "mean_AbsSHAP_slope", "mean_AbsSHAP_aspect", "mean_AbsSHAP_landcover_triclass")) |>
  cols_hide(c(is_min_test_RMSE, AOI, resampling_method, num_folds)) |>
  opt_horizontal_padding(scale = 3) |>
  opt_vertical_padding(scale = 0.5)

CONUSModelMetrics_gt
### saving table
CONUSModelMetrics_gt |>
  gtsave(
    "CONUSModelMetrics.png", expand = 5,
    path = here("Outputs", "ChartsGraphsTables", "ModelMetrics")
  )


### Final Model table
FinalModelMetrics_gt <- FinalModels |>
  gt() |> # use 'gt' to make an awesome table...
  gt_theme_guardian() |>
  # tab_header(
  #   title = paste(), # ...with this title
  #   subtitle = "")  |>  # and this subtitle
  ##tab_style(style = cell_fill("bisque"),
  ##        locations = cells_body()) |>  # add fill color to table
  fmt_number( # A column (numeric data)
    columns = c(best_model_R2, testing_RMSE, testing_MAE), # What column variable? FinalVoATop25$VoA_Rating
    decimals = 2 # With four decimal places
  ) |>
  # fmt_number( # A column (numeric data)
  #   columns = c(), # What column variable? FinalVoATop25$VoA_Rating
  #   decimals = 0 # With four decimal places
  # ) |>
  # data_color( # Update cell colors, testing different color palettes
  #   columns = c(best_model_R2),
  #   fn = scales::col_numeric( # <- bc it's numeric
  #     palette = brewer.pal(11, "RdYlGn"), # A color scheme (gradient)
  #     domain = c(), # Column scale endpoints
  #     reverse = FALSE
  #   )
  # ) |>
  # data_color( # Update cell colors, testing different color palettes
  #   columns = c(testing_RMSE, testing_MAE), # ...for dose column
  #   # direction = "row",
  #   fn = scales::col_numeric( # <- bc it's numeric
  #     palette = brewer.pal(11, "RdYlGn"), # A color scheme (gradient)
  #     domain = c(), # Column scale endpoints
  #     reverse = TRUE
  #   )
  # ) |>
  cols_label(AOI = "Region", vars_included = "Covariates", resampling_method = "Resampling", num_folds = "# Folds", best_model_R2 = html("R<sup>2</sup>"), testing_RMSE = "RMSE", testing_MAE = "MAE") |> # Update labels
  # cols_move_to_end(columns = c("mean_AbsSHAP_DecFeb_prcpSumCDMSum", "mean_AbsSHAP_DecFeb_tmeanCDMSum", "mean_AbsSHAP_elevation", "mean_AbsSHAP_slope", "mean_AbsSHAP_aspect", "mean_AbsSHAP_landcover_triclass")) |>
  cols_hide(c(is_min_test_RMSE)) |>
  opt_vertical_padding(scale = 2)

FinalModelMetrics_gt
### saving table
FinalModelMetrics_gt |>
  gtsave(
    "FinalModelMetrics.png", expand = 5,
    path = here("Outputs", "ChartsGraphsTables", "ModelMetrics")
  )
