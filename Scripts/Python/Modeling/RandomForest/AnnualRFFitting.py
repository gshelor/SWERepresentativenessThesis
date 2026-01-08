##### Script for testing Python version of RF model #####
import pandas as pd
import geopandas as gpd
from pathlib import Path
from sklearn.model_selection import train_test_split, GridSearchCV, KFold, LeaveOneOut, ShuffleSplit
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, make_scorer
import numpy as np
import joblib
from tqdm import tqdm

## Configuration
# Set random seed for reproducibility
np.random.seed(802)

## File Paths
DATA_DIR = Path("Data")
SNOTEL_DATA_PATH = DATA_DIR / "SNOTEL" / "Combined" / "GIS" / "SnotelCombined_AnnualCovars_AnnualZeroPts.gpkg"
FITTED_MODELS_DIR = DATA_DIR / "FittedModels" / "Caret" / "RF"
FITTED_MODELS_DIR.mkdir(parents=True, exist_ok=True)

## Loading Data
print("Loading data...")
snotel_sf = gpd.read_file(SNOTEL_DATA_PATH)
snotel_df = pd.DataFrame(snotel_sf.drop(columns='geometry'))

columns_to_check_na = [
    "OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect",
    "OctApr_tminMean", "OctApr_tmaxMean", "OctApr_sradMean", "landcover_triclass",
    "OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "OctMay_tminMean",
    "OctMay_tmaxMean", "OctMay_sradMean", "SepMay_prcpSumCDMSum",
    "SepMay_tmeanCDMSum", "SepMay_tminMean", "SepMay_tmaxMean", "SepMay_sradMean"
]
snotel_df.dropna(subset=columns_to_check_na, inplace=True)

snotel_df['landcover'] = snotel_df['landcover'].astype(int)
snotel_df['landcover_triclass'] = snotel_df['landcover_triclass'].astype(int)
print("Data loaded and preprocessed.")

## Splitting data into training and testing datasets
print("Splitting data into training and testing sets...")
snotel_train, snotel_test = train_test_split(snotel_df, test_size=0.25, random_state=802)

train_y = snotel_train['peak_swe']
test_y = snotel_test['peak_swe']

feature_sets = {
    "OctApr_x": ["OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope",
                 "OctApr_tminMean", "OctApr_tmaxMean", "OctApr_sradMean", "landcover_triclass"],
    "OctAprAspect_x": ["OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect",
                       "OctApr_tminMean", "OctApr_tmaxMean", "OctApr_sradMean", "landcover_triclass"],
    "OctAprAspect_Nosrad_x": ["OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect",
                              "OctApr_tminMean", "OctApr_tmaxMean", "landcover_triclass"],
    "OctAprAspect_Nosradtmin_x": ["OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect",
                                  "OctApr_tmaxMean", "landcover_triclass"],
    "OctAprAspect_Nosradtmintmax_x": ["OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect",
                                      "landcover_triclass"],
    "OctAprAspect_Notmin_x": ["OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect",
                              "OctApr_tmaxMean", "OctApr_sradMean", "landcover_triclass"],
    "OctAprAspect_Notmintmax_x": ["OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect",
                                  "OctApr_sradMean", "landcover_triclass"],
    "OctAprAspect_Notmax_x": ["OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "elevation", "slope", "aspect",
                              "OctApr_tminMean", "OctApr_sradMean", "landcover_triclass"],
    "OctAprDecFeb_x": ["OctApr_prcpSumCDMSum", "OctApr_tmeanCDMSum", "DecFeb_prcpSumCDMSum",
                       "DecFeb_tmeanCDMSum", "elevation", "slope", "OctApr_tminMean", "OctApr_tmaxMean",
                       "OctApr_sradMean", "landcover_triclass"],
    "OctMay_x": ["OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "elevation", "slope",
                 "OctMay_tminMean", "OctMay_tmaxMean", "OctMay_sradMean", "landcover_triclass"],
    "OctMayAspect_Nosradtmin_x": ["OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "elevation", "slope", "aspect",
                                  "OctMay_tmaxMean", "landcover_triclass"],
    "OctMayAspect_Nosradtmintmax_x": ["OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "elevation", "slope", "aspect",
                                      "landcover_triclass"],
    "OctMayAspect_Notmintmax_x": ["OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "elevation", "slope", "aspect",
                                  "OctMay_sradMean", "landcover_triclass"],
    "OctMayDecFeb_x": ["OctMay_prcpSumCDMSum", "OctMay_tmeanCDMSum", "DecFeb_prcpSumCDMSum",
                       "DecFeb_tmeanCDMSum", "elevation", "slope", "OctMay_tminMean", "OctMay_tmaxMean",
                       "OctMay_sradMean", "landcover_triclass"],
    "SepMay_x": ["SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "elevation", "slope",
                 "SepMay_tminMean", "SepMay_tmaxMean", "SepMay_sradMean", "landcover_triclass"],
    "SepMayDecFeb_x": ["SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "DecFeb_prcpSumCDMSum",
                       "DecFeb_tmeanCDMSum", "elevation", "slope", "SepMay_tminMean", "SepMay_tmaxMean",
                       "SepMay_sradMean", "landcover_triclass"],
    "SepMayAspect_Nosradtmin_x": ["SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "elevation", "slope", "aspect",
                                  "SepMay_tmaxMean", "landcover_triclass"],
    "SepMayAspect_Nosradtmintmax_x": ["SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "elevation", "slope", "aspect",
                                      "landcover_triclass"],
    "SepMayAspect_Notmintmax_x": ["SepMay_prcpSumCDMSum", "SepMay_tmeanCDMSum", "elevation", "slope", "aspect",
                                  "SepMay_sradMean", "landcover_triclass"]
}

train_df_dict = {name: snotel_train[cols] for name, cols in feature_sets.items()}
test_df_dict = {name: snotel_test[cols] for name, cols in feature_sets.items()}

train_df_list = list(train_df_dict.values())
test_df_list = list(test_df_dict.values())
model_names = list(feature_sets.keys())
print("Data split into training and testing sets.")



## Hyperparameter Tuning with GridSearchCV


# Custom RMSE scorer for GridSearchCV
def rmse_scorer(y_true, y_pred):
    # Ensure predictions are non-negative
    y_pred[y_pred < 0] = 0
    return np.sqrt(mean_squared_error(y_true, y_pred))

rmse_scorer = make_scorer(rmse_scorer, greater_is_better=False) # 'greater_is_better=False' means lower RMSE is better

def tune_and_train_rf_model_gridsearch(X_train, y_train, cv_strategy, model_name):
    # Determine the mtry (max_features) range
    num_features = X_train.shape[1]
    param_grid = {'max_features': list(range(2, num_features)) if num_features > 2 else [1.0]}

    # Base Random Forest model with n_jobs=12
    rf = RandomForestRegressor(n_estimators=2500, random_state=802, n_jobs=12, oob_score=True)

    # Define the cross-validation strategy
    if cv_strategy == "oob":
        print(f"Manually tuning {model_name} with OOB error...")
        best_oob_rmse = float('inf')
        best_params = None
        best_model = None
        results_df = pd.DataFrame(columns=['max_features', 'RMSE', 'Rsquared'])

        for max_features in param_grid['max_features']:
            # RandomForestRegressor with n_jobs=12
            model = RandomForestRegressor(n_estimators=2500, max_features=max_features,
                                          random_state=802, oob_score=True, n_jobs=12)
            model.fit(X_train, y_train)
            
            # OOB predictions for RMSE (ensure non-negative)
            oob_preds = model.oob_prediction_
            oob_preds[oob_preds < 0] = 0
            
            oob_rmse = np.sqrt(mean_squared_error(y_train, oob_preds))
            oob_r2 = model.oob_score_

            current_results = pd.DataFrame([{'max_features': max_features, 'RMSE': oob_rmse, 'Rsquared': oob_r2}])
            results_df = pd.concat([results_df, current_results], ignore_index=True)

            if oob_rmse < best_oob_rmse:
                best_oob_rmse = oob_rmse
                best_params = {'max_features': max_features}
                best_model = model

        return {'bestTune': best_params, 'results': results_df, 'model': best_model}

    else:
        print(f"Tuning {model_name} with GridSearchCV ({cv_strategy})...")
        if cv_strategy == "cv":
            cv_obj = KFold(n_splits=10, shuffle=True, random_state=802)
        elif cv_strategy == "LOOCV":
            cv_obj = LeaveOneOut()
        elif cv_strategy == "boot": # Simulate bootstrapping with ShuffleSplit for CV
            cv_obj = ShuffleSplit(n_splits=10, test_size=0.25, random_state=802)

        # GridSearchCV with n_jobs=12
        grid_search = GridSearchCV(
            estimator=rf,
            param_grid=param_grid,
            cv=cv_obj,
            scoring=rmse_scorer,
            n_jobs=12, # Changed n_jobs to 12
            verbose=1,
            return_train_score=True
        )

        grid_search.fit(X_train, y_train)

        best_params = grid_search.best_params_
        best_score = -grid_search.best_score_
        
        results_df = pd.DataFrame(grid_search.cv_results_)
        results_df = results_df[['param_max_features', 'mean_test_score', 'std_test_score']]
        results_df.rename(columns={
            'param_max_features': 'max_features',
            'mean_test_score': 'RMSE',
            'std_test_score': 'RMSE_SD'
        }, inplace=True)
        results_df['RMSE'] = -results_df['RMSE']
        # results_df['Rsquared'] = grid_search.cv_results_['mean_train_score']
        # results_df['Rsquared'] = results_df['mean_train_score']
        
        best_model = grid_search.best_estimator_

        return {'bestTune': best_params, 'results': results_df, 'model': best_model}

## Training Models

# Training models resampled using OOB error
print("\n--- Training models with OOB error ---")
OOBModels = []
for i, train_data in enumerate(tqdm(train_df_list, desc="Fitting OOB Models")):
    model_info = tune_and_train_rf_model_gridsearch(train_data, train_y, "oob", model_names[i])
    OOBModels.append(model_info)

# Evaluating test metrics for OOBModels
print("\n--- Evaluating test metrics for OOBModels ---")
for i, model_info in enumerate(OOBModels):
    print(f"\nModel {i+1}: {model_names[i]}")
    temp_model = model_info['model']
    print("Best Tune:", model_info['bestTune'])
    print("Results:\n", model_info['results'])
    
    RF_model_preds = temp_model.predict(test_df_list[i])
    RF_model_preds[RF_model_preds < 0] = 0
    test_rmse = np.sqrt(mean_squared_error(test_y, RF_model_preds))
    print(f"Test RMSE: {test_rmse}")

# Training models resampled using bootstrapping (simulated with ShuffleSplit)
# print("\n--- Training models with bootstrapping (simulated) ---")
# BootModels = []
# for i, train_data in enumerate(tqdm(train_df_list, desc="Fitting Bootstrapped Models")):
#     model_info = tune_and_train_rf_model_gridsearch(train_data, train_y, "boot", model_names[i])
#     BootModels.append(model_info)

# # Evaluating test metrics for BootModels
# print("\n--- Evaluating test metrics for BootModels ---")
# for i, model_info in enumerate(BootModels):
#     print(f"\nModel {i+1}: {model_names[i]}")
#     temp_model = model_info['model']
#     print("Best Tune:", model_info['bestTune'])
#     print("Results:\n", model_info['results'])
    
#     RF_model_preds = temp_model.predict(test_df_list[i])
#     RF_model_preds[RF_model_preds < 0] = 0
#     test_rmse = np.sqrt(mean_squared_error(test_y, RF_model_preds))
#     print(f"Test RMSE: {test_rmse}")

# Training models resampled using cross-validation (KFold)
print("\n--- Training models with cross-validation ---")
CVModels = []
for i, train_data in enumerate(tqdm(train_df_list, desc="Fitting CV Models")):
    model_info = tune_and_train_rf_model_gridsearch(train_data, train_y, "cv", model_names[i])
    CVModels.append(model_info)

# Evaluating test metrics for CVModels
print("\n--- Evaluating test metrics for CVModels ---")
for i, model_info in enumerate(CVModels):
    print(f"\nModel {i+1}: {model_names[i]}")
    temp_model = model_info['model']
    print("Best Tune:", model_info['bestTune'])
    print("Results:\n", model_info['results'])
    
    RF_model_preds = temp_model.predict(test_df_list[i])
    RF_model_preds[RF_model_preds < 0] = 0
    test_rmse = np.sqrt(mean_squared_error(test_y, RF_model_preds))
    print(f"Test RMSE: {test_rmse}")

## Saving Models
print("\nSaving fitted models...")
#joblib.dump(OOBModels, FITTED_MODELS_DIR / "OOBModelList.pkl", compress=True)
#joblib.dump(BootModels, FITTED_MODELS_DIR / "BootModelList.pkl", compress=True)
#joblib.dump(CVModels, FITTED_MODELS_DIR / "CVModelList.pkl", compress=True)
print("Models saved.")