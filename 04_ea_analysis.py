# -*- coding: utf-8 -*-
"""
EA MCF Study

@author: ying
"""

# Earnings Announcement Modified Causal Forest Analsyis


# import modules
import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plot
from mcf import ModifiedCausalForest

# PATH = "/Users/ ying/Desktop/master-thesis-code/data/"
PATH = "C:/Users/dfd20/Desktop/y/data/"

# choose the sue measure dataset for estimate
DATANAME = 'ea_prepared_bycqtr.csv'  # SUE based on time-series model
DATANAME = 'ea_forecast_bycqtr.csv' # SUE based on analyst forecast

# load in data using pandas
data = pd.read_csv(PATH + DATANAME)

#%% 
# Prepare CPI and GDP data
CPI = 'cpi_gdp.csv'
cpi_data = pd.read_csv(PATH + CPI)

data['date'] = pd.to_datetime(data['date'])
cpi_data['date'] = pd.to_datetime(cpi_data['date'])

data['month'] = data['date'].dt.to_period('M')
cpi_data['month'] = cpi_data['date'].dt.to_period('M')

# Merge to the original data using the 'month' column
merged_df = pd.merge(data, cpi_data.drop(columns='date'), on='month', how='left')

# Optional: drop 'month' column or move it
# merged_df = merged_df.drop(columns=['month', 'PERMNO'])
merged_df = merged_df.drop(columns=['month'])

#%%
# Part 1: Prepare the target variable

conditions = [
    (data['earning_qtile'] == 5),
    (data['earning_qtile'] == 1),
    (data['earning_qtile'] == 3)
   ]

# create a list of the values we want to assign for each condition
values = [2, 0, 1]

# create a new column and use np.select to assign values to it using our lists as arguments
merged_df['sue_qtile'] = np.select(conditions, values, default=np.nan)
# merged_df['sue_qtile'] = np.select(conditions, values)
# have a look at the distribution of treatment variable
merged_df['sue_qtile'].value_counts()
# check for the NA values 
merged_df.isna().sum()

# drop the rest qtiles
merged_df.dropna(subset=['sue_qtile'], inplace = True)

# factorized the date and symbol variable
merged_df['cqtr_index'], qtr_mapping = pd.factorize(merged_df['datacqtr'])
merged_df['symbol_index'], symbol_mapping = pd.factorize(merged_df['tsymbol'])


# prepare the seed for reproducebility
np.random.seed(72)
# randomly split the dataset into two parts
df_train = merged_df.sample(frac = 0.5)
df_test = merged_df.drop(df_train.index)

# specify input arguments for the modified random forest
D_NAME = ['sue_qtile']
X_NAME_ORD = ['pre_ann_rh_sahc_d', 'pre_ann_rh_sahc_w', 
              'pre_ann_taq_sanb_d', 'pre_ann_taq_sanb_w', 
              'pre_ann_car_d', 'pre_ann_car_w', 
              'pre_ann_std_d', 'pre_ann_std_w', 
              'pre_ann_turnover_d', 'pre_ann_turnover_w',
              'pre_ann_mkt_std_d', 'pre_ann_mkt_std_w', 
              'pre_ann_mcap', 'pre_ann_BtM', 'pre_ann_gp', 'pre_ann_inv', 'pre_ann_mom', 
              'CPI', 'GDP', 
              ]

C_NAME = ['cqtr_index']



#%%
## Robinhood's Abnormal Purchase react to Earnings Announcement 
# specify Robinhood input arguments
Y_NAME_D = ['post_ann_rh_sahc_d']
Y_NAME_W = ['post_ann_rh_sahc_w']

# Create an instance of the Modified Causal Forest model
rh_mcf_d = ModifiedCausalForest(
    var_y_name=Y_NAME_D,  # Outcome variable
    var_d_name=D_NAME,    # Treatment variable
    var_x_name_ord=X_NAME_ORD,  # Ordered covariates
    # var_x_name_unord=X_NAME_UNORD,  # Unordered covariate
    var_cluster_name=C_NAME,
    gen_panel_data=True,
    # gen_panel_in_rf=True,
    _int_show_plots=False  # Disable plots for faster performance
)

# Train the Modified Causal Forest on the training data
rh_mcf_d.train(df_train)
# Predict treatment effects using the model on prediction data
results_rh_d = rh_mcf_d.predict(df_test)
# Extract the dictionary of estimates
results_rh_d_dict = results_rh_d[0]
# Access the average treatenet effect (ATE)
ate_rh_d = results_rh_d_dict.get('ate')
print("Average Treatment Effect (ATE) for Robinhood daily:\n", ate_rh_d)

# Create an instance of the Modified Causal Forest model
rh_mcf_w = ModifiedCausalForest(
    var_y_name=Y_NAME_W,  # Outcome variable
    var_d_name=D_NAME,    # Treatment variable
    var_x_name_ord=X_NAME_ORD,  # Ordered covariates
    # var_x_name_unord=["x_unord0"],  # Unordered covariate
    var_cluster_name=C_NAME,
    gen_panel_data=True,
    _int_show_plots=False  # Disable plots for faster performance
)

# Train the Modified Causal Forest on the training data
rh_mcf_w.train(df_train)
# Predict treatment effects using the model on prediction data
results_rh_w = rh_mcf_w.predict(df_test)
# Extract the dictionary of estimates
results_rh_w_dict = results_rh_w[0]
# Access the average treatenet effect (ATE)
ate_rh_w = results_rh_w_dict.get('ate')
print("Average Treatment Effect (ATE) for Robinhood weekly:\n", ate_rh_w)

#%%
## TAQ Retail's Abnormal Purchase react to Earnings Announcement 
# specify taq input arguments
Y_NAME_D = ['post_ann_taq_sanb_d']
Y_NAME_W = ['post_ann_taq_sanb_w']

# Create an instance of the Modified Causal Forest model
taq_mcf_d = ModifiedCausalForest(
    var_y_name=Y_NAME_D,  # Outcome variable
    var_d_name=D_NAME,    # Treatment variable
    var_x_name_ord=X_NAME_ORD,  # Ordered covariates
    # var_x_name_unord=["x_unord0"],  # Unordered covariate
    var_cluster_name=C_NAME,
    gen_panel_data=True,
    # gen_panel_in_rf=True,
    _int_show_plots=False  # Disable plots for faster performance
)

# Train the Modified Causal Forest on the training data
taq_mcf_d.train(df_train)
# Predict treatment effects using the model on prediction data
results_taq_d = taq_mcf_d.predict(df_test)
# Extract the dictionary of estimates
results_taq_d_dict = results_taq_d[0]
# Access the average treatenet effect (ATE)
ate_taq_d = results_taq_d_dict.get('ate')
print("Average Treatment Effect (ATE) for TAQ daily:\n", ate_taq_d)

# Create an instance of the Modified Causal Forest model
taq_mcf_w = ModifiedCausalForest(
    var_y_name=Y_NAME_W,  # Outcome variable
    var_d_name=D_NAME,    # Treatment variable
    var_x_name_ord=X_NAME_ORD,  # Ordered covariates
    # var_x_name_unord=["x_unord0"],  # Unordered covariate
    var_cluster_name=C_NAME,
    gen_panel_data=True,
    # gen_panel_in_rf=True,
    _int_show_plots=False  # Disable plots for faster performance
)

# Train the Modified Causal Forest on the training data
taq_mcf_w.train(df_train)
# Predict treatment effects using the model on prediction data
results_taq_w = taq_mcf_w.predict(df_test)
# Extract the dictionary of estimates
results_taq_w_dict = results_taq_w[0]
# Access the average treatenet effect (ATE)
ate_taq_w = results_taq_w_dict.get('ate')
print("Average Treatment Effect (ATE) for TAQ weekly:\n", ate_taq_w)

#%%
## The Return react to Earnings Announcement 
# specify taq input arguments
Y_NAME_D = ['post_ann_car_d']
Y_NAME_W = ['post_ann_car_w']

# Create an instance of the Modified Causal Forest model
ret_mcf_d = ModifiedCausalForest(
    var_y_name=Y_NAME_D,  # Outcome variable
    var_d_name=D_NAME,    # Treatment variable
    var_x_name_ord=X_NAME_ORD,  # Ordered covariates
    # var_x_name_unord=["x_unord0"],  # Unordered covariate
    var_cluster_name=C_NAME,
    gen_panel_data=True,
    # gen_panel_in_rf=True,
    _int_show_plots=False  # Disable plots for faster performance
)

# Train the Modified Causal Forest on the training data
ret_mcf_d.train(df_train)
# Predict treatment effects using the model on prediction data
results_ret_d = ret_mcf_d.predict(df_test)
# Extract the dictionary of estimates
results_ret_d_dict = results_ret_d[0]
# Access the average treatenet effect (ATE)
ate_ret_d = results_ret_d_dict.get('ate')
print("Average Treatment Effect (ATE) for Return daily:\n", ate_ret_d)

# Create an instance of the Modified Causal Forest model
ret_mcf_w = ModifiedCausalForest(
    var_y_name=Y_NAME_W,  # Outcome variable
    var_d_name=D_NAME,    # Treatment variable
    var_x_name_ord=X_NAME_ORD,  # Ordered covariates
    # var_x_name_unord=["x_unord0"],  # Unordered covariate
    var_cluster_name=C_NAME,
    gen_panel_data=True,
    _int_show_plots=False  # Disable plots for faster performance
)

# Train the Modified Causal Forest on the training data
ret_mcf_w.train(df_train)
# Predict treatment effects using the model on prediction data
results_ret_w = ret_mcf_w.predict(df_test)
# Extract the dictionary of estimates
results_ret_w_dict = results_ret_w[0]
# Access the average treatenet effect (ATE)
ate_ret_w = results_ret_w_dict.get('ate')
print("Average Treatment Effect (ATE) for Return weekly:\n", ate_ret_w)