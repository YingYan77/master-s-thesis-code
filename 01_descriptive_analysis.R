################################################################################
## author:    Ying Yan
## file name: 01_Descriptive_Analysis.R
## summary:   this file aims to study generate the descriptive analysis table 
##            used in the thesis. 
## date:      21 April 2025
################################################################################


## empty memory (!)
rm(list=ls())

## load packages
library(data.table)
library(stargazer)

## read the taq daily aggregated dataset, cleaned from 00_TAQ_RH_CRSP_Merged.R
rh_taq <- fread("data/rh_taq.csv")
crsp <- fread("data/crsp_volatility_compustat_merged.csv")

rh_taq_crsp_fundas_daily <- merge(rh_taq, crsp, by = c("date", "tsymbol"))

################################################################################
######################## STEP 1: DESCRIPRIVE STATISTICS ########################
################################################################################

names(rh_taq_crsp_fundas_daily)

rh_taq_crsp_fundas_daily[, mcap_mil := mcap / 1000000]

rh_vars <- c("users_open", "users_close", "holding_change", "sa_holding_change")

taq_vars <- c("trd_sell", "trd_buy", "diff_trd", "sa_diff_trd")

sym_vars <- c("ret", "volatility", "mcap_mil", "mkt_ret", "mkt_std_dev", "momentum", "BtM", "gp", "inv")

stargazer(rh_taq_crsp_fundas_daily[, ..rh_vars], type = "text", median = TRUE, iqr = TRUE)
stargazer(rh_taq_crsp_fundas_daily[, ..rh_vars], type = "latex", median = TRUE, iqr = TRUE)

stargazer(rh_taq_crsp_fundas_daily[, ..taq_vars], type = "text", median = TRUE, iqr = TRUE)
stargazer(rh_taq_crsp_fundas_daily[, ..taq_vars], type = "latex", median = TRUE, iqr = TRUE)

stargazer(rh_taq_crsp_fundas_daily[, ..sym_vars], type = "text", median = TRUE, iqr = TRUE)
stargazer(rh_taq_crsp_fundas_daily[, ..sym_vars], type = "latex", median = TRUE, iqr = TRUE)

