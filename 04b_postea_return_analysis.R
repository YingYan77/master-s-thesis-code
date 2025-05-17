################################################################################
## author:    Ying Yan
## file name: 05b_PostEA_Return_Analysis.R
## summary:   this file aims to study the Earnings Accountment effect on Robinhood 
##            investors' trading behaviour using IBES Detail, CRSP and RH data. 
## date:      28 Feb 2024
################################################################################

## empty memory (!)
rm(list=ls())

## load packages
library(dplyr)
library(data.table)
library(tidyverse)

## regression
library("sandwich")
library("lmtest")
library(fixest)

## latex
library(knitr)
library(kableExtra)

## read the function library
source("00_functions.R")

################################################################################
############## STEP 1: COMPUTE STANDARDIZED UNEXPECTED EARNINGS ################
################################################################################
## Load the quarterly earning estimation and EA IBES dataset 
ibes_actual <- fread("data/ibes_compustat_merged.csv")

## check on the number of estimations/forecasts for each EA event  
## we are interested on the EA date effect on the stock return & changes in holding 
## so we use the actual earning announcement report date as the aggregation date 
ea_event <- ibes_actual[, .(length = .N), by = .(anndats, tsymbol)]

summary(ea_event$length) # the average estimates for each event are 10
table(ea_event$length) # 10% of the events have only one estimate. 
## check on the number of EA event for each stock
EAs <- ea_event[, .(dates = .N), by = .(tsymbol)]

table(EAs$dates) # most stocks have 9-10 quarterly EAs

## check the number of stocks left
length(unique(ibes_actual$tsymbol))
## in order to match the earning event date, we have to adjust the event date in 
## IBES dataset, i.e., the overnight changes from t-1 (16:00) to t (9:00) belongs 
## to t, which means that any event happens between t-1 16:00 to t-1 24:00 is
## assigned to date t. 
ibes_actual$anntims <- hms::hms(seconds = ibes_actual$anntims)

## load the NYSE business day calendar
bizdays::load_quantlib_calendars(c("UnitedStates/NYSE"), from = "2016-01-01", to = "2020-12-30" )

# Adjust date for time >= 16:00
ibes_actual[, anndats_adj := anndats]  # Start by copying the original date
ibes_actual[anntims > hms::as_hms("16:00:00"), anndats_adj := bizdays::add.bizdays(anndats, 1, "QuantLib/UnitedStates/NYSE")] # Adjust to next business day

# 1123 EAs happen during the trading hour, which account for 1.8% of the cases
nrow(ibes_actual[anntims < hms::as_hms("16:00:00") & anntims > hms::as_hms("9:30:00"), ]) / nrow(ibes_actual)
# focus on off-trading hour earnings announcement events
ibes_overnight <- ibes_actual[anntims > hms::as_hms("16:00:00") | anntims < hms::as_hms("9:30:00"), ] 

## next, merge with crsp dataset to compute the SUE
crsp <- fread("data/crsp_volatility_compustat_merged.csv")

# crsp[, paste0("ret_lead_", 0:n_lags) := shift(ret, n = 0:n_lags, type = "lead"), by = .(tsymbol)]
# crsp[,  paste0("ret_lag_", 0:n_lags)  := shift(ret, n = 0:n_lags, type = "lag"), by = .(tsymbol)]

ibes_overnight <- ibes_overnight[!is.na(sue)]
## merge robinhood and sue data 
eas_all <- merge(crsp, ibes_overnight, by.x = c("date", "tsymbol"), by.y = c("anndats_adj","tsymbol"), all.x = TRUE)
length(unique(eas_all[!is.na(sue)]$tsymbol))

## create the dummy variable that indicate whether this is an EA event
eas_all[, ea := fifelse(is.na(value), 0, 1)]
table(eas_all$ea) # after merge, there are 20107 events left.

################################################################################
################### STEP 2: COMPUTE CONTROL VARIABLES FOR MCF ##################
################################################################################

## read the taq daily aggregated dataset, cleaned from 00_TAQ_RH_CRSP_Merged.R
rh_taq <- fread("data/rh_taq.csv")

## merge Robinhood, TAQ and CRSP
rh_taq_crsp_ea <- merge(rh_taq, eas_all, by = c("date", "tsymbol"))


## calculate the cumulative return from t+1 to t+5
for (k in c(1,5,10,20)) {
  rh_taq_crsp_ea[, paste0("post_ann_ret_", k) := 
                   frollapply(shift(1 + ret, 1, type = "lead"), k, prod, fill = NA, align = "left"),
                 by = .(tsymbol)]
  
  rh_taq_crsp_ea[, paste0("post_ann_mktret_", k) := 
                   frollapply(shift(1 + mkt_ret, 1, type = "lead"), k, prod, fill = NA, align = "left"),
                 by = .(tsymbol)]
  
  rh_taq_crsp_ea[, paste0("post_ann_car_", k) := 
                   get(paste0("post_ann_ret_", k)) - get(paste0("post_ann_mktret_", k)),
                 by = .(tsymbol)]
}

## get the earnings announcement event dates 
ea <- rh_taq_crsp_ea[ea == 1,]
## another subset for the rest of dates
nea <- rh_taq_crsp_ea[ea == 0,]

## find the quintile interval for the changes in eps
ea[, earning_qtile := findInterval(sue, 
                                   quantile(sue, prob=c(.2,.4,.6,.8), na.rm = T), left.open = T) + 1, 
   by = .(datacqtr)] 


names(ea)

################################################################################
#################### STEP 3: CHECK ROBINHOOD TRADES AND EPS ####################
################################################################################

## create the list of CAR time period
car_vars <- c("post_ann_car_1", "post_ann_car_5", "post_ann_car_10", "post_ann_car_20")

## get the tested average CAR table for each quintile 
sue_ret <- get_sue_holding_ret_result(varlist = car_vars, data = ea, invlist = c("sa_holding_change", "sa_diff_trd"))
sue_ret_sapdiff <- get_sue_qtile_diff_ret(varlist = car_vars, data = ea, invlist = c("sa_holding_change", "sa_diff_trd"))

sue_ret <- sue_ret |>
  select(-p.value) |> 
  bind_rows(sue_ret_sapdiff)

## print out the result
rh_sue_ret <- sue_ret |> filter(investor == "sa_holding_change")
kable(reshape_sue_timelag(rh_sue_ret), format = "latex", booktabs = TRUE, align = "c")

taq_sue_ret <- sue_ret |> filter(investor == "sa_diff_trd")
kable(reshape_sue_timelag(taq_sue_ret), format = "latex", booktabs = TRUE, align = "c")


