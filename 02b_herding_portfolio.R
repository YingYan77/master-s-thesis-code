## empty memory (!)
rm(list=ls())

## load packages
library(dplyr)
library(data.table)
library(tidyverse)
library(zoo)

library("sandwich")
library("lmtest")

## latex
library(knitr)
library(kableExtra)

## read the function library
source("00_functions.R")

################################################################################
############################## DATASET PREPARTION ##############################
################################################################################

## read data
rh_taq <- fread("data/rh_taq.csv")

## CRSP data
crsp <- fread("data/crsp_volatility_compustat_merged.csv")
## calculate the t+1 to t+5 return for portfolio computation
n_lags = 5
crsp[, paste0("ret_lead_", 0:n_lags) := shift(ret, n = 0:n_lags, type = "lead"), by = .(tsymbol)]
## calculate the cumulative return from t+1 to t+5
crsp[, cumret_5d := frollapply(shift(1 + ret, 1, type = "lead"), 5, prod, fill = NA, align = "left"), by = .(tsymbol)]
## calculate the cumulative return from t+1 to t+10
crsp[, cumret_10d := frollapply(shift(1 + ret, 1, type = "lead"), 10, prod, fill = NA, align = "left"), by = .(tsymbol)]
## calculate the cumulative return from t+1 to t+20
crsp[, cumret_20d := frollapply(shift(1 + ret, 1, type = "lead"), 20, prod, fill = NA, align = "left"), by = .(tsymbol)]
# for (i in 1:10) { 
#   crsp[, paste0("cumret_",i, "d") := frollapply(shift(1 + ret, 1, type = "lead"), i, prod, fill = NA, align = "left"), by = .(tsymbol)]
#   }

## merge Robinhood and CRSP
rh_taq_crsp <- merge(rh_taq, crsp, by = c("date", "tsymbol"))

## FF5 daily factors 
ff5 <- read.csv("data/ff5.csv", skip=3)
ff5$Date  <- as.Date(as.character(ff5$X), format = "%Y%m%d")
ff5 <- dplyr::select(ff5, -c(X))

## Momontum factor 
mom <- read.csv("data/mom.csv", skip=13)
mom$Date  <- as.Date(mom$X, format = "%Y%m%d")
mom <- dplyr::select(mom, -c(X))

ff6 <- merge(ff5, mom, by = c("Date"))

################################################################################
############################ ASSET PRICING ANALYSIS ############################
################################################################################

# create the 5% interval
rh_taq_crsp <- get_qtile_data(rh_taq_crsp, c("sa_holding_change", "sa_diff_trd"), extreme = .05)
# run the whole anlysis
alphas.5 <- get_ret_alpha(
  rh_taq_crsp,
  varlists = c("sa_holding_change", "sa_diff_trd")
)


# create 0.5% extreme qtile variable for target measure
rh_taq_crsp <- get_qtile_data(rh_taq_crsp, c("sa_holding_change", "sa_diff_trd"), extreme = .005)
# run the whole anlysis
alphas.05 <- get_ret_alpha(rh_taq_crsp, varlists = c("sa_holding_change", "sa_diff_trd"))

################################################################################
########################## ADJUST FORMATTING TO PRINT ##########################
################################################################################

alphas.5.equal <- alphas.5 |> 
  filter(weight == 'equal' & day %in% c(1,5,10) & model != 'Excess') |> 
  mutate(
    col_id = paste(type, model, sep = "."),
    row_id = paste0("day", day, ".decile", decile)
  ) |> 
  select(row_id, col_id, intercept) |> 
  pivot_wider(names_from = col_id, values_from = intercept)

## formatting the standard error
se.5.equal <- alphas.5 |> 
  filter(weight == 'equal' & day %in% c(1,5,10) & model != 'Excess') |> 
  mutate(
    se = paste0("(", round(std.error, 3), ")"),
    col_id = paste(type, model, sep = "."),
    row_id = paste0("day", day, ".decile", decile)
  ) |> 
  select(row_id, col_id, se) |> 
  pivot_wider(names_from = col_id, values_from = se)


## print out 
table.5.equal <- rbind(alphas.5.equal, se.5.equal)
table.5.equal <- table.5.equal |> arrange(row_id)
kable(table.5.equal, format = "latex", booktabs = TRUE, align = "c")


alphas.5.weight <- alphas.5 |> 
  filter(weight == 'holding_changes' & day %in% c(1,5,10) & model != 'Excess') |> 
  mutate(
    col_id = paste(type, model, sep = "."),
    row_id = paste0("day", day, ".decile", decile)
  ) |> 
  select(row_id, col_id, intercept) |> 
  pivot_wider(names_from = col_id, values_from = intercept)

# kable(alphas.5.equal, format = "latex", booktabs = TRUE, align = "c")
# kable(alphas.5.weight, format = "latex", booktabs = TRUE, align = "c")


se.5.weight <- alphas.5 |> 
  filter(weight == 'holding_changes' & day %in% c(1,5,10) & model != 'Excess') |> 
  mutate(
    se = paste0("(", round(std.error, 3), ")"),
    col_id = paste(type, model, sep = "."),
    row_id = paste0("day", day, ".decile", decile)
  ) |> 
  select(row_id, col_id, se) |> 
  pivot_wider(names_from = col_id, values_from = se)

# kable(se.5.equal, format = "latex", booktabs = TRUE, align = "c")
# kable(se.5.weight, format = "latex", booktabs = TRUE, align = "c")

## print out 
table.5.weight <- rbind(alphas.5.weight, se.5.weight)
table.5.weight <- table.5.weight |> arrange(row_id)
kable(table.5.weight, format = "latex", booktabs = TRUE, align = "c")

################################################################################
########################### ADJUST THE TOP INTERVALS ###########################
################################################################################

## 0.5% extreme positive return
alphas.05.equal <- alphas.05 |> 
  filter(weight == 'equal' & day %in% c(1,5,10) & model != 'Excess') |> 
  mutate(
    col_id = paste(type, model, sep = "."),
    row_id = paste0("day", day, ".decile", decile)
  ) |> 
  select(row_id, col_id, intercept) |> 
  pivot_wider(names_from = col_id, values_from = intercept)

se.05.equal <- alphas.05 |> 
  filter(weight == 'equal' & day %in% c(1,5,10) & model != 'Excess') |> 
  mutate(
    se = paste0("(", round(std.error, 3), ")"),
    col_id = paste(type, model, sep = "."),
    row_id = paste0("day", day, ".decile", decile)
  ) |> 
  select(row_id, col_id, se) |> 
  pivot_wider(names_from = col_id, values_from = se)


kable(alphas.05.equal, format = "latex", booktabs = TRUE, align = "c")
kable(se.05.equal, format = "latex", booktabs = TRUE, align = "c")


alphas.05.weight <- alphas.05 |> 
  filter(weight == 'holding_changes' & day %in% c(1,5,10) & model != 'Excess') |> 
  mutate(
    col_id = paste(type, model, sep = "."),
    row_id = paste0("day", day, ".decile", decile)
  ) |> 
  select(row_id, col_id, intercept) |> 
  pivot_wider(names_from = col_id, values_from = intercept)

se.05.weight <- alphas.05 |> 
  filter(weight == 'holding_changes' & day %in% c(1,5,10) & model != 'Excess') |> 
  mutate(
    se = paste0("(", round(std.error, 3), ")"),
    col_id = paste(type, model, sep = "."),
    row_id = paste0("day", day, ".decile", decile)
  ) |> 
  select(row_id, col_id, se) |> 
  pivot_wider(names_from = col_id, values_from = se)


kable(alphas.05.weight, format = "latex", booktabs = TRUE, align = "c")
kable(se.05.weight, format = "latex", booktabs = TRUE, align = "c")

################################################################################
## STEP 4b: 10% extreme interval and drop top 0.5%
rh_taq_10 <- rh_taq_crsp[(sa_holding_change %between% quantile(sa_holding_change, c(.005, .995), na.rm = TRUE)) | 
                           (sa_diff_trd %between% quantile(sa_diff_trd, c(.005, .995), na.rm = TRUE))]

summary(rh_taq_10[, .(sa_holding_change, sa_diff_trd) ])

# create qtile variable for target measure
rh_taq_10 <- get_qtile_data(rh_taq_10, c("sa_holding_change", "sa_diff_trd"), extreme = .10)
# run the whole anlysis
alphas.10 <- get_ret_alpha(rh_taq_10, varlists = c("sa_holding_change", "sa_diff_trd"))

alphas.10.equal <- alphas.10 |> 
  filter(weight == 'equal' & day %in% c(1,5,10) & model != 'Excess') |> 
  mutate(
    col_id = paste(type, model, sep = "."),
    row_id = paste0("day", day, ".decile", decile)
  ) |> 
  select(row_id, col_id, intercept) |> 
  pivot_wider(names_from = col_id, values_from = intercept)

se.10.equal <- alphas.10 |> 
  filter(weight == 'equal' & day %in% c(1,5,10) & model != 'Excess') |> 
  mutate(
    se = paste0("(", round(std.error, 3), ")"),
    col_id = paste(type, model, sep = "."),
    row_id = paste0("day", day, ".decile", decile)
  ) |> 
  select(row_id, col_id, se) |> 
  pivot_wider(names_from = col_id, values_from = se)

kable(alphas.10.equal, format = "latex", booktabs = TRUE, align = "c")
kable(se.10.equal, format = "latex", booktabs = TRUE, align = "c")

table.10.equal <- rbind(alphas.10.equal, se.10.equal)
table.10.equal <- table.10.equal |> arrange(row_id)
kable(table.10.equal, format = "latex", booktabs = TRUE, align = "c")


alphas.10.weight <- alphas.10 |> 
  filter(weight == 'holding_changes' & day %in% c(1,5,10) & model != 'Excess') |> 
  mutate(
    col_id = paste(type, model, sep = "."),
    row_id = paste0("day", day, ".decile", decile)
  ) |> 
  select(row_id, col_id, intercept) |> 
  pivot_wider(names_from = col_id, values_from = intercept)

se.10.weight <- alphas.10 |> 
  filter(weight == 'holding_changes' & day %in% c(1,5,10) & model != 'Excess') |> 
  mutate(
    se = paste0("(", round(std.error, 3), ")"),
    col_id = paste(type, model, sep = "."),
    row_id = paste0("day", day, ".decile", decile)
  ) |> 
  select(row_id, col_id, se) |> 
  pivot_wider(names_from = col_id, values_from = se)

kable(alphas.10.weight, format = "latex", booktabs = TRUE, align = "c")
kable(se.10.weight, format = "latex", booktabs = TRUE, align = "c")

alphas.10.weight <- rbind(alphas.10.weight, se.10.weight)
alphas.10.weight <- alphas.10.weight |> arrange(row_id)
kable(alphas.10.weight, format = "latex", booktabs = TRUE, align = "c")

