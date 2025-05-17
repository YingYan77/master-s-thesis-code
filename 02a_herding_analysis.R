################################################################################
## author:    Ying Yan
## file name: 02_Herding_Analysis.R
## summary:   this file aims to study generate the descriptive analysis table 
##            used in the thesis. 
## date:      24 April 2025
################################################################################

## empty memory (!)
rm(list=ls())

library(data.table)
library(sandwich)
library(lmtest)
library(car)
library(fixest)

library(fastDummies)
library(zoo)
library(dplyr)
library(stringr)

## graph's libraries
library(ggplot2)
library(tidyr)
library(patchwork)

## latex
library(knitr)
library(kableExtra)
library(stargazer)

## read the function library
source("00_functions.R")

## read the taq daily aggregated dataset, cleaned from 00_TAQ_RH_CRSP_Merged.R
rh_taq <- fread("data/rh_taq.csv")
crsp <- fread("data/crsp_volatility_compustat_merged.csv")

## merge Robinhood, TAQ and CRSP
rh_taq_crsp <- merge(rh_taq, crsp, by = c("date", "tsymbol"))

################################################################################
## STEP 1: VARIABLE PREPRATON

## create the lags of log_holding_change for 5 days
n_lags = 5
lag_vars <- c("sa_holding_change", "sa_diff_trd", "ret", "volatility", "mkt_ret", "mkt_std_dev")

for (var in lag_vars) {
  rh_taq_crsp[, paste0(var, "_lag_", 0:n_lags) := shift(get(var), n = 0:n_lags, type = "lag"), by = .(tsymbol)]
}

################################################################################
## STEP 2.A: REGRESSION MODEL FOR ROBINHOOD

rh_taq_5 <- rh_taq_crsp
y_vars <- c("sa_holding_change", "sa_diff_trd")
for (var in y_vars){
  rh_taq_5 <- get_lags_qtile_data(rh_taq_5, var, extreme = .05)
}

## run regression from lag 0 to lag 5 data
rh_betas <- data.frame()
rh_beta.table <- data.frame()
rh_waldtest <- data.frame()

for (l in (1:5)) {
  output <- get_daily_qtile_estimates(rh_taq_5, y = "sa_holding_change", var = "sa_holding_change", control = "sa_diff_trd", nlag = l)
  betas_lag <- output[[1]][1:5,c("coefficient","estimate", "Lag")]|> 
    mutate(Interval = str_extract(coefficient, "\\d+$")) |> 
    select(estimate, Lag, Interval)
  beta.table_lag <- output[[1]][1:5,] |> 
    mutate(Interval = str_extract(coefficient, "\\d+$"), 
           sd = paste0("(",round(std.error,3), ")"))
  rh_betas <- rbind(rh_betas, betas_lag)
  rh_beta.table <- rbind(rh_beta.table, beta.table_lag)
  rh_waldtest <- rbind(rh_waldtest, output[[2]])
}

rh_beta.table.arrange <- rh_beta.table[order(rh_beta.table$Interval),]


# GET THE WALT TEST RESULT
rh_waldtest.arrange <- rh_waldtest[order(rh_waldtest$interval),] |> 
  mutate(f_stat = paste0("(",round(F_stat,3), ")"))

kable(rh_waldtest.arrange$difference, format = "latex", booktabs = TRUE, align = "c")
kable(rh_waldtest.arrange$f_stat, format = "latex", booktabs = TRUE, align = "c")


################################################################################
## STEP 2.B: REGRESSION MODEL FOR TAQ

## run regression from lag 0 to lag 5 data
taq_betas <- data.frame()
taq_beta.table <- data.frame()
taq_waldtest <- data.frame()

for (l in (1:5)) {
  output <- get_daily_qtile_estimates(rh_taq_5, y = "sa_diff_trd", var = "sa_diff_trd", control = "sa_holding_change", nlag = l)
  betas_lag <- output[[1]][1:5,c("coefficient","estimate", "Lag")]|> 
    mutate(Interval = str_extract(coefficient, "\\d+$")) |> 
    select(estimate, Lag, Interval)
  # betas_lag$partial_effect[2:6] <- betas_lag$estimate[2:6] + betas_lag$estimate[1]
  # betas_lag$partial_effect[1] <- betas_lag$estimate[1]
  
  beta.table_lag <- output[[1]][1:5,] |> 
    mutate(Interval = str_extract(coefficient, "\\d+$"), 
           sd = paste0("(",round(std.error,3), ")"))
  
  taq_betas <- rbind(taq_betas, betas_lag)
  taq_beta.table <- rbind(taq_beta.table, beta.table_lag)
  taq_waldtest <- rbind(taq_waldtest, output[[2]])
}

taq_beta.table.arrange <- taq_beta.table[order(taq_beta.table$Interval),]

## GET THE WALD TEST REUSLT
taq_waldtest.arrange <- taq_waldtest[order(taq_waldtest$interval),] |> 
  mutate(f_stat = paste0("(",round(F_stat,3), ")"))

kable(taq_waldtest.arrange$difference, format = "latex", booktabs = TRUE, align = "c")
kable(taq_waldtest.arrange$f_stat, format = "latex", booktabs = TRUE, align = "c")


################################################################################
## STEP 2.C: PRINT THE RESULT INTO TABLE 

rh_beta.table.arrange$investor <- "rh"
taq_beta.table.arrange$investor <- "taq"
herding.beta.table <- rbind(rh_beta.table.arrange, taq_beta.table.arrange)

kable(reshape_betas(herding.beta.table), format = "latex", booktabs = TRUE, align = "c")
# R-SQR NUMBER
kable(taq_beta.table.arrange$adj.r.sqrd, format = "latex", booktabs = TRUE, align = "c")
kable(rh_beta.table.arrange$adj.r.sqrd, format = "latex", booktabs = TRUE, align = "c")

################################################################################
## STEP 3: VISUALIZE REGRESSOR FOR RH AND TAQ

rh_lag_plot <- plot_daily_return_beta(rh_betas, xaxis = "Lag", clas = "Interval", ylim_d = -1, ylim_u = 2)
taq_lag_plot <- plot_daily_return_beta(taq_betas, xaxis = "Lag", clas = "Interval", ylim_d = -1, ylim_u = 2)

herding_plot <- rh_lag_plot + taq_lag_plot + 
  plot_layout(guides = "collect", axes = "collect") 

ggsave("plots/rh_herding.pdf", rh_lag_plot, 
       width = 8,
       height = 4.5,
       # units = "in", 
       limitsize = FALSE)

ggsave("plots/taq_herding.pdf", taq_lag_plot, 
       width = 8,
       height = 4.5,
       # units = "in", 
       limitsize = FALSE)


################################################################################
## STEP 4.A: ADJUST FOR QUINTILES 0.5% extreme interval

rh_taq_05 <- rh_taq_crsp

y_vars <- c("sa_holding_change", "sa_diff_trd")
for (var in y_vars){
  rh_taq_05 <- get_lags_qtile_data(rh_taq_05, var, extreme = .005)
}

## run Robinhood regression from lag 0 to lag 5 data
rh_betas <- data.frame()
rh_beta.table <- data.frame()
rh_waldtest <- data.frame()

for (l in (1:5)) {
  output <- get_daily_qtile_estimates(rh_taq_05, y = "sa_holding_change", var = "sa_holding_change", control = "sa_diff_trd", nlag = l)
  betas_lag <- output[[1]][1:5,c("coefficient","estimate", "Lag")]|> 
    mutate(Interval = str_extract(coefficient, "\\d+$")) |> 
    select(estimate, Lag, Interval)
  beta.table_lag <- output[[1]][1:5,] |> 
    mutate(Interval = str_extract(coefficient, "\\d+$"), 
           sd = paste0("(",round(std.error,3), ")"))
  rh_betas <- rbind(rh_betas, betas_lag)
  rh_beta.table <- rbind(rh_beta.table, beta.table_lag)
  rh_waldtest <- rbind(rh_waldtest, output[[2]])
}

rh_beta.table.arrange <- rh_beta.table[order(rh_beta.table$Interval),]



## run TAQ regression from lag 1 to lag 5 data
taq_betas <- data.frame()
taq_beta.table <- data.frame()
taq_waldtest <- data.frame()

for (l in (1:5)) {
  output <- get_daily_qtile_estimates(rh_taq_05, y = "sa_diff_trd", var = "sa_diff_trd", control = "sa_holding_change", nlag = l)
  betas_lag <- output[[1]][1:5,c("coefficient","estimate", "Lag")]|> 
    mutate(Interval = str_extract(coefficient, "\\d+$")) |> 
    select(estimate, Lag, Interval)

  beta.table_lag <- output[[1]][1:5,] |> 
    mutate(Interval = str_extract(coefficient, "\\d+$"), 
           sd = paste0("(",round(std.error,3), ")"))
  
  taq_betas <- rbind(taq_betas, betas_lag)
  taq_beta.table <- rbind(taq_beta.table, beta.table_lag)
  taq_waldtest <- rbind(taq_waldtest, output[[2]])
}

taq_beta.table.arrange <- taq_beta.table[order(taq_beta.table$Interval),]

## PRINT OUT THE RESULT
rh_beta.table.arrange$investor <- "rh"
taq_beta.table.arrange$investor <- "taq"
herding.beta.table <- rbind(rh_beta.table.arrange, taq_beta.table.arrange)

kable(reshape_betas(herding.beta.table), format = "latex", booktabs = TRUE, align = "c")
# R-SQR NUMBER
kable(taq_beta.table.arrange$adj.r.sqrd, format = "latex", booktabs = TRUE, align = "c")
kable(rh_beta.table.arrange$adj.r.sqrd, format = "latex", booktabs = TRUE, align = "c")

################################################################################
## STEP 4b: ADJUST FOR DECILES 10% extreme interval and drop top 0.5%
rh_taq_10 <- rh_taq_crsp[(sa_holding_change %between% quantile(sa_holding_change, c(.005, .995), na.rm = TRUE)) | 
                           (sa_diff_trd %between% quantile(sa_diff_trd, c(.005, .995), na.rm = TRUE))]

summary(rh_taq_10[, .(sa_holding_change, sa_diff_trd) ])

## create the lags of log_holding_change for 5 days
n_lags = 5
lag_vars <- c("sa_holding_change", "sa_diff_trd", "ret", "volatility", "mkt_ret", "mkt_std_dev")

for (var in lag_vars) {
  rh_taq_10[, paste0(var, "_lag_", 0:n_lags) := shift(get(var), n = 0:n_lags, type = "lag"), by = .(tsymbol)]
}

y_vars <- c("sa_holding_change", "sa_diff_trd")
for (var in y_vars){
  rh_taq_10 <- get_lags_qtile_data(rh_taq_10, var, extreme = .10)
}

## run Robinhood regression from lag 0 to lag 5 data
rh_betas <- data.frame()
rh_beta.table <- data.frame()
rh_waldtest <- data.frame()

for (l in (1:5)) {
  output <- get_daily_dcile_estimates(rh_taq_10, y = "sa_holding_change", var = "sa_holding_change", control = "sa_diff_trd", nlag = l)
  betas_lag <- output[[1]][1:9,c("coefficient", "estimate", "Lag")]|> 
    mutate(Interval = str_extract(coefficient, "\\d+$")) |> 
    select(estimate, Lag, Interval)
  beta.table_lag <- output[[1]][1:9,] |> 
    mutate(Interval = str_extract(coefficient, "\\d+$"), sd = paste0("(",round(std.error,3), ")"))
  rh_betas <- rbind(rh_betas, betas_lag)
  rh_beta.table <- rbind(rh_beta.table, beta.table_lag)
  rh_waldtest <- rbind(rh_waldtest, output[[2]])
}

rh_beta.table.arrange <- rh_beta.table[order(rh_beta.table$Interval),]



## run TAQ regression from lag 1 to lag 5 data
taq_betas <- data.frame()
taq_beta.table <- data.frame()
taq_waldtest <- data.frame()

for (l in (1:5)) {
  output <- get_daily_dcile_estimates(rh_taq_10, y = "sa_diff_trd", var = "sa_diff_trd", control = "sa_holding_change", nlag = l)
  betas_lag <- output[[1]][1:9,c("coefficient", "estimate", "Lag")]|> 
    mutate(Interval = str_extract(coefficient, "\\d+$")) |> 
    select(estimate, Lag, Interval)
  
  beta.table_lag <- output[[1]][1:9,] |> 
    mutate(Interval = str_extract(coefficient, "\\d+$"), sd = paste0("(",round(std.error,3), ")"))
  
  taq_betas <- rbind(taq_betas, betas_lag)
  taq_beta.table <- rbind(taq_beta.table, beta.table_lag)
  taq_waldtest <- rbind(taq_waldtest, output[[2]])
}

taq_beta.table.arrange <- taq_beta.table[order(taq_beta.table$Interval),]

## PRINT OUT THE RESULT
rh_beta.table.arrange$investor <- "rh"
taq_beta.table.arrange$investor <- "taq"
herding.beta.table <- rbind(rh_beta.table.arrange, taq_beta.table.arrange)

kable(reshape_betas(herding.beta.table), format = "latex", booktabs = TRUE, align = "c")
kable(rh_beta.table.arrange$adj.r.sqrd, format = "latex", booktabs = TRUE, align = "c")
kable(taq_beta.table.arrange$adj.r.sqrd, format = "latex", booktabs = TRUE, align = "c")