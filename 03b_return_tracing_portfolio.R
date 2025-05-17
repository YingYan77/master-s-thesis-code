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
########################## FUNCTIONS for further use ###########################
################################################################################

get_qtile_data <- function(dt, varlists, extreme){
  dt_clean <- dt[complete.cases(dt[, ..varlists])]
  
  dt_clean[, (paste0("qtile_", varlists)) := 
             lapply(.SD, function(var) findInterval(var, 
                                                    quantile(var, prob = c(extreme,.25,.5,.75,(1-extreme)), na.rm = TRUE), 
                                                    left.open = TRUE) + 1),
           by = .(date),
           .SDcols = varlists]
  
  return(dt_clean)
}

get_buyhold_ret_portfolio <- function(dt, investor, weight_var, weight = c("equal", "holding_changes"), 
                                  portfolio = c("1", "10", "101")){
  # # get weight variables 
  # weight_var <- investor
  
  # get the top and bottom deciles 
  bottom5 <- dt[get(paste0("qtile_", investor))==1, ]
  top5 <- dt[get(paste0("qtile_", investor))==6, ]
  extreme5 <- dt[get(paste0("qtile_", investor)) %in% c(1,6), ]
  
  # calculate the daily weighted average return 
  if (weight == "equal") {
    bottom_portfolio <- bottom5[, .(retp_lead_0 = mean(ret_lead_0, na.rm = TRUE)*100, 
                                    retp_lead_1 = mean(ret_lead_1, na.rm = TRUE)*100, 
                                    retp_lead_5 = mean((cumret_5d - 1)/5, na.rm = TRUE)*100, 
                                    retp_lead_10 = mean((cumret_10d - 1)/10, na.rm = TRUE)*100), 
                                by = .(date)]
    
    top_portfolio <- top5[, .(retp_lead_0 = mean(ret_lead_0, na.rm = TRUE)*100, 
                              retp_lead_1 = mean(ret_lead_1, na.rm = TRUE)*100, 
                              retp_lead_5 = mean((cumret_5d - 1)/5, na.rm = TRUE)*100, 
                              retp_lead_10 = mean((cumret_10d - 1)/10, na.rm = TRUE)*100), 
                          by = .(date)]
    
    
    extreme_portfolio <- extreme5[, .(retp_lead_0 = mean(ret_lead_0, na.rm = TRUE)*100, 
                                      retp_lead_1 = mean(ret_lead_1, na.rm = TRUE)*100, 
                                      retp_lead_5 = mean((cumret_5d - 1)/5, na.rm = TRUE)*100, 
                                      retp_lead_10 = mean((cumret_10d - 1)/10, na.rm = TRUE)*100), 
                                  by = .(date)]
    
  } else if (weight == "holding_changes") {
    bottom5[, change_weight := abs(get(weight_var)) / sum(abs(get(weight_var))), 
            by = .(date)]
    top5[, change_weight := abs(get(weight_var)) / sum(abs(get(weight_var))), 
         by = .(date)]
    extreme5[, change_weight := abs(get(weight_var)) / sum(abs(get(weight_var))), 
            by = .(date)]
    
    bottom_portfolio <- bottom5[, .(retp_lead_0 = sum(change_weight*ret_lead_0, na.rm = TRUE)*100, 
                                    retp_lead_1 = sum(change_weight*ret_lead_1, na.rm = TRUE)*100, 
                                    retp_lead_5 = sum(change_weight*(cumret_5d - 1)/5, na.rm = TRUE)*100, 
                                    retp_lead_10 = sum(change_weight*(cumret_10d - 1)/10, na.rm = TRUE)*100), 
                                by = .(date)]
    
    top_portfolio <- top5[, .(retp_lead_0 = sum(change_weight*ret_lead_0, na.rm = TRUE)*100, 
                              retp_lead_1 = sum(change_weight*ret_lead_1, na.rm = TRUE)*100, 
                              retp_lead_5 = sum(change_weight*(cumret_5d - 1)/5, na.rm = TRUE)*100, 
                              retp_lead_10 = sum(change_weight*(cumret_10d - 1)/10, na.rm = TRUE)*100), 
                          by = .(date)]
    
    
    extreme_portfolio <- extreme5[, .(retp_lead_0 = sum(change_weight*ret_lead_0, na.rm = TRUE)*100, 
                                      retp_lead_1 = sum(change_weight*ret_lead_1, na.rm = TRUE)*100, 
                                      retp_lead_5 = sum(change_weight*(cumret_5d - 1)/5, na.rm = TRUE)*100, 
                                      retp_lead_10 = sum(change_weight*(cumret_10d - 1)/10, na.rm = TRUE)*100), 
                                  by = .(date)]
  } 
  
  if (portfolio == "10") {
    return(top_portfolio)
  } else if (portfolio == "1") {
    return(bottom_portfolio)
  } else if (portfolio == "101") {
    return(extreme_portfolio)
  }
}

get_alpha <- function(dt, ret_choice, hperiod, investor, dcile, weight) {
  
  Excess <- lm(get(ret_choice) ~ 1, data = dt)
  NW_Excess <- NeweyWest(Excess, lag = 1, prewhite = FALSE, adjust = TRUE)
  se_Excess <- sqrt(diag(NW_Excess))[1]
  # se_Excess <- sqrt(diag(vcovHAC(Excess, prewhite = FALSE)))[1]
  alpha_Excess <- Excess$coefficients[1]
  t_Excess <- alpha_Excess / se_Excess
  
  CAPM <- lm(get(ret_choice) ~ Mkt.RF, data = dt)
  # Newey-West adjusted t-test
  NW_CAPM <- NeweyWest(CAPM, lag = 1, prewhite = FALSE, adjust = TRUE)
  se_CAPM <- sqrt(diag(NW_CAPM))[1]
  # se_CAPM <- sqrt(diag(vcovHAC(CAPM, prewhite = FALSE)))[1]
  alpha_CAPM <- CAPM$coefficients[1]
  t_CAPM <- alpha_CAPM / se_CAPM 
  
  FF3 <- lm(get(ret_choice) ~ Mkt.RF + SMB + HML, data = dt)
  NW_FF3 <- NeweyWest(FF3, lag = 1, prewhite = FALSE, adjust = TRUE)
  se_FF3 <- sqrt(diag(NW_FF3))[1]
  # se_FF3 <- sqrt(diag(vcovHAC(FF3, prewhite = FALSE)))[1]
  alpha_FF3 <- FF3$coefficients[1]
  t_FF3 <- alpha_FF3 / se_FF3 
  
  FF6 <- lm(get(ret_choice) ~ Mkt.RF + SMB + HML + RMW + CMA + Mom, data = dt)
  NW_FF6 <- NeweyWest(FF6, lag = 1, prewhite = FALSE, adjust = TRUE)
  se_FF6 <- sqrt(diag(NW_FF6))[1]
  # se_FF6 <- sqrt(diag(vcovHAC(FF6, prewhite = FALSE)))[1]
  alpha_FF6 <- FF6$coefficients[1]
  t_FF6 <- alpha_FF6 / se_FF6 
  
  output <- rbind(c(alpha_Excess, se_Excess, t_Excess), 
                  c(alpha_CAPM, se_CAPM, t_CAPM), 
                  c(alpha_FF3, se_FF3, t_FF3), 
                  c(alpha_FF6, se_FF6, t_FF6))
  colnames(output) <- c('estimate', 'std.error', 't.value')
  output <- as.data.frame(output)
  output <- output |> 
    mutate(intercept = case_when(abs(t.value) >= 2.576 ~ paste0(round(estimate,3), "***"), 
                                 abs(t.value) >= 1.96 ~ paste0(round(estimate,3), "**"), 
                                 abs(t.value) >= 1.645 ~ paste0(round(estimate,3), "*"), 
                                 abs(t.value) < 1.645 ~ paste0(round(estimate,3))
    ))
  
  output$model <- c("Excess", "CAPM", "FF3", "FF6")
  output$type <- investor
  output$day <- hperiod
  output$decile <- dcile
  output$weight <- weight
  
  return(output)
}

get_diff_alpha <- function(dt1, dt10, ret_choice, 
                           model = c("Excess", "CAPM", "FF3", "FF6"), 
                           hperiod, investor, weight) {
  n <- max(nrow(dt1), nrow(dt10))
  if (model == "Excess") {
    model_1 <- lm(get(ret_choice) ~ 1, data = dt1)
    model_10 <- lm(get(ret_choice) ~ 1, data = dt10)
  } else if (model == "CAPM") {
    model_1 <- lm(get(ret_choice) ~ Mkt.RF, data = dt1)
    model_10 <- lm(get(ret_choice) ~ Mkt.RF, data = dt10)
  } else if (model == "FF3") {
    model_1 <- lm(get(ret_choice) ~ Mkt.RF + SMB + HML, data = dt1)
    model_10 <- lm(get(ret_choice) ~ Mkt.RF + SMB + HML, data = dt10)
  } else if (model == "FF6") {
    model_1 <- lm(get(ret_choice) ~ Mkt.RF + SMB + HML + CMA + Mom, data = dt1)
    model_10 <- lm(get(ret_choice) ~ Mkt.RF + SMB + HML + CMA + Mom, data = dt10)
  }
  
  vcov_1 <- NeweyWest(model_1, lag = 1, prewhite = FALSE, adjust = TRUE)
  se_1 <- sqrt(diag(vcov_1))[1]
  alpha_1 <- coef(model_1)[1]
  residuals_1 <- residuals(model_1)
  
  vcov_10 <- NeweyWest(model_10, lag = 1, prewhite = FALSE, adjust = TRUE)
  se_10 <- sqrt(diag(vcov_10))[1]
  alpha_10 <- coef(model_10)[1]
  residuals_10 <- residuals(model_10)
  
  cov_residuals <- sum(residuals_1 * residuals_10) / n - 1 
  # Estimate the covariance of the intercepts (using similar scaling as in Newey-West)
  # Scale by the inverse of the sum of the squared weights (1/n) in each model
  X1 <- model.matrix(model_1)
  X10 <- model.matrix(model_10)
  cov_alpha1_alpha2 <- cov_residuals * solve(t(X1) %*% X1)[1, 1] * solve(t(X10) %*% X10)[1, 1]
  
  var_diff_alpha <- vcov_1[1,1] + vcov_10[1,1] - 2*cov_alpha1_alpha2
  se_diff <- sqrt(var_diff_alpha)
  alpha_diff <- alpha_10 - alpha_1
  t_diff <- alpha_diff/se_diff
  
  output <- as.data.frame(c(alpha_diff, se_diff, t_diff))
  output <- as.data.frame(t(output))
  
  colnames(output) <- c('estimate', 'std.error', 't.value')
  rownames(output) <- NULL
  output <- output |> 
    mutate(intercept = case_when(abs(t.value) >= 2.576 ~ paste0(round(estimate,3), "***"), 
                                 abs(t.value) >= 1.96 ~ paste0(round(estimate,3), "**"), 
                                 abs(t.value) >= 1.645 ~ paste0(round(estimate,3), "*"), 
                                 abs(t.value) < 1.645 ~ paste0(round(estimate,3))
    ))
  
  output$model <- model
  output$type <- investor
  output$day <- hperiod
  output$decile <- "longshort"
  output$weight <- weight
  
  return(output)
}

# get_alpha(rh_portfolio, "Retp.Rf", hperiod = 1, "holding_ratio", "10", weight= 'equal')

get_ret_alpha_extreme <- function(dt, varlists, retvar) {
  alphas <- data.frame()
  
  for (i in varlists) {
    for (j in c(0, 1, 5, 10)) {
      for (l in c("equal", "holding_changes")) {
        for (k in c("10", "1", "101")) {
          ret_portfolio <- get_buyhold_ret_portfolio(dt, investor = retvar, weight_var = i, weight = l, portfolio = k)
          ret_portfolio <- merge(ret_portfolio, ff6, by.x = c("date"), by.y = c("Date"))
          ret_portfolio[, Retp.Rf := get(paste0("retp_lead_",j)) - RF] 
          result <- get_alpha(ret_portfolio, "Retp.Rf", hperiod = j, investor = i, dcile = k, weight = l)
          
          alphas <- rbind(alphas, result)
        }
        for (n in c("Excess", "CAPM", "FF3", "FF6")) {
          ret_portfolio_1 <- get_buyhold_ret_portfolio(dt, investor = retvar, weight_var = i, weight = l, portfolio = 1)
          ret_portfolio_1 <- merge(ret_portfolio_1, ff6, by.x = c("date"), by.y = c("Date"))
          ret_portfolio_1[, Retp.Rf := get(paste0("retp_lead_",j)) - RF] 
          
          ret_portfolio_10 <- get_buyhold_ret_portfolio(dt, investor = retvar, weight_var = i, weight = l, portfolio = 10)
          ret_portfolio_10 <- merge(ret_portfolio_10, ff6, by.x = c("date"), by.y = c("Date"))
          ret_portfolio_10[, Retp.Rf := get(paste0("retp_lead_",j)) - RF] 
          
          result <- get_diff_alpha(ret_portfolio_1, ret_portfolio_10, 
                                   model = n, "Retp.Rf", hperiod = j, investor = i, weight = l)
          
          alphas <- rbind(alphas, result)
        }
      }
    }
  }
  return(alphas)
  
}

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
crsp[,  paste0("ret_lag_", 0:n_lags)  := shift(ret, n = 0:n_lags, type = "lag"), by = .(tsymbol)]
## calculate the cumulative return from t+1 to t+5
crsp[, cumret_5d := frollapply(shift(1 + ret, 1, type = "lead"), 5, prod, fill = NA, align = "left"), by = .(tsymbol)]
## calculate the cumulative return from t+1 to t+10
crsp[, cumret_10d := frollapply(shift(1 + ret, 1, type = "lead"), 10, prod, fill = NA, align = "left"), by = .(tsymbol)]

# Apply rolling product for each lag and store it as a new column
for (k in 1:5) {
  crsp[, paste0("cumret_lag_", k) := 
         frollapply(shift(1 + ret, 1, type = "lag"), k, prod, fill = NA),
       by = .(tsymbol)]
  
  crsp[, paste0("cum_mktret_lag_", k) := 
         frollapply(shift(1 + mkt_ret, 1, type = "lag"), k, prod, fill = NA),
       by = .(tsymbol)]
  
  crsp[, paste0("car_lag_", k) := 
         get(paste0("cumret_lag_", k)) - get(paste0("cum_mktret_lag_", k)),
       by = .(tsymbol)]
}
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

## the main analysis for return
rh_taq_crsp <- get_qtile_data(rh_taq_crsp, c("ret", "sa_diff_trd", "sa_holding_change"), extreme = .05)

# run the whole anlysis
alphas.5 <- get_ret_alpha_extreme(rh_taq_crsp, 
                                  varlists = c("sa_holding_change", "sa_diff_trd"), 
                                  retvar = "ret")

## alternative analysis
rh_taq_crsp <- get_qtile_data(rh_taq_crsp, c("ret_lag_2", "sa_diff_trd", "sa_holding_change"), extreme = .05)

# run the whole anlysis
alphas.5 <- get_ret_alpha_extreme(rh_taq_crsp, 
                                  varlists = c("sa_holding_change", "sa_diff_trd"), 
                                  retvar = "ret_lag_2")


# additional analysos for CAR
rh_taq_crsp <- get_qtile_data(rh_taq_crsp, c("car_lag_5", "sa_diff_trd", "sa_holding_change"), extreme = .05)

# run the whole anlysis
alphas.5 <- get_ret_alpha_extreme(rh_taq_crsp,
                                  varlists = c("sa_holding_change", "sa_diff_trd"),
                                  retvar = "car_lag_5")

# create qtile variable for target measure
rh_taq_crsp <- get_qtile_data(rh_taq_crsp, c("car_lag_1", "sa_diff_trd", "sa_holding_change"), extreme = .05)

# run the whole anlysis
alphas.5 <- get_ret_alpha_extreme(rh_taq_crsp,
                                  varlists = c("sa_holding_change", "sa_diff_trd"),
                                  retvar = "car_lag_1")



################################################################################
########################## ADJUST FORMATTING TO PRINT ##########################
################################################################################

alphas.5.equal <- alphas.5 |> 
  filter(weight == 'equal' & day!=0 & model != 'Excess') |> 
  mutate(
    col_id = paste(type, model, sep = "."),
    row_id = paste0("day", day, ".decile", decile)
  ) |> 
  select(row_id, col_id, intercept) |> 
  pivot_wider(names_from = col_id, values_from = intercept)

alphas.5.weight <- alphas.5 |> 
  filter(weight == 'holding_changes' & day!=0 & model != 'Excess') |> 
  mutate(
    col_id = paste(type, model, sep = "."),
    row_id = paste0("day", day, ".decile", decile)
  ) |> 
  select(row_id, col_id, intercept) |> 
  pivot_wider(names_from = col_id, values_from = intercept)

# kable(alphas.5.equal, format = "latex", booktabs = TRUE, align = "c")
# kable(alphas.5.weight, format = "latex", booktabs = TRUE, align = "c")

## formatting the standard error
se.5.equal <- alphas.5 |> 
  filter(weight == 'equal' & day!=0 & model != 'Excess') |> 
  mutate(
    se = paste0("(", round(std.error, 3), ")"),
    col_id = paste(type, model, sep = "."),
    row_id = paste0("day", day, ".decile", decile)
  ) |> 
  select(row_id, col_id, se) |> 
  pivot_wider(names_from = col_id, values_from = se)

se.5.weight <- alphas.5 |> 
  filter(weight == 'holding_changes' & day!=0 & model != 'Excess') |> 
  mutate(
    se = paste0("(", round(std.error, 3), ")"),
    col_id = paste(type, model, sep = "."),
    row_id = paste0("day", day, ".decile", decile)
  ) |> 
  select(row_id, col_id, se) |> 
  pivot_wider(names_from = col_id, values_from = se)

# kable(se.5.equal, format = "latex", booktabs = TRUE, align = "c")
# kable(se.5.weight, format = "latex", booktabs = TRUE, align = "c")


table.5.equal <- rbind(alphas.5.equal, se.5.equal)
table.5.equal <- table.5.equal |> arrange(row_id)
kable(table.5.equal, format = "latex", booktabs = TRUE, align = "c")

table.5.weight <- rbind(alphas.5.weight, se.5.weight)
table.5.weight <- table.5.weight |> arrange(row_id)
kable(table.5.weight, format = "latex", booktabs = TRUE, align = "c")


