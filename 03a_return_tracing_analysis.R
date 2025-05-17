################################################################################
## author:    Ying Yan
## file name: 04a_Representativeness_Analysis.R
## summary:   this file aims to study generate the descriptive analysis table 
##            used in the thesis. 
## date:      12 May 2025
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

## compute the cumulative return and CAR
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

summary(crsp$cum_mktret_lag_5)
summary(crsp$cumret_lag_5)
summary(crsp$car_lag_5)

## merge Robinhood, TAQ and CRSP
rh_taq_crsp <- merge(rh_taq, crsp, by = c("date", "tsymbol"))

## create the lags of log_holding_change for 5 days
n_lags = 5
lag_vars <- c("sa_holding_change", "sa_diff_trd", "ret", "volatility", "mkt_ret", "mkt_std_dev")

for (var in lag_vars) {
  rh_taq_crsp[, paste0(var, "_lag_", 0:n_lags) := shift(get(var), n = 0:n_lags, type = "lag"), by = .(tsymbol)]
}

################################################################################
## STEP 1: DEFINE FUNCTIONS
get_lags_qtile_data <- function(dt, varname, extreme){
  var_lags <- grep(paste0("^", varname, "_lag_"), names(dt), value = TRUE)
  dt_clean <- dt[complete.cases(dt[, ..var_lags])]
  
  dt_clean[, (paste0("qtile_", var_lags)) := 
             lapply(.SD, function(var) findInterval(var, 
                                                    quantile(var, prob = c(extreme,.25,.5,.75,(1-extreme)), na.rm = TRUE), 
                                                    left.open = TRUE) + 1),
           by = .(date),
           .SDcols = var_lags]
  
  return(dt_clean)
}

get_daily_qtile_estimates <- function(dt, y, var, control, nlag) {
  
  # get dataframe
  lag.data <- dt[!is.na(get(paste0("qtile_", var, "_lag_", nlag)))]
  lag.data <- dummy_cols(lag.data, select_columns = paste0("qtile_", var, "_lag_", nlag))
  
  # dt <- get(paste0("robintrack_lag_", 1))
  
  # get lag quintle variables
  qitle.vars <- paste0("qtile_", var, "_lag_", nlag, "_", c(1,2,3,5,6), collapse = " + ")
  
  # add control variables if it exist
  mkt.ret.vars <- paste0("mkt_ret_lag_", 1:5, collapse = " + ")
  mkt.sd.vars <- paste0("mkt_std_dev_lag_", 1:5, collapse = " + ")
  
  # get return and volatility variables
  if (var == "ret" | var == "cumret" | var == "car") {
    ret.vars <- paste0("ret_lag_", setdiff(1:5, nlag), collapse = " + ")
    sd.vars <- paste0("volatility_lag_", setdiff(1:5, nlag), collapse = " + ")
    rh.vars <- paste0("sa_holding_change_lag_",1:5, collapse = " + ")
    taq.vars <- paste0("sa_diff_trd_lag_", 1:5, collapse = " + ")
  } else {
    ret.vars <- paste0("ret_lag_", 1:5, collapse = " + ")
    sd.vars <- paste0("volatility_lag_", 1:5, collapse = " + ")
    rh.vars <- paste0("sa_holding_change_lag_", setdiff(1:5, nlag), collapse = " + ")
    taq.vars <- paste0("sa_diff_trd_lag_", setdiff(1:5, nlag), collapse = " + ")
  }
  
  # get regressors ready
  if (missing(control)) {
    # without control variable
    regressors.herding <- paste(qitle.vars, 
                                ret.vars, 
                                sd.vars, 
                                # rh.vars, 
                                # taq.vars,
                                mkt.ret.vars, 
                                mkt.sd.vars, 
                                "mcap", "BtM", "gp", "inv", "momentum", 
                                sep = " + "
    )
    
    regressors.represent <- paste(qitle.vars, 
                                  ret.vars,
                                  sd.vars, 
                                  rh.vars,
                                  taq.vars,
                                  mkt.ret.vars,
                                  mkt.sd.vars, 
                                  "mcap", "BtM", "gp", "inv", "momentum", 
                                  sep = " + "
    )
  } else {
    lag.data <- lag.data[!is.na(get(paste0("qtile_", control, "_lag_", nlag)))]
    lag.data <- dummy_cols(lag.data, select_columns = paste0("qtile_", control, "_lag_", nlag))
    control.vars <- paste0("qtile_", control, "_lag_", nlag, "_", c(1,2,3,5,6), collapse = " + ")
    regressors.herding <- paste(qitle.vars, 
                                control.vars, 
                                ret.vars, 
                                sd.vars, 
                                # rh.vars, 
                                # taq.vars,
                                mkt.ret.vars, 
                                mkt.sd.vars,  
                                "mcap", "BtM", "gp", "inv", "momentum", 
                                sep = " + "
    )
    
    regressors.represent <- paste(qitle.vars, 
                                  control.vars, 
                                  # ret.vars,
                                  sd.vars, 
                                  rh.vars,
                                  taq.vars,
                                  mkt.ret.vars, 
                                  mkt.sd.vars,  
                                  "mcap", "BtM", "gp", "inv", "momentum", 
                                  sep = " + "
    )
  }
  
  if (var == "ret" | var == "cumret" | var == "car") {
    formula <- as.formula(paste(y,"~", regressors.represent, "| date"))
  } else {
    formula <- as.formula(paste(y,"~", regressors.herding, "| date"))
  }
  
  
  # run regression
  # model <- lm(formula, data = lag.data)
  model <- feols(formula, data = lag.data)
  
  # save betas and clustered standard error
  coefs <- model$coefficients
  # cluster_varv <- vcovCL(model, cluster = ~ date + tsymbol)
  cluster_varv <- vcov(model, cluster = ~ tsymbol + date)
  model_sum <- summary(model,  cluster = ~ tsymbol + date)
  # cluster_se <- sqrt(diag(cluster_varv))
  cluster_se <- model_sum$se
  n.obs <- model_sum$nobs
  adj.r.squared <- fitstat(model, "ar2")$ar2
  
  # save into output
  output.table <- data.frame(
    coefficient = names(coefs),
    estimate = coefs,
    std.error = cluster_se,
    t.value = coefs / cluster_se,
    Lag = nlag,
    nobs = n.obs,
    adj.r.sqrd = adj.r.squared
  )
  
  rownames(output.table) <- NULL
  
  output.table <- output.table |>
    mutate(intercept =  case_when(abs(t.value) >= 2.576 ~ paste0(round(estimate,3), "***"),
                                  abs(t.value) >= 1.96 ~ paste0(round(estimate,3), "**"),
                                  abs(t.value) >= 1.645 ~ paste0(round(estimate,3), "*"),
                                  abs(t.value) < 1.645 ~ paste0(round(estimate,3))))
  
  if (var == "ret" | var == "cumret" | var == "car") {
    return(output.table)
  } else {
    ## calculate the difference
    library(car)
    
    # Define the variable pairs to test
    lags <- c(1,2,3,5,6)
    results <- list()
    
    # Loop through each lag to compare the holding_change and diff_trd coefficients
    for (lag in c(1,2,3,5,6)) {
      var1 <- paste0("qtile_", var, "_lag_" , nlag, "_", lag)
      var2 <- paste0("qtile_", control, "_lag_" , nlag, "_", lag)
      hypothesis <- paste(var1, "=", var2)
      
      # Perform Wald test and store the result
      test_result <- linearHypothesis(model, hypothesis, vcov = cluster_varv)
      results[[paste0("lag_", lag)]] <- test_result
    }
    
    
    wald_test_table <- data.frame(
      interval = c(1,2,3,5,6),
      lag = nlag,
      X_variable = paste0("qtile_", var, "_lag_" , nlag, "_", lags),
      ctrl_variable = paste0("qtile_", control, "_lag_" , nlag, "_", lags),
      # F_stat = sapply(results, function(res) res[2, "F"]),
      # p_value = sapply(results, function(res) res[2, "Pr(>F)"])
      F_stat = sapply(results, function(res) res[2, "Chisq"]),
      p_value = sapply(results, function(res) res[2, "Pr(>Chisq)"])
    )
    
    
    wald_test_table <- wald_test_table |>
      mutate(difference = case_when(abs(p_value) <= 0.01 ~ paste0(round(coefs[X_variable] - coefs[ctrl_variable],3), "***"),
                                    abs(p_value) <= 0.05 & abs(p_value) >= 0.01 ~ paste0(round(coefs[X_variable] - coefs[ctrl_variable],3), "**"),
                                    abs(p_value) <= 0.1  & abs(p_value) >= 0.05 ~ paste0(round(coefs[X_variable] - coefs[ctrl_variable],3), "*"),
                                    abs(p_value) > 0.1 ~ paste0(round(coefs[X_variable] - coefs[ctrl_variable],3))))
    
    output <- list(output.table, wald_test_table)
    
    return(output)
  }
  
  
}

plot_daily_return_beta <- function(df, xaxis, clas, ylim_d, ylim_u) {
  
  plot <- ggplot(df, aes(x = get(xaxis), y = estimate, color = factor(get(clas)), group = factor(get(clas)))) +
    geom_line() +
    geom_point() + 
    ylim(ylim_d, ylim_u) + 
    scale_color_discrete(name = clas) +
    labs(
      x = xaxis,
      y = "Beta Estimates"
    ) +
    theme_bw() + 
    scale_x_continuous(breaks = seq(min(df[[xaxis]]), max(df[[xaxis]]), by = 1))
  
  return(plot)
}

################################################################################
## STEP 2.A: REGRESSION MODEL FOR ROBINHOOD

robintrack <- rh_taq_crsp
y_vars <- c("sa_holding_change", "cumret", "car", "ret")
for (name in y_vars){
  robintrack <- get_lags_qtile_data(robintrack, name, extreme = .05)
}

## run regression from lag 0 to lag 5 data
rh_betas <- data.frame()
rh_beta.table <- data.frame()

for (l in (1:5)) {
  output <- get_daily_qtile_estimates(robintrack, y = "sa_holding_change", var = "ret", nlag = l)
  betas_lag <- output[1:5,c("coefficient","estimate", "Lag")]|> 
    mutate(Interval = str_extract(coefficient, "\\d+$")) |> 
    select(estimate, Lag, Interval)
  beta.table_lag <- output[1:5,] |> 
    mutate(Interval = str_extract(coefficient, "\\d+$"), 
           sd = paste0("(",round(std.error,3), ")"))
  rh_betas <- rbind(rh_betas, betas_lag)
  rh_beta.table <- rbind(rh_beta.table, beta.table_lag)
}

summary(rh_betas)

rh_beta.table.arrange <- rh_beta.table[order(rh_beta.table$Interval),]

################################################################################
## STEP 2.B: REGRESSION MODEL FOR TAQ

taq <- rh_taq_crsp
y_vars <- c("sa_diff_trd", "cumret", "car", "ret")
for (name in y_vars){
  taq <- get_lags_qtile_data(taq, name, extreme = .05)
}

## run regression from lag 0 to lag 5 data
taq_betas <- data.frame()
taq_beta.table <- data.frame()

for (l in (1:5)) {
  output <- get_daily_qtile_estimates(taq, y = "sa_diff_trd", var = "ret", nlag = l)
  betas_lag <- output[1:5,c("coefficient","estimate", "Lag")]|> 
    mutate(Interval = str_extract(coefficient, "\\d+$")) |> 
    select(estimate, Lag, Interval)
  
  beta.table_lag <- output[1:5,] |> 
    mutate(Interval = str_extract(coefficient, "\\d+$"), 
           sd = paste0("(",round(std.error,3), ")"))
  
  taq_betas <- rbind(taq_betas, betas_lag)
  taq_beta.table <- rbind(taq_beta.table, beta.table_lag)
}

summary(taq_betas)

taq_beta.table.arrange <- taq_beta.table[order(taq_beta.table$Interval),]


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
## STEP 3: VISUALIZE REGRESSOR FOR ROBINHOOD AND TAQ

# log_holding changes
rh_lag_plot <- plot_daily_return_beta(rh_betas, xaxis = "Lag", clas = "Interval", ylim_d = -0.15, ylim_u = 0.45)
taq_lag_plot <- plot_daily_return_beta(taq_betas, xaxis = "Lag", clas = "Interval", ylim_d = -0.15, ylim_u = 0.45)


attention_plot <- rh_lag_plot + taq_lag_plot + 
  plot_layout(guides = "collect", axes = "collect") 

ggsave("plots/rh_ret.pdf", rh_lag_plot, 
       width = 8,
       height = 4.5,
       # units = "in", 
       limitsize = FALSE)

ggsave("plots/taq_ret.pdf", taq_lag_plot, 
       width = 8,
       height = 4.5,
       # units = "in", 
       limitsize = FALSE)


################################################################################
## STEP 4: ADJUST FOR CAR

## run regression from lag 1 to lag 5 data
rh_betas <- data.frame()
rh_beta.table <- data.frame()

for (l in (1:5)) {
  output <- get_daily_qtile_estimates(robintrack, y = "sa_holding_change", var = "car", nlag = l)
  betas_lag <- output[1:5,c("coefficient","estimate", "Lag")]|> 
    mutate(Interval = str_extract(coefficient, "\\d+$")) |> 
    select(estimate, Lag, Interval)
  beta.table_lag <- output[1:5,] |> 
    mutate(Interval = str_extract(coefficient, "\\d+$"), 
           sd = paste0("(",round(std.error,3), ")"))
  rh_betas <- rbind(rh_betas, betas_lag)
  rh_beta.table <- rbind(rh_beta.table, beta.table_lag)
}

summary(rh_betas)

rh_beta.table.arrange <- rh_beta.table[order(rh_beta.table$Interval),]


## run TAQ regression from lag 0 to lag 5 data
taq_betas <- data.frame()
taq_beta.table <- data.frame()

for (l in (1:5)) {
  output <- get_daily_qtile_estimates(taq, y = "sa_diff_trd", var = "car", nlag = l)
  betas_lag <- output[1:5,c("coefficient","estimate", "Lag")]|> 
    mutate(Interval = str_extract(coefficient, "\\d+$")) |> 
    select(estimate, Lag, Interval)
  
  beta.table_lag <- output[1:5,] |> 
    mutate(Interval = str_extract(coefficient, "\\d+$"), 
           sd = paste0("(",round(std.error,3), ")"))
  
  taq_betas <- rbind(taq_betas, betas_lag)
  taq_beta.table <- rbind(taq_beta.table, beta.table_lag)
}

taq_beta.table.arrange <- taq_beta.table[order(taq_beta.table$Interval),]

## print the table together
rh_beta.table.arrange$investor <- "rh"
taq_beta.table.arrange$investor <- "taq"
herding.beta.table <- rbind(rh_beta.table.arrange, taq_beta.table.arrange)

kable(reshape_betas(herding.beta.table), format = "latex", booktabs = TRUE, align = "c")
# R-SQR NUMBER
kable(taq_beta.table.arrange$adj.r.sqrd, format = "latex", booktabs = TRUE, align = "c")
kable(rh_beta.table.arrange$adj.r.sqrd, format = "latex", booktabs = TRUE, align = "c")
