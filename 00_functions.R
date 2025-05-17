################################################################################
## author:    Ying Yan
## file name: 000_functions.R
## summary:   this file contains all functions that are used in all coding files
## date:      6 Nov 2024
################################################################################



######################### funtions used in 02a and 03a #########################
get_lags_qtile_data <- function(dt, varname, extreme){
  var_lags <- grep(paste0("^", varname, "_lag_"), names(dt), value = TRUE)
  dt_clean <- dt[complete.cases(dt[, ..var_lags])]
  
  if (extreme == 0.10) {
    dt_clean[, (paste0("qtile_", var_lags)) := 
               lapply(.SD, function(var) findInterval(var, 
                                                      quantile(var, prob = c(.1,.2,.3,.4,.5,.6,.7,.8,.9), na.rm = TRUE), 
                                                      left.open = TRUE) + 1),
             by = .(date),
             .SDcols = var_lags]
  } else {
    dt_clean[, (paste0("qtile_", var_lags)) := 
               lapply(.SD, function(var) findInterval(var, 
                                                      quantile(var, prob = c(extreme,.25,.5,.75,(1-extreme)), na.rm = TRUE), 
                                                      left.open = TRUE) + 1),
             by = .(date),
             .SDcols = var_lags]
  }
  
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

get_daily_dcile_estimates <- function(dt, y, var, control, nlag) {
  
  # get dataframe
  lag.data <- dt[!is.na(get(paste0("qtile_", var, "_lag_", nlag)))]
  lag.data <- dummy_cols(lag.data, select_columns = paste0("qtile_", var, "_lag_", nlag))
  
  # dt <- get(paste0("robintrack_lag_", 1))
  
  # get lag quintle variables
  qitle.vars <- paste0("qtile_", var, "_lag_", nlag, "_", c(1,2,3,4,5,7,8,9,10), collapse = " + ")
  
  # add control variables if it exist
  mkt.ret.vars <- paste0("mkt_ret_lag_", 1:5, collapse = " + ")
  mkt.sd.vars <- paste0("mkt_std_dev_lag_", 1:5, collapse = " + ")
  
  # get return and volatility variables
  if (var == "ret" | var == "cumret" | var == "car") {
    ret.vars <- paste0("ret_lag_", setdiff(1:5, nlag), collapse = " + ")
    sd.vars <- paste0("volatility_lag_", setdiff(1:5, nlag), collapse = " + ")
    rh.vars <- paste0("sa_holding_change_lag_", 1:5, collapse = " + ")
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
                                  # ret.vars,
                                  sd.vars, 
                                  rh.vars,
                                  taq.vars,
                                  # mkt.ret.vars, 
                                  mkt.sd.vars, 
                                  "mcap", "BtM", "gp", "inv", "momentum", 
                                  sep = " + "
    )
  } else {
    lag.data <- lag.data[!is.na(get(paste0("qtile_", control, "_lag_", nlag)))]
    lag.data <- dummy_cols(lag.data, select_columns = paste0("qtile_", control, "_lag_", nlag))
    control.vars <- paste0("qtile_", control, "_lag_", nlag, "_", c(1,2,3,4,5,7,8,9,10), collapse = " + ")
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
    mutate(intercept = case_when(abs(t.value) >= 2.576 ~ paste0(round(estimate,3), "***"), 
                                 abs(t.value) >= 1.96 ~ paste0(round(estimate,3), "**"), 
                                 abs(t.value) >= 1.645 ~ paste0(round(estimate,3), "*"), 
                                 abs(t.value) < 1.645 ~ paste0(round(estimate,3))))
  
  if (var == "ret" | var == "cumret" | var == "car") {
    return(output.table)
  } else {
    ## calculate the difference 
    library(car)
    
    # Define the variable pairs to test
    lags <- c(1,2,3,4,5,7,8,9,10)
    results <- list()
    
    # Loop through each lag to compare the holding_change and diff_trd coefficients
    for (lag in c(1,2,3,4,5,7,8,9,10)) {
      var1 <- paste0("qtile_", var, "_lag_" , nlag, "_", lag)
      var2 <- paste0("qtile_", control, "_lag_" , nlag, "_", lag)
      hypothesis <- paste(var1, "=", var2)
      
      # Perform Wald test and store the result
      test_result <- linearHypothesis(model, hypothesis, vcov = cluster_varv)
      results[[paste0("lag_", lag)]] <- test_result
    }
    
    
    wald_test_table <- data.frame(
      interval = c(1,2,3,4,5,7,8,9,10),
      lag = nlag,
      X_variable = paste0("qtile_", var, "_lag_" , nlag, "_", lags),
      ctrl_variable = paste0("qtile_", control, "_lag_" , nlag, "_", lags),
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

reshape_betas <- function(dt) {
  library(dplyr)
  library(tidyr)
  
  # Ensure sue_qtile is treated as a factor or integer
  dt <- dt %>%
    mutate(Interval = as.integer(Interval), 
           Lag = as.integer(Lag), 
           col_name = paste0(investor, "_", Lag))
  
  # Reshape estimates
  estimates_wide <- dt %>%
    select(Interval, col_name, intercept) %>%
    pivot_wider(names_from = col_name, values_from = intercept)
  
  # Reshape standard errors
  stderr_wide <- dt %>%
    select(Interval, col_name, sd) %>%
    pivot_wider(names_from = col_name, values_from = sd)
  
  # list(estimates = estimates_wide, std_errors = stderr_wide)
  
  combined_result <- rbind(estimates_wide, stderr_wide)
  combined_result <- combined_result |> arrange(Interval)
  
  return(combined_result)
}


######################### funtions used in 02b and 03b #########################

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

# portfolio function for herding 
get_buyhold_portfolio <- function(dt, investor, weight = c("equal", "holding_changes"), 
                                  portfolio = c("1", "10", "10-1")){
  # get weight variables 
  weight_var <- investor
  
  # get the top and bottom deciles 
  bottom5 <- dt[get(paste0("qtile_", investor))==1, ]
  top5 <- dt[get(paste0("qtile_", investor))==6, ]
  
  # calculate the daily weighted average return 
  if (weight == "equal") {
    bottom_portfolio <- bottom5[, .(retp_lead_0 = mean(ret_lead_0, na.rm = TRUE)*100, 
                                    retp_lead_1 = mean(ret_lead_1, na.rm = TRUE)*100, 
                                    retp_lead_5 = mean((cumret_5d - 1)/5, na.rm = TRUE)*100, 
                                    retp_lead_10 = mean((cumret_10d - 1)/10, na.rm = TRUE)*100, 
                                    retp_lead_20 = mean((cumret_10d - 1)/20, na.rm = TRUE)*100), 
                                by = .(date)]
    
    top_portfolio <- top5[, .(retp_lead_0 = mean(ret_lead_0, na.rm = TRUE)*100, 
                              retp_lead_1 = mean(ret_lead_1, na.rm = TRUE)*100, 
                              retp_lead_5 = mean((cumret_5d - 1)/5, na.rm = TRUE)*100, 
                              retp_lead_10 = mean((cumret_10d - 1)/10, na.rm = TRUE)*100, 
                              retp_lead_20 = mean((cumret_10d - 1)/20, na.rm = TRUE)*100), 
                          by = .(date)]
    
    
  } else if (weight == "holding_changes") {
    bottom5[, change_weight := abs(get(weight_var)) / sum(abs(get(weight_var))), 
            by = .(date)]
    top5[, change_weight := abs(get(weight_var)) / sum(abs(get(weight_var))), 
         by = .(date)]
    
    bottom_portfolio <- bottom5[, .(retp_lead_0 = sum(change_weight*ret_lead_0, na.rm = TRUE)*100, 
                                    retp_lead_1 = sum(change_weight*ret_lead_1, na.rm = TRUE)*100, 
                                    retp_lead_5 = sum(change_weight*(cumret_5d - 1)/5, na.rm = TRUE)*100, 
                                    retp_lead_10 = sum(change_weight*(cumret_10d - 1)/10, na.rm = TRUE)*100,
                                    retp_lead_20 = sum(change_weight*(cumret_20d - 1)/20, na.rm = TRUE)*100), 
                                by = .(date)]
    
    top_portfolio <- top5[, .(retp_lead_0 = sum(change_weight*ret_lead_0, na.rm = TRUE)*100, 
                              retp_lead_1 = sum(change_weight*ret_lead_1, na.rm = TRUE)*100, 
                              retp_lead_5 = sum(change_weight*(cumret_5d - 1)/5, na.rm = TRUE)*100, 
                              retp_lead_10 = sum(change_weight*(cumret_10d - 1)/10, na.rm = TRUE)*100,
                              retp_lead_20 = sum(change_weight*(cumret_20d - 1)/20, na.rm = TRUE)*100), 
                          by = .(date)]
    
    
  } 
  
  long_short_portfolio <- top_portfolio[, .(date = date, 
                                            retp_lead_0 = retp_lead_0 - bottom_portfolio$retp_lead_0, 
                                            retp_lead_1 = retp_lead_1 - bottom_portfolio$retp_lead_1,
                                            retp_lead_5 = retp_lead_5 - bottom_portfolio$retp_lead_5, 
                                            retp_lead_10 = retp_lead_10 - bottom_portfolio$retp_lead_10, 
                                            retp_lead_20 = retp_lead_20 - bottom_portfolio$retp_lead_20)]
  
  if (portfolio == "10") {
    return(top_portfolio)
  } else if (portfolio == "1") {
    return(bottom_portfolio)
  } else if (portfolio == "10-1") {
    return(long_short_portfolio)
  }
}

# portfolio function for return tracing
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
  ## Newey-West adjusted t-test
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

# for herding
get_ret_alpha <- function(dt, varlists) {
  alphas <- data.frame()
  
  for (i in varlists) {
    for (j in c(0,1,5,10,20)) {
      for (l in c("equal", "holding_changes")) {
        for (k in c("10", "1")) {
          ret_portfolio <- get_buyhold_portfolio(dt, investor = i, weight = l, portfolio = k)
          ret_portfolio <- merge(ret_portfolio, ff6, by.x = c("date"), by.y = c("Date"))
          ret_portfolio[, Retp.Rf := get(paste0("retp_lead_",j)) - RF] 
          result <- get_alpha(ret_portfolio, "Retp.Rf", hperiod = j, investor = i, dcile = k, weight = l)
          
          alphas <- rbind(alphas, result)
        }
        for (n in c("Excess", "CAPM", "FF3", "FF6")) {
          ret_portfolio_1 <- get_buyhold_portfolio(dt, i, weight = l, portfolio = 1)
          ret_portfolio_1 <- merge(ret_portfolio_1, ff6, by.x = c("date"), by.y = c("Date"))
          ret_portfolio_1[, Retp.Rf := get(paste0("retp_lead_",j)) - RF] 
          
          ret_portfolio_10 <- get_buyhold_portfolio(dt, i, weight = l, portfolio = 10)
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
# for return tracing
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

########################### funtions used in 04b ###############################

## create the function to calculate the standard error clustered test average value
get_sue_holding_ret_result <- function(varlist, data, invlist) {
  test_result <- data.frame()
  
  for (investor in invlist) { 
    data[, sap_qtile := findInterval(get(investor), 
                                     quantile(get(investor), prob=c(.2,.4,.6,.8), na.rm = T), left.open = T) + 1, 
         by = .(datacqtr, earning_qtile)] 
    
    for (i in 1:5) { 
      for (j in 1:5) {
        subset_data <- data[earning_qtile == i & sap_qtile == j]
        for (var in varlist) { 
          formula <- as.formula(paste(var,"~ 1"))
          model <- lm(formula, data = subset_data)
          clustered_se <- vcovCL(model, cluster = subset_data$datacqtr)
          result <- coeftest(model, vcov. = clustered_se)
          
          output <- data.frame(
            var_name = var,
            estimate = result[1],
            std.error = result[2],
            t.value = result[3],
            p.value = result[4], 
            sue_qtile = i, 
            sap_qtile = j, 
            investor = investor
          )
          
          test_result <- rbind(test_result, output)
          
        }
      }
    }}
  
  test_result <- test_result |>
    mutate(mean = case_when(abs(p.value) <= 0.01 ~ paste0(round(estimate*100,3), "***"),
                            abs(p.value) <= 0.05 & abs(p.value) >= 0.01 ~ paste0(round(estimate*100,3), "**"),
                            abs(p.value) <= 0.1  & abs(p.value) >= 0.05 ~ paste0(round(estimate*100,3), "*"),
                            abs(p.value) > 0.1 ~ paste0(round(estimate*100,3))), 
           sd = paste0("(",round(std.error*100,3), ")"))
  
  return(test_result)
}

get_sue_qtile_diff_ret <- function(varlist, data, invlist) {
  test_result <- data.frame()
  
  for (investor in invlist) { 
    data[, sap_qtile := findInterval(get(investor), 
                                     quantile(get(investor), prob=c(.2,.4,.6,.8), na.rm = T), left.open = T) + 1, 
         by = .(datacqtr, earning_qtile)] 
    
    for (i in 1:5) { 
      
      subset_data <- data[earning_qtile == i & (sap_qtile == 1 | sap_qtile == 5)]
      subset_data$sap_indicator <- ifelse(subset_data$sap_qtile == 5, 1, 0)
      
      for (var in varlist) {
        formula <- as.formula(paste(var,"~ sap_indicator | datacqtr"))
        
        model <- feols(formula, data = subset_data, cluster = ~datacqtr)
        coefs <- model$coefficients
        cluster_se <- model$se
        
        output <- data.frame(
          var_name = var,
          estimate = coefs[1],
          std.error = cluster_se[1],
          t.value = coefs / cluster_se,
          sue_qtile = i, 
          sap_qtile = 6, 
          investor = investor
        )
        
        test_result <- rbind(test_result, output)
      }
    }
    
    for (j in 1:5) { 
      
      subset_data <- data[(earning_qtile == 1 | earning_qtile == 5) & sap_qtile == j]
      subset_data$sue_indicator <- ifelse(subset_data$earning_qtile == 5, 1, 0)
      
      for (var in varlist) {
        formula <- as.formula(paste(var,"~ sue_indicator | datacqtr"))
        
        model <- feols(formula, data = subset_data, cluster = ~datacqtr)
        coefs <- model$coefficients
        cluster_se <- model$se
        
        output <- data.frame(
          var_name = var,
          estimate = coefs[1],
          std.error = cluster_se[1],
          t.value = coefs / cluster_se,
          sue_qtile = 6, 
          sap_qtile = j, 
          investor = investor
        )
        
        test_result <- rbind(test_result, output)
      }
    }}
  
  
  test_result <- test_result |> 
    mutate(mean = case_when(abs(t.value) >= 2.576 ~ paste0(round(estimate*100,3), "***"), 
                            abs(t.value) >= 1.96 ~ paste0(round(estimate*100,3), "**"), 
                            abs(t.value) >= 1.645 ~ paste0(round(estimate*100,3), "*"), 
                            abs(t.value) < 1.645 ~ paste0(round(estimate*100,3))), 
           sd = paste0("(",round(std.error*100,3), ")"))
  
  return(test_result)
}

reshape_sue_timelag <- function(dt) {
  library(dplyr)
  library(tidyr)
  
  # Ensure sue_qtile is treated as a factor or integer
  dt <- dt %>%
    mutate(sue_qtile = as.integer(sue_qtile), 
           sap_qtile = as.integer(sap_qtile),
           row_name = paste0(var_name, "_", sue_qtile))
  
  # Reshape estimates
  estimates_wide <- dt %>%
    select(row_name, sap_qtile, mean) %>%
    pivot_wider(names_from = sap_qtile, values_from = mean)
  
  # Reshape standard errors
  stderr_wide <- dt %>%
    select(row_name, sap_qtile, sd) %>%
    pivot_wider(names_from = sap_qtile, values_from = sd)
  
  # list(estimates = estimates_wide, std_errors = stderr_wide)
  
  combined_result <- rbind(estimates_wide, stderr_wide)
  combined_result <- combined_result |> arrange(row_name)
  
  return(combined_result)
}
