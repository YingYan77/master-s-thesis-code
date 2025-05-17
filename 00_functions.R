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



