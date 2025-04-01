require(timeROC);require(pec);require(survcompare);require(SurvMetrics)

# Function to calculate the inverse probability of censoring weights
ipc_weights <- function(time, km_fit) {
  km_prob <- summary(km_fit, times = time)$surv
  return(1 / km_prob)
}

# Function to calculate the Brier score at a specific time considering censoring
brier_score_at_time <- function(time, predicted_probs, end, status, km_fit) {
  n <- length(end)
  brier_score <- 0
  for (i in 1:n) {
    surv_prob <- predicted_probs[i]
    weight <- if (end[i] >= time) ipc_weights(time, km_fit) else 0
    if (status[i] == 1 && end[i] <= time) {
      brier_score <- brier_score + (1 - surv_prob)^2 * weight
    } else if (end[i] >= time) {
      brier_score <- brier_score + (surv_prob)^2 * weight
    }
  }
  return(brier_score / n)
}

# General function to calculate time-dependent Brier score
calculate_brier_score <- function(times, predicted_probs, end, status) {
  actual_data <- data.frame(end = end, status = status)
  km_fit <- survfit(Surv(end, 1 - status) ~ 1, data = actual_data)
  brier_scores <- sapply(times, function(t) brier_score_at_time(t, predicted_probs, end, status, km_fit))
  brier_scores_df <- data.frame(time = times, brier_score = brier_scores)
  return(brier_scores_df)
}

predict_evaluation <- function(sim_data, sampled_data, rep_times, models, eval_time, epvs){
  eval = list()
  col1 = rgb(0.4, 0.6, 1, 0.5)   # 淡い青色の透過度50%
  col2 = rgb(0.4, 1, 0.4, 0.5)   # 淡い緑色の透過度50%
  col3 = rgb(1, 0.6, 0.2, 0.5)   # 淡いオレンジ色の透過度50%
  col4 = rgb(1, 0.4, 0.4, 0.5)   # 淡い赤色の透過度50%
  
  for (epv in epvs){
    for (model in models) {
      
      if (model %in% c("naive_ivc","naive_tvc")){
      eval[[as.character(epv)]]$models[[model]]$calp <- calPlot(object = c(sampled_data[[as.character(epv)]]$models[[model]]),
                                                      time = eval_time,
                                                      data = sampled_data[["val_sim"]],
                                                      formula = Surv(end, failed) ~ -1,
                                                      type = "survival",
                                                      method = "nne",
                                                      #q = 20,
                                                      pseudo = TRUE,
                                                      showPseudo = FALSE,
                                                      pseudo.col = "gray",
                                                      bars = FALSE,
                                                      add =FALSE,
                                                      legend = FALSE,
                                                      col = c(rep(col1, rep_times),
                                                              rep(col2, rep_times)),
                                                      lwd = c(rep(2, rep_times),
                                                              rep(2, rep_times)),
                                                      lty = c(rep("twodash", rep_times),
                                                              rep("dashed", rep_times)),
                                                      xlab = "Predicted survival probabilities",
                                                      ylab = "Observed survival frequencies"
                                                      )
      calpredict <- eval[[as.character(epv)]]$models[[model]]$calp$predictions
      }

      for (r in 1:rep_times) {
          if (model == "proxy_tvc"){
            sampled_data[["val_sim"]]$visited <- 1
            tmp_proxy <- predict(object=sampled_data[[as.character(epv)]]$models[["proxy_model"]][[r]],
                                 newdata=sampled_data[["val_sim"]],
                                 type="expected",
                                 se.fit=FALSE,
                                 na.action=na.pass)
            sampled_data[["val_sim"]]$proxy <- tmp_proxy
          }
          tmp_predict <- predict(object=sampled_data[[as.character(epv)]]$models[[model]][[r]],
                                 newdata=sampled_data[["val_sim"]],
                                 type="survival",
                                 se.fit=FALSE,
                                 na.action=na.pass)

          if (model %in% c("naive_ivc","naive_tvc")){
            cal_formula <- as.formula(paste(calpredict[1]," ~ 1 + ", calpredict[1+r]))
            cal_model_summary <- summary(lm(cal_formula,
                                            data=eval[[as.character(epv)]]$models[[model]]$calp$predictions))
          } else {
            eval[[as.character(epv)]]$models[[model]]$predictions[[r]] <- tmp_predict
            cal_model_summary <- summary(lm(unlist(calpredict[1])~ 1 + tmp_predict))
          }
          
          eval[[as.character(epv)]]$models[[model]]$oos_calIntercept[r] <- cal_model_summary$coefficients[1]
          eval[[as.character(epv)]]$models[[model]]$oos_calSlope[r] <- cal_model_summary$coefficients[2]
  
          eval[[as.character(epv)]]$models[[model]]$oos_auc[r] <- timeROC(T=sampled_data[["val_sim"]]$end,
                                                                            delta=sampled_data[["val_sim"]]$failed,
                                                                            marker=tmp_predict,
                                                                            cause=1,
                                                                            weighting="marginal",
                                                                            times=c(eval_time))$AUC[2]
          
          #eval[[as.character(epv)]]$models[[model]]$oos_bs[r] <- calculate_brier_score(end = sampled_data[["val_sim"]]$end,
          #                                                                                     status = sampled_data[["val_sim"]]$failed,
          #                                                                                     times = c(eval_time),
          #                                                                                     predicted_probs = tmp_predict)
          
          #eval[[as.character(epv)]]$models[[model]]$oos_bs[r] <- Brier(Surv(sampled_data[["val_sim"]]$end,
          #                                                                           sampled_data[["val_sim"]]$failed),
          #                                                                      pre_sp=tmp_predict,
          #                                                                      t_star=c(eval_time))
          #eval[[as.character(epv)]]$models[[model]]$oos_bs[r] <- surv_brierscore(y_predicted_newdata=tmp_predict,
          #                                                                   #df_brier_train,
          #                                                                   df_newdata=sampled_data[["val_sim"]]$failed,
          #                                                                   time_points=c(eval_time),
          #                                                                   weighted = FALSE)
      }
    }
  }
  return(eval)
}
