require(survival);require(dplyr)


invariant_data_transform <- function(data){
  data <- data %>% group_by(id) %>%
    summarise(start = min(start), # initial observed time
              end = max(end), # last observed time
              across(starts_with("X"), first),
              failed = last(failed))
  return(data)
}

fixed_weights_data_transform <- function(data){
  data <- data %>% group_by(id) %>% mutate(fixedweights=n())
  data$fixedweights <- (1/data$fixedweights)
  return(data)
}

dynamic_weights_data_transform <- function(data, trimming){
  #data$visited_start <- 0
  data$visited <- 1
  tmp_data <- data[data$failed!=1,]
  tmp_model <- coxph(Surv(start, end, visited, type = "counting") ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
                     data = tmp_data)
  tmp_predict <- predict(object=tmp_model,
                         newdata=data,
                         type="expected",
                         se.fit=FALSE,
                         na.action=na.pass)
  data$dynamicweights <- (1/pmax(tmp_predict, trimming))
  return(data)
}

proxy_data_transform <- function(data){
  data$visited <- 1
  tmp_data <- data[data$failed!=1,]
  tmp_model <- coxph(Surv(start, end, visited, type = "counting") ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
                     data = data)
  tmp_predict <- predict(object=tmp_model,
                         newdata=data,
                         type="expected",
                         se.fit=FALSE,
                         na.action=na.pass)
  data$proxy <- tmp_predict
  return(list(data, tmp_model))
}

estimation <- function(models, sampled_data, rep_times, epvs, trimming){
  model_functions <- list(
    naive_ivc = function(data) {
      tmp_iv_data <- invariant_data_transform(data)
      naive_ivc <- coxph(Surv(start, end, failed) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
                         data = tmp_iv_data,
                         x = TRUE, y = TRUE)
      return(naive_ivc)
    },
    naive_tvc = function(data) {
      coxph(Surv(start, end, failed) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
            data = data,
            x = TRUE, y = TRUE)
    },
    fixed_weighted_tvc = function(data) {
      tmp_fixed_weights_data <- fixed_weights_data_transform(data)
      fixed_weighted_tvc <- coxph(Surv(start, end, failed) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
                                  data = tmp_fixed_weights_data,
                                  x = TRUE, y = TRUE,
                                  weights = fixedweights)
      return(fixed_weighted_tvc)
    },
    dynamic_weighted_tvc = function(data) {
      tmp_dynamic_weights_data <- dynamic_weights_data_transform(data, trimming)
      #tmp_dynamic_weights_data <- tmp_dynamic_weights_data[tmp_dynamic_weights_data$dynamicweights!=0,]
      dynamic_weighted_tvc <- coxph(Surv(start, end, failed) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
                                  data = tmp_dynamic_weights_data,
                                  x = TRUE, y = TRUE,
                                  weights = dynamicweights)
      return(dynamic_weighted_tvc)
    },
    proxy_tvc = function(data) {
      tmp_proxy_data <- proxy_data_transform(data)
      proxy_data <- tmp_proxy_data[[1]]
      proxy_model <- tmp_proxy_data[[2]]
      proxy_tvc <- coxph(Surv(start, end, failed) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + proxy,
                                    data = proxy_data,
                                    x = TRUE, y = TRUE)
      return(list(proxy_tvc, proxy_model))
    }
  )
  for (epv in epvs) {
    for (r in 1:rep_times) {
      for (model in models) {
        if (model %in% names(model_functions)) {
          model_result <- model_functions[[model]](sampled_data[[as.character(epv)]]$data[[r]])
          if (model!="proxy_tvc"){
            sampled_data[[as.character(epv)]]$models[[model]][[r]] <- model_result  
          } else {
            sampled_data[[as.character(epv)]]$models[[model]][[r]] <- model_result[[1]]
            sampled_data[[as.character(epv)]]$models[["proxy_model"]][[r]] <- model_result[[2]]
          }
        }
      }
    }
  }
  return(sampled_data)
}
