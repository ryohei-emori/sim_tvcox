
## source data: tv-covariates
require(coxed)
require(MASS)

Generate_X <- function(n = 10000, t_max = 100, p = 10) {
  
  random_matrix <- matrix(rnorm(p^2), ncol = p)
  Sigma <- cov(random_matrix)
  
  X_continuous <- matrix(0, nrow = n * t_max, ncol = 6)
  X_binary_tv <- matrix(0, nrow = n * t_max, ncol = 2)
  X_binary_ti <- matrix(0, nrow = n, ncol = 2)
  
  for (i in 1:n) {
    X_continuous[(i-1)*t_max+1, ] <- mvrnorm(1, mu = rep(0, 6), Sigma = Sigma[1:6, 1:6])
    for (t in 2:t_max) {
      X_continuous[(i-1)*t_max+t, ] <- X_continuous[(i-1)*t_max+t-1, ] + mvrnorm(1, mu = rep(0, 6), Sigma = Sigma[1:6, 1:6])
    }
    
    X_binary_tv[(i-1)*t_max+1:t_max, ] <- mvrnorm(t_max, mu = rep(0, 2), Sigma = Sigma[7:8, 7:8]) > 0
    
    X_initial <- mvrnorm(1, mu = rep(0, 2), Sigma = Sigma[9:10, 9:10])
    X_binary_ti[i, ] <- X_initial > 0
  }
  
  X <- data.frame(
    id = rep(1:n, each = t_max),
    time = rep(1:t_max, times = n),
    X1 = X_continuous[, 1],
    X2 = X_continuous[, 2],
    X3 = X_continuous[, 3],
    X4 = X_continuous[, 4],
    X5 = X_continuous[, 5],
    X6 = X_continuous[, 6],
    X7 = X_binary_tv[, 1],
    X8 = X_binary_tv[, 2],
    X9 = rep(X_binary_ti[, 1], each = t_max),
    X10 = rep(X_binary_ti[, 2], each = t_max)
  )
  
  return(X)
}

source_data <- function(N = 10000, T = 100, betas, censoring_rate = 0.0) {
  
  source <- sim.survdata(N = N,
                         T = T,
                         type = "tvc",
                         num.data.frames = 1,
                         fixed.hazard = TRUE,
                         knots = 15,
                         beta = betas,
                         xvars = 10,
                         censor = censoring_rate
  )
  
  X <- Generate_X(n = N, t_max = T, p = 10)
  baseline_hazard <- source$baseline$hazard
  
  X$linear_predictor <- as.vector(as.matrix(X[, c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10")]) %*% betas)
  X$hazard <- baseline_hazard * exp(X$linear_predictor)
  
  U <- runif(N)
  X <- X %>%
    group_by(id) %>%
    mutate(
      cumulative_hazard = cumsum(hazard),
      survp = 1 - exp(-cumulative_hazard),
      event_time = time[which(cumulative_hazard >= -log(U[which(id == first(id))]))[1]],
      event_time = ifelse(is.na(event_time), T + 1, event_time),  # T + 1に変更
      event = as.numeric(time == event_time)
    ) %>%
    filter(time <= event_time) %>%
    ungroup()
  
  result_data <- X %>%
    group_by(id) %>%
    mutate(
      start = time,
      end = ifelse(time == event_time, time + 1, lead(time, default = T)),
      end = pmin(end, T),
      failed = event
    ) %>%
    filter(time <= event_time & time != end) %>%
    filter(!(start == end & end == T)) %>%  # start = end = T のレコードを除外
    ungroup()
  
  validation_data <- X %>%
    group_by(id) %>%
    summarise(
      start = min(time),
      end = ifelse(max(event_time) > T, T, max(event_time) + 1),  # 修正
      failed = ifelse(max(event_time) > T, 0, 1),  # event_time > Tならfailed = 0
      event_time = max(event_time),
      across(starts_with("X"), first)
    ) %>%
    ungroup()
  
  return(list(result_data = result_data, validation_data = validation_data))
}

## EPV-based sampling
EPV_sample <- function(source, epvs, p = 10) {
  sampled_data <- list()
  source_data <- source$result_data
  
  print(str(source_data))
  failed_counts <- sum(source_data$failed)
  failed_ids <- unique(source_data$id[source_data$failed == 1])
  events_rate <- if (length(unique(source_data$id)) > 0) failed_counts / length(unique(source_data$id)) else 0
  print(events_rate)
  
  for (epv in epvs) {
  
    required_events <- as.numeric(epv) * p
    sampled_failed_ids <- sample(failed_ids, min(length(failed_ids), required_events))
    non_failed_epv <- (length(sampled_failed_ids) / events_rate) - length(sampled_failed_ids)
    non_failed_ids <- setdiff(unique(source_data$id), failed_ids)
    
    if (length(non_failed_ids) > 0) {
      sampled_non_failed_ids <- sample(non_failed_ids, min(length(non_failed_ids), round(non_failed_epv)))
    } else {
      sampled_non_failed_ids <- character(0)
    }
    
    sampled_ids <- c(sampled_failed_ids, sampled_non_failed_ids)
    sampled_records <- source_data[source_data$id %in% sampled_ids, ]
    sampled_data[[epv]] <- sampled_records
  }
  return(sampled_data)
}

## informative sampling
update_end_values <- function(data) {
  ids <- unique(data$id)
  
  for (i in ids) {
    ind_data <- data[data$id == i, ]
    ind_data <- ind_data[order(ind_data$start), ]
    
    ind_data$end[-nrow(ind_data)] <- ind_data$start[-1]
    data[data$id == i, ] <- ind_data
  }
  return(data)
}

IO_sample <- function(EPV_source, epvs, dist="exp") {
  
  IO_sample <- list()
  for (epv in epvs){
    data <- EPV_source[[epv]]
    ids <- unique(data$id)
    all_obs <- list()
    
    for (i in ids) {
      ind_data <- data[data$id == i, ]
      time <- 1
      ind_end <- max(ind_data$end)
      obs_list <- list()
      matched <- ind_data[ind_data$start == time, ]
      obs_list[[length(obs_list)+1]] <- matched 
      
      while (time < ind_end) {
        ind_intensity <- as.numeric(ind_data[ind_data$start == time, "survp"])
        
        if (is.na(ind_intensity) || is.infinite(ind_intensity)) {
          print(paste("Invalid ind_intensity detected for id:", i, "time:", time, "value:", ind_intensity))
          break
        }
        
        if (dist == "exp") {
          lambda <- 5 + (1 / ind_intensity)
          
          if (is.na(lambda) || is.infinite(lambda)) {
            print(paste("Invalid lambda detected for id:", i, "time:", time, "value:", lambda))
            break
          }
          
          interval <- ceiling(rexp(1, rate = lambda))
        } else if (dist == "lnorm") {
          interval <- ceiling(rlnorm(1, meanlog = log(ind_intensity), sdlog = 1))
        }
        time <- time + interval
        if (time <= ind_end) {
          matched <- ind_data[ind_data$start == time, ]
          obs_list[[length(obs_list) + 1]] <- matched 
        } else if (time > ind_end) {
          #obs_list[[length(obs_list) + 1]] <- tail(ind_data, n = 1)
          break
        }
      }
      all_obs[[i]] <- do.call(rbind, obs_list)
    }
    data <- do.call(rbind, all_obs)
    #data <- update_end_values(data)
    IO_sample[[epv]] <- data
  }
  return(IO_sample)
}


## estimation 
invariant_transform <- function(data){
  data <- data %>% group_by(id) %>%
    reframe(start = min(start), # initial observed time
              end = max(end), # last observed time
              across(starts_with("X"), first),
              failed = last(failed))
  return(data)
}

weights_transform <- function(data, trimming){
  tmp_data <- data[data$failed!=1,]
  tmp_data$visited <- 1
  tmp_model <- coxph(Surv(start, end, visited, type = "counting") ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
                     data = tmp_data)
  tmp_predict <- predict(object=tmp_model,
                         newdata=data %>% mutate(visited=1),
                         type="expected",
                         se.fit=FALSE,
                         na.action=na.pass)
  data$weights <- pmin(1/pmax(tmp_predict, trimming))
  return(data)
}

model_list <- list(
IV_Cox = function(data) {
  tmp_iv_data <- invariant_transform(data)
  coxph(Surv(start, end, failed) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
        data = tmp_iv_data,
        x = TRUE, y = TRUE)
},
TV_Cox = function(data) {
  coxph(Surv(start, end, failed) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
        data = data,
        x = TRUE, y = TRUE)
},
WTV_Cox = function(data) {
  tmp_weights <- weights_transform(data, trimming=1e-3)
  #tmp_weights <- tmp_weights[tmp_weights$dynamicweights!=0,]
  coxph(Surv(start, end, failed) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
        data = tmp_weights,
        x = TRUE, y = TRUE,
        weights = weights)
}
)

## Prediction
Est_Pred <- function(epvs, models, IO_data, r){
  
  Pred_res <- IO_data[["validation"]]
  for (epv in epvs){
    data <- IO_data[[epv]]
    print(epv)
  for (model in models){
    ## Estimation
    Est_model <- model_list[[model]](data)
    print(model)
    ## Predicition
    if (model %in% c("IV_Cox","TV_Cox", "WTV_Cox")){
      tmp_predict <- predict(object = Est_model,
                             newdata = Pred_res,
                             type="survival",
                             se.fit=FALSE,
                             na.action=na.pass)
      Pred_res[[model]] <- tmp_predict
    } else if (model %in% c("IV_RSF","TV_RSF","WTV_RSF")){
      tmp_predict <- 
      Pred_res[[model]] <- tmp_predict
    }
  }}
  write.csv(Pred_res, paste0(r,"_Pred_res.csv"))
}


## replication wrapper
rep_wrapper <- function(R, epvs){
  
  for (r in R){
    print(paste0("Replication time: ", r))
    print("Sampling source data")
    source <- source_data(N=10000, T=100, betas, censoring_rate=0.0)
    print("Sampling based on each EPVs")
    EPV_source <- EPV_sample(source = source, epvs = epvs, p=10)
    print("Informative Sampling")
    IO_data <- IO_sample(EPV_source = EPV_source, epvs = epvs, dist="exp")
    IO_data[["validation"]] <- source$validation_data
    #print(str(IO_data))
    print("Start estimating and Predicting")
    Est_Pred(epvs = epvs, models = models, IO_data = IO_data, r = r)
  }
}

## drawing results
AUCPlot <- function(epvs, models, R, Scores_path, type="boxplot"){

  result_df <- read_csv(Scores_path)
  if (type=="boxplot"){
    p <- ggplot(result_df, aes(x = factor(EPV), y = oos_auc, fill = Model)) +
      geom_boxplot(width = 0.2, position = position_dodge(width = 0.8)) +
      scale_fill_brewer(palette = "Set2") +
      labs(#title = "",
        x = "EPV",
        y = "AUC with Out-of-sample") +
      theme_minimal() +
      theme(legend.title = element_blank(),
            axis.text.x = element_text(size = rel(1.2), colour = "black"),
            axis.text.y = element_text(size = rel(1.2), colour = "black"))
  }
  return(list(p, result_df))
}

## execution 
betas <- rnorm(n = 10, mean = -2, sd = 0.1)
epvs <- c("15","25","35")
models <- c("IV_Cox",
            "TV_Cox",
            "WTV_Cox")
R <- 1
rep_wrapper(R, epvs)
