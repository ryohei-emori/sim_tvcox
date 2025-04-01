setwd("/Users/emo_157/Desktop/simulation_coxph/sim_coxphtv")
#link : https://search.r-project.org/CRAN/refmans/coxed/html/sim.survdata.html
library(coxed);library(dplyr);library(timeROC);library(ggplot2);library(pec)
#library(SurvMetrics);

#### A simulated data set with time-varying covariates ####
set.seed(123)
betas <- rnorm(n = 10, mean = 0, sd = 0.1)
rep_times <- 25
covs_num <- 10
simdata <- sim.survdata(N = 10000,
                        T = 84,
                        type = "tvc",
                        #hazard.fun = NULL,
                        num.data.frames = rep_times+1,
                        fixed.hazard = TRUE,
                        knots = 15,
                        #spline = TRUE,
                        #X = NULL,
                        beta = betas,
                        xvars = covs_num,
                        #mu = 0,
                        #sd = 0.5,
                        #covariate = 10,
                        #low = 0,
                        #high = 1,
                        #compare = median,
                        censor = 0.85,
                        #censor.cond = FALSE
                        )

View(simdata[[1]]$data)
plot(simdata[[1]]$baseline$hazard)
sum(simdata[[1]]$data$failed)

simdata[[1]]$data %>% group_by(id) %>%
  summarise(end = max(end)) %>%
  summary(end)

eval_time <- 32

## Resampling for adjusting sample size in each dataset
sample_sizes <- c(3000, 5000, 7000, 10000)
sampled_data <- list()
#Surv_obj <- with(simdata[[rep_times+1]]$data, Surv(end, failed))
sampled_data[["val_sim"]] <- simdata[[rep_times+1]]$data %>%
  group_by(id) %>%
  summarise(
    start = min(start),  # 最初の観測時間
    end = max(end),    # 最後の観測時間
    across(starts_with("X"), first),
    failed = last(failed)
  )
#true_intensity <- sum(sampled_data[["val_sim"]]$failed) / nrow(sampled_data[["val_sim"]])

for (sample_size in sample_sizes) {
  # sampled_data[[as.character(sample_size)]] <- list(ids = list(),
  #                                                   data = list(),
  #                                                   tvcov_model = list(),
  #                                                   oos_auc = list(),
  #                                                   oos_ppv = list()) # initialization

  for (r in 1:rep_times) {
    unique_ids <- unique(simdata[[r]]$data$id)
    sampled_ids <- as.numeric(sample(unique_ids, size = sample_size))

    sampled_data[[as.character(sample_size)]]$ids[[r]] <- sampled_ids
    sampled_data[[as.character(sample_size)]]$data[[r]] <- simdata[[r]]$data %>% 
      filter(id %in% sampled_ids)
    sampled_data[[as.character(sample_size)]]$data[[r]]$events_num <- sum(sampled_data[[as.character(sample_size)]]$data[[r]]$failed)

    # estimation
    model <- coxph(Surv(start, end, failed) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10,
                   data = sampled_data[[as.character(sample_size)]]$data[[r]],
                   x = TRUE,
                   y = TRUE)
    sampled_data[[as.character(sample_size)]]$tvcov_model[[r]] <- model

    # prediction, score
    predict_tmp <- predict(object=model,
                           newdata=sampled_data[["val_sim"]],
                           type="expected",
                           se.fit=FALSE,
                           na.action=na.pass)

    sampled_data[[as.character(sample_size)]]$oos_auc[[r]] <- timeROC(T=sampled_data[["val_sim"]]$end,
                                                                      delta=sampled_data[["val_sim"]]$failed,
                                                                      marker=exp(-predict_tmp),
                                                                      cause=1,
                                                                      weighting="marginal",
                                                                      times=c(eval_time))$AUC
    
    # sampled_data[[as.character(sample_size)]]$oos_ppv[[r]] <- SeSpPPVNPV(cutpoint=as.numeric(quantile(predict_tmp,0.9)),
    #                                                                      T=sampled_data[["val_sim"]]$end,
    #                                                                      delta=sampled_data[["val_sim"]]$failed,
    #                                                                      marker=predict_tmp,
    #                                                                      cause=1,
    #                                                                      weighting = "marginal",
    #                                                                      times=c(eval_time),
    #                                                                      iid = FALSE)$PPV
    
    #sampled_data[[as.character(sample_size)]]$oos_bs[[r]] <- Brier(Surv_obj, predict_tmp, 60)
  }
}




#### evaluation ####
eval <- list()
for (sample_size in sample_sizes) {
  eval[[as.character(sample_size)]]$est_betas <- vector("list", covs_num)
  
  for (i in 1:covs_num) {
    for (r in 1:rep_times) {
      eval[[as.character(sample_size)]]$est_betas[[i]] <- c(eval[[as.character(sample_size)]]$est_betas[[i]], 
                                                           as.numeric(sampled_data[[as.character(sample_size)]]$tvcov_model[[r]]$coefficients[i])
                                                       )
    }
  }
}


for (sample_size in sample_sizes) {
  for (i in 1:covs_num) {
    difftmp <- betas[i] - eval[[as.character(sample_size)]]$est_betas[[i]]
    eval[[as.character(sample_size)]]$mse[i] <- mean(difftmp^2)
    eval[[as.character(sample_size)]]$variance[i] <- var(eval[[as.character(sample_size)]]$est_betas[[i]])
    eval[[as.character(sample_size)]]$bias[i] <- mean(difftmp)
  }
  for (r in 1:rep_times){
    eval[[as.character(sample_size)]]$oos_auc[r] <- as.numeric(sampled_data[[as.character(sample_size)]]$oos_auc[[r]][2])
    #eval[[as.character(sample_size)]]$oos_ppv[r] <- as.numeric(sampled_data[[as.character(sample_size)]]$oos_ppv[[r]][2])
  }
}



#### Right-censored Calibration ####
col1 = rgb(0.4, 0.6, 1, 0.5)   # 淡い青色の透過度50%
col2 = rgb(0.4, 1, 0.4, 0.5)   # 淡い緑色の透過度50%
col3 = rgb(1, 0.6, 0.2, 0.5)   # 淡いオレンジ色の透過度50%
col4 = rgb(1, 0.4, 0.4, 0.5)   # 淡い赤色の透過度50%

for (sample_size in sample_sizes) {
  eval[[as.character(sample_size)]]$calp <- calPlot(object = c(sampled_data[[as.character(sample_size)]]$tvcov_model),
                       time = 40,
                       data = sampled_data$val_sim,
                       formula = Surv(end, failed) ~ -1,
                       type = "risk",
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
                       xlab = "Predicted risk probabilities",
                       ylab = "Observed risk frequencies"
  )
  col_predictnames<- colnames(eval[[as.character(sample_size)]]$calp$predictions)
  
  for (r in 1:(rep_times)) {
    formula <- as.formula(paste(col_predictnames[1]," ~ 1 + ", col_predictnames[1+r]))
    model_summary <- summary(lm(formula, data=eval[[as.character(sample_size)]]$calp$predictions))
    eval[[as.character(sample_size)]]$oos_calIntercept[r] <- model_summary$coefficients[1]
    eval[[as.character(sample_size)]]$oos_calSlope[r] <- model_summary$coefficients[2]
  }
}


#### Right-censored AUC ####

# transforming into long-format
group_names <- setdiff(names(eval), c("calp","scores"))
data <- lapply(group_names, function(x) {
  data.frame(group = x, score = eval[[x]][["oos_auc"]])
})

eval[["scores"]]$oos_auc[["data"]] <- do.call(rbind, data)
eval[["scores"]]$oos_auc[["data"]]$group <- factor(eval[["scores"]]$oos_auc[["data"]]$group, levels = group_names)

eval[["scores"]]$oos_auc[["plot"]] <- ggplot(eval[["scores"]]$oos_auc[["data"]], aes(x = group, y = score, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +  # 透明度を追加
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  # ボックスプロットの調整
  labs(x = "Group", y = "Out-of-Sample AUC", title = "Comparison of Groups") +  # タイトルを追加
  theme_minimal(base_size = 14) +  # ベースフォントサイズの変更
  theme(legend.position = "right") +  # 凡例の位置調整
  theme(plot.title = element_text(hjust = 0.5)) +  # タイトルのスタイル調整
  scale_fill_brewer(palette = "Set2")  # 洗練された色のパレットに変更

plot(eval[["scores"]]$oos_auc[["plot"]])






## A simulated data set with time-varying coefficients
simdata <- sim.survdata(N=1000, T=100, type="tvbeta", num.data.frames = 1)
simdata$betas
