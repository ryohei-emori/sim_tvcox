setwd("/Users/emo_157/Desktop/24_CPM_EHR/simulation_coxph/sim_coxphtv")
source("src/01_simulation_data.R")
source("src/02_compared_estimation.R")
source("src/03_prediction_evaluation.R")
source("src/04_drawing_figures.R")

#### requirments ####
set.seed(123)

betas <- rnorm(n = 10, mean = 1, sd = 0.3)
rep_times <- 25
covs_num <- 10
censoring_rate <- 0.85
epvs <- c("10","15","20","25","30","35")
models <- c("naive_ivc",
            "naive_tvc",
            "dynamic_weighted_tvc")


#### simdata ####
simdata <- gen_simdata(betas,
                       covs_num,
                       rep_times,
                       censoring_rate)
#plot(simdata[[1]]$baseline$hazard, type="l") #randomly chosen baseline hazard

epv_data <- epv_sampling(simdata, epvs, rep_times, covs_num)


## informative observation sampling: visits conditionally at random
io_data <- io_sampling_repeating(sampled_data=epv_data,
                                 epvs=epvs,
                                 rep_times=rep_times,
                                 dist="exp")
io_data[["val_sim"]] <- io_validation(simdata, rep_times)
#io_plot(io_data[["30"]]$data[[1]], kde=FALSE)
#summary(io_data$val_sim$end)

eval_time <- summary(io_data$val_sim$end)[[5]] #3rd quantile


#### estimation ####
io_ests <- estimation(models,
                      sampled_data = io_data,
                      rep_times = rep_times,
                      epvs = epvs,
                      trimming = 0.01)


#### prediction & evaluation ####
io_eval <- predict_evaluation(sim_data, 
                                  sampled_data=io_ests,
                                  rep_times=rep_times,
                                  models=models,
                                  eval_time,
                                  epvs = epvs)


#### drawing figures ####
auc <- aucplot(epvs=epvs,
               models=c("naive_ivc",
                        "naive_tvc",
                        "dynamic_weighted_tvc"),
               rep_times, eval=io_eval, type="boxplot")
auc[1]





