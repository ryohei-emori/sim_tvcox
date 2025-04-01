require(coxed)

gen_simdata <- function(betas, covs_num, rep_times, censoring_rate) {
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
                          censor = censoring_rate,
                          #censor.cond = FALSE
  )
  return(simdata)
}

val_simdata <- function(simdata, rep_times){
  val_simdata <- simdata[[rep_times+1]]$data %>%
                                  group_by(id) %>%
                                  summarise(
                                    start = min(start),  # initial observed time
                                    end = max(end),    # last observed time
                                    across(starts_with("X"), first),
                                    failed = last(failed))
  return(val_simdata)
}

# simple proportional sampling (not informative sampling)
sampled_simdata <- function(simdata, epvs, rep_times){
  sampled_data <- list()
  for (epv in epvs) {
    for (r in 1:rep_times) {
      unique_ids <- unique(simdata[[r]]$data$id)
      sampled_ids <- as.numeric(sample(unique_ids, size = epv))
      id_counts <- table(simdata[[r]]$data$id)
      
      sampled_data[[as.character(epv)]]$data[[r]] <- simdata[[r]]$data %>% 
        filter(id %in% sampled_ids)
      sampled_data[[as.character(epv)]]$data[[r]]$survp <- unlist(lapply(sampled_ids, function(id) {
        num_records <- id_counts[as.character(id)]
        survp_row <- simdata[[r]]$ind.survive[id, 1:num_records]
        return(survp_row)
      }))
      sampled_data[[as.character(epv)]]$events_num[[r]] <- sum(sampled_data[[as.character(epv)]]$data[[r]]$failed)
    }
  }
  return(sampled_data)
}

# epv-based sampling : proportional and epv sampling
epv_sampling <- function(simdata, epvs, rep_times, covs_num) {
  sampled_data <- list()
  for (epv in epvs) {
    for (r in 1:rep_times) {
      final_failed <- simdata[[r]]$data %>%
        group_by(id) %>%
        summarize(failed = last(failed))
      
      # failed sampling
      failed_counts <- sum(final_failed$failed == 1)
      failed_ids <- unique(simdata[[r]]$data$id[simdata[[r]]$data$failed == 1])
      events_rate <- failed_counts / length(unique(simdata[[r]]$data$id))
      required_events <- epv * covs_num
      sampled_failed_ids <- sample(failed_ids, min(length(failed_ids), required_events))
      
      # non-failed sampling
      non_failed_epv <- (length(sampled_failed_ids) / events_rate) - length(sampled_failed_ids)
      non_failed_ids <- setdiff(unique(simdata[[r]]$data$id), failed_ids)
      sampled_non_failed_ids <- sample(non_failed_ids, min(length(non_failed_ids), round(non_failed_epv)))
      
      # gathering failed & non-failed sample
      sampled_ids <- c(sampled_failed_ids, sampled_non_failed_ids)
      sampled_records <- simdata[[r]]$data[simdata[[r]]$data$id %in% sampled_ids, ]
      
      sampled_data[[as.character(epv)]]$data[[r]] <- sampled_records
      # saved survival probability
      id_counts <- table(simdata[[r]]$data$id)
      sampled_data[[as.character(epv)]]$data[[r]]$survp <- unlist(lapply(sampled_ids, function(id) {
        num_records <- id_counts[as.character(id)]
        survp_row <- simdata[[r]]$ind.survive[id, 1:num_records]
        return(survp_row)
      }))
    }
    }
    return(sampled_data)
}


# informative observation sampling
io_sampling_repeating <- function(sampled_data, epvs, rep_times, dist="exp"){
  io_sampled_data <- list()
  
  for (epv in epvs) {
    print(paste0("Now sampling epv ",epv))
    for (r in 1:rep_times){
      # generating from individual visits at random observed process
      io_sampled_data[[as.character(epv)]]$data[[r]] <- io_sampling(sampled_data[[as.character(epv)]]$data[[r]],
                                                                            dist=dist)
    }
  }
  return(io_sampled_data)
}

io_sampling <- function(data, dist="exp", ind.survf=TRUE) {
  ids <- unique(data$id)
  all_obs <- list()
  
  for (i in ids) {
    ind_data <- data[data$id == i, ]
    time <- 0
    ind_end <- max(ind_data$end)
    obs_list <- list()
    
    while (time + 1 < ind_end) {
      ind_intensity <- 7 * ind_data[ind_data$start == time, "survp"]
      if (dist=="exp"){
        lambda <- (1 / ind_intensity)
        interval <- ceiling(rexp(1, rate = lambda))
      } else if (dist=="lnorm"){
        interval <- ceiling(rlnorm(1, meanlog=log(ind_intensity), sdlog=1))
      }
      time <- time + interval
      
      if (time + 1 < ind_end) {
        matched <- ind_data[ind_data$start == time, ]
        if (nrow(matched) > 0) {
          obs_list[[length(obs_list) + 1]] <- matched
        }
      } else {
        obs_list[[length(obs_list) + 1]] <- tail(ind_data, n = 1)
        break
      }
    }
    all_obs[[i]] <- do.call(rbind, obs_list)
  }
  data <- do.call(rbind, all_obs)
  data <- update_end_values(data)
  return(data)
}

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

io_validation <- function(simdata, rep_times){
  
  id_counts <- table(simdata[[rep_times+1]]$data$id)
  ids <- unique(simdata[[rep_times+1]]$data$id)
  simdata[[rep_times+1]]$data$survp <- unlist(lapply(ids, function(id) {
    num_records <- id_counts[as.character(id)]
    survp_row <- simdata[[rep_times+1]]$ind.survive[id, 1:num_records]
    return(survp_row)}))
  
  init_obs_time_data <- simdata[[rep_times+1]]$data %>%
    filter(start == 0) %>%
    group_by(id) %>%
    summarise(
      ind.intensity = 5 * survp,
      lambda = 1/3 + (1 / ind.intensity),
      interval = ceiling(rexp(1, rate = lambda)),
      init_obs_time = 0 + interval
    ) %>%
    ungroup()

  updated_simdata <- simdata[[rep_times+1]]$data %>%
    left_join(init_obs_time_data, by = "id") %>%
    mutate(start = init_obs_time) %>%
    arrange(id, start)
  
  io_val_simdata <- updated_simdata %>%
    group_by(id) %>%
    summarise(
      start = first(start),  # init_obs_time
      end = max(end),        # last observed time
      across(starts_with("X"), first),
      failed = last(failed)
    )
  
  return(io_val_simdata)
}


## Illustration of informative observation sample
io_plot <- function(data, kde=TRUE){

  data <- data %>% 
      group_by(id) %>%
      summarize(
        id = first(id),
        followup_time = max(end) - min(start),
        avg_interval= mean(end - start),
        failed = as.integer(any(failed==1))
      )

  if (kde){
    p <- ggplot(data, aes(x=avg_interval, fill = as.factor(failed))) +
      geom_density(color="black", alpha=0.3) +
      scale_fill_manual(values=c("0"="blue", "1"="red"),
                        labels=c("Non Event", "Event")) +
      labs(#title="Avg. Interval by Individuals",
           x="Avg. Interval",
           y="Density",
           fill="Event Status") +
      theme_minimal() +
      theme(axis.text.x = element_text(size = rel(1.2), colour = "black"),
                    axis.text.y = element_text(size = rel(1.2), colour = "black"))

  } else {
    p <- ggplot(data, aes(x=avg_interval, fill = as.factor(failed))) +
      geom_histogram(binwidth = 0.1, color="black", alpha=0.3) +
      scale_fill_manual(values=c("0"="blue", "1"="red"),
                        labels=c("Non Event", "Event")) +
      labs(#title="Avg. Interval by Individuals",
           x="Avg. Interval",
           y="Count",
           fill="Event Status") +
      theme_minimal() +
      theme(axis.text.x = element_text(size = rel(1.2), colour = "black"),
                    axis.text.y = element_text(size = rel(1.2), colour = "black"))
  }
  return(p)
}
