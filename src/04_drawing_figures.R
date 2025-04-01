require(survival)
require(ggplot2)
require(reshape2)

aucplot <- function(epvs, models, rep_times, eval, type="boxplot"){
  
  result_df <- data.frame()
  
  for (epv in epvs){
    for (model in models){
      tmp_data <- data.frame(EPV = rep(epv, each=rep_times),
                             Model = rep(model, each=rep_times),
                             oos_auc = eval[[as.character(epv)]]$models[[model]]$oos_auc,
                             Rep = 1:rep_times)
      result_df <- rbind(result_df, tmp_data)
    }
  }
  if (type=="boxplot"){
    p <- ggplot(result_df, aes(x = factor(EPV), y = oos_auc, fill = Model)) +
      geom_boxplot(width = 0.2, position = position_dodge(width = 0.8)) +
      scale_fill_brewer(palette = "Set2") +
      labs(#title = "AUC by Sample Size and Model with Adjusted Box Width",
        x = "EPV",
        y = "AUC with Out-of-sample") +
      theme_minimal() +
      theme(legend.title = element_blank(),
            axis.text.x = element_text(size = rel(1.2), colour = "black"),
            axis.text.y = element_text(size = rel(1.2), colour = "black"))
  } else if (type=="violin") {
    p <- ggplot(result_df, aes(x = factor(EPV), y = oos_auc, fill = Model)) +
      geom_violin(trim = FALSE, adjust = 1.5, position = position_dodge(0.7)) +
      scale_fill_brewer(palette = "Set2") +
      labs(#title = "AUC by Sample Size and Model (Violin Plot)",
        x = "EPV",
        y = "AUC with Out-of-sample") +
      theme_minimal() +
      theme(legend.title = element_blank(),
            axis.text.x = element_text(size = rel(1.2), colour = "black"),
            axis.text.y = element_text(size = rel(1.2), colour = "black"))
  }
  return(list(p, result_df))
}


# pseudo_calplot <- function(models, sampled_data, rep_times, pseudo_values){
# 
#     all_predictions <- lapply(models, function(model_list) {
#       sapply(model_list, function(model) predict(model, newdata = sampled_data[["val_sim"]], type = "expected"))})
#     
#     all_predictions_df <- do.call(rbind, all_predictions)
#     all_predictions_melted <- melt(all_predictions, variable.name = "Model", value.name = "Predicted")
#     colnames(all_predictions_melted)[colnames(all_predictions_melted) == "L1"] <- "Model"
#     all_predictions_melted$Model <- factor(all_predictions_melted$Model)
#     
#     col1 = rgb(0.4, 0.6, 1, 0.5)   # 淡い青色の透過度50%
#     col2 = rgb(0.4, 1, 0.4, 0.5)   # 淡い緑色の透過度50%
#     col3 = rgb(1, 0.6, 0.2, 0.5)   # 淡いオレンジ色の透過度50%
#     col4 = rgb(1, 0.4, 0.4, 0.5)   # 淡い赤色の透過度50%
#     colors <- c(col1, col2, col3, col4)
#     linetypes <- c("solid", "dashed", "dotdash", "longdash")  # モデルごとの線の形状
#     
#     all_predictions_melted$Pseudo <- rep(as.numeric(pseudo_values), times = length(models)*rep_times)
#     
#     p <- ggplot(all_predictions_melted, aes(x = Predicted, y = Pseudo, group = Model)) +
#       geom_point(aes(color = Model), alpha = 0.1) +
#       geom_smooth(aes(linetype = Model, color = Model), method = "loess", se = FALSE) +  # スムージングを適用
#       scale_color_manual(values = colors) +
#       scale_linetype_manual(values = linetypes) +
#       geom_abline(intercept = 0, slope = 1, color = "gray") +
#       coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
#       theme_minimal() +
#       labs(x = "Predicted Risk Probabilities", y = "Observed Risk Frequencies")
#     
#     return(p)
# }

#pseudo_calplot(io_sampled_data[["3000"]]$models,
#               io_sampled_data, rep_times=2, pseudo_values=unlist(io_eval[["3000"]]$models$naive_ivc$calp$predictions[1]))


# Step 1: Group predictions into quantiles
# q <- 10  # number of quantiles
# groups <- quantile(predictions, probs=seq(0, 1, length.out=q+1))
# cut_points <- cut(predictions, breaks=groups, include.lowest=TRUE, labels=FALSE)
# 
# # Step 2: Calculate mean predicted probabilities for each quantile group
# mean_predictions <- tapply(predictions, cut_points, mean)
# 
# # Step 3: Calculate mean pseudo-observations for each quantile group
# mean_pseudo_obs <- tapply(pseudo_obs, cut_points, mean)
# 
# # Step 4: Create plotFrames data frame
# plotFrames <- data.frame(Predicted=mean_predictions, Observed=mean_pseudo_obs)
# 
# 
# 
# ## example
# cl <- calPlot(object = c(io_sampled_data[["10000"]]$models[["naive_ivc"]]),
#               time = eval_time,
#               data = sampled_data[["val_sim"]],
#               formula = Surv(end, failed) ~ -1,
#               type = "survival",
#               method = "nne",
#               #q = 20,
#               pseudo = TRUE,
#               showPseudo = TRUE,
#               pseudo.col = "gray",
#               bars = FALSE,
#               add =FALSE,
#               legend = FALSE,
#               col = c(rep(col1, rep_times),
#                       rep(col2, rep_times)),
#               lwd = c(rep(2, rep_times),
#                       rep(2, rep_times)),
#               lty = c(rep("twodash", rep_times),
#                       rep("dashed", rep_times)),
#               xlab = "Predicted survival probabilities",
#               ylab = "Observed survival frequencies"
# )
# 
# groups <- quantile(cl$predictions$Model.1, probs=seq(0, 1, length.out=20+1))
# cut_points <- cut(cl$predictions$Model.1, breaks=groups, include.lowest=TRUE, labels=FALSE)
# mean_pseudo_obs <- tapply(cl$predictions$t.33, cut_points, mean)
# mean_predictions <- tapply(cl$predictions$Model.1, cut_points, mean)
# plotFrames <- data.frame(Predicted=mean_predictions, Observed=mean_pseudo_obs)
# 
# plot(mean_predictions, mean_pseudo_obs)

## 
# groups <- quantile(io_eval[["7000"]]$models$dynamic_weighted_tvc$predictions[[1]], probs=seq(0, 1, length.out=20+1))
# cut_points <- cut(io_eval[["7000"]]$models$dynamic_weighted_tvc$predictions[[1]], breaks=groups, include.lowest=TRUE, labels=FALSE)
# mean_pseudo_obs <- tapply(io_eval[["7000"]]$models$naive_ivc$calp$predictions$t.33, cut_points, mean)
# mean_predictions <- tapply(io_eval[["7000"]]$models$dynamic_weighted_tvc$predictions[[1]], cut_points, mean)
# plotFrames <- data.frame(Predicted=mean_predictions, Observed=mean_pseudo_obs)
# 
# plot(mean_predictions, mean_pseudo_obs)

