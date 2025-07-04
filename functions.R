capitalize_first <- function(x) {
  sapply(x, function(name) {
    paste0(toupper(substring(name, 1, 1)), substring(name, 2))
  })
}

generate_train_test <- function(data, train_ratio = 0.7) {
  train_rows <- round(train_ratio * nrow(data))
  train_indices <- sample(seq_len(nrow(data)), size = train_rows)
  test_indices <- setdiff(seq_len(nrow(data)), train_indices)
  
  train_data <- data[train_indices, , drop = FALSE]
  test_data <- data[test_indices, , drop = FALSE]
  
  return(list(
    train_df = train_data,
    test_df = test_data
  ))
}


robust_feature_elimination <- function(data, 
                                       repeat_runs = 20, 
                                       directory = ".", 
                                       file_prefix = "output",
                                       min_features = 1) {
  
  n_cores <- 16 # Assuming 16 cores
  # Setup parallel backend
  cl <- makeCluster(n_cores)
  clusterEvalQ(cl, {
    .libPaths("/path/to/Rlib")
  })
  registerDoParallel(cl)
  
  # Run RFE repeats in parallel
  all_results <- foreach(i = 1:repeat_runs, .packages = c("survival", "ranger", "dplyr"),
                         .export = c("generate_train_test", "perform_feature_elimination")) %dopar% {
                           cat("RFE run:", i, "\n")
                           
                           # Generate train/test split
                           train_test_data <- generate_train_test(data)
                           train_df <- train_test_data$train_df
                           test_df <- train_test_data$test_df
                           
                           feature_c_index_list <- list()
                           count <- 1
                           
                           # Run your existing perform_feature_elimination function
                           res <- perform_feature_elimination(train_df, test_df, count, feature_c_index_list)
                           res
                         }
  
  stopCluster(cl)
  
  # Aggregate c-index for each number of features across runs
  num_features_all <- sort(unique(unlist(lapply(all_results, function(res) sapply(res, function(x) length(x$variables))))))
  
  agg_df <- data.frame()
  
  for (nf in num_features_all) {
    c_indices_nf <- c()
    for (run in all_results) {
      c_index_for_nf <- NA
      for (step in run) {
        if (length(step$variables) == nf) {
          c_index_for_nf <- step$c_index
          break
        }
      }
      c_indices_nf <- c(c_indices_nf, c_index_for_nf)
    }
    c_indices_nf <- c_indices_nf[!is.na(c_indices_nf)]
    if(length(c_indices_nf) < 2) next
    
    agg_df <- rbind(agg_df, data.frame(
      Num_Features = nf,
      Mean_CIndex = mean(c_indices_nf),
      SD_CIndex = sd(c_indices_nf),
      Lower_CI = quantile(c_indices_nf, 0.025),
      Upper_CI = quantile(c_indices_nf, 0.975)
    ))
  }
  
  # Filter out small feature sets (optional)
  agg_df <- agg_df %>% filter(Num_Features >= min_features)
  
  # Find feature count with highest mean c-index
  best_nf <- agg_df$Num_Features[which.max(agg_df$Mean_CIndex)]
  
  # Stability selection: get most frequent features in best subsets across runs
  feature_counts <- list()
  for (run in all_results) {
    for (step in run) {
      if (length(step$variables) == best_nf) {
        feature_counts <- c(feature_counts, step$variables)
        break
      }
    }
  }
  feature_freq <- sort(table(unlist(feature_counts)), decreasing = TRUE)
  final_features <- names(feature_freq)[1:best_nf]
  
  # Plot mean c-index with error bars
  p <- ggplot(agg_df, aes(x = Num_Features, y = Mean_CIndex)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2) +
    labs(x = "Number of Features", y = "Mean C-Index", 
         title = paste0(file_prefix, ": RFE Performance Across Runs")) +
    theme_minimal()
  
  # Ensure output directory exists
  if (!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE)
  }
  
  # Save results
  ggsave(file.path(directory, paste0(file_prefix, "_robust_RFE_plot.png")), p, width = 10, height = 6)
  write.csv(agg_df, file.path(directory, paste0(file_prefix, "_robust_RFE_summary.csv")), row.names = FALSE)
  write.csv(data.frame(Feature = final_features), 
            file.path(directory, paste0(file_prefix, "_selected_features.csv")), row.names = FALSE)
  
}


evaluate_performance <- function(train_df, test_df, n_bootstrap = 30, conf_level = 0.95) {
  # Extract survival information
  train_survival <- Surv(train_df$overall_survival, train_df$deceased)
  test_survival  <- Surv(test_df$overall_survival, test_df$deceased)
  
  # Remove survival columns from predictors
  predictors_train <- train_df[, !(colnames(train_df) %in% c("overall_survival", "deceased"))]
  predictors_test  <- test_df[,  !(colnames(test_df) %in%  c("overall_survival", "deceased"))]
  
  train_x_matrix <- as.matrix(predictors_train)
  test_x_matrix  <- as.matrix(predictors_test)
  
  # Bootstrap CI function
  bootstrap_ci <- function(predictions, true_values, n_bootstrap, conf_level) {
    c_index_values <- numeric(n_bootstrap)
    for (i in 1:n_bootstrap) {
      sample_indices <- sample(1:length(true_values), length(true_values), replace = TRUE)
      boot_true_values <- true_values[sample_indices]
      boot_predictions <- predictions[sample_indices]
      c_index_values[i] <- concordance(boot_true_values ~ boot_predictions)$concordance
    }
    
    lower_ci <- quantile(c_index_values, (1 - conf_level) / 2)
    upper_ci <- quantile(c_index_values, 1 - (1 - conf_level) / 2)
    mean_c_index <- mean(c_index_values)
    std_c_index <- sd(c_index_values)
    
    p_value <- 2 * (1 - pnorm(mean_c_index, mean = 0.5, sd = std_c_index))  # assuming comparison to 0.5 baseline
    
    return(list(
      mean_c_index = mean_c_index,
      std_c_index = std_c_index,
      lower_ci = lower_ci,
      upper_ci = upper_ci,
      p_value = p_value
    ))
  }
  
  # Cox PH
  eval_cox_model <- coxph(train_survival ~ ., data = predictors_train)
  cox_predictions <- predict(eval_cox_model, newdata = predictors_test)
  cox_results <- bootstrap_ci(cox_predictions, test_survival, n_bootstrap, conf_level)
  
  # Ranger
  mtry_value <- floor(sqrt(ncol(predictors_train)))
  ranger_df <- predictors_train
  ranger_df$overall_survival <- train_df$overall_survival
  ranger_df$deceased <- train_df$deceased
  
  eval_r_fit <- ranger(
    formula = Surv(overall_survival, deceased) ~ .,
    data = ranger_df,
    mtry = mtry_value,
    importance = "permutation",
    splitrule = "maxstat",
    min.node.size = 50,
    num.trees = 1000
  )
  
  ranger_predictions <- predict(eval_r_fit, data = predictors_test)$survival
  ranger_results <- bootstrap_ci(ranger_predictions, test_survival, n_bootstrap, conf_level)
  
  # GLMBoost
  glmboost_df <- predictors_train
  glmboost_df$overall_survival <- train_df$overall_survival
  glmboost_df$deceased <- train_df$deceased
  
  eval_boosted_cox_model <- glmboost(Surv(overall_survival, deceased) ~ ., data = glmboost_df, family = CoxPH())
  glmboost_predictions <- predict(eval_boosted_cox_model, newdata = predictors_test)
  glmboost_results <- bootstrap_ci(glmboost_predictions, test_survival, n_bootstrap, conf_level)
  
  # Elastic Net
  elastic_net_model <- cv.glmnet(x = train_x_matrix, y = train_survival, family = "cox", alpha = 0.5)
  elastic_net_predictions <- predict(elastic_net_model, newx = test_x_matrix, s = "lambda.min")
  elastic_net_results <- bootstrap_ci(elastic_net_predictions, test_survival, n_bootstrap, conf_level)
  
  # Results
  cindex_df <- data.frame(
    CIndexCox = cox_results$mean_c_index,
    CIndexRanger = ranger_results$mean_c_index,
    CIndexGLMBoost = glmboost_results$mean_c_index,
    CIndexElasticNet = elastic_net_results$mean_c_index,
    StdCox = cox_results$std_c_index,
    StdRanger = ranger_results$std_c_index,
    StdGLMBoost = glmboost_results$std_c_index,
    StdElasticNet = elastic_net_results$std_c_index,
    LowerCI_Cox = cox_results$lower_ci,
    LowerCI_Ranger = ranger_results$lower_ci,
    LowerCI_GLMBoost = glmboost_results$lower_ci,
    LowerCI_ElasticNet = elastic_net_results$lower_ci,
    UpperCI_Cox = cox_results$upper_ci,
    UpperCI_Ranger = ranger_results$upper_ci,
    UpperCI_GLMBoost = glmboost_results$upper_ci,
    UpperCI_ElasticNet = elastic_net_results$upper_ci,
    PValueCox = cox_results$p_value,
    PValueRanger = ranger_results$p_value,
    PValueGLMBoost = glmboost_results$p_value,
    PValueElasticNet = elastic_net_results$p_value
  )
  
  return(cindex_df)
}


combine_results <- function(results) {
  # Extract results
  feature_cox_results <- results$feature_cox_results
  feature_ranger_results <- results$feature_ranger_results
  feature_glmboost_results <- results$feature_glmboost_results
  feature_elasticnet_results <- results$feature_elasticnet_results
  
  c_index_cox_results <- results$c_index_cox_results
  c_index_ranger_results <- results$c_index_ranger_results
  c_index_glmboost_results <- results$c_index_glmboost_results
  c_index_elasticnet_results <- results$c_index_elasticnet_results
  
  # Get minimum length across all result vectors
  all_lengths <- c(
    length(feature_cox_results), length(feature_ranger_results),
    length(feature_glmboost_results), length(feature_elasticnet_results),
    length(c_index_cox_results), length(c_index_ranger_results),
    length(c_index_glmboost_results), length(c_index_elasticnet_results)
  )
  
  if (min(all_lengths) == 0) {
    stop("At least one of the result vectors is empty. Cannot compute summary.")
  }
  
  min_length <- min(all_lengths)
  
  # Truncate all to common min length
  feature_selection_df <- data.frame(
    FeatureCox = feature_cox_results[1:min_length],
    FeatureRanger = feature_ranger_results[1:min_length],
    FeatureGLMBoost = feature_glmboost_results[1:min_length],
    FeatureElasticNet = feature_elasticnet_results[1:min_length]
  )
  
  feature_cindex_df <- data.frame(
    CIndexCox = c_index_cox_results[1:min_length],
    CIndexRanger = c_index_ranger_results[1:min_length],
    CIndexGLMBoost = c_index_glmboost_results[1:min_length],
    CIndexElasticNet = c_index_elasticnet_results[1:min_length]
  )
  
  # Calculate standard deviations
  cox_feature_sd <- sd(feature_selection_df$FeatureCox)
  ranger_feature_sd <- sd(feature_selection_df$FeatureRanger)
  boost_feature_sd <- sd(feature_selection_df$FeatureGLMBoost)
  penalised_feature_sd <- sd(feature_selection_df$FeatureElasticNet)
  
  cox_Cindex_sd <- sd(feature_cindex_df$CIndexCox)
  ranger_Cindex_sd <- sd(feature_cindex_df$CIndexRanger)
  boost_Cindex_sd <- sd(feature_cindex_df$CIndexGLMBoost)
  penalised_Cindex_sd <- sd(feature_cindex_df$CIndexElasticNet)
  
  # Calculate means
  mean_feature_values <- colMeans(feature_selection_df)
  mean_cindex_values <- colMeans(feature_cindex_df)
  
  # Sample sizes
  n_feature <- sapply(feature_selection_df, length)
  n_cindex <- sapply(feature_cindex_df, length)
  
  # 95% Confidence Intervals
  ci_feature <- 1.96 * c(cox_feature_sd, ranger_feature_sd, boost_feature_sd, penalised_feature_sd) / sqrt(n_feature)
  ci_cindex <- 1.96 * c(cox_Cindex_sd, ranger_Cindex_sd, boost_Cindex_sd, penalised_Cindex_sd) / sqrt(n_cindex)
  
  # Compute bounds
  lower_ci_feature <- mean_feature_values - ci_feature
  upper_ci_feature <- mean_feature_values + ci_feature
  
  lower_ci_cindex <- mean_cindex_values - ci_cindex
  upper_ci_cindex <- mean_cindex_values + ci_cindex
  
  # Summary data frames
  feature_selected_result_summary <- data.frame(
    Model = names(mean_feature_values),
    MeanFeatureSelected = mean_feature_values,
    StandardDeviation = c(cox_feature_sd, ranger_feature_sd, boost_feature_sd, penalised_feature_sd),
    LowerCI = lower_ci_feature,
    UpperCI = upper_ci_feature
  )
  
  feature_cindex_result_summary <- data.frame(
    Model = names(mean_cindex_values),
    MeanCIndex = mean_cindex_values,
    StandardDeviation = c(cox_Cindex_sd, ranger_Cindex_sd, boost_Cindex_sd, penalised_Cindex_sd),
    LowerCI = lower_ci_cindex,
    UpperCI = upper_ci_cindex
  )
  
  return(list(
    FeatureSelected = feature_selected_result_summary,
    CIndex = feature_cindex_result_summary
  ))
}


plot_results <- function(no_fs_results, multi_results, uni_results, rsf_vi_results, rsf_md_results, rsf_vh_results, custom_name) {
  
  # Create directory if it doesn't exist
  dir_name <- paste0(custom_name, "_plots")
  if (!dir.exists(dir_name)) {
    dir.create(dir_name)
  }
  
  # Extract C-Index and feature selection results
  extract_cindex <- function(results, source_name) {
    data.frame(
      Source = source_name,
      MeanCIndex = results$CIndex$MeanCIndex,
      SDCIndex = results$CIndex$StandardDeviation,
      LowerCI = results$CIndex$LowerCI,
      UpperCI = results$CIndex$UpperCI
    )
  }
  
  extract_features <- function(results, source_name) {
    data.frame(
      Source = source_name,
      MeanFeatureSelected = results$FeatureSelected$MeanFeatureSelected,
      SDFeatureSelected = results$FeatureSelected$StandardDeviation,
      LowerCI = results$FeatureSelected$LowerCI,
      UpperCI = results$FeatureSelected$UpperCI
    )
  }
  # Extract and merge results, adding the source name to each
  c_index_results <- rbind(
    extract_cindex(no_fs_results, "No FS"),
    extract_cindex(multi_results, "Multi"),
    extract_cindex(uni_results, "Uni"),
    extract_cindex(rsf_vi_results, "RSF VI"),
    extract_cindex(rsf_md_results, "RSF MD"),
    extract_cindex(rsf_vh_results, "RSF VH")
  )
  feature_results <- rbind(
    extract_features(no_fs_results, "No FS"),
    extract_features(multi_results, "Multi"),
    extract_features(uni_results, "Uni"),
    extract_features(rsf_vi_results, "RSF VI"),
    extract_features(rsf_md_results, "RSF MD"),
    extract_features(rsf_vh_results, "RSF VH")
  )
  # Save merged results to CSV files
  write.csv(c_index_results, file.path(dir_name, paste0(custom_name, "_cindex_results.csv")), row.names = FALSE)
  write.csv(feature_results, file.path(dir_name, paste0(custom_name,"_feature_results.csv")), row.names = FALSE)
  
  # Extract results from the input lists
  no_fs_feature_selected_result_summary <- no_fs_results$FeatureSelected
  no_fs_feature_cindex_result_summary <- no_fs_results$CIndex
  multi_feature_selected_result_summary <- multi_results$FeatureSelected
  multi_feature_cindex_result_summary <- multi_results$CIndex
  uni_feature_selected_result_summary <- uni_results$FeatureSelected
  uni_feature_cindex_result_summary <- uni_results$CIndex
  rsf_cindex_result_summary <- rsf_vi_results$CIndex
  rsf_feature_selected_result_summary <- rsf_vi_results$FeatureSelected
  rsfmd_cindex_result_summary <- rsf_md_results$CIndex
  rsfmd_feature_selected_result_summary <- rsf_md_results$FeatureSelected
  rsfvh_cindex_result_summary <- rsf_vh_results$CIndex
  rsfvh_feature_selected_result_summary <- rsf_vh_results$FeatureSelected
  
  # Combine c-index results from the four data frames
  c_index_combined <- cbind(
    no_fs_feature_cindex_result_summary$MeanCIndex, 
    multi_feature_cindex_result_summary$MeanCIndex, 
    uni_feature_cindex_result_summary$MeanCIndex, 
    rsf_cindex_result_summary$MeanCIndex, 
    rsfmd_cindex_result_summary$MeanCIndex, 
    rsfvh_cindex_result_summary$MeanCIndex
  )
  
  c_index_df <- as.data.frame(c_index_combined)
  rownames(c_index_df) <- c("CoxPH", "Random Forest", "GLMBoost", "ElasticNet")
  colnames(c_index_df) <- c("No FS", "Multivariate","Univariate", "RF Var Imp", "RF Min Depth", "RF Max Stat")
  c_index_df$Model <- rownames(c_index_df)
  rownames(c_index_df) <- NULL
  g3 <- melt(c_index_df)
  gg <- ggplot(g3, aes(variable, Model, fill = value)) + 
    geom_tile(color = "white") +
    geom_text(aes(label = round(value, 2)), size = 3, color = "black") + 
    scale_fill_gradient2(low = "yellow", high = "blue", limits = c(min(g3$value), 1)) +
    theme_minimal() +
    labs(title = "Mean Performance of Models + Feature Selection", x = "Feature Selection Method", y = "Machine Learning Algorithm", fill = NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(hjust = 1),
          axis.title = element_text(size = 10),
          axis.ticks = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = "right",
          panel.grid = element_blank(),
          panel.border = element_blank(),
          plot.margin = margin(30, 30, 30, 30, unit = "pt")) +
    coord_equal(ratio = 0.7)
  
  ggsave(file.path(dir_name, paste0(custom_name, "_cindex_heatmap.pdf")), gg, width = 8, height = 6, device = "pdf")
  # Combine feature selection results from the four data frames
  feature_combined <- cbind(
    no_fs_feature_selected_result_summary$MeanFeatureSelected, 
    multi_feature_selected_result_summary$MeanFeatureSelected, 
    uni_feature_selected_result_summary$MeanFeatureSelected, 
    rsf_feature_selected_result_summary$MeanFeatureSelected, 
    rsfmd_feature_selected_result_summary$MeanFeatureSelected, 
    rsfvh_feature_selected_result_summary$MeanFeatureSelected
  )
  
  feature_df <- as.data.frame(feature_combined)
  rownames(feature_df) <- c("CoxPH", "Random Forest", "GLMBoost", "ElasticNet")
  colnames(feature_df) <- c("No FS", "Multivariate","Univariate", "RF Var Imp", "RF Min Depth", "RF Max Stat")
  feature_df$Model <- rownames(feature_df)
  rownames(feature_df) <- NULL
  g4 <- melt(feature_df)
  
  gg2 <- ggplot(g4, aes(variable, Model, fill = value)) + 
    geom_tile(color = "white") +
    geom_text(aes(label = round(value, 2)), size = 3, color = "black") + 
    scale_fill_gradient2(low = "yellow", high = "purple", limits = c(min(g4$value), max(g4$value))) + 
    theme_minimal() + 
    labs(title = "Mean Features Selected by Models + Feature Selection", x = "Feature Selection Method", y = "Machine Learning Algorithm", fill = NULL) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(hjust = 1),
          axis.title = element_text(size = 10),
          axis.ticks = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = "right",
          panel.grid = element_blank(),
          panel.border = element_blank(),
          plot.margin = margin(30, 30, 30, 30, unit = "pt")) +
    coord_equal(ratio = 0.7)
  
  ggsave(file.path(dir_name, paste0(custom_name, "_feature_heatmap.pdf")), gg2, width = 8, height = 6, device = "pdf")
}


perform_parallel_feature_selection <- function(data, survival_data, repeats, folds, feature_selection, fs_name = "unknown") {

  num_cores <- 16 # Assuming 16 cores
  cl <- makeCluster(num_cores)
  # Export .libPaths to each worker
  clusterEvalQ(cl, {
    .libPaths("/path/to/Rlib")
  })
  registerDoParallel(cl)
  all_indices <- expand.grid(rep = 1:repeats, fold = 1:folds)
  # Precompute folds only once per repeat
  folds_list <- lapply(1:repeats, function(rep) {
    createFolds(factor(data$deceased), k = folds)
  })
  
  all_results <- foreach(i = 1:nrow(all_indices), .combine = list, .multicombine = TRUE,
                         .packages = c('caret', 'survival', 'glmnet', 'mboost', 'ranger', 'randomForestSRC'),
                         .export = c("feature_selection")) %dopar% {
                           
                           rep <- all_indices$rep[i]
                           fold <- all_indices$fold[i]
                           
                           fold_indices <- folds_list[[rep]]
                           train_indices <- unlist(fold_indices[-fold])
                           test_indices <- unlist(fold_indices[fold])
                           
                           train_data <- data[train_indices, ]
                           test_data <- data[test_indices, ]
                           
                           train_survival <- Surv(train_data$overall_survival, train_data$deceased)
                           test_survival <- Surv(test_data$overall_survival, test_data$deceased)
                           
                           # his includes survival info — needed for feature selection functions like mboost, glmnet (with Surv)
                           train_df_for_selection <- train_data[, !(names(train_data) %in% c("sampleID"))]
                           
                           # is drops survival columns — used for modeling
                           train_df_for_model <- train_data[, !(names(train_data) %in% c("sampleID", "overall_survival", "deceased"))]
                           test_df_for_model  <- test_data[,  !(names(test_data)  %in% c("sampleID", "overall_survival", "deceased"))]
                           
                           feature_selected <- feature_selection(train_df_for_selection)
                           selected_features <- feature_selected$feature
                           
                           
                           train_x_matrix <- as.matrix(train_df_for_model[, selected_features, drop = FALSE])
                           test_x_matrix  <- as.matrix(test_df_for_model[,  selected_features, drop = FALSE])
                           train_df <- as.data.frame(train_x_matrix)
                           test_df  <- as.data.frame(test_x_matrix)
                           
                           # COXPH
                           train_df$overall_survival <- train_data$overall_survival
                           train_df$deceased <- train_data$deceased
                           cox_formula <- as.formula("Surv(overall_survival, deceased) ~ .")
                           eval_cox_model <- coxph(cox_formula, data = train_df)
                           cox_predictions <- predict(eval_cox_model, newdata = test_df)
                           c_index_cox <- concordance(Surv(test_data$overall_survival, test_data$deceased) ~ cox_predictions)$concordance
                           selected_features_cox <- names(coef(eval_cox_model))[coef(eval_cox_model) != 0]
                           
                           # RANGER
                           ranger_df <- train_df[, selected_features, drop = FALSE]
                           ranger_df$overall_survival <- train_data$overall_survival
                           ranger_df$deceased <- train_data$deceased
                           
                           eval_r_fit <- ranger(
                             formula = Surv(overall_survival, deceased) ~ .,
                             data = ranger_df,
                             mtry = floor(sqrt(ncol(train_df))),
                             importance = "permutation",
                             splitrule = "maxstat",
                             min.node.size = 50,
                             num.trees = 1000
                           )
                           
                           ranger_predictions <- predict(eval_r_fit, data = test_df)$survival
                           c_index_ranger_result <- concordance(Surv(test_data$overall_survival, test_data$deceased) ~ ranger_predictions)
                           c_index_ranger <- if (length(c_index_ranger_result$concordance) > 1) {
                             mean(c_index_ranger_result$concordance)
                           } else {
                             c_index_ranger_result$concordance
                           }
                           selected_features_ranger_df <- data.frame(
                             Variable = names(eval_r_fit$variable.importance),
                             Importance = eval_r_fit$variable.importance
                           )
                           selected_features_ranger <- selected_features_ranger_df[selected_features_ranger_df$Importance > 0, ]
                           
                           # GLMBOOST
                           glmboost_model <- glmboost(Surv(overall_survival, deceased) ~ ., data = train_df, family = CoxPH())
                           glmboost_predictions <- predict(glmboost_model, newdata = test_df)
                           c_index_glmboost <- concordance(Surv(test_data$overall_survival, test_data$deceased) ~ glmboost_predictions)$concordance
                           selected_features_boosted_cox <- names(coef(glmboost_model))[coef(glmboost_model) != 0]
                           
                           # ELASTICNET
                           elastic_net_model <- cv.glmnet(x = train_x_matrix, y = train_survival, family = "cox", alpha = 0.5)
                           predictions <- predict(elastic_net_model, newx = test_x_matrix, s = "lambda.min")
                           c_index <- concordance(test_survival ~ predictions)$concordance
                           selected_features_elasticnet <- coef(elastic_net_model, s = "lambda.min")
                           selected_features_elasticnet <- as.matrix(selected_features_elasticnet)
                           selected_features_elasticnet_df <- data.frame(
                             feature = rownames(selected_features_elasticnet),
                             coefficient = selected_features_elasticnet
                           )
                           names(selected_features_elasticnet_df)[2] <- "coefficient"
                           selected_features_elasticnet <- selected_features_elasticnet_df[selected_features_elasticnet_df$coefficient != 0, ]
                           
                           list(
                             selected_features = feature_selected$feature,
                             c_index_elasticnet = c_index,
                             c_index_cox = c_index_cox,
                             c_index_ranger = c_index_ranger,
                             c_index_glmboost = c_index_glmboost,
                             feature_elasticnet = nrow(selected_features_elasticnet),
                             feature_cox = length(selected_features_cox),
                             feature_ranger = nrow(selected_features_ranger),
                             feature_glmboost = length(selected_features_boosted_cox)
                           )
                         }
  
  stopCluster(cl)
  
  flat_results <- Filter(Negate(is.null), all_results)
  
  selected_features_list <- lapply(flat_results, function(res) res$selected_features)
  c_index_elasticnet_results <- sapply(flat_results, function(res) res$c_index_elasticnet)
  c_index_cox_results <- sapply(flat_results, function(res) res$c_index_cox)
  c_index_ranger_results <- sapply(flat_results, function(res) res$c_index_ranger)
  c_index_glmboost_results <- sapply(flat_results, function(res) res$c_index_glmboost)
  feature_elasticnet_results <- sapply(flat_results, function(res) res$feature_elasticnet)
  feature_cox_results <- sapply(flat_results, function(res) res$feature_cox)
  feature_ranger_results <- sapply(flat_results, function(res) res$feature_ranger)
  feature_glmboost_results <- sapply(flat_results, function(res) res$feature_glmboost)
  
  return(list(
    selected_features_list = selected_features_list,
    c_index_elasticnet_results = c_index_elasticnet_results,
    c_index_cox_results = c_index_cox_results,
    c_index_ranger_results = c_index_ranger_results,
    c_index_glmboost_results = c_index_glmboost_results,
    feature_elasticnet_results = feature_elasticnet_results,
    feature_cox_results = feature_cox_results,
    feature_ranger_results = feature_ranger_results,
    feature_glmboost_results = feature_glmboost_results
  ))
}

univariate <- function(data) {
  # Exclude survival columns from features
  feature_columns <- setdiff(colnames(data), c("overall_survival", "deceased"))
  uni_cox_df <- data.frame(feature = character(),
                           coefficients = numeric(),
                           p_value = numeric(),
                           c_index = numeric(),
                           stringsAsFactors = FALSE)
  for (feature in feature_columns) {
    # Create proper formula and supply data argument
    formula <- reformulate(termlabel = feature, response = "Surv(overall_survival, deceased)")
    cox <- try(coxph(formula, data = data))
    if (inherits(cox, "try-error")) {
      next  
    }
    p_value <- summary(cox)$logtest["pvalue"]
    c_index <- summary(cox)$concordance["C"]
    coefficients <- coef(cox)
    result_row <- data.frame(feature = feature,
                             coefficients = coefficients,
                             p_value = p_value,
                             c_index = c_index)
    
    uni_cox_df <- rbind(uni_cox_df, result_row)
  }
  # Apply FDR correction 
  uni_cox_df$adj_p_value <- p.adjust(uni_cox_df$p_value, method = "fdr")
  # Return only significant, non-zero features
  significant_uni_cox_df <- subset(uni_cox_df, adj_p_value < 0.05 & coefficients != 0)
  return(significant_uni_cox_df)
}


multivariate <- function(data) {
  # Remove survival columns from predictors
  predictor_data <- data[, !(colnames(data) %in% c("overall_survival", "deceased"))]
  # Build Cox model
  eval_cox_model <- tryCatch({
    coxph(Surv(data$overall_survival, data$deceased) ~ ., data = predictor_data)
  }, error = function(e) {
    message("Cox model did not converge: ", e$message)
    return(NULL)
  })
  if (is.null(eval_cox_model)) {
    return(NULL)
  }
  cox_summary <- summary(eval_cox_model)
  cox_coefficients <- coef(eval_cox_model)
  cox_p_values <- cox_summary$coefficients[, "Pr(>|z|)"]
  # Select features
  significant_features <- names(cox_coefficients[cox_p_values < 0.05 & cox_coefficients != 0])
  significant_cox_df <- data.frame(feature = significant_features)
  return(significant_cox_df)
}


rsfvh <- function(data) {
  # Exclude survival columns from predictors
  predictor_data <- data[, !(colnames(data) %in% c("overall_survival", "deceased"))]
  # Construct formula explicitly
  f <- as.formula("Surv(overall_survival, deceased) ~ .")
  # Run var.select on data that includes survival columns plus only predictor columns
  rsfvh_model <- var.select(f,
                            data = cbind(data[, c("overall_survival", "deceased")], predictor_data),
                            ntree = 1000,
                            method = "vh",
                            nodesize = 5,
                            nsplit = 20,
                            splitrule = "logrank",
                            nrep = 3,
                            K = 10,
                            nstep = 1)
  rsfvhfeature_Selection <- data.frame(feature = unlist(rsfvh_model$topvars))
  return(rsfvhfeature_Selection)
}


rsfmd <- function(data) {
  # Exclude survival columns from predictors
  predictor_data <- data[, !(colnames(data) %in% c("overall_survival", "deceased"))]
  # Combine survival columns and predictors explicitly
  model_data <- cbind(data[, c("overall_survival", "deceased")], predictor_data)
  # Define formula
  f <- as.formula("Surv(overall_survival, deceased) ~ .")
  # Run variable selection with method 'md'
  rsfmd_model <- var.select(f,
                            data = model_data,
                            ntree = 1000,
                            method = "md",
                            nodesize = 5,
                            nsplit = 20,
                            splitrule = "logrank")
  rsfmdfeature_Selection <- data.frame(feature = unlist(rsfmd_model$topvars))
  return(rsfmdfeature_Selection)
}


rsfvi <- function(data) {
  # Exclude survival columns from predictors
  predictor_data <- data[, !(colnames(data) %in% c("overall_survival", "deceased"))]
  # Combine survival columns and predictors explicitly
  model_data <- cbind(data[, c("overall_survival", "deceased")], predictor_data)
  # Set mtry as sqrt of number of predictor variables (excluding survival)
  mtry_value <- floor(sqrt(ncol(predictor_data)))
  # Fit the random survival forest model
  rsf_model <- rfsrc(Surv(overall_survival, deceased) ~ ., 
                     data = model_data,
                     ntree = 1000,
                     mtry = mtry_value,
                     nodesize = 3,
                     nsplit = 10)
  # Extract variable importance
  var_importance_rsf <- vimp(rsf_model)
  var_importance_df <- data.frame(
    feature = names(var_importance_rsf$importance),
    Importance = var_importance_rsf$importance
  )
  # Select features with positive importance
  rsffeature_selection <- var_importance_df[var_importance_df$Importance > 0, ]
  return(rsffeature_selection)
}


all_features <- function(data) {
  # Remove sampleID and survival columns
  feature_cols <- colnames(data)[-(1:2)]
  
  all_features_df <- data.frame(feature = feature_cols, stringsAsFactors = FALSE)
  return(all_features_df)
}


perform_and_save_results <- function(data, survival_data, repeats, folds, custom_name) {
  
  # Perform feature selection and combine results
  results <- perform_parallel_feature_selection(data, survival_data, repeats, folds, all_features)
  no_fs_results <- combine_results(results)
  results <- perform_parallel_feature_selection(data, survival_data, repeats, folds, multivariate)
  multi_results <- combine_results(results)
  results <- perform_parallel_feature_selection(data, survival_data, repeats, folds, univariate)
  uni_results <- combine_results(results)
  results <- perform_parallel_feature_selection(data, survival_data, repeats, folds, rsfvi)
  rsf_vi_results <- combine_results(results)
  results <- perform_parallel_feature_selection(data, survival_data, repeats, folds, rsfmd)
  rsf_md_results <- combine_results(results)
  results <- perform_parallel_feature_selection(data, survival_data, repeats, folds, rsfvh)
  rsf_vh_results <- combine_results(results)
  # Plot the results
  plot_results(no_fs_results,multi_results, uni_results, rsf_vi_results, rsf_md_results, rsf_vh_results, custom_name)
}


perform_feature_selection_CV <- function(data,feature_counts_df, output_dir) {
  
  # Perform feature selection and update counts
  updated_feature_counts_df <- perform_feature_parallel_CV(data, multivariate, feature_counts_df)
  updated_feature_counts_df <- perform_feature_parallel_CV(data, rsfvi, updated_feature_counts_df)
  updated_feature_counts_df <- perform_feature_parallel_CV(data, rsfmd, updated_feature_counts_df)
  updated_feature_counts_df <- perform_feature_parallel_CV(data, rsfvh, updated_feature_counts_df)
  
  # Sort the feature_counts_df by Count column in descending order
  sorted_feature_counts <- updated_feature_counts_df[order(-updated_feature_counts_df$Count), ]
  # Filter features based on the threshold
  sorted_feature_counts <- sorted_feature_counts %>%
    filter(Count >= 0.5 * 100)
  
  # Define the directory and file path
  dir_path <- output_dir
  file_path <- file.path(dir_path, "features_cv.csv")
  
  # Create the directory if it doesn't exist
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  # Write the data frame to a CSV file
  write.csv(sorted_feature_counts, file = file_path, row.names = FALSE)
  
  return(sorted_feature_counts)
}

perform_feature_parallel_CV <- function(data, feature_selection, feature_counts_df, repeats = 5, folds = 5) {
  
  # Identify survival columns (adjust if different)
  time_col <- "overall_survival"
  event_col <- "deceased"
  
  num_cores <- 16 # Assume 16 cores
  cl <- makeCluster(num_cores)
  # Export R library path (for HPC systems like Katana)
  clusterEvalQ(cl, {
    .libPaths("/path/to/Rlib")
  })
  registerDoParallel(cl)
  on.exit(stopCluster(cl))  # Ensure cluster is stopped even if an error occurs
  
  # Precompute folds per repeat
  folds_list <- lapply(1:repeats, function(rep) {
    createFolds(factor(data[[event_col]]), k = folds)
  })
  
  all_indices <- expand.grid(rep = 1:repeats, fold = 1:folds)
  
  all_results <- foreach(i = 1:nrow(all_indices), .combine = 'list', .multicombine = TRUE,
                         .packages = c('caret', 'survival', 'glmnet', 'mboost', 'ranger', 'randomForestSRC'),
                         .export = c("feature_selection")) %dopar% {
                           
                           rep <- all_indices$rep[i]
                           fold <- all_indices$fold[i]
                           
                           fold_indices <- folds_list[[rep]]
                           train_indices <- unlist(fold_indices[-fold])
                           
                           train_data <- data[train_indices, ]
                           
                           train_features <- train_data[, !(colnames(train_data) %in% c("sampleID"))]
                           feature_selected <- feature_selection(train_features)
                           feature_selected$feature
                         }
  
  selected_features <- unlist(all_results)  
  feature_counts <- table(selected_features)
  
  # Update the passed-in feature_counts_df
  for (feature in names(feature_counts)) {
    row_index <- which(feature_counts_df$Feature == feature)
    if (length(row_index) == 1) {
      feature_counts_df$Count[row_index] <- feature_counts_df$Count[row_index] + feature_counts[[feature]]
    }
  }
  
  return(feature_counts_df)
}
                     
perform_bootstrapping <- function(data, n_bootstrap = 100) {
  # Define a function for feature selection
  feature_selection <- function(data) {
    multi_features <- multivariate(data)
    rsfvi_features <- rsfvi(data)
    rsfvh_features <- rsfvh(data)
    rsfmd_features <- rsfmd(data)
    
    common_values <- Reduce(intersect, list(
      multi_features$feature,
      rsfvi_features$feature,
      rsfvh_features$feature,
      rsfmd_features$feature
    ))
    return(common_values)
  }
  # Create a cluster with all available cores
  num_cores <- 16 # Assume 16 cores
  cl <- makeCluster(num_cores)
  clusterEvalQ(cl, {
    .libPaths("/path/to/Rlib")
  })
  # Register the cluster
  registerDoParallel(cl)
  selected_features_list <- foreach(i = 1:n_bootstrap,
                                    .combine = 'list',
                                    .multicombine = TRUE,
                                    .packages = c("survival", "randomForestSRC", "glmnet", "mboost", "ranger"),
                                    .export = c("multivariate", "rsfvi", "rsfvh", "rsfmd", "feature_selection")) %dopar% {
                                      # Bootstrap 70% of the data with replacement
                                      boot_data <- data[sample(nrow(data), size = 0.7 * nrow(data), replace = TRUE), ]
                                      
                                      # Run survival-aware feature selection
                                      selected_features <- feature_selection(boot_data)
                                      return(selected_features)
                                    }
  
  stopCluster(cl)
  
  # Count feature occurrences across all bootstraps
  feature_table <- table(unlist(selected_features_list))
  sorted_feature_df <- data.frame(
    feature = names(feature_table),
    count = as.integer(feature_table)
  )
  sorted_feature_df <- sorted_feature_df[order(-sorted_feature_df$count), ]
  
  return(sorted_feature_df)
}
