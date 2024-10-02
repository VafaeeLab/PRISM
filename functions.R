# Impute missing values with row mean
impute_row_mean <- function(row) {
  row_mean <- mean(row, na.rm = TRUE)
  row[is.na(row)] <- row_mean
  return(row)
}
generate_train_test <- function(data, survival_data, train_ratio = 0.7) {
  train_rows <- round(train_ratio * nrow(data))
  train_indices <- sample(seq_len(nrow(data)), size = train_rows)
  train_data <- data[train_indices, ]
  train_x_matrix <- as.matrix(train_data[, -c(1:2)])
  test_indices <- setdiff(seq_len(nrow(data)), train_indices)
  test_data <- data[test_indices, ]
  test_x_matrix <- as.matrix(test_data[, -c(1:2)])
  
  train_df <- as.data.frame(train_x_matrix)
  test_df <- as.data.frame(test_x_matrix)
  
  train_survival <- survival_data[train_indices]
  test_survival <- survival_data[test_indices]
  
  return(list(train_df = train_df,
              test_df = test_df,
              train_survival = train_survival,
              test_survival = test_survival))
}

# Recursive function for feature elimination
perform_feature_elimination <- function(train_data, test_data, train_survival_data, test_survival_data, count, feature_c_index_list) {
  
  while (ncol(train_data) >= 1 ) {
    print(ncol(train_data))
    
    # Random Forest (Ranger)
    mtry_value <- floor(sqrt(ncol(train_data)))
    eval_r_fit <- ranger(
      formula = train_survival_data ~ .,
      data = train_data,
      mtry = mtry_value,
      importance = "permutation",
      splitrule = "maxstat",
      min.node.size = 50,
      num.trees = 1000
    )
    
    # Store feature importance scores
    selected_features_ranger_df <- data.frame(
      Variable = names(eval_r_fit$variable.importance),
      Importance = eval_r_fit$variable.importance
    )
    
    # Rank features based on importance scores
    selected_features_ranger_df <- selected_features_ranger_df %>%
      arrange(desc(Importance))
    
    ranger_predictions <- predict(eval_r_fit, data = test_data)$survival
    RFE_c_index_ranger <- concordance(test_survival_data ~ ranger_predictions)
    # Check if there are multiple c-index values
    if (length(RFE_c_index_ranger$concordance) > 1) {
      # Compute the mean of multiple c-index values
      mean_c_index <- mean(RFE_c_index_ranger$concordance)
    } else {
      # Only one c-index value, do nothing
      mean_c_index <- RFE_c_index_ranger$concordance
    }
    
    feature_c_index_list[[count]] <- list(variables = selected_features_ranger_df$Variable, c_index = mean_c_index)
    count <- count + 1
    
    if (ncol(train_data) == 1) {
      break
    }
    # Remove the least important feature
    least_important_feature <- tail(selected_features_ranger_df$Variable, 1)
    # Find the index of the least important feature in column names
    feature_index <- which(colnames(train_data) == least_important_feature)
    train_data <- train_data[, -feature_index, drop = FALSE]
    test_data <- test_data[, -feature_index, drop = FALSE]
    
  }
  
  return(feature_c_index_list)
}


get_max_c_index_index <- function(result_list) {
  # Initialize variables to store the maximum c_index value and its corresponding index
  max_c_index <- -Inf
  max_index <- NULL
  
  # Iterate through each element in the list
  for (i in seq_along(result_list)) {
    # Extract the c_index value from the current element
    current_c_index <- result_list[[i]]$c_index
    
    # Check if the current c_index value is greater than the current maximum
    if (current_c_index > max_c_index) {
      max_c_index <- current_c_index
      max_index <- i
    }
  }
  
  # Return the index with the highest c_index value
  return(max_index)
}

# Function to perform min-max normalization ignoring the first three columns
min_max_normalize <- function(x) {
  if (ncol(x) > 3) {
    x[, -c(1:3)] <- apply(x[, -c(1:3)], 2, function(y) (y - min(y)) / (max(y) - min(y)))
  }
  return(x)
}


# Main function to process data and save plot and dataframe
fusion <- function(data, survival_data, directory = ".", file_prefix = "output") {
  
  data <- data[,-1]
  # Generate train and test data
  train_test_data <- generate_train_test(data, survival_data)
  train_df <- train_test_data$train_df
  test_df <- train_test_data$test_df
  train_survival <- train_test_data$train_survival
  test_survival <- train_test_data$test_survival
  
  feature_c_index_list <- list()
  count <- 1
  
  # Perform feature elimination
  result <- perform_feature_elimination(train_df, test_df, train_survival, test_survival, count, feature_c_index_list)
  
  # Extract the number of features and c-index values
  num_features <- sapply(result, function(x) length(x$variables))
  c_index_values <- sapply(result, function(x) x$c_index)
  
  # Create dataframe for plotting
  result_df <- data.frame(num_features = num_features, c_index_values = c_index_values)
  
  # Create plot
  result_plot <- ggplot(result_df, aes(x = num_features, y = c_index_values)) +
    geom_point() +   # Scatter plot
    geom_line() +    # Connect points with a line
    labs(x = "Number of Features", y = "C-Index", title = paste(file_prefix, "Performance"))
  
  # Ensure directory exists
  if (!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE)
  }
  
  # Get best results
  max_index <- get_max_c_index_index(result)

  # Save plot
  plot_file_path <- file.path(directory, paste0(file_prefix, "_plot.png"))
  ggsave(plot_file_path, plot = result_plot, width = 10, height = 6, units = "in")
  
  # Save dataframe
  data_file_path <- file.path(directory, paste0(file_prefix, "_data.csv"))
  write.csv(result[[max_index]], data_file_path, row.names = FALSE)
}

evaluate_performance <- function(train_df, train_survival, test_df, test_survival) {
  
  train_x_matrix <- as.matrix(train_df)
  test_x_matrix <- as.matrix(test_df)
  # Cox Model
  eval_cox_model <- coxph(train_survival ~ ., data = train_df)
  cox_predictions <- predict(eval_cox_model, newdata = test_df)
  c_index_cox <- concordance(test_survival ~ cox_predictions)
  c_index_cox_results <- c_index_cox$concordance
  
  # Random Forest (Ranger)
  mtry_value <- floor(sqrt(ncol(train_df)))
  eval_r_fit <- ranger(
    formula = train_survival ~ .,
    data = train_df,
    mtry = mtry_value,
    importance = "permutation",
    splitrule = "maxstat",
    min.node.size = 50,
    num.trees = 1000
  )
  
  ranger_predictions <- predict(eval_r_fit, data = test_df)$survival
  c_index_ranger <- concordance(test_survival ~ ranger_predictions)
  # Check if there are multiple c-index values
  if (length(c_index_ranger$concordance) > 1) {
    # Compute the mean of multiple c-index values
    mean_c_index <- mean(c_index_ranger$concordance)
  } else {
    # Only one c-index value, do nothing
    mean_c_index <- c_index_ranger$concordance
  }
  c_index_ranger_results <- mean_c_index
  
  # Boosted Cox Model (GLMBoost)
  eval_boosted_cox_model <- glmboost(train_survival ~ ., data = train_df, family = CoxPH())
  glmboost_predictions <- predict(eval_boosted_cox_model, newdata = test_df)
  c_index_glmboost <- concordance(test_survival ~ glmboost_predictions)
  c_index_glmboost_results <- c_index_glmboost$concordance
  
  # Elastic Net Cox Model
  elastic_net_model <- cv.glmnet(x = train_x_matrix, y = train_survival, family = "cox", alpha = 0.5)
  predictions <- predict(elastic_net_model, newx = test_x_matrix, s = "lambda.min")
  c_index <- concordance(test_survival ~ predictions)
  c_index_elasticnet_results <- c_index$concordance
  
  cindex_df <- data.frame(
    CIndexCox = c_index_cox_results,
    CIndexRanger = c_index_ranger_results,
    CIndexGLMBoost = c_index_glmboost_results,
    CIndexElasticNet = c_index_elasticnet_results
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
  
  features_selected <- results$selected_features_list
  
  # Combine results into data frames
  feature_selection_df <- data.frame(
    FeatureCox = feature_cox_results,
    FeatureRanger = feature_ranger_results,
    FeatureGLMBoost = feature_glmboost_results,
    FeatureElasticNet = feature_elasticnet_results
  )
  
  feature_cindex_df <- data.frame(
    CIndexCox = c_index_cox_results,
    CIndexRanger = c_index_ranger_results,
    CIndexGLMBoost = c_index_glmboost_results,
    CIndexElasticNet = c_index_elasticnet_results
  )
  
  # Calculate standard deviations
  cox_feature_sd <- sd(feature_cox_results)
  boost_feature_sd <- sd(feature_glmboost_results)
  ranger_feature_sd <- sd(feature_ranger_results)
  penalised_feature_sd <- sd(feature_elasticnet_results)
  
  cox_Cindex_sd <- sd(c_index_cox_results)
  boost_Cindex_sd <- sd(c_index_glmboost_results)
  ranger_Cindex_sd <- sd(c_index_ranger_results)
  penalised_Cindex_sd <- sd(c_index_elasticnet_results)
  
  # Calculate mean values
  mean_feature_values <- apply(feature_selection_df, 2, mean)
  mean_cindex_values <- apply(feature_cindex_df, 2, mean)
  
  # Create summary data frames
  feature_selected_result_summary <- data.frame(
    Model = names(mean_feature_values),
    MeanFeatureSelected = mean_feature_values,
    StandardDeviation = c(cox_feature_sd, ranger_feature_sd, boost_feature_sd, penalised_feature_sd)
  )
  
  feature_cindex_result_summary <- data.frame(
    Model = names(mean_cindex_values),
    MeanCIndex = mean_cindex_values,
    StandardDeviation = c(cox_Cindex_sd, ranger_Cindex_sd, boost_Cindex_sd, penalised_Cindex_sd)
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
  
  # Merge results
  no_fs_feature_merged <- cbind(no_fs_feature_cindex_result_summary, no_fs_feature_selected_result_summary)
  multi_feature_merged <- cbind(multi_feature_cindex_result_summary, multi_feature_selected_result_summary)
  uni_feature_merged <- cbind(uni_feature_cindex_result_summary, uni_feature_selected_result_summary)
  rsf_feature_merged <- cbind(rsf_cindex_result_summary, rsf_feature_selected_result_summary)
  rsfmd_feature_merged <- cbind(rsfmd_cindex_result_summary, rsfmd_feature_selected_result_summary)
  rsfvh_feature_merged <- cbind(rsfvh_cindex_result_summary, rsfvh_feature_selected_result_summary)
  
  # Save merged results to CSV files
  write.csv(no_fs_feature_merged, file.path(dir_name, paste0(custom_name, "_no_fs_results.csv")), row.names = FALSE)
  write.csv(multi_feature_merged, file.path(dir_name, paste0(custom_name, "_multi_results.csv")), row.names = FALSE)
  write.csv(uni_feature_merged, file.path(dir_name, paste0(custom_name, "_uni_results.csv")), row.names = FALSE)
  write.csv(rsf_feature_merged, file.path(dir_name, paste0(custom_name, "_rsf_results.csv")), row.names = FALSE)
  write.csv(rsfmd_feature_merged, file.path(dir_name, paste0(custom_name, "_rsfmd_results.csv")), row.names = FALSE)
  write.csv(rsfvh_feature_merged, file.path(dir_name, paste0(custom_name, "_rsfvh_results.csv")), row.names = FALSE)
  
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
  
  ggsave(file.path(dir_name, paste0(custom_name, "_mean_performance.png")), gg, width = 8, height = 6)
  
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
  
  ggsave(file.path(dir_name, paste0(custom_name, "_mean_features.png")), gg2, width = 8, height = 6)
  
  # Prepare data for the boxplots
  boxplot_data <- list(
    None = no_fs_feature_cindex_result_summary$MeanCIndex, 
    Multivariate = multi_feature_cindex_result_summary$MeanCIndex, 
    RFVI = rsf_cindex_result_summary$MeanCIndex,
    Univariate = uni_feature_cindex_result_summary$MeanCIndex,
    RFMD = rsfmd_cindex_result_summary$MeanCIndex,
    RFMS = rsfvh_cindex_result_summary$MeanCIndex
  )
  
  boxplot_df <- stack(boxplot_data)
  colnames(boxplot_df) <- c("CIndex", "Model")
  
  # Boxplot using ggplot2
  gg3 <- ggplot(boxplot_df, aes(x = Model, y = CIndex, fill = Model)) + 
    geom_boxplot() +
    scale_fill_manual(values = c("#a6cee3", "#fdbf6f","#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3")) +
    theme_minimal() +
    labs(title = "Performance of each FS measured by the predictive models", x = "Model", y = "C-Index") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(hjust = 1),
          axis.title = element_text(size = 10),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          panel.grid = element_blank(),
          plot.margin = margin(30, 30, 30, 30, unit = "pt"))
  
  ggsave(file.path(dir_name, paste0(custom_name, "_cindex_boxplot.png")), gg3, width = 8, height = 6)
}

perform_feature_parallel_selection <- function(data, survival_data, repeats, folds, feature_selection) {
  
  # Set up parallel backend
  num_cores <- 16
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  # Perform feature selection in parallel
  results <- foreach(rep = 1:repeats, .combine = 'c', .packages = c('caret', 'survival', 'randomForestSRC', 'glmnet', 'ranger', 'mboost')) %:%
    foreach(fold = 1:folds, .combine = 'c') %dopar% {
      
      fold_indices <- createFolds(survival_data, k = folds)
      train_indices <- unlist(fold_indices[-fold])
      test_indices <- unlist(fold_indices[fold])
      
      train_data <- data[train_indices, ]
      train_survival <- survival_data[train_indices]
      train_x_matrix <- as.matrix(train_data[, -1])
      
      test_data <- data[test_indices, ]
      test_survival <- survival_data[test_indices]
      test_x_matrix <- as.matrix(test_data[, -1])
      
      train_df <- as.data.frame(train_x_matrix)
      test_df <- as.data.frame(test_x_matrix)
      
      feature_selected <- feature_selection(train_df)
      
      train_x_matrix <- as.matrix(train_x_matrix[,feature_selected$feature])
      test_x_matrix <- as.matrix(test_x_matrix[, feature_selected$feature])
      train_df <- as.data.frame(train_df[,feature_selected$feature])
      test_df <- as.data.frame(test_df[,feature_selected$feature])
      
      # Cox Model
      eval_cox_model <- coxph(train_survival ~ ., data = train_df)
      cox_predictions <- predict(eval_cox_model, newdata = test_df)
      c_index_cox <- concordance(test_survival ~ cox_predictions)$concordance
      selected_features_cox <- names(coef(eval_cox_model))[coef(eval_cox_model) != 0]
      
      # Random Forest (Ranger)
      mtry_value <- floor(sqrt(ncol(train_df)))
      eval_r_fit <- ranger(
        formula = train_survival ~ .,
        data = train_df,
        mtry = mtry_value,
        importance = "permutation",
        splitrule = "maxstat",
        min.node.size = 50,
        num.trees = 1000
      )
      ranger_predictions <- predict(eval_r_fit, data = test_df)$survival
      c_index_ranger_results <- concordance(test_survival ~ ranger_predictions)
      
      # Check if there are multiple c-index values
      if (length(c_index_ranger_results$concordance) > 1) {
        # Compute the mean of multiple c-index values
        mean_c_index <- mean(c_index_ranger_results$concordance)
      } else {
        # Only one c-index value, do nothing
        mean_c_index <- c_index_ranger_results$concordance
      }
      c_index_ranger <- mean_c_index
      
      selected_features_ranger_df <- data.frame(
        Variable = names(eval_r_fit$variable.importance),
        Importance = eval_r_fit$variable.importance
      )
      selected_features_ranger <- selected_features_ranger_df[selected_features_ranger_df$Importance > 0, ]
      
      # Boosted Cox Model (GLMBoost)
      eval_boosted_cox_model <- glmboost(train_survival ~ ., data = train_df, family = CoxPH())
      glmboost_predictions <- predict(eval_boosted_cox_model, newdata = test_df)
      c_index_glmboost <- concordance(test_survival ~ glmboost_predictions)$concordance
      selected_features_boosted_cox <- names(coef(eval_boosted_cox_model))[coef(eval_boosted_cox_model) != 0]
      
      # Elastic Net Cox Model
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
      selected_features_elasticnet <-  selected_features_elasticnet_df[selected_features_elasticnet_df$coefficient != 0, ]
      
      list(
        selected_features = feature_selected$feature,
        c_index_elasticnet = c_index,
        c_index_cox = c_index_cox,
        c_index_ranger = c_index_ranger,
        c_index_glmboost = c_index_glmboost,
        feature_elasticnet = nrow(selected_features_elasticnet),
        feature_cox = length(selected_features_cox),
        feature_ranger = nrow(selected_features_ranger),
        feature_glmboost = length(selected_features_boosted_cox) - 1
      )
    }
  
  # Stop the parallel backend
  stopCluster(cl)
  
  # Combine results from parallel execution
  results <- matrix(results, nrow = repeats * folds, byrow = TRUE)
  
  selected_features_list <- lapply(results[, 1], unlist)
  c_index_elasticnet_results <- unlist(results[, 2])
  c_index_cox_results <- unlist(results[, 3])
  c_index_ranger_results <- unlist(results[, 4])
  c_index_glmboost_results <- unlist(results[, 5])
  feature_elasticnet_results <- unlist(results[, 6])
  feature_cox_results <- unlist(results[, 7])
  feature_ranger_results <- unlist(results[, 8])
  feature_glmboost_results <- unlist(results[, 9])
  
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

# Univariate --------------------------------------------------------------------
univariate <- function(data) {
  uni_cox_df <- data.frame(feature = character(),
                           coefficients = numeric(),
                           p_value = numeric(),
                           c_index = numeric(),
                           stringsAsFactors = FALSE)
  
  for (j in 3:ncol(data)) {
    feature <- colnames(data)[j]
    cox <- try(coxph(Surv(overall_survival, deceased) ~ get(feature), data = data))
    
    if (inherits(cox, "try-error")) {
      next  
    }
    
    p_value <- summary(cox)$logtest["pvalue"] # likelihood ratio test
    c_index <- summary(cox)$concordance["C"]  # C-index
    coefficients <- coef(cox) # Coefficient
    result_row <- data.frame(feature = feature, coefficients = coefficients, p_value = p_value, c_index = c_index)
    
    uni_cox_df <- rbind(uni_cox_df, result_row)
  }
  
  significant_uni_cox_df <- subset(uni_cox_df, p_value < 0.05 & coefficients != 0)
  return(significant_uni_cox_df)
}

# Multi-variate --------------------------------------------------------------------
multivariate <- function(data) {
  # Cox Model
  eval_cox_model <- coxph(Surv(overall_survival, deceased) ~ ., data = data)
  # Extract coefficients and corresponding p-values
  cox_summary <- summary(eval_cox_model)
  cox_coefficients <- coef(eval_cox_model)
  cox_p_values <- cox_summary$coefficients[, "Pr(>|z|)"]
  
  # Extract significant features based on p-values < 0.05 and non-zero coefficients
  significant_features <- names(cox_coefficients[cox_p_values < 0.05 & cox_coefficients != 0])
  
  # Create a data frame of significant features
  significant_cox_df <- data.frame(feature = significant_features)
  
  return(significant_cox_df)
}


# Random Forest Variable Hunting ---------------------------------------------
rsfvh <- function(data) {
  rsfvh_model <- var.select(Surv(overall_survival, deceased) ~ ., 
                            data = data,
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

# Random Forest Min Depth ---------------------------------------------
rsfmd <- function(data) {
  rsfmd_model <- var.select(Surv(overall_survival, deceased) ~ ., 
                            data = data,
                            ntree = 1000,
                            method = "md",
                            nodesize = 5,
                            nsplit = 20,
                            splitrule = "logrank")
  
  rsfmdfeature_Selection <- data.frame(feature = unlist(rsfmd_model$topvars))
  
  return(rsfmdfeature_Selection)
}

# Random Forest Variable Importance ---------------------------------------------
rsfvi <- function(data) {
  mtry_value <- floor(sqrt(ncol(data)))
  rsf_model <- rfsrc(Surv(overall_survival, deceased) ~ ., 
                     data = data,
                     ntree = 1000,
                     mtry = mtry_value,
                     nodesize = 3,
                     nsplit = 10)
  
  var_importance_rsf <- vimp(rsf_model)
  
  var_importance_df <- data.frame(
    feature = names(var_importance_rsf$importance),
    Importance = var_importance_rsf$importance
  )
  
  rsffeature_selection <- var_importance_df[var_importance_df$Importance > 0, ]
  return(rsffeature_selection)
}


all_features <- function(data) {
  all_features <- colnames(data)[-c(1,2)]
  all_features_df <- data.frame(features = all_features)
  return(all_features_df)
}

perform_and_save_results <- function(data, survival_data, repeats, folds, custom_name) {
  
  # Perform feature selection and combine results
  results <- perform_feature_parallel_selection(data, survival_data, repeats, folds, all_features)
  no_fs_results <- combine_results(results)
  results <- perform_feature_parallel_selection(data, survival_data, repeats, folds, multivariate)
  multi_results <- combine_results(results)
  results <- perform_feature_parallel_selection(data, survival_data, repeats, folds, univariate)
  uni_results <- combine_results(results)
  results <- perform_feature_parallel_selection(data, survival_data, repeats, folds, rsfvi)
  rsf_vi_results <- combine_results(results)
  results <- perform_feature_parallel_selection(data, survival_data, repeats, folds, rsfmd)
  rsf_md_results <- combine_results(results)
  results <- perform_feature_parallel_selection(data, survival_data, repeats, folds, rsfvh)
  rsf_vh_results <- combine_results(results)
  # Plot the results
  plot_results(no_fs_results,multi_results, uni_results, rsf_vi_results, rsf_md_results, rsf_vh_results, custom_name)
}

perform_feature_selection_CV <- function(data, survival_data,feature_counts_df, output_dir) {
  # Perform feature selection and update counts
  updated_feature_counts_df <- perform_feature_parallel_CV(data, survival_data, multivariate, feature_counts_df)
  updated_feature_counts_df <- perform_feature_parallel_CV(data, survival_data, rsfvi, updated_feature_counts_df)
  updated_feature_counts_df <- perform_feature_parallel_CV(data, survival_data, rsfmd, updated_feature_counts_df)
  updated_feature_counts_df <- perform_feature_parallel_CV(data, survival_data, rsfvh, updated_feature_counts_df)
  
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

perform_feature_parallel_CV <- function(data, survival_data, feature_selection, feature_counts_df, repeats = 5, folds = 5) {
  
  # Set up parallel backend
  num_cores <- 16
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  # Perform feature selection in parallel
  results <- foreach(rep = 1:repeats, .combine = 'list', .packages = c('caret', 'survival', 'randomForestSRC', 'glmnet', 'ranger', 'mboost')) %:%
    foreach(fold = 1:folds, .combine = 'list') %dopar% {
      
      fold_indices <- createFolds(survival_data, k = folds)
      train_indices <- unlist(fold_indices[-fold])
      train_data <- data[train_indices, ]
      train_survival <- survival_data[train_indices]
      train_x_matrix <- as.matrix(train_data)
      train_df <- as.data.frame(train_x_matrix)
      
      feature_selected <- feature_selection(train_df)
      selected_features_list <- list()
      selected_features_list <- c(selected_features_list, feature_selected$feature)
      return(selected_features_list)
    }
  
  # Stop the parallel backend
  stopCluster(cl)
  
  # Flatten the list of results and update feature_counts_df
  all_selected_features <- unlist(results)
  feature_counts <- table(all_selected_features)
  
  for (feature in names(feature_counts)) {
    row_index <- which(feature_counts_df$Feature == feature)
    if (length(row_index) > 0) {
      feature_counts_df$Count[row_index] <- feature_counts_df$Count[row_index] + feature_counts[feature]
    } 
  }
  
  return(feature_counts_df)
}


perform_bootstrapping <- function(data) {
  # Define a function for feature selection
  feature_selection <- function(data, multivariate, rsfvi, rsfvh, rsfmd) {
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
  cl <- makeCluster(16)
  # Register the cluster
  registerDoParallel(cl)
  
  # Bootstrap and run feature selection 100 times in parallel
  selected_features_list <- foreach(i = 1:100, 
                                    .packages = c("survival", "randomForestSRC"),
                                    .export = c("multivariate", "rsfvi", "rsfvh", "rsfmd", "feature_selection")) %dopar% {
                                      # Bootstrap 70% of the data
                                      boot_data <- data[sample(nrow(data), size = 0.7 * nrow(data), replace = TRUE), ]
                                      
                                      # Run feature selection
                                      selected_features <- feature_selection(boot_data, multivariate, rsfvi, rsfvh, rsfmd)
                                      
                                      # Return selected features
                                      selected_features
                                    }
  
  # Stop the cluster
  stopCluster(cl)
  
  # Initialize a list to store the count of each feature
  feature_count <- list()
  
  # Count the occurrence of each feature across all selected features lists
  for (i in 1:length(selected_features_list)) {
    for (feature in selected_features_list[[i]]) {
      if (feature %in% names(feature_count)) {
        feature_count[[feature]] <- feature_count[[feature]] + 1
      } else {
        feature_count[[feature]] <- 1
      }
    }
  }
  
  # Sort feature_count based on count
  sorted_feature_count <- feature_count[order(unlist(feature_count), decreasing = TRUE)]
  
  # Convert the sorted_feature_count list to a dataframe
  sorted_feature_count_df <- data.frame(
    feature = names(sorted_feature_count),
    count = unlist(sorted_feature_count)
  )
  
  # Return the selected features dataframe
  return(sorted_feature_count_df)
}

# Define the categorization function
categorize_cnv <- function(copy_number) {
  case_when(
    copy_number < 0.5 ~ -2,       # Homozygous deletion
    copy_number >= 0.5 & copy_number < 1.5 ~ -1,  # Hemizygous deletion
    copy_number >= 1.5 & copy_number <= 2.5 ~ 0,  # Neutral / No change
    copy_number > 2.5 & copy_number <= 4 ~ 1,     # Low-level amplification
    copy_number > 4 ~ 2,           # High-level amplification
    TRUE ~ NA_real_  # Handle missing or unexpected values
  )
}

# Example function to create KM plots for each gene
create_cnv_km_plots <- function(data, gene, cancer) {
  # Select relevant columns including the specific gene column
  data <- data[, c("case_id", "overall_survival", "deceased", gene)]
  
  # Ensure the specific gene column is numeric
  data[[gene]] <- as.numeric(data[[gene]])
  
  # Categorize CNV into Gain, Loss, and Neutral
  data$strata <- ifelse(data[[gene]] > 0, "Gain",
                        ifelse(data[[gene]] < 0, "Loss", "Neutral"))
  
  # Create a Surv object
  surv_object <- Surv(time = data$overall_survival, event = data$deceased)
  # Fit survival model
  fit <- survfit(surv_object ~ strata, data = data)
  
  # Plot Kaplan-Meier curve
  ggsurv <- ggsurvplot(fit,
                       data = data,
                       risk.table = TRUE,
                       pval = TRUE,
                       title = paste("Kaplan-Meier Curve", cancer, "for Gene:", gene),
                       xlab = "Time (days)",
                       ylab = "Survival Probability",
                       palette = "Dark2")
  
  # Print the plot explicitly
  print(ggsurv)
}


create_ge_km_plots <- function(data, gene, cancer) {
  # Select relevant columns including the specific gene column
  data <- data[, c("case_id", "overall_survival", "deceased", gene)]
  
  # Get median value
  column_value <- data[[gene]]
  median_value <- median(column_value, na.rm = TRUE)
  
  # Create strata based on the median value
  data <- data %>%
    mutate(strata = case_when(
      data[[gene]] >= median_value ~ "High",
      data[[gene]] < median_value ~ "Low",
      TRUE ~ NA_character_
    ))
  
  # Remove NA strata
  combined_data <- na.omit(data)
  
  # Create a Surv object
  surv_object <- Surv(time = combined_data$overall_survival, event = combined_data$deceased)
  
  # Fit survival model
  fit <- survfit(surv_object ~ strata, data = combined_data)
  
  # Plot Kaplan-Meier curve
  ggsurv <- ggsurvplot(fit,
                       data = combined_data,
                       risk.table = TRUE,
                       pval = TRUE,
                       title = paste("Kaplan-Meier Curve", cancer, "for miRNA:", gene),
                       xlab = "Time (days)",
                       ylab = "Survival Probability",
                       palette = "Dark2")
  
  # Print the plot explicitly
  print(ggsurv)
}