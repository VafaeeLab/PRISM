# ============================================================
# Function: impute_row_mean
# ============================================================

# Description:
# This function imputes missing values in a row of data by replacing them with the mean of the non-missing values in the same row.

# Parameters:
# - row: A numeric vector (representing a single row of data) with potentially missing values (NA).

# Returns:
# - A numeric vector with missing values replaced by the mean of the non-missing values in the row.

# Procedure:
# 1. The function calculates the mean of the non-missing values in the input row using `mean()` with `na.rm = TRUE`.
# 2. It then replaces all missing values (`NA`) in the row with the calculated mean.
# 3. The function returns the row with imputed values.
impute_row_mean <- function(row) {
  row_mean <- mean(row, na.rm = TRUE)
  row[is.na(row)] <- row_mean
  return(row)
}
# ============================================================
# Function: generate_train_test
# ============================================================

# Description:
# This function splits the input data into training and testing sets based on a specified ratio. 
# It also separates the survival data corresponding to each set. The training data will be used to train models, 
# and the testing data will be used to evaluate model performance.

# Parameters:
# - data: A data frame containing the feature data. The first two columns are assumed to be non-feature columns 
#         (e.g., survival time and event indicator).
# - survival_data: A vector containing survival data corresponding to each row in the `data` frame.
# - train_ratio: A numeric value between 0 and 1, specifying the proportion of data to be used for training (default is 0.7, i.e., 70%).

# Returns:
# - A list with the following components:
#   - `train_df`: A data frame containing the features for the training set.
#   - `test_df`: A data frame containing the features for the test set.
#   - `train_survival`: A vector containing the survival data for the training set.
#   - `test_survival`: A vector containing the survival data for the test set.

# Procedure:
# 1. The function computes the number of rows to include in the training set based on the `train_ratio`.
# 2. It randomly selects rows for the training set using a random sampling method.
# 3. The test set is derived by excluding the training set rows from the original data.
# 4. The training and testing sets are separated into `train_df` and `test_df`, respectively.
# 5. The corresponding survival data for both the training and testing sets is extracted.
# 6. The function returns a list containing the training and testing sets along with their corresponding survival data.
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

# ============================================================
# Function: perform_feature_elimination
# ============================================================

# Description:
# This function performs feature elimination recursively by evaluating a Random Forest model at each iteration. 
# At each step, the feature with the lowest importance score (as calculated by Random Forest) is removed, 
# and the model's performance (C-index) is evaluated. The process continues until only one feature remains in the data.
# The function stores the selected features and C-index values at each step in a list and returns the list.

# Parameters:
# - train_data: A data frame containing the training data (features) used for model training.
# - test_data: A data frame containing the test data (features) used for model evaluation.
# - train_survival_data: A vector or data frame containing the survival data for the training set.
# - test_survival_data: A vector or data frame containing the survival data for the test set.
# - count: An integer counter to track the iteration number and store results in a list.
# - feature_c_index_list: A list to store the selected features and the corresponding C-index values at each iteration.

# Returns:
# - A list of results where each element is a list containing:
#   - `variables`: The names of the selected features.
#   - `c_index`: The C-index value computed from the model's predictions.

# Procedure:
# 1. The function iterates while there are more than one feature in the training data.
# 2. At each iteration, it trains a Random Forest model (`ranger` package) on the current features.
# 3. The model's feature importance scores are calculated, and features are ranked in descending order of importance.
# 4. The model's C-index is calculated on the test data using the predictions.
# 5. The least important feature (the one with the lowest importance score) is removed from both the training and test data.
# 6. The process continues until only one feature remains in the data.
# 7. The function returns a list of the selected features and corresponding C-index values at each iteration.

# Notes:
# - The Random Forest model uses a maxstat split rule, with a minimum node size of 50 and 1000 trees.
# - The function uses the concordance function from survival analysis to compute the C-index.
# - The function stops when only one feature remains in the dataset.
# - Feature importance is evaluated using permutation importance.
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

# ============================================================
# Function: get_max_c_index_index
# ============================================================

# Description:
# This function takes a list of results and finds the index of the element that has the highest C-index value. 
# The function compares the C-index value in each element and returns the index of the element with the maximum value.

# Parameters:
# - result_list: A list where each element is expected to be a list or data frame containing a field `c_index` representing the C-index value.

# Returns:
# - An integer index corresponding to the element in the list that has the highest C-index value.

# Procedure:
# 1. Initialize `max_c_index` to a very low value (negative infinity) and `max_index` as NULL.
# 2. Iterate through each element in the `result_list` and extract the C-index value.
# 3. Compare each C-index value to the current `max_c_index` and update the `max_c_index` and `max_index` if a higher value is found.
# 4. Return the index of the element with the highest C-index value.

# Notes:
# - The function assumes that each element in the list contains a `c_index` field or attribute.
# - If the list is empty or no C-index values are found, the function will return `NULL`.
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

# ============================================================
# Function: min_max_normalize
# ============================================================

# Description:
# This function normalizes the numeric columns of a data frame (excluding the first three columns) using the Min-Max scaling method. 
# It transforms the values in each column to a range between 0 and 1 by subtracting the minimum value and dividing by the range 
# (max - min) of that column.

# Parameters:
# - x: A data frame where columns beyond the first three are numeric and will be normalized using Min-Max scaling.

# Returns:
# - A data frame with the same structure as the input, but with normalized values in columns beyond the first three.

# Procedure:
# 1. Check if the number of columns in `x` is greater than 3.
# 2. Apply Min-Max normalization on all columns excluding the first three. The normalization formula is:
#    (value - min(value)) / (max(value) - min(value)).
# 3. Return the modified data frame.

# Notes:
# - The first three columns are not normalized, assuming they are identifiers and survival data.
# - The function assumes that the columns to be normalized are numeric.
min_max_normalize <- function(x) {
  if (ncol(x) > 3) {
    x[, -c(1:3)] <- apply(x[, -c(1:3)], 2, function(y) (y - min(y)) / (max(y) - min(y)))
  }
  return(x)
}


# ============================================================
# Function: fusion
# ============================================================

# Description:
# This function performs a feature elimination process for survival analysis models, evaluates performance 
# using the concordance index (C-index), and generates a plot of the relationship between the number of features 
# and the C-index values. It also saves the best results and associated plot to a specified directory.

# Parameters:
# - data: A data frame containing the features for the survival analysis.
# - survival_data: A vector or data frame of survival times and censoring indicators.
# - directory (optional): The directory where the results (plot and data) will be saved. Default is the current working directory (".").
# - file_prefix (optional): A prefix for the output filenames. Default is "output".

# Returns:
# - Saves a plot showing the relationship between the number of features and the C-index values as a PNG file.
# - Saves the best feature set and corresponding C-index value to a CSV file.

# Procedure:
# 1. Prepare the data by removing the first column (assumed to be an identifier column).
# 2. Generate the train and test datasets using the `generate_train_test` function.
# 3. Perform feature elimination on the training and testing data using the `perform_feature_elimination` function.
# 4. Extract the number of selected features and corresponding C-index values from the results of feature elimination.
# 5. Create a scatter plot with the number of features on the x-axis and the C-index values on the y-axis, and connect the points with a line.
# 6. Ensure the directory exists or create it if needed.
# 7. Find the index of the best result based on the highest C-index value using the `get_max_c_index_index` function.
# 8. Save the plot and the best results (feature set and C-index) to the specified directory.

# Notes:
# - The function assumes the existence of other helper functions: `generate_train_test`, `perform_feature_elimination`, and `get_max_c_index_index`.
# - The plot saved will be named based on the `file_prefix` argument, followed by "_plot.png".
# - The best results will be saved in a CSV file, named based on the `file_prefix`, followed by "_data.csv".
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
                           
# ============================================================
# Function: evaluate_performance
# ============================================================

# Description:
# This function evaluates the performance of four survival analysis models: Cox Proportional Hazards (CoxPH), 
# Random Forest (Ranger), Boosted Cox (GLMBoost), and Elastic Net (Cox model). It calculates the concordance 
# index (C-index) for each model on a test set and returns the results in a data frame.

# Parameters:
# - train_df: A data frame containing the training data (features).
# - train_survival: A vector or data frame of survival times and censoring indicators for the training set.
# - test_df: A data frame containing the test data (features).
# - test_survival: A vector or data frame of survival times and censoring indicators for the test set.

# Returns:
# - A data frame (`cindex_df`) containing the C-index results for each of the four models (CoxPH, Random Forest,
#   GLMBoost, and ElasticNet) on the test set.

# Procedure:
# 1. Prepare the training and testing matrices for the survival analysis models.
# 2. Fit and evaluate the following models:
#    - Cox Proportional Hazards model (CoxPH) and calculate the C-index.
#    - Random Forest model (Ranger) and calculate the C-index.
#    - Boosted Cox Proportional Hazards model (GLMBoost) and calculate the C-index.
#    - Elastic Net Cox model (ElasticNet) and calculate the C-index.
# 3. Combine the C-index results from all models into a data frame.
# 4. Return the data frame containing the C-index results for each model.

# Notes:
# - The C-index is used as a measure of the model's discriminatory ability (the ability to correctly rank pairs of 
#   subjects in terms of their predicted survival times).
evaluate_performance <- function(train_df, train_survival, test_df, test_survival, n_bootstrap = 30, conf_level = 0.95) {
  # Initialize function to calculate bootstrap confidence intervals and p-values
  bootstrap_ci <- function(predictions, true_values, n_bootstrap, conf_level) {
    c_index_values <- numeric(n_bootstrap)
    for (i in 1:n_bootstrap) {
      sample_indices <- sample(1:length(true_values), length(true_values), replace = TRUE)
      boot_true_values <- true_values[sample_indices]
      boot_predictions <- predictions[sample_indices]
      c_index_values[i] <- concordance(boot_true_values ~ boot_predictions)$concordance
    }
    
    # Compute the confidence interval
    lower_ci <- quantile(c_index_values, (1 - conf_level) / 2)
    upper_ci <- quantile(c_index_values, 1 - (1 - conf_level) / 2)
    mean_c_index <- mean(c_index_values)
    std_c_index <- sd(c_index_values)
    
    # Assuming p-value calculation is based on comparing mean c-index to a reference or baseline model
    # For simplicity, let's just calculate a one-sided p-value assuming a normal distribution
    p_value <- 2 * (1 - pnorm(mean_c_index, mean = 0.5, sd = std_c_index))  # Example: comparing to 0.5 baseline
    
    return(list(mean_c_index = mean_c_index, std_c_index = std_c_index, lower_ci = lower_ci, upper_ci = upper_ci, p_value = p_value))
  }
  
  # Prepare matrices
  train_x_matrix <- as.matrix(train_df)
  test_x_matrix <- as.matrix(test_df)
  
  # Cox Model
  eval_cox_model <- coxph(train_survival ~ ., data = train_df)
  cox_predictions <- predict(eval_cox_model, newdata = test_df)
  c_index_cox <- concordance(test_survival ~ cox_predictions)
  cox_results <- bootstrap_ci(cox_predictions, test_survival, n_bootstrap, conf_level)
  
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
  ranger_results <- bootstrap_ci(ranger_predictions, test_survival, n_bootstrap, conf_level)
  
  # Boosted Cox Model (GLMBoost)
  eval_boosted_cox_model <- glmboost(train_survival ~ ., data = train_df, family = CoxPH())
  glmboost_predictions <- predict(eval_boosted_cox_model, newdata = test_df)
  glmboost_results <- bootstrap_ci(glmboost_predictions, test_survival, n_bootstrap, conf_level)
  
  # Elastic Net Cox Model
  elastic_net_model <- cv.glmnet(x = train_x_matrix, y = train_survival, family = "cox", alpha = 0.5)
  elastic_net_predictions <- predict(elastic_net_model, newx = test_x_matrix, s = "lambda.min")
  elastic_net_results <- bootstrap_ci(elastic_net_predictions, test_survival, n_bootstrap, conf_level)
  
  # Create data frame with results
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

# ============================================================
# Function: combine_results
# ============================================================

# Description:
# This function combines the results of feature selection and model performance 
# (C-index) across multiple methods, calculates the mean, standard deviation, 
# and 95% confidence intervals for the features selected and the C-index values.

# Parameters:
# - results: A list containing the following elements:
#     - feature_cox_results: Feature selection results for Cox regression.
#     - feature_ranger_results: Feature selection results for Random Forest.
#     - feature_glmboost_results: Feature selection results for GLMBoost.
#     - feature_elasticnet_results: Feature selection results for ElasticNet.
#     - c_index_cox_results: C-index results for Cox regression.
#     - c_index_ranger_results: C-index results for Random Forest.
#     - c_index_glmboost_results: C-index results for GLMBoost.
#     - c_index_elasticnet_results: C-index results for ElasticNet.
#     - selected_features_list: List of selected features from each method.

# Returns:
# - A list containing two dataframes:
#     - FeatureSelected: A dataframe summarizing feature selection results 
#       including mean values, standard deviation, and 95% confidence intervals.
#     - CIndex: A dataframe summarizing C-index results including mean values, 
#       standard deviation, and 95% confidence intervals.

# Procedure:
# 1. Extracts feature selection results and C-index results from `results`.
# 2. Combines results into data frames (`feature_selection_df`, `feature_cindex_df`).
# 3. Calculates standard deviation for each method.
# 4. Calculates mean values for feature selection and C-index results.
# 5. Calculates 95% confidence intervals for each method.
# 6. Creates two summary dataframes: one for feature selection and one for C-index.
# 7. Returns a list containing the summary dataframes.
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
  
  # Sample sizes
  n_feature <- sapply(feature_selection_df, length)
  n_cindex <- sapply(feature_cindex_df, length)
  
  # Calculate 95% confidence intervals
  ci_feature <- 1.96 * c(cox_feature_sd, ranger_feature_sd, boost_feature_sd, penalised_feature_sd) / sqrt(n_feature)
  ci_cindex <- 1.96 * c(cox_Cindex_sd, ranger_Cindex_sd, boost_Cindex_sd, penalised_Cindex_sd) / sqrt(n_cindex)
  
  # Compute lower and upper CI bounds
  lower_ci_feature <- mean_feature_values - ci_feature
  upper_ci_feature <- mean_feature_values + ci_feature
  
  lower_ci_cindex <- mean_cindex_values - ci_cindex
  upper_ci_cindex <- mean_cindex_values + ci_cindex
  
  # Create summary data frames
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

# ============================================================
# Function: plot_results
# ============================================================

# Description:
# This function generates heatmaps and saves the plots as PNG files. It compares the performance (C-index) 
# and the number of features selected by different machine learning models across various feature selection 
# methods. It also saves the combined results into CSV files.

# Parameters:
# - no_fs_results: A list containing the results of the feature selection with no feature selection (FS).
#   Must include 'FeatureSelected' and 'CIndex' components.
# - multi_results: A list containing the results of the feature selection using multivariate selection.
#   Must include 'FeatureSelected' and 'CIndex' components.
# - uni_results: A list containing the results of the feature selection using univariate selection.
#   Must include 'FeatureSelected' and 'CIndex' components.
# - rsf_vi_results: A list containing the results of the Random Forest with variable importance.
#   Must include 'FeatureSelected' and 'CIndex' components.
# - rsf_md_results: A list containing the results of the Random Forest with minimum depth.
#   Must include 'FeatureSelected' and 'CIndex' components.
# - rsf_vh_results: A list containing the results of the Random Forest with maximum stat.
#   Must include 'FeatureSelected' and 'CIndex' components.
# - custom_name: A custom string to prefix the generated plot and result filenames.

# Returns:
# - Saves the following plots and result files in a directory named `custom_name_plots`:
#     - Mean performance plot (`custom_name_mean_performance.png`): Heatmap of the C-index values for each feature selection method and model.
#     - Mean features plot (`custom_name_mean_features.png`): Heatmap of the number of features selected for each feature selection method and model.
#     - CSV files containing the combined results of feature selection and C-index for each method (`_no_fs_results.csv`, `_multi_results.csv`, etc.).

# Procedure:
# 1. Create a directory to store the plots if it doesn't already exist.
# 2. Extract the summary of feature selection results and C-index from each input result.
# 3. Combine the C-index results from different methods and models into a single data frame.
# 4. Create a heatmap of the C-index values and save it as a PNG file.
# 5. Combine the feature selection results from different methods and models into a single data frame.
# 6. Create a heatmap of the number of features selected and save it as a PNG file.
# 7. Save the combined results to CSV files.
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
}

# ============================================================
# Function: perform_feature_parallel_selection
# ============================================================

# Description:
# This function performs parallelized feature selection using different machine learning models 
# (Cox regression, Random Forest, GLMBoost, and Elastic Net) with repeated k-fold cross-validation. 
# It returns the selected features and performance metrics (C-index) for each model.

# Parameters:
# - data: The dataset to use for feature selection. Rows are samples, columns are features.
# - survival_data: The survival data containing overall survival and deceased information.
# - repeats: The number of repeats for cross-validation.
# - folds: The number of folds for cross-validation.
# - feature_selection: A function that selects features from the training data.

# Returns:
# - A list containing:
#     - selected_features_list: A list of selected features for each fold and repeat.
#     - c_index_elasticnet_results: C-index results for Elastic Net.
#     - c_index_cox_results: C-index results for Cox regression.
#     - c_index_ranger_results: C-index results for Random Forest.
#     - c_index_glmboost_results: C-index results for GLMBoost.
#     - feature_elasticnet_results: Number of features selected by Elastic Net.
#     - feature_cox_results: Number of features selected by Cox regression.
#     - feature_ranger_results: Number of features selected by Random Forest.
#     - feature_glmboost_results: Number of features selected by GLMBoost.

# Procedure:
# 1. Set up a parallel backend using `doParallel` with the specified number of cores.
# 2. For each repeat and fold, divide the data into training and test sets using `createFolds`.
# 3. Perform feature selection using the specified `feature_selection` function.
# 4. For each model (Cox, Random Forest, GLMBoost, Elastic Net), train the model on the training data.
#    - Evaluate each model using C-index on the test data.
#    - Store the selected features and performance metrics (C-index).
# 5. Stop the parallel backend after computation.
# 6. Return a list containing all the results from each repeat and fold.
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
                           
# ============================================================
# Function: univariate
# ============================================================

# Description:
# This function performs univariate Cox proportional hazards regression 
# for each feature in the dataset to assess its association with survival.

# Parameters:
# - data: A dataframe containing survival data (overall_survival, deceased)
#         and feature variables.

# Returns:
# - significant_uni_cox_df: A dataframe containing features with significant 
#   CoxPH results (p-value < 0.05 and non-zero coefficients), including:
#     - feature: Feature name
#     - coefficients: Cox regression coefficient
#     - p_value: P-value from the likelihood ratio test
#     - c_index: Concordance index (C-index) measuring predictive ability

# Procedure:
# 1. Initializes an empty dataframe (`uni_cox_df`) to store CoxPH results.
# 2. Iterates over each feature (starting from the 3rd column of `data`).
# 3. Runs a Cox proportional hazards model for each feature:
#    - Dependent variable: overall_survival, deceased
#    - Independent variable: Single feature
# 4. Extracts:
#    - P-value from the likelihood ratio test
#    - C-index (concordance index)
#    - Coefficient of the feature
# 5. Filters features with p-value < 0.05 and non-zero coefficients.
# 6. Returns the filtered dataframe (`significant_uni_cox_df`).
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

# ============================================================
# Function: multivariate
# ============================================================

# Description:
# This function performs multivariate Cox proportional hazards regression 
# to identify significant features associated with survival.

# Parameters:
# - data: A dataframe containing survival data (overall_survival, deceased) 
#         and feature variables.

# Returns:
# - significant_cox_df: A dataframe containing significant features 
#   (p-value < 0.05 and non-zero coefficients), including:
#     - feature: Feature name

# Procedure:
# 1. Fits a multivariate Cox proportional hazards model using all features.
# 2. Extracts:
#    - Coefficients of the Cox model.
#    - P-values associated with each feature.
# 3. Selects features with:
#    - P-value < 0.05 (statistical significance).
#    - Non-zero coefficients (contributing to the model).
# 4. Returns a dataframe (`significant_cox_df`) with the selected features.
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

# ============================================================
# Function: rsfvh
# ============================================================

# Description:
# This function performs variable selection using the "Variable Hunting" (VH) 
# method from Random Survival Forests (RSF) to identify important survival-related features.

# Parameters:
# - data: A dataframe containing survival data (overall_survival, deceased) 
#         and feature variables.

# Returns:
# - rsfvhfeature_Selection: A dataframe containing selected features, including:
#     - feature: Feature name

# Procedure:
# 1. Fits a Random Survival Forest model using the `var.select` function with:
#    - ntree = 1000 (number of trees)
#    - method = "vh" (Variable Hunting method)
#    - nodesize = 5 (minimum terminal node size)
#    - nsplit = 20 (number of split points per variable)
#    - splitrule = "logrank" (splitting criterion)
#    - nrep = 3 (number of repetitions for stability)
#    - K = 10 (number of variables randomly selected per split)
#    - nstep = 1 (step size for forward selection)
# 2. Extracts the top selected variables from `rsfvh_model$topvars`.
# 3. Returns a dataframe (`rsfvhfeature_Selection`) with the selected features.
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

# ============================================================
# Function: rsfmd
# ============================================================

# Description:
# This function performs variable selection using the "Minimal Depth" (MD) 
# method from Random Survival Forests (RSF) to identify important survival-related features.

# Parameters:
# - data: A dataframe containing survival data (overall_survival, deceased) 
#         and feature variables.

# Returns:
# - rsfmdfeature_Selection: A dataframe containing selected features, including:
#     - feature: Feature name

# Procedure:
# 1. Fits a Random Survival Forest model using the `var.select` function with:
#    - ntree = 1000 (number of trees)
#    - method = "md" (Minimal Depth method)
#    - nodesize = 5 (minimum terminal node size)
#    - nsplit = 20 (number of split points per variable)
#    - splitrule = "logrank" (splitting criterion)
# 2. Extracts the top selected variables from `rsfmd_model$topvars`.
# 3. Returns a dataframe (`rsfmdfeature_Selection`) with the selected features.
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

# ============================================================
# Function: rsfvi
# ============================================================

# Description:
# This function performs variable selection using the Variable Importance (VI) 
# method from Random Survival Forests (RSF) to identify important survival-related features.

# Parameters:
# - data: A dataframe containing survival data (overall_survival, deceased) 
#         and feature variables.

# Returns:
# - rsffeature_selection: A dataframe containing selected features based on 
#   positive variable importance scores, including:
#     - feature: Feature name
#     - Importance: Variable importance score

# Procedure:
# 1. Computes `mtry_value` as the square root of the total number of features.
# 2. Fits a Random Survival Forest model using `rfsrc` with:
#    - ntree = 1000 (number of trees)
#    - mtry = mtry_value (number of variables randomly selected at each split)
#    - nodesize = 3 (minimum terminal node size)
#    - nsplit = 10 (number of split points per variable)
# 3. Computes variable importance using `vimp(rsf_model)`.
# 4. Extracts features with positive importance scores.
# 5. Returns a dataframe (`rsffeature_selection`) with the selected features.
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

# ============================================================
# Function: all_features
# ============================================================

# Description:
# This function extracts all feature names from the dataset, excluding 
# the first two columns (assumed to contain survival-related data).

# Parameters:
# - data: A dataframe containing survival data (overall_survival, deceased) 
#         and feature variables.

# Returns:
# - all_features_df: A dataframe containing all feature names, including:
#     - features: Feature name

# Procedure:
# 1. Extracts column names from `data`, excluding the first two columns.
# 2. Stores the feature names in a dataframe (`all_features_df`).
# 3. Returns the dataframe with the list of all feature names.
all_features <- function(data) {
  all_features <- colnames(data)[-c(1,2)]
  all_features_df <- data.frame(features = all_features)
  return(all_features_df)
}
                           
# ============================================================
# Function: perform_and_save_results
# ============================================================

# Description:
# This function performs feature selection using multiple methods, 
# combines the results, and plots the comparisons.

# Parameters:
# - data: A dataframe containing feature data.
# - survival_data: A dataframe containing survival information.
# - repeats: Number of times cross-validation should be repeated.
# - folds: Number of folds for cross-validation.
# - custom_name: A string used for naming the plot.

# Procedure:
# 1. Runs feature selection without filtering (`all_features`).
# 2. Runs feature selection using:
#    - Multivariate Cox regression (`multivariate`).
#    - Univariate Cox regression (`univariate`).
#    - Random Survival Forest Variable Importance (`rsfvi`).
#    - Random Survival Forest Minimal Depth (`rsfmd`).
#    - Random Survival Forest Variable Hunting (`rsfvh`).
# 3. Combines the results for each method.
# 4. Plots the performance comparison across feature selection methods.

# Dependencies:
# - Requires `perform_feature_parallel_selection`, `combine_results`, 
#   and `plot_results` functions.
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


# ============================================================
# Function: perform_feature_selection_CV & perform_feature_parallel_CV
# ============================================================

# Description:
# Performs feature selection using multiple filter methods with cross-validation.

# Parameters:
# - data: Input dataset containing omics features and survival data.
# - labels: Corresponding survival labels.
# - n_folds: Number of folds for cross-validation (default: 5).
# - n_repeats: Number of repetitions for cross-validation (default: 5).
# - methods: List of feature selection methods to apply.

# Returns:
# - selected_features: A list of features that appear ≥50 times across CV runs.

# Procedure:
# 1. Runs each filter method (CoxPH, Random Forest, etc.).
# 2. Extracts important features based on:
#    - Non-zero coefficients
#    - P-values (CoxPH)
#    - Importance scores (Random Forest)
# 3. Repeats this process 100 times (5-fold, 5-repeats for each method).
# 4. Tracks feature occurrences and selects those appearing ≥50 times.
# 5. Outputs final feature set for model evaluation.                   
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

# ============================================================
# Function: perform_bootstrapping
# ============================================================

# Description:
# This function performs bootstrapping and parallelized feature selection.
# It selects features based on multiple methods and tracks their occurrences
# across 100 bootstrap iterations.

# Parameters:
# - data: A dataframe containing omics features and survival data.

# Returns:
# - sorted_feature_count_df: A dataframe containing selected features and 
#   their occurrence counts across bootstrap iterations.

# Procedure:
# 1. Defines an internal function (feature_selection) to select features using:
#    - Multivariate CoxPH (`multivariate`)
#    - Random Survival Forest Variable Importance (`rsfvi`)
#    - Random Survival Forest Variable Hunting (`rsfvh`)
#    - Random Survival Forest Minimal Depth (`rsfmd`)
#    The function identifies common features across all four methods.
#
# 2. Initializes parallel computing using all available cores (16).
#
# 3. Runs feature selection on 100 bootstrap samples (70% of the original data)
#    in parallel:
#    - Each iteration draws a bootstrap sample.
#    - Feature selection is performed.
#    - The selected features are stored.
#
# 4. Stops parallel processing after execution.
#
# 5. Counts occurrences of each selected feature across bootstrap iterations.
#
# 6. Sorts features by occurrence count in descending order.
#
# 7. Returns a dataframe (`sorted_feature_count_df`) with feature names and their counts.                      
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

                           
# ============================================================
# Function: categorize_cnv
# ============================================================

# Description:
# This function categorizes copy number variation (CNV) values into predefined categories based on the copy number.
# It assigns numeric values to each category for further analysis.

# Parameters:
# - copy_number: A numeric value representing the copy number for a particular gene or genomic region.

# Returns:
# - A numeric value representing the CNV category:
#     - -2: Homozygous deletion (copy number < 0.5)
#     - -1: Hemizygous deletion (0.5 <= copy number < 1.5)
#     - 0: Neutral / No change (1.5 <= copy number <= 2.5)
#     - 1: Low-level amplification (2.5 < copy number <= 4)
#     - 2: High-level amplification (copy number > 4)
#     - NA: For missing or unexpected copy number values

# Procedure:
# 1. The function checks the value of `copy_number` and categorizes it based on predefined thresholds using the `case_when` function.
# 2. The function returns an integer representing the category corresponding to the provided `copy_number` value.
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

# ============================================================
# Function: create_cnv_km_plots
# ============================================================

# Description:
# This function creates a Kaplan-Meier (KM) survival plot for a specific gene's copy number variation (CNV) values.
# It categorizes the CNV values into Gain, Loss, or Neutral, then fits a survival model and plots the Kaplan-Meier curve.

# Parameters:
# - data: A data frame containing the following columns:
#     - `case_id`: Unique identifier for each case.
#     - `overall_survival`: Survival time in days.
#     - `deceased`: Binary variable (1 if deceased, 0 if alive).
#     - `gene`: A column representing the CNV values for a specific gene.
# - gene: The name of the gene column in the data frame to be used for CNV analysis.
# - cancer: The name of the cancer type used for the plot title.

# Returns:
# - A Kaplan-Meier plot showing survival curves based on CNV categories (Gain, Loss, Neutral) for the specified gene.
#   The plot includes:
#     - The survival probability over time.
#     - A risk table indicating the number of cases at risk at various time points.
#     - A p-value indicating the statistical significance of the survival difference between the CNV categories.

# Procedure:
# 1. The function selects relevant columns (`case_id`, `overall_survival`, `deceased`, and the specified `gene` column).
# 2. It ensures that the CNV values for the gene are numeric.
# 3. The CNV values are categorized into three groups: "Gain" for positive values, "Loss" for negative values, and "Neutral" for zero.
# 4. The survival data is then used to create a Kaplan-Meier survival object.
# 5. The function fits a survival model based on the categorized CNV values and generates the Kaplan-Meier curve.
# 6. The plot includes a risk table, p-value, and customized plot labels for better interpretation.
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

# ============================================================
# Function: create_ge_km_plots
# ============================================================

# Description:
# This function creates a Kaplan-Meier (KM) survival plot for a specific gene's expression (GE) values.
# The gene expression values are categorized into "High" and "Low" based on the median value, then the Kaplan-Meier curve is plotted.

# Parameters:
# - data: A data frame containing the following columns:
#     - `case_id`: Unique identifier for each case.
#     - `overall_survival`: Survival time in days.
#     - `deceased`: Binary variable (1 if deceased, 0 if alive).
#     - `gene`: A column representing the gene expression values for a specific gene.
# - gene: The name of the gene column in the data frame to be used for expression analysis.
# - cancer: The name of the cancer type used for the plot title.

# Returns:
# - A Kaplan-Meier plot showing survival curves based on the high and low gene expression levels.
#   The plot includes:
#     - The survival probability over time.
#     - A risk table indicating the number of cases at risk at various time points.
#     - A p-value indicating the statistical significance of the survival difference between high and low expression groups.

# Procedure:
# 1. The function selects relevant columns (`case_id`, `overall_survival`, `deceased`, and the specified `gene` column).
# 2. It calculates the median gene expression value to create the "High" and "Low" categories.
# 3. The data is then split into "High" and "Low" expression categories based on whether the gene expression is greater than or equal to the median or less than it.
# 4. NA values in the strata are removed.
# 5. A Kaplan-Meier survival object is created based on the survival data.
# 6. The function fits a survival model using the categorized gene expression levels and generates the Kaplan-Meier curve.
# 7. The plot includes a risk table, p-value, and customized plot labels for better interpretation.
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
