library("survival")
library("survminer")
library(foreach)
library(glmnet)
library(ggplot2)
library(ranger)
library(mboost)
library(caret)
library(randomForestSRC)
library(viridis)
library(tidyverse)
library(reshape)
library(parallel)
library(doParallel)
.libPaths("/path/to/Rlib")
setwd("/path/to/pwd")
source("functions.R")
# Set seed for reproducibility
set.seed(123)
# Number of repeats and folds
repeats <- 5
folds <- 5

####################################################################
#                           BRCA
#
####################################################################
# ME ------------------------------------------------------------
BRCA_ME=read.csv("BRCA/BRCA_ME_clean.csv", row.names = NULL)
BRCA_ME$overall_survival[BRCA_ME$overall_survival <= 0] <- 0.001
BRCA_ME <- BRCA_ME[!is.na(BRCA_ME$overall_survival), ]
get_features <- colnames(BRCA_ME)[-c(1:3)]
train_test_data <- generate_train_test(BRCA_ME)
BRCA_ME_train_df <- train_test_data$train_df
BRCA_ME_test_df <- train_test_data$test_df
feature_counts_df <- data.frame(Feature = get_features, Count = 0)
perform_feature_selection_CV(BRCA_ME_train_df,feature_counts_df, "BRCA/ME")
train_survival_cols <- BRCA_ME_train_df[, c("overall_survival", "deceased")]
test_survival_cols  <- BRCA_ME_test_df[, c("overall_survival", "deceased")]
sorted_feature_counts=read.csv("BRCA/ME/features_cv.csv", row.names = NULL)
BRCA_ME_test <- BRCA_ME_test_df[, intersect(colnames(BRCA_ME_test_df), sorted_feature_counts$Feature)]
BRCA_ME_train <- BRCA_ME_train_df[, intersect(colnames(BRCA_ME_train_df), sorted_feature_counts$Feature)]
BRCA_ME_train <- cbind(train_survival_cols, BRCA_ME_train)
BRCA_ME_test  <- cbind(test_survival_cols, BRCA_ME_test)
cv_result <- evaluate_performance(BRCA_ME_train, BRCA_ME_test)
boostrap_feature_counts <- perform_bootstrapping(BRCA_ME_train_df[,-1])
# Compute the 70th percentile threshold
threshold <- quantile(boostrap_feature_counts$count, probs = 0.70)
# Filter features that meet or exceed this threshold
top_features <- boostrap_feature_counts[boostrap_feature_counts$count >= threshold, ]
write.csv(top_features, file = "BRCA/ME/features_bs", row.names = FALSE)
BRCA_ME_test <- BRCA_ME_test_df[, intersect(colnames(BRCA_ME_test_df), top_features$feature)]
BRCA_ME_train <- BRCA_ME_train_df[, intersect(colnames(BRCA_ME_train_df), top_features$feature)]
BRCA_ME_train <- cbind(train_survival_cols, BRCA_ME_train)
BRCA_ME_test  <- cbind(test_survival_cols, BRCA_ME_test)
boostrap_result <- evaluate_performance(BRCA_ME_train, BRCA_ME_test)
no_feature_selection <- evaluate_performance(BRCA_ME_train_df[,-1], BRCA_ME_test_df[,-1]) 
write.csv(cv_result, file = "BRCA/ME/results_with_cv", row.names = FALSE)
write.csv(boostrap_result, file = "BRCA/ME/results_with_bs", row.names = FALSE)
write.csv(no_feature_selection, file = "BRCA/ME/results_without_fs", row.names = FALSE)

# GE ------------------------------------------------------------
BRCA_GE=read.csv("BRCA/BRCA_GE_clean.csv", row.names = NULL)
BRCA_GE$overall_survival[BRCA_GE$overall_survival <= 0] <- 0.001
BRCA_GE <- BRCA_GE[!is.na(BRCA_GE$overall_survival), ]
get_features <- colnames(BRCA_GE)[-c(1:3)]
train_test_data <- generate_train_test(BRCA_GE)
BRCA_GE_train_df <- train_test_data$train_df
BRCA_GE_test_df <- train_test_data$test_df
feature_counts_df <- data.frame(Feature = get_features, Count = 0)
perform_feature_selection_CV(BRCA_GE_train_df,feature_counts_df, "BRCA/GE")
train_survival_cols <- BRCA_GE_train_df[, c("overall_survival", "deceased")]
test_survival_cols  <- BRCA_GE_test_df[, c("overall_survival", "deceased")]
sorted_feature_counts=read.csv("BRCA/GE/features_cv.csv", row.names = NULL)
BRCA_GE_test <- BRCA_GE_test_df[, intersect(colnames(BRCA_GE_test_df), sorted_feature_counts$Feature)]
BRCA_GE_train <- BRCA_GE_train_df[, intersect(colnames(BRCA_GE_train_df), sorted_feature_counts$Feature)]
BRCA_GE_train <- cbind(train_survival_cols, BRCA_GE_train)
BRCA_GE_test  <- cbind(test_survival_cols, BRCA_GE_test)
cv_result <- evaluate_performance(BRCA_GE_train, BRCA_GE_test)
boostrap_feature_counts <- perform_bootstrapping(BRCA_GE_train_df[,-1])
# Compute the 70th percentile threshold
threshold <- quantile(boostrap_feature_counts$count, probs = 0.70)
# Filter features that meet or exceed this threshold
top_features <- boostrap_feature_counts[boostrap_feature_counts$count >= threshold, ]
write.csv(top_features, file = "BRCA/GE/features_bs", row.names = FALSE)
BRCA_GE_test <- BRCA_GE_test_df[, intersect(colnames(BRCA_GE_test_df), top_features$feature)]
BRCA_GE_train <- BRCA_GE_train_df[, intersect(colnames(BRCA_GE_train_df), top_features$feature)]
BRCA_GE_train <- cbind(train_survival_cols, BRCA_GE_train)
BRCA_GE_test  <- cbind(test_survival_cols, BRCA_GE_test)
boostrap_result <- evaluate_performance(BRCA_GE_train, BRCA_GE_test)
no_feature_selection <- evaluate_performance(BRCA_GE_train_df[,-1], BRCA_GE_test_df[,-1]) 
write.csv(cv_result, file = "BRCA/GE/results_with_cv", row.names = FALSE)
write.csv(boostrap_result, file = "BRCA/GE/results_with_bs", row.names = FALSE)
write.csv(no_feature_selection, file = "BRCA/GE/results_without_fs", row.names = FALSE)

# DM ------------------------------------------------------------
BRCA_DM=read.csv("BRCA/BRCA_DM_clean.csv", row.names = NULL)
BRCA_DM$overall_survival[BRCA_DM$overall_survival <= 0] <- 0.001
BRCA_DM <- BRCA_DM[!is.na(BRCA_DM$overall_survival), ]
get_features <- colnames(BRCA_DM)[-c(1:3)]
train_test_data <- generate_train_test(BRCA_DM)
BRCA_DM_train_df <- train_test_data$train_df
BRCA_DM_test_df <- train_test_data$test_df
feature_counts_df <- data.frame(Feature = get_features, Count = 0)
perform_feature_selection_CV(BRCA_DM_train_df,feature_counts_df, "BRCA/DM")
train_survival_cols <- BRCA_DM_train_df[, c("overall_survival", "deceased")]
test_survival_cols  <- BRCA_DM_test_df[, c("overall_survival", "deceased")]
sorted_feature_counts=read.csv("BRCA/DM/features_cv.csv", row.names = NULL)
BRCA_DM_test <- BRCA_DM_test_df[, intersect(colnames(BRCA_DM_test_df), sorted_feature_counts$Feature)]
BRCA_DM_train <- BRCA_DM_train_df[, intersect(colnames(BRCA_DM_train_df), sorted_feature_counts$Feature)]
BRCA_DM_train <- cbind(train_survival_cols, BRCA_DM_train)
BRCA_DM_test  <- cbind(test_survival_cols, BRCA_DM_test)
cv_result <- evaluate_performance(BRCA_DM_train, BRCA_DM_test)
boostrap_feature_counts <- perform_bootstrapping(BRCA_DM_train_df[,-1])
# Compute the 70th percentile threshold
threshold <- quantile(boostrap_feature_counts$count, probs = 0.70)
# Filter features that meet or exceed this threshold
top_features <- boostrap_feature_counts[boostrap_feature_counts$count >= threshold, ]
write.csv(top_features, file = "BRCA/DM/features_bs", row.names = FALSE)
BRCA_DM_test <- BRCA_DM_test_df[, intersect(colnames(BRCA_DM_test_df), top_features$feature)]
BRCA_DM_train <- BRCA_DM_train_df[, intersect(colnames(BRCA_DM_train_df), top_features$feature)]
BRCA_DM_train <- cbind(train_survival_cols, BRCA_DM_train)
BRCA_DM_test  <- cbind(test_survival_cols, BRCA_DM_test)
boostrap_result <- evaluate_performance(BRCA_DM_train, BRCA_DM_test)
no_feature_selection <- evaluate_performance(BRCA_DM_train_df[,-1], BRCA_DM_test_df[,-1]) 
write.csv(cv_result, file = "BRCA/DM/results_with_cv", row.names = FALSE)
write.csv(boostrap_result, file = "BRCA/DM/results_with_bs", row.names = FALSE)
write.csv(no_feature_selection, file = "BRCA/DM/results_without_fs", row.names = FALSE)

# CNV ------------------------------------------------------------
BRCA_CNV=read.csv("BRCA/BRCA_CNV_clean.csv", row.names = NULL)
BRCA_CNV$overall_survival[BRCA_CNV$overall_survival <= 0] <- 0.001
BRCA_CNV <- BRCA_CNV[!is.na(BRCA_CNV$overall_survival), ]
get_features <- colnames(BRCA_CNV)[-c(1:3)]
train_test_data <- generate_train_test(BRCA_CNV)
BRCA_CNV_train_df <- train_test_data$train_df
BRCA_CNV_test_df <- train_test_data$test_df
feature_counts_df <- data.frame(Feature = get_features, Count = 0)
perform_feature_selection_CV(BRCA_CNV_train_df,feature_counts_df, "BRCA/CNV")
train_survival_cols <- BRCA_CNV_train_df[, c("overall_survival", "deceased")]
test_survival_cols  <- BRCA_CNV_test_df[, c("overall_survival", "deceased")]
sorted_feature_counts=read.csv("BRCA/CNV/features_cv.csv", row.names = NULL)
BRCA_CNV_test <- BRCA_CNV_test_df[, intersect(colnames(BRCA_CNV_test_df), sorted_feature_counts$Feature)]
BRCA_CNV_train <- BRCA_CNV_train_df[, intersect(colnames(BRCA_CNV_train_df), sorted_feature_counts$Feature)]
BRCA_CNV_train <- cbind(train_survival_cols, BRCA_CNV_train)
BRCA_CNV_test  <- cbind(test_survival_cols, BRCA_CNV_test)
cv_result <- evaluate_performance(BRCA_CNV_train, BRCA_CNV_test)
boostrap_feature_counts <- perform_bootstrapping(BRCA_CNV_train_df[,-1])
# Compute the 70th percentile threshold
threshold <- quantile(boostrap_feature_counts$count, probs = 0.70)
# Filter features that meet or exceed this threshold
top_features <- boostrap_feature_counts[boostrap_feature_counts$count >= threshold, ]
write.csv(top_features, file = "BRCA/CNV/features_bs", row.names = FALSE)
BRCA_CNV_test <- BRCA_CNV_test_df[, intersect(colnames(BRCA_CNV_test_df), top_features$feature)]
BRCA_CNV_train <- BRCA_CNV_train_df[, intersect(colnames(BRCA_CNV_train_df), top_features$feature)]
BRCA_CNV_train <- cbind(train_survival_cols, BRCA_CNV_train)
BRCA_CNV_test  <- cbind(test_survival_cols, BRCA_CNV_test)
boostrap_result <- evaluate_performance(BRCA_CNV_train, BRCA_CNV_test)
no_feature_selection <- evaluate_performance(BRCA_CNV_train_df[,-1], BRCA_CNV_test_df[,-1]) 
write.csv(cv_result, file = "BRCA/CNV/results_with_cv", row.names = FALSE)
write.csv(boostrap_result, file = "BRCA/CNV/results_with_bs", row.names = FALSE)
write.csv(no_feature_selection, file = "BRCA/CNV/results_without_fs", row.names = FALSE)


####################################################################
#                           CESC
#
####################################################################
# ME ------------------------------------------------------------
CESC_ME=read.csv("CESC/CESC_ME_clean.csv", row.names = NULL)
CESC_ME$overall_survival[CESC_ME$overall_survival <= 0] <- 0.001
CESC_ME <- CESC_ME[!is.na(CESC_ME$overall_survival), ]
get_features <- colnames(CESC_ME)[-c(1:3)]
train_test_data <- generate_train_test(CESC_ME)
CESC_ME_train_df <- train_test_data$train_df
CESC_ME_test_df <- train_test_data$test_df
feature_counts_df <- data.frame(Feature = get_features, Count = 0)
perform_feature_selection_CV(CESC_ME_train_df,feature_counts_df, "CESC/ME")
train_survival_cols <- CESC_ME_train_df[, c("overall_survival", "deceased")]
test_survival_cols  <- CESC_ME_test_df[, c("overall_survival", "deceased")]
sorted_feature_counts=read.csv("CESC/ME/features_cv.csv", row.names = NULL)
CESC_ME_test <- CESC_ME_test_df[, intersect(colnames(CESC_ME_test_df), sorted_feature_counts$Feature)]
CESC_ME_train <- CESC_ME_train_df[, intersect(colnames(CESC_ME_train_df), sorted_feature_counts$Feature)]
CESC_ME_train <- cbind(train_survival_cols, CESC_ME_train)
CESC_ME_test  <- cbind(test_survival_cols, CESC_ME_test)
cv_result <- evaluate_performance(CESC_ME_train, CESC_ME_test)
boostrap_feature_counts <- perform_bootstrapping(CESC_ME_train_df[,-1])
# Compute the 70th percentile threshold
threshold <- quantile(boostrap_feature_counts$count, probs = 0.70)
# Filter features that meet or exceed this threshold
top_features <- boostrap_feature_counts[boostrap_feature_counts$count >= threshold, ]
write.csv(top_features, file = "CESC/ME/features_bs", row.names = FALSE)
CESC_ME_test <- CESC_ME_test_df[, intersect(colnames(CESC_ME_test_df), top_features$feature)]
CESC_ME_train <- CESC_ME_train_df[, intersect(colnames(CESC_ME_train_df), top_features$feature)]
CESC_ME_train <- cbind(train_survival_cols, CESC_ME_train)
CESC_ME_test  <- cbind(test_survival_cols, CESC_ME_test)
boostrap_result <- evaluate_performance(CESC_ME_train, CESC_ME_test)
no_feature_selection <- evaluate_performance(CESC_ME_train_df[,-1], CESC_ME_test_df[,-1]) 
write.csv(cv_result, file = "CESC/ME/results_with_cv", row.names = FALSE)
write.csv(boostrap_result, file = "CESC/ME/results_with_bs", row.names = FALSE)
write.csv(no_feature_selection, file = "CESC/ME/results_without_fs", row.names = FALSE)

# GE ------------------------------------------------------------
CESC_GE=read.csv("CESC/CESC_GE_clean.csv", row.names = NULL)
CESC_GE$overall_survival[CESC_GE$overall_survival <= 0] <- 0.001
CESC_GE <- CESC_GE[!is.na(CESC_GE$overall_survival), ]
get_features <- colnames(CESC_GE)[-c(1:3)]
train_test_data <- generate_train_test(CESC_GE)
CESC_GE_train_df <- train_test_data$train_df
CESC_GE_test_df <- train_test_data$test_df
feature_counts_df <- data.frame(Feature = get_features, Count = 0)
perform_feature_selection_CV(CESC_GE_train_df,feature_counts_df, "CESC/GE")
train_survival_cols <- CESC_GE_train_df[, c("overall_survival", "deceased")]
test_survival_cols  <- CESC_GE_test_df[, c("overall_survival", "deceased")]
sorted_feature_counts=read.csv("CESC/GE/features_cv.csv", row.names = NULL)
CESC_GE_test <- CESC_GE_test_df[, intersect(colnames(CESC_GE_test_df), sorted_feature_counts$Feature)]
CESC_GE_train <- CESC_GE_train_df[, intersect(colnames(CESC_GE_train_df), sorted_feature_counts$Feature)]
CESC_GE_train <- cbind(train_survival_cols, CESC_GE_train)
CESC_GE_test  <- cbind(test_survival_cols, CESC_GE_test)
cv_result <- evaluate_performance(CESC_GE_train, CESC_GE_test)
boostrap_feature_counts <- perform_bootstrapping(CESC_GE_train_df[,-1])
# Compute the 70th percentile threshold
threshold <- quantile(boostrap_feature_counts$count, probs = 0.70)
# Filter features that meet or exceed this threshold
top_features <- boostrap_feature_counts[boostrap_feature_counts$count >= threshold, ]
write.csv(top_features, file = "CESC/GE/features_bs", row.names = FALSE)
CESC_GE_test <- CESC_GE_test_df[, intersect(colnames(CESC_GE_test_df), top_features$feature)]
CESC_GE_train <- CESC_GE_train_df[, intersect(colnames(CESC_GE_train_df), top_features$feature)]
CESC_GE_train <- cbind(train_survival_cols, CESC_GE_train)
CESC_GE_test  <- cbind(test_survival_cols, CESC_GE_test)
boostrap_result <- evaluate_performance(CESC_GE_train, CESC_GE_test)
no_feature_selection <- evaluate_performance(CESC_GE_train_df[,-1], CESC_GE_test_df[,-1]) 
write.csv(cv_result, file = "CESC/GE/results_with_cv", row.names = FALSE)
write.csv(boostrap_result, file = "CESC/GE/results_with_bs", row.names = FALSE)
write.csv(no_feature_selection, file = "CESC/GE/results_without_fs", row.names = FALSE)

# DM ------------------------------------------------------------
CESC_DM=read.csv("CESC/CESC_DM_clean.csv", row.names = NULL)
CESC_DM$overall_survival[CESC_DM$overall_survival <= 0] <- 0.001
CESC_DM <- CESC_DM[!is.na(CESC_DM$overall_survival), ]
get_features <- colnames(CESC_DM)[-c(1:3)]
train_test_data <- generate_train_test(CESC_DM)
CESC_DM_train_df <- train_test_data$train_df
CESC_DM_test_df <- train_test_data$test_df
feature_counts_df <- data.frame(Feature = get_features, Count = 0)
perform_feature_selection_CV(CESC_DM_train_df,feature_counts_df, "CESC/DM")
train_survival_cols <- CESC_DM_train_df[, c("overall_survival", "deceased")]
test_survival_cols  <- CESC_DM_test_df[, c("overall_survival", "deceased")]
sorted_feature_counts=read.csv("CESC/DM/features_cv.csv", row.names = NULL)
CESC_DM_test <- CESC_DM_test_df[, intersect(colnames(CESC_DM_test_df), sorted_feature_counts$Feature)]
CESC_DM_train <- CESC_DM_train_df[, intersect(colnames(CESC_DM_train_df), sorted_feature_counts$Feature)]
CESC_DM_train <- cbind(train_survival_cols, CESC_DM_train)
CESC_DM_test  <- cbind(test_survival_cols, CESC_DM_test)
cv_result <- evaluate_performance(CESC_DM_train, CESC_DM_test)
boostrap_feature_counts <- perform_bootstrapping(CESC_DM_train_df[,-1])
# Compute the 70th percentile threshold
threshold <- quantile(boostrap_feature_counts$count, probs = 0.70)
# Filter features that meet or exceed this threshold
top_features <- boostrap_feature_counts[boostrap_feature_counts$count >= threshold, ]
write.csv(top_features, file = "CESC/DM/features_bs", row.names = FALSE)
CESC_DM_test <- CESC_DM_test_df[, intersect(colnames(CESC_DM_test_df), top_features$feature)]
CESC_DM_train <- CESC_DM_train_df[, intersect(colnames(CESC_DM_train_df), top_features$feature)]
CESC_DM_train <- cbind(train_survival_cols, CESC_DM_train)
CESC_DM_test  <- cbind(test_survival_cols, CESC_DM_test)
boostrap_result <- evaluate_performance(CESC_DM_train, CESC_DM_test)
no_feature_selection <- evaluate_performance(CESC_DM_train_df[,-1], CESC_DM_test_df[,-1]) 
write.csv(cv_result, file = "CESC/DM/results_with_cv", row.names = FALSE)
write.csv(boostrap_result, file = "CESC/DM/results_with_bs", row.names = FALSE)
write.csv(no_feature_selection, file = "CESC/DM/results_without_fs", row.names = FALSE)

# CNV ------------------------------------------------------------
CESC_CNV=read.csv("CESC/CESC_CNV_clean.csv", row.names = NULL)
CESC_CNV$overall_survival[CESC_CNV$overall_survival <= 0] <- 0.001
CESC_CNV <- CESC_CNV[!is.na(CESC_CNV$overall_survival), ]
get_features <- colnames(CESC_CNV)[-c(1:3)]
train_test_data <- generate_train_test(CESC_CNV)
CESC_CNV_train_df <- train_test_data$train_df
CESC_CNV_test_df <- train_test_data$test_df
feature_counts_df <- data.frame(Feature = get_features, Count = 0)
perform_feature_selection_CV(CESC_CNV_train_df,feature_counts_df, "CESC/CNV")
train_survival_cols <- CESC_CNV_train_df[, c("overall_survival", "deceased")]
test_survival_cols  <- CESC_CNV_test_df[, c("overall_survival", "deceased")]
sorted_feature_counts=read.csv("CESC/CNV/features_cv.csv", row.names = NULL)
CESC_CNV_test <- CESC_CNV_test_df[, intersect(colnames(CESC_CNV_test_df), sorted_feature_counts$Feature)]
CESC_CNV_train <- CESC_CNV_train_df[, intersect(colnames(CESC_CNV_train_df), sorted_feature_counts$Feature)]
CESC_CNV_train <- cbind(train_survival_cols, CESC_CNV_train)
CESC_CNV_test  <- cbind(test_survival_cols, CESC_CNV_test)
cv_result <- evaluate_performance(CESC_CNV_train, CESC_CNV_test)
boostrap_feature_counts <- perform_bootstrapping(CESC_CNV_train_df[,-1])
# Compute the 70th percentile threshold
threshold <- quantile(boostrap_feature_counts$count, probs = 0.70)
# Filter features that meet or exceed this threshold
top_features <- boostrap_feature_counts[boostrap_feature_counts$count >= threshold, ]
write.csv(top_features, file = "CESC/CNV/features_bs", row.names = FALSE)
CESC_CNV_test <- CESC_CNV_test_df[, intersect(colnames(CESC_CNV_test_df), top_features$feature)]
CESC_CNV_train <- CESC_CNV_train_df[, intersect(colnames(CESC_CNV_train_df), top_features$feature)]
CESC_CNV_train <- cbind(train_survival_cols, CESC_CNV_train)
CESC_CNV_test  <- cbind(test_survival_cols, CESC_CNV_test)
boostrap_result <- evaluate_performance(CESC_CNV_train, CESC_CNV_test)
no_feature_selection <- evaluate_performance(CESC_CNV_train_df[,-1], CESC_CNV_test_df[,-1]) 
write.csv(cv_result, file = "CESC/CNV/results_with_cv", row.names = FALSE)
write.csv(boostrap_result, file = "CESC/CNV/results_with_bs", row.names = FALSE)
write.csv(no_feature_selection, file = "CESC/CNV/results_without_fs", row.names = FALSE)


####################################################################
#                           OV
#
####################################################################
# ME ------------------------------------------------------------
OV_ME=read.csv("OV/OV_ME_clean.csv", row.names = NULL)
OV_ME$overall_survival[OV_ME$overall_survival <= 0] <- 0.001
OV_ME <- OV_ME[!is.na(OV_ME$overall_survival), ]
get_features <- colnames(OV_ME)[-c(1:3)]
train_test_data <- generate_train_test(OV_ME)
OV_ME_train_df <- train_test_data$train_df
OV_ME_test_df <- train_test_data$test_df
feature_counts_df <- data.frame(Feature = get_features, Count = 0)
perform_feature_selection_CV(OV_ME_train_df,feature_counts_df, "OV/ME")
train_survival_cols <- OV_ME_train_df[, c("overall_survival", "deceased")]
test_survival_cols  <- OV_ME_test_df[, c("overall_survival", "deceased")]
sorted_feature_counts=read.csv("OV/ME/features_cv.csv", row.names = NULL)
OV_ME_test <- OV_ME_test_df[, intersect(colnames(OV_ME_test_df), sorted_feature_counts$Feature)]
OV_ME_train <- OV_ME_train_df[, intersect(colnames(OV_ME_train_df), sorted_feature_counts$Feature)]
OV_ME_train <- cbind(train_survival_cols, OV_ME_train)
OV_ME_test  <- cbind(test_survival_cols, OV_ME_test)
cv_result <- evaluate_performance(OV_ME_train, OV_ME_test)
boostrap_feature_counts <- perform_bootstrapping(OV_ME_train_df[,-1])
# Compute the 70th percentile threshold
threshold <- quantile(boostrap_feature_counts$count, probs = 0.70)
# Filter features that meet or exceed this threshold
top_features <- boostrap_feature_counts[boostrap_feature_counts$count >= threshold, ]
write.csv(top_features, file = "OV/ME/features_bs", row.names = FALSE)
OV_ME_test <- OV_ME_test_df[, intersect(colnames(OV_ME_test_df), top_features$feature)]
OV_ME_train <- OV_ME_train_df[, intersect(colnames(OV_ME_train_df), top_features$feature)]
OV_ME_train <- cbind(train_survival_cols, OV_ME_train)
OV_ME_test  <- cbind(test_survival_cols, OV_ME_test)
boostrap_result <- evaluate_performance(OV_ME_train, OV_ME_test)
no_feature_selection <- evaluate_performance(OV_ME_train_df[,-1], OV_ME_test_df[,-1]) 
write.csv(cv_result, file = "OV/ME/results_with_cv", row.names = FALSE)
write.csv(boostrap_result, file = "OV/ME/results_with_bs", row.names = FALSE)
write.csv(no_feature_selection, file = "OV/ME/results_without_fs", row.names = FALSE)

# GE ------------------------------------------------------------
OV_GE=read.csv("OV/OV_GE_clean.csv", row.names = NULL)
OV_GE$overall_survival[OV_GE$overall_survival <= 0] <- 0.001
OV_GE <- OV_GE[!is.na(OV_GE$overall_survival), ]
get_features <- colnames(OV_GE)[-c(1:3)]
train_test_data <- generate_train_test(OV_GE)
OV_GE_train_df <- train_test_data$train_df
OV_GE_test_df <- train_test_data$test_df
feature_counts_df <- data.frame(Feature = get_features, Count = 0)
perform_feature_selection_CV(OV_GE_train_df,feature_counts_df, "OV/GE")
train_survival_cols <- OV_GE_train_df[, c("overall_survival", "deceased")]
test_survival_cols  <- OV_GE_test_df[, c("overall_survival", "deceased")]
sorted_feature_counts=read.csv("OV/GE/features_cv.csv", row.names = NULL)
OV_GE_test <- OV_GE_test_df[, intersect(colnames(OV_GE_test_df), sorted_feature_counts$Feature)]
OV_GE_train <- OV_GE_train_df[, intersect(colnames(OV_GE_train_df), sorted_feature_counts$Feature)]
OV_GE_train <- cbind(train_survival_cols, OV_GE_train)
OV_GE_test  <- cbind(test_survival_cols, OV_GE_test)
cv_result <- evaluate_performance(OV_GE_train, OV_GE_test)
boostrap_feature_counts <- perform_bootstrapping(OV_GE_train_df[,-1])
# Compute the 70th percentile threshold
threshold <- quantile(boostrap_feature_counts$count, probs = 0.70)
# Filter features that meet or exceed this threshold
top_features <- boostrap_feature_counts[boostrap_feature_counts$count >= threshold, ]
write.csv(top_features, file = "OV/GE/features_bs", row.names = FALSE)
OV_GE_test <- OV_GE_test_df[, intersect(colnames(OV_GE_test_df), top_features$feature)]
OV_GE_train <- OV_GE_train_df[, intersect(colnames(OV_GE_train_df), top_features$feature)]
OV_GE_train <- cbind(train_survival_cols, OV_GE_train)
OV_GE_test  <- cbind(test_survival_cols, OV_GE_test)
boostrap_result <- evaluate_performance(OV_GE_train, OV_GE_test)
no_feature_selection <- evaluate_performance(OV_GE_train_df[,-1], OV_GE_test_df[,-1]) 
write.csv(cv_result, file = "OV/GE/results_with_cv", row.names = FALSE)
write.csv(boostrap_result, file = "OV/GE/results_with_bs", row.names = FALSE)
write.csv(no_feature_selection, file = "OV/GE/results_without_fs", row.names = FALSE)

# DM ------------------------------------------------------------
OV_DM=read.csv("OV/OV_DM_clean.csv", row.names = NULL)
OV_DM$overall_survival[OV_DM$overall_survival <= 0] <- 0.001
OV_DM <- OV_DM[!is.na(OV_DM$overall_survival), ]
get_features <- colnames(OV_DM)[-c(1:3)]
train_test_data <- generate_train_test(OV_DM)
OV_DM_train_df <- train_test_data$train_df
OV_DM_test_df <- train_test_data$test_df
feature_counts_df <- data.frame(Feature = get_features, Count = 0)
perform_feature_selection_CV(OV_DM_train_df,feature_counts_df, "OV/DM")
train_survival_cols <- OV_DM_train_df[, c("overall_survival", "deceased")]
test_survival_cols  <- OV_DM_test_df[, c("overall_survival", "deceased")]
sorted_feature_counts=read.csv("OV/DM/features_cv.csv", row.names = NULL)
OV_DM_test <- OV_DM_test_df[, intersect(colnames(OV_DM_test_df), sorted_feature_counts$Feature)]
OV_DM_train <- OV_DM_train_df[, intersect(colnames(OV_DM_train_df), sorted_feature_counts$Feature)]
OV_DM_train <- cbind(train_survival_cols, OV_DM_train)
OV_DM_test  <- cbind(test_survival_cols, OV_DM_test)
cv_result <- evaluate_performance(OV_DM_train, OV_DM_test)
boostrap_feature_counts <- perform_bootstrapping(OV_DM_train_df[,-1])
# Compute the 70th percentile threshold
threshold <- quantile(boostrap_feature_counts$count, probs = 0.70)
# Filter features that meet or exceed this threshold
top_features <- boostrap_feature_counts[boostrap_feature_counts$count >= threshold, ]
write.csv(top_features, file = "OV/DM/features_bs", row.names = FALSE)
OV_DM_test <- OV_DM_test_df[, intersect(colnames(OV_DM_test_df), top_features$feature)]
OV_DM_train <- OV_DM_train_df[, intersect(colnames(OV_DM_train_df), top_features$feature)]
OV_DM_train <- cbind(train_survival_cols, OV_DM_train)
OV_DM_test  <- cbind(test_survival_cols, OV_DM_test)
boostrap_result <- evaluate_performance(OV_DM_train, OV_DM_test)
no_feature_selection <- evaluate_performance(OV_DM_train_df[,-1], OV_DM_test_df[,-1]) 
write.csv(cv_result, file = "OV/DM/results_with_cv", row.names = FALSE)
write.csv(boostrap_result, file = "OV/DM/results_with_bs", row.names = FALSE)
write.csv(no_feature_selection, file = "OV/DM/results_without_fs", row.names = FALSE)

# CNV ------------------------------------------------------------
OV_CNV=read.csv("OV/OV_CNV_clean.csv", row.names = NULL)
OV_CNV$overall_survival[OV_CNV$overall_survival <= 0] <- 0.001
OV_CNV <- OV_CNV[!is.na(OV_CNV$overall_survival), ]
get_features <- colnames(OV_CNV)[-c(1:3)]
train_test_data <- generate_train_test(OV_CNV)
OV_CNV_train_df <- train_test_data$train_df
OV_CNV_test_df <- train_test_data$test_df
feature_counts_df <- data.frame(Feature = get_features, Count = 0)
perform_feature_selection_CV(OV_CNV_train_df,feature_counts_df, "OV/CNV")
train_survival_cols <- OV_CNV_train_df[, c("overall_survival", "deceased")]
test_survival_cols  <- OV_CNV_test_df[, c("overall_survival", "deceased")]
sorted_feature_counts=read.csv("OV/CNV/features_cv.csv", row.names = NULL)
OV_CNV_test <- OV_CNV_test_df[, intersect(colnames(OV_CNV_test_df), sorted_feature_counts$Feature)]
OV_CNV_train <- OV_CNV_train_df[, intersect(colnames(OV_CNV_train_df), sorted_feature_counts$Feature)]
OV_CNV_train <- cbind(train_survival_cols, OV_CNV_train)
OV_CNV_test  <- cbind(test_survival_cols, OV_CNV_test)
cv_result <- evaluate_performance(OV_CNV_train, OV_CNV_test)
boostrap_feature_counts <- perform_bootstrapping(OV_CNV_train_df[,-1])
# Compute the 70th percentile threshold
threshold <- quantile(boostrap_feature_counts$count, probs = 0.70)
# Filter features that meet or exceed this threshold
top_features <- boostrap_feature_counts[boostrap_feature_counts$count >= threshold, ]
write.csv(top_features, file = "OV/CNV/features_bs", row.names = FALSE)
OV_CNV_test <- OV_CNV_test_df[, intersect(colnames(OV_CNV_test_df), top_features$feature)]
OV_CNV_train <- OV_CNV_train_df[, intersect(colnames(OV_CNV_train_df), top_features$feature)]
OV_CNV_train <- cbind(train_survival_cols, OV_CNV_train)
OV_CNV_test  <- cbind(test_survival_cols, OV_CNV_test)
boostrap_result <- evaluate_performance(OV_CNV_train, OV_CNV_test)
no_feature_selection <- evaluate_performance(OV_CNV_train_df[,-1], OV_CNV_test_df[,-1]) 
write.csv(cv_result, file = "OV/CNV/results_with_cv", row.names = FALSE)
write.csv(boostrap_result, file = "OV/CNV/results_with_bs", row.names = FALSE)
write.csv(no_feature_selection, file = "OV/CNV/results_without_fs", row.names = FALSE)


####################################################################
#                           UCEC
#
####################################################################
# ME ------------------------------------------------------------
UCEC_ME=read.csv("UCEC/UCEC_ME_clean.csv", row.names = NULL)
UCEC_ME$overall_survival[UCEC_ME$overall_survival <= 0] <- 0.001
UCEC_ME <- UCEC_ME[!is.na(UCEC_ME$overall_survival), ]
get_features <- colnames(UCEC_ME)[-c(1:3)]
train_test_data <- generate_train_test(UCEC_ME)
UCEC_ME_train_df <- train_test_data$train_df
UCEC_ME_test_df <- train_test_data$test_df
feature_counts_df <- data.frame(Feature = get_features, Count = 0)
perform_feature_selection_CV(UCEC_ME_train_df,feature_counts_df, "UCEC/ME")
train_survival_cols <- UCEC_ME_train_df[, c("overall_survival", "deceased")]
test_survival_cols  <- UCEC_ME_test_df[, c("overall_survival", "deceased")]
sorted_feature_counts=read.csv("UCEC/ME/features_cv.csv", row.names = NULL)
UCEC_ME_test <- UCEC_ME_test_df[, intersect(colnames(UCEC_ME_test_df), sorted_feature_counts$Feature)]
UCEC_ME_train <- UCEC_ME_train_df[, intersect(colnames(UCEC_ME_train_df), sorted_feature_counts$Feature)]
UCEC_ME_train <- cbind(train_survival_cols, UCEC_ME_train)
UCEC_ME_test  <- cbind(test_survival_cols, UCEC_ME_test)
cv_result <- evaluate_performance(UCEC_ME_train, UCEC_ME_test)
boostrap_feature_counts <- perform_bootstrapping(UCEC_ME_train_df[,-1])
# Compute the 70th percentile threshold
threshold <- quantile(boostrap_feature_counts$count, probs = 0.70)
# Filter features that meet or exceed this threshold
top_features <- boostrap_feature_counts[boostrap_feature_counts$count >= threshold, ]
write.csv(top_features, file = "UCEC/ME/features_bs", row.names = FALSE)
UCEC_ME_test <- UCEC_ME_test_df[, intersect(colnames(UCEC_ME_test_df), top_features$feature)]
UCEC_ME_train <- UCEC_ME_train_df[, intersect(colnames(UCEC_ME_train_df), top_features$feature)]
UCEC_ME_train <- cbind(train_survival_cols, UCEC_ME_train)
UCEC_ME_test  <- cbind(test_survival_cols, UCEC_ME_test)
boostrap_result <- evaluate_performance(UCEC_ME_train, UCEC_ME_test)
no_feature_selection <- evaluate_performance(UCEC_ME_train_df[,-1], UCEC_ME_test_df[,-1]) 
write.csv(cv_result, file = "UCEC/ME/results_with_cv", row.names = FALSE)
write.csv(boostrap_result, file = "UCEC/ME/results_with_bs", row.names = FALSE)
write.csv(no_feature_selection, file = "UCEC/ME/results_without_fs", row.names = FALSE)

# GE ------------------------------------------------------------
UCEC_GE=read.csv("UCEC/UCEC_GE_clean.csv", row.names = NULL)
UCEC_GE$overall_survival[UCEC_GE$overall_survival <= 0] <- 0.001
UCEC_GE <- UCEC_GE[!is.na(UCEC_GE$overall_survival), ]
get_features <- colnames(UCEC_GE)[-c(1:3)]
train_test_data <- generate_train_test(UCEC_GE)
UCEC_GE_train_df <- train_test_data$train_df
UCEC_GE_test_df <- train_test_data$test_df
feature_counts_df <- data.frame(Feature = get_features, Count = 0)
perform_feature_selection_CV(UCEC_GE_train_df,feature_counts_df, "UCEC/GE")
train_survival_cols <- UCEC_GE_train_df[, c("overall_survival", "deceased")]
test_survival_cols  <- UCEC_GE_test_df[, c("overall_survival", "deceased")]
sorted_feature_counts=read.csv("UCEC/GE/features_cv.csv", row.names = NULL)
UCEC_GE_test <- UCEC_GE_test_df[, intersect(colnames(UCEC_GE_test_df), sorted_feature_counts$Feature)]
UCEC_GE_train <- UCEC_GE_train_df[, intersect(colnames(UCEC_GE_train_df), sorted_feature_counts$Feature)]
UCEC_GE_train <- cbind(train_survival_cols, UCEC_GE_train)
UCEC_GE_test  <- cbind(test_survival_cols, UCEC_GE_test)
cv_result <- evaluate_performance(UCEC_GE_train, UCEC_GE_test)
boostrap_feature_counts <- perform_bootstrapping(UCEC_GE_train_df[,-1])
# Compute the 70th percentile threshold
threshold <- quantile(boostrap_feature_counts$count, probs = 0.70)
# Filter features that meet or exceed this threshold
top_features <- boostrap_feature_counts[boostrap_feature_counts$count >= threshold, ]
write.csv(top_features, file = "UCEC/GE/features_bs", row.names = FALSE)
UCEC_GE_test <- UCEC_GE_test_df[, intersect(colnames(UCEC_GE_test_df), top_features$feature)]
UCEC_GE_train <- UCEC_GE_train_df[, intersect(colnames(UCEC_GE_train_df), top_features$feature)]
UCEC_GE_train <- cbind(train_survival_cols, UCEC_GE_train)
UCEC_GE_test  <- cbind(test_survival_cols, UCEC_GE_test)
boostrap_result <- evaluate_performance(UCEC_GE_train, UCEC_GE_test)
no_feature_selection <- evaluate_performance(UCEC_GE_train_df[,-1], UCEC_GE_test_df[,-1]) 
write.csv(cv_result, file = "UCEC/GE/results_with_cv", row.names = FALSE)
write.csv(boostrap_result, file = "UCEC/GE/results_with_bs", row.names = FALSE)
write.csv(no_feature_selection, file = "UCEC/GE/results_without_fs", row.names = FALSE)

# DM ------------------------------------------------------------
UCEC_DM=read.csv("UCEC/UCEC_DM_clean.csv", row.names = NULL)
UCEC_DM$overall_survival[UCEC_DM$overall_survival <= 0] <- 0.001
UCEC_DM <- UCEC_DM[!is.na(UCEC_DM$overall_survival), ]
get_features <- colnames(UCEC_DM)[-c(1:3)]
train_test_data <- generate_train_test(UCEC_DM)
UCEC_DM_train_df <- train_test_data$train_df
UCEC_DM_test_df <- train_test_data$test_df
feature_counts_df <- data.frame(Feature = get_features, Count = 0)
perform_feature_selection_CV(UCEC_DM_train_df,feature_counts_df, "UCEC/DM")
train_survival_cols <- UCEC_DM_train_df[, c("overall_survival", "deceased")]
test_survival_cols  <- UCEC_DM_test_df[, c("overall_survival", "deceased")]
sorted_feature_counts=read.csv("UCEC/DM/features_cv.csv", row.names = NULL)
UCEC_DM_test <- UCEC_DM_test_df[, intersect(colnames(UCEC_DM_test_df), sorted_feature_counts$Feature)]
UCEC_DM_train <- UCEC_DM_train_df[, intersect(colnames(UCEC_DM_train_df), sorted_feature_counts$Feature)]
UCEC_DM_train <- cbind(train_survival_cols, UCEC_DM_train)
UCEC_DM_test  <- cbind(test_survival_cols, UCEC_DM_test)
cv_result <- evaluate_performance(UCEC_DM_train, UCEC_DM_test)
boostrap_feature_counts <- perform_bootstrapping(UCEC_DM_train_df[,-1])
# Compute the 70th percentile threshold
threshold <- quantile(boostrap_feature_counts$count, probs = 0.70)
# Filter features that meet or exceed this threshold
top_features <- boostrap_feature_counts[boostrap_feature_counts$count >= threshold, ]
write.csv(top_features, file = "UCEC/DM/features_bs", row.names = FALSE)
UCEC_DM_test <- UCEC_DM_test_df[, intersect(colnames(UCEC_DM_test_df), top_features$feature)]
UCEC_DM_train <- UCEC_DM_train_df[, intersect(colnames(UCEC_DM_train_df), top_features$feature)]
UCEC_DM_train <- cbind(train_survival_cols, UCEC_DM_train)
UCEC_DM_test  <- cbind(test_survival_cols, UCEC_DM_test)
boostrap_result <- evaluate_performance(UCEC_DM_train, UCEC_DM_test)
no_feature_selection <- evaluate_performance(UCEC_DM_train_df[,-1], UCEC_DM_test_df[,-1]) 
write.csv(cv_result, file = "UCEC/DM/results_with_cv", row.names = FALSE)
write.csv(boostrap_result, file = "UCEC/DM/results_with_bs", row.names = FALSE)
write.csv(no_feature_selection, file = "UCEC/DM/results_without_fs", row.names = FALSE)

# CNV ------------------------------------------------------------
UCEC_CNV=read.csv("UCEC/UCEC_CNV_clean.csv", row.names = NULL)
UCEC_CNV$overall_survival[UCEC_CNV$overall_survival <= 0] <- 0.001
UCEC_CNV <- UCEC_CNV[!is.na(UCEC_CNV$overall_survival), ]
get_features <- colnames(UCEC_CNV)[-c(1:3)]
train_test_data <- generate_train_test(UCEC_CNV)
UCEC_CNV_train_df <- train_test_data$train_df
UCEC_CNV_test_df <- train_test_data$test_df
feature_counts_df <- data.frame(Feature = get_features, Count = 0)
perform_feature_selection_CV(UCEC_CNV_train_df,feature_counts_df, "UCEC/CNV")
train_survival_cols <- UCEC_CNV_train_df[, c("overall_survival", "deceased")]
test_survival_cols  <- UCEC_CNV_test_df[, c("overall_survival", "deceased")]
sorted_feature_counts=read.csv("UCEC/CNV/features_cv.csv", row.names = NULL)
UCEC_CNV_test <- UCEC_CNV_test_df[, intersect(colnames(UCEC_CNV_test_df), sorted_feature_counts$Feature)]
UCEC_CNV_train <- UCEC_CNV_train_df[, intersect(colnames(UCEC_CNV_train_df), sorted_feature_counts$Feature)]
UCEC_CNV_train <- cbind(train_survival_cols, UCEC_CNV_train)
UCEC_CNV_test  <- cbind(test_survival_cols, UCEC_CNV_test)
cv_result <- evaluate_performance(UCEC_CNV_train, UCEC_CNV_test)
boostrap_feature_counts <- perform_bootstrapping(UCEC_CNV_train_df[,-1])
# Compute the 70th percentile threshold
threshold <- quantile(boostrap_feature_counts$count, probs = 0.70)
# Filter features that meet or exceed this threshold
top_features <- boostrap_feature_counts[boostrap_feature_counts$count >= threshold, ]
write.csv(top_features, file = "UCEC/CNV/features_bs", row.names = FALSE)
UCEC_CNV_test <- UCEC_CNV_test_df[, intersect(colnames(UCEC_CNV_test_df), top_features$feature)]
UCEC_CNV_train <- UCEC_CNV_train_df[, intersect(colnames(UCEC_CNV_train_df), top_features$feature)]
UCEC_CNV_train <- cbind(train_survival_cols, UCEC_CNV_train)
UCEC_CNV_test  <- cbind(test_survival_cols, UCEC_CNV_test)
boostrap_result <- evaluate_performance(UCEC_CNV_train, UCEC_CNV_test)
no_feature_selection <- evaluate_performance(UCEC_CNV_train_df[,-1], UCEC_CNV_test_df[,-1]) 
write.csv(cv_result, file = "UCEC/CNV/results_with_cv", row.names = FALSE)
write.csv(boostrap_result, file = "UCEC/CNV/results_with_bs", row.names = FALSE)
write.csv(no_feature_selection, file = "UCEC/CNV/results_without_fs", row.names = FALSE)
