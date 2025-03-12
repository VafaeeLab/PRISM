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

source("functions.R")
# Set seed for reproducibility
set.seed(123)
# Number of repeats and folds
repeats <- 5
folds <- 5

# ============================================================
# Cross-validation Pipeline
# ============================================================

# This script processes multi-omics data (ME, DM, GE, CNV) for four cancer types: 
# BRCA, OV, CESC, and UCEC. The main objectives are:
# 1. Extract survival data and gene-level features.
# 2. Create train and test datasets for each omics type.
# 3. Track occurrences of gene-level features.
# 4. Perform feature selection using cross-validation (CV).
# 5. Evaluate selected features against models without feature selection.

# --------------------
# Workflow
# --------------------
# 1. Data Loading: Reads omics datasets for each cancer type.
# 2. Feature Extraction: Extracts survival data and gene-level features.
# 3. Data Splitting: Creates train and test sets for each omics dataset.
# 4. Feature Selection (perform_feature_selection_CV function):
#    - Applies four filter methods: 
#      - Multivariate CoxPH
#      - Random Forest and its extensions
#    - Features are selected based on:
#      - Non-zero coefficients
#      - P-values (CoxPH)
#      - Importance scores (Random Forest)
#    - Each method runs 5-fold, 5-repeat CV (100 runs total).
#    - Features appearing ≥50 times are retained.
# 5. Model Evaluation: 
#    - Selected features are evaluated on the training set.
#    - Performance is compared against models without feature selection.




####################################################################
#                           BRCA
#
####################################################################
# ME ------------------------------------------------------------
BRCA_ME=read.csv("BRCA/BRCA_ME_data.csv", row.names = NULL)
BRCA_ME$overall_survival[BRCA_ME$overall_survival <= 0] <- 0.001
survival_data <- Surv(BRCA_ME$overall_survival, BRCA_ME$deceased)
get_features <- colnames(BRCA_ME)[-c(1:4)]
train_test_data <- generate_train_test(BRCA_ME, survival_data)
BRCA_ME_train_df <- train_test_data$train_df
BRCA_ME_test_df <- train_test_data$test_df
BRCA_ME_train_survival <- train_test_data$train_survival
BRCA_ME_test_survival <- train_test_data$test_survival
feature_counts_df <- data.frame(Feature = get_features, Count = 0)
perform_feature_selection_CV(BRCA_ME_train_df,BRCA_ME_train_survival, feature_counts_df, "BRCA/ME")
sorted_feature_counts=read.csv("BRCA/ME/features_cv.csv", row.names = NULL)
BRCA_ME_test <- BRCA_ME_test_df[, intersect(colnames(BRCA_ME_test_df), sorted_feature_counts$Feature)]
BRCA_ME_train <- BRCA_ME_train_df[, intersect(colnames(BRCA_ME_train_df), sorted_feature_counts$Feature)]
result <- evaluate_performance(BRCA_ME_train, BRCA_ME_train_survival, BRCA_ME_test, BRCA_ME_test_survival)
no_feature_selection <- evaluate_performance(BRCA_ME_train_df[,-c(1:2)], BRCA_ME_train_survival, BRCA_ME_test_df[,-c(1:2)], BRCA_ME_test_survival) 
write.csv(result, file = "BRCA/ME/results_with_fs", row.names = FALSE)
write.csv(no_feature_selection, file = "BRCA/ME/results_without_fs", row.names = FALSE)

# DM ------------------------------------------------------------
BRCA_METH=read.csv("BRCA/BRCA_METH_data.csv", row.names = NULL)
BRCA_METH$overall_survival[BRCA_METH$overall_survival <= 0] <- 0.001
survival_data <- Surv(BRCA_METH$overall_survival, BRCA_METH$deceased)
get_features <- colnames(BRCA_METH)[-c(1:4)]
train_test_data <- generate_train_test(BRCA_METH, survival_data)
BRCA_METH_train_df <- train_test_data$train_df
BRCA_METH_test_df <- train_test_data$test_df
BRCA_METH_train_survival <- train_test_data$train_survival
BRCA_METH_test_survival <- train_test_data$test_survival
feature_counts_df <- data.frame(Feature = get_features, Count = 0)
perform_feature_selection_CV(BRCA_METH_train_df,BRCA_METH_train_survival, feature_counts_df, "BRCA/METH")
sorted_feature_counts=read.csv("BRCA/METH/features_cv.csv", row.names = NULL)
BRCA_METH_test <- BRCA_METH_test_df[, intersect(colnames(BRCA_METH_test_df), sorted_feature_counts$Feature)]
BRCA_METH_train <- BRCA_METH_train_df[, intersect(colnames(BRCA_METH_train_df), sorted_feature_counts$Feature)]
result <- evaluate_performance(BRCA_METH_train, BRCA_METH_train_survival, BRCA_METH_test, BRCA_METH_test_survival)
no_feature_selection <- evaluate_performance(BRCA_METH_train_df[,-c(1:2)], BRCA_METH_train_survival, BRCA_METH_test_df[,-c(1:2)], BRCA_METH_test_survival) 
write.csv(result, file = "BRCA/METH/results_with_fs", row.names = FALSE)
write.csv(no_feature_selection, file = "BRCA/METH/results_without_fs", row.names = FALSE)

# GE ------------------------------------------------------------
BRCA_GE=read.csv("BRCA/BRCA_GE_data.csv", row.names = NULL)
BRCA_GE$overall_survival[BRCA_GE$overall_survival <= 0] <- 0.001
survival_data <- Surv(BRCA_GE$overall_survival, BRCA_GE$deceased)
get_features <- colnames(BRCA_GE)[-c(1:4)]
train_test_data <- generate_train_test(BRCA_GE, survival_data)
BRCA_GE_train_df <- train_test_data$train_df
BRCA_GE_test_df <- train_test_data$test_df
BRCA_GE_train_survival <- train_test_data$train_survival
BRCA_GE_test_survival <- train_test_data$test_survival
feature_counts_df <- data.frame(Feature = get_features, Count = 0)
perform_feature_selection_CV(BRCA_GE_train_df,BRCA_GE_train_survival, feature_counts_df, "BRCA/GE")
sorted_feature_counts=read.csv("BRCA/GE/features_cv.csv", row.names = NULL)
BRCA_GE_test <- BRCA_GE_test_df[, intersect(colnames(BRCA_GE_test_df), sorted_feature_counts$Feature)]
BRCA_GE_train <- BRCA_GE_train_df[, intersect(colnames(BRCA_GE_train_df), sorted_feature_counts$Feature)]
result <- evaluate_performance(BRCA_GE_train, BRCA_GE_train_survival, BRCA_GE_test, BRCA_GE_test_survival)
no_feature_selection <- evaluate_performance(BRCA_GE_train_df[,-c(1:2)], BRCA_GE_train_survival, BRCA_GE_test_df[,-c(1:2)], BRCA_GE_test_survival) 
write.csv(result, file = "BRCA/GE/results_with_fs", row.names = FALSE)
write.csv(no_feature_selection, file = "BRCA/GE/results_without_fs", row.names = FALSE)

# CNV ------------------------------------------------------------
BRCA_CNV=read.csv("BRCA/BRCA_CNV_data.csv", row.names = NULL)
BRCA_CNV$overall_survival[BRCA_CNV$overall_survival <= 0] <- 0.001
survival_data <- Surv(BRCA_CNV$overall_survival, BRCA_CNV$deceased)
get_features <- colnames(BRCA_CNV)[-c(1:4)]
train_test_data <- generate_train_test(BRCA_CNV, survival_data)
BRCA_CNV_train_df <- train_test_data$train_df
BRCA_CNV_test_df <- train_test_data$test_df
BRCA_CNV_train_survival <- train_test_data$train_survival
BRCA_CNV_test_survival <- train_test_data$test_survival
feature_counts_df <- data.frame(Feature = get_features, Count = 0)
perform_feature_selection_CV(BRCA_CNV_train_df,BRCA_CNV_train_survival, feature_counts_df, "BRCA/CNV")
sorted_feature_counts=read.csv("BRCA/CNV/features_cv.csv", row.names = NULL)
BRCA_CNV_test <- BRCA_CNV_test_df[, intersect(colnames(BRCA_CNV_test_df), sorted_feature_counts$Feature)]
BRCA_CNV_train <- BRCA_CNV_train_df[, intersect(colnames(BRCA_CNV_train_df), sorted_feature_counts$Feature)]
result <- evaluate_performance(BRCA_CNV_train, BRCA_CNV_train_survival, BRCA_CNV_test, BRCA_CNV_test_survival)
no_feature_selection <- evaluate_performance(BRCA_CNV_train_df[,-c(1:2)], BRCA_CNV_train_survival, BRCA_CNV_test_df[,-c(1:2)], BRCA_CNV_test_survival) 
write.csv(result, file = "BRCA/CNV/results_with_fs", row.names = FALSE)
write.csv(no_feature_selection, file = "BRCA/CNV/results_without_fs", row.names = FALSE)

####################################################################
#                           OV
#
####################################################################

# ME ------------------------------------------------------------
OV_ME=read.csv("OV/OV_ME_data.csv", row.names = NULL)
OV_ME$overall_survival[OV_ME$overall_survival <= 0] <- 0.001
survival_data <- Surv(OV_ME$overall_survival, OV_ME$deceased)
get_features <- colnames(OV_ME)[-c(1:4)]
train_test_data <- generate_train_test(OV_ME, survival_data)
OV_ME_train_df <- train_test_data$train_df
OV_ME_test_df <- train_test_data$test_df
OV_ME_train_survival <- train_test_data$train_survival
OV_ME_test_survival <- train_test_data$test_survival
feature_counts_df <- data.frame(Feature = get_features, Count = 0)
perform_feature_selection_CV(OV_ME_train_df,OV_ME_train_survival, feature_counts_df, "OV/ME")
sorted_feature_counts=read.csv("OV/ME/features_cv.csv", row.names = NULL)
OV_ME_test <- OV_ME_test_df[, intersect(colnames(OV_ME_test_df), sorted_feature_counts$Feature)]
OV_ME_train <- OV_ME_train_df[, intersect(colnames(OV_ME_train_df), sorted_feature_counts$Feature)]
result <- evaluate_performance(OV_ME_train, OV_ME_train_survival, OV_ME_test, OV_ME_test_survival)
no_feature_selection <- evaluate_performance(OV_ME_train_df[,-c(1:2)], OV_ME_train_survival, OV_ME_test_df[,-c(1:2)], OV_ME_test_survival) 
write.csv(result, file = "OV/ME/results_with_fs", row.names = FALSE)
write.csv(no_feature_selection, file = "OV/ME/results_without_fs", row.names = FALSE)

# DM ------------------------------------------------------------
OV_METH=read.csv("OV/OV_METH_data.csv", row.names = NULL)
OV_METH$overall_survival[OV_METH$overall_survival <= 0] <- 0.001
survival_data <- Surv(OV_METH$overall_survival, OV_METH$deceased)
get_features <- colnames(OV_METH)[-c(1:4)]
train_test_data <- generate_train_test(OV_METH, survival_data)
OV_METH_train_df <- train_test_data$train_df
OV_METH_test_df <- train_test_data$test_df
OV_METH_train_survival <- train_test_data$train_survival
OV_METH_test_survival <- train_test_data$test_survival
feature_counts_df <- data.frame(Feature = get_features, Count = 0)
perform_feature_selection_CV(OV_METH_train_df,OV_METH_train_survival, feature_counts_df, "OV/METH")
sorted_feature_counts=read.csv("OV/METH/features_cv.csv", row.names = NULL)
OV_METH_test <- OV_METH_test_df[, intersect(colnames(OV_METH_test_df), sorted_feature_counts$Feature)]
OV_METH_train <- OV_METH_train_df[, intersect(colnames(OV_METH_train_df), sorted_feature_counts$Feature)]
result <- evaluate_performance(OV_METH_train, OV_METH_train_survival, OV_METH_test, OV_METH_test_survival)
no_feature_selection <- evaluate_performance(OV_METH_train_df[,-c(1:2)], OV_METH_train_survival, OV_METH_test_df[,-c(1:2)], OV_METH_test_survival) 
write.csv(result, file = "OV/METH/results_with_fs", row.names = FALSE)
write.csv(no_feature_selection, file = "OV/METH/results_without_fs", row.names = FALSE)

# GE ------------------------------------------------------------
OV_GE=read.csv("OV/OV_GE_data.csv", row.names = NULL)
OV_GE$overall_survival[OV_GE$overall_survival <= 0] <- 0.001
survival_data <- Surv(OV_GE$overall_survival, OV_GE$deceased)
get_features <- colnames(OV_GE)[-c(1:4)]
train_test_data <- generate_train_test(OV_GE, survival_data)
OV_GE_train_df <- train_test_data$train_df
OV_GE_test_df <- train_test_data$test_df
OV_GE_train_survival <- train_test_data$train_survival
OV_GE_test_survival <- train_test_data$test_survival
feature_counts_df <- data.frame(Feature = get_features, Count = 0)
perform_feature_selection_CV(OV_GE_train_df,OV_GE_train_survival, feature_counts_df, "OV/GE")
sorted_feature_counts=read.csv("OV/GE/features_cv.csv", row.names = NULL)
OV_GE_test <- OV_GE_test_df[, intersect(colnames(OV_GE_test_df), sorted_feature_counts$Feature)]
OV_GE_train <- OV_GE_train_df[, intersect(colnames(OV_GE_train_df), sorted_feature_counts$Feature)]
result <- evaluate_performance(OV_GE_train, OV_GE_train_survival, OV_GE_test, OV_GE_test_survival)
no_feature_selection <- evaluate_performance(OV_GE_train_df[,-c(1:2)], OV_GE_train_survival, OV_GE_test_df[,-c(1:2)], OV_GE_test_survival) 
write.csv(result, file = "OV/GE/results_with_fs", row.names = FALSE)
write.csv(no_feature_selection, file = "OV/GE/results_without_fs", row.names = FALSE)

# CNV ------------------------------------------------------------
OV_CNV=read.csv("OV/OV_CNV_data.csv", row.names = NULL)
OV_CNV$overall_survival[OV_CNV$overall_survival <= 0] <- 0.001
survival_data <- Surv(OV_CNV$overall_survival, OV_CNV$deceased)
get_features <- colnames(OV_CNV)[-c(1:4)]
train_test_data <- generate_train_test(OV_CNV, survival_data)
OV_CNV_train_df <- train_test_data$train_df
OV_CNV_test_df <- train_test_data$test_df
OV_CNV_train_survival <- train_test_data$train_survival
OV_CNV_test_survival <- train_test_data$test_survival
feature_counts_df <- data.frame(Feature = get_features, Count = 0)
perform_feature_selection_CV(OV_CNV_train_df,OV_CNV_train_survival, feature_counts_df, "OV/CNV")
sorted_feature_counts=read.csv("OV/CNV/features_cv.csv", row.names = NULL)
OV_CNV_test <- OV_CNV_test_df[, intersect(colnames(OV_CNV_test_df), sorted_feature_counts$Feature)]
OV_CNV_train <- OV_CNV_train_df[, intersect(colnames(OV_CNV_train_df), sorted_feature_counts$Feature)]
result <- evaluate_performance(OV_CNV_train, OV_CNV_train_survival, OV_CNV_test, OV_CNV_test_survival)
no_feature_selection <- evaluate_performance(OV_CNV_train_df[,-c(1:2)], OV_CNV_train_survival, OV_CNV_test_df[,-c(1:2)], OV_CNV_test_survival) 
write.csv(result, file = "OV/CNV/results_with_fs", row.names = FALSE)
write.csv(no_feature_selection, file = "OV/CNV/results_without_fs", row.names = FALSE)


####################################################################
#                           CESC
#
####################################################################
# ME ------------------------------------------------------------
CESC_ME=read.csv("CESC/CESC_ME_data.csv", row.names = NULL)
CESC_ME$overall_survival[CESC_ME$overall_survival <= 0] <- 0.001
survival_data <- Surv(CESC_ME$overall_survival, CESC_ME$deceased)
get_features <- colnames(CESC_ME)[-c(1:4)]
train_test_data <- generate_train_test(CESC_ME, survival_data)
CESC_ME_train_df <- train_test_data$train_df
CESC_ME_test_df <- train_test_data$test_df
CESC_ME_train_survival <- train_test_data$train_survival
CESC_ME_test_survival <- train_test_data$test_survival
feature_counts_df <- data.frame(Feature = get_features, Count = 0)
perform_feature_selection_CV(CESC_ME_train_df,CESC_ME_train_survival, feature_counts_df, "CESC/ME")
sorted_feature_counts=read.csv("CESC/ME/features_cv.csv", row.names = NULL)
CESC_ME_test <- CESC_ME_test_df[, intersect(colnames(CESC_ME_test_df), sorted_feature_counts$Feature)]
CESC_ME_train <- CESC_ME_train_df[, intersect(colnames(CESC_ME_train_df), sorted_feature_counts$Feature)]
result <- evaluate_performance(CESC_ME_train, CESC_ME_train_survival, CESC_ME_test, CESC_ME_test_survival)
no_feature_selection <- evaluate_performance(CESC_ME_train_df[,-c(1:2)], CESC_ME_train_survival, CESC_ME_test_df[,-c(1:2)], CESC_ME_test_survival) 
write.csv(result, file = "CESC/ME/results_with_fs", row.names = FALSE)
write.csv(no_feature_selection, file = "CESC/ME/results_without_fs", row.names = FALSE)

# DM ------------------------------------------------------------
CESC_METH=read.csv("CESC/CESC_METH_data.csv", row.names = NULL)
CESC_METH$overall_survival[CESC_METH$overall_survival <= 0] <- 0.001
survival_data <- Surv(CESC_METH$overall_survival, CESC_METH$deceased)
get_features <- colnames(CESC_METH)[-c(1:4)]
train_test_data <- generate_train_test(CESC_METH, survival_data)
CESC_METH_train_df <- train_test_data$train_df
CESC_METH_test_df <- train_test_data$test_df
CESC_METH_train_survival <- train_test_data$train_survival
CESC_METH_test_survival <- train_test_data$test_survival
feature_counts_df <- data.frame(Feature = get_features, Count = 0)
perform_feature_selection_CV(CESC_METH_train_df,CESC_METH_train_survival, feature_counts_df, "CESC/METH")
sorted_feature_counts=read.csv("CESC/METH/features_cv.csv", row.names = NULL)
CESC_METH_test <- CESC_METH_test_df[, intersect(colnames(CESC_METH_test_df), sorted_feature_counts$Feature)]
CESC_METH_train <- CESC_METH_train_df[, intersect(colnames(CESC_METH_train_df), sorted_feature_counts$Feature)]
result <- evaluate_performance(CESC_METH_train, CESC_METH_train_survival, CESC_METH_test, CESC_METH_test_survival)
no_feature_selection <- evaluate_performance(CESC_METH_train_df[,-c(1:2)], CESC_METH_train_survival, CESC_METH_test_df[,-c(1:2)], CESC_METH_test_survival) 
write.csv(result, file = "CESC/METH/results_with_fs", row.names = FALSE)
write.csv(no_feature_selection, file = "CESC/METH/results_without_fs", row.names = FALSE)

# GE ------------------------------------------------------------
CESC_GE=read.csv("CESC/CESC_GE_data.csv", row.names = NULL)
CESC_GE$overall_survival[CESC_GE$overall_survival <= 0] <- 0.001
survival_data <- Surv(CESC_GE$overall_survival, CESC_GE$deceased)
get_features <- colnames(CESC_GE)[-c(1:4)]
train_test_data <- generate_train_test(CESC_GE, survival_data)
CESC_GE_train_df <- train_test_data$train_df
CESC_GE_test_df <- train_test_data$test_df
CESC_GE_train_survival <- train_test_data$train_survival
CESC_GE_test_survival <- train_test_data$test_survival
feature_counts_df <- data.frame(Feature = get_features, Count = 0)
perform_feature_selection_CV(CESC_GE_train_df,CESC_GE_train_survival, feature_counts_df, "CESC/GE")
sorted_feature_counts=read.csv("CESC/GE/features_cv.csv", row.names = NULL)
CESC_GE_test <- CESC_GE_test_df[, intersect(colnames(CESC_GE_test_df), sorted_feature_counts$Feature)]
CESC_GE_train <- CESC_GE_train_df[, intersect(colnames(CESC_GE_train_df), sorted_feature_counts$Feature)]
result <- evaluate_performance(CESC_GE_train, CESC_GE_train_survival, CESC_GE_test, CESC_GE_test_survival)
no_feature_selection <- evaluate_performance(CESC_GE_train_df[,-c(1:2)], CESC_GE_train_survival, CESC_GE_test_df[,-c(1:2)], CESC_GE_test_survival) 
write.csv(result, file = "CESC/GE/results_with_fs", row.names = FALSE)
write.csv(no_feature_selection, file = "CESC/GE/results_without_fs", row.names = FALSE)

# CNV ------------------------------------------------------------
CESC_CNV=read.csv("CESC/CESC_CNV_data.csv", row.names = NULL)
CESC_CNV$overall_survival[CESC_CNV$overall_survival <= 0] <- 0.001
survival_data <- Surv(CESC_CNV$overall_survival, CESC_CNV$deceased)
get_features <- colnames(CESC_CNV)[-c(1:4)]
train_test_data <- generate_train_test(CESC_CNV, survival_data)
CESC_CNV_train_df <- train_test_data$train_df
CESC_CNV_test_df <- train_test_data$test_df
CESC_CNV_train_survival <- train_test_data$train_survival
CESC_CNV_test_survival <- train_test_data$test_survival
feature_counts_df <- data.frame(Feature = get_features, Count = 0)
perform_feature_selection_CV(CESC_CNV_train_df,CESC_CNV_train_survival, feature_counts_df, "CESC/CNV")
sorted_feature_counts=read.csv("CESC/CNV/features_cv.csv", row.names = NULL)
CESC_CNV_test <- CESC_CNV_test_df[, intersect(colnames(CESC_CNV_test_df), sorted_feature_counts$Feature)]
CESC_CNV_train <- CESC_CNV_train_df[, intersect(colnames(CESC_CNV_train_df), sorted_feature_counts$Feature)]
result <- evaluate_performance(CESC_CNV_train, CESC_CNV_train_survival, CESC_CNV_test, CESC_CNV_test_survival)
no_feature_selection <- evaluate_performance(CESC_CNV_train_df[,-c(1:2)], CESC_CNV_train_survival, CESC_CNV_test_df[,-c(1:2)], CESC_CNV_test_survival) 
write.csv(result, file = "CESC/CNV/results_with_fs", row.names = FALSE)
write.csv(no_feature_selection, file = "CESC/CNV/results_without_fs", row.names = FALSE)



####################################################################
#                           UCEC
#
####################################################################
# ME ------------------------------------------------------------
UCEC_ME=read.csv("UCEC/UCEC_ME_data.csv", row.names = NULL)
UCEC_ME$overall_survival[UCEC_ME$overall_survival <= 0] <- 0.001
survival_data <- Surv(UCEC_ME$overall_survival, UCEC_ME$deceased)
get_features <- colnames(UCEC_ME)[-c(1:4)]
train_test_data <- generate_train_test(UCEC_ME, survival_data)
UCEC_ME_train_df <- train_test_data$train_df
UCEC_ME_test_df <- train_test_data$test_df
UCEC_ME_train_survival <- train_test_data$train_survival
UCEC_ME_test_survival <- train_test_data$test_survival
feature_counts_df <- data.frame(Feature = get_features, Count = 0)
perform_feature_selection_CV(UCEC_ME_train_df,UCEC_ME_train_survival, feature_counts_df, "UCEC/ME")
sorted_feature_counts=read.csv("UCEC/ME/features_cv.csv", row.names = NULL)
UCEC_ME_test <- UCEC_ME_test_df[, intersect(colnames(UCEC_ME_test_df), sorted_feature_counts$Feature)]
UCEC_ME_train <- UCEC_ME_train_df[, intersect(colnames(UCEC_ME_train_df), sorted_feature_counts$Feature)]
result <- evaluate_performance(UCEC_ME_train, UCEC_ME_train_survival, UCEC_ME_test, UCEC_ME_test_survival)
no_feature_selection <- evaluate_performance(UCEC_ME_train_df[,-c(1:2)], UCEC_ME_train_survival, UCEC_ME_test_df[,-c(1:2)], UCEC_ME_test_survival) 
write.csv(result, file = "UCEC/ME/results_with_fs", row.names = FALSE)
write.csv(no_feature_selection, file = "UCEC/ME/results_without_fs", row.names = FALSE)

# DM ------------------------------------------------------------
UCEC_METH=read.csv("UCEC/UCEC_METH_data.csv", row.names = NULL)
UCEC_METH$overall_survival[UCEC_METH$overall_survival <= 0] <- 0.001
survival_data <- Surv(UCEC_METH$overall_survival, UCEC_METH$deceased)
get_features <- colnames(UCEC_METH)[-c(1:4)]
train_test_data <- generate_train_test(UCEC_METH, survival_data)
UCEC_METH_train_df <- train_test_data$train_df
UCEC_METH_test_df <- train_test_data$test_df
UCEC_METH_train_survival <- train_test_data$train_survival
UCEC_METH_test_survival <- train_test_data$test_survival
feature_counts_df <- data.frame(Feature = get_features, Count = 0)
perform_feature_selection_CV(UCEC_METH_train_df,UCEC_METH_train_survival, feature_counts_df, "UCEC/METH")
sorted_feature_counts=read.csv("UCEC/METH/features_cv.csv", row.names = NULL)
UCEC_METH_test <- UCEC_METH_test_df[, intersect(colnames(UCEC_METH_test_df), sorted_feature_counts$Feature)]
UCEC_METH_train <- UCEC_METH_train_df[, intersect(colnames(UCEC_METH_train_df), sorted_feature_counts$Feature)]
result <- evaluate_performance(UCEC_METH_train, UCEC_METH_train_survival, UCEC_METH_test, UCEC_METH_test_survival)
no_feature_selection <- evaluate_performance(UCEC_METH_train_df[,-c(1:2)], UCEC_METH_train_survival, UCEC_METH_test_df[,-c(1:2)], UCEC_METH_test_survival) 
write.csv(result, file = "UCEC/METH/results_with_fs", row.names = FALSE)
write.csv(no_feature_selection, file = "UCEC/METH/results_without_fs", row.names = FALSE)

# GE ------------------------------------------------------------
UCEC_GE=read.csv("UCEC/UCEC_GE_data.csv", row.names = NULL)
UCEC_GE$overall_survival[UCEC_GE$overall_survival <= 0] <- 0.001
survival_data <- Surv(UCEC_GE$overall_survival, UCEC_GE$deceased)
get_features <- colnames(UCEC_GE)[-c(1:4)]
train_test_data <- generate_train_test(UCEC_GE, survival_data)
UCEC_GE_train_df <- train_test_data$train_df
UCEC_GE_test_df <- train_test_data$test_df
UCEC_GE_train_survival <- train_test_data$train_survival
UCEC_GE_test_survival <- train_test_data$test_survival
feature_counts_df <- data.frame(Feature = get_features, Count = 0)
perform_feature_selection_CV(UCEC_GE_train_df,UCEC_GE_train_survival, feature_counts_df, "UCEC/GE")
sorted_feature_counts=read.csv("UCEC/GE/features_cv.csv", row.names = NULL)
UCEC_GE_test <- UCEC_GE_test_df[, intersect(colnames(UCEC_GE_test_df), sorted_feature_counts$Feature)]
UCEC_GE_train <- UCEC_GE_train_df[, intersect(colnames(UCEC_GE_train_df), sorted_feature_counts$Feature)]
result <- evaluate_performance(UCEC_GE_train, UCEC_GE_train_survival, UCEC_GE_test, UCEC_GE_test_survival)
no_feature_selection <- evaluate_performance(UCEC_GE_train_df[,-c(1:2)], UCEC_GE_train_survival, UCEC_GE_test_df[,-c(1:2)], UCEC_GE_test_survival) 
write.csv(result, file = "UCEC/GE/results_with_fs", row.names = FALSE)
write.csv(no_feature_selection, file = "UCEC/GE/results_without_fs", row.names = FALSE)

# CNV ------------------------------------------------------------
UCEC_CNV=read.csv("UCEC/UCEC_CNV_data.csv", row.names = NULL)
UCEC_CNV$overall_survival[UCEC_CNV$overall_survival <= 0] <- 0.001
survival_data <- Surv(UCEC_CNV$overall_survival, UCEC_CNV$deceased)
get_features <- colnames(UCEC_CNV)[-c(1:4)]
train_test_data <- generate_train_test(UCEC_CNV, survival_data)
UCEC_CNV_train_df <- train_test_data$train_df
UCEC_CNV_test_df <- train_test_data$test_df
UCEC_CNV_train_survival <- train_test_data$train_survival
UCEC_CNV_test_survival <- train_test_data$test_survival
feature_counts_df <- data.frame(Feature = get_features, Count = 0)
perform_feature_selection_CV(UCEC_CNV_train_df,UCEC_CNV_train_survival, feature_counts_df, "UCEC/CNV")
sorted_feature_counts=read.csv("UCEC/CNV/features_cv.csv", row.names = NULL)
UCEC_CNV_test <- UCEC_CNV_test_df[, intersect(colnames(UCEC_CNV_test_df), sorted_feature_counts$Feature)]
UCEC_CNV_train <- UCEC_CNV_train_df[, intersect(colnames(UCEC_CNV_train_df), sorted_feature_counts$Feature)]
result <- evaluate_performance(UCEC_CNV_train, UCEC_CNV_train_survival, UCEC_CNV_test, UCEC_CNV_test_survival)
no_feature_selection <- evaluate_performance(UCEC_CNV_train_df[,-c(1:2)], UCEC_CNV_train_survival, UCEC_CNV_test_df[,-c(1:2)], UCEC_CNV_test_survival) 
write.csv(result, file = "UCEC/CNV/results_with_fs", row.names = FALSE)
write.csv(no_feature_selection, file = "UCEC/CNV/results_without_fs", row.names = FALSE)



