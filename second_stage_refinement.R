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
library(dplyr)


source("functions.R")
set.seed(123)

# *********************************************************************
# Second Stage Refinement R Script Documentation
# *********************************************************************

# Overview:
# The Second Stage Refinement script processes the omics data by:
# 1. Reading the omics features selected by the cross-validation (CV) pipeline.
# 2. Min-Max normalizing the omics data.
# 3. Performing RFE on each of the single omics data - First Stage Refinement 
# 4. Once we have a refined set, we concat each possible omics combination and run it through RFE again - Second Stage Refinement
# This script prepares multi-omics datasets for further downstream analysis or modeling.

# *********************************************************************

####################################################################
#                           BRCA
#
####################################################################

ME_features=read.csv("BRCA/ME/features_cv.csv", row.names = NULL)
GE_features=read.csv("BRCA/GE/features_cv.csv", row.names = NULL)
METH_features=read.csv("BRCA/METH/features_cv.csv", row.names = NULL)
CNV_features=read.csv("BRCA/CNV/features_cv.csv", row.names = NULL)

BRCA_ME=read.csv("BRCA/BRCA_ME_data.csv", row.names = NULL)
BRCA_GE=read.csv("BRCA/BRCA_GE_data.csv", row.names = NULL)
BRCA_METH=read.csv("BRCA/BRCA_METH_data.csv", row.names = NULL)
BRCA_CNV=read.csv("BRCA/BRCA_CNV_data.csv", row.names = NULL)

BRCA_METH <- BRCA_METH[,-1]
BRCA_ME <- BRCA_ME[,-1]
BRCA_GE <- BRCA_GE[,-1]
BRCA_CNV <- BRCA_CNV[,-1]

# Preprocess-----------------------------------------------------------------------

BRCA_ME_filtered <- BRCA_ME[, intersect(colnames(BRCA_ME), ME_features$Feature)]
BRCA_ME_filtered <- cbind(BRCA_ME[, 1:3], BRCA_ME_filtered)
colnames(BRCA_ME_filtered)[-c(1:3)] <- paste("ME_", colnames(BRCA_ME_filtered)[-c(1:3)], sep = "")

BRCA_GE_filtered <- BRCA_GE[, intersect(colnames(BRCA_GE), GE_features$Feature)]
BRCA_GE_filtered <- cbind(BRCA_GE[, 1:3], BRCA_GE_filtered)
colnames(BRCA_GE_filtered)[-c(1:3)] <- paste("GE_", colnames(BRCA_GE_filtered)[-c(1:3)], sep = "")

BRCA_METH_filtered <- BRCA_METH[, intersect(colnames(BRCA_METH), METH_features$Feature)]
BRCA_METH_filtered <- cbind(BRCA_METH[, 1:3], BRCA_METH_filtered)
colnames(BRCA_METH_filtered)[-c(1:3)] <- paste("METH_", colnames(BRCA_METH_filtered)[-c(1:3)], sep = "")

BRCA_CNV_filtered <- BRCA_CNV[, intersect(colnames(BRCA_CNV), CNV_features$Feature)]
BRCA_CNV_filtered <- cbind(BRCA_CNV[, 1:3], BRCA_CNV_filtered)
colnames(BRCA_CNV_filtered)[-c(1:3)] <- paste("CNV_", colnames(BRCA_CNV_filtered)[-c(1:3)], sep = "")

# Apply min-max normalization to features (excluding first three columns)
BRCA_ME_normalized <- min_max_normalize(BRCA_ME_filtered)
BRCA_GE_normalized <- min_max_normalize(BRCA_GE_filtered)
BRCA_CNV_normalized <- min_max_normalize(BRCA_CNV_filtered)

# First Stage Refinement-------------------------------------------------------------------------------
# ME 
BRCA_ME_normalized$overall_survival[BRCA_ME_normalized$overall_survival <= 0] <- 0.001
BRCA_ME_normalized_survival_data <- Surv(BRCA_ME_normalized$overall_survival, BRCA_ME_normalized$deceased)
fusion(BRCA_ME_normalized, BRCA_ME_normalized_survival_data, directory = "BRCA/LF", file_prefix = "BRCA_ME_Fusion")

#-------------------------------------------------------------------------------
# GE 
BRCA_GE_normalized$overall_survival[BRCA_GE_normalized$overall_survival <= 0] <- 0.001
BRCA_GE_normalized_survival_data <- Surv(BRCA_GE_normalized$overall_survival, BRCA_GE_normalized$deceased)
fusion(BRCA_GE_normalized, BRCA_GE_normalized_survival_data, directory = "BRCA/LF", file_prefix = "BRCA_GE_Fusion")

#-------------------------------------------------------------------------------
# METH 
BRCA_METH_filtered$overall_survival[BRCA_METH_filtered$overall_survival <= 0] <- 0.001
BRCA_METH_filtered_survival_data <- Surv(BRCA_METH_filtered$overall_survival, BRCA_METH_filtered$deceased)
fusion(BRCA_METH_filtered, BRCA_METH_filtered_survival_data, directory = "BRCA/LF", file_prefix = "BRCA_METH_Fusion")

#-------------------------------------------------------------------------------
# CNV 
BRCA_CNV_normalized$overall_survival[BRCA_CNV_normalized$overall_survival <= 0] <- 0.001
BRCA_CNV_normalized_survival_data <- Surv(BRCA_CNV_normalized$overall_survival, BRCA_CNV_normalized$deceased)
fusion(BRCA_CNV_normalized, BRCA_CNV_normalized_survival_data, directory = "BRCA/LF", file_prefix = "BRCA_CNV_Fusion")


ME_fusion=read.csv("BRCA/LF/BRCA_ME_Fusion_data.csv", row.names = NULL)
METH_fusion=read.csv("BRCA/LF/BRCA_METH_Fusion_data.csv", row.names = NULL)
GE_fusion=read.csv("BRCA/LF/BRCA_GE_Fusion_data.csv", row.names = NULL)
CNV_fusion=read.csv("BRCA/LF/BRCA_CNV_Fusion_data.csv", row.names = NULL)

# Preprocess-----------------------------------------------------------------------

BRCA_ME_filtered <- BRCA_ME_filtered[, names(BRCA_ME_filtered) %in% ME_fusion$variables, drop = FALSE]
BRCA_ME_filtered <- cbind(BRCA_ME[, 1:3], BRCA_ME_filtered)
BRCA_GE_filtered <- BRCA_GE_filtered[, names(BRCA_GE_filtered) %in% GE_fusion$variables, drop = FALSE]
BRCA_GE_filtered <- cbind(BRCA_GE[, 1:3], BRCA_GE_filtered)
BRCA_METH_filtered <- BRCA_METH_filtered[, names(BRCA_METH_filtered) %in% METH_fusion$variables, drop = FALSE]
BRCA_METH_filtered <- cbind(BRCA_METH[, 1:3], BRCA_METH_filtered)
BRCA_CNV_filtered <- BRCA_CNV_filtered[, names(BRCA_CNV_filtered) %in% CNV_fusion$variables, drop = FALSE]
BRCA_CNV_filtered <- cbind(BRCA_CNV[, 1:3], BRCA_CNV_filtered)

# Apply min-max normalization to features (excluding first three columns)
BRCA_ME_normalized <- min_max_normalize(BRCA_ME_filtered)
BRCA_GE_normalized <- min_max_normalize(BRCA_GE_filtered)
BRCA_CNV_normalized <- min_max_normalize(BRCA_CNV_filtered)

# Stage 2 Refinement-------------------------------------------------------------------
# ME AND GE
BRCA_ME_GE <- merge(BRCA_ME_normalized, BRCA_GE_normalized, by = c("case_id", "deceased", "overall_survival"))
BRCA_ME_GE$overall_survival[BRCA_ME_GE$overall_survival <= 0] <- 0.001
BRCA_ME_GE_survival_data <- Surv(BRCA_ME_GE$overall_survival, BRCA_ME_GE$deceased)
fusion(BRCA_ME_GE, BRCA_ME_GE_survival_data, directory = "BRCA/LF", file_prefix = "BRCA_ME_GE")
# ME AND METH
BRCA_ME_METH <- merge(BRCA_ME_normalized, BRCA_METH_filtered, by = c("case_id", "deceased", "overall_survival"))
BRCA_ME_METH$overall_survival[BRCA_ME_METH$overall_survival <= 0] <- 0.001
BRCA_ME_METH_survival_data <- Surv(BRCA_ME_METH$overall_survival, BRCA_ME_METH$deceased)
fusion(BRCA_ME_METH, BRCA_ME_METH_survival_data, directory = "BRCA/LF", file_prefix = "BRCA_ME_METH")
# METH AND GE
BRCA_METH_GE <- merge(BRCA_METH_filtered, BRCA_GE_normalized, by = c("case_id", "deceased", "overall_survival"))
BRCA_METH_GE$overall_survival[BRCA_METH_GE$overall_survival <= 0] <- 0.001
BRCA_METH_GE_survival_data <- Surv(BRCA_METH_GE$overall_survival, BRCA_METH_GE$deceased)
fusion(BRCA_METH_GE, BRCA_METH_GE_survival_data, directory = "BRCA/LF", file_prefix = "BRCA_METH_GE")
# GE AND CNV
BRCA_GE_CNV <- merge(BRCA_GE_normalized, BRCA_CNV_normalized, by = c("case_id", "deceased", "overall_survival"))
BRCA_GE_CNV$overall_survival[BRCA_GE_CNV$overall_survival <= 0] <- 0.001
BRCA_GE_CNV_survival_data <- Surv(BRCA_GE_CNV$overall_survival, BRCA_GE_CNV$deceased)
fusion(BRCA_GE_CNV, BRCA_GE_CNV_survival_data, directory = "BRCA/LF", file_prefix = "BRCA_GE_CNV")
# ME AND CNV
BRCA_ME_CNV <- merge(BRCA_ME_normalized, BRCA_CNV_normalized, by = c("case_id", "deceased", "overall_survival"))
BRCA_ME_CNV$overall_survival[BRCA_ME_CNV$overall_survival <= 0] <- 0.001
BRCA_ME_CNV_survival_data <- Surv(BRCA_ME_CNV$overall_survival, BRCA_ME_CNV$deceased)
fusion(BRCA_ME_CNV, BRCA_ME_CNV_survival_data, directory = "BRCA/LF", file_prefix = "BRCA_ME_CNV")
# METH AND CNV ------------------------
BRCA_METH_CNV <- merge(BRCA_METH_filtered, BRCA_CNV_normalized, by = c("case_id", "deceased", "overall_survival"))
BRCA_METH_CNV$overall_survival[BRCA_METH_CNV$overall_survival <= 0] <- 0.001
BRCA_METH_CNV_survival_data <- Surv(BRCA_METH_CNV$overall_survival, BRCA_METH_CNV$deceased)
fusion(BRCA_METH_CNV, BRCA_METH_CNV_survival_data, directory = "BRCA/LF", file_prefix = "BRCA_METH_CNV")
# METH AND GE AND ME
BRCA_METH_GE_ME <- merge(BRCA_METH_filtered, BRCA_GE_normalized, by = c("case_id", "deceased", "overall_survival"))
BRCA_METH_GE_ME <- merge(BRCA_METH_GE_ME, BRCA_ME_normalized, by = c("case_id", "deceased", "overall_survival"))
BRCA_METH_GE_ME $overall_survival[BRCA_METH_GE_ME$overall_survival <= 0] <- 0.001
BRCA_METH_GE_ME_survival_data <- Surv(BRCA_METH_GE_ME$overall_survival, BRCA_METH_GE_ME$deceased)
fusion(BRCA_METH_GE_ME, BRCA_METH_GE_ME_survival_data, directory = "BRCA/LF", file_prefix = "BRCA_METH_GE_ME")
# METH AND GE AND CNV
BRCA_METH_GE_CNV <- merge(BRCA_METH_filtered, BRCA_GE_normalized, by = c("case_id", "deceased", "overall_survival"))
BRCA_METH_GE_CNV <- merge(BRCA_METH_GE_CNV, BRCA_CNV_normalized, by = c("case_id", "deceased", "overall_survival"))
BRCA_METH_GE_CNV $overall_survival[BRCA_METH_GE_CNV$overall_survival <= 0] <- 0.001
BRCA_METH_GE_CNV_survival_data <- Surv(BRCA_METH_GE_CNV$overall_survival, BRCA_METH_GE_CNV$deceased)
fusion(BRCA_METH_GE_CNV, BRCA_METH_GE_CNV_survival_data, directory = "BRCA/LF", file_prefix = "BRCA_METH_GE_CNV")
# METH AND ME AND CNV
BRCA_METH_ME_CNV <- merge(BRCA_METH_filtered, BRCA_ME_normalized, by = c("case_id", "deceased", "overall_survival"))
BRCA_METH_ME_CNV <- merge(BRCA_METH_ME_CNV, BRCA_CNV_normalized, by = c("case_id", "deceased", "overall_survival"))
BRCA_METH_ME_CNV $overall_survival[BRCA_METH_ME_CNV$overall_survival <= 0] <- 0.001
BRCA_METH_ME_CNV_survival_data <- Surv(BRCA_METH_ME_CNV$overall_survival, BRCA_METH_ME_CNV$deceased)
fusion(BRCA_METH_ME_CNV, BRCA_METH_ME_CNV_survival_data, directory = "BRCA/LF", file_prefix = "BRCA_METH_ME_CNV")
# ME AND GE AND CNV
BRCA_ME_GE_CNV <- merge(BRCA_ME_normalized, BRCA_GE_normalized, by = c("case_id", "deceased", "overall_survival"))
BRCA_ME_GE_CNV <- merge(BRCA_ME_GE_CNV, BRCA_CNV_normalized, by = c("case_id", "deceased", "overall_survival"))
BRCA_ME_GE_CNV $overall_survival[BRCA_ME_GE_CNV$overall_survival <= 0] <- 0.001
BRCA_ME_GE_CNV_survival_data <- Surv(BRCA_ME_GE_CNV$overall_survival, BRCA_ME_GE_CNV$deceased)
fusion(BRCA_ME_GE_CNV, BRCA_ME_GE_CNV_survival_data, directory = "BRCA/LF", file_prefix = "BRCA_ME_GE_CNV")
# ME AND GE AND CNV AND METH
BRCA_ME_GE_CNV_METH <- merge(BRCA_ME_normalized, BRCA_GE_normalized, by = c("case_id", "deceased", "overall_survival"))
BRCA_ME_GE_CNV_METH <- merge(BRCA_ME_GE_CNV_METH, BRCA_CNV_normalized, by = c("case_id", "deceased", "overall_survival"))
BRCA_ME_GE_CNV_METH <- merge(BRCA_ME_GE_CNV_METH, BRCA_METH_filtered, by = c("case_id", "deceased", "overall_survival"))
BRCA_ME_GE_CNV_METH $overall_survival[BRCA_ME_GE_CNV_METH$overall_survival <= 0] <- 0.001
BRCA_ME_GE_CNV_METH_survival_data <- Surv(BRCA_ME_GE_CNV_METH$overall_survival, BRCA_ME_GE_CNV_METH$deceased)
fusion(BRCA_ME_GE_CNV_METH, BRCA_ME_GE_CNV_METH_survival_data, directory = "BRCA/LF", file_prefix = "BRCA_ME_GE_CNV_METH")



####################################################################
#                           OV
#
####################################################################

ME_features=read.csv("OV/ME/features_cv.csv", row.names = NULL)
GE_features=read.csv("OV/GE/features_cv.csv", row.names = NULL)
METH_features=read.csv("OV/METH/features_cv.csv", row.names = NULL)
CNV_features=read.csv("OV/CNV/features_cv.csv", row.names = NULL)

OV_ME=read.csv("OV/OV_ME_data.csv", row.names = NULL)
OV_GE=read.csv("OV/OV_GE_data.csv", row.names = NULL)
OV_METH=read.csv("OV/OV_METH_data.csv", row.names = NULL)
OV_CNV=read.csv("OV/OV_CNV_data.csv", row.names = NULL)

OV_METH <- OV_METH[,-1]
OV_ME <- OV_ME[,-1]
OV_GE <- OV_GE[,-1]
OV_CNV <- OV_CNV[,-1]

# Preprocess-----------------------------------------------------------------------

OV_ME_filtered <- OV_ME[, intersect(colnames(OV_ME), ME_features$Feature)]
OV_ME_filtered <- cbind(OV_ME[, 1:3], OV_ME_filtered)
colnames(OV_ME_filtered)[-c(1:3)] <- paste("ME_", colnames(OV_ME_filtered)[-c(1:3)], sep = "")

OV_GE_filtered <- OV_GE[, intersect(colnames(OV_GE), GE_features$Feature)]
OV_GE_filtered <- cbind(OV_GE[, 1:3], OV_GE_filtered)
colnames(OV_GE_filtered)[-c(1:3)] <- paste("GE_", colnames(OV_GE_filtered)[-c(1:3)], sep = "")

OV_METH_filtered <- OV_METH[, intersect(colnames(OV_METH), METH_features$Feature)]
OV_METH_filtered <- cbind(OV_METH[, 1:3], OV_METH_filtered)
colnames(OV_METH_filtered)[-c(1:3)] <- paste("METH_", colnames(OV_METH_filtered)[-c(1:3)], sep = "")

OV_CNV_filtered <- OV_CNV[, intersect(colnames(OV_CNV), CNV_features$Feature)]
OV_CNV_filtered <- cbind(OV_CNV[, 1:3], OV_CNV_filtered)
colnames(OV_CNV_filtered)[-c(1:3)] <- paste("CNV_", colnames(OV_CNV_filtered)[-c(1:3)], sep = "")

# Apply min-max normalization to features (excluding first three columns)
OV_ME_normalized <- min_max_normalize(OV_ME_filtered)
OV_GE_normalized <- min_max_normalize(OV_GE_filtered)
OV_CNV_normalized <- min_max_normalize(OV_CNV_filtered)

#-------------------------------------------------------------------------------
# ME 
OV_ME_normalized$overall_survival[OV_ME_normalized$overall_survival <= 0] <- 0.001
OV_ME_normalized_survival_data <- Surv(OV_ME_normalized$overall_survival, OV_ME_normalized$deceased)
fusion(OV_ME_normalized, OV_ME_normalized_survival_data, directory = "OV/LF", file_prefix = "OV_ME_Fusion")

#-------------------------------------------------------------------------------
# GE 
OV_GE_normalized$overall_survival[OV_GE_normalized$overall_survival <= 0] <- 0.001
OV_GE_normalized_survival_data <- Surv(OV_GE_normalized$overall_survival, OV_GE_normalized$deceased)
fusion(OV_GE_normalized, OV_GE_normalized_survival_data, directory = "OV/LF", file_prefix = "OV_GE_Fusion")

#-------------------------------------------------------------------------------
# METH 
OV_METH_filtered$overall_survival[OV_METH_filtered$overall_survival <= 0] <- 0.001
OV_METH_filtered_survival_data <- Surv(OV_METH_filtered$overall_survival, OV_METH_filtered$deceased)
fusion(OV_METH_filtered, OV_METH_filtered_survival_data, directory = "OV/LF", file_prefix = "OV_METH_Fusion")

#-------------------------------------------------------------------------------
# CNV 
OV_CNV_normalized$overall_survival[OV_CNV_normalized$overall_survival <= 0] <- 0.001
OV_CNV_normalized_survival_data <- Surv(OV_CNV_normalized$overall_survival, OV_CNV_normalized$deceased)
fusion(OV_CNV_normalized, OV_CNV_normalized_survival_data, directory = "OV/LF", file_prefix = "OV_CNV_Fusion")


ME_fusion=read.csv("OV/LF/OV_ME_Fusion_data.csv", row.names = NULL)
METH_fusion=read.csv("OV/LF/OV_METH_Fusion_data.csv", row.names = NULL)
GE_fusion=read.csv("OV/LF/OV_GE_Fusion_data.csv", row.names = NULL)
CNV_fusion=read.csv("OV/LF/OV_CNV_Fusion_data.csv", row.names = NULL)

# Preprocess-----------------------------------------------------------------------

OV_ME_filtered <- OV_ME_filtered[, names(OV_ME_filtered) %in% ME_fusion$variables, drop = FALSE]
OV_ME_filtered <- cbind(OV_ME[, 1:3], OV_ME_filtered)
OV_GE_filtered <- OV_GE_filtered[, names(OV_GE_filtered) %in% GE_fusion$variables, drop = FALSE]
OV_GE_filtered <- cbind(OV_GE[, 1:3], OV_GE_filtered)
OV_METH_filtered <- OV_METH_filtered[, names(OV_METH_filtered) %in% METH_fusion$variables, drop = FALSE]
OV_METH_filtered <- cbind(OV_METH[, 1:3], OV_METH_filtered)
OV_CNV_filtered <- OV_CNV_filtered[, names(OV_CNV_filtered) %in% CNV_fusion$variables, drop = FALSE]
OV_CNV_filtered <- cbind(OV_CNV[, 1:3], OV_CNV_filtered)

# Apply min-max normalization to features (excluding first three columns)
OV_ME_normalized <- min_max_normalize(OV_ME_filtered)
OV_GE_normalized <- min_max_normalize(OV_GE_filtered)
OV_CNV_normalized <- min_max_normalize(OV_CNV_filtered)

# Late Fusion -------------------------------------------------------------------
# ME AND GE
OV_ME_GE <- merge(OV_ME_normalized, OV_GE_normalized, by = c("case_id", "deceased", "overall_survival"))
OV_ME_GE$overall_survival[OV_ME_GE$overall_survival <= 0] <- 0.001
OV_ME_GE_survival_data <- Surv(OV_ME_GE$overall_survival, OV_ME_GE$deceased)
fusion(OV_ME_GE, OV_ME_GE_survival_data, directory = "OV/LF", file_prefix = "OV_ME_GE")
# ME AND METH
OV_ME_METH <- merge(OV_ME_normalized, OV_METH_filtered, by = c("case_id", "deceased", "overall_survival"))
OV_ME_METH$overall_survival[OV_ME_METH$overall_survival <= 0] <- 0.001
OV_ME_METH_survival_data <- Surv(OV_ME_METH$overall_survival, OV_ME_METH$deceased)
fusion(OV_ME_METH, OV_ME_METH_survival_data, directory = "OV/LF", file_prefix = "OV_ME_METH")
# METH AND GE
OV_METH_GE <- merge(OV_METH_filtered, OV_GE_normalized, by = c("case_id", "deceased", "overall_survival"))
OV_METH_GE$overall_survival[OV_METH_GE$overall_survival <= 0] <- 0.001
OV_METH_GE_survival_data <- Surv(OV_METH_GE$overall_survival, OV_METH_GE$deceased)
fusion(OV_METH_GE, OV_METH_GE_survival_data, directory = "OV/LF", file_prefix = "OV_METH_GE")
# GE AND CNV
OV_GE_CNV <- merge(OV_GE_normalized, OV_CNV_normalized, by = c("case_id", "deceased", "overall_survival"))
OV_GE_CNV$overall_survival[OV_GE_CNV$overall_survival <= 0] <- 0.001
OV_GE_CNV_survival_data <- Surv(OV_GE_CNV$overall_survival, OV_GE_CNV$deceased)
fusion(OV_GE_CNV, OV_GE_CNV_survival_data, directory = "OV/LF", file_prefix = "OV_GE_CNV")
# ME AND CNV
OV_ME_CNV <- merge(OV_ME_normalized, OV_CNV_normalized, by = c("case_id", "deceased", "overall_survival"))
OV_ME_CNV$overall_survival[OV_ME_CNV$overall_survival <= 0] <- 0.001
OV_ME_CNV_survival_data <- Surv(OV_ME_CNV$overall_survival, OV_ME_CNV$deceased)
fusion(OV_ME_CNV, OV_ME_CNV_survival_data, directory = "OV/LF", file_prefix = "OV_ME_CNV")
# METH AND CNV ------------------------
OV_METH_CNV <- merge(OV_METH_filtered, OV_CNV_normalized, by = c("case_id", "deceased", "overall_survival"))
OV_METH_CNV$overall_survival[OV_METH_CNV$overall_survival <= 0] <- 0.001
OV_METH_CNV_survival_data <- Surv(OV_METH_CNV$overall_survival, OV_METH_CNV$deceased)
fusion(OV_METH_CNV, OV_METH_CNV_survival_data, directory = "OV/LF", file_prefix = "OV_METH_CNV")
# METH AND GE AND ME
OV_METH_GE_ME <- merge(OV_METH_filtered, OV_GE_normalized, by = c("case_id", "deceased", "overall_survival"))
OV_METH_GE_ME <- merge(OV_METH_GE_ME, OV_ME_normalized, by = c("case_id", "deceased", "overall_survival"))
OV_METH_GE_ME $overall_survival[OV_METH_GE_ME$overall_survival <= 0] <- 0.001
OV_METH_GE_ME_survival_data <- Surv(OV_METH_GE_ME$overall_survival, OV_METH_GE_ME$deceased)
fusion(OV_METH_GE_ME, OV_METH_GE_ME_survival_data, directory = "OV/LF", file_prefix = "OV_METH_GE_ME")
# METH AND GE AND CNV
OV_METH_GE_CNV <- merge(OV_METH_filtered, OV_GE_normalized, by = c("case_id", "deceased", "overall_survival"))
OV_METH_GE_CNV <- merge(OV_METH_GE_CNV, OV_CNV_normalized, by = c("case_id", "deceased", "overall_survival"))
OV_METH_GE_CNV $overall_survival[OV_METH_GE_CNV$overall_survival <= 0] <- 0.001
OV_METH_GE_CNV_survival_data <- Surv(OV_METH_GE_CNV$overall_survival, OV_METH_GE_CNV$deceased)
fusion(OV_METH_GE_CNV, OV_METH_GE_CNV_survival_data, directory = "OV/LF", file_prefix = "OV_METH_GE_CNV")
# METH AND ME AND CNV
OV_METH_ME_CNV <- merge(OV_METH_filtered, OV_ME_normalized, by = c("case_id", "deceased", "overall_survival"))
OV_METH_ME_CNV <- merge(OV_METH_ME_CNV, OV_CNV_normalized, by = c("case_id", "deceased", "overall_survival"))
OV_METH_ME_CNV $overall_survival[OV_METH_ME_CNV$overall_survival <= 0] <- 0.001
OV_METH_ME_CNV_survival_data <- Surv(OV_METH_ME_CNV$overall_survival, OV_METH_ME_CNV$deceased)
fusion(OV_METH_ME_CNV, OV_METH_ME_CNV_survival_data, directory = "OV/LF", file_prefix = "OV_METH_ME_CNV")
# ME AND GE AND CNV
OV_ME_GE_CNV <- merge(OV_ME_normalized, OV_GE_normalized, by = c("case_id", "deceased", "overall_survival"))
OV_ME_GE_CNV <- merge(OV_ME_GE_CNV, OV_CNV_normalized, by = c("case_id", "deceased", "overall_survival"))
OV_ME_GE_CNV $overall_survival[OV_ME_GE_CNV$overall_survival <= 0] <- 0.001
OV_ME_GE_CNV_survival_data <- Surv(OV_ME_GE_CNV$overall_survival, OV_ME_GE_CNV$deceased)
fusion(OV_ME_GE_CNV, OV_ME_GE_CNV_survival_data, directory = "OV/LF", file_prefix = "OV_ME_GE_CNV")
# ME AND GE AND CNV AND METH
OV_ME_GE_CNV_METH <- merge(OV_ME_normalized, OV_GE_normalized, by = c("case_id", "deceased", "overall_survival"))
OV_ME_GE_CNV_METH <- merge(OV_ME_GE_CNV_METH, OV_CNV_normalized, by = c("case_id", "deceased", "overall_survival"))
OV_ME_GE_CNV_METH <- merge(OV_ME_GE_CNV_METH, OV_METH_filtered, by = c("case_id", "deceased", "overall_survival"))
OV_ME_GE_CNV_METH $overall_survival[OV_ME_GE_CNV_METH$overall_survival <= 0] <- 0.001
OV_ME_GE_CNV_METH_survival_data <- Surv(OV_ME_GE_CNV_METH$overall_survival, OV_ME_GE_CNV_METH$deceased)
fusion(OV_ME_GE_CNV_METH, OV_ME_GE_CNV_METH_survival_data, directory = "OV/LF", file_prefix = "OV_ME_GE_CNV_METH")


####################################################################
#                           CESC
#
####################################################################

ME_features=read.csv("CESC/ME/features_cv.csv", row.names = NULL)
GE_features=read.csv("CESC/GE/features_cv.csv", row.names = NULL)
METH_features=read.csv("CESC/METH/features_cv.csv", row.names = NULL)
CNV_features=read.csv("CESC/CNV/features_cv.csv", row.names = NULL)

CESC_ME=read.csv("CESC/CESC_ME_data.csv", row.names = NULL)
CESC_GE=read.csv("CESC/CESC_GE_data.csv", row.names = NULL)
CESC_METH=read.csv("CESC/CESC_METH_data.csv", row.names = NULL)
CESC_CNV=read.csv("CESC/CESC_CNV_data.csv", row.names = NULL)

CESC_METH <- CESC_METH[,-1]
CESC_ME <- CESC_ME[,-1]
CESC_GE <- CESC_GE[,-1]
CESC_CNV <- CESC_CNV[,-1]

# Preprocess-----------------------------------------------------------------------

CESC_ME_filtered <- CESC_ME[, intersect(colnames(CESC_ME), ME_features$Feature)]
CESC_ME_filtered <- cbind(CESC_ME[, 1:3], CESC_ME_filtered)
colnames(CESC_ME_filtered)[-c(1:3)] <- paste("ME_", colnames(CESC_ME_filtered)[-c(1:3)], sep = "")

CESC_GE_filtered <- CESC_GE[, intersect(colnames(CESC_GE), GE_features$Feature)]
CESC_GE_filtered <- cbind(CESC_GE[, 1:3], CESC_GE_filtered)
colnames(CESC_GE_filtered)[-c(1:3)] <- paste("GE_", colnames(CESC_GE_filtered)[-c(1:3)], sep = "")

CESC_METH_filtered <- CESC_METH[, intersect(colnames(CESC_METH), METH_features$Feature)]
CESC_METH_filtered <- cbind(CESC_METH[, 1:3], CESC_METH_filtered)
colnames(CESC_METH_filtered)[-c(1:3)] <- paste("METH_", colnames(CESC_METH_filtered)[-c(1:3)], sep = "")

CESC_CNV_filtered <- CESC_CNV[, intersect(colnames(CESC_CNV), CNV_features$Feature)]
CESC_CNV_filtered <- cbind(CESC_CNV[, 1:3], CESC_CNV_filtered)
colnames(CESC_CNV_filtered)[-c(1:3)] <- paste("CNV_", colnames(CESC_CNV_filtered)[-c(1:3)], sep = "")

# Apply min-max normalization to features (excluding first three columns)
CESC_ME_normalized <- min_max_normalize(CESC_ME_filtered)
CESC_GE_normalized <- min_max_normalize(CESC_GE_filtered)
CESC_CNV_normalized <- min_max_normalize(CESC_CNV_filtered)

#-------------------------------------------------------------------------------
# ME 
CESC_ME_normalized$overall_survival[CESC_ME_normalized$overall_survival <= 0] <- 0.001
CESC_ME_normalized_survival_data <- Surv(CESC_ME_normalized$overall_survival, CESC_ME_normalized$deceased)
fusion(CESC_ME_normalized, CESC_ME_normalized_survival_data, directory = "CESC/LF", file_prefix = "CESC_ME_Fusion")

#-------------------------------------------------------------------------------
# GE 
CESC_GE_normalized$overall_survival[CESC_GE_normalized$overall_survival <= 0] <- 0.001
CESC_GE_normalized_survival_data <- Surv(CESC_GE_normalized$overall_survival, CESC_GE_normalized$deceased)
fusion(CESC_GE_normalized, CESC_GE_normalized_survival_data, directory = "CESC/LF", file_prefix = "CESC_GE_Fusion")

#-------------------------------------------------------------------------------
# METH 
CESC_METH_filtered$overall_survival[CESC_METH_filtered$overall_survival <= 0] <- 0.001
CESC_METH_filtered_survival_data <- Surv(CESC_METH_filtered$overall_survival, CESC_METH_filtered$deceased)
fusion(CESC_METH_filtered, CESC_METH_filtered_survival_data, directory = "CESC/LF", file_prefix = "CESC_METH_Fusion")

#-------------------------------------------------------------------------------
# CNV 
CESC_CNV_normalized$overall_survival[CESC_CNV_normalized$overall_survival <= 0] <- 0.001
CESC_CNV_normalized_survival_data <- Surv(CESC_CNV_normalized$overall_survival, CESC_CNV_normalized$deceased)
fusion(CESC_CNV_normalized, CESC_CNV_normalized_survival_data, directory = "CESC/LF", file_prefix = "CESC_CNV_Fusion")


ME_fusion=read.csv("CESC/LF/CESC_ME_Fusion_data.csv", row.names = NULL)
METH_fusion=read.csv("CESC/LF/CESC_METH_Fusion_data.csv", row.names = NULL)
GE_fusion=read.csv("CESC/LF/CESC_GE_Fusion_data.csv", row.names = NULL)
CNV_fusion=read.csv("CESC/LF/CESC_CNV_Fusion_data.csv", row.names = NULL)

# Preprocess-----------------------------------------------------------------------

CESC_ME_filtered <- CESC_ME_filtered[, names(CESC_ME_filtered) %in% ME_fusion$variables, drop = FALSE]
CESC_ME_filtered <- cbind(CESC_ME[, 1:3], CESC_ME_filtered)
CESC_GE_filtered <- CESC_GE_filtered[, names(CESC_GE_filtered) %in% GE_fusion$variables, drop = FALSE]
CESC_GE_filtered <- cbind(CESC_GE[, 1:3], CESC_GE_filtered)
CESC_METH_filtered <- CESC_METH_filtered[, names(CESC_METH_filtered) %in% METH_fusion$variables, drop = FALSE]
CESC_METH_filtered <- cbind(CESC_METH[, 1:3], CESC_METH_filtered)
CESC_CNV_filtered <- CESC_CNV_filtered[, names(CESC_CNV_filtered) %in% CNV_fusion$variables, drop = FALSE]
CESC_CNV_filtered <- cbind(CESC_CNV[, 1:3], CESC_CNV_filtered)

# Apply min-max normalization to features (excluding first three columns)
CESC_ME_normalized <- min_max_normalize(CESC_ME_filtered)
CESC_GE_normalized <- min_max_normalize(CESC_GE_filtered)
CESC_CNV_normalized <- min_max_normalize(CESC_CNV_filtered)

# Late Fusion -------------------------------------------------------------------
# ME AND GE
CESC_ME_GE <- merge(CESC_ME_normalized, CESC_GE_normalized, by = c("case_id", "deceased", "overall_survival"))
CESC_ME_GE$overall_survival[CESC_ME_GE$overall_survival <= 0] <- 0.001
CESC_ME_GE_survival_data <- Surv(CESC_ME_GE$overall_survival, CESC_ME_GE$deceased)
fusion(CESC_ME_GE, CESC_ME_GE_survival_data, directory = "CESC/LF", file_prefix = "CESC_ME_GE")
# ME AND METH
CESC_ME_METH <- merge(CESC_ME_normalized, CESC_METH_filtered, by = c("case_id", "deceased", "overall_survival"))
CESC_ME_METH$overall_survival[CESC_ME_METH$overall_survival <= 0] <- 0.001
CESC_ME_METH_survival_data <- Surv(CESC_ME_METH$overall_survival, CESC_ME_METH$deceased)
fusion(CESC_ME_METH, CESC_ME_METH_survival_data, directory = "CESC/LF", file_prefix = "CESC_ME_METH")
# METH AND GE
CESC_METH_GE <- merge(CESC_METH_filtered, CESC_GE_normalized, by = c("case_id", "deceased", "overall_survival"))
CESC_METH_GE$overall_survival[CESC_METH_GE$overall_survival <= 0] <- 0.001
CESC_METH_GE_survival_data <- Surv(CESC_METH_GE$overall_survival, CESC_METH_GE$deceased)
fusion(CESC_METH_GE, CESC_METH_GE_survival_data, directory = "CESC/LF", file_prefix = "CESC_METH_GE")
# GE AND CNV
CESC_GE_CNV <- merge(CESC_GE_normalized, CESC_CNV_normalized, by = c("case_id", "deceased", "overall_survival"))
CESC_GE_CNV$overall_survival[CESC_GE_CNV$overall_survival <= 0] <- 0.001
CESC_GE_CNV_survival_data <- Surv(CESC_GE_CNV$overall_survival, CESC_GE_CNV$deceased)
fusion(CESC_GE_CNV, CESC_GE_CNV_survival_data, directory = "CESC/LF", file_prefix = "CESC_GE_CNV")
# ME AND CNV
CESC_ME_CNV <- merge(CESC_ME_normalized, CESC_CNV_normalized, by = c("case_id", "deceased", "overall_survival"))
CESC_ME_CNV$overall_survival[CESC_ME_CNV$overall_survival <= 0] <- 0.001
CESC_ME_CNV_survival_data <- Surv(CESC_ME_CNV$overall_survival, CESC_ME_CNV$deceased)
fusion(CESC_ME_CNV, CESC_ME_CNV_survival_data, directory = "CESC/LF", file_prefix = "CESC_ME_CNV")
# METH AND CNV ------------------------
CESC_METH_CNV <- merge(CESC_METH_filtered, CESC_CNV_normalized, by = c("case_id", "deceased", "overall_survival"))
CESC_METH_CNV$overall_survival[CESC_METH_CNV$overall_survival <= 0] <- 0.001
CESC_METH_CNV_survival_data <- Surv(CESC_METH_CNV$overall_survival, CESC_METH_CNV$deceased)
fusion(CESC_METH_CNV, CESC_METH_CNV_survival_data, directory = "CESC/LF", file_prefix = "CESC_METH_CNV")
# METH AND GE AND ME
CESC_METH_GE_ME <- merge(CESC_METH_filtered, CESC_GE_normalized, by = c("case_id", "deceased", "overall_survival"))
CESC_METH_GE_ME <- merge(CESC_METH_GE_ME, CESC_ME_normalized, by = c("case_id", "deceased", "overall_survival"))
CESC_METH_GE_ME $overall_survival[CESC_METH_GE_ME$overall_survival <= 0] <- 0.001
CESC_METH_GE_ME_survival_data <- Surv(CESC_METH_GE_ME$overall_survival, CESC_METH_GE_ME$deceased)
fusion(CESC_METH_GE_ME, CESC_METH_GE_ME_survival_data, directory = "CESC/LF", file_prefix = "CESC_METH_GE_ME")
# METH AND GE AND CNV
CESC_METH_GE_CNV <- merge(CESC_METH_filtered, CESC_GE_normalized, by = c("case_id", "deceased", "overall_survival"))
CESC_METH_GE_CNV <- merge(CESC_METH_GE_CNV, CESC_CNV_normalized, by = c("case_id", "deceased", "overall_survival"))
CESC_METH_GE_CNV $overall_survival[CESC_METH_GE_CNV$overall_survival <= 0] <- 0.001
CESC_METH_GE_CNV_survival_data <- Surv(CESC_METH_GE_CNV$overall_survival, CESC_METH_GE_CNV$deceased)
fusion(CESC_METH_GE_CNV, CESC_METH_GE_CNV_survival_data, directory = "CESC/LF", file_prefix = "CESC_METH_GE_CNV")
# METH AND ME AND CNV
CESC_METH_ME_CNV <- merge(CESC_METH_filtered, CESC_ME_normalized, by = c("case_id", "deceased", "overall_survival"))
CESC_METH_ME_CNV <- merge(CESC_METH_ME_CNV, CESC_CNV_normalized, by = c("case_id", "deceased", "overall_survival"))
CESC_METH_ME_CNV $overall_survival[CESC_METH_ME_CNV$overall_survival <= 0] <- 0.001
CESC_METH_ME_CNV_survival_data <- Surv(CESC_METH_ME_CNV$overall_survival, CESC_METH_ME_CNV$deceased)
fusion(CESC_METH_ME_CNV, CESC_METH_ME_CNV_survival_data, directory = "CESC/LF", file_prefix = "CESC_METH_ME_CNV")
# ME AND GE AND CNV
CESC_ME_GE_CNV <- merge(CESC_ME_normalized, CESC_GE_normalized, by = c("case_id", "deceased", "overall_survival"))
CESC_ME_GE_CNV <- merge(CESC_ME_GE_CNV, CESC_CNV_normalized, by = c("case_id", "deceased", "overall_survival"))
CESC_ME_GE_CNV $overall_survival[CESC_ME_GE_CNV$overall_survival <= 0] <- 0.001
CESC_ME_GE_CNV_survival_data <- Surv(CESC_ME_GE_CNV$overall_survival, CESC_ME_GE_CNV$deceased)
fusion(CESC_ME_GE_CNV, CESC_ME_GE_CNV_survival_data, directory = "CESC/LF", file_prefix = "CESC_ME_GE_CNV")
# ME AND GE AND CNV AND METH
CESC_ME_GE_CNV_METH <- merge(CESC_ME_normalized, CESC_GE_normalized, by = c("case_id", "deceased", "overall_survival"))
CESC_ME_GE_CNV_METH <- merge(CESC_ME_GE_CNV_METH, CESC_CNV_normalized, by = c("case_id", "deceased", "overall_survival"))
CESC_ME_GE_CNV_METH <- merge(CESC_ME_GE_CNV_METH, CESC_METH_filtered, by = c("case_id", "deceased", "overall_survival"))
CESC_ME_GE_CNV_METH $overall_survival[CESC_ME_GE_CNV_METH$overall_survival <= 0] <- 0.001
CESC_ME_GE_CNV_METH_survival_data <- Surv(CESC_ME_GE_CNV_METH$overall_survival, CESC_ME_GE_CNV_METH$deceased)
fusion(CESC_ME_GE_CNV_METH, CESC_ME_GE_CNV_METH_survival_data, directory = "CESC/LF", file_prefix = "CESC_ME_GE_CNV_METH")


####################################################################
#                           UCEC
#
####################################################################

ME_features=read.csv("UCEC/ME/features_cv.csv", row.names = NULL)
GE_features=read.csv("UCEC/GE/features_cv.csv", row.names = NULL)
METH_features=read.csv("UCEC/METH/features_cv.csv", row.names = NULL)
CNV_features=read.csv("UCEC/CNV/features_cv.csv", row.names = NULL)

UCEC_ME=read.csv("UCEC/UCEC_ME_data.csv", row.names = NULL)
UCEC_GE=read.csv("UCEC/UCEC_GE_data.csv", row.names = NULL)
UCEC_METH=read.csv("UCEC/UCEC_METH_data.csv", row.names = NULL)
UCEC_CNV=read.csv("UCEC/UCEC_CNV_data.csv", row.names = NULL)

UCEC_METH <- UCEC_METH[,-1]
UCEC_ME <- UCEC_ME[,-1]
UCEC_GE <- UCEC_GE[,-1]
UCEC_CNV <- UCEC_CNV[,-1]

# Preprocess-----------------------------------------------------------------------

UCEC_ME_filtered <- UCEC_ME[, intersect(colnames(UCEC_ME), ME_features$Feature)]
UCEC_ME_filtered <- cbind(UCEC_ME[, 1:3], UCEC_ME_filtered)
colnames(UCEC_ME_filtered)[-c(1:3)] <- paste("ME_", colnames(UCEC_ME_filtered)[-c(1:3)], sep = "")

UCEC_GE_filtered <- UCEC_GE[, intersect(colnames(UCEC_GE), GE_features$Feature)]
UCEC_GE_filtered <- cbind(UCEC_GE[, 1:3], UCEC_GE_filtered)
colnames(UCEC_GE_filtered)[-c(1:3)] <- paste("GE_", colnames(UCEC_GE_filtered)[-c(1:3)], sep = "")

UCEC_METH_filtered <- UCEC_METH[, intersect(colnames(UCEC_METH), METH_features$Feature)]
UCEC_METH_filtered <- cbind(UCEC_METH[, 1:3], UCEC_METH_filtered)
colnames(UCEC_METH_filtered)[-c(1:3)] <- paste("METH_", colnames(UCEC_METH_filtered)[-c(1:3)], sep = "")

UCEC_CNV_filtered <- UCEC_CNV[, intersect(colnames(UCEC_CNV), CNV_features$Feature)]
UCEC_CNV_filtered <- cbind(UCEC_CNV[, 1:3], UCEC_CNV_filtered)
colnames(UCEC_CNV_filtered)[-c(1:3)] <- paste("CNV_", colnames(UCEC_CNV_filtered)[-c(1:3)], sep = "")

# Apply min-max normalization to features (excluding first three columns)
UCEC_ME_normalized <- min_max_normalize(UCEC_ME_filtered)
UCEC_GE_normalized <- min_max_normalize(UCEC_GE_filtered)
UCEC_CNV_normalized <- min_max_normalize(UCEC_CNV_filtered)

#-------------------------------------------------------------------------------
# ME 
UCEC_ME_normalized$overall_survival[UCEC_ME_normalized$overall_survival <= 0] <- 0.001
UCEC_ME_normalized_survival_data <- Surv(UCEC_ME_normalized$overall_survival, UCEC_ME_normalized$deceased)
fusion(UCEC_ME_normalized, UCEC_ME_normalized_survival_data, directory = "UCEC/LF", file_prefix = "UCEC_ME_Fusion")

#-------------------------------------------------------------------------------
# GE 
UCEC_GE_normalized$overall_survival[UCEC_GE_normalized$overall_survival <= 0] <- 0.001
UCEC_GE_normalized_survival_data <- Surv(UCEC_GE_normalized$overall_survival, UCEC_GE_normalized$deceased)
fusion(UCEC_GE_normalized, UCEC_GE_normalized_survival_data, directory = "UCEC/LF", file_prefix = "UCEC_GE_Fusion")

#-------------------------------------------------------------------------------
# METH 
UCEC_METH_filtered$overall_survival[UCEC_METH_filtered$overall_survival <= 0] <- 0.001
UCEC_METH_filtered_survival_data <- Surv(UCEC_METH_filtered$overall_survival, UCEC_METH_filtered$deceased)
fusion(UCEC_METH_filtered, UCEC_METH_filtered_survival_data, directory = "UCEC/LF", file_prefix = "UCEC_METH_Fusion")

#-------------------------------------------------------------------------------
# CNV 
UCEC_CNV_normalized$overall_survival[UCEC_CNV_normalized$overall_survival <= 0] <- 0.001
UCEC_CNV_normalized_survival_data <- Surv(UCEC_CNV_normalized$overall_survival, UCEC_CNV_normalized$deceased)
fusion(UCEC_CNV_normalized, UCEC_CNV_normalized_survival_data, directory = "UCEC/LF", file_prefix = "UCEC_CNV_Fusion")


ME_fusion=read.csv("UCEC/LF/UCEC_ME_Fusion_data.csv", row.names = NULL)
METH_fusion=read.csv("UCEC/LF/UCEC_METH_Fusion_data.csv", row.names = NULL)
GE_fusion=read.csv("UCEC/LF/UCEC_GE_Fusion_data.csv", row.names = NULL)
CNV_fusion=read.csv("UCEC/LF/UCEC_CNV_Fusion_data.csv", row.names = NULL)

# Preprocess-----------------------------------------------------------------------

UCEC_ME_filtered <- UCEC_ME_filtered[, names(UCEC_ME_filtered) %in% ME_fusion$variables, drop = FALSE]
UCEC_ME_filtered <- cbind(UCEC_ME[, 1:3], UCEC_ME_filtered)
UCEC_GE_filtered <- UCEC_GE_filtered[, names(UCEC_GE_filtered) %in% GE_fusion$variables, drop = FALSE]
UCEC_GE_filtered <- cbind(UCEC_GE[, 1:3], UCEC_GE_filtered)
UCEC_METH_filtered <- UCEC_METH_filtered[, names(UCEC_METH_filtered) %in% METH_fusion$variables, drop = FALSE]
UCEC_METH_filtered <- cbind(UCEC_METH[, 1:3], UCEC_METH_filtered)
UCEC_CNV_filtered <- UCEC_CNV_filtered[, names(UCEC_CNV_filtered) %in% CNV_fusion$variables, drop = FALSE]
UCEC_CNV_filtered <- cbind(UCEC_CNV[, 1:3], UCEC_CNV_filtered)

# Apply min-max normalization to features (excluding first three columns)
UCEC_ME_normalized <- min_max_normalize(UCEC_ME_filtered)
UCEC_GE_normalized <- min_max_normalize(UCEC_GE_filtered)
UCEC_CNV_normalized <- min_max_normalize(UCEC_CNV_filtered)

# Late Fusion -------------------------------------------------------------------
# ME AND GE
UCEC_ME_GE <- merge(UCEC_ME_normalized, UCEC_GE_normalized, by = c("case_id", "deceased", "overall_survival"))
UCEC_ME_GE$overall_survival[UCEC_ME_GE$overall_survival <= 0] <- 0.001
UCEC_ME_GE_survival_data <- Surv(UCEC_ME_GE$overall_survival, UCEC_ME_GE$deceased)
fusion(UCEC_ME_GE, UCEC_ME_GE_survival_data, directory = "UCEC/LF", file_prefix = "UCEC_ME_GE")
# ME AND METH
UCEC_ME_METH <- merge(UCEC_ME_normalized, UCEC_METH_filtered, by = c("case_id", "deceased", "overall_survival"))
UCEC_ME_METH$overall_survival[UCEC_ME_METH$overall_survival <= 0] <- 0.001
UCEC_ME_METH_survival_data <- Surv(UCEC_ME_METH$overall_survival, UCEC_ME_METH$deceased)
fusion(UCEC_ME_METH, UCEC_ME_METH_survival_data, directory = "UCEC/LF", file_prefix = "UCEC_ME_METH")
# METH AND GE
UCEC_METH_GE <- merge(UCEC_METH_filtered, UCEC_GE_normalized, by = c("case_id", "deceased", "overall_survival"))
UCEC_METH_GE$overall_survival[UCEC_METH_GE$overall_survival <= 0] <- 0.001
UCEC_METH_GE_survival_data <- Surv(UCEC_METH_GE$overall_survival, UCEC_METH_GE$deceased)
fusion(UCEC_METH_GE, UCEC_METH_GE_survival_data, directory = "UCEC/LF", file_prefix = "UCEC_METH_GE")
# GE AND CNV
UCEC_GE_CNV <- merge(UCEC_GE_normalized, UCEC_CNV_normalized, by = c("case_id", "deceased", "overall_survival"))
UCEC_GE_CNV$overall_survival[UCEC_GE_CNV$overall_survival <= 0] <- 0.001
UCEC_GE_CNV_survival_data <- Surv(UCEC_GE_CNV$overall_survival, UCEC_GE_CNV$deceased)
fusion(UCEC_GE_CNV, UCEC_GE_CNV_survival_data, directory = "UCEC/LF", file_prefix = "UCEC_GE_CNV")
# ME AND CNV
UCEC_ME_CNV <- merge(UCEC_ME_normalized, UCEC_CNV_normalized, by = c("case_id", "deceased", "overall_survival"))
UCEC_ME_CNV$overall_survival[UCEC_ME_CNV$overall_survival <= 0] <- 0.001
UCEC_ME_CNV_survival_data <- Surv(UCEC_ME_CNV$overall_survival, UCEC_ME_CNV$deceased)
fusion(UCEC_ME_CNV, UCEC_ME_CNV_survival_data, directory = "UCEC/LF", file_prefix = "UCEC_ME_CNV")
# METH AND CNV ------------------------
UCEC_METH_CNV <- merge(UCEC_METH_filtered, UCEC_CNV_normalized, by = c("case_id", "deceased", "overall_survival"))
UCEC_METH_CNV$overall_survival[UCEC_METH_CNV$overall_survival <= 0] <- 0.001
UCEC_METH_CNV_survival_data <- Surv(UCEC_METH_CNV$overall_survival, UCEC_METH_CNV$deceased)
fusion(UCEC_METH_CNV, UCEC_METH_CNV_survival_data, directory = "UCEC/LF", file_prefix = "UCEC_METH_CNV")
# METH AND GE AND ME
UCEC_METH_GE_ME <- merge(UCEC_METH_filtered, UCEC_GE_normalized, by = c("case_id", "deceased", "overall_survival"))
UCEC_METH_GE_ME <- merge(UCEC_METH_GE_ME, UCEC_ME_normalized, by = c("case_id", "deceased", "overall_survival"))
UCEC_METH_GE_ME $overall_survival[UCEC_METH_GE_ME$overall_survival <= 0] <- 0.001
UCEC_METH_GE_ME_survival_data <- Surv(UCEC_METH_GE_ME$overall_survival, UCEC_METH_GE_ME$deceased)
fusion(UCEC_METH_GE_ME, UCEC_METH_GE_ME_survival_data, directory = "UCEC/LF", file_prefix = "UCEC_METH_GE_ME")
# METH AND GE AND CNV
UCEC_METH_GE_CNV <- merge(UCEC_METH_filtered, UCEC_GE_normalized, by = c("case_id", "deceased", "overall_survival"))
UCEC_METH_GE_CNV <- merge(UCEC_METH_GE_CNV, UCEC_CNV_normalized, by = c("case_id", "deceased", "overall_survival"))
UCEC_METH_GE_CNV $overall_survival[UCEC_METH_GE_CNV$overall_survival <= 0] <- 0.001
UCEC_METH_GE_CNV_survival_data <- Surv(UCEC_METH_GE_CNV$overall_survival, UCEC_METH_GE_CNV$deceased)
fusion(UCEC_METH_GE_CNV, UCEC_METH_GE_CNV_survival_data, directory = "UCEC/LF", file_prefix = "UCEC_METH_GE_CNV")
# METH AND ME AND CNV
UCEC_METH_ME_CNV <- merge(UCEC_METH_filtered, UCEC_ME_normalized, by = c("case_id", "deceased", "overall_survival"))
UCEC_METH_ME_CNV <- merge(UCEC_METH_ME_CNV, UCEC_CNV_normalized, by = c("case_id", "deceased", "overall_survival"))
UCEC_METH_ME_CNV $overall_survival[UCEC_METH_ME_CNV$overall_survival <= 0] <- 0.001
UCEC_METH_ME_CNV_survival_data <- Surv(UCEC_METH_ME_CNV$overall_survival, UCEC_METH_ME_CNV$deceased)
fusion(UCEC_METH_ME_CNV, UCEC_METH_ME_CNV_survival_data, directory = "UCEC/LF", file_prefix = "UCEC_METH_ME_CNV")
# ME AND GE AND CNV
UCEC_ME_GE_CNV <- merge(UCEC_ME_normalized, UCEC_GE_normalized, by = c("case_id", "deceased", "overall_survival"))
UCEC_ME_GE_CNV <- merge(UCEC_ME_GE_CNV, UCEC_CNV_normalized, by = c("case_id", "deceased", "overall_survival"))
UCEC_ME_GE_CNV $overall_survival[UCEC_ME_GE_CNV$overall_survival <= 0] <- 0.001
UCEC_ME_GE_CNV_survival_data <- Surv(UCEC_ME_GE_CNV$overall_survival, UCEC_ME_GE_CNV$deceased)
fusion(UCEC_ME_GE_CNV, UCEC_ME_GE_CNV_survival_data, directory = "UCEC/LF", file_prefix = "UCEC_ME_GE_CNV")
# ME AND GE AND CNV AND METH
UCEC_ME_GE_CNV_METH <- merge(UCEC_ME_normalized, UCEC_GE_normalized, by = c("case_id", "deceased", "overall_survival"))
UCEC_ME_GE_CNV_METH <- merge(UCEC_ME_GE_CNV_METH, UCEC_CNV_normalized, by = c("case_id", "deceased", "overall_survival"))
UCEC_ME_GE_CNV_METH <- merge(UCEC_ME_GE_CNV_METH, UCEC_METH_filtered, by = c("case_id", "deceased", "overall_survival"))
UCEC_ME_GE_CNV_METH $overall_survival[UCEC_ME_GE_CNV_METH$overall_survival <= 0] <- 0.001
UCEC_ME_GE_CNV_METH_survival_data <- Surv(UCEC_ME_GE_CNV_METH$overall_survival, UCEC_ME_GE_CNV_METH$deceased)
fusion(UCEC_ME_GE_CNV_METH, UCEC_ME_GE_CNV_METH_survival_data, directory = "UCEC/LF", file_prefix = "UCEC_ME_GE_CNV_METH")



