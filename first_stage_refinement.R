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
library(parallel)
library(doParallel)

.libPaths("/gpath/to/Rlib")
setwd("/path/to/pwd")
source("functions.R")
set.seed(123)

# *********************************************************************
# First Stage Refinement R Script Documentation
# *********************************************************************

# Overview:
# The First Stage Refinement script processes the omics data by:
# 1. Reading the omics features selected by the Feature Selection pipeline.
# 2. Performing feature-level fusion on all possible combinations of omics data.


# *********************************************************************

####################################################################
#                           BRCA
#
####################################################################

ME_features=read.csv("BRCA/ME/features_cv.csv", row.names = NULL)
colnames(ME_features) <- capitalize_first(colnames(ME_features))
GE_features=read.csv("BRCA/GE/features_bs", row.names = NULL)
colnames(GE_features) <- capitalize_first(colnames(GE_features))
DM_features=read.csv("BRCA/DM/features_bs", row.names = NULL)
colnames(DM_features) <- capitalize_first(colnames(DM_features))
CNV_features=read.csv("BRCA/CNV/features_bs", row.names = NULL)
colnames(CNV_features) <- capitalize_first(colnames(CNV_features))

BRCA_ME=read.csv("BRCA/BRCA_ME_clean.csv", row.names = NULL)
BRCA_GE=read.csv("BRCA/BRCA_GE_clean.csv", row.names = NULL)
BRCA_DM=read.csv("BRCA/BRCA_DM_clean.csv", row.names = NULL)
BRCA_CNV=read.csv("BRCA/BRCA_CNV_clean.csv", row.names = NULL)


# Preprocess-----------------------------------------------------------------------
BRCA_ME_filtered <- BRCA_ME[, intersect(colnames(BRCA_ME), ME_features$Feature)]
BRCA_ME_filtered <- cbind(BRCA_ME[, 1:3], BRCA_ME_filtered)
colnames(BRCA_ME_filtered)[-c(1:3)] <- paste("ME_", colnames(BRCA_ME_filtered)[-c(1:3)], sep = "")
BRCA_GE_filtered <- BRCA_GE[, intersect(colnames(BRCA_GE), GE_features$Feature)]
BRCA_GE_filtered <- cbind(BRCA_GE[, 1:3], BRCA_GE_filtered)
colnames(BRCA_GE_filtered)[-c(1:3)] <- paste("GE_", colnames(BRCA_GE_filtered)[-c(1:3)], sep = "")
BRCA_DM_filtered <- BRCA_DM[, intersect(colnames(BRCA_DM), DM_features$Feature)]
BRCA_DM_filtered <- cbind(BRCA_DM[, 1:3], BRCA_DM_filtered)
colnames(BRCA_DM_filtered)[-c(1:3)] <- paste("DM_", colnames(BRCA_DM_filtered)[-c(1:3)], sep = "")
BRCA_CNV_filtered <- BRCA_CNV[, intersect(colnames(BRCA_CNV), CNV_features$Feature)]
BRCA_CNV_filtered <- cbind(BRCA_CNV[, 1:3], BRCA_CNV_filtered)
colnames(BRCA_CNV_filtered)[-c(1:3)] <- paste("CNV_", colnames(BRCA_CNV_filtered)[-c(1:3)], sep = "")


BRCA_ME_normalized <- BRCA_ME_filtered
BRCA_GE_normalized <- BRCA_GE_filtered
BRCA_CNV_normalized <- BRCA_CNV_filtered

# First Stage Refinement -------------------------------------------------------------------
# ME AND GE
BRCA_ME_GE <- merge(BRCA_ME_normalized, BRCA_GE_normalized, by = c("sampleID", "deceased", "overall_survival"))
BRCA_ME_GE$overall_survival[BRCA_ME_GE$overall_survival <= 0] <- 0.001
BRCA_ME_GE <- BRCA_ME_GE[!is.na(BRCA_ME_GE$overall_survival), ]
robust_feature_elimination(BRCA_ME_GE, directory = "BRCA/1S", file_prefix = "BRCA_ME_GE")

# ME AND DM
BRCA_ME_DM <- merge(BRCA_ME_normalized, BRCA_DM_filtered, by = c("sampleID", "deceased", "overall_survival"))
BRCA_ME_DM$overall_survival[BRCA_ME_DM$overall_survival <= 0] <- 0.001
BRCA_ME_DM <- BRCA_ME_DM[!is.na(BRCA_ME_DM$overall_survival), ]
robust_feature_elimination(BRCA_ME_DM, directory = "BRCA/1S", file_prefix = "BRCA_ME_DM")

# DM AND GE
BRCA_DM_GE <- merge(BRCA_DM_filtered, BRCA_GE_normalized, by = c("sampleID", "deceased", "overall_survival"))
BRCA_DM_GE$overall_survival[BRCA_DM_GE$overall_survival <= 0] <- 0.001
BRCA_DM_GE <- BRCA_DM_GE[!is.na(BRCA_DM_GE$overall_survival), ]
robust_feature_elimination(BRCA_DM_GE, directory = "BRCA/1S", file_prefix = "BRCA_DM_GE")

# GE AND CNV
BRCA_GE_CNV <- merge(BRCA_GE_normalized, BRCA_CNV_normalized, by = c("sampleID", "deceased", "overall_survival"))
BRCA_GE_CNV$overall_survival[BRCA_GE_CNV$overall_survival <= 0] <- 0.001
BRCA_GE_CNV <- BRCA_GE_CNV[!is.na(BRCA_GE_CNV$overall_survival), ]
robust_feature_elimination(BRCA_GE_CNV, directory = "BRCA/1S", file_prefix = "BRCA_GE_CNV")

# ME AND CNV
BRCA_ME_CNV <- merge(BRCA_ME_normalized, BRCA_CNV_normalized, by = c("sampleID", "deceased", "overall_survival"))
BRCA_ME_CNV$overall_survival[BRCA_ME_CNV$overall_survival <= 0] <- 0.001
BRCA_ME_CNV <- BRCA_ME_CNV[!is.na(BRCA_ME_CNV$overall_survival), ]
robust_feature_elimination(BRCA_ME_CNV, directory = "BRCA/1S", file_prefix = "BRCA_ME_CNV")

# DM AND CNV ------------------------
BRCA_DM_CNV <- merge(BRCA_DM_filtered, BRCA_CNV_normalized, by = c("sampleID", "deceased", "overall_survival"))
BRCA_DM_CNV$overall_survival[BRCA_DM_CNV$overall_survival <= 0] <- 0.001
BRCA_DM_CNV <- BRCA_DM_CNV[!is.na(BRCA_DM_CNV$overall_survival), ]
robust_feature_elimination(BRCA_DM_CNV, directory = "BRCA/1S", file_prefix = "BRCA_DM_CNV")

# DM AND GE AND ME
BRCA_DM_GE_ME <- merge(BRCA_DM_filtered, BRCA_GE_normalized, by = c("sampleID", "deceased", "overall_survival"))
BRCA_DM_GE_ME <- merge(BRCA_DM_GE_ME, BRCA_ME_normalized, by = c("sampleID", "deceased", "overall_survival"))
BRCA_DM_GE_ME $overall_survival[BRCA_DM_GE_ME$overall_survival <= 0] <- 0.001
BRCA_DM_GE_ME <- BRCA_DM_GE_ME[!is.na(BRCA_DM_GE_ME$overall_survival), ]
robust_feature_elimination(BRCA_DM_GE_ME, directory = "BRCA/1S", file_prefix = "BRCA_DM_GE_ME")

# DM AND GE AND CNV
BRCA_DM_GE_CNV <- merge(BRCA_DM_filtered, BRCA_GE_normalized, by = c("sampleID", "deceased", "overall_survival"))
BRCA_DM_GE_CNV <- merge(BRCA_DM_GE_CNV, BRCA_CNV_normalized, by = c("sampleID", "deceased", "overall_survival"))
BRCA_DM_GE_CNV $overall_survival[BRCA_DM_GE_CNV$overall_survival <= 0] <- 0.001
BRCA_DM_GE_CNV <- BRCA_DM_GE_CNV[!is.na(BRCA_DM_GE_CNV$overall_survival), ]
robust_feature_elimination(BRCA_DM_GE_CNV, directory = "BRCA/1S", file_prefix = "BRCA_DM_GE_CNV")

# DM AND ME AND CNV
BRCA_DM_ME_CNV <- merge(BRCA_DM_filtered, BRCA_ME_normalized, by = c("sampleID", "deceased", "overall_survival"))
BRCA_DM_ME_CNV <- merge(BRCA_DM_ME_CNV, BRCA_CNV_normalized, by = c("sampleID", "deceased", "overall_survival"))
BRCA_DM_ME_CNV $overall_survival[BRCA_DM_ME_CNV$overall_survival <= 0] <- 0.001
BRCA_DM_ME_CNV <- BRCA_DM_ME_CNV[!is.na(BRCA_DM_ME_CNV$overall_survival), ]
robust_feature_elimination(BRCA_DM_ME_CNV, directory = "BRCA/1S", file_prefix = "BRCA_DM_ME_CNV")

# ME AND GE AND CNV
BRCA_ME_GE_CNV <- merge(BRCA_ME_normalized, BRCA_GE_normalized, by = c("sampleID", "deceased", "overall_survival"))
BRCA_ME_GE_CNV <- merge(BRCA_ME_GE_CNV, BRCA_CNV_normalized, by = c("sampleID", "deceased", "overall_survival"))
BRCA_ME_GE_CNV $overall_survival[BRCA_ME_GE_CNV$overall_survival <= 0] <- 0.001
BRCA_ME_GE_CNV <- BRCA_ME_GE_CNV[!is.na(BRCA_ME_GE_CNV$overall_survival), ]
robust_feature_elimination(BRCA_ME_GE_CNV, directory = "BRCA/1S", file_prefix = "BRCA_ME_GE_CNV")

# ME AND GE AND CNV AND DM
BRCA_ME_GE_CNV_DM <- merge(BRCA_ME_normalized, BRCA_GE_normalized, by = c("sampleID", "deceased", "overall_survival"))
BRCA_ME_GE_CNV_DM <- merge(BRCA_ME_GE_CNV_DM, BRCA_CNV_normalized, by = c("sampleID", "deceased", "overall_survival"))
BRCA_ME_GE_CNV_DM <- merge(BRCA_ME_GE_CNV_DM, BRCA_DM_filtered, by = c("sampleID", "deceased", "overall_survival"))
BRCA_ME_GE_CNV_DM $overall_survival[BRCA_ME_GE_CNV_DM$overall_survival <= 0] <- 0.001
BRCA_ME_GE_CNV_DM <- BRCA_ME_GE_CNV_DM[!is.na(BRCA_ME_GE_CNV_DM$overall_survival), ]
robust_feature_elimination(BRCA_ME_GE_CNV_DM, directory = "BRCA/1S", file_prefix = "BRCA_ME_GE_CNV_DM")


####################################################################
#                           CESC
#
####################################################################

ME_features=read.csv("CESC/ME/features_cv.csv", row.names = NULL)
colnames(ME_features) <- capitalize_first(colnames(ME_features))
GE_features=read.csv("CESC/GE/features_bs", row.names = NULL)
colnames(GE_features) <- capitalize_first(colnames(GE_features))
DM_features=read.csv("CESC/DM/features_cv.csv", row.names = NULL)
colnames(DM_features) <- capitalize_first(colnames(DM_features))
CNV_features=read.csv("CESC/CNV/features_bs", row.names = NULL)
colnames(CNV_features) <- capitalize_first(colnames(CNV_features))

CESC_ME=read.csv("CESC/CESC_ME_clean.csv", row.names = NULL)
CESC_GE=read.csv("CESC/CESC_GE_clean.csv", row.names = NULL)
CESC_DM=read.csv("CESC/CESC_DM_clean.csv", row.names = NULL)
CESC_CNV=read.csv("CESC/CESC_CNV_clean.csv", row.names = NULL)


# Preprocess-----------------------------------------------------------------------
CESC_ME_filtered <- CESC_ME[, intersect(colnames(CESC_ME), ME_features$Feature)]
CESC_ME_filtered <- cbind(CESC_ME[, 1:3], CESC_ME_filtered)
colnames(CESC_ME_filtered)[-c(1:3)] <- paste("ME_", colnames(CESC_ME_filtered)[-c(1:3)], sep = "")
CESC_GE_filtered <- CESC_GE[, intersect(colnames(CESC_GE), GE_features$Feature)]
CESC_GE_filtered <- cbind(CESC_GE[, 1:3], CESC_GE_filtered)
colnames(CESC_GE_filtered)[-c(1:3)] <- paste("GE_", colnames(CESC_GE_filtered)[-c(1:3)], sep = "")
CESC_DM_filtered <- CESC_DM[, intersect(colnames(CESC_DM), DM_features$Feature)]
CESC_DM_filtered <- cbind(CESC_DM[, 1:3], CESC_DM_filtered)
colnames(CESC_DM_filtered)[-c(1:3)] <- paste("DM_", colnames(CESC_DM_filtered)[-c(1:3)], sep = "")
CESC_CNV_filtered <- CESC_CNV[, intersect(colnames(CESC_CNV), CNV_features$Feature)]
CESC_CNV_filtered <- cbind(CESC_CNV[, 1:3], CESC_CNV_filtered)
colnames(CESC_CNV_filtered)[-c(1:3)] <- paste("CNV_", colnames(CESC_CNV_filtered)[-c(1:3)], sep = "")


CESC_ME_normalized <- CESC_ME_filtered
CESC_GE_normalized <- CESC_GE_filtered
CESC_CNV_normalized <- CESC_CNV_filtered

# First Stage Refinement -------------------------------------------------------------------
# ME AND GE
CESC_ME_GE <- merge(CESC_ME_normalized, CESC_GE_normalized, by = c("sampleID", "deceased", "overall_survival"))
CESC_ME_GE$overall_survival[CESC_ME_GE$overall_survival <= 0] <- 0.001
CESC_ME_GE <- CESC_ME_GE[!is.na(CESC_ME_GE$overall_survival), ]
robust_feature_elimination(CESC_ME_GE, directory = "CESC/1S", file_prefix = "CESC_ME_GE")

# ME AND DM
CESC_ME_DM <- merge(CESC_ME_normalized, CESC_DM_filtered, by = c("sampleID", "deceased", "overall_survival"))
CESC_ME_DM$overall_survival[CESC_ME_DM$overall_survival <= 0] <- 0.001
CESC_ME_DM <- CESC_ME_DM[!is.na(CESC_ME_DM$overall_survival), ]
robust_feature_elimination(CESC_ME_DM, directory = "CESC/1S", file_prefix = "CESC_ME_DM")

# DM AND GE
CESC_DM_GE <- merge(CESC_DM_filtered, CESC_GE_normalized, by = c("sampleID", "deceased", "overall_survival"))
CESC_DM_GE$overall_survival[CESC_DM_GE$overall_survival <= 0] <- 0.001
CESC_DM_GE <- CESC_DM_GE[!is.na(CESC_DM_GE$overall_survival), ]
robust_feature_elimination(CESC_DM_GE, directory = "CESC/1S", file_prefix = "CESC_DM_GE")

# GE AND CNV
CESC_GE_CNV <- merge(CESC_GE_normalized, CESC_CNV_normalized, by = c("sampleID", "deceased", "overall_survival"))
CESC_GE_CNV$overall_survival[CESC_GE_CNV$overall_survival <= 0] <- 0.001
CESC_GE_CNV <- CESC_GE_CNV[!is.na(CESC_GE_CNV$overall_survival), ]
robust_feature_elimination(CESC_GE_CNV, directory = "CESC/1S", file_prefix = "CESC_GE_CNV")

# ME AND CNV
CESC_ME_CNV <- merge(CESC_ME_normalized, CESC_CNV_normalized, by = c("sampleID", "deceased", "overall_survival"))
CESC_ME_CNV$overall_survival[CESC_ME_CNV$overall_survival <= 0] <- 0.001
CESC_ME_CNV <- CESC_ME_CNV[!is.na(CESC_ME_CNV$overall_survival), ]
robust_feature_elimination(CESC_ME_CNV, directory = "CESC/1S", file_prefix = "CESC_ME_CNV")

# DM AND CNV ------------------------
CESC_DM_CNV <- merge(CESC_DM_filtered, CESC_CNV_normalized, by = c("sampleID", "deceased", "overall_survival"))
CESC_DM_CNV$overall_survival[CESC_DM_CNV$overall_survival <= 0] <- 0.001
CESC_DM_CNV <- CESC_DM_CNV[!is.na(CESC_DM_CNV$overall_survival), ]
robust_feature_elimination(CESC_DM_CNV, directory = "CESC/1S", file_prefix = "CESC_DM_CNV")

# DM AND GE AND ME
CESC_DM_GE_ME <- merge(CESC_DM_filtered, CESC_GE_normalized, by = c("sampleID", "deceased", "overall_survival"))
CESC_DM_GE_ME <- merge(CESC_DM_GE_ME, CESC_ME_normalized, by = c("sampleID", "deceased", "overall_survival"))
CESC_DM_GE_ME $overall_survival[CESC_DM_GE_ME$overall_survival <= 0] <- 0.001
CESC_DM_GE_ME <- CESC_DM_GE_ME[!is.na(CESC_DM_GE_ME$overall_survival), ]
robust_feature_elimination(CESC_DM_GE_ME, directory = "CESC/1S", file_prefix = "CESC_DM_GE_ME")

# DM AND GE AND CNV
CESC_DM_GE_CNV <- merge(CESC_DM_filtered, CESC_GE_normalized, by = c("sampleID", "deceased", "overall_survival"))
CESC_DM_GE_CNV <- merge(CESC_DM_GE_CNV, CESC_CNV_normalized, by = c("sampleID", "deceased", "overall_survival"))
CESC_DM_GE_CNV $overall_survival[CESC_DM_GE_CNV$overall_survival <= 0] <- 0.001
CESC_DM_GE_CNV <- CESC_DM_GE_CNV[!is.na(CESC_DM_GE_CNV$overall_survival), ]
robust_feature_elimination(CESC_DM_GE_CNV, directory = "CESC/1S", file_prefix = "CESC_DM_GE_CNV")

# DM AND ME AND CNV
CESC_DM_ME_CNV <- merge(CESC_DM_filtered, CESC_ME_normalized, by = c("sampleID", "deceased", "overall_survival"))
CESC_DM_ME_CNV <- merge(CESC_DM_ME_CNV, CESC_CNV_normalized, by = c("sampleID", "deceased", "overall_survival"))
CESC_DM_ME_CNV $overall_survival[CESC_DM_ME_CNV$overall_survival <= 0] <- 0.001
CESC_DM_ME_CNV <- CESC_DM_ME_CNV[!is.na(CESC_DM_ME_CNV$overall_survival), ]
robust_feature_elimination(CESC_DM_ME_CNV, directory = "CESC/1S", file_prefix = "CESC_DM_ME_CNV")

# ME AND GE AND CNV
CESC_ME_GE_CNV <- merge(CESC_ME_normalized, CESC_GE_normalized, by = c("sampleID", "deceased", "overall_survival"))
CESC_ME_GE_CNV <- merge(CESC_ME_GE_CNV, CESC_CNV_normalized, by = c("sampleID", "deceased", "overall_survival"))
CESC_ME_GE_CNV $overall_survival[CESC_ME_GE_CNV$overall_survival <= 0] <- 0.001
CESC_ME_GE_CNV <- CESC_ME_GE_CNV[!is.na(CESC_ME_GE_CNV$overall_survival), ]
robust_feature_elimination(CESC_ME_GE_CNV, directory = "CESC/1S", file_prefix = "CESC_ME_GE_CNV")

# ME AND GE AND CNV AND DM
CESC_ME_GE_CNV_DM <- merge(CESC_ME_normalized, CESC_GE_normalized, by = c("sampleID", "deceased", "overall_survival"))
CESC_ME_GE_CNV_DM <- merge(CESC_ME_GE_CNV_DM, CESC_CNV_normalized, by = c("sampleID", "deceased", "overall_survival"))
CESC_ME_GE_CNV_DM <- merge(CESC_ME_GE_CNV_DM, CESC_DM_filtered, by = c("sampleID", "deceased", "overall_survival"))
CESC_ME_GE_CNV_DM $overall_survival[CESC_ME_GE_CNV_DM$overall_survival <= 0] <- 0.001
CESC_ME_GE_CNV_DM <- CESC_ME_GE_CNV_DM[!is.na(CESC_ME_GE_CNV_DM$overall_survival), ]
robust_feature_elimination(CESC_ME_GE_CNV_DM, directory = "CESC/1S", file_prefix = "CESC_ME_GE_CNV_DM")


####################################################################
#                           OV
#
####################################################################

ME_features=read.csv("OV/ME/features_cv.csv", row.names = NULL)
colnames(ME_features) <- capitalize_first(colnames(ME_features))
GE_features=read.csv("OV/GE/features_bs", row.names = NULL)
colnames(GE_features) <- capitalize_first(colnames(GE_features))
DM_features=read.csv("OV/DM/features_cv.csv", row.names = NULL)
colnames(DM_features) <- capitalize_first(colnames(DM_features))
CNV_features=read.csv("OV/CNV/features_bs", row.names = NULL)
colnames(CNV_features) <- capitalize_first(colnames(CNV_features))

OV_ME=read.csv("OV/OV_ME_clean.csv", row.names = NULL)
OV_GE=read.csv("OV/OV_GE_clean.csv", row.names = NULL)
OV_DM=read.csv("OV/OV_DM_clean.csv", row.names = NULL)
OV_CNV=read.csv("OV/OV_CNV_clean.csv", row.names = NULL)


# Preprocess-----------------------------------------------------------------------
OV_ME_filtered <- OV_ME[, intersect(colnames(OV_ME), ME_features$Feature)]
OV_ME_filtered <- cbind(OV_ME[, 1:3], OV_ME_filtered)
colnames(OV_ME_filtered)[-c(1:3)] <- paste("ME_", colnames(OV_ME_filtered)[-c(1:3)], sep = "")
OV_GE_filtered <- OV_GE[, intersect(colnames(OV_GE), GE_features$Feature)]
OV_GE_filtered <- cbind(OV_GE[, 1:3], OV_GE_filtered)
colnames(OV_GE_filtered)[-c(1:3)] <- paste("GE_", colnames(OV_GE_filtered)[-c(1:3)], sep = "")
OV_DM_filtered <- OV_DM[, intersect(colnames(OV_DM), DM_features$Feature)]
OV_DM_filtered <- cbind(OV_DM[, 1:3], OV_DM_filtered)
colnames(OV_DM_filtered)[-c(1:3)] <- paste("DM_", colnames(OV_DM_filtered)[-c(1:3)], sep = "")
OV_CNV_filtered <- OV_CNV[, intersect(colnames(OV_CNV), CNV_features$Feature)]
OV_CNV_filtered <- cbind(OV_CNV[, 1:3], OV_CNV_filtered)
colnames(OV_CNV_filtered)[-c(1:3)] <- paste("CNV_", colnames(OV_CNV_filtered)[-c(1:3)], sep = "")


OV_ME_normalized <- OV_ME_filtered
OV_GE_normalized <- OV_GE_filtered
OV_CNV_normalized <- OV_CNV_filtered

# First Stage Refinement -------------------------------------------------------------------
# ME AND GE
OV_ME_GE <- merge(OV_ME_normalized, OV_GE_normalized, by = c("sampleID", "deceased", "overall_survival"))
OV_ME_GE$overall_survival[OV_ME_GE$overall_survival <= 0] <- 0.001
OV_ME_GE <- OV_ME_GE[!is.na(OV_ME_GE$overall_survival), ]
robust_feature_elimination(OV_ME_GE, directory = "OV/1S", file_prefix = "OV_ME_GE")

# ME AND DM
OV_ME_DM <- merge(OV_ME_normalized, OV_DM_filtered, by = c("sampleID", "deceased", "overall_survival"))
OV_ME_DM$overall_survival[OV_ME_DM$overall_survival <= 0] <- 0.001
OV_ME_DM <- OV_ME_DM[!is.na(OV_ME_DM$overall_survival), ]
robust_feature_elimination(OV_ME_DM, directory = "OV/1S", file_prefix = "OV_ME_DM")

# DM AND GE
OV_DM_GE <- merge(OV_DM_filtered, OV_GE_normalized, by = c("sampleID", "deceased", "overall_survival"))
OV_DM_GE$overall_survival[OV_DM_GE$overall_survival <= 0] <- 0.001
OV_DM_GE <- OV_DM_GE[!is.na(OV_DM_GE$overall_survival), ]
robust_feature_elimination(OV_DM_GE, directory = "OV/1S", file_prefix = "OV_DM_GE")

# GE AND CNV
OV_GE_CNV <- merge(OV_GE_normalized, OV_CNV_normalized, by = c("sampleID", "deceased", "overall_survival"))
OV_GE_CNV$overall_survival[OV_GE_CNV$overall_survival <= 0] <- 0.001
OV_GE_CNV <- OV_GE_CNV[!is.na(OV_GE_CNV$overall_survival), ]
robust_feature_elimination(OV_GE_CNV, directory = "OV/1S", file_prefix = "OV_GE_CNV")

# ME AND CNV
OV_ME_CNV <- merge(OV_ME_normalized, OV_CNV_normalized, by = c("sampleID", "deceased", "overall_survival"))
OV_ME_CNV$overall_survival[OV_ME_CNV$overall_survival <= 0] <- 0.001
OV_ME_CNV <- OV_ME_CNV[!is.na(OV_ME_CNV$overall_survival), ]
robust_feature_elimination(OV_ME_CNV, directory = "OV/1S", file_prefix = "OV_ME_CNV")

# DM AND CNV ------------------------
OV_DM_CNV <- merge(OV_DM_filtered, OV_CNV_normalized, by = c("sampleID", "deceased", "overall_survival"))
OV_DM_CNV$overall_survival[OV_DM_CNV$overall_survival <= 0] <- 0.001
OV_DM_CNV <- OV_DM_CNV[!is.na(OV_DM_CNV$overall_survival), ]
robust_feature_elimination(OV_DM_CNV, directory = "OV/1S", file_prefix = "OV_DM_CNV")

# DM AND GE AND ME
OV_DM_GE_ME <- merge(OV_DM_filtered, OV_GE_normalized, by = c("sampleID", "deceased", "overall_survival"))
OV_DM_GE_ME <- merge(OV_DM_GE_ME, OV_ME_normalized, by = c("sampleID", "deceased", "overall_survival"))
OV_DM_GE_ME $overall_survival[OV_DM_GE_ME$overall_survival <= 0] <- 0.001
OV_DM_GE_ME <- OV_DM_GE_ME[!is.na(OV_DM_GE_ME$overall_survival), ]
robust_feature_elimination(OV_DM_GE_ME, directory = "OV/1S", file_prefix = "OV_DM_GE_ME")

# DM AND GE AND CNV
OV_DM_GE_CNV <- merge(OV_DM_filtered, OV_GE_normalized, by = c("sampleID", "deceased", "overall_survival"))
OV_DM_GE_CNV <- merge(OV_DM_GE_CNV, OV_CNV_normalized, by = c("sampleID", "deceased", "overall_survival"))
OV_DM_GE_CNV $overall_survival[OV_DM_GE_CNV$overall_survival <= 0] <- 0.001
OV_DM_GE_CNV <- OV_DM_GE_CNV[!is.na(OV_DM_GE_CNV$overall_survival), ]
robust_feature_elimination(OV_DM_GE_CNV, directory = "OV/1S", file_prefix = "OV_DM_GE_CNV")

# DM AND ME AND CNV
OV_DM_ME_CNV <- merge(OV_DM_filtered, OV_ME_normalized, by = c("sampleID", "deceased", "overall_survival"))
OV_DM_ME_CNV <- merge(OV_DM_ME_CNV, OV_CNV_normalized, by = c("sampleID", "deceased", "overall_survival"))
OV_DM_ME_CNV $overall_survival[OV_DM_ME_CNV$overall_survival <= 0] <- 0.001
OV_DM_ME_CNV <- OV_DM_ME_CNV[!is.na(OV_DM_ME_CNV$overall_survival), ]
robust_feature_elimination(OV_DM_ME_CNV, directory = "OV/1S", file_prefix = "OV_DM_ME_CNV")

# ME AND GE AND CNV
OV_ME_GE_CNV <- merge(OV_ME_normalized, OV_GE_normalized, by = c("sampleID", "deceased", "overall_survival"))
OV_ME_GE_CNV <- merge(OV_ME_GE_CNV, OV_CNV_normalized, by = c("sampleID", "deceased", "overall_survival"))
OV_ME_GE_CNV $overall_survival[OV_ME_GE_CNV$overall_survival <= 0] <- 0.001
OV_ME_GE_CNV <- OV_ME_GE_CNV[!is.na(OV_ME_GE_CNV$overall_survival), ]
robust_feature_elimination(OV_ME_GE_CNV, directory = "OV/1S", file_prefix = "OV_ME_GE_CNV")

# ME AND GE AND CNV AND DM
OV_ME_GE_CNV_DM <- merge(OV_ME_normalized, OV_GE_normalized, by = c("sampleID", "deceased", "overall_survival"))
OV_ME_GE_CNV_DM <- merge(OV_ME_GE_CNV_DM, OV_CNV_normalized, by = c("sampleID", "deceased", "overall_survival"))
OV_ME_GE_CNV_DM <- merge(OV_ME_GE_CNV_DM, OV_DM_filtered, by = c("sampleID", "deceased", "overall_survival"))
OV_ME_GE_CNV_DM $overall_survival[OV_ME_GE_CNV_DM$overall_survival <= 0] <- 0.001
OV_ME_GE_CNV_DM <- OV_ME_GE_CNV_DM[!is.na(OV_ME_GE_CNV_DM$overall_survival), ]
robust_feature_elimination(OV_ME_GE_CNV_DM, directory = "OV/1S", file_prefix = "OV_ME_GE_CNV_DM")

####################################################################
#                           UCEC
#
####################################################################

ME_features=read.csv("UCEC/ME/features_cv.csv", row.names = NULL)
colnames(ME_features) <- capitalize_first(colnames(ME_features))
GE_features=read.csv("UCEC/GE/features_bs", row.names = NULL)
colnames(GE_features) <- capitalize_first(colnames(GE_features))
DM_features=read.csv("UCEC/DM/features_cv.csv", row.names = NULL)
colnames(DM_features) <- capitalize_first(colnames(DM_features))
CNV_features=read.csv("UCEC/CNV/features_bs", row.names = NULL)
colnames(CNV_features) <- capitalize_first(colnames(CNV_features))

UCEC_ME=read.csv("UCEC/UCEC_ME_clean.csv", row.names = NULL)
UCEC_GE=read.csv("UCEC/UCEC_GE_clean.csv", row.names = NULL)
UCEC_DM=read.csv("UCEC/UCEC_DM_clean.csv", row.names = NULL)
UCEC_CNV=read.csv("UCEC/UCEC_CNV_clean.csv", row.names = NULL)


# Preprocess-----------------------------------------------------------------------
UCEC_ME_filtered <- UCEC_ME[, intersect(colnames(UCEC_ME), ME_features$Feature)]
UCEC_ME_filtered <- cbind(UCEC_ME[, 1:3], UCEC_ME_filtered)
colnames(UCEC_ME_filtered)[-c(1:3)] <- paste("ME_", colnames(UCEC_ME_filtered)[-c(1:3)], sep = "")
UCEC_GE_filtered <- UCEC_GE[, intersect(colnames(UCEC_GE), GE_features$Feature)]
UCEC_GE_filtered <- cbind(UCEC_GE[, 1:3], UCEC_GE_filtered)
colnames(UCEC_GE_filtered)[-c(1:3)] <- paste("GE_", colnames(UCEC_GE_filtered)[-c(1:3)], sep = "")
UCEC_DM_filtered <- UCEC_DM[, intersect(colnames(UCEC_DM), DM_features$Feature)]
UCEC_DM_filtered <- cbind(UCEC_DM[, 1:3], UCEC_DM_filtered)
colnames(UCEC_DM_filtered)[-c(1:3)] <- paste("DM_", colnames(UCEC_DM_filtered)[-c(1:3)], sep = "")
UCEC_CNV_filtered <- UCEC_CNV[, intersect(colnames(UCEC_CNV), CNV_features$Feature)]
UCEC_CNV_filtered <- cbind(UCEC_CNV[, 1:3], UCEC_CNV_filtered)
colnames(UCEC_CNV_filtered)[-c(1:3)] <- paste("CNV_", colnames(UCEC_CNV_filtered)[-c(1:3)], sep = "")


UCEC_ME_normalized <- UCEC_ME_filtered
UCEC_GE_normalized <- UCEC_GE_filtered
UCEC_CNV_normalized <- UCEC_CNV_filtered

# First Stage Refinement -------------------------------------------------------------------
# ME AND GE
UCEC_ME_GE <- merge(UCEC_ME_normalized, UCEC_GE_normalized, by = c("sampleID", "deceased", "overall_survival"))
UCEC_ME_GE$overall_survival[UCEC_ME_GE$overall_survival <= 0] <- 0.001
UCEC_ME_GE <- UCEC_ME_GE[!is.na(UCEC_ME_GE$overall_survival), ]
robust_feature_elimination(UCEC_ME_GE, directory = "UCEC/1S", file_prefix = "UCEC_ME_GE")

# ME AND DM
UCEC_ME_DM <- merge(UCEC_ME_normalized, UCEC_DM_filtered, by = c("sampleID", "deceased", "overall_survival"))
UCEC_ME_DM$overall_survival[UCEC_ME_DM$overall_survival <= 0] <- 0.001
UCEC_ME_DM <- UCEC_ME_DM[!is.na(UCEC_ME_DM$overall_survival), ]
robust_feature_elimination(UCEC_ME_DM, directory = "UCEC/1S", file_prefix = "UCEC_ME_DM")

# DM AND GE
UCEC_DM_GE <- merge(UCEC_DM_filtered, UCEC_GE_normalized, by = c("sampleID", "deceased", "overall_survival"))
UCEC_DM_GE$overall_survival[UCEC_DM_GE$overall_survival <= 0] <- 0.001
UCEC_DM_GE <- UCEC_DM_GE[!is.na(UCEC_DM_GE$overall_survival), ]
robust_feature_elimination(UCEC_DM_GE, directory = "UCEC/1S", file_prefix = "UCEC_DM_GE")

# GE AND CNV
UCEC_GE_CNV <- merge(UCEC_GE_normalized, UCEC_CNV_normalized, by = c("sampleID", "deceased", "overall_survival"))
UCEC_GE_CNV$overall_survival[UCEC_GE_CNV$overall_survival <= 0] <- 0.001
UCEC_GE_CNV <- UCEC_GE_CNV[!is.na(UCEC_GE_CNV$overall_survival), ]
robust_feature_elimination(UCEC_GE_CNV, directory = "UCEC/1S", file_prefix = "UCEC_GE_CNV")

# ME AND CNV
UCEC_ME_CNV <- merge(UCEC_ME_normalized, UCEC_CNV_normalized, by = c("sampleID", "deceased", "overall_survival"))
UCEC_ME_CNV$overall_survival[UCEC_ME_CNV$overall_survival <= 0] <- 0.001
UCEC_ME_CNV <- UCEC_ME_CNV[!is.na(UCEC_ME_CNV$overall_survival), ]
robust_feature_elimination(UCEC_ME_CNV, directory = "UCEC/1S", file_prefix = "UCEC_ME_CNV")

# DM AND CNV ------------------------
UCEC_DM_CNV <- merge(UCEC_DM_filtered, UCEC_CNV_normalized, by = c("sampleID", "deceased", "overall_survival"))
UCEC_DM_CNV$overall_survival[UCEC_DM_CNV$overall_survival <= 0] <- 0.001
UCEC_DM_CNV <- UCEC_DM_CNV[!is.na(UCEC_DM_CNV$overall_survival), ]
robust_feature_elimination(UCEC_DM_CNV, directory = "UCEC/1S", file_prefix = "UCEC_DM_CNV")

# DM AND GE AND ME
UCEC_DM_GE_ME <- merge(UCEC_DM_filtered, UCEC_GE_normalized, by = c("sampleID", "deceased", "overall_survival"))
UCEC_DM_GE_ME <- merge(UCEC_DM_GE_ME, UCEC_ME_normalized, by = c("sampleID", "deceased", "overall_survival"))
UCEC_DM_GE_ME $overall_survival[UCEC_DM_GE_ME$overall_survival <= 0] <- 0.001
UCEC_DM_GE_ME <- UCEC_DM_GE_ME[!is.na(UCEC_DM_GE_ME$overall_survival), ]
robust_feature_elimination(UCEC_DM_GE_ME, directory = "UCEC/1S", file_prefix = "UCEC_DM_GE_ME")

# DM AND GE AND CNV
UCEC_DM_GE_CNV <- merge(UCEC_DM_filtered, UCEC_GE_normalized, by = c("sampleID", "deceased", "overall_survival"))
UCEC_DM_GE_CNV <- merge(UCEC_DM_GE_CNV, UCEC_CNV_normalized, by = c("sampleID", "deceased", "overall_survival"))
UCEC_DM_GE_CNV $overall_survival[UCEC_DM_GE_CNV$overall_survival <= 0] <- 0.001
UCEC_DM_GE_CNV <- UCEC_DM_GE_CNV[!is.na(UCEC_DM_GE_CNV$overall_survival), ]
robust_feature_elimination(UCEC_DM_GE_CNV, directory = "UCEC/1S", file_prefix = "UCEC_DM_GE_CNV")

# DM AND ME AND CNV
UCEC_DM_ME_CNV <- merge(UCEC_DM_filtered, UCEC_ME_normalized, by = c("sampleID", "deceased", "overall_survival"))
UCEC_DM_ME_CNV <- merge(UCEC_DM_ME_CNV, UCEC_CNV_normalized, by = c("sampleID", "deceased", "overall_survival"))
UCEC_DM_ME_CNV $overall_survival[UCEC_DM_ME_CNV$overall_survival <= 0] <- 0.001
UCEC_DM_ME_CNV <- UCEC_DM_ME_CNV[!is.na(UCEC_DM_ME_CNV$overall_survival), ]
robust_feature_elimination(UCEC_DM_ME_CNV, directory = "UCEC/1S", file_prefix = "UCEC_DM_ME_CNV")

# ME AND GE AND CNV
UCEC_ME_GE_CNV <- merge(UCEC_ME_normalized, UCEC_GE_normalized, by = c("sampleID", "deceased", "overall_survival"))
UCEC_ME_GE_CNV <- merge(UCEC_ME_GE_CNV, UCEC_CNV_normalized, by = c("sampleID", "deceased", "overall_survival"))
UCEC_ME_GE_CNV $overall_survival[UCEC_ME_GE_CNV$overall_survival <= 0] <- 0.001
UCEC_ME_GE_CNV <- UCEC_ME_GE_CNV[!is.na(UCEC_ME_GE_CNV$overall_survival), ]
robust_feature_elimination(UCEC_ME_GE_CNV, directory = "UCEC/1S", file_prefix = "UCEC_ME_GE_CNV")

# ME AND GE AND CNV AND DM
UCEC_ME_GE_CNV_DM <- merge(UCEC_ME_normalized, UCEC_GE_normalized, by = c("sampleID", "deceased", "overall_survival"))
UCEC_ME_GE_CNV_DM <- merge(UCEC_ME_GE_CNV_DM, UCEC_CNV_normalized, by = c("sampleID", "deceased", "overall_survival"))
UCEC_ME_GE_CNV_DM <- merge(UCEC_ME_GE_CNV_DM, UCEC_DM_filtered, by = c("sampleID", "deceased", "overall_survival"))
UCEC_ME_GE_CNV_DM $overall_survival[UCEC_ME_GE_CNV_DM$overall_survival <= 0] <- 0.001
UCEC_ME_GE_CNV_DM <- UCEC_ME_GE_CNV_DM[!is.na(UCEC_ME_GE_CNV_DM$overall_survival), ]
robust_feature_elimination(UCEC_ME_GE_CNV_DM, directory = "UCEC/1S", file_prefix = "UCEC_ME_GE_CNV_DM")




