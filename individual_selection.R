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

setwd("/srv/scratch/z5309282") 
source("functions.R")

# Set seed for reproducibility
repeats <- 5
folds <- 5
set.seed(123)

####################################################################
#                           BRCA
#
####################################################################
# ME ------------------------------------------------------------
BRCA_ME=read.csv("BRCA/BRCA_ME_data.csv", row.names = NULL)
BRCA_ME <- BRCA_ME[,-1]
BRCA_ME$overall_survival[BRCA_ME$overall_survival <= 0] <- 0.001
BRCA_ME_survival_data <- Surv(BRCA_ME$overall_survival, BRCA_ME$deceased)
perform_and_save_results(BRCA_ME, BRCA_ME_survival_data, repeats, folds, "BRCA_ME")
# GE ------------------------------------------------------------
BRCA_GE=read.csv("BRCA/BRCA_GE_data.csv", row.names = NULL)
BRCA_GE <- BRCA_GE[,-1]
BRCA_GE$overall_survival[BRCA_GE$overall_survival <= 0] <- 0.001
BRCA_GE_survival_data <- Surv(BRCA_GE$overall_survival, BRCA_GE$deceased)
perform_and_save_results(BRCA_GE, BRCA_GE_survival_data, repeats, folds, "BRCA_GE")
# METH ------------------------------------------------------------
BRCA_METH=read.csv("BRCA/BRCA_METH_data.csv", row.names = NULL)
BRCA_METH <- BRCA_METH[,-1]
BRCA_METH$overall_survival[BRCA_METH$overall_survival <= 0] <- 0.001
BRCA_METH_survival_data <- Surv(BRCA_METH$overall_survival, BRCA_METH$deceased)
perform_and_save_results(BRCA_METH, BRCA_METH_survival_data, repeats, folds, "BRCA_METH")
# CNV ------------------------------------------------------------
BRCA_CNV=read.csv("BRCA/BRCA_CNV_data.csv", row.names = NULL)
BRCA_CNV <- BRCA_CNV[,-1]
BRCA_CNV$overall_survival[BRCA_CNV$overall_survival <= 0] <- 0.001
BRCA_CNV_survival_data <- Surv(BRCA_CNV$overall_survival, BRCA_CNV$deceased)
perform_and_save_results(BRCA_CNV, BRCA_CNV_survival_data, repeats, folds, "BRCA_CNV")

####################################################################
#                           OV
#
####################################################################
# ME ------------------------------------------------------------
OV_ME=read.csv("OV/OV_ME_data.csv", row.names = NULL)
OV_ME <- OV_ME[,-1]
OV_ME$overall_survival[OV_ME$overall_survival <= 0] <- 0.001
OV_ME_survival_data <- Surv(OV_ME$overall_survival, OV_ME$deceased)
perform_and_save_results(OV_ME, OV_ME_survival_data, repeats, folds, "OV_ME")
# GE ------------------------------------------------------------
OV_GE=read.csv("OV/OV_GE_data.csv", row.names = NULL)
OV_GE <- OV_GE[,-1]
OV_GE$overall_survival[OV_GE$overall_survival <= 0] <- 0.001
OV_GE_survival_data <- Surv(OV_GE$overall_survival, OV_GE$deceased)
perform_and_save_results(OV_GE, OV_GE_survival_data, repeats, folds, "OV_GE")
# METH ------------------------------------------------------------
OV_METH=read.csv("OV/OV_METH_data.csv", row.names = NULL)
OV_METH <- OV_METH[,-1]
OV_METH$overall_survival[OV_METH$overall_survival <= 0] <- 0.001
OV_METH_survival_data <- Surv(OV_METH$overall_survival, OV_METH$deceased)
perform_and_save_results(OV_METH, OV_METH_survival_data, repeats, folds, "OV_METH")
# CNV ------------------------------------------------------------
OV_CNV=read.csv("OV/OV_CNV_data.csv", row.names = NULL)
OV_CNV <- OV_CNV[,-1]
OV_CNV$overall_survival[OV_CNV$overall_survival <= 0] <- 0.001
OV_CNV_survival_data <- Surv(OV_CNV$overall_survival, OV_CNV$deceased)
perform_and_save_results(OV_CNV, OV_CNV_survival_data, repeats, folds, "OV_CNV")

####################################################################
#                           CESC
#
####################################################################
# ME ------------------------------------------------------------
CESC_ME=read.csv("CESC/CESC_ME_data.csv", row.names = NULL)
CESC_ME <- CESC_ME[,-1]
CESC_ME$overall_survival[CESC_ME$overall_survival <= 0] <- 0.001
CESC_ME_survival_data <- Surv(CESC_ME$overall_survival, CESC_ME$deceased)
perform_and_save_results(CESC_ME, CESC_ME_survival_data, repeats, folds, "CESC_ME")
# GE ------------------------------------------------------------
CESC_GE=read.csv("CESC/CESC_GE_data.csv", row.names = NULL)
CESC_GE <- CESC_GE[,-1]
CESC_GE$overall_survival[CESC_GE$overall_survival <= 0] <- 0.001
CESC_GE_survival_data <- Surv(CESC_GE$overall_survival, CESC_GE$deceased)
perform_and_save_results(CESC_GE, CESC_GE_survival_data, repeats, folds, "CESC_GE")
# METH ------------------------------------------------------------
CESC_METH=read.csv("CESC/CESC_METH_data.csv", row.names = NULL)
CESC_METH <- CESC_METH[,-1]
CESC_METH$overall_survival[CESC_METH$overall_survival <= 0] <- 0.001
CESC_METH_survival_data <- Surv(CESC_METH$overall_survival, CESC_METH$deceased)
perform_and_save_results(CESC_METH, CESC_METH_survival_data, repeats, folds, "CESC_METH")
# CNV ------------------------------------------------------------
CESC_CNV=read.csv("CESC/CESC_CNV_data.csv", row.names = NULL)
CESC_CNV <- CESC_CNV[,-1]
CESC_CNV$overall_survival[CESC_CNV$overall_survival <= 0] <- 0.001
CESC_CNV_survival_data <- Surv(CESC_CNV$overall_survival, CESC_CNV$deceased)
perform_and_save_results(CESC_CNV, CESC_CNV_survival_data, repeats, folds, "CESC_CNV")

####################################################################
#                           UCEC
#
####################################################################
# ME ------------------------------------------------------------
UCEC_ME=read.csv("UCEC/UCEC_ME_data.csv", row.names = NULL)
UCEC_ME <- UCEC_ME[,-1]
UCEC_ME$overall_survival[UCEC_ME$overall_survival <= 0] <- 0.001
UCEC_ME_survival_data <- Surv(UCEC_ME$overall_survival, UCEC_ME$deceased)
perform_and_save_results(UCEC_ME, UCEC_ME_survival_data, repeats, folds, "UCEC_ME")
# GE ------------------------------------------------------------
UCEC_GE=read.csv("UCEC/UCEC_GE_data.csv", row.names = NULL)
UCEC_GE <- UCEC_GE[,-1]
UCEC_GE$overall_survival[UCEC_GE$overall_survival <= 0] <- 0.001
UCEC_GE_survival_data <- Surv(UCEC_GE$overall_survival, UCEC_GE$deceased)
perform_and_save_results(UCEC_GE, UCEC_GE_survival_data, repeats, folds, "UCEC_GE")
# METH ------------------------------------------------------------
UCEC_METH=read.csv("UCEC/UCEC_METH_data.csv", row.names = NULL)
UCEC_METH <- UCEC_METH[,-1]
UCEC_METH$overall_survival[UCEC_METH$overall_survival <= 0] <- 0.001
UCEC_METH_survival_data <- Surv(UCEC_METH$overall_survival, UCEC_METH$deceased)
perform_and_save_results(UCEC_METH, UCEC_METH_survival_data, repeats, folds, "UCEC_METH")
# CNV ------------------------------------------------------------
UCEC_CNV=read.csv("UCEC/UCEC_CNV_data.csv", row.names = NULL)
UCEC_CNV <- UCEC_CNV[,-1]
UCEC_CNV$overall_survival[UCEC_CNV$overall_survival <= 0] <- 0.001
UCEC_CNV_survival_data <- Surv(UCEC_CNV$overall_survival, UCEC_CNV$deceased)
perform_and_save_results(UCEC_CNV, UCEC_CNV_survival_data, repeats, folds, "UCEC_CNV")