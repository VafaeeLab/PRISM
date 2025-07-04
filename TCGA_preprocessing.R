library(dplyr)
library(tidyverse)
library(readr)
library(dplyr)
library(tibble)
library(miRBaseConverter)
library(caret)
.libPaths("/path/to/Rlib")
setwd("/path/to/pwd")
# ---------- Filter Gene Expression Data ----------
BRCA_GE <- read_csv("/path/to/BRCA_primary_tumor_GE_with_clinical.csv") 
CESC_GE <- read_csv("/path/to/CESC_primary_tumor_GE_with_clinical.csv") 
UCEC_GE <- read_csv("/path/to/UCEC_primary_tumor_GE_with_clinical.csv") 
OV_GE <- read_csv("/path/to/OV_primary_tumor_GE_with_clinical.csv") 

filter_omics_GE_matrix <- function(raw_df, clinical_cols = 3, na_threshold = 0.2, var_percentile = 0.90) {
  # 1. Separate clinical and omics data
  clinical_data <- raw_df[, 1:clinical_cols]
  omics_data <- raw_df[, -(1:clinical_cols)]
  
  # 2. Remove features (columns) with > na_threshold missing values
  feature_na_pct <- colMeans(is.na(omics_data))
  omics_data <- omics_data[, feature_na_pct <= na_threshold]
  
  # 3. Remove samples (rows) with > na_threshold missing values
  sample_na_pct <- rowMeans(is.na(omics_data))
  omics_data <- omics_data[sample_na_pct <= na_threshold, ]
  
  # 4. Match clinical data to filtered omics samples
  clinical_data <- clinical_data[rownames(omics_data), ]
  
  # 5. Mean imputation (column-wise)
  for (col in colnames(omics_data)) {
    col_mean <- mean(omics_data[[col]], na.rm = TRUE)
    omics_data[[col]][is.na(omics_data[[col]])] <- col_mean
  }
  
  # 6. Remove low-variance features using percentile threshold
  feature_vars <- apply(omics_data, 2, var)
  var_threshold <- quantile(feature_vars, probs = var_percentile, na.rm = TRUE)
  omics_data <- omics_data[, feature_vars >= var_threshold]
  
  # 7. Combine clinical and omics data
  filtered_df <- cbind(clinical_data, omics_data)
  
  return(filtered_df)
}

BRCA_GE_clean <- filter_omics_GE_matrix(BRCA_GE)
CESC_GE_clean <- filter_omics_GE_matrix(CESC_GE)
UCEC_GE_clean <- filter_omics_GE_matrix(UCEC_GE)
OV_GE_clean <- filter_omics_GE_matrix(OV_GE)

# ---------- Filter miRNA Expression Data ----------
BRCA_ME <- read_csv("/path/to/BRCA_primary_tumor_ME_with_clinical.csv") 
CESC_ME <- read_csv("/path/to/CESC_primary_tumor_ME_with_clinical.csv") 
UCEC_ME <- read_csv("/path/to/UCEC_primary_tumor_ME_with_clinical.csv") 
OV_ME <- read_csv("/path/to/OV_primary_tumor_ME_with_clinical.csv") 


filter_omics_ME_matrix <- function(raw_df, clinical_cols = 3, na_threshold = 0.2, var_threshold = 0, apply_expression_filter = TRUE) {
  # 1. Separate clinical and omics data
  clinical_data <- raw_df[, 1:clinical_cols]
  omics_data <- raw_df[, -(1:clinical_cols)]
  
  # 2. Remove features (columns) with > na_threshold missing values
  feature_na_pct <- colMeans(is.na(omics_data))
  omics_data <- omics_data[, feature_na_pct <= na_threshold]
  
  # 3. Remove samples (rows) with > na_threshold missing values
  sample_na_pct <- rowMeans(is.na(omics_data))
  omics_data <- omics_data[sample_na_pct <= na_threshold, ]
  
  # 4. Match clinical data to filtered omics samples
  clinical_data <- clinical_data[rownames(omics_data), ]
  
  # 5. Mean imputation (column-wise)
  for (col in colnames(omics_data)) {
    col_mean <- mean(omics_data[[col]], na.rm = TRUE)
    omics_data[[col]][is.na(omics_data[[col]])] <- col_mean
  }
  
  # 6. miRNA expression-level filtering
  if (apply_expression_filter) {
    keep_miRNA <- sapply(omics_data, function(x) {
      expr_gt0 <- sum(x > 0, na.rm = TRUE) / length(x)
      expr_gt1 <- sum(x > 1, na.rm = TRUE) / length(x)
      return(expr_gt0 > 0.5 && expr_gt1 > 0.1)
    })
    omics_data <- omics_data[, keep_miRNA]
  }
  
  # 7. Remove low-variance features
  feature_vars <- apply(omics_data, 2, var)
  omics_data <- omics_data[, feature_vars > var_threshold]
  
  # 8. Combine clinical and omics data
  filtered_df <- cbind(clinical_data, omics_data)
  
  return(filtered_df)
}


BRCA_ME_clean <- filter_omics_ME_matrix(BRCA_ME)
CESC_ME_clean <- filter_omics_ME_matrix(CESC_ME)
UCEC_ME_clean <- filter_omics_ME_matrix(UCEC_ME)
OV_ME_clean <- filter_omics_ME_matrix(OV_ME)


# ---------- Filter CNV Expression Data ----------
BRCA_CNV <- read_csv("/path/to/BRCA_primary_tumor_CNV_with_clinical.csv") 
CESC_CNV <- read_csv("/path/to/CESC_primary_tumor_CNV_with_clinical.csv") 
UCEC_CNV <- read_csv("/path/to/UCEC_primary_tumor_CNV_with_clinical.csv") 
OV_CNV <- read_csv("/path/to/OV_primary_tumor_CNV_with_clinical.csv") 

filter_omics_CNV_matrix<- function(raw_df, clinical_cols = 3, na_threshold = 0.2, var_percentile = 0.90) {
  # 1. Separate clinical and CNV data
  clinical_data <- raw_df[, 1:clinical_cols]
  cnv_data <- raw_df[, -(1:clinical_cols)]
  
  # 2. Remove features (columns) with > na_threshold missing values
  feature_na_pct <- colMeans(is.na(cnv_data))
  cnv_data <- cnv_data[, feature_na_pct <= na_threshold]
  
  # 3. Impute using the mode (most common value per column)
  impute_mode <- function(x) {
    ux <- na.omit(unique(x))
    ux[which.max(tabulate(match(x, ux)))]
  }
  
  for (col in colnames(cnv_data)) {
    mode_val <- impute_mode(cnv_data[[col]])
    cnv_data[[col]][is.na(cnv_data[[col]])] <- mode_val
  }
  
  # 4. Remove low-variance features using percentile threshold
  feature_vars <- apply(cnv_data, 2, var)
  var_threshold <- quantile(feature_vars, probs = var_percentile, na.rm = TRUE)
  cnv_data <- cnv_data[, feature_vars >= var_threshold]
  
  # 5. Recombine clinical and CNV data
  filtered_df <- cbind(clinical_data, cnv_data)
  
  return(filtered_df)
}

filter_omics_CNV_BRCA_matrix <- function(raw_df, clinical_cols = 3,
                                         na_threshold = 0.2,
                                         var_percentile = 0.85) {
  # 1. Separate clinical and CNV data
  clinical_data <- raw_df[, 1:clinical_cols]
  cnv_data <- raw_df[, -(1:clinical_cols)]
  
  # 2. Remove features with too many missing values
  feature_na_pct <- colMeans(is.na(cnv_data))
  cnv_data <- cnv_data[, feature_na_pct <= na_threshold, drop = FALSE]
  
  # 3. Mode imputation (works well for discrete CNV values)
  impute_mode <- function(x) {
    ux <- na.omit(unique(x))
    ux[which.max(tabulate(match(x, ux)))]
  }
  
  for (col in colnames(cnv_data)) {
    mode_val <- impute_mode(cnv_data[[col]])
    cnv_data[[col]][is.na(cnv_data[[col]])] <- mode_val
  }
  
  # 4. Remove low-variance features using percentile threshold
  feature_vars <- apply(cnv_data, 2, var)
  var_threshold <- quantile(feature_vars, probs = var_percentile, na.rm = TRUE)
  cnv_data <- cnv_data[, feature_vars >= var_threshold, drop = FALSE]
  
  # 5. Combine with clinical data
  filtered_df <- cbind(clinical_data, cnv_data)
  
  return(filtered_df)
}

BRCA_CNV_clean <- filter_omics_CNV_BRCA_matrix(BRCA_CNV)
CESC_CNV_clean <- filter_omics_CNV_matrix(CESC_CNV)
UCEC_CNV_clean <- filter_omics_CNV_matrix(UCEC_CNV)
OV_CNV_clean <- filter_omics_CNV_matrix(OV_CNV)


# ---------- Filter DM Expression Data ----------
BRCA_DM27 <- read_csv("/path/to/BRCA_primary_tumor_DM_27_with_clinical.csv") 
BRCA_DM450 <- read_csv("/path/to/BRCA_primary_tumor_DM_450_with_clinical.csv") 
CESC_DM <- read_csv("/path/to/CESC_primary_tumor_DM_450_with_clinical.csv") 
UCEC_DM450 <- read_csv("/path/to/UCEC_primary_tumor_DM_450_with_clinical.csv") 
UCEC_DM27 <- read_csv("/path/to/UCEC_primary_tumor_DM_27_with_clinical.csv")
OV_DM <- read_csv("/path/to/OV_primary_tumor_DM_27_with_clinical.csv") 

identical(colnames(UCEC_DM27), colnames(BRCA_DM27))

filter_omics_DM_matrix <- function(dm27_df, dm450_df, clinical_cols = 3, na_threshold = 0.5, var_threshold = 0.001) {
  # 1. Extract omics column names from 27k (excluding clinical columns)
  dm27_features <- colnames(dm27_df)[-(1:clinical_cols)]
  
  # 2. Separate clinical and omics data from 450k
  clinical_data <- dm450_df[, 1:clinical_cols]
  omics_data <- dm450_df[, -(1:clinical_cols)]
  
  # 3. Restrict 450k data to features that overlap with 27k
  common_features <- intersect(colnames(omics_data), dm27_features)
  omics_data <- omics_data[, common_features, drop = FALSE]
  
  # 4. Remove features with > na_threshold missing values
  feature_na_pct <- colMeans(is.na(omics_data))
  omics_data <- omics_data[, feature_na_pct <= na_threshold]
  
  # 5. Remove samples (rows) with > na_threshold missing values
  sample_na_pct <- rowMeans(is.na(omics_data))
  omics_data <- omics_data[sample_na_pct <= na_threshold, , drop = FALSE]
  
  # 6. Subset clinical data to match filtered omics data
  clinical_data <- clinical_data[rownames(omics_data), , drop = FALSE]
  
  # 7. Mean imputation (column-wise)
  for (col in colnames(omics_data)) {
    col_mean <- mean(omics_data[[col]], na.rm = TRUE)
    omics_data[[col]][is.na(omics_data[[col]])] <- col_mean
  }
  
  # 8. Remove low-variance features
  feature_vars <- apply(omics_data, 2, var)
  omics_data <- omics_data[, feature_vars > var_threshold, drop = FALSE]
  
  # 9. Combine clinical and omics data
  filtered_df <- cbind(clinical_data, omics_data)
  
  return(filtered_df)
}

filter_omics_DM_OV_UCEC_matrix <- function(dm27_df, dm450_df, clinical_cols = 3, na_threshold = 0.5, var_threshold = 0.001, var_percentile = 0.90) {
  # 1. Extract omics column names from 27k (excluding clinical columns)
  dm27_features <- colnames(dm27_df)[-(1:clinical_cols)]
  
  # 2. Separate clinical and omics data from 450k
  clinical_data <- dm450_df[, 1:clinical_cols]
  omics_data <- dm450_df[, -(1:clinical_cols)]
  
  # 3. Restrict 450k data to features that overlap with 27k
  common_features <- intersect(colnames(omics_data), dm27_features)
  omics_data <- omics_data[, common_features, drop = FALSE]
  
  # 4. Remove features with > na_threshold missing values
  feature_na_pct <- colMeans(is.na(omics_data))
  omics_data <- omics_data[, feature_na_pct <= na_threshold]
  
  # 5. Remove samples (rows) with > na_threshold missing values
  sample_na_pct <- rowMeans(is.na(omics_data))
  omics_data <- omics_data[sample_na_pct <= na_threshold, , drop = FALSE]
  
  # 6. Subset clinical data to match filtered omics data
  clinical_data <- clinical_data[rownames(omics_data), , drop = FALSE]
  
  # 7. Mean imputation (column-wise)
  for (col in colnames(omics_data)) {
    col_mean <- mean(omics_data[[col]], na.rm = TRUE)
    omics_data[[col]][is.na(omics_data[[col]])] <- col_mean
  }
  
  # 8. Remove low-variance features using percentile threshold
  feature_vars <- apply(omics_data, 2, var)
  var_percentile_cutoff <- quantile(feature_vars, probs = var_percentile, na.rm = TRUE)
  final_threshold <- max(var_threshold, var_percentile_cutoff)
  valid_features <- names(feature_vars[feature_vars >= final_threshold])
  omics_data <- omics_data[, valid_features, drop = FALSE]
  # 9. Combine clinical and omics data
  filtered_df <- cbind(clinical_data, omics_data)
  
  return(filtered_df)
}

BRCA_DM_clean <- filter_omics_DM_matrix(BRCA_DM27, BRCA_DM450)
CESC_DM_clean <- filter_omics_DM_matrix(BRCA_DM27, CESC_DM)
UCEC_DM_clean <- filter_omics_DM_OV_UCEC_matrix(UCEC_DM27, UCEC_DM450)
OV_DM_clean <- filter_omics_DM_OV_UCEC_matrix(BRCA_DM27, OV_DM)

min_max_scale_omics <- function(df, clinical_cols = 3) {
  # Extract clinical data (keep as is)
  clinical_data <- df[, 1:clinical_cols]
  
  # Extract omics features
  omics_data <- df[, (clinical_cols + 1):ncol(df)]
  
  # Min-max scale each feature column: (x - min) / (max - min)
  scaled_omics <- as.data.frame(
    apply(omics_data, 2, function(x) {
      rng <- range(x, na.rm = TRUE)
      if(diff(rng) == 0) {
        # If no variation, return zeros (or original)
        return(rep(0, length(x)))
      } else {
        return((x - rng[1]) / diff(rng))
      }
    })
  )
  
  # Combine clinical and scaled omics data
  scaled_df <- cbind(clinical_data, scaled_omics)
  
  return(scaled_df)
}

# ---------- ME Data Rename MIMAT IDs to miRNA names ----------
convert_mimat_to_mirnames <- function(df, clinical_cols = 3, version = "v22") {
  # Extract MIMAT IDs (miRNA accession numbers)
  mimat_ids <- colnames(df)[(clinical_cols + 1):ncol(df)]
  
  # Convert to miRNA names
  converted <- miRNA_AccessionToName(mimat_ids, targetVersion = version)
  
  # Extract converted names (handle NAs if any)
  mir_names <- converted$TargetName
  mir_names[is.na(mir_names)] <- mimat_ids[is.na(mir_names)]  # fallback to original MIMAT ID
  
  # Replace feature column names
  new_colnames <- colnames(df)
  new_colnames[(clinical_cols + 1):ncol(df)] <- mir_names
  colnames(df) <- new_colnames
  
  return(df)
}

BRCA_ME_clean <- convert_mimat_to_mirnames(BRCA_ME_clean)
UCEC_ME_clean <- convert_mimat_to_mirnames(UCEC_ME_clean)
CESC_ME_clean <- convert_mimat_to_mirnames(CESC_ME_clean)
OV_ME_clean   <- convert_mimat_to_mirnames(OV_ME_clean)

# ---------- Min-max Normalise the data GE and ME ----------

BRCA_GE_scaled <- min_max_scale_omics(BRCA_GE_clean, clinical_cols = 3)
BRCA_ME_scaled <- min_max_scale_omics(BRCA_ME_clean, clinical_cols = 3)
UCEC_GE_scaled <- min_max_scale_omics(UCEC_GE_clean, clinical_cols = 3)
UCEC_ME_scaled <- min_max_scale_omics(UCEC_ME_clean, clinical_cols = 3)
CESC_GE_scaled <- min_max_scale_omics(CESC_GE_clean, clinical_cols = 3)
CESC_ME_scaled <- min_max_scale_omics(CESC_ME_clean, clinical_cols = 3)
OV_GE_scaled <- min_max_scale_omics(OV_GE_clean, clinical_cols = 3)
OV_ME_scaled <- min_max_scale_omics(OV_ME_clean, clinical_cols = 3)


save_filtered_datasets <- function(common_samples,
                                   ge_df, cnv_df, me_df, dm_df,
                                   dir_path, prefix = NULL) {
  # Create directory if it doesn't exist
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  # Use last part of directory name as prefix if not provided
  if (is.null(prefix)) {
    prefix <- basename(normalizePath(dir_path, mustWork = FALSE))
  }
  
  # Filter datasets
  ge_filtered  <- ge_df[ge_df$sampleID %in% common_samples, ]
  cnv_filtered <- cnv_df[cnv_df$sampleID %in% common_samples, ]
  me_filtered  <- me_df[me_df$sampleID %in% common_samples, ]
  dm_filtered  <- dm_df[dm_df$sampleID %in% common_samples, ]
  
  # Save filtered datasets with dynamic file names
  write.csv(ge_filtered,  file.path(dir_path, paste0(prefix, "_GE_clean.csv")),  row.names = FALSE)
  write.csv(cnv_filtered, file.path(dir_path, paste0(prefix, "_CNV_clean.csv")), row.names = FALSE)
  write.csv(me_filtered,  file.path(dir_path, paste0(prefix, "_ME_clean.csv")),  row.names = FALSE)
  write.csv(dm_filtered,  file.path(dir_path, paste0(prefix, "_DM_clean.csv")),  row.names = FALSE)
  
}

# ---------- BRCA Data ----------

common_samples <- Reduce(intersect, list(
  BRCA_GE_scaled$sampleID,
  BRCA_CNV_clean$sampleID,
  BRCA_ME_scaled$sampleID,
  BRCA_DM_clean$sampleID
))

# View the number of common samples
length(common_samples)

save_filtered_datasets(
  common_samples = common_samples,
  ge_df = BRCA_GE_scaled,
  cnv_df = BRCA_CNV_clean,
  me_df = BRCA_ME_scaled,
  dm_df = BRCA_DM_clean,
  dir_path = "BRCA"
)

# ---------- CESC Data ----------

common_samples <- Reduce(intersect, list(
  CESC_GE_scaled$sampleID,
  CESC_CNV_clean$sampleID,
  CESC_ME_scaled$sampleID,
  CESC_DM_clean$sampleID
))

# View the number of common samples
length(common_samples)

save_filtered_datasets(
  common_samples = common_samples,
  ge_df = CESC_GE_scaled,
  cnv_df = CESC_CNV_clean,
  me_df = CESC_ME_scaled,
  dm_df = CESC_DM_clean,
  dir_path = "CESC"
)

# ---------- UCEC Data ----------

common_samples <- Reduce(intersect, list(
  UCEC_GE_scaled$sampleID,
  UCEC_CNV_clean$sampleID,
  UCEC_ME_scaled$sampleID,
  UCEC_DM_clean$sampleID
))

# View the number of common samples
length(common_samples)

save_filtered_datasets(
  common_samples = common_samples,
  ge_df = UCEC_GE_scaled,
  cnv_df = UCEC_CNV_clean,
  me_df = UCEC_ME_scaled,
  dm_df = UCEC_DM_clean,
  dir_path = "UCEC"
)

# ---------- OV Data ----------

common_samples <- Reduce(intersect, list(
  OV_GE_scaled$sampleID,
  OV_CNV_clean$sampleID,
  OV_ME_scaled$sampleID,
  OV_DM_clean$sampleID
))

# View the number of common samples
length(common_samples)

save_filtered_datasets(
  common_samples = common_samples,
  ge_df = OV_GE_scaled,
  cnv_df = OV_CNV_clean,
  me_df = OV_ME_scaled,
  dm_df = OV_DM_clean,
  dir_path = "OV"
)
