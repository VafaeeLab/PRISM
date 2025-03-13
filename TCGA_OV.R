library(TCGAbiolinks)
library(dplyr)
library(DT)
library(sesameData)
library(tidyverse)
library(SummarizedExperiment)
library(dbplyr)
library(sesame)
library(ggVennDiagram)
library(parallel)
library(doParallel)


source("functions.R")
missing_threshold <- 0.2

# ============================================================
# Pre-process Clinical and Omics Data (OV)
# ============================================================

# Description:
# This section prepares the clinical and omics data for TCGA OV. It processes clinical data, including the creation 
# of a "deceased" variable and an "overall survival" variable. Additionally, it builds queries to retrieve gene expression, 
# copy number variation, DNA methylation, and miRNA expression data, and identifies common samples across all datasets.
# The common samples are used for further downstream analysis.

# Key Steps:
# 1. Clinical Data Preprocessing:
#    - Creates a "deceased" variable based on the 'vital_status'.
#    - Creates an "overall survival" variable, using either 'days_to_death' for deceased patients 
#      or 'days_to_last_follow_up' for alive patients.
#    - Selects relevant clinical data for merging with omics data.
# 2. Build Queries for Omics Data Retrieval:
#    - Retrieves gene expression data using RNA-Seq (STAR - Counts) from the GDC.
#    - Retrieves copy number variation (CNV) data at the gene level.
#    - Retrieves DNA methylation data using the Illumina Human Methylation 450 platform.
#    - Retrieves miRNA expression data using miRNA-Seq (BCGSC miRNA Profiling).
# 3. Identify Common Samples Across All Omics Data:
#    - Finds the common samples across all omics datasets to be used in downstream analysis.

# ============================================================

# Clinical Data Preprocessing
clin_OV <- GDCquery_clinic("TCGA-OV", "Clinical")
clin_OV$deceased <- ifelse(clin_OV$vital_status == "Alive", 0, 1)

# Create overall survival variable
clin_OV$overall_survival <- ifelse(clin_OV$vital_status == "Alive",
                                     clin_OV$days_to_last_follow_up,
                                     clin_OV$days_to_death)

# Select relevant columns for clinical data
clin_OV_data <- clin_OV[c('submitter_id', 'deceased', 'overall_survival')]
dim(clin_OV_data)  # Check dimensions of the clinical data

# Build a query to retrieve gene expression data (RNA-Seq)
query_OV_RNA_TP <- GDCquery(project = 'TCGA-OV',
                              data.category = 'Transcriptome Profiling',
                              experimental.strategy = 'RNA-Seq',
                              workflow.type = 'STAR - Counts',
                              sample.type = 'Primary Tumor',
                              access  = 'open')

# Get results for gene expression data
Output_query_OV_RNA_TP <- getResults(query_OV_RNA_TP)
OV_GE_sample <- Output_query_OV_RNA_TP[c('cases.submitter_id')]

# Build a query to retrieve Copy Number Variation (CNV) data
query_OV_CNV <- GDCquery(project = 'TCGA-OV',
                           data.category = 'Copy Number Variation',
                           sample.type = 'Primary Tumor',
                           data.type = 'Gene Level Copy Number',
                           access = 'open')

# Get results for CNV data
Output_query_OV_CNV <- getResults(query_OV_CNV)
OV_CNV_sample <- Output_query_OV_CNV[c('cases.submitter_id')]

# Build a query to retrieve DNA Methylation data
query_OV_Meth <- GDCquery(project = 'TCGA-OV',
                            data.category = 'DNA Methylation',
                            platform = 'Illumina Human Methylation 450',
                            sample.type = 'Primary Tumor',
                            data.type = 'Methylation Beta Value',
                            access = 'open')

# Get results for DNA methylation data
Output_query_OV_Meth <- getResults(query_OV_Meth)
OV_Meth_sample <- Output_query_OV_Meth[c('cases.submitter_id')]

# Build a query to retrieve miRNA expression data
query_OV_ME <- GDCquery(project = 'TCGA-OV',
                          data.category = 'Transcriptome Profiling',
                          experimental.strategy = 'miRNA-Seq',
                          workflow.type = 'BCGSC miRNA Profiling',
                          data.type = 'miRNA Expression Quantification',
                          sample.type = 'Primary Tumor',
                          access = 'open')

# Get results for miRNA expression data
Output_query_OV_ME <- getResults(query_OV_ME)
OV_ME_sample <- Output_query_OV_ME[c('cases.submitter_id')]

# Identify common samples across all omics datasets
common_samples <- Reduce(intersect, list(OV_GE_sample[[1]], 
                                         OV_CNV_sample[[1]], 
                                         OV_Meth_sample[[1]], 
                                         OV_ME_sample[[1]]))

# ============================================================
# Pre-process Gene Expression Data (OV)
# ============================================================

# Description:
# This code handles the pre-processing of TCGA OV gene expression data (FPKM values), 
# filtering out features (genes) and samples with a high percentage of missing data. 
# It also applies imputation for missing values and focuses on high-variance genes that exhibit more dynamic expression. 
# The data is transformed using log-transformation and reshaped into a wide format, which is then saved to a CSV file.

# Key Steps:
# 1. Data Download & Preparation: Downloads and prepares TCGA OV gene expression data.
# 2. Missing Data Handling: Filters out genes and samples with excessive missing values.
# 3. Imputation: Imputes missing gene expression values with the row mean.
# 4. Variance Filtering: Selects genes with high variance to prioritize biologically relevant genes.
# 5. Clinical Data Merging: Merges gene expression data with clinical data (e.g., survival status).
# 6. Log Transformation: Applies log transformation to the expression counts to normalize the data.
# 7. Reshaping: Reshapes the data into a wide format where genes are columns, and samples are rows.
# 8. File Saving: Saves the processed gene expression data to a CSV file.

# ==========================================================
# Download and prepare data from GDC (Gene Expression)
# The data includes gene expression (FPKM) and clinical information.
GDCdownload(query_OV_RNA_TP) 
tcga_OV_GE <- GDCprepare(query_OV_RNA_TP)  # Prepare OV gene expression data
dim(tcga_OV_GE)  # Check dimensions (number of features and samples)
colnames(colData(tcga_OV_GE))  # Check column names (includes both clinical and expression data)

# Extract gene expression data (FPKM values) from the prepared data
OV_matrix_GE <- assay(tcga_OV_GE, 'fpkm_unstrand')
# Pre-process the gene expression data: 
# - Remove features (genes) with a high percentage of missing values.
# - Remove samples with a high percentage of missing values.

# Calculate the percentage of missing values per gene and filter genes with missing data above a threshold
missing_percentage_features <- rowMeans(is.na(OV_matrix_GE)) 
selected_features <- which(missing_percentage_features <= missing_threshold)  # Only include features with missing data <= threshold
OV_matrix_GE_filtered_features <- OV_matrix_GE[selected_features, ]

# Calculate the percentage of missing values per sample and filter samples with missing data above a threshold
missing_percentage_samples <- colMeans(is.na(OV_matrix_GE_filtered_features))
selected_samples <- which(missing_percentage_samples <= missing_threshold)  # Only include samples with missing data <= threshold
OV_matrix_GE_filtered <- OV_matrix_GE_filtered_features[, selected_samples]
# Impute missing values (if any) by filling with the row mean
OV_matrix_GE_filtered <- t(apply(OV_matrix_GE_filtered, 1, impute_row_mean))
# Filter genes with low variance (focus on genes with higher variance, indicating more dynamic expression patterns)
gene_variances <- apply(OV_matrix_GE_filtered, 1, var)
variance_threshold <- quantile(gene_variances, probs = 0.95)  # Use the 95th percentile as the variance threshold
high_var_genes_indices <- which(gene_variances >= variance_threshold)  # Identify genes with high variance
high_var_genes <- OV_matrix_GE_filtered[gene_variances >= variance_threshold, ]  # Select high variance genes
# Check how many genes remain after filtering by variance
dim(high_var_genes)
# Retrieve gene metadata (gene names)
OV_gene_metadata <- as.data.frame(rowData(tcga_OV_GE))
OV_gene_data <- OV_gene_metadata[c('gene_id', 'gene_name')]
# Merge high variance genes with gene names and prepare the dataset for further analysis
OV_GE <- high_var_genes %>%
  as.data.frame() %>%
  rownames_to_column(var = 'gene_id') %>%
  gather(key = 'case_id', value = 'counts', -gene_id) %>%
  left_join(., OV_gene_data, by = "gene_id")
# Pre-process: Take only samples that are shared by all omics data (common samples)
OV_GE$case_id <- gsub('-01.*', '', OV_GE$case_id)  # Remove sample suffix for uniformity
# Filter for samples that are in both the gene expression and clinical data
OV_GE_data <- OV_GE %>%
  filter(case_id %in% common_samples)
# Add clinical information (e.g., survival status) to the gene expression data
OV_GE_data <- merge(OV_GE_data, clin_OV_data, by.x = 'case_id', by.y = 'submitter_id')
# Log-transformation of counts (FPKM) to normalize gene expression data
OV_GE_data$counts <- log(OV_GE_data$counts + 1)
# Clean and reshape the data: Remove unnecessary columns and handle missing values
OV_RNA <- OV_GE_data[,-1]  # Remove gene_id column
OV_RNA <- OV_GE_data[,-2]  # Remove case_id column
OV_RNA <- na.omit(OV_RNA)  # Remove rows with missing values
# Reshape data: Convert long format to wide format (genes as columns)
OV_RNA <- pivot_wider(OV_RNA, names_from = gene_name, values_from = counts, values_fn = mean)
# Create a directory to save the results (if it doesn't exist)
dir_name <- 'OV'
if (!dir.exists(dir_name)) {
  dir.create(dir_name)
}
# Save the processed gene expression data to a CSV file
write.csv(OV_RNA, 'OV/OV_GE_data.csv')

# ============================================================
# Pre-process Copy Number Variation Data (OV)
# ============================================================

# Description:
# This code handles the pre-processing of TCGA OV copy number variation (CNV) data. 
# It filters out genes and samples with a high percentage of missing data, imputes missing values, 
# and focuses on high-variance genes. The CNV data is then categorized, merged with clinical data, 
# reshaped into a wide format, and saved to a CSV file for downstream analysis.

# Key Steps:
# 1. Data Download & Preparation: Downloads and prepares TCGA OV CNV data.
# 2. Missing Data Handling: Filters out genes and samples with excessive missing values.
# 3. Imputation: Imputes missing CNV values by filling them with the row mean.
# 4. Variance Filtering: Selects genes with high variance to retain biologically relevant features.
# 5. CNV Categorization: Categorizes CNV data into relevant categories (e.g., deletion, amplification).
# 6. Parallel Processing: Implements parallel processing to speed up CNV categorization using multiple CPU cores.
# 7. Gene Metadata Retrieval: Retrieves gene metadata (gene IDs and gene names) for further analysis.
# 8. Clinical Data Merging: Merges CNV data with clinical information (e.g., survival status).
# 9. Data Reshaping: Reshapes data into a wide format with genes as columns and samples as rows.
# 10. Handle List Columns: Converts list columns into character vectors for consistency.
# 11. File Saving: Saves the processed CNV data to a CSV file for downstream analysis.

# ============================================================

# Download and prepare data from GDC (Copy Number Variation)
GDCdownload(query_OV_CNV) 
tcga_OV_CNV <- GDCprepare(query_OV_CNV)  # Prepare OV CNV data
dim(tcga_OV_CNV)  # Check dimensions (number of features and samples)

# Filter out the features (genes) with more than 20% missing values across all samples.
# Filter the samples that have more than 20% missing values across all features.

# Extract CNV data (copy number values) from the prepared data
OV_matrix_CNV <- assay(tcga_OV_CNV, 'copy_number')
# Calculate the percentage of missing values per gene (feature)
missing_percentage_features <- rowMeans(is.na(OV_matrix_CNV))
selected_features <- which(missing_percentage_features <= missing_threshold)  # Only include features with missing data <= threshold
OV_matrix_CNV_filtered_features <- OV_matrix_CNV[selected_features, ]
# Calculate the percentage of missing values per sample
missing_percentage_samples <- colMeans(is.na(OV_matrix_CNV_filtered_features))
selected_samples <- which(missing_percentage_samples <= missing_threshold)  # Only include samples with missing data <= threshold
OV_matrix_CNV_filtered <- OV_matrix_CNV_filtered_features[, selected_samples]

# Impute missing values (if any) by filling with the row mean
OV_matrix_CNV_filtered <- t(apply(OV_matrix_CNV_filtered, 1, impute_row_mean))
# Filter genes with low variance (focus on genes with higher variance, indicating more dynamic expression patterns)
gene_variances <- apply(OV_matrix_CNV_filtered, 1, var)
variance_threshold <- quantile(gene_variances, probs = 0.95)  # Use the 95th percentile as the variance threshold
high_var_genes_indices <- which(gene_variances >= variance_threshold)  # Identify genes with high variance
OV_matrix_CNV_filtered <- OV_matrix_CNV_filtered[gene_variances >= variance_threshold, ]  # Select high variance genes

# Setup parallel processing to apply categorization to CNV data
num_cores <- detectCores() - 1  # Use one less than the total cores available
cl <- makeCluster(num_cores)
# Export the necessary objects and functions to the cluster
clusterExport(cl, varlist = c("categorize_cnv", "OV_matrix_CNV_filtered"))
# Load the dplyr package on each worker
clusterEvalQ(cl, library(dplyr))
# Function to apply categorization to a matrix chunk
apply_categorization <- function(rows) {
  apply(OV_matrix_CNV_filtered[rows, , drop = FALSE], c(1, 2), categorize_cnv)
}
# Export the apply_categorization function to the cluster
clusterExport(cl, varlist = "apply_categorization")
# Split the matrix into chunks for parallel processing
split_indices <- split(seq_len(nrow(OV_matrix_CNV_filtered)), rep(1:num_cores, length.out = nrow(OV_matrix_CNV_filtered)))
# Apply categorization in parallel
OV_matrix_CNV_categorized_chunks <- parLapply(cl, split_indices, apply_categorization)
# Combine the results back into a single matrix
OV_matrix_CNV_categorized <- do.call(rbind, OV_matrix_CNV_categorized_chunks)
# Stop the cluster
stopCluster(cl)
# Convert the resulting matrix to a data frame if needed
OV_matrix_CNV_filtered <- as.matrix(OV_matrix_CNV_categorized)
# Retrieve gene metadata (gene names)
OV_gene_metadata <- as.data.frame(rowData(tcga_OV_CNV))
OV_gene_data <- OV_gene_metadata[c('gene_id', 'gene_name')]

# Merge categorized CNV data with gene names and prepare the dataset for further analysis
OV_CNV <- OV_matrix_CNV_filtered %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'copy_number', -gene_id) %>% 
  left_join(., OV_gene_data, by = "gene_id")

# Pre-process: Take only samples that are shared by all omics data (common samples)
OV_CNV$case_id <- gsub('-01.*', '', OV_CNV$case_id)  # Remove sample suffix for uniformity
# Filter for samples that are in both the CNV and clinical data
OV_CNV_data <- OV_CNV %>%
  filter(case_id %in% common_samples)
# Add clinical information (e.g., survival status) to the CNV data
OV_CNV_data <- merge(OV_CNV_data, clin_OV_data, by.x = 'case_id', by.y = 'submitter_id')
# Clean and reshape the data: Remove unnecessary columns and handle missing values
OV_CNV <- OV_CNV_data[,-1]  # Remove gene_id column
OV_CNV <- OV_CNV_data[,-2]  # Remove case_id column
OV_CNV <- na.omit(OV_CNV)  # Remove rows with missing values
# Reshape data: Convert long format to wide format (genes as columns)
OV_CNV <- pivot_wider(OV_CNV, names_from = gene_name, values_from = copy_number, values_fn = mean)
# Check if any columns are lists and convert them to character vectors
list_columns <- sapply(OV_CNV, function(x) is.list(x))
OV_CNV[list_columns] <- lapply(OV_CNV[list_columns], function(x) sapply(x, function(y) toString(y)))
# Save the processed CNV data to a CSV file
write.table(OV_CNV, 'OV/OV_CNV_data.csv', col.names = NA, sep = ",")

# ============================================================
# Pre-process DNA Methylation Data (OV)
# ============================================================

# Description:
# This code handles the pre-processing of TCGA OV DNA methylation data using the Illumina Human Methylation 27 platform. 
# It filters out probes and samples with excessive missing data, applies imputation to missing values, 
# and focuses on high-variance probes. The data is reshaped into a wide format, merged with clinical data, 
# and saved to a CSV file for further analysis.

# Key Steps:
# 1. Data Download & Preparation: Downloads and prepares TCGA OV DNA methylation data (27k probes).
# 2. Data Filtering: Filters out probes and samples with more than 50% missing values.
# 3. Imputation: Imputes missing methylation beta values by filling them with the row mean.
# 4. Variance Filtering: Selects probes with high variance to focus on biologically relevant methylation changes.
# 5. Data Reshaping: Converts the data into a wide format with probes as columns and samples as rows.
# 6. Clinical Data Merging: Merges methylation data with clinical data (e.g., survival status).
# 7. File Saving: Saves the processed methylation data to a CSV file for downstream analysis.

# ============================================================
                                                                        
# Download and prepare the DNA methylation data for ovarian cancer (OV)
GDCdownload(query_OV_Meth)
tcga_OV_Meth <- GDCprepare(query_OV_Meth)

# Get the dimensions of the data (number of features and samples)
dim(tcga_OV_Meth)
# Extract the methylation matrix (beta values)
OV_matrix_METH = assay(tcga_OV_Meth)
# Filter out features with more than 50% missing values
missing_percentage_features <- rowMeans(is.na(OV_matrix_METH))
selected_features <- which(missing_percentage_features <= 0.5)
OV_matrix_METH_filtered_features <- OV_matrix_METH[selected_features, ]
# Filter out samples with more than 50% missing values
missing_percentage_samples <- colMeans(is.na(OV_matrix_METH_filtered_features))
selected_samples <- which(missing_percentage_samples <= 0.5)
OV_matrix_METH_filtered <- OV_matrix_METH_filtered_features[, selected_samples]
# Apply mean imputation for missing data across rows
OV_matrix_METH_filtered <- t(apply(OV_matrix_METH_filtered, 1, impute_row_mean))
# Calculate gene variance and filter genes with the top 5% variance
gene_variances <- apply(OV_matrix_METH_filtered, 1, var)
variance_threshold <- quantile(gene_variances, probs = 0.95)
high_var_genes_indices <- which(gene_variances >= variance_threshold)
# Filter out the high variance genes
high_var_genes <- OV_matrix_METH_filtered[gene_variances >= variance_threshold, ]
dim(high_var_genes)
# Reshape the methylation data and merge with clinical data
OV_METH <- high_var_genes %>%
  as.data.frame() %>%
  rownames_to_column(var = 'CpG') %>%
  gather(key = 'case_id', value = 'beta value', -CpG) %>%
  pivot_wider(names_from = CpG, values_from = 'beta value', values_fn = mean)
# Clean up the case ID format to match common sample IDs
OV_METH$case_id <- gsub('-01.*', '', OV_METH$case_id)
# Filter data to include only samples common across all omics
OV_METH_data <- OV_METH %>% 
  filter(case_id %in% common_samples)
# Merge with clinical data
OV_METH_data <- merge(OV_METH_data, clin_OV_data, by.x = 'case_id', by.y = 'submitter_id')
# Reorder columns: case_id, clinical data, and methylation data
num_cols <- ncol(OV_METH_data)
OV_METH_data <- cbind(
  OV_METH_data[, 1], 
  OV_METH_data[, (num_cols - 1):num_cols], 
  OV_METH_data[, -c(1, (num_cols - 1):num_cols)]  
)
# Rename first column to "case_id"
colnames(OV_METH_data)[1] <- "case_id" 
# Remove any remaining missing values
OV_METH_data <- na.omit(OV_METH_data)
# Write the processed data to a CSV file
write.csv(OV_METH_data, 'OV/OV_METH_data.csv')

# ============================================================
# Pre-process MiRNA Data (OV)
# ============================================================

# Description:
# This code handles the pre-processing of TCGA OV miRNA expression data. It filters out miRNAs and samples 
# with excessive missing data, applies imputation to missing values, and log-transforms the miRNA counts. 
# The data is reshaped into a wide format, retaining miRNAs as features and samples as rows. Finally, it filters 
# miRNAs based on their expression levels (keeping only those expressed in more than 50% of samples with a count > 0 
# and >10% of samples with a count > 1), merges the data with clinical information, and saves the processed data to a CSV file.

# Key Steps:
# 1. Data Download & Preparation: Downloads and prepares TCGA OV miRNA expression data (reads per million).
# 2. Data Filtering: Filters out miRNAs and samples with excessive missing values.
# 3. Imputation: Imputes missing miRNA counts by filling them with the row mean.
# 4. Log Transformation: Applies a log transformation (log(RPM + 1)) to normalize miRNA counts.
# 5. Data Reshaping: Converts the data into a wide format with miRNAs as columns and samples as rows.
# 6. Expression Threshold Filtering: Retains miRNAs expressed in more than 50% of samples with counts > 0 
#    and >10% of samples with counts > 1.
# 7. Clinical Data Merging: Merges miRNA expression data with clinical data (e.g., survival status).
# 8. File Saving: Saves the processed miRNA expression data to a CSV file for further analysis.

# ============================================================
                                                                        
# Download and prepare data from GDC (miRNA expression)
GDCdownload(query_OV_ME) 
tcga_OV_ME <- GDCprepare(query_OV_ME)
dim(tcga_OV_ME)  # Check dimensions (features and samples)
# Obtain only reads per million (RPM) data
rpm <- startsWith(names(tcga_OV_ME), "reads_per_")
tcga_OV_ME_data <- tcga_OV_ME[, c(1, which(rpm))]
# Clean column names
tcga_OV_ME_columns <- names(tcga_OV_ME_data)
tcga_OV_ME_columns <- ifelse(seq_along(tcga_OV_ME_columns) == 1, tcga_OV_ME_columns, 
                               gsub('^reads_per_million_miRNA_mapped_', "", tcga_OV_ME_columns))
names(tcga_OV_ME_data) <- tcga_OV_ME_columns
# Convert data to matrix
OV_matrix_ME <- as.matrix(tcga_OV_ME_data)
dim(OV_matrix_ME)  # Check dimensions of the data
# Filter out features (miRNAs) and samples 
missing_percentage_features <- rowMeans(is.na(OV_matrix_ME))
selected_features <- which(missing_percentage_features <= missing_threshold)
OV_matrix_ME_filtered_features <- OV_matrix_ME[selected_features, ]
missing_percentage_samples <- colMeans(is.na(OV_matrix_ME_filtered_features))
selected_samples <- which(missing_percentage_samples <= missing_threshold)
OV_matrix_ME_filtered <- OV_matrix_ME_filtered_features[, selected_samples]

# Impute missing values (if any) by filling with the row mean
OV_matrix_ME_filtered <- t(apply(OV_matrix_ME_filtered, 1, impute_row_mean))
dim(OV_matrix_ME_filtered)  # Check dimensions after imputation

# Reshape the data into a long format (miRNAs as features, case_id as the sample identifier)
OV_ME <- OV_matrix_ME_filtered %>%
  as.data.frame() %>%
  gather(key = 'case_id', value = 'counts', -miRNA_ID)

# Clean up case IDs and filter for common samples
OV_ME$case_id <- gsub('-01.*', '', OV_ME$case_id)
OV_ME_data <- OV_ME %>% filter(case_id %in% common_samples)
# Merge miRNA data with clinical data
OV_ME_data <- merge(OV_ME_data, clin_OV_data, by.x = 'case_id', by.y = 'submitter_id')
OV_ME_data$counts <- as.numeric(OV_ME_data$counts)
# Log (RPM + 1) the miRNA counts
OV_ME_data$counts <- log(OV_ME_data$counts + 1)
# Remove rows with missing values
OV_ME <- na.omit(OV_ME_data)
# Reshape data into a wide format with miRNAs as columns and samples as rows
OV_ME <- pivot_wider(OV_ME, names_from = miRNA_ID, values_from = counts, values_fn = mean)
# Filter miRNAs based on expression thresholds
filter_ME <- OV_ME[, -c(1:3)]  # Exclude the first three columns (case_id, clinical data)
percent_gt_0 <- colMeans(filter_ME > 0, na.rm = TRUE) * 100  # Percentage of samples with counts > 0
percent_gt_1 <- colMeans(filter_ME > 1, na.rm = TRUE) * 100  # Percentage of samples with counts > 1

# Retain miRNAs expressed in more than 50% of samples with counts > 0 and >10% with counts > 1
selected_columns <- names(filter_ME)[percent_gt_0 > 50 & percent_gt_1 > 10]

# Subset the dataset with selected miRNAs
OV_ME_filtered <- OV_ME[, selected_columns]
dim(OV_ME_filtered)  # Check dimensions after filtering

# Append the first three columns back (case_id and clinical data)
appended_data <- cbind(OV_ME[, 1:3], OV_ME_filtered)
# Save the processed miRNA expression data to a CSV file
write.csv(appended_data, 'OV/OV_ME_data.csv')
