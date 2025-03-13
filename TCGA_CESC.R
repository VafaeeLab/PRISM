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
# Pre-process Clinical and Omics Data (CESC)
# ============================================================

# Description:
# This section prepares the clinical and omics data for TCGA CESC (Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma). 
# It processes clinical data by creating a "deceased" variable and an "overall survival" variable. Additionally, it builds queries 
# to retrieve gene expression, copy number variation, DNA methylation, and miRNA expression data, and identifies common samples 
# across all datasets. These common samples will be used for further analysis.

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
clin_CESC <- GDCquery_clinic("TCGA-CESC", "Clinical")
clin_CESC$deceased <- ifelse(clin_CESC$vital_status == "Alive", 0, 1)

# Create overall survival variable
clin_CESC$overall_survival <- ifelse(clin_CESC$vital_status == "Alive",
                                      clin_CESC$days_to_last_follow_up,
                                      clin_CESC$days_to_death)

# Select relevant columns for clinical data
clin_CESC_data <- clin_CESC[c('submitter_id', 'deceased', 'overall_survival')]
dim(clin_CESC_data)  # Check dimensions of the clinical data

# Build a query to retrieve gene expression data (RNA-Seq)
query_CESC_RNA_TP <- GDCquery(project = 'TCGA-CESC',
                              data.category = 'Transcriptome Profiling',
                              experimental.strategy = 'RNA-Seq',
                              workflow.type = 'STAR - Counts',
                              sample.type = 'Primary Tumor',
                              access  = 'open')

# Get results for gene expression data
Output_query_CESC_RNA_TP <- getResults(query_CESC_RNA_TP)
CESC_GE_sample <- Output_query_CESC_RNA_TP[c('cases.submitter_id')]

# Build a query to retrieve Copy Number Variation (CNV) data
query_CESC_CNV <- GDCquery(project = 'TCGA-CESC',
                           data.category = 'Copy Number Variation',
                           sample.type = 'Primary Tumor',
                           data.type = 'Gene Level Copy Number',
                           access = 'open')

# Get results for CNV data
Output_query_CESC_CNV <- getResults(query_CESC_CNV)
CESC_CNV_sample <- Output_query_CESC_CNV[c('cases.submitter_id')]

# Build a query to retrieve DNA Methylation data
query_CESC_Meth <- GDCquery(project = 'TCGA-CESC',
                            data.category = 'DNA Methylation',
                            platform = 'Illumina Human Methylation 450',
                            sample.type = 'Primary Tumor',
                            data.type = 'Methylation Beta Value',
                            access = 'open')

# Get results for DNA methylation data
Output_query_CESC_Meth <- getResults(query_CESC_Meth)
CESC_Meth_sample <- Output_query_CESC_Meth[c('cases.submitter_id')]

# Build a query to retrieve miRNA expression data
query_CESC_ME <- GDCquery(project = 'TCGA-CESC',
                          data.category = 'Transcriptome Profiling',
                          experimental.strategy = 'miRNA-Seq',
                          workflow.type = 'BCGSC miRNA Profiling',
                          data.type = 'miRNA Expression Quantification',
                          sample.type = 'Primary Tumor',
                          access = 'open')

# Get results for miRNA expression data
Output_query_CESC_ME <- getResults(query_CESC_ME)
CESC_ME_sample <- Output_query_CESC_ME[c('cases.submitter_id')]

# Identify common samples across all omics datasets
common_samples <- Reduce(intersect, list(CESC_GE_sample[[1]], 
                                         CESC_CNV_sample[[1]], 
                                         CESC_Meth_sample[[1]], 
                                         CESC_ME_sample[[1]]))

# ============================================================
# Pre-process Gene Expression Data (CESC)
# ============================================================

# Description:
# This section focuses on preprocessing the gene expression data for TCGA CESC (Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma).
# It retrieves FPKM-based gene expression data, filters out genes and samples with 20% missing values, selects genes with high variance, 
# and merges the expression data with clinical information. The processed gene expression data is then saved as a CSV file.

# Key Steps:
# 1. Gene Expression Data Retrieval and Preprocessing:
#    - Retrieves gene expression data using RNA-Seq (FPKM) from the GDC.
#    - Filters out genes and samples with excessive missing values.
#    - Selects genes with higher variance (top 5% of variance).
# 2. Merges Gene Expression Data with Clinical Information:
#    - Adds gene metadata and clinical data.
#    - Log-transforms the expression counts for normalization.
# 3. Filtering and Saving the Processed Data:
#    - Filters out samples that are not shared across all omics data.
#    - Saves the processed gene expression data in a CSV file.

# ============================================================

# Gene Expression Data Preprocessing
# Download and prepare gene expression data for TCGA-CESC
GDCdownload(query_CESC_RNA_TP)
tcga_CESC_GE <- GDCprepare(query_CESC_RNA_TP)
# Check the dimensions of the gene expression data (features and samples)
dim(tcga_CESC_GE) 
# Check column names to ensure clinical and expression data are included
colnames(colData(tcga_CESC_GE))
# Retrieve FPKM (unstranded) gene expression matrix
CESC_matrix_GE <- assay(tcga_CESC_GE, 'fpkm_unstrand')

# Preprocess: Filter features and samples with excessive missing values
missing_percentage_features <- rowMeans(is.na(CESC_matrix_GE))
selected_features <- which(missing_percentage_features <= missing_threshold)
CESC_matrix_GE_filtered_features <- CESC_matrix_GE[selected_features, ]

missing_percentage_samples <- colMeans(is.na(CESC_matrix_GE_filtered_features))
selected_samples <- which(missing_percentage_samples <= missing_threshold)
CESC_matrix_GE_filtered <- CESC_matrix_GE_filtered_features[, selected_samples]
# Select genes with high variance (top 5% of variance)
CESC_matrix_GE_filtered <- t(apply(CESC_matrix_GE_filtered, 1, impute_row_mean))
gene_variances <- apply(CESC_matrix_GE_filtered, 1, var)
variance_threshold <- quantile(gene_variances, probs = 0.95)
high_var_genes_indices <- which(gene_variances >= variance_threshold)
high_var_genes <- CESC_matrix_GE_filtered[gene_variances >= variance_threshold, ]

# Check how many genes remain after filtering
dim(high_var_genes)
# Retrieve gene metadata (gene_id and gene_name)
CESC_gene_metadata <- as.data.frame(rowData(tcga_CESC_GE))
CESC_gene_data <- CESC_gene_metadata[c('gene_id', 'gene_name')]
# Combine high variance genes with gene metadata
CESC_GE <- high_var_genes %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
  left_join(., CESC_gene_data, by = "gene_id")
# Filter out samples that are not shared across all omics data
CESC_GE$case_id <- gsub('-01.*', '', CESC_GE$case_id)
CESC_GE_data <- CESC_GE %>% 
  filter(case_id %in% common_samples)

# Add clinical information to CESC_GE_data
CESC_GE_data <- merge(CESC_GE_data, clin_CESC_data, by.x = 'case_id', by.y = 'submitter_id')
# Log-transform the expression counts (FPKM + 1)
CESC_GE_data$counts <- log(CESC_GE_data$counts + 1)
# Filter out unnecessary columns and handle missing values
CESC_RNA <- CESC_GE_data[,-1]
CESC_RNA <- CESC_GE_data[,-2]
CESC_RNA <- na.omit(CESC_RNA)
# Pivot the data to wide format (gene names as columns)
CESC_RNA <- pivot_wider(CESC_RNA, names_from = gene_name, values_from = counts, values_fn = mean)
# Create a directory for the output file if it doesn't exist
dir_name <- 'CESC'
if (!dir.exists(dir_name)) {
  dir.create(dir_name)
}
# Save the processed gene expression data as a CSV file
write.csv(CESC_RNA, 'CESC/CESC_GE_data.csv')

# ============================================================
# Pre-process Copy Number Variation Data (CESC)
# ============================================================

# Description:
# This code handles the pre-processing of TCGA CESC copy number variation (CNV) data. 
# It filters out genes and samples with a high percentage of missing data, imputes missing values, 
# and focuses on high-variance genes. The CNV data is then categorized, merged with clinical data, 
# reshaped into a wide format, and saved to a CSV file for downstream analysis.

# Key Steps:
# 1. Data Download & Preparation: Downloads and prepares TCGA CESC CNV data.
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
GDCdownload(query_CESC_CNV) 
tcga_CESC_CNV <- GDCprepare(query_CESC_CNV)  # Prepare CESC CNV data
dim(tcga_CESC_CNV)  # Check dimensions (number of features and samples)

# Filter out the features (genes) with more than 20% missing values across all samples.
# Filter the samples that have more than 20% missing values across all features.

# Extract CNV data (copy number values) from the prepared data
CESC_matrix_CNV <- assay(tcga_CESC_CNV, 'copy_number')
# Calculate the percentage of missing values per gene (feature)
missing_percentage_features <- rowMeans(is.na(CESC_matrix_CNV))
selected_features <- which(missing_percentage_features <= missing_threshold)  # Only include features with missing data <= threshold
CESC_matrix_CNV_filtered_features <- CESC_matrix_CNV[selected_features, ]
# Calculate the percentage of missing values per sample
missing_percentage_samples <- colMeans(is.na(CESC_matrix_CNV_filtered_features))
selected_samples <- which(missing_percentage_samples <= missing_threshold)  # Only include samples with missing data <= threshold
CESC_matrix_CNV_filtered <- CESC_matrix_CNV_filtered_features[, selected_samples]

# Impute missing values (if any) by filling with the row mean
CESC_matrix_CNV_filtered <- t(apply(CESC_matrix_CNV_filtered, 1, impute_row_mean))
# Filter genes with low variance (focus on genes with higher variance, indicating more dynamic expression patterns)
gene_variances <- apply(CESC_matrix_CNV_filtered, 1, var)
variance_threshold <- quantile(gene_variances, probs = 0.95)  # Use the 95th percentile as the variance threshold
high_var_genes_indices <- which(gene_variances >= variance_threshold)  # Identify genes with high variance
CESC_matrix_CNV_filtered <- CESC_matrix_CNV_filtered[gene_variances >= variance_threshold, ]  # Select high variance genes

# Setup parallel processing to apply categorization to CNV data
num_cores <- detectCores() - 1  # Use one less than the total cores available
cl <- makeCluster(num_cores)
# Export the necessary objects and functions to the cluster
clusterExport(cl, varlist = c("categorize_cnv", "CESC_matrix_CNV_filtered"))
# Load the dplyr package on each worker
clusterEvalQ(cl, library(dplyr))
# Function to apply categorization to a matrix chunk
apply_categorization <- function(rows) {
  apply(CESC_matrix_CNV_filtered[rows, , drop = FALSE], c(1, 2), categorize_cnv)
}
# Export the apply_categorization function to the cluster
clusterExport(cl, varlist = "apply_categorization")
# Split the matrix into chunks for parallel processing
split_indices <- split(seq_len(nrow(CESC_matrix_CNV_filtered)), rep(1:num_cores, length.out = nrow(CESC_matrix_CNV_filtered)))
# Apply categorization in parallel
CESC_matrix_CNV_categorized_chunks <- parLapply(cl, split_indices, apply_categorization)
# Combine the results back into a single matrix
CESC_matrix_CNV_categorized <- do.call(rbind, CESC_matrix_CNV_categorized_chunks)
# Stop the cluster
stopCluster(cl)
# Convert the resulting matrix to a data frame if needed
CESC_matrix_CNV_filtered <- as.matrix(CESC_matrix_CNV_categorized)
# Retrieve gene metadata (gene names)
CESC_gene_metadata <- as.data.frame(rowData(tcga_CESC_CNV))
CESC_gene_data <- CESC_gene_metadata[c('gene_id', 'gene_name')]

# Merge categorized CNV data with gene names and prepare the dataset for further analysis
CESC_CNV <- CESC_matrix_CNV_filtered %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'copy_number', -gene_id) %>% 
  left_join(., CESC_gene_data, by = "gene_id")

# Pre-process: Take only samples that are shared by all omics data (common samples)
CESC_CNV$case_id <- gsub('-01.*', '', CESC_CNV$case_id)  # Remove sample suffix for uniformity
# Filter for samples that are in both the CNV and clinical data
CESC_CNV_data <- CESC_CNV %>%
  filter(case_id %in% common_samples)
# Add clinical information (e.g., survival status) to the CNV data
CESC_CNV_data <- merge(CESC_CNV_data, clin_CESC_data, by.x = 'case_id', by.y = 'submitter_id')
# Clean and reshape the data: Remove unnecessary columns and handle missing values
CESC_CNV <- CESC_CNV_data[,-1]  # Remove gene_id column
CESC_CNV <- CESC_CNV_data[,-2]  # Remove case_id column
CESC_CNV <- na.omit(CESC_CNV)  # Remove rows with missing values
# Reshape data: Convert long format to wide format (genes as columns)
CESC_CNV <- pivot_wider(CESC_CNV, names_from = gene_name, values_from = copy_number, values_fn = mean)
# Check if any columns are lists and convert them to character vectors
list_columns <- sapply(CESC_CNV, function(x) is.list(x))
CESC_CNV[list_columns] <- lapply(CESC_CNV[list_columns], function(x) sapply(x, function(y) toString(y)))
# Save the processed CNV data to a CSV file
write.table(CESC_CNV, 'CESC/CESC_CNV_data.csv', col.names = NA, sep = ",")


# ============================================================
# Pre-process DNA Methylation Data (CESC)
# ============================================================

# Description:
# This code handles the pre-processing of TCGA CESC DNA methylation data using the Illumina Human Methylation 27 platform. 
# It filters out probes and samples with excessive missing data, applies imputation to missing values, 
# and focuses on high-variance probes. The data is reshaped into a wide format, merged with clinical data, 
# and saved to a CSV file for further analysis.

# Key Steps:
# 1. Data Download & Preparation: Downloads and prepares TCGA CESC DNA methylation data (27k probes).
# 2. Data Filtering: Filters out probes and samples with more than 50% missing values.
# 3. Imputation: Imputes missing methylation beta values by filling them with the row mean.
# 4. Variance Filtering: Selects probes with high variance to focus on biologically relevant methylation changes.
# 5. Data Reshaping: Converts the data into a wide format with probes as columns and samples as rows.
# 6. Clinical Data Merging: Merges methylation data with clinical data (e.g., survival status).
# 7. File Saving: Saves the processed methylation data to a CSV file for downstream analysis.

# ============================================================

# Download and prepare data from GDC (DNA Methylation 27)
query_CESC_Meth27 <- GDCquery(project = 'TCGA-CESC',
                              data.category = 'DNA Methylation',
                              platform = 'Illumina Human Methylation 27',
                              sample.type = 'Primary Tumor',
                              data.type = 'Methylation Beta Value',
                              access = 'open')
GDCdownload(query_CESC_Meth27) 
GDCdownload(query_CESC_Meth) 
# Prepare data for TCGA CESC methylation (both 27k and 450k platforms)                                                                          
tcga_CESC_Meth27 <- GDCprepare(query_CESC_Meth27)
tcga_CESC_Meth <- GDCprepare(query_CESC_Meth)
dim(tcga_CESC_Meth)  # Check dimensions (features and samples)
# Extract methylation data (beta values) for both datasets
CESC_matrix_METH <- assay(tcga_CESC_Meth)
CESC_matrix_METH27 <- assay(tcga_CESC_Meth27)
# Filter 450k data to include only 27k probes
rownames_to_keep <- rownames(CESC_matrix_METH27)
CESC_matrix_METH <- CESC_matrix_METH[rownames(CESC_matrix_METH) %in% rownames_to_keep, ]
# Handle missing data by filtering out features (probes) and samples with excessive missing values
missing_percentage_features <- rowMeans(is.na(CESC_matrix_METH))
selected_features <- which(missing_percentage_features <= 0.5)  # Filter features with >50% missing data
CESC_matrix_METH_filtered_features <- CESC_matrix_METH[selected_features, ]
missing_percentage_samples <- colMeans(is.na(CESC_matrix_METH_filtered_features))
selected_samples <- which(missing_percentage_samples <= 0.5)  # Filter samples with >50% missing data
CESC_matrix_METH_filtered <- CESC_matrix_METH_filtered_features[, selected_samples]
dim(CESC_matrix_METH_filtered)  # Check the size of the filtered data
# Impute missing values (if any) by filling with the row mean
CESC_matrix_METH_filtered <- t(apply(CESC_matrix_METH_filtered, 1, impute_row_mean))
# Focus on probes with high variance to retain biologically relevant methylation features
gene_variances <- apply(CESC_matrix_METH_filtered, 1, var)
variance_threshold <- quantile(gene_variances, probs = 0.95)  # Use 95th percentile as the variance threshold
high_var_genes_indices <- which(gene_variances >= variance_threshold)  # Identify high-variance probes
high_var_genes <- CESC_matrix_METH_filtered[gene_variances >= variance_threshold, ]  # Filter high-variance probes
dim(high_var_genes)  # Check how many probes remain after filtering

# Reshape the data into a wide format (probes as columns, samples as rows)
CESC_METH <- high_var_genes %>%
  as.data.frame() %>%
  rownames_to_column(var = 'CpG') %>%
  gather(key = 'case_id', value = 'beta value', -CpG) %>%
  pivot_wider(names_from = CpG, values_from = 'beta value', values_fn = mean)
# Clean up case IDs (remove suffix) and filter for common samples
CESC_METH$case_id <- gsub('-01.*', '', CESC_METH$case_id)
CESC_METH_data <- CESC_METH %>% filter(case_id %in% common_samples)
# Merge methylation data with clinical data
CESC_METH_data <- merge(CESC_METH_data, clin_CESC_data, by.x = 'case_id', by.y = 'submitter_id')
# Rearrange columns: Move clinical data to the beginning of the data frame
num_cols <- ncol(CESC_METH_data)
CESC_METH_data <- cbind(
  CESC_METH_data[, 1],  # Keep case_id in the first column
  CESC_METH_data[, (num_cols - 1):num_cols],  # Move the last two columns (clinical data) to the front
  CESC_METH_data[, -c(1, (num_cols - 1):num_cols)]  # Keep the remaining columns (methylation data)
)
# Assign column names for clarity
colnames(CESC_METH_data)[1] <- "case_id" 
# Remove rows with missing values
CESC_METH_data <- na.omit(CESC_METH_data)
# Save the processed methylation data to a CSV file
write.csv(CESC_METH_data, 'CESC/CESC_METH_data.csv')

# ============================================================
# Pre-process MiRNA Data (CESC)
# ============================================================

# Description:
# This code handles the pre-processing of TCGA CESC miRNA expression data. It filters out miRNAs and samples 
# with excessive missing data, applies imputation to missing values, and log-transforms the miRNA counts. 
# The data is reshaped into a wide format, retaining miRNAs as features and samples as rows. Finally, it filters 
# miRNAs based on their expression levels (keeping only those expressed in more than 50% of samples with a count > 0 
# and >10% of samples with a count > 1), merges the data with clinical information, and saves the processed data to a CSV file.

# Key Steps:
# 1. Data Download & Preparation: Downloads and prepares TCGA CESC miRNA expression data (reads per million).
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
GDCdownload(query_CESC_ME) 
tcga_CESC_ME <- GDCprepare(query_CESC_ME)
dim(tcga_CESC_ME)  # Check dimensions (features and samples)
# Obtain only reads per million (RPM) data
rpm <- startsWith(names(tcga_CESC_ME), "reads_per_")
tcga_CESC_ME_data <- tcga_CESC_ME[, c(1, which(rpm))]
# Clean column names
tcga_CESC_ME_columns <- names(tcga_CESC_ME_data)
tcga_CESC_ME_columns <- ifelse(seq_along(tcga_CESC_ME_columns) == 1, tcga_CESC_ME_columns, 
                               gsub('^reads_per_million_miRNA_mapped_', "", tcga_CESC_ME_columns))
names(tcga_CESC_ME_data) <- tcga_CESC_ME_columns
# Convert data to matrix
CESC_matrix_ME <- as.matrix(tcga_CESC_ME_data)
dim(CESC_matrix_ME)  # Check dimensions of the data
# Filter out features (miRNAs) and samples 
missing_percentage_features <- rowMeans(is.na(CESC_matrix_ME))
selected_features <- which(missing_percentage_features <= missing_threshold)
CESC_matrix_ME_filtered_features <- CESC_matrix_ME[selected_features, ]
missing_percentage_samples <- colMeans(is.na(CESC_matrix_ME_filtered_features))
selected_samples <- which(missing_percentage_samples <= missing_threshold)
CESC_matrix_ME_filtered <- CESC_matrix_ME_filtered_features[, selected_samples]

# Impute missing values (if any) by filling with the row mean
CESC_matrix_ME_filtered <- t(apply(CESC_matrix_ME_filtered, 1, impute_row_mean))
dim(CESC_matrix_ME_filtered)  # Check dimensions after imputation

# Reshape the data into a long format (miRNAs as features, case_id as the sample identifier)
CESC_ME <- CESC_matrix_ME_filtered %>%
  as.data.frame() %>%
  gather(key = 'case_id', value = 'counts', -miRNA_ID)

# Clean up case IDs and filter for common samples
CESC_ME$case_id <- gsub('-01.*', '', CESC_ME$case_id)
CESC_ME_data <- CESC_ME %>% filter(case_id %in% common_samples)
# Merge miRNA data with clinical data
CESC_ME_data <- merge(CESC_ME_data, clin_CESC_data, by.x = 'case_id', by.y = 'submitter_id')
CESC_ME_data$counts <- as.numeric(CESC_ME_data$counts)
# Log (RPM + 1) the miRNA counts
CESC_ME_data$counts <- log(CESC_ME_data$counts + 1)
# Remove rows with missing values
CESC_ME <- na.omit(CESC_ME_data)
# Reshape data into a wide format with miRNAs as columns and samples as rows
CESC_ME <- pivot_wider(CESC_ME, names_from = miRNA_ID, values_from = counts, values_fn = mean)
# Filter miRNAs based on expression thresholds
filter_ME <- CESC_ME[, -c(1:3)]  # Exclude the first three columns (case_id, clinical data)
percent_gt_0 <- colMeans(filter_ME > 0, na.rm = TRUE) * 100  # Percentage of samples with counts > 0
percent_gt_1 <- colMeans(filter_ME > 1, na.rm = TRUE) * 100  # Percentage of samples with counts > 1

# Retain miRNAs expressed in more than 50% of samples with counts > 0 and >10% with counts > 1
selected_columns <- names(filter_ME)[percent_gt_0 > 50 & percent_gt_1 > 10]

# Subset the dataset with selected miRNAs
CESC_ME_filtered <- CESC_ME[, selected_columns]
dim(CESC_ME_filtered)  # Check dimensions after filtering

# Append the first three columns back (case_id and clinical data)
appended_data <- cbind(CESC_ME[, 1:3], CESC_ME_filtered)
# Save the processed miRNA expression data to a CSV file
write.csv(appended_data, 'CESC/CESC_ME_data.csv')
