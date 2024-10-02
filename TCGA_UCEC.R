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

setwd("/srv/scratch/z5309282") 
source("functions.R")

missing_threshold <- 0.2

clin_UCEC <- GDCquery_clinic("TCGA-UCEC", "Clinical")

clin_UCEC$deceased <- ifelse(clin_UCEC$vital_status == "Alive", 0, 1)
# create an "overall survival" variable that is equal to days_to_death for dead patients
# and to days_to_last_follow_up for patients who are still alive

clin_UCEC$overall_survival <- ifelse(clin_UCEC$vital_status == "Alive",
                                     clin_UCEC$days_to_last_follow_up,
                                     clin_UCEC$days_to_death)

clin_UCEC_data <- clin_UCEC[c('submitter_id', 'deceased', 'overall_survival')]

dim(clin_UCEC_data)


# build a query to retrieve gene expression data ------------
query_UCEC_RNA_TP <- GDCquery(project = 'TCGA-UCEC',
                              data.category = 'Transcriptome Profiling',
                              experimental.strategy = 'RNA-Seq',
                              workflow.type = 'STAR - Counts',
                              sample.type = 'Primary Tumor',
                              access  = 'open')

Output_query_UCEC_RNA_TP <- getResults(query_UCEC_RNA_TP)
UCEC_GE_sample <- Output_query_UCEC_RNA_TP[c('cases.submitter_id')]

# build a query to retrieve Copy Number Variation ------------

query_UCEC_CNV <- GDCquery(project = 'TCGA-UCEC',
                           data.category = 'Copy Number Variation',
                           sample.type = 'Primary Tumor',
                           data.type = 'Gene Level Copy Number',
                           access = 'open')
Output_query_UCEC_CNV <- getResults(query_UCEC_CNV)
UCEC_CNV_sample <- Output_query_UCEC_CNV[c('cases.submitter_id')]


# build a query to retrieve DNA Methylation data ------------

query_UCEC_Meth <- GDCquery(project = 'TCGA-UCEC',
                            data.category = 'DNA Methylation',
                            platform = 'Illumina Human Methylation 450',
                            sample.type = 'Primary Tumor',
                            data.type = 'Methylation Beta Value',
                            access = 'open')

Output_query_UCEC_Meth <- getResults(query_UCEC_Meth)
UCEC_Meth_sample <- Output_query_UCEC_Meth[c('cases.submitter_id')]


# build a query to retrieve miRNA expression data ------------

query_UCEC_ME <- GDCquery(project = 'TCGA-UCEC',
                          data.category = 'Transcriptome Profiling',
                          experimental.strategy = 'miRNA-Seq',
                          workflow.type = 'BCGSC miRNA Profiling',
                          data.type = 'miRNA Expression Quantification',
                          sample.type = 'Primary Tumor',
                          access = 'open')

Output_query_UCEC_ME <- getResults(query_UCEC_ME)
UCEC_ME_sample <- Output_query_UCEC_ME[c('cases.submitter_id')]

# Get COMMON SAMPLES ACROSS ALL OMICS DATA
common_samples <- Reduce(intersect, list(UCEC_GE_sample[[1]], UCEC_CNV_sample[[1]], UCEC_Meth_sample[[1]], UCEC_ME_sample[[1]]))


# Pre-process gene expression data ------------------------------------


#GDCdownload(query_UCEC_RNA_TP)
tcga_UCEC_GE <- GDCprepare(query_UCEC_RNA_TP)
dim(tcga_UCEC_GE) # Gets HOW MANY FEATURES THEN SAMPLES
colnames(colData(tcga_UCEC_GE)) # GETS COLUMN NAMES -> BOTH CLINCIAL AND EXPRESSION DATA IS PRESENT
UCEC_matrix_GE <- assay(tcga_UCEC_GE, 'fpkm_unstrand')

# Pre-process -> remove missing values + get high gene variance

missing_percentage_features <- rowMeans(is.na(UCEC_matrix_GE))
selected_features <- which(missing_percentage_features <= missing_threshold)
UCEC_matrix_GE_filtered_features <- UCEC_matrix_GE[selected_features, ]
missing_percentage_samples <- colMeans(is.na(UCEC_matrix_GE_filtered_features))
selected_samples <- which(missing_percentage_samples <= missing_threshold)
UCEC_matrix_GE_filtered <- UCEC_matrix_GE_filtered_features[, selected_samples]


#  By focusing on genes with higher variance, you prioritize those that exhibit more dynamic expression patterns across samples. 
#  This can help identify genes that are likely to be biologically relevant or associated with specific conditions.
UCEC_matrix_GE_filtered <- t(apply(UCEC_matrix_GE_filtered, 1, impute_row_mean))
gene_variances <- apply(UCEC_matrix_GE_filtered, 1, var)
variance_threshold <- quantile(gene_variances, probs = 0.95)
high_var_genes_indices <- which(gene_variances >= variance_threshold)
high_var_genes <- UCEC_matrix_GE_filtered[gene_variances >= variance_threshold, ]


dim(high_var_genes) # GETS HOW MANY GENES REMAIN AFTER FILTERING

UCEC_gene_metadata <- as.data.frame(rowData(tcga_UCEC_GE)) # To get gene name
UCEC_gene_data <- UCEC_gene_metadata[c('gene_id', 'gene_name')]


UCEC_GE <- high_var_genes %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
  left_join(., UCEC_gene_data, by = "gene_id")

# Pre-process -> Take only samples that are shared by all omics data
UCEC_GE$case_id <- gsub('-01.*', '', UCEC_GE$case_id)

UCEC_GE_data<- UCEC_GE %>% 
  filter(case_id %in% common_samples)

# Add clinical information to UCEC_GE_data
UCEC_GE_data <- merge(UCEC_GE_data, clin_UCEC_data, by.x = 'case_id', by.y = 'submitter_id')

# Log (FPKM + 1) the counts
UCEC_GE_data$counts <- log(UCEC_GE_data$counts + 1)


UCEC_RNA <- UCEC_GE_data[,-1]
UCEC_RNA <- UCEC_GE_data[,-2]
UCEC_RNA <- na.omit(UCEC_RNA)

UCEC_RNA <- pivot_wider(UCEC_RNA, names_from = gene_name, values_from = counts, values_fn = mean)

dir_name <- 'UCEC'
if (!dir.exists(dir_name)) {
  dir.create(dir_name)
}

write.csv(UCEC_RNA,'UCEC/UCEC_GE_data.csv')

# Pre-process Copy Number Variation data ------------------------------------

#GDCdownload(query_UCEC_CNV)
tcga_UCEC_CNV <- GDCprepare(query_UCEC_CNV)
dim(tcga_UCEC_CNV) # Gets HOW MANY FEATURES THEN SAMPLES

# Filter out the features that had more than 20% missing values across all samples 
# and filtered the samples that had more than 20% missing values across all features

UCEC_matrix_CNV <- assay(tcga_UCEC_CNV, 'copy_number')

missing_percentage_features <- rowMeans(is.na(UCEC_matrix_CNV))
selected_features <- which(missing_percentage_features <= missing_threshold)
UCEC_matrix_CNV_filtered_features <- UCEC_matrix_CNV[selected_features, ]
missing_percentage_samples <- colMeans(is.na(UCEC_matrix_CNV_filtered_features))
selected_samples <- which(missing_percentage_samples <= missing_threshold)
UCEC_matrix_CNV_filtered <- UCEC_matrix_CNV_filtered_features[, selected_samples]

UCEC_matrix_CNV_filtered <- t(apply(UCEC_matrix_CNV_filtered, 1, impute_row_mean))
gene_variances <- apply(UCEC_matrix_CNV_filtered, 1, var)
variance_threshold <- quantile(gene_variances, probs = 0.95)
high_var_genes_indices <- which(gene_variances >= variance_threshold)
UCEC_matrix_CNV_filtered <- UCEC_matrix_CNV_filtered[gene_variances >= variance_threshold, ]


-------------------
  
# Setup parallel processing
num_cores <- detectCores() - 1  # Use one less than the total cores available
cl <- makeCluster(num_cores)
# Export the necessary objects and functions to the cluster
clusterExport(cl, varlist = c("categorize_cnv", "UCEC_matrix_CNV_filtered"))
# Load the dplyr package on each worker
clusterEvalQ(cl, library(dplyr))
# Function to apply categorization to a matrix chunk
apply_categorization <- function(rows) {
  apply(UCEC_matrix_CNV_filtered[rows, , drop = FALSE], c(1, 2), categorize_cnv)
}
# Export the apply_categorization function to the cluster
clusterExport(cl, varlist = "apply_categorization")
# Split the matrix into chunks for parallel processing
split_indices <- split(seq_len(nrow(UCEC_matrix_CNV_filtered)), rep(1:num_cores, length.out = nrow(UCEC_matrix_CNV_filtered)))
# Apply categorization in parallel
UCEC_matrix_CNV_categorized_chunks <- parLapply(cl, split_indices, apply_categorization)
# Combine the results back into a single matrix
UCEC_matrix_CNV_categorized <- do.call(rbind, UCEC_matrix_CNV_categorized_chunks)
# Stop the cluster
stopCluster(cl)
# Convert the resulting matrix to a data frame if needed
UCEC_matrix_CNV_filtered <- as.matrix(UCEC_matrix_CNV_categorized)

UCEC_gene_metadata <- as.data.frame(rowData(tcga_UCEC_CNV))
UCEC_gene_data <- UCEC_gene_metadata[c('gene_id', 'gene_name')]


UCEC_CNV <- UCEC_matrix_CNV_filtered %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'copy_number', -gene_id) %>% 
  left_join(., UCEC_gene_data, by = "gene_id")

UCEC_CNV$case_id <- gsub('-01.*', '', UCEC_CNV$case_id)
UCEC_CNV_data<- UCEC_CNV %>% 
  filter(case_id %in% common_samples)

UCEC_CNV_data <- merge(UCEC_CNV_data, clin_UCEC_data, by.x = 'case_id', by.y = 'submitter_id')

UCEC_CNV <- UCEC_CNV_data[,-1]
UCEC_CNV <- UCEC_CNV_data[,-2]
UCEC_CNV <- na.omit(UCEC_CNV)
UCEC_CNV <- pivot_wider(UCEC_CNV, names_from = gene_name, values_from = copy_number, values_fn = mean)

list_columns <- sapply(UCEC_CNV, function(x) is.list(x))

# Convert each list column to character vectors
UCEC_CNV[list_columns] <- lapply(UCEC_CNV[list_columns], function(x) sapply(x, function(y) toString(y)))


# Write the modified data frame to a CSV file
write.table(UCEC_CNV, 'UCEC/UCEC_CNV_data.csv', col.names = NA, sep = ",")


# Pre-process DNA Methylation data ------------------------------------

# build a query to retrieve DNA Methylation data 27 ------------

query_UCEC_Meth27 <- GDCquery(project = 'TCGA-UCEC',
                              data.category = 'DNA Methylation',
                              platform = 'Illumina Human Methylation 27',
                              sample.type = 'Primary Tumor',
                              data.type = 'Methylation Beta Value',
                              access = 'open')
#GDCdownload(query_UCEC_Meth27)
#GDCdownload(query_UCEC_Meth)
tcga_UCEC_Meth27 <- GDCprepare(query_UCEC_Meth27)
tcga_UCEC_Meth <- GDCprepare(query_UCEC_Meth)
dim(tcga_UCEC_Meth)
UCEC_matrix_METH=assay(tcga_UCEC_Meth) 
UCEC_matrix_METH27=assay(tcga_UCEC_Meth27) 

# Get only 27K probes for the 450K data
rownames_to_keep <- rownames(UCEC_matrix_METH27)
UCEC_matrix_METH <- UCEC_matrix_METH[rownames(UCEC_matrix_METH) %in% rownames_to_keep, ]

missing_percentage_features <- rowMeans(is.na(UCEC_matrix_METH))
selected_features <- which(missing_percentage_features <= 0.5)
UCEC_matrix_METH_filtered_features <- UCEC_matrix_METH[selected_features, ]
missing_percentage_samples <- colMeans(is.na(UCEC_matrix_METH_filtered_features))
selected_samples <- which(missing_percentage_samples <= 0.5)
UCEC_matrix_METH_filtered <- UCEC_matrix_METH_filtered_features[, selected_samples]

dim(UCEC_matrix_METH_filtered)

# Apply the custom function to each row of the dataframe
UCEC_matrix_METH_filtered <- t(apply(UCEC_matrix_METH_filtered, 1, impute_row_mean))

gene_variances <- apply(UCEC_matrix_METH_filtered, 1, var)
variance_threshold <- quantile(gene_variances, probs = 0.95)
high_var_genes_indices <- which(gene_variances >= variance_threshold)

high_var_genes <- UCEC_matrix_METH_filtered[gene_variances >= variance_threshold, ]
dim(high_var_genes) # GETS HOW MANY GENES REMAIN AFTER FILTERING


UCEC_METH <- high_var_genes %>%
  as.data.frame() %>%
  rownames_to_column(var = 'CpG') %>%
  gather(key = 'case_id', value = 'beta value', -CpG) %>%
  pivot_wider(names_from = CpG, values_from = 'beta value', values_fn = mean)

UCEC_METH$case_id <- gsub('-01.*', '', UCEC_METH$case_id)
UCEC_METH_data<- UCEC_METH %>% 
  filter(case_id %in% common_samples)

UCEC_METH_data <- merge(UCEC_METH_data, clin_UCEC_data, by.x = 'case_id', by.y = 'submitter_id')
# Get the number of columns in the matrix
num_cols <- ncol(UCEC_METH_data)

# Create a new matrix with the desired column order
UCEC_METH_data <- cbind(
  UCEC_METH_data[, 1],  # First column remains unchanged
  UCEC_METH_data[, (num_cols - 1):num_cols],  # Last two columns moved to positions 2 and 3
  UCEC_METH_data[, -c(1, (num_cols - 1):num_cols)]  # Remaining columns
)
# Assign column names
colnames(UCEC_METH_data)[1] <- "case_id" 

UCEC_METH_data <- na.omit(UCEC_METH_data)

write.csv(UCEC_METH_data,'UCEC/UCEC_METH_data.csv')

# Pre-process MiRNA data ---------------------------------------------

#GDCdownload(query_UCEC_ME)
tcga_UCEC_ME <- GDCprepare(query_UCEC_ME)
dim(tcga_UCEC_ME) # Gets HOW MANY FEATURES THEN SAMPLES
# Obtain only reads per million sample data
rpm <- startsWith(names(tcga_UCEC_ME), "reads_per_")
tcga_UCEC_ME_data <- tcga_UCEC_ME[, c(1, which(rpm))]
tcga_UCEC_ME_columns <- names(tcga_UCEC_ME_data)
tcga_UCEC_ME_columns <- ifelse(seq_along(tcga_UCEC_ME_columns) == 1, tcga_UCEC_ME_columns, gsub('^reads_per_million_miRNA_mapped_', "", tcga_UCEC_ME_columns))
names(tcga_UCEC_ME_data) <- tcga_UCEC_ME_columns
UCEC_matrix_ME <- as.matrix(tcga_UCEC_ME_data)
dim(UCEC_matrix_ME)

missing_percentage_features <- rowMeans(is.na(UCEC_matrix_ME))
selected_features <- which(missing_percentage_features <= missing_threshold)
UCEC_matrix_ME_filtered_features <- UCEC_matrix_ME[selected_features, ]
missing_percentage_samples <- colMeans(is.na(UCEC_matrix_ME_filtered_features))
selected_samples <- which(missing_percentage_samples <= missing_threshold)
UCEC_matrix_ME_filtered <- UCEC_matrix_ME_filtered_features[, selected_samples]

UCEC_matrix_ME_filtered <- t(apply(UCEC_matrix_ME_filtered, 1, impute_row_mean))

dim(UCEC_matrix_ME_filtered)
UCEC_ME <- UCEC_matrix_ME_filtered %>% 
  as.data.frame() %>% 
  gather(key = 'case_id', value = 'counts', -miRNA_ID)

UCEC_ME$case_id <- gsub('-01.*', '', UCEC_ME$case_id)
UCEC_ME_data<- UCEC_ME %>% 
  filter(case_id %in% common_samples)

UCEC_ME_data <- merge(UCEC_ME_data, clin_UCEC_data, by.x = 'case_id', by.y = 'submitter_id')
UCEC_ME_data$counts <- as.numeric(UCEC_ME_data$counts)

# Log (RPM + 1) the counts
UCEC_ME_data$counts <- log(UCEC_ME_data$counts + 1)
UCEC_ME <- na.omit(UCEC_ME_data)
UCEC_ME <- pivot_wider(UCEC_ME, names_from = miRNA_ID, values_from = counts, values_fn = mean)

filter_ME <- UCEC_ME[, -c(1:3)]
# Calculate the number of samples
num_samples <- nrow(filter_ME)

# Calculate the percentage of samples where values are greater than 0 and 1 for each column
percent_gt_0 <- colMeans(filter_ME > 0, na.rm = TRUE) * 100
percent_gt_1 <- colMeans(filter_ME > 1, na.rm = TRUE) * 100

# Retain columns that satisfy the conditions
selected_columns <- names(filter_ME)[percent_gt_0 > 50 & percent_gt_1 > 10]

# Subset the dataset with selected columns
UCEC_ME_filtered <- UCEC_ME[, selected_columns]
dim(UCEC_ME_filtered)
# Append the first three columns back
appended_data <- cbind(UCEC_ME[, 1:3], UCEC_ME_filtered)

write.csv(appended_data,'UCEC/UCEC_ME_data.csv')