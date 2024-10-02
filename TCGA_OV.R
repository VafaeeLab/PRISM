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

clin_OV <- GDCquery_clinic("TCGA-OV", "Clinical")

clin_OV$deceased <- ifelse(clin_OV$vital_status == "Alive", 0, 1)
# create an "overall survival" variable that is equal to days_to_death for dead patients
# and to days_to_last_follow_up for patients who are still alive

clin_OV$overall_survival <- ifelse(clin_OV$vital_status == "Alive",
                                     clin_OV$days_to_last_follow_up,
                                     clin_OV$days_to_death)

clin_OV_data <- clin_OV[c('submitter_id', 'deceased', 'overall_survival')]

dim(clin_OV_data)


# build a query to retrieve gene expression data ------------
query_OV_RNA_TP <- GDCquery(project = 'TCGA-OV',
                              data.category = 'Transcriptome Profiling',
                              experimental.strategy = 'RNA-Seq',
                              workflow.type = 'STAR - Counts',
                              sample.type = 'Primary Tumor',
                              access  = 'open')

Output_query_OV_RNA_TP <- getResults(query_OV_RNA_TP)
OV_GE_sample <- Output_query_OV_RNA_TP[c('cases.submitter_id')]

# build a query to retrieve Copy Number Variation ------------

query_OV_CNV <- GDCquery(project = 'TCGA-OV',
                           data.category = 'Copy Number Variation',
                           sample.type = 'Primary Tumor',
                           data.type = 'Gene Level Copy Number',
                           access = 'open')
Output_query_OV_CNV <- getResults(query_OV_CNV)
OV_CNV_sample <- Output_query_OV_CNV[c('cases.submitter_id')]


# build a query to retrieve DNA Methylation data ------------

query_OV_Meth <- GDCquery(project = 'TCGA-OV',
                            data.category = 'DNA Methylation',
                            platform = 'Illumina Human Methylation 27',
                            sample.type = 'Primary Tumor',
                            data.type = 'Methylation Beta Value',
                            access = 'open')

Output_query_OV_Meth <- getResults(query_OV_Meth)
OV_Meth_sample <- Output_query_OV_Meth[c('cases.submitter_id')]


# build a query to retrieve miRNA expression data ------------

query_OV_ME <- GDCquery(project = 'TCGA-OV',
                          data.category = 'Transcriptome Profiling',
                          experimental.strategy = 'miRNA-Seq',
                          workflow.type = 'BCGSC miRNA Profiling',
                          data.type = 'miRNA Expression Quantification',
                          sample.type = 'Primary Tumor',
                          access = 'open')

Output_query_OV_ME <- getResults(query_OV_ME)
OV_ME_sample <- Output_query_OV_ME[c('cases.submitter_id')]


# Get COMMON SAMPLES ACROSS ALL OMICS DATA
common_samples <- Reduce(intersect, list(OV_GE_sample[[1]], OV_CNV_sample[[1]], OV_Meth_sample[[1]], OV_ME_sample[[1]]))


# Pre-process gene expression data ------------------------------------


# GDCdownload(query_OV_RNA_TP)
tcga_OV_GE <- GDCprepare(query_OV_RNA_TP)
dim(tcga_OV_GE) # Gets HOW MANY FEATURES THEN SAMPLES
colnames(colData(tcga_OV_GE)) # GETS COLUMN NAMES -> BOTH CLINCIAL AND EXPRESSION DATA IS PRESENT
OV_matrix_GE <- assay(tcga_OV_GE, 'fpkm_unstrand')

# Pre-process -> impute missing values + get high gene variance

missing_percentage_features <- rowMeans(is.na(OV_matrix_GE))
selected_features <- which(missing_percentage_features <= missing_threshold)
OV_matrix_GE_filtered_features <- OV_matrix_GE[selected_features, ]
missing_percentage_samples <- colMeans(is.na(OV_matrix_GE_filtered_features))
selected_samples <- which(missing_percentage_samples <= missing_threshold)
OV_matrix_GE_filtered <- OV_matrix_GE_filtered_features[, selected_samples]


#  By focusing on genes with higher variance, you prioritize those that exhibit more dynamic expression patterns across samples. 
#  This can help identify genes that are likely to be biologically relevant or associated with specific conditions.
# Apply the custom function to each row of the dataframe
OV_matrix_GE_filtered <- t(apply(OV_matrix_GE_filtered, 1, impute_row_mean))
gene_variances <- apply(OV_matrix_GE_filtered, 1, var)
variance_threshold <- quantile(gene_variances, probs = 0.95)
high_var_genes_indices <- which(gene_variances >= variance_threshold)
high_var_genes <- OV_matrix_GE_filtered[gene_variances >= variance_threshold, ]


dim(high_var_genes) # GETS HOW MANY GENES REMAIN AFTER FILTERING

OV_gene_metadata <- as.data.frame(rowData(tcga_OV_GE)) # To get gene name
OV_gene_data <- OV_gene_metadata[c('gene_id', 'gene_name')]


OV_GE <- high_var_genes %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
  left_join(., OV_gene_data, by = "gene_id")

# Pre-process -> Take only samples that are shared by all omics data
OV_GE$case_id <- gsub('-01.*', '', OV_GE$case_id)

OV_GE_data<- OV_GE %>% 
  filter(case_id %in% common_samples)

# Add clinical information to OV_GE_data
OV_GE_data <- merge(OV_GE_data, clin_OV_data, by.x = 'case_id', by.y = 'submitter_id')

# Log (FPKM + 1) the counts
OV_GE_data$counts <- log(OV_GE_data$counts + 1)

OV_RNA <- OV_GE_data[,-1]
OV_RNA <- OV_GE_data[,-2]
OV_RNA <- na.omit(OV_RNA)

OV_RNA <- pivot_wider(OV_RNA, names_from = gene_name, values_from = counts, values_fn = mean)

dir_name <- 'OV'
if (!dir.exists(dir_name)) {
  dir.create(dir_name)
}

write.csv(OV_RNA,'OV/OV_GE_data.csv')

# Pre-process Copy Number Variation data ------------------------------------

#GDCdownload(query_OV_CNV)
tcga_OV_CNV <- GDCprepare(query_OV_CNV)
dim(tcga_OV_CNV) # Gets HOW MANY FEATURES THEN SAMPLES

# Filter out the features that had more than 20% missing values across all samples 
# and filtered the samples that had more than 20% missing values across all features

OV_matrix_CNV <- assay(tcga_OV_CNV, 'copy_number')

missing_percentage_features <- rowMeans(is.na(OV_matrix_CNV))
selected_features <- which(missing_percentage_features <= missing_threshold)
OV_matrix_CNV_filtered_features <- OV_matrix_CNV[selected_features, ]
missing_percentage_samples <- colMeans(is.na(OV_matrix_CNV_filtered_features))
selected_samples <- which(missing_percentage_samples <= missing_threshold)
OV_matrix_CNV_filtered <- OV_matrix_CNV_filtered_features[, selected_samples]

OV_matrix_CNV_filtered <- t(apply(OV_matrix_CNV_filtered, 1, impute_row_mean))
gene_variances <- apply(OV_matrix_CNV_filtered, 1, var)
variance_threshold <- quantile(gene_variances, probs = 0.95)
high_var_genes_indices <- which(gene_variances >= variance_threshold)
OV_matrix_CNV_filtered <- OV_matrix_CNV_filtered[gene_variances >= variance_threshold, ]

-------------------
  
# Setup parallel processing
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

OV_gene_metadata <- as.data.frame(rowData(tcga_OV_CNV))
OV_gene_data <- OV_gene_metadata[c('gene_id', 'gene_name')]

OV_CNV <- OV_matrix_CNV_filtered %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'copy_number', -gene_id) %>% 
  left_join(., OV_gene_data, by = "gene_id")

OV_CNV$case_id <- gsub('-01.*', '', OV_CNV$case_id)
OV_CNV_data<- OV_CNV %>% 
  filter(case_id %in% common_samples)

OV_CNV_data <- merge(OV_CNV_data, clin_OV_data, by.x = 'case_id', by.y = 'submitter_id')

OV_CNV <- OV_CNV_data[,-1]
OV_CNV <- OV_CNV_data[,-2]
OV_CNV <- na.omit(OV_CNV)
OV_CNV <- pivot_wider(OV_CNV, names_from = gene_name, values_from = copy_number, values_fn = mean)

list_columns <- sapply(OV_CNV, function(x) is.list(x))

# Convert each list column to character vectors
OV_CNV[list_columns] <- lapply(OV_CNV[list_columns], function(x) sapply(x, function(y) toString(y)))


# Write the modified data frame to a CSV file
write.table(OV_CNV, 'OV/OV_CNV_data.csv', col.names = NA, sep = ",")


# Pre-process DNA Methylation data ------------------------------------

#GDCdownload(query_OV_Meth)
tcga_OV_Meth <- GDCprepare(query_OV_Meth)
dim(tcga_OV_Meth)
OV_matrix_METH=assay(tcga_OV_Meth) 

missing_percentage_features <- rowMeans(is.na(OV_matrix_METH))
selected_features <- which(missing_percentage_features <= 0.5)
OV_matrix_METH_filtered_features <- OV_matrix_METH[selected_features, ]
missing_percentage_samples <- colMeans(is.na(OV_matrix_METH_filtered_features))
selected_samples <- which(missing_percentage_samples <= 0.5)
OV_matrix_METH_filtered <- OV_matrix_METH_filtered_features[, selected_samples]

dim(OV_matrix_METH_filtered)

# Apply the custom function to each row of the dataframe
OV_matrix_METH_filtered <- t(apply(OV_matrix_METH_filtered, 1, impute_row_mean))

gene_variances <- apply(OV_matrix_METH_filtered, 1, var)
variance_threshold <- quantile(gene_variances, probs = 0.95)
high_var_genes_indices <- which(gene_variances >= variance_threshold)

high_var_genes <- OV_matrix_METH_filtered[gene_variances >= variance_threshold, ]
dim(high_var_genes) # GETS HOW MANY GENES REMAIN AFTER FILTERING


OV_METH <- high_var_genes %>%
  as.data.frame() %>%
  rownames_to_column(var = 'CpG') %>%
  gather(key = 'case_id', value = 'beta value', -CpG) %>%
  pivot_wider(names_from = CpG, values_from = 'beta value', values_fn = mean)

OV_METH$case_id <- gsub('-01.*', '', OV_METH$case_id)
OV_METH_data<- OV_METH %>% 
  filter(case_id %in% common_samples)

OV_METH_data <- merge(OV_METH_data, clin_OV_data, by.x = 'case_id', by.y = 'submitter_id')
# Get the number of columns in the matrix
num_cols <- ncol(OV_METH_data)

# Create a new matrix with the desired column order
OV_METH_data <- cbind(
  OV_METH_data[, 1],  # First column remains unchanged
  OV_METH_data[, (num_cols - 1):num_cols],  # Last two columns moved to positions 2 and 3
  OV_METH_data[, -c(1, (num_cols - 1):num_cols)]  # Remaining columns
)
# Assign column names
colnames(OV_METH_data)[1] <- "case_id" 

OV_METH_data <- na.omit(OV_METH_data)

write.csv(OV_METH_data,'OV/OV_METH_data.csv')

# Pre-process MiRNA data ---------------------------------------------

#GDCdownload(query_OV_ME)
tcga_OV_ME <- GDCprepare(query_OV_ME)
dim(tcga_OV_ME) # Gets HOW MANY FEATURES THEN SAMPLES
# Obtain only reads per million sample data
rpm <- startsWith(names(tcga_OV_ME), "reads_per_")
tcga_OV_ME_data <- tcga_OV_ME[, c(1, which(rpm))]
tcga_OV_ME_columns <- names(tcga_OV_ME_data)
tcga_OV_ME_columns <- ifelse(seq_along(tcga_OV_ME_columns) == 1, tcga_OV_ME_columns, gsub('^reads_per_million_miRNA_mapped_', "", tcga_OV_ME_columns))
names(tcga_OV_ME_data) <- tcga_OV_ME_columns
OV_matrix_ME <- as.matrix(tcga_OV_ME_data)
dim(OV_matrix_ME)

missing_percentage_features <- rowMeans(is.na(OV_matrix_ME))
selected_features <- which(missing_percentage_features <= missing_threshold)
OV_matrix_ME_filtered_features <- OV_matrix_ME[selected_features, ]
missing_percentage_samples <- colMeans(is.na(OV_matrix_ME_filtered_features))
selected_samples <- which(missing_percentage_samples <= missing_threshold)
OV_matrix_ME_filtered <- OV_matrix_ME_filtered_features[, selected_samples]

OV_matrix_ME_filtered <- t(apply(OV_matrix_ME_filtered, 1, impute_row_mean))

dim(OV_matrix_ME_filtered)
OV_ME <- OV_matrix_ME_filtered %>% 
  as.data.frame() %>% 
  gather(key = 'case_id', value = 'counts', -miRNA_ID)

OV_ME$case_id <- gsub('-01.*', '', OV_ME$case_id)
OV_ME_data<- OV_ME %>% 
  filter(case_id %in% common_samples)

OV_ME_data <- merge(OV_ME_data, clin_OV_data, by.x = 'case_id', by.y = 'submitter_id')
OV_ME_data$counts <- as.numeric(OV_ME_data$counts)

# Log (RPM + 1) the counts
OV_ME_data$counts <- log(OV_ME_data$counts + 1)
OV_ME <- na.omit(OV_ME_data)
OV_ME <- pivot_wider(OV_ME, names_from = miRNA_ID, values_from = counts, values_fn = mean)

filter_ME <- OV_ME[, -c(1:3)]
# Calculate the number of samples
num_samples <- nrow(filter_ME)

# Calculate the percentage of samples where values are greater than 0 and 1 for each column
percent_gt_0 <- colMeans(filter_ME > 0, na.rm = TRUE) * 100
percent_gt_1 <- colMeans(filter_ME > 1, na.rm = TRUE) * 100

# Retain columns that satisfy the conditions
selected_columns <- names(filter_ME)[percent_gt_0 > 50 & percent_gt_1 > 10]

# Subset the dataset with selected columns
OV_ME_filtered <- OV_ME[, selected_columns]
dim(OV_ME_filtered)
# Append the first three columns back
appended_data <- cbind(OV_ME[, 1:3], OV_ME_filtered)

write.csv(appended_data,'OV/OV_ME_data.csv')

