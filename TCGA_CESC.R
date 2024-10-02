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

clin_CESC <- GDCquery_clinic("TCGA-CESC", "Clinical")

clin_CESC$deceased <- ifelse(clin_CESC$vital_status == "Alive", 0, 1)
# create an "CESCerall survival" variable that is equal to days_to_death for dead patients
# and to days_to_last_follow_up for patients who are still alive

clin_CESC$overall_survival <- ifelse(clin_CESC$vital_status == "Alive",
                                   clin_CESC$days_to_last_follow_up,
                                   clin_CESC$days_to_death)

clin_CESC_data <- clin_CESC[c('submitter_id', 'deceased', 'overall_survival')]

dim(clin_CESC_data)


# build a query to retrieve gene expression data ------------
query_CESC_RNA_TP <- GDCquery(project = 'TCGA-CESC',
                            data.category = 'Transcriptome Profiling',
                            experimental.strategy = 'RNA-Seq',
                            workflow.type = 'STAR - Counts',
                            sample.type = 'Primary Tumor',
                            access  = 'open')

Output_query_CESC_RNA_TP <- getResults(query_CESC_RNA_TP)
CESC_GE_sample <- Output_query_CESC_RNA_TP[c('cases.submitter_id')]

# build a query to retrieve Copy Number Variation ------------

query_CESC_CNV <- GDCquery(project = 'TCGA-CESC',
                         data.category = 'Copy Number Variation',
                         sample.type = 'Primary Tumor',
                         data.type = 'Gene Level Copy Number',
                         access = 'open')
Output_query_CESC_CNV <- getResults(query_CESC_CNV)
CESC_CNV_sample <- Output_query_CESC_CNV[c('cases.submitter_id')]


# build a query to retrieve DNA Methylation data ------------

query_CESC_Meth <- GDCquery(project = 'TCGA-CESC',
                          data.category = 'DNA Methylation',
                          platform = 'Illumina Human Methylation 450',
                          sample.type = 'Primary Tumor',
                          data.type = 'Methylation Beta Value',
                          access = 'open')

Output_query_CESC_Meth <- getResults(query_CESC_Meth)
CESC_Meth_sample <- Output_query_CESC_Meth[c('cases.submitter_id')]


# build a query to retrieve miRNA expression data ------------

query_CESC_ME <- GDCquery(project = 'TCGA-CESC',
                        data.category = 'Transcriptome Profiling',
                        experimental.strategy = 'miRNA-Seq',
                        workflow.type = 'BCGSC miRNA Profiling',
                        data.type = 'miRNA Expression Quantification',
                        sample.type = 'Primary Tumor',
                        access = 'open')

Output_query_CESC_ME <- getResults(query_CESC_ME)
CESC_ME_sample <- Output_query_CESC_ME[c('cases.submitter_id')]


# Get COMMON SAMPLES ACROSS ALL OMICS DATA
common_samples <- Reduce(intersect, list(CESC_GE_sample[[1]], CESC_CNV_sample[[1]], CESC_Meth_sample[[1]], CESC_ME_sample[[1]]))


# Pre-process gene expression data ------------------------------------

# GDCdownload(query_CESC_RNA_TP)
tcga_CESC_GE <- GDCprepare(query_CESC_RNA_TP)
dim(tcga_CESC_GE) # Gets HOW MANY FEATURES THEN SAMPLES
colnames(colData(tcga_CESC_GE)) # GETS COLUMN NAMES -> BOTH CLINCIAL AND EXPRESSION DATA IS PRESENT
CESC_matrix_GE <- assay(tcga_CESC_GE, 'fpkm_unstrand')

# Pre-process -> remCESCe missing values + get high gene variance

missing_percentage_features <- rowMeans(is.na(CESC_matrix_GE))
selected_features <- which(missing_percentage_features <= missing_threshold)
CESC_matrix_GE_filtered_features <- CESC_matrix_GE[selected_features, ]
missing_percentage_samples <- colMeans(is.na(CESC_matrix_GE_filtered_features))
selected_samples <- which(missing_percentage_samples <= missing_threshold)
CESC_matrix_GE_filtered <- CESC_matrix_GE_filtered_features[, selected_samples]


#  By focusing on genes with higher variance, you prioritize those that exhibit more dynamic expression patterns across samples. 
#  This can help identify genes that are likely to be biologically relevant or associated with specific conditions.
CESC_matrix_GE_filtered <- t(apply(CESC_matrix_GE_filtered, 1, impute_row_mean))
gene_variances <- apply(CESC_matrix_GE_filtered, 1, var)
variance_threshold <- quantile(gene_variances, probs = 0.95)
high_var_genes_indices <- which(gene_variances >= variance_threshold)
high_var_genes <- CESC_matrix_GE_filtered[gene_variances >= variance_threshold, ]


dim(high_var_genes) # GETS HOW MANY GENES REMAIN AFTER FILTERING

CESC_gene_metadata <- as.data.frame(rowData(tcga_CESC_GE)) # To get gene name
CESC_gene_data <- CESC_gene_metadata[c('gene_id', 'gene_name')]


CESC_GE <- high_var_genes %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
  left_join(., CESC_gene_data, by = "gene_id")

# Pre-process -> Take only samples that are shared by all omics data
CESC_GE$case_id <- gsub('-01.*', '', CESC_GE$case_id)

CESC_GE_data<- CESC_GE %>% 
  filter(case_id %in% common_samples)

# Add clinical information to CESC_GE_data
CESC_GE_data <- merge(CESC_GE_data, clin_CESC_data, by.x = 'case_id', by.y = 'submitter_id')

# Log (FPKM + 1) the counts
CESC_GE_data$counts <- log(CESC_GE_data$counts + 1)


CESC_RNA <- CESC_GE_data[,-1]
CESC_RNA <- CESC_GE_data[,-2]
CESC_RNA <- na.omit(CESC_RNA)

CESC_RNA <- pivot_wider(CESC_RNA, names_from = gene_name, values_from = counts, values_fn = mean)

dir_name <- 'CESC'
if (!dir.exists(dir_name)) {
  dir.create(dir_name)
}

write.csv(CESC_RNA,'CESC/CESC_GE_data.csv')

# Pre-process Copy Number Variation data ------------------------------------

#GDCdownload(query_CESC_CNV)
tcga_CESC_CNV <- GDCprepare(query_CESC_CNV)
dim(tcga_CESC_CNV) # Gets HOW MANY FEATURES THEN SAMPLES

# Filter out the features that had more than 20% missing values across all samples 
# and filtered the samples that had more than 20% missing values across all features

CESC_matrix_CNV <- assay(tcga_CESC_CNV, 'copy_number')

missing_percentage_features <- rowMeans(is.na(CESC_matrix_CNV))
selected_features <- which(missing_percentage_features <= missing_threshold)
CESC_matrix_CNV_filtered_features <- CESC_matrix_CNV[selected_features, ]
missing_percentage_samples <- colMeans(is.na(CESC_matrix_CNV_filtered_features))
selected_samples <- which(missing_percentage_samples <= missing_threshold)
CESC_matrix_CNV_filtered <- CESC_matrix_CNV_filtered_features[, selected_samples]

CESC_matrix_CNV_filtered <- t(apply(CESC_matrix_CNV_filtered, 1, impute_row_mean))
gene_variances <- apply(CESC_matrix_CNV_filtered, 1, var)
variance_threshold <- quantile(gene_variances, probs = 0.95)
high_var_genes_indices <- which(gene_variances >= variance_threshold)
CESC_matrix_CNV_filtered <- CESC_matrix_CNV_filtered[gene_variances >= variance_threshold, ]

-------------------
  
# Setup parallel processing
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

CESC_gene_metadata <- as.data.frame(rowData(tcga_CESC_CNV))
CESC_gene_data <- CESC_gene_metadata[c('gene_id', 'gene_name')]

CESC_CNV <- CESC_matrix_CNV_filtered %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'copy_number', -gene_id) %>% 
  left_join(., CESC_gene_data, by = "gene_id")

CESC_CNV$case_id <- gsub('-01.*', '', CESC_CNV$case_id)
CESC_CNV_data<- CESC_CNV %>% 
  filter(case_id %in% common_samples)

CESC_CNV_data <- merge(CESC_CNV_data, clin_CESC_data, by.x = 'case_id', by.y = 'submitter_id')

CESC_CNV <- CESC_CNV_data[,-1]
CESC_CNV <- CESC_CNV_data[,-2]
CESC_CNV <- na.omit(CESC_CNV)
CESC_CNV <- pivot_wider(CESC_CNV, names_from = gene_name, values_from = copy_number, values_fn = mean)

list_columns <- sapply(CESC_CNV, function(x) is.list(x))

# Convert each list column to character vectors
CESC_CNV[list_columns] <- lapply(CESC_CNV[list_columns], function(x) sapply(x, function(y) toString(y)))


# Write the modified data frame to a CSV file
write.table(CESC_CNV, 'CESC/CESC_CNV_data.csv', col.names = NA, sep = ",")


# Pre-process DNA Methylation data ------------------------------------

#GDCdownload(query_CESC_Meth)
tcga_CESC_Meth <- GDCprepare(query_CESC_Meth)
dim(tcga_CESC_Meth)
CESC_matrix_METH=assay(tcga_CESC_Meth) 

# build a query to retrieve DNA Methylation data 27 ------------

query_BRCA_Meth27 <- GDCquery(project = 'TCGA-BRCA',
                              data.category = 'DNA Methylation',
                              platform = 'Illumina Human Methylation 27',
                              sample.type = 'Primary Tumor',
                              data.type = 'Methylation Beta Value',
                              access = 'open')
#GDCdownload(query_BRCA_Meth27)
tcga_BRCA_Meth27 <- GDCprepare(query_BRCA_Meth27)
BRCA_matrix_METH27=assay(tcga_BRCA_Meth27) 

# Get only 27K probes for the 450K data
rownames_to_keep <- rownames(BRCA_matrix_METH27)
CESC_matrix_METH <- CESC_matrix_METH[rownames(CESC_matrix_METH) %in% rownames_to_keep, ]

missing_percentage_features <- rowMeans(is.na(CESC_matrix_METH))
selected_features <- which(missing_percentage_features <= 0.5)
CESC_matrix_METH_filtered_features <- CESC_matrix_METH[selected_features, ]
missing_percentage_samples <- colMeans(is.na(CESC_matrix_METH_filtered_features))
selected_samples <- which(missing_percentage_samples <= 0.5)
CESC_matrix_METH_filtered <- CESC_matrix_METH_filtered_features[, selected_samples]

dim(CESC_matrix_METH_filtered)
# Apply the custom function to each row of the dataframe
CESC_matrix_METH_filtered <- t(apply(CESC_matrix_METH_filtered, 1, impute_row_mean))

gene_variances <- apply(CESC_matrix_METH_filtered, 1, var)
variance_threshold <- quantile(gene_variances, probs = 0.95)
high_var_genes_indices <- which(gene_variances >= variance_threshold)

high_var_genes <- CESC_matrix_METH_filtered[gene_variances >= variance_threshold, ]
dim(high_var_genes) # GETS HOW MANY GENES REMAIN AFTER FILTERING


CESC_METH <- high_var_genes %>%
  as.data.frame() %>%
  rownames_to_column(var = 'CpG') %>%
  gather(key = 'case_id', value = 'beta value', -CpG) %>%
  pivot_wider(names_from = CpG, values_from = 'beta value', values_fn = mean)

CESC_METH$case_id <- gsub('-01.*', '', CESC_METH$case_id)
CESC_METH_data<- CESC_METH %>% 
  filter(case_id %in% common_samples)

CESC_METH_data <- merge(CESC_METH_data, clin_CESC_data, by.x = 'case_id', by.y = 'submitter_id')
# Get the number of columns in the matrix
num_cols <- ncol(CESC_METH_data)

# Create a new matrix with the desired column order
CESC_METH_data <- cbind(
  CESC_METH_data[, 1],  # First column remains unchanged
  CESC_METH_data[, (num_cols - 1):num_cols],  # Last two columns mCESCed to positions 2 and 3
  CESC_METH_data[, -c(1, (num_cols - 1):num_cols)]  # Remaining columns
)
# Assign column names
colnames(CESC_METH_data)[1] <- "case_id" 

CESC_METH_data <- na.omit(CESC_METH_data)

write.csv(CESC_METH_data,'CESC/CESC_METH_data.csv')

# Pre-process MiRNA data ---------------------------------------------

#GDCdownload(query_CESC_ME)
tcga_CESC_ME <- GDCprepare(query_CESC_ME)
dim(tcga_CESC_ME) # Gets HOW MANY FEATURES THEN SAMPLES
# Obtain only reads per million sample data
rpm <- startsWith(names(tcga_CESC_ME), "reads_per_")
tcga_CESC_ME_data <- tcga_CESC_ME[, c(1, which(rpm))]
tcga_CESC_ME_columns <- names(tcga_CESC_ME_data)
tcga_CESC_ME_columns <- ifelse(seq_along(tcga_CESC_ME_columns) == 1, tcga_CESC_ME_columns, gsub('^reads_per_million_miRNA_mapped_', "", tcga_CESC_ME_columns))
names(tcga_CESC_ME_data) <- tcga_CESC_ME_columns
CESC_matrix_ME <- as.matrix(tcga_CESC_ME_data)
dim(CESC_matrix_ME)

missing_percentage_features <- rowMeans(is.na(CESC_matrix_ME))
selected_features <- which(missing_percentage_features <= missing_threshold)
CESC_matrix_ME_filtered_features <- CESC_matrix_ME[selected_features, ]
missing_percentage_samples <- colMeans(is.na(CESC_matrix_ME_filtered_features))
selected_samples <- which(missing_percentage_samples <= missing_threshold)
CESC_matrix_ME_filtered <- CESC_matrix_ME_filtered_features[, selected_samples]

CESC_matrix_ME_filtered <- t(apply(CESC_matrix_ME_filtered, 1, impute_row_mean))

dim(CESC_matrix_ME_filtered)
CESC_ME <- CESC_matrix_ME_filtered %>% 
  as.data.frame() %>% 
  gather(key = 'case_id', value = 'counts', -miRNA_ID)

CESC_ME$case_id <- gsub('-01.*', '', CESC_ME$case_id)
CESC_ME_data<- CESC_ME %>% 
  filter(case_id %in% common_samples)

CESC_ME_data <- merge(CESC_ME_data, clin_CESC_data, by.x = 'case_id', by.y = 'submitter_id')
CESC_ME_data$counts <- as.numeric(CESC_ME_data$counts)

# Log (RPM + 1) the counts
CESC_ME_data$counts <- log(CESC_ME_data$counts + 1)
CESC_ME <- na.omit(CESC_ME_data)
CESC_ME <- pivot_wider(CESC_ME, names_from = miRNA_ID, values_from = counts, values_fn = mean)

filter_ME <- CESC_ME[, -c(1:3)]
# Calculate the number of samples
num_samples <- nrow(filter_ME)

# Calculate the percentage of samples where values are greater than 0 and 1 for each column
percent_gt_0 <- colMeans(filter_ME > 0, na.rm = TRUE) * 100
percent_gt_1 <- colMeans(filter_ME > 1, na.rm = TRUE) * 100

# Retain columns that satisfy the conditions
selected_columns <- names(filter_ME)[percent_gt_0 > 50 & percent_gt_1 > 10]

# Subset the dataset with selected columns
CESC_ME_filtered <- CESC_ME[, selected_columns]
dim(CESC_ME_filtered)
# Append the first three columns back
appended_data <- cbind(CESC_ME[, 1:3], CESC_ME_filtered)

write.csv(appended_data,'CESC/CESC_ME_data.csv')