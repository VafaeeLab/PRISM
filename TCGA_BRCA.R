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

clin_BRCA <- GDCquery_clinic("TCGA-BRCA", "Clinical")

clin_BRCA$deceased <- ifelse(clin_BRCA$vital_status == "Alive", 0, 1)
# create an "overall survival" variable that is equal to days_to_death for dead patients
# and to days_to_last_follow_up for patients who are still alive

clin_BRCA$overall_survival <- ifelse(clin_BRCA$vital_status == "Alive",
                                     clin_BRCA$days_to_last_follow_up,
                                     clin_BRCA$days_to_death)

clin_BRCA_data <- clin_BRCA[c('submitter_id', 'deceased', 'overall_survival')]

dim(clin_BRCA_data)


# build a query to retrieve gene expression data ------------
query_BRCA_RNA_TP <- GDCquery(project = 'TCGA-BRCA',
                              data.category = 'Transcriptome Profiling',
                              experimental.strategy = 'RNA-Seq',
                              workflow.type = 'STAR - Counts',
                              sample.type = 'Primary Tumor',
                              access  = 'open')

Output_query_BRCA_RNA_TP <- getResults(query_BRCA_RNA_TP)
BRCA_GE_sample <- Output_query_BRCA_RNA_TP[c('cases.submitter_id')]

# build a query to retrieve Copy Number Variation ------------

query_BRCA_CNV <- GDCquery(project = 'TCGA-BRCA',
                           data.category = 'Copy Number Variation',
                           sample.type = 'Primary Tumor',
                           data.type = 'Gene Level Copy Number',
                           access = 'open')
Output_query_BRCA_CNV <- getResults(query_BRCA_CNV)
BRCA_CNV_sample <- Output_query_BRCA_CNV[c('cases.submitter_id')]


# build a query to retrieve DNA Methylation data ------------

query_BRCA_Meth <- GDCquery(project = 'TCGA-BRCA',
                            data.category = 'DNA Methylation',
                            platform = 'Illumina Human Methylation 450',
                            sample.type = 'Primary Tumor',
                            data.type = 'Methylation Beta Value',
                            access = 'open')

Output_query_BRCA_Meth <- getResults(query_BRCA_Meth)
BRCA_Meth_sample <- Output_query_BRCA_Meth[c('cases.submitter_id')]


# build a query to retrieve miRNA expression data ------------

query_BRCA_ME <- GDCquery(project = 'TCGA-BRCA',
                          data.category = 'Transcriptome Profiling',
                          experimental.strategy = 'miRNA-Seq',
                          workflow.type = 'BCGSC miRNA Profiling',
                          data.type = 'miRNA Expression Quantification',
                          sample.type = 'Primary Tumor',
                          access = 'open')

Output_query_BRCA_ME <- getResults(query_BRCA_ME)
BRCA_ME_sample <- Output_query_BRCA_ME[c('cases.submitter_id')]

# Get COMMON SAMPLES ACROSS ALL OMICS DATA
common_samples <- Reduce(intersect, list(BRCA_GE_sample[[1]], BRCA_CNV_sample[[1]], BRCA_Meth_sample[[1]], BRCA_ME_sample[[1]]))


# Pre-process gene expression data ------------------------------------


#GDCdownload(query_BRCA_RNA_TP)
tcga_BRCA_GE <- GDCprepare(query_BRCA_RNA_TP)
dim(tcga_BRCA_GE) # Gets HOW MANY FEATURES THEN SAMPLES
colnames(colData(tcga_BRCA_GE)) # GETS COLUMN NAMES -> BOTH CLINCIAL AND EXPRESSION DATA IS PRESENT
BRCA_matrix_GE <- assay(tcga_BRCA_GE, 'fpkm_unstrand')

# Pre-process -> remove missing values + get high gene variance

missing_percentage_features <- rowMeans(is.na(BRCA_matrix_GE))
selected_features <- which(missing_percentage_features <= missing_threshold)
BRCA_matrix_GE_filtered_features <- BRCA_matrix_GE[selected_features, ]
missing_percentage_samples <- colMeans(is.na(BRCA_matrix_GE_filtered_features))
selected_samples <- which(missing_percentage_samples <= missing_threshold)
BRCA_matrix_GE_filtered <- BRCA_matrix_GE_filtered_features[, selected_samples]


#  By focusing on genes with higher variance, you prioritize those that exhibit more dynamic expression patterns across samples. 
#  This can help identify genes that are likely to be biologically relevant or associated with specific conditions.
BRCA_matrix_GE_filtered <- t(apply(BRCA_matrix_GE_filtered, 1, impute_row_mean))
gene_variances <- apply(BRCA_matrix_GE_filtered, 1, var)
variance_threshold <- quantile(gene_variances, probs = 0.95)
high_var_genes_indices <- which(gene_variances >= variance_threshold)
high_var_genes <- BRCA_matrix_GE_filtered[gene_variances >= variance_threshold, ]


dim(high_var_genes) # GETS HOW MANY GENES REMAIN AFTER FILTERING

BRCA_gene_metadata <- as.data.frame(rowData(tcga_BRCA_GE)) # To get gene name
BRCA_gene_data <- BRCA_gene_metadata[c('gene_id', 'gene_name')]


BRCA_GE <- high_var_genes %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
  left_join(., BRCA_gene_data, by = "gene_id")

# Pre-process -> Take only samples that are shared by all omics data
BRCA_GE$case_id <- gsub('-01.*', '', BRCA_GE$case_id)

BRCA_GE_data<- BRCA_GE %>% 
  filter(case_id %in% common_samples)

# Add clinical information to BRCA_GE_data
BRCA_GE_data <- merge(BRCA_GE_data, clin_BRCA_data, by.x = 'case_id', by.y = 'submitter_id')

# Log (FPKM + 1) the counts
BRCA_GE_data$counts <- log(BRCA_GE_data$counts + 1)


BRCA_RNA <- BRCA_GE_data[,-1]
BRCA_RNA <- BRCA_GE_data[,-2]
BRCA_RNA <- na.omit(BRCA_RNA)

BRCA_RNA <- pivot_wider(BRCA_RNA, names_from = gene_name, values_from = counts, values_fn = mean)

dir_name <- 'BRCA'
if (!dir.exists(dir_name)) {
  dir.create(dir_name)
}

write.csv(BRCA_RNA,'BRCA/BRCA_GE_data.csv')

# Pre-process Copy Number Variation data ------------------------------------

#GDCdownload(query_BRCA_CNV)
tcga_BRCA_CNV <- GDCprepare(query_BRCA_CNV)
dim(tcga_BRCA_CNV) # Gets HOW MANY FEATURES THEN SAMPLES

# Filter out the features that had more than 20% missing values across all samples 
# and filtered the samples that had more than 20% missing values across all features

BRCA_matrix_CNV <- assay(tcga_BRCA_CNV, 'copy_number')

missing_percentage_features <- rowMeans(is.na(BRCA_matrix_CNV))
selected_features <- which(missing_percentage_features <= missing_threshold)
BRCA_matrix_CNV_filtered_features <- BRCA_matrix_CNV[selected_features, ]
missing_percentage_samples <- colMeans(is.na(BRCA_matrix_CNV_filtered_features))
selected_samples <- which(missing_percentage_samples <= missing_threshold)
BRCA_matrix_CNV_filtered <- BRCA_matrix_CNV_filtered_features[, selected_samples]

BRCA_matrix_CNV_filtered <- t(apply(BRCA_matrix_CNV_filtered, 1, impute_row_mean))
gene_variances <- apply(BRCA_matrix_CNV_filtered, 1, var)
variance_threshold <- quantile(gene_variances, probs = 0.95)
high_var_genes_indices <- which(gene_variances >= variance_threshold)
BRCA_matrix_CNV_filtered <- BRCA_matrix_CNV_filtered[gene_variances >= variance_threshold, ]


-------------------

# Setup parallel processing
num_cores <- detectCores() - 1  # Use one less than the total cores available
cl <- makeCluster(num_cores)
# Export the necessary objects and functions to the cluster
clusterExport(cl, varlist = c("categorize_cnv", "BRCA_matrix_CNV_filtered"))
# Load the dplyr package on each worker
clusterEvalQ(cl, library(dplyr))
# Function to apply categorization to a matrix chunk
apply_categorization <- function(rows) {
  apply(BRCA_matrix_CNV_filtered[rows, , drop = FALSE], c(1, 2), categorize_cnv)
}
# Export the apply_categorization function to the cluster
clusterExport(cl, varlist = "apply_categorization")
# Split the matrix into chunks for parallel processing
split_indices <- split(seq_len(nrow(BRCA_matrix_CNV_filtered)), rep(1:num_cores, length.out = nrow(BRCA_matrix_CNV_filtered)))
# Apply categorization in parallel
BRCA_matrix_CNV_categorized_chunks <- parLapply(cl, split_indices, apply_categorization)
# Combine the results back into a single matrix
BRCA_matrix_CNV_categorized <- do.call(rbind, BRCA_matrix_CNV_categorized_chunks)
# Stop the cluster
stopCluster(cl)
# Convert the resulting matrix to a data frame if needed
BRCA_matrix_CNV_filtered <- as.matrix(BRCA_matrix_CNV_categorized)

BRCA_gene_metadata <- as.data.frame(rowData(tcga_BRCA_CNV))
BRCA_gene_data <- BRCA_gene_metadata[c('gene_id', 'gene_name')]


BRCA_CNV <- BRCA_matrix_CNV_filtered %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'copy_number', -gene_id) %>% 
  left_join(., BRCA_gene_data, by = "gene_id")

BRCA_CNV$case_id <- gsub('-01.*', '', BRCA_CNV$case_id)
BRCA_CNV_data<- BRCA_CNV %>% 
  filter(case_id %in% common_samples)

BRCA_CNV_data <- merge(BRCA_CNV_data, clin_BRCA_data, by.x = 'case_id', by.y = 'submitter_id')

BRCA_CNV <- BRCA_CNV_data[,-1]
BRCA_CNV <- BRCA_CNV_data[,-2]
BRCA_CNV <- na.omit(BRCA_CNV)
BRCA_CNV <- pivot_wider(BRCA_CNV, names_from = gene_name, values_from = copy_number, values_fn = mean)

list_columns <- sapply(BRCA_CNV, function(x) is.list(x))

# Convert each list column to character vectors
BRCA_CNV[list_columns] <- lapply(BRCA_CNV[list_columns], function(x) sapply(x, function(y) toString(y)))


# Write the modified data frame to a CSV file
write.table(BRCA_CNV, 'BRCA/BRCA_CNV_data.csv', col.names = NA, sep = ",")


# Pre-process DNA Methylation data ------------------------------------

# build a query to retrieve DNA Methylation data 27 ------------

query_BRCA_Meth27 <- GDCquery(project = 'TCGA-BRCA',
                              data.category = 'DNA Methylation',
                              platform = 'Illumina Human Methylation 27',
                              sample.type = 'Primary Tumor',
                              data.type = 'Methylation Beta Value',
                              access = 'open')
#GDCdownload(query_BRCA_Meth27)
#GDCdownload(query_BRCA_Meth)
tcga_BRCA_Meth27 <- GDCprepare(query_BRCA_Meth27)
tcga_BRCA_Meth <- GDCprepare(query_BRCA_Meth)
dim(tcga_BRCA_Meth)
BRCA_matrix_METH=assay(tcga_BRCA_Meth) 
BRCA_matrix_METH27=assay(tcga_BRCA_Meth27) 

# Get only 27K probes for the 450K data
rownames_to_keep <- rownames(BRCA_matrix_METH27)
BRCA_matrix_METH <- BRCA_matrix_METH[rownames(BRCA_matrix_METH) %in% rownames_to_keep, ]

missing_percentage_features <- rowMeans(is.na(BRCA_matrix_METH))
selected_features <- which(missing_percentage_features <= 0.5)
BRCA_matrix_METH_filtered_features <- BRCA_matrix_METH[selected_features, ]
missing_percentage_samples <- colMeans(is.na(BRCA_matrix_METH_filtered_features))
selected_samples <- which(missing_percentage_samples <= 0.5)
BRCA_matrix_METH_filtered <- BRCA_matrix_METH_filtered_features[, selected_samples]

dim(BRCA_matrix_METH_filtered)

# Apply the custom function to each row of the dataframe
BRCA_matrix_METH_filtered <- t(apply(BRCA_matrix_METH_filtered, 1, impute_row_mean))

gene_variances <- apply(BRCA_matrix_METH_filtered, 1, var)
variance_threshold <- quantile(gene_variances, probs = 0.95)
high_var_genes_indices <- which(gene_variances >= variance_threshold)

high_var_genes <- BRCA_matrix_METH_filtered[gene_variances >= variance_threshold, ]
dim(high_var_genes) # GETS HOW MANY GENES REMAIN AFTER FILTERING


BRCA_METH <- high_var_genes %>%
  as.data.frame() %>%
  rownames_to_column(var = 'CpG') %>%
  gather(key = 'case_id', value = 'beta value', -CpG) %>%
  pivot_wider(names_from = CpG, values_from = 'beta value', values_fn = mean)

BRCA_METH$case_id <- gsub('-01.*', '', BRCA_METH$case_id)
BRCA_METH_data<- BRCA_METH %>% 
  filter(case_id %in% common_samples)

BRCA_METH_data <- merge(BRCA_METH_data, clin_BRCA_data, by.x = 'case_id', by.y = 'submitter_id')
# Get the number of columns in the matrix
num_cols <- ncol(BRCA_METH_data)

# Create a new matrix with the desired column order
BRCA_METH_data <- cbind(
  BRCA_METH_data[, 1],  # First column remains unchanged
  BRCA_METH_data[, (num_cols - 1):num_cols],  # Last two columns moved to positions 2 and 3
  BRCA_METH_data[, -c(1, (num_cols - 1):num_cols)]  # Remaining columns
)
# Assign column names
colnames(BRCA_METH_data)[1] <- "case_id" 

BRCA_METH_data <- na.omit(BRCA_METH_data)

write.csv(BRCA_METH_data,'BRCA/BRCA_METH_data.csv')

# Pre-process MiRNA data ---------------------------------------------

#GDCdownload(query_BRCA_ME)
tcga_BRCA_ME <- GDCprepare(query_BRCA_ME)
dim(tcga_BRCA_ME) # Gets HOW MANY FEATURES THEN SAMPLES
# Obtain only reads per million sample data
rpm <- startsWith(names(tcga_BRCA_ME), "reads_per_")
tcga_BRCA_ME_data <- tcga_BRCA_ME[, c(1, which(rpm))]
tcga_BRCA_ME_columns <- names(tcga_BRCA_ME_data)
tcga_BRCA_ME_columns <- ifelse(seq_along(tcga_BRCA_ME_columns) == 1, tcga_BRCA_ME_columns, gsub('^reads_per_million_miRNA_mapped_', "", tcga_BRCA_ME_columns))
names(tcga_BRCA_ME_data) <- tcga_BRCA_ME_columns
BRCA_matrix_ME <- as.matrix(tcga_BRCA_ME_data)
dim(BRCA_matrix_ME)

missing_percentage_features <- rowMeans(is.na(BRCA_matrix_ME))
selected_features <- which(missing_percentage_features <= missing_threshold)
BRCA_matrix_ME_filtered_features <- BRCA_matrix_ME[selected_features, ]
missing_percentage_samples <- colMeans(is.na(BRCA_matrix_ME_filtered_features))
selected_samples <- which(missing_percentage_samples <= missing_threshold)
BRCA_matrix_ME_filtered <- BRCA_matrix_ME_filtered_features[, selected_samples]

BRCA_matrix_ME_filtered <- t(apply(BRCA_matrix_ME_filtered, 1, impute_row_mean))
dim(BRCA_matrix_ME_filtered)
BRCA_ME <- BRCA_matrix_ME_filtered %>% 
  as.data.frame() %>% 
  gather(key = 'case_id', value = 'counts', -miRNA_ID)

BRCA_ME$case_id <- gsub('-01.*', '', BRCA_ME$case_id)
BRCA_ME_data<- BRCA_ME %>% 
  filter(case_id %in% common_samples)

BRCA_ME_data <- merge(BRCA_ME_data, clin_BRCA_data, by.x = 'case_id', by.y = 'submitter_id')
BRCA_ME_data$counts <- as.numeric(BRCA_ME_data$counts)

# Log (RPM + 1) the counts
BRCA_ME_data$counts <- log(BRCA_ME_data$counts + 1)
BRCA_ME <- na.omit(BRCA_ME_data)
BRCA_ME <- pivot_wider(BRCA_ME, names_from = miRNA_ID, values_from = counts, values_fn = mean)

filter_ME <- BRCA_ME[, -c(1:3)]
# Calculate the number of samples
num_samples <- nrow(filter_ME)

# Calculate the percentage of samples where values are greater than 0 and 1 for each column
percent_gt_0 <- colMeans(filter_ME > 0, na.rm = TRUE) * 100
percent_gt_1 <- colMeans(filter_ME > 1, na.rm = TRUE) * 100

# Retain columns that satisfy the conditions
selected_columns <- names(filter_ME)[percent_gt_0 > 50 & percent_gt_1 > 10]

# Subset the dataset with selected columns
BRCA_ME_filtered <- BRCA_ME[, selected_columns]
dim(BRCA_ME_filtered)
# Append the first three columns back
appended_data <- cbind(BRCA_ME[, 1:3], BRCA_ME_filtered)

write.csv(appended_data,'BRCA/BRCA_ME_data.csv')



# BRCA ----------------------------------------------------------------------

names_GE <- BRCA_GE$case_id
names_CNV <- BRCA_CNV$case_id
names_METH <- BRCA_METH$case_id
names_ME <- BRCA_ME$case_id
names_Clinical <- clin_BRCA$submitter_id

# Collect the names of each variable in a list
all_names <- list(
  RNA = names_GE,
  miRNA = names_ME,
  DNA_Methylation = names_METH,
  CNV = names_CNV,
  Clinical =names_Clinical
)
# Create the Venn diagram

category.names=c("GE","ME","DM","CNV", "Clinical")
venn.plot <- ggVennDiagram(
  x = all_names,
  category.names = category.names,
  label.col = "white",  # Label color
  cat.col = c("red", "blue", "green", "yellow", "orange"),  # Category colors
  margin = 0.1,
  category.names.fontface = "bold",
)
plot(venn.plot)
