library(UCSCXenaTools)
library(dplyr)
library(tidyverse)
library(readr)
library(dplyr)
library(tibble)
library(miRBaseConverter)

setwd("/path/to/pwd")
.libPaths("/path/to/Rlib")

# ---------- Step 1: Load BRCA clinical data ----------
brca_cohort = XenaData %>% 
  filter(XenaHostNames == "tcgaHub") %>% # select TCGA Hub
  XenaScan("TCGA Breast Cancer")   # select BRCA cohort

cli_query = brca_cohort %>%
  filter(DataSubtype == "phenotype") %>%  # select clinical dataset
  XenaGenerate() %>%  # generate a XenaHub object
  XenaQuery() %>%
  XenaDownload()

cli = XenaPrepare(cli_query)
# Only get relevant columns
brca_clinical <- cli$BRCA_clinicalMatrix[, c("sampleID", "vital_status", "days_to_last_followup", "days_to_death")]
# Create overall survival variable
brca_clinical$overall_survival <- ifelse(brca_clinical$vital_status == "LIVING",
                                         brca_clinical$days_to_last_followup,
                                         brca_clinical$days_to_death)
clin_BRCA <- brca_clinical



# ---------- Step 1: Load CESC clinical data ----------
cesc_cohort = XenaData %>% 
  filter(XenaHostNames == "tcgaHub") %>% # select TCGA Hub
  XenaScan("TCGA Cervical Cancer")   # select BRCA cohort

cli_query = cesc_cohort %>%
  filter(DataSubtype == "phenotype") %>%  # select clinical dataset
  XenaGenerate() %>%  # generate a XenaHub object
  XenaQuery() %>%
  XenaDownload()

cli = XenaPrepare(cli_query)
# Only get relevant columns
cesc_clinical <- cli$CESC_clinicalMatrix[, c("sampleID", "vital_status", "days_to_last_followup", "days_to_death")]
# Create overall survival variable
cesc_clinical$overall_survival <- ifelse(cesc_clinical$vital_status == "LIVING",
                                         cesc_clinical$days_to_last_followup,
                                         cesc_clinical$days_to_death)
clin_CESC <- cesc_clinical


# ---------- Step 1: Load CESC clinical data ----------
cesc_cohort = XenaData %>% 
  filter(XenaHostNames == "tcgaHub") %>% # select TCGA Hub
  XenaScan("TCGA Cervical Cancer")   # select BRCA cohort

cli_query = cesc_cohort %>%
  filter(DataSubtype == "phenotype") %>%  # select clinical dataset
  XenaGenerate() %>%  # generate a XenaHub object
  XenaQuery() %>%
  XenaDownload()

cli = XenaPrepare(cli_query)
# Only get relevant columns
cesc_clinical <- cli$CESC_clinicalMatrix[, c("sampleID", "vital_status", "days_to_last_followup", "days_to_death")]
# Create overall survival variable
cesc_clinical$overall_survival <- ifelse(cesc_clinical$vital_status == "LIVING",
                                         cesc_clinical$days_to_last_followup,
                                         cesc_clinical$days_to_death)
clin_CESC <- cesc_clinical


# ---------- Step 1: Load UCEC clinical data ----------
ucec_cohort = XenaData %>% 
  filter(XenaHostNames == "tcgaHub") %>% # select TCGA Hub
  XenaScan("TCGA Endometrioid Cancer")   # select BRCA cohort

cli_query = ucec_cohort %>%
  filter(DataSubtype == "phenotype") %>%  # select clinical dataset
  XenaGenerate() %>%  # generate a XenaHub object
  XenaQuery() %>%
  XenaDownload()

cli = XenaPrepare(cli_query)
# Only get relevant columns
ucec_clinical <- cli$UCEC_clinicalMatrix[, c("sampleID", "vital_status", "days_to_last_followup", "days_to_death")]
# Create overall survival variable
ucec_clinical$overall_survival <- ifelse(ucec_clinical$vital_status == "LIVING",
                                         ucec_clinical$days_to_last_followup,
                                         ucec_clinical$days_to_death)
clin_UCEC <- ucec_clinical


# ---------- Step 1: Load OV clinical data ----------
ov_cohort = XenaData %>% 
  filter(XenaHostNames == "tcgaHub") %>% # select TCGA Hub
  XenaScan("TCGA Ovarian Cancer")   # select BRCA cohort

cli_query = ov_cohort %>%
  filter(DataSubtype == "phenotype") %>%  # select clinical dataset
  XenaGenerate() %>%  # generate a XenaHub object
  XenaQuery() %>%
  XenaDownload()

cli = XenaPrepare(cli_query)
# Only get relevant columns
ov_clinical <- cli$OV_clinicalMatrix[, c("sampleID", "vital_status", "days_to_last_followup", "days_to_death")]
# Create overall survival variable
ov_clinical$overall_survival <- ifelse(ov_clinical$vital_status == "LIVING",
                                         ov_clinical$days_to_last_followup,
                                         ov_clinical$days_to_death)
clin_OV <- ov_clinical

# ---------- Step 2: Convert vital status to binary ----------
clin_BRCA$deceased <- ifelse(clin_BRCA$vital_status == "LIVING", 0, 1)
clin_UCEC$deceased <- ifelse(clin_UCEC$vital_status == "LIVING", 0, 1)
clin_CESC$deceased <- ifelse(clin_CESC$vital_status == "LIVING", 0, 1)
clin_OV$deceased <- ifelse(clin_OV$vital_status == "LIVING", 0, 1)

# ---------- Step 3: Get only tumor samples ----------
clin_BRCA_filtered <- clin_BRCA[grepl("-01$", clin_BRCA$sampleID), ]
clin_UCEC_filtered <- clin_UCEC[grepl("-01$", clin_UCEC$sampleID), ]
clin_CESC_filtered <- clin_CESC[grepl("-01$", clin_CESC$sampleID), ]
clin_OV_filtered <- clin_OV[grepl("-01$", clin_OV$sampleID), ]

# ---------- Step 4: Get the necessary survival columns ----------
clin_BRCA <- clin_BRCA_filtered[c('sampleID', 'deceased', 'overall_survival')]
clin_UCEC <- clin_UCEC_filtered[c('sampleID', 'deceased', 'overall_survival')]
clin_CESC <- clin_CESC_filtered[c('sampleID', 'deceased', 'overall_survival')]
clin_OV <- clin_OV_filtered[c('sampleID', 'deceased', 'overall_survival')]


process_omics_matrix <- function(expr_matrix, clinical_df, cancer_type, omics_name, output_dir) {

  # Transpose and clean
  expr_clean <- expr_matrix %>%
    column_to_rownames(var = colnames(expr_matrix)[1]) %>%
    t() %>%
    as.data.frame()
  
  expr_clean$sampleID <- rownames(expr_clean)
  expr_clean$sample_type_code <- substr(expr_clean$sampleID, 14, 15)
  
  # Join with clinical data
  merged <- expr_clean %>%
    left_join(clinical_df, by = "sampleID")
  
  # Reorder clinical columns first
  move_clinical_first <- function(df) {
    clinical_cols <- c('sampleID', 'deceased', 'overall_survival')
    other_cols <- setdiff(colnames(df), clinical_cols)
    df[, c(clinical_cols, other_cols)]
  }
  merged <- move_clinical_first(merged)
  
  # Filter for primary tumor samples (sample type code "01")
  primary_tumor <- merged %>%
    filter(sample_type_code == "01") %>%
    select(-sample_type_code)
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Construct output file name
  output_file <- file.path(output_dir, paste0(cancer_type, "_primary_tumor_", omics_name, "_with_clinical.csv"))
  
  # Write to CSV
  write_csv(primary_tumor, output_file)
  
  message("File saved to: ", output_file)
}


get_sample_type_distribution <- function(expr_matrix) {
  library(dplyr)
  library(tibble)
  
  # Transpose and clean
  expr_clean <- expr_matrix %>%
    column_to_rownames(var = colnames(expr_matrix)[1]) %>%
    t() %>%
    as.data.frame()
  
  # Extract sample ID and sample type code
  expr_clean$sampleID <- rownames(expr_clean)
  expr_clean$sample_type_code <- substr(expr_clean$sampleID, 14, 15)
  
  # Return table of sample_type_code
  return(table(expr_clean$sample_type_code))
}

# ---------- Step 5.1: Get the Gene Expression Data ----------
downloadTCGA(project = "BRCA", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = '/path/to')
expr_path <- "/path/to/TCGA.BRCA.sampleMap/HiSeqV2.gz"
BRCA_matrix_GE <- read_tsv(expr_path)

process_omics_matrix(
  expr_matrix = BRCA_matrix_GE,
  clinical_df = clin_BRCA,
  cancer_type = "BRCA",
  omics_name = "GE",
  output_dir = "mRNA_expression"
)


downloadTCGA(project = "CESC", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = '/path/to')
expr_path <- "/path/to/TCGA.CESC.sampleMap/HiSeqV2.gz"
CESC_matrix_GE <- read_tsv(expr_path)

process_omics_matrix(
  expr_matrix = CESC_matrix_GE,
  clinical_df = clin_CESC,
  cancer_type = "CESC",
  omics_name = "GE",
  output_dir = "mRNA_expression"
)

downloadTCGA(project = "UCEC", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = '/path/to')
expr_path <- "/path/to/TCGA.UCEC.sampleMap/HiSeqV2.gz"
UCEC_matrix_GE <- read_tsv(expr_path)

process_omics_matrix(
  expr_matrix = UCEC_matrix_GE,
  clinical_df = clin_UCEC,
  cancer_type = "UCEC",
  omics_name = "GE",
  output_dir = "mRNA_expression"
)

downloadTCGA(project = "OV", data_type = "Gene Expression RNASeq", file_type = "IlluminaHiSeq RNASeqV2", destdir = '/path/to')
expr_path <- "/path/to/TCGA.OV.sampleMap/HiSeqV2.gz"
OV_matrix_GE <- read_tsv(expr_path)

process_omics_matrix(
  expr_matrix = OV_matrix_GE,
  clinical_df = clin_OV,
  cancer_type = "OV",
  omics_name = "GE",
  output_dir = "mRNA_expression"
)


# ---------- Step 5.2: Get the CNV Data ----------
downloadTCGA(project = "BRCA", data_type = "Gene Level Copy Number", file_type = "Gistic2 thresholded", destdir = '/path/to')
expr_path <- "/path/to/TCGA.BRCA.sampleMap/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz"
BRCA_matrix_CNV <- read_tsv(expr_path)

process_omics_matrix(
  expr_matrix = BRCA_matrix_CNV,
  clinical_df = clin_BRCA,
  cancer_type = "BRCA",
  omics_name = "CNV",
  output_dir = "copy_number_variation"
)

downloadTCGA(project = "CESC", data_type = "Gene Level Copy Number", file_type = "Gistic2 thresholded", destdir = '/path/to')
expr_path <- "/path/to/TCGA.CESC.sampleMap/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz"
CESC_matrix_CNV <- read_tsv(expr_path)

process_omics_matrix(
  expr_matrix = CESC_matrix_CNV,
  clinical_df = clin_CESC,
  cancer_type = "CESC",
  omics_name = "CNV",
  output_dir = "copy_number_variation"
)

downloadTCGA(project = "UCEC", data_type = "Gene Level Copy Number", file_type = "Gistic2 thresholded", destdir = '/path/to')
expr_path <- "/path/to/PRISM/TCGA.UCEC.sampleMap/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz"
UCEC_matrix_CNV <- read_tsv(expr_path)

process_omics_matrix(
  expr_matrix = UCEC_matrix_CNV,
  clinical_df = clin_UCEC,
  cancer_type = "UCEC",
  omics_name = "CNV",
  output_dir = "copy_number_variation"
)

downloadTCGA(project = "OV", data_type = "Gene Level Copy Number", file_type = "Gistic2 thresholded", destdir = '/path/to')
expr_path <- "/path/to/TCGA.OV.sampleMap/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz"
OV_matrix_CNV <- read_tsv(expr_path)

process_omics_matrix(
  expr_matrix = OV_matrix_CNV,
  clinical_df = clin_OV,
  cancer_type = "OV",
  omics_name = "CNV",
  output_dir = "copy_number_variation"
)


# ---------- Step 5.3: Get the Methylation Data ----------
downloadTCGA(project = "BRCA", data_type = "DNA Methylation", file_type = "Methylation27K", destdir = '/path/to')
expr_path <- "/path/to/TCGA.BRCA.sampleMap/HumanMethylation27.gz"
BRCA_matrix_DM_27 <- read_tsv(expr_path)

process_omics_matrix(
  expr_matrix = BRCA_matrix_DM_27,
  clinical_df = clin_BRCA,
  cancer_type = "BRCA",
  omics_name = "DM_27",
  output_dir = "methylation"
)

downloadTCGA(project = "BRCA", data_type = "DNA Methylation", file_type = "Methylation450K", destdir = '/path/to')
expr_path <- "/path/to/TCGA.BRCA.sampleMap/HumanMethylation450.gz"
BRCA_matrix_DM_450 <- read_tsv(expr_path)

process_omics_matrix(
  expr_matrix = BRCA_matrix_DM_450,
  clinical_df = clin_BRCA,
  cancer_type = "BRCA",
  omics_name = "DM_450",
  output_dir = "methylation"
)

downloadTCGA(project = "CESC", data_type = "DNA Methylation", file_type = "Methylation450K", destdir = '/path/to')
expr_path <- "/path/to/TCGA.CESC.sampleMap/HumanMethylation450.gz"
CESC_matrix_DM_450 <- read_tsv(expr_path)

process_omics_matrix(
  expr_matrix = CESC_matrix_DM_450,
  clinical_df = clin_CESC,
  cancer_type = "CESC",
  omics_name = "DM_450",
  output_dir = "methylation"
)

downloadTCGA(project = "UCEC", data_type = "DNA Methylation", file_type = "Methylation27K", destdir = '/path/to')
expr_path <- "/path/to/TCGA.UCEC.sampleMap/HumanMethylation27.gz"
UCEC_matrix_DM_27 <- read_tsv(expr_path)

process_omics_matrix(
  expr_matrix = UCEC_matrix_DM_27,
  clinical_df = clin_UCEC,
  cancer_type = "UCEC",
  omics_name = "DM_27",
  output_dir = "methylation"
)

downloadTCGA(project = "UCEC", data_type = "DNA Methylation", file_type = "Methylation450K", destdir = '/path/to')
expr_path <- "/path/to/TCGA.UCEC.sampleMap/HumanMethylation450.gz"
UCEC_matrix_DM_450 <- read_tsv(expr_path)

process_omics_matrix(
  expr_matrix = UCEC_matrix_DM_450,
  clinical_df = clin_UCEC,
  cancer_type = "UCEC",
  omics_name = "DM_450",
  output_dir = "methylation"
)

downloadTCGA(project = "OV", data_type = "DNA Methylation", file_type = "Methylation27K", destdir = '/path/to')
expr_path <- "/path/to/TCGA.OV.sampleMap/HumanMethylation27.gz"
OV_matrix_DM_27 <- read_tsv(expr_path)

process_omics_matrix(
  expr_matrix = OV_matrix_DM_27,
  clinical_df = clin_OV,
  cancer_type = "OV",
  omics_name = "DM_27",
  output_dir = "methylation"
)


# ---------- Step 5.3: Get the miRNA Data ----------
downloadTCGA(project = "BRCA", data_type = "miRNA Mature Strand Expression RNASeq", file_type = "IlluminaHiSeq RNASeq", destdir = '/path/to')
expr_path <- "/path/to/TCGA.BRCA.sampleMap/miRNA_HiSeq_gene.gz"
BRCA_matrix_ME <- read_tsv(expr_path)

process_omics_matrix(
  expr_matrix = BRCA_matrix_ME,
  clinical_df = clin_BRCA,
  cancer_type = "BRCA",
  omics_name = "ME",
  output_dir = "miRNA_expression"
)

downloadTCGA(project = "CESC", data_type = "miRNA Mature Strand Expression RNASeq", file_type = "IlluminaHiSeq RNASeq", destdir = '/path/to')
expr_path <- "/path/to/TCGA.CESC.sampleMap/miRNA_HiSeq_gene.gz"
CESC_matrix_ME <- read_tsv(expr_path)

process_omics_matrix(
  expr_matrix = CESC_matrix_ME,
  clinical_df = clin_CESC,
  cancer_type = "CESC",
  omics_name = "ME",
  output_dir = "miRNA_expression"
)

downloadTCGA(project = "UCEC", data_type = "miRNA Mature Strand Expression RNASeq", file_type = "IlluminaHiSeq RNASeq", destdir = '/path/to')
expr_path <- "/path/to/TCGA.UCEC.sampleMap/miRNA_HiSeq_gene.gz"
UCEC_matrix_ME <- read_tsv(expr_path)

process_omics_matrix(
  expr_matrix = UCEC_matrix_ME,
  clinical_df = clin_UCEC,
  cancer_type = "UCEC",
  omics_name = "ME",
  output_dir = "miRNA_expression"
)

downloadTCGA(project = "OV", data_type = "miRNA Mature Strand Expression RNASeq", file_type = "IlluminaHiSeq RNASeq", destdir = '/path/to')
expr_path <- "/path/to/TCGA.OV.sampleMap/miRNA_HiSeq_gene.gz"
OV_matrix_ME <- read_tsv(expr_path)

process_omics_matrix(
  expr_matrix = OV_matrix_ME,
  clinical_df = clin_OV,
  cancer_type = "OV",
  omics_name = "ME",
  output_dir = "miRNA_expression"
)
