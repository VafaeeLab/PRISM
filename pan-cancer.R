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
library(UpSetR)
library(enrichR)
library(grid)
library(dplyr)
library(multiMiR)
library(igraph)
library(ggraph)
library(tidygraph)
library(stringr)
library(IlluminaHumanMethylation27k.db)
library(enrichplot)
library(clusterProfiler)
library(DOSE)


source("functions.R")
set.seed(123)

# *********************************************************************
# Pan-cancer Analysis R Script Documentation
# *********************************************************************

# Overview:
# This script conducts a pan-cancer analysis for four cancer types (BRCA, OV, CESC, and UCEC) to uncover 
# shared multi-modal signatures across these cancers. The steps in the analysis are as follows:
# 1. **MiRA Targets**: Identifies unique & common miRNA gene targets across the four cancers 
# 2. **Gene Targets Disease Associations**: Using the gene targets, identifies further overlap
#    via diease associations
# 3. **GSEA**: Extracts the GO Terms and KEGG Pathways from the gene set
# 4. **GO TERMS AND KEGG PATHWAYS**: Filters for the top 10 GO terms & KEGG Pathways with the highest aggregated scores.
#    and plots it
# 5. **GO & KEGG Overlaps**: Plots the  top 10 GO & KEGG overlaps that appear in all four cancer, sorted by highest aggregated scores
# 6. **KM plots**: Creates KM plots for survival analysis
# This analysis aims to identify shared biological themes, potential overlaps in miRNA targets, and key pathways across 
# the cancers under study.

# *********************************************************************

# Load multi-omics feature-level fusion datasets for each cancer type
BRCA_data = read.csv("BRCA/LF/BRCA_METH_ME_CNV_data.csv", row.names = NULL)
OV_data = read.csv("OV/LF/OV_ME_METH_data.csv", row.names = NULL)
CESC_data = read.csv("CESC/LF/CESC_METH_ME_CNV_data.csv", row.names = NULL)
UCEC_data = read.csv("UCEC/LF/UCEC_METH_ME_CNV_data.csv", row.names = NULL)

# Load gene expression datasets for each cancer type
BRCA_gene_data = read.csv("BRCA/BRCA_GE_data.csv", row.names = NULL)
OV_gene_data = read.csv("OV/OV_GE_data.csv", row.names = NULL)
CESC_gene_data = read.csv("CESC/CESC_GE_data.csv", row.names = NULL)
UCEC_gene_data = read.csv("UCEC/UCEC_GE_data.csv", row.names = NULL)

# Remove the first column from gene expression datasets 
BRCA_gene_data <- BRCA_gene_data[, -1]
OV_gene_data <- OV_gene_data[, -1]
CESC_gene_data <- CESC_gene_data[,-1]
UCEC_gene_data <- UCEC_gene_data[, -1]


# Extract the 'variables' column
BRCA_vars <- BRCA_data$variables
OV_vars <- OV_data$variables
CESC_vars <- CESC_data$variables
UCEC_vars <- UCEC_data$variables

# Generate and save an UpSet plot to visualize the overlap of features between cancer types
listInput <- list(BRCA = BRCA_vars, OV = OV_vars, CESC = CESC_vars, UCEC = UCEC_vars)
png("upset_plot.png", width = 1600, height = 1200, res = 300)
upset(fromList(listInput), order.by = "freq")
grid.text("Pan-cancer Signature Overlap", x = 0.7, y = 0.98, gp = gpar(fontsize = 10, fontface = "bold"))
dev.off()
dev.new()

# Function to extract gene names for CNV features
extract_cnv_names <- function(vars) {
  gene_names <- sub("^CNV_", "", vars[grepl("^CNV_", vars)])
  return(gene_names)
}
# Function to extract gene names for ME features
extract_me_names <- function(vars) {
  gene_names <- sub("^ME_", "", vars[grepl("^ME_", vars)])
  return(gene_names)
}
# Function to extract gene names for METH features
extract_meth_names <- function(vars) {
  gene_names <- sub("^METH_", "", vars[grepl("^METH_", vars)])
  return(gene_names)
}

# Extract GE and CNV gene names for each cancer type
BRCA_cnv_genes <- extract_cnv_names(BRCA_vars)
OV_cnv_genes <- extract_cnv_names(OV_vars)
CESC_cnv_genes <- extract_cnv_names(CESC_vars)
UCEC_cnv_genes <- extract_cnv_names(UCEC_vars)
BRCA_me_genes <- extract_me_names(BRCA_vars)
OV_me_genes <- extract_me_names(OV_vars)
CESC_me_genes <- extract_me_names(CESC_vars)
UCEC_me_genes <- extract_me_names(UCEC_vars)
BRCA_meth_genes <- extract_meth_names(BRCA_vars)
OV_meth_genes <- extract_meth_names(OV_vars)
CESC_meth_genes <- extract_meth_names(CESC_vars)
UCEC_meth_genes <- extract_meth_names(UCEC_vars)



#################################################################################
#                             MiRA Targets
#################################################################################

BRCA_me_genes <- extract_me_names(BRCA_vars)
OV_me_genes <- extract_me_names(OV_vars)
CESC_me_genes <- extract_me_names(CESC_vars)
UCEC_me_genes <- extract_me_names(UCEC_vars)

# Function to convert miRNA IDs to hsa-miR-XX-5p or hsa-miR-XX-3p format
convert_miRNA_ids <- function(mirna_ids) {
  mirna_ids <- gsub("\\.", "-", mirna_ids)
  converted_ids <- c()
  for (id in mirna_ids) {
    converted_id_5p <- paste0(id, "-5p")
    converted_id_3p <- paste0(id, "-3p")
    converted_ids <- c(converted_ids, converted_id_5p, converted_id_3p)
  }
  return(converted_ids)
}

# Convert miRNA IDs for BRCA and retrieve target genes using multiMiR
BRCA_me_genes <- convert_miRNA_ids(BRCA_me_genes)
BRCA_gene_targets <- get_multimir(mirna = BRCA_me_genes,summary = TRUE)
# Extract mature miRNA IDs and their corresponding target gene symbols
mature_mirna_ids <- BRCA_gene_targets@summary$mature_mirna_id
target_symbols <- BRCA_gene_targets@summary$target_symbol
# Create a dataframe of miRNA-gene target interactions for BRCA
BRCA_gene_target_df <- data.frame(mature_mirna_id = mature_mirna_ids,
                 target_symbol = target_symbols)
# Filter out empty entries
BRCA_gene_target_df <- BRCA_gene_target_df %>%
  filter(mature_mirna_id != "" & target_symbol != "")

# Filter target genes to retain only those that are expressed in BRCA data
expressed_genes <- colnames(BRCA_gene_data)[4:ncol(BRCA_gene_data)]
BRCA_gene_target_df <- BRCA_gene_target_df[BRCA_gene_target_df$target_symbol %in% expressed_genes, ]

# Convert miRNA IDs for OV and retrieve target genes using multiMiR
OV_me_genes <- convert_miRNA_ids(OV_me_genes)
OV_gene_targets <- get_multimir(mirna = OV_me_genes, summary = TRUE)
# Extract mature miRNA IDs and their corresponding target gene symbols
mature_mirna_ids <- OV_gene_targets@summary$mature_mirna_id
target_symbols <- OV_gene_targets@summary$target_symbol
# Create a dataframe of miRNA-gene target interactions for OV
OV_gene_target_df <- data.frame(
  mature_mirna_id = mature_mirna_ids,
  target_symbol = target_symbols
)
# Filter out empty entries
OV_gene_target_df <- OV_gene_target_df %>%
  filter(mature_mirna_id != "" & target_symbol != "")
# Filter target genes to retain only those that are expressed in OV data
expressed_genes <- colnames(OV_gene_data)[4:ncol(OV_gene_data)]
OV_gene_target_df <- OV_gene_target_df[OV_gene_target_df$target_symbol %in% expressed_genes, ]

# Convert miRNA IDs for CESC and retrieve target genes using multiMiR
CESC_me_genes <- convert_miRNA_ids(CESC_me_genes)
CESC_gene_targets <- get_multimir(mirna = CESC_me_genes, summary = TRUE)
# Extract mature miRNA IDs and their corresponding target gene symbols
mature_mirna_ids <- CESC_gene_targets@summary$mature_mirna_id
target_symbols <- CESC_gene_targets@summary$target_symbol
# Create a dataframe of miRNA-gene target interactions for CESC
CESC_gene_target_df <- data.frame(
  mature_mirna_id = mature_mirna_ids,
  target_symbol = target_symbols
)
# Filter out empty entries
CESC_gene_target_df <- CESC_gene_target_df %>%
  filter(mature_mirna_id != "" & target_symbol != "")
# Filter target genes to retain only those that are expressed in CESC data
expressed_genes <- colnames(CESC_gene_data)[4:ncol(CESC_gene_data)]
CESC_gene_target_df <- CESC_gene_target_df[CESC_gene_target_df$target_symbol %in% expressed_genes, ]
# Convert miRNA IDs for UCEC and retrieve target genes using multiMiR
UCEC_me_genes <- convert_miRNA_ids(UCEC_me_genes)
UCEC_gene_targets <- get_multimir(mirna = UCEC_me_genes, summary = TRUE)
# Extract mature miRNA IDs and their corresponding target gene symbols
mature_mirna_ids <- UCEC_gene_targets@summary$mature_mirna_id
target_symbols <- UCEC_gene_targets@summary$target_symbol
# Create a dataframe of miRNA-gene target interactions for UCEC
UCEC_gene_target_df <- data.frame(
  mature_mirna_id = mature_mirna_ids,
  target_symbol = target_symbols
)
# Filter out empty entries
UCEC_gene_target_df <- UCEC_gene_target_df %>%
  filter(mature_mirna_id != "" & target_symbol != "")
# Filter target genes to retain only those that are expressed in UCEC data
expressed_genes <- colnames(UCEC_gene_data)[4:ncol(UCEC_gene_data)]
UCEC_gene_target_df <- UCEC_gene_target_df[UCEC_gene_target_df$target_symbol %in% expressed_genes, ]

# Create a list of target genes for each cancer type
listInput <- list(
  BRCA = BRCA_gene_target_df$target_symbol,
  OV = OV_gene_target_df$target_symbol,
  CESC = CESC_gene_target_df$target_symbol,
  UCEC = UCEC_gene_target_df$target_symbol
)

# Generate an UpSet plot to visualize shared and unique miRNA-target interactions across cancer types
upset_data <- fromList(listInput)
upset(upset_data, order.by = "freq")

# Combine the gene targets from all four data frames
combined_targets <- bind_rows(
  BRCA_gene_target_df %>% dplyr::select(target_symbol),
  OV_gene_target_df %>% dplyr::select(target_symbol),
  CESC_gene_target_df %>% dplyr::select(target_symbol),
  UCEC_gene_target_df %>% dplyr::select(target_symbol)
)

# Count the frequency of each gene target
gene_freq <- combined_targets %>%
  group_by(target_symbol) %>%
  summarise(frequency = n()) %>%
  arrange(desc(frequency))

# Identify common genes in all four data frames
common_genes <- Reduce(intersect, list(
  BRCA_gene_target_df$target_symbol,
  OV_gene_target_df$target_symbol,
  CESC_gene_target_df$target_symbol,
  UCEC_gene_target_df$target_symbol
))

# Combine the gene targets from all four data frames
combined_targets <- bind_rows(
  BRCA_gene_target_df %>% filter(target_symbol %in% common_genes) %>% dplyr::select(target_symbol),
  OV_gene_target_df %>% filter(target_symbol %in% common_genes) %>% dplyr::select(target_symbol),
  CESC_gene_target_df %>% filter(target_symbol %in% common_genes) %>% dplyr::select(target_symbol),
  UCEC_gene_target_df %>% filter(target_symbol %in% common_genes) %>% dplyr::select(target_symbol)
)

# Count the frequency of each common gene target
common_gene_freq <- combined_targets %>%
  group_by(target_symbol) %>%
  summarise(frequency = n()) %>%
  arrange(desc(frequency))



#################################################################################
#                      Gene Targets Disease Associations
#################################################################################

# Retrieve miRNA-disease associations for BRCA using multiMiR's 'mir2disease' table
BRCA_DD_targets <- get_multimir(mirna = BRCA_me_genes, table = 'mir2disease', summary = TRUE)
# Extract relevant columns: miRNA IDs, associated diseases/drugs, and PubMed IDs
mature_mirna_ids <- BRCA_DD_targets@data$mature_mirna_id
disease_drug <- BRCA_DD_targets@data$disease_drug
paper_pubmedID <- BRCA_DD_targets@data$paper_pubmedID
# Create a dataframe for BRCA miRNA-disease associations
BRCA_DD_target_df <- data.frame(
  mature_mirna_id = mature_mirna_ids,
  disease_drug = disease_drug,
  paper_pubmedID = paper_pubmedID
)
# Filter out empty entries
BRCA_DD_target_df <- BRCA_DD_target_df %>%
  filter(mature_mirna_id != "" & disease_drug != "" & paper_pubmedID != "")
# Retrieve miRNA-disease associations for OV
OV_DD_targets <- get_multimir(mirna = OV_me_genes, table = 'mir2disease', summary = TRUE)
# Extract relevant columns
mature_mirna_ids <- OV_DD_targets@data$mature_mirna_id
disease_drug <- OV_DD_targets@data$disease_drug
paper_pubmedID <- OV_DD_targets@data$paper_pubmedID
# Create a dataframe for OV miRNA-disease associations
OV_DD_target_df <- data.frame(
  mature_mirna_id = mature_mirna_ids,
  disease_drug = disease_drug,
  paper_pubmedID = paper_pubmedID
)
# Filter out empty entries
OV_DD_target_df <- OV_DD_target_df %>%
  filter(mature_mirna_id != "" & disease_drug != "" & paper_pubmedID != "")

# Retrieve miRNA-disease associations for CESC
CESC_DD_targets <- get_multimir(mirna = CESC_me_genes, table = 'mir2disease', summary = TRUE)
# Extract relevant columns
mature_mirna_ids <- CESC_DD_targets@data$mature_mirna_id
disease_drug <- CESC_DD_targets@data$disease_drug
paper_pubmedID <- CESC_DD_targets@data$paper_pubmedID
# Create a dataframe for CESC miRNA-disease associations
CESC_DD_target_df <- data.frame(
  mature_mirna_id = mature_mirna_ids,
  disease_drug = disease_drug,
  paper_pubmedID = paper_pubmedID
)
# Filter out empty entries
CESC_DD_target_df <- CESC_DD_target_df %>%
  filter(mature_mirna_id != "" & disease_drug != "" & paper_pubmedID != "")

# Retrieve miRNA-disease associations for UCEC
UCEC_DD_targets <- get_multimir(mirna = UCEC_me_genes, table = 'mir2disease', summary = TRUE)
# Extract relevant columns
mature_mirna_ids <- UCEC_DD_targets@data$mature_mirna_id
disease_drug <- UCEC_DD_targets@data$disease_drug
paper_pubmedID <- UCEC_DD_targets@data$paper_pubmedID
# Create a dataframe for UCEC miRNA-disease associations
UCEC_DD_target_df <- data.frame(
  mature_mirna_id = mature_mirna_ids,
  disease_drug = disease_drug,
  paper_pubmedID = paper_pubmedID
)
# Filter out empty entries
UCEC_DD_target_df <- UCEC_DD_target_df %>%
  filter(mature_mirna_id != "" & disease_drug != "" & paper_pubmedID != "")

# Combine all cancer-type miRNA-disease association data frames into a single dataset
all_disease_associations <- bind_rows(
  mutate(BRCA_DD_target_df, cancer_type = "BRCA"),
  mutate(OV_DD_target_df, cancer_type = "OV"),
  mutate(CESC_DD_target_df, cancer_type = "CESC"),
  mutate(UCEC_DD_target_df, cancer_type = "UCEC")
)


# Get the overlapping miRNA and 
# identify and compare diseases associated with overlapping miRNAs across different cancer types.

overlapping_miRNAs <- c("hsa.mir.150", "hsa.mir.22", "hsa.mir.30a", "hsa.mir.31", 
                        "hsa.mir.7a.1", "hsa.mir.7a.3", "hsa.mir.135b", "hsa.mir.140", 
                        "hsa.mir.142", "hsa.mir.144", "hsa.mir.148a", "hsa.mir.155", 
                        "hsa.mir.196a.1", "hsa.mir.335")

overlapping_miRNAs <- convert_miRNA_ids(overlapping_miRNAs)
overlapping_miRNAs <- gsub("mir", "miR", overlapping_miRNAs)

filtered_associations <- all_disease_associations %>%
  filter(mature_mirna_id %in% overlapping_miRNAs)


# Function to revert miRNA IDs to their original form
revert_miRNA_ids <- function(mirna_ids) {
  # Initialize an empty vector to store results
  original_ids <- character(length(mirna_ids))
  
  # Loop through each miRNA ID and revert format
  for (i in seq_along(mirna_ids)) {
    # Remove the "-5p" and "-3p" suffixes
    id <- gsub("-5p|-3p", "", mirna_ids[i])
    # Replace "-" with "." to revert to the original format
    id <- gsub("-", ".", id)
    # Store the result
    original_ids[i] <- id
  }
  
  return(original_ids)
}

# Apply the revert_miRNA_ids function to mature_mirna_id column
filtered_associations <- filtered_associations %>%
  mutate(original_mirna_id = revert_miRNA_ids(mature_mirna_id))

# Group by disease and count the number of unique miRNAs associated with each disease
common_diseases <- filtered_associations %>%
  group_by(disease_drug) %>%
  summarise(num_miRNAs = n_distinct(mature_mirna_id)) %>%
  filter(num_miRNAs > 1)  # Keep only diseases associated with more than one miRNA

# Filter the original dataset to include only these common diseases
filtered_common_diseases <- filtered_associations %>%
  filter(disease_drug %in% common_diseases$disease_drug)

# Example edges data frame
edges <- filtered_common_diseases %>%
  select(original_mirna_id, disease_drug, cancer_type) %>%
  rename(from = original_mirna_id, to = disease_drug) %>%
  group_by(from, to) %>%
  summarise(cancer_type = paste(sort(unique(cancer_type)), collapse = ", "))  # Combine cancer types into one string
# Create a graph object
graph <- graph_from_data_frame(edges, directed = FALSE)
# Filter nodes and edges
nodes_to_keep <- which(degree(graph) > 1)  # Nodes with more than one edge
edges_filtered <- edges %>%
  filter(from %in% V(graph)[nodes_to_keep]$name & to %in% V(graph)[nodes_to_keep]$name)
# Create a new graph object with filtered edges
graph_filtered <- graph_from_data_frame(edges_filtered, directed = FALSE)
# Remove isolated nodes (degree 0 or 1)
graph_filtered <- induced_subgraph(graph_filtered, which(degree(graph_filtered) > 1))
# Convert to tidygraph object
tidy_graph <- as_tbl_graph(graph_filtered)
# Assuming `tidy_graph` is already prepared and `edges_filtered` is your filtered data
p <- ggraph(tidy_graph, layout = "fr") + 
  geom_edge_link(aes(color = cancer_type), alpha = 0.6) +
  geom_node_point(aes(color = name %in% edges_filtered$from), size = 5) +
  geom_node_text(aes(label = name), repel = TRUE, size = 2, max.overlaps = 100, box.padding = 0.5) +  # Reduced size
  theme_void() +  # Remove all axes and grid lines
  labs(title = "miRNA-Disease Associations Across Cancers",
       edge_color = "Cancer Type",
       node_color = "Type") +
  scale_color_manual(values = c("blue", "red", "purple"),  # Define colors based on cancer type combinations
                     labels = c("Disease", "miRNA"),
                     name = "Association Type") +
  theme(legend.position = "bottom")  # Adjust legend position for clarity
# Export the plot to a PDF
ggsave("miRNA_Disease_Associations.pdf", plot = p, width = 10, height = 8)


# Extract gene symbols for each cancer type from their respective data frames
BRCA_genes <- BRCA_gene_target_df$target_symbol
CESC_genes <- CESC_gene_target_df$target_symbol
UCEC_genes <- UCEC_gene_target_df$target_symbol
OV_genes <- OV_gene_target_df$target_symbol
# Find the intersection of gene symbols across all four cancer types (BRCA, CESC, UCEC, OV)
geneList <- Reduce(intersect, list(BRCA_genes, CESC_genes, UCEC_genes, OV_genes))
# Convert the list of gene symbols to ENTREZ IDs using the bitr function from the clusterProfiler package
# This is necessary if you need to perform enrichment analysis based on ENTREZ IDs
geneList_entrez <- bitr(geneList, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# Perform enrichment analysis using the Disease Gene Network (DisGeNET) to identify disease-related genes
edo <- enrichDGN(geneList_entrez$ENTREZID)
# Perform pairwise similarity analysis on the enriched terms to improve the visualization of the enrichment map
edo <- pairwise_termsim(edo)
# Open a PDF device to save the generated enrichment map plot
pdf("emapplot_output.pdf", width = 8, height = 6)  # Set the dimensions of the output plot
# Generate and display the enrichment map plot
emapplot(edo)
# Close the PDF device, saving the plot to the file
dev.off()

#################################################################################
#                                     GSEA
#################################################################################
# Retrieve available Enrichr databases
dbs <- listEnrichrDbs()
# Function to filter databases based on specific criteria (KEGG 2021 and GO BP 2023)
filter_dbs <- function(name) {
  if (grepl("^KEGG", name)) {
    parts <- strsplit(name, "_")[[1]]
    if (length(parts) >= 2 && parts[2] == "2021") {
      return(TRUE)
    }
  } else if (grepl("^GO_Biological_Process_2023", name)) {  # Removed `$` to allow additional text
    return(TRUE)
  }
  return(FALSE)
}

# Apply filtering function to Enrichr databases
filtered_dbs <- dbs[sapply(dbs$libraryName, filter_dbs), ]
dbs <- filtered_dbs$libraryName

# Function to perform GSEA using EnrichR and store results in a data frame
perform_gsea_df <- function(gene_list, cancer_type, data_type) {
  if (length(gene_list) > 0) {
    enrichr_results <- enrichr(genes = gene_list, dbs)
    results_list <- lapply(enrichr_results, function(res) {
      if (nrow(res) > 0) {
        res_df <- as.data.frame(res)
        res_df$cancer_type <- cancer_type
        res_df$data_type <- data_type
        return(res_df)
      } else {
        print("No data")
        return(NULL)
      }
    })
    results_df <- do.call(rbind, results_list)
    return(results_df)
  } else {
    return(NULL)
  }
}

# Perform GSEA for each cancer type using gene targets
BRCA_GSEA_results <- perform_gsea_df(BRCA_gene_target_df$target_symbol, "BRCA", "GE")
OV_GSEA_results <- perform_gsea_df(OV_gene_target_df$target_symbol, "OV", "GE")
CESC_GSEA_results <- perform_gsea_df(CESC_gene_target_df$target_symbol, "CESC", "GE")
UCEC_GSEA_results <- perform_gsea_df(UCEC_gene_target_df$target_symbol, "UCEC", "GE")

# Extract significant GO and KEGG pathway results (Adjusted P-value ≤ 0.05)
BRCA_significant_results_GO <- BRCA_GSEA_results[grepl("GO", rownames(BRCA_GSEA_results)) & BRCA_GSEA_results$Adjusted.P.value <= 0.05, ]
BRCA_significant_results_KEGG <- BRCA_GSEA_results[grepl("KEGG", rownames(BRCA_GSEA_results)) & BRCA_GSEA_results$Adjusted.P.value <= 0.05, ]

OV_significant_results_GO <- OV_GSEA_results[grepl("GO", rownames(OV_GSEA_results)) & OV_GSEA_results$Adjusted.P.value <= 0.05, ]
OV_significant_results_KEGG <- OV_GSEA_results[grepl("KEGG", rownames(OV_GSEA_results)) & OV_GSEA_results$Adjusted.P.value <= 0.05, ]

CESC_significant_results_GO <- CESC_GSEA_results[grepl("GO", rownames(CESC_GSEA_results)) & CESC_GSEA_results$Adjusted.P.value <= 0.05, ]
CESC_significant_results_KEGG <- CESC_GSEA_results[grepl("KEGG", rownames(CESC_GSEA_results)) & CESC_GSEA_results$Adjusted.P.value <= 0.05, ]

UCEC_significant_results_GO <- UCEC_GSEA_results[grepl("GO", rownames(UCEC_GSEA_results)) & UCEC_GSEA_results$Adjusted.P.value <= 0.05, ]
UCEC_significant_results_KEGG <- UCEC_GSEA_results[grepl("KEGG", rownames(UCEC_GSEA_results)) & UCEC_GSEA_results$Adjusted.P.value <= 0.05, ]


#################################################################################
#                            GO TERMS AND KEGG PATHWAYS
#################################################################################

# Extract GO term IDs from the 'Term' column in the BRCA, OV, CESC, and UCEC data frames
# and clean up the IDs by removing parentheses
BRCA_significant_results_GO$ID <- str_extract(BRCA_significant_results_GO$Term, "\\(GO:\\d+\\)")
BRCA_significant_results_GO <- BRCA_significant_results_GO %>%
  mutate(ID = gsub("[()]", "", ID)) 
# Remove duplicate GO terms within the same cancer type (BRCA)
BRCA_significant_results_GO <- BRCA_significant_results_GO %>%
  distinct(ID, cancer_type, .keep_all = TRUE)

OV_significant_results_GO$ID <- str_extract(OV_significant_results_GO$Term, "\\(GO:\\d+\\)")
OV_significant_results_GO <- OV_significant_results_GO %>%
  mutate(ID = gsub("[()]", "", ID)) 
# Remove duplicate GO terms within the same cancer type (OV)
OV_significant_results_GO <- OV_significant_results_GO %>%
  distinct(ID, cancer_type, .keep_all = TRUE)

CESC_significant_results_GO$ID <- str_extract(CESC_significant_results_GO$Term, "\\(GO:\\d+\\)")
CESC_significant_results_GO <- CESC_significant_results_GO %>%
  mutate(ID = gsub("[()]", "", ID)) 
# Remove duplicate GO terms within the same cancer type (CESC)
CESC_significant_results_GO <- CESC_significant_results_GO %>%
  distinct(ID, cancer_type, .keep_all = TRUE)

UCEC_significant_results_GO$ID <- str_extract(UCEC_significant_results_GO$Term, "\\(GO:\\d+\\)")
UCEC_significant_results_GO <- UCEC_significant_results_GO %>%
  mutate(ID = gsub("[()]", "", ID)) 
# Remove duplicate GO terms within the same cancer type (UCEC)
UCEC_significant_results_GO <- UCEC_significant_results_GO %>%
  distinct(ID, cancer_type, .keep_all = TRUE)

# Calculate the overlap fraction (numerator/denominator) for each gene set in the GO results
BRCA_significant_results_GO$Overlaps <- sapply(strsplit(BRCA_significant_results_GO$Overlap, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
OV_significant_results_GO$Overlaps <- sapply(strsplit(OV_significant_results_GO$Overlap, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
CESC_significant_results_GO$Overlaps <- sapply(strsplit(CESC_significant_results_GO$Overlap, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
UCEC_significant_results_GO$Overlaps <- sapply(strsplit(UCEC_significant_results_GO$Overlap, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))

# Plot enrichment results for GO terms for each cancer type (BRCA, OV, CESC, UCEC)
plotEnrich(BRCA_significant_results_GO, showTerms = 10, numChar = 100, y = "Count",
           orderBy = "Combined.Score", xlab = 'GO Terms', title = 'Enrichment plot for BRCA')

plotEnrich(OV_significant_results_GO, showTerms = 10, numChar = 100, y = "Count",
           orderBy = "Combined.Score", xlab = 'GO Terms', title = 'Enrichment plot for OV')

plotEnrich(CESC_significant_results_GO, showTerms = 10, numChar = 100, y = "Count",
           orderBy = "Combined.Score", xlab = 'GO Terms', title = 'Enrichment plot for CESC')

plotEnrich(UCEC_significant_results_GO, showTerms = 10, numChar = 100, y = "Count",
           orderBy = "Combined.Score", xlab = 'GO Terms', title = 'Enrichment plot for UCEC')

# Plot enrichment results for KEGG pathways for each cancer type (BRCA, OV, CESC, UCEC)
plotEnrich(BRCA_significant_results_KEGG, showTerms = 10, numChar = 100, y = "Count",
           orderBy = "Combined.Score", xlab = 'KEGG Pathways', title = 'Enrichment plot for BRCA')

plotEnrich(OV_significant_results_KEGG, showTerms = 10, numChar = 100, y = "Count",
           orderBy = "Combined.Score", xlab = 'KEGG Pathways', title = 'Enrichment plot for OV')

plotEnrich(CESC_significant_results_KEGG, showTerms = 10, numChar = 100, y = "Count",
           orderBy = "Combined.Score", xlab = 'KEGG Pathways', title = 'Enrichment plot for CESC')

plotEnrich(UCEC_significant_results_KEGG, showTerms = 10, numChar = 100, y = "Count",
           orderBy = "Combined.Score", xlab = 'KEGG Pathways', title = 'Enrichment plot for UCEC')



#################################################################################
#                       KEGG Pathway Overlap Across All Cancers 
#################################################################################
# Extract KEGG terms from the significant results of each cancer type (BRCA, OV, CESC, UCEC)
BRCA_terms <- BRCA_significant_results_KEGG$Term
OV_terms <- OV_significant_results_KEGG$Term
CESC_terms <- CESC_significant_results_KEGG$Term
UCEC_terms <- UCEC_significant_results_KEGG$Term

# Combine all unique KEGG terms from the four cancer types
all_terms <- unique(c(BRCA_terms, OV_terms, CESC_terms, UCEC_terms))

# Create a binary matrix indicating the presence of each term in each cancer type
binary_matrix <- data.frame(
  Term = all_terms,
  BRCA = as.numeric(all_terms %in% BRCA_terms),
  OV = as.numeric(all_terms %in% OV_terms),
  CESC = as.numeric(all_terms %in% CESC_terms),
  UCEC = as.numeric(all_terms %in% UCEC_terms)
)

# Plot an upset plot to visualize the intersections of KEGG pathways across cancer types
upset(binary_matrix, sets = c("BRCA", "OV", "CESC", "UCEC"), order.by = "freq")

# Filter terms that are shared across all four cancer types
shared_terms <- binary_matrix %>%
  filter(BRCA == 1 & OV == 1 & CESC == 1 & UCEC == 1)

# Define a function to count the number of genes in a given gene string
count_genes <- function(gene_string) {
  if (is.na(gene_string)) {
    return(0)
  } else {
    genes <- strsplit(gene_string, ";")[[1]]
    return(length(genes))
  }
}

# Apply the gene count function to each cancer type's KEGG results
BRCA_significant_results_KEGG$num_genes <- sapply(BRCA_significant_results_KEGG$Genes, count_genes)
OV_significant_results_KEGG$num_genes <- sapply(OV_significant_results_KEGG$Genes, count_genes)
CESC_significant_results_KEGG$num_genes <- sapply(CESC_significant_results_KEGG$Genes, count_genes)
UCEC_significant_results_KEGG$num_genes <- sapply(UCEC_significant_results_KEGG$Genes, count_genes)

# Extract the Combined Score and the number of genes for each shared KEGG term in each cancer type
BRCA_scores <- BRCA_significant_results_KEGG %>%
  dplyr::filter(Term %in% shared_terms$Term) %>%
  dplyr::select(Term, Combined.Score, num_genes) %>%
  dplyr::rename(BRCA_score = Combined.Score, BRCA_num_genes = num_genes)

OV_scores <- OV_significant_results_KEGG %>%
  dplyr::filter(Term %in% shared_terms$Term) %>%
  dplyr::select(Term, Combined.Score, num_genes) %>%
  dplyr::rename(OV_score = Combined.Score, OV_num_genes = num_genes)

CESC_scores <- CESC_significant_results_KEGG %>%
  dplyr::filter(Term %in% shared_terms$Term) %>%
  dplyr::select(Term, Combined.Score, num_genes) %>%
  dplyr::rename(CESC_score = Combined.Score, CESC_num_genes = num_genes)

UCEC_scores <- UCEC_significant_results_KEGG %>%
  dplyr::filter(Term %in% shared_terms$Term) %>%
  dplyr::select(Term, Combined.Score, num_genes) %>%
  dplyr::rename(UCEC_score = Combined.Score, UCEC_num_genes = num_genes)

# Merge the scores and number of genes into a single data frame
shared_terms_enrichment <- shared_terms %>%
  dplyr::select(Term) %>%
  left_join(BRCA_scores, by = "Term") %>%
  left_join(OV_scores, by = "Term") %>%
  left_join(CESC_scores, by = "Term") %>%
  left_join(UCEC_scores, by = "Term")

# Reshape the data to long format for plotting the scores
shared_terms_long_score <- shared_terms_enrichment %>%
  pivot_longer(cols = -Term, names_to = "Cancer", values_to = "Score") %>%
  mutate(Cancer = gsub("_score", "", Cancer))

# Reshape the data to long format for plotting the number of genes
shared_terms_long_genes <- shared_terms_enrichment %>%
  pivot_longer(cols = -Term, names_to = "Cancer", values_to = "Genes") %>%
  mutate(Cancer = gsub("_num_genes", "", Cancer))

# Combine both data frames (scores and genes) into one
shared_terms_long <- full_join(shared_terms_long_score, shared_terms_long_genes, 
                               by = c("Term", "Cancer"))

# Filter out the rows for scores and genes separately
shared_terms_filtered <- shared_terms_long %>%
  filter(!grepl("_num_genes$", Cancer)) %>%
  filter(!grepl("_score$", Cancer))

# Aggregate the scores and the number of genes across cancers for each term
shared_terms_aggregated <- shared_terms_filtered %>%
  group_by(Term) %>%
  summarise(
    Total_Score = sum(Score, na.rm = TRUE),
    Total_Genes = sum(Genes, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  arrange(desc(Total_Score)) %>%
  slice_head(n = 10)

# Filter the shared terms based on the top 10 terms by total score
top_10_terms <- shared_terms_aggregated$Term

filtered_shared_terms <- shared_terms_filtered %>%
  filter(Term %in% top_10_terms)

# Plotting the data: Visualize the top 10 KEGG pathways with points sized by number of genes
# and colored by combined score across the four cancer types
ggplot(filtered_shared_terms, aes(x = Cancer, y = Term)) +
  geom_point(aes(size = Genes, color = Score)) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(x = "Cancer Type", y = "KEGG Pathway", size = "Num of Genes", color = "Combined Score",
       title = "Top 10 KEGG Pathway Enrichment in All Cancers") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())



#################################################################################
#                         GO TERMS Overlap Across All Cancers 
#################################################################################
BRCA_terms <- BRCA_significant_results_GO$Term
OV_terms <- OV_significant_results_GO$Term
CESC_terms <- CESC_significant_results_GO$Term
UCEC_terms <- UCEC_significant_results_GO$Term
all_terms <- unique(c(BRCA_terms, OV_terms, CESC_terms, UCEC_terms))
binary_matrix <- data.frame(
  Term = all_terms,
  BRCA = as.numeric(all_terms %in% BRCA_terms),
  OV = as.numeric(all_terms %in% OV_terms),
  CESC = as.numeric(all_terms %in% CESC_terms),
  UCEC = as.numeric(all_terms %in% UCEC_terms)
)
upset(binary_matrix, sets = c("BRCA", "OV", "CESC", "UCEC"), order.by = "freq")

shared_terms <- binary_matrix %>%
  filter(BRCA == 1 & OV == 1 & CESC == 1 & UCEC == 1)


# Define a function to count the number of genes in a string
count_genes <- function(gene_string) {
  if (is.na(gene_string)) {
    return(0)
  } else {
    genes <- strsplit(gene_string, ";")[[1]]
    return(length(genes))
  }
}

BRCA_significant_results_GO$num_genes <- sapply(BRCA_significant_results_GO$Genes, count_genes)
OV_significant_results_GO$num_genes <- sapply(OV_significant_results_GO$Genes, count_genes)
CESC_significant_results_GO$num_genes <- sapply(CESC_significant_results_GO$Genes, count_genes)
UCEC_significant_results_GO$num_genes <- sapply(UCEC_significant_results_GO$Genes, count_genes)

# Extract Scores and num_genes for each cancer type
BRCA_scores <- BRCA_significant_results_GO %>%
  dplyr::filter(Term %in% shared_terms$Term) %>%
  dplyr::select(Term, Combined.Score, num_genes) %>%
  dplyr::rename(BRCA_score = Combined.Score, BRCA_num_genes = num_genes)

OV_scores <- OV_significant_results_GO %>%
  dplyr::filter(Term %in% shared_terms$Term) %>%
  dplyr::select(Term, Combined.Score, num_genes) %>%
  dplyr::rename(OV_score = Combined.Score, OV_num_genes = num_genes)

CESC_scores <- CESC_significant_results_GO %>%
  dplyr::filter(Term %in% shared_terms$Term) %>%
  dplyr::select(Term, Combined.Score, num_genes) %>%
  dplyr::rename(CESC_score = Combined.Score, CESC_num_genes = num_genes)

UCEC_scores <- UCEC_significant_results_GO %>%
  dplyr::filter(Term %in% shared_terms$Term) %>%
  dplyr::select(Term, Combined.Score, num_genes) %>%
  dplyr::rename(UCEC_score = Combined.Score, UCEC_num_genes = num_genes)

# Merge scores into one data frame
shared_terms_enrichment <- shared_terms %>%
  dplyr::select(Term) %>%
  left_join(BRCA_scores, by = "Term") %>%
  left_join(OV_scores, by = "Term") %>%
  left_join(CESC_scores, by = "Term") %>%
  left_join(UCEC_scores, by = "Term")

# Reshape the data to long format for Score
shared_terms_long_score <- shared_terms_enrichment %>%
  pivot_longer(cols = -Term, names_to = "Cancer", values_to = "Score") %>%
  mutate(Cancer = gsub("_score", "", Cancer))

shared_terms_long_genes <- shared_terms_enrichment %>%
  pivot_longer(cols = -Term, names_to = "Cancer", values_to = "Genes") %>%
  mutate(Cancer = gsub("_num_genes", "", Cancer))

shared_terms_long <- full_join(shared_terms_long_score, shared_terms_long_genes, 
                               by = c("Term", "Cancer"))

shared_terms_filtered <- shared_terms_long %>%
  filter(!grepl("_num_genes$", Cancer)) %>%
  filter(!grepl("_score$", Cancer))

# Aggregate scores across cancers for each term
shared_terms_aggregated <- shared_terms_filtered %>%
  group_by(Term) %>%
  summarise(
    Total_Score = sum(Score, na.rm = TRUE),
    Total_Genes = sum(Genes, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  arrange(desc(Total_Score)) %>%
  slice_head(n = 10)

# Filter shared_terms_filtered based on shared_terms_aggregated
top_10_terms <- shared_terms_aggregated$Term

filtered_shared_terms <- shared_terms_filtered %>%
  filter(Term %in% top_20_terms)

# Plotting
ggplot(filtered_shared_terms, aes(x = Cancer, y = Term)) +
  geom_point(aes(size = Genes, color = Score)) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(x = "Cancer Type", y = "GO Term", size = "Num of Genes", color = "Combined Score",
       title = "Top 10 GO Term Enrichment in All Cancers") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

#################################################################################
#                           KM PLOTS
#################################################################################

BRCA_CNV_data = read.csv("BRCA/BRCA_CNV_data.csv", row.names = NULL)
CESC_CNV_data = read.csv("CESC/CESC_CNV_data.csv", row.names = NULL)
UCEC_CNV_data = read.csv("UCEC/UCEC_CNV_data.csv", row.names = NULL)

BRCA_ME_data <- read.csv("BRCA/BRCA_ME_data.csv", row.names = NULL)
CESC_ME_data <- read.csv("CESC/CESC_ME_data.csv", row.names = NULL)
UCEC_ME_data <- read.csv("UCEC/UCEC_ME_data.csv", row.names = NULL)
OV_ME_data <- read.csv("OV/OV_ME_data.csv", row.names = NULL)

BRCA_METH_data <- read.csv("BRCA/BRCA_METH_data.csv", row.names = NULL)
CESC_METH_data <- read.csv("CESC/CESC_METH_data.csv", row.names = NULL)
UCEC_METH_data <- read.csv("UCEC/UCEC_METH_data.csv", row.names = NULL)
OV_METH_data <- read.csv("OV/OV_METH_data.csv", row.names = NULL)

BRCA_CNV_data <- BRCA_CNV_data[,-1]
CESC_CNV_data <- CESC_CNV_data[,-1]
UCEC_CNV_data <- UCEC_CNV_data[,-1]

BRCA_ME_data <- BRCA_ME_data[,-1]
CESC_ME_data <- CESC_ME_data[,-1]
UCEC_ME_data <- UCEC_ME_data[,-1]
OV_ME_data <- OV_ME_data[,-1]

BRCA_METH_data <- BRCA_METH_data[,-1]
CESC_METH_data <- CESC_METH_data[,-1]
UCEC_METH_data <- UCEC_METH_data[,-1]
OV_METH_data <- OV_METH_data[,-1]


# Adapted function for creating KM plots for gene expression/miRNA data - example
gene <- "TIMP2"
data <- BRCA_gene_data
cancer <- "BRCA"
data <- data[, c("case_id", "overall_survival", "deceased", gene)]
column_value <- data[[gene]]
median_value <- median(column_value, na.rm = TRUE)
data <- data %>%
  mutate(strata = case_when(
    data[[gene]] >= median_value ~ "High",
    data[[gene]] < median_value ~ "Low",
    TRUE ~ NA_character_
  ))
combined_data <- na.omit(data)
surv_object <- Surv(time = combined_data$overall_survival, event = combined_data$deceased)
fit <- survfit(surv_object ~ strata, data = combined_data)
ggsurv <- ggsurvplot(fit,
                     data = combined_data,
                     risk.table = TRUE,
                     pval = TRUE,
                     title = paste("Kaplan-Meier Curve", cancer, "for gene:", gene),
                     xlab = "Time (days)",
                     ylab = "Survival Probability",
                     palette = "Dark2")
print(ggsurv)

# Adapted function for creating KM plots for methylation data - example
cpg_probe <-"cg17525406"
data <- UCEC_METH_data
cancer <- "UCEC"
data <- data[, c("case_id", "overall_survival", "deceased", cpg_probe)]
data[[cpg_probe]] <- as.numeric(data[[cpg_probe]])
# Categorize methylation levels
data$strata <- cut(data[[cpg_probe]], 
                   breaks = c(-Inf, 0.2, 0.8, Inf), 
                   labels = c("Low", "Partial", "High"),
                   right = FALSE)
surv_object <- Surv(time = data$overall_survival, event = data$deceased)
fit <- survfit(surv_object ~ strata, data = data)
ggsurv <- ggsurvplot(fit,
                     data = data,
                     risk.table = TRUE,
                     pval = TRUE,
                     title = paste("Kaplan-Meier Curve", cancer, "for CpG Probe:", cpg_probe),
                     xlab = "Time (days)",
                     ylab = "Survival Probability",
                     palette = "Dark2")
print(ggsurv)

# Adapted function for creating KM plots for CNV data - example
gene <-"AHNAK"
data <- BRCA_CNV_data
cancer <- "BRCA"
# Select relevant columns including the specific gene column
data <- data[, c("case_id", "overall_survival", "deceased", gene)]
# Ensure the specific gene column is numeric
data[[gene]] <- as.numeric(data[[gene]])
# Categorize CNV into Gain, Loss, and Neutral
data$strata <- ifelse(data[[gene]] > 0, "Gain",
                      ifelse(data[[gene]] < 0, "Loss", "Neutral"))
surv_object <- Surv(time = data$overall_survival, event = data$deceased)
fit <- survfit(surv_object ~ strata, data = data)
ggsurv <- ggsurvplot(fit,
                     data = data,
                     risk.table = TRUE,
                     pval = TRUE,
                     title = paste("Kaplan-Meier Curve", cancer, "for Gene:", gene),
                     xlab = "Time (days)",
                     ylab = "Survival Probability",
                     palette = "Dark2")
print(ggsurv)
