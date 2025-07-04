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

.libPaths("/path/to/Rlib")
setwd("/path/to/pwd")

# Load multi-omics feature-level fusion datasets for each cancer type
BRCA_data = read.csv("BRCA/1S/BRCA_ME_GE_selected_features.csv", row.names = NULL)
OV_data = read.csv("OV/2S/OV_ME_DM_selected_features.csv", row.names = NULL)
CESC_data = read.csv("CESC/2S/CESC_ME_DM_selected_features.csv", row.names = NULL)
UCEC_data = read.csv("UCEC/2S/UCEC_ME_GE_CNV_DM_selected_features.csv", row.names = NULL)


# Load gene expression datasets for each cancer type
BRCA_gene_data = read.csv("BRCA/BRCA_GE_clean.csv", row.names = NULL)
OV_gene_data = read.csv("OV/OV_GE_clean.csv", row.names = NULL)
CESC_gene_data = read.csv("CESC/CESC_GE_clean.csv", row.names = NULL)
UCEC_gene_data = read.csv("UCEC/UCEC_GE_clean.csv", row.names = NULL)

# Remove the sampleID from gene expression datasets 
BRCA_gene_data <- BRCA_gene_data[, -1]
OV_gene_data <- OV_gene_data[, -1]
CESC_gene_data <- CESC_gene_data[,-1]
UCEC_gene_data <- UCEC_gene_data[, -1]

# Remove rows with NA in the overall_survival column
BRCA_gene_data  <- BRCA_gene_data[!is.na(BRCA_gene_data[[2]]), ]
OV_gene_data  <- OV_gene_data[!is.na(OV_gene_data[[2]]), ]
CESC_gene_data  <- CESC_gene_data[!is.na(CESC_gene_data[[2]]), ]
UCEC_gene_data  <- UCEC_gene_data [!is.na(UCEC_gene_data [[2]]), ]

# Extract the 'Feature' column
BRCA_vars <- BRCA_data$Feature
OV_vars <- OV_data$Feature
CESC_vars <- CESC_data$Feature
UCEC_vars <- UCEC_data$Feature

# Generate and save an UpSet plot to visualize the overlap of features between cancer types
listInput <- list(BRCA = BRCA_vars, OV = OV_vars, CESC = CESC_vars, UCEC = UCEC_vars)
png("upset_plot.png", width = 1600, height = 1200, res = 300)
upset(fromList(listInput), order.by = "freq")
grid.text("Pan-cancer Signature Overlap", x = 0.7, y = 0.98, gp = gpar(fontsize = 10, fontface = "bold"))
dev.off()
dev.new()

# Get all combinations of cancer types (2 or more)
cancer_combinations <- unlist(lapply(2:length(listInput), function(k) {
  combn(names(listInput), k, simplify = FALSE)
}), recursive = FALSE)

# Initialize empty list to hold results
overlap_table <- list()

# Find overlaps for each combination
for (combo in cancer_combinations) {
  # Get the intersection of features across the selected cancer types
  overlapping_features <- Reduce(intersect, listInput[combo])
  
  # If there's any overlap, store the result
  if (length(overlapping_features) > 0) {
    combo_name <- paste(combo, collapse = ", ")
    overlap_table[[combo_name]] <- overlapping_features
  }
}

# Convert to data frame for better viewing
overlap_df <- do.call(rbind, lapply(names(overlap_table), function(name) {
  data.frame(
    `Cancer Types` = name,
    `Overlapping Features` = paste(overlap_table[[name]], collapse = ", "),
    stringsAsFactors = FALSE
  )
}))

write.csv(overlap_df, "pan_cancer_feature_overlap_table.csv", row.names = FALSE)


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
# Function to extract gene names for DM features
extract_dm_names <- function(vars) {
  gene_names <- sub("^DM_", "", vars[grepl("^DM_", vars)])
  return(gene_names)
}
# Function to extract gene names for GE features
extract_ge_names <- function(vars) {
  gene_names <- sub("^GE_", "", vars[grepl("^GE_", vars)])
  return(gene_names)
}


# Extract GE and CNV gene names for each cancer type
BRCA_ge_genes <- extract_ge_names(BRCA_vars)
UCEC_ge_genes <- extract_ge_names(UCEC_vars)

UCEC_cnv_genes <- extract_cnv_names(UCEC_vars)
BRCA_cnv_genes <- extract_cnv_names(BRCA_vars)

OV_me_genes <- extract_me_names(OV_vars)
CESC_me_genes <- extract_me_names(CESC_vars)
UCEC_me_genes <- extract_me_names(UCEC_vars)
BRCA_me_genes <- extract_me_names(BRCA_vars)

OV_dm_genes <- extract_dm_names(OV_vars)
CESC_dm_genes <- extract_dm_names(CESC_vars)
#################################################################################
#                             MiRA Targets
#################################################################################

# Function to convert miRNA IDs to hsa-miR-XX-5p or hsa-miR-XX-3p format
convert_miRNA_ids <- function(mirna_ids) {
  mirna_ids <- gsub("\\.", "-", mirna_ids)
  return(mirna_ids)
}
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
expressed_genes <- colnames(BRCA_gene_data)[3:ncol(BRCA_gene_data)]
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
expressed_genes <- colnames(OV_gene_data)[3:ncol(OV_gene_data)]
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
expressed_genes <- colnames(CESC_gene_data)[3:ncol(CESC_gene_data)]
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
expressed_genes <- colnames(UCEC_gene_data)[3:ncol(UCEC_gene_data)]
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
# Save the plot to a PDF file
pdf("pan-cancer_gene_target_overlaps.pdf", width = 8, height = 6)
upset(upset_data, order.by = "freq")
dev.off()


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

# Combine into a named list
miRNA_lists <- list(
  BRCA = BRCA_me_genes,
  OV = OV_me_genes,
  CESC = CESC_me_genes,
  UCEC = UCEC_me_genes
)
# Flatten into one long vector with cancer labels
miRNA_all <- unlist(miRNA_lists)
# Count how many times each miRNA appears across cancers
miRNA_counts <- table(miRNA_all)
# Get miRNAs that appear in 2 or more cancers
miRNA_overlap <- names(miRNA_counts[miRNA_counts >= 2])

overlapping_miRNAs<- convert_miRNA_ids(miRNA_overlap)
filtered_associations <- all_disease_associations %>%
  filter(mature_mirna_id %in% overlapping_miRNAs)

# Group by disease and count the number of unique miRNAs associated with each disease
common_diseases <- filtered_associations %>%
  group_by(disease_drug) %>%
  summarise(num_miRNAs = n_distinct(mature_mirna_id))


# Filter the original dataset to include only these common diseases
filtered_common_diseases <- filtered_associations %>%
  filter(disease_drug %in% common_diseases$disease_drug)

# Example edges data frame
edges <- filtered_common_diseases %>%
  select(mature_mirna_id, disease_drug, cancer_type) %>%
  rename(from = mature_mirna_id, to = disease_drug) %>%
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
# Plot
p <- ggraph(tidy_graph, layout = "fr") + 
  geom_edge_link(aes(color = cancer_type), alpha = 0.6, show.legend = TRUE) +
  geom_node_point(aes(color = name %in% edges_filtered$from), size = 5, show.legend = TRUE) +
  geom_node_text(aes(label = name), repel = TRUE, size = 4, max.overlaps = 100, box.padding = 0.5) +
  theme_void() +
  guides(
    edge_color = guide_legend(
      title = "Cancer Type",
      override.aes = list(size = 1, linetype = 1, shape = NA)  # Removes dot, keeps line
    ),
    color = guide_legend(
      title = "Node Type",
      override.aes = list(size = 4, shape = 16)  # Keeps circle for node legend
    )
  ) +
  scale_color_manual(
    values = c("FALSE" = "blue", "TRUE" = "red"),
    labels = c("Disease", "miRNA"),
    name = "Node Type"
  ) +
  theme(
    legend.justification = c("right", "top"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    plot.title = element_blank()
  )
# Save the plot
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
pdf("emapplot_output_miRNA.pdf", width = 10, height = 8)  # Set the dimensions of the output plot
# Generate and display the enrichment map plot
emapplot(
  edo,
  layout = "nicely",          # More spaced layout than default
  showCategory = 30,          # Max number of categories shown (adjust as needed)
)
# Close the PDF device, saving the plot to the file
dev.off()



#################################################################################
#                                   GSEA
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
BRCA_GSEA_results <- perform_gsea_df(BRCA_gene_target_df$target_symbol, "BRCA", "ME")
OV_GSEA_results <- perform_gsea_df(OV_gene_target_df$target_symbol, "OV", "ME")
CESC_GSEA_results <- perform_gsea_df(CESC_gene_target_df$target_symbol, "CESC", "ME")
UCEC_GSEA_results <- perform_gsea_df(UCEC_gene_target_df$target_symbol, "UCEC", "ME")

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

save_plotEnrich_pdf <- function(enrich_result, filename, xlab_title, plot_title) {
  pdf(filename, width = 8, height = 6)     # 1. open PDF device
  p <- plotEnrich(enrich_result, showTerms = 10, numChar = 100, y = "Count",
                  orderBy = "Combined.Score", xlab = xlab_title, title = plot_title)
  print(p)                                # 2. explicitly print the plot if it's a ggplot object
  dev.off()                               # 3. close the device to save the file
}

# GO term enrichment plots
save_plotEnrich_pdf(BRCA_significant_results_GO, "miRNA_BRCA_GO_enrichment.pdf", "GO Terms", "Enrichment plot for BRCA")
save_plotEnrich_pdf(OV_significant_results_GO, "miRNA_OV_GO_enrichment.pdf", "GO Terms", "Enrichment plot for OV")
save_plotEnrich_pdf(CESC_significant_results_GO, "miRNA_CESC_GO_enrichment.pdf", "GO Terms", "Enrichment plot for CESC")
save_plotEnrich_pdf(UCEC_significant_results_GO, "miRNA_UCEC_GO_enrichment.pdf", "GO Terms", "Enrichment plot for UCEC")

# KEGG pathway enrichment plots
save_plotEnrich_pdf(BRCA_significant_results_KEGG, "miRNA_BRCA_KEGG_enrichment.pdf", "KEGG Pathways", "Enrichment plot for BRCA")
save_plotEnrich_pdf(OV_significant_results_KEGG, "miRNA_OV_KEGG_enrichment.pdf", "KEGG Pathways", "Enrichment plot for OV")
save_plotEnrich_pdf(CESC_significant_results_KEGG, "miRNA_CESC_KEGG_enrichment.pdf", "KEGG Pathways", "Enrichment plot for CESC")
save_plotEnrich_pdf(UCEC_significant_results_KEGG, "miRNA_UCEC_KEGG_enrichment.pdf", "KEGG Pathways", "Enrichment plot for UCEC")


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

# Open PDF device to save plot
pdf("KEGG_pathways_upset_plot.pdf", width = 8, height = 6)
# Create the UpSet plot
upset(binary_matrix, sets = c("BRCA", "OV", "CESC", "UCEC"), order.by = "freq")
# Close the PDF device to finalize the file
dev.off()

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
p <- ggplot(filtered_shared_terms, aes(x = Cancer, y = Term)) +
  geom_point(aes(size = Genes, color = Score)) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(x = "Cancer Type", y = "KEGG Pathway", size = "Num of Genes", color = "Combined Score",
       title = "Top 10 KEGG Pathway Enrichment in All Cancers") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave("KEGG_Pathway_Enrichment_All_Cancers.pdf", plot = p, width = 8, height = 6)

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
# Open PDF device to save plot
pdf("GO_Terms_upset_plot.pdf", width = 8, height = 6)
# Create the UpSet plot
upset(binary_matrix, sets = c("BRCA", "OV", "CESC", "UCEC"), order.by = "freq")
# Close the PDF device to finalize the file
dev.off()
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
  filter(Term %in% top_10_terms)

# Plotting
p <-ggplot(filtered_shared_terms, aes(x = Cancer, y = Term)) +
  geom_point(aes(size = Genes, color = Score)) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(x = "Cancer Type", y = "GO Term", size = "Num of Genes", color = "Combined Score",
       title = "Top 10 GO Term Enrichment in All Cancers") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggsave("GO_Term_Enrichment_All_Cancers.pdf", plot = p, width = 8, height = 6)

