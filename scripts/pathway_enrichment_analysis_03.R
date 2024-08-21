#clean up your environment
rm(list=ls())

# For data management
install.packages('tidyverse')
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
install.packages("msigdbr")

# For visualisation
install.packages('pheatmap')
install.packages("DOSE")
install.packages("enrichplot")
install.packages("ggupset")

BiocManager::install("DOSE")
BiocManager::install("enrichplot")

# Load required libraries
library(patchwork)
library(Matrix)
library(ggplot2)
library(glmGamPoi)
library(dplyr)
library(Seurat)
library(patchwork)



# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(pheatmap)
library(clusterProfiler) # for PEA analysis
library('org.Hs.eg.db')
library(DOSE)
library(enrichplot) # for visualisations
library(ggupset) # for visualisations
library(msigdbr) # this package provides access to the MSigDB gene sets


# Load the Seurat object and assign it to a new variable name
seurat_object <- readRDS("/home/cyang40/chingyao/long_covid_project/LCvsCT/merged_LN_PC_UMAP.rds")

# Check metadata of merged seurat_objec
colnames(seurat_object@meta.data)
head(seurat_object@meta.data)
seurat_object
str(seurat_object)


## Perform DE analysis for each cluster, including positive and negative markers
seurat_object.markers <- FindAllMarkers(seurat_object)

# Save DE results in the Seurat object
seurat_object[["RNA"]]@misc$markers <- seurat_object.markers

# Export DE results to a CSV file
write.csv(seurat_object.markers, "de_markers.csv", row.names = FALSE)

# Save the Seurat object
saveRDS(seurat_object, file = "/home/cyang40/chingyao/long_covid_project/LCvsCT/merged_LN_PC_UMAP_02.rds")

# Load the Seurat object and assign it to a new variable name
seurat_object <- readRDS("/home/cyang40/chingyao/long_covid_project/LCvsCT/merged_LN_PC_UMAP_02.rds")

colnames(seurat_object@meta.data)


# Check if DE results are present in the Seurat object
if (!is.null(seurat_object[["RNA"]]@misc$markers)) {
  head(seurat_object[["RNA"]]@misc$markers)
} else {
  print("No DE markers found in the object.")
}

# Load DE results from CSV file
de_markers_path <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/de_markers.csv"
seurat_object.markers <- read.csv(de_markers_path)

# Check the structure of the loaded data
head(seurat_object.markers)

# Count the number of genes in seurat_object.markers
number_of_genes <- nrow(seurat_object.markers)

# Print the results
print(paste("Number of genes in seurat_object.markers:", number_of_genes))

# Count the number of unique genes in the DE results
unique_genes <- length(unique(seurat_object.markers$gene))

# Print the results
print(paste("Number of unique genes in seurat_object.markers:", unique_genes))

# Assume seurat_object.markers is your data frame with marker results
# Filter for significant markers (example criteria)
significant_markers <- seurat_object.markers %>%
  filter(p_val_adj <= 0.05, abs(avg_log2FC) >= 1)

# Extract gene names
gene_list <- significant_markers$gene

# Count the number of genes in gene_list
number_of_genes_in_gene_list <- length(gene_list)

# Print the results
print(paste("Number of genes in gene_list:", number_of_genes_in_gene_list))

# Convert gene symbols to Entrez IDs
gene_list_map <- clusterProfiler::bitr(
  geneID = gene_list,  # Use the list of genes from your single-cell analysis
  fromType = "SYMBOL", # Input type is gene symbols
  toType = "ENTREZID", # Convert to Entrez IDs
  OrgDb = "org.Hs.eg.db" # Use the human gene annotation database
)

head(gene_list_map)


# Extract all genes that are expressed in at least 3 cells
background_genes <- rownames(seurat_object[["RNA"]]@data)

# Count the number of genes in background_genes
number_of_genes_in_background_genes <- length(background_genes)

# Print the results
print(paste("Number of genes in background_genes:", number_of_genes_in_background_genes))


# Define the path for saving the background genes
background_genes_path <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/background_genes.txt"

# Save the list of background genes to a text file
write.table(background_genes, background_genes_path, row.names = FALSE, col.names = FALSE, quote = FALSE)

# Load the list of background genes from the text file
background_genes <- read.table(background_genes_path, stringsAsFactors = FALSE)[, 1]


# Convert gene symbols to Entrez IDs for background
background_genes_map <- clusterProfiler::bitr(
  geneID = background_genes,  # Use all expressed genes as the background
  fromType = "SYMBOL",        # Input type is gene symbols
  toType = "ENTREZID",        # Convert to Entrez IDs
  OrgDb = "org.Hs.eg.db"      # Use the human gene annotation database
)

# Check the first few mapped background genes
head(background_genes_map)


### Perform GO enrichment analysis
# Perform GO enrichment analysis for Biological Processes
ego <- enrichGO(
  gene          = gene_list_map$ENTREZID,
  universe      = background_genes_map$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",  # Biological Process
  pAdjustMethod = "BH",
  #pvalueCutoff  = 0.01,
  #qvalueCutoff  = 0.05,
   pvalueCutoff  = 0.05,
   qvalueCutoff  = 0.2,
  readable      = TRUE
)

# View the top enriched GO terms
head(ego)

# Save the 'ego' object to a file
saveRDS(ego, file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/GO/ego_results_all_02.rds")

# Load the 'ego' object from the file
ego <- readRDS("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/GO/ego_results_all_02.rds")

# View the top enriched GO terms to confirm it loaded correctly
head(ego)

# Load the enrichplot package for visualization
library(enrichplot)

# Bar plot of enriched GO terms
#barplot(ego, showCategory = 20)

# Define the file path for saving the bar plot
barplot_file_path <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/GO/barplot_ego_all_02.png"

# Export the bar plot to the specified file path
png(filename = barplot_file_path, width = 1200, height = 1200) # You can adjust width and height as needed
#barplot(ego, showCategory = 20)
print(barplot(ego, showCategory = 20) +
           ggtitle(paste("GO Enrichment for All")))
dev.off()


# Dot plot of enriched GO terms
#dotplot(ego)

# Define the file path for saving the dot plot
dotplot_file_path <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/GO/dotplot_ego_all_02.png"

# Export the dot plot to the specified file path
png(filename = dotplot_file_path, width = 1200, height = 1200) # Adjust size if needed
#dotplot(ego, showCategory = 20)
print(dotplot(ego, showCategory = 20) +
           ggtitle(paste("GO Enrichment for All")))
dev.off()


### Perform enrichment analysis for each cluster

# Load DE results from CSV file
de_markers_path <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/de_markers.csv"
seurat_object.markers <- read.csv(de_markers_path)

# Check the structure of the loaded data
head(seurat_object.markers)

# Filter for significant markers
significant_markers <- seurat_object.markers %>%
  filter(p_val_adj <= 0.05, abs(avg_log2FC) >= 1)

head(significant_markers)

# Extract gene names
gene_list <- significant_markers$gene

# Count the number of genes in gene_list
number_of_genes_in_gene_list <- length(gene_list)

# Print the results
print(paste("Number of genes in gene_list:", number_of_genes_in_gene_list))

# Define the path for saving the background genes
background_genes_path <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/background_genes.txt"

# Load the list of background genes from the text file
background_genes <- read.table(background_genes_path, stringsAsFactors = FALSE)[, 1]

# Convert gene symbols to Entrez IDs for background
background_genes_map <- clusterProfiler::bitr(
  geneID = background_genes,  # Use all expressed genes as the background
  fromType = "SYMBOL",        # Input type is gene symbols
  toType = "ENTREZID",        # Convert to Entrez IDs
  OrgDb = "org.Hs.eg.db"      # Use the human gene annotation database
)

# Check the first few mapped background genes
head(background_genes_map)


# Prepare list to store pathway results
pathway_results <- list()

# Perform enrichment analysis for each cluster
for (cluster in unique(significant_markers$cluster)) {
  # Get the DE genes for the current cluster
  cluster_genes <- significant_markers$gene[significant_markers$cluster == cluster]
  
  # Convert gene symbols to Entrez IDs for cluster genes
  entrez_cluster_genes <- clusterProfiler::bitr(
    geneID = cluster_genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Hs.eg.db"
  )

  # Perform GO enrichment analysis using the pre-processed background
  go_enrichment <- enrichGO(
    gene = entrez_cluster_genes$ENTREZID,
    universe = background_genes_map$ENTREZID,
    OrgDb = org.Hs.eg.db,
    keyType = 'ENTREZID',
    ont = "BP",  # Biological Process
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
  
  # Store results in a list
  pathway_results[[as.character(cluster)]] <- go_enrichment
}

# Print pathway results for each cluster
for (cluster in names(pathway_results)) {
  cat("Cluster:", cluster, "\n")
  print(head(pathway_results[[cluster]]@result))
  cat("\n")
}


# Check the results
pathway_results


# Save pathway results to a file
saveRDS(pathway_results, file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/GO/GO_pathway_results_each_cluster_02.rds")

# Load pathway results from the file for review
pathway_results <- readRDS("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/GO/GO_pathway_results_each_cluster_02.rds")

# Check the results
pathway_results

# Define the directory for saving plots
plot_dir <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/GO/GO_pathway_plots_02"
dir.create(plot_dir, showWarnings = FALSE)

# Visualize results for all clusters
for (cluster in names(pathway_results)) {
  print(paste("Generating plots for Cluster", cluster))
  
  # Define file names for the plots
  barplot_file <- file.path(plot_dir, paste0("GO_Enrichment_Cluster_", cluster, "_barplot_02.png"))
  dotplot_file <- file.path(plot_dir, paste0("GO_Enrichment_Cluster_", cluster, "_dotplot_02.png"))
  
  # Generate and save the barplot for the current cluster
  png(filename = barplot_file, width = 800, height = 600)
  #print(barplot(pathway_results[[cluster]], showCategory = 10, main = paste("GO Enrichment for Cluster", cluster)))
  print(barplot(pathway_results[[cluster]], showCategory = 10) +
           ggtitle(paste("GO Enrichment for Cluster", cluster)))
  dev.off()
  
  # Generate and save the dotplot for the current cluster
  png(filename = dotplot_file, width = 800, height = 600)
  print(dotplot(pathway_results[[cluster]], showCategory = 10) + 
          ggtitle(paste("GO Enrichment for Cluster", cluster)))
  dev.off()
  
 }


# Define the directory for saving plots
plot_dir <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/GO/GO_pathway_plots_02"
dir.create(plot_dir, showWarnings = FALSE)

# Visualize results for all clusters, excluding Cluster 9
for (cluster in names(pathway_results)) {
  # Check if the cluster is 9 and skip it
  if (cluster == "9") {
    print(paste("Skipping Cluster", cluster, "due to zero enriched terms."))
    next
  }
  
  # Proceed with plotting if the cluster is not 9
  print(paste("Generating plots for Cluster", cluster))
  
  # Define file names for the plots
  barplot_file <- file.path(plot_dir, paste0("GO_Enrichment_Cluster_", cluster, "_barplot_02.png"))
  dotplot_file <- file.path(plot_dir, paste0("GO_Enrichment_Cluster_", cluster, "_dotplot_02.png"))
  
  # Generate and save the barplot for the current cluster
  png(filename = barplot_file, width = 800, height = 600)
  #print(barplot(pathway_results[[cluster]], showCategory = 10, main = paste("GO Enrichment for Cluster", cluster)))
  print(barplot(pathway_results[[cluster]], showCategory = 10) +
           ggtitle(paste("GO Enrichment for Cluster", cluster)))
  dev.off()
  
  # Generate and save the dotplot for the current cluster
  png(filename = dotplot_file, width = 800, height = 600)
  print(dotplot(pathway_results[[cluster]], showCategory = 10) + 
          ggtitle(paste("GO Enrichment for Cluster", cluster)))
  dev.off()
}


### using the Hallmark category from the MSigDB database

# Load DE results from CSV file
de_markers_path <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/de_markers.csv"
seurat_object.markers <- read.csv(de_markers_path)

# Check the structure of the loaded data
head(seurat_object.markers)

# Assume seurat_object.markers is your data frame with marker results
# Filter for significant markers (example criteria)
significant_markers <- seurat_object.markers %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 1)

# Extract gene names
gene_list <- significant_markers$gene

# Count the number of genes in gene_list
number_of_genes_in_gene_list <- length(gene_list)

# Print the results
print(paste("Number of genes in gene_list:", number_of_genes_in_gene_list))

# Convert gene symbols to Entrez IDs
gene_list_map <- clusterProfiler::bitr(
  geneID = gene_list,  # Use the list of genes from your single-cell analysis
  fromType = "SYMBOL", # Input type is gene symbols
  toType = "ENTREZID", # Convert to Entrez IDs
  OrgDb = "org.Hs.eg.db" # Use the human gene annotation database
)

head(gene_list_map)

# Define the path for saving the background genes
background_genes_path <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/background_genes.txt"

# Load the list of background genes from the text file
background_genes <- read.table(background_genes_path, stringsAsFactors = FALSE)[, 1]

# Convert gene symbols to Entrez IDs for background
background_genes_map <- clusterProfiler::bitr(
  geneID = background_genes,  # Use all expressed genes as the background
  fromType = "SYMBOL",        # Input type is gene symbols
  toType = "ENTREZID",        # Convert to Entrez IDs
  OrgDb = "org.Hs.eg.db"      # Use the human gene annotation database
)

# Check the first few mapped background genes
head(background_genes_map)

# Retrieve MSigDB Gene Sets for Homo Sapiens
m_df <- msigdbr(species = "Homo sapiens")
head(m_df)

# Retrieve hallmark gene sets for Homo sapiens
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

# View the first few entries to verify
head(m_t2g)

# Perform enrichment analysis using hallmark gene sets with additional parameters
em <- enricher(
  gene          = gene_list_map$ENTREZID,
  TERM2GENE     = m_t2g,
  universe      = background_genes_map$ENTREZID,
  pAdjustMethod = "BH",         # Benjamini-Hochberg adjustment
  pvalueCutoff  = 0.05,         # P-value cutoff for significance
  qvalueCutoff  = 0.1,         # Q-value cutoff for significance
  )


summary(em)

# View the top enriched hallmark gene sets
head(em)

# Save the 'em' object to a file
saveRDS(em, file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/hallmark/em_results_all_02.rds")

# Load the 'em' object from the file
em <- readRDS("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/hallmark/em_results_all_02.rds")

# View the top enriched hallmark terms to confirm it loaded correctly
head(em)


# Visualize the results
# Bar plot of enriched hallmark gene sets
#barplot(em, showCategory = 20)

# Define the file path for saving the bar plot
barplot_file_path <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/hallmark/barplot_em_all_02.png"

# Export the bar plot to the specified file path
png(filename = barplot_file_path, width = 1200, height = 1200) # You can adjust width and height as needed
#barplot(em, showCategory = 20)
print(barplot(em, showCategory = 20) +
           ggtitle(paste("Hallmark Enrichment for All")))
dev.off()

# Dot plot of enriched hallmark gene sets
#dotplot(em)

# Define the file path for saving the dot plot
dotplot_file_path <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/hallmark/dotplot_em_all_02.png"

# Export the dot plot to the specified file path
png(filename = dotplot_file_path, width = 1200, height = 1200) # Adjust size if needed
#dotplot(em, showCategory = 20)
print(dotplot(em, showCategory = 20) +
           ggtitle(paste("Hallmark Enrichment for All")))
dev.off()


### Perform enrichment analysis for each cluster

# Load DE results from CSV file
de_markers_path <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/de_markers.csv"
seurat_object.markers <- read.csv(de_markers_path)

# Check the structure of the loaded data
head(seurat_object.markers)

# Filter for significant markers
significant_markers <- seurat_object.markers %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 1)

head(significant_markers)

# Extract gene names
gene_list <- significant_markers$gene

# Count the number of genes in gene_list
number_of_genes_in_gene_list <- length(gene_list)

# Print the results
print(paste("Number of genes in gene_list:", number_of_genes_in_gene_list))

## using the Hallmark category from the MSigDB database
# Retrieve MSigDB Gene Sets for Homo Sapiens
m_df <- msigdbr(species = "Homo sapiens")
head(m_df)

# Retrieve hallmark gene sets for Homo sapiens
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

# View the first few entries to verify
head(m_t2g)

# Define the path for saving the background genes
background_genes_path <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/background_genes.txt"

# Load the list of background genes from the text file
background_genes <- read.table(background_genes_path, stringsAsFactors = FALSE)[, 1]

# Convert gene symbols to Entrez IDs for background
background_genes_map <- clusterProfiler::bitr(
  geneID = background_genes,  # Use all expressed genes as the background
  fromType = "SYMBOL",        # Input type is gene symbols
  toType = "ENTREZID",        # Convert to Entrez IDs
  OrgDb = "org.Hs.eg.db"      # Use the human gene annotation database
)

# Check the first few mapped background genes
head(background_genes_map)


# Prepare list to store pathway results
pathway_results <- list()

# Perform enrichment analysis for each cluster
for (cluster in unique(significant_markers$cluster)) {
  # Get the DE genes for the current cluster
  cluster_genes <- significant_markers$gene[significant_markers$cluster == cluster]
  
  # Convert gene symbols to Entrez IDs for cluster genes
  entrez_cluster_genes <- clusterProfiler::bitr(
    geneID = cluster_genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Hs.eg.db"
  )

  
  # Perform enrichment analysis using hallmark gene sets with additional parameters
em <- enricher(
  gene          = entrez_cluster_genes$ENTREZID,
  TERM2GENE     = m_t2g,
  universe      = background_genes_map$ENTREZID,
  pAdjustMethod = "BH",         # Benjamini-Hochberg adjustment
  pvalueCutoff  = 0.05,         # P-value cutoff for significance
  qvalueCutoff  = 0.1,         # Q-value cutoff for significance
  )
  
  # Store results in a list
  pathway_results[[as.character(cluster)]] <- em
}

# Print pathway results for each cluster
for (cluster in names(pathway_results)) {
  cat("Cluster:", cluster, "\n")
  print(head(pathway_results[[cluster]]@result))
  cat("\n")
}

# Check the results
pathway_results


# Save pathway results to a file
saveRDS(pathway_results, file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/hallmark/hallmark_pathway_results_each_cluster_02.rds")

# Load pathway results from the file for review
pathway_results <- readRDS("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/hallmark/hallmark_pathway_results_each_cluster_02.rds")

# Check the results
pathway_results


# Define the directory for saving plots
plot_dir <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/hallmark/hallmark_pathway_plots_02"
dir.create(plot_dir, showWarnings = FALSE)

# Visualize results for all clusters
for (cluster in names(pathway_results)) {
  print(paste("Generating plots for Cluster", cluster))
  
  # Define file names for the plots
  barplot_file <- file.path(plot_dir, paste0("hallmark_Enrichment_Cluster_", cluster, "_barplot_02.png"))
  dotplot_file <- file.path(plot_dir, paste0("hallmark_Enrichment_Cluster_", cluster, "_dotplot_02.png"))
  
  # Generate and save the barplot for the current cluster
  png(filename = barplot_file, width = 800, height = 600)
  #print(barplot(pathway_results[[cluster]], showCategory = 10, main = paste("Hallmark Enrichment for Cluster", cluster)))
  print(barplot(pathway_results[[cluster]], showCategory = 10) +
           ggtitle(paste("Hallmark Enrichment for Cluster", cluster)))
  dev.off()
  
  # Generate and save the dotplot for the current cluster
  png(filename = dotplot_file, width = 800, height = 600)
  print(dotplot(pathway_results[[cluster]], showCategory = 10) + 
          ggtitle(paste("Hallmark Enrichment for Cluster", cluster)))
  dev.off()
  
 }



### using the C2:CP (curated pathways) category from the MSigDB database

# Assume seurat_object.markers is your data frame with marker results
# Load DE results from CSV file
de_markers_path <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/de_markers.csv"
seurat_object.markers <- read.csv(de_markers_path)

# Check the structure of the loaded data
head(seurat_object.markers)

# Filter for significant markers (example criteria)
significant_markers <- seurat_object.markers

# Extract gene names
gene_list <- significant_markers$gene

# Count the number of genes in gene_list
number_of_genes_in_gene_list <- length(gene_list)

# Print the results
print(paste("Number of genes in gene_list:", number_of_genes_in_gene_list))

# Convert gene symbols to Entrez IDs
gene_list_map <- clusterProfiler::bitr(
  geneID = gene_list,  # Use the list of genes from your single-cell analysis
  fromType = "SYMBOL", # Input type is gene symbols
  toType = "ENTREZID", # Convert to Entrez IDs
  OrgDb = "org.Hs.eg.db" # Use the human gene annotation database
)

head(gene_list_map)


# Retrieve C2:CP (curated pathways) gene sets for Homo sapiens
c2_cp_t2g <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP") %>% 
  dplyr::select(gs_name, entrez_gene)

# View the first few entries to verify
head(c2_cp_t2g)

# Check for common IDs between your gene list and the TERM2GENE data frame
common_ids <- intersect(gene_list_map$ENTREZID, c2_cp_t2g$entrez_gene)
length(common_ids)  # This should be greater than zero for potential enrichment results

## Examine distribution of gene set sizes
table(dplyr::count(c2_cp_t2g, gs_name))

# First, create a summary of counts per gene set
gene_set_counts <- dplyr::count(c2_cp_t2g, gs_name)

# Display the table of counts
table(gene_set_counts$n)

# Or you can view the number of genes in each set directly
print(gene_set_counts, n=29)


# Define the path for saving the background genes
background_genes_path <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/background_genes.txt"

# Load the list of background genes from the text file
background_genes <- read.table(background_genes_path, stringsAsFactors = FALSE)[, 1]

# Convert gene symbols to Entrez IDs for background
background_genes_map <- clusterProfiler::bitr(
  geneID = background_genes,  # Use all expressed genes as the background
  fromType = "SYMBOL",        # Input type is gene symbols
  toType = "ENTREZID",        # Convert to Entrez IDs
  OrgDb = "org.Hs.eg.db"      # Use the human gene annotation database
)

# Check the first few mapped background genes
head(background_genes_map)


# Perform enrichment analysis using C2:CP gene sets
em_c2_cp <- enricher(
  gene          = gene_list_map$ENTREZID,
  TERM2GENE     = c2_cp_t2g,
  universe      = background_genes_map$ENTREZID,
  #pAdjustMethod = "BH",         # Benjamini-Hochberg adjustment
  #pvalueCutoff  = 0.01,         # P-value cutoff for significance
  #qvalueCutoff  = 0.05,         # Q-value cutoff for significance
  )

summary(em_c2_cp)

# View the top enriched c2_cp gene sets
head(em_c2_cp)

# Save the 'em_c2_cp' object to a file
saveRDS(em_c2_cp, file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/c2_cp/em_c2_cp_results_all_02.rds")

# Load the 'em_c2_cp' object from the file
em_c2_cp <- readRDS("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/c2_cp/em_c2_cp_results_all_02.rds")

# View the top enriched C2:CP gene sets
head(em_c2_cp)

# Visualize the results
# Bar plot of enriched C2:CP gene sets
#barplot(em_c2_cp, showCategory = 20)

# Define the file path for saving the bar plot
barplot_file_path <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/c2_cp/barplot_c2_cp_all_02.png"

# Export the bar plot to the specified file path
png(filename = barplot_file_path, width = 1200, height = 1200) # You can adjust width and height as needed
#barplot(em_c2_cp, showCategory = 20)
print(barplot(em_c2_cp, showCategory = 20) +
           ggtitle(paste("C2:CP Enrichment for All")))
dev.off()

# Dot plot of enriched C2:CP gene sets
#dotplot(em_c2_cp)

# Define the file path for saving the dot plot
dotplot_file_path <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/c2_cp/dotplot_c2_cp_all_02.png"

# Export the dot plot to the specified file path
png(filename = dotplot_file_path, width = 1200, height = 1200) # Adjust size if needed
#dotplot(em_c2_cp, showCategory = 20)
print(dotplot(em_c2_cp, showCategory = 20) +
           ggtitle(paste("C2:CP Enrichment for All")))
dev.off()


### Perform enrichment analysis for each cluster

# Load DE results from CSV file
de_markers_path <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/de_markers.csv"
seurat_object.markers <- read.csv(de_markers_path)

# Check the structure of the loaded data
head(seurat_object.markers)

# Filter for significant markers (example criteria)
significant_markers <- seurat_object.markers

# Extract gene names
gene_list <- significant_markers$gene

# Count the number of genes in gene_list
number_of_genes_in_gene_list <- length(gene_list)

# Print the results
print(paste("Number of genes in gene_list:", number_of_genes_in_gene_list))

# using the Hallmark category from the MSigDB database
# Retrieve MSigDB Gene Sets for Homo Sapiens
m_df <- msigdbr(species = "Homo sapiens")
head(m_df)

# Retrieve C2:CP (curated pathways) gene sets for Homo sapiens
c2_cp_t2g <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP") %>% 
  dplyr::select(gs_name, entrez_gene)

# View the first few entries to verify
head(c2_cp_t2g)


# Define the path for saving the background genes
background_genes_path <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/background_genes.txt"

# Load the list of background genes from the text file
background_genes <- read.table(background_genes_path, stringsAsFactors = FALSE)[, 1]

# Convert gene symbols to Entrez IDs for background
background_genes_map <- clusterProfiler::bitr(
  geneID = background_genes,  # Use all expressed genes as the background
  fromType = "SYMBOL",        # Input type is gene symbols
  toType = "ENTREZID",        # Convert to Entrez IDs
  OrgDb = "org.Hs.eg.db"      # Use the human gene annotation database
)

# Check the first few mapped background genes
head(background_genes_map)


# Prepare list to store pathway results
pathway_results <- list()

# Perform enrichment analysis for each cluster
for (cluster in unique(significant_markers$cluster)) {
  # Get the DE genes for the current cluster
  cluster_genes <- significant_markers$gene[significant_markers$cluster == cluster]
  
  # Convert gene symbols to Entrez IDs for cluster genes
  entrez_cluster_genes <- clusterProfiler::bitr(
    geneID = cluster_genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Hs.eg.db"
  )

  
  # Perform enrichment analysis using hallmark gene sets with additional parameters
em_c2_cp <- enricher(
  gene          = entrez_cluster_genes$ENTREZID,
  TERM2GENE     = c2_cp_t2g,
  universe      = background_genes_map$ENTREZID,
  #pAdjustMethod = "BH",         # Benjamini-Hochberg adjustment
  #pvalueCutoff  = 0.05,         # P-value cutoff for significance
  #qvalueCutoff  = 0.1,         # Q-value cutoff for significance
  )
  
  # Store results in a list
  pathway_results[[as.character(cluster)]] <- em_c2_cp
}

# Print pathway results for each cluster
for (cluster in names(pathway_results)) {
  cat("Cluster:", cluster, "\n")
  print(head(pathway_results[[cluster]]@result))
  cat("\n")
}

# Check the results
pathway_results


# Save pathway results to a file
saveRDS(pathway_results, file = "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/c2_cp/c2_cp_pathway_results_each_cluster_02.rds")

# Load pathway results from the file for review
pathway_results <- readRDS("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/c2_cp/c2_cp_pathway_results_each_cluster_02.rds")

# Check the results
pathway_results


# Define the directory for saving plots
plot_dir <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/c2_cp/c2_cp_pathway_plots_02"
dir.create(plot_dir, showWarnings = FALSE)

# Define the clusters to skip
clusters_to_skip <- c("9", "13", "14")

# Visualize results for all clusters, excluding specified clusters
for (cluster in names(pathway_results)) {
  # Check if the cluster is in the list to skip
  if (cluster %in% clusters_to_skip) {
    print(paste("Skipping Cluster", cluster, "due to zero enriched terms."))
    next
  }
  
  # Proceed with plotting if the cluster is not in the skip list
  print(paste("Generating plots for Cluster", cluster))
  
  # Define file names for the plots
  barplot_file <- file.path(plot_dir, paste0("c2_cp_Enrichment_Cluster_", cluster, "_barplot_02.png"))
  dotplot_file <- file.path(plot_dir, paste0("c2_cp_Enrichment_Cluster_", cluster, "_dotplot_02.png"))
  
  # Generate and save the barplot for the current cluster
  png(filename = barplot_file, width = 800, height = 600)
  #print(barplot(pathway_results[[cluster]], showCategory = 10, main = paste("C2_CP Enrichment for Cluster", cluster)))
  print(barplot(pathway_results[[cluster]], showCategory = 10) +
           ggtitle(paste("C2_CP Enrichment for Cluster", cluster)))
  dev.off()
  
  # Generate and save the dotplot for the current cluster
  png(filename = dotplot_file, width = 800, height = 600)
  print(dotplot(pathway_results[[cluster]], showCategory = 10) + 
          ggtitle(paste("C2_CP Enrichment for Cluster", cluster)))
  dev.off()
}


# Define the directory for saving plots
plot_dir <- "/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/c2_cp/c2_cp_pathway_plots_02"
dir.create(plot_dir, showWarnings = FALSE)

# Visualize results for all clusters
for (cluster in names(pathway_results)) {
  print(paste("Generating plots for Cluster", cluster))
  
  # Define file names for the plots
  barplot_file <- file.path(plot_dir, paste0("c2_cp_Enrichment_Cluster_", cluster, "_barplot_02.png"))
  dotplot_file <- file.path(plot_dir, paste0("c2_cp_Enrichment_Cluster_", cluster, "_dotplot_02.png"))
  
  # Generate and save the barplot for the current cluster
  png(filename = barplot_file, width = 800, height = 600)
  #print(barplot(pathway_results[[cluster]], showCategory = 10, main = paste("Hallmark Enrichment for Cluster", cluster)))
  print(barplot(pathway_results[[cluster]], showCategory = 10) +
           ggtitle(paste("C2_CP Enrichment for Cluster", cluster)))
  dev.off()
  
  # Generate and save the dotplot for the current cluster
  png(filename = dotplot_file, width = 800, height = 600)
  print(dotplot(pathway_results[[cluster]], showCategory = 10) + 
          ggtitle(paste("C2_CP Enrichment for Cluster", cluster)))
  dev.off()
  
 }

