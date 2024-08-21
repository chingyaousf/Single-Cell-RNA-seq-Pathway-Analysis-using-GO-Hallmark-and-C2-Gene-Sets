# Load required libraries
library(patchwork)
library(Matrix)
library(ggplot2)
library(SingleR)
library(celldex)
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


go_results <- readRDS("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/GO/GO_pathway_results_each_cluster_02.rds")
hallmark_results <- readRDS("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/hallmark/hallmark_pathway_results_each_cluster_02.rds")
c2_cp_results <- readRDS("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/c2_cp/c2_cp_pathway_results_each_cluster_02.rds")


c2_cp_results 

# Define clusters associated with the NK/T-cell group
nk_t_cell_clusters <- c("0", "1", "2", "4", "9", "12")

# Extract c2_cp pathway results for each NK/T-cell cluster
c2_cp_nk_t_cell_results <- setNames(
  lapply(nk_t_cell_clusters, function(cl) c2_cp_results[[cl]]), 
  nk_t_cell_clusters
)

# Check the structure of the cluster 0 result
str(c2_cp_nk_t_cell_results[["0"]])

# View the first few rows of the data
head(c2_cp_nk_t_cell_results[["0"]])

# Check the names of the columns (if applicable)
names(c2_cp_nk_t_cell_results[["0"]])

# Define a function to extract pathways from a cluster result
extract_pathways <- function(cluster_result) {
  cluster_result@result$Description
}

# Extract c2_cp pathway descriptions for each NK/T-cell cluster
c2_cp_nk_t_cell_pathways <- lapply(c2_cp_nk_t_cell_results, extract_pathways)

# View the first few extracted pathways for each cluster
head(c2_cp_nk_t_cell_pathways)

# Show the first few pathways for each of the first few clusters
concise_output <- lapply(c2_cp_nk_t_cell_pathways, function(paths) head(paths, 5))
concise_output

# Identify common pathways across all NK/T-cell clusters
common_pathways <- Reduce(intersect, c2_cp_nk_t_cell_pathways)
head(common_pathways)

# Identify unique pathways for each cluster
unique_pathways <- lapply(c2_cp_nk_t_cell_pathways, function(pathways) {
  setdiff(pathways, common_pathways)
})
head(unique_pathways, n=10)

# Define a function to extract the top N pathways from a cluster result
extract_top_pathways <- function(cluster_result, n = 10) {
  top_pathways <- cluster_result@result$Description[1:n]
  return(top_pathways)
}

# Extract the top 10 or 20 pathways for each NK/T-cell cluster
top_10_pathways <- lapply(c2_cp_nk_t_cell_results, extract_top_pathways, n = 10)
top_20_pathways <- lapply(c2_cp_nk_t_cell_results, extract_top_pathways, n = 20)

# Identify common pathways across all NK/T-cell clusters for the top 10 pathways
common_top_10_pathways <- Reduce(intersect, top_10_pathways)
common_top_10_pathways

# Identify unique pathways for each cluster for the top 10 pathways
unique_top_10_pathways <- lapply(top_10_pathways, function(pathways) {
  setdiff(pathways, common_top_10_pathways)
})
unique_top_10_pathways 

# Prepare data for heatmap
top_10_pathways_df <- data.frame(
  Cluster = rep(names(top_10_pathways), each = 10),
  Pathway = unlist(top_10_pathways),
  Rank = rep(1:10, times = length(top_10_pathways))
)

# Plot a heatmap for the top 10 pathways for each cluster
ggplot(top_10_pathways_df, aes(x = factor(Rank), y = Cluster, fill = Pathway)) +
  geom_tile(color = "white") +
  scale_fill_viridis_d(option = "C") +
  labs(x = "Rank", y = "Cluster", fill = "Pathway") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top 10 C2:CP Pathways Across Clusters")

# Save the plot as a file in the specified path
ggsave("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/c2_cp/Top_10_c2_cp_Pathways_heatmap_nk_t_cell.png", width = 10, height = 8)

# Prepare data for bar chart
top_10_pathways_df <- data.frame(
  Cluster = factor(rep(names(top_10_pathways), each = 10)),
  Pathway = unlist(top_10_pathways),
  Rank = rep(1:10, times = length(top_10_pathways))
)

# Plot a bar chart for the top 10 pathways for each cluster
ggplot(top_10_pathways_df, aes(x = Pathway, y = Rank, fill = Cluster)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis_d(option = "C") +
  labs(x = "Pathway", y = "Rank", fill = "Cluster") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top 10 C2:CP Pathways Across Clusters") +
  coord_flip()  # Flip coordinates for better readability

# Save the plot
ggsave("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/c2_cp/top_10_c2_cp_pathways_barplot_nk_t_cell.png", width = 12, height = 8)



## mono_mac_clusters
# Load the C2:CP pathway results
c2_cp_results <- readRDS("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/c2_cp/c2_cp_pathway_results_each_cluster_02.rds")

# Define clusters associated with the Mono/Mac group
mono_mac_clusters <- c("3", "5", "7")

# Extract C2:CP pathway results for each Mono/Mac cluster
c2_cp_mono_mac_results <- setNames(
  lapply(mono_mac_clusters, function(cl) c2_cp_results[[cl]]), 
  mono_mac_clusters
)

# Check the structure of the cluster 3 result
str(c2_cp_mono_mac_results[["3"]])

# View the first few rows of the data
head(c2_cp_mono_mac_results[["3"]])

# Check the names of the columns (if applicable)
names(c2_cp_mono_mac_results[["3"]])

# Define a function to extract pathways from a cluster result
extract_pathways <- function(cluster_result) {
  cluster_result@result$Description
}

# Extract C2:CP pathway descriptions for each Mono/Mac cluster
c2_cp_mono_mac_pathways <- lapply(c2_cp_mono_mac_results, extract_pathways)

# View the first few extracted pathways for each cluster
head(c2_cp_mono_mac_pathways)

# Show the first few pathways for each of the first few clusters
concise_output <- lapply(c2_cp_mono_mac_pathways, function(paths) head(paths, 5))
concise_output

# Identify common pathways across all Mono/Mac clusters
common_pathways <- Reduce(intersect, c2_cp_mono_mac_pathways)
head(common_pathways)

# Identify unique pathways for each cluster
unique_pathways <- lapply(c2_cp_mono_mac_pathways, function(pathways) {
  setdiff(pathways, common_pathways)
})
head(unique_pathways, n=10)

# Define a function to extract the top N pathways from a cluster result
extract_top_pathways <- function(cluster_result, n = 10) {
  top_pathways <- cluster_result@result$Description[1:n]
  return(top_pathways)
}

# Extract the top 10 or 20 pathways for each Mono/Mac cluster
top_10_pathways <- lapply(c2_cp_mono_mac_results, extract_top_pathways, n = 10)
top_20_pathways <- lapply(c2_cp_mono_mac_results, extract_top_pathways, n = 20)

# Identify common pathways across all Mono/Mac clusters for the top 10 pathways
common_top_10_pathways <- Reduce(intersect, top_10_pathways)
common_top_10_pathways

# Identify unique pathways for each cluster for the top 10 pathways
unique_top_10_pathways <- lapply(top_10_pathways, function(pathways) {
  setdiff(pathways, common_top_10_pathways)
})
unique_top_10_pathways 

# Prepare data for heatmap
top_10_pathways_df <- data.frame(
  Cluster = rep(names(top_10_pathways), each = 10),
  Pathway = unlist(top_10_pathways),
  Rank = rep(1:10, times = length(top_10_pathways))
)

# Plot a heatmap for the top 10 pathways for each cluster
ggplot(top_10_pathways_df, aes(x = factor(Rank), y = Cluster, fill = Pathway)) +
  geom_tile(color = "white") +
  scale_fill_viridis_d(option = "F") +
  labs(x = "Rank", y = "Cluster", fill = "Pathway") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top 10 C2:CP Pathways Across Mono/Mac Clusters")

# Save the plot as a file in the specified path
ggsave("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/c2_cp/Top_10_c2_cp_Pathways_heatmap_mono_mac.png", width = 10, height = 8)

# Prepare data for bar chart
top_10_pathways_df <- data.frame(
  Cluster = factor(rep(names(top_10_pathways), each = 10)),
  Pathway = unlist(top_10_pathways),
  Rank = rep(1:10, times = length(top_10_pathways))
)

# Plot a bar chart for the top 10 pathways for each cluster
ggplot(top_10_pathways_df, aes(x = Pathway, y = Rank, fill = Cluster)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis_d(option = "F") +
  labs(x = "Pathway", y = "Rank", fill = "Cluster") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top 10 C2:CP Pathways Across Mono/Mac Clusters") +
  coord_flip()  # Flip coordinates for better readability

# Save the plot
ggsave("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/c2_cp/top_10_c2_cp_pathways_barplot_mono_mac.png", width = 12, height = 8)




# Load the C2:CP pathway results
c2_cp_results <- readRDS("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/c2_cp/c2_cp_pathway_results_each_cluster_02.rds")

# Define clusters associated with the B-cell group
b_cell_clusters <- c("6", "8", "11")

# Extract C2:CP pathway results for each B-cell cluster
c2_cp_b_cell_results <- setNames(
  lapply(b_cell_clusters, function(cl) c2_cp_results[[cl]]), 
  b_cell_clusters
)

# Check the structure of the cluster 6 result
str(c2_cp_b_cell_results[["6"]])

# View the first few rows of the data
head(c2_cp_b_cell_results[["6"]])

# Check the names of the columns (if applicable)
names(c2_cp_b_cell_results[["6"]])

# Define a function to extract pathways from a cluster result
extract_pathways <- function(cluster_result) {
  cluster_result@result$Description
}

# Extract C2:CP pathway descriptions for each B-cell cluster
c2_cp_b_cell_pathways <- lapply(c2_cp_b_cell_results, extract_pathways)

# View the first few extracted pathways for each cluster
head(c2_cp_b_cell_pathways)

# Show the first few pathways for each of the first few clusters
concise_output <- lapply(c2_cp_b_cell_pathways, function(paths) head(paths, 5))
concise_output

# Identify common pathways across all B-cell clusters
common_pathways <- Reduce(intersect, c2_cp_b_cell_pathways)
head(common_pathways)

# Identify unique pathways for each cluster
unique_pathways <- lapply(c2_cp_b_cell_pathways, function(pathways) {
  setdiff(pathways, common_pathways)
})
head(unique_pathways, n=10)

# Define a function to extract the top N pathways from a cluster result
extract_top_pathways <- function(cluster_result, n = 10) {
  top_pathways <- cluster_result@result$Description[1:n]
  return(top_pathways)
}

# Extract the top 10 or 20 pathways for each B-cell cluster
top_10_pathways <- lapply(c2_cp_b_cell_results, extract_top_pathways, n = 10)
top_20_pathways <- lapply(c2_cp_b_cell_results, extract_top_pathways, n = 20)

# Identify common pathways across all B-cell clusters for the top 10 pathways
common_top_10_pathways <- Reduce(intersect, top_10_pathways)
common_top_10_pathways

# Identify unique pathways for each cluster for the top 10 pathways
unique_top_10_pathways <- lapply(top_10_pathways, function(pathways) {
  setdiff(pathways, common_top_10_pathways)
})
unique_top_10_pathways 

# Prepare data for heatmap
top_10_pathways_df <- data.frame(
  Cluster = rep(names(top_10_pathways), each = 10),
  Pathway = unlist(top_10_pathways),
  Rank = rep(1:10, times = length(top_10_pathways))
)

# Plot a heatmap for the top 10 pathways for each cluster
ggplot(top_10_pathways_df, aes(x = factor(Rank), y = Cluster, fill = Pathway)) +
  geom_tile(color = "white") +
  scale_fill_viridis_d(option = "C") +
  labs(x = "Rank", y = "Cluster", fill = "Pathway") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top 10 C2:CP Pathways Across B-cell Clusters")

# Save the plot as a file in the specified path
ggsave("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/c2_cp/Top_10_c2_cp_Pathways_heatmap_b_cell.png", width = 10, height = 8)

# Prepare data for bar chart
top_10_pathways_df <- data.frame(
  Cluster = factor(rep(names(top_10_pathways), each = 10)),
  Pathway = unlist(top_10_pathways),
  Rank = rep(1:10, times = length(top_10_pathways))
)

# Plot a bar chart for the top 10 pathways for each cluster
ggplot(top_10_pathways_df, aes(x = Pathway, y = Rank, fill = Cluster)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis_d(option = "C") +
  labs(x = "Pathway", y = "Rank", fill = "Cluster") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top 10 C2:CP Pathways Across B-cell Clusters") +
  coord_flip()  # Flip coordinates for better readability

# Save the plot
ggsave("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/c2_cp/top_10_c2_cp_pathways_barplot_b_cell.png", width = 12, height = 8)






