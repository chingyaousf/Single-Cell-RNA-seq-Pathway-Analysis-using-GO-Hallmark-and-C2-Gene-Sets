
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


go_results

# Define clusters associated with the NK/T-cell group
nk_t_cell_clusters <- c("0", "1", "2", "4", "9", "12")


# Extract GO pathway results for each NK/T-cell cluster
go_nk_t_cell_results <- setNames(
  lapply(nk_t_cell_clusters, function(cl) go_results[[cl]]), 
  nk_t_cell_clusters
)

go_nk_t_cell_results


# Check the structure of the cluster 0 result
str(go_nk_t_cell_results[["0"]])

# View the first few rows of the data
head(go_nk_t_cell_results[["0"]])

# Check the names of the columns (if applicable)
names(go_nk_t_cell_results[["0"]])


# Define a function to extract pathways from a cluster result
extract_pathways <- function(cluster_result) {
  cluster_result@result$Description
}

# Extract GO pathway descriptions for each NK/T-cell cluster
go_nk_t_cell_pathways <- lapply(go_nk_t_cell_results, extract_pathways)

# The list already has names corresponding to the cluster numbers, no need to rename it
# View the first few extracted pathways for each cluster
head(go_nk_t_cell_pathways)


# Show the first few pathways for each of the first few clusters
concise_output <- lapply(go_nk_t_cell_pathways, function(paths) head(paths, 5))
concise_output


# View the first few clusters and their first few pathways
concise_output <- head(go_nk_t_cell_pathways, 2)  # Adjust the number '2' to show more or fewer clusters
concise_output <- lapply(concise_output, function(paths) head(paths, 5))
concise_output

# Identify common pathways across all NK/T-cell clusters
common_pathways <- Reduce(intersect, go_nk_t_cell_pathways)

head(common_pathways)


# Identify unique pathways for each cluster
unique_pathways <- lapply(go_nk_t_cell_pathways, function(pathways) {
  setdiff(pathways, common_pathways)
})

head(unique_pathways, n=10)

# Limit the output to the first 2 clusters and first 5 unique pathways within each cluster
concise_unique_pathways <- head(unique_pathways, n = 3)
concise_unique_pathways <- lapply(concise_unique_pathways, function(paths) head(paths, 5))

concise_unique_pathways



# Define a function to extract the top N pathways from a cluster result
extract_top_pathways <- function(cluster_result, n = 10) {
  top_pathways <- cluster_result@result$Description[1:n]
  return(top_pathways)
}

# Extract the top 10 or 20 pathways for each NK/T-cell cluster
top_10_pathways <- lapply(go_nk_t_cell_results, extract_top_pathways, n = 10)
top_20_pathways <- lapply(go_nk_t_cell_results, extract_top_pathways, n = 20)

# View the top 10 pathways for each cluster
top_10_pathways


# Identify common pathways across all NK/T-cell clusters for the top 10 pathways
common_top_10_pathways <- Reduce(intersect, top_10_pathways)

common_top_10_pathways

# Identify unique pathways for each cluster for the top 10 pathways
unique_top_10_pathways <- lapply(top_10_pathways, function(pathways) {
  setdiff(pathways, common_top_10_pathways)
})

unique_top_10_pathways 


# Similarly for the top 20 pathways
common_top_20_pathways <- Reduce(intersect, top_20_pathways)
unique_top_20_pathways <- lapply(top_20_pathways, function(pathways) {
  setdiff(pathways, common_top_20_pathways)
})


common_top_20_pathways
unique_top_20_pathways



# Load required libraries
library(ggplot2)
library(reshape2)

# Prepare data for bar chart
top_10_pathways_df <- data.frame(
  Cluster = rep(names(top_10_pathways), each = 10),
  Pathway = unlist(top_10_pathways),
  Rank = rep(1:10, times = length(top_10_pathways))
)

# Plot a heapmap for the top 10 pathways for each cluster
ggplot(top_10_pathways_df, aes(x = factor(Rank), y = Cluster, fill = Pathway)) +
  geom_tile(color = "white") +
  scale_fill_viridis_d(option = "C") +
  labs(x = "Rank", y = "Cluster", fill = "Pathway") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top 10 Pathways Across Clusters")

# Save the plot as a file in the specified path
ggsave("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/GO/Top_10_Pathways_heapmap_nk_t_cell.png", width = 10, height = 8)


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
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  axis.text.y = element_text(size = 15) ) +
  ggtitle("Top 10 Pathways Across Clusters") +
  coord_flip()  # Flip coordinates for better readability

# Save the plot
ggsave("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/GO/top_10_pathways_barplot_nk_t_cell.png", width = 12, height = 8)



go_results <- readRDS("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/GO/GO_pathway_results_each_cluster_02.rds")
go_results

# Define clusters associated with the Mono/Mac group
mono_mac_clusters <- c("3", "5", "7")

# Extract GO pathway results for each Mono/Mac cluster
go_mono_mac_results <- setNames(
  lapply(mono_mac_clusters, function(cl) go_results[[cl]]), 
  mono_mac_clusters
)

# Check the structure of the cluster 3 result
str(go_mono_mac_results[["3"]])

# View the first few rows of the data
head(go_mono_mac_results[["3"]])

# Check the names of the columns (if applicable)
names(go_mono_mac_results[["3"]])

# Define a function to extract pathways from a cluster result
extract_pathways <- function(cluster_result) {
  cluster_result@result$Description
}

# Extract GO pathway descriptions for each Mono/Mac cluster
go_mono_mac_pathways <- lapply(go_mono_mac_results, extract_pathways)

# View the first few extracted pathways for each cluster
head(go_mono_mac_pathways)

# Show the first few pathways for each of the first few clusters
concise_output <- lapply(go_mono_mac_pathways, function(paths) head(paths, 5))
concise_output

# Identify common pathways across all Mono/Mac clusters
common_pathways <- Reduce(intersect, go_mono_mac_pathways)
head(common_pathways)

# Identify unique pathways for each cluster
unique_pathways <- lapply(go_mono_mac_pathways, function(pathways) {
  setdiff(pathways, common_pathways)
})
head(unique_pathways, n=10)

# Limit the output to the first 2 clusters and first 5 unique pathways within each cluster
concise_unique_pathways <- head(unique_pathways, n = 3)
concise_unique_pathways <- lapply(concise_unique_pathways, function(paths) head(paths, 5))
concise_unique_pathways

# Define a function to extract the top N pathways from a cluster result
extract_top_pathways <- function(cluster_result, n = 10) {
  top_pathways <- cluster_result@result$Description[1:n]
  return(top_pathways)
}

# Extract the top 10 or 20 pathways for each Mono/Mac cluster
top_10_pathways <- lapply(go_mono_mac_results, extract_top_pathways, n = 10)
top_20_pathways <- lapply(go_mono_mac_results, extract_top_pathways, n = 20)

# Identify common pathways across all Mono/Mac clusters for the top 10 pathways
common_top_10_pathways <- Reduce(intersect, top_10_pathways)
common_top_10_pathways

# Identify unique pathways for each cluster for the top 10 pathways
unique_top_10_pathways <- lapply(top_10_pathways, function(pathways) {
  setdiff(pathways, common_top_10_pathways)
})
unique_top_10_pathways 

# Similarly for the top 20 pathways
common_top_20_pathways <- Reduce(intersect, top_20_pathways)
unique_top_20_pathways <- lapply(top_20_pathways, function(pathways) {
  setdiff(pathways, common_top_20_pathways)
})

common_top_20_pathways
unique_top_20_pathways

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
  ggtitle("Top GO 10 Pathways Across Mono/Mac Clusters")

# Save the plot as a file in the specified path
ggsave("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/GO/Top_10_GO_Pathways_heatmap_mono_mac.png", width = 10, height = 8)

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
  ggtitle("Top 10 GO Pathways Across Mono/Mac Clusters") +
  coord_flip()  # Flip coordinates for better readability

# Save the plot
ggsave("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/GO/top_10_GO_pathways_barplot_mono_mac.png", width = 12, height = 8)



# Load the GO pathway results
go_results <- readRDS("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/GO/GO_pathway_results_each_cluster_02.rds")

# Define clusters associated with the B-cell group
b_cell_clusters <- c("6", "8", "11")

# Extract GO pathway results for each B-cell cluster
go_b_cell_results <- setNames(
  lapply(b_cell_clusters, function(cl) go_results[[cl]]), 
  b_cell_clusters
)

# Check the structure of the cluster 6 result
str(go_b_cell_results[["6"]])

# View the first few rows of the data
head(go_b_cell_results[["6"]])

# Check the names of the columns (if applicable)
names(go_b_cell_results[["6"]])

# Define a function to extract pathways from a cluster result
extract_pathways <- function(cluster_result) {
  cluster_result@result$Description
}

# Extract GO pathway descriptions for each B-cell cluster
go_b_cell_pathways <- lapply(go_b_cell_results, extract_pathways)

# View the first few extracted pathways for each cluster
head(go_b_cell_pathways)

# Show the first few pathways for each of the first few clusters
concise_output <- lapply(go_b_cell_pathways, function(paths) head(paths, 5))
concise_output

# Identify common pathways across all B-cell clusters
common_pathways <- Reduce(intersect, go_b_cell_pathways)
head(common_pathways)

# Identify unique pathways for each cluster
unique_pathways <- lapply(go_b_cell_pathways, function(pathways) {
  setdiff(pathways, common_pathways)
})
head(unique_pathways, n=10)

# Limit the output to the first 2 clusters and first 5 unique pathways within each cluster
concise_unique_pathways <- head(unique_pathways, n = 3)
concise_unique_pathways <- lapply(concise_unique_pathways, function(paths) head(paths, 5))
concise_unique_pathways

# Define a function to extract the top N pathways from a cluster result
extract_top_pathways <- function(cluster_result, n = 10) {
  top_pathways <- cluster_result@result$Description[1:n]
  return(top_pathways)
}

# Extract the top 10 or 20 pathways for each B-cell cluster
top_10_pathways <- lapply(go_b_cell_results, extract_top_pathways, n = 10)
top_20_pathways <- lapply(go_b_cell_results, extract_top_pathways, n = 20)

# Identify common pathways across all B-cell clusters for the top 10 pathways
common_top_10_pathways <- Reduce(intersect, top_10_pathways)
common_top_10_pathways

# Identify unique pathways for each cluster for the top 10 pathways
unique_top_10_pathways <- lapply(top_10_pathways, function(pathways) {
  setdiff(pathways, common_top_10_pathways)
})
unique_top_10_pathways 

# Similarly for the top 20 pathways
common_top_20_pathways <- Reduce(intersect, top_20_pathways)
unique_top_20_pathways <- lapply(top_20_pathways, function(pathways) {
  setdiff(pathways, common_top_20_pathways)
})

common_top_20_pathways
unique_top_20_pathways

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
  ggtitle("Top 10 GO Pathways Across B-cell Clusters")

# Save the plot as a file in the specified path
ggsave("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/GO/Top_10_GO_Pathways_heatmap_b_cell.png", width = 10, height = 8)

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
  ggtitle("Top 10 GO Pathways Across B-cell Clusters") +
  coord_flip()  # Flip coordinates for better readability

# Save the plot
ggsave("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/GO/top_10_GO_pathways_barplot_b_cell.png", width = 12, height = 8)
