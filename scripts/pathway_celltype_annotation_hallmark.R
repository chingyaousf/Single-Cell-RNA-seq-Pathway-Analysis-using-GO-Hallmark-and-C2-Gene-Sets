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


hallmark_results

# Define clusters associated with the NK/T-cell group
nk_t_cell_clusters <- c("0", "1", "2", "4", "9", "12")


# Extract hallmark pathway results for each NK/T-cell cluster
hallmark_nk_t_cell_results <- setNames(
  lapply(nk_t_cell_clusters, function(cl) hallmark_results[[cl]]), 
  nk_t_cell_clusters
)

hallmark_nk_t_cell_results


# Check the structure of the cluster 0 result
str(hallmark_nk_t_cell_results[["0"]])

# View the first few rows of the data
head(hallmark_nk_t_cell_results[["0"]])

# Check the names of the columns (if applicable)
names(hallmark_nk_t_cell_results[["0"]])


# Define a function to extract pathways from a cluster result
extract_pathways <- function(cluster_result) {
  cluster_result@result$Description
}

# Extract hallmark pathway descriptions for each NK/T-cell cluster
hallmark_nk_t_cell_pathways <- lapply(hallmark_nk_t_cell_results, extract_pathways)

# The list already has names corresponding to the cluster numbers, no need to rename it
# View the first few extracted pathways for each cluster
head(hallmark_nk_t_cell_pathways)


# Show the first few pathways for each of the first few clusters
concise_output <- lapply(hallmark_nk_t_cell_pathways, function(paths) head(paths, 10))
concise_output


# View the first few clusters and their first few pathways
concise_output <- head(hallmark_nk_t_cell_pathways, 2)  # Adjust the number '2' to show more or fewer clusters
concise_output <- lapply(concise_output, function(paths) head(paths, 5))
concise_output

# Identify common pathways across all NK/T-cell clusters
common_pathways <- Reduce(intersect, hallmark_nk_t_cell_pathways)

head(common_pathways)


# Identify unique pathways for each cluster
unique_pathways <- lapply(hallmark_nk_t_cell_pathways, function(pathways) {
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
top_10_pathways <- lapply(hallmark_nk_t_cell_results, extract_top_pathways, n = 10)
top_20_pathways <- lapply(hallmark_nk_t_cell_results, extract_top_pathways, n = 20)

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

# Plot a bar chart for the top 10 pathways for each cluster
ggplot(top_10_pathways_df, aes(x = factor(Rank), y = Cluster, fill = Pathway)) +
  geom_tile(color = "white") +
  scale_fill_viridis_d(option = "H") +
  labs(x = "Rank", y = "Cluster", fill = "Pathway") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top 10 hallmark Pathways Across Clusters")

# Save the plot as a file in the specified path
ggsave("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/hallmark/Top_10_hallmark_Pathways_heapmap_nk_t_cell.png", width = 10, height = 8)


# Plot the heatmap with improved contrast and legend positioning
ggplot(top_10_pathways_df, aes(x = factor(Rank), y = Cluster, fill = Pathway)) +
  geom_tile(color = "white") +
  scale_fill_viridis_d(option = "B", direction = 1) +  # Improved contrast with a different color option
  labs(x = "Rank", y = "Cluster", fill = "Pathway") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",  # Move legend to the right side
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) +
  ggtitle("Top 10 Hallmark Pathways Across Clusters")

# Save the plot as a file in the specified path
ggsave("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/hallmark/Top_10_hallmark_Pathways_heatmap_nk_t_cell_improved.png", width = 10, height = 8)


# Plot the heatmap with a different color palette
ggplot(top_10_pathways_df, aes(x = factor(Rank), y = Cluster, fill = Pathway)) +
  geom_tile(color = "white") +
  scale_fill_brewer(palette = "Set3") +  # Using a qualitative palette for distinct colors
  labs(x = "Rank", y = "Cluster", fill = "Pathway") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",  # Legend on the right side
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) +
  ggtitle("Top 10 Hallmark Pathways Across Clusters")

# Save the plot as a file in the specified path
ggsave("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/hallmark/Top_10_hallmark_Pathways_heatmap_nk_t_cell_brewer.png", width = 10, height = 8)



# Plot the heatmap with a Viridis color palette that supports a large number of colors
ggplot(top_10_pathways_df, aes(x = factor(Rank), y = Cluster, fill = Pathway)) +
  geom_tile(color = "white") +
  scale_fill_viridis_d(option = "C") +  # Using the 'C' option from viridis, which supports more colors
  labs(x = "Rank", y = "Cluster", fill = "Pathway") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) +
  ggtitle("Top 10 Hallmark Pathways Across Clusters")

# Save the plot as a file in the specified path
ggsave("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/hallmark/Top_10_hallmark_Pathways_heatmap_nk_t_cell_viridis.png", width = 10, height = 8)





# Generate a custom color palette with enough colors
custom_colors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(unique(top_10_pathways_df$Pathway)))

# Plot the heatmap with the custom color palette
ggplot(top_10_pathways_df, aes(x = factor(Rank), y = Cluster, fill = Pathway)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = custom_colors) +
  labs(x = "Rank", y = "Cluster", fill = "Pathway") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) +
  ggtitle("Top 10 Hallmark Pathways Across Clusters")

# Save the plot as a file in the specified path
ggsave("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/hallmark/Top_10_hallmark_Pathways_heatmap_nk_t_cell_custom.png", width = 10, height = 8)


# Load required library
library(RColorBrewer)

# Generate a custom color palette with enough colors for your data
custom_colors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(unique(top_10_pathways_df$Pathway)))

# Plot the heatmap with the custom color palette
ggplot(top_10_pathways_df, aes(x = factor(Rank), y = Cluster, fill = Pathway)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = custom_colors) +
  labs(x = "Rank", y = "Cluster", fill = "Pathway") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) +
  ggtitle("Top 10 Hallmark Pathways Across Clusters")

# Save the plot as a file in the specified path
ggsave("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/hallmark/Top_10_hallmark_Pathways_heatmap_nk_t_cell_custom.png", width = 10, height = 8)


# Prepare data for bar chart
top_10_pathways_df <- data.frame(
  Cluster = factor(rep(names(top_10_pathways), each = 10)),
  Pathway = unlist(top_10_pathways),
  Rank = rep(1:10, times = length(top_10_pathways))
)

# Plot a bar chart for the top 10 pathways for each cluster
ggplot(top_10_pathways_df, aes(x = Pathway, y = Rank, fill = Cluster)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis_d(option = "E") +
  labs(x = "Pathway", y = "Rank", fill = "Cluster") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top 10 hallmark Pathways Across Clusters") +
  coord_flip()  # Flip coordinates for better readability

# Save the plot
ggsave("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/hallmark/top_10_hallmark_pathways_barplot_nk_t_cell.png", width = 12, height = 8)





# Load the hallmark pathway results
hallmark_results <- readRDS("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/hallmark/hallmark_pathway_results_each_cluster_02.rds")

hallmark_results

# Define clusters associated with the Mono/Mac group
mono_mac_clusters <- c("3", "5", "7")

# Extract hallmark pathway results for each Mono/Mac cluster
hallmark_mono_mac_results <- setNames(
  lapply(mono_mac_clusters, function(cl) hallmark_results[[cl]]), 
  mono_mac_clusters
)

# Check the structure of the cluster 3 result
str(hallmark_mono_mac_results[["3"]])

# View the first few rows of the data
head(hallmark_mono_mac_results[["3"]])

# Check the names of the columns (if applicable)
names(hallmark_mono_mac_results[["3"]])

# Define a function to extract pathways from a cluster result
extract_pathways <- function(cluster_result) {
  cluster_result@result$Description
}

# Extract hallmark pathway descriptions for each Mono/Mac cluster
hallmark_mono_mac_pathways <- lapply(hallmark_mono_mac_results, extract_pathways)

# View the first few extracted pathways for each cluster
head(hallmark_mono_mac_pathways)

# Show the first few pathways for each of the first few clusters
concise_output <- lapply(hallmark_mono_mac_pathways, function(paths) head(paths, 10))
concise_output

# View the first few clusters and their first few pathways
concise_output <- head(hallmark_mono_mac_pathways, 2)  # Adjust the number '2' to show more or fewer clusters
concise_output <- lapply(concise_output, function(paths) head(paths, 5))
concise_output

# Identify common pathways across all Mono/Mac clusters
common_pathways <- Reduce(intersect, hallmark_mono_mac_pathways)
head(common_pathways)

# Identify unique pathways for each cluster
unique_pathways <- lapply(hallmark_mono_mac_pathways, function(pathways) {
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
top_10_pathways <- lapply(hallmark_mono_mac_results, extract_top_pathways, n = 10)
top_20_pathways <- lapply(hallmark_mono_mac_results, extract_top_pathways, n = 20)

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

# Prepare data for bar chart
top_10_pathways_df <- data.frame(
  Cluster = rep(names(top_10_pathways), each = 10),
  Pathway = unlist(top_10_pathways),
  Rank = rep(1:10, times = length(top_10_pathways))
)

# Plot a bar chart for the top 10 pathways for each cluster
ggplot(top_10_pathways_df, aes(x = factor(Rank), y = Cluster, fill = Pathway)) +
  geom_tile(color = "white") +
  scale_fill_viridis_d(option = "H") +
  labs(x = "Rank", y = "Cluster", fill = "Pathway") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top 10 Hallmark Pathways Across Mono/Mac Clusters")

# Save the plot as a file in the specified path
ggsave("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/hallmark/Top_10_hallmark_Pathways_heatmap_mono_mac.png", width = 10, height = 8)

# Prepare data for bar chart
top_10_pathways_df <- data.frame(
  Cluster = factor(rep(names(top_10_pathways), each = 10)),
  Pathway = unlist(top_10_pathways),
  Rank = rep(1:10, times = length(top_10_pathways))
)

# Plot a bar chart for the top 10 pathways for each cluster
ggplot(top_10_pathways_df, aes(x = Pathway, y = Rank, fill = Cluster)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis_d(option = "E") +
  labs(x = "Pathway", y = "Rank", fill = "Cluster") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top 10 Hallmark Pathways Across Mono/Mac Clusters") +
  coord_flip()  # Flip coordinates for better readability

# Save the plot
ggsave("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/hallmark/top_10_hallmark_pathways_barplot_mono_mac.png", width = 12, height = 8)



## b_cell_clusters
# Load the hallmark pathway results
hallmark_results <- readRDS("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/hallmark/hallmark_pathway_results_each_cluster_02.rds")

# Define clusters associated with the B-cell group
b_cell_clusters <- c("6", "8", "11")

# Extract hallmark pathway results for each B-cell cluster
hallmark_b_cell_results <- setNames(
  lapply(b_cell_clusters, function(cl) hallmark_results[[cl]]), 
  b_cell_clusters
)

# Check the structure of the cluster 6 result
str(hallmark_b_cell_results[["6"]])

# View the first few rows of the data
head(hallmark_b_cell_results[["6"]])

# Check the names of the columns (if applicable)
names(hallmark_b_cell_results[["6"]])

# Define a function to extract pathways from a cluster result
extract_pathways <- function(cluster_result) {
  cluster_result@result$Description
}

# Extract hallmark pathway descriptions for each B-cell cluster
hallmark_b_cell_pathways <- lapply(hallmark_b_cell_results, extract_pathways)

# View the first few extracted pathways for each cluster
head(hallmark_b_cell_pathways)

# Show the first few pathways for each of the first few clusters
concise_output <- lapply(hallmark_b_cell_pathways, function(paths) head(paths, 10))
concise_output

# View the first few clusters and their first few pathways
concise_output <- head(hallmark_b_cell_pathways, 2)  # Adjust the number '2' to show more or fewer clusters
concise_output <- lapply(concise_output, function(paths) head(paths, 5))
concise_output

# Identify common pathways across all B-cell clusters
common_pathways <- Reduce(intersect, hallmark_b_cell_pathways)
head(common_pathways)

# Identify unique pathways for each cluster
unique_pathways <- lapply(hallmark_b_cell_pathways, function(pathways) {
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
top_10_pathways <- lapply(hallmark_b_cell_results, extract_top_pathways, n = 10)
top_20_pathways <- lapply(hallmark_b_cell_results, extract_top_pathways, n = 20)

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

# Prepare data for bar chart
top_10_pathways_df <- data.frame(
  Cluster = rep(names(top_10_pathways), each = 10),
  Pathway = unlist(top_10_pathways),
  Rank = rep(1:10, times = length(top_10_pathways))
)

# Plot a bar chart for the top 10 pathways for each cluster
ggplot(top_10_pathways_df, aes(x = factor(Rank), y = Cluster, fill = Pathway)) +
  geom_tile(color = "white") +
  scale_fill_viridis_d(option = "C") +
  labs(x = "Rank", y = "Cluster", fill = "Pathway") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Top 10 Hallmark Pathways Across B-cell Clusters")

# Save the plot as a file in the specified path
ggsave("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/hallmark/Top_10_hallmark_Pathways_heatmap_b_cell.png", width = 10, height = 8)

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
  ggtitle("Top 10 Hallmark Pathways Across B-cell Clusters") +
  coord_flip()  # Flip coordinates for better readability

# Save the plot
ggsave("/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis/hallmark/top_10_hallmark_pathways_barplot_b_cell.png", width = 12, height = 8)
