# Single-Cell-RNA-seq-Pathway-Analysis-using-GO-Hallmark-and-C2-Gene-Sets

### Pathway Analysis Overview

**Pathway analysis** is a crucial bioinformatics approach used to interpret gene expression data by identifying biological pathways that are significantly enriched in a given set of genes. By understanding these pathways, researchers can uncover the biological mechanisms underlying various conditions or phenotypes. This repository contains the code and results for a pathway analysis of single-cell RNA-seq data, utilizing three key gene sets: GO (Gene Ontology), Hallmark, and C2

(Curated Pathways). The analysis was performed at both the whole cluster level and for individual clusters, with pathways annotated to define subset cell types within each cluster.

### Data and Methods

-   **Data Sources**: Single-cell RNA-seq data.

-   **Gene Sets**: GO (Biological Process), Hallmark, and C2:CP

    (Curated Pathways).

-   **Tools Used**:

    -   `Seurat` for data processing and differential expression analysis.

    -   `clusterProfiler` for enrichment analysis.

    -   `msigdbr` for accessing MSigDB gene sets.

    -   `ggplot2`, `pheatmap`, `enrichplot` for visualization.

-   **Analysis Flow**:

    1.  Load and preprocess the Seurat object.

    2.  Perform differential expression analysis (DEA) across clusters.

    3.  Conduct enrichment analysis for GO, Hallmark, and C2

        gene sets.

    4.  Visualize top pathways and annotate subset cell types.

### Reproducing the Analysis

1.  **Set Up the Environment**: Install necessary R packages.

    ```         
    install.packages('tidyverse')
    BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"))
    install.packages("msigdbr")
    ```

2.  **Run the Analysis**:

    Execute the R scripts

3.  **Visualize Results**: Plots are saved in the respective subdirectories under `pathway_analysis/`.

### Results

-   **GO Enrichment**:

    -   Identified key biological processes across all clusters and specific clusters.

-   **Hallmark Gene Sets**:

    -   Highlighted hallmark pathways relevant to the conditions studied.

-   **C2:CP**

    **Pathways**:

    -   Provided insight into curated pathways involved in each cluster.

### Future Directions

-   Extend the analysis to additional gene sets or custom pathway databases.

-   Refine cell type annotations based on pathway enrichment patterns.
