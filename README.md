# scRNA-seq

*You may need to refresh a couple of times if you see the message "Unable to render code block" when loading the (pdf of the) Jupyter notebook above*

I wrote these R scripts to analyze a large, single-cell RNA sequencing (scRNA-seq) dataset (~900,000 cells) using R and various Bioconductor bioinformatics packages. 
 
### Overall Functionality:
The scripts accomplish the following tasks:
1. Load necessary R packages for analysis (Seurat, DESeq2, pheatmap, ggplot2, etc.).
2. Read scRNA-seq data from an h5 file and converts it into a Seurat object.
3. Conduct quality control (QC) analysis, filtering out low-quality cells based on gene counts, mitochondrial content, etc.
4. Generate QC plots to visualize the distribution of various QC metrics.
5. Perform normalization, variable gene selection, scaling, dimension reduction (UMAP), clustering, and visualization using Seurat's workflows.
6. Use Azimuth for cell annotation based on the Human Lung Cell Atlas.
7. Conduct pseudobulk analysis, differential expression analysis (DESeq2), and visualization of significant genes.
8. Print session information, clears workspace, and manages parallelization using the future package.

### Input Files:
- Input: Single-cell RNA sequencing data in an h5 file (`16plex_900k_32_NSCLC_multiplex_count_filtered_feature_bc_matrix.h5`).
- Additional CSV files for sample information and annotation.

### Required Packages and Tools:
- R Packages: Seurat, SeuratDisk, Azimuth, DESeq2, pheatmap, ggplot2, EnhancedVolcano, BPCells, future, dplyr.

### Outputs:
- Seurat objects at different analysis stages.
- QC plots (Violin plots, heatmaps).
- Intermediate data files for normalization, clustering, and differential expression analysis.
- Visualizations (UMAP plots, volcano plots) highlighting various aspects of the data.

[The dataset](https://www.10xgenomics.com/resources/datasets/aggregate-of-900k-human-non-small-cell-lung-cancer-and-normal-adjacent-cells-multiplexed-samples-16-probe-barcodes-1-standard)
 is from 10X Genomics. 
 
