# Suppress warnings
options(warn = -1)

# load packages
library(dplyr)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(BPCells)
library(Azimuth)
library(future)

library(DESeq2)

library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(EnhancedVolcano)

# change the current plan to access parallelization
plan("multisession", workers = availableCores())
plan()

# increase size of plots from default
options(repr.plot.width = 14, 
        repr.plot.height = 14) # from 7, 7

# load data from h5 file
file_path <- "data/16plex_900k_32_NSCLC_multiplex_count_filtered_feature_bc_matrix.h5"
nsclc_data <- open_matrix_10x_hdf5(
    path = file_path
)

# Write the matrix to a directory
mat_dir <- "data/nsclc_counts"
write_matrix_dir(
    mat = nsclc_data, 
    dir = mat_dir)

# Now that we have the matrix on disk, we can load it
nsclc_mat <- open_matrix_dir(dir = mat_dir)
nsclc_mat <- Azimuth:::ConvertEnsembleToSymbol(mat = nsclc_mat, 
                                               species = "human")

# Create Seurat Object
nsclc <- CreateSeuratObject(counts = nsclc_mat)

# Now that we have the matrix on disk, we can load it
mat_dir <- "data/nsclc_counts"
nsclc_mat <- open_matrix_dir(dir = mat_dir)
nsclc_mat <- Azimuth:::ConvertEnsembleToSymbol(mat = nsclc_mat, 
                                               species = "human")

# Create Seurat Object
(nsclc <- CreateSeuratObject(counts = nsclc_mat))

nsclc[['sample']] <- as.integer(sub('.*-', '', colnames(nsclc)))
sample_table <- read.csv('data/16plex_900k_32_NSCLC_multiplex_aggregation.csv')
sample_table$sample <- (1:nrow(sample_table))
nsclc[[]] <- merge(nsclc[[]], sample_table, by = "sample", all.x = TRUE)
nsclc[[]]

# mitochondrial counts to metadata
nsclc[["percent.mt"]] <- PercentageFeatureSet(nsclc, pattern = "^MT-")
summary(nsclc[["percent.mt"]])

# Visualize QC metrics as a violin plot
VlnPlot(nsclc, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)

summary(nsclc[[c("nFeature_RNA", "nCount_RNA", "percent.mt")]])

# filter out cells with too low/high counts, too high mitochondrial counts, too many/few genes
nsclc
nsclc <- subset(nsclc, 
                subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & 
                            percent.mt < 5 & 
                            nCount_RNA > 600 & nCount_RNA < 6000)
nsclc
summary(nsclc[[c("nFeature_RNA", "nCount_RNA", "percent.mt")]])

# visualize QC metrics post-filter
VlnPlot(nsclc, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        alpha = 0,
        ncol = 3)

# SLURM submission script
file_to_read <- "~/tools/slurm_scripts/Rscript_submit_himem_1node_12cpus_mem500_231125_slurm.sh"
file_content <- readLines(file_to_read)
cat(file_content, sep = "\n")

# R script for data normalization, selection of highly variable genes, data scaling, dimension reduction, clustering, and visualization
file_to_read <- "~/tools/R_scripts/scrnaseq_sctransform_to_umap_231127.r"
file_content <- readLines(file_to_read)
cat(file_content, sep = "\n")

# load Seurat object output by SLURM submission
filename <- 'nsclc_2023_11_27_08_03_25.rds'
nsclc <- readRDS(filename)
nsclc

# visualize UMAP and group by the 32 samples of origin
print(
DimPlot(nsclc, reduction = "umap", 
        group.by = c('description'),
       raster = FALSE)
    )

# visualize UMAP and group by the 5 biological states
print(
DimPlot(nsclc, reduction = "umap", 
        group.by = c('disease.state'), 
       raster = FALSE)
    )

# R script for Azimuth cell annotation using the Human Lung Cell Atlas as a reference
file_to_read <- "~/tools/R_scripts/scrnaseq_run_azimuth_231128.r"
file_content <- readLines(file_to_read)
cat(file_content, sep = "\n")

# read Azimuth annotated rds
nsclc <- readRDS('scrnaseq_run_azimuth_231128_2023_11_28_11_01_02.rds')
nsclc

# visualization of Azimuth-annotated data
print(
DimPlot(nsclc, group.by = "predicted.ann_level_3", 
        label = TRUE, label.size = 5, 
        raster = FALSE)  + NoLegend()
    )

# visualization of Azimuth-annotated data
print(
DimPlot(nsclc, group.by = "predicted.ann_level_4", 
        label = TRUE, label.size = 4, 
        raster = FALSE)  + NoLegend()
    )

# pseudobulk by donor
bulk_by_donor <- AggregateExpression(nsclc, return.seurat = TRUE, 
                                     assays = "RNA", 
                                     group.by = "description")
Cells(bulk_by_donor)

# extract the five biological groups (conditions = LUAD LUAD-Tx LUSC LUSC-Tx NADJ) from sample ids
sample_ids <- Cells(bulk_by_donor)

# Function to extract the desired parts
extract_parts <- function(input) {
  # Split the string by '-'
  parts <- strsplit(input, "-")[[1]]
  
  # Extract the parts based on positions
  after_last_hyphen <- tail(parts, 1)
  after_second_last_hyphen <- tail(head(parts, -1), 1)
  remainder <- paste(head(parts, -2), collapse = "-")
  
  return(list(after_last_hyphen, after_second_last_hyphen, remainder))
}

# Apply the function to the input vector
result <- lapply(sample_ids, extract_parts)

# Extracting the three resulting vectors
after_last_hyphen <- sapply(result, `[[`, 1)
after_second_last_hyphen <- sapply(result, `[[`, 2)
condition <- sapply(result, `[[`, 3)

print(condition)

# create coldata for DESeq dataset
coldata <- data.frame(condition = as.factor(condition), 
                     row.names = colnames(bulk_by_donor[['RNA']]$counts))
coldata

# create DESeq object
countdata <- bulk_by_donor[['RNA']]$counts
dds <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = coldata,
                                 design = ~ condition)
dds

# drop non-informative (low/no count) genes
# at least 4 samples (columns) with a count of 10 or higher
keep <- rowSums(counts(dds) >= 10) >= 4 
dds <- dds[keep, ]
dds

# heatmap of sample distances
rld <- rlog(dds)

sampleDists <- dist(t(assay(rld)))

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- colnames(rld)
colnames(sampleDistMatrix) <- NULL 
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

# PCA plot
object <- rld
ntop <- 500 # number of variable genes to use for PCA
intgroup <- 'condition'
pcsToUse = 1:2

rv <- rowVars(assay(object))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
pca <- prcomp(t(assay(object)[select, ]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)

intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
group <- colData(object)[[intgroup]]

pcs <- paste0("PC", pcsToUse)
d <- data.frame(V1 = pca$x[, pcsToUse[1]], V2 = pca$x[, pcsToUse[2]], 
                group = group, intgroup.df, name = colnames(object))
colnames(d)[1:2] <- pcs

ggplot(data = d, aes(x = PC1, y = PC2, color = group)) + 
    geom_point(size = 4) + 
    xlab(paste0(pcs[1], ": ", round(percentVar[pcsToUse[1]] * 100), "% variance")) + 
    ylab(paste0(pcs[2], ": ", round(percentVar[pcsToUse[2]] * 100), "% variance")) + 
    coord_fixed() + 
    theme_gray(base_size = 24)

# differential expression analysis
dds <- DESeq(dds)
dds

res <- results(dds)
summary(res)
res

# using more stringent thresholds: BH FDR adjusted p-value 0.05 vs 0.1
res <- results(dds, alpha = 0.05)
summary(res)
res

# significant genes with strongest down-regulation in NADJ
resSig <- subset(res, padj < 0.05)
head(resSig[ order(resSig$log2FoldChange), ])

# significant genes with strongest up-regulation in NADJ
head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ])

EnhancedVolcano(as.data.frame(res), x = 'log2FoldChange', y = 'padj', 
                lab = rownames(res), 
               pCutoff = 0.05, FCcutoff = 1.5, 
               xlim = c(-10, 10))

res[c('AGER', 'EMP2'),]

# compare NADJ to LUSC
res <- results(dds, alpha = 0.05, 
               contrast = c("condition", "NADJ", "LUSC"))
summary(res)
res

EnhancedVolcano(as.data.frame(res), x = 'log2FoldChange', y = 'padj', 
                lab = rownames(res), 
               pCutoff = 0.05, FCcutoff = 1.5, 
               xlim = c(-10, 10))

res[c('AGER', 'EMP2'),]

# Explicitly close multisession workers by switching plan
plan(sequential)
plan()

sessionInfo()

