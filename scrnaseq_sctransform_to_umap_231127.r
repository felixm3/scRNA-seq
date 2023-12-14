# clear workspace
rm(list = ls())
gc()

args <- commandArgs(TRUE) # to access command line arguments
print(args)

############################################################################################################################################################


work_dir <- args[1]                                                                                                  ### WORKING/OUTPUT DIRECTORY
filename <- paste('scrnaseq_sctransform_to_umap_231127', system("date '+_%Y_%m_%d_%H_%M_%S'", intern=TRUE), '.rds', sep = '')   ### OUTPUT FILENAME

# load packages
library(BPCells)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(Azimuth)
library(future)
library(dplyr)

options(future.globals.maxSize = 16e+09) # 16G

# change the current plan to access parallelization
plan("multisession", workers = availableCores())
plan()

# load data
cat(paste('\n', system("date", intern=TRUE), '---> loading data...', '\n'))
setwd(work_dir)
filename <- "nsclc_2023_11_27_08_03_25.rds"
nsclc <- readRDS(filename)
cat(paste('\n', system("date", intern=TRUE), '---> loading data COMPLETE.', '\n'))

# normalize & cluster
cat(paste('\n', system("date", intern=TRUE), '---> starting sctransform to umap...', '\n'))
nsclc <- SCTransform(nsclc, ncells = 100000, conserve.memory = TRUE) %>% 
            RunPCA() %>% 
                FindNeighbors(dims = 1:30) %>% 
                    FindClusters() %>% 
                        RunUMAP(dims = 1:30)
cat(paste('\n', system("date", intern=TRUE), '---> sctransform to umap COMPLETE.', '\n'))

# save objects
cat(paste('\n', system("date", intern=TRUE), '---> saving output...', '\n'))
saveRDS(nsclc, file = filename)
cat(paste('\n', system("date", intern=TRUE), '---> save objects COMPLETE.', '\n'))

# cleanup
rm(nsclc)
gc()

# Explicitly close multisession workers by switching plan
plan(sequential)
plan()

############################################################################################################################################################

## print session info ##
cat("\nSession Info below: \n")
sessionInfo()

# clear workspace
rm(list = ls())
gc()


