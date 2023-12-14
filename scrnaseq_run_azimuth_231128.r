# clear workspace
rm(list = ls())
gc()

args <- commandArgs(TRUE) # to access command line arguments
print(args)

############################################################################################################################################################


work_dir <- args[1]                                                                                                            ### WORKING/OUTPUT DIRECTORY
input_filename <- "nsclc_2023_11_27_08_03_25.rds"                                                                              ### INPUT FILENAME
output_filename <- paste('scrnaseq_run_azimuth_231128', system("date '+_%Y_%m_%d_%H_%M_%S'", intern=TRUE), '.rds', sep = '')   ### OUTPUT FILENAME

# load packages
library(BPCells)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(Azimuth)
library(SeuratData)
library(future)
library(dplyr)

options(future.globals.maxSize = 16e+09) # 16G

# change the current plan to access parallelization
plan("multisession", workers = availableCores())
plan()

# load data
cat(paste('\n', system("date", intern=TRUE), '---> loading data...', '\n'))
setwd(work_dir)
nsclc <- readRDS(input_filename)
cat(paste('\n', system("date", intern=TRUE), '---> loading data COMPLETE.', '\n'))

# annotate with Azimuth
cat(paste('\n', system("date", intern=TRUE), '---> starting annotate with Azimuth...', '\n'))
options(timeout = 1200)
nsclc <- RunAzimuth(nsclc, reference = "lungref")
cat(paste('\n', system("date", intern=TRUE), '---> annotate with Azimuth COMPLETE.', '\n'))

# save objects
cat(paste('\n', system("date", intern=TRUE), '---> saving output...', '\n'))
saveRDS(nsclc, file = output_filename)
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


