#!/usr/bin/Rscript

#### inferCNV analysis of individual samples to identify tumor cells (pt2: after tumor vs non-tumor separation)
#### Author: Jana Biermann, PhD

print(paste('Start:', Sys.time()))

library(dplyr)
library(Seurat)
library(infercnv)
library(stringr)

pat <- commandArgs()[6]

# Read in Seurat object
patient <- readRDS(file = paste0('inferCNV_subcluster_', pat, '/', pat, '_cnv.rds'))

# Set up counts matrix
counts_matrix <- as.data.frame(GetAssayData(object = patient, slot = 'counts'))

# Annotate tumor cells vs non-tumor cells based on inferCNV pt1
annot <- patient@meta.data$tumor
annot <- as.data.frame(annot)
colnames(annot) <- c('V1')
annot$V1 <- as.character(annot$V1)
rownames(annot) <- as.character(colnames(counts_matrix))
annot <- na.omit(annot)

# Load gene order file
gene_order <- read.table('misc/gen_pos_gene_name.txt', header = F, row.names = 1)

# Create inferCNV object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix = counts_matrix,
                                    annotations_file = annot,
                                    gene_order_file = gene_order,
                                    ref_group_names = 'immune')

# Run inferCNV
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff = 0.1,
                             out_dir = paste0('inferCNV_final_', pat),
                             cluster_by_groups = F,
                             denoise = T,
                             HMM = T,
                             output_format = 'pdf',
                             num_threads = 16)

# Save object
saveRDS(patient, paste0('inferCNV_final_', pat, '/', pat, '_final.rds'))

print(paste('End:', Sys.time()))