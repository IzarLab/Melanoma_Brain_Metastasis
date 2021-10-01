#!/usr/bin/Rscript

#### inferCNV analysis of individual samples to identify tumor cells
#### Author: Jana Biermann, PhD

print(paste('Start:', Sys.time()))

library(dplyr)
library(Seurat)
library(infercnv)
library(stringr)


pat <- commandArgs()[6]
percentile <- commandArgs()[7]

# Read in Seurat object
patient <- readRDS(file = paste0('data/', pat, '/data_', pat, '_cb.rds'))

# Set up counts matrix
counts_matrix <- as.data.frame(GetAssayData(object = patient, slot = 'counts'))

# Annotate immune cells based on SingleR output
annot <- as.data.frame(patient@meta.data$celltype_bped_fine)
colnames(annot) <- c('V1')
annot$V1 <- as.character(annot$V1)
rownames(annot) <- as.character(colnames(counts_matrix))
immune <- c('CD4+ T-cells', 'CD4+ Tcm', 'CD4+ Tem', 'CD8+ T-cells', 'CD8+ Tcm', 'CD8+ Tem', 'Class-switched memory B-cells', 'DC', 'Eosinophils', 'Macrophages', 'Macrophages M1', 'Macrophages M2', 'Memory B-cells', 'Monocytes', 'naive B-cells', 'Neutrophils', 'NK cells', 'Plasma cells', 'Tregs')
annot$cell <- ifelse(annot$V1 %in% immune, 'immune', 'non-immune')
annot$V1 <- NULL

# Load gene order file
gene_order <- read.table('misc/gen_pos_gene_name.txt', header = F, row.names = 1)

# Create inferCNV object
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = counts_matrix,
                                     annotations_file = annot,
                                     gene_order_file = gene_order,
                                     ref_group_names = 'immune')

# Run inferCNV
infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff = 0.1,
                              out_dir = paste0('inferCNV_subcluster_', pat),
                              cluster_by_groups = F,
                              denoise = T,
                              HMM = T,
                              analysis_mode = 'subclusters',
                              output_format = 'pdf',
                              num_threads = 4)

# Identify malignant cells
immune <- c('CD4+ T-cells', 'CD4+ Tcm', 'CD4+ Tem', 'CD8+ T-cells', 'CD8+ Tcm', 'CD8+ Tem', 'Class-switched memory B-cells', 'DC', 'Eosinophils', 'Macrophages', 'Macrophages M1', 'Macrophages M2', 'Memory B-cells', 'Monocytes', 'naive B-cells', 'Neutrophils', 'NK cells', 'Plasma cells', 'Tregs')
patient[['immune']] <- ifelse(patient@meta.data$celltype_bped_fine %in% immune, 'immune', 'non-immune')
patient <- add_to_seurat(patient, paste0('inferCNV_subcluster_', pat))
cnv_cols <- grep('proportion_scaled_cnv_chr', names(patient@meta.data), value = T)
cnvs <- patient@meta.data[, cnv_cols]
patient@meta.data$cnv_avg <- rowMeans(cnvs)

# Step 1: Patient-specific percentile as first threshold (based on cell composition of sample)
quan <- quantile(patient@meta.data[rownames(patient@meta.data[patient@meta.data$immune == 'non-immune', ]), 'cnv_avg'], percentile)
patient@meta.data$malignant <- ifelse(patient@meta.data$cnv_avg > quan, 'malignant', 'non-malignant')

# Step 2: Using SingleR annotation of cell as second threshold to identifying tumor cells
patient[['tumor']] <- ifelse(patient@meta.data$malignant == 'malignant' & patient@meta.data$celltype_bped_fine == 'Melanocytes', 'tumor',
                             ifelse(patient@meta.data$immune == 'immune', 'immune', NA))

# Save object
saveRDS(patient, paste0('inferCNV_subcluster_', pat, '/', pat, '_cnv.rds'))


print(paste('End:', Sys.time()))