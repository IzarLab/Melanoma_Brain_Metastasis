#!/usr/bin/env Rscript

#### Estimate CNAs using inferCNV with patient argument provided
#### Author: Jana Biermann, PhD

print(paste('Start:', Sys.time()))

library(dplyr)
library(Seurat)
library(infercnv)
library(stringr)
library(gplots)
library(ggplot2)
library(viridis)

# Get patient argument
pat <- commandArgs()[6]

# Read-in Seurat object
system(paste0("aws s3 sync s3://snrna-seq/Seurat/", pat, "/ data/", pat, "/ ",
              "--exclude '*' --include 'data_", pat, "_cb.rds'"))
seu <- readRDS(file = paste0('data/', pat, '/data_', pat, '_cb_DF.rds'))

# Set up counts matrix
counts_matrix <- as.data.frame(GetAssayData(object = seu, slot = 'counts'))

# Annotate immune cells based on SingleR output
annot <- as.data.frame(seu$celltype_bped_fine)
annot[is.na(annot)] <- 'unknown'
immune <- c('CD4+ T-cells', 'CD4+ Tcm', 'CD4+ Tem', 'CD8+ T-cells', 'CD8+ Tcm',
            'CD8+ Tem', 'Class-switched memory B-cells', 'DC', 'Eosinophils', 'Macrophages',
            'Macrophages M1', 'Macrophages M2', 'Memory B-cells', 'Monocytes', 'naive B-cells',
            'Neutrophils', 'NK cells', 'Plasma cells', 'Tregs')
annot$cell <- ifelse(annot$`seu$celltype_bped_fine` %in% immune, 'immune', 'non-immune')
annot$`seu$celltype_bped_fine` <- NULL

# Load gene order file
gene_order <- read.table('misc/refdata-gex-GRCh38-2020-A_gen_pos.txt', header = F,row.names = 1)


# Create infercnv object
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
seu[['immune']] <- ifelse(seu$celltype_bped_fine %in% immune, 'immune', 'non-immune')
seu <- add_to_seurat(seu, paste0('inferCNV_subcluster_', pat))
cnv_cols <- grep('proportion_scaled_cnv_chr', names(seu@meta.data), value = T)
cnvs <- seu@meta.data[, cnv_cols]
seu$cnv_avg <- rowMeans(cnvs)

# For the majority of samples cut-off >0.1 (exceptions below)
seu$malignant <- ifelse(seu$cnv_avg > 0.1, 'malignant', 'non-malignant')

## MBM16_sn, MPM04_sn
#seu$malignant<-ifelse(seu$cnv_avg > 0.1 & seu@reductions$pca@cell.embeddings[,1] < 5, 'malignant','non-malignant')

## MBM07_sn
#seu$malignant<-ifelse(seu$cnv_avg > 0.1 & seu@reductions$pca@cell.embeddings[,2] > (-25), 'malignant','non-malignant')

## MBM11_sn, MBM19_sn, MBM04_sc
#seu$malignant<-ifelse(seu$cnv_avg > 0.1 & seu$celltype_bped_fine=='Melanocytes', 'malignant','non-malignant')
#seu$malignant<-ifelse(is.na(seu$malignant), 'malignant',seu$malignant)


# Add CNV metrics
cnv_cols <- grep('proportion_scaled_cnv_chr', names(seu@meta.data), value = T)
cnvs <- seu@meta.data[, cnv_cols]
seu$proportion_scaled_cnv_avg <- rowMeans(cnvs)

cnv_cols <- grep('proportion_cnv_chr', names(seu@meta.data), value = T)
cnvs <- seu@meta.data[, cnv_cols]
seu$proportion_cnv_avg <- rowMeans(cnvs)

cnv_cols <- grep('has_cnv_chr', names(seu@meta.data), value = T)
cnvs <- seu@meta.data[, cnv_cols]
seu$has_cnv_avg <- rowMeans(cnvs)


# Save object
saveRDS(seu, paste0('inferCNV_subcluster_', pat, '/data_', pat, '_cnv.rds'))

# Plots
pdf(paste0('inferCNV_subcluster_', pat, '/plots_', pat, '_cnv.pdf'))
textplot(addmargins(table(seu$immune, seu$malignant)), cex = 1.2, halign = 'left')
textplot(table(seu$celltype_bped_fine, seu$malignant), cex = 0.9, halign = 'left')
hist(seu$cnv_avg, breaks = 100, main = 'Average absolute CNV level; all cells',
     xlab = 'Average absolute CNV proportion')
abline(v = q10, col = 'red')
abline(v = 0.1, col = 'blue')

DimPlot(seu, reduction = 'pca', group.by = 'malignant')
DimPlot(seu, reduction = 'umap', group.by = 'malignant')
DimPlot(seu, reduction = 'umap', group.by = 'immune')
DimPlot(seu, reduction = 'umap', label = T, group.by = 'celltype_bped_main', repel = T,
        label.size = 2.5, raster = T, shuffle = T) + 
  guides(col = guide_legend(nrow = 30,override.aes = list(size = 5))) + 
  theme(legend.text = element_text(size = 5))

FeaturePlot(seu, features = c('proportion_scaled_cnv_avg'), min.cutoff = 'q05',
            max.cutoff = 'q95', order = T, raster = T) + scale_color_viridis(direction = -1)
FeaturePlot(seu, features = c('MLANA', 'PMEL', 'MET', 'proportion_scaled_cnv_avg'),
            min.cutoff = 'q05', max.cutoff = 'q95', order = T, raster = T)
VlnPlot(seu, features = c('proportion_scaled_cnv_avg'), group.by = 'celltype_bped_main',
        pt.size = 0)
VlnPlot(seu, features = c('proportion_scaled_cnv_avg'), pt.size = 0)
DimPlot(seu, reduction = 'umap', shuffle = T, raster = T, label = T)
dev.off()

system(paste0("aws s3 sync inferCNV_subcluster_", pat, "/ s3://snrna-seq/inferCNV/", pat,
              "/inferCNV_subcluster_", pat, "/ --exclude '*DS_S*' --exclude '._*' "))

print(paste('End:', Sys.time()))