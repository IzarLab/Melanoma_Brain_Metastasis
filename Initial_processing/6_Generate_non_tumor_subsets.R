#!/usr/bin/env Rscript

#### Generate non-tumor subsets for cell type annotation
#### Author: Jana Biermann, PhD

library(dplyr)
library(Seurat)
library(ggplot2)
library(gplots)

cohort <- commandArgs()[6]

# Read-in object
seu <- readRDS(paste0('data/MBPM/', cohort, '/data_', cohort, '_anchor2000_dims30.rds'))
seu <- subset(seu, malignant == 'non-malignant')

# Rerun workflow
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 75)
ElbowPlot(seu, ndims = 75)

if (cohort == 'MBM_sc') {
  seu <- RunUMAP(seu, dims = 1:35)
  seu <- FindNeighbors(seu, dims = 1:35)
  seu <- FindClusters(seu, resolution = 0.3)
}

if (cohort == 'MBM_sn') {
  seu <- RunUMAP(seu, dims = 1:40)
  seu <- FindNeighbors(seu, dims = 1:40)
  seu <- FindClusters(seu, resolution = 0.3)
}

if (cohort == 'MPM_sn') {
  seu <- RunUMAP(seu, dims = 1:40)
  seu <- FindNeighbors(seu, dims = 1:40)
  seu <- FindClusters(seu, resolution = 0.5)
}

# Save object
saveRDS(seu, file = paste0('data/MBPM/', cohort, '/data_', cohort, '_notum.rds'))

### stats
stats <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 11))
colnames(stats) <- c('sample', 'n_features', 'n_cells', 'median_features', 'median_counts')
rownames(stats) <- cohort
stats$sample <- cohort
stats$n_features <- dim(seu@assays$integrated@data)[1]
stats$n_cells <- dim(seu@assays$integrated@data)[2]
stats$median_features <- round(median(seu@meta.data$nFeature_RNA))
stats$median_counts <- round(median(seu@meta.data$nCount_RNA))

pdf(file = 'data/MBPM/', cohort, '/plots_', cohort, '_notum.pdf')
textplot(t(stats), cex = 1.2, halign = 'left')
DimPlot(seu, reduction = 'pca', group.by = 'ident', shuffle = T, raster = T)
DimPlot(seu, reduction = 'pca', group.by = 'patient', shuffle = T, raster = T)
DimPlot(seu, reduction = 'pca', group.by = 'malignant', shuffle = T, raster = T)
ElbowPlot(seu, ndims = 75)
DimPlot(seu, reduction = 'umap', label = T, group.by = 'ident', shuffle = T, raster = T)
DimPlot(seu, reduction = 'umap', label = F, group.by = 'patient', shuffle = T, raster = T)
FeatureScatter(seu, 'G2M.Score', 'S.Score', group.by = 'cell_cycle')
DimPlot(seu, reduction = 'umap', label = T, group.by = 'cell_cycle', shuffle = T, raster = T)
DimPlot(seu, reduction = 'umap', label = T, group.by = 'celltype_bped_fine', repel = T,
        label.size = 2.5, shuffle = T, raster = T) + 
  guides(col = guide_legend(nrow = 20,override.aes = list(size = 5))) + 
  theme(legend.text = element_text(size = 6))
DimPlot(seu, reduction = 'umap', label = T, group.by = 'celltype_hpca_main', repel = T,
        label.size = 2.5, shuffle = T, raster = T) + 
  guides(col = guide_legend(nrow = 20,override.aes = list(size = 5))) + 
  theme(legend.text = element_text(size = 6))

FeaturePlot(seu, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'doublet_scores'),
                  order = T, raster = T)
FeaturePlot(seu, features = c('proportion_scaled_cnv_avg', 'proportion_cnv_avg', 'has_cnv_avg'), 
                  min.cutoff = 'q05', max.cutoff = 'q95', order = T, raster = T)

FeaturePlot(seu, features = c('MLANA', 'PMEL', 'MET', 'TYR'), min.cutoff = 'q05',
            max.cutoff = 'q95', order = T, raster = T)
FeaturePlot(seu, features = c('PTPRC', 'CD8A', 'CD68', 'MKI67'), min.cutoff = 'q05',
                  max.cutoff = 'q95', order = T, raster = T)
dev.off()
