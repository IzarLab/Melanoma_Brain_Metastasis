#!/usr/bin/env Rscript

### title: Global analysis of integrated object 
### author: Jana Biermann, PhD

print(Sys.time())

library(Seurat)
library(dplyr)
library(ggplot2)
library(gplots)
'%notin%' <- Negate('%in%')


colBP <- c('#A80D11', '#008DB8')
colSCSN <- c('#E1AC24', '#288F56')

#### UMAPs
seu <- readRDS('data/MBPM/data_MBPM_scn.rds')
seu <- subset(seu, cell_type_int %notin% c('Low-quality cells', 'Doublets', 'Contamination', 'Undetermined'))
seu <- ScaleData(seu) %>% RunPCA(npcs = 85)
seu <- RunUMAP(seu, dims = 1:85, spread = 1.5, min.dist = 0.1)


# Plots
ifelse(!dir.exists(file.path(paste0('data/MBPM/global'))), 
       dir.create(file.path(paste0('data/MBPM/global')), recursive = T), FALSE)
pdf(file = 'data/MBPM/global/plots_MBPM_global.pdf')
DimPlot(seu, reduction = 'umap', label = F, group.by = 'sequencing', shuffle = T, raster = T, cols = colSCSN) + NoLegend()
DimPlot(seu, reduction = 'umap', label = F, group.by = 'sequencing', shuffle = T, raster = T, cols = colSCSN)

DimPlot(seu, reduction = 'umap', label = F, group.by = 'organ', shuffle = T, raster = T, cols = colBP) + NoLegend()
DimPlot(seu, reduction = 'umap', label = F, group.by = 'organ', shuffle = T, raster = T, cols = colBP)

DimPlot(seu, reduction = 'umap', label = F, group.by = 'cell_type_main', shuffle = T, raster = T, repel = T) + NoLegend()
DimPlot(seu, reduction = 'umap', label = T, group.by = 'cell_type_main', shuffle = T, raster = T, repel = T)

DimPlot(seu, reduction = 'umap', label = F, group.by = 'cell_type_int', shuffle = T, raster = T, repel = T) + NoLegend()
DimPlot(seu, reduction = 'umap', label = T, group.by = 'cell_type_int', shuffle = T, raster = T, repel = T) + 
  guides(col = guide_legend(nrow = 21, override.aes = list(size = 3))) + 
  theme(legend.text = element_text(size = 8))

DimPlot(seu, reduction = 'umap', label = F, group.by = 'cell_type_fine', shuffle = T, raster = T, repel = T) + NoLegend()
DimPlot(seu, reduction = 'umap', label = T, group.by = 'cell_type_fine', shuffle = T, raster = T, repel = T, label.size = 3) + 
  guides(color = guide_legend(ncol = 1), override.aes = list(size = 1)) + 
  theme(legend.text = element_text(size = 5), legend.key.size = unit(0.2, 'cm'))

FeaturePlot(seu, features = 'proportion_scaled_cnv_avg', min.cutoff = 'q05', max.cutoff = 'q95', 
            order = T, raster = T, pt.size = 0.7) + scale_color_viridis(direction = -1)
dev.off()

print(Sys.time())