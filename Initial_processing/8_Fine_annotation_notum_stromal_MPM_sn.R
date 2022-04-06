#!/usr/bin/env Rscript

#### Fine cell type annotation of stromal cell subset for MPM_sn
#### Author: Jana Biermann, PhD

library(dplyr)
library(Seurat)
library(ggplot2)
library(gplots)
library(viridis)

path.ct <- 'data/cell_type_DEG/MPM_sn/stromal/'
filename <- 'MPM_sn_stromal'

seu <- readRDS('data/cell_type_DEG/MPM_sn/main/data_MPM_sn_notum.rds')

# Subset to stromal cells
seu <- subset(seu, cell_type_main %in% c('Stromal cells', 'Low-quality cells'))
seu <- ScaleData(object = seu)
seu <- RunPCA(object = seu, npcs = 75)
ElbowPlot(seu, ndims = 75)
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.3)
seu <- RunUMAP(object = seu, dims = 1:30)
DimPlot(seu, reduction = 'umap', label = T, group.by = 'ident', shuffle = T, raster = T)

# Cell cycle assignment
DefaultAssay(seu) <- 'RNA'
seu <- CellCycleScoring(seu, s.features = cc.genes.updated.2019$s.genes,
                        g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
DefaultAssay(seu) <- 'integrated'
seu$cell_cycle <- ifelse(seu$G2M.Score > 0.1 | seu$S.Score > 0.1, 'cycling',
                         'non-cycling')

# Differential gene expression (DGE) based on clusters
markers <- FindAllMarkers(seu, only.pos = TRUE, assay = 'RNA', min.pct = 0.25,
                          logfc.threshold = 0.25, test.use = 'MAST')
write.csv(markers, paste0(path.ct, 'markers_', filename, '.csv'), row.names = F)
markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  top_n(15, wt = avg_log2FC) -> tops

DefaultAssay(seu) <- 'RNA'
seu <- NormalizeData(seu) %>%
  ScaleData()
DefaultAssay(seu) <- 'integrated'

pdf(file = paste0(path.ct, 'markers_', filename, '_heatmap.pdf'), width = 10, height = 18)
DoHeatmap(seu, features = tops$gene, group.by = 'ident', raster = T, assay = 'RNA')
dev.off()

# Fine cell type assignment after manual annotation based on DGE
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.3 %in% c(3, 5, 6, 7),
                             'Low-quality tumor/immune', 'NA')
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.3 %in% c(4, 8), 'Epithelial cells',
                             seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.3 %in% c(10), 'Endothelial cells',
                             seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.3 %in% c(0, 1, 9), 'CAFs',
                             seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.3 %in% c(2), 'Pericytes',
                             seu$cell_type_fine)

# Differential gene expression (DGE) based on cell types
Idents(seu) <- seu$cell_type_fine
markers <- FindAllMarkers(seu, only.pos = TRUE, assay = 'RNA', min.pct = 0.25,
                          logfc.threshold = 0.25, test.use = 'MAST')
write.csv(markers, paste0(path.ct, 'markers_', filename, '_celltype.csv'), row.names = F)
markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  top_n(15, wt = avg_log2FC) -> tops
DefaultAssay(seu) <- 'RNA'
seu <- NormalizeData(seu) %>%
  ScaleData()
DefaultAssay(seu) <- 'integrated'
pdf(file = paste0(path.ct, 'markers_', filename, '_heatmap_celltype.pdf'), width = 10,
    height = 18)
DoHeatmap(seu, features = tops$gene, group.by = 'cell_type_fine', raster = T,
          assay = 'RNA')
dev.off()

# Save object and labels
saveRDS(seu, paste0(path.ct, 'data_', filename, '.rds'))
ct <- seu@meta.data %>%
  select('barcode_pat', 'cell_cycle', 'cell_type_fine')
write.csv(ct, paste0(path.ct, 'data_', filename, '_celltype.csv'))
