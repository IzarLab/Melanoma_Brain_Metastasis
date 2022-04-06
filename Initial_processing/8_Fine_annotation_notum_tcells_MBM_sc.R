#!/usr/bin/env Rscript

#### Fine cell type annotation of T cell subset for MBM_sc
#### Author: Jana Biermann, PhD

library(dplyr)
library(Seurat)
library(ggplot2)
library(gplots)
library(viridis)

path.ct <- 'data/cell_type_DEG/MBM_sc/tcells/'
filename <- 'MBM_sc_tcells'

seu <- readRDS('data/cell_type_DEG/MBM_sc/main/data_MBM_sc_notum.rds')

# Subset to T cells
seu <- subset(seu, cell_type_main == 'T/NK cells')
seu <- ScaleData(object = seu)
seu <- RunPCA(object = seu, npcs = 75)
ElbowPlot(seu, ndims = 75)
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.4)
seu <- RunUMAP(object = seu, dims = 1:30)
DimPlot(seu, reduction = 'umap', label = T, group.by = 'ident', shuffle = T, raster = T)

# Get module scores for signatures
sigs <- read.csv('~/signatures/Tcellsignatures.csv', na.strings = '')
for (c in 1:ncol(sigs)) {
  seu <- AddModuleScore(object = seu, features = list(na.omit(sigs[, c])),
                        name = colnames(sigs)[c], assay = 'RNA', search = F)
}

# Cell cycle assignment
DefaultAssay(seu) <- 'RNA'
seu <- CellCycleScoring(seu, s.features = cc.genes.updated.2019$s.genes,
                        g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
DefaultAssay(seu) <- 'integrated'
seu$cell_cycle <- ifelse(seu$G2M.Score > 0.12 | seu$S.Score > 0.12, 'cycling',
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
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(0), 'CD4+ T cells', 'NA')
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(6), 'Tfh-like cells',
                             seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(2, 7), 'CD8+ T cells TOX+',
                             seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(1), 'CD8+ T cells TCF7+',
                             seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(3), 'Tregs',
                             seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(4), 'Cycling cells',
                             seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(5), 'NK cells',
                             seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(8, 9), 'Myeloid doublets',
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
