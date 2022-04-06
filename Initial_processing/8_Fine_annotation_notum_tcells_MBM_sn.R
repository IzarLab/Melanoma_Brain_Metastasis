#!/usr/bin/env Rscript

#### Fine cell type annotation of T cell subset for MBM_sn
#### Author: Jana Biermann, PhD

library(dplyr)
library(Seurat)
library(ggplot2)
library(gplots)
library(viridis)
'%notin%' <- Negate('%in%')

path.ct <- 'data/cell_type_DEG/MBM_sn/tcells/'
filename <- 'MBM_sn_tcells'

seu <- readRDS('data/cell_type_DEG/MBM_sn/main/data_MBM_sn_notum.rds')

# Subset to T cells
seu <- subset(seu, cell_type_main %in% c('T/NK cells'))

## Reintegrate
DefaultAssay(seu) <- 'RNA'

# Split 
obj.list <- SplitObject(seu, split.by = 'patient')
# Keep only patient with >100 cells
for (i in 1:length(obj.list)) {
  if (dim(obj.list[[i]]@assays$RNA@counts)[2] < 100) {
    obj.list[i] <- NA
  } else {
    obj.list[[i]] <- NormalizeData(obj.list[[i]])
    obj.list[[i]] <- FindVariableFeatures(obj.list[[i]])
  }
}
# Remove NAs
obj.list <- obj.list[!is.na(obj.list)]

# Integration
anchors <- FindIntegrationAnchors(object.list = obj.list)
seu <- IntegrateData(anchorset = anchors)

# Rerun workflow
seu <- ScaleData(object = seu)
seu <- RunPCA(object = seu, npcs = 75)
ElbowPlot(seu, ndims = 75)
seu <- FindNeighbors(seu, dims = 1:20)
seu <- FindClusters(seu, resolution = 0.4)
seu <- RunUMAP(object = seu, dims = 1:20)
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
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(1), 'CD8+ T cells TCF7+',
                             'NA')
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(0, 5, 9, 11),
                             'CD8+ T cells TOX+', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(3), 'Tfh-like cells',
                             seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(2), 'CD4+ T cells',
                             seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(4, 8, 10), 'Tregs',
                             seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(7), 'NK cells',
                             seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(6, 12), 'Contamination',
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
