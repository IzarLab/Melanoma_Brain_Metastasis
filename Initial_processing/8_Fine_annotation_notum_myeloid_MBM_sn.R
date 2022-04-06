#!/usr/bin/env Rscript

#### Fine cell type annotation of myeloid cell subset for MBM_sn
#### Author: Jana Biermann, PhD

library(dplyr)
library(Seurat)
library(ggplot2)
library(gplots)
library(viridis)
'%notin%' <- Negate('%in%')

path.ct <- 'data/cell_type_DEG/MBM_sn/myeloid/'
filename <- 'MBM_sn_myeloid'

seu <- readRDS('data/cell_type_DEG/MBM_sn/main/data_MBM_sn_notum.rds')

# Subset to myeloid cells
seu <- subset(seu, cell_type_main == 'Myeloid cells')
seu <- ScaleData(object = seu)
seu <- RunPCA(object = seu, npcs = 75)
ElbowPlot(seu, ndims = 75)
seu <- FindNeighbors(seu, dims = 1:20)
seu <- FindClusters(seu, resolution = 0.4)
seu <- RunUMAP(object = seu, dims = 1:20)
DimPlot(seu, reduction = 'umap', label = T, group.by = 'ident', shuffle = T, raster = T)

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
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.6)
seu <- RunUMAP(object = seu, dims = 1:30)
DimPlot(seu, reduction = 'umap', label = T, group.by = 'ident', shuffle = T, raster = T)
FeaturePlot(seu, features = c('CD3D', 'CD3E', 'TOX'), min.cutoff = 'q05',
            max.cutoff = 'q95', order = T, raster = T)

seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.6 == 8, 'T-cell doublets', NA)
ct <- seu@meta.data %>%
  select('barcode_pat', 'cell_cycle', 'cell_type_fine')
ct <- subset(ct, is.na(ct$cell_type_fine) == F)
write.csv(ct, paste0(path.ct, 'data_', filename, '_celltype_pt1.csv'), row.names = F)

# Remove T-cell doublets
seu <- subset(seu, integrated_snn_res.0.6 != 8)
seu <- ScaleData(object = seu)
seu <- RunPCA(object = seu, npcs = 75)
ElbowPlot(seu, ndims = 75)
seu <- FindNeighbors(seu, dims = 1:20)
seu <- FindClusters(seu, resolution = 0.3)
seu <- RunUMAP(object = seu, dims = 1:20)
DimPlot(seu, reduction = 'umap', label = T, group.by = 'ident', shuffle = T, raster = T)

# Get module scores for signatures
sigs <- read.csv('~/signatures/myeloid_signatures.csv', na.strings = '')
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

# Fine cell type assignment after manual annotation based on DGE reintegrated
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.3 %in% c(0, 1, 2, 5),
                             'Monocyte-derived macrophages', 'NA')
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.3 %in% c(3), 'Microglia',
                             seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.3 %in% c(4), 'DCs',
                             seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.3 %in% c(6, 7, 8), 'unknown',
                             seu$cell_type_fine)


#### Subset to DCs ####
dc <- subset(seu, integrated_snn_res.0.3 %in% c(4))
dc <- ScaleData(dc)
dc <- RunPCA(dc, npcs = 75)
ElbowPlot(dc, ndims = 75)
dc <- FindNeighbors(dc, dims = 1:20)
dc <- FindClusters(dc, resolution = 0.4)
dc <- RunUMAP(dc, dims = 1:20)
DimPlot(dc, reduction = 'umap', label = T, group.by = 'ident', shuffle = T, raster = T)

# Get module scores for signatures
sigs <- read.csv('~/signatures/myeloid_signatures.csv', na.strings = '')
for (c in 1:ncol(sigs)) {
  dc <- AddModuleScore(object = dc, features = list(na.omit(sigs[, c])),
                       name = colnames(sigs)[c], assay = 'RNA', search = F)
}

# Differential gene expression (DGE) based on clusters
Idents(dc) <- dc$integrated_snn_res.0.4
markers <- FindAllMarkers(dc, only.pos = TRUE, assay = 'RNA', min.pct = 0.25,
                          logfc.threshold = 0.25, test.use = 'MAST')
write.csv(markers, paste0(path.ct, 'markers_', filename, '_DC.csv'), row.names = F)
markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  top_n(20, wt = avg_log2FC) -> tops

DefaultAssay(dc) <- 'RNA'
dc <- NormalizeData(dc) %>%
  ScaleData()
DefaultAssay(dc) <- 'integrated'

pdf(file = paste0(path.ct, 'markers_', filename, '_DC_heatmap.pdf'), width = 10,
    height = 18)
DoHeatmap(dc, features = tops$gene, group.by = 'ident', raster = T, assay = 'RNA')
dev.off()

# Fine cell type assignment after manual annotation based on DGE
dc$cell_type_fine <- ifelse(dc$integrated_snn_res.0.4 %in% c(0), 'Monocytes', 'NA')
dc$cell_type_fine <- ifelse(dc$integrated_snn_res.0.4 %in% c(2), 'cDC1',
                            dc$cell_type_fine)
dc$cell_type_fine <- ifelse(dc$integrated_snn_res.0.4 %in% c(3), 'DC3', dc$cell_type_fine)
dc$cell_type_fine <- ifelse(dc$integrated_snn_res.0.4 %in% c(1), 'cDC2',
                            dc$cell_type_fine)

# Save object and labels
saveRDS(dc, paste0(path.ct, 'data_', filename, '_DC.rds'))
celltype <- dc@meta.data %>%
  select('barcode_pat', 'cell_cycle', 'cell_type_fine')
write.csv(celltype, paste0(path.ct, 'data_', filename, '_DC_celltype.csv'))


#### Add DC annotation ####
ct_dc <- read.csv(paste0(path.ct, 'data_', filename, '_DC_celltype.csv'))
ct_dc$X <- NULL
seu$barcode_all <- rownames(seu@meta.data)
seu@meta.data <- left_join(seu@meta.data, ct_dc, by = 'barcode_pat')
rownames(seu@meta.data) <- seu$barcode_all
seu$cell_type_fine <- ifelse(is.na(seu$cell_type_fine.y) == T, seu$cell_type_fine.x,
                             seu$cell_type_fine.y)
seu$cell_cycle <- seu$cell_cycle.x
seu$cell_cycle.y <- NULL
seu$cell_cycle.x <- NULL
seu$cell_type_fine.x <- NULL
seu$cell_type_fine.y <- NULL


#### Subset to MDMs ####
mdm <- subset(seu, cell_type_fine %in% c('Monocyte-derived macrophages', 'Monocytes',
                                         'unknown'))
mdm$original_clusters <- mdm$integrated_snn_res.0.3
mdm <- ScaleData(mdm)
mdm <- RunPCA(mdm, npcs = 75)
ElbowPlot(mdm, ndims = 75)
mdm <- FindNeighbors(mdm, dims = 1:20)
mdm <- FindClusters(mdm, resolution = 0.7)
mdm <- RunUMAP(mdm, dims = 1:20)
DimPlot(mdm, reduction = 'umap', label = T, group.by = 'ident', shuffle = T, raster = T)

# Get module scores for signatures
sigs <- read.csv('~/signatures/myeloid_signatures.csv', na.strings = '')
for (c in 1:ncol(sigs)) {
  mdm <- AddModuleScore(object = mdm, features = list(na.omit(sigs[, c])),
                        name = colnames(sigs)[c], assay = 'RNA', search = F)
}

# Differential gene expression (DGE) based on clusters
markers <- FindAllMarkers(mdm, only.pos = TRUE, assay = 'RNA', min.pct = 0.25,
                          logfc.threshold = 0.25, test.use = 'MAST')
write.csv(markers, paste0(path.ct, 'markers_', filename, '_MDM.csv'), row.names = F)
markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  top_n(-20, wt = p_val_adj) -> tops

DefaultAssay(mdm) <- 'RNA'
mdm <- NormalizeData(mdm) %>%
  ScaleData()
DefaultAssay(mdm) <- 'integrated'

pdf(file = paste0(path.ct, 'markers_', filename, '_MDM_heatmap.pdf'), width = 10,
    height = 18)
DoHeatmap(mdm, features = tops$gene, group.by = 'ident', raster = T, assay = 'RNA')
dev.off()

# Fine cell type assignment after manual annotation based on DGE
mdm$cell_type_fine <- ifelse(mdm$integrated_snn_res.0.7 %in% c(2, 3, 5, 9), 'MDM', 'NA')
mdm$cell_type_fine <- ifelse(mdm$integrated_snn_res.0.7 %in% c(0), 'MDM M2-like',
                             mdm$cell_type_fine)
mdm$cell_type_fine <- ifelse(mdm$integrated_snn_res.0.7 %in% c(7), 'MDM M1-like',
                             mdm$cell_type_fine)
mdm$cell_type_fine <- ifelse(mdm$integrated_snn_res.0.7 %in% c(4), 'Proinflammatory MDM',
                             mdm$cell_type_fine)
mdm$cell_type_fine <- ifelse(mdm$integrated_snn_res.0.7 %in% c(1, 6), 'MDM FTL+',
                             mdm$cell_type_fine)
mdm$cell_type_fine <- ifelse(mdm$integrated_snn_res.0.7 %in% c(8), 'Monocytes',
                             mdm$cell_type_fine)
mdm$cell_type_fine <- ifelse(mdm$integrated_snn_res.0.7 %in% c(10, 11, 12, 13, 14),
                             'Contamination', mdm$cell_type_fine)

# Save object and labels
saveRDS(mdm, paste0(path.ct, 'data_', filename, '_MDM.rds'))
celltype <- mdm@meta.data %>%
  select('barcode_pat', 'cell_cycle', 'cell_type_fine')
write.csv(celltype, paste0(path.ct, 'data_', filename, '_MDM_celltype.csv'))


#### Add annotation from MDM ####
ct_mdm <- read.csv(paste0(path.ct, 'data_', filename, '_MDM_celltype.csv'))
ct_mdm$X <- NULL
seu$barcode_all <- rownames(seu@meta.data)
seu@meta.data <- left_join(seu@meta.data, ct_mdm, by = 'barcode_pat')
rownames(seu@meta.data) <- seu$barcode_all
seu$cell_type_fine <- ifelse(is.na(seu$cell_type_fine.y) == T, seu$cell_type_fine.x,
                             seu$cell_type_fine.y)
seu$cell_cycle <- seu$cell_cycle.x
seu$cell_cycle.y <- NULL
seu$cell_cycle.x <- NULL
seu$cell_type_fine.x <- NULL
seu$cell_type_fine.y <- NULL

# Differential gene expression (DGE) based on cell type
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
DoHeatmap(seu, features = tops$gene, group.by = 'ident', raster = T, assay = 'RNA')
dev.off()

# Save object and labels
saveRDS(seu, paste0(path.ct, 'data_', filename, '.rds'))
ct <- seu@meta.data %>%
  select('barcode_pat', 'cell_cycle', 'cell_type_fine')
write.csv(ct, paste0(path.ct, 'data_', filename, '_celltype.csv'))
