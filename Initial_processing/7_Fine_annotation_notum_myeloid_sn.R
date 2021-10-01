#!/usr/bin/env Rscript

#### Fine cell type annotation of myeloid cell subset in sn 
#### Author: Jana Biermann, PhD

library(dplyr)
library(Seurat)


path.ct <- 'data/cell_type_assignment/sn/myeloid/'
filename <- 'MBPM_sn_myeloid'

# Read-in object
seu <- readRDS('data/MBPM/sn/data_MBPM_sn_notum.rds')

# Subset to myeloid cells
seu <- subset(seu, cell_type_main == 'Myeloid cells')
DefaultAssay(seu) <- 'RNA'

# Split for re-integration
obj.list <- SplitObject(seu, split.by = 'patient')
# Keep only samples with >100 cells
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

# Re-integration
anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = 'LogNormalize')
seu <- IntegrateData(anchorset = anchors, normalization.method = 'LogNormalize')
rm(anchors, obj.list)

# Seurat workflow
seu <- ScaleData(object = seu)
seu <- RunPCA(object = seu, npcs = 75)
ElbowPlot(seu, ndims = 75)
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.4)
seu <- RunUMAP(object = seu, dims = 1:30)

# Get module scores for signatures
sigs <- read.csv('data/signatures/myeloid_signatures.csv', na.strings = '')
for (c in 1:ncol(sigs)) {
  seu <- AddModuleScore(object = seu, features = list(na.omit(sigs[, c])), name = colnames(sigs)[c], assay = 'RNA', search = F)
}

# Cell cycle
DefaultAssay(seu) <- 'RNA'
seu <- CellCycleScoring(seu, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
DefaultAssay(seu) <- 'integrated'
seu$cell_cycle <- ifelse(seu$G2M.Score > 0.08 | seu$S.Score > 0.08, 'cycling', 'non-cycling')

# Differential gene expression (DGE) based on clusters
markers <- FindAllMarkers(seu, only.pos = TRUE, assay = 'RNA', min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, 'data/cell_type_assignment/sn/myeloid/markers_MBPM_sn_myeloid.csv', row.names = F)

# Fine cell type assignment after manual annotation based on DGE
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(0, 1, 2, 4, 7), 'Monocyte-derived macrophages', 'NA')
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(5), 'Microglia', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(6), 'DCs', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(3), 'Monocytes', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(9), 'Cycling cells', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(8), 'Contamination', seu$cell_type_fine)

# Save intermidate labels
celltype <- seu@meta.data %>% dplyr::select('barcode_all', 'cell_cycle', 'cell_type_fine')
write.csv(celltype, 'data/cell_type_assignment/sn/myeloid/data_MBPM_sn_myeloid_celltype.csv')


#### Subset to MDM and monocytes and annotate separately #### 
### pt1
mm <- subset(seu, cell_type_fine %in% c('Monocyte-derived macrophages', 'Monocytes'))

# Rerun Seurat workflow
mm <- ScaleData(mm)
mm <- RunPCA(mm)
ElbowPlot(mm, ndims = 75)
mm <- RunUMAP(mm, dims = 1:20)
mm <- FindNeighbors(mm, dims = 1:20)
mm <- FindClusters(mm, resolution = 0.2)

# Differential gene expression (DGE) based on clusters
ifelse(!dir.exists(file.path(path.ct)), dir.create(file.path(path.ct), recursive = T), FALSE)
markers <- FindAllMarkers(mm, only.pos = TRUE, assay = 'RNA', min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, paste0(path.ct, 'markers_', filename, '.csv'), row.names = F)

# Fine cell type assignment after manual annotation based on DGE
mm$cell_type_fine <- ifelse(mm$integrated_snn_res.0.2 %in% c(3), 'Patient-specific MDM', mm$cell_type_fine)

# Save labels pt1
celltype <- mm@meta.data %>% dplyr::select('barcode_all', 'cell_cycle', 'cell_type_fine')
write.csv(celltype, paste0(path.ct, 'data_', filename, '_celltype_pt1.csv'))


### pt2
mm <- subset(mm, cell_type_fine %in% c('Monocyte-derived macrophages', 'Monocytes'))

# Rerun Seurat workflow
mm <- ScaleData(mm)
mm <- RunPCA(mm)
ElbowPlot(mm, ndims = 75)
mm <- RunUMAP(mm, dims = 1:20)
mm <- FindNeighbors(mm, dims = 1:20)
mm <- FindClusters(mm, resolution = 0.2)

# Fine cell type assignment after manual annotation based on DGE
mm$cell_type_fine <- ifelse(mm$integrated_snn_res.0.2 %in% c(0, 2), 'TAM MDM', mm$cell_type_fine)
mm$cell_type_fine <- ifelse(mm$integrated_snn_res.0.2 %in% c(1), 'TAM MDM FTL+', mm$cell_type_fine)
mm$cell_type_fine <- ifelse(mm$integrated_snn_res.0.2 %in% c(3), 'Monocytes', mm$cell_type_fine)


# Save labels pt2
celltype <- mm@meta.data %>%
  dplyr::select('barcode_all', 'cell_cycle', 'cell_type_fine')
write.csv(celltype, paste0(path.ct, 'data_', filename, '_celltype_pt2.csv'))


#### Subset to DCs and annotate separately ####
dc <- subset(seu, cell_type_fine %in% c('DCs'))
dc <- ScaleData(dc)
dc <- RunPCA(dc, npcs = 75)
ElbowPlot(dc, ndims = 75)
dc <- FindNeighbors(dc, dims = 1:40)
dc <- FindClusters(dc, resolution = 0.5)
dc <- RunUMAP(dc, dims = 1:40)

# Get module scores for signatures
sigs <- read.csv('data/signatures/myeloid_signatures.csv', na.strings = '')
for (c in 1:ncol(sigs)) {
  dc <- AddModuleScore(object = dc, features = list(na.omit(sigs[, c])), name = colnames(sigs)[c], assay = 'RNA', search = F)
}

# Differential gene expression (DGE) based on clusters
markers <- FindAllMarkers(dc, only.pos = TRUE, assay = 'RNA', min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, paste0(path.ct, 'markers_', filename, '_DC.csv'), row.names = F)

# Fine cell type assignment after manual annotation based on DGE
dc$cell_type_fine <- ifelse(dc$integrated_snn_res.0.5 %in% c(3), 'cDC1', 'NA')
dc$cell_type_fine <- ifelse(dc$integrated_snn_res.0.5 %in% c(1), 'cDC2', dc$cell_type_fine)
dc$cell_type_fine <- ifelse(dc$integrated_snn_res.0.5 %in% c(2), 'DC3', dc$cell_type_fine)
dc$cell_type_fine <- ifelse(dc$integrated_snn_res.0.5 %in% c(0), 'T-cell doublets', dc$cell_type_fine)

# Save intermidate labels
celltype <- dc@meta.data %>% dplyr::select('barcode_all', 'cell_cycle', 'cell_type_fine')
write.csv(celltype, paste0(path.ct, 'data_', filename, '_DC_celltype.csv'))


#### Combine all annotaions for final myeloid object #### 
# Add annotaions
mdm1 <- read.csv(paste0(path.ct, 'data_', filename, '_celltype_pt1.csv'))
mdm2 <- read.csv(paste0(path.ct, 'data_', filename, '_celltype_pt2.csv'))
mdm <- left_join(mdm1, mdm2, by = 'barcode_all')
mdm$cell_cycle <- mdm$cell_cycle.x
mdm$cell_type_fine <- ifelse(is.na(mdm$cell_type_fine.y) == T, mdm$cell_type_fine.x, mdm$cell_type_fine.y)
mdm <- mdm[, c('barcode_all', 'cell_cycle', 'cell_type_fine')]

ct <- read.csv('data/cell_type_assignment/sn/myeloid/data_MBPM_sn_myeloid_celltype.csv')
ct <- left_join(ct, mdm, by = 'barcode_all')
ct$cell_cycle <- ifelse(is.na(ct$cell_cycle.y) == T, ct$cell_cycle.x, ct$cell_cycle.y)
ct$cell_type_fine <- ifelse(is.na(ct$cell_type_fine.y) == T, ct$cell_type_fine.x, ct$cell_type_fine.y)
ct <- ct[, c('barcode_all', 'cell_cycle', 'cell_type_fine')]

dc <- read.csv(paste0(path.ct, 'data_', filename, '_DC_celltype.csv'))
ct <- left_join(ct, dc, by = 'barcode_all')
ct$cell_cycle <- ifelse(is.na(ct$cell_cycle.y) == T, ct$cell_cycle.x, ct$cell_cycle.y)
ct$cell_type_fine <- ifelse(is.na(ct$cell_type_fine.y) == T, ct$cell_type_fine.x, ct$cell_type_fine.y)
ct <- ct[, c('barcode_all', 'cell_cycle', 'cell_type_fine')]

seu@meta.data <- left_join(seu@meta.data, ct, by = 'barcode_all')
rownames(seu@meta.data) <- seu$barcode_all
seu$cell_cycle <- ifelse(is.na(seu$cell_cycle.y) == T, seu$cell_cycle.x, seu$cell_cycle.y)
seu$cell_type_fine <- ifelse(is.na(seu$cell_type_fine.y) == T, seu$cell_type_fine.x, seu$cell_type_fine.y)


# Save object and labels
saveRDS(seu, 'data/cell_type_assignment/sn/myeloid/data_MBPM_sn_myeloid.rds')
celltype <- seu@meta.data %>% dplyr::select('barcode_all', 'cell_cycle', 'cell_type_fine')
write.csv(celltype, 'data/cell_type_assignment/sn/myeloid/data_MBPM_sn_myeloid_celltype.csv')

# DGE based on cell type
Idents(seu) <- seu$cell_type_fine
markers <- FindAllMarkers(seu, only.pos = TRUE, assay = 'RNA', min.pct = 0.25, logfc.threshold = 0.25)
ifelse(!dir.exists(file.path('data/cell_type_DEG/sn/myeloid')), dir.create(file.path('data/cell_type_DEG/sn/myeloid'), recursive = T), FALSE)
write.csv(markers, 'data/cell_type_DEG/sn/myeloid/markers_MBPM_sn_myeloid_type.csv', row.names = F)
markers$cluster <- as.character(markers$cluster)
tops <- markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(-20, wt = p_val_adj) %>%
  group_by(cluster) %>% slice_head(n = 20)

DefaultAssay(seu) <- 'RNA'
seu <- NormalizeData(seu) %>% ScaleData()
DefaultAssay(seu) <- 'integrated'

pdf(file = 'data/cell_type_DEG/sn/myeloid/plots_MBPM_sn_myeloid_heatmap_type.pdf', width = 15, height = 18)
DoHeatmap(seu, features = tops$gene, group.by = 'cell_type_fine', raster = T, assay = 'RNA')
dev.off()
