#!/usr/bin/env Rscript

#### Fine cell type annotation of myeloid cell subset in sc 
#### Author: Jana Biermann, PhD

library(dplyr)
library(Seurat)


path.ct <- 'data/cell_type_assignment/sc/myeloid/'
filename <- 'MBPM_sc_myeloid'

# Read-in object
seu <- readRDS('data/MBPM/sc/data_MBPM_sc_notum.rds')

# Subset to myeloid cells
seu <- subset(seu, cell_type_main == 'Myeloid cells')
seu <- ScaleData(object = seu)
seu <- RunPCA(object = seu, npcs = 75)
ElbowPlot(seu, ndims = 75)
seu <- FindNeighbors(seu, dims = 1:40)
seu <- FindClusters(seu, resolution = 0.4)
seu <- RunUMAP(object = seu, dims = 1:40)

# Get module scores for sigs
sigs <- read.csv('data/signatures/myeloid_signatures.csv', na.strings = '')
for (c in 1:ncol(sigs)) {
  seu <- AddModuleScore(object = seu, features = list(na.omit(sigs[, c])), name = colnames(sigs)[c], assay = 'RNA', search = F)
}

# Cell cycle
DefaultAssay(seu) <- 'RNA'
seu <- CellCycleScoring(seu, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
DefaultAssay(seu) <- 'integrated'
seu$cell_cycle <- ifelse(seu$G2M.Score > 0.12 | seu$S.Score > 0.12, 'cycling', 'non-cycling')

# Differential gene expression (DGE) based on clusters
markers <- FindAllMarkers(seu, only.pos = TRUE, assay = 'RNA', min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, paste0(path.ct, 'markers_', filename, '.csv'), row.names = F)

# Fine cell type assignment after manual annotation based on DGE
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(0, 1, 3, 6, 7), 'Monocyte-derived macrophages', 'NA')
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(2), 'Microglia', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(9, 10, 11, 12), 'DCs', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(11), 'cDC1', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(9), 'cDC2', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(10), 'mo-DCs', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(12), 'DC3', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(5), 'Monocytes', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(4), 'T-cell doublets', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(8), 'Cycling cells', seu$cell_type_fine)


#### Subset to MDMs for separate annotation ####
mdm <- subset(seu, integrated_snn_res.0.4 %in% c(0, 1, 3, 6, 7))
mdm$original_clusters <- mdm$integrated_snn_res.0.4
mdm <- ScaleData(mdm)
mdm <- RunPCA(mdm, npcs = 75)
ElbowPlot(mdm, ndims = 75)
mdm <- FindNeighbors(mdm, dims = 1:40)
mdm <- FindClusters(mdm, resolution = 0.3)
mdm <- RunUMAP(mdm, dims = 1:40)

# Get module scores for signatures
sigs <- read.csv('data/signatures/myeloid_signatures.csv', na.strings = '')
for (c in 1:ncol(sigs)) {
  mdm <- AddModuleScore(object = mdm, features = list(na.omit(sigs[, c])), name = colnames(sigs)[c], assay = 'RNA', search = F)
}

# Differential gene expression (DGE) based on MDM clusters
Idents(mdm) <- mdm$integrated_snn_res.0.3
markers <- FindAllMarkers(mdm, only.pos = TRUE, assay = 'RNA', min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, paste0(path.ct, 'MDM/markers_', filename, '_MDM.csv'), row.names = F)

# Fine cell type assignment after manual annotation based on DGE for MDM
mdm$cell_type_fine <- ifelse(mdm$integrated_snn_res.0.3 %in% c(0, 2, 4), 'TAM MDM', 'NA')
mdm$cell_type_fine <- ifelse(mdm$integrated_snn_res.0.3 %in% c(1), 'TAM MDM M2-like', mdm$cell_type_fine)
mdm$cell_type_fine <- ifelse(mdm$integrated_snn_res.0.3 %in% c(5), 'TAM MDM M1-like', mdm$cell_type_fine)
mdm$cell_type_fine <- ifelse(mdm$integrated_snn_res.0.3 %in% c(0) & mdm$Monocytes1 > 0.4, 'Monocytes', mdm$cell_type_fine)
mdm$cell_type_fine <- ifelse(mdm$integrated_snn_res.0.3 %in% c(3), 'TAM MDM FTL+', mdm$cell_type_fine)


# Save intermediate object and labels
saveRDS(mdm, paste0(path.ct, 'MDM/data_', filename, '_MDM.rds'))
celltype <- mdm@meta.data %>% select('barcode_all', 'cell_cycle', 'cell_type_fine')
write.csv(celltype, paste0(path.ct, 'MDM/data_', filename, '_MDM_celltype.csv'))


#### Combine all annotaions for final myeloid object #### 
# Add annotation from MDM
celltype <- read.csv(paste0(path.ct, 'MDM/data_', filename, '_MDM_celltype.csv'))
seu@meta.data <- left_join(seu@meta.data, celltype, by = 'barcode_all')
rownames(seu@meta.data) <- seu$barcode_all

# Save object and labels
saveRDS(seu, paste0(path.ct, 'data_', filename, '.rds'))
celltype <- seu@meta.data %>% select('barcode_all', 'cell_cycle', 'cell_type_fine')
write.csv(celltype, paste0(path.ct, 'data_', filename, '_celltype.csv'))

# DGE based on cell type
Idents(seu) <- seu$cell_type_fine
markers <- FindAllMarkers(seu, only.pos = TRUE, assay = 'RNA', min.pct = 0.25, logfc.threshold = 0.25)
ifelse(!dir.exists(file.path('data/cell_type_DEG/sc/myeloid')), dir.create(file.path('data/cell_type_DEG/sc/myeloid'), recursive = T), FALSE)
write.csv(markers, paste0('data/cell_type_DEG/sc/myeloid/markers_', filename, '_types.csv'), row.names = F)
tops <- markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(-20, wt = p_val_adj) %>%
  group_by(cluster) %>% slice_head(n = 20)

DefaultAssay(seu) <- 'RNA'
seu <- NormalizeData(seu) %>% ScaleData()
DefaultAssay(seu) <- 'integrated'

pdf(file = paste0('data/cell_type_DEG/sc/myeloid/markers_', filename, '_heatmap_type.pdf'), width = 10, height = 18)
DoHeatmap(seu, features = tops$gene, group.by = 'cell_type_fine', raster = T, assay = 'RNA')
dev.off()
