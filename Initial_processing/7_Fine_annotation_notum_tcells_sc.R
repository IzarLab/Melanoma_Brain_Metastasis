#!/usr/bin/env Rscript

#### Fine cell type annotation of T cell subset in sc 
#### Author: Jana Biermann, PhD

library(dplyr)
library(Seurat)


# Read-in object
seu <- readRDS('data/MBPM/sc/data_MBPM_sc_notum.rds')

# Subset to T cells
seu <- subset(seu, cell_type_main == 'T/NK cells')
seu <- ScaleData(object = seu)
seu <- RunPCA(object = seu, npcs = 75)
ElbowPlot(seu, ndims = 75)
seu <- FindNeighbors(seu, dims = 1:50)
seu <- FindClusters(seu, resolution = 0.6)
seu <- RunUMAP(object = seu, dims = 1:50)

# Get module scores for sigs
sigs <- read.csv('data/signatures/Tcellsignatures.csv', na.strings = '')
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
write.csv(markers, 'data/cell_type_assignment/sc/tcells/markers_MBPM_sc_tcells.csv', row.names = F)

# Fine cell type assignment after manual annotation based on DGE
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.5 %in% c(2, 4, 7), 'CD4+ T cells', 'NA')
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.5 %in% c(7), 'Tfh-like cells', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.5 %in% c(0), 'CD8+ T cells TOX+', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.5 %in% c(1), 'CD8+ T cells TCF7+', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.5 %in% c(3), 'Tregs', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.5 %in% c(5), 'Cycling cells', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.5 %in% c(6), 'NK cells', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.5 %in% c(8), 'Myeloid doublets', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.5 %in% c(9), 'Contamination', seu$cell_type_fine)


# Save object and labels
saveRDS(seu, 'data/cell_type_assignment/sc/tcells/data_MBPM_sc_tcells.rds')
celltype <- seu@meta.data %>%
  select('barcode_all', 'cell_cycle', 'cell_type_fine')
write.csv(celltype, 'data/cell_type_assignment/sc/tcells/data_MBPM_sc_tcells_celltype.csv')

# DGE based on cell type
Idents(seu) <- seu$cell_type_fine
markers <- FindAllMarkers(seu, only.pos = TRUE, assay = 'RNA', min.pct = 0.25, logfc.threshold = 0.25)
ifelse(!dir.exists(file.path('data/cell_type_DEG/sc/tcells')), dir.create(file.path('data/cell_type_DEG/sc/tcells'), recursive = T), FALSE)
write.csv(markers, 'data/cell_type_DEG/sc/tcells/markers_MBPM_sc_tcells_type.csv', row.names = F)
tops <- markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(-20, wt = p_val_adj) %>%
  group_by(cluster) %>% slice_head(n = 20)

DefaultAssay(seu) <- 'RNA'
seu <- NormalizeData(seu) %>% ScaleData()
DefaultAssay(seu) <- 'integrated'

pdf(file = 'data/cell_type_DEG/sc/tcells/plots_MBPM_sc_tcells_heatmap_type.pdf', width = 10, height = 17)
DoHeatmap(seu, features = tops$gene, group.by = 'cell_type_fine', raster = T, assay = 'RNA')
dev.off()
