#!/usr/bin/env Rscript

#### Fine cell type annotation of B cell subset in sc 
####Author: Jana Biermann, PhD

library(dplyr)
library(Seurat)


# Read-in object
seu <- readRDS('data/MBPM/sc/data_MBPM_sc_notum.rds')

# Subset to B cells
seu <- subset(seu, cell_type_main == 'B/Plasma cells')
seu <- ScaleData(object = seu)
seu <- RunPCA(object = seu, npcs = 75)
ElbowPlot(seu, ndims = 75)
seu <- FindNeighbors(seu, dims = 1:40)
seu <- FindClusters(seu, resolution = 0.2)
seu <- RunUMAP(object = seu, dims = 1:40)

# Cell cycle
DefaultAssay(seu) <- 'RNA'
seu <- CellCycleScoring(seu, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
DefaultAssay(seu) <- 'integrated'
seu$cell_cycle <- ifelse(seu$G2M.Score > 0.1 | seu$S.Score > 0.1, 'cycling', 'non-cycling')

# Differential gene expression (DGE) based on clusters
markers <- FindAllMarkers(seu, only.pos = TRUE, assay = 'RNA', min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, 'data/cell_type_assignment/sc/bcells/markers_MBPM_sc_bcells.csv', row.names = F)

# Fine cell type assignment after manual annotation based on DGE
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.2 %in% c(1, 2), 'Activated B cells', 'NA')
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.2 %in% c(3), 'Naïve B cells', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.2 %in% c(0), 'Plasma cells', seu$cell_type_fine)


# Save object and labels
saveRDS(seu, 'data/cell_type_assignment/sc/bcells/data_MBPM_sc_bcells.rds')
celltype <- seu@meta.data %>% select('barcode_all', 'cell_cycle', 'cell_type_fine')
write.csv(celltype, 'data/cell_type_assignment/sc/bcells/data_MBPM_sc_bcells_celltype.csv')

# DGE based on cell type
Idents(seu) <- seu$cell_type_fine
markers <- FindAllMarkers(seu, only.pos = TRUE, assay = 'RNA', min.pct = 0.25, logfc.threshold = 0.25)
ifelse(!dir.exists(file.path('data/cell_type_DEG/sc/bcells')), dir.create(file.path('data/cell_type_DEG/sc/bcells'), recursive = T), FALSE)
write.csv(markers, 'data/cell_type_DEG/sc/bcells/markers_MBPM_sc_bcells_type.csv', row.names = F)
tops <- markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(-20, wt = p_val_adj) %>%
  group_by(cluster) %>% slice_head(n = 20)

DefaultAssay(seu) <- 'RNA'
seu <- NormalizeData(seu) %>% ScaleData()
DefaultAssay(seu) <- 'integrated'

pdf(file = 'data/cell_type_DEG/sc/bcells/plots_MBPM_sc_bcells_heatmap_type.pdf', width = 10, height = 17)
DoHeatmap(seu, features = tops$gene, group.by = 'cell_type_fine', raster = T, assay = 'RNA')
dev.off()
