#!/usr/bin/env Rscript

#### Fine cell type annotation of stromal subset in sn 
#### Author: Jana Biermann, PhD

library(dplyr)
library(Seurat)


# Read-in object
seu <- readRDS('data/MBPM/sn/data_MBPM_sn_notum.rds')

# Subset to stromal cells
seu <- subset(seu, cell_type_main == 'Stromal cells')
seu <- ScaleData(object = seu)
seu <- RunPCA(object = seu, npcs = 75)
ElbowPlot(seu, ndims = 75)
seu <- FindNeighbors(seu, dims = 1:40)
seu <- FindClusters(seu, resolution = 0.15)
seu <- RunUMAP(object = seu, dims = 1:40)

# Cell cycle
DefaultAssay(seu) <- 'RNA'
seu <- CellCycleScoring(seu, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
DefaultAssay(seu) <- 'integrated'
seu$cell_cycle <- ifelse(seu$G2M.Score > 0.08 | seu$S.Score > 0.08, 'cycling', 'non-cycling')

# Differential gene expression (DGE) based on clusters
markers <- FindAllMarkers(seu, only.pos = TRUE, assay = 'RNA', min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, 'data/cell_type_assignment/sn/stromal/markers_MBPM_sn_stromal.csv', row.names = F)

# Fine cell type assignment after manual annotation based on DGE
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.15 %in% c(0, 2), 'CAFs', 'NA')
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.15 %in% c(1), 'Pericytes', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.15 %in% c(3), 'Cycling CAFs', seu$cell_type_fine)

# Save object and labels
saveRDS(seu, 'data/cell_type_assignment/sn/stromal/data_MBPM_sn_stromal.rds')
celltype <- seu@meta.data %>% dplyr::select('barcode_all', 'cell_cycle', 'cell_type_fine')
write.csv(celltype, 'data/cell_type_assignment/sn/stromal/data_MBPM_sn_stromal_celltype.csv')

# DGE based on cell type
Idents(seu) <- seu$cell_type_fine
markers <- FindAllMarkers(seu, only.pos = TRUE, assay = 'RNA', min.pct = 0.25, logfc.threshold = 0.25)
ifelse(!dir.exists(file.path('data/cell_type_DEG/sn/stromal')), dir.create(file.path('data/cell_type_DEG/sn/stromal'), recursive = T), FALSE)
write.csv(markers, 'data/cell_type_DEG/sn/stromal/markers_MBPM_sn_stromal_type.csv', row.names = F)
tops <- markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(-20, wt = p_val_adj) %>%
  group_by(cluster) %>% slice_head(n = 20)

DefaultAssay(seu) <- 'RNA'
seu <- NormalizeData(seu) %>% ScaleData()
DefaultAssay(seu) <- 'integrated'

pdf(file = 'data/cell_type_DEG/sn/stromal/plots_MBPM_sn_stromal_heatmap_type.pdf', width = 10, height = 17)
DoHeatmap(seu, features = tops$gene, group.by = 'cell_type_fine', raster = T, assay = 'RNA')
dev.off()
