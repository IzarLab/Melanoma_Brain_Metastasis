#!/usr/bin/env Rscript

#### Fine cell type annotation of CNS cell subset in sn 
#### Author: Jana Biermann, PhD

library(dplyr)
library(Seurat)


# Read-in object
seu <- readRDS('data/MBPM/sn/data_MBPM_sn_notum.rds')

# Subset to CNS, low-quality and stromal cells
seu <- subset(seu, cell_type_main %in% c('CNS cells', 'Low-quality cells', 'Stromal cells'))
seu <- ScaleData(object = seu)
seu <- RunPCA(object = seu, npcs = 75)
ElbowPlot(seu, ndims = 75)
seu <- FindNeighbors(seu, dims = 1:20)
seu <- FindClusters(seu, resolution = 0.4)
seu <- RunUMAP(object = seu, dims = 1:20)

# Get module scores for sigs
sigs <- read.csv('data/signatures/brain.csv', na.strings = '')
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
write.csv(markers, 'data/cell_type_assignment/sn/cns/markers_MBPM_sn_cns.csv', row.names = F)

# Fine cell type assignment after manual annotation based on DGE
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(5), 'Oligodendrocytes', 'NA')
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(2, 6, 7, 9, 11), 'Low-quality tumor/immune', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(3, 12, 13), 'Neurons', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(8), 'Astrocytes', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(0, 1, 4, 14), 'CAFs', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(10), 'Epithelial cells', seu$cell_type_fine)

# Add stromal cell-type assignment
fb <- read.csv('data/cell_type_assignment/sn/stromal/data_MBPM_sn_stromal_celltype.csv')
seu@meta.data <- left_join(seu@meta.data, fb, by = 'barcode_all')
rownames(seu@meta.data) <- seu$barcode_all
seu$cell_cycle <- ifelse(is.na(seu$cell_cycle.y) == T, seu$cell_cycle.x, seu$cell_cycle.y)
seu$cell_type_fine <- ifelse(is.na(seu$cell_type_fine.y) == T, seu$cell_type_fine.x, seu$cell_type_fine.y)


# Save object and labels
saveRDS(seu, 'data/cell_type_assignment/sn/cns/data_MBPM_sn_cns.rds')
celltype <- seu@meta.data %>%
  select('barcode_all', 'cell_cycle', 'cell_type_fine')
write.csv(celltype, 'data/cell_type_assignment/sn/cns/data_MBPM_sn_cns_celltype.csv')

# DGE based on cell type
Idents(seu) <- seu$cell_type_fine
markers <- FindAllMarkers(seu, only.pos = TRUE, assay = 'RNA', min.pct = 0.25, logfc.threshold = 0.25)
ifelse(!dir.exists(file.path('data/cell_type_DEG/sn/cns')), dir.create(file.path('data/cell_type_DEG/sn/cns'), recursive = T), FALSE)
write.csv(markers, 'data/cell_type_DEG/sn/cns/markers_MBPM_sn_cns_type.csv', row.names = F)
markers$cluster <- as.character(markers$cluster)
tops <- markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(-20, wt = p_val_adj) %>%
  group_by(cluster) %>% slice_head(n = 20)

DefaultAssay(seu) <- 'RNA'
seu <- NormalizeData(seu) %>% ScaleData()
DefaultAssay(seu) <- 'integrated'

pdf(file = 'data/cell_type_DEG/sn/cns/plots_MBPM_sn_cns_heatmap_type.pdf', width = 10, height = 17)
DoHeatmap(seu, features = tops$gene, group.by = 'cell_type_fine', raster = T, assay = 'RNA')
dev.off()
