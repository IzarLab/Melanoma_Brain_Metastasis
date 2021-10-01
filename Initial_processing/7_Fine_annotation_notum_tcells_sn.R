#!/usr/bin/env Rscript

#### Fine cell type annotation of T cell subset in sn 
#### Author: Jana Biermann, PhD

library(dplyr)
library(Seurat)


# Read-in object
seu <- readRDS('data/MBPM/sn/data_MBPM_sn_notum.rds')

# Subset to T cells
seu <- subset(seu, cell_type_main == 'T/NK cells')
DefaultAssay(seu) <- 'RNA'

# Split for re-integration
obj.list <- SplitObject(seu, split.by = 'patient')
# Keep only patient with >50 cells
for (i in 1:length(obj.list)) {
  if (dim(obj.list[[i]]@assays$RNA@counts)[2] < 50) {
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
seu <- IntegrateData(anchorset = anchors, normalization.method = 'LogNormalize', k.weight = 50)
rm(anchors, obj.list)

# Seurat workflow
seu <- ScaleData(object = seu)
seu <- RunPCA(object = seu, npcs = 75)
ElbowPlot(seu, ndims = 75)
seu <- FindNeighbors(seu, dims = 1:20)
seu <- FindClusters(seu, resolution = 0.4)
seu <- RunUMAP(object = seu, dims = 1:20)

# Get module scores for signatures (include CD8_TCF7_sc_top20exclusive and Tfh_sc_top20exclusive markers from sc)
sigs <- read.csv('data/signatures/Tcellsignatures.csv', na.strings = '')
for (c in 1:ncol(sigs)) {
  seu <- AddModuleScore(object = seu, features = list(na.omit(sigs[, c])), name = colnames(sigs)[c], assay = 'RNA', search = F)
}

# Cell cycle
DefaultAssay(seu) <- 'RNA'
seu <- CellCycleScoring(seu, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
DefaultAssay(seu) <- 'integrated'
seu$cell_cycle <- ifelse(seu$G2M.Score > 0.1 | seu$S.Score > 0.1, 'cycling', 'non-cycling')

# Differential gene expression (DGE) based on clusters
markers <- FindAllMarkers(seu, only.pos = TRUE, assay = 'RNA', min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, 'data/cell_type_assignment/sn/tcells/markers_MBPM_sn_tcells.csv', row.names = F)

# Fine cell type assignment after manual annotation based on DGE
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(1, 8) & seu$CD8_TCF7_sc_top20exclusive1 > 0.4, 'CD8+ T cells TCF7+', 'CD4+ T cells')
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(0, 2), 'CD8+ T cells TOX+', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(6, 9), 'Tfh-like cells', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(0, 2, 1) & seu$Tfh_sc_top20exclusive1 > 0.25, 'Tfh-like cells', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(4), 'CD4+ T cells', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(3), 'Tregs', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(5), 'Cycling cells', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(7), 'NK cells', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$integrated_snn_res.0.4 %in% c(10, 11), 'Contamination', seu$cell_type_fine)


# Save object and labels
saveRDS(seu, 'data/cell_type_assignment/sn/tcells/data_MBPM_sn_tcells.rds')
celltype <- seu@meta.data %>%
  select('barcode_all', 'cell_cycle', 'cell_type_fine')
write.csv(celltype, 'data/cell_type_assignment/sn/tcells/data_MBPM_sn_tcells_celltype.csv')

# DGE based on cell type
Idents(seu) <- seu$cell_type_fine
markers <- FindAllMarkers(seu, only.pos = TRUE, assay = 'RNA', min.pct = 0.25, logfc.threshold = 0.25)
ifelse(!dir.exists(file.path('data/cell_type_DEG/sn/tcells')), dir.create(file.path('data/cell_type_DEG/sn/tcells'), recursive = T), FALSE)
write.csv(markers, 'data/cell_type_DEG/sn/tcells/markers_MBPM_sn_tcells_type.csv', row.names = F)
tops <- markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% top_n(-20, wt = p_val_adj) %>%
  group_by(cluster) %>% slice_head(n = 20)

DefaultAssay(seu) <- 'RNA'
seu <- NormalizeData(seu) %>% ScaleData()
DefaultAssay(seu) <- 'integrated'

pdf(file = 'data/cell_type_DEG/sn/tcells/plots_MBPM_sn_tcells_heatmap_type.pdf', width = 10, height = 17)
DoHeatmap(seu, features = tops$gene, group.by = 'cell_type_fine', raster = T, assay = 'RNA')
dev.off()
