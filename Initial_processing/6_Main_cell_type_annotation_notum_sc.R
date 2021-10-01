#!/usr/bin/env Rscript

#### Main cell type annotation of integrated Seurat object in sc 
#### Author: Jana Biermann, PhD

library(dplyr)
library(Seurat)


# Read-in object
seu <- readRDS('data/MBPM/sc/data_MBPM_sc_notum.rds')

# Differential gene expression (DGE) based on clusters
markers <- FindAllMarkers(seu, assay = 'RNA', only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, 'data/MBPM/sc/markers_MBPM_sc_notum.csv', row.names = F)

# Main cell type assignment after manual annotation based on DGE
seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.6 %in% c(0, 3, 6, 7, 8, 9, 10, 11, 14), 'Myeloid cells', 'NA')
seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.6 %in% c(19), 'CNS cells', seu$cell_type_main)
seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.6 %in% c(23), 'Endothelial cells', seu$cell_type_main)
seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.6 %in% c(15, 21), 'Stromal cells', seu$cell_type_main)
seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.6 %in% c(13, 16, 20, 22), 'B/Plasma cells', seu$cell_type_main)
seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.6 %in% c(1, 2, 4, 5, 12, 17), 'T/NK cells', seu$cell_type_main)
seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.6 %in% c(18), 'Low-quality cells', seu$cell_type_main)


# Save object and cell-type annotation
saveRDS(seu, 'data/MBPM/sc/data_MBPM_sc_notum.rds')
celltype <- seu@meta.data %>% dplyr::select('barcode_all', 'cell_type_main')
write.csv(celltype, 'data/MBPM/sc/data_MBPM_sc_celltype_main.csv')

# DGE based on cell types
DefaultAssay(seu) <- 'RNA'
seu <- NormalizeData(seu) %>% ScaleData()
DefaultAssay(seu) <- 'integrated'

Idents(seu) <- seu$cell_type_main
markers <- FindAllMarkers(seu, only.pos = TRUE, assay = 'RNA', min.pct = 0.25, logfc.threshold = 0.25)
ifelse(!dir.exists(file.path('data/cell_type_DEG/sc/global')), dir.create(file.path('data/cell_type_DEG/sc/global'), recursive = T), FALSE)
write.csv(markers, 'data/cell_type_DEG/sc/global/markers_MBPM_sc_global_type.csv', row.names = F)
