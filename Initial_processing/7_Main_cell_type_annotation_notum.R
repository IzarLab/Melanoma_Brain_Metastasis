#!/usr/bin/env Rscript

#### Main cell type annotation of integrated Seurat objects
#### Author: Jana Biermann, PhD

library(dplyr)
library(Seurat)
library(ggplot2)
library(gplots)

# Select one cohort
cohort <- 'MBM_sc'
cohort <- 'MBM_sn'
cohort <- 'MPM_sn'

# Set up folders
celltype <- 'main'
folder <- paste0('data/cell_type_DEG/', cohort, '/', celltype, '/')
ifelse(!dir.exists(file.path(folder)), 
       dir.create(file.path(folder), recursive = T), FALSE)

# Read in object
seu <- readRDS(paste0('data/MBPM/', cohort, '/data_', cohort, '_notum.rds'))

# Differential gene expression (DGE) based on clusters
markers <- FindAllMarkers(seu, assay = 'RNA', only.pos = T, min.pct = 0.25,
                          logfc.threshold = 0.25, test.use = 'MAST')
write.csv(markers, paste0(folder, 'markers_', cohort, '_notum.csv'), row.names = F)
DimPlot(seu)

# Main cell type assignment after manual annotation based on DGE
if (cohort == 'MBM_sc') {
  seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.3 %in% c(0, 1, 5, 7, 11, 12, 17), 
                               'Myeloid cells', 'NA')
  seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.3 %in% c(14, 15),
                               'CNS cells', seu$cell_type_main)
  seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.3 %in% c(20),
                               'Endothelial cells', seu$cell_type_main)
  seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.3 %in% c(19),
                               'Stromal cells', seu$cell_type_main)
  seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.3 %in% c(8, 10, 18), 
                               'B/Plasma cells', seu$cell_type_main)
  seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.3 %in% c(2, 3, 4, 6, 9, 13, 16), 
                               'T/NK cells', seu$cell_type_main)
}

if (cohort == 'MBM_sn') {
  seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.3 %in% c(0, 7, 14),
                               'Myeloid cells', 'NA')
  seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.3 %in% c(5, 10, 13), 
                               'CNS cells', seu$cell_type_main)
  seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.3 %in% c(12),
                               'Endothelial cells', seu$cell_type_main)
  seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.3 %in% c(6), 
                               'Stromal cells', seu$cell_type_main)
  seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.3 %in% c(9, 11),
                               'B/Plasma cells', seu$cell_type_main)
  seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.3 %in% c(1, 2, 4), 
                               'T/NK cells', seu$cell_type_main)
  seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.3 %in% c(3, 8),
                               'Low-quality cells', seu$cell_type_main)
}

if (cohort == 'MPM_sn') {
  seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.5 %in% c(1, 2, 5, 14, 17, 18), 
                               'Myeloid cells', 'NA')
  seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.5 %in% c(10),
                               'Endothelial cells', seu$cell_type_main)
  seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.5 %in% c(0, 12, 13, 19), 
                               'Stromal cells', seu$cell_type_main)
  seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.5 %in% c(4, 15, 16), 
                               'B/Plasma cells', seu$cell_type_main)
  seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.5 %in% c(3, 6, 7, 8, 11), 
                               'T/NK cells', seu$cell_type_main)
  seu$cell_type_main <- ifelse(seu$integrated_snn_res.0.5 %in% c(9),
                               'Low-quality cells', seu$cell_type_main)
}

DimPlot(seu, group.by = 'cell_type_main')

# Save object and annotation info
saveRDS(seu, paste0(folder, 'data_', cohort, '_notum.rds'))
ct <- seu@meta.data %>%
  dplyr::select('barcode_pat', 'cell_type_main')
write.csv(ct, paste0(folder, 'data_', cohort, '_celltype_main.csv'), row.names = F)

# Differential gene expression (DGE) based on cell types
DefaultAssay(seu) <- 'RNA'
seu <- NormalizeData(seu) %>%
  ScaleData()
DefaultAssay(seu) <- 'integrated'

Idents(seu) <- seu$cell_type_main
markers <- FindAllMarkers(seu, only.pos = TRUE, assay = 'RNA', min.pct = 0.25,
                          logfc.threshold = 0.25, test.use = 'MAST')
write.csv(markers, paste0(folder, 'markers_', cohort, '_main_celltype.csv'),
          row.names = F)
markers$cluster <- as.character(markers$cluster)
markers %>%
  filter(p_val_adj < 0.05) %>%
  arrange(cluster) %>%
  group_by(cluster) %>%
  top_n(-15, wt = p_val_adj) %>%
  group_by(cluster) %>%
  slice_head(n = 15) -> tops

pdf(file = paste0(folder, 'plots_', cohort, '_main_heatmap_celltype.pdf'), width = 10,
    height = 17)
DoHeatmap(seu, features = tops$gene, group.by = 'cell_type_main', raster = T,
          assay = 'RNA')
dev.off()

