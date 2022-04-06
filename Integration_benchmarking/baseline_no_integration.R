#!/usr/bin/env Rscript

#### Baseline for integration analysis (raw) and LISI score
#### Author: Jana Biermann, PhD

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrastr)
library(gplots)
library(lisi)
library(tidyr)
library(magrittr)
library(viridis)
library(scales)

colBP <- c('#A80D11', '#008DB8')
colSCSN <- c('#E1AC24', '#288F56')

### Select one label
label <- 'MBPM_scn'
label <- 'MBPM_sn'
label <- 'MBPM_scn_tcell'
label <- 'MBPM_sn_tcell'


# Raw PCA coordinates
seu <- readRDS('data/MBPM/data_MBPM_scn.rds')
if (label == 'MBPM_sn') {
  seu <- subset(seu, sequencing == 'Single nuclei')
}

if (label == 'MBPM_scn_tcell') {
  seu <- subset(seu, cell_type_main == 'T/NK cells' & cell_type_fine != 'NK cells')
}

if (label == 'MBPM_sn_tcell') {
  seu <- subset(seu, sequencing == 'Single nuclei' & cell_type_main == 'T/NK cells' &
                  cell_type_fine != 'NK cells')
}

DefaultAssay(seu) <- 'RNA'
seu <- DietSeurat(seu, assays = 'RNA')
seu <- NormalizeData(seu) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(npcs = 50)
seu <- RunUMAP(seu, dims = 1:30)

DimPlot(seu, reduction = 'umap', group.by = 'patient', label = T, repel = T, shuffle = T,
        raster = T)

# LISI with raw UMAP coordinates
vars <- data.frame(patient = seu$patient, cell_type_main = seu$cell_type_main,
                   cell_type_int = seu$cell_type_int, cell_type_fine = seu$cell_type_fine,
                   sequencing = seu$sequencing, organ = seu$organ, 
                   row.names = rownames(seu@meta.data))

coords_raw_umap <- as.matrix(seu@reductions$umap@cell.embeddings[, 1:2])
res_raw_umap <- compute_lisi(coords_raw_umap, vars, names(vars))
write.csv(res_raw_umap, paste0('data/integration_benchmarking/lisi_', label,
                               '_raw_umap.csv'))

# Add LISI score to Seurat object
res_raw_umap <- read.csv(paste0('data/integration_benchmarking/lisi_', label,
                                '_raw_umap.csv'))
colnames(res_raw_umap)[1] <- 'barcode_all'
colnames(res_raw_umap)[2:ncol(res_raw_umap)] <- paste0('raw_umap_',
                                                       colnames(res_raw_umap[2:ncol(res_raw_umap)]))

seu@meta.data <- left_join(seu@meta.data, res_raw_umap, by = 'barcode_all')
rownames(seu@meta.data) <- seu$barcode_all


# Plots raw
pdf(paste0('data/integration_benchmarking/plots_', label, '_raw.pdf'))
DimPlot(seu, reduction = 'umap', group.by = 'patient', label = T, repel = T, shuffle = T,
        raster = T)
DimPlot(seu, reduction = 'umap', group.by = 'patient', label = F, repel = T, shuffle = T,
        raster = T) + NoLegend()
DimPlot(seu, reduction = 'umap', group.by = 'sequencing', label = T, repel = T,
        shuffle = T, raster = T, cols = colSCSN)
DimPlot(seu, reduction = 'umap', group.by = 'sequencing', label = F, repel = T,
        shuffle = T, raster = T, cols = colSCSN) + NoLegend()
DimPlot(seu, reduction = 'umap', group.by = 'organ', label = T, repel = T, shuffle = T,
        raster = T, cols = colBP)
DimPlot(seu, reduction = 'umap', group.by = 'organ', label = F, repel = T, shuffle = T,
        raster = T, cols = colBP) + NoLegend()
DimPlot(seu, reduction = 'umap', group.by = 'cell_type_main', shuffle = T, raster = T) +
  NoLegend()
DimPlot(seu, reduction = 'umap', group.by = 'cell_type_main', label = T, repel = T,
        shuffle = T, raster = T)

DimPlot(seu, reduction = 'umap', group.by = 'cell_type_int', shuffle = T, raster = T) +
  NoLegend()
DimPlot(seu, reduction = 'umap', group.by = 'cell_type_int', label = T, repel = T,
        shuffle = T, raster = T)

DimPlot(seu, reduction = 'umap', group.by = 'cell_type_fine', shuffle = T, raster = T) +
  NoLegend()
DimPlot(seu, reduction = 'umap', group.by = 'cell_type_fine', label = T, repel = T,
        shuffle = T, raster = T)

FeaturePlot(seu, reduction = 'umap', order = T, raster = T,
            features = c('raw_umap_patient', 'raw_umap_cell_type_fine',
                         'raw_umap_organ', 'raw_umap_sequencing'))
FeaturePlot(seu, reduction = 'umap', order = T, raster = T,
            features = c('raw_umap_cell_type_main', 'raw_umap_cell_type_int',
                         'raw_umap_cell_type_fine'))
dev.off()
