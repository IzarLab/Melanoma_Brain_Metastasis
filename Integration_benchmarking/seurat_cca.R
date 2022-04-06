#!/usr/bin/env Rscript

#### Seurat CCA integration and LISI score
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

# Seurat PCA coordinates
seu <- readRDS('data/MBPM/data_MBPM_scn.rds')

if (label == 'MBPM_sn') {
  seu <- readRDS('data/MBPM/MBPM_sn/data_MBPM_sn_anchor2000_dims30.rds')
  ct <- read.csv('data/cell_type_DEG/cell_types_MBPM_sn.csv')
  seu$barcode_all <- rownames(seu@meta.data)
  seu@meta.data <- left_join(seu@meta.data, ct, by = 'barcode_pat')
  rownames(seu@meta.data) <- seu$barcode_all
}

if (label == 'MBPM_scn_tcell') {
  seu <- subset(seu, cell_type_main == 'T/NK cells' & cell_type_fine != 'NK cells')
  table(seu$patient) %>%
    sort
  
  DefaultAssay(seu) <- 'RNA'
  obj.list <- SplitObject(seu, split.by = 'patient')
  
  # Keep only patients with >50 cells
  for (i in 1:length(obj.list)) {
    obj.list[[i]] <- NormalizeData(obj.list[[i]])
    if (dim(obj.list[[i]]@assays$RNA@counts)[2] < 50) {
      obj.list[i] <- NA
    }
  }
  # Remove NAs
  obj.list <- obj.list[!is.na(obj.list)]
  
  # Find integration anchors
  anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:30, k.filter = 50)
  
  # Integrate data sets
  seu <- IntegrateData(anchorset = anchors, dims = 1:30, k.weight = 50)
  
  # Rerun Seurat workflow
  seu <- ScaleData(seu) %>%
    RunPCA() %>%
    RunUMAP(dims = 1:30)
  
  saveRDS(seu, paste0('data/integration_benchmarking/data_', label, '_reint.rds'))
}

if (label == 'MBPM_sn_tcell') {
  seu <- subset(seu, sequencing == 'Single nuclei' & cell_type_main == 'T/NK cells' &
                  cell_type_fine != 'NK cells')
  
  DefaultAssay(seu) <- 'RNA'
  obj.list <- SplitObject(seu, split.by = 'patient')
  
  # Keep only patients with >50 cells
  for (i in 1:length(obj.list)) {
    obj.list[[i]] <- NormalizeData(obj.list[[i]])
    if (dim(obj.list[[i]]@assays$RNA@counts)[2] < 50) {
      obj.list[i] <- NA
    }
  }
  # Remove NAs
  obj.list <- obj.list[!is.na(obj.list)]
  
  # Find integration anchors
  anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:30, k.filter = 50)
  
  # Integrate data sets
  seu <- IntegrateData(anchorset = anchors, dims = 1:30, k.weight = 50)
  
  # Rerun Seurat workflow
  seu <- ScaleData(seu) %>%
    RunPCA() %>%
    RunUMAP(dims = 1:30)
  
  saveRDS(seu, paste0('data/integration_benchmarking/data_', label, '_reint.rds'))
}


DimPlot(seu, reduction = 'umap', group.by = 'patient', label = T, repel = T, shuffle = T,
        raster = T)


# LISI with UMAP coordinates
vars <- data.frame(patient = seu$patient, cell_type_main = seu$cell_type_main,
                   cell_type_int = seu$cell_type_int, cell_type_fine = seu$cell_type_fine,
                   sequencing = seu$sequencing, organ = seu$organ, row.names = rownames(seu@meta.data))

coords_seurat_umap <- as.matrix(seu@reductions$umap@cell.embeddings[, 1:2])
res_seurat_umap <- compute_lisi(coords_seurat_umap, vars, names(vars))
write.csv(res_seurat_umap, paste0('data/integration_benchmarking/lisi_', label,
                                  '_seurat_umap.csv'))


# Add LISI score to Seurat object
res_seurat_umap <- read.csv(paste0('data/integration_benchmarking/lisi_', label,
                                   '_seurat_umap.csv'))
colnames(res_seurat_umap)[1] <- 'barcode_all'
colnames(res_seurat_umap)[2:ncol(res_seurat_umap)] <- paste0('seurat_umap_',
                                                             colnames(res_seurat_umap[2:ncol(res_seurat_umap)]))

seu$barcode_all <- rownames(seu@meta.data)
seu@meta.data <- left_join(seu@meta.data, res_seurat_umap, by = 'barcode_all')
rownames(seu@meta.data) <- seu$barcode_all


# Plots 
pdf(paste0('data/integration_benchmarking/plots_', label, '_seurat.pdf'))
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
            features = c('seurat_umap_patient', 'seurat_umap_cell_type_fine',
                         'seurat_umap_organ', 'seurat_umap_sequencing'))
FeaturePlot(seu, reduction = 'umap', order = T, raster = T,
            features = c('seurat_umap_cell_type_main', 'seurat_umap_cell_type_int',
                         'seurat_umap_cell_type_fine'))
dev.off()
