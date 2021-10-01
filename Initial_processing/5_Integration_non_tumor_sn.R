#!/usr/bin/env Rscript

#### Seurat re-intergation of single-nuclei samples 
#### Author: Jana Biermann, PhD

library(dplyr)
library(Seurat)

# Read-in non-tumor object
seu <- readRDS('data/MBPM/data_MBPM_notum.rds')
DefaultAssay(seu) <- 'RNA'
seu <- DietSeurat(seu, assays = 'RNA')
seu <- subset(seu, sequencing == 'Single nuclei')

# Integration
obj.list <- SplitObject(seu, split.by = 'patient')
anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = 'LogNormalize', verbose = T)
seu <- IntegrateData(anchorset = anchors, normalization.method = 'LogNormalize', verbose = T)
rm(anchors, obj.list)

# Seurat workflow
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = 75)
seu <- RunUMAP(seu, dims = 1:75)
seu <- FindNeighbors(seu, dims = 1:75)
seu <- FindClusters(seu, resolution = 0.6)

# Save object
ifelse(!dir.exists(file.path('data/MBPM/sn/')), dir.create(file.path('data/MBPM/sn')), FALSE)
saveRDS(seu, file = 'data/MBPM/sn/data_MBPM_sn_notum.rds')


print(Sys.time())