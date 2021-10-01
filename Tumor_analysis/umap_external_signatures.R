#!/usr/bin/env Rscript

#### title: Tumor cell UMAPs and external signautres applied to tumor cells 
#### author: Jana Biermann, PhD

library(dplyr)
library(Seurat)
library(ggplot2)
library(gplots)

colBP <- c('#A80D11', '#008DB8')
colSCSN <- c('#E1AC24', '#288F56')
cycling <- c('#5AAE61', '#762A83')
colMITF <- c('#C7B122', '#C70E7B', '#22BAC7')


# Read-in integrated object
seu <- readRDS('data/MBPM/data_MBPM.rds')

# Subset to sn tumor cells
seu <- subset(seu, sequencing == 'Single nuclei' & cell_type_main == 'Tumor cells')
seu <- ScaleData(seu) %>% RunPCA() %>% RunUMAP(dims = 1:50)

# Apply signautres
sigs <- read.csv('signatures/melanoma_sigs.csv', na.strings = c('', 'NA'))
for (c in 1:ncol(sigs)) {
  sig <- as.character(na.omit(sigs[, c]))
  seu <- AddModuleScore(object = seu, features = list(sig), name = colnames(sigs[c]), assay = 'RNA', search = F)
}


pdf('data/tumor/plots_umaps_external_sigs.pdf')
# UMAPs
plt <- FeaturePlot(seu, features = c('AXL_sig1', 'MITF_sig1'), blend = T, cols = colMITF[2:3], raster = T, 
                   ncol = 2, combine = F)
CombinePlots(plt)
DimPlot(seu, group.by = 'cell_cycle', shuffle = T, raster = T, cols = cycling)
DimPlot(seu, group.by = 'cell_cycle', shuffle = T, raster = T, cols = cycling) + NoLegend()
DimPlot(seu, group.by = 'organ', shuffle = T, raster = T, cols = colBP)
DimPlot(seu, group.by = 'organ', shuffle = T, raster = T, cols = colBP) + NoLegend()
DimPlot(seu, group.by = 'patient', shuffle = T, raster = T)
DimPlot(seu, group.by = 'patient', shuffle = T, raster = T) + NoLegend()

# MITF + AXL
VlnPlot(seu, features = c('MITF_sig1', 'AXL_sig1'), pt.size = 0, group.by = 'organ', cols = colBP)

# Bulk signatures for MBM and MPM
VlnPlot(seu, features = c('fischer_DEG_mbm_100', 'fischer_DEG_mpm_100'), pt.size = 0, group.by = 'organ', 
        cols = colBP)
FeatureScatter(seu, feature1 = 'fischer_DEG_mbm_100', feature2 = 'fischer_DEG_mpm_100', group.by = 'organ', 
               raster = T, shuffle = T, cols = colBP)
dev.off()
