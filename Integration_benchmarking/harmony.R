#!/usr/bin/env Rscript

#### Integration using Harmony and LISI score
#### Author: Jana Biermann, PhD

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrastr)
library(harmony)
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


#### Run Harmony ####
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

seu <- RunHarmony(object = seu, group.by.vars = 'patient', dims.use = 1:30,
                  assay.use = 'RNA', plot_convergence = TRUE)
seu <- RunUMAP(seu, reduction = 'harmony', dims = 1:30)


# Save object
ifelse(!dir.exists(file.path('data/integration_benchmarking')),
       dir.create(file.path('data/integration_benchmarking'), recursive = T),
       FALSE)
saveRDS(seu, paste0('data/integration_benchmarking/data_', label, '_harmony.rds'))


#### LISI harmony ####
# coords: matrix of cells (rows) and coordinates (PC scores, tSNE or UMAP dimensions, etc.)
# vars: data frame with categorical variables (one row for each cell)
# The score indicates the effective number of different categories represented in the local neighborhood of each cell.
# If the cells are well-mixed, then we might expect the LISI score to be near 2 for a categorical variable with 2 categories.

seu <- readRDS(paste0('data/integration_benchmarking/data_', label, '_harmony.rds'))

# Harmony UMAP
vars <- data.frame(patient = seu$patient, cell_type_main = seu$cell_type_main,
                   cell_type_int = seu$cell_type_int, cell_type_fine = seu$cell_type_fine,
                   sequencing = seu$sequencing, organ = seu$organ, row.names = rownames(seu@meta.data))

coords_harmony_umap <- as.matrix(seu@reductions$umap@cell.embeddings[, 1:2])
res_harmony_umap <- compute_lisi(coords_harmony_umap, vars, names(vars))
write.csv(res_harmony_umap, paste0('data/integration_benchmarking/lisi_', label,
                                   '_harmony_umap.csv'))


# Add LISI score to Seurat object
res_harmony_umap <- read.csv(paste0('data/integration_benchmarking/lisi_', label,
                                    '_harmony_umap.csv'))
colnames(res_harmony_umap)[1] <- 'barcode_all'
colnames(res_harmony_umap)[2:ncol(res_harmony_umap)] <- paste0('harmony_umap_',
                                                               colnames(res_harmony_umap[2:ncol(res_harmony_umap)]))

seu@meta.data <- left_join(seu@meta.data, res_harmony_umap, by = 'barcode_all')
rownames(seu@meta.data) <- seu$barcode_all


# Plots harmony
pdf(paste0('data/integration_benchmarking/plots_', label, '_harmony.pdf'))
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

DimPlot(seu, reduction = 'harmony', group.by = 'patient', label = T, repel = T,
        shuffle = T, raster = T)
DimPlot(seu, reduction = 'harmony', group.by = 'sequencing', label = T, repel = T,
        shuffle = T, raster = T)
DimPlot(seu, reduction = 'harmony', group.by = 'organ', label = T, repel = T, shuffle = T,
        raster = T)

FeaturePlot(seu, reduction = 'umap', order = T, raster = T,
            features = c('harmony_umap_patient', 'harmony_umap_cell_type_fine',
                         'harmony_umap_organ', 'harmony_umap_sequencing'))
FeaturePlot(seu, reduction = 'umap', order = T, raster = T,
            features = c('harmony_umap_cell_type_main', 'harmony_umap_cell_type_int',
                         'harmony_umap_cell_type_fine'))
dev.off()
