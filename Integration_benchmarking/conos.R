#!/usr/bin/env Rscript

#### Integration using Conos and LISI score
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
library(igraph)
library(leidenAlg)
library(conos)

colBP <- c('#A80D11', '#008DB8')
colSCSN <- c('#E1AC24', '#288F56')

### Select one label
label <- 'MBPM_scn'
label <- 'MBPM_sn'
label <- 'MBPM_scn_tcell'
label <- 'MBPM_sn_tcell'


# Read-in object
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


#### Run Conos ####
DefaultAssay(seu) <- 'RNA'
seu <- DietSeurat(seu, assays = 'RNA')

# Split object
obj.list <- SplitObject(seu, split.by = 'patient')
for (i in 1:length(obj.list)) {
  if (dim(obj.list[[i]]@assays$RNA@counts)[2] < 50) {
    obj.list[i] <- NA
  } else {
    obj.list[[i]] <- NormalizeData(obj.list[[i]]) %>%
      FindVariableFeatures() %>%
      ScaleData() %>%
      RunPCA(npcs = 50)
  }
}

# Remove NAs
obj.list <- obj.list[!is.na(obj.list)]

# Prepare Conos object
con <- Conos$new(obj.list, n.cores = 1)

# Build the Joint Graph
con$buildGraph(space = 'PCA', ncomps = 50, verbose = TRUE,
               score.component.variance = TRUE)

ifelse(!dir.exists(file.path('data/integration_benchmarking')),
       dir.create(file.path('data/integration_benchmarking'), recursive = T),
       FALSE)
saveRDS(con, paste0('data/integration_benchmarking/data_', label, '_conos.rds'))

# Get Conos UMAP coordinates
con$embedGraph(method = 'UMAP')
write.csv(con$embeddings$UMAP, paste0('data/integration_benchmarking/conos_', label,
                                      '_umap_embedding.csv'))


#### Run LISI ####
# Conos UMAP coordinates
vars <- data.frame(patient = seu$patient, cell_type_main = seu$cell_type_main,
                   cell_type_int = seu$cell_type_int, cell_type_fine = seu$cell_type_fine,
                   sequencing = seu$sequencing, organ = seu$organ, row.names = rownames(seu@meta.data))

coords_conos_umap <- read.csv(paste0('data/integration_benchmarking/conos_', label,
                                     '_umap_embedding.csv'), row.names = 1)
colnames(coords_conos_umap) <- c('Conos_umap_1', 'Conos_umap_2')
coords_conos_umap <- as.matrix(coords_conos_umap)
vars <- vars[rownames(coords_conos_umap), ]
seu <- subset(seu, cells = rownames(coords_conos_umap))
res_conos_umap <- compute_lisi(coords_conos_umap, vars, names(vars))
write.csv(res_conos_umap, paste0('data/integration_benchmarking/lisi_', label,
                                 '_conos_umap.csv'))

# Add LISI score to Seurat object
res_conos_umap <- read.csv(paste0('data/integration_benchmarking/lisi_', label,
                                  '_conos_umap.csv'))
colnames(res_conos_umap)[1] <- 'barcode_all'
colnames(res_conos_umap)[2:ncol(res_conos_umap)] <- paste0('conos_umap_',
                                                           colnames(res_conos_umap[2:ncol(res_conos_umap)]))

seu$barcode_all <- rownames(seu@meta.data)
seu@meta.data <- left_join(seu@meta.data, res_conos_umap, by = 'barcode_all')
rownames(seu@meta.data) <- seu$barcode_all


# Plots
pdf(paste0('data/integration_benchmarking/plots_', label, '_conos.pdf'))
DimPlot(seu, reduction = 'conos_umap', group.by = 'patient', label = T, repel = T,
        shuffle = T, raster = T)
DimPlot(seu, reduction = 'conos_umap', group.by = 'patient', label = F, repel = T,
        shuffle = T, raster = T) + NoLegend()
DimPlot(seu, reduction = 'conos_umap', group.by = 'sequencing', label = T, repel = T,
        shuffle = T, raster = T, cols = colSCSN)
DimPlot(seu, reduction = 'conos_umap', group.by = 'sequencing', label = F, repel = T,
        shuffle = T, raster = T, cols = colSCSN) + NoLegend()
DimPlot(seu, reduction = 'conos_umap', group.by = 'organ', label = T, repel = T,
        shuffle = T, raster = T, cols = colBP)
DimPlot(seu, reduction = 'conos_umap', group.by = 'organ', label = F, repel = T,
        shuffle = T, raster = T, cols = colBP) + NoLegend()

DimPlot(seu, reduction = 'conos_umap', group.by = 'cell_type_main', label = T, repel = T,
        shuffle = T, raster = T)
DimPlot(seu, reduction = 'conos_umap', group.by = 'cell_type_main', label = F, repel = T,
        shuffle = T, raster = T) + NoLegend()
DimPlot(seu, reduction = 'conos_umap', group.by = 'cell_type_int', label = T, repel = T,
        shuffle = T, raster = T)
DimPlot(seu, reduction = 'conos_umap', group.by = 'cell_type_int', label = F, repel = T,
        shuffle = T, raster = T) + NoLegend()
DimPlot(seu, reduction = 'conos_umap', group.by = 'cell_type_fine', label = T, repel = T,
        shuffle = T, raster = T)
DimPlot(seu, reduction = 'conos_umap', group.by = 'cell_type_fine', label = F, repel = T,
        shuffle = T, raster = T) + NoLegend()

FeaturePlot(seu, reduction = 'conos_umap', order = T, raster = T,
            features = c('conos_umap_patient', 'conos_umap_cell_type_fine',
                         'conos_umap_organ', 'conos_umap_sequencing'))
FeaturePlot(seu, reduction = 'conos_umap', order = T, raster = T,
            features = c('conos_umap_cell_type_main', 'conos_umap_cell_type_int',
                         'conos_umap_cell_type_fine'))
dev.off()
