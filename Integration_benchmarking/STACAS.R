#!/usr/bin/env Rscript

#### Title: Integration of T cells using STACAS and LISI score
#### Authors: Jana Biermann, PhD; Massimo Andreatta

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
library(STACAS)
library(TILPRED)
library(AUCell)
library(SingleCellExperiment)
library(reshape2)
library(grid)
library(gridExtra)

colBP <- c('#A80D11', '#008DB8')
colSCSN <- c('#E1AC24', '#288F56')

### Select one label
label <- 'MBPM_scn_tcell'
label <- 'MBPM_sn_tcell'


#### Run STACAS ####
seu <- readRDS('data/MBPM/data_MBPM_scn.rds')
if (label == 'MBPM_scn_tcell') {
  seu <- subset(seu, cell_type_main == 'T/NK cells' & cell_type_fine != 'NK cells')
}

if (label == 'MBPM_sn_tcell') {
  seu <- subset(seu, sequencing == 'Single nuclei' & cell_type_main == 'T/NK cells' &
                  cell_type_fine != 'NK cells')
}

cellCycle.symbol <- read.csv('misc/stacas_ccgenes.csv', as.is = T)$x

DefaultAssay(seu) <- 'RNA'
seu <- DietSeurat(seu, assays = 'RNA')
data.list <- SplitObject(seu, split.by = 'patient')

var.genes.n <- 800
var.genes.integrated.n <- 500

for (i in 1:length(data.list)) {
  if (dim(data.list[[i]]@assays$RNA@counts)[2] < 50) {
    data.list[i] <- NA
  } else {
    data.list[[i]] <- NormalizeData(data.list[[i]], verbose = FALSE)
    data.list[[i]] <- FindVariableFeatures(data.list[[i]], selection.method = 'vst', 
                                           nfeatures = var.genes.n * 2, verbose = FALSE)
    mito.genes <- grep(pattern = '^MT-', rownames(data.list[[i]]), value = TRUE)
    ribo.genes <- grep(pattern = '^RP[LS]', rownames(data.list[[i]]), value = TRUE)
    
    data.list[[i]]@assays$RNA@var.features <- setdiff(data.list[[i]]@assays$RNA@var.features,
                                                      cellCycle.symbol)
    data.list[[i]]@assays$RNA@var.features <- setdiff(data.list[[i]]@assays$RNA@var.features,
                                                      mito.genes)
    data.list[[i]]@assays$RNA@var.features <- setdiff(data.list[[i]]@assays$RNA@var.features,
                                                      ribo.genes)
    data.list[[i]]@assays$RNA@var.features <- head(data.list[[i]]@assays$RNA@var.features,
                                                   var.genes.n)
  }
}
# Remove NAs
data.list <- data.list[!is.na(data.list)]

# Get anchors
ndim <- 20
ref.anchors <- FindAnchors.STACAS(data.list, dims = 1:ndim,
                                  anchor.features = var.genes.integrated.n)

# Filter anchors with default threshold 
ref.anchors.filtered <- ref.anchors  # no filtering

all.genes <- row.names(data.list[[1]])
for (i in 2:length(data.list)) {
  all.genes <- intersect(all.genes, row.names(data.list[[i]]))
}

mySampleTree <- SampleTree.STACAS(ref.anchors.filtered)
print(mySampleTree)
names(data.list)

ref.integrated <- IntegrateData(anchorset = ref.anchors.filtered, dims = 1:20,
                                features.to.integrate = all.genes, sample.tree = mySampleTree, preserve.order = T,
                                k.weight = 50)

set.seed(1234)
ndim <- 75
length(ref.integrated@assays$integrated@var.features)

ref.integrated <- ScaleData(ref.integrated, verbose = TRUE)
ref.integrated <- RunPCA(ref.integrated,
                         features = ref.integrated@assays$integrated@var.features,
                         ndims.print = 1:5, nfeatures.print = 5)

ndim <- 20  #how many PCA components to retain
ref.integrated <- RunUMAP(ref.integrated, reduction = 'pca', dims = 1:ndim,
                          seed.use = 123, n.neighbors = 30, min.dist = 0.3)
saveRDS(ref.integrated, paste0('data/integration_benchmarking/data_', label,
                               '_stacas.rds'))


#### Run LISI ####
seu <- readRDS(paste0('data/integration_benchmarking/data_', label, '_stacas.rds'))

# LISI with UMAP coordinates
vars <- data.frame(patient = seu$patient, cell_type_main = seu$cell_type_main,
                   cell_type_int = seu$cell_type_int, cell_type_fine = seu$cell_type_fine,
                   sequencing = seu$sequencing, organ = seu$organ, row.names = rownames(seu@meta.data))

coords_stacas_umap <- as.matrix(seu@reductions$umap@cell.embeddings[, 1:2])
res_stacas_umap <- compute_lisi(coords_stacas_umap, vars, names(vars))
write.csv(res_stacas_umap, paste0('data/integration_benchmarking/lisi_', label,
                                  '_stacas_umap.csv'))

# Add LISI score to Seurat object
res_stacas_umap <- read.csv(paste0('data/integration_benchmarking/lisi_', label,
                                   '_stacas_umap.csv'))
colnames(res_stacas_umap)[1] <- 'barcode_all'
colnames(res_stacas_umap)[2:ncol(res_stacas_umap)] <- paste0('stacas_umap_',
                                                             colnames(res_stacas_umap[2:ncol(res_stacas_umap)]))

seu$barcode_all <- rownames(seu@meta.data)
seu@meta.data <- left_join(seu@meta.data, res_stacas_umap, by = 'barcode_all')
rownames(seu@meta.data) <- seu$barcode_all


# Plots 
pdf(paste0('data/integration_benchmarking/plots_', label, '_stacas.pdf'))
DimPlot(seu, reduction = 'pca', group.by = 'patient', label = F, repel = T, shuffle = T,
        raster = T)
DimPlot(seu, reduction = 'pca', group.by = 'sequencing', label = T, repel = T,
        shuffle = T, raster = T)
DimPlot(seu, reduction = 'pca', group.by = 'organ', label = T, repel = T, shuffle = T,
        raster = T)
DimPlot(seu, reduction = 'pca', group.by = 'cell_type_main', label = F, repel = T,
        shuffle = T, raster = T)
DimPlot(seu, reduction = 'pca', group.by = 'cell_type_int', label = F, repel = T,
        shuffle = T, raster = T)
DimPlot(seu, reduction = 'pca', group.by = 'cell_type_fine', label = F, repel = T,
        shuffle = T, raster = T)

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
            features = c('stacas_umap_patient', 'stacas_umap_cell_type_fine',
                         'stacas_umap_organ', 'stacas_umap_sequencing'))
FeaturePlot(seu, reduction = 'umap', order = T, raster = T,
            features = c('stacas_umap_cell_type_main', 'stacas_umap_cell_type_int',
                         'stacas_umap_cell_type_fine'))
dev.off()
