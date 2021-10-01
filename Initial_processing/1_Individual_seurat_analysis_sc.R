#!/usr/bin/env Rscript

#### Seurat analysis for cellbender output of scRNA-seq (sample name provided as argument)
#### Author: Jana Biermann, PhD

print(paste('Start:', Sys.time()))

library(dplyr)
library(Seurat)
library(rscrublet)
library(SingleR)
library(SingleCellExperiment)


pat <- commandArgs()[6]
pat_num <- substr(pat, 8, 15)


### Load dataset
CD45neg.data <- Read10X(data.dir = paste0('data/', pat, '/CD45negGEXBI', pat_num, '/filtered_feature_bc_matrix'))
CD45pos.data <- Read10X(data.dir = paste0('data/', pat, '/CD45posGEXBI', pat_num, '/filtered_feature_bc_matrix'))

### Initialize Seurat object with raw data
CD45neg_raw <- CreateSeuratObject(counts = CD45neg.data, project = paste0(pat, '_CD45neg'), min.cells = 3, min.features = 200)
CD45pos_raw <- CreateSeuratObject(counts = CD45pos.data, project = paste0(pat, '_CD45pos'), min.cells = 3, min.features = 200)

# Annotate MT genes
CD45neg_raw[['percent.mt']] <- PercentageFeatureSet(CD45neg_raw, pattern = '^MT-')
CD45pos_raw[['percent.mt']] <- PercentageFeatureSet(CD45pos_raw, pattern = '^MT-')

# Annotate patient info
CD45neg_raw[['patient']] <- pat
CD45pos_raw[['patient']] <- pat

# Identify doublets using rscrublet
scrub_CD45neg <- scrubDoublets(as.matrix(CD45neg_raw@assays$RNA@counts), expected_doublet_rate = 0.037)
scrub_CD45pos <- scrubDoublets(as.matrix(CD45pos_raw@assays$RNA@counts), expected_doublet_rate = 0.037)
CD45neg_raw[['predicted_doublets']] <- scrub_CD45neg$scrubDoublets
CD45pos_raw[['predicted_doublets']] <- scrub_CD45pos$scrubDoublets
CD45neg_raw[['doublet_scores']] <- scrub_CD45neg$doublet_scores_obs
CD45pos_raw[['doublet_scores']] <- scrub_CD45pos$doublet_scores_obs

### QC filtering
CD45neg <- subset(CD45neg_raw, subset = nFeature_RNA > 750 & nFeature_RNA < 6000 & nCount_RNA > 1000 & nCount_RNA < 100000 & percent.mt < 10 & ScrubDoublet == F)
CD45pos <- subset(CD45pos_raw, subset = nFeature_RNA > 750 & nFeature_RNA < 6000 & nCount_RNA > 1000 & nCount_RNA < 100000 & percent.mt < 10 & ScrubDoublet == F)

### Normalization
CD45neg <- NormalizeData(CD45neg, normalization.method = 'LogNormalize', scale.factor = 10000)
CD45pos <- NormalizeData(CD45pos, normalization.method = 'LogNormalize', scale.factor = 10000)

### Identification of highly variable features
CD45neg <- FindVariableFeatures(CD45neg, selection.method = 'vst', nfeatures = 2000)
CD45pos <- FindVariableFeatures(CD45pos, selection.method = 'vst', nfeatures = 2000)

### Scaling the data
CD45neg <- ScaleData(CD45neg, features = rownames(CD45neg))
CD45pos <- ScaleData(CD45pos, features = rownames(CD45pos))

### Perform linear dimensional reduction
# PCA
CD45neg <- RunPCA(CD45neg, features = VariableFeatures(object = CD45neg))
CD45pos <- RunPCA(CD45pos, features = VariableFeatures(object = CD45pos))

### Determine the dimensionality of the dataset
CD45neg <- JackStraw(CD45neg, num.replicate = 100)
CD45neg <- ScoreJackStraw(CD45neg, reduction = 'pca', dims = 1:20)
CD45pos <- JackStraw(CD45pos, num.replicate = 100)
CD45pos <- ScoreJackStraw(CD45pos, reduction = 'pca', dims = 1:20)

### Cluster the cells
CD45neg <- FindNeighbors(CD45neg, dims = 1:15)
CD45neg <- FindClusters(CD45neg, resolution = 1.2)
CD45pos <- FindNeighbors(CD45pos, dims = 1:15)
CD45pos <- FindClusters(CD45pos, resolution = 1.2)

### Run non-linear dimensional reduction (UMAP/tSNE)
CD45neg <- RunUMAP(CD45neg, dims = 1:30)
CD45neg <- RunTSNE(CD45neg, dims = 1:30)
CD45pos <- RunUMAP(CD45pos, dims = 1:30)
CD45pos <- RunTSNE(CD45pos, dims = 1:30)

### Save objects
saveRDS(CD45neg, file = paste0('data/', pat, '/', CD45neg@meta.data$orig.ident[1], '.rds'))
saveRDS(CD45pos, file = paste0('data/', pat, '/', CD45pos@meta.data$orig.ident[1], '.rds'))


##### Merge CD45neg and CD45pos #####
patient <- merge(CD45neg, y = CD45pos, add.cell.ids = c('CD45neg', 'CD45pos'), project = pat)

### Normalization
patient <- NormalizeData(patient, normalization.method = 'LogNormalize', scale.factor = 10000)

### Identification of highly variable features
patient <- FindVariableFeatures(patient, selection.method = 'vst', nfeatures = 2000)

### Scaling the data
patient <- ScaleData(patient, features = rownames(patient))

### Perform linear dimensional reduction
# PCA
patient <- RunPCA(patient, features = VariableFeatures(object = patient))

### Determine the dimensionality of the dataset
patient <- JackStraw(patient, num.replicate = 100)
patient <- ScoreJackStraw(patient, reduction = 'pca', dims = 1:20)

### Cluster the cells
patient <- FindNeighbors(patient, dims = 1:15)
patient <- FindClusters(patient, resolution = 1.2)

### Run non-linear dimensional reduction (UMAP/tSNE)
patient <- RunUMAP(patient, dims = 1:30)
patient <- RunTSNE(patient, dims = 1:30)

### Preliminary cell type identification
patient_sce <- as.SingleCellExperiment(patient)
patient_pred_bped_fine <- SingleR(test = patient_sce, ref = bped, labels = bped$label.fine)
pruneScores(patient_pred_bped_fine)
patient[['celltype_bped_fine']] <- patient_pred_bped_fine$pruned.labels

### Save objects
saveRDS(patient, file = paste0('data/', pat, '/', pat, '.rds'))

print(paste('End:', Sys.time()))