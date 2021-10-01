#!/usr/bin/env Rscript

#### Seurat analysis for cellbender output of snRNA-seq (sample name and doublet rate provided as arguments)
#### Author: Jana Biermann, PhD

print(paste('Start:', Sys.time()))

library(dplyr)
library(Seurat)
library(SingleR)
library(SingleCellExperiment)
library(Matrix)


pat <- commandArgs()[6]
doublet_rate <- doublet_rate


### Load dataset
seu.data <- Read10X_h5(paste0('data/', pat, '/', pat, '_filtered.h5'), use.names = TRUE, unique.features = TRUE)

### Initialize Seurat object with raw data
seu_raw <- CreateSeuratObject(counts = seu.data, project = pat, min.cells = 1, min.features = 1)

# Annotate MT genes
seu_raw[['percent.mt']] <- PercentageFeatureSet(seu_raw, pattern = '^MT-')

# Annotate patient info
seu_raw[['patient']] <- df[pat, ]

# Identify doublets using scrublet
writeMM(seu_raw@assays$RNA@counts, paste0('data/', pat, '/matrix_', pat, '_raw.mtx'))
system(paste('python3 misc/scrublet_code.py', pat, doublet_rate))
doublets <- read.table(paste0('data/', pat, '/doublets_', pat, '_raw.txt'), header = T)
seu_raw[['predicted_doublets']] <- doublets$predicted_doublets
seu_raw[['doublet_scores']] <- doublets$doublet_scores
system(paste0('rm data/', pat, '/matrix_', pat, '_raw.mtx'))
system(paste0('rm data/', pat, '/doublets_', pat, '_raw.txt'))


### QC filtering
minFeature <- 500
maxFeature <- 7500
minCount <- 1000
maxCount <- 40000
maxMT <- 20
seu <- subset(seu_raw, subset = nFeature_RNA > minFeature & nFeature_RNA < maxFeature & nCount_RNA > minCount &
  nCount_RNA < maxCount & percent.mt < maxMT & predicted_doublets == F)

### Standard workflow
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu, features = rownames(seu))
seu <- RunPCA(seu)
seu <- JackStraw(seu, num.replicate = 100)
seu <- ScoreJackStraw(seu, reduction = 'pca', dims = 1:20)
seu <- FindNeighbors(seu, dims = 1:20)
seu <- FindClusters(seu)
seu <- RunUMAP(seu, dims = 1:20)

### Preliminary cell type identification
seu_sce <- as.SingleCellExperiment(seu)

bped <- BlueprintEncodeData()
pred_bped_main <- SingleR(test = seu_sce, ref = bped, labels = bped$label.main)
pruneScores(pred_bped_main)
seu[['celltype_bped_main']] <- pred_bped_main$pruned.labels
pred_bped_fine <- SingleR(test = seu_sce, ref = bped, labels = bped$label.fine)
pruneScores(pred_bped_fine)
seu[['celltype_bped_fine']] <- pred_bped_fine$pruned.labels

### Save object
saveRDS(seu, file = paste0('data/', pat, '/data_', seu$patient[1], '_cb.rds'))

print(paste('End:', Sys.time()))