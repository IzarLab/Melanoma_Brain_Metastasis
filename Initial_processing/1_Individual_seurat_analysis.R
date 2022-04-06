#!/usr/bin/env Rscript

#### Seurat analysis using CellBender output with sample name provided as argument
#### Author: Jana Biermann, PhD

print(paste('Start:', Sys.time()))

library(dplyr)
library(Seurat)
library(purrr)
library(DropletUtils)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(scater)
library(Matrix)
library(DoubletFinder)

pat <- commandArgs()[6]
doublet_rate <- read.csv('misc/mbpm_doublet_rate.csv')
print(pat)


### Load filtered.h5 files from AWS
if (grepl("sc", pat) == T) {
  system(paste0("aws s3 sync s3://scrna-seq/cellbender/CRv6.1.1_CBv0.2.0/", pat, "/ ~/data/", pat,
                "/ --exclude '*' --include '", pat, "_filtered.h5' --quiet"))
}
if (grepl("sn", pat) == T) {
  system(paste0("aws s3 sync s3://snrna-seq/cellbender/CRv6.1.1_CBv0.2.0/", pat, "/ ~/data/", pat,
                "/ --exclude '*' --include '", pat, "_filtered.h5' --quiet"))
}
seu.data <- Read10X_h5(paste0('data/', pat, '/', pat, '_filtered.h5'), use.names = TRUE,
                       unique.features = TRUE)

### Initialize the Seurat object with raw data
seu_raw <- CreateSeuratObject(counts = seu.data, project = pat, min.cells = 1, 
                              min.features = 1)

### Meta data annotations
seu_raw <- PercentageFeatureSet(seu_raw, pattern = '^MT-', col.name = 'percent.mt')
seu_raw <- PercentageFeatureSet(seu_raw, pattern = '^RPS', col.name = 'percent.rps')
seu_raw <- PercentageFeatureSet(seu_raw, pattern = '^RPL', col.name = 'percent.rpl')
seu_raw$percent.rp <- seu_raw$percent.rps + seu_raw$percent.rpl
seu_raw$barcode <- rownames(seu_raw@meta.data)
seu_raw$barcode_pat <- paste0(rownames(seu_raw@meta.data), '_', pat)

if (grepl('sc', pat)) {
  seu_raw[['patient']] <- substr(pat, 1, 8)
  seu_raw[['ID']] <- strsplit(pat, '_')[[1]][1]
  seu_raw[['batch']] <- strsplit(pat, '_')[[1]][3]
  seu_raw[['sequencing']] <- 'Single cells'
  seu_raw[['organ']] <- 'Brain'
  seu_raw[['cancer']] <- 'MBM'
}
if (grepl('MPM', pat)) {
  seu_raw[['patient']] <- pat
  seu_raw[['ID']] <- strsplit(pat, '_')[[1]][1]
  seu_raw[['sequencing']] <- 'Single nuclei'
  seu_raw[['organ']] <- 'Peripheral'
  seu_raw[['cancer']] <- 'MPM'
}
if (grepl('MBM', pat) & grepl('sn', pat)) {
  seu_raw[['patient']] <- pat
  seu_raw[['ID']] <- strsplit(pat, '_')[[1]][1]
  seu_raw[['sequencing']] <- 'Single nuclei'
  seu_raw[['organ']] <- 'Brain'
  seu_raw[['cancer']] <- 'MBM'
}

### Identify doublets using scrublet
doublet_rate_tmp <- doublet_rate[doublet_rate$sample == pat, 2]
writeMM(seu_raw@assays$RNA@counts, paste0('data/', pat, '/matrix_', pat, '_raw.mtx'))
system(paste('python3 brain_mets/scrublet/scrublet_code.py', pat, doublet_rate_tmp))
doublets <- read.table(paste0('data/', pat, '/doublets_', pat, '_raw.txt'), header = T)
seu_raw[['predicted_doublets']] <- doublets$predicted_doublets
seu_raw[['doublet_scores']] <- doublets$doublet_scores
system(paste0('rm data/', pat, '/matrix_', pat, '_raw.mtx'))
system(paste0('rm data/', pat, '/doublets_', pat, '_raw.txt'))

### Seurat workflow
seu_raw <- NormalizeData(seu_raw)
seu_raw <- FindVariableFeatures(seu_raw)
seu_raw <- ScaleData(seu_raw)
seu_raw <- RunPCA(seu_raw)
seu_raw <- RunUMAP(seu_raw, dims = 1:40)
seu_raw <- FindNeighbors(seu_raw, dims = 1:40)
seu_raw <- FindClusters(seu_raw, resolution = 0.1)

### Identify doublets using DF
sweep.list <- paramSweep_v3(seu_raw, PCs = 1:40, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK <- bcmvn %>%
  arrange(desc(BCmetric))
pK <- pK[1, 2]
pK <- as.numeric(levels(pK[[1]]))[pK[[1]]]
nExp <- round(doublet_rate_tmp * dim(seu_raw@assays$RNA@counts)[2])
seu_raw <- doubletFinder_v3(seu_raw, PCs = 1:40, pK = pK, nExp = nExp)
seu_raw$doublet <- seu_raw@meta.data[, paste0('DF.classifications_0.25_', pK, '_', nExp)]
seu_raw$DF_score <- seu_raw@meta.data[, paste0('pANN_0.25_', pK, '_', nExp)]

### Filtering steps 
minFeature <- 500
maxFeature <- 10000
minCount <- 1000
maxCount <- 60000
maxMT <- 10

seu <- subset(seu_raw, subset = nFeature_RNA > minFeature & nFeature_RNA < maxFeature &
                nCount_RNA > minCount & nCount_RNA < maxCount & percent.mt < maxMT & 
                doublet == 'Singlet' & predicted_doublets == F)

### Seurat workflow
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = 50)
seu <- RunUMAP(seu, dims = 1:40)
seu <- FindNeighbors(seu, dims = 1:40)
seu <- FindClusters(seu)

### Preliminary cell type identification using SingleR
seu_sce <- as.SingleCellExperiment(seu)

bped <- BlueprintEncodeData()
pred_bped_main <- SingleR(test = seu_sce, ref = bped, labels = bped$label.main)
pruneScores(pred_bped_main)
seu[['celltype_bped_main']] <- pred_bped_main$pruned.labels
pred_bped_fine <- SingleR(test = seu_sce, ref = bped, labels = bped$label.fine)
pruneScores(pred_bped_fine)
seu[['celltype_bped_fine']] <- pred_bped_fine$pruned.labels

iced <- DatabaseImmuneCellExpressionData()
pred_iced_main <- SingleR(test = seu_sce, ref = iced, labels = iced$label.main)
pruneScores(pred_iced_main)
seu[['celltype_iced_main']] <- pred_iced_main$pruned.labels
pred_iced_fine <- SingleR(test = seu_sce, ref = iced, labels = iced$label.fine)
pruneScores(pred_iced_fine)
seu[['celltype_iced_fine']] <- pred_iced_fine$pruned.labels

hpca <- HumanPrimaryCellAtlasData()
pred_hpca_main <- SingleR(test = seu_sce, ref = hpca, labels = hpca$label.main)
pruneScores(pred_hpca_main)
seu[['celltype_hpca_main']] <- pred_hpca_main$pruned.labels
pred_hpca_fine <- SingleR(test = seu_sce, ref = hpca, labels = hpca$label.fine)
pruneScores(pred_hpca_fine)
seu[['celltype_hpca_fine']] <- pred_hpca_fine$pruned.labels

mid <- MonacoImmuneData()
pred_mid_main <- SingleR(test = seu_sce, ref = mid, labels = mid$label.main)
pruneScores(pred_mid_main)
seu[['celltype_mid_main']] <- pred_mid_main$pruned.labels
pred_mid_fine <- SingleR(test = seu_sce, ref = mid, labels = mid$label.fine)
pruneScores(pred_mid_fine)
seu[['celltype_mid_fine']] <- pred_mid_fine$pruned.labels


### Save object
ifelse(!dir.exists(file.path(paste0('data/', pat, '/'))), 
       dir.create(file.path(paste0('data/', pat, '/'))), FALSE)
seu@meta.data[, paste0('DF.classifications_0.25_', pK, '_', nExp)] <- NULL
seu@meta.data[, paste0('pANN_0.25_', pK, '_', nExp)] <- NULL
saveRDS(seu, file = paste0('data/', pat, '/data_', pat, '_cb_DF.rds'))

### Sync out
system(paste0("aws s3 sync data/", pat, "/ s3://snrna-seq/Seurat/", pat,
              "/ --exclude '*' --include '*_cb*' --exclude '.*' --quiet"))

print(paste('End:', Sys.time()))