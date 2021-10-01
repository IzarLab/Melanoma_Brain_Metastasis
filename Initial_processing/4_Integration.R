#!/usr/bin/env Rscript

#### Seurat intergation of individual objects 
#### Author: Jana Biermann, PhD

library(dplyr)
library(Seurat)


# Sample list
pats <- c('MBM01_sc', 'MBM02_sc', 'MBM03_sc', 'MBM04_sc', 'MBM05_sc', 
          'MBM05_sn', 'MBM06_sn', 'MBM07_sn', 'MBM08_sn', 'MBM09_sn', 
          'MBM10_sn', 'MBM11_sn', 'MBM12_sn', 'MBM13_sn', 'MBM14_sn', 
          'MBM15_sn', 'MBM16_sn', 'MBM17_sn', 'MBM18_sn', 'MBM19_sn', 
          'MBM20_sn', 'MBM21_sn', 'MPM01_sn', 'MPM02_sn', 'MPM04_sn', 
          'MPM05_sn')

# Load objects
for (pat in pats) {
  assign(pat, readRDS(paste0('data/', pat, '/data_', pat, '_cb.rds')))
}

# Prepare obj.list for integration
obj.list <- NULL
for (obj in grep(paste(c('MBM', 'MPM'), collapse = '|'), ls(), value = T)) {
  obj.list <- c(object.list, eval(parse(text = obj)))
}

# Integration
anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = 'LogNormalize', verbose = T, 
                                  dims = 1:50, anchor.features = 2000)
seu <- IntegrateData(anchorset = anchors, normalization.method = 'LogNormalize', verbose = T, dims = 1:50)
rm(anchors, obj.list)

# Seurat workflow
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = 85)
seu <- RunUMAP(seu, dims = 1:85, spread = 1.5, min.dist = 0.1)

# Annotation
seu[['percent.rps']] <- PercentageFeatureSet(object = seu, pattern = '^RPS', assay = 'RNA')
seu[['percent.rpl']] <- PercentageFeatureSet(object = seu, pattern = '^RPL', assay = 'RNA')
seu$barcode_all <- rownames(seu@meta.data)
seu$tumor <- ifelse(is.na(seu$tumor) == T, 'normal', seu$tumor)

# Save object
ifelse(!dir.exists(file.path('data/MBPM/')), dir.create(file.path('data/MBPM/')), FALSE)
saveRDS(seu, file = 'data/MBPM/data_MBPM.rds')


### Subset to non-tumor cells for cell-type annotation
seu <- subset(seu, tumor != 'tumor')

# Redo workflow
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = 75)
seu <- RunUMAP(seu, dims = 1:75)
seu <- FindNeighbors(seu, dims = 1:75)
seu <- FindClusters(seu, resolution = 0.8)

# Save object
saveRDS(seu, file = 'data/MBPM/data_MBPM_notum.rds')


print(Sys.time())