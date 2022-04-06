#!/usr/bin/env Rscript

#### Merge CD45pos and CD45neg batches for scRNA-seq samples
#### Author: Jana Biermann, PhD

library(dplyr)
library(Seurat)
library(gplots)
library(ggplot2)
library(ggrastr)

pat <- commandArgs()[6]  #provide patient name, e.g. 'MBM02_sc'
print(pat)

# Read in
scCD45neg <- readRDS(paste0('data/', pat, '_CD45neg/data_', pat, '_CD45neg_cb_DF.rds'))
scCD45pos <- readRDS(paste0('data/', pat, '_CD45pos/data_', pat, '_CD45pos_cb_DF.rds'))

# Merge and rerun workflow
merged <- merge(scCD45neg, y = scCD45pos, project = pat)
merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged)
merged <- ScaleData(merged)
merged <- RunPCA(merged)
merged <- FindNeighbors(merged, dims = 1:40)
merged <- FindClusters(merged)
merged <- RunUMAP(merged, dims = 1:40)

# Save object
ifelse(!dir.exists(file.path(paste0('data/', pat))), 
       dir.create(file.path(paste0('data/',pat)), recursive = T), FALSE)
saveRDS(merged, file = paste0('data/', pat, '/data_', pat, '_cb_DF.rds'))

# Sync out
system(paste0("aws s3 sync data/", pat, "/ s3://snrna-seq/Seurat/", pat,
              "/ --exclude '*' --include '*_cb*' --exclude '.*' --quiet"))
