#!/usr/bin/env Rscript

### title: Analysis of non-immune and non-tumor cell types (UMAPs, DEG)
### author: Jana Biermann, PhD

print(Sys.time())

library(Seurat)
library(dplyr)
library(ggplot2)
library(gplots)
'%notin%' <- Negate('%in%')

colBP <- c('#A80D11', '#008DB8')
colSCSN <- c('#E1AC24', '#288F56')
colCNS <- c('#92D84F', '#7473A6', '#FF8C8D')


#### Central nervous system (CNS) cells ####
seu <- readRDS('data/MBPM/data_MBPM.rds')
seu <- subset(seu, cell_type_int %notin% c('Low-quality cells', 'Doublets', 'Contamination', 'Not determined'))
seu <- subset(seu, cell_type_main == 'CNS cells')
seu <- ScaleData(seu) %>% RunPCA(npcs = 50)
ElbowPlot(seu, ndims = 50)
seu <- RunUMAP(seu, dims = 1:20)

### Re-integrate
DefaultAssay(seu) <- 'RNA'
obj.list <- SplitObject(seu, split.by = 'patient')

# Keep only patients with >40 cells
for (i in 1:length(obj.list)) {
  if (dim(obj.list[[i]]@assays$RNA@counts)[2] < 40) {
    obj.list[i] <- NA
  }
}
# Remove NAs
obj.list <- obj.list[!is.na(obj.list)]

# Integration
anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:20)
seu <- IntegrateData(anchorset = anchors, dims = 1:20, k.weight = 40)
seu <- ScaleData(seu) %>% RunPCA() %>% RunUMAP(dims = 1:20)

# Plots
ifelse(!dir.exists(file.path(paste0('data/MBPM/global'))), dir.create(file.path(paste0('data/MBPM/global')), recursive = T), FALSE)
pdf(file = 'data/MBPM/global/plots_MBPM_cns_reintegrated.pdf')
DimPlot(seu, reduction = 'umap', label = F, group.by = 'sequencing', shuffle = T, raster = T, cols = colSCSN) + NoLegend()
DimPlot(seu, reduction = 'umap', label = F, group.by = 'sequencing', shuffle = T, raster = T, cols = colSCSN)

DimPlot(seu, reduction = 'umap', label = F, group.by = 'organ', shuffle = T, raster = T, cols = colBP) + NoLegend()
DimPlot(seu, reduction = 'umap', label = F, group.by = 'organ', shuffle = T, raster = T, cols = colBP)

DimPlot(seu, reduction = 'umap', label = F, group.by = 'patient', shuffle = T, raster = T)

DimPlot(seu, reduction = 'umap', label = F, group.by = 'cell_type_fine', shuffle = T, raster = T, repel = T, cols = colCNS) + NoLegend()
DimPlot(seu, reduction = 'umap', label = T, group.by = 'cell_type_fine', shuffle = T, raster = T, repel = T, cols = colCNS)
dev.off()


#### Stromal and epithelial cells #####
seu <- readRDS('data/MBPM/data_MBPM.rds')
seu <- subset(seu, cell_type_int %notin% c('Low-quality cells', 'Doublets', 'Contamination', 'Not determined'))
seu <- subset(seu, cell_type_main %in% c('Stromal cells', 'Epithelial cells', 'Endothelial cells'))
seu <- ScaleData(seu) %>% RunPCA(npcs = 50)
ElbowPlot(seu, ndims = 50)
seu <- RunUMAP(seu, dims = 1:30)

# Plot
pdf(file = 'data/MBPM/global/plots_MBPM_stromal_epithelial.pdf')
DimPlot(seu, reduction = 'umap', label = F, group.by = 'sequencing', shuffle = T, raster = T, cols = colSCSN) + NoLegend()
DimPlot(seu, reduction = 'umap', label = F, group.by = 'sequencing', shuffle = T, raster = T, cols = colSCSN)

DimPlot(seu, reduction = 'umap', label = F, group.by = 'organ', shuffle = T, raster = T, cols = colBP) + NoLegend()
DimPlot(seu, reduction = 'umap', label = F, group.by = 'organ', shuffle = T, raster = T, cols = colBP)

DimPlot(seu, reduction = 'umap', label = F, group.by = 'patient', shuffle = T, raster = T)

DimPlot(seu, reduction = 'umap', label = F, group.by = 'cell_type_fine', shuffle = T, raster = T, repel = T) + NoLegend()
DimPlot(seu, reduction = 'umap', label = T, group.by = 'cell_type_fine', shuffle = T, raster = T, repel = T)
dev.off()


#### DEGs non-immune non-tumor ####
seu <- readRDS('data/MBPM/data_MBPM.rds')
seu <- subset(seu, cell_type_int %notin% c('Low-quality cells', 'Doublets', 'Contamination', 'Not determined'))
seu <- subset(seu, cell_type_main %in% c('Stromal cells', 'CNS cells', 'Epithelial cells', 'Endothelial cells'))
sc <- subset(seu, sequencing == 'Single cells')
sn <- subset(seu, sequencing == 'Single nuclei')

# sc
Idents(sc) <- sc$cell_type_fine
markers_sc <- FindAllMarkers(sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = 'RNA')
ifelse(!dir.exists(file.path('data/cell_type_DEG/sc/nonimmune_nontumor/')), dir.create(file.path('data/cell_type_DEG/sc/nonimmune_nontumor/'), recursive = T), FALSE)
write.csv(markers_sc, 'data/cell_type_DEG/sc/nonimmune_nontumor/markers_MBPM_sc_type.csv', row.names = F)

# sn
Idents(sn) <- sn$cell_type_fine
markers_sn <- FindAllMarkers(sn, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay = 'RNA')
ifelse(!dir.exists(file.path('data/cell_type_DEG/sn/nonimmune_nontumor/')), dir.create(file.path('data/cell_type_DEG/sn/nonimmune_nontumor/'), recursive = T), FALSE)
write.csv(markers_sn, 'data/cell_type_DEG/sn/nonimmune_nontumor/markers_MBPM_sn_type.csv', row.names = F)


print(Sys.time())