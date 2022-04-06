#!/usr/bin/env Rscript

### title: Analysis of non-immune and non-tumor cell types (UMAPs, DEG)
### author: Jana Biermann, PhD

library(Seurat)
library(dplyr)
library(ggplot2)
library(gplots)
'%notin%' <- Negate('%in%')

colBP <- c('#A80D11', '#008DB8')
colSCSN <- c('#E1AC24', '#288F56')
colCNS <- c('#92D84F', '#7473A6', '#FF8C8D')
colStromal <- c('#800080', '#008080', '#FFD700', '#F08080')

#### Central nervous system (CNS) cells ####
seu <- readRDS('data/MBPM/data_MBPM_scn.rds')
seu <- subset(seu, cell_type_int %notin% c('Low-quality cells', 'Doublets',
                                           'Contamination', 'Undetermined'))
seu <- subset(seu, cell_type_main == 'CNS cells')
seu <- ScaleData(seu) %>%
  RunPCA(npcs = 50)
ElbowPlot(seu, ndims = 50)
seu <- RunUMAP(seu, dims = 1:20)

# Plots
ifelse(!dir.exists(file.path(paste0('data/MBPM/global'))),
       dir.create(file.path(paste0('data/MBPM/global')), recursive = T),
       FALSE)
pdf(file = 'data/MBPM/global/plots_MBM_scn_cns.pdf')
DimPlot(seu, reduction = 'umap', label = F, group.by = 'sequencing', shuffle = T,
        raster = T, cols = colSCSN) + NoLegend()
DimPlot(seu, reduction = 'umap', label = F, group.by = 'sequencing', shuffle = T,
        raster = T, cols = colSCSN)
DimPlot(seu, reduction = 'umap', label = F, group.by = 'organ', shuffle = T, raster = T,
        cols = colBP) + NoLegend()
DimPlot(seu, reduction = 'umap', label = F, group.by = 'organ', shuffle = T, raster = T,
        cols = colBP)
DimPlot(seu, reduction = 'umap', label = F, group.by = 'patient', shuffle = T, raster = T)
DimPlot(seu, reduction = 'umap', label = F, group.by = 'cell_type_fine', shuffle = T,
        raster = T, repel = T, cols = colCNS) + NoLegend()
DimPlot(seu, reduction = 'umap', label = T, group.by = 'cell_type_fine', shuffle = T,
        raster = T, repel = T, cols = colCNS)
dev.off()


#### Stromal and epithelial cells #####
seu <- readRDS('data/MBPM/data_MBPM_scn.rds')
seu <- subset(seu, cell_type_int %notin% c('Low-quality cells', 'Doublets',
                                           'Contamination', 'Undetermined'))
seu <- subset(seu, cell_type_main %in% c('Stromal cells', 'Epithelial cells',
                                         'Endothelial cells'))
seu <- ScaleData(seu) %>%
  RunPCA(npcs = 50)
ElbowPlot(seu, ndims = 50)
seu <- RunUMAP(seu, dims = 1:20)

# Plot
pdf(file = 'data/MBPM/global/plots_MBPM_scn_stromal_epithelial.pdf')
DimPlot(seu, reduction = 'umap', label = F, group.by = 'sequencing', shuffle = T,
        raster = T, cols = colSCSN) + NoLegend()
DimPlot(seu, reduction = 'umap', label = F, group.by = 'sequencing', shuffle = T,
        raster = T, cols = colSCSN)
DimPlot(seu, reduction = 'umap', label = F, group.by = 'organ', shuffle = T, raster = T,
        cols = colBP) + NoLegend()
DimPlot(seu, reduction = 'umap', label = F, group.by = 'organ', shuffle = T, raster = T,
        cols = colBP)
DimPlot(seu, reduction = 'umap', label = F, group.by = 'patient', shuffle = T, raster = T)
DimPlot(seu, reduction = 'umap', label = F, group.by = 'cell_type_fine', shuffle = T,
        raster = T, repel = T, cols = colStromal) + NoLegend()
DimPlot(seu, reduction = 'umap', label = T, group.by = 'cell_type_fine', shuffle = T,
        raster = T, repel = T, cols = colStromal)
dev.off()


#### DEGs non-immune non-tumor ####
seu <- readRDS('data/MBPM/data_MBPM_scn.rds')
seu <- subset(seu, cell_type_int %notin% c('Low-quality cells', 'Doublets',
                                           'Contamination', 'Not determined'))
seu <- subset(seu, cell_type_main %in% c('Stromal cells', 'CNS cells', 'Epithelial cells',
                                         'Endothelial cells'))
sc <- subset(seu, sequencing == 'Single cells')
sn <- subset(seu, sequencing == 'Single nuclei' & organ == 'Brain')
snP <- subset(seu, sequencing == 'Single nuclei' & organ == 'Peripheral')

# MBM_sc
Idents(sc) <- sc$cell_type_fine
markers_sc <- FindAllMarkers(sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,
                             assay = 'RNA')
ifelse(!dir.exists(file.path('data/cell_type_DEG/MBM_sc/nonimmune_nontumor/')),
       dir.create(file.path('data/cell_type_DEG/MBM_sc/nonimmune_nontumor/'), recursive = T),
       FALSE)
write.csv(markers_sc,
          'data/cell_type_DEG/MBM_sc/nonimmune_nontumor/markers_MBM_sc_celltype.csv',
          row.names = F)

# MBM_sn
Idents(sn) <- sn$cell_type_fine
markers_sn <- FindAllMarkers(sn, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,
                             assay = 'RNA')
ifelse(!dir.exists(file.path('data/cell_type_DEG/MBM_sn/nonimmune_nontumor/')),
       dir.create(file.path('data/cell_type_DEG/MBM_sn/nonimmune_nontumor/'), recursive = T),
       FALSE)
write.csv(markers_sn,
          'data/cell_type_DEG/MBM_sn/nonimmune_nontumor/markers_MBM_sn_celltype.csv',
          row.names = F)

# MPM_sn
Idents(snP) <- snP$cell_type_fine
markers_snP <- FindAllMarkers(snP, only.pos = TRUE, min.pct = 0.25,
                              logfc.threshold = 0.25, assay = 'RNA')
ifelse(!dir.exists(file.path('data/cell_type_DEG/MPM_sn/nonimmune_nontumor/')),
       dir.create(file.path('data/cell_type_DEG/MPM_sn/nonimmune_nontumor/'), recursive = T),
       FALSE)
write.csv(markers_snP,
          'data/cell_type_DEG/MPM_sn/nonimmune_nontumor/markers_MPM_sn_celltype.csv',
          row.names = F)
