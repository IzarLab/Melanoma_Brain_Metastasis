#!/usr/bin/env Rscript

### title: T cell analysis (UMAPs, violins, DC) 
### author: Jana Biermann, PhD

print(Sys.time())

library(Seurat)
library(destiny)
library(dplyr)
library(ggplot2)
library(gplots)
library(viridis)
library(scales)
library(patchwork)
library(reshape2)
library(DropletUtils)
'%notin%' <- Negate('%in%')

colBP <- c('#A80D11', '#008DB8')
colSCSN <- c('#E1AC24', '#288F56')


#### UMAPs #####
seu <- readRDS('data/MBPM/data_MBPM.rds')
seu <- subset(seu, cell_type_main == 'T/NK cells' & 
                cell_type_int %notin% c('Low-quality cells', 'Doublets', 'Contamination', 'Not determined'))

ifelse(!dir.exists(file.path('data/MBPM/tcells')), dir.create(file.path('data/MBPM/tcells'), recursive = T),FALSE)
pdf(file = 'data/MBPM/tcells/plots_MBPM_tcells.pdf')
DimPlot(seu, reduction = 'umap', label = F, group.by = 'patient', shuffle = T, raster = T)
DimPlot(seu, reduction = 'umap', label = F, group.by = 'sequencing', shuffle = T, raster = T, cols = colSCSN)
DimPlot(seu, reduction = 'umap', label = F, group.by = 'sequencing', shuffle = T, raster = T, cols = colSCSN) + NoLegend()
DimPlot(seu, reduction = 'umap', label = F, group.by = 'organ', shuffle = T, raster = T, cols = colBP)
DimPlot(seu, reduction = 'umap', label = F, group.by = 'organ', shuffle = T, raster = T, cols = colBP) + NoLegend()
FeaturePlot(seu, features = c('frequency'), order = T, raster = T) + 
  scale_color_viridis(option = 'A', direction = -1, na.value = 'lightgray')
DimPlot(seu, reduction = 'umap', label = F, group.by = 'cell_type_fine', shuffle = T, raster = T) + NoLegend()
DimPlot(seu, reduction = 'umap', label = T, group.by = 'cell_type_fine', shuffle = T, raster = T, repel = T)
dev.off()


#### Violin plots #####
genes <- c('CD4', 'CD8A', 'MKI67', 'NCAM1', 'GNLY', 'KLRF1', 'KLRG1', 'TOX', 'TOX2', 'CD2', 'CD28', 'CTLA4', 
           'HAVCR2', 'LAG3', 'ICOS', 'PDCD1', 'TIGIT', 'GZMA', 'GZMB', 'GZMK', 'CCL3', 'CCL5', 'IFNG', 'TCF7', 
           'SELL', 'CCR7', 'LEF1', 'IL7R', 'BATF', 'CXCL13', 'FOXP3', 'IL2RA', 'CD69', 'CXCR6', 'ITGAE', 'FABP5')

pdf('data/MBPM/tcells/plots_tcell_markers.pdf', height = 10, width = 5)
VlnPlot(seu, features = genes, group.by = 'cell_type_fine', pt.size = 0, assay = 'RNA', stack = T, flip = T) + 
  NoLegend() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

# Subset to sn tox
sn_tox <- subset(seu, cell_type_fine == 'CD8+ T cells TOX+' & sequencing == 'Single nuclei')
genes_tox <- c('ITGAE', 'CD69', 'ENTPD1', 'LAG3', 'SLAMF6')

pdf('data/MBPM/tcells/plots_tcell_markers_tox.pdf')
VlnPlot(sn_tox, features = genes_tox, group.by = 'cell_type_fine', pt.size = 0, assay = 'RNA', stack = T, 
        flip = T, split.by = 'organ', split.plot = F, cols = colBP) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()


#### DC analysis of CD8+ T cells ####
seu <- subset(seu, cell_type_fine %in% c('CD8+ T cells TOX+', 'CD8+ T cells TCF7+'))
seu <- ScaleData(seu)
seu$clone_size <- ifelse(seu$frequency < 2, 'Non-expanded', 'Expanded')

# Apply signatures
sigs <- read.csv('signatures/Tcellsignatures.csv', na.strings = c('', 'NA'))
for (c in 1:ncol(sigs)) {
  sig <- as.character(na.omit(sigs[, c]))
  seu <- AddModuleScore(object = seu, features = list(sig), name = colnames(sigs[c]), assay = 'RNA', search = F)
}

# Generate diffusion map
es <- as.ExpressionSet(as.data.frame(t(seu@assays$integrated@data)))
es@phenoData@data <- seu@meta.data
dm <- DiffusionMap(es, verbose = T, n_pcs = 30)

# Plots
ifelse(!dir.exists(file.path('data/diffusion/tcells')), dir.create(file.path('data/diffusion/tcells'), recursive = T), FALSE)
pdf('data/diffusion/tcells/plots_MBPM_tcells_dm_cd8.pdf')
par(mar = c(5.1, 4.1, 4.1, 2.1), xpd = TRUE)
palette(colSCSN)
plot(dm, col = as.factor(es@phenoData@data$sequencing), main = 'sequencing', pch = 20)
legend('bottomright', inset = c(-0.05, 0), legend = levels(as.factor(es@phenoData@data$sequencing)), 
       pch = 16, col = as.factor(levels(as.factor(es@phenoData@data$sequencing))), bty = 'n')

palette(hue_pal()(7)[2:3])
plot(dm, col = as.factor(es@phenoData@data$cell_type_fine), main = 'cell_type_fine', pch = 20)
legend('bottomright', inset = c(-0.05, -0.1), legend = levels(as.factor(es@phenoData@data$cell_type_fine)), 
       pch = 16, col = as.factor(levels(as.factor(es@phenoData@data$cell_type_fine))), bty = 'n')

palette(colBP)
plot(dm, col = as.factor(es@phenoData@data$organ), main = 'organ', pch = 20)
legend('bottomright', inset = c(-0.05, 0), legend = levels(as.factor(es@phenoData@data$organ)), pch = 16, 
       col = as.factor(levels(as.factor(es@phenoData@data$organ))), bty = 'n')

palette('ggplot2')
plot(dm, col = as.factor(seu$clone_size), main = 'clone_size', pch = 20)
legend('bottomright', inset = c(-0.01, -0.1), legend = levels(as.factor(seu$clone_size)), pch = 16, 
       col = as.factor(levels(as.factor(seu$clone_size))), bty = 'n')

par(mar = c(5.1, 4.1, 4.1, 7), xpd = TRUE)
palette(viridis(100))
plot(dm, col = es@phenoData@data$Exhaustion_Terminal_differentiation1, 
     main = 'Exhaustion_Terminal_differentiation1', pch = 20)
colorlegend(es@phenoData@data$Exhaustion_Terminal_differentiation1, 
            viridis(length(es@phenoData@data$Exhaustion_Terminal_differentiation1)), 
            posx = c(0.88, 0.9), posy = c(0, 0.65))
dev.off()


print(Sys.time())