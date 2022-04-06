#!/usr/bin/env Rscript

### title: Diffusion component (D.C.) analysis and violin plots of B cells 
### author: Jana Biermann, PhD

library(Seurat)
library(destiny)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(gplots)
library(viridis)
library(scales)

colBP <- c('#A80D11', '#008DB8')
colSCSN <- c('#E1AC24', '#288F56')
colBcell <- c('#C7B122', '#C70E7B', '#22BAC7')

# Read-in object and subset to B cells
seu <- readRDS('data/MBPM/data_MBPM_scn.rds')
seu <- subset(seu, cell_type_fine %in% c('Activated B cells', 'Naïve B cells',
                                         'Plasma cells'))
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = 50)
ElbowPlot(seu, ndims = 50)
seu <- RunUMAP(seu, dims = 1:30)

#### Diffusion maps ####
# Generate diffusion map
es <- as.ExpressionSet(as.data.frame(t(seu@assays$integrated@data)))
es@phenoData@data <- seu@meta.data
dm <- DiffusionMap(es, verbose = T, n_pcs = 20)

# Diffusion pseudotime calculation
dpt <- DPT(dm)
es$pseudotime_dpt <- rank(dpt$dpt)
seu$pseudotime_dpt <- rank(dpt$dpt)
seu <- AddMetaData(seu, dm@eigenvectors, col.name = colnames(dm@eigenvectors))

saveRDS(dm, 'data/diffusion/bcells/data_MBPM_scn_bcells_dm.rds')

# Plots
ifelse(!dir.exists(file.path('data/diffusion/bcells')),
       dir.create(file.path('data/diffusion/bcells'), recursive = T),
       FALSE)
pdf('data/diffusion/bcells/plots_MBPM_scn_bcells_dm.pdf')

# UMAPs
DimPlot(seu, group.by = 'cell_type_fine', label = T, raster = T, shuffle = T, repel = T,
        cols = colBcell)
DimPlot(seu, group.by = 'cell_type_fine', raster = T, shuffle = T, cols = colBcell) +
  NoLegend()
DimPlot(seu, group.by = 'sequencing', raster = T, shuffle = T, cols = colSCSN)
DimPlot(seu, group.by = 'sequencing', raster = T, shuffle = T, cols = colSCSN) +
  NoLegend()
DimPlot(seu, group.by = 'organ', raster = T, shuffle = T, cols = colBP)
DimPlot(seu, group.by = 'organ', raster = T, shuffle = T, cols = colBP) + NoLegend()

# 3D diffusion maps
par(mar = c(5.1, 4.1, 4.1, 2.1), xpd = TRUE)
palette(colSCSN)
tmp_var <- seu$sequencing
plot(dm, col = as.factor(tmp_var), main = 'sequencing', pch = 20)
legend('bottomright', inset = c(0.1, -0.2), legend = levels(as.factor(tmp_var)), pch = 16,
       col = as.factor(levels(as.factor(tmp_var))), bty = 'n')

palette(colBP)
tmp_var <- seu$organ
plot(dm, col = as.factor(tmp_var), main = 'organ', pch = 20)
legend('bottomright', inset = c(0.1, -0.2), legend = levels(as.factor(tmp_var)), pch = 16,
       col = as.factor(levels(as.factor(tmp_var))), bty = 'n')

palette(colBcell)
tmp_var <- seu$cell_type_fine
plot(dm, col = as.factor(tmp_var), main = 'cell_type_fine', pch = 20)
legend('bottomright', inset = c(-0.05, -0.2), legend = levels(as.factor(tmp_var)),
       pch = 16, col = as.factor(levels(as.factor(tmp_var))), bty = 'n')

# 2D diffusion maps
set.seed(1)
dm_df <- as.data.frame(dm@eigenvectors)[sample(nrow(as.data.frame(dm@eigenvectors))), 1:2]
ggplot(dm_df, aes(DC1, DC2, col = seu@meta.data[rownames(dm_df), 'cell_type_fine'])) +
  geom_point(alpha = 1, shape = 16, size = 1.5) + theme_classic() +
  ggtitle('cell_type_fine') + theme(legend.title = element_blank()) +
  scale_color_manual(values = colBcell)

ggplot(dm_df, aes(DC1, DC2, col = seu@meta.data[rownames(dm_df), 'sequencing'])) +
  geom_point(alpha = 1, shape = 16, size = 1.5) + theme_classic() +
  ggtitle('sequencing') + theme(legend.title = element_blank()) +
  scale_color_manual(values = colSCSN)

ggplot(dm_df, aes(DC1, DC2, col = seu@meta.data[rownames(dm_df), 'organ'])) +
  geom_point(alpha = 1, shape = 16, size = 1.5) + theme_classic() + ggtitle('organ') +
  theme(legend.title = element_blank()) + scale_color_manual(values = colBP)

ggplot(dm_df, aes(DC1, DC2, col = seu@meta.data[rownames(dm_df), 'patient'])) +
  geom_point(alpha = 1, shape = 16, size = 1.5) + theme_classic() + ggtitle('patient') +
  theme(legend.title = element_blank())

ggplot(dm_df, aes(DC1, DC2, col = seu@meta.data[rownames(dm_df), 'pseudotime_dpt'])) +
  geom_point(alpha = 1, shape = 16, size = 1.5) + theme_classic() +
  ggtitle('pseudotime_dpt') + theme(legend.title = element_blank()) +
  scale_color_viridis()

dev.off()


#### Violin plots #####
seu <- readRDS('data/MBPM/data_MBPM_scn.rds')
seu <- subset(seu, cell_type_fine %in% c('Activated B cells', 'Naïve B cells',
                                         'Plasma cells'))
sn <- subset(seu, sequencing == 'Single nuclei')

genes <- c('BANK1', 'BLK', # activated B
           'LAIR1', 'TCL1A', #naive
           'SDC1', 'CD38', 'PRDM1') # plasma

ifelse(!dir.exists(file.path('data/MBPM/bcells')),
       dir.create(file.path('data/MBPM/bcells'), recursive = T), FALSE)
pdf('data/MBPM/bcells/plots_MBPM_bcells_violin.pdf')
VlnPlot(sn, features = genes, group.by = 'cell_type_fine', pt.size = 0,
        assay = 'RNA', stack = T, flip = T, split.by = 'organ', cols = colBP) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = 'bottom') + ggtitle('sn')
dev.off()
