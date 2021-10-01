#!/usr/bin/env Rscript

### title: Diffusion Component (D.C.) analysis of B cells 
### author: Jana Biermann, PhD

print(Sys.time())

library(Seurat)
library(destiny)
library(dplyr)
library(ggplot2)
library(gplots)
library(viridis)
library(scales)


# Read-in object and subset to B cells
seu <- readRDS('data/MBPM/data_MBPM.rds')
seu <- subset(seu, cell_type_fine %in% c('Activated B cells', 'Naïve B cells', 'Plasma cells'))
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = 50)
ElbowPlot(seu, ndims = 50)
seu <- RunUMAP(seu, dims = 1:30)

# Generate diffusion map
es <- as.ExpressionSet(as.data.frame(t(seu@assays$integrated@data)))
es@phenoData@data <- seu@meta.data
dm <- DiffusionMap(es, verbose = T, n_pcs = 30)

# Diffusion pseudotime calculation
dpt <- DPT(dm)
es$pseudotime_dpt <- rank(dpt$dpt)
seu$pseudotime_dpt <- rank(dpt$dpt)
seu <- AddMetaData(seu, dm@eigenvectors, col.name = colnames(dm@eigenvectors))

# Plots
ifelse(!dir.exists(file.path('data/diffusion/bcells')), dir.create(file.path('data/diffusion/bcells'), recursive = T), FALSE)
pdf('data/diffusion/bcells/plots_MBPM_bcells_dm.pdf')
set.seed(1)
dm_df <- as.data.frame(dm@eigenvectors)[sample(nrow(as.data.frame(dm@eigenvectors))), 1:2]
ggplot(dm_df, aes(DC1, DC2, col = seu@meta.data[rownames(dm_df), 'cell_type_fine'])) + 
  geom_point(alpha = 1, shape = 16, size = 1.5) + theme_classic() + ggtitle('cell_type_fine') + 
  theme(legend.title = element_blank())

ggplot(dm_df, aes(DC1, DC2, col = seu@meta.data[rownames(dm_df), 'sequencing'])) + 
  geom_point(alpha = 1, shape = 16, size = 1.5) + theme_classic() + ggtitle('sequencing') + 
  theme(legend.title = element_blank())

ggplot(dm_df, aes(DC1, DC2, col = seu@meta.data[rownames(dm_df), 'organ'])) + geom_point(alpha = 1, shape = 16, 
                                                                                         size = 1.5) + theme_classic() + ggtitle('organ') + theme(legend.title = element_blank())

ggplot(dm_df, aes(DC1, DC2, col = seu@meta.data[rownames(dm_df), 'patient'])) + 
  geom_point(alpha = 1, shape = 16, size = 1.5) + theme_classic() + ggtitle('patient') + 
  theme(legend.title = element_blank())

ggplot(dm_df, aes(DC1, DC2, col = -seu@meta.data[rownames(dm_df), 'pseudotime_dpt'])) + 
  geom_point(alpha = 1, shape = 16, size = 1.5) + theme_classic() + ggtitle('pseudotime_dpt') + 
  theme(legend.title = element_blank()) + scale_color_viridis()
dev.off()


print(Sys.time())