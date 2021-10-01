#!/usr/bin/env Rscript

### title: Generating averaged heatmap with Viper output 
### author: Jana Biermann, PhD

library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
'%notin%' <- Negate('%in%')

colBP <- c('#A80D11', '#008DB8')
colSCSN <- c('#E1AC24', '#288F56')


# Read-in integrated object
seu <- readRDS('data/MBPM/data_MBPM.rds')

# Read-in viper output and modify
viper <- read.csv('data/viper/MBM_MPM_sn_all_internal_reference_vpmat.csv')
viper <- as.data.frame(t(viper))
colnames(viper) <- viper[1, ]
viper <- viper[-1, ]

vp_mat <- as.matrix(viper)
vp3 <- apply(vp_mat, 2, as.numeric)
rownames(vp3) <- rownames(vp_mat)

# Generate seu assay obj with Viper matrix
seu <- subset(seu, barcode_all %in% colnames(viper))
seu@assays$viper <- CreateAssayObject(data = vp3)
seu@assays$viper@scale.data <- vp3
DefaultAssay(seu) <- 'viper'

# Generate heatmap
avg <- AverageExpression(seu, assays = 'viper', return.seurat = T, group.by = 'organ', slot = 'scale.data')
avg$organ <- rownames(avg@meta.data)
tops <- rownames(vp3)

pdf(paste0('data/viper/plots_viper_MBPM_sn_tum_organ.pdf'))
DoHeatmap(avg, features = tops, group.by = 'organ', raster = T, assay = 'viper', draw.lines = F, group.colors = colBP) + 
  scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name = 'RdBu')))(100))
DoHeatmap(avg, features = tops, group.by = 'organ', raster = T, assay = 'viper', draw.lines = F, size = 0, group.colors = colBP) + 
  scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name = 'RdBu')))(100)) + 
  theme(plot.margin = margin(0, 9, 0, 0, 'cm'))
dev.off()
