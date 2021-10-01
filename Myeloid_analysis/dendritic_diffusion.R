#!/usr/bin/env Rscript

### title: Diffusion component analysis of dendritic cells (DCs)
### author: Jana Biermann, PhD

print(Sys.time())

library(Seurat)
library(destiny)
library(dplyr)
library(ggplot2)
library(gplots)
library(viridis)
library(scales)

colBP <- c('#A80D11', '#008DB8')
colSCSN <- c('#E1AC24', '#288F56')
colDC <- c('#DE8C00', '#F564E3', '#7CAE00', '#00B4F0', '#00C08B')


# Read-in integrated object and subset to monocytes and DCs
seu <- readRDS('data/MBPM/data_MBPM.rds')
seu <- subset(seu, cell_type_fine %in% c('Monocytes', 'cDC1', 'cDC2', 'DC3', 'mo-DCs'))

# Apply signatures
sigs <- read.csv('signatures/myeloid_signatures.csv', na.strings = '')
for (c in 1:ncol(sigs)) {
  sig <- as.character(na.omit(sigs[, c]))
  seu <- AddModuleScore(object = seu, features = list(sig), name = colnames(sigs[c]), assay = 'RNA', search = F)
}

# Generate diffusion map
es <- as.ExpressionSet(as.data.frame(t(seu@assays$integrated@data)))
es@phenoData@data <- seu@meta.data
dm <- DiffusionMap(es, verbose = T, n_pcs = 30)

# Plots
pdf('data/diffusion/myeloid/DC/plots_MBPM_myeloid_dc_dm.pdf')
par(mar = c(5.1, 4.1, 4.1, 2.1), xpd = TRUE)
palette(colSCSN)
plot(dm, col = as.factor(es@phenoData@data$sequencing), main = 'sequencing', pch = 20)
legend('bottomright', inset = c(-0.05, 0), legend = levels(as.factor(es@phenoData@data$sequencing)), 
       pch = 16, col = as.factor(levels(as.factor(es@phenoData@data$sequencing))), bty = 'n')

palette(colBP)
plot(dm, col = as.factor(es@phenoData@data$organ), main = 'organ', pch = 20)
legend('bottomright', inset = c(-0.05, 0), legend = levels(as.factor(es@phenoData@data$organ)), pch = 16, 
       col = as.factor(levels(as.factor(es@phenoData@data$organ))), bty = 'n')

palette(colDC)
plot(dm, col = as.factor(seu$cell_type_fine), main = 'cell_type_fine', pch = 20)
legend('bottomright', inset = c(-0.05, -0.2), legend = levels(as.factor(seu$cell_type_fine)), pch = 16, 
       col = as.factor(levels(as.factor(seu$cell_type_fine))), bty = 'n')

par(mar = c(5.1, 4.1, 4.1, 7), xpd = TRUE)
palette(viridis(100))
plot(dm, col = es@phenoData@data$DC3_zilionis1, main = 'DC3_zilionis1', pch = 20)
colorlegend(es@phenoData@data$DC3_zilionis1, viridis(length(es@phenoData@data$DC3_zilionis1)), 
            posx = c(0.88, 0.9), posy = c(0, 0.65))
dev.off()


print(Sys.time())