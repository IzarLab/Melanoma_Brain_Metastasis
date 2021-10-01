#!/usr/bin/env Rscript

### title: Analysis myeloid cells (UMAP, violin plots, diffusion map of monocytes + TAM MDMs) 
### author: Jana Biermann, PhD

print(Sys.time())

library(Seurat)
library(destiny)
library(dplyr)
library(ggplot2)
library(gplots)
library(viridis)
library(scales)
'%notin%' <- Negate('%in%')


colBP <- c('#A80D11', '#008DB8')
colSCSN <- c('#E1AC24', '#288F56')
colMM <- c('#00C08B', '#619CFF', '#FF64B0', '#B79F00', '#F8766D', '#C77CFF', '#00BFC4', '#00B4F0', '#DE8C00', 
           '#F564E3', '#7CAE00', '#00BA38')

# Read-in integrated object, subset and rerun workflow
seu_all <- readRDS('data/MBPM/data_MBPM.rds')
seu <- subset(seu_all, cell_type_main == 'Myeloid cells')
seu <- subset(seu, cell_type_int %notin% c('Low-quality cells', 'Doublets', 'Contamination', 'Not determined'))
seu <- ScaleData(seu) %>% RunPCA() %>% RunUMAP(dims = 1:50)


#### UMAPs ####
pdf('data/MBPM/myeloid/plots_myeloid_umap.pdf')
DimPlot(seu, group.by = 'cell_type_fine', label = T, raster = T, shuffle = T, cols = colMM)
DimPlot(seu, group.by = 'cell_type_fine', raster = T, shuffle = T, cols = colMM) + NoLegend()
DimPlot(seu, group.by = 'sequencing', raster = T, shuffle = T, cols = colSCSN)
DimPlot(seu, group.by = 'sequencing', raster = T, shuffle = T, cols = colSCSN) + NoLegend()
DimPlot(seu, group.by = 'organ', raster = T, shuffle = T, cols = colBP)
DimPlot(seu, group.by = 'organ', raster = T, shuffle = T, cols = colBP) + NoLegend()
dev.off()


#### Violin plots ####
genes_suppl <- c('S100A8', 'S100A9', 'S100A12', 'VCAN', 'FCN1', 'LYZ', 'MX1', 'IFI6', 'LAP3', 'OAS1', 
                 'CXCL10', 'LY6E', 'RNASE1', 'F13A1', 'STAB1', 'SLC40A1', 'MRC1', 'CD163', 'SELENOP', 'DAB2', 'SIGLEC1', 
                 'GPNMB', 'APOC1', 'IQGAP2', 'LIPA', 'FTL', 'FTH1', 'CD5L', 'TMSB10', 'SPP1', 'F13A1', 'RBPJ', 'NAMPT', 
                 'SNX9', 'DNM2', 'SLC16A10', 'SPRED1', 'TREM2', 'MERTK', 'NAV3', 'C3', 'TMEM119', 'P2RY12', 'THBD', 
                 'CD83', 'ITGAX', 'PPARG', 'CLEC9A', 'XCR1', 'CLNK', 'ENPP1', 'BCL11A', 'PLAC8', 'IL18R1', 'FLT3', 
                 'SLC38A1', 'RTN1', 'IDO1', 'IDO2', 'CD1C', 'FCER1A', 'IRF4', 'CLEC10A', 'CLEC4A', 'LAMP3', 'CCL19', 
                 'IDO1', 'BATF3', 'IRF8', 'CD80', 'CD274', 'PDCD1LG2', 'CD200', 'CCR7', 'FSCN1', 'CD40', 'MKI67', 
                 'TOP2A')

genes_main <- c('VCAN', 'FCN1', 'LYZ', 'MX1', 'IFI6', 'LY6E', 'RNASE1', 'MRC1', 'SELENOP', 'FTL', 'FTH1', 
                'DNM2', 'SLC16A10', 'TMEM119', 'P2RY12', 'THBD', 'CD83', 'CLEC9A', 'XCR1', 'CLNK', 'PLAC8', 'FLT3', 
                'IDO1', 'LAMP3', 'CCL19', 'CD200', 'MKI67', 'TOP2A')

pdf('data/MBPM/myeloid/plots_myeloid_markers_main.pdf', height = 9, width = 4)
VlnPlot(seu, features = genes_main, group.by = 'cell_type_fine', pt.size = 0, assay = 'RNA', stack = T, 
        flip = T) + NoLegend() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

pdf('data/MBPM/myeloid/plots_myeloid_markers_suppl.pdf', height = 14, width = 4)
VlnPlot(seu, features = genes13, group.by = 'cell_type_fine', pt.size = 0, assay = 'RNA', stack = T, 
        flip = T) + NoLegend() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()


#### Diffusion component analysis of monocytes + TAM MDMs ####
# Subset to monoctes + TAM MDMs
seu <- subset(seu, cell_type_fine %in% c('Monocytes', 'TAM MDM', 'TAM MDM FTL+', 'TAM MDM M1-like', 'TAM MDM M2-like'))

# Apply signatures
sigs <- read.csv('signatures/myeloid_signatures.csv', na.strings = '')
for (c in 1:ncol(sigs)) {
  sig <- as.character(na.omit(sigs[, c]))
  seu <- AddModuleScore(object = seu, features = list(sig), name = colnames(sigs[c]), assay = 'RNA', search = F)
}

# Generate diffusion map
es <- as.ExpressionSet(as.data.frame(t(seu@assays$integrated@data)))
es@phenoData@data <- seu@meta.data
dm <- DiffusionMap(es, verbose = T, n_pcs = 20)

# Plots
pdf('data/diffusion/myeloid/plots_MBPM_myeloid_dm.pdf')
par(mar = c(5.1, 4.1, 4.1, 2.1), xpd = TRUE)
palette(colSCSN)
plot(dm, col = as.factor(es@phenoData@data$sequencing), main = 'sequencing', pch = 20)
legend('bottomright', inset = c(-0.05, 0), legend = levels(as.factor(es@phenoData@data$sequencing)), 
       pch = 16, col = as.factor(levels(as.factor(es@phenoData@data$sequencing))), bty = 'n')

palette(colBP)
plot(dm, col = as.factor(es@phenoData@data$organ), main = 'organ', pch = 20)
legend('bottomright', inset = c(-0.05, 0), legend = levels(as.factor(es@phenoData@data$organ)), pch = 16, 
       col = as.factor(levels(as.factor(es@phenoData@data$organ))), bty = 'n')

palette(colMM)
plot(dm, col = as.factor(seu$cell_type_fine), main = 'cell_type_fine', pch = 20)
legend('bottomright', inset = c(-0.05, -0.2), legend = levels(as.factor(seu$cell_type_fine)), pch = 16, 
       col = as.factor(levels(as.factor(seu$cell_type_fine))), bty = 'n')
dev.off()

print(Sys.time())