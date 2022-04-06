#!/usr/bin/env Rscript

#### Integration analysis: LISI plots
#### Author: Jana Biermann, PhD

library(Seurat)
library(dplyr)
library(ggplot2)
library(gplots)
library(ggrastr)
library(viridis)
library(scales)

colBP <- c('#A80D11', '#008DB8')
colSCSN <- c('#E1AC24', '#288F56')

### Select one label
label <- 'MBPM_scn'
label <- 'MBPM_sn'
label <- 'MBPM_scn_tcell'
label <- 'MBPM_sn_tcell'


seu <- readRDS('data/MBPM/data_MBPM_scn.rds')
if (label == 'MBPM_sn') {
  seu <- subset(seu, sequencing == 'Single nuclei')
}

if (label == 'MBPM_scn_tcell') {
  seu <- subset(seu, cell_type_main == 'T/NK cells' & cell_type_fine != 'NK cells')
}

if (label == 'MBPM_sn_tcell') {
  seu <- subset(seu, sequencing == 'Single nuclei' & cell_type_main == 'T/NK cells' &
                  cell_type_fine != 'NK cells')
}

DefaultAssay(seu) <- 'RNA'
seu <- DietSeurat(seu, assays = 'RNA')

vars <- data.frame(patient = seu$patient, cell_type_main = seu$cell_type_main,
                   cell_type_int = seu$cell_type_int, cell_type_fine = seu$cell_type_fine,
                   sequencing = seu$sequencing, organ = seu$organ, 
                   row.names = rownames(seu@meta.data))


#### Combine LISI scores ####
res_raw_umap <- read.csv(paste0('data/integration_benchmarking/lisi_', label,
                                '_raw_umap.csv'), row.names = 1)
res_seurat_umap <- read.csv(paste0('data/integration_benchmarking/lisi_', label,
                                   '_seurat_umap.csv'), row.names = 1)
res_harmony_umap <- read.csv(paste0('data/integration_benchmarking/lisi_', label,
                                    '_harmony_umap.csv'), row.names = 1)
res_conos_umap <- read.csv(paste0('data/integration_benchmarking/lisi_', label,
                                  '_conos_umap.csv'), row.names = 1)

colnames(res_raw_umap) <- paste0('raw_umap_', colnames(res_raw_umap))
colnames(res_seurat_umap) <- paste0('seurat_umap_', colnames(res_seurat_umap))
colnames(res_harmony_umap) <- paste0('harmony_umap_', colnames(res_harmony_umap))
colnames(res_conos_umap) <- paste0('conos_umap_', colnames(res_conos_umap))

if (grepl('tcell', label)) {
  res_stacas_umap <- read.csv(paste0('data/integration_benchmarking/lisi_', label,
                                     '_stacas_umap.csv'), row.names = 1)
  colnames(res_stacas_umap) <- paste0('stacas_umap_', colnames(res_stacas_umap))
}

if (grepl('tcell', label)) {
  bc <- rownames(res_stacas_umap)
  vars <- vars[bc, ]
  df <- cbind.data.frame(res_raw_umap[bc, ], res_seurat_umap, res_harmony_umap[bc, ], 
                         res_conos_umap, res_stacas_umap)
  seu <- subset(seu, cells = bc)
} else {
  df <- cbind.data.frame(res_raw_umap, res_seurat_umap, res_harmony_umap,
                         res_conos_umap)
}

# Convert to long format
df_long <- reshape2::melt(df, value.name = 'LISI score')
temp <- unlist(lapply(names(vars), function(x) {
  rep(x, dim(seu@assays$RNA@counts)[2])
}))
df_long$category <- temp
df_long$tool <- substr(df_long$variable, 1, 3)
write.csv(df_long, paste0('data/integration_benchmarking/lisi_', label, '_combined.csv'),
          row.names = F)

# Plots
pdf(paste0('data/integration_benchmarking/plots_', label, '_lisi.pdf'))
ggplot(df_long[df_long$category == 'patient', ], 
       aes(x = reorder(variable, `LISI score`, FUN = median), 
           y = `LISI score`, fill = tool)) + 
  geom_boxplot() + 
  theme_classic() +
  ggtitle('iLISI score of UMAP coordinates (patient)') +
  xlab('Integration method') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        legend.key.size = unit(0.3, 'cm'))

ggplot(df_long[grepl('cell_type_fine', df_long$category), ], 
       aes(x = reorder(variable, `LISI score`, FUN = median), 
           y = `LISI score`, fill = tool)) + 
  geom_boxplot() +
  theme_classic() + 
  ggtitle('cLISI score of UMAP coordinates (cell_type_fine)') +
  xlab('Integration method') + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        legend.key.size = unit(0.3, 'cm'))
dev.off()
