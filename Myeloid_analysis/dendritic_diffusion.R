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
seu <- readRDS('data/MBPM/data_MBPM_scn.rds')
seu <- subset(seu, cell_type_fine %in% c('Monocytes', 'cDC1', 'cDC2', 'DC3'))

# Apply signatures
sigs <- read.csv('signatures/myeloid_signatures.csv', na.strings = '')
for (c in 1:ncol(sigs)) {
  sig <- as.character(na.omit(sigs[, c]))
  seu <- AddModuleScore(object = seu, features = list(sig), name = colnames(sigs[c]), 
                        assay = 'RNA', search = F)
}

# Generate diffusion map
es <- as.ExpressionSet(as.data.frame(t(seu@assays$integrated@data)))
es@phenoData@data <- seu@meta.data
dm <- DiffusionMap(es, verbose = T, n_pcs = 30)

# Plots
pdf('data/diffusion/myeloid/DC/plots_MBPM_myeloid_dc_dm.pdf')
set.seed(1)
dm_df<-as.data.frame(dm@eigenvectors)[sample(nrow(as.data.frame(dm@eigenvectors))),1:2]
ggplot(dm_df, aes(DC1, DC2, col = seu@meta.data[rownames(dm_df), 'cell_type_fine']))+
  geom_point(alpha = 1,shape = 16,size = 1.5) + theme_classic() + ggtitle('cell_type_fine')+
  theme(legend.title = element_blank())

ggplot(dm_df,aes(DC1,DC2,col=seu@meta.data[rownames(dm_df),'sequencing']))+
  geom_point(alpha = 1,shape=16,size=1.5) + theme_classic() +ggtitle('sequencing')+
  theme(legend.title = element_blank())

ggplot(dm_df,aes(DC1,DC2,col=seu@meta.data[rownames(dm_df),'organ']))+
  geom_point(alpha = 1,shape=16,size=1.5) + theme_classic() +ggtitle('organ')+
  theme(legend.title = element_blank())

ggplot(dm_df,aes(DC1,DC2,col=seu@meta.data[rownames(dm_df),'cDC11']))+
  geom_point(alpha = 1,shape=16,size=1.5) + theme_classic() +ggtitle('cDC11')+
  theme(legend.title = element_blank())+ scale_color_viridis()

ggplot(dm_df,aes(DC1,DC2,col=seu@meta.data[rownames(dm_df),'cDC21']))+
  geom_point(alpha = 1,shape=16,size=1.5) + theme_classic() +ggtitle('cDC21')+
  theme(legend.title = element_blank())+ scale_color_viridis()

ggplot(dm_df,aes(DC1,DC2,col=seu@meta.data[rownames(dm_df),'DC3_zilionis1']))+
  geom_point(alpha = 1,shape=16,size=1.5) + theme_classic() +ggtitle('DC3_zilionis1')+
  theme(legend.title = element_blank())+ scale_color_viridis()

ggplot(dm_df,aes(DC1,DC2,col=seu@meta.data[rownames(dm_df),'Monocytes1']))+
  geom_point(alpha = 1,shape=16,size=1.5) + theme_classic() +ggtitle('Monocytes1')+
  theme(legend.title = element_blank())+ scale_color_viridis()
dev.off()


print(Sys.time())