#!/usr/bin/env Rscript

### title: T cell analysis (UMAPs, violins, diffusion maps) 
### author: Jana Biermann, PhD

library(Seurat)
library(destiny)
library(dplyr)
library(ggplot2)
library(gplots)
library(viridis)
library(scales)
library(patchwork)
library(reshape2)
'%notin%' <- Negate('%in%')

colBP <- c('#A80D11', '#008DB8')
colSCSN <- c('#E1AC24', '#288F56')
colExp<-c('#C82F71','#7473A6')


#### UMAPs MBPM_sn #####
seu_all<-readRDS('data/MBPM/data_MBPM_scn.rds')
seu<-subset(seu_all,cell_type_main=='T/NK cells' &
              sequencing=='Single nuclei' &
              cell_type_int %notin% c('Low-quality cells','Doublets','Contamination',
                                      'Undetermined','Cycling cells'))

# Reintegrate
DefaultAssay(seu)<-'RNA'
# Split 
obj.list <- SplitObject(seu, split.by = "patient")
# Keep only patient with >50 cells
for(i in 1:length(obj.list)){
  if(dim(obj.list[[i]]@assays$RNA@counts)[2]<50){
    obj.list[i]<-NA
  }else{
    obj.list[[i]]<- NormalizeData(obj.list[[i]])
    obj.list[[i]]<- FindVariableFeatures(obj.list[[i]])
  }
}
# Remove NAs
obj.list<-obj.list[!is.na(obj.list)]

# Integration
anchors <- FindIntegrationAnchors(object.list = obj.list)
seu <- IntegrateData(anchorset = anchors, k.weight = 50)

seu <- ScaleData(seu) %>% RunPCA() 
ElbowPlot(seu,ndims = 50)
seu <- RunUMAP(seu,dims = 1:20)


# Get module scores for sigs
sigs<-read.csv('~/signatures/Tcellsignatures.csv',na.strings = '')
for(c in 1:ncol(sigs)){
  seu<-AddModuleScore(object = seu,features = list(na.omit(sigs[,c])),
                      name = colnames(sigs)[c],assay = 'RNA',search=F)
}


ifelse(!dir.exists(file.path(paste0('data/MBPM/tcells'))), 
       dir.create(file.path(paste0('data/MBPM/tcells')),recursive = T), FALSE)
saveRDS(seu,"data/MBPM/tcells/data_MBPM_sn_tcells_reint.rds")

pdf(file = 'data/MBPM/tcells/plots_MBPM_sn_tcells_reint.pdf')
DimPlot(seu, reduction = 'umap', label = F, group.by = 'patient', shuffle = T, raster = T)
DimPlot(seu, reduction = 'umap', label = F, group.by = 'patient', shuffle = T, raster = T)+NoLegend()
DimPlot(seu, reduction = 'umap', label = F, group.by = 'sequencing', shuffle = T, raster = T, cols = colSCSN)
DimPlot(seu, reduction = 'umap', label = F, group.by = 'sequencing', shuffle = T, raster = T, cols = colSCSN) + NoLegend()
DimPlot(seu, reduction = 'umap', label = F, group.by = 'organ', shuffle = T, raster = T, cols = colBP)
DimPlot(seu, reduction = 'umap', label = F, group.by = 'organ', shuffle = T, raster = T, cols = colBP) + NoLegend()
DimPlot(seu, reduction = 'umap', label = F, group.by = 'cell_type_fine', shuffle = T, raster = T) + NoLegend()
DimPlot(seu, reduction = 'umap', label = T, group.by = 'cell_type_fine', shuffle = T, raster = T, repel = T)
dev.off()


#### UMAPs MBM_sc #####
seu<-readRDS('data/cell_type_DEG/MBM_sc/tcells/data_MBM_sc_tcells.rds')
seu<-subset(seu,cell_type_fine %notin% c('Low-quality cells','Doublets','Contamination',
                                         'Undetermined','Cycling cells','Myeloid doublets'))

seu <- ScaleData(seu)
seu <- RunPCA(seu,npcs = 50)
ElbowPlot(seu,ndims = 50)
seu <- RunUMAP(seu, dims = 1:30)


ifelse(!dir.exists(file.path('data/MBPM/tcells')), 
       dir.create(file.path('data/MBPM/tcells'), recursive = T),FALSE)
pdf(file = 'data/MBPM/tcells/plots_MBM_sc_tcells.pdf')
DimPlot(seu, reduction = 'umap', label = F, group.by = 'patient', shuffle = T, raster = T)
DimPlot(seu, reduction = 'umap', label = F, group.by = 'sequencing', shuffle = T, raster = T, cols = colSCSN)
DimPlot(seu, reduction = 'umap', label = F, group.by = 'sequencing', shuffle = T, raster = T, cols = colSCSN) + NoLegend()
DimPlot(seu, reduction = 'umap', label = F, group.by = 'organ', shuffle = T, raster = T, cols = colBP)
DimPlot(seu, reduction = 'umap', label = F, group.by = 'organ', shuffle = T, raster = T, cols = colBP) + NoLegend()
DimPlot(seu, reduction = 'umap', label = F, group.by = 'cell_type_fine', shuffle = T, raster = T) + NoLegend()
DimPlot(seu, reduction = 'umap', label = T, group.by = 'cell_type_fine', shuffle = T, raster = T, repel = T)
dev.off()


#### UMAPs MBM_sn #####
seu<-readRDS('data/cell_type_DEG/MBM_sn/tcells/data_MBM_sn_tcells.rds')
seu<-subset(seu,cell_type_fine %notin% c('Low-quality cells','Doublets',
                                         'Contamination','Undetermined'))

seu <- ScaleData(seu)
seu <- RunPCA(seu,npcs = 50)
ElbowPlot(seu,ndims = 50)
seu <- RunUMAP(seu, dims = 1:20)


ifelse(!dir.exists(file.path('data/MBPM/tcells')), 
       dir.create(file.path('data/MBPM/tcells'), recursive = T),FALSE)
pdf(file = 'data/MBPM/tcells/plots_MBM_sn_tcells.pdf')
DimPlot(seu, reduction = 'umap', label = F, group.by = 'patient', shuffle = T, raster = T)
DimPlot(seu, reduction = 'umap', label = F, group.by = 'sequencing', shuffle = T, raster = T, cols = colSCSN)
DimPlot(seu, reduction = 'umap', label = F, group.by = 'sequencing', shuffle = T, raster = T, cols = colSCSN) + NoLegend()
DimPlot(seu, reduction = 'umap', label = F, group.by = 'organ', shuffle = T, raster = T, cols = colBP)
DimPlot(seu, reduction = 'umap', label = F, group.by = 'organ', shuffle = T, raster = T, cols = colBP) + NoLegend()
DimPlot(seu, reduction = 'umap', label = F, group.by = 'cell_type_fine', shuffle = T, raster = T) + NoLegend()
DimPlot(seu, reduction = 'umap', label = T, group.by = 'cell_type_fine', shuffle = T, raster = T, repel = T)
dev.off()


#### UMAPs MPM_sn #####
seu<-readRDS('data/cell_type_DEG/MPM_sn/tcells/data_MPM_sn_tcells.rds')
seu<-subset(seu,cell_type_fine %notin% c('Low-quality cells','Doublets','Contamination',
                                         'Undetermined','Cycling cells'))

seu <- ScaleData(seu)
seu <- RunPCA(seu,npcs = 50)
ElbowPlot(seu,ndims = 50)
seu <- RunUMAP(seu, dims = 1:20)


ifelse(!dir.exists(file.path('data/MBPM/tcells')), 
       dir.create(file.path('data/MBPM/tcells'), recursive = T),FALSE)
pdf(file = 'data/MBPM/tcells/plots_MPM_sn_tcells.pdf')
DimPlot(seu, reduction = 'umap', label = F, group.by = 'patient', shuffle = T, raster = T)
DimPlot(seu, reduction = 'umap', label = F, group.by = 'sequencing', shuffle = T, raster = T, cols = colSCSN)
DimPlot(seu, reduction = 'umap', label = F, group.by = 'sequencing', shuffle = T, raster = T, cols = colSCSN) + NoLegend()
DimPlot(seu, reduction = 'umap', label = F, group.by = 'organ', shuffle = T, raster = T, cols = colBP)
DimPlot(seu, reduction = 'umap', label = F, group.by = 'organ', shuffle = T, raster = T, cols = colBP) + NoLegend()
DimPlot(seu, reduction = 'umap', label = F, group.by = 'cell_type_fine', shuffle = T, raster = T) + NoLegend()
DimPlot(seu, reduction = 'umap', label = T, group.by = 'cell_type_fine', shuffle = T, raster = T, repel = T)
dev.off()


#### Violin plots #####
seu_all<-readRDS('data/MBPM/data_MBPM_scn.rds')
seu<-subset(seu_all,cell_type_main=='T/NK cells' &
              cell_type_int %notin% c('Low-quality cells','Doublets','Contamination',
                                      'Undetermined','Cycling cells'))

sn<-subset(seu,sequencing=='Single nuclei')
sc<-subset(seu,sequencing=='Single cells')

genes <- c('CD4','IL7R','SELL', 'CCR7', #CD4
           'TCF7','LEF1', # cd8 TCF
           'CD8A', 'TOX', # cd8 tox
           'LAG3','HAVCR2','PDCD1','ENTPD1', 'TIGIT', #immune check points
           'GZMA', # effector genes
           'CCL4','CCL5', # chemokines
           'NCAM1', 'GNLY', 'KLRF1', # NK
           'CD28', 'CTLA4',  
           'BATF', 'TOX2', #tfh
           'FOXP3', 'IL2RA') #treg


pdf('data/MBPM/tcells/plots_tcell_markers.pdf', height = 10, width = 5)
VlnPlot(sc, features = genes, group.by = 'cell_type_fine', pt.size = 0, assay = 'RNA', 
        stack = T, flip = T) + 
  NoLegend() + 
  ggtitle('sc')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
VlnPlot(sn, features = genes, group.by = 'cell_type_fine', pt.size = 0, assay = 'RNA', 
        stack = T, flip = T,split.by = 'organ',cols = colBP) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = 'bottom')+
  ggtitle('sn')
dev.off()

# Subset to sn TOX+
sn_tox <- subset(seu, cell_type_fine == 'CD8+ T cells TOX+' & sequencing == 'Single nuclei')
genes_BvsP<-c('HLA-A','HLA-B','HLA-C','B2M')
genes_BvsP2<-c('SLAMF6','TIGIT','LAG3','HAVCR2','PDCD1')

pdf('data/MBPM/tcells/plots_tcell_markers_tox.pdf')
VlnPlot(sn_tox, features = genes_BvsP, group.by = 'cell_type_fine', 
        pt.size = 0, assay = 'RNA', stack = T, 
        flip = T, split.by = 'organ', split.plot = F, cols = colBP) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
VlnPlot(sn_tox, features = genes_BvsP2, group.by = 'cell_type_fine', 
        pt.size = 0, assay = 'RNA', stack = T, 
        flip = T, split.by = 'organ', split.plot = F, cols = colBP) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()


#### D.C. analysis of CD8+ T cells MBM_sc ####
seu_all<-readRDS('data/cell_type_DEG/MBM_sc/tcells/data_MBM_sc_tcells.rds')
seu <- subset(seu_all, cell_type_fine %in% c('CD8+ T cells TOX+', 'CD8+ T cells TCF7+'))
seu$clone_size <- ifelse(seu$frequency < 2, 'Non-expanded', 'Expanded')

# Apply signatures
sigs <- read.csv('~/signatures/Tcellsignatures.csv', na.strings = c('', 'NA'))
for (c in 1:ncol(sigs)) {
  sig <- as.character(na.omit(sigs[, c]))
  seu <- AddModuleScore(object = seu, features = list(sig), 
                        name = colnames(sigs[c]), assay = 'RNA', search = F)
}

# Generate diffusion map
es <- as.ExpressionSet(as.data.frame(t(seu@assays$integrated@data)))
es@phenoData@data <- seu@meta.data
dm <- DiffusionMap(es, verbose = T, n_pcs = 20)


ifelse(!dir.exists(file.path('data/diffusion/MBM_sc/tcells')), 
       dir.create(file.path('data/diffusion/MBM_sc/tcells'), recursive = T), FALSE)
saveRDS(dm,'data/diffusion/MBM_sc/tcells/data_dm_MBM_sc_tcells_cd8.rds')

# Plots
pdf('data/diffusion/MBM_sc/tcells/plots_MBM_sc_tcells_dm_cd8.pdf')
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

palette(colExp)
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

set.seed(1)
dm_df<-as.data.frame(dm@eigenvectors)[sample(nrow(as.data.frame(dm@eigenvectors))),1:2]
ggplot(dm_df,aes(DC1,DC2,col=seu@meta.data[rownames(dm_df),'patient']))+
  geom_point(alpha = 1,shape=16,size=1.5) + theme_classic() +ggtitle('patient')+
  theme(legend.title = element_blank())
ggplot(dm_df,aes(DC1,DC2,col=seu@meta.data[rownames(dm_df),'sequencing']))+
  geom_point(alpha = 1,shape=16,size=1.5) + theme_classic() +ggtitle('sequencing')+
  theme(legend.title = element_blank())+scale_color_manual(values = colSCSN)
ggplot(dm_df,aes(DC1,DC2,col=seu@meta.data[rownames(dm_df),'organ']))+
  geom_point(alpha = 1,shape=16,size=1.5) + theme_classic() +ggtitle('organ')+
  theme(legend.title = element_blank())+scale_color_manual(values = colBP)
ggplot(dm_df,aes(DC1,DC2,col=seu@meta.data[rownames(dm_df),'cell_type_fine']))+
  geom_point(alpha = 1,shape=16,size=1.5) + theme_classic() +ggtitle('cell_type_fine')+
  theme(legend.title = element_blank())
ggplot(dm_df,aes(DC1,DC2,col=seu@meta.data[rownames(dm_df),'clone_size']))+
  geom_point(alpha = 1,shape=16,size=1.5) + theme_classic() +ggtitle('clone_size')+
  theme(legend.title = element_blank())+scale_color_manual(values = colExp)
ggplot(dm_df,aes(DC1,DC2,col=seu@meta.data[rownames(dm_df),'Exhaustion_Terminal_differentiation1']))+
  geom_point(alpha = 1,shape=16,size=1.5) + theme_classic() +ggtitle('Exhaustion_Terminal_differentiation1')+
  theme(legend.title = element_blank())+ scale_color_viridis()
dev.off()

