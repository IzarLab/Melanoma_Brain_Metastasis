#!/usr/bin/env Rscript

### title: Analysis of myeloid cells (UMAP, violin plots, diffusion maps, volcano plots) 
### author: Jana Biermann, PhD

library(Seurat)
library(destiny)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(gplots)
library(viridis)
library(scales)
library(plyr)
library(ggrastr)
library(ggrepel)
library(ggpubr)
library(patchwork)
library(hypeR)
'%notin%' <- Negate('%in%')

directory<-'data/MBPM/myeloid/'

HALLMARK <- msigdb_gsets(species='Homo sapiens', category='H')
KEGG     <- msigdb_gsets(species='Homo sapiens', category='C2', subcategory='CP:KEGG')
REACTOME <- msigdb_gsets(species='Homo sapiens', category='C2', subcategory='CP:REACTOME')
C5_BP <- msigdb_gsets(species='Homo sapiens', category='C5', subcategory='BP')

# Colors for UMAPs
colBP <- c('#A80D11', '#008DB8')
colSCSN <- c('#E1AC24', '#288F56')
colMye<-c('#316395','#ff9900','#109618','#dc3912','#0099c6','#990099',
          '#b82e2e','#66aa00','#3366cc')
colMye_MBM_sn<-c('#316395','#ff9900','#109618','#0099c6','#990099',
          '#b82e2e','#66aa00','#3366cc')
colMye_MPM_sn<-c('#316395','#ff9900','#109618','#dc3912','#0099c6','#990099',
                 '#66aa00','#3366cc')
colMye_MBM_sc<-c('#316395','#dd4477','#109618','#ff9900','#0099c6','#990099',
                '#b82e2e','#66aa00','#3366cc')
colMM<-c('#0099c6','#990099','#66aa00')


#### UMAPs MBPM_sn ####
filename<-'MBPM_sn_myeloid'

# Read-in integrated object, subset and rerun workflow
seu_all <- readRDS('data/MBPM/data_MBPM_scn.rds')
seu <- subset(seu_all, cell_type_main == 'Myeloid cells' & 
                sequencing == 'Single nuclei' &
                cell_type_int %notin% c('Low-quality cells', 'Doublets', 
                                        'Contamination', 'Undetermined',
                                        'Cycling cells'))
DefaultAssay(seu)<-'RNA'
seu<-NormalizeData(seu) %>% FindVariableFeatures()
seu <- ScaleData(seu) %>% RunPCA() 
ElbowPlot(seu,ndims = 50)
seu <- RunUMAP(seu,dims = 1:20)

# Adjust to updated cell types
seu$cell_type_myeloid<- ifelse(seu$cell_type_fine %in% 
                                 c('MDM M1-like', 'MDM M2-like','Proinflammatory MDM'),
                               'MDM',seu$cell_type_fine)

# Apply signatures
sigs <- read.csv('~/signatures/myeloid_signatures.csv', na.strings = '')
for (c in 1:ncol(sigs)) {
  sig <- as.character(na.omit(sigs[, c]))
  seu <- AddModuleScore(object = seu, features = list(sig), name = colnames(sigs[c]), 
                        assay = 'RNA', search = F)
}

saveRDS(seu,'data/MBPM/myeloid/data_MBPM_sn_myeloid_nonreint.rds')

# Plot UMAPs MBPM_sn
pdf(paste0(directory,'plots_',filename,'_umap.pdf'))
DimPlot(seu, group.by = 'cell_type_int', label = T, raster = T, shuffle = T)
DimPlot(seu, group.by = 'cell_type_int', raster = T, shuffle = T) + NoLegend()
DimPlot(seu, group.by = 'cell_type_myeloid', label = T, raster = T, shuffle = T,
        repel = T,cols = colMye)
DimPlot(seu, group.by = 'cell_type_myeloid', raster = T, shuffle = T,cols = colMye) + 
  NoLegend()
DimPlot(seu, group.by = 'cell_type_fine', label = T, raster = T, shuffle = T,repel = T) 
DimPlot(seu, group.by = 'cell_type_fine', raster = T, shuffle = T,cols = colMM) + NoLegend()
DimPlot(seu, group.by = 'sequencing', raster = T, shuffle = T, cols = colSCSN[2])
DimPlot(seu, group.by = 'sequencing', raster = T, shuffle = T, cols = colSCSN[2]) + NoLegend()
DimPlot(seu, group.by = 'organ', raster = T, shuffle = T, cols = colBP)
DimPlot(seu, group.by = 'organ', raster = T, shuffle = T, cols = colBP) + NoLegend()
DimPlot(seu, group.by = 'patient', raster = T, shuffle = T)
DimPlot(seu, group.by = 'patient', raster = T, shuffle = T) + NoLegend()
FeaturePlot(seu, features = c('M1_ada1','M2_ada1'),raster = T,order = T,
            min.cutoff = 'q05',max.cutoff = 'q95')
FeatureScatter(seu,feature1 = 'M1_ada1',feature2 = 'M2_ada1',shuffle = T,
               raster = T,group.by = 'cell_type_myeloid')+ s
cale_color_manual(values = colMye)
dev.off()


#### UMAPs MBM_sc ####
filename<-'MBM_sc_myeloid'
seu<-readRDS('data/cell_type_DEG/MBM_sc/myeloid/data_MBM_sc_myeloid.rds')
seu <- subset(seu, cell_type_fine %notin% c('Low-quality cells', 'Doublets', 
                                            'Contamination', 'Undetermined',
                                            'Cycling cells'))
seu <- ScaleData(seu) %>% RunPCA() 
ElbowPlot(seu,ndims = 50)
seu <- RunUMAP(seu,dims = 1:25)

# Adjust to updated cell types
seu$cell_type_myeloid<- ifelse(seu$cell_type_fine %in% 
                                 c('MDM M1-like', 'MDM M2-like','Proinflammatory MDM',
                                   'Proinflammatory MDM (glycolysis)',
                                   'Proinflammatory MDM (zinc)'),
                               'MDM',seu$cell_type_fine)

# Apply signatures
sigs <- read.csv('~/signatures/myeloid_signatures.csv', na.strings = '')
for (c in 1:ncol(sigs)) {
  sig <- as.character(na.omit(sigs[, c]))
  seu <- AddModuleScore(object = seu, features = list(sig), 
                        name = colnames(sigs[c]), assay = 'RNA', search = F)
}

# Plot UMAPs MBM_sc
pdf(paste0(directory,'plots_',filename,'_umap.pdf'))
DimPlot(seu, group.by = 'cell_type_myeloid', label = T, raster = T, shuffle = T,
        repel = T)+ scale_color_manual(values = colMye_MBM_sc)
DimPlot(seu, group.by = 'cell_type_myeloid', raster = T, shuffle = T) + 
  NoLegend()+ scale_color_manual(values = colMye_MBM_sc)
DimPlot(seu, group.by = 'cell_type_fine', label = T, raster = T, shuffle = T,repel = T) 
DimPlot(seu, group.by = 'cell_type_fine', raster = T, shuffle = T) + NoLegend()
DimPlot(seu, group.by = 'sequencing', raster = T, shuffle = T, cols = colSCSN[2])
DimPlot(seu, group.by = 'sequencing', raster = T, shuffle = T, cols = colSCSN[2]) + NoLegend()
DimPlot(seu, group.by = 'organ', raster = T, shuffle = T, cols = colBP)
DimPlot(seu, group.by = 'organ', raster = T, shuffle = T, cols = colBP) + NoLegend()
DimPlot(seu, group.by = 'patient', raster = T, shuffle = T)
DimPlot(seu, group.by = 'patient', raster = T, shuffle = T) + NoLegend()
FeaturePlot(seu, features = c('M1_ada1','M2_ada1'),raster = T,order = T,
            min.cutoff = 'q05',max.cutoff = 'q95')
FeatureScatter(seu,feature1 = 'M1_ada1',feature2 = 'M2_ada1',shuffle = T,
               raster = T,group.by = 'cell_type_myeloid')+ 
  scale_color_manual(values = colMye_MBM_sc)
dev.off()


#### UMAPs MBM_sn ####
filename<-'MBM_sn_myeloid'
seu<-readRDS('data/cell_type_DEG/MBM_sn/myeloid/data_MBM_sn_myeloid.rds')
seu <- subset(seu, cell_type_fine %notin% c('Low-quality cells', 'Doublets', 
                                            'Contamination', 'Undetermined',
                                            'Cycling cells'))
seu <- ScaleData(seu) %>% RunPCA() 
ElbowPlot(seu,ndims = 50)
seu <- RunUMAP(seu,dims = 1:20)

# Adjust to updated cell types
seu$cell_type_myeloid<- ifelse(seu$cell_type_fine %in% 
                                 c('MDM M1-like', 'MDM M2-like','Proinflammatory MDM'),
                               'MDM',seu$cell_type_fine)

# Apply signatures
sigs <- read.csv('~/signatures/myeloid_signatures.csv', na.strings = '')
for (c in 1:ncol(sigs)) {
  sig <- as.character(na.omit(sigs[, c]))
  seu <- AddModuleScore(object = seu, features = list(sig), 
                        name = colnames(sigs[c]), assay = 'RNA', search = F)
}

# Plot UMAPs MBM_sn
pdf(paste0(directory,'plots_',filename,'_umap.pdf'))
DimPlot(seu, group.by = 'cell_type_myeloid', label = T, raster = T, shuffle = T,
        repel = T)+ scale_color_manual(values = colMye_MBM_sn)
DimPlot(seu, group.by = 'cell_type_myeloid', raster = T, shuffle = T) + 
  NoLegend()+ scale_color_manual(values = colMye_MBM_sn)
DimPlot(seu, group.by = 'cell_type_fine', label = T, raster = T, shuffle = T,repel = T) 
DimPlot(seu, group.by = 'cell_type_fine', raster = T, shuffle = T) + NoLegend()
DimPlot(seu, group.by = 'sequencing', raster = T, shuffle = T, cols = colSCSN[2])
DimPlot(seu, group.by = 'sequencing', raster = T, shuffle = T, cols = colSCSN[2]) + NoLegend()
DimPlot(seu, group.by = 'organ', raster = T, shuffle = T, cols = colBP)
DimPlot(seu, group.by = 'organ', raster = T, shuffle = T, cols = colBP) + NoLegend()
DimPlot(seu, group.by = 'patient', raster = T, shuffle = T)
DimPlot(seu, group.by = 'patient', raster = T, shuffle = T) + NoLegend()
FeaturePlot(seu, features = c('M1_ada1','M2_ada1'),raster = T,order = T,
            min.cutoff = 'q05',max.cutoff = 'q95')
FeatureScatter(seu,feature1 = 'M1_ada1',feature2 = 'M2_ada1',shuffle = T,
               raster = T,group.by = 'cell_type_myeloid')+ 
  scale_color_manual(values = colMye_MBM_sn)
dev.off()

#### UMAPs MPM_sn ####
filename<-'MPM_sn_myeloid'

seu<-readRDS('data/cell_type_DEG/MPM_sn/myeloid/data_MPM_sn_myeloid.rds')
seu <- subset(seu, cell_type_fine %notin% c('Low-quality cells', 'Doublets', 
                                            'Contamination', 'Undetermined',
                                            'Cycling cells'))
seu <- ScaleData(seu) %>% RunPCA() 
ElbowPlot(seu,ndims = 50)
seu <- RunUMAP(seu,dims = 1:20)

# Adjust to updated cell types
seu$cell_type_myeloid<- ifelse(seu$cell_type_fine %in% 
                                 c('MDM M1-like', 'MDM M2-like','Proinflammatory MDM'),
                               'MDM',seu$cell_type_fine)

# Apply signatures
sigs <- read.csv('~/signatures/myeloid_signatures.csv', na.strings = '')
for (c in 1:ncol(sigs)) {
  sig <- as.character(na.omit(sigs[, c]))
  seu <- AddModuleScore(object = seu, features = list(sig), 
                        name = colnames(sigs[c]), assay = 'RNA', search = F)
}

# Plot UMAPs MPM_sn
pdf(paste0(directory,'plots_',filename,'_umap.pdf'))
DimPlot(seu, group.by = 'cell_type_myeloid', label = T, raster = T, shuffle = T,
        repel = T)+ scale_color_manual(values = colMye_MPM_sn)
DimPlot(seu, group.by = 'cell_type_myeloid', raster = T, shuffle = T) + 
  NoLegend()+ scale_color_manual(values = colMye_MPM_sn)
DimPlot(seu, group.by = 'cell_type_fine', label = T, raster = T, shuffle = T,repel = T) 
DimPlot(seu, group.by = 'cell_type_fine', raster = T, shuffle = T) + NoLegend()
DimPlot(seu, group.by = 'sequencing', raster = T, shuffle = T, cols = colSCSN[2])
DimPlot(seu, group.by = 'sequencing', raster = T, shuffle = T, cols = colSCSN[2]) + NoLegend()
DimPlot(seu, group.by = 'organ', raster = T, shuffle = T, cols = colBP[2])
DimPlot(seu, group.by = 'organ', raster = T, shuffle = T, cols = colBP[2]) + NoLegend()
DimPlot(seu, group.by = 'patient', raster = T, shuffle = T)
DimPlot(seu, group.by = 'patient', raster = T, shuffle = T) + NoLegend()
FeaturePlot(seu, features = c('M1_ada1','M2_ada1'),raster = T,order = T,
            min.cutoff = 'q05',max.cutoff = 'q95')
FeatureScatter(seu,feature1 = 'M1_ada1',feature2 = 'M2_ada1',shuffle = T,
               raster = T,group.by = 'cell_type_myeloid')+ 
  scale_color_manual(values = colMye_MPM_sn)
dev.off()


#### Violin plots ####
filename<-'MBPM_sn_myeloid'

seu_all <- readRDS('data/MBPM/data_MBPM_scn.rds')
seu <- subset(seu_all, cell_type_main == 'Myeloid cells' & 
                sequencing == 'Single nuclei' &
                cell_type_int %notin% c('Low-quality cells', 'Doublets', 'Contamination',
                                        'Undetermined','Cycling cells'))

mdm <- subset(seu, cell_type_fine %in% c('MDM','MDM FTL+'))

genes_main <- c('CLEC9A', 'XCR1','CLNK', # cDC1
                'SLC38A1', 'RTN1','SIRPA', #cDC2
                #'MKI67', 'TOP2A', # cycling
                'IDO1', 'LAMP3','CD200', #DC3
                'MS4A2','KIT', #mast
                'MRC1','MERTK', #MDM
                'FCGR1A','CD38', #M1
                'CD163L1','SELENOP','F13A1','DAB2','SIGLEC1', #M2
                'FTL','FTH1', #FTL
                'NAV3','C3', 'P2RY12', # microglia
                'VCAN', 'FCN1','LYZ'  # Monocytes
                )

genes_dc <- c('CLEC9A', 'XCR1', 'CLNK', 'ENPP1', 'BCL11A','BATF3', #cDC1
              'SLC38A1', 'RTN1',  'CD1D', 'FCER1A', 'CLEC10A', #cDC2
              'LAMP3', 'IDO1', 'IDO2', 'CD80', 'CD274', 
              'PDCD1LG2', 'CD200', 'CCR7', 'FSCN1', 'CD40' #DC3
              )


genes_BvsP<-c('CD163','MRC1','TLR2','MERTK', 'AXL','IL2RA','HLA-A',
              'HLA-B','HLA-C','B2M','CIITA','CD74','HLA-DRA','HLA-DPA1',
              'HLA-DQB1','HLA-DQA1','HLA-DPB1','MM2','MMP4')

pdf(paste0(directory,'plots_',filename,'_markers_BvsP.pdf'), height = 9, width = 4)
VlnPlot(mdm, features = genes_BvsP, group.by = 'cell_type_myeloid', pt.size = 0, 
        assay = 'RNA', stack = T, flip = T,split.by = 'organ',cols = colBP) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = 'bottom',
        axis.title.x = element_blank())
dev.off()

pdf(paste0(directory,'plots_',filename,'_markers_main.pdf'), height = 9, width = 4)
VlnPlot(seu, features = genes_main, group.by = 'cell_type_myeloid', pt.size = 0, 
        assay = 'RNA', stack = T, flip = T,split.by = 'organ',cols = colBP) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = 'bottom',
        axis.title.x = element_blank())

DotPlot(seu,assay = 'RNA',features = genes_main,
        group.by = 'cell_type_myeloid', dot.scale = 5)+ 
  scale_color_viridis() + coord_flip()+RotatedAxis()+
  theme(plot.margin = margin(0, 0, 8, 0, 'cm'),
        text = element_text(size=5),
        axis.text= element_text(size=5))

VlnPlot(subset(seu,cell_type_int=='Dendritic cells'), features = genes_dc, 
        group.by = 'cell_type_myeloid', pt.size = 0, 
        assay = 'RNA', stack = T, flip = T,split.by = 'organ',cols = colBP) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = 'bottom',
        axis.title.x = element_blank())
dev.off()


#### Diffusion component analysis of monocytes + MDMs for MBPM_sn ####
# Subset to monocytes + MDMs
seu<-readRDS('data/MBPM/myeloid/data_MBPM_sn_myeloid_nonreint.rds')

seu$cell_type_myeloid<- ifelse(seu$cell_type_fine %in% 
                                 c('MDM M1-like', 'MDM M2-like','Proinflammatory MDM'),
                               'MDM',seu$cell_type_fine)

seu <- subset(seu, cell_type_myeloid %in% c('Monocytes', 'MDM', 'MDM FTL+'))
seu <- ScaleData(seu)


# Apply signatures
sigs <- read.csv('~/signatures/myeloid_signatures.csv', na.strings = '')
for (c in 1:ncol(sigs)) {
  sig <- as.character(na.omit(sigs[, c]))
  seu <- AddModuleScore(object = seu, features = list(sig), name = colnames(sigs[c]), 
                        assay = 'RNA', search = F)
}

# Generate diffusion map
es <- as.ExpressionSet(as.data.frame(t(seu@assays$integrated@data)))
es@phenoData@data <- seu@meta.data
dm <- DiffusionMap(es, verbose = T, n_pcs = 20)

# Diffusion pseudotime calculation
dpt <- DPT(dm)
es$pseudotime_dpt <- rank(dpt$dpt)
seu$pseudotime_dpt <- rank(dpt$dpt)
seu <- AddMetaData(seu, dm@eigenvectors, col.name = colnames(dm@eigenvectors))

# Plots
pdf(paste0(directory,'plots_',filename,'_MM_dm.pdf'))
par(mar = c(5.1, 4.1, 4.1, 2.1), xpd = TRUE)
palette(colMM)
tmp_var<-seu$cell_type_myeloid
plot(dm, col = as.factor(tmp_var), main = 'cell_type_myeloid', pch = 20)
legend('bottomright', inset = c(0.1, -0.2), legend = levels(as.factor(tmp_var)), 
       pch = 16, col = as.factor(levels(as.factor(tmp_var))), bty = 'n')

palette(colBP)
tmp_var<-seu$organ
plot(dm, col = as.factor(tmp_var), main = 'organ', pch = 20)
legend('bottomright', inset = c(0.1, -0.2), legend = levels(as.factor(tmp_var)), 
       pch = 16, col = as.factor(levels(as.factor(tmp_var))), bty = 'n')

palette(hue_pal()(length(unique(seu$cell_type_fine))))
tmp_var<-seu$cell_type_fine
plot(dm, col = as.factor(tmp_var), main = 'cell_type_fine', pch = 20)
legend('bottomright', inset = c(-0.05, -0.2), legend = levels(as.factor(tmp_var)), 
       pch = 16, col = as.factor(levels(as.factor(tmp_var))), bty = 'n')

par(mar = c(5.1, 4.1, 4.1, 7), xpd = TRUE)
palette(viridis(100))
sig_col<-seu$M1_ada1
plot(dm, col = sig_col, main = 'M1_ada1', pch = 20)
colorlegend(sig_col, viridis(length(sig_col)), 
            posx = c(0.88, 0.9), posy = c(0, 0.65))

palette(viridis(100))
sig_col<-seu$M2_ada1
plot(dm, col = sig_col, main = 'M2_ada1', pch = 20)
colorlegend(sig_col, viridis(length(sig_col)), 
            posx = c(0.88, 0.9), posy = c(0, 0.65))

palette(viridis(100))
sig_col<-seu$pseudotime_dpt
plot(dm, col = sig_col, main = 'pseudotime_dpt', pch = 20)
colorlegend(sig_col, viridis(length(sig_col)), 
            posx = c(0.88, 0.9), posy = c(0, 0.65))

plot(dm, col = seu@assays$RNA@data['FTL',], main = 'FTL', pch = 20)
colorlegend(seu@assays$RNA@data['FTL',], viridis(length(seu@assays$RNA@data['FTL',])), 
            posx = c(0.88, 0.9), posy = c(0, 0.65))

plot(dm, col = seu@assays$RNA@data['SPP1',], main = 'SPP1', pch = 20)
colorlegend(seu@assays$RNA@data['SPP1',], viridis(length(seu@assays$RNA@data['SPP1',])), 
            posx = c(0.88, 0.9), posy = c(0, 0.65))

set.seed(1)
dm_df<-as.data.frame(dm@eigenvectors)[sample(nrow(as.data.frame(dm@eigenvectors))),1:2]
ggplot(dm_df,aes(DC1,DC2,col=seu@meta.data[rownames(dm_df),'cell_type_myeloid']))+
  geom_point(alpha = 1,shape=16,size=1.5) + theme_classic() +ggtitle('cell_type_myeloid')+
  theme(legend.title = element_blank())+scale_color_manual(values = colMM)

ggplot(dm_df,aes(DC1,DC2,col=seu@meta.data[rownames(dm_df),'cell_type_fine']))+
  geom_point(alpha = 1,shape=16,size=1.5) + theme_classic() +ggtitle('cell_type_fine')+
  theme(legend.title = element_blank())

ggplot(dm_df,aes(DC1,DC2,col=seu@meta.data[rownames(dm_df),'organ']))+
  geom_point(alpha = 1,shape=16,size=1.5) + theme_classic() +ggtitle('organ')+
  theme(legend.title = element_blank())+scale_color_manual(values = colBP)

ggplot(dm_df,aes(DC1,DC2,col=seu@meta.data[rownames(dm_df),'M1_ada1']))+
  geom_point(alpha = 1,shape=16,size=1.5) + theme_classic() +ggtitle('M1_ada1')+
  theme(legend.title = element_blank())+scale_color_viridis()

ggplot(dm_df,aes(DC1,DC2,col=seu@meta.data[rownames(dm_df),'M2_ada1']))+
  geom_point(alpha = 1,shape=16,size=1.5) + theme_classic() +ggtitle('M2_ada1')+
  theme(legend.title = element_blank())+scale_color_viridis()

ggplot(dm_df,aes(DC1,DC2,col=seu@assays$RNA@data['FTL',rownames(dm_df)]))+
  geom_point(alpha = 1,shape=16,size=1.5) + theme_classic() +ggtitle('FTL')+
  theme(legend.title = element_blank())+scale_color_viridis()

plot(dpt)
dev.off()

saveRDS(dm,paste0(directory,'plots_',filename,'_MM_dm.rds'))

#### Diffusion component analysis of monocytes + MDMs for MBM_sc ####
filename<-'MBM_sc_myeloid'

# Subset to monocytes + MDMs
seu<-readRDS('data/cell_type_DEG/MBM_sc/myeloid/data_MBM_sc_myeloid.rds')

seu$cell_type_myeloid<- ifelse(seu$cell_type_fine %in% 
                                 c('MDM M1-like', 'MDM M2-like','Proinflammatory MDM',
                                   'Proinflammatory MDM (glycolysis)',
                                   'Proinflammatory MDM (zinc)'),
                               'MDM',seu$cell_type_fine)

seu <- subset(seu, cell_type_myeloid %in% c('MDM', 'MDM FTL+', 'Monocytes'))
seu <- ScaleData(seu)

# Apply signatures
sigs <- read.csv('~/signatures/myeloid_signatures.csv', na.strings = '')
for (c in 1:ncol(sigs)) {
  sig <- as.character(na.omit(sigs[, c]))
  seu <- AddModuleScore(object = seu, features = list(sig), name = colnames(sigs[c]), 
                        assay = 'RNA', search = F)
}

# Generate diffusion map
es <- as.ExpressionSet(as.data.frame(t(seu@assays$integrated@data)))
es@phenoData@data <- seu@meta.data
dm <- DiffusionMap(es, verbose = T, n_pcs = 20)

# Plots
pdf(paste0(directory,'plots_',filename,'_MM_dm.pdf'))
par(mar = c(5.1, 4.1, 4.1, 2.1), xpd = TRUE)
palette(colMM)
tmp_var<-seu$cell_type_myeloid
plot(dm, col = as.factor(tmp_var), main = 'cell_type_myeloid', pch = 20)
legend('bottomright', inset = c(0.1, -0.2), legend = levels(as.factor(tmp_var)), 
       pch = 16, col = as.factor(levels(as.factor(tmp_var))), bty = 'n')

palette(colBP)
tmp_var<-seu$organ
plot(dm, col = as.factor(tmp_var), main = 'organ', pch = 20)
legend('bottomright', inset = c(0.1, -0.2), legend = levels(as.factor(tmp_var)), 
       pch = 16, col = as.factor(levels(as.factor(tmp_var))), bty = 'n')

palette(hue_pal()(length(unique(seu$cell_type_fine))))
tmp_var<-seu$cell_type_fine
plot(dm, col = as.factor(tmp_var), main = 'cell_type_fine', pch = 20)
legend('bottomright', inset = c(-0.05, -0.2), legend = levels(as.factor(tmp_var)), 
       pch = 16, col = as.factor(levels(as.factor(tmp_var))), bty = 'n')

par(mar = c(5.1, 4.1, 4.1, 7), xpd = TRUE)
palette(viridis(100))
sig_col<-seu$M1_ada1
plot(dm, col = sig_col, main = 'M1_ada1', pch = 20)
colorlegend(sig_col, viridis(length(sig_col)), 
            posx = c(0.88, 0.9), posy = c(0, 0.65))

palette(viridis(100))
sig_col<-seu$M2_ada1
plot(dm, col = sig_col, main = 'M2_ada1', pch = 20)
colorlegend(sig_col, viridis(length(sig_col)), 
            posx = c(0.88, 0.9), posy = c(0, 0.65))

plot(dm, col = seu@assays$RNA@data['FTL',], main = 'FTL', pch = 20)
colorlegend(seu@assays$RNA@data['FTL',], viridis(length(seu@assays$RNA@data['FTL',])), 
            posx = c(0.88, 0.9), posy = c(0, 0.65))

plot(dm, col = seu@assays$RNA@data['SPP1',], main = 'SPP1', pch = 20)
colorlegend(seu@assays$RNA@data['SPP1',], viridis(length(seu@assays$RNA@data['SPP1',])), 
            posx = c(0.88, 0.9), posy = c(0, 0.65))

set.seed(1)
dm_df<-as.data.frame(dm@eigenvectors)[sample(nrow(as.data.frame(dm@eigenvectors))),1:2]
ggplot(dm_df,aes(DC1,DC2,col=seu@meta.data[rownames(dm_df),'cell_type_myeloid']))+
  geom_point(alpha = 1,shape=16,size=1.5) + theme_classic() +ggtitle('cell_type_myeloid')+
  theme(legend.title = element_blank())+scale_color_manual(values = colMM)

ggplot(dm_df,aes(DC1,DC2,col=seu@meta.data[rownames(dm_df),'cell_type_fine']))+
  geom_point(alpha = 1,shape=16,size=1.5) + theme_classic() +ggtitle('cell_type_fine')+
  theme(legend.title = element_blank())

ggplot(dm_df,aes(DC1,DC2,col=seu@meta.data[rownames(dm_df),'organ']))+
  geom_point(alpha = 1,shape=16,size=1.5) + theme_classic() +ggtitle('organ')+
  theme(legend.title = element_blank())+scale_color_manual(values = colBP)

ggplot(dm_df,aes(DC1,DC2,col=seu@meta.data[rownames(dm_df),'M1_ada1']))+
  geom_point(alpha = 1,shape=16,size=1.5) + theme_classic() +ggtitle('M1_ada1')+
  theme(legend.title = element_blank())+scale_color_viridis()

ggplot(dm_df,aes(DC1,DC2,col=seu@meta.data[rownames(dm_df),'M2_ada1']))+
  geom_point(alpha = 1,shape=16,size=1.5) + theme_classic() +ggtitle('M2_ada1')+
  theme(legend.title = element_blank())+scale_color_viridis()

ggplot(dm_df,aes(DC1,DC2,col=seu@assays$RNA@data['FTL',rownames(dm_df)]))+
  geom_point(alpha = 1,shape=16,size=1.5) + theme_classic() +ggtitle('FTL')+
  theme(legend.title = element_blank())+scale_color_viridis()
dev.off()

saveRDS(dm,paste0(directory,'plots_',filename,'_MM_dm.rds'))


#### Diffusion component analysis of monocytes + MDMs for MBM_sn ####
filename<-'MBM_sn_myeloid'

# Subset to monocytes + MDMs
seu<-readRDS('data/cell_type_DEG/MBM_sn/myeloid/data_MBM_sn_myeloid.rds')
seu <- subset(seu, cell_type_myeloid %in% c('MDM', 'MDM FTL+', 'Monocytes'))
seu <- ScaleData(seu)

# Apply signatures
sigs <- read.csv('~/signatures/myeloid_signatures.csv', na.strings = '')
for (c in 1:ncol(sigs)) {
  sig <- as.character(na.omit(sigs[, c]))
  seu <- AddModuleScore(object = seu, features = list(sig), name = colnames(sigs[c]), 
                        assay = 'RNA', search = F)
}

# Generate diffusion map
es <- as.ExpressionSet(as.data.frame(t(seu@assays$integrated@data)))
es@phenoData@data <- seu@meta.data
dm <- DiffusionMap(es, verbose = T, n_pcs = 20)

# Plots
pdf(paste0(directory,'plots_',filename,'_MM_dm.pdf'))
par(mar = c(5.1, 4.1, 4.1, 2.1), xpd = TRUE)
palette(colMM)
tmp_var<-seu$cell_type_myeloid
plot(dm, col = as.factor(tmp_var), main = 'cell_type_myeloid', pch = 20)
legend('bottomright', inset = c(0.1, -0.2), legend = levels(as.factor(tmp_var)), 
       pch = 16, col = as.factor(levels(as.factor(tmp_var))), bty = 'n')

palette(colBP)
tmp_var<-seu$organ
plot(dm, col = as.factor(tmp_var), main = 'organ', pch = 20)
legend('bottomright', inset = c(0.1, -0.2), legend = levels(as.factor(tmp_var)), 
       pch = 16, col = as.factor(levels(as.factor(tmp_var))), bty = 'n')

palette(hue_pal()(length(unique(seu$cell_type_fine))))
tmp_var<-seu$cell_type_fine
plot(dm, col = as.factor(tmp_var), main = 'cell_type_fine', pch = 20)
legend('bottomright', inset = c(-0.05, -0.2), legend = levels(as.factor(tmp_var)), 
       pch = 16, col = as.factor(levels(as.factor(tmp_var))), bty = 'n')

par(mar = c(5.1, 4.1, 4.1, 7), xpd = TRUE)
palette(viridis(100))
sig_col<-seu$M1_ada1
plot(dm, col = sig_col, main = 'M1_ada1', pch = 20)
colorlegend(sig_col, viridis(length(sig_col)), 
            posx = c(0.88, 0.9), posy = c(0, 0.65))

palette(viridis(100))
sig_col<-seu$M2_ada1
plot(dm, col = sig_col, main = 'M2_ada1', pch = 20)
colorlegend(sig_col, viridis(length(sig_col)), 
            posx = c(0.88, 0.9), posy = c(0, 0.65))

plot(dm, col = seu@assays$RNA@data['FTL',], main = 'FTL', pch = 20)
colorlegend(seu@assays$RNA@data['FTL',], viridis(length(seu@assays$RNA@data['FTL',])), 
            posx = c(0.88, 0.9), posy = c(0, 0.65))

plot(dm, col = seu@assays$RNA@data['SPP1',], main = 'SPP1', pch = 20)
colorlegend(seu@assays$RNA@data['SPP1',], viridis(length(seu@assays$RNA@data['SPP1',])), 
            posx = c(0.88, 0.9), posy = c(0, 0.65))

set.seed(1)
dm_df<-as.data.frame(dm@eigenvectors)[sample(nrow(as.data.frame(dm@eigenvectors))),1:2]
ggplot(dm_df,aes(DC1,DC2,col=seu@meta.data[rownames(dm_df),'cell_type_myeloid']))+
  geom_point(alpha = 1,shape=16,size=1.5) + theme_classic() +ggtitle('cell_type_myeloid')+
  theme(legend.title = element_blank())+scale_color_manual(values = colMM)

ggplot(dm_df,aes(DC1,DC2,col=seu@meta.data[rownames(dm_df),'cell_type_fine']))+
  geom_point(alpha = 1,shape=16,size=1.5) + theme_classic() +ggtitle('cell_type_fine')+
  theme(legend.title = element_blank())

ggplot(dm_df,aes(DC1,DC2,col=seu@meta.data[rownames(dm_df),'organ']))+
  geom_point(alpha = 1,shape=16,size=1.5) + theme_classic() +ggtitle('organ')+
  theme(legend.title = element_blank())+scale_color_manual(values = colBP)

ggplot(dm_df,aes(DC1,DC2,col=seu@meta.data[rownames(dm_df),'M1_ada1']))+
  geom_point(alpha = 1,shape=16,size=1.5) + theme_classic() +ggtitle('M1_ada1')+
  theme(legend.title = element_blank())+scale_color_viridis()

ggplot(dm_df,aes(DC1,DC2,col=seu@meta.data[rownames(dm_df),'M2_ada1']))+
  geom_point(alpha = 1,shape=16,size=1.5) + theme_classic() +ggtitle('M2_ada1')+
  theme(legend.title = element_blank())+scale_color_viridis()

ggplot(dm_df,aes(DC1,DC2,col=seu@assays$RNA@data['FTL',rownames(dm_df)]))+
  geom_point(alpha = 1,shape=16,size=1.5) + theme_classic() +ggtitle('FTL')+
  theme(legend.title = element_blank())+scale_color_viridis()
dev.off()

saveRDS(dm,paste0(directory,'plots_',filename,'_MM_dm.rds'))


#### Diffusion component analysis of monocytes + MDMs for MPM_sn ####
filename<-'MPM_sn_myeloid'

# Subset to monocytes + MDMs
seu<-readRDS('data/cell_type_DEG/MPM_sn/myeloid/data_MPM_sn_myeloid.rds')

seu$cell_type_myeloid<- ifelse(seu$cell_type_fine %in% 
                                 c('MDM M1-like', 'MDM M2-like','Proinflammatory MDM'),
                               'MDM',seu$cell_type_fine)

seu <- subset(seu, cell_type_myeloid %in% c('MDM', 'MDM FTL+', 'Monocytes'))
seu <- ScaleData(seu)

# Apply signatures
sigs <- read.csv('~/signatures/myeloid_signatures.csv', na.strings = '')
for (c in 1:ncol(sigs)) {
  sig <- as.character(na.omit(sigs[, c]))
  seu <- AddModuleScore(object = seu, features = list(sig), 
                        name = colnames(sigs[c]), assay = 'RNA', search = F)
}

# Generate diffusion map
es <- as.ExpressionSet(as.data.frame(t(seu@assays$integrated@data)))
es@phenoData@data <- seu@meta.data
dm <- DiffusionMap(es, verbose = T, n_pcs = 20)

# Plots
pdf(paste0(directory,'plots_',filename,'_MM_dm.pdf'))
par(mar = c(5.1, 4.1, 4.1, 2.1), xpd = TRUE)
palette(colMM)
tmp_var<-seu$cell_type_myeloid
plot(dm, col = as.factor(tmp_var), main = 'cell_type_myeloid', pch = 20)
legend('bottomright', inset = c(0.1, -0.2), legend = levels(as.factor(tmp_var)), 
       pch = 16, col = as.factor(levels(as.factor(tmp_var))), bty = 'n')

palette(colBP)
tmp_var<-seu$organ
plot(dm, col = as.factor(tmp_var), main = 'organ', pch = 20)
legend('bottomright', inset = c(0.1, -0.2), legend = levels(as.factor(tmp_var)), 
       pch = 16, col = as.factor(levels(as.factor(tmp_var))), bty = 'n')

palette(hue_pal()(length(unique(seu$cell_type_fine))))
tmp_var<-seu$cell_type_fine
plot(dm, col = as.factor(tmp_var), main = 'cell_type_fine', pch = 20)
legend('bottomright', inset = c(-0.05, -0.2), legend = levels(as.factor(tmp_var)), 
       pch = 16, col = as.factor(levels(as.factor(tmp_var))), bty = 'n')

par(mar = c(5.1, 4.1, 4.1, 7), xpd = TRUE)
palette(viridis(100))
sig_col<-seu$M1_ada1
plot(dm, col = sig_col, main = 'M1_ada1', pch = 20)
colorlegend(sig_col, viridis(length(sig_col)), 
            posx = c(0.88, 0.9), posy = c(0, 0.65))

palette(viridis(100))
sig_col<-seu$M2_ada1
plot(dm, col = sig_col, main = 'M2_ada1', pch = 20)
colorlegend(sig_col, viridis(length(sig_col)), 
            posx = c(0.88, 0.9), posy = c(0, 0.65))

plot(dm, col = seu@assays$RNA@data['FTL',], main = 'FTL', pch = 20)
colorlegend(seu@assays$RNA@data['FTL',], viridis(length(seu@assays$RNA@data['FTL',])), 
            posx = c(0.88, 0.9), posy = c(0, 0.65))

plot(dm, col = seu@assays$RNA@data['SPP1',], main = 'SPP1', pch = 20)
colorlegend(seu@assays$RNA@data['SPP1',], viridis(length(seu@assays$RNA@data['SPP1',])), 
            posx = c(0.88, 0.9), posy = c(0, 0.65))

set.seed(1)
dm_df<-as.data.frame(dm@eigenvectors)[sample(nrow(as.data.frame(dm@eigenvectors))),1:2]
ggplot(dm_df,aes(DC1,DC2,col=seu@meta.data[rownames(dm_df),'cell_type_myeloid']))+
  geom_point(alpha = 1,shape=16,size=1.5) + theme_classic() +ggtitle('cell_type_myeloid')+
  theme(legend.title = element_blank())+scale_color_manual(values = colMM)

ggplot(dm_df,aes(DC1,DC2,col=seu@meta.data[rownames(dm_df),'cell_type_fine']))+
  geom_point(alpha = 1,shape=16,size=1.5) + theme_classic() +ggtitle('cell_type_fine')+
  theme(legend.title = element_blank())

ggplot(dm_df,aes(DC1,DC2,col=seu@meta.data[rownames(dm_df),'organ']))+
  geom_point(alpha = 1,shape=16,size=1.5) + theme_classic() +ggtitle('organ')+
  theme(legend.title = element_blank())+scale_color_manual(values = colBP)

ggplot(dm_df,aes(DC1,DC2,col=seu@meta.data[rownames(dm_df),'M1_ada1']))+
  geom_point(alpha = 1,shape=16,size=1.5) + theme_classic() +ggtitle('M1_ada1')+
  theme(legend.title = element_blank())+scale_color_viridis()

ggplot(dm_df,aes(DC1,DC2,col=seu@meta.data[rownames(dm_df),'M2_ada1']))+
  geom_point(alpha = 1,shape=16,size=1.5) + theme_classic() +ggtitle('M2_ada1')+
  theme(legend.title = element_blank())+scale_color_viridis()

ggplot(dm_df,aes(DC1,DC2,col=seu@assays$RNA@data['FTL',rownames(dm_df)]))+
  geom_point(alpha = 1,shape=16,size=1.5) + theme_classic() +ggtitle('FTL')+
  theme(legend.title = element_blank())+scale_color_viridis()
dev.off()

saveRDS(dm,paste0(directory,'plots_',filename,'_MM_dm.rds'))


#### DEG and pathway analysis of MDMs MBM_sc ####
filename<-'MBM_sc_myeloid'

# Subset to MDMs
seu<-readRDS('data/cell_type_DEG/MBM_sc/myeloid/data_MBM_sc_myeloid.rds')
seu$cell_type_myeloid<- ifelse(seu$cell_type_fine %in% 
                                 c('MDM M1-like', 'MDM M2-like','Proinflammatory MDM',
                                   'Proinflammatory MDM (glycolysis)',
                                   'Proinflammatory MDM (zinc)'),
                               'MDM',seu$cell_type_fine)
seu <- subset(seu, cell_type_myeloid %in% c('MDM', 'MDM FTL+'))
seu <- ScaleData(seu)

# Differential gene expression (DGE) based on cell types
Idents(seu)<-seu$cell_type_myeloid
markers <- FindAllMarkers(seu, only.pos = TRUE,assay='RNA',test.use = 'MAST',
                          min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, paste0(directory,'markers_',filename,'.csv'), row.names = F)
markers<- markers %>% filter(p_val_adj < 0.05)

background<-rownames(seu@assays$RNA@data)
pdf(paste0(directory,'markers_',filename,'_pathway.pdf'),width = 10,height = 7)
for(pw in c(C5_BP,HALLMARK,KEGG,REACTOME)){
  for(i in unique(markers$cluster)){
    subs<-markers %>% dplyr::filter(cluster==i) %>% magrittr::use_series(gene)
    print(hyp_dots(hypeR(as.character(subs), genesets=pw, background=background),
                   title = paste0(pw$name,' cluster ',i),abrv=80)+ theme_bw() + 
            theme(axis.title.y = element_blank(), 
                  axis.text.y = element_text(color = 'black'), 
                  axis.text.x = element_text(color = 'black'), 
                  panel.grid.minor = element_blank()) + 
            scale_color_gradient(high = 'black', low = '#C51B8A'))
  }
}
dev.off()


# Markers for volcano plot
Idents(seu)<-seu$cell_type_myeloid
markers <- FindMarkers(seu, ident.1 = 'MDM', ident.2 ='MDM FTL+', min.pct = 0, 
                       logfc.threshold = 0,test.use = 'MAST',
                       assay = 'RNA')
markers$gene<-rownames(markers)
tops<-markers %>% filter(p_val_adj<0.05) %>% arrange(desc(avg_log2FC))%>% head(n=20)
tops<-rbind.data.frame(tops,markers %>% filter(p_val_adj<0.05) %>% arrange(avg_log2FC)%>% head(n=20))

markers$diffexpressed<-ifelse(markers$avg_log2FC>0 & markers$p_val_adj<0.05,'Up',
                              ifelse(markers$avg_log2FC < 0 & markers$p_val_adj<0.05,'Down','NS'))
markers$label<-ifelse(markers$gene %in% tops$gene,markers$gene,NA)

write.csv(markers,paste0('data/MBPM/myeloid/markers_',filename,'_FTL_full.csv'),row.names = F)

pdf(paste0('data/MBPM/myeloid/plots_volcano_',filename,'_FTL.pdf'))
ggplot(data=markers, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=label)) +
  geom_point() + 
  geom_text_repel(size=3) +
  scale_color_manual(values=c('darkblue', 'grey', 'darkred')) +
  labs(color='Log2FC\n(MDM vs\nMDM FTL+)') +
  theme_bw() +
  ggtitle(paste0('MDM vs MDM FTL+'))+
  theme(panel.grid.major = element_blank(),title = element_text(size=10))

ggplot(data=markers, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=label)) +
  geom_point() + 
  scale_color_manual(values=c('darkblue', 'grey', 'darkred')) +
  labs(color='Log2FC\n(MDM vs\nMDM FTL+)') +
  theme_bw() +
  ggtitle(paste0('MDM vs MDM FTL+'))+
  theme(panel.grid.major = element_blank(),title = element_text(size=10))
dev.off()


#### DEG and pathway analysis of MDMs MBM_sn ####
filename<-'MBM_sn_myeloid'

# Subset to monocytes + MDMs
seu<-readRDS('data/cell_type_DEG/MBM_sn/myeloid/data_MBM_sn_myeloid.rds')
seu$cell_type_myeloid<- ifelse(seu$cell_type_fine %in% 
                                 c('MDM M1-like', 'MDM M2-like','Proinflammatory MDM'),
                               'MDM',seu$cell_type_fine)
seu <- subset(seu, cell_type_myeloid %in% c('MDM', 'MDM FTL+'))
seu <- ScaleData(seu)

# Differential gene expression (DGE) based on cell types
Idents(seu)<-seu$cell_type_myeloid
markers <- FindAllMarkers(seu, only.pos = TRUE,assay='RNA',test.use = 'MAST',
                          min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, paste0(directory,'markers_',filename,'.csv'), row.names = F)
markers<- markers %>% filter(p_val_adj < 0.05)


background<-rownames(seu@assays$RNA@data)
pdf(paste0(directory,'markers_',filename,'_pathway.pdf'),width = 10,height = 7)
for(pw in c(C5_BP,HALLMARK,KEGG,REACTOME)){
  for(i in unique(markers$cluster)){
    subs<-markers %>% dplyr::filter(cluster==i) %>% magrittr::use_series(gene)
    print(hyp_dots(hypeR(as.character(subs), genesets=pw, background=background),
                   title = paste0(pw$name,' cluster ',i),abrv=80)+ theme_bw() + 
            theme(axis.title.y = element_blank(), 
                  axis.text.y = element_text(color = 'black'), 
                  axis.text.x = element_text(color = 'black'), 
                  panel.grid.minor = element_blank()) + 
            scale_color_gradient(high = 'black', low = '#C51B8A'))
  }
}
dev.off()


# Markers for volcano plot
Idents(seu)<-seu$cell_type_myeloid
markers <- FindMarkers(seu, ident.1 = 'MDM', ident.2 ='MDM FTL+', min.pct = 0, 
                       logfc.threshold = 0,test.use = 'MAST',
                       assay = 'RNA')
markers$gene<-rownames(markers)
tops<-markers %>% filter(p_val_adj<0.05) %>% arrange(desc(avg_log2FC))%>% head(n=20)
tops<-rbind.data.frame(tops,markers %>% filter(p_val_adj<0.05) %>% arrange(avg_log2FC)%>% head(n=20))

markers$diffexpressed<-ifelse(markers$avg_log2FC>0 & markers$p_val_adj<0.05,'Up',
                              ifelse(markers$avg_log2FC < 0 & markers$p_val_adj<0.05,'Down','NS'))
markers$label<-ifelse(markers$gene %in% tops$gene,markers$gene,NA)

write.csv(markers,paste0('data/MBPM/myeloid/markers_',filename,'_FTL_full.csv'),row.names = F)
#markers$avg_log2FC<-markers$avg_log2FC * (-1) # change order for nicer plot

pdf(paste0('data/MBPM/myeloid/plots_volcano_',filename,'_FTL.pdf'))
ggplot(data=markers, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=label)) +
  geom_point() + 
  geom_text_repel(size=3) +
  scale_color_manual(values=c(colMM[2], 'grey', colMM[1])) +
  labs(color='Log2FC\n(MDM FTL+ vs\nMDM)') +
  theme_bw() +
  ggtitle(paste0('MDM FTL+ vs MDM'))+
  theme(panel.grid.major = element_blank(),title = element_text(size=10))

ggplot(data=markers, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=label)) +
  geom_point() + 
  scale_color_manual(values=c(colMM[2], 'grey', colMM[1])) +
  labs(color='Log2FC\n(MDM FTL+ vs\nMDM)') +
  theme_bw() +
  ggtitle(paste0('MDM FTL+ vs MDM'))+
  theme(panel.grid.major = element_blank(),title = element_text(size=10))
dev.off()


#### DEG and pathway analysis of MDMs MPM_sn ####
filename<-'MPM_sn_myeloid'

# Subset to monocytes + MDMs
seu<-readRDS('data/cell_type_DEG/MPM_sn/myeloid/data_MPM_sn_myeloid.rds')
seu$cell_type_myeloid<- ifelse(seu$cell_type_fine %in% 
                                 c('MDM M1-like', 'MDM M2-like','Proinflammatory MDM'),
                               'MDM',seu$cell_type_fine)
seu <- subset(seu, cell_type_myeloid %in% c('MDM', 'MDM FTL+'))
seu <- ScaleData(seu)

# Differential gene expression (DGE) based on cell types
Idents(seu)<-seu$cell_type_myeloid
markers <- FindAllMarkers(seu, only.pos = TRUE,assay='RNA',test.use = 'MAST',
                          min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, paste0(directory,'markers_',filename,'.csv'), row.names = F)
markers<- markers %>% filter(p_val_adj < 0.05)

background<-rownames(seu@assays$RNA@data)
pdf(paste0(directory,'markers_',filename,'_pathway.pdf'),width = 10,height = 7)
for(pw in c(C5_BP,HALLMARK,KEGG,REACTOME)){
  for(i in unique(markers$cluster)){
    subs<-markers %>% dplyr::filter(cluster==i) %>% magrittr::use_series(gene)
    print(hyp_dots(hypeR(as.character(subs), genesets=pw, background=background),
                   title = paste0(pw$name,' cluster ',i),abrv=80)+ theme_bw() + 
            theme(axis.title.y = element_blank(), 
                  axis.text.y = element_text(color = 'black'), 
                  axis.text.x = element_text(color = 'black'), 
                  panel.grid.minor = element_blank()) + 
            scale_color_gradient(high = 'black', low = '#C51B8A'))
  }
}
dev.off()


# Markers for volcano plot
Idents(seu)<-seu$cell_type_myeloid
markers <- FindMarkers(seu, ident.1 = 'MDM', ident.2 ='MDM FTL+', min.pct = 0, 
                       logfc.threshold = 0,test.use = 'MAST',
                       assay = 'RNA')
markers$gene<-rownames(markers)
tops<-markers %>% filter(p_val_adj<0.05) %>% arrange(desc(avg_log2FC))%>% head(n=20)
tops<-rbind.data.frame(tops,markers %>% filter(p_val_adj<0.05) %>% arrange(avg_log2FC)%>% head(n=20))

markers$diffexpressed<-ifelse(markers$avg_log2FC>0 & markers$p_val_adj<0.05,'Up',
                              ifelse(markers$avg_log2FC < 0 & markers$p_val_adj<0.05,'Down','NS'))
markers$label<-ifelse(markers$gene %in% tops$gene,markers$gene,NA)
write.csv(markers,paste0('data/MBPM/myeloid/markers_',filename,'_FTL_full.csv'),row.names = F)
#markers$avg_log2FC<-markers$avg_log2FC* (-1) # change order for nicer plot

pdf(paste0('data/MBPM/myeloid/plots_volcano_',filename,'_FTL.pdf'))
ggplot(data=markers, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=label)) +
  geom_point() + 
  geom_text_repel(size=3) +
  scale_color_manual(values=c(colMM[2], 'grey', colMM[1])) +
  labs(color='Log2FC\n(MDM FTL+ vs\nMDM)') +
  theme_bw() +
  ggtitle(paste0('MDM FTL+ vs MDM'))+
  theme(panel.grid.major = element_blank(),title = element_text(size=10))

ggplot(data=markers, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=label)) +
  geom_point() + 
  scale_color_manual(values=c(colMM[2], 'grey', colMM[1])) +
  labs(color='Log2FC\n(MDM FTL+ vs\nMDM)') +
  theme_bw() +
  ggtitle(paste0('MDM FTL+ vs MDM'))+
  theme(panel.grid.major = element_blank(),title = element_text(size=10))
dev.off()


