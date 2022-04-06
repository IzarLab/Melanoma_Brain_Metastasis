#!/usr/bin/env Rscript

#### title: Tumor cell integration, UMAPs and pathway signatures applied to tumor cells 
#### author: Jana Biermann, PhD

library(dplyr)
library(Seurat)
library(ggplot2)
library(gplots)
'%notin%' <- Negate('%in%')

colBP <- c('#A80D11', '#008DB8')
colSCSN <- c('#E1AC24', '#288F56')
cycling <- c('#5AAE61', '#762A83')
colMITF <- c('#C7B122', '#C70E7B', '#22BAC7')


# Read-in and re-integrate
seu<-readRDS('data/MBPM/data_MBPM_scn.rds')
DefaultAssay(seu)<-'RNA'
seu<-DietSeurat(seu,assays='RNA')
seu<-subset(seu,malignant=='malignant')
seu<-subset(seu,sequencing=='Single nuclei')

# Remove samples with <200 cells and normalize
obj.list <- SplitObject(seu, split.by = 'patient')
for(i in 1:length(obj.list)){
  if(dim(obj.list[[i]]@assays$RNA@counts)[2]<200){
    obj.list[i]<-NA
  }else{
    p=obj.list[[i]]
    p <- NormalizeData(p)
    obj.list[[i]]=p
  }
}
# Remove NAs
obj.list<-obj.list[!is.na(obj.list)]

# Integration
anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:30)
seu <- IntegrateData(anchorset = anchors, dims = 1:30)

# Processing
seu<- ScaleData(seu)
seu <- RunPCA(seu,npcs = 100)
seu <- RunUMAP(seu, dims = 1:25)
seu <- FindNeighbors(seu, dims = 1:25)
seu <- FindClusters(seu, resolution = 0.2)

# Get module scores for sigs
sigs<-read.csv('~/signatures/melanoma_sigs.csv',na.strings = '')
for(c in 1:ncol(sigs)){
  seu<-AddModuleScore(object = seu,features = list(na.omit(sigs[,c])),
                      name = colnames(sigs)[c],assay = 'RNA',search=F)
}

# Cell cycle assignment
DefaultAssay(seu) <- 'RNA'
seu <- CellCycleScoring(seu, s.features = cc.genes.updated.2019$s.genes, 
                        g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
DefaultAssay(seu) <- 'integrated'
seu$cell_cycle <- ifelse(seu$G2M.Score > 0.1 | seu$S.Score > 0.1, 'cycling', 'non-cycling')

# Remove doublets
seu<-subset(seu,integrated_snn_res.0.2 %notin% c(5,6))
seu <- ScaleData(seu) %>% RunPCA() %>% RunUMAP(dims = 1:50)

# Apply signatures
sigs <- read.csv('~/signatures/melanoma_sigs.csv', na.strings = c('', 'NA'))
sigs <- read.csv('~/signatures/differential_sigs_MBPM.csv', 
                 na.strings = c('', 'NA'))
for (c in 1:ncol(sigs)) {
  sig <- as.character(na.omit(sigs[, c]))
  seu <- AddModuleScore(object = seu, features = list(sig), name = colnames(sigs[c]), assay = 'RNA', search = F)
}

final<-c('KEGG_OXIDATIVE_PHOSPHORYLATION','KEGG_INSULIN_SIGNALING_PATHWAY',
         'WP_PI3KAKT_SIGNALING_PATHWAY','KEGG_ERBB_SIGNALING_PATHWAY',
         'WP_EGFR_TYROSINE_KINASE_INHIBITOR_RESISTANCE','PID_PDGFRB_PATHWAY',
         'PID_KIT_PATHWAY','HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
         'KEGG_CELL_ADHESION_MOLECULES_CAMS','HALLMARK_MTORC1_SIGNALING')

msigdbgenes <- read.csv('~/brain_mets/signatures/relevant_msigdb_pathways_hypeR.csv', na.strings = '')
for (s in final) {
  seu <- AddModuleScore(object = seu, features = list(msigdbgenes[, s]), 
                        name = s, assay = 'RNA', search = F)
}

fb<-read.csv('data/NMF/gsea.metagenes.0.3.csv',row.names = 1)
for (c in 1:ncol(fb)) {
  sig <- as.character(na.omit(fb[, c]))
  seu <- AddModuleScore(object = seu, features = list(sig), name = colnames(fb[c]), 
                        assay = 'RNA', search = F)
}


pdf('data/tumor/plots_tumor_umaps.pdf')
# UMAPs
plt <- FeaturePlot(seu, features = c('AXL_sig1', 'MITF_sig1'), blend = T, 
                   cols = colMITF[2:3], raster = T, ncol = 2, combine = F)
CombinePlots(plt)
DimPlot(seu, group.by = 'cell_cycle', shuffle = T, raster = T, cols = cycling)
DimPlot(seu, group.by = 'cell_cycle', shuffle = T, raster = T, cols = cycling) + NoLegend()
DimPlot(seu, group.by = 'organ', shuffle = T, raster = T, cols = colBP)
DimPlot(seu, group.by = 'organ', shuffle = T, raster = T, cols = colBP) + NoLegend()
DimPlot(seu, group.by = 'patient', shuffle = T, raster = T)
DimPlot(seu, group.by = 'patient', shuffle = T, raster = T) + NoLegend()

# MITF + AXL
VlnPlot(seu, features = c('MITF_sig1', 'AXL_sig1'), pt.size = 0, group.by = 'organ',
        cols = colBP)
dev.off()


pdf('data/tumor/plots_umaps_final_sigs_and_kimono.pdf',width = 10,height = 10)
p1<-FeaturePlot(seu, features = c(paste0(final,'1')),min.cutoff = 'q05', 
               max.cutoff = 'q95',raster = T,order = T,combine = F)
p1 <- lapply(X = p1,FUN = function(p) p + theme(title = element_text(size=8)))
CombinePlots(plots = p1,ncol=3)

p2<-FeaturePlot(seu, features = c(paste0(colnames(fb),'1')),min.cutoff = 'q05', 
                max.cutoff = 'q95',raster = T,order = T,combine = F)
p2 <- lapply(X = p2,FUN = function(p) p + theme(title = element_text(size=6)))
CombinePlots(plots = p2, ncol=3)

p3<- VlnPlot(seu, features = c(paste0(final,'1')), pt.size = 0, 
             group.by = 'patient',combine=F)
p3 <- lapply(X = p3,FUN = function(p) p + theme(axis.text.x =element_blank(),
                                                axis.title.x=element_blank(),
                                                title = element_text(size=8)))
CombinePlots(plots = p3, ncol=3,legend='bottom')

p4<- VlnPlot(seu, features = c(paste0(colnames(fb),'1')), pt.size = 0, 
             group.by = 'patient',combine=F)
p4 <- lapply(X = p4,FUN = function(p) p + theme(axis.text.x =element_blank(),
                                                axis.title.x=element_blank(),
                                                title = element_text(size=8)))
CombinePlots(plots = p4, ncol=3,legend='bottom')
dev.off()

