#!/usr/bin/env Rscript

### title: Microglia (MG) analysis (DEG, diffusion component analysis, pathway analysis) 
### author: Jana Biermann, PhD

print(Sys.time())

library(Seurat)
library(destiny)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(gplots)
library(ggpubr)
library(scales)
library(viridis)
library(hypeR)
library(ggpubr)
library(ggrepel)
'%notin%' <- Negate('%in%')


colBP <- c('#A80D11', '#008DB8')
colSCSN <- c('#E1AC24', '#288F56')

# Read-in Seurat object
seu_all <- readRDS('data/MBPM/data_MBPM_scn.rds')

#### Single nuclei workflow ####
seu <- subset(seu_all, cell_type_fine %in% c('Microglia') & sequencing == 'Single nuclei')

# Reintegrate
table(seu$patient) %>% sort
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

# Integrate
anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:20, k.filter = 50)
seu <- IntegrateData(anchorset = anchors, dims = 1:20, k.weight = 50)

# Rerun Seurat workflow
seu <- ScaleData(seu) %>% RunPCA()
ElbowPlot(seu, ndims = 50)
seu <- FindNeighbors(seu, dims = 1:20)
seu <- FindClusters(seu, resolution = 0.5)
seu<-RunUMAP(seu, dims = 1:20)
DimPlot(seu, reduction = 'umap', label = F,shuffle = T, raster=T)
DimPlot(seu, reduction = 'umap', group.by = 'patient',label = F,shuffle = T, raster=T)

# Identify markers for the two MG clusters
markers <- FindAllMarkers(seu, only.pos = TRUE, assay = 'RNA', 
                          min.pct = 0.25, logfc.threshold = 0.25,test.use = 'MAST')
ifelse(!dir.exists(file.path(paste0('data/DEG/microglia/'))), 
       dir.create(file.path(paste0('data/DEG/microglia/')), recursive = T), FALSE)
write.csv(markers, paste0('data/DEG/microglia/markers_MBM_sn_microglia_cluster.csv'), 
          row.names = F)

# Prepare volcano plots
Idents(seu)<-seu$mg_cluster
markers <- FindMarkers(seu, ident.1 = 'Activated MG', ident.2 ='Homeostatic MG', 
                       min.pct = 0, logfc.threshold = 0,assay = 'RNA')

markers$gene<-rownames(markers)
tops<-markers %>% filter(p_val_adj<0.05) %>% arrange(desc(avg_log2FC))%>% head(n=20)
tops<-rbind.data.frame(tops,markers %>% filter(p_val_adj<0.05) %>% 
                         arrange(avg_log2FC)%>% head(n=20))

markers$diffexpressed<-ifelse(markers$avg_log2FC>0 & markers$p_val_adj<0.05,'Up',
                              ifelse(markers$avg_log2FC < 0 & markers$p_val_adj<0.05,'Down','NS'))
markers$label<-ifelse(markers$gene %in% tops$gene,markers$gene,NA)
write.csv(markers,paste0('data/DEG/microglia/markers_MBM_sn_microglia_cluster_full.csv'),
          row.names = F)

# Prepare heatmap
DefaultAssay(seu)<-'RNA'
seu<-NormalizeData(seu)
seu<-ScaleData(seu,assay='RNA')
avg<-AverageExpression(seu,assays = 'RNA',return.seurat = T,
                       group.by = 'mg_cluster',slot = 'scale.data')
avg$mg_cluster<-rownames(avg@meta.data)

### Plots heatmap and volcano
pdf(paste0('data/DEG/microglia/plots_volcano_MBM_sn_microglia.pdf'))
DoHeatmap(avg, features = tops$gene, group.by = 'mg_cluster',raster = T,assay = 'RNA',
          draw.lines =F,group.colors = c('deeppink', 'midnightblue'),
          size = 3,angle = 0,hjust = 0.5)+
  scale_fill_gradient2(low = "#FF00FF",high = "#FFFF00",mid = "#000000",midpoint = 0)+
  theme(plot.margin = margin(0, 8, 0, 0, "cm"),
        axis.text.y = element_text(colour = 'black',face = 'italic'))

ggplot(data=markers, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=label)) +
  geom_point() + 
  geom_text_repel(size=3) +
  scale_color_manual(values=c("darkblue", "grey", "darkred")) +
  labs(color="Log2FC\n(Activated vs\nHomeostatic)") +
  theme_bw() +
  coord_cartesian(xlim = c(-2.5, 2.5))+
  ggtitle(paste0('MBM_sn: Activated vs Homeostatic MG'))+
  theme(panel.grid.major = element_blank(),title = element_text(size=10))

ggplot(data=markers, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=label)) +
  geom_point() + 
  scale_color_manual(values=c("darkblue", "grey", "darkred")) +
  labs(color="Log2FC\n(Activated vs\nHomeostatic)") +
  theme_bw() +
  coord_cartesian(xlim = c(-2.5, 2.5))+
  ggtitle(paste0('MBM_sn: Activated vs Homeostatic MG'))+
  theme(panel.grid.major = element_blank(),title = element_text(size=10))

dev.off()

# Rename clusters
seu$mg_cluster <- seu$integrated_snn_res.0.5
seu$mg_cluster <- ifelse(seu$mg_cluster == '0', 'Homeostatic MG', 'Activated MG')

# Save object
ifelse(!dir.exists(file.path(paste0('data/MBPM/myeloid/microglia'))), 
       dir.create(file.path(paste0('data/MBPM/myeloid/microglia')), recursive = T), FALSE)
saveRDS(seu,'data/MBPM/myeloid/microglia/data_mg_sn.rds')


### Diffusion component analysis
set.seed(1)
es <- as.ExpressionSet(as.data.frame(t(seu@assays$integrated@data)))
es@phenoData@data <- seu@meta.data
dm <- DiffusionMap(es, verbose = T, n_pcs = 20)

### Plots diffusion maps
pdf('data/MBPM/myeloid/microglia/plots_dm_microglia_sn.pdf')
par(mar = c(5.1, 4.1, 4.1, 2.1), xpd = TRUE)
palette(c('deeppink', 'midnightblue'))
tmp_var<-seu$mg_cluster
plot(dm, col = as.factor(tmp_var), main = 'mg_cluster', pch = 20)
legend('bottomright',  legend = levels(as.factor(tmp_var)), 
       pch = 16, col = as.factor(levels(as.factor(tmp_var))), bty = 'n',cex=0.5)

tmp_var<-seu$patient
palette(hue_pal()(length(unique(tmp_var))))
plot(dm, col = as.factor(tmp_var), main = 'patient', pch = 20)
legend('bottomright',  legend = levels(as.factor(tmp_var)), 
       pch = 16, col = as.factor(levels(as.factor(tmp_var))), bty = 'n',cex=0.5)

par(mar = c(5.1, 4.1, 4.1, 7), xpd = TRUE)
palette(viridis(100))
plot(dm, col = seu@assays$RNA@data['SPP1',], main = 'SPP1', pch = 20)
colorlegend(seu@assays$RNA@data['SPP1',], viridis(length(seu@assays$RNA@data['SPP1',])), 
            posx = c(0.88, 0.9), posy = c(0, 0.65))

set.seed(1)
dm_df<-as.data.frame(dm@eigenvectors)[sample(nrow(as.data.frame(dm@eigenvectors))),1:2]
ggplot(dm_df, aes(DC1, DC2, col = seu@meta.data[rownames(dm_df), 'mg_cluster'])) + 
  geom_point(alpha = 1, shape = 16, size = 1.5) + theme_classic() + 
  ggtitle('mg_cluster') + theme(legend.title = element_blank())

ggplot(dm_df, aes(DC1, DC2, col = seu@meta.data[rownames(dm_df), 'patient'])) + 
  geom_point(alpha = 1, shape = 16, size = 1.5) + theme_classic() + 
  ggtitle('patient') + theme(legend.title = element_blank())

ggplot(dm_df, aes(DC1, DC2, col = seu@assays$integrated@data['SPP1', rownames(dm_df)])) + 
  geom_point(alpha = 1, shape = 16, size = 1.5) + theme_classic() + 
  ggtitle('SPP1') + theme(legend.title = element_blank()) + scale_color_viridis()

dev.off()


# Pathway analysis
markers<-read.csv(paste0('data/DEG/microglia/markers_MBM_sn_microglia_cluster.csv'))
markers<- markers %>% filter(p_val_adj < 0.05)

HALLMARK <- msigdb_gsets(species="Homo sapiens", category="H")
KEGG     <- msigdb_gsets(species="Homo sapiens", category="C2", subcategory="CP:KEGG")
REACTOME <- msigdb_gsets(species="Homo sapiens", category="C2", subcategory="CP:REACTOME")
C5_BP <- msigdb_gsets(species="Homo sapiens", category="C5", subcategory="BP")

background<-rownames(seu@assays$RNA@data)

pdf(paste0('data/MBPM/myeloid/microglia/markers_microglia_sn_pathway.pdf'),
    width = 10,height = 7)
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



#### Single cell workflow #####
seu <- subset(seu_all, cell_type_fine == 'Microglia' & sequencing == 'Single cells')

# Re-integrate
table(seu$patient) %>% sort
DefaultAssay(seu) <- 'RNA'
obj.list <- SplitObject(seu, split.by = 'patient')
anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:20, k.filter = 50)
seu <- IntegrateData(anchorset = anchors, dims = 1:20,k.weight = 50)

# Rerun Seurat workflow
seu <- ScaleData(seu) %>% RunPCA() %>%RunUMAP(dims = 1:20) %>% 
  FindNeighbors(dims = 1:20) %>% FindClusters(resolution = 0.2)
DimPlot(seu, reduction = 'umap', label = F,shuffle = T, raster=T)
DimPlot(seu, reduction = 'umap', group.by = 'patient',label = F,shuffle = T, raster=T)

# Identify markers for the two MG clusters
markers <- FindAllMarkers(seu, only.pos = TRUE, assay = 'RNA', 
                          min.pct = 0.25, logfc.threshold = 0.25,test.use = 'MAST')
ifelse(!dir.exists(file.path(paste0('data/DEG/microglia/'))), 
       dir.create(file.path(paste0('data/DEG/microglia/')), recursive = T), FALSE)
write.csv(markers, paste0('data/DEG/microglia/markers_MBM_sc_microglia_cluster.csv'), 
          row.names = F)

# Rename clusters
seu$mg_cluster <- Idents(seu)
seu$mg_cluster <- ifelse(seu$mg_cluster == '1', 'Homeostatic MG', 'Activated MG')

# Save object
saveRDS(seu,'data/MBPM/myeloid/microglia/data_mg_sc.rds')


# Prepare volcano plots
Idents(seu)<-seu$mg_cluster
markers <- FindMarkers(seu, ident.1 = 'Activated MG', ident.2 ='Homeostatic MG', 
                       min.pct = 0, logfc.threshold = 0,assay = 'RNA')
markers$gene<-rownames(markers)
tops<-markers %>% filter(p_val_adj<0.05) %>% arrange(desc(avg_log2FC))%>% head(n=20)
tops<-rbind.data.frame(tops,markers %>% filter(p_val_adj<0.05) %>% 
                         arrange(avg_log2FC)%>% head(n=20))

markers$diffexpressed<-ifelse(markers$avg_log2FC>0 & markers$p_val_adj<0.05,'Up',
                              ifelse(markers$avg_log2FC < 0 & markers$p_val_adj<0.05,'Down','NS'))
markers$label<-ifelse(markers$gene %in% tops$gene,markers$gene,NA)
write.csv(markers,paste0('data/DEG/microglia/markers_MBM_sc_microglia_cluster_full.csv'),
          row.names = F)
markers<-read.csv(paste0('data/DEG/microglia/markers_MBM_sc_microglia_cluster_full.csv'))

# Prepare heatmap
DefaultAssay(seu)<-'RNA'
seu<-NormalizeData(seu)
seu<-ScaleData(seu,assay='RNA')
avg<-AverageExpression(seu,assays = 'RNA',return.seurat = T,
                       group.by = 'mg_cluster',slot = 'scale.data')
avg$mg_cluster<-rownames(avg@meta.data)

### Plots heatmap and volcano
pdf(paste0('data/DEG/microglia/plots_volcano_MBM_sc_microglia.pdf'))
DoHeatmap(avg, features = tops$gene, group.by = 'mg_cluster',raster = T,assay = 'RNA',
          draw.lines =F,group.colors = c('deeppink', 'midnightblue'),
          size = 3,angle = 0,hjust = 0.5)+
  scale_fill_gradient2(low = "#FF00FF",high = "#FFFF00",mid = "#000000",midpoint = 0)+
  theme(plot.margin = margin(0, 8, 0, 0, "cm"),
        axis.text.y = element_text(colour = 'black',face = 'italic'))

ggplot(data=markers, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=label)) +
  geom_point() + 
  geom_text_repel(size=3) +
  scale_color_manual(values=c("darkblue", "grey", "darkred")) +
  labs(color="Log2FC\n(Activated vs\nHomeostatic)") +
  theme_bw() +
  ggtitle(paste0('MBM_sc: Activated vs Homeostatic MG'))+
  theme(panel.grid.major = element_blank(),title = element_text(size=10))

ggplot(data=markers, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=label)) +
  geom_point() + 
  scale_color_manual(values=c("darkblue", "grey", "darkred")) +
  labs(color="Log2FC\n(Activated vs\nHomeostatic)") +
  theme_bw() +
  ggtitle(paste0('MBM_sc: Activated vs Homeostatic MG'))+
  theme(panel.grid.major = element_blank(),title = element_text(size=10))

dev.off()


### Diffusion component analysis
set.seed(1)
es <- as.ExpressionSet(as.data.frame(t(seu@assays$integrated@data)))
es@phenoData@data <- seu@meta.data
dm <- DiffusionMap(es, verbose = T, n_pcs = 20)

### Plots diffusion maps
pdf('data/MBPM/myeloid/microglia/plots_dm_microglia_sc.pdf')
par(mar = c(5.1, 4.1, 4.1, 2.1), xpd = TRUE)
palette(c('deeppink', 'midnightblue'))
tmp_var<-seu$mg_cluster
plot(dm, col = as.factor(tmp_var), main = 'mg_cluster', pch = 20)
legend('bottomright',  legend = levels(as.factor(tmp_var)), 
       pch = 16, col = as.factor(levels(as.factor(tmp_var))), bty = 'n',cex=0.5)

tmp_var<-seu$patient
palette(hue_pal()(length(unique(tmp_var))))
plot(dm, col = as.factor(tmp_var), main = 'patient', pch = 20)
legend('bottomright',  legend = levels(as.factor(tmp_var)), 
       pch = 16, col = as.factor(levels(as.factor(tmp_var))), bty = 'n',cex=0.5)

par(mar = c(5.1, 4.1, 4.1, 7), xpd = TRUE)
palette(viridis(100))
plot(dm, col = seu@assays$RNA@data['SPP1',], main = 'SPP1', pch = 20)
colorlegend(seu@assays$RNA@data['SPP1',], viridis(length(seu@assays$RNA@data['SPP1',])), 
            posx = c(0.88, 0.9), posy = c(0, 0.65))

set.seed(1)
dm_df<-as.data.frame(dm@eigenvectors)[sample(nrow(as.data.frame(dm@eigenvectors))),1:2]
ggplot(dm_df, aes(DC1, DC2, col = seu@meta.data[rownames(dm_df), 'mg_cluster'])) + 
  geom_point(alpha = 1, shape = 16, size = 1.5) + theme_classic() + 
  ggtitle('mg_cluster') + theme(legend.title = element_blank())

ggplot(dm_df, aes(DC1, DC2, col = seu@meta.data[rownames(dm_df), 'patient'])) + 
  geom_point(alpha = 1, shape = 16, size = 1.5) + theme_classic() + 
  ggtitle('patient') + theme(legend.title = element_blank())

ggplot(dm_df, aes(DC1, DC2, col = seu@assays$integrated@data['SPP1', rownames(dm_df)])) + 
  geom_point(alpha = 1, shape = 16, size = 1.5) + theme_classic() + 
  ggtitle('SPP1') + theme(legend.title = element_blank()) + scale_color_viridis()
dev.off()


### Pathway analysis
HALLMARK <- msigdb_gsets(species="Homo sapiens", category="H")
KEGG     <- msigdb_gsets(species="Homo sapiens", category="C2", subcategory="CP:KEGG")
REACTOME <- msigdb_gsets(species="Homo sapiens", category="C2", subcategory="CP:REACTOME")
C5_BP <- msigdb_gsets(species="Homo sapiens", category="C5", subcategory="BP")

background<-rownames(seu@assays$RNA@data)

pdf('data/MBPM/myeloid/microglia/plots_microglia_sc_pathways.pdf',width = 10,height = 7)
for(pw in c(C5_BP,HALLMARK,KEGG,REACTOME)){
  for(i in unique(markers$cluster)){
    subs<-markers %>% dplyr::filter(cluster==i) %>% magrittr::use_series(gene)
    print(hyp_dots(hypeR(as.character(subs), genesets=pw, background=background),
                   title = paste0(pw$name,' cluster ',i),abrv=80))
  }
}
dev.off()
