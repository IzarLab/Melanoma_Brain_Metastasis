#!/usr/bin/env Rscript

### title: SCENIC output integration into Seurat object and downstream analysis
### author: Jana Biermann, PhD

library(dplyr)
library(Seurat)
library(ggplot2)
library(gplots)
library(reshape2)
library(viridis)
library(ggrastr)
library(ggpubr)
library(ggrepel)

colBP <- c('#A80D11', '#008DB8')
colSCSN <- c('#E1AC24', '#288F56')

seu<-readRDS('data/tumor/MBPM_sn_tum/data_MBPM_sn_tum.rds')

# Combine individual SCENIC results (output from pySCENIC)
auc<-read.csv(paste0('data/tumor/scenic/auc_mtx_',unique(seu$patient)[1],'.csv'),
              check.names = F)
for(pat in unique(seu$patient)[2:length(unique(seu$patient))]){
  tmp<-read.csv(paste0('data/tumor/scenic/auc_mtx_',pat,'.csv'),check.names = F)
  auc<-full_join(auc,tmp,by='Regulon')
} 
auc[is.na(auc)] <- 0
colnames(auc)<-gsub('.', '-', colnames(auc),fixed = T)
auc$X<-NULL
rownames(auc)<-auc$Regulon
auc$Regulon<-NULL
all(colnames(auc)==rownames(seu@meta.data))
saveRDS(auc,'data/tumor/scenic/data_MBPM_sn_auc_combined.rds')

# Add to seurat object of tumor cells
seu@assays$scenic<-CreateAssayObject(data=as.matrix(auc))
seu@assays$scenic@scale.data<-as.matrix(auc)
DefaultAssay(seu) <- "scenic"
seu@assays$scenic@key <- "scenic_"
seu <- RunPCA(seu, features = rownames(seu@assays$scenic@scale.data),npcs = 20)
ElbowPlot(seu)
seu <- FindNeighbors(seu, dims = 1:20)
seu <- FindClusters(seu, resolution=0.2) 
seu <- RunUMAP(seu, dims = 1:20)

saveRDS(seu,'data/tumor/scenic/data_MBPM_sn_tum_scenic.rds')


# Find differentially expressed TFs
Idents(seu)<-seu$organ
markers <- FindMarkers(seu, ident.1 = 'Brain', ident.2 ='Peripheral',
                          min.pct = 0, logfc.threshold = 0,
                          assay = 'scenic')
markers$gene<-rownames(markers)
tops<-markers %>% filter(p_val_adj<0.05) %>% arrange(desc(avg_log2FC))%>% head(n=20)
tops<-rbind.data.frame(tops,markers %>% filter(p_val_adj<0.05) %>% arrange(avg_log2FC)%>% head(n=20))
markers$diffexpressed<-ifelse(markers$avg_log2FC>0 & markers$p_val_adj<0.05,'Up',
                              ifelse(markers$avg_log2FC < 0 & markers$p_val_adj<0.05,'Down','NS'))
markers$label<-ifelse((markers$gene %in% tops$gene) |
                        grepl('MITF',markers$gene),markers$gene,NA)
write.csv(markers,paste0('data/tumor/scenic/markers_MBPM_sn_scenic.csv'),row.names = F)
markers$p_val_adj<-ifelse(markers$p_val_adj ==0,1e-306,markers$p_val_adj)

# Averaged heatmap
avg<-AverageExpression(seu,assays = 'scenic',return.seurat = T,group.by = 'organ',
                       slot = 'scale.data')
avg$organ<-rownames(avg@meta.data)


### Plots
pdf(paste0('data/tumor/scenic/plots_MBPM_sn_scenic.pdf'),width = 10,height = 10)
# Averaged heatmap
DoHeatmap(avg, features = tops$gene, group.by = 'organ',raster = T,
          assay = 'scenic',draw.lines =F,group.colors = colBP)+
  scale_fill_gradientn(colors = c("white", "black",'black','black'))

# Volcano plot
ggplot(data=markers, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=label)) +
  geom_point() + 
  geom_text_repel(size=2,max.overlaps = 20) +
  scale_color_manual(values=c("darkblue", "grey", "darkred")) +
  labs(color="Log2FC\n(Brain vs\nPeripheral)") +
  theme_bw() +
  ggtitle('Brain vs Peripheral: SCENIC TFs')+
  theme(panel.grid.major = element_blank(),title = element_text(size=10))

# Violin plot
VlnPlot(seu,features = c('MITF(+)','SOX4(+)','RELB(+)','TEAD1(+)','NF1(+)'),
        group.by = 'organ', pt.size = 0,cols = colBP)

df2<-data.frame(patient=seu$patient,
                organ=seu$organ,
                MITF=seu@assays$scenic@data['MITF(+)',],
                SOX4=seu@assays$scenic@data['SOX4(+)',],
                RELB=seu@assays$scenic@data['RELB(+)',],
                TEAD1=seu@assays$scenic@data['TEAD1(+)',],
                NF1=seu@assays$scenic@data['NF1(+)',])
df2_long<-melt(df2[,-1])

ggplot(df2_long, aes(x=organ, y=value,fill=organ)) + 
  geom_violin()+
  theme_bw()+
  scale_fill_manual(values = colBP)+
  ylab('TF expression')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid =element_blank(),
        axis.title.x = element_blank())+
  facet_grid(~variable)

# Box plot
ggplot(df2_long, aes(x=organ, y=value,fill=organ)) + 
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values = colBP)+
  ylab('TF expression')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.grid =element_blank(),
        axis.title.x = element_blank())+
  stat_compare_means(aes(group = organ,label='p.format'), 
                     method = "wilcox.test",size=2) + 
  facet_grid(~variable)

# UMAPs
print(DimPlot(seu, reduction = "umap",label = T,group.by = 'patient',raster = T,shuffle = T))
print(DimPlot(seu, reduction = "umap",label = T,group.by = 'organ',raster = T,shuffle = T))
dev.off()

