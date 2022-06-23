#!/usr/bin/env Rscript

### title: Differential gene expression (DGE) of tumor cells comparing MBM and MPM
### author: Jana Biermann, PhD

print(Sys.time())

library(dplyr)
library(Seurat)
library(pheatmap)
library(ggplot2)
library(gplots)
library(scales)
library(viridis)
library(ggrastr)
library(ggrepel)

colBP <- c('#A80D11', '#008DB8')
colSCSN <- c('#E1AC24', '#288F56')
type<- 'tumor'

# Read-in integrated dataset
seu <- readRDS('data/MBPM/data_MBPM_scn.rds')

# Subset to sn tumor cells
seu <- subset(seu, sequencing == 'Single nuclei' & cell_type_main == 'Tumor cells')


#### DGE
# markers
Idents(seu)<-seu$organ
markers <- FindMarkers(seu, ident.1 = 'Brain', ident.2 ='Peripheral', min.pct = 0, 
                       logfc.threshold = 0,test.use = 'MAST',
                       assay = 'RNA', max.cells.per.ident = min(table(seu$organ)))
markers$gene<-rownames(markers)
tops<-markers %>% filter(p_val_adj<0.05) %>% arrange(desc(avg_log2FC))%>% head(n=20)
tops<-rbind.data.frame(tops,markers %>% filter(p_val_adj<0.05) %>% arrange(avg_log2FC)%>% head(n=20))

markers$diffexpressed<-ifelse(markers$avg_log2FC>0 & markers$p_val_adj<0.05,'Up',
                              ifelse(markers$avg_log2FC < 0 & markers$p_val_adj<0.05,'Down','NS'))
markers$label<-ifelse(markers$gene %in% tops$gene,markers$gene,NA)

ifelse(!dir.exists(file.path(paste0('data/DEG/organ/',type))), 
       dir.create(file.path(paste0('data/DEG/organ/',type)),recursive = T), FALSE)
write.csv(markers,paste0('data/DEG/organ/',type,'/markers_MBPM_sn_',type,'_organ_full.csv'),row.names = F)
markers$p_val_adj<-ifelse(markers$p_val_adj<1e-399,1e-399,markers$p_val_adj)

# Apply markers as signatures
mbm_sig<-markers %>% 
  arrange(desc(avg_log2FC))%>% 
  head(n=100)
seu <- AddModuleScore(object = seu, features = list(mbm_sig$gene), name = 'mbm_sig', 
                      assay = 'RNA', search = F)

mpm_sig<-markers %>% 
  arrange(desc(avg_log2FC))%>% 
  tail(n=100)
seu <- AddModuleScore(object = seu, features = list(mpm_sig$gene), name = 'mpm_sig', 
                      assay = 'RNA', search = F)

# Averaged heatmap
DefaultAssay(seu) <- 'RNA'
seu <- NormalizeData(seu)
seu <- ScaleData(seu, assay = 'RNA')
avg <- AverageExpression(seu, assays = 'RNA', return.seurat = T, group.by = 'organ', slot = 'scale.data')
avg$organ <- rownames(avg@meta.data)

# Plot volcano
pdf('data/DEG/organ/tumor/plots_avg_heatmap_MBPM_sn_tumor_organ.pdf')
DoHeatmap(avg, features = tops$gene, group.by = 'organ', raster = T, 
          assay = 'RNA', draw.lines = F, group.colors = colBP)

ggplot(data=markers, aes(x=avg_log2FC, y=-log10(p_val_adj), col=diffexpressed, label=label)) +
  geom_point() + 
  geom_text_repel(size=3) +
  scale_color_manual(values=c("darkblue", "grey", "darkred")) +
  labs(color="Log2FC\n(Brain vs\nPeripheral)") +
  theme_bw() +
  ggtitle(paste0('Brain vs Peripheral: Tumor cells'))+
  theme(panel.grid.major = element_blank(),title = element_text(size=10))
dev.off()


## Large heatmaps
# Top 20,30,40,50 DEGs plus additional genes
for(ngenes in c(200)){
  # Top DEGs
  tops<-markers %>% 
    filter(p_val_adj<0.05) %>% 
    arrange(desc(avg_log2FC))%>% 
    head(n=ngenes/2)
  
  additionals<-c('MET','PMEL','ST6GALNAC3','DAAM2','SAMHD1','NRG3','NRCAM',
                 'NCAM1', 'NELL1','LRRC1', 'AKAP12', 'CADPS2','TBXAS1', 'WNK4')
  for(gene in additionals){
    if(gene %notin% tops$gene){
      tops<-rbind.data.frame(tops,markers[markers$gene==gene,])
    }
  }
  
  tops<-rbind.data.frame(tops,markers %>% 
                           filter(p_val_adj<0.05) %>% 
                           arrange(avg_log2FC)%>% 
                           head(n=ngenes/2))
  
  # Cluster anno
  anno <- data.frame(organ=seu$organ,
                     patient=seu$patient)
  rownames(anno)<-rownames(seu@meta.data)
  ann_colors = list(organ = c(`Brain`=colBP[1],`Peripheral`= colBP[2]))
  
  # Matrix
  mat<- as.matrix(seu@assays$RNA@scale.data[tops$gene,])
  mat[mat > 2] <- 2 
  mat<-mat[,names(sort(seu$mbm_sig1,decreasing = T))]
  anno <- anno[colnames(mat),]
  
  tiff(filename =paste0('data/tumor/MBPM_sn_tum/plot_heatmap_MBPM_sn_BvsP_',ngenes,'genes_plusAdditional.tiff'),
       width=600,height=1400)
  print(pheatmap(mat,cluster_rows = F, cluster_cols = F, annotation_col = anno,
                 show_colnames = F,col=PurpleAndYellow(),
                 main = paste0('BvsP ',ngenes,' DEGs + additional genes\n(Tumor cells, MBPM_sn)'),
                 annotation_colors =ann_colors,scale = 'none'))
  dev.off()
}


print(Sys.time())
