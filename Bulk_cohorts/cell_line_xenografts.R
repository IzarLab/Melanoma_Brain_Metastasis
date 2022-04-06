#!/usr/bin/env Rscript

#### Title: DEG, heatmap and signature scores of our MBM/ECM tumor signatures on xenografts 
#### Author: Jana Biermann, PhD

library(dplyr)
library(ggplot2)
library(gplots)
library(singscore)
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(viridis)
library(reshape2)
library(ggpubr)
'%notin%' <- Negate('%in%')

colBP <- c('#A80D11', '#008DB8')
colSCSN <- c('#E1AC24', '#288F56')
colBP_light<-c('#FF2D34','#018FFA')


#### DESeq2
gex <- read.csv('data/NYU_data/5b1.4l.counts.raw.csv')
rownames(gex) <- gex$gene
gex$gene <- NULL
cts <- as.matrix(gex)

meta <- read.csv('data/NYU_data/metadata_xenograft.csv')
meta<-subset(meta,injection=='IC')
coldata <- meta
rownames(coldata)<-coldata$sample_name
coldata$organ <- factor(coldata$organ)
coldata$cell_line <- factor(coldata$cell_line)
coldata$group <- factor(coldata$group)
cts <- cts[,rownames(coldata)]

all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ cell_line+organ)
dds
dds$organ <- relevel(dds$organ, ref = "Peripheral")

# DEG
dds <- DESeq(dds)
res <- results(dds)
res
resOrdered <- res[order(res$pvalue),]
summary(res)
sum(res$padj < 0.05, na.rm=TRUE)
write.csv(as.data.frame(resOrdered), file=paste0("data/NYU_data/data_BvsP_DESeq2.csv"))
plotMA(res)

# Sample distance heatmap
vsd <- vst(dds, blind=FALSE)
distmat<-as.matrix(dist(t(assay(vsd))))
colnames(distmat)<-meta$sample
rownames(distmat)<-meta$sample
distmat_melt<-reshape2::melt(distmat)
distmat_melt$organ<-ifelse(grepl('brain',distmat_melt$Var1),'Brain','ECM')
distmat_melt$tissue<-ifelse(grepl('lung',distmat_melt$Var1),'Lung',
                            ifelse(grepl('liver',distmat_melt$Var1),'Liver','Brain'))


### Signatures
# Read-in DEGs from our cohort
markers_mbm <- read.csv('data/DEG/organ/tumor/markers_MBPM_sn_tumor_organ_full.csv')
mbm_markers <- data.frame(gene = markers_mbm$gene, mbm_logFC = markers_mbm$avg_log2FC, 
                          mbm_adjp = markers_mbm$p_val_adj)
mbm_markers <- mbm_markers %>% arrange(desc(mbm_logFC))%>% filter(mbm_adjp<0.05)

# Apply SingScore on bulk data
rankData <- rankGenes(counts(dds))
scoredf <- simpleScore(rankData, upSet = mbm_markers$gene[1:100], 
                       downSet = mbm_markers$gene[(nrow(mbm_markers)-100) :nrow(mbm_markers)])
scoredf$sample_name <- rownames(scoredf)
scoredf <- left_join(scoredf, meta[, c('sample_name', 'organ')], by = 'sample_name')
scoredf$organ <- factor(scoredf$organ, levels = c("Brain", "Peripheral"))



### Plots
pdf(paste0('data/NYU_data/plots_cellline_xenograft_BvsP.pdf'))
# Heatmap
anno <- data.frame(organ=meta$organ,
                   tissue=meta$group,
                   cell_line=meta$cell_line)
rownames(anno)<-meta$sample

hc<-hclust(dist(distmat),method = 'complete')
anno <- anno[hc$order,]
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

ann_colors = list(
  organ = c(`Brain`=colBP[1],`Peripheral`= colBP[2]),
  tissue = c(`brain`=colBP[1],`lung`= '#8214A0',`liver`= '#0AB45A'),
  cell_line = c(`X4L`='#F2B701',`X5B1`= '#E73F74'))

pheatmap(distmat,cluster_rows = T, cluster_cols = T, annotation_col = anno,
         show_colnames = F,col=colors,main = 'Sample distance matrix',
         annotation_colors =ann_colors)

# Signatures
box1<-ggboxplot(scoredf, x = 'organ', y = 'UpScore', color = 'organ', add = 'jitter') + 
  stat_compare_means(aes(group = organ, label = sprintf('p = %5.7f', as.numeric(..p.format..))), 
                     method = 'wilcox.test', size = 4) + color_palette(palette = colBP_light) + 
  theme(legend.position = 'none',axis.title.x = element_blank())+ 
  ylab('MBM signature up-score')
box2<-ggboxplot(scoredf, x = 'organ', y = 'DownScore', color = 'organ', add = 'jitter') + 
  stat_compare_means(aes(group = organ, label = sprintf('p = %5.7f', as.numeric(..p.format..))), 
                     method = 'wilcox.test', size = 4) + color_palette(palette = colBP_light) + 
  theme(legend.position = 'none',axis.title.x = element_blank())+
  ylab('MPM signature down-score')
box1+box2

dev.off()
