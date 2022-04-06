#!/usr/bin/env Rscript

### title: Correlation of cell line RNA-seq and scRNA-seq data for matched patients
### author: Jana Biermann, PhD

library(dplyr)
library(Seurat)
library(ggplot2)
library(gplots)
library(ggrepel)
library(viridis)
library(DESeq2)
library(ggrepel)
library(ggVennDiagram)
'%notin%' <- Negate('%in%')

#### Choose one patient ####
# MBM03_sc == BI3
pat<- 'MBM03_sc'
BI<- 'BI3'

# MBM05_sc == BI5
pat<- 'MBM05_sc'
BI<- 'BI5'
#############################


# Read-in scRNA-seq data
seu<-readRDS('data/MBPM/data_MBPM_scn.rds')
sc<-subset(seu, patient == pat)
DefaultAssay(sc)<-'RNA'
sc_avg<-data.frame(sc_avg=apply(sc@assays$RNA@counts, 1, mean))
sc_avg$gene<-rownames(sc_avg)

# Read-in cell line bulk RNA-seq data
bulk<-read.csv('data/celllines_RNAseq/data_counts_combined.csv')
rownames(bulk)<-bulk$gene
bulk<-bulk[,grepl(BI,colnames(bulk))]
bulk_avg<-data.frame(bulk_avg=apply(bulk, 1, mean))
bulk_avg$gene<-rownames(bulk_avg)

# Combine
comb_sc<-inner_join(sc_avg,bulk_avg,by='gene')
comb_sc$label<-ifelse((log2(comb_sc$sc_avg+1)>3) | (log2(comb_sc$bulk_avg+1) >15), 
                      comb_sc$gene,NA)
corr_sc <- cor.test(comb_sc$sc_avg, comb_sc$bulk_avg, method = 'spearman')


pdf(paste0('data/celllines_RNAseq/plots_',pat,'_comparisons.pdf'))
ggplot(data=comb_sc, aes(x=log2(sc_avg+1), y=log2(bulk_avg+1),label=label)) +
  geom_point() + 
  ggrepel::geom_label_repel(aes(label = label), fill = 'white', size = 2, box.padding = 0.1, 
                            label.padding = 0.1, max.overlaps = 17) + 
  theme_bw() +
  ggtitle(paste0(pat,' scRNA-seq vs cell line bulk RNA-seq\nRho: ', 
                 round(corr_sc$estimate[[1]], 3), ', cor P = ', round(corr_sc$p.value,5)))+
  theme(panel.grid.major = element_blank(),title = element_text(size=10))

ggplot(data=comb_sc, aes(x=log2(sc_avg+1), y=log2(bulk_avg+1),label=label)) +
  geom_point() + 
  theme_bw() +
  ggtitle(paste0(pat,' scRNA-seq vs cell line bulk RNA-seq\nRho: ', 
                 round(corr_sc$estimate[[1]], 3), ', cor P = ', round(corr_sc$p.value,5)))+
  theme(panel.grid.major = element_blank(),title = element_text(size=10))
dev.off()

