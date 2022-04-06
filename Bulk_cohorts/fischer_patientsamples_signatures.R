#!/usr/bin/env Rscript

#### Title: Scoring markers from our MBM/ECM tumor signatures on bulk patient samples 
#### Author: Jana Biermann, PhD

library(dplyr)
library(ggplot2)
library(gplots)
library(singscore)
library(DESeq2)
'%notin%' <- Negate('%in%')

colBP <- c('#A80D11', '#008DB8')
colSCSN <- c('#E1AC24', '#288F56')

folder<-'Fischer_patientsamples'


#### DESeq2  ####
gex <- read.csv(paste0('data/',folder,'/davies_2019_htseq_counts.csv'))
gex<-gex[gex$hgnc_symbol != '',]
rownames(gex) <- gex$hgnc_symbol
gex <- gex[, - c(1:4)]
cts <- as.matrix(gex)

meta <- read.csv(paste0('data/',folder,'/davies_2019_sample_info.csv'))
meta<-meta[meta$rna_seq_id %in% colnames(gex),]
meta$organ<-ifelse(meta$accession_site == 'Brain','Brain','ECM')
table(meta$patient_id,meta$organ)
coldata <- meta[,c('rna_seq_id','organ','sample_id','patient_id','accession_site')]
rownames(coldata)<-coldata$rna_seq_id
coldata$organ <- factor(coldata$organ)
coldata$sample_id <- factor(coldata$sample_id)
coldata$patient_id <- factor(coldata$patient_id)

head(cts,2)
coldata
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ patient_id+organ)
dds
dds$organ <- relevel(dds$organ, ref = 'ECM')


# DEG
dds <- DESeq(dds)
saveRDS(dds,paste0('data/',folder,'/data_dds_fischer_patientsamples.rds'))

res <- results(dds)
res
resOrdered <- res[order(res$pvalue),]
summary(res)
sum(res$pvalue < 0.05, na.rm=TRUE)
sum(res$padj < 0.05, na.rm=TRUE)
write.csv(as.data.frame(resOrdered), 
          file=paste0('data/',folder,'/data_fischer_patientsamples_BvsP_DESeq2.csv'))
plotMA(res)
plotCounts(dds, gene=which.min(res$padj), intgroup='organ')

# Read-in DEGs from our cohort
markers_mbm <- read.csv('data/DEG/organ/tumor/markers_MBPM_sn_tumor_organ_full.csv')
mbm_markers <- data.frame(gene = markers_mbm$gene, mbm_logFC = markers_mbm$avg_log2FC, 
                          mbm_adjp = markers_mbm$p_val_adj)
mbm_markers <- mbm_markers %>% arrange(desc(mbm_logFC))%>% filter(mbm_adjp<0.05)
mbm_up<-mbm_markers %>% filter(mbm_logFC>0) %>% 
  top_n(n = 100,mbm_logFC) %>% dplyr::select(gene)
mbm_dn<-mbm_markers %>% filter(mbm_logFC< 0) %>% 
  top_n(n = -100,mbm_logFC) %>% arrange(mbm_logFC) %>% dplyr::select(gene)


# Apply SingScore on bulk data
rankData <- rankGenes(counts(dds))
scoredf <- simpleScore(rankData, upSet = mbm_up$gene[1:100], downSet = mbm_dn$gene[1:100])
scoredf$rna_seq_id <- rownames(scoredf)
scoredf <- left_join(scoredf, as.data.frame(colData(dds))[, c('rna_seq_id', 'organ')], 
                     by = 'rna_seq_id')
scoredf$organ <- factor(scoredf$organ, levels = c('Brain', 'ECM'))


# Plots
pdf(paste0('data/',folder,'/plots_fischer_patientsamples_signature_scores_BvsP_DESeq2.pdf'))
box1<-ggboxplot(scoredf, x = 'organ', y = 'UpScore', color = 'organ', add = 'jitter') + 
  stat_compare_means(aes(group = organ, label = sprintf('p = %5.7f', as.numeric(..p.format..))), 
                     method = 'wilcox.test', size = 4) + color_palette(palette = colBP) + 
  theme(legend.position = 'none',axis.title.x = element_blank())+ 
  ylab('MBM signature up-score')
box2<-ggboxplot(scoredf, x = 'organ', y = 'DownScore', color = 'organ', add = 'jitter') + 
  stat_compare_means(aes(group = organ, label = sprintf('p = %5.7f', as.numeric(..p.format..))), 
                     method = 'wilcox.test', size = 4) + color_palette(palette = colBP) + 
  theme(legend.position = 'none',axis.title.x = element_blank())+
  ylab('MPM signature down-score')
box1+box2
dev.off()

