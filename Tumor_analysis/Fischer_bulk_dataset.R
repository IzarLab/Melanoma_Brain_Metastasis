#!/usr/bin/env Rscript

#### Title: Scoring markers from our MBM/MPM tumor signatures on Fischer bulk signatures 
#### Author: Jana Biermann, PhD

library(dplyr)
library(Seurat)
library(ggplot2)
library(gplots)
library(singscore)
library(tidyr)
library(reshape2)
library(ggpubr)
library(ggVennDiagram)
'%notin%' <- Negate('%in%')

colBP <- c('#A80D11', '#008DB8')
colSCSN <- c('#E1AC24', '#288F56')


#### Score MBM/MPM signatures from our cohort on external dataset ####
# Reformat gene expression table from Fischer et al.
gex <- read.table('data/Davies/voom.logCPM_138.txt', header = T, sep = '\t')
gex <- gex[-1, ]
gex <- gex[, -1]
colnames(gex)[1] <- 'gene'
gex <- gex[!is.na(gex$gene), ]
gex <- gex[gex$gene != '43161', ]
gex <- gex[gex$gene != '43160', ]
rownames(gex) <- gex$gene
gex$gene <- NULL

# Reformat clinical table from Fischer et al.
clin <- read.csv('data/Davies/clinical_data.csv')
clin$name <- paste0('MC.', clin$RNA.seq.ID)
clin <- clin[clin$name %in% colnames(gex), ]
table(clin$Accession.Site)
clin$extracranial <- ifelse(clin$Accession.Site == 'Brain', 'BM', 'ECM')
table(clin$extracranial)

# Read-in top 100 markers
markers <- read.csv('signatures/melanoma_makers.csv')
markers_b <- markers$MBM_top100
markers_p <- markers$MPM_top100

# Apply SingScore on bulk data
rankData <- rankGenes(gex)
scoredf <- simpleScore(rankData, upSet = markers_b, downSet = markers_p)
scoredf$name <- rownames(scoredf)
scoredf <- left_join(scoredf, clin[, c('name', 'extracranial')], by = 'name')

# Plot
pdf('data/tumor/plots_score_MBPM_sigs_on_bulk.pdf')
ggboxplot(scoredf, x = 'extracranial', y = 'UpScore', color = 'extracranial', add = 'jitter') + 
  stat_compare_means(aes(group = extracranial, label = sprintf('p = %5.7f', as.numeric(..p.format..))), 
                     method = 'wilcox.test', size = 4) + color_palette(palette = colBP) + 
  theme(legend.position = 'none')
ggboxplot(scoredf, x = 'extracranial', y = 'DownScore', color = 'extracranial', add = 'jitter') + 
  stat_compare_means(aes(group = extracranial, label = sprintf('p = %5.7f', as.numeric(..p.format..))), 
                     method = 'wilcox.test', size = 4) + color_palette(palette = colBP) + 
  theme(legend.position = 'none')
dev.off()


#### Compare log2FCs in DEGs of our cohort and Fischer et al. ####
# Read-in Fischer DEGs
dav_deg <- read.csv('data/DEG/data_DEG_DESeq2_100_logFC.csv')
davies <- data.frame(gene = dav_deg$gene, davies_logFC = dav_deg$avg_log2FC, 
                     davies_adjp = dav_deg$p_val_adj)
davies <- davies %>% arrange(desc(davies_logFC))

# Read-in DEGs from our cohort
markers_mbm <- read.csv('data/DEG/organ/tumor/markers_MBPM_sn_tumor_organ_full.csv')
mbm_markers <- data.frame(gene = markers_mbm$gene, mbm_logFC = markers_mbm$avg_log2FC, 
                          mbm_adjp = markers_mbm$p_val_adj)
mbm_markers_dn <- mbm_markers %>% arrange(desc(mbm_logFC)) %>% top_n(-100, wt = mbm_logFC)
mbm_markers_up <- mbm_markers %>% arrange(desc(mbm_logFC)) %>% top_n(100, wt = mbm_logFC)
mbm_markers <- rbind.data.frame(mbm_markers_up, mbm_markers_dn)

# Combine DEG tables
comb_full <- full_join(davies, mbm_markers, by = 'gene')
comb_full$sig <- ifelse(comb_full$davies_logFC > 0 & comb_full$mbm_logFC > 0, 'both_up', NA)
comb_full$sig <- ifelse(comb_full$davies_logFC > 0 & (is.na(comb_full$mbm_logFC) == T), 'Davies_up', comb_full$sig)
comb_full$sig <- ifelse(comb_full$mbm_logFC > 0 & (is.na(comb_full$davies_logFC) == T), 'MBM_up', comb_full$sig)
comb_full$sig <- ifelse(comb_full$davies_logFC < 0 & comb_full$mbm_logFC < 0, 'both_down', comb_full$sig)
comb_full$sig <- ifelse(comb_full$davies_logFC < 0 & (is.na(comb_full$mbm_logFC) == T), 'Davies_down', comb_full$sig)
comb_full$sig <- ifelse(comb_full$mbm_logFC < 0 & (is.na(comb_full$davies_logFC) == T), 'MBM_down', comb_full$sig)
table(comb_full$sig, useNA = 'always')
write.csv(comb_full, 'data/Davies/table_combined_topbottom.csv', row.names = F)

# Prepare data for plotting
comb_full$davies_logFC[is.na(comb_full$davies_logFC)] <- 0
comb_full$mbm_logFC[is.na(comb_full$mbm_logFC)] <- 0

comb_up <- subset(comb_full, davies_logFC > 0 & mbm_logFC > 0)
corr_up <- cor.test(comb_up$davies_logFC, comb_up$mbm_logFC, method = 'spearman')

comb_dn <- subset(comb_full, davies_logFC < 0 & mbm_logFC < 0)
corr_dn <- cor.test(comb_dn$davies_logFC, comb_dn$mbm_logFC, method = 'spearman')


# Plots
pdf('data/tumor/plots_scatter_MBPM_sn_tumor_log2FC_100topbottom.pdf')
# MBM
ggplot(comb_up, aes(x = mbm_logFC, y = davies_logFC)) + 
  geom_point() + theme_bw() + 
  ggtitle(paste0('Shared DEG MBM vs MBMbulk\nTop100 genes, Log2FC >0, P <0.05\nRho: ', 
                 round(corr_up$estimate[[1]], 3), ', cor P = ', corr_up$p.value)) + 
  ggrepel::geom_label_repel(aes(label = gene), fill = 'white', size = 2, box.padding = 0.1, 
                            label.padding = 0.1, max.overlaps = 17) + 
  xlab('Log2FC of DEGs in MBM_sn') + ylab('Log2FC of DEGs in MBMbulk')

ggplot(comb_up, aes(x = mbm_logFC, y = davies_logFC)) + 
  geom_point() + theme_bw() + 
  ggtitle(paste0('Shared DEG MBM vs MBMbulk\nTop100 genes, Log2FC >0, P <0.05\nRho: ', 
                              round(corr_up$estimate[[1]], 3), ', cor P = ', corr_up$p.value)) + 
  xlab('Log2FC of DEGs in MBM_sn') + ylab('Log2FC of DEGs in MBMbulk')

# MPM
ggplot(comb_dn, aes(x = mbm_logFC, y = davies_logFC)) + 
  geom_point() + theme_bw() + 
  ggtitle(paste0('Shared DEG MPM vs MPMbulk\nTop100 genes, Log2FC <0, P <0.05\nRho: ', 
                              round(corr_dn$estimate[[1]], 3), ', cor P = ', corr_dn$p.value)) + 
  ggrepel::geom_label_repel(aes(label = gene), fill = 'white', size = 2, box.padding = 0.1, 
                            label.padding = 0.1, max.overlaps = 17) + 
  xlab('Log2FC of DEGs in MPM_sn') + ylab('Log2FC of DEGs in MPMbulk')

ggplot(comb_dn, aes(x = mbm_logFC, y = davies_logFC)) + 
  geom_point() + theme_bw() + 
  ggtitle(paste0('Shared DEG MPM vs MPMbulk\nTop100 genes, Log2FC <0, P <0.05\nRho: ', 
                              round(corr_dn$estimate[[1]], 3), ', cor P = ', corr_dn$p.value)) + 
  xlab('Log2FC of DEGs in MPM_sn') + ylab('Log2FC of DEGs in MPMbulk')
dev.off()


#### Venn diagrams ####
pdf('data/tumor/plots_venn_marker_overlap.pdf')
lis <- list(MBM = mbm_markers_up$gene, Davies = davies$gene[1:100])
ggVennDiagram(lis) + ggtitle('Overlapping genes in MBM sig and MBMbulk sig') + 
  scale_fill_gradient(low = '#F4FAFE', high = '#4981BF')

lis2 <- list(MBM = mbm_markers_dn$gene, Davies = davies$gene[101:200])
ggVennDiagram(lis2) + ggtitle('Overlapping genes in MPM sig and MPMbulk sig') + 
  scale_fill_gradient(low = '#F4FAFE', high = '#4981BF')
dev.off()
