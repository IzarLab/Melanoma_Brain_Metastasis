#!/usr/bin/env Rscript

### title: Overlap of MBM enriched genes from cell line RNA-seq, ATAC-seq, TF and motif enrichment
### author: Jana Biermann, PhD

library(dplyr)
library(ggplot2)
library(gplots)
library(viridis)
library(scales)
library(ggrastr)
library(ggrepel)
library(ggpubr)
library(patchwork)
library(ggrastr)
library(ggvenn)
'%notin%' <- Negate('%in%')

colBP <- c('#A80D11', '#008DB8')
directory<-'data/celllines_ATACseq/'


#### TOBIAS plot ####
df <- read.table('Melanoma_ATACseq Analysis/All brain vs All peripheral Analysis/Footprinting Analysis/bindetect_results.txt',header = T)
df$sig<-ifelse(df$Brain_Peripheral_change>0.2 & df$Brain_Peripheral_pvalue<0.05,'Brain',
                 ifelse(df$Brain_Peripheral_change< -0.2 & df$Brain_Peripheral_pvalue<0.05,'Peripheral',NA))
table(df$sig,useNA = 'always')
df$gene<-ifelse(is.na(df$sig)==F,df$name,NA)


pdf(paste0(directory,'TF_footprints_all.pdf'))
ggplot(df, aes(x =Brain_Peripheral_change, y = -log10(Brain_Peripheral_pvalue),col=sig)) + 
  geom_point() + theme_bw() + 
  ggtitle(paste0('BINDetect all BvsP ungrouped')) + 
  ggrepel::geom_text_repel(aes(label = gene),col='black',size=2.5,max.overlaps =40) + 
  xlab('Differential binding score') + ylab('-log10(pvalue)')+
  scale_color_manual(values = c(colBP[1],colBP[2]))+
  xlim(-0.6,0.6)

ggplot(df, aes(x =Brain_Peripheral_change, y = -log10(Brain_Peripheral_pvalue),col=sig)) + 
  geom_point() + theme_bw() + 
  ggtitle(paste0('BINDetect all BvsP ungrouped')) + 
  xlab('Differential binding score') + ylab('-log10(pvalue)')+
  scale_color_manual(values = c(colBP[1],colBP[2]))+
  xlim(-0.6,0.6)
dev.off()


#### Find overlap TOBIAS, motif enrichment, SCENIC & VIPER ####
## TOBIAS
tobias <- read.table('Melanoma_ATACseq Analysis/All brain vs All peripheral Analysis/Footprinting Analysis/bindetect_results.txt',header = T)
tobias$sig<-ifelse(tobias$Brain_Peripheral_change>0 & 
                     tobias$Brain_Peripheral_pvalue<0.05,'Brain',
               ifelse(tobias$Brain_Peripheral_change< 0 & 
                        tobias$Brain_Peripheral_pvalue<0.05,'Peripheral',NA))
tobias <- tobias %>% 
  arrange(desc(Brain_Peripheral_change)) %>% 
  filter(Brain_Peripheral_pvalue<0.05)
table(tobias$sig,useNA = 'always')
intersect(tobias[tobias$sig=='Brain',2],tobias[tobias$sig=='Peripheral',2])
tobias$name[which(duplicated(tobias$name))]

# Remove duplicates
tob<-aggregate(Brain_Peripheral_change ~ name, mean, data = tobias)
tob<-data.frame(gene=tob$name,
                tobias_diffBind=tob$Brain_Peripheral_change)
tob <- tob %>% arrange(desc(tobias_diffBind))
tob$gene<-toupper(tob$gene)

## Motif enrichment
# Motif MBM
motif_up <- read.delim('\Melanoma_ATACseq Analysis/All brain vs All peripheral Analysis/Motif Analysis/BrainCL_over_PeriCL.deseq.Padj0.01.LG2FC.1.up.bed_motif/knownResults.txt',header = T)
motif_up <- motif_up %>% arrange(P.value) %>% filter(P.value<0.05)

# Motif ECM
motif_dn <- read.delim('Melanoma_ATACseq Analysis/All brain vs All peripheral Analysis/Motif Analysis/BrainCL_over_PeriCL.deseq.Padj0.01.LG2FC.-1.down.bed_motif/knownResults.txt',header = T)
motif_dn <- motif_dn %>% arrange(P.value) %>% filter(P.value<0.05)

# Combine motif enrichments and reformat
both<-intersect(motif_up$Motif.Name,motif_dn$Motif.Name)
motif_up <- motif_up[motif_up$Motif.Name %notin% both,]
motif_dn <- motif_dn[motif_dn$Motif.Name %notin% both,]

motif_up$gene<-toupper(lapply(motif_up$Motif.Name, 
                              function(x){strsplit(x,"(",fixed = T)[[1]][1]}))
motif_up$direction<-1
motif_up<-motif_up[,-c(2,4:9)]
motif_dn$gene<-toupper(lapply(motif_dn$Motif.Name, 
                              function(x){strsplit(x,"(",fixed = T)[[1]][1]}))
motif_dn$direction<- -1
motif_dn<-motif_dn[,-c(2,4:9)]
motif<-rbind.data.frame(motif_up,motif_dn)
colnames(motif)[2]<-'motif_pval'
colnames(motif)[4]<-'motif_direction'

## SCENIC
markers<-read.csv('data/tumor/scenic/markers_MBPM_sn_scenic.csv')
scenic<-data.frame(gene=markers$gene,
                   scenic_logFC=markers$avg_log2FC,
                   scenic_adjp=markers$p_val_adj)
scenic$gene<-gsub('(+)','',scenic$gene,fixed = T)
scenic <- scenic %>% arrange(desc(scenic_logFC)) %>% filter(scenic_adjp<0.05)

## VIPER
viper_res<-read.csv('data/viper/MBP_vs_MPM_MRs_t_test.csv')
viper<-data.frame(gene=viper_res$genes,
                  viper_logFC=viper_res$avg_log2FC,
                  viper_adjp=viper_res$p.value)
viper <- viper %>% arrange(desc(viper_logFC)) %>% filter(viper_adjp<0.05)


## Combine all tables
comb_full <- full_join(viper, scenic, by = 'gene')
intersect(toupper(comb_full$gene),toupper(motif$gene))
comb_full<-full_join(comb_full,motif,by='gene')
comb_full<-full_join(comb_full,tob,by='gene')
table(comb_full$sig, useNA = 'always')
comb_up <- subset(comb_full, viper_logFC > 0 & scenic_logFC > 0 & 
                    motif_direction>0 & tobias_diffBind >0)
comb_dn <- subset(comb_full, viper_logFC < 0 & scenic_logFC < 0 & 
                    motif_direction< 0 & tobias_diffBind < 0)


# Prepare venn diagrams
scenic_up<-scenic %>% filter(scenic_logFC>0) %>% select(gene)
viper_up<-viper %>% filter(viper_logFC>0)%>% select(gene)
motif_up<-motif %>% filter(motif_direction>0)%>% select(gene)
tobias_up<-tob %>% filter(tobias_diffBind>0)%>% select(gene)
intersect(scenic_up$gene,
          intersect(viper_up$gene,
                    intersect(motif_up$gene,tob$gene)))

scenic_dn<-scenic %>% filter(scenic_logFC< 0) %>% select(gene)
viper_dn<-viper %>% filter(viper_logFC< 0) %>% select(gene)
motif_dn<-motif %>% filter(motif_direction< 0)%>% select(gene)
tobias_dn<-tob %>% filter(tobias_diffBind< 0)%>% select(gene)
intersect(scenic_dn$gene,
          intersect(viper_dn$gene,
                    intersect(motif_dn$gene,tob$gene)))

lis_all <- list(
  scenic_MBM = scenic_up$gene, 
  viper_MBM = viper_up$gene,
  motif_MBM = motif_up$gene,
  tobias_MBM = tobias_up$gene,
  scenic_ECM = scenic_dn$gene, 
  viper_ECM = viper_dn$gene,
  motif_ECM = motif_dn$gene,
  tobias_ECM = tobias_dn$gene)


## Plots
pdf(paste0(directory,'plots_motif_tobias_scenic_viper_comparison.pdf'))
# MBM
ggplot(comb_up, aes(x =viper_logFC, y = scenic_logFC,col=tobias_diffBind)) + 
  geom_point() + theme_bw()+scale_color_gradient()+
  ggrepel::geom_text_repel(aes(label = gene),size=2.5,col='black')+
  ggtitle('MBM (Delta>0 & P<0.05 in scenic, viper, motif, tobias)')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# ECM
ggplot(comb_dn, aes(x =viper_logFC, y = scenic_logFC,col=tobias_diffBind)) + 
  geom_point() + theme_bw()+scale_color_gradient()+
  ggrepel::geom_text_repel(aes(label = gene),size=2.5,col='black')+
  ggtitle('ECM (Delta<0 & P<0.05 in scenic, viper, motif, tobias)')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# MBM
ggvenn(lis_all[1:4],fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF",'forestgreen'))+
  ggtitle('Overlap of differential TFs/proteins in MBM\nusing motifs, TOBIAS, SCENIC and VIPER\n(Adj P < 0.05; log2FC > 0)\n')

# ECM
ggvenn(lis_all[5:8],fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF",'forestgreen'))+
  ggtitle('Overlap of differential TFs/proteins in ECM\nusing motifs, TOBIAS, SCENIC and VIPER\n(Adj P < 0.05; log2FC < 0)\n')

# MBM
ggvenn(lis_all[c(3,4)],fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF",'forestgreen'))+
  ggtitle('Overlap of differential TFs in MBM\nusing motifs and TOBIAS\n(Adj P < 0.05; log2FC > 0)\n')

# ECM
ggvenn(lis_all[c(7,8)],fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF",'forestgreen'))+
  ggtitle('Overlap of differential TFs in ECM\nusing motifs and TOBIAS\n(Adj P < 0.05; log2FC < 0)\n')
dev.off()

