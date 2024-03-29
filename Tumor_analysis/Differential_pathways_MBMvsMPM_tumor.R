#!/usr/bin/env Rscript

### title: Differential pathways compared between MBM and MPM tumor cells 
### author: Jana Biermann, PhD

library(hypeR)
library(dplyr)
library(Seurat)
library(ggplot2)
library(gplots)
library(DropletUtils)
library(msigdbr)
library(scales)
library(viridis)
library(plyr)
library(ggrastr)
library(patchwork)
colBP <- c('#A80D11', '#008DB8')
colSCSN <- c('#E1AC24', '#288F56')


# Subset object to sn tumor cells
seu <- readRDS('data/MBPM/data_MBPM_scn.rds')
seu <- subset(seu, sequencing == 'Single nuclei' & cell_type_fine == 'Tumor cells')

# DEG
Idents(seu) <- seu$organ
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, 
                          assay = 'RNA', max.cells.per.ident = min(table(seu$organ)),
                          test.use = 'MAST')
ifelse(!dir.exists(file.path('data/DEG/organ/tumor')), 
       dir.create(file.path('data/DEG/organ/tumor'), recursive = T), FALSE)
write.csv(markers, 'data/DEG/organ/tumor/markers_MBPM_sn_tumor_organ.csv', row.names = F)


#### Get significant pathways for marker sets ####
# Read-in genesets from hypeR
HALLMARK <- msigdb_gsets(species = 'Homo sapiens', category = 'H')
KEGG <- msigdb_gsets(species = 'Homo sapiens', category = 'C2', subcategory = 'CP:KEGG')
REACTOME <- msigdb_gsets(species = 'Homo sapiens', category = 'C2', subcategory = 'CP:REACTOME')
PID <- msigdb_gsets(species = 'Homo sapiens', category = 'C2', subcategory = 'CP:PID')
CP <- msigdb_gsets(species = 'Homo sapiens', category = 'C2', subcategory = 'CP')
C5_BP <- msigdb_gsets(species = 'Homo sapiens', category = 'C5', subcategory = 'BP')
C5_MF <- msigdb_gsets(species = 'Homo sapiens', category = 'C5', subcategory = 'MF')
C6 <- msigdb_gsets(species = 'Homo sapiens', category = 'C6')
C7 <- msigdb_gsets(species = 'Homo sapiens', category = 'C7', subcategory = 'IMMUNESIGDB')
sigs_wiki <- read.csv('signatures/msigdb_wiki.csv', na.strings = c('', NA))
genesets_wiki <- as.list(as.data.frame(sigs_wiki))

background <- rownames(seu@assays$RNA@data)
markers <- markers %>% group_by(cluster) %>% dplyr::filter(p_val_adj < 0.05)

# Loop through pathways and collect significant ones
sig.pw <- c()
for (pw in c(C5_BP, C5_MF, CP, HALLMARK, KEGG, REACTOME, C6, C7, PID)) {
  for (i in unique(markers$cluster)) {
    subs <- markers %>% dplyr::filter(cluster == i) %>% magrittr::use_series(gene)
    hyp <- hypeR(as.character(subs), genesets = pw, background = background)
    sig.pw <- c(sig.pw, hyp$data %>% filter(fdr < 1e-04) %>% select(label))
  }
}
for (i in unique(markers$cluster)) {
  subs <- markers %>% dplyr::filter(cluster == i) %>% magrittr::use_series(gene)
  hyp <- hypeR(as.character(subs), genesets = genesets_wiki, background = background)
  sig.pw <- c(sig.pw, hyp$data %>% filter(fdr < 1e-04) %>% select(label))
}
sig.pw <- unlist(sig.pw) %>% unique
length(sig.pw)


#### Apply significant pathways as module scores ####
msigdbgenes <- read.csv('~/brain_mets/signatures/relevant_msigdb_pathways_hypeR.csv', na.strings = '')

for (s in sig.pw) {
  seu <- AddModuleScore(object = seu, features = list(msigdbgenes[, s]), 
                        name = s, assay = 'RNA', search = F)
}


### Significance permutation using GmAMisc package (modified version)
source('Tumor_analysis/perm.t.test.mod.R')

# Set up result data frame
perm.stats <- data.frame(pathway = NA, sample1.lab = NA, n1 = NA, sample1_lci = NA, mean1 = NA, sample1_uci = NA, 
                         sample2.lab = NA, n2 = NA, sample2_lci = NA, mean2 = NA, sample2_uci = NA, meanDiff.1. = NA, p.lowertail = NA, 
                         p.uppertail = NA, two.sided.p = NA, p.equal.var = NA, p.unequal.var = NA, p.wilcox = NA)

# Loop through significant pathways, remove outliers and permute
for (s in sig.pw) {
  q10 <- quantile(seu@meta.data[, paste0(s, '1')], probs = 0.1)
  q90 <- quantile(seu@meta.data[, paste0(s, '1')], probs = 0.9)
  tmp <- data.frame(seu@meta.data[, paste0(s, '1')], seu$organ)
  tmp <- subset(tmp, tmp[, 1] > q10 & tmp[, 1] < q90)
  
  perm.t.test.mod(tmp, format = 'long', B = 10000, sample1.lab = 'Brain', sample2.lab = 'Peripheral', 
                  pathway = s, plot = F)
  perm.stats <- rbind.data.frame(perm.stats, res_perm)
}
perm.stats <- perm.stats[-1, ]
write.csv(perm.stats, 'data/DEG/organ/downsampled/tumor/permpathways_MBPM_sn_tumor_organ_down.csv', row.names = F)


#### Plot violins ####
final<-c('KEGG_OXIDATIVE_PHOSPHORYLATION', 'WP_PI3KAKT_SIGNALING_PATHWAY',
         'KEGG_INSULIN_SIGNALING_PATHWAY', 'KEGG_ERBB_SIGNALING_PATHWAY', 
         'WP_EGFR_TYROSINE_KINASE_INHIBITOR_RESISTANCE', 'PID_PDGFRB_PATHWAY',
         'PID_KIT_PATHWAY', 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
         'HALLMARK_MTORC1_SIGNALING', 'KEGG_CELL_ADHESION_MOLECULES_CAMS',
         'KEGG_MELANOGENESIS')

msigdbgenes <- read.csv('~/brain_mets/signatures/relevant_msigdb_pathways_hypeR.csv', na.strings = '')
for (s in final) {
  seu <- AddModuleScore(object = seu, features = list(msigdbgenes[, s]), 
                        name = s, assay = 'RNA', search = F)
}


pdf('data/DEG/organ/tumor/plots_vln_MBPM_sn_tumor_organ_final.pdf',width = 10,height = 10)
plots <- VlnPlot(object = seu, features = paste0(final, 1), pt.size = 0, group.by = 'sequencing', 
                 split.by = 'organ', split.plot = T, ncol = 3, combine = F, cols = colBP)
plots <- lapply(X = plots, FUN = function(p) p + geom_boxplot(width = 0.1, outlier.color = NA) + 
                  theme(plot.title = element_text(size = 4), 
                        axis.title.x = element_blank(), 
                        axis.text.x = element_blank(), 
                        legend.position = 'none', 
                        axis.ticks.x = element_blank()))
print(CombinePlots(plots = plots, ncol = 4, legend = 'bottom'))

plots2 <- VlnPlot(object = seu, features = paste0(final, 1), pt.size = 0, group.by = 'patient', 
                  combine = F)
plots2 <- lapply(X = plots2, FUN = function(p) p +
                  theme(plot.title = element_text(size = 4), axis.title.x = element_blank(), 
                        axis.text.x = element_blank(), 
                        legend.position = 'none', 
                        axis.ticks.x = element_blank(),
                        legend.key.size = unit(0.2,'cm'),
                        legend.text = element_text(size=5)))
print(CombinePlots(plots = plots2, ncol = 3, legend = 'bottom'))
dev.off()

