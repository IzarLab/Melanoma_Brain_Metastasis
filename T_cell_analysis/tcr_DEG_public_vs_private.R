#!/usr/bin/env Rscript

### title: DGE (Differential gene expression) in expanded CD8+ T cells comparing public vs private TCRs 
### author: Jana Biermann, PhD

print(Sys.time())

library(dplyr)
library(Seurat)
library(ggplot2)
library(gplots)
library(DropletUtils)
library(ggrepel)
'%notin%' <- Negate('%in%')

# Read-in integrated object, add TCR info and subset
seu1 <- readRDS('data/MBPM/data_MBPM.rds')
all_tcr <- read.csv('data/TCR/data_all_TCR_v611.csv')
seu1@meta.data <- left_join(seu1@meta.data, all_tcr, by = 'barcode_all')
rownames(seu1@meta.data) <- seu1@meta.data$barcode_all
seu1$clone_size <- ifelse(seu1$frequency < 2, 'Non-expanded', 'Expanded')
seu <- subset(seu1, sequencing == 'Single cells' & cell_type_int == 'CD8+ T cells' & 
                clone_size %in% c('Expanded'))
seu$both_chains <- ifelse(grepl('TRB:', seu$cdr3s_aa) == T & 
                            grepl('TRA:', seu$cdr3s_aa) == T, 'both', 'not_both')
seu$mait <- ifelse(seu$mait_evidence %in% c('TRA:gene', 'TRA:gene;TRB:gene', 'TRA:gene+junction', 
                                            'TRA:gene+junction;TRB:gene', 'TRB:gene'), 'MAIT', 'no_MAIT')
seu$inkt <- ifelse(seu$inkt_evidence %in% c('TRB:gene'), 'iNKT', 'no_iNKT')
seu$public_info <- paste0(seu$disease_TRB1, ';', seu$disease_TRB2)
seu <- subset(seu, both_chains == 'both' & public_info %notin% c('EBV;NA', 'flu;NA') & 
                mait == 'no_MAIT' & inkt == 'no_iNKT')
seu$public <- ifelse(seu$public_info %in% c('melanoma_mets;NA', 'melanoma_primary;NA', 
                                            'NA;melanoma_mets'), 'public_melanoma', 'private')
seu <- ScaleData(seu)


#### Remove TRAV/TRBV genes ####
DefaultAssay(seu) <- 'RNA'
seu <- seu[grep('^TRAV|^TRBV', rownames(seu), value = T, invert = T), ]
seu <- NormalizeData(seu)


#### DEG ####
Idents(seu) <- seu$public
markers <- FindMarkers(seu, ident.1 = 'private', ident.2 = 'public_melanoma', min.pct = 0, logfc.threshold = 0, 
                       assay = 'RNA', max.cells.per.ident = min(table(seu$public)))
markers$gene <- rownames(markers)

tops <- markers %>% filter(p_val_adj < 0.05) %>% top_n(n = 10, wt = avg_log2FC)
tops <- rbind.data.frame(tops, markers %>% filter(p_val_adj < 0.05) %>% top_n(n = -10, wt = avg_log2FC))
markers$diffexpressed <- ifelse(markers$avg_log2FC > 0 & markers$p_val_adj < 0.05, 'Up', 
                                ifelse(markers$avg_log2FC < 0 & markers$p_val_adj < 0.05, 'Down', 'NS'))
markers$label <- ifelse(markers$gene %in% tops$gene, markers$gene, NA)

ifelse(!dir.exists(file.path('data/DEG/public_TCRs/')), dir.create(file.path('data/DEG/public_TCRs/'), recursive = T), FALSE)
write.csv(markers, 'data/DEG/public_TCRs/markers_MBM_sc_public.csv', row.names = F)

pdf('data/DEG/public_TCRs/plots_volcano_MBM_sc_public.pdf')
ggplot(data = markers, aes(x = avg_log2FC, y = -log10(p_val_adj), col = diffexpressed, label = label)) + 
  geom_point() + geom_text_repel(size = 3) + 
  scale_color_manual(values = c('darkblue', 'grey', 'darkred')) + 
  labs(color = 'Log2FC\n(Private vs\nPublic melanoma)') + 
  theme_bw() + coord_cartesian(xlim = c(-3, 3)) + ggtitle('Private vs public') + 
  theme(panel.grid.major = element_blank(), title = element_text(size = 10))

ggplot(data = markers, aes(x = avg_log2FC, y = -log10(p_val_adj), col = diffexpressed, label = label)) + 
  geom_point() + scale_color_manual(values = c('darkblue', 'grey', 'darkred')) + 
  labs(color = 'Log2FC\n(Private vs\nPublic melanoma)') + 
  theme_bw() + coord_cartesian(xlim = c(-3, 3)) + ggtitle('Private vs public') + 
  theme(panel.grid.major = element_blank(), title = element_text(size = 10))
dev.off()


#### Distribution of TCRs ####
seu2 <- subset(seu1, sequencing == 'Single cells' & cell_type_int == 'CD8+ T cells')
seu2$both_chains <- ifelse(grepl('TRB:', seu2$cdr3s_aa) == T & 
                             grepl('TRA:', seu2$cdr3s_aa) == T, 'both', 'not_both')
seu2$mait <- ifelse(seu2$mait_evidence %in% c('TRA:gene', 'TRA:gene;TRB:gene', 'TRA:gene+junction', 
                                              'TRA:gene+junction;TRB:gene', 'TRB:gene'), 'MAIT', 'no_MAIT')
seu2$inkt <- ifelse(seu2$inkt_evidence %in% c('TRB:gene'), 'iNKT', 'no_iNKT')
seu2 <- subset(seu2, both_chains == 'both' & mait == 'no_MAIT' & inkt == 'no_iNKT')
seu2$public_info <- paste0(seu2$disease_TRB1, ';', seu2$disease_TRB2)
seu2$public <- ifelse(seu2$public_info == 'EBV;NA', 'EBV', 
                      ifelse(seu2$public_info %in% c('melanoma_mets;NA', 'melanoma_primary;NA', 'NA;melanoma_mets'), 
                             'public_melanoma', 'private'))
seu2 <- subset(seu2, public != 'EBV')
seu2$clone_pub <- paste0(seu2$clone_size, '_', seu2$public)
df <- seu2@meta.data %>%
  select('clone_size', 'public', 'patient', 'clone_pub') %>%
  group_by(patient, clone_pub) %>%
  tally()

pdf('data/DEG/public_TCRs/plots_MBM_sc_TCR_distribution.pdf')
ggplot(data = df, aes(x = patient, y = n, fill = clone_pub)) + 
  geom_bar(stat = 'identity', position = position_dodge(), color = 'black', size = 0.3) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  ggtitle('Distribution of TCRs per patient') + xlab('Patient') + ylab('# Cells')
dev.off()


print(Sys.time())