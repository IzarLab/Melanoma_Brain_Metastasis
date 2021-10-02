#!/usr/bin/env Rscript

### title: DGE (Differential gene expression) in CD8+ T cells comparing expanded vs non-expanded TCRs
### author: Jana Biermann, PhD

print(Sys.time())

library(dplyr)
library(Seurat)
library(ggplot2)
library(gplots)
library(ggrepel)
library(hypeR)
'%notin%' <- Negate('%in%')


# Read-in integrated object and add TCR info
seu1 <- readRDS('data/MBPM/data_MBPM.rds')
all_tcr <- read.csv('data/TCR/data_all_TCR_v611.csv')
seu1@meta.data <- left_join(seu1@meta.data, all_tcr, by = 'barcode_all')
rownames(seu1@meta.data) <- seu1@meta.data$barcode_all
seu1$both_chains <- ifelse(grepl('TRB:', seu1$cdr3s_aa) == T & grepl('TRA:', seu1$cdr3s_aa) == T, 'both', 'not_both')
seu1$mait <- ifelse(seu1$mait_evidence %in% c('TRA:gene', 'TRA:gene;TRB:gene', 'TRA:gene+junction',
                                              'TRA:gene+junction;TRB:gene', 'TRB:gene'), 'MAIT', 'no_MAIT')
seu1$inkt <- ifelse(seu1$inkt_evidence %in% c('TRB:gene'), 'iNKT', 'no_iNKT')
seu1$public_info <- paste0(seu1$disease_TRB1, ';', seu1$disease_TRB2)


#### Choose one comparison ####

### Compare all expanded CD8+ T cells vs non-expanded CD8+ T cells
type <- 'expanded_withPublicPrivateMaitInkt'

seu <- subset(seu1, sequencing == 'Single cells' & 
                cell_type_int == 'CD8+ T cells' & 
                clone_size %in% c('Non-expanded', 'Expanded') & 
                both_chains == 'both')
seu <- ScaleData(seu)


### Compare expanded CD8+ T cells with public TCRs vs non-expanded CD8+ T cells with public TCRs
type <- 'expanded_withPublic_withoutPrivateMaitInkt'

seu <- subset(seu1, sequencing == 'Single cells' & 
                cell_type_int == 'CD8+ T cells' & 
                clone_size %in% c('Non-expanded', 'Expanded') & 
                both_chains == 'both' & 
                mait == 'no_MAIT' & 
                inkt == 'no_iNKT')

# Exclude private TCRs (NA;NA)
seu <- subset(seu, public_info %notin% c('EBV;NA', 'flu;NA', 'NA;NA'))
seu <- ScaleData(seu)


### Compare expanded CD8+ T cells with private TCRs vs non-expanded CD8+ T cells with private TCRs
type <- 'expanded_withPrivate_withoutPublicMaitInkt'

seu <- subset(seu1, sequencing == 'Single cells' & 
                cell_type_int == 'CD8+ T cells' & 
                clone_size %in% c('Non-expanded', 'Expanded') & 
                both_chains == 'both' & 
                mait == 'no_MAIT' & 
                inkt == 'no_iNKT')

# Exclude public TCRs (labeled as melanoma)
seu <- subset(seu, public_info %notin% c('EBV;NA', 'flu;NA', 'melanoma_mets;NA', 'melanoma_primary;NA', 
                                         'NA;melanoma_mets'))
seu <- ScaleData(seu)


#### Remove TRAV/TRBV genes ####
DefaultAssay(seu) <- 'RNA'
seu <- seu[grep('^TRAV|^TRBV', rownames(seu), value = T, invert = T), ]
seu <- NormalizeData(seu)


#### DEG ####
Idents(seu) <- seu$clone_size
markers <- FindMarkers(seu, ident.1 = 'Expanded', ident.2 = 'Non-expanded', min.pct = 0, logfc.threshold = 0, 
                       assay = 'RNA', max.cells.per.ident = min(table(seu$clone_size)))
markers$gene <- rownames(markers)

tops <- markers %>%
  filter(p_val_adj < 0.05) %>% top_n(n = 10, wt = avg_log2FC)
tops <- rbind.data.frame(tops, markers %>% filter(p_val_adj < 0.05) %>% top_n(n = -10, wt = avg_log2FC))
markers$diffexpressed <- ifelse(markers$avg_log2FC > 0 & markers$p_val_adj < 0.05, 'Up', 
                                ifelse(markers$avg_log2FC < 0 & markers$p_val_adj < 0.05, 'Down', 'NS'))
markers$label <- ifelse(markers$gene %in% tops$gene | markers$gene %in% c('TOX', 'TCF7'), markers$gene, NA)

ifelse(!dir.exists(file.path(paste0('data/DEG/clone_size/', type))), 
       dir.create(file.path(paste0('data/DEG/clone_size/', type)), recursive = T), FALSE)
write.csv(markers, paste0('data/DEG/clone_size/', type, '/markers_MBM_sc_', type, '_clone_size.csv'), row.names = F)

# Plot
pdf(paste0('data/DEG/clone_size/', type, '/plots_volcano_MBM_sc_', type, '_clone_size.pdf'))
ggplot(data = markers, aes(x = avg_log2FC, y = -log10(p_val_adj), col = diffexpressed, label = label)) + 
  geom_point() + geom_text_repel() + scale_color_manual(values = c('darkblue', 'grey', 'darkred')) + 
  labs(color = 'Log2FC\n(Expanded vs\nNon-expanded)') + theme_bw() + coord_cartesian(xlim = c(-3, 3)) + 
  ggtitle(paste0('Expanded vs Non-expanded ', type)) + 
  theme(panel.grid.major = element_blank(), title = element_text(size = 10))

ggplot(data = markers, aes(x = avg_log2FC, y = -log10(p_val_adj), col = diffexpressed, label = label)) + 
  geom_point() + scale_color_manual(values = c('darkblue', 'grey', 'darkred')) + 
  labs(color = 'Log2FC\n(Expanded vs\nNon-expanded)') + 
  theme_bw() + coord_cartesian(xlim = c(-3, 3)) + ggtitle(paste0('Expanded vs Non-expanded ', type)) + 
  theme(panel.grid.major = element_blank(), title = element_text(size = 10))
dev.off()


#### Pathway enrichment ####
HALLMARK <- msigdb_gsets(species='Homo sapiens', category='H')

# DEG
markers <- FindAllMarkers(seu, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25, assay = 'RNA',
                          max.cells.per.ident = min(table(seu$clone_size)))

# Subset markers to 'expanded' and run hypeR
markers <- markers %>% group_by(cluster) %>% dplyr::filter(p_val_adj<0.05)
background <- rownames(seu@assays$RNA@data)
subs <- markers %>% dplyr::filter(cluster=='Expanded') %>% magrittr::use_series(gene)
hyp <- hypeR(as.character(subs), genesets=HALLMARK, background = background)

pdf(paste0('data/DEG/clone_size/', type, '/plots_pathway_MBM_sc_', type, '_clone_size.pdf'))
hyp_dots(hyp, title = 'Hallmarks cluster Expanded', abrv=80, top = 12) +
  theme_bw()+ 
  theme(axis.title.y = element_blank(), axis.text.y=element_text(color='black'),
        axis.text.x=element_text(color='black'), panel.grid.minor = element_blank()) +
  scale_color_gradient(high = 'black', low = '#C51B8A')
dev.off()


print(Sys.time())