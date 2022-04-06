#!/usr/bin/env Rscript

### title: DGE (Differential gene expression) in CD8+ T cells comparing expanded vs non-expanded TCRs
### author: Jana Biermann, PhD

library(dplyr)
library(Seurat)
library(ggplot2)
library(gplots)
library(ggrepel)
library(hypeR)
library(DropletUtils)
'%notin%' <- Negate('%in%')

colExp<-c('#C82F71','#7473A6')


### Compare all expanded CD8+ T cells vs non-expanded CD8+ T cells
type <- 'expanded_vs_nonexpanded'

# Read-in integrated object
seu <- readRDS('data/MBPM/data_MBPM_scn.rds')
seu <- subset(seu, sequencing == 'Single cells' & 
                cell_type_int == 'CD8+ T cells' & 
                clone_size %in% c('Non-expanded', 'Expanded') & 
                both_chains == 'both')
seu <- ScaleData(seu)


#### Remove TRAV/TRBV genes ####
DefaultAssay(seu) <- 'RNA'
seu <- seu[grep('^TRAV|^TRBV', rownames(seu), value = T, invert = T), ]
seu <- NormalizeData(seu)


#### DEG ####
Idents(seu) <- seu$clone_size
markers <- FindMarkers(seu, ident.1 = 'Expanded', ident.2 = 'Non-expanded', 
                       min.pct = 0, logfc.threshold = 0,test.use = 'MAST',
                       assay = 'RNA', max.cells.per.ident = min(table(seu$clone_size)))
markers$gene <- rownames(markers)

tops <- markers %>%
  filter(p_val_adj < 0.05) %>% 
  top_n(n = 10, wt = avg_log2FC)
tops <- rbind.data.frame(tops, markers %>% 
                           filter(p_val_adj < 0.05) %>% 
                           top_n(n = -10, wt = avg_log2FC))
markers$diffexpressed <- ifelse(markers$avg_log2FC > 0 & markers$p_val_adj < 0.05, 'Up', 
                                ifelse(markers$avg_log2FC < 0 & markers$p_val_adj < 0.05, 
                                       'Down', 'NS'))
markers$label <- ifelse(markers$gene %in% tops$gene | markers$gene %in% c('TOX', 'TCF7'), 
                        markers$gene, NA)

ifelse(!dir.exists(file.path(paste0('data/DEG/clone_size/', type))), 
       dir.create(file.path(paste0('data/DEG/clone_size/', type)), recursive = T), FALSE)
write.csv(markers, paste0('data/DEG/clone_size/', type, '/markers_MBM_sc_', type, '_clone_size.csv'), row.names = F)


# Plot
pdf(paste0('data/DEG/clone_size/', type, '/plots_volcano_MBM_sc_', type, '_clone_size.pdf'))
ggplot(data = markers, aes(x = avg_log2FC, y = -log10(p_val_adj), col = diffexpressed, label = label)) + 
  geom_point() + 
  geom_text_repel() + 
  scale_color_manual(values = c(colExp[2], 'grey', colExp[1])) + 
  labs(color = 'Log2FC\n(Expanded vs\nNon-expanded)') + 
  theme_bw() + 
  coord_cartesian(xlim = c(-3, 3)) + 
  ggtitle(paste0('Expanded vs Non-expanded ', type)) + 
  theme(panel.grid.major = element_blank(), title = element_text(size = 10))

ggplot(data = markers, aes(x = avg_log2FC, y = -log10(p_val_adj), col = diffexpressed, label = label)) + 
  geom_point() + 
  scale_color_manual(values = c(colExp[2], 'grey', colExp[1])) + 
  labs(color = 'Log2FC\n(Expanded vs\nNon-expanded)') + 
  theme_bw() + 
  coord_cartesian(xlim = c(-3, 3)) + 
  ggtitle(paste0('Expanded vs Non-expanded ', type)) + 
  theme(panel.grid.major = element_blank(), title = element_text(size = 10))
dev.off()
