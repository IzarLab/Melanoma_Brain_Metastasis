#!/usr/bin/env Rscript

### title: Differential gene expression (DGE) in TOX+ CD8+ T cells
### comparing single-cell (sc) and single-nuclei (sn) RNA-seq 
### authors: Jana Biermann, PhD; Yiping Wang, PhD

library(dplyr)
library(Seurat)
library(ggplot2)
library(gplots)
library(ggrepel)
library(DropletUtils)
library(scales)
'%notin%' <- Negate('%in%')

colBP <- c('#A80D11', '#008DB8')
colSCSN <- c('#E1AC24', '#288F56')

# Read-in integrated object and subset
seu <- readRDS('data/MBPM/data_MBPM_scn.rds')
seu <- subset(seu, cell_type_fine == 'CD8+ T cells TOX+')

# Remove RPL/S and MT genes
DefaultAssay(seu) <- 'RNA'
seu <- seu[grep('^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA|^MT-', rownames(seu), 
                value = T, invert = T), ]

# Downsampling to adjust for differences in quality between the sequencing techniques
downsampled <- downsampleBatches(as.matrix(seu@assays$RNA@counts), batch = seu$sequencing)
seu[['downsampled']] <- CreateAssayObject(counts = downsampled)
DefaultAssay(seu) <- 'downsampled'
seu <- NormalizeData(seu)

# DEG
Idents(seu) <- seu$sequencing
markers <- FindMarkers(seu, ident.1 = 'Single cells', ident.2 = 'Single nuclei', 
                       assay = 'downsampled', test.use = 'MAST',
                       logfc.threshold = 0, min.pct = 0, 
                       max.cells.per.ident = min(table(seu$sequencing)))
markers$gene <- rownames(markers)

# Label different log2FC categories
markers$diffexpressed <- ifelse(markers$p_val_adj < 0.05 & markers$avg_log2FC > 2, 
                                'Upregulated', 
                                ifelse(markers$p_val_adj < 0.05 & markers$avg_log2FC < (-2), 
                                       'Downregulated', 'NS'))

# Add top 20 DEGs as labels and save file
tops <- markers %>% 
  filter(p_val_adj < 0.05) %>% 
  top_n(n = 20, wt = avg_log2FC)
tops <- rbind.data.frame(tops, markers %>% 
                           filter(p_val_adj < 0.05) %>% 
                           top_n(n = -20, wt = avg_log2FC))
markers$label <- ifelse(markers$gene %in% tops$gene, markers$gene, NA)
ifelse(!dir.exists(file.path('data/DEG/sequencing/downsampled_noRPMT/cd8_tox/')), 
       dir.create(file.path('data/DEG/sequencing/downsampled_noRPMT/cd8_tox/'), recursive = T), FALSE)
write.csv(markers, 'data/DEG/sequencing/downsampled_noRPMT/cd8_tox/markers_MBPM_cd8_tox_down_full.csv', row.names = F)

# Average expression for sequencing groups
sn_avglogexp <- log2(rowMeans(expm1(seu@assays$downsampled@data[, seu$sequencing == 'Single nuclei'])) + 1)
sc_avglogexp <- log2(rowMeans(expm1(seu@assays$downsampled@data[, seu$sequencing == 'Single cells'])) + 1)
all(names(sc_avglogexp) == names(sn_avglogexp))
avlogexp <- cbind.data.frame(sn_avglogexp, sc_avglogexp)
avlogexp$gene <- rownames(avlogexp)
avlogexp <- left_join(avlogexp, markers, by = 'gene')

# Plot
pdf('data/DEG/sequencing/downsampled_noRPMT/cd8_tox/plots_heatmap_MBPM_cd8_tox_down.pdf')
ggplot(avlogexp, aes(x = sc_avglogexp, y = sn_avglogexp)) + 
  geom_point(aes(color = diffexpressed), size = 0.2) + 
  scale_color_manual(values = c(colSCSN[2], 'grey', colSCSN[1])) + 
  theme_minimal() + 
  theme(legend.position = 'right') + 
  ggrepel::geom_text_repel(aes(label = label), size = 2, max.overlaps = 20) + 
  labs(title = paste0('TOX+ CD8+ T cells (R^2: ', 
                      round(cor.test(avlogexp$sc_avglogexp, avlogexp$sn_avglogexp, 
                                     method = 'pearson')$estimate, 2), ')')) + 
  xlab('Log2 of average normalized expression in single cells') + 
  ylab('Log2 of average normalized expression in single nuclei') + 
  guides(color = guide_legend(title = 'Differential expression')) + 
  geom_abline(intercept = 0, slope = 1, col = 'black', size = 0.2)

ggplot(avlogexp, aes(x = sc_avglogexp, y = sn_avglogexp)) + 
  geom_point(aes(color = diffexpressed), size = 0.2) + 
  scale_color_manual(values = c(colSCSN[2], 'grey', colSCSN[1])) + 
  theme_minimal() + 
  theme(legend.position = 'right') + 
  labs(title = paste0('TOX+ CD8+ T cells (R^2: ', 
                      round(cor.test(avlogexp$sc_avglogexp, avlogexp$sn_avglogexp, 
                                     method = 'pearson')$estimate, 2), ')')) + 
  xlab('Log2 of average normalized expression in single cells') + 
  ylab('Log2 of average normalized expression in single nuclei') + 
  guides(color = guide_legend(title = 'Differential expression')) + 
  geom_abline(intercept = 0, slope = 1, col = 'black', size = 0.2)
dev.off()

