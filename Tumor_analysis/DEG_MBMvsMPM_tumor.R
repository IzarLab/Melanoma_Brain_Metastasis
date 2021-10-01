#!/usr/bin/env Rscript

### title: Differential gene expression (DGE) of tumor cells comparing MBM and MPM
### author: Jana Biermann, PhD

print(Sys.time())

library(dplyr)
library(Seurat)


colBP <- c('#A80D11', '#008DB8')
colSCSN <- c('#E1AC24', '#288F56')

# Read-in integrated dataset
seu <- readRDS('data/MBPM/data_MBPM.rds')

# Subset to sn tumor cells
seu <- subset(seu1, sequencing == 'Single nuclei' & cell_type_main == 'Tumor cells')


#### DGE
Idents(seu) <- seu$organ
markers <- FindMarkers(seu, ident.1 = 'Brain', ident.2 = 'Peripheral', min.pct = 0, logfc.threshold = 0, 
                       assay = 'RNA', max.cells.per.ident = min(table(seu$organ)))
markers$gene <- rownames(markers)
tops <- markers %>% filter(p_val_adj < 0.05) %>% arrange(desc(avg_log2FC)) %>% top_n(n = 20, wt = avg_log2FC)
tops <- rbind.data.frame(tops, markers %>% filter(p_val_adj < 0.05) %>% arrange(avg_log2FC) %>% top_n(n = -20, wt = avg_log2FC))
ifelse(!dir.exists(file.path('data/DEG/organ/tumor')), dir.create(file.path('data/DEG/organ/tumor'), recursive = T), FALSE)
write.csv(markers, 'data/DEG/organ/tumor/markers_MBPM_sn_tumor_organ_full.csv', row.names = F)

# Heatmap
DefaultAssay(seu) <- 'RNA'
seu <- NormalizeData(seu)
seu <- ScaleData(seu, assay = 'RNA')
avg <- AverageExpression(seu, assays = 'RNA', return.seurat = T, group.by = 'organ', slot = 'scale.data')
avg$organ <- rownames(avg@meta.data)

# Plot
pdf('data/DEG/organ/tumor/plots_avg_heatmap_MBPM_sn_tumor_organ.pdf')
DoHeatmap(avg, features = tops$gene, group.by = 'organ', raster = T, assay = 'RNA', draw.lines = F, group.colors = colBP)
dev.off()


print(Sys.time())