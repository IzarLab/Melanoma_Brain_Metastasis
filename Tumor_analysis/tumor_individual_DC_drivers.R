#!/usr/bin/env Rscript

### title: Individual diffusion component (DC) analysis of tumor cells
### and overlap of drivers in DC1-3 between samples 
### author: Jana Biermann, PhD

print(Sys.time())

library(dplyr)
library(Seurat)
library(ggplot2)
library(gplots)
library(reshape2)
library(viridis)
library(destiny)
library(ggthemes)
library(colorRamps)
library(hypeR)

colBP <- c('#A80D11', '#008DB8')
colSCSN <- c('#E1AC24', '#288F56')


# Read-in integrated object and subset to non-cycling tumor cells
seu <- readRDS('data/MBPM/data_MBPM.rds')
seu <- subset(seu, cell_type_main == 'Tumor cells' & cell_cycle == 'non-cycling')
DefaultAssay(seu) <- 'RNA'
seu <- DietSeurat(seu, assays = 'RNA', counts = T)

# Split object and remove samples with <500 cells
obj.list <- SplitObject(seu, split.by = 'patient')
for (i in 1:length(obj.list)) {
  if (dim(obj.list[[i]]@assays$RNA@counts)[2] < 500) {
    obj.list[i] <- NA
  }
}
# Remove NAs
obj.list <- obj.list[!is.na(obj.list)]


# Loop through samples and perform DC analysis
for (seu in obj.list) {
  pat <- seu$patient[1]
  seu <- NormalizeData(seu) %>% ScaleData()
  
  # Run DC analysis
  es <- as.ExpressionSet(as.data.frame(t(seu@assays$RNA@data)))
  es@phenoData@data <- seu@meta.data
  dm <- DiffusionMap(es, verbose = T, n_pcs = 20)
  ifelse(!dir.exists(file.path(paste0('data/tumor/diffusion/single/', pat))), 
         dir.create(file.path(paste0('data/tumor/diffusion/single/', pat)), recursive = T), FALSE)
  saveRDS(dm, paste0('data/tumor/diffusion/single/', pat, '/data_dm_', pat, '.rds'))
  
  # Add DC info to Seurat and save loadings table
  seu <- AddMetaData(seu, dm@eigenvectors, col.name = colnames(dm@eigenvectors))
  seu[['DC']] <- CreateDimReducObject(embeddings = dm@eigenvectors, key = 'DC_', assay = 'RNA')
  seu <- ProjectDim(seu, reduction = 'DC')
  write.csv(as.data.frame(seu@reductions$DC@feature.loadings.projected), 
            paste0('data/tumor/diffusion/single/', pat, '/table_dm_loadings_', pat, '.csv'))
}

# Combine genes from individual patients
res <- matrix(data = NA, nrow = 150, ncol = length(names(obj.list)))
colnames(res) <- names(obj.list)
for (pat in names(obj.list)) {
  loading <- read.csv(paste0('data/tumor/diffusion/single/', pat, '/table_dm_loadings_', pat, '.csv'))
  markers_pat <- c()
  for (DC in c('DC_1', 'DC_2', 'DC_3')) {
    subs <- as.data.frame(abs(loading[, DC]))
    subs$genes <- loading[, 'X']
    subs <- subs %>% arrange(desc(`abs(loading[, DC])`)) %>% dplyr::slice_head(n = 50) %>% select(genes)
    markers_pat <- c(markers_pat, subs)
  }
  markers_pat <- unique(unlist(markers_pat))
  res[1:length(markers_pat), pat] <- markers_pat
}
write.csv(res, 'data/tumor/diffusion/table_single_DC1-3_drivers.csv', row.names = F)


# Plot MBM bar plots and pathway enrichments with genes detected in >25% of samples
res_mbm <- res[, grep('MBM', colnames(res))]
df <- as.data.frame(sort(table(res_mbm)/ncol(res_mbm), decreasing = F))
df_plot <- df[df$Freq > 0.25, ]
pdf('data/tumor/diffusion/single/plot_MBM_single_DC1-3_drivers_barplot.pdf', height = 14)
ggplot(df_plot, aes(x = Freq, y = res_mbm)) + ggtitle('MBM - Top overlapping drivers of DC1-3') + 
  geom_bar(stat = 'identity', color = 'white', fill = colBP[1], size = 0.25) + theme_pubr() + 
  theme(legend.position = 'none', plot.margin = margin(0, 5, 0, 0, 'cm')) + 
  xlab('Percentage of overlap') + ylab('Gene') + scale_x_continuous(expand = c(0, 0))
dev.off()

genesets_w <- as.list(as.data.frame(read.csv('signatures/wouters_mel_sigs.csv', na.strings = c('', NA))))
df_pw <- df[df$Freq > 0.25, ]
pdf('data/tumor/diffusion/single/plot_MBM_single_DC1-3_drivers_pathways.pdf')
hyp_dots(hypeR(as.character(df_pw$res), genesets = genesets_w), title = paste0('Wouters et al. pathway selection'), 
         abrv = 80) + theme_bw() + 
  theme(axis.title.y = element_blank(), axis.text.y = element_text(color = 'black'), 
        axis.text.x = element_text(color = 'black'), panel.grid.minor = element_blank()) + 
  scale_color_gradient(high = 'black', low = '#C51B8A')
dev.off()


# Plot MPM bar plots with genes detected in >25% of samples
res_mpm <- res[, grep('MPM', colnames(res))]
df <- as.data.frame(sort(table(res_mpm)/ncol(res_mpm), decreasing = F))
df_plot <- df[df$Freq > 0.25, ]
pdf('data/tumor/diffusion/single/plot_MPM_single_DC1-3_drivers_barplot.pdf', height = 14)
ggplot(df_plot, aes(x = Freq, y = res_mpm)) + ggtitle('MPM - Top overlapping drivers of DC1-3') + 
  geom_bar(stat = 'identity', color = 'white', fill = colBP[2], size = 0.25) + 
  theme_pubr() + theme(legend.position = 'none', plot.margin = margin(0, 5, 0, 0, 'cm')) + 
  xlab('Percentage of overlap') + ylab('Gene') + scale_x_continuous(expand = c(0, 0))
dev.off()

df_pw <- df[df$Freq > 0.25, ]
genesets_w <- as.list(as.data.frame(read.csv('signatures/wouters_mel_sigs.csv', na.strings = c('', NA))))
pdf('data/tumor/diffusion/single/plot_MPM_single_DC1-3_drivers_pathways.pdf')
hyp_dots(hypeR(as.character(df_pw$res), genesets = genesets_w), title = paste0('Wouters et al. pathway selection'), 
         abrv = 80) + theme_bw() + 
  theme(axis.title.y = element_blank(), axis.text.y = element_text(color = 'black'), 
        axis.text.x = element_text(color = 'black'), panel.grid.minor = element_blank()) + 
  scale_color_gradient(high = 'black', low = '#C51B8A')
dev.off()


print(Sys.time())