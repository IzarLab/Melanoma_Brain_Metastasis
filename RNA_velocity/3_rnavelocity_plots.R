#!/usr/bin/env Rscript

### title: RNA velocity plots in R
### author: Jana Biermann, PhD

library(plyr)
library(dplyr)
library(ggplot2)
library(ggrastr)
library(scales)
library(viridis)
library(patchwork)
library(reticulate)

### Choose one label
label <- 'MBM_sn'
label <- 'MPM_sn'


# Load python libraries
scv <- import('scvelo')
sc <- import('scanpy')

colMM <- c('#0099c6', '#990099', '#66aa00')
colMM <- c(colMM, sc$pl$palettes$default_102)

# Load adata (Harmony-integrated)
if (label == 'MBM_sn') {
  adata <- scv$read('data/RNA_velocity/Harmony_MDM_MBM_velocity.h5ad')
  markers <- c('CD163L1', 'SELENOP', 'F13A1', 'DAB2', 'SIGLEC1', 'FTL',
               'FTH1', 'SPP1', 'MERTK', 'MRC1', 'CD83', 'CD38', 'CD80', 'FCGR1A',
               'LYZ', 'VCAN', 'FCN1', 'S100A9')
}
if (label == 'MPM_sn') {
  adata <- scv$read('data/RNA_velocity/Harmony_MDM_MPM_velocity.h5ad')
  markers <- c('CD163L1', 'SELENOP', 'F13A1', 'DAB2', 'SIGLEC1', 'FTL',
               'FTH1', 'SPP1', 'MERTK', 'MRC1', 'CD83', 'CD38', 'CD80', 'FCGR1A',
               'LYZ', 'VCAN', 'FCN1')
}
adata


# Apply signatures
sigs <- read.csv('~/signatures/myeloid_signatures.csv', na.strings = '')
sc$tl$score_genes(adata, gene_list = as.character(na.omit(sigs[, 'M1_ada'])),
                  score_name = 'M1_ada2')
sc$tl$score_genes(adata, gene_list = as.character(na.omit(sigs[, 'M2_ada'])),
                  score_name = 'M2_ada2')

# UMAP
scv$pl$scatter(adata, legend_loc = 'right margin', size = 60,
               color = c('latent_time', 'velocity_pseudotime', 'cell_type_fine',
                         'M1_ada2', 'M2_ada2', 'patient', 'root_cells', 'end_points'),
               color_map = 'Blues', perc = c(5, 95), palette = colMM, dpi = 300,
               ncols = as.integer(3), save = paste0('data/RNA_velocity/plots_MDM_',
                                                    label, '_sigs_umap.pdf'))

# Stream
scv$pl$velocity_embedding_stream(adata, basis = 'umap',
                                 legend_loc = 'right margin', 
                                 color = c('latent_time', 'velocity_pseudotime', 
                                           'cell_type_fine', 'M1_ada2', 'M2_ada2', 
                                           'patient', 'root_cells', 'end_points'), 
                                 color_map = 'Blues', palette = colMM,
                                 ncols = as.integer(3), alpha = 0.5, size = 80,
                                 perc = c(5, 95), 
                                 save = paste0('data/RNA_velocity/plots_MDM_',
                                                                label, '_stream.svg'))

# Arrows
scv$pl$velocity_embedding_grid(adata, basis = 'umap', legend_loc = 'right margin', 
                               color = c('latent_time', 'velocity_pseudotime', 
                                         'cell_type_fine', 'M1_ada2', 'M2_ada2', 
                                         'patient', 'root_cells', 'end_points'),
                               color_map = 'Blues', palette = colMM, ncols = as.integer(3),
                               arrow_size = 3, arrow_length = 3, density = 0.5, 
                               perc = c(5, 95), dpi = 300, linewidth = 0.7, 
                               arrow_color = 'black',
                               save = paste0('data/RNA_velocity/plots_MDM_', label,
                                             '_arrows.pdf'))

# UMAP markers
scv$pl$scatter(adata, legend_loc = 'right margin', size = 60, color = markers,
               color_map = 'Blues', perc = c(5, 95), ncols = as.integer(5), dpi = 300,
               save = paste0('data/RNA_velocity/plots_MDM_', label, '_markers_umap.pdf'))

# PAGA
scv$pl$paga(adata, basis = 'umap', size = 50, alpha = 0.8,
            legend_loc = 'right margin', color = 'cell_type_fine',
            min_edge_width = 2, node_size_scale = 1.5, dpi = 300,
            save = paste0('data/RNA_velocity/plots_MDM_', label, '_paga.pdf'))


# Speed and coherence
scv$pl$scatter(adata, color = c('velocity_length',
                                'velocity_confidence'), legend_loc = 'right margin',
               size = 40, color_map = 'coolwarm',
               perc = c(5, 95), dpi = 300, alpha = 0.7,
               save = paste0('data/RNA_velocity/plots_MDM_',
                             label, '_speed_coherence.pdf'))

# Pseudotime histogram
df <- data.frame(barcode = rownames(adata$obs),
                 cell_type_fine = adata$obs['cell_type_fine'],
                 M2_ada2 = adata$obs['M2_ada2'], M1_ada2 = adata$obs['M1_ada2'],
                 velocity_pseudotime = adata$obs['velocity_pseudotime'],
                 latent_time = adata$obs['latent_time'],
                 root_cells = adata$obs['root_cells'], end_points = adata$obs['end_points'])

pdf(paste0('data/RNA_velocity/plots_MDM_', label, '_histogram.pdf'))
for (i in names(df)[3:ncol(df)]) {
  df[, i] <- rescale(df[, i])
  wilc <- wilcox.test(df[, i] ~ df$cell_type_fine)
  p1 <- ggplot(df, aes(x = df[, i], fill = cell_type_fine)) +
    geom_histogram(aes(y = ..density..),
                   position = 'identity', bins = 300,
                   alpha = 0.4) + geom_density(alpha = 0.2) +
    scale_fill_manual(values = colMM) + theme_classic() +
    geom_boxplot(width = 0.5, outlier.color = NA) +
    ggtitle(paste0(i, '  (Wilcoxon P: ', (scientific(wilc$p.value, 3)), ')')) +
    theme(axis.title.x = element_blank())
  p1
  
  heat <- df[, -2]
  heat[2:ncol(heat)] <- apply(heat[2:ncol(heat)], 2, rescale)
  heat_long <- reshape2::melt(heat, var.name = 'barcode')
  order_x <- heat %>%
    dplyr::arrange(eval(parse(text = i))) %>%
    dplyr::select(barcode)
  order_y <- c(paste0(names(heat)[which(names(heat) != i)]), i)
  p2 <- ggplot(heat_long, 
               aes(x = factor(barcode, level = order_x$barcode), 
                   y = factor(variable, level = order_y), 
                   fill = value)) + 
    rasterise(geom_tile(),dpi = 300) + 
    scale_fill_viridis() + 
    theme_classic() +
    xlab('Cells') + 
    labs(fill = 'Score') +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_blank())
  p2
  
  print(p1/p2)
}
dev.off()
