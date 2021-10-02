#!/usr/bin/env Rscript

### title: Plot pie charts of CellPhoneDB results 
### author: Jana Biermann, PhD

library(ggplot2)
library(dplyr)
library(reshape2)
library(poolr)
library(scales)


#### Brain ####

# Collect p values for selected interactions and cell types
pats_b <- c('MBM05_sn', 'MBM06_sn', 'MBM07_sn', 'MBM08_sn', 'MBM09_sn', 'MBM10_sn', 'MBM11_sn', 'MBM12_sn', 
            'MBM13_sn', 'MBM14_sn', 'MBM15_sn', 'MBM16_sn', 'MBM17_sn', 'MBM18_sn', 'MBM19_sn', 'MBM20_sn')
sel_intr <- c('CD74_APP', 'CD74_COPA', 'COPA_SORT1', 'EGFR_COPA', 'NRG2_ERBB4', 'NRG3_ERBB4', 'SPP1_CD44', 
              'LGALS9_MET')
sel_cell <- c('Astrocytes', 'Oligodendrocytes', 'Microglia', 'Neurons')

df_pval_b <- NULL
for (pat in pats_b) {
  temp_pval = read.table(paste0('data/CellPhoneDB/', pat, '/cell_type_fine/pvalues.txt'), header = T, 
                         stringsAsFactors = F, sep = '\t', comment.char = '', check.names = F)
  tmp_df <- temp_pval[, colnames(temp_pval)[c(2, 12:ncol(temp_pval))]]
  
  # Select cell types and interactions
  tmp_df <- tmp_df[, c('interacting_pair', grep('Tumor', colnames(tmp_df), value = T))]
  tmp_df <- tmp_df[tmp_df$interacting_pair %in% sel_intr, ]
  tmp_long <- reshape2::melt(tmp_df, id.vars = 'interacting_pair', value.name = 'pval', variable.name = 'cell')
  tmp_long$pval_label <- ifelse(tmp_long$pval < 0.05, 'sig_sample', 'non_sig_sample')
  tmp_long$patient <- pat
  
  # Combine
  df_pval_b <- rbind.data.frame(df_pval_b, tmp_long)
}

df_pval_b$intr_cell <- paste0(df_pval_b$interacting_pair, '_', df_pval_b$cell)
df_pval_b <- subset(df_pval_b, cell %in% unique(grep(paste(sel_cell, collapse = '|'), df_pval_b$cell, 
                                                     value = T)) | cell == 'Tumor cells|Tumor cells')
df_pval_b$cell <- droplevels(df_pval_b$cell)


# Read-in median means
meds_b <- read.csv(paste0('data/CellPhoneDB/cell_type_fine/MedianMeans_brain.csv'), check.names = FALSE, 
                   row.names = 1)

# Select cell types and interactions and prepare plot
meds_b <- meds_b[sel_intr, ]
meds_b <- meds_b[, grep('Tumor', colnames(meds_b), value = T)]
meds_b <- meds_b[, grep(paste(unique(df_pval_b$cell), collapse = '|'), colnames(meds_b), value = T)]
meds_b$intr <- rownames(meds_b)
meds_b_melt <- reshape2::melt(meds_b, id.vars = 'intr', value.name = 'median_means', variable.name = 'cell')

meds_b_melt$intr_cell <- paste0(meds_b_melt$intr, '_', meds_b_melt$cell)
plot_b <- left_join(df_pval_b, meds_b_melt[, c('median_means', 'intr_cell')], by = 'intr_cell')

# Remove NAs
for (intr_cel in unique(plot_b$intr_cell)) {
  sub <- subset(plot_b, plot_b$intr_cell == intr_cel)
  if (all(sub$pval_label == 'non_sig_sample') == T) {
    plot_b[plot_b$intr_cell == intr_cel, 'median_means'] <- rep(0, length(plot_b[plot_b$intr_cell == 
                                                                                   intr_cel, 'median_means']))
  }
}

# Plot
pdf('data/CellPhoneDB/cell_type_fine/plots_pies_brain.pdf')
ggplot(data = plot_b, aes(x = '', y = '', fill = pval_label, alpha = median_means)) + 
  geom_bar(stat = 'identity', position = 'stack') + coord_polar('y', start = 0) + 
  scale_fill_manual(values = c('#54BCD1', '#C70E7B')) + facet_grid(cell ~ interacting_pair, drop = F) + 
  theme_void() + theme(strip.text.y = element_text(vjust = 0.5, hjust = 0), 
                       strip.text.x = element_text(vjust = 0.5, hjust = 0, angle = 90))
dev.off()



#### Peripheral #####

# Collect p values for selected interactions and cell types
pats_p <- c('MPM01_sn', 'MPM02_sn', 'MPM04_sn', 'MPM05_sn')
sel_intr_p <- c('CADM1_CADM1', 'CADM1_NECTIN3', 'CD226_NECTIN2', 'COL16A1_a11b1 complex', 'COL16A1_a1b1 complex', 
                'COL16A1_a2b1 complex', 'COL18A1_a1b1 complex', 'COL19A1_a11b1 complex', 'COL19A1_a1b1 complex', 
                'COL19A1_a2b1 complex', 'COL6A1_a1b1 complex', 'LAMC1_a2b1 complex', 'NECTIN2_NECTIN3', 
                'TIGIT_NECTIN2', 'TIGIT_NECTIN3')

df_pval_p <- NULL
for (pat in pats_p) {
  temp_pval = read.table(paste0('data/CellPhoneDB/', pat, '/cell_type_fine/pvalues.txt'), header = T, 
                         stringsAsFactors = F, sep = '\t', comment.char = '', check.names = F)
  
  # Select cell types and interactions
  tmp_df <- temp_pval[, colnames(temp_pval)[c(2, 12:ncol(temp_pval))]]
  tmp_df <- tmp_df[, c('interacting_pair', grep('Tumor', colnames(tmp_df), value = T))]
  tmp_df <- tmp_df[tmp_df$interacting_pair %in% sel_intr_p, ]
  tmp_long <- reshape2::melt(tmp_df, id.vars = 'interacting_pair', value.name = 'pval', variable.name = 'cell')
  tmp_long$pval_label <- ifelse(tmp_long$pval < 0.05, 'sig_sample', 'non_sig_sample')
  tmp_long$patient <- pat
  
  # Combine
  df_pval_p <- rbind.data.frame(df_pval_p, tmp_long)
}
df_pval_p$intr_cell <- paste0(df_pval_p$interacting_pair, '_', df_pval_p$cell)

# Read-in median means
meds_p <- read.csv(paste0('data/CellPhoneDB/cell_type_fine/MedianMeans_peripheral.csv'), check.names = FALSE, 
                   row.names = 1)

# Select cell types and interactions and prepare plot
meds_p <- meds_p[sel_intr_p, ]
meds_p <- meds_p[, grep('Tumor', colnames(meds_p), value = T)]

# Remove columns without signifificant interactions
meds_p <- meds_p[, colSums(meds_p, na.rm = T) > 0]
meds_p$intr <- rownames(meds_p)
meds_p_melt <- reshape2::melt(meds_p, id.vars = 'intr', value.name = 'median_means', variable.name = 'cell')

meds_p_melt$intr_cell <- paste0(meds_p_melt$intr, '_', meds_p_melt$cell)
df_pval_p <- subset(df_pval_p, intr_cell %in% meds_p_melt$intr_cell)
plot_p <- left_join(df_pval_p, meds_p_melt[, c('median_means', 'intr_cell')], by = 'intr_cell')

# Remove NAs
for (intr_cel in unique(plot_p$intr_cell)) {
  sub <- subset(plot_p, plot_p$intr_cell == intr_cel)
  if (all(sub$pval_label == 'non_sig_sample') == T) {
    plot_p[plot_p$intr_cell == intr_cel, 'median_means'] <- rep(0, length(plot_p[plot_p$intr_cell == 
                                                                                   intr_cel, 'median_means']))
  }
}
plot_p$cell <- droplevels(plot_p$cell)

# Plot
pdf('data/CellPhoneDB/cell_type_fine/plots_pies_peripheral.pdf', height = 11)
ggplot(data = plot_p, aes(x = '', y = '', fill = pval_label, alpha = median_means)) + 
  geom_bar(stat = 'identity', position = 'stack') + coord_polar('y', start = 0) + 
  scale_fill_manual(values = c('#54BCD1', '#C70E7B')) + facet_grid(cell ~ interacting_pair, drop = F) + 
  theme_void() + theme(strip.text.y = element_text(vjust = 0.5, hjust = 0), 
                       strip.text.x = element_text(vjust = 0.5, hjust = 0, angle = 90))
dev.off()

