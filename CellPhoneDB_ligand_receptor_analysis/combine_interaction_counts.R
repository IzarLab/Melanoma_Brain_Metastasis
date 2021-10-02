#!/usr/bin/env Rscript

### title: Analysis of interaction counts between cell types of individual samples obtained through CellPhoneDB
### author: Jana Biermann, PhD 
### Based on: https://github.com/IzarLab/CUIMC-NYP_COVID_autopsy_lung/blob/main/code/CellPhoneDB_ligand_receptor_analysis/CPDB_interactionCountsCombined.R

library(pheatmap)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(Seurat)
'%notin%' <- Negate('%in%')


#### Getting peripheral interaction counts ####

# Read-in count_network from each individual sample and combine into a single file
pats <- c('MPM01_sn', 'MPM02_sn', 'MPM04_sn', 'MPM05_sn')
counts_mat <- NULL
for (pat in pats) {
  print(pat)
  count_network <- read.delim(paste0('data/CellPhoneDB/', pat, '/cell_type_fine/count_network.txt'), 
                              row.names = NULL, stringsAsFactors = FALSE)
  
  # Combining the source and target names for sorting
  count_network$source_target <- paste0(count_network$SOURCE, '_', count_network$TARGET)
  
  # Removing spaces in the names for sorting and comparison
  count_network$source_target <- gsub(' ', '-', count_network$source_target)
  
  # Sorting
  count_network <- count_network[order(count_network$source_target), ]
  
  # Getting only the required columns
  temp_counts <- count_network[, c('count', 'source_target')]
  colnames(temp_counts) <- paste0(pat, '-', colnames(temp_counts))
  temp_counts$source_target <- temp_counts[, paste0(pat, '-source_target')]
  
  if (is.null(counts_mat)) {
    counts_mat <- cbind(count_network$source_target, temp_counts)
    colnames(counts_mat)[1] <- 'source_target'
    counts_mat <- counts_mat[, 1:3]
  } else {
    counts_mat <- full_join(counts_mat, temp_counts, by = 'source_target')
  }
}
counts_mat <- cbind.data.frame(colsplit(string = counts_mat$source_target, 
                                        pattern = '_', names = c('SOURCE', 'TARGET')), 
                               counts_mat[, 2:ncol(counts_mat)])

# Keeping only the counts and other required columns
counts_mat <- counts_mat[, c(1, 2, grep('count', colnames(counts_mat)))]

# Save counts matrix
ifelse(!dir.exists(file.path('data/CellPhoneDB/cell_type_fine')), 
       dir.create(file.path('data/CellPhoneDB/cell_type_fine'), recursive = T), FALSE)
write.csv(counts_mat, 'data/CellPhoneDB/cell_type_fine/ligandReceptors_combinedCountsMat_peripheral.csv', row.names = F)


#### Getting brain interaction counts ####

pats <- c('MBM05_sn', 'MBM06_sn', 'MBM07_sn', 'MBM08_sn', 'MBM09_sn', 'MBM10_sn', 'MBM11_sn', 'MBM12_sn', 
          'MBM13_sn', 'MBM14_sn', 'MBM15_sn', 'MBM16_sn', 'MBM17_sn', 'MBM18_sn', 'MBM19_sn', 'MBM20_sn', 'MBM21_sn')

counts_mat <- NULL
for (pat in pats) {
  print(pat)
  count_network <- read.delim(paste0('data/CellPhoneDB/', pat, '/cell_type_fine/count_network.txt'), 
                              row.names = NULL, stringsAsFactors = FALSE)
  
  # Combining the source and target names for sorting
  count_network$source_target <- paste0(count_network$SOURCE, '_', count_network$TARGET)
  
  # Removing spaces in the names for sorting and comparison
  count_network$source_target <- gsub(' ', '-', count_network$source_target)
  
  # Sorting
  count_network <- count_network[order(count_network$source_target), ]
  
  # Getting only the required columns
  temp_counts <- count_network[, c('count', 'source_target')]
  colnames(temp_counts) <- paste0(pat, '-', colnames(temp_counts))
  temp_counts$source_target <- temp_counts[, paste0(pat, '-source_target')]
  
  if (is.null(counts_mat)) {
    counts_mat <- cbind(count_network$source_target, temp_counts)
    colnames(counts_mat)[1] <- 'source_target'
    counts_mat <- counts_mat[, 1:3]
  } else {
    counts_mat <- full_join(counts_mat, temp_counts, by = 'source_target')
  }
}
counts_mat <- cbind.data.frame(colsplit(string = counts_mat$source_target, 
                                        pattern = '_', names = c('SOURCE', 'TARGET')), 
                               counts_mat[, 2:ncol(counts_mat)])

# Keeping only the counts and other required columns
counts_mat <- counts_mat[, c(1, 2, grep('count', colnames(counts_mat)))]

# Save counts matrix
write.csv(counts_mat, 'data/CellPhoneDB/cell_type_fine/ligandReceptors_combinedCountsMat_brain.csv', 
          row.names = F)


#### Get median counts for brain and peripheral samples separately #####

# peripheral
p_counts <- read.csv(file = 'data/CellPhoneDB/cell_type_fine/ligandReceptors_combinedCountsMat_peripheral.csv')
temp <- apply(p_counts[, 3:ncol(p_counts)], MARGIN = 1, FUN = median, na.rm = T)
p_counts_med <- cbind(p_counts[, 1:2], temp)
colnames(p_counts_med)[3] <- 'count'
head(p_counts_med)
p_matrix <- dcast(p_counts_med, p_counts_med$SOURCE ~ p_counts_med$TARGET, fill = 0)
colnames(p_matrix)[1] <- ''
write.csv(p_matrix, 'data/CellPhoneDB/cell_type_fine/ligandReceptors_counts_median_peripheral.csv', row.names = F)

# brain
b_counts <- read.csv(file = 'data/CellPhoneDB/cell_type_fine/ligandReceptors_combinedCountsMat_brain.csv')
temp <- apply(b_counts[, 3:ncol(b_counts)], MARGIN = 1, FUN = median, na.rm = T)
b_counts_med <- cbind(b_counts[, 1:2], temp)
colnames(b_counts_med)[3] <- 'count'
head(b_counts_med)
b_matrix <- dcast(b_counts_med, b_counts_med$SOURCE ~ b_counts_med$TARGET, fill = 0)
colnames(b_matrix)[1] <- ''
write.csv(b_matrix, 'data/CellPhoneDB/cell_type_fine/ligandReceptors_counts_median_brain.csv', row.names = F)


#### Median combined counts heatmap ####

# peripheral
colnames(p_matrix) <- rownames(p_matrix)
pdf('data/CellPhoneDB/cell_type_fine/heatmap_median_combined_counts_peripheral.pdf')
pheatmap(p_matrix, show_rownames = T, show_colnames = T, scale = 'none', cluster_cols = T, border_color = 'white', 
         cluster_rows = T, clustering_method = 'ward.D2', main = 'Peripheral median combined counts cell_type_fine', 
         color = colorRampPalette(c('dodgerblue4', 'peachpuff', 'deeppink4'))(1000))
dev.off()

# brain
colnames(b_matrix) <- rownames(b_matrix)
pdf('data/CellPhoneDB/cell_type_fine/heatmap_median_combined_counts_brain.pdf')
pheatmap(b_matrix, show_rownames = T, show_colnames = T, scale = 'none', cluster_cols = T, border_color = 'white', 
         cluster_rows = T, clustering_method = 'ward.D2', main = 'Brain median combined counts cell_type_fine', 
         color = colorRampPalette(c('dodgerblue4', 'peachpuff', 'deeppink4'))(1000))
dev.off()

