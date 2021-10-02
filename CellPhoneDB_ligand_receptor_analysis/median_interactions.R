#!/usr/bin/env Rscript

### title: Generating medians of mean interaction scores and Stouffer-integrated p values across
### individual samples obtained through CellPhoneDB 
### author: Jana Biermann, PhD 
### Based on: https://github.com/IzarLab/CUIMC-NYP_COVID_autopsy_lung/blob/main/code/CellPhoneDB_ligand_receptor_analysis/CPDB_topInteractions_CtrlCovidCompare.R

library(ggplot2)
library(dplyr)
library(reshape2)
library(poolr)


#### Read-in p values from individual CellPhoneDB samples and integrate ####

# peripheral
pats_p <- c('MPM01_sn', 'MPM02_sn', 'MPM04_sn', 'MPM05_sn')

# Getting the shared interactions and cell-cells
intr_pairs <- NULL
cell_cells <- NULL
pval_p <- list()
for (pat in pats_p) {
  # Reading the cpdb data
  temp_all_pval = read.table(paste0('data/CellPhoneDB/', pat, '/cell_type_fine/pvalues.txt'), header = T, 
                             stringsAsFactors = F, sep = '\t', comment.char = '', check.names = F)
  
  if (!is.null(intr_pairs)) {
    intr_pairs = union(intr_pairs, temp_all_pval$interacting_pair)
    cell_cells = union(cell_cells, colnames(temp_all_pval)[12:ncol(temp_all_pval)])
  } else {
    intr_pairs = temp_all_pval$interacting_pair
    cell_cells = colnames(temp_all_pval)[12:ncol(temp_all_pval)]
  }
  temp_all_pval <- cbind.data.frame(interacting_pair = temp_all_pval[, 2], temp_all_pval[, 12:ncol(temp_all_pval)])
  pval_p[[pat]] = temp_all_pval
}
intr_pairs_p <- intr_pairs
cell_cells_p <- cell_cells

# Integrate p values and save file
vals_p <- matrix(data = NA, nrow = length(intr_pairs_p), ncol = length(cell_cells_p))
colnames(vals_p) <- cell_cells_p
rownames(vals_p) <- intr_pairs_p

for (ro in intr_pairs_p) {
  for (co in cell_cells_p) {
    collect_vals <- list()
    collect_vals <- lapply(pval_p, function(x) {
      c(collect_vals, x[x$interacting_pair == ro, co])
    })
    if (is.null(stouffer(unlist(collect_vals), na.rm = T)$p)) {
      vals_p[ro, co] <- 1
    } else {
      vals_p[ro, co] <- stouffer(unlist(collect_vals), na.rm = T)$p
    }
  }
}
vals_p[is.na(vals_p)] <- 1
write.csv(vals_p, 'data/CellPhoneDB/cell_type_fine/StoufferPvals_peripheral.csv', row.names = T)


# brain
pats_b <- c('MBM05_sn', 'MBM06_sn', 'MBM07_sn', 'MBM08_sn', 'MBM09_sn', 'MBM10_sn', 'MBM11_sn', 'MBM12_sn', 
            'MBM13_sn', 'MBM14_sn', 'MBM15_sn', 'MBM16_sn', 'MBM17_sn', 'MBM18_sn', 'MBM19_sn', 'MBM20_sn', 'MBM21_sn')

# Getting the common interactions and cell-cells
intr_pairs <- NULL
cell_cells <- NULL
pval_b <- list()
for (pat in pats_b) {
  # Reading the cpdb data
  temp_all_pval = read.table(paste0('data/CellPhoneDB/', pat, '/cell_type_fine/pvalues.txt'), header = T, 
                             stringsAsFactors = F, sep = '\t', comment.char = '', check.names = F)
  
  if (!is.null(intr_pairs)) {
    intr_pairs = union(intr_pairs, temp_all_pval$interacting_pair)
    cell_cells = union(cell_cells, colnames(temp_all_pval)[12:ncol(temp_all_pval)])
  } else {
    intr_pairs = temp_all_pval$interacting_pair
    cell_cells = colnames(temp_all_pval)[12:ncol(temp_all_pval)]
  }
  temp_all_pval <- cbind.data.frame(interacting_pair = temp_all_pval[, 2], temp_all_pval[, 12:ncol(temp_all_pval)])
  pval_b[[pat]] = temp_all_pval
}
intr_pairs_b <- intr_pairs
cell_cells_b <- cell_cells

# Integrate pvals and save file
vals_b <- matrix(data = NA, nrow = length(intr_pairs_b), ncol = length(cell_cells_b))
colnames(vals_b) <- cell_cells_b
rownames(vals_b) <- intr_pairs_b

for (ro in intr_pairs_b) {
  for (co in cell_cells_b) {
    collect_vals <- NULL
    collect_vals <- lapply(pval_b, function(x) {
      c(collect_vals, x[x$interacting_pair == ro, co])
    })
    if (is.null(unlist(collect_vals))) {
      vals_b[ro, co] <- 1
    } else {
      vals_b[ro, co] <- stouffer(unlist(collect_vals))$p
    }
  }
}
vals_b[is.na(vals_b)] <- 1
write.csv(vals_b, 'data/CellPhoneDB/cell_type_fine/StoufferPvals_brain.csv', row.names = T)


#### Caluclate median of means of interaction scores ####

# peripheral
means_p <- list()
for (pat in pats_p) {
  temp_all_means = read.table(paste0('data/CellPhoneDB/', pat, '/cell_type_fine/significant_means.txt'), 
                              header = T, stringsAsFactors = F, sep = '\t', comment.char = '', check.names = F)
  temp_all_means <- cbind.data.frame(interacting_pair = temp_all_means[, 2], temp_all_means[, 12:ncol(temp_all_means)])
  means_p[[pat]] = temp_all_means
}

# Compute the median of all means
meds_p <- matrix(data = NA, nrow = length(intr_pairs_p), ncol = length(cell_cells_p))
colnames(meds_p) <- cell_cells_p
rownames(meds_p) <- intr_pairs_p

for (ro in intr_pairs_p) {
  for (co in cell_cells_p) {
    collect_means <- list()
    collect_means <- lapply(means_p, function(x) {
      c(collect_means, x[x$interacting_pair == ro, co])
    })
    if (is.null(median(unlist(collect_means), na.rm = T))) {
      meds_p[ro, co] <- 0
    } else {
      meds_p[ro, co] <- median(unlist(collect_means), na.rm = T)
    }
  }
}
write.csv(meds_p, 'data/CellPhoneDB/cell_type_fine/MedianMeans_peripheral.csv', row.names = T)


# brain
means_b <- list()
for (pat in pats_b) {
  temp_all_means = read.table(paste0('data/CellPhoneDB/', pat, '/cell_type_fine/significant_means.txt'), 
                              header = T, stringsAsFactors = F, sep = '\t', comment.char = '', check.names = F)
  temp_all_means <- cbind.data.frame(interacting_pair = temp_all_means[, 2], temp_all_means[, 12:ncol(temp_all_means)])
  means_b[[pat]] = temp_all_means
}

# Compute the median of all means
meds_b <- matrix(data = NA, nrow = length(intr_pairs_b), ncol = length(cell_cells_b))
colnames(meds_b) <- cell_cells_b
rownames(meds_b) <- intr_pairs_b

for (ro in intr_pairs_b) {
  for (co in cell_cells_b) {
    collect_means <- list()
    collect_means <- lapply(means_b, function(x) {
      c(collect_means, x[x$interacting_pair == ro, co])
    })
    if (is.null(median(unlist(collect_means), na.rm = T))) {
      meds_b[ro, co] <- 0
    } else {
      meds_b[ro, co] <- median(unlist(collect_means), na.rm = T)
    }
  }
}
write.csv(meds_b, 'data/CellPhoneDB/cell_type_fine/MedianMeans_brain.csv', row.names = T)

