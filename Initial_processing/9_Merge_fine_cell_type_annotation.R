#!/usr/bin/env Rscript

#### Final fine cell type annotation of integrated Seurat object
#### Author: Jana Biermann, PhD

library(dplyr)
library(Seurat)
library(ggplot2)
library(gplots)
library(viridis)

# Read-in integrated object
seu <- readRDS('data/MBPM/MBPM_scn/data_MBPM_scn_anchor2000_dims50.rds')

# Read in cell type annotations from first round of subsets
sc_my_pt1 <- read.csv('data/cell_type_DEG/MBM_sc/myeloid/data_MBM_sc_myeloid_celltype_reint_pt1.csv')
sc_my_pt2 <- read.csv('data/cell_type_DEG/MBM_sc/myeloid/data_MBM_sc_myeloid_celltype.csv')
sn_my_pt1 <- read.csv('data/cell_type_DEG/MBM_sn/myeloid/data_MBM_sn_myeloid_celltype_pt1.csv')
sn_my_pt2 <- read.csv('data/cell_type_DEG/MBM_sn/myeloid/data_MBM_sn_myeloid_celltype.csv')
snP_my <- read.csv('data/cell_type_DEG/MPM_sn/myeloid/data_MPM_sn_myeloid_celltype.csv')
sc_b <- read.csv('data/cell_type_DEG/MBM_sc/bcells/data_MBM_sc_bcells_celltype.csv')
sn_b <- read.csv('data/cell_type_DEG/MBM_sn/bcells/data_MBM_sn_bcells_celltype.csv')
snP_b <- read.csv('data/cell_type_DEG/MPM_sn/bcells/data_MPM_sn_bcells_celltype.csv')
sc_t <- read.csv('data/cell_type_DEG/MBM_sc/tcells/data_MBM_sc_tcells_celltype.csv')
sn_t <- read.csv('data/cell_type_DEG/MBM_sn/tcells/data_MBM_sn_tcells_celltype.csv')
snP_t <- read.csv('data/cell_type_DEG/MPM_sn/tcells/data_MPM_sn_tcells_celltype.csv')
sc_neu <- read.csv('data/cell_type_DEG/MBM_sc/cns/data_MBM_sc_cns_celltype.csv')
sn_neu <- read.csv('data/cell_type_DEG/MBM_sn/cns/data_MBM_sn_cns_celltype.csv')
snP_neu <- read.csv('data/cell_type_DEG/MPM_sn/stromal/data_MPM_sn_stromal_celltype.csv')

sc_glob <- read.csv('data/cell_type_DEG/MBM_sc/main/data_MBPM_sc_celltype_main.csv')
sn_glob <- read.csv('data/cell_type_DEG/MBM_sn/main/data_MBPM_sn_celltype_main.csv')
snP_glob <- read.csv('data/cell_type_DEG/MPM_sn/main/data_MPM_sn_celltype_main.csv')

### Add cell type and cell cycle annotations from subsets and tidy up annotation
# cell_type_main
celltypes <- rbind.data.frame(sc_glob, sn_glob, snP_glob)
seu$barcode_all <- rownames(seu@meta.data)
seu@meta.data <- left_join(seu@meta.data, celltypes, by = 'barcode_pat')
rownames(seu@meta.data) <- seu$barcode_all
seu$cell_type_main <- ifelse(seu$malignant == 'malignant', 'Tumor cells',
                             seu$cell_type_main)
print(DimPlot(seu, label = T, group.by = 'cell_type_main', shuffle = T, raster = T,
              repel = T))

# cell_type_fine
ct_fine <- rbind.data.frame(sc_my_pt1, sc_my_pt2, sn_my_pt1, sn_my_pt2, snP_my, sc_b,
                            sn_b, snP_b, sc_t, sn_t, snP_t, sc_neu, sn_neu, snP_neu)
ct_fine$X <- NULL
seu$barcode_all <- rownames(seu@meta.data)
seu@meta.data <- left_join(seu@meta.data, ct_fine, by = 'barcode_pat')
rownames(seu@meta.data) <- seu$barcode_all

seu$cell_type_fine <- ifelse(seu$malignant == 'malignant', 'Tumor cells',
                             seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$cell_type_main == 'Endothelial cells',
                             'Endothelial cells', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$cell_type_main == 'B/Plasma cells' &
                               is.na(seu$cell_type_fine) == T, 'Doublets', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(is.na(seu$cell_type_fine) == T, 'NA', seu$cell_type_fine)
seu$cell_type_main <- ifelse(seu$cell_type_fine %in% c('CAFs', 'Pericytes'),
                             'Stromal cells', seu$cell_type_main)
seu$cell_type_main <- ifelse(seu$cell_type_fine == 'Epithelial cells', 'Epithelial cells',
                             seu$cell_type_main)
seu$cell_type_main <- ifelse(seu$cell_type_fine == 'Endothelial cells',
                             'Endothelial cells', seu$cell_type_main)
seu$cell_type_main <- ifelse(seu$cell_type_fine == 'Endothelial cells',
                             'Endothelial cells', seu$cell_type_main)
seu$cell_type_main <- ifelse(seu$cell_type_fine %in% c('Neurons', 'Astrocytes', 'Oligodendrocytes'), 
                             'CNS cells', seu$cell_type_main)
seu$cell_type_fine <- ifelse(seu$cell_type_fine == 'NA', 'Undetermined',
                             seu$cell_type_fine)
seu$cell_type_main <- ifelse(seu$cell_type_fine %in% c('Low-quality tumor', 'Low-quality tumor/immune', 
                                                       'Cycling low-quality cells'), 
                             'Low-quality cells', seu$cell_type_main)
print(DimPlot(seu, label = T, group.by = 'cell_type_fine', shuffle = T, raster = T, repel = T))


### Add cell_type_int
# T cells
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('T/NK cells') &
                              seu$cell_type_fine %in% c('CD8+ T cells TOX+', 'CD8+ T cells TCF7+'),
                            'CD8+ T cells', NA)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('T/NK cells') &
                              seu$cell_type_fine %in% c('CD4+ T cells', 'Tfh-like cells'), 'CD4+ T cells',
                            seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('T/NK cells') &
                              seu$cell_type_fine %in% c('Tregs'), 'Tregs', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('T/NK cells') &
                              seu$cell_type_fine %in% c('NK cells'), 'NK cells', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('T/NK cells') &
                              seu$cell_type_fine %in% c('Cycling cells'), 'Cycling cells', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_fine %in% c('Contamination'), 'Contamination',
                            seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_fine %in% c('Myeloid doublets',
                                                      'T-cell doublets'), 'Doublets', seu$cell_type_int)

# Myeloid cells
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('Myeloid cells') &
                              seu$cell_type_fine %in% c('Monocytes'), 'Monocytes', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('Myeloid cells') &
                              seu$cell_type_fine %in% c('cDC1', 'cDC2', 'DC3'), 'Dendritic cells',
                            seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('Myeloid cells') &
                              seu$cell_type_fine %in% c('Microglia'), 'Microglia', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('Myeloid cells') &
                              seu$cell_type_fine %in% c('MDM', 'MDM FTL+', 'MDM M1-like', 'MDM M2-like',
                                                        'Proinflammatory MDM (glycolysis)', 'Proinflammatory MDM (zinc)',
                                                        'Proinflammatory MDM'), 'MDM', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('Myeloid cells') &
                              seu$cell_type_fine %in% c('Cycling cells'), 'Cycling cells', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('Myeloid cells') &
                              seu$cell_type_fine %in% c('Mast cells'), 'Mast cells', seu$cell_type_int)

# B cells
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('B/Plasma cells') &
                              seu$cell_type_fine %in% c('Plasma cells'), 'Plasma cells', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('B/Plasma cells') &
                              seu$cell_type_fine %in% c('Activated B cells', 'Naïve B cells'), 'B cells',
                            seu$cell_type_int)

# Stromal, CNS, endothelial, tumor
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('CNS cells') &
                              seu$cell_type_fine %in% c('Astrocytes', 'Oligodendrocytes', 'Neurons'),
                            'CNS cells', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('Stromal cells') &
                              seu$cell_type_fine %in% c('CAFs', 'Pericytes'), 'Stromal cells', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('Epithelial cells') &
                              seu$cell_type_fine %in% c('Epithelial cells'), 'Epithelial cells', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_fine %in% c('Low-quality tumor',
                                                      'Low-quality tumor/immune', 'Undetermined'), 'Low-quality cells', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('Endothelial cells') &
                              seu$cell_type_fine %in% c('Endothelial cells'), 'Endothelial cells',
                            seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('Tumor cells') &
                              seu$cell_type_fine %in% c('Tumor cells'), 'Tumor cells', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_fine %in% c('Doublets', 'Myeloid doublets',
                                                      'T-cell doublets'), 'Doublets', seu$cell_type_int)

print(DimPlot(seu, label = T, group.by = 'cell_type_int', shuffle = T, raster = T,
              repel = T))


### Update cell cycle
seu$cell_cycle <- ifelse(is.na(seu$cell_cycle.y) == T, seu$cell_cycle.x, seu$cell_cycle.y)
seu$cell_cycle.x <- NULL
seu$cell_cycle.y <- NULL
print(DimPlot(seu, label = T, group.by = 'cell_cycle', shuffle = T, raster = T,
              repel = T))


### Summarize myeloid cells
seu$cell_type_myeloid_granular <- seu$cell_type_fine
seu$cell_type_myeloid <- ifelse(seu$cell_type_fine %in% c('MDM M1-like', 'MDM M2-like', 'Proinflammatory MDM', 
                                                          'Proinflammatory MDM (glycolysis)', 'Proinflammatory MDM (zinc)'), 
                                'MDM', seu$cell_type_fine)
seu$cell_type_fine <- seu$cell_type_myeloid


### Save object
saveRDS(seu, 'data/MBPM/data_MBPM_scn.rds')
ct <- seu@meta.data %>%
  dplyr::select('barcode_pat', 'cell_type_main', 'cell_type_int', 'cell_type_fine',
                'cell_cycle')
write.csv(ct, 'data/MBPM/data_MBPM_scn_celltype_assignment.csv', row.names = F)
