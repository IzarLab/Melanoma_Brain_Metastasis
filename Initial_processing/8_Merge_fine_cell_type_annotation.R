#!/usr/bin/env Rscript

#### Final fine cell type annotation for integrated Seurat object 
#### Author: Jana Biermann, PhD

library(dplyr)
library(Seurat)


# Read-in integrated object
seu <- readRDS('data/MBPM/data_MBPM.rds')

# Read-in cell-type annotations from first round of subsets
sc_my <- read.csv('data/cell_type_assignment/sc/myeloid/data_MBPM_sc_myeloid_celltype.csv')
sn_my <- read.csv('data/cell_type_assignment/sn/myeloid/data_MBPM_sn_myeloid_celltype.csv')
sc_b <- read.csv('data/cell_type_assignment/sc/bcells/data_MBPM_sc_bcells_celltype.csv')
sn_b <- read.csv('data/cell_type_assignment/sn/bcells/data_MBPM_sn_bcells_celltype.csv')
sc_t <- read.csv('data/cell_type_assignment/sc/tcells/data_MBPM_sc_tcells_celltype.csv')
sn_t <- read.csv('data/cell_type_assignment/sn/tcells/data_MBPM_sn_tcells_celltype.csv')
sc_neu <- read.csv('data/cell_type_assignment/sc/cns/data_MBPM_sc_cns_celltype.csv')
sn_neu <- read.csv('data/cell_type_assignment/sn/cns/data_MBPM_sn_cns_celltype.csv')
sc_glob <- read.csv('data/MBPM/sc/data_MBPM_sc_celltype_main.csv')
sn_glob <- read.csv('data/MBPM/sn/data_MBPM_sn_celltype_main.csv')


### Add cell type and cell cycle annotations from subsets and tidy up annotation cell_type_main
celltypes <- rbind.data.frame(sc_glob, sn_glob)
seu$barcode_all <- rownames(seu@meta.data)
seu@meta.data <- left_join(seu@meta.data, celltypes, by = 'barcode_all')
rownames(seu@meta.data) <- seu$barcode_all
seu$cell_type_main <- ifelse(seu$tumor == 'tumor', 'Tumor cells', seu$cell_type_main)

# cell_type_fine
celltypes <- rbind.data.frame(sc_my, sn_my, sc_b, sn_b, sc_t, sn_t, sc_neu, sn_neu)
seu$barcode_all <- rownames(seu@meta.data)
seu@meta.data <- left_join(seu@meta.data, celltypes, by = 'barcode_all')
rownames(seu@meta.data) <- seu$barcode_all

seu$cell_type_fine <- ifelse(seu$tumor == 'tumor', 'Tumor cells', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$cell_type_main == 'Endothelial cells', 'Endothelial cells', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$cell_type_main == 'T/NK cells' & seu$cell_type_fine == 'Cycling cells', 'Cycling T/NK cells', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$cell_type_main == 'Myeloid cells' & seu$cell_type_fine == 'Cycling cells', 'Cycling myeloid cells', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$cell_type_fine == 'Cycling cells', 'Cycling low-quality cells', seu$cell_type_fine)


### Updating cell_type_main
seu$cell_type_main <- ifelse(seu$cell_type_fine %in% c('CAFs', 'Cycling CAFs', 'Pericytes'), 'Stromal cells', seu$cell_type_main)
seu$cell_type_main <- ifelse(seu$cell_type_fine %in% c('Low-quality tumor', 'Low-quality tumor/immune', 'Cycling low-quality cells'), 'Low-quality cells', seu$cell_type_main)
seu$cell_type_main <- ifelse(seu$cell_type_fine %in% c('Neurons', 'Astrocytes', 'Oligodendrocytes'), 'CNS cells', seu$cell_type_main)
seu$cell_type_main <- ifelse(seu$cell_type_fine %in% c('Epithelial cells'), 'Epithelial cells', seu$cell_type_main)


### Creating cell_type_int 
# T cells
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('T/NK cells') & seu$cell_type_fine %in% c('CD8+ T cells TOX+', 'CD8+ T cells TCF7+'), 'CD8+ T cells', NA)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('T/NK cells') & seu$cell_type_fine %in% c('CD4+ T cells', 'Tfh-like cells'), 'CD4+ T cells', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('T/NK cells') & seu$cell_type_fine %in% c('Tregs'), 'Tregs', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('T/NK cells') & seu$cell_type_fine %in% c('NK cells'), 'NK cells', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('T/NK cells') & seu$cell_type_fine %in% c('Cycling T/NK cells'), 'Cycling T/NK cells', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_fine %in% c('Contamination'), 'Contamination', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_fine %in% c('Myeloid doublets', 'T-cell doublets'), 'Doublets', seu$cell_type_int)

# Myeloid cells
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('Myeloid cells') & seu$cell_type_fine %in% c('Monocytes'), 'Monocytes', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('Myeloid cells') & seu$cell_type_fine %in% c('cDC1', 'cDC2', 'DC3', 'mo-DCs'), 'Dendritic cells', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('Myeloid cells') & seu$cell_type_fine %in% c('Microglia'), 'Microglia', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('Myeloid cells') & seu$cell_type_fine %in% c('TAM MDM', 'TAM MDM FTL+', 'TAM MDM M1-like', 'TAM MDM M2-like', 'Patient-specific MDM'), 'TAM MDM', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('Myeloid cells') & seu$cell_type_fine %in% c('Cycling myeloid cells'), 'Cycling myeloid cells', seu$cell_type_int)

# B cells
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('B/Plasma cells') & seu$cell_type_fine %in% c('Plasma cells'), 'Plasma cells', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('B/Plasma cells') & seu$cell_type_fine %in% c('Activated B cells', 'Naïve B cells'), 'B cells', seu$cell_type_int)

# Stromal, CNS, endothelial, tumor cells
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('CNS cells') & seu$cell_type_fine %in% c('Astrocytes', 'Oligodendrocytes', 'Neurons'), 'CNS cells', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('Stromal cells') & seu$cell_type_fine %in% c('CAFs', 'Cycling CAFs', 'Pericytes'), 'Stromal cells', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('Epithelial cells') & seu$cell_type_fine %in% c('Epithelial cells'), 'Epithelial cells', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('Low-quality cells') & seu$cell_type_fine %in% c('Low-quality tumor', 'Low-quality tumor/immune', 'Cycling low-quality cells'), 'Low-quality cells', 
                            seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('Endothelial cells') & seu$cell_type_fine %in% c('Endothelial cells'), 'Endothelial cells', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_main %in% c('Tumor cells') & seu$cell_type_fine %in% c('Tumor cells'), 'Tumor cells', seu$cell_type_int)

# Not determined
seu$cell_type_int <- ifelse(is.na(seu$cell_type_int) == T, 'Not determined', seu$cell_type_int)
seu$cell_type_fine <- ifelse(is.na(seu$cell_type_fine) == T, 'Not determined', seu$cell_type_fine)
seu$cell_type_int <- ifelse(seu$cell_type_int %in% c('CNS cells') & seu$organ == 'Peripheral', 'Not determined', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_int %in% c('Microglia') & seu$organ == 'Peripheral', 'Not determined', seu$cell_type_int)
seu$cell_type_int <- ifelse(seu$cell_type_int %in% c('Epithelial cells') & seu$organ == 'Brain', 'Not determined', seu$cell_type_int)
seu$cell_type_fine <- ifelse(seu$cell_type_fine %in% c('Astrocytes', 'Oligodendrocytes', 'Neurons', 'Microglia') & seu$organ == 'Peripheral', 'Not determined', seu$cell_type_fine)
seu$cell_type_fine <- ifelse(seu$cell_type_fine %in% c('Epithelial cells') & seu$organ == 'Brain', 'Not determined', seu$cell_type_fine)

### Update cell cycle
seu$cell_cycle <- ifelse(is.na(seu$cell_cycle.y) == T, seu$cell_cycle.x, seu$cell_cycle.y)
seu$cell_cycle.x <- NULL
seu$cell_cycle.y <- NULL


### Save final integrated object
saveRDS(seu, 'data/MBPM/data_MBPM.rds')
