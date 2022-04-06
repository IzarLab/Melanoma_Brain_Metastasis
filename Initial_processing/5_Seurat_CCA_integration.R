#!/usr/bin/env Rscript

#### Seurat CCA integration of MBPM_scn/MBM_sc/MBM_sn/MPM_sn with anchors, dims and cohort provided
#### Author: Jana Biermann, PhD

library(dplyr)
library(Seurat)
print(Sys.time())

anchor <- as.numeric(commandArgs()[6])  #2000 anchors
dims <- as.numeric(commandArgs()[7])  #MBPM_scn:50; MBM_sc: 30; MBM_sn: 30; MPM_sn: 30 
cohort <- commandArgs()[8]

# Load objects
pat_list <- c('MBM01_sc', 'MBM02_sc', 'MBM03_sc', 'MBM04_sc', 'MBM05_sc', 'MBM05_sn',
              'MBM06_sn', 'MBM07_sn', 'MBM08_sn', 'MBM09_sn', 'MBM10_sn', 'MBM11_sn', 'MBM12_sn',
              'MBM13_sn', 'MBM14_sn', 'MBM15_sn', 'MBM16_sn', 'MBM17_sn', 'MBM18_sn', 'MBM19_sn',
              'MBM20_sn', 'MBM21_sn', 'MPM01_sn', 'MPM02_sn', 'MPM04_sn', 'MPM05_sn', 'MPM06_sn',
              'MPM07_sn', 'MPM08_sn', 'MPM09_sn', 'MPM10_sn', 'MPM11_sn')

if (cohort == 'MBM_sc') {
  pat_list <- grep('sc', pat_list, value = T)
}

if (cohort == 'MBM_sn') {
  pat_list <- grep('sc', pat_list, value = T, invert = T)
  pat_list <- grep('MBM', pat_list, value = T)
}

if (cohort == 'MPM_sn') {
  pat_list <- grep('MPM', pat_list, value = T)
}

for (pat in pat_list) {
  assign(pat, readRDS(paste0('data/', pat, '/data_', pat, '_cnv.rds')))
}

# Create obj.list
obj.list <- NULL
for (obj in grep(paste(c('_sc', '_sn'), collapse = '|'), ls(), value = T)) {
  obj.list <- c(obj.list, eval(parse(text = obj)))
}

# Integration
anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = anchor,
                                  dims = 1:dims)
rm(seu)
seu <- IntegrateData(anchorset = anchors, dims = 1:dims)
rm(anchors, obj.list)

# Seurat workflow
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 100)
seu <- RunUMAP(seu, dims = 1:dims)
seu <- FindNeighbors(seu, dims = 1:dims)
seu <- FindClusters(seu)

# Cell cycle
seu <- CellCycleScoring(seu, s.features = cc.genes.updated.2019$s.genes,
                        g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)
seu$cell_cycle <- ifelse(seu$G2M.Score > 0.1 | seu$S.Score > 0.25, 'cycling', 'non-cycling')

# Save object
ifelse(!dir.exists(file.path(paste0('data/MBPM/', cohort))),
       dir.create(file.path(paste0('data/MBPM/', cohort)), recursive = T), FALSE)
saveRDS(seu, file = paste0('data/MBPM/', cohort, '/data_', cohort, '_anchor', anchor,
                           '_dims', dims, '.rds'))

print(Sys.time())