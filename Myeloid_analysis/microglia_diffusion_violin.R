#!/usr/bin/env Rscript

### title: Microglia (MG) analysis (DEG, diffusion component analysis, violin plots) 
### author: Jana Biermann, PhD

print(Sys.time())

library(Seurat)
library(destiny)
library(dplyr)
library(ggplot2)
library(gplots)
'%notin%' <- Negate('%in%')


colBP <- c('#A80D11', '#008DB8')
colSCSN <- c('#E1AC24', '#288F56')

# Read-in Seurat object
seu_all <- readRDS('data/MBPM/data_MBPM.rds')

#### Single nuclei workflow ####
seu <- subset(seu_all, cell_type_fine %in% c('Microglia') & sequencing == 'Single nuclei' & organ == 'Brain')
seu <- ScaleData(seu) %>% RunPCA() %>% RunUMAP(dims = 1:30)
ElbowPlot(seu, ndims = 50)
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.3)

# Save object
seu <- readRDS('data/MBPM/myeloid/data_mg_sn.rds')

# Identify markers for the two MG clusters
markers <- FindMarkers(seu, ident.1 = '0', ident.2 = '1', min.pct = 0, logfc.threshold = 0, assay = 'RNA')
markers$gene <- rownames(markers)
ifelse(!dir.exists(file.path(paste0('data/DEG/microglia/'))), dir.create(file.path(paste0('data/DEG/microglia/')), recursive = T), FALSE)
write.csv(markers, paste0('data/DEG/microglia/markers_MBM_sn_microglia_cluster.csv'), row.names = F)

# Rename clusters
seu$mg_cluster <- Idents(seu)
seu$mg_cluster <- ifelse(seu$mg_cluster == '1', 'Homeostatic MG', 'Activated MG')

# Apply MG activated signautre to Seurat object (top 100 genes from DEG)
sigs <- read.csv('signatures/myeloid_signatures.csv', na.strings = c('', 'NA'))
for (c in 1:ncol(sigs)) {
  sig <- as.character(na.omit(sigs[, c]))
  seu <- AddModuleScore(object = seu, features = list(sig), name = colnames(sigs[c]), assay = 'RNA', search = F)
}


### Diffusion component analysis
es <- as.ExpressionSet(as.data.frame(t(seu@assays$integrated@data)))
es@phenoData@data <- seu@meta.data
dm <- DiffusionMap(es, verbose = T, n_pcs = 20)

# Add DCs to Seurat object
seu <- AddMetaData(seu, dm@eigenvectors, col.name = colnames(dm@eigenvectors))

# Remove M2 contamination cluster
seu$mg_cluster2 <- ifelse(seu$DC1 < 0 & seu$DC2 > 0, 'M2', seu$mg_cluster)
seu <- subset(seu, mg_cluster2 != 'M2')
seu <- ScaleData(seu) %>% RunPCA() %>% RunUMAP(dims = 1:30)
ElbowPlot(seu, ndims = 50)
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.2)
seu$mg_cluster_new <- Idents(seu)
seu$mg_cluster_new <- ifelse(seu$mg_cluster_new == '1', 'Homeostatic MG', 'Activated MG')


### Plots
ifelse(!dir.exists(file.path(paste0('data/diffusion/myeloid/microglia'))), dir.create(file.path(paste0('data/diffusion/myeloid/microglia')), recursive = T), FALSE)
pdf('data/diffusion/myeloid/microglia/plots_dm_microglia_sn.pdf')
# Diffusion maps
set.seed(1)
ggplot(dm_df, aes(DC1, DC2, col = seu@meta.data[rownames(dm_df), 'mg_cluster_new'])) + 
  geom_point(alpha = 1, shape = 16, size = 1.5) + theme_classic() + 
  ggtitle('mg_cluster_new') + theme(legend.title = element_blank())

ggplot(dm_df, aes(DC1, DC2, col = seu@assays$integrated@data['SPP1', rownames(dm_df)])) + 
  geom_point(alpha = 1, shape = 16, size = 1.5) + theme_classic() + 
  ggtitle('SPP1') + theme(legend.title = element_blank()) + scale_color_viridis()

ggplot(dm_df, aes(DC1, DC2, col = seu@assays$integrated@data['HAVCR2', rownames(dm_df)])) + 
  geom_point(alpha = 1, shape = 16, size = 1.5) + theme_classic() + 
  ggtitle('HAVCR2') + theme(legend.title = element_blank()) + scale_color_viridis()

# Violin plots
genes <- c('APOE', 'CST3', 'FTL', 'SPP1', 'B2M', 'PLD3', 'CD74', 'CD86')
plts <- VlnPlot(seu, features = genes, cols = c('deeppink', 'midnightblue'), pt.size = 0, 
                group.by = 'mg_cluster_new', combine = F)
plts <- lapply(plts, function(x) {
  x + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), legend.position = 'none')
})
CombinePlots(plts, legend = 'bottom')
dev.off()



#### Single cell workflow #####
seu <- subset(seu_all, cell_type_fine == 'Microglia' & sequencing == 'Single cells' & organ == 'Brain')
seu <- ScaleData(seu) %>% RunPCA() %>% RunUMAP(dims = 1:30)
ElbowPlot(seu)
seu <- FindNeighbors(seu, dims = 1:15)
seu <- FindClusters(seu, resolution = 0.15)
seu <- subset(seu, integrated_snn_res.0.15 != 2)

# Re-integrate
DefaultAssay(seu) <- 'RNA'
obj.list <- SplitObject(seu, split.by = 'patient')
anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:30, k.filter = 50)
seu <- IntegrateData(anchorset = anchors, dims = 1:30)

# Rerun Seurat workflow
seu <- ScaleData(seu) %>% RunPCA() %>%RunUMAP(dims = 1:30) %>% 
  FindNeighbors(dims = 1:30) %>% FindClusters(resolution = 0.2)
seu <- FindClusters(seu, resolution = 0.2)
seu <- subset(seu, integrated_snn_res.0.2 != 2)
seu <- ScaleData(seu) %>% RunPCA() %>% RunUMAP(dims = 1:30) %>%
  FindNeighbors(dims = 1:30) %>% FindClusters(resolution = 0.2)

# Save object
seu <- readRDS('data/MBPM/myeloid/data_mg_sc.rds')

# Identify markers for the two MG clusters
markers <- FindMarkers(seu, ident.1 = '0', ident.2 = '1', min.pct = 0, logfc.threshold = 0, assay = 'RNA')
markers$gene <- rownames(markers)
ifelse(!dir.exists(file.path(paste0('data/DEG/microglia/'))), dir.create(file.path(paste0('data/DEG/microglia/')), recursive = T), FALSE)
write.csv(markers, paste0('data/DEG/microglia/markers_MBM_sc_microglia_cluster.csv'), row.names = F)

# Rename clusters
seu$mg_cluster <- Idents(seu)
seu$mg_cluster <- ifelse(seu$mg_cluster == '1', 'Homeostatic MG', 'Activated MG')

# Apply MG activated signautre to Seurat object (top 100 genes from DEG)
sigs <- read.csv('signatures/myeloid_signatures.csv', na.strings = c('', 'NA'))
for (c in 1:ncol(sigs)) {
  sig <- as.character(na.omit(sigs[, c]))
  seu <- AddModuleScore(object = seu, features = list(sig), name = colnames(sigs[c]), assay = 'RNA', search = F)
}

pdf('data/MBPM/myeloid/plots_vln_mg_cluster_sc.pdf')
genes <- c('APOE', 'CST3', 'FTL', 'SPP1', 'B2M', 'PLD3', 'CD74', 'CD86')
plts <- VlnPlot(seu, features = genes, cols = c('deeppink', 'midnightblue'), pt.size = 0, assay = 'RNA', 
                group.by = 'mg_cluster', combine = F)
plts <- lapply(plts, function(x) {
  x + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), legend.position = 'none')
})
CombinePlots(plts, legend = 'bottom')
dev.off()


print(Sys.time())