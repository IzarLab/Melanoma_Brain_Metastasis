#!/usr/bin/env Rscript

### title: Integration of MBM05_sc and MBM05_sn for direct comparison 
### author: Jana Biermann, PhD

library(dplyr)
library(Seurat)
library(gplots)
library(ggplot2)
library(scales)
library(viridis)
library(ggpubr)
library(DropletUtils)
library(SingleR)
'%notin%' <- Negate('%in%')

colSCSN <- c('#E1AC24', '#288F56')
colSCSN2 <- c('#FADA8B', '#E1AC24', '#80CFA3', '#288F56')


#### Re-integration ####
# Split integrated object
seu <- readRDS('data/MBPM/data_MBPM_scn.rds')
sc <- subset(seu, patient == 'MBM05_sc')
sn <- subset(seu, patient == 'MBM05_sn')
rm(seu)

# Integration
anchors <- FindIntegrationAnchors(object.list = c(sc, sn), dims = 1:30)
seu <- IntegrateData(anchorset = anchors, dims = 1:30)
seu <- ScaleData(seu)
seu <- RunPCA(seu, verbose = FALSE)
seu <- RunUMAP(seu, reduction = 'pca', dims = 1:30)

# Annotation
seu <- PercentageFeatureSet(seu, pattern = '^RPS', col.name = 'percent.rps', assay = 'RNA')
seu <- PercentageFeatureSet(seu, pattern = '^RPL', col.name = 'percent.rpl', assay = 'RNA')

# Apply signatures
sigs <- read.csv('~/signatures/stress_mito_signatures.csv', na.strings = c(''))
for (c in 1:ncol(sigs)) {
  sig <- as.character(na.omit(sigs[, c]))
  seu <- AddModuleScore(object = seu, features = list(sig), name = colnames(sigs[c]), assay = 'RNA', search = F)
}

# Save object
ifelse(!dir.exists(file.path('data/seq_comparison/sc_sn')), 
       dir.create(file.path('data/seq_comparison/sc_sn'), recursive = T), FALSE)
saveRDS(seu, file = paste0('data/seq_comparison/sc_sn/integrated_sc_sn.rds'))

# Plots for cell-type-specific signautres
df_sigs <- seu@meta.data %>%
  select(sequencing, cell_type_int, 'stress_module1', 'mito_module1', 'IFN_module1') %>%
  group_by(sequencing, cell_type_int) %>%
  mutate(mean_stress = mean(stress_module1, na.rm = T)) %>%
  mutate(mean_mito = mean(mito_module1, na.rm = T)) %>%
  mutate(mean_inf = mean(IFN_module1, na.rm = T))

g_stress <- ggplot(df_sigs, aes(sequencing, cell_type_int, fill = mean_stress)) + 
  geom_tile(color = 'white') + scale_fill_viridis() + theme_minimal() + 
  ggtitle('stress_module1') + coord_fixed() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 7, hjust = 1), 
                                    axis.title.x = element_blank(), axis.title.y = element_blank()) 
g_mito <- ggplot(df_sigs, aes(sequencing, cell_type_int, fill = mean_mito)) + 
  geom_tile(color = 'white') + coord_fixed() + scale_fill_viridis() + theme_minimal() + ggtitle('mito_module1') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 7, hjust = 1), 
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank()) 
g_inf <- ggplot(df_sigs, aes(sequencing, cell_type_int, fill = mean_inf)) + 
  geom_tile(color = 'white') + scale_fill_viridis() + coord_fixed() + theme_minimal() + ggtitle('IFN_module1') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 7, hjust = 1), 
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank()) 
g_marsh<-ggplot(df_sigs, aes(sequencing, cell_type_int, fill = mean_marsh))+
  geom_tile(color = 'white')+
  scale_fill_viridis()+
  theme_minimal()+ ggtitle('stress_marsh1')+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 7, hjust = 1),axis.title.x=element_blank(),axis.title.y=element_blank())+coord_fixed()
g_vanhove<-ggplot(df_sigs, aes(sequencing, cell_type_int, fill = mean_vanhove))+
  geom_tile(color = 'white')+
  scale_fill_viridis()+
  theme_minimal()+ ggtitle('stress_vanhove1')+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 7, hjust = 1),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_fixed()
g_brink<-ggplot(df_sigs, aes(sequencing, cell_type_int, fill = mean_brink))+
  geom_tile(color = 'white')+
  scale_fill_viridis()+
  theme_minimal()+ ggtitle('stress_brink1')+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 7, hjust = 1),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_fixed()
g_denisenko<-ggplot(df_sigs, aes(sequencing, cell_type_int, fill = mean_denisenko))+
  geom_tile(color = 'white')+
  scale_fill_viridis()+
  theme_minimal()+ ggtitle('stress_denisenko1')+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 7, hjust = 1),axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.y = element_blank())+coord_fixed()


# Plots part 1
pdf(file = 'data/seq_comparison/sc_sn/integrated_sc_sn.pdf')
VlnPlot(seu, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'doublet_scores',
                          'stress_module1', 'IFN_module1'), pt.size = 0, group.by = 'sequencing', cols = colSCSN)
VlnPlot(seu, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'doublet_scores', 
                          'stress_module1', 'IFN_module1'), pt.size = 0, group.by = 'sequencing', cols = colSCSN) + 
  theme(legend.position = 'right')

DimPlot(seu, reduction = 'umap', label = F, group.by = 'sequencing', shuffle = T, raster = T, cols = colSCSN) + NoLegend()
DimPlot(seu, reduction = 'umap', label = F, group.by = 'sequencing', shuffle = T, raster = T, cols = colSCSN)

DimPlot(seu, reduction = 'umap', label = F, group.by = 'cell_type_main') + NoLegend()
DimPlot(seu, reduction = 'umap', label = T, group.by = 'cell_type_main', repel = T)

DimPlot(seu, reduction = 'umap', label = F, group.by = 'cell_type_int') + NoLegend()
DimPlot(seu, reduction = 'umap', label = T, group.by = 'cell_type_int', repel = T, label.size = 3)

DimPlot(seu, reduction = 'umap', label = F, group.by = 'cell_type_fine') + NoLegend()
DimPlot(seu, reduction = 'umap', label = T, group.by = 'cell_type_fine', repel = T, label.size = 3) + guides(color = guide_legend(ncol = 1), override.aes = list(size = 1)) + theme(legend.text = element_text(size = 5), 
                                                                                                                                                                                    legend.key.size = unit(0.2, 'cm'))

FeaturePlot(seu, features = 'proportion_scaled_cnv_avg') + scale_color_viridis(direction = -1)
FeaturePlot(seu, features = 'proportion_scaled_cnv_avg') + scale_color_viridis(direction = -1) + NoLegend()

VlnPlot(seu, features = 'proportion_scaled_cnv_avg', group.by = 'sequencing', pt.size = 0, cols = colSCSN)
VlnPlot(seu, features = 'proportion_scaled_cnv_avg', group.by = 'sequencing', pt.size = 0, cols = colSCSN) + NoLegend()

cti <- data.frame(seu$sequencing, seu$cell_type_main)
ggplot(cti, aes(x = seu.sequencing, fill = seu.cell_type_main)) + 
  geom_bar(colour = 'black', position = 'fill', size = 0.25) + ggtitle('cell_type_main') + 
  theme_classic() + xlab('Sequencing technique') + 
  ylab('Fraction (%)') + labs(fill = 'Cell type') + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

cti <- data.frame(seu$sequencing, seu$cell_type_int)
ggplot(cti, aes(x = seu.sequencing, fill = seu.cell_type_int)) + 
  geom_bar(colour = 'black', position = 'fill', size = 0.25) + ggtitle('cell_type_int') + 
  theme_classic() + xlab('Sequencing technique') + 
  ylab('Fraction (%)') + labs(fill = 'Cell type') + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

cti <- data.frame(seu$sequencing, seu$cell_type_fine)
ggplot(cti, aes(x = seu.sequencing, fill = seu.cell_type_fine)) + 
  geom_bar(colour = 'black', position = 'fill', size = 0.25) + ggtitle('cell_type_fine') + 
  theme_classic() + xlab('Sequencing technique') + 
  ylab('Fraction (%)') + labs(fill = 'Cell type') + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

VlnPlot(seu, features = c('stress_module1','mito_module1','IFN_module1','Ig_module1'),
        group.by = 'sequencing',pt.size = 0,cols = colSCSN)
VlnPlot(seu, features = c('stress_marsh1','stress_vanhove1','stress_brink1','stress_denisenko1'),
        group.by = 'sequencing',pt.size = 0,cols = colSCSN,ncol = 2)

VlnPlot(seu, features = c('stress_marsh1','stress_vanhove1','stress_brink1','stress_denisenko1'),
        group.by = 'cell_type_int',split.by = 'sequencing',pt.size = 0,cols = colSCSN,ncol = 2)+
  theme(legend.position = 'bottom')
VlnPlot(seu, features = c('stress_module1','mito_module1','IFN_module1','Ig_module1'),
        group.by = 'cell_type_int',split.by = 'sequencing',pt.size = 0,cols = colSCSN,ncol = 2)+
  theme(legend.position = 'bottom')

VlnPlot(seu, features = c('mito_module1','IFN_module1','stress_module1',
                          'stress_marsh1','stress_vanhove1','stress_brink1','stress_denisenko1'),
        group.by = 'cell_type_int',split.by = 'sequencing',pt.size = 0,cols = colSCSN,
        stack = T,flip = T)+
  theme(legend.position = 'bottom',axis.title.x = element_blank())

FeaturePlot(seu, features = 'proportion_scaled_cnv_avg') + scale_color_viridis(direction = -1)
VlnPlot(seu, features = 'proportion_scaled_cnv_avg', group.by = 'sequencing', pt.size = 0, cols = colSCSN)

print(g_stress + g_mito + g_inf)
print(g_marsh+g_vanhove+g_brink+g_denisenko)
dev.off()


#### in-silico separation of MBM05_sn ####
seu$batch <- ifelse(seu$orig.ident == 'MBM05_sn' & 
                      seu$cell_type_main %in% c('Myeloid cells', 'B/Plasma cells', 'T/NK cells'), 'MBM05_sn_CD45+', 
                    ifelse(seu$orig.ident == 'MBM05_sn' & 
                             seu$cell_type_main %notin% c('Myeloid cells', 'B/Plasma cells', 'T/NK cells'), 
                           'MBM05_sn_CD45-', seu$orig.ident))
seu$batch <- ifelse(seu$batch == 'MBM05_sc_CD45neg', 'MBM05_sc_CD45-', 
                    ifelse(seu$batch == 'MBM05_sc_CD45pos', 'MBM05_sc_CD45+', seu$batch))

# Plots part 2
pdf(file = 'data/seq_comparison/sc_sn/integrated_sc_sn_batch.pdf')
FeatureScatter(seu, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA', group.by = 'batch', 
               raster = T, cols = colSCSN2)

VlnPlot(seu, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), pt.size = 0, 
        group.by = 'batch', cols = colSCSN2)

DimPlot(seu, reduction = 'umap', label = F, group.by = 'batch', shuffle = T, raster = T, cols = colSCSN2) + NoLegend()
DimPlot(seu, reduction = 'umap', label = F, group.by = 'batch', shuffle = T, raster = T, cols = colSCSN2)
dev.off()


#### Raw files (pre-QC) #### 
# Raw MBM05_sn
patient.data <- Read10X(data.dir = paste0('data/MBM05_sn/filtered_feature_bc_matrix'))

# Initialize the Seurat object with the raw (non-normalized data)
patient_raw <- CreateSeuratObject(counts = patient.data, project = 'MBM05_sn', min.cells = 1, min.features = 1)

# Annotate MT genes
patient_raw[['percent.mt']] <- PercentageFeatureSet(patient_raw, pattern = '^MT-')

# Annotate pat and facs info
patient_raw[['patient']] <- 'MBM05_sn'
patient_raw[['sequencing']] <- 'Single nuclei'

patient_raw <- NormalizeData(patient_raw)
patient_raw <- FindVariableFeatures(patient_raw)
patient_raw <- ScaleData(patient_raw, features = rownames(patient_raw))
patient_raw <- RunPCA(patient_raw)
patient_raw <- RunUMAP(patient_raw, dims = 1:30)
patient_raw <- FindNeighbors(patient_raw, dims = 1:30)
patient_raw <- FindClusters(patient_raw, resolution = 0.2)

# Cell-type annotation using SingleR for in-silico separation of batches
patient_sce <- as.SingleCellExperiment(patient_raw)
bped <- BlueprintEncodeData()
pred_bped_main <- SingleR(test = patient_sce, ref = bped, labels = bped$label.main)
pruneScores(pred_bped_main)
patient_raw[['celltype_bped_main']] <- pred_bped_main$pruned.labels

DimPlot(patient_raw, reduction = 'umap', label = T, group.by = 'ident')
DimPlot(patient_raw, reduction = 'umap', label = T, group.by = 'celltype_bped_main')
FeaturePlot(patient_raw, features = 'PTPRC')

patient_raw$batch <- ifelse(patient_raw$RNA_snn_res.0.2 %in% c(2, 3, 5), 'raw_MBM05_sn_CD45+', 'raw_MBM05_sn_CD45-')
DimPlot(patient_raw, reduction = 'umap', label = T, group.by = 'batch')


# Raw MBM05_sc CD45neg
raw_neg <- Read10X(data.dir = paste0('data/MBM05_sc/CD45negGEXBI5/outs/filtered_feature_bc_matrix/'))
raw_neg <- CreateSeuratObject(counts = raw_neg, project = 'MBM05_sc_CD45-', min.cells = 1, min.features = 1)
raw_neg[['percent.mt']] <- PercentageFeatureSet(raw_neg, pattern = '^MT-')
raw_neg[['patient']] <- 'MBM05_sc'
raw_neg[['sequencing']] <- 'Single cells'
raw_neg[['batch']] <- 'raw_MBM05_sc_CD45-'

# Raw MBM05_sc CD45pos
raw_pos <- Read10X(data.dir = paste0('data/MBM05_sc/CD45posGEXBI5/outs/filtered_feature_bc_matrix/'))
raw_pos <- CreateSeuratObject(counts = raw_pos, project = 'MBM05_sc_CD45+', min.cells = 1, min.features = 1)
raw_pos[['percent.mt']] <- PercentageFeatureSet(raw_pos, pattern = '^MT-')
raw_pos[['patient']] <- 'MBM05_sc'
raw_pos[['sequencing']] <- 'Single cells'
raw_pos[['batch']] <- 'raw_MBM05_sc_CD45+'

# Merge MBM05_sc CD45neg, CD45pos and MBM05_sn
mer <- merge(patient_raw, y = c(raw_neg, raw_pos))
saveRDS(mer, file = paste0('data/seq_comparison/sc_sn/raw_merged_sc_sn.rds'))

# Merge merged object (MBM05_sc CD45neg, CD45pos and MBM05_sn) with re-integrated version from above
mer2 <- merge(seu, y = mer)
order_x <- c('raw_MBM05_sc_CD45-', 'raw_MBM05_sc_CD45+', 'raw_MBM05_sn_CD45-', 'raw_MBM05_sn_CD45+', 'MBM05_sc_CD45-', 'MBM05_sc_CD45+', 'MBM05_sn_CD45-', 'MBM05_sn_CD45+')
mer2$batch <- factor(x = mer2$batch, levels = order_x)

# Plots part 3
pdf(file = 'data/seq_comparison/sc_sn/plots_raw_merged_sc_sn_batch.pdf')
plts <- VlnPlot(mer2, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'), 
                pt.size = 0, group.by = 'batch', combine = F, cols = c(colSCSN2, colSCSN2))
plts <- lapply(plts, function(x) {
  x + theme(axis.title.x = element_blank(), 
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))
})
CombinePlots(plots = plts, legend = 'none', ncol = 3)
dev.off()
