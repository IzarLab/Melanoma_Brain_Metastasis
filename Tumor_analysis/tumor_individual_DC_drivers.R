#!/usr/bin/env Rscript

### title: Individual diffusion component (D.C.) analysis of tumor cells
### and overlap of drivers in DC1-3 between samples 
### author: Jana Biermann, PhD

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
library(SingleCellExperiment)
library(cowplot)
library(ggpubr)

colBP <- c('#A80D11', '#008DB8')
colSCSN <- c('#E1AC24', '#288F56')

#### Read-in integrated object and subset to non-cycling tumor cells ####
seu <- readRDS('data/MBPM/data_MBPM_scn.rds')
seu <- subset(seu, sequencing =='Single nuclei' & 
                cell_type_fine == 'Tumor cells' & 
                cell_cycle == 'non-cycling')
DefaultAssay(seu) <- 'RNA'
seu <- DietSeurat(seu, assays = 'RNA', counts = T)
# Remove RP and MT genes
seu <- seu[grep("^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA|^MT-", rownames(seu), value = T, invert = T), ]

# Split object and remove samples with <200 cells
obj.list <- SplitObject(seu, split.by = 'patient')
for (i in 1:length(obj.list)) {
  if (dim(obj.list[[i]]@assays$RNA@counts)[2] < 200) {
    obj.list[i] <- NA
  }
}
# Remove NAs
obj.list <- obj.list[!is.na(obj.list)]


#### Loop through samples and perform DC analysis ####
for (seu in obj.list) {
  pat <- seu$patient[1]
  seu <- NormalizeData(seu) %>% ScaleData()
  
  # signatures
  sigs<-read.csv('~/signatures/melanoma_sigs.csv',na.strings = c('','NA'))
  for (c in 1:ncol(sigs)) {
    sig <- as.character(na.omit(sigs[, c]))
    seu <- AddModuleScore(object = seu, features = list(sig), 
                          name = colnames(sigs[c]), assay = 'RNA', search = F)
  }
  
  # Run D.C. analysis
  es <- as.ExpressionSet(as.data.frame(t(seu@assays$RNA@data)))
  es@phenoData@data <- seu@meta.data
  set.seed(1)
  dm <- DiffusionMap(es, verbose = T, n_pcs = 20)
  ifelse(!dir.exists(file.path(paste0('data/tumor/diffusion/single_noMTRP/', pat))), 
         dir.create(file.path(paste0('data/tumor/diffusion/single_noMTRP/', pat)), recursive = T), FALSE)
  #saveRDS(dm, paste0('data/tumor/diffusion/single_noMTRP/', pat, '/data_dm_', pat, '.rds'))
  
  pdf(paste0('data/tumor/diffusion/single_noMTRP/', pat,'/plots_dm_',pat,'.pdf'))
  par(mar=c(5.1, 4.1, 4.1, 6.5),xpd=TRUE)
  print(plot(dm,col=es@phenoData@data$MITF_sig1,main='MITF_sig1',pch=20))
  colorlegend(es@phenoData@data$MITF_sig1,viridis(length(es@phenoData@data$MITF_sig1)),posx = c(0.88, 0.9),posy = c(0,0.65))
  
  print(plot(dm,col=es@phenoData@data$AXL_sig1,main='AXL_sig1',pch=20))
  colorlegend(es@phenoData@data$AXL_sig1,viridis(length(es@phenoData@data$AXL_sig1)),posx = c(0.88, 0.9),posy = c(0,0.65))
  
  print(plot_grid(plot(dm, 1:2, col_by='MITF_sig1')+ theme_classic() + ggtitle('MITF_sig1')+ theme(legend.position = "none"),
                  plot(dm, 1:2, col_by='AXL_sig1')+ theme_classic() + ggtitle('AXL_sig1')+ theme(legend.position = "none"),
                  plot(dm, 1:2, col_by='S.Score')+ theme_classic() + ggtitle('S.Score')+ theme(legend.position = "none"),
                  plot(dm, 1:2, col_by='G2M.Score')+ theme_classic() + ggtitle('G2M.Score')+ theme(legend.position = "none"),
                  labels = c("A", "B", "C","D"),ncol = 2, nrow = 2))
  dev.off()
  
  # Add DC info to Seurat and save loadings table
  seu <- AddMetaData(seu, dm@eigenvectors, col.name = colnames(dm@eigenvectors))
  seu[['DC']] <- CreateDimReducObject(embeddings = dm@eigenvectors, key = 'DC_', assay = 'RNA')
  seu <- ProjectDim(seu, reduction = 'DC')
  write.csv(as.data.frame(seu@reductions$DC@feature.loadings.projected), 
            paste0('data/tumor/diffusion/single_noMTRP/', pat, '/table_dm_loadings_', pat, '.csv'))
}



#### Combine genes from individual patients ####
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
write.csv(res, 'data/tumor/diffusion/single/table_single_DC1-3_drivers_noMTRP.csv', row.names = F)



#### Count and check overlap ####
res<-read.csv('data/tumor/diffusion/single/table_single_DC1-3_drivers_noMTRP.csv')

res_mbm <- res[, grep('MBM', colnames(res))]
res_mbm_long <- tidyr::gather(res_mbm,key='patient',value = 'driver')
mbm_counts <- res_mbm_long %>% group_by(driver) %>% tally() %>% na.omit()
mbm_counts$percentage<- mbm_counts$n / length(grep('MBM', colnames(res)))
length(which(mbm_counts$percentage > 0.3))

res_mpm <- res[, grep('MPM', colnames(res))]
res_mpm_long <- tidyr::gather(res_mpm,key='patient',value = 'driver')
mpm_counts <- res_mpm_long %>% group_by(driver) %>% tally() %>% na.omit()
mpm_counts$percentage<- mpm_counts$n / length(grep('MPM', colnames(res)))
length(which(mpm_counts$percentage > 0.3))

intersect(mbm_counts$driver[which(mbm_counts$percentage > 0.3)],
          mpm_counts$driver[which(mpm_counts$percentage > 0.3)])

length(intersect(mbm_counts$driver[which(mbm_counts$percentage > 0.3)],
          mpm_counts$driver[which(mpm_counts$percentage > 0.3)]))

# unique to MBM
setdiff(mbm_counts$driver[which(mbm_counts$percentage > 0.3)],
          mpm_counts$driver[which(mpm_counts$percentage > 0.3)])

length(setdiff(mbm_counts$driver[which(mbm_counts$percentage > 0.3)],
        mpm_counts$driver[which(mpm_counts$percentage > 0.3)]))

# unique to MPM
setdiff(mpm_counts$driver[which(mpm_counts$percentage > 0.3)],
        mbm_counts$driver[which(mbm_counts$percentage > 0.3)])

length(setdiff(mpm_counts$driver[which(mpm_counts$percentage > 0.3)],
               mbm_counts$driver[which(mbm_counts$percentage > 0.3)]))


#### Top shared genes + enrichment ####
library(hypeR)
HALLMARK <- msigdb_gsets(species="Homo sapiens", category="H")
genesets_w <- as.list(as.data.frame(read.csv('~/signatures/wouters_mel_sigs.csv', na.strings = c('', NA))))

# Plot MBM bar plots and pathway enrichments with genes detected in >30% of samples
plot_mbm_30 <- mbm_counts[which(mbm_counts$percentage > 0.3),] %>% arrange(desc(percentage))
write.csv(plot_mbm_30, 'data/tumor/diffusion/single/table_MBM_drivers_noMTRP_0.3.csv', row.names = F)

pdf('data/tumor/diffusion/single/plot_MBM_single_DC1-3_drivers_barplot_noMTRP.pdf', 
    height = 14)
ggplot(plot_mbm_30, aes(x = percentage, y = reorder(driver,percentage))) + 
  ggtitle('MBM - Top overlapping drivers of DC1-3') + 
  geom_bar(stat = 'identity', color = 'white', fill = colBP[1], size = 0.25) + 
  theme_pubr() + 
  theme(legend.position = 'none') + 
  xlab('Percentage of overlap') + 
  ylab('Genes') + 
  scale_x_continuous(expand = c(0, 0))

ggplot(plot_mbm_30, aes(x = n, y = reorder(driver,n))) + 
  ggtitle('MBM - Top overlapping drivers of DC1-3') + 
  geom_bar(stat = 'identity', color = 'white', fill = colBP[1], size = 0.25) + 
  theme_pubr() + 
  theme(legend.position = 'none') + 
  xlab('# samples') + 
  ylab('Genes') +
  scale_x_continuous(expand = c(0, 0),breaks = c(seq(0,max(plot_mbm_30$n)+2,2)))
dev.off()

# MBM pathways with >25% shared genes
pw_mbm_25 <- mbm_counts[which(mbm_counts$percentage > 0.25),] 

pdf('data/tumor/diffusion/single/plot_MBM_single_DC1-3_drivers_pathways_noMTRP.pdf')
hyp_dots(hypeR(as.character(pw_mbm_25$driver), genesets = genesets_w), 
         title = paste0('Wouters et al. pathway selection'), abrv = 80) + 
  theme_bw() + 
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_text(color = 'black'), 
        axis.text.x = element_text(color = 'black'), 
        panel.grid.minor = element_blank()) + 
  scale_color_gradient(high = 'black', low = '#C51B8A')
hyp_dots(hypeR(as.character(pw_mbm_25$driver), genesets = HALLMARK), 
         title = paste0('Hallmarks'),abrv = 80) + 
  theme_bw() + 
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_text(color = 'black'), 
        axis.text.x = element_text(color = 'black'), 
        panel.grid.minor = element_blank()) + 
  scale_color_gradient(high = 'black', low = '#C51B8A')
dev.off()


# Plot MPM bar plots with genes detected in >30% of samples
plot_mpm_30 <- mpm_counts[which(mpm_counts$percentage > 0.3),] %>% arrange(desc(percentage))
write.csv(plot_mpm_30, 'data/tumor/diffusion/single/table_MPM_drivers_noMTRP_0.3.csv', row.names = F)

pdf('data/tumor/diffusion/single/plot_MPM_single_DC1-3_drivers_barplot_noMTRP.pdf', 
    height = 14)
ggplot(plot_mpm_30, aes(x = percentage, y = reorder(driver,percentage))) + 
  ggtitle('MPM - Top overlapping drivers of DC1-3') + 
  geom_bar(stat = 'identity', color = 'white', fill = colBP[2], size = 0.25) + 
  theme_pubr() + 
  theme(legend.position = 'none') + 
  xlab('Percentage of overlap') + 
  ylab('Genes') + 
  scale_x_continuous(expand = c(0, 0))

ggplot(plot_mpm_30, aes(x = n, y = reorder(driver,n))) + 
  ggtitle('MPM - Top overlapping drivers of DC1-3') + 
  geom_bar(stat = 'identity', color = 'white', fill = colBP[2], size = 0.25) + 
  theme_pubr() + 
  theme(legend.position = 'none') + 
  xlab('# samples') + 
  ylab('Genes') +
  scale_x_continuous(expand = c(0, 0),breaks = c(seq(0,max(plot_mpm_30$n)+2,2)))
dev.off()

# MPM pathways with >25% shared genes
pw_mpm_25 <- mpm_counts[which(mpm_counts$percentage > 0.25),] 

pdf('data/tumor/diffusion/single/plot_MPM_single_DC1-3_drivers_pathways_noMTRP.pdf')
hyp_dots(hypeR(as.character(pw_mpm_25$driver), genesets = genesets_w), 
         title = paste0('Wouters et al. pathway selection'), abrv = 80) + 
  theme_bw() + 
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_text(color = 'black'), 
        axis.text.x = element_text(color = 'black'), 
        panel.grid.minor = element_blank()) + 
  scale_color_gradient(high = 'black', low = '#C51B8A')

hyp_dots(hypeR(as.character(pw_mpm_25$driver), genesets = HALLMARK), 
         title = paste0('Hallmarks'),abrv = 80) + 
  theme_bw() + 
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_text(color = 'black'), 
        axis.text.x = element_text(color = 'black'), 
        panel.grid.minor = element_blank()) + 
  scale_color_gradient(high = 'black', low = '#C51B8A')
dev.off()

