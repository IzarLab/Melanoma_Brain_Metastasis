#!/usr/bin/env Rscript

### title: Boxplots and bar plots of cell-type frequencies 
### author: Jana Biermann, PhD

library(dplyr)
library(Seurat)
library(ggplot2)
library(gplots)
library(viridis)
library(ggpubr)
'%notin%' <- Negate('%in%')

colBP <- c('#A80D11', '#008DB8')
colSCSN <- c('#E1AC24', '#288F56')


#### Boxplots of cell-type frequencies ####
seu <- readRDS('data/MBPM/data_MBPM.rds')
seu <- subset(seu, sequencing == 'Single nuclei')
seu <- subset(seu, cell_type_int %notin% c('Low-quality cells', 'Doublets', 'Contamination', 'Not determined'))
seu$immune_status <- ifelse(seu$cell_type_main %in% c('Myeloid cells', 'T/NK cells', 'B/Plasma cells'),
                            'Immune', 'Non-immune')

df <- data.frame(patient = seu$patient, organ = seu$organ, cell_type_main = seu$cell_type_main, 
                 cell_type_int = seu$cell_type_int, cell_type_fine = seu$cell_type_fine, 
                 immune_status = seu$immune_status)

# cell_type_main
df_main = df %>%
    group_by(patient, cell_type_main, organ) %>%
    tally() %>%
    group_by(patient) %>%
    mutate(freq = n/sum(n))

# cell_type_int
df_int = df %>%
    group_by(patient, cell_type_int, organ) %>%
    tally() %>%
    group_by(patient) %>%
    mutate(freq = n/sum(n))

# cell_type_fine
df_fine = df %>%
    group_by(patient, cell_type_fine, organ) %>%
    tally() %>%
    group_by(patient) %>%
    mutate(freq = n/sum(n))


## Generate boxplots of cell-type frequencies
pdf(file = 'data/cell_type_assignment/plots_MBPM_sn_boxplot_freq.pdf', width = 10)
# cell_type_main
ggboxplot(df_main, x = 'cell_type_main', y = 'freq', color = 'organ', add = 'jitter', 
          order = sort(unique(df_main$cell_type_main))) + ylim(0, 1) + 
    stat_compare_means(aes(group = organ, label = sprintf('p = %5.3f', as.numeric(..p.format..))),
                       method = 'wilcox.test', size = 3) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle('cell_type_main') + color_palette(palette = colBP)

# cell_type_int
ggboxplot(df_int, x = 'cell_type_int', y = 'freq', color = 'organ', add = 'jitter', 
          order = sort(unique(df_int$cell_type_int))) + ylim(0, 1) + 
    stat_compare_means(aes(group = organ, label = sprintf('p = %5.3f', as.numeric(..p.format..))), 
                       method = 'wilcox.test', size = 3) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    ggtitle('cell_type_int') + color_palette(palette = colBP)

# cell_type_fine
ggboxplot(df_fine, x = 'cell_type_fine', y = 'freq', color = 'organ', add = 'jitter', 
          order = sort(unique(df_fine$cell_type_fine))) + ylim(0, 1) + 
    stat_compare_means(aes(group = organ, label = sprintf('p = %5.3f', as.numeric(..p.format..))),
                       method = 'wilcox.test', size = 1.8) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    ggtitle('cell_type_fine') + color_palette(palette = colBP)
dev.off()


#### Bar plots of cell-type frequencies (single nuclei) ####
seu <- readRDS('data/MBPM/data_MBPM.rds')
seu <- subset(seu, sequencing == 'Single nuclei')
seu <- subset(seu, cell_type_int %notin% c('Low-quality cells', 'Doublets', 'Contamination', 'Not determined'))
seu$immune_status <- ifelse(seu$cell_type_main %in% c('Myeloid cells', 'T/NK cells', 'B/Plasma cells'), 'Immune', 'Non-immune')

df <- data.frame(patient = seu$patient, organ = seu$organ, sequencing = seu$sequencing, 
                 cell_type_main = seu$cell_type_main, cell_type_int = seu$cell_type_int, 
                 cell_type_fine = seu$cell_type_fine, immune_status = seu$immune_status)

immune <- df %>% filter(immune_status == 'Immune')
tcell <- df %>% filter(cell_type_main == 'T/NK cells')
myeloid <- df %>% filter(cell_type_main == 'Myeloid cells')

pdf(file = 'data/cell_type_assignment/plots_MBPM_sn_barplots.pdf')
## All sn cell types 
# Stacked bar plots by patient
ggplot(df, aes(x = patient, fill = cell_type_main)) + 
    geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('Single nuclei') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(df, aes(x = patient, fill = cell_type_int)) + 
    geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('Single nuclei') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Stacked bar plots by tissue of origin
ggplot(df, aes(x = organ, fill = cell_type_main)) + 
    geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('Single nuclei') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), plot.margin = margin(0, 8, 0, 0, 'cm'))

ggplot(df, aes(x = organ, fill = cell_type_int)) + 
    geom_bar(position = 'fill', size = 0.25, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('Single nuclei') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), plot.margin = margin(0, 8, 0, 0, 'cm'))

## All sn immune cells 
# Stacked bar plots by patient
ggplot(immune, aes(x = patient, fill = cell_type_main)) + 
    geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('Single nuclei (only immune cells)') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(immune, aes(x = patient, fill = cell_type_int)) + 
    geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('Single nuclei (only immune cells)') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(immune, aes(x = patient, fill = cell_type_fine)) + 
    geom_bar(position = 'fill', size = 0.25, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('Single nuclei (only immune cells)') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Stacked bar plots by tissue of origin
ggplot(immune, aes(x = organ, fill = cell_type_main)) + 
    geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('Single nuclei (only immune cells)') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), plot.margin = margin(0, 8, 0, 0, 'cm'))

ggplot(immune, aes(x = organ, fill = cell_type_int)) + 
    geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('Single nuclei (only immune cells)') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), plot.margin = margin(0, 8, 0, 0, 'cm'))

ggplot(immune, aes(x = organ, fill = cell_type_fine)) + 
    geom_bar(position = 'fill', size = 0.25, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('Single nuclei (only immune cells)') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), plot.margin = margin(0, 8, 0, 0, 'cm'))

## Only sn T cells 
# Stacked bar plots by patient
ggplot(tcell, aes(x = patient, fill = cell_type_int)) + 
    geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('Single nuclei (only T/NK cells)') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(tcell, aes(x = patient, fill = cell_type_fine)) + 
    geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('Single nuclei (only T/NK cells)') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Stacked bar plots by tissue of origin
ggplot(tcell, aes(x = organ, fill = cell_type_int)) + 
    geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('Single nuclei (only T/NK cells)') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), plot.margin = margin(0, 8, 0, 0, 'cm'))

ggplot(tcell, aes(x = organ, fill = cell_type_fine)) + 
    geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('Single nuclei (only T/NK cells)') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), plot.margin = margin(0, 8, 0, 0, 'cm'))

## Only sn myeloid cells 
# Stacked bar plots by patient
ggplot(myeloid, aes(x = patient, fill = cell_type_int)) + 
    geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('Single nuclei (only myeloid cells)') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(myeloid, aes(x = patient, fill = cell_type_fine)) + 
    geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('Single nuclei (only myeloid cells)') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Stacked bar plots by tissue of origin
ggplot(myeloid, aes(x = organ, fill = cell_type_int)) + 
    geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('Single nuclei (only myeloid cells)') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), plot.margin = margin(0, 8, 0, 0, 'cm'))

ggplot(myeloid, aes(x = organ, fill = cell_type_fine)) + 
    geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('Single nuclei (only myeloid cells)') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), plot.margin = margin(0, 8, 0, 0, 'cm'))
dev.off()


#### Bar plots of cell-type frequencies (MBM) ####
seu <- readRDS('data/MBPM/data_MBPM.rds')
seu <- subset(seu, organ == 'Brain')
seu <- subset(seu, cell_type_int %notin% c('Low-quality cells', 'Doublets', 'Contamination', 'Not determined'))
seu$immune_status <- ifelse(seu$cell_type_main %in% c('Myeloid cells', 'T/NK cells', 'B/Plasma cells'), 'Immune', 'Non-immune')

df <- data.frame(patient = seu$patient, organ = seu$organ, sequencing = seu$sequencing, 
                 cell_type_main = seu$cell_type_main, cell_type_int = seu$cell_type_int, 
                 cell_type_fine = seu$cell_type_fine, immune_status = seu$immune_status)

nonimmune <- df %>% filter(immune_status == 'Non-immune')
immune <- df %>% filter(immune_status == 'Immune')
tcell <- df %>% filter(cell_type_main == 'T/NK cells')
myeloid <- df %>% filter(cell_type_main == 'Myeloid cells')

pdf(file = 'data/cell_type_assignment/plots_MBM_scn_barplots.pdf')
## All immune cells 
# Stacked bar plots by sequencing
ggplot(immune, aes(x = sequencing, fill = cell_type_main)) + 
    geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('MBM') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(immune, aes(x = sequencing, fill = cell_type_int)) + 
    geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('MBM') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(immune, aes(x = sequencing, fill = cell_type_fine)) + 
    geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('MBM') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

## All non-immune cells 
# Stacked bar plots by sequencing
ggplot(nonimmune, aes(x = sequencing, fill = cell_type_main)) + 
    geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('MBM') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(nonimmune, aes(x = sequencing, fill = cell_type_int)) + 
    geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('MBM') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(nonimmune, aes(x = sequencing, fill = cell_type_fine)) + 
    geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('MBM') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

## Only myeloid cells 
# Stacked bar plots by sequencing
ggplot(myeloid, aes(x = sequencing, fill = cell_type_int)) + 
    geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('MBM') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(myeloid, aes(x = sequencing, fill = cell_type_fine)) + 
    geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('MBM') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

## Only T cells 
# Stacked bar plots by sequencing
ggplot(tcell, aes(x = sequencing, fill = cell_type_int)) + 
    geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('MBM') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(tcell, aes(x = sequencing, fill = cell_type_fine)) + 
    geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('MBM') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()
#### FACS plot ####
facs <- read.csv('data/MBPM/facs_output.csv')
facs_melt <- reshape2::melt(facs, variable.name = 'CD45')

order_facs <- c('CD45pos', 'CD45neg')
facs_melt$CD45 <- factor(x = facs_melt$CD45, levels = order_facs)

pdf('data/MBPM/plots_facs.pdf')
ggplot(facs_melt, aes(x = patient, y = value, fill = CD45)) + 
    geom_bar(position = 'fill', stat = 'identity', size = 0.3, color = 'black') + 
    theme_classic() + ylab('CD45 fraction') + ggtitle('FACS') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    scale_fill_manual(values = c('#B25D91', '#F4B95A'))
dev.off()
