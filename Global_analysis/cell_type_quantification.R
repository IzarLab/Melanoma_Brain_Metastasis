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
colBP_scn <- c('#762A83','#A80D11', '#008DB8')
colSCSN <- c('#E1AC24', '#288F56')


#### Boxplots of cell-type frequencies ####
seu <- readRDS('data/MBPM/data_MBPM_scn.rds')
seu <- subset(seu, sequencing == 'Single nuclei')
seu <- subset(seu, cell_type_int %notin% c('Low-quality cells', 'Doublets', 'Contamination',
                                           'Undetermined', 'Cycling cells'))
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
pdf(file = 'data/cell_type_DEG/plots_MBPM_sn_boxplot_freq.pdf', width = 10)
# cell_type_main
ggboxplot(df_main, x = 'cell_type_main', y = 'freq', color = 'organ', add = 'jitter', 
          order = sort(unique(df_main$cell_type_main))) + ylim(0, 1) + 
    stat_compare_means(aes(group = organ),label = 'p.format',
                       method = 'wilcox.test', size = 3) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle('cell_type_main') + color_palette(palette = colBP)

# cell_type_int
ggboxplot(df_int, x = 'cell_type_int', y = 'freq', color = 'organ', add = 'jitter', 
          order = sort(unique(df_int$cell_type_int))) + ylim(0, 1) + 
    stat_compare_means(aes(group = organ),label = 'p.format',
                       method = 'wilcox.test', size = 2) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    ggtitle('cell_type_int') + color_palette(palette = colBP)

# cell_type_fine
ggboxplot(df_fine, x = 'cell_type_fine', y = 'freq', color = 'organ', add = 'jitter', 
          order = sort(unique(df_fine$cell_type_fine))) + ylim(0, 1) + 
    stat_compare_means(aes(group = organ),label = 'p.format',
                       method = 'wilcox.test', size = 1.8) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    ggtitle('cell_type_fine') + color_palette(palette = colBP)
dev.off()


#### Boxplots of cell-type frequencies MBPM_scn immune ####
seu <- readRDS('data/MBPM/data_MBPM_scn.rds')
seu <- subset(seu, cell_type_int %notin% c('Low-quality cells', 'Doublets', 'Contamination',
                                           'Undetermined','Cycling cells'))
seu$immune_status <- ifelse(seu$cell_type_main %in% c('Myeloid cells', 'T/NK cells', 'B/Plasma cells'),
                            'Immune', 'Non-immune')
seu <- subset(seu, immune_status =='Immune')

seu$organ_seq<-ifelse(seu$organ=='Brain' & seu$sequencing =='Single cells', 'MBM_sc',
                      ifelse(seu$organ=='Brain' & seu$sequencing =='Single nuclei', 'MBM_sn',
                             ifelse(seu$organ=='Peripheral' & seu$sequencing =='Single nuclei', 
                                    'MPM_sn',NA)))

df <- data.frame(patient = seu$patient, 
                 organ = seu$organ, 
                 sequencing = seu$sequencing, 
                 cell_type_main = seu$cell_type_main, 
                 cell_type_int = seu$cell_type_int, 
                 cell_type_fine = seu$cell_type_fine, 
                 organ_seq = seu$organ_seq)

# cell_type_main
df_main = df %>%
    group_by(patient, cell_type_main, organ_seq) %>%
    tally() %>%
    group_by(patient) %>%
    mutate(freq = n/sum(n))

df_main_organ = df %>%
    filter(sequencing == 'Single nuclei') %>%
    group_by(patient, cell_type_main, organ) %>%
    tally() %>%
    group_by(patient) %>%
    mutate(freq = n/sum(n))

df_main_seq = df %>%
    filter(organ == 'Brain') %>%
    group_by(patient, cell_type_main, sequencing) %>%
    tally() %>%
    group_by(patient) %>%
    mutate(freq = n/sum(n))


# cell_type_int
df_int = df %>%
    group_by(patient, cell_type_int, organ_seq) %>%
    tally() %>%
    group_by(patient) %>%
    mutate(freq = n/sum(n))

df_int_organ = df %>%
    filter(sequencing == 'Single nuclei') %>%
    group_by(patient, cell_type_int, organ) %>%
    tally() %>%
    group_by(patient) %>%
    mutate(freq = n/sum(n))

df_int_seq = df %>%
    filter(organ == 'Brain') %>%
    group_by(patient, cell_type_int, sequencing) %>%
    tally() %>%
    group_by(patient) %>%
    mutate(freq = n/sum(n))


# cell_type_fine
df_fine = df %>%
    group_by(patient, cell_type_fine, organ_seq) %>%
    tally() %>%
    group_by(patient) %>%
    mutate(freq = n/sum(n))

df_fine_organ = df %>%
    filter(sequencing == 'Single nuclei') %>%
    group_by(patient, cell_type_fine, organ) %>%
    tally() %>%
    group_by(patient) %>%
    mutate(freq = n/sum(n))

df_fine_seq = df %>%
    filter(organ == 'Brain') %>%
    group_by(patient, cell_type_fine, sequencing) %>%
    tally() %>%
    group_by(patient) %>%
    mutate(freq = n/sum(n))



## Generate boxplots of cell-type frequencies
pdf(file = 'data/cell_type_DEG/plots_MBPM_scn_boxplot_freq_immune.pdf', width = 10)
# cell_type_main
ggboxplot(df_main, x = 'cell_type_main', y = 'freq', color = 'organ_seq', add = 'jitter', 
          order = sort(unique(df_main$cell_type_main))) + ylim(0, 1) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle('cell_type_main') + color_palette(palette = colBP_scn)

ggboxplot(df_main_organ, x = 'cell_type_main', y = 'freq', color = 'organ', add = 'jitter', 
          order = sort(unique(df_main_organ$cell_type_main))) + ylim(0, 1) + 
    stat_compare_means(aes(group = organ),label = 'p.format',
                       method = 'wilcox.test', size = 3) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle('cell_type_main Single nuclei') + color_palette(palette = colBP)

ggboxplot(df_main_seq, x = 'cell_type_main', y = 'freq', color = 'sequencing', add = 'jitter', 
          order = sort(unique(df_main_seq$cell_type_main))) + ylim(0, 1) + 
    stat_compare_means(aes(group = sequencing),label = 'p.format',
                       method = 'wilcox.test', size = 3) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle('cell_type_main MBM') + color_palette(palette = colSCSN)

# cell_type_int
ggboxplot(df_int, x = 'cell_type_int', y = 'freq', color = 'organ_seq', add = 'jitter', 
          order = sort(unique(df_int$cell_type_int))) + ylim(0, 1) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    ggtitle('cell_type_int') + color_palette(palette = colBP_scn)

ggboxplot(df_int_organ, x = 'cell_type_int', y = 'freq', color = 'organ', add = 'jitter', 
          order = sort(unique(df_int_organ$cell_type_int))) + ylim(0, 1) + 
    stat_compare_means(aes(group = organ),label = 'p.format',
                       method = 'wilcox.test', size = 3) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle('cell_type_int Single nuclei') + color_palette(palette = colBP)

ggboxplot(df_int_seq, x = 'cell_type_int', y = 'freq', color = 'sequencing', add = 'jitter', 
          order = sort(unique(df_int_seq$cell_type_int))) + ylim(0, 1) + 
    stat_compare_means(aes(group = sequencing),label = 'p.format',
                       method = 'wilcox.test', size = 3) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle('cell_type_int MBM') + color_palette(palette = colSCSN)


# cell_type_fine
ggboxplot(df_fine, x = 'cell_type_fine', y = 'freq', color = 'organ_seq', add = 'jitter', 
          order = sort(unique(df_fine$cell_type_fine))) + ylim(0, 1) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    ggtitle('cell_type_fine') + color_palette(palette = colBP_scn)

ggboxplot(df_fine_organ, x = 'cell_type_fine', y = 'freq', color = 'organ', add = 'jitter', 
          order = sort(unique(df_fine_organ$cell_type_fine))) + ylim(0, 1) + 
    stat_compare_means(aes(group = organ),label = 'p.format',
                       method = 'wilcox.test', size = 2) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle('cell_type_fine Single nuclei') + color_palette(palette = colBP)

ggboxplot(df_fine_seq, x = 'cell_type_fine', y = 'freq', color = 'sequencing', add = 'jitter', 
          order = sort(unique(df_fine_seq$cell_type_fine))) + ylim(0, 1) + 
    stat_compare_means(aes(group = sequencing),label = 'p.format',
                       method = 'wilcox.test', size = 1.8) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle('cell_type_fine MBM') + color_palette(palette = colSCSN)

dev.off()



#### Boxplots of cell-type frequencies MBPM_scn nonimmune ####
seu <- readRDS('data/MBPM/data_MBPM_scn.rds')
seu <- subset(seu, cell_type_int %notin% c('Low-quality cells', 'Doublets', 'Contamination',
                                           'Undetermined','Cycling cells'))
seu$immune_status <- ifelse(seu$cell_type_main %in% c('Myeloid cells', 'T/NK cells', 'B/Plasma cells'),
                            'Immune', 'Non-immune')
seu <- subset(seu, immune_status =='Non-immune')

seu$organ_seq<-ifelse(seu$organ=='Brain' & seu$sequencing =='Single cells', 'MBM_sc',
                      ifelse(seu$organ=='Brain' & seu$sequencing =='Single nuclei', 'MBM_sn',
                             ifelse(seu$organ=='Peripheral' & seu$sequencing =='Single nuclei', 
                                    'MPM_sn',NA)))

df <- data.frame(patient = seu$patient, 
                 organ = seu$organ, 
                 sequencing = seu$sequencing, 
                 cell_type_main = seu$cell_type_main, 
                 cell_type_int = seu$cell_type_int, 
                 cell_type_fine = seu$cell_type_fine, 
                 organ_seq = seu$organ_seq)

# cell_type_main
df_main = df %>%
    group_by(patient, cell_type_main, organ_seq) %>%
    tally() %>%
    group_by(patient) %>%
    mutate(freq = n/sum(n))

df_main_organ = df %>%
    filter(sequencing == 'Single nuclei') %>%
    group_by(patient, cell_type_main, organ) %>%
    tally() %>%
    group_by(patient) %>%
    mutate(freq = n/sum(n))

df_main_seq = df %>%
    filter(organ == 'Brain') %>%
    group_by(patient, cell_type_main, sequencing) %>%
    tally() %>%
    group_by(patient) %>%
    mutate(freq = n/sum(n))


# cell_type_int
df_int = df %>%
    group_by(patient, cell_type_int, organ_seq) %>%
    tally() %>%
    group_by(patient) %>%
    mutate(freq = n/sum(n))

df_int_organ = df %>%
    filter(sequencing == 'Single nuclei') %>%
    group_by(patient, cell_type_int, organ) %>%
    tally() %>%
    group_by(patient) %>%
    mutate(freq = n/sum(n))

df_int_seq = df %>%
    filter(organ == 'Brain') %>%
    group_by(patient, cell_type_int, sequencing) %>%
    tally() %>%
    group_by(patient) %>%
    mutate(freq = n/sum(n))


# cell_type_fine
df_fine = df %>%
    group_by(patient, cell_type_fine, organ_seq) %>%
    tally() %>%
    group_by(patient) %>%
    mutate(freq = n/sum(n))

df_fine_organ = df %>%
    filter(sequencing == 'Single nuclei') %>%
    group_by(patient, cell_type_fine, organ) %>%
    tally() %>%
    group_by(patient) %>%
    mutate(freq = n/sum(n))

df_fine_seq = df %>%
    filter(organ == 'Brain') %>%
    group_by(patient, cell_type_fine, sequencing) %>%
    tally() %>%
    group_by(patient) %>%
    mutate(freq = n/sum(n))



## Generate boxplots of cell-type frequencies
pdf(file = 'data/cell_type_DEG/plots_MBPM_scn_boxplot_freq_nonimmune.pdf', width = 10)
# cell_type_main
ggboxplot(df_main, x = 'cell_type_main', y = 'freq', color = 'organ_seq', add = 'jitter', 
          order = sort(unique(df_main$cell_type_main))) + ylim(0, 1) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle('cell_type_main') + color_palette(palette = colBP_scn)

ggboxplot(df_main_organ, x = 'cell_type_main', y = 'freq', color = 'organ', add = 'jitter', 
          order = sort(unique(df_main_organ$cell_type_main))) + ylim(0, 1) + 
    stat_compare_means(aes(group = organ),label = 'p.format',
                       method = 'wilcox.test', size = 3) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle('cell_type_main Single nuclei') + color_palette(palette = colBP)

ggboxplot(df_main_seq, x = 'cell_type_main', y = 'freq', color = 'sequencing', add = 'jitter', 
          order = sort(unique(df_main_seq$cell_type_main))) + ylim(0, 1) + 
    stat_compare_means(aes(group = sequencing),label = 'p.format',
                       method = 'wilcox.test', size = 3) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle('cell_type_main MBM') + color_palette(palette = colSCSN)

# cell_type_int
ggboxplot(df_int, x = 'cell_type_int', y = 'freq', color = 'organ_seq', add = 'jitter', 
          order = sort(unique(df_int$cell_type_int))) + ylim(0, 1) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    ggtitle('cell_type_int') + color_palette(palette = colBP_scn)

ggboxplot(df_int_organ, x = 'cell_type_int', y = 'freq', color = 'organ', add = 'jitter', 
          order = sort(unique(df_int_organ$cell_type_int))) + ylim(0, 1) + 
    stat_compare_means(aes(group = organ),label = 'p.format',
                       method = 'wilcox.test', size = 3) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle('cell_type_int Single nuclei') + color_palette(palette = colBP)

ggboxplot(df_int_seq, x = 'cell_type_int', y = 'freq', color = 'sequencing', add = 'jitter', 
          order = sort(unique(df_int_seq$cell_type_int))) + ylim(0, 1) + 
    stat_compare_means(aes(group = sequencing),label = 'p.format',
                       method = 'wilcox.test', size = 3) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle('cell_type_int MBM') + color_palette(palette = colSCSN)


# cell_type_fine
ggboxplot(df_fine, x = 'cell_type_fine', y = 'freq', color = 'organ_seq', add = 'jitter', 
          order = sort(unique(df_fine$cell_type_fine))) + ylim(0, 1) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    ggtitle('cell_type_fine') + color_palette(palette = colBP_scn)

ggboxplot(df_fine_organ, x = 'cell_type_fine', y = 'freq', color = 'organ', add = 'jitter', 
          order = sort(unique(df_fine_organ$cell_type_fine))) + ylim(0, 1) + 
    stat_compare_means(aes(group = organ),label = 'p.format',
                       method = 'wilcox.test', size = 2) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle('cell_type_fine Single nuclei') + color_palette(palette = colBP)

ggboxplot(df_fine_seq, x = 'cell_type_fine', y = 'freq', color = 'sequencing', add = 'jitter', 
          order = sort(unique(df_fine_seq$cell_type_fine))) + ylim(0, 1) + 
    stat_compare_means(aes(group = sequencing),label = 'p.format',
                       method = 'wilcox.test', size = 1.8) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle('cell_type_fine MBM') + color_palette(palette = colSCSN)

dev.off()


#### Bar plots of cell-type frequencies (single nuclei) ####
seu <- readRDS('data/MBPM/data_MBPM_scn.rds')
seu <- subset(seu, sequencing == 'Single nuclei')
seu <- subset(seu, cell_type_int %notin% c('Low-quality cells', 'Doublets', 'Contamination', 
                                           'Undetermined','Cycling cells'))
seu$immune_status <- ifelse(seu$cell_type_main %in% c('Myeloid cells', 'T/NK cells', 
                                                      'B/Plasma cells'), 
                            'Immune', 'Non-immune')

df <- data.frame(patient = seu$patient, organ = seu$organ, sequencing = seu$sequencing, 
                 cell_type_main = seu$cell_type_main, cell_type_int = seu$cell_type_int, 
                 cell_type_fine = seu$cell_type_fine, immune_status = seu$immune_status)

immune <- df %>% filter(immune_status == 'Immune')
tcell <- df %>% filter(cell_type_main == 'T/NK cells')
myeloid <- df %>% filter(cell_type_main == 'Myeloid cells')

pdf(file = 'data/cell_type_DEG/plots_MBPM_sn_barplots.pdf')
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
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
          plot.margin = margin(0, 8, 0, 0, 'cm'))

ggplot(df, aes(x = organ, fill = cell_type_int)) + 
    geom_bar(position = 'fill', size = 0.25, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('Single nuclei') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
          plot.margin = margin(0, 8, 0, 0, 'cm'))

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
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
          plot.margin = margin(0, 8, 0, 0, 'cm'))

ggplot(immune, aes(x = organ, fill = cell_type_int)) + 
    geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('Single nuclei (only immune cells)') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
          plot.margin = margin(0, 8, 0, 0, 'cm'))

ggplot(immune, aes(x = organ, fill = cell_type_fine)) + 
    geom_bar(position = 'fill', size = 0.25, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('Single nuclei (only immune cells)') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
          plot.margin = margin(0, 5, 0, 0, 'cm'))

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
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.margin = margin(0, 8, 0, 0, 'cm'))

ggplot(tcell, aes(x = organ, fill = cell_type_fine)) + 
    geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('Single nuclei (only T/NK cells)') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
          plot.margin = margin(0, 8, 0, 0, 'cm'))

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
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
          plot.margin = margin(0, 8, 0, 0, 'cm'))

ggplot(myeloid, aes(x = organ, fill = cell_type_fine)) + 
    geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('Single nuclei (only myeloid cells)') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
          plot.margin = margin(0, 8, 0, 0, 'cm'))
dev.off()


#### Bar plots of cell-type frequencies (MBM) ####
seu <- readRDS('data/MBPM/data_MBPM_scn.rds')
seu <- subset(seu, organ == 'Brain')
seu <- subset(seu, cell_type_int %notin% c('Low-quality cells', 'Doublets', 'Contamination',
                                           'Undetermined','Cycling cells'))
seu$immune_status <- ifelse(seu$cell_type_main %in% c('Myeloid cells', 'T/NK cells', 
                                                      'B/Plasma cells'), 
                            'Immune', 'Non-immune')

df <- data.frame(patient = seu$patient, organ = seu$organ, sequencing = seu$sequencing, 
                 cell_type_main = seu$cell_type_main, cell_type_int = seu$cell_type_int, 
                 cell_type_fine = seu$cell_type_fine, immune_status = seu$immune_status)

nonimmune <- df %>% filter(immune_status == 'Non-immune')
immune <- df %>% filter(immune_status == 'Immune')
tcell <- df %>% filter(cell_type_main == 'T/NK cells')
myeloid <- df %>% filter(cell_type_main == 'Myeloid cells')

pdf(file = 'data/cell_type_DEG/plots_MBM_scn_barplots.pdf')
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


#### Bar plots of cell-type frequencies in sc (MBM_sc) ####
seu <- readRDS('data/MBPM/data_MBPM_scn.rds')
seu <- subset(seu, organ == 'Brain')
seu <- subset(seu, cell_type_int %notin% c('Low-quality cells', 'Doublets', 'Contamination', 
                                           'Undetermined','Cycling cells')
              & sequencing == 'Single cells')
seu$immune_status <- ifelse(seu$cell_type_main %in% c('Myeloid cells', 'T/NK cells',
                                                      'B/Plasma cells'), 
                            'Immune', 'Non-immune')

df <- data.frame(patient = seu$patient, organ = seu$organ, sequencing = seu$sequencing, 
                 cell_type_main = seu$cell_type_main, cell_type_int = seu$cell_type_int, 
                 cell_type_fine = seu$cell_type_fine, immune_status = seu$immune_status)

nonimmune <- df %>% filter(immune_status == 'Non-immune')
immune <- df %>% filter(immune_status == 'Immune')
tcell <- df %>% filter(cell_type_main == 'T/NK cells')
myeloid <- df %>% filter(cell_type_main == 'Myeloid cells')

pdf(file = 'data/cell_type_DEG/plots_MBM_sc_barplots.pdf')
## All immune cells 
# Stacked bar plots by patient
ggplot(immune, aes(x = patient, fill = cell_type_main)) + 
    geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('MBM_sc') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(immune, aes(x = patient, fill = cell_type_int)) + 
    geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('MBM_sc') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(immune, aes(x = patient, fill = cell_type_fine)) + 
    geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('MBM_sc') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

## All non-immune cells 
# Stacked bar plots by sequencing
ggplot(nonimmune, aes(x = patient, fill = cell_type_main)) + 
    geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('MBM_sc') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(nonimmune, aes(x = patient, fill = cell_type_int)) + 
    geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('MBM_sc') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(nonimmune, aes(x = patient, fill = cell_type_fine)) + 
    geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('MBM_sc') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

## Only myeloid cells 
# Stacked bar plots by patient
ggplot(myeloid, aes(x = patient, fill = cell_type_int)) + 
    geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('MBM_sc') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(myeloid, aes(x = patient, fill = cell_type_fine)) + 
    geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('MBM_sc') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

## Only T cells 
# Stacked bar plots by patient
ggplot(tcell, aes(x = patient, fill = cell_type_int)) + 
    geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('MBM_sc') + 
    theme(legend.position = 'right', axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(tcell, aes(x = patient, fill = cell_type_fine)) + 
    geom_bar(position = 'fill', size = 0.3, color = 'black') + theme_classic() + 
    ylab('Cell type fraction') + ggtitle('MBM_sc') + 
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
