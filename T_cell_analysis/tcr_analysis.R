#!/usr/bin/env Rscript

### title: Integrating TCR output from Cell Ranger vdj (v6.1.1) and 
### visualization of TCR distribution
### author: Jana Biermann, PhD

print(Sys.time())

library(Seurat)
library(dplyr)
library(ggplot2)
library(gplots)
library(viridis)
library(scales)
library(stringr)
library(reshape2)
library(patchwork)
'%notin%' <- Negate('%in%')

colExp<-c('#C82F71','#7473A6')
colSha<-c('#a67473','#2fc886','#2f71c8')
colMAIT<-c('#009F3F','#E0C126','#F08080','#2fbec8')
colExpPub<- c('#d85a90','#a8285f','#807fae','#504f7c')

#### Add TCR-seq data ####

# Read-in integrated object
seu <- readRDS('data/MBPM/data_MBPM_scn.rds')

# Combine Cell Ranger vdj output
all_tcr <- NULL
full_seu <- seu@meta.data %>%
  select(barcode_all, barcode_backup, orig.ident, patient, sequencing) %>%
  filter(sequencing == 'Single cells')

for (s in c('MBM01_sc', 'MBM02_sc', 'MBM03_sc', 'MBM04_sc', 'MBM05_sc')) {
  filt <- read.csv(paste0('data/', s, '/TCR', s, '/v611/filtered_contig_annotations.csv'))
  clon <- read.csv(paste0('data/', s, '/TCR', s, '/v611/clonotypes.csv'))
  filt$clonotype_id <- filt$raw_clonotype_id
  mer <- left_join(filt, clon, by = 'clonotype_id')
  mer <- mer %>% select(barcode, frequency, proportion, cdr3s_aa, cdr3s_nt, inkt_evidence, mait_evidence)
  mer$BI <- s
  deduped.mer <- unique(mer[, 1:8])
  tmp_seu <- full_seu %>% filter(pat == s)
  tmp_seu <- tmp_seu[grep('_CD45pos', tmp_seu$orig.ident), ]
  tmp_seu <- cbind.data.frame(tmp_seu, 
                              colsplit(string = tmp_seu$barcode_backup, pattern = '_', names = c('barcode', 'rest')))
  final <- tmp_seu %>% select(barcode_all, barcode)
  final <- left_join(final, deduped.mer, by = 'barcode')
  all_tcr <- rbind.data.frame(all_tcr, final)
}
all_tcr$barcode <- NULL
write.csv(all_tcr, 'data/TCR/data_all_TCR_v611.csv', row.names = F)

# Create data frame with unique TCRs
uniq_tcr <- data.frame(cdr3s_aa = all_tcr[!duplicated(all_tcr$cdr3s_aa), ]$cdr3s_aa)
uniq_tcr <- na.omit(uniq_tcr)

# Check how many TR combined
which(str_count(uniq_tcr$cdr3s_aa, pattern = 'TRA:') == 3)
which(str_count(uniq_tcr$cdr3s_aa, pattern = 'TRB:') == 3)

# Generate columns in uniq_tcr to prepare for splitting
uniq_tcr$TRB1 <- NA
uniq_tcr$TRB2 <- NA
uniq_tcr$TRB3 <- NA
uniq_tcr$TRA1 <- NA
uniq_tcr$TRA2 <- NA

# Fill columns
for (r in 1:nrow(uniq_tcr)) {
  tmp <- uniq_tcr[r, 'cdr3s_aa']
  if (str_count(tmp, 'TRB:') == 1 & str_count(tmp, 'TRA:') == 1) {
    tmp2 <- str_split(tmp, ';TRA:')
    tra <- tmp2[[1]][2]
    trb <- tmp2[[1]][1]
    trb <- str_remove(trb, 'TRB:')
    uniq_tcr[r, 'TRB1'] <- trb
    uniq_tcr[r, 'TRA1'] <- tra
  }
  
  if (str_count(tmp, 'TRB:') == 2 & str_count(tmp, 'TRA:') == 1) {
    tmp2 <- str_split(tmp, ';TRA:')
    tra <- tmp2[[1]][2]
    trb <- str_split(tmp2[[1]][1], ';TRB:')
    trb2 <- trb[[1]][2]
    trb1 <- str_remove(trb[[1]][1], 'TRB:')
    uniq_tcr[r, 'TRB1'] <- trb1
    uniq_tcr[r, 'TRB2'] <- trb2
    uniq_tcr[r, 'TRA1'] <- tra
  }
  
  if (str_count(tmp, 'TRB:') == 3 & str_count(tmp, 'TRA:') == 1) {
    tmp2 <- str_split(tmp, ';TRA:')
    tra <- tmp2[[1]][2]
    trb <- str_split(tmp2[[1]][1], ';TRB:')
    trb2 <- trb[[1]][2]
    trb3 <- trb[[1]][3]
    trb1 <- str_remove(trb[[1]][1], 'TRB:')
    uniq_tcr[r, 'TRB1'] <- trb1
    uniq_tcr[r, 'TRB2'] <- trb2
    uniq_tcr[r, 'TRB3'] <- trb3
    uniq_tcr[r, 'TRA1'] <- tra
  }
  
  if (str_count(tmp, 'TRB:') == 1 & str_count(tmp, 'TRA:') == 2) {
    tmp2 <- str_split(tmp, ';TRA:')
    trb <- str_split(tmp2[[1]][1], 'TRB:')[[1]][2]
    tra1 <- tmp2[[1]][2]
    tra2 <- tmp2[[1]][3]
    uniq_tcr[r, 'TRB1'] <- trb1
    uniq_tcr[r, 'TRA1'] <- tra1
    uniq_tcr[r, 'TRA2'] <- tra2
  }
  
  if (str_count(tmp, 'TRB:') == 2 & str_count(tmp, 'TRA:') == 2) {
    tmp2 <- str_split(tmp, ';TRA:')
    trb2 <- str_split(tmp2[[1]][1], ';TRB:')[[1]][2]
    trb1 <- str_split(str_split(tmp2[[1]][1], ';TRB:')[[1]][1], 'TRB:')[[1]][2]
    tra1 <- tmp2[[1]][2]
    tra2 <- tmp2[[1]][3]
    uniq_tcr[r, 'TRB1'] <- trb1
    uniq_tcr[r, 'TRB2'] <- trb2
    uniq_tcr[r, 'TRA1'] <- tra1
    uniq_tcr[r, 'TRA2'] <- tra2
  }
  
  if (str_count(tmp, 'TRB:') == 1 & str_count(tmp, 'TRA:') == 0) {
    trb <- str_split(tmp, 'TRB:')[[1]][2]
    uniq_tcr[r, 'TRB1'] <- trb
  }
  if (str_count(tmp, 'TRB:') == 0 & str_count(tmp, 'TRA:') == 1) {
    tra <- str_split(tmp, 'TRA:')[[1]][2]
    uniq_tcr[r, 'TRA1'] <- tra
  }
}


#### Add public TCRs ####

# Read-in public TCR cohorts and combine
val <- read.csv('data/TCR/TCR-cdr3/01_Valpione.csv')
yuk <- read.csv('data/TCR/TCR-cdr3/02_Yusko.csv')
pru <- read.csv('data/TCR/TCR-cdr3/03_Pruessmann.csv')
vir <- read.csv('data/TCR/TCR-cdr3/viral_davis.csv')
public <- rbind.data.frame(val, yuk, pru, vir)
public <- public[!duplicated(public$amino_acid), ]

# Add info for TRB1
uniq_tcr$amino_acid <- NULL
uniq_tcr$amino_acid <- uniq_tcr$TRB1
uniq_tcr <- left_join(uniq_tcr, public, by = 'amino_acid')
colnames(uniq_tcr)[8:11] <- c('sample_TRB1', 'nt_TRB1', 'cohort_TRB1', 'disease_TRB1')

# Add info for TRB2
uniq_tcr$amino_acid <- uniq_tcr$TRB2
uniq_tcr <- left_join(uniq_tcr, public, by = 'amino_acid')
colnames(uniq_tcr)[12:15] <- c('sample_TRB2', 'nt_TRB2', 'cohort_TRB2', 'disease_TRB2')

# Add info for TRB3
uniq_tcr$amino_acid <- uniq_tcr$TRB3
uniq_tcr <- left_join(uniq_tcr, public, by = 'amino_acid')
colnames(uniq_tcr)[16:19] <- c('sample_TRB3', 'nt_TRB3', 'cohort_TRB3', 'disease_TRB3')
uniq_tcr$amino_acid <- NULL

# Merge unique TCRs back into all_tcr
all_tcr <- left_join(all_tcr, uniq_tcr, by = 'cdr3s_aa')
write.csv(all_tcr, 'data/TCR/data_all_TCR_v611.csv', row.names = F)

# Add to Seurat and subset
seu@meta.data <- left_join(seu@meta.data, all_tcr, by = 'barcode_all')
rownames(seu@meta.data) <- seu@meta.data$barcode_all
seu$tcr <- ifelse(is.na(seu$cdr3s_aa) == T, 'no_TCR', 'TCR')

seu <- subset(seu, cell_type_fine %in% c('CD8+ T cells TOX+', 'CD8+ T cells TCF7+', 'CD4+ T cells', 
                                         'Tfh-like cells', 'Tregs', 'Cycling T/NK cells') & 
                sequencing == 'Single cells' & tcr == 'TCR')
seu$mait <- ifelse(seu$mait_evidence %in% c('TRA:gene', 'TRA:gene;TRB:gene', 'TRA:gene+junction', 
                                            'TRA:gene+junction;TRB:gene', 'TRB:gene'), 'MAIT', 'no_MAIT')
seu$inkt <- ifelse(seu$inkt_evidence %in% c('TRB:gene'), 'iNKT', 'no_iNKT')
seu$public_info <- paste0(seu$disease_TRB1, ';', seu$mait, ';', seu$inkt)

seu$both_chains <- ifelse(grepl('TRB:', seu$cdr3s_aa) == T & grepl('TRA:', seu$cdr3s_aa) == T, 'both', 'not_both')
seu$clone_size <- ifelse(seu$frequency > 1 & seu$both_chains == 'both', 'Expanded', 
                         ifelse(seu$frequency == 1 & seu$both_chains == 'both', 'Non-expanded', NA))
seu$mait_inkt <- paste0(seu$mait, '_', seu$inkt)


#### Plot TCR data ####

### Count TCR overlap in CD8+
seu <- readRDS('data/MBPM/data_MBPM_scn.rds')
seu$public_info <- paste0(seu$disease_TRB1, ';', seu$mait, ';', seu$inkt)
seu$public<-ifelse(grepl('melanoma',seu$disease_TRB1) |
                     grepl('melanoma',seu$disease_TRB2),'public',
                   ifelse(seu$tcr=='TCR','private',NA))
seu$exp_pub<-paste0(seu$clone_size,'_',seu$public)

seu <- subset(seu, sequencing == 'Single cells' & 
                cell_type_main=='T/NK cells' &
                cell_type_int %notin% c('Low-quality cells','Doublets','Contamination',
                                        'Undetermined','NK cells'))

seu_both <- subset(seu, both_chains == 'both' & cell_type_int == 'CD8+ T cells')

res <- data.frame(patient = unique(seu_both$patient))
for (pat in unique(seu_both$patient)) {
  res$shared_TCR[res$patient == pat] <- 
    length(intersect(subset(seu_both, patient == pat & 
                              cell_type_fine == 'CD8+ T cells TCF7+')$cdr3s_aa, 
                     subset(seu_both, patient == pat & 
                              cell_type_fine == 'CD8+ T cells TOX+')$cdr3s_aa))
  res$unique_TOX[res$patient == pat] <- 
    length(unique(subset(seu_both, patient == pat & 
                           cell_type_fine == 'CD8+ T cells TOX+')$cdr3s_aa)) - res$shared_TCR[res$patient == pat]
  res$unique_TCF7[res$patient == pat] <- 
    length(unique(subset(seu_both, patient == pat & 
                           cell_type_fine == 'CD8+ T cells TCF7+')$cdr3s_aa)) - res$shared_TCR[res$patient == pat]
}

# Plots
df <- seu@meta.data[, c('barcode_all', 'patient', 'frequency', 'clone_size', 'cdr3s_aa', 'cell_type_fine', 
                        'mait', 'mait_evidence', 'inkt', 'inkt_evidence', 'both_chains', 'disease_TRB1', 
                        'disease_TRB2', 'exp_pub', 'mait_inkt')]
df_both <- df %>% filter(both_chains == 'both')

pdf(file = 'data/MBPM/tcells/plots_MBM_sc_tcell_TCR.pdf')
# Stacked bar plot both_chains
ggplot(df, aes(x = patient, fill = both_chains)) + 
  geom_bar(colour = 'black', position = 'fill', size = 0.25) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  ggtitle('TCR chains recovered across patients') + 
  xlab('Patient') + ylab('Fraction')

# Stacked bar plot mait_inkt
ggplot(df_both, aes(x = patient, fill = mait_inkt)) + 
  geom_bar(colour = 'black', position = 'fill', size = 0.25) + 
  ggtitle('MAIT + iNKT in sc T cells') + 
  theme_classic() + 
  xlab('') + 
  ylab('Fraction of T cells') + 
  scale_fill_manual(values = colMAIT)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0))

# Stacked bar plot exp_pub
df_exp_pub <- df_both %>% group_by(patient,exp_pub) %>% summarise(n = n())
ggplot(df_exp_pub, aes(x = patient,y=n, fill = exp_pub)) + 
  geom_bar(stat="identity",position=position_dodge(), colour = 'black',size = 0.25) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  xlab('') + ylab('# T cells') + scale_fill_manual(values = colExpPub)

# Stacked bar plot exp_pub CD8
df_exp_pub_cd8 <- df_both %>% filter(grepl('CD8',cell_type_fine)) %>% group_by(patient,exp_pub) %>% summarise(n = n())
ggplot(df_exp_pub_cd8, aes(x = patient,y=n, fill = exp_pub)) + 
  geom_bar(stat="identity",position=position_dodge(), colour = 'black',size = 0.25) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  xlab('') + ylab('# CD8+ T cells') + scale_fill_manual(values = colExpPub)

# Bar plot clones
df_unique <- df_both %>% group_by(patient) %>% summarise(nClones = n_distinct(cdr3s_aa))
p1<-ggplot(df_unique, aes(x = patient, y = nClones)) + 
  geom_bar(stat = 'identity') + 
  theme_classic() + 
  theme(axis.text.x = element_blank(),axis.title.x = element_blank()) + 
  ylab('# Clones')

# Bar plot cells
df_cells <- df_both %>% group_by(patient) %>% summarise(nCells = n_distinct(barcode_all))
p2<- ggplot(df_cells, aes(x = patient, y = nCells)) + 
  geom_bar(stat = 'identity') + theme_classic() + 
  theme(axis.text.x = element_blank(),axis.title.x = element_blank()) + 
  ylab('# Cells')

# Stacked bar plot clone_size
p3<-ggplot(df_both, aes(x = patient, fill = clone_size)) + 
  geom_bar(colour = 'black', position = 'fill', size = 0.25) + 
  theme_classic() + scale_fill_manual(values = colExp)+
  theme(axis.title.x = element_blank()) + 
  ylab('Fraction of T cells')

p1/p2/p3

# Circle stacked bar plot stratified clonotype + pat
df_ct_strat <- df_both %>%
  group_by(patient, clone_size, cell_type_fine) %>%
  summarise(n = n()) %>%
  mutate(freq = n/sum(n))

ggplot(df_ct_strat, aes(x = '', y = freq, group = cell_type_fine, 
                        fill = cell_type_fine, color = cell_type_fine)) + 
  geom_bar(width = 1, stat = 'identity', color = 'white', size = 0.1) + 
  coord_polar('y', start = 0) + 
  theme_void() + 
  facet_grid(clone_size ~ patient)

# Circle stacked bar plot stratified cell type + pat
df_clonesize_strat <- df_both %>%
  group_by(patient, cell_type_fine, clone_size) %>%
  summarise(n = n()) %>%
  mutate(freq = n/sum(n))

ggplot(df_clonesize_strat, aes(x = '', y = freq, group = clone_size, fill = clone_size, color = clone_size)) + 
  geom_bar(width = 1, stat = 'identity', color = 'white', size = 0.1) + 
  coord_polar('y', start = 0) + 
  theme_void() + scale_fill_manual(values = colExp)+
  facet_grid(cell_type_fine ~ patient)

# Stacked bar plot TCR overlap
df_overlap <- melt(res) %>% group_by(patient, variable)
ggplot(df_overlap, aes(x = patient, y = value, fill = variable)) + 
  geom_bar(colour = 'black', position = 'fill', stat = 'identity') + theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ggtitle('Distribution of TCRs between CD8+ T cells') + 
  xlab('Patient') + 
  ylab('Fraction of CD8+ T cells')+
  scale_fill_manual(values = colSha)
dev.off()

