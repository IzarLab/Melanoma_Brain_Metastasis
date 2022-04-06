#!/usr/bin/env Rscript

### title: Counting number of IG combinations in B cells 
### authors: Jana Biermann,PhD; Yiping Wang, PhD
library(Seurat)
library(destiny)
library(dplyr)
library(ggplot2)
library(cowplot)
library(scater)
library(SingleCellExperiment)
library(gplots)
library(viridis)
library(ggvenn)
library(ggalluvial)


seu <- readRDS("data/MBPM/data_MBPM_scn.rds")
seu <- subset(seu, cell_type_int %in%c('Plasma cells','B cells'))
ig <- grep("^IG", rownames(seu@assays$RNA@counts), value = T)

# Create gene expression dataframe with only IG genes, delete empty rows and columns
tab <- as.data.frame(seu@assays$RNA@data[ig, ])
tab <- tab[, colSums(tab) > 0]
tab <- tab[rowSums(tab) > 0, ]

# List heavy and light variable regions
ighv <- grep("^IGHV", rownames(tab), value = T)
iglv <- grep("^(IGLV|IGKV)", rownames(tab), value = T)

# Set up empty dataframe for results
res <- matrix(data = 0, nrow = length(iglv), ncol = length(ighv))
colnames(res) <- ighv
rownames(res) <- iglv
res <- as.data.frame(res)

# List constant regions, manually remove the ones that don't follow the right
# pattern, and set up result dataframe for constant chains
'%notin%' <- Negate('%in%')
ig_con <- grep("^(IGHG|IGHA|IGHM|IGHE)", rownames(tab), value = T)
ig_con <- ig_con[ig_con %notin% c('IGHEP2',"IGHGP",'IGHMBP2')]
res_con <- as.data.frame(matrix(data = NA, nrow = 2000, ncol = 8))
colnames(res_con) <- c("barcode", "orig.ident",'patient','sequencing','organ', "heavy", "light", "constant")

# Loop to identify cells with variable and constant chains expressed
i <- 1
for (c in 1:ncol(tab)) {
  # Remove values from previous loop
  rm("igl")
  rm("igh")
  # Go through gene matrix cell by cell
  cell <- c(tab[, c])
  names(cell) <- rownames(tab)
  
  # Find cell's highest expressed heavy chain
  heav <- cell[ighv]
  igh <- names(sort(heav, decreasing = T)[1])
  
  # Find cell's highest expressed light chain
  lig <- cell[iglv]
  igl <- names(sort(lig, decreasing = T)[1])
  
  # Find cell's highest expressed constant region
  con <- cell[ig_con]
  con <- names(sort(con, decreasing = T)[1])
  
  # Only add cell to variable regions result dataframe if it has both heavy and
  # light variable chains expressed
  # Also add cell to patient dataframe only if groups!="", to avoid doublecounting
  if (cell[igh] > 0 & cell[igl] > 0) {
    res[igl, igh] <- res[igl, igh] + 1
  }
  
  # Only add cell to constant regions result dataframe if it has heavy and light
  # variable chains and a constant chain expressed
  if (cell[igh] > 0 & cell[igl] > 0 & cell[con] > 0) {
    res_con$barcode[i] <- colnames(tab)[c]
    res_con$heavy[i] <- igh
    res_con$light[i] <- igl
    res_con$constant[i] <- con
    res_con$orig.ident[i] <- seu$orig.ident[colnames(seu) == colnames(tab)[c]]
    res_con$patient[i] <- seu$patient[colnames(seu) == colnames(tab)[c]]
    res_con$sequencing[i] <- seu$sequencing[colnames(seu) == colnames(tab)[c]]
    res_con$organ[i] <- seu$organ[colnames(seu) == colnames(tab)[c]]
    i <- i + 1
  }
}

# Remove columns and rows that are empty
res <- res[, colSums(res) > 0]
res <- res[rowSums(res) > 0, ]

# Add IG subtype info
res_con$subtype <- substr(res_con$constant, 1, 4)
res_con <- res_con[!is.na(res_con$light) & !is.na(res_con$heavy), ]
ifelse(!dir.exists(file.path('data/MBPM/IG_analysis')), 
       dir.create(file.path('data/MBPM/IG_analysis'),recursive = T), FALSE)
write.csv(res_con, "data/MBPM/IG_analysis/table_IG_MBPM_bcells_chains.csv", row.names = F, quote = F)

# Save heatmap of combination frequencies, change height dimensions so light labels on y axis don't have
# some dropout if figure too short
pdf(paste0("data/MBPM/IG_analysis/plots_IG_combination_frequency.pdf"), height = 12, width = 12)
out = heatmap.2(as.matrix(res), trace = "none", scale = "none", margins = c(6,6), cexRow = 0.7, 
                cexCol = 0.7, hclustfun = function(d) {hclust(d, method = "average")}, 
                xlab = "Heavy chains", ylab = "Light chains", col = rev(heat.colors(100)), 
                keysize = 1.5, lhei = c(1.7, 5), density.info = "none")
dev.off()

# write out the order of light and heavy chains as plotted in combination frequency heatmap,
# for use in IG_combination_plots.R
lightorder = rev(rownames(res)[out$rowInd])
write.table(lightorder, "data/MBPM/IG_analysis/lightorder.txt", col.names = F, row.names = F, quote = F, sep = "\t")
heavyorder = colnames(res)[out$colInd]
write.table(heavyorder, "data/MBPM/IG_analysis/heavyorder.txt", col.names = F, row.names = F, quote = F, sep = "\t")


# write out overall frequency of light and heavy chains
lightoveralldf = data.frame(lightchain = rownames(res), frequency = rowSums(res))
heavyoveralldf = data.frame(heavychain = colnames(res), frequency = colSums(res))
write.table(lightoveralldf, "data/MBPM/IG_analysis/lightoverallfreq.txt", row.names = F, col.names = T, quote = F, sep = "\t")
write.table(heavyoveralldf, "data/MBPM/IG_analysis/heavyoverallfreq.txt", row.names = F, col.names = T, quote = F, sep = "\t")


# count number of light chain occurrences in each patient
theme_set(theme_bw())
res_patient_tally_light = res_con %>% group_by(light, patient) %>% tally()
res_seq_tally_light = res_con %>% group_by(light, sequencing) %>% tally()
res_organ_tally_light = res_con %>% group_by(light, organ) %>% tally()

pdf("data/MBPM/IG_analysis/plots_usage_light_per_patient.pdf",width = 12)
# print barplot of light chain usage in each patient, sorted by overall light
# chain usage across all patients, using information in lightoveralldf
ggplot(res_patient_tally_light, aes(fill = patient, x = n, y = factor(light, level = lightoveralldf$lightchain[order(lightoveralldf$frequency)]))) + 
  geom_bar(position = "stack", stat = "identity",colour="black",size=0.25) + 
  ylab("Light Chain") + xlab("Usage")+theme_bw()

ggplot(res_seq_tally_light, aes(fill = sequencing, x = n, y = factor(light, level = lightoveralldf$lightchain[order(lightoveralldf$frequency)]))) + 
  geom_bar(position = "stack", stat = "identity",colour="black",size=0.25) + 
  ylab("Light Chain") + xlab("Usage")+theme_bw()

ggplot(res_organ_tally_light, aes(fill = organ, x = n, y = factor(light, level = lightoveralldf$lightchain[order(lightoveralldf$frequency)]))) + 
  geom_bar(position = "stack", stat = "identity",colour="black",size=0.25) + 
  ylab("Light Chain") + xlab("Usage")+theme_bw()
dev.off()


# similar analyses as above, for heavy chains
res_patient_tally_heavy = res_con %>% group_by(heavy, patient) %>% tally()
res_seq_tally_heavy = res_con %>% group_by(heavy, sequencing) %>% tally()
res_organ_tally_heavy = res_con %>% group_by(heavy, organ) %>% tally()

pdf("data/MBPM/IG_analysis/plots_usage_heavy_per_patient.pdf",width = 12)
ggplot(res_patient_tally_heavy, aes(fill = patient, x = n, y = factor(heavy, level = heavyoveralldf$heavychain[order(heavyoveralldf$frequency)]))) + 
  geom_bar(position = "stack", stat = "identity",colour="black",size=0.25) + 
  ylab("Light Chain") + xlab("Usage")+theme_bw()

ggplot(res_seq_tally_heavy, aes(fill = sequencing, x = n, y = factor(heavy, level = heavyoveralldf$heavychain[order(heavyoveralldf$frequency)]))) + 
  geom_bar(position = "stack", stat = "identity",colour="black",size=0.25) + 
  ylab("Light Chain") + xlab("Usage")+theme_bw()

ggplot(res_organ_tally_heavy, aes(fill = organ, x = n, y = factor(heavy, level = heavyoveralldf$heavychain[order(heavyoveralldf$frequency)]))) + 
  geom_bar(position = "stack", stat = "identity",colour="black",size=0.25) + 
  ylab("Light Chain") + xlab("Usage")+theme_bw()
dev.off()


#### Overlap of heavy/light between cell types per patient ####
tmp<-data.frame(cell_type_fine=seu@meta.data[,c("cell_type_fine")],
                barcode=rownames(seu@meta.data))
res_con2<-left_join(res_con,tmp,by='barcode')
res_con2$heavy_light<-paste0(res_con2$heavy,'_',res_con2$light)

for(pat in unique(res_con2$patient)){
  tab <- res_con2 %>% 
    filter(patient == pat) %>% 
    select(cell_type_fine,heavy_light)
  tab$cell_type_fine <- factor(tab$cell_type_fine,
                               levels = c('Naïve B cells','Activated B cells',  'Plasma cells'))
  
  lis_all <- list(
    naive = tab[tab$cell_type_fine=='Naïve B cells',2],
    activated = tab[tab$cell_type_fine=='Activated B cells',2],
    plasma = tab[tab$cell_type_fine=='Plasma cells',2])
  
  ## Plots per patient
  pdf(paste0('data/MBPM/IG_analysis/plots_HL_chains_',pat,'.pdf'))
  # Venn diagram
  print(ggvenn(lis_all,fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"))+
          ggtitle(paste0(pat,': Overlap of B cell heavy and light chains')))
  
  # Alluvial plot
  tab2<- tab %>% 
    select(cell_type_fine,heavy_light) %>%
    group_by(cell_type_fine,heavy_light) %>% 
    tally()
  tab2$cell_type_fine<-factor(tab2$cell_type_fine, 
                              levels = c('Naïve B cells','Activated B cells',  'Plasma cells'))
  tab2$label<-ifelse(tab2$n>= sort(tab2$n,decreasing = T)[5],
                     tab2$heavy_light,NA)
  
  print(ggplot(tab2,aes(y = n, x=cell_type_fine,stratum=heavy_light,
                        alluvium=heavy_light,fill = heavy_light)) +
          geom_alluvium(aes(fill = heavy_light)) +
          guides(fill = FALSE) +
          geom_stratum() +
          geom_text(stat = "stratum", aes(label = label),size=3) +
          ggtitle(paste0(pat,': Overlap of B cell heavy and light chains'))+
          theme_classic()+ylab('Frequency'))
  dev.off()
}
