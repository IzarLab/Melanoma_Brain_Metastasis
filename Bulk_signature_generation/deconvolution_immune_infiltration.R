# Author: Somnath Tagore, Ph.D. Title: Deconvolution and Immune infiltration
# Script Name: deconvolution_immune_infiltration.R 
# Last Updated: 09/24/2021

# Packages required for this analysis
formatR::tidy_app()

#Install packages
install.packages("remotes")
remotes::install_github("icbi-lab/immunedeconv")

library(reshape2)
library(ggplot2)

source('CIBERSORT.R')
library(org.Hs.eg.db)
library(annotate)

library(dplyr)
library(tidyr)
library(immunedeconv)
library(tibble)

#Load TCGA-SKCM data
load('skcm-tissue_rawcounts.rda')
dset.skcm <- rawcounts
saveRDS(dset.skcm,file="dset.skcm.rds")

#Load Davies data
dset.davies <- read.csv("davies.csv")
rownames(dset.davies) <- dset.davies[,1]
dset.davies<-dset.davies[,-1]
saveRDS(dset.davies,file="~/Documents/Ben_Izar_project/Melanoma_Jana/August_13_2021/dset.davies.rawcounts.rds")

tcga.tumor<-rawcounts
tcga.tumor<-dset.davies

#removing #N/A from samples
tcga.tumor <- tcga.tumor[is.finite(rowSums(tcga.tumor)),]
dim(tcga.tumor)

#any(is.na(tcga.normal))
any(is.na(tcga.tumor))

#log2cpm+1 normalization
mixed.tcga.tumor.normal <- tcga.tumor

mixed.tcga.tumor.normal.log2cpm <- apply(mixed.tcga.tumor.normal,2,function(x){
  y <- 1E6*x/sum(x) + 1
  z <- log(y,2)
  return(z)
})
# 
mixed.tcga.tumor.normal.log2cpm.ges <- t(apply(mixed.tcga.tumor.normal.log2cpm,1,function(x){
  y <- (x - mean(x))/sd(x)
  return(y)
}))

dset.log2cpm <- mixed.tcga.tumor.normal.log2cpm.ges

#Convert  Entrez ID to Gene Symbol
rownames(dset.log2cpm)<-getSYMBOL(rownames(dset.log2cpm),data='org.Hs.eg')

#Run de-convolution
res = immunedeconv::deconvolute(dset.log2cpm, "quantiseq", tumor = TRUE)
res_quantiseq = deconvolute(dataset_racle$expr_mat, "quantiseq", tumor = TRUE)
saveRDS(res,file="davies.quantiseq.all.samples.with.tumor.comp.gt.0.7.rds")

#Displaying result
pdf(file="cell.type.bar.plot.tumor.pdf",height=30,width=20)

res %>%
  gather(sample, fraction, -cell_type) %>%
  # plot as stacked bar chart
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
  geom_bar(stat='identity') +
  coord_flip() +
  scale_fill_brewer(palette="Paired") +theme(text = element_text(size=10))+#geom_text(aes(label=round(fraction, digits = 2)))+
  scale_x_discrete(limits = rev(levels(res_quantiseq)))
dev.off()

pdf(file="cell.type.tile.tumor.pdf",height=30,width=20)
ggplot(melt(res), aes(x=cell_type, y=variable, fill=value)) + 
  geom_tile() +
  scale_fill_distiller(palette = "RdBu") +
  geom_text(aes(label=round(value, digits = 2)),size=3) +theme(text = element_text(size = 7))
dev.off()

pdf(file="cell.type.dot.plot.tumor.pdf",height=30,width=20)
res_mcp_counter = deconvolute(dset.log2cpm, "mcp_counter")
res_mcp_counter %>%
  gather(sample, score, -cell_type) %>%
  ggplot(aes(x=sample, y=score, color=cell_type)) +
  geom_point(size=1) +
  facet_wrap(~cell_type, scales="free_x", ncol=3) +
  scale_color_brewer(palette="Paired", guide=FALSE) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text = element_text(size = 2))
dev.off()

