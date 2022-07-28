# Author: Somnath Tagore, Ph.D. Title: Running KINOMO for performing Non-negative Matrix Factorization using gene expression data 
# Script Name: kinomo_run.R
# Last Updated: 01/24/2022

#Instructions
#The KINOMO repsository can be accessed via https://github.com/IzarLab/KINOMO.git
#Convert the gene expression data (raw counts) to Seurat object and perform necessary normalization, scaling etc.
#Run the following R script

#!/usr/bin/env Rscript

## NMF

print(Sys.time())

library(dplyr)
library(Seurat)
library(ggplot2)
library(gplots)
library(NMF)
library(registry)
library(rngtools)
library(purrr)
library(cowplot)
library(stringr)
library(pkgmaker)
library(cluster)

#source("/KINOMO/KINOMO.R")

pat <- 'Filename'

integrated<-readRDS(paste0(Filename,".rds"))

#Seurat Object creation and necessary pre-filteration, normalization. Please customize if necessary.
DefaultAssay(object = integrated) <- "RNA"

seu.1<-integrated

mito.genes <- grep(pattern = "^MT-", x = rownames(seu.1), value = TRUE)
rbl.genes <- grep(pattern = "^RB-", x = rownames(seu.1), value = TRUE)
rsl.genes <- grep(pattern = "^RS-", x = rownames(seu.1), value = TRUE)
rpl.genes <- grep(pattern = "^RPL-", x = rownames(seu.1), value = TRUE)
rbl.genes <- grep(pattern = "^RBL-", x = rownames(seu.1), value = TRUE)
rps.genes <- grep(pattern = "^RPS-", x = rownames(seu.1), value = TRUE)
rbs.genes <- grep(pattern = "^RBS-", x = rownames(seu.1), value = TRUE)
rbl1.genes <- grep(pattern = "^RB", x = rownames(seu.1), value = TRUE)
rsl1.genes <- grep(pattern = "^RS", x = rownames(seu.1), value = TRUE)
rpl1.genes <- grep(pattern = "^RPL", x = rownames(seu.1), value = TRUE)
rbl1.genes <- grep(pattern = "^RBL", x = rownames(seu.1), value = TRUE)
rps1.genes <- grep(pattern = "^RPS", x = rownames(seu.1), value = TRUE)
rbs1.genes <- grep(pattern = "^RBS", x = rownames(seu.1), value = TRUE)

counts <- GetAssayData(seu.1, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% c(mito.genes,rbl.genes,rsl.genes,rpl.genes,rbl.genes,rps.genes,rbs.genes,rbl1.genes,rsl1.genes,rpl1.genes,rbl1.genes,rps1.genes,rbs1.genes))),]
seu <- subset(seu.1, features = rownames(counts))

seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu, features = rownames(seu))
seu <- RunPCA(seu, features = VariableFeatures(object = seu))
seu <- FindNeighbors(seu, dims = 1:15)
seu <- FindClusters(seu,resolution=0.4)
seu <- RunUMAP(seu, dims = 1:20)

mat<-as.matrix(seu@assays$RNA@data)
mat<-mat[rowSums(mat)>0,]

#Estimating the factorization rank
#A critical parameter in KINOMO-NMF is the factorization rank r. It defines the number of metagenes used to approximate the target matrix. Given a 
#NMF method and 
#the target matrix, a common way of deciding on r is to try different values, compute some quality measure of the results, and choose the best value 
#according to this quality criteria.
#Several approaches have then been proposed to choose the optimal value of r. For example, (Brunet2004) proposed to take the first value of r 
#for which the cophenetic coefficient #starts decreasing, (Hutchins2008) suggested to choose the first value where the RSS curve presents an 
#inflection point, and (Frigyesi2008) considered the smallest value at which #the decrease in the RSS is lower than the decrease of the RSS obtained 
#from random data.

pdf(file = paste0(pat,"_KINOMO_nmf.pdf"),height = 20,width=20)
#nmf.estimate.rank <- KINOMO(mat, 2:10, nrun=10)
nmf.estimate.rank <- nmf(mat, 2:10, nrun=10)
plot(nmf.estimate.rank)
dev.off()

#test
seu1<-seu
seu2<-seu
seu3<-seu
seu4<-seu
seu5<-seu

#Perform factorization using estimated rank
myrank<-4
#nmf.estimate.rank <- KINOMO(mat, rank=myrank)
nmf.estimate.rank <- nmf(mat, rank=myrank)
saveRDS(nmf.estimate.rank,paste0(pat,"/",pat,"_KINOMO_nmf_rank_",myrank,".rds"))
res<-nmf.estimate.rank

#top30 meta-genes
tab<-matrix(NA,30,myrank)
colnames(tab)<-paste0(pat,'_factor',seq(1:myrank))
for(c in 1:myrank){
  genes<-order(as.data.frame(res@fit@W)[,c],decreasing = T) %>% head(n=30)
  tab[,c]<-rownames(as.data.frame(res@fit@W))[genes]
  colnames(tab)[c]<-paste0(pat,'_factor',c)
}
write.csv(tab,paste0(pat,"/",pat, "_KINOMO_nmf_rank_",myrank,"_top30_W.csv"))
for(c in 1:ncol(tab)){
  seu1<-AddModuleScore(object = seu1,features = list(tab[,c]),name = colnames(tab)[c],assay = 'RNA',search=T)
}
genelist <- c(tab[,1],tab[,2],tab[,3])
genelist <- genelist[!duplicated(genelist)]


#top50 meta-genes
tab1<-matrix(NA,50,myrank)
colnames(tab1)<-paste0(pat,'_factor',seq(1:myrank))
for(c in 1:myrank){
  genes<-order(as.data.frame(res@fit@W)[,c],decreasing = T) %>% head(n=50)
  tab1[,c]<-rownames(as.data.frame(res@fit@W))[genes]
  colnames(tab1)[c]<-paste0(pat,'_factor',c)
}
write.csv(tab1,paste0(pat,"/",pat, "_KINOMO_nmf_rank_",myrank,"_top50_W.csv"))
for(c in 1:ncol(tab1)){
  seu2<-AddModuleScore(object = seu2,features = list(tab1[,c]),name = colnames(tab1)[c],assay = 'RNA',search=T)
}
genelist1 <- c(tab1[,1],tab1[,2],tab1[,3])
genelist1 <- genelist1[!duplicated(genelist1)]


#top100 meta-genes
tab2<-matrix(NA,100,myrank)
colnames(tab2)<-paste0(pat,'_factor',seq(1:myrank))
for(c in 1:myrank){
  genes<-order(as.data.frame(res@fit@W)[,c],decreasing = T) %>% head(n=100)
  tab2[,c]<-rownames(as.data.frame(res@fit@W))[genes]
  colnames(tab2)[c]<-paste0(pat,'_factor',c)
}
write.csv(tab2,paste0(pat,"/",pat, "_KINOMO_nmf_rank_",myrank,"_top100_W.csv"))
for(c in 1:ncol(tab2)){
  seu3<-AddModuleScore(object = seu2,features = list(tab2[,c]),name = colnames(tab2)[c],assay = 'RNA',search=T)
}
genelist2 <- c(tab2[,1],tab2[,2],tab2[,3])
genelist2 <- genelist2[!duplicated(genelist2)]


#top200 meta-genes
tab3<-matrix(NA,200,myrank)
colnames(tab3)<-paste0(pat,'_factor',seq(1:myrank))
for(c in 1:myrank){
  genes<-order(as.data.frame(res@fit@W)[,c],decreasing = T) %>% head(n=200)
  tab3[,c]<-rownames(as.data.frame(res@fit@W))[genes]
  colnames(tab3)[c]<-paste0(pat,'_factor',c)
}
write.csv(tab3,paste0(pat,"/",pat, "_KINOMO_nmf_rank_",myrank,"_top200_W.csv"))
for(c in 1:ncol(tab3)){
  seu4<-AddModuleScore(object = seu3,features = list(tab3[,c]),name = colnames(tab3)[c],assay = 'RNA',search=T)
}
genelist3 <- c(tab3[,1],tab3[,2],tab3[,3])
genelist3 <- genelist3[!duplicated(genelist3)]


#top300 meta-genes
tab4<-matrix(NA,300,myrank)
colnames(tab4)<-paste0(pat,'_factor',seq(1:myrank))
for(c in 1:myrank){
  genes<-order(as.data.frame(res@fit@W)[,c],decreasing = T) %>% head(n=300)
  tab4[,c]<-rownames(as.data.frame(res@fit@W))[genes]
  colnames(tab4)[c]<-paste0(pat,'_factor',c)
}
write.csv(tab4,paste0(pat,"/",pat, "_KINOMO_nmf_rank_",myrank,"_top300_W.csv"))
for(c in 1:ncol(tab4)){
  seu5<-AddModuleScore(object = seu5,features = list(tab4[,c]),name = colnames(tab4)[c],assay = 'RNA',search=T)
}
genelist4 <- c(tab4[,1],tab4[,2],tab4[,3])
genelist4 <- genelist4[!duplicated(genelist4)]

#colnames(seu@meta.data)


#Plot meta-gene signatures
pdf(file = paste0(pat,"/",pat, "_KINOMO_nmf_rank_",myrank,".pdf"),height = 20,width=20)

print(FeaturePlot(seu1, features = c(paste0(pat,'_factor',seq(1,myrank),'1')),min.cutoff = "q10", max.cutoff = "q99"))
print(FeaturePlot(seu2, features = c(paste0(pat,'_factor',seq(1,myrank),'1')),min.cutoff = "q10", max.cutoff = "q99"))
print(FeaturePlot(seu3, features = c(paste0(pat,'_factor',seq(1,myrank),'1')),min.cutoff = "q10", max.cutoff = "q99"))
print(FeaturePlot(seu4, features = c(paste0(pat,'_factor',seq(1,myrank),'1')),min.cutoff = "q10", max.cutoff = "q99"))
print(FeaturePlot(seu5, features = c(paste0(pat,'_factor',seq(1,myrank),'1')),min.cutoff = "q10", max.cutoff = "q99"))

dev.off()
write.csv(res@fit@W,paste0(pat,"/",pat, "_KINOMO_nmf_rank_",myrank,"_W.csv"))
write.csv(res@fit@H,paste0(pat,"/",pat, "_KINOMO_nmf_rank_",myrank,"_H.csv"))

#Plot meta-genes heatmaps
pdf(file=paste0(pat,"/",pat, "_KINOMO_nmf_rank_",myrank,"_basismap_30.pdf"),width = 10, height = 10)
basismap(res[rownames(res)%in%genelist,])
dev.off()

pdf(file=paste0(pat,"/",pat, "_KINOMO_nmf_rank_",myrank,"_basismap_50.pdf"),width = 20, height = 20)
#basismap(res[rownames(res)%in%c("FAM87B","LINC01128", "LINC00115","FAM41C","AL645608.6"),],color = "-RdYlBu")
basismap(res[rownames(res)%in%genelist1,])
dev.off()

pdf(file=paste0(pat,"/",pat, "_KINOMO_nmf_rank_",myrank,"_basismap_100.pdf"),width = 30, height = 30)
basismap(res[rownames(res)%in%genelist2,])
dev.off()

pdf(file=paste0(pat,"/",pat, "_KINOMO_nmf_rank_",myrank,"_basismap_200.pdf"),width = 40, height = 40)
basismap(res[rownames(res)%in%genelist3,])
dev.off()

pdf(file=paste0(pat,"/",pat, "_KINOMO_nmf_rank_",myrank,"_basismap_300.pdf"),width = 50, height = 50)
basismap(res[rownames(res)%in%genelist4,])
dev.off()

#After the script successfully runs, the following files would be created in the current directory:
# 1. Rank plot
# 2. Meta-gene signature UMAPS
# 3. Meta-gene heatmaps
# 4. Associated .csv files
# 5. Associated .rds files
