
# Author: Somnath Tagore, Ph.D. Title: Master Regulator Analysis of Melanoma data 
# Script Name: protein_activity_melanoma.R 
# Last Updated: 09/24/2021

# Packages required for this analysis
formatR::tidy_app()
library(cluster)
library(ggplot2)
library(viper)
library(annotate)
library(dplyr)
library(Seurat)
library(pheatmap)
library(RColorBrewer)
library(fpc)
library(ggrepel)
library(pheatmap)
library(org.Hs.eg.db)
library(ReactomePA)

save_pheatmap_pdf <- function(x, filename, width=10, height=10) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

quantile_breaks <- function(xs, n = 10){
  breaks <- quantile(xs, probs = seq(from = 0, to = 1, length.out = n))
  unique(breaks)
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

ClusterColors <- function(k, offset = 0) {
  hues <- seq(15, 375, length = k + 1) + offset
  return(hcl(h = hues, l = 65, c = 100)[1:k])
}

#Read tumor data for each sample
#Here, each sample correspond to categoris MBM/MPM
#An example snippet is given below. Please change the sample names ('sample_id') as required.

sample_id <- 'MPM05_sn.rds'
data_tumor <- readRDS(file=sample_id)

#Prepare data for network reconstruction using ARACNe
#Normalization using log2(TPM+1)
data_tumor <- t(t(as.matrix(data_tumor)) / (colSums(as.matrix(data_tumor)) / 1e6))
data_tumor <- log2(data_tumor + 1)
saveRDS(data_tumor,file="data_tumor_tpm_count.rds")

#Run ARACNe (https://github.com/califano-lab/ARACNe-AP)
#Read the network files as prepared by ARACNe
#Read all network files corresponding to each sample or each sub-cluster (Sample-wise)
sample_id_regulon_1 <- readRDS(file="sample_id_pruned_1.rds")
sample_id_regulon_2 <- readRDS(file="sample_id_pruned_2.rds") 

#Prepare data for protein activity
#Internal signature generation
data_tumor <- t(t(as.matrix(data_tumor)) / (colSums(as.matrix(data_tumor)) / 1e6))
data_tumor <- log2(data_tumor + 1)
mean_data<-apply(data_tumor,1,mean)
sd_data<-apply(data_tumor,1,sd)
signature_int<-(data_tumor - mean_data)/sd_data
signature_int[1:50,1:4]
dim(signature_int)
signature_int <- na.omit(signature_int)
saveRDS(signature_int,file="data_tumor_all_internal_reference_signature.rds")

#Protein activity prediction using VIPER
data_tumor_all_internal_reference_vpmat <- viper(eset = signature_int, regulon = as.list(sample_id_regulon_1,sample_id_regulon_2), method = "none")
saveRDS(data_tumor_all_internal_reference_vpmat,file="data_tumor_all_internal_reference_vpmat.rds")

#Clustering based on protein activity
# pamk clustering based on vipersimilarity
input<- data_tumor_all_internal_reference_vpmat
viperSimilarity_dist <- as.dist(viperSimilarity(input))
dim(viperSimilarity_dist)
viperSimilarity_dist[1:5]
saveRDS(viperSimilarity_dist,file="viperSimilarity_dist.rds")

silhouette_score_mean<- c()
silhouette_score_total<-list()
for (k in k_min:k_max) {
  print(k)
  clustering <- pam(viperSimilarity_dist, k)
  sil <- silhouette(clustering$cluster, viperSimilarity_dist)
  silhouette_score_mean <- c(silhouette_score_mean, mean(sil[,3]))
  silhouette_score_total[[k-1]]<-sil
}
silhouette_score_summary <- data.frame('k' = k_min:k_max, 'Silhouette_Score' = silhouette_score_mean)
saveRDS(silhouette_score_summary,file="silhouette_score_summary.rds")

#Plot silhouette score summary
ggplot(silhouette_score_summary, aes(x=k, y=Silhouette_Score)) +
  geom_point() + geom_line()+ theme(axis.text = element_text(size=15), text = element_text(size=20)) +
  labs(y="Silhouette Score")
clustering_opt<-silhouette_score_total[[which.max(silhouette_score_summary[,2])]]
rownames(clustering_opt)<-colnames(input)
p1<-fviz_silhouette(clustering_opt)
print(p1)
saveRDS(clustering_opt,file="clustering_opt.rds")
saveRDS(silhouette_score_total,file="silhouette_score_total.rds")

#Generate PCA and UMAP visualization for VIPER
input<-data_tumor_all_internal_reference_vpmat
silhouette_score_summary<-silhouette_score_summary
clustering_opt<-clustering_opt

dim(silhouette_score_summary)
silhouette_score_summary
vpmat<-input[,match(colnames(input),rownames(clustering_opt))]
dim(vpmat)

vpmat <- t(vpmat)

vpmat <- as.data.frame(vpmat)

#Extracting cluster-wise cells
vpmat$Cluster <- clustering_opt[,colnames(clustering_opt)=='cluster']
dim(vpmat)

#Assigning unique colors to clusters
clust.colors <- ClusterColors(length(unique(vpmat$Cluster)))
names(clust.colors) <- unique(vpmat$Cluster)

#PCA
pca.obj <- prcomp(vpmat)

#Visualize PCA
dtp <- data.frame("Cluster" = vpmat$Cluster, pca.obj$x[,1:2]) # the first two componets are selected (NB: you can also select 3 for 3D plottings or 3+)
ggplot(data = dtp) +
  geom_point(aes(x = PC1, y = PC2), color = vpmat$Cluster) +
  theme_bw()
ggsave(file="PCA_protein_activity_tumor_updated.pdf")

#UMAP
tumor.umap = umap(vpmat[,1:622])

#Creating a dataframe for visualization
df <- data.frame(x = tumor.umap$layout[,1],
                 y = tumor.umap$layout[,2],
                 Cluster = vpmat$Cluster)

#Visualize UMAP
ggplot(data = df) +
  geom_point(aes(x = x, y = y), color = vpmat$Cluster) +
  theme_bw()
ggsave(file="UMAP_protein_activity_tumor_updated.pdf")

#Differential Protein activity across VIPER-inferred protein activity clusters
#Here, protein_activity_cluster_1 = cells assigned to cluster 1 based on vpmat$Cluster 
#and protein_activity_cluster_2 = cells assigned to cluster 2 based on vpmat$Cluster

int.ref.pamk <- rowTtest(protein_activity_cluster_1, protein_activity_cluster_2)
write.csv(int.ref.pamk,file="int.ref.pamk.csv")

#Display top/bottom ('top.bottom') proteins based on 'int.ref.pamk'
vpmat.sorted <- vpmat[top.bottom,]

#Assign cluster ids as per vpmat$Cluster
Cluster <- vpmat$Cluster

#Visualize Master Regulator Heatmap
mat.breaks <- quantile_breaks(vpmat.sorted,n = 100)
test1 <-  pheatmap(vpmat.sorted, main = "Viper Clustering: Master Regulators: Differentially Active Proteins", 
                   fontsize = 16, annotation_col = Cluster,
                   cluster_cols = FALSE, show_colnames = FALSE, cluster_rows = FALSE, 
                   show_rownames = TRUE, fontsize_row = 10, fontsize_col = 12, breaks = mat.breaks, 
                   color = colorRampPalette(rev(brewer.pal(n = 7,  name = "RdBu")))(length(mat.breaks)))
save_pheatmap_pdf(test1, "VIPER_MRs_internal_reference.pdf")

