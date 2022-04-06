library(Seurat)
#library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SingleR)
library(SingleCellExperiment)
library(scater)
library(pheatmap)
library(rlist)
library(stringr)
library(grid)
library(ggrastr)
library(metR)
library(akima)

homedir = "/data"

#read in table of signature genes, manually assigned by grouping spatial DEG
sigs_to_test_orig = read.table("Moran_i_genes_grouped_jca_v2.csv",header=T,sep=",",quote=NULL)
sigs_to_test_orig = list.append(sigs_to_test_orig,c("VIM"))
names(sigs_to_test_orig)[length(sigs_to_test_orig)] = "VIM"
sigs_to_test_orig = list.append(sigs_to_test_orig,c("TIMP1"))
names(sigs_to_test_orig)[length(sigs_to_test_orig)] = "TIMP1"
patslist = list(c("MBM05_rep1_slide"),c("MBM06_slide"),c("MBM07_slide"),c("MBM08_slide"),c("MBM11_rep1_slide"),c("MBM18_slide"),c("MBM13_slide"),c("MPM08_pre_slide"),c("MPM10_slide"),c("MPM06_slide"),c("MBM05_rep2_slide"),c("MBM11_rep2_slide"),c("puck5"),c("puck6final"),c("puck7_20_feature_threshold"),c("puck8_20_feature_threshold"))

#read in gene signatures to be used only for MBM13 scoring
mbm13_only = read.table("MBPM_signatures_spatial_scoring.csv",header=T,sep=",",quote=NULL)
#nonsense = nonsense+1
mbm13_only = mbm13_only[,c("AXL_sig","MITF_sig","GOBP_OXIDATIVE_PHOSPHORYLATION","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")]

for (pats in patslist)
{
  #determine which signatures to plot, depending on whether MBM13 sample is being processed
  if (pats[1]=="MBM13_slide")
  {
    sigs_to_test = sigs_to_test_orig
    for (asig in names(mbm13_only)) {
      sigs_to_test = list.append(sigs_to_test, mbm13_only[[asig]])
      names(sigs_to_test)[length(sigs_to_test)] = asig
    }
  }
  else
  {
    sigs_to_test = sigs_to_test_orig
  }
  pdf(paste0("/data/mbpm_puckdata_add_rctd_and_sigs/Moran_i_genes_grouped_jca_",pats[1],".pdf"))
  for (patidx in 1:length(pats)) {
    #read in puck rds objects
    pat = pats[patidx]
    if (str_starts(pats[1],"puck"))
    {
      puck = readRDS(paste0(homedir,"/",pats[1],"_with_spatial.rds"))
    }
    else
    {
      system(paste0("aws s3 cp s3://snrna-seq/MBPM/newpuckdata/",pat,".rds ",homedir,"/",pat,".rds"))
      puck = readRDS(paste0(homedir,"/",pat,".rds"))
      system(paste0("rm ",homedir,"/",pat,".rds"))
    }

    #read in rctd_cell_type information, add to puck rds object
    rctd = readRDS(paste0("/data/",pat,"_rctd_main.rds"))
    with_doublets_list = as.vector(rctd@results$results_df$first_type)
    doublets_idxs = (rctd@results$results_df$spot_class %in% c("doublet_certain","doublet_uncertain"))
    with_doublets_list[doublets_idxs] = paste0(rctd@results$results_df$first_type[doublets_idxs],"_",rctd@results$results_df$second_type[doublets_idxs])
    rare_doublet_types = table(with_doublets_list)
    rare_doublet_types = rare_doublet_types[rare_doublet_types<=10]
    #with_doublets_list[with_doublets_list %in% names(rare_doublet_types)] = "Rare doublets"
    plotdf = data.frame(x = rctd@spatialRNA@coords$x, y = rctd@spatialRNA@coords$y, cell_type = with_doublets_list)
    puck$rctd_cell_type = ""
    puck$rctd_cell_type[match(rownames(rctd@results$results_df),colnames(puck))] = as.vector(rctd@results$results_df$first_type)
    puck$rctd_cell_type_with_doublets = ""
    puck$rctd_cell_type_with_doublets[match(rownames(rctd@results$results_df),colnames(puck))] = with_doublets_list
    puck = subset(puck, rctd_cell_type!="")

    #calculate signature scores for each puck, keep track of maximum and minimum score values
    min_score = 100
    max_score = -100
    for (j in 1:length(sigs_to_test)) {
      if (sum(na.omit(sigs_to_test[[j]]) %in% rownames(puck@assays$SCT)) > 0) {
	asig = sigs_to_test[[j]]
	asig = asig[asig %in% rownames(puck)]
	if (pats[1]=="puck7_20_feature_threshold" || pats[1]=="puck8_20_feature_threshold" || pats[1]=="MBM06_slide" || pats[1]=="MPM06_slide")
	{
	  puck = AddModuleScore(puck, features = list(na.omit(asig)), name = names(sigs_to_test)[j], assay = "SCT", search = T, nbin = 3)
	}
	else
	{
	  puck = AddModuleScore(puck, features = list(na.omit(asig)), name = names(sigs_to_test)[j], assay = "SCT", search = T)
	}

	current_min = quantile(puck[[paste0(names(sigs_to_test)[j],"1")]][[1]],.1)
	if (current_min < min_score)
	{
	  min_score = current_min
	}
	current_max = quantile(puck[[paste0(names(sigs_to_test)[j],"1")]][[1]],.9)
	if (current_max > max_score)
	{
	  max_score = current_max
	}
      }
    }

    #plot signature strength across the puck, using a uniform intensity scale for all signatures
    for (j in 1:length(sigs_to_test)) {
      if (sum(na.omit(sigs_to_test[[j]]) %in% rownames(puck@assays$SCT)) > 0) {
	tplot = SpatialPlot(puck, features = c(paste0(names(sigs_to_test)[j], "1")), alpha = c(0.1, 1), max.cutoff = "q90", min.cutoff = "q10", pt.size.factor = 0.4, stroke = 0.1) + scale_color_gradient(low = "grey", high = "black", limits = c(min_score, max_score)) + scale_fill_gradient(low = "grey", high = "black", limits = c(min_score, max_score)) + ggtitle(paste0(pat, " ", names(sigs_to_test)[j]))
	tplot = rasterise(tplot)
	AugmentPlot(tplot, dpi = 300)
	print(tplot)
      }
    }

    #save rds object with rctd assignments and signature strengths
    saveRDS(puck,paste0("/data/mbpm_puckdata_add_rctd_and_sigs/",pats[1],"_with_rctd_sigs.rds"))
    system(paste0("aws s3 cp /data/mbpm_puckdata_add_rctd_and_sigs/",pats[1],"_with_rctd_sigs.rds s3://uveal-melanoma/figurefolder/mbpm_puckdata_add_rctd_and_sigs/",pats[1],"_with_rctd_sigs.rds"))
    system(paste0("rm /data/mbpm_puckdata_add_rctd_and_sigs/",pats[1],"_with_rctd_sigs.rds"))
  }
  dev.off()
}
