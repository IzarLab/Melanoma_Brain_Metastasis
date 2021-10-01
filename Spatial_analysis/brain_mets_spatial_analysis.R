#!/usr/bin/env Rscript

### title: Analyze gene expression signatures and individual genes with spatial
### patterns in four slide-seq pucks author: Yiping Wang date: 09/27/2021

library(Seurat)
library(ggplot2)
library(rlist)
library(grid)

pucks = c("puck5", "puck6", "puck7", "puck8")
puck_store_folder = "/data"
filenames = c(paste0(puck_store_folder, "/puck5.rds"), paste0(puck_store_folder, 
    "/puck6final.rds"), paste0(puck_store_folder, "/puck7_20_feature_threshold.rds"), 
    paste0(puck_store_folder, "/puck8_20_feature_threshold.rds"))

pdf("puck_spatial_sigs_final.pdf")
for (i in 1:length(pucks)) {
    puck = readRDS(filenames[i])

    DefaultAssay(puck) = "SCT"

    # load in gene signatures to test on spatial data mostly testing plasma cell
    # signature, plus top-20, 50 and 100 plasma cell signatures as well as plasma
    # cell signatures without IGHG1-3 genes
    sigs_to_test = list()

    markers_bcells = read.table("markers_MBPM_sn_bcells_type.csv", sep = ",", header = T, 
        quote = "\"")
    sigs_to_test = list.append(sigs_to_test, markers_bcells$gene[markers_bcells$cluster == 
        "Plasma cells" & markers_bcells$p_val_adj < 0.05])
    sigs_to_test = list.append(sigs_to_test, markers_bcells$gene[markers_bcells$cluster == 
        "Plasma cells" & markers_bcells$p_val_adj < 0.05][1:20])
    sigs_to_test = list.append(sigs_to_test, markers_bcells$gene[markers_bcells$cluster == 
        "Plasma cells" & markers_bcells$p_val_adj < 0.05][1:50])
    sigs_to_test = list.append(sigs_to_test, markers_bcells$gene[markers_bcells$cluster == 
        "Plasma cells" & markers_bcells$p_val_adj < 0.05][1:100])
    sigs_to_test = list.append(sigs_to_test, setdiff(markers_bcells$gene[markers_bcells$cluster == 
        "Plasma cells" & markers_bcells$p_val_adj < 0.05], c("IGHG1", "IGHG2", "IGHG3")))
    sigs_to_test = list.append(sigs_to_test, setdiff(markers_bcells$gene[markers_bcells$cluster == 
        "Plasma cells" & markers_bcells$p_val_adj < 0.05][1:20], c("IGHG1", "IGHG2", 
        "IGHG3")))
    sigs_to_test = list.append(sigs_to_test, setdiff(markers_bcells$gene[markers_bcells$cluster == 
        "Plasma cells" & markers_bcells$p_val_adj < 0.05][1:50], c("IGHG1", "IGHG2", 
        "IGHG3")))
    sigs_to_test = list.append(sigs_to_test, setdiff(markers_bcells$gene[markers_bcells$cluster == 
        "Plasma cells" & markers_bcells$p_val_adj < 0.05][1:100], c("IGHG1", "IGHG2", 
        "IGHG3")))

    names(sigs_to_test) = c("plasma_cell", "plasma_cell_top20", "plasma_cell_top50", 
        "plasma_cell_top100", "plasma_cell_without_ighg123", "plasma_cell_without_ighg123_top20", 
        "plasma_cell_without_ighg123_top50", "plasma_cell_without_ighg123_top100")

    # calculate gene signature scores for each puck if signature is for top 100
    # plasma cell genes, print out plot of top 10% of cells in signature strength
    for (j in 1:length(sigs_to_test)) {
        if (sum(na.omit(sigs_to_test[[j]]) %in% rownames(puck@assays$SCT)) > 0) {
            asig = sigs_to_test[[j]]
            asig = asig[asig %in% rownames(puck)]
            if (pucks[i] == "puck7" || pucks[i] == "puck8") {
                puck = AddModuleScore(puck, features = list(na.omit(asig)), name = names(sigs_to_test)[j], 
                  assay = "SCT", search = T, nbin = 3)
            } else {
                puck = AddModuleScore(puck, features = list(na.omit(asig)), name = names(sigs_to_test)[j], 
                  assay = "SCT", search = T)
            }

            if (make_binary_plots) {
                eval(parse(text = paste0("puck$placeholder = puck$", names(sigs_to_test)[j], 
                  "1")))
                sig_cutoff = quantile(puck$placeholder, 0.9)
                if (sig_cutoff != max(puck$placeholder)) {
                  puck2 = subset(puck, placeholder > sig_cutoff)
                } else {
                  puck2 = puck
                }
                if (make_final_plots) {
                  if (names(sigs_to_test)[j] == "plasma_cell_top100" || names(sigs_to_test)[j] == 
                    "plasma_cell_without_ighg123_top100") {
                    print(SpatialPlot(puck2, features = c(paste0(names(sigs_to_test)[j], 
                      "1")), alpha = c(0.1, 1), max.cutoff = "q90", min.cutoff = "q10") + 
                      scale_color_gradient(low = "grey", high = "black") + scale_fill_gradient(low = "grey", 
                      high = "black") + ggtitle(paste0(pucks[i], " ", names(sigs_to_test)[j])))
                  }
                }
            }
        }
    }

    # print plots of top 10% highest expressing cells for various melanocyte and
    # immune cell genes
    genes_to_test = c("TIMP1", "B2M", "HLA-A", "HLA-B", "HLA-C", "FOXP3", "SOX4", 
        "CD58", "APP", "MET", "PMEL")
    DefaultAssay(puck) = "Spatial"

    for (j in 1:length(genes_to_test)) {
        puck$placeholder = t(as.matrix(puck@assays$Spatial[genes_to_test[j]]))
        sig_cutoff = quantile(puck$placeholder, 0.9)
        puck1 = subset(puck, placeholder > sig_cutoff)
        if (dim(puck1)[2] > 1) {
            print(SpatialPlot(puck1, features = c(genes_to_test[j])) + ggtitle(paste0(pucks[i], 
                " ", genes_to_test[j])))
        }
    }
}
dev.off()

# determine spatially variable genes for each puck, and store in rds files
spatialpuckfilenames = c("puck5_with_spatial.rds", "puck6final_with_spatial.rds", 
    "puck7_20_feature_threshold_with_spatial.rds", "puck8_20_feature_threshold_with_spatial.rds")
spatialpucks = list()
for (i in 1:length(pucks)) {
    puck = readRDS(filenames[i])
    DefaultAssay(puck) = "SCT"
    puck = FindSpatiallyVariableFeatures(puck, assay = "SCT", slot = "scale.data", 
        features = VariableFeatures(puck)[1:1000], selection.method = "moransi", 
        x.cuts = 100, y.cuts = 100)
    spatialpucks = list.append(spatialpucks, puck)
    saveRDS(paste0(puck_store_folder, "/", spatialpuckfilenames[i]))
    spatialgenenames = SpatiallyVariableFeatures(puck, selection.method = "moransi")
    tempdf = data.frame(gene = spatialgenenames, moransi_pval = puck[["SCT"]][["MoransI_p.value"]][spatialgenenames, 
        ])
    write.table(tempdf, paste0(pucks[i], "_spatially_variable_features.csv"), sep = ",", 
        quote = F)
}

# for puck 6, find genes that have significant Spearman correlation with TIMP1,
# with Bonferroni-corrected p-val < .05
puck6 = readRDS(paste0(puck_store_folder, "/", spatialpuckfilenames[2]))
DefaultAssay(puck6) = "SCT"
puck6 = ScaleData(puck6)
puck6_pos_corr_genes = rownames(puck6@assays$SCT@scale.data)[unlist(lapply(1:3000, 
    function(x) {
        res = cor.test(puck6@assays$SCT@scale.data[x, ], puck6@assays$SCT@scale.data["TIMP1", 
            ], method = "spearman")
        if (res$estimate > 0) {
            res$p.val
        } else {
            1
        }
    })) < 0.05/31758]
puck6_neg_corr_genes = rownames(puck6@assays$SCT@scale.data)[unlist(lapply(1:3000, 
    function(x) {
        res = cor.test(puck6@assays$SCT@scale.data[x, ], puck6@assays$SCT@scale.data["TIMP1", 
            ], method = "spearman")
        if (res$estimate < 0) {
            res$p.val
        } else {
            1
        }
    })) < 0.05/31758]

# print plots of gene expression for positive and negatively TIMP1-correlated
# genes
puck6_pos_corr_genes = intersect(SpatiallyVariableFeatures(puck6, selection.method = "moransi"), 
    puck6_pos_corr_genes)
puck6_neg_corr_genes = intersect(SpatiallyVariableFeatures(puck6, selection.method = "moransi"), 
    puck6_neg_corr_genes)
pdf("puck6_TIMP1_pos_corr_genes.pdf", width = 7, height = 35)
print(SpatialFeaturePlot(puck6, features = puck6_pos_corr_genes, ncol = 5, alpha = c(0.1, 
    1), max.cutoff = "q95"))
dev.off()
pdf("puck6_TIMP1_neg_corr_genes.pdf", width = 7, height = 35)
print(SpatialFeaturePlot(puck6, features = puck6_neg_corr_genes, ncol = 5, alpha = c(0.1, 
    1), max.cutoff = "q95"))
dev.off()

# print individual plots of genes that by manual inspection are positively or
# negatively correlated with TIMP1
manual_pos = c("APOD", "LAMA4", "USH1C", "NBEA", "ITPR2", "COL19A1", "AKAP6", "HMGA2")
manual_neg = c("CXCL9", "WARS", "B2M", "HLA-A", "IGHG1", "CD74", "IGHG3", "HLA-DPB1", 
    "IGHG2", "FTL", "CTSD", "C1QA")
pdf("puck6_pos_corr_genes_vln.pdf", width = 10, height = 7)
print(VlnPlot(puck6, features = manual_pos, group.by = "celltype_bped_main", pt.size = 0, 
    assay = "SCT", stack = T, flip = T) + NoLegend() + theme(axis.text.x = element_text(angle = 90, 
    vjust = 0.5, hjust = 1)))
dev.off()
pdf("puck6_neg_corr_genes_vln.pdf", width = 10, height = 7)
print(VlnPlot(puck6, features = manual_neg, group.by = "celltype_bped_main", pt.size = 0, 
    assay = "SCT", stack = T, flip = T) + NoLegend() + theme(axis.text.x = element_text(angle = 90, 
    vjust = 0.5, hjust = 1)))
dev.off()

for (i in 1:length(manual_pos)) {
    pdf(paste0("puck6_individual_gene_plots/puck6_", manual_pos[i], ".pdf"))
    eval(parse(text = paste0("puck6$placeholder = puck6@assays$SCT@counts[\"", manual_pos[i], 
        "\",]")))
    sig_cutoff = quantile(puck6$placeholder, 0.9)
    if (sig_cutoff != max(puck6$placeholder)) {
        puck6temp = subset(puck6, placeholder > sig_cutoff)
    } else {
        puck6temp = puck6
    }
    print(SpatialFeaturePlot(puck6temp, features = manual_pos[i], ncol = 1, alpha = c(0.1, 
        1), max.cutoff = "q95"))
    dev.off()
}

for (i in 1:length(manual_neg)) {
    pdf(paste0("puck6_individual_gene_plots/puck6_", manual_neg[i], ".pdf"))
    eval(parse(text = paste0("puck6$placeholder = puck6@assays$SCT@counts[\"", manual_neg[i], 
        "\",]")))
    sig_cutoff = quantile(puck6$placeholder, 0.9)
    if (sig_cutoff != max(puck6$placeholder)) {
        puck6temp = subset(puck6, placeholder > sig_cutoff)
    } else {
        puck6temp = puck6
    }
    print(SpatialFeaturePlot(puck6, features = manual_neg[i], ncol = 1, alpha = c(0.1, 
        1), max.cutoff = "q95"))
    dev.off()
}
