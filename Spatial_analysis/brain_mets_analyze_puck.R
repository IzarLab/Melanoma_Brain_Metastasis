#!/usr/bin/env Rscript

### title: Load results of slideseq-tools pipeline for pucks 5-8, and convert to
### Seurat rds format author: Yiping Wang date: 09/27/2021

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SingleR)
library(SingleCellExperiment)
library(scater)
library(pheatmap)
library(gplots)

pats = c("puck5", "puck6", "puck7", "puck8")
puck_nrs = c("puck5", "puck6", "puck7", "puck8")
capitalpucks = c("Puck5", "Puck6", "Puck7", "Puck8")
for (i in 1:length(pats)) {
    pat <- pats[i]
    puck_nr <- puck_nrs[i]
    capitalpuck <- capitalpucks[i]

    # read in puck data processed by slideseq-tools pipeline
    expr <- data.table::fread(input = paste0(capitalpuck, ".AllIllumina.digital_expression.txt"), 
        data.table = FALSE)

    rownames(expr) <- expr[, 1]
    expr <- expr[, -1]

    # load number of reads per cell, calculate saturation per cell
    umi <- as.data.frame(colSums(expr))
    colnames(umi) <- "umi"
    umi$barcode <- colnames(expr)
    reads <- read.delim(paste0(capitalpuck, ".numReads_perCell_XC_mq_10.txt"))
    reads <- reads[, 1:2]
    colnames(reads) <- c("nReads", "barcode")
    umi <- left_join(umi, reads, "barcode")
    umi$saturation <- 1 - (umi$umi/umi$nReads)

    # read in spatial positions
    positions <- read.csv(paste0(pat, "_BeadLocationsForR_rough.csv"))
    write.table(colnames(expr), paste0("Illumina_barcode_", pat, ".txt"), sep = "\t", 
        quote = F, row.names = F, col.names = F)
    rownames(positions) = positions$barcodes
    write.table(data.frame(rownames(positions), positions$xcoord, positions$ycoord), 
        paste0("bead_barcode_", pat, ".txt"), sep = "\t", quote = F, row.names = F, 
        col.names = F)
    source(paste0("cmatcher bead_barcode_", pat, ".txt Illumina_barcode_", pat, ".txt cmatcher_", 
        pat, ".txt cmatcher_details_", pat, ".txt 180402 1"))

    details_test = read.table(paste0("/mnt/vdb/home/ubuntu2/cmatcher_details_", puck_nr, 
        ".txt"), header = T, sep = "\t", quote = NULL, fill = T)
    details_test$IlluminaBarcodes = substr(details_test$IlluminaBarcodes, 1, 14)
    details_test = details_test[match(unique(details_test$IlluminaBarcodes), details_test$IlluminaBarcodes), 
        ]

    colnames(expr) = substr(colnames(expr), 1, 14)
    umi$barcode = substr(umi$barcode, 1, 14)

    exprmatched = expr[, na.omit(match(details_test$IlluminaBarcodes, colnames(expr)))]
    umimatched = umi[na.omit(match(details_test$IlluminaBarcodes, colnames(expr))), 
        ]

    # create seurat obj
    puck <- CreateSeuratObject(counts = exprmatched, project = puck_nr, assay = "Spatial")
    puck[["patient"]] <- pat
    puck[["puck"]] <- puck_nr

    # add positional info
    positions = details_test
    rownames(x = positions) <- positions$IlluminaBarcodes
    positions <- positions[, 5:6]
    positions = positions[match(colnames(exprmatched), rownames(positions)), ]
    puck[["image"]] <- new(Class = "SlideSeq", assay = "Spatial", coordinates = positions)

    puck$log_nCount_Spatial <- log(puck$nCount_Spatial)

    # standard workflow
    puck <- SCTransform(puck, assay = "Spatial", ncells = 3000, verbose = FALSE)
    puck <- RunPCA(puck)
    puck <- RunUMAP(puck, dims = 1:30)
    puck <- FindNeighbors(puck, dims = 1:30)
    puck <- FindClusters(puck, resolution = 0.3, verbose = FALSE)

    # SingleR
    puck_sce <- as.SingleCellExperiment(puck)
    bped <- BlueprintEncodeData()
    pred_bped_fine <- SingleR(test = puck_sce, ref = bped, labels = bped$label.fine)
    pruneScores(pred_bped_fine)
    puck[["celltype_bped_fine"]] <- pred_bped_fine$pruned.labels
    pred_bped_main <- SingleR(test = puck_sce, ref = bped, labels = bped$label.main)
    pruneScores(pred_bped_main)
    puck[["celltype_bped_main"]] <- pred_bped_main$pruned.labels

    saveRDS(puck, paste0(pat, ".rds"))

    puck[["nReads"]] <- umimatched$nReads
    puck[["saturation"]] <- umimatched$saturation

    stats <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 10))
    colnames(stats) <- c("sample", "puck", "saturation", "n_reads", "n_UMIs", "n_features", 
        "n_cells", "median_reads", "median_features", "median_counts")
    rownames(stats) <- pat
    stats$sample <- pat
    stats$puck <- puck_nr
    stats$saturation <- round(mean(puck@meta.data$saturation), digits = 2)
    stats$n_reads <- format(sum(puck@meta.data$nReads), big.mark = ",")
    stats$n_UMIs <- format(sum(puck@meta.data$nCount_Spatial), big.mark = ",")
    stats$n_features <- format(dim(puck@assays$Spatial@counts)[1], big.mark = ",")
    stats$n_cells <- format(dim(puck@assays$Spatial@counts)[2], big.mark = ",")
    stats$median_reads <- format(round(median(puck@meta.data$nReads)), big.mark = ",")
    stats$median_features <- round(median(puck@meta.data$nFeature_Spatial))
    stats$median_counts <- round(median(puck@meta.data$nCount_Spatial))

    ## plots
    pdf(file = paste0(puck_nr, ".pdf"))
    textplot(t(stats), cex = 1.2, halign = "left")
    print(VlnPlot(puck, features = c("nFeature_Spatial", "nCount_Spatial"), pt.size = 0, 
        log = TRUE, ncol = 2, split.by = NULL))
    print(VlnPlot(puck, features = "nReads", log = T))
    print(VlnPlot(puck, features = "saturation"))
    print(SpatialFeaturePlot(puck, features = "log_nCount_Spatial") + theme(legend.position = "right"))
    Idents(puck) <- puck[["seurat_clusters"]]

    print(DimPlot(puck, reduction = "umap", label = TRUE))

    print(SpatialDimPlot(puck, stroke = 0))

    # SingleR
    print(plotScoreHeatmap(pred_bped_main, clusters = puck@meta.data$seurat_clusters, 
        fontsize = 6, main = "pred_bped_main"))
    print(DimPlot(puck, reduction = "umap", label = T, group.by = "celltype_bped_main", 
        repel = T, label.size = 2.5) + ggtitle("Cell type identification using SingleR (celltype_bped_main)") + 
        guides(col = guide_legend(nrow = 30, override.aes = list(size = 5))) + theme(legend.text = element_text(size = 6)))
    Idents(puck) <- puck[["celltype_bped_main"]]
    dev.off()
    print(SpatialDimPlot(puck, stroke = 0) + ggtitle("Cell type identification using SingleR (celltype_bped_main)") + 
        guides(col = guide_legend(nrow = 30, override.aes = list(size = 5))) + theme(legend.text = element_text(size = 6)))
    SpatialFeaturePlot(puck, features = c("MLANA", "PTPRC", "CD4", "CD8A"), alpha = c(0.1, 
        1))
    SpatialFeaturePlot(puck, features = c("MKI67", "CD19", "CD68", "CD79A"), alpha = c(0.1, 
        1))

    dev.off()
}
