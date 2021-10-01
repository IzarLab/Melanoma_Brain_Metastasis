# Author: Somnath Tagore, Ph.D. 
# Title: DGE using deseq2
# Script Name: decseq2.run.R 
# Last Updated: 09/24/2021

# Packages required for this analysis
formatR::tidy_app()

library("DESeq2")
library(ggplot2)

# Data import
cts <- as.matrix(read.csv("countdata.tsv", sep = "\t", row.names = "gene_id"))
coldata <- read.csv("annotation.csv", row.names = 1)
coldata <- coldata[, c("condition", "type")]
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)

# Display the count data
head(cts, 2)

# Construct the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~condition)

# Pre-filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Setting the factor levels
dds$condition <- factor(dds$condition, levels = c("control", "brainmets"))

# DGE
dds <- DESeq(dds)
res <- results(dds)
res

# Results table
res <- results(dds, name = "condition_brainmets_vs_control")
write.csv(res, file = "deseq2.res.csv")

