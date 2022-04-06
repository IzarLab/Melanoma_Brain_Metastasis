#!/usr/bin/env Rscript

#### Combining counts after running Kallisto
#### Authors: Somnath Tagore, PhD; Jana Biermann, PhD

library(dplyr)
library(tximport)
library(biomaRt)
library(readr)
library(rhdf5)
library(AnnotationDbi)

sample_dir <- 'data/celllines_RNAseq'
sample_names <- grep('-', x = list.files(sample_dir), value = T)

# Create the sample metadata object
clin <- read.csv('data/celllines_RNAseq/clinical.csv')

# Specify the abundance file paths 
files <- file.path(sample_dir, clin$Sample_Name, 'abundance.h5')
names(files) <- clin$Sample_Name

# Construct the data frame to convert transcript IDs to gene symbols 
mart <- biomaRt::useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl',
                         host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c('ensembl_transcript_id', 'hgnc_symbol'),
                      mart = mart)

# Load kallisto abundance files into R and convert them
txi <- tximport(files, type = 'kallisto', tx2gene = t2g, txIdCol = 'ensembl_transcript_id',
                ignoreTxVersion = T, ignoreAfterBar = T, txOut = FALSE)

# Get count matrix and round to integer
counts <- as.data.frame(txi$counts)
counts <- counts[-which(rownames(counts) == ''), ]
counts <- round(counts, 0)
counts <- data.frame(gene = rownames(counts), counts)

# Save counts
write.csv(counts, 'data/celllines_RNAseq/data_counts_combined.csv', row.names = F)
