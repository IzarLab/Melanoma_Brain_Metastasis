#!/usr/bin/env Rscript

### title: Compare copy number alterations among TCGA, Davies, and MBPM datasets
### author: Yiping Wang date: 09/27/2021

library(copynumber)
library(stringr)
library(infercnv)
library(grid)
library(rlist)
library(matrixStats)
library(dplyr)

tcga_clinical = read.table("gdc_download_skcm_clinical.tsv", sep = "\t", header = T, 
    quote = NULL)
tcga_aliquot = read.table("gdc_download_skcm_aliquot.tsv", sep = "\t", header = T, 
    quote = NULL)
davies_keyfile = read.table("Davies_processed_WES/key.file.to.share.txt", sep = "\t", 
    header = T, quote = "\"")
human_genome_length = 3.2e+09

# create dataframe of parameters that indicate which dataset to process in
# comparing cnv's
parameter_df = data.frame(useDavies = F, useTCGA = F, useMBPM = T, useMBPM_BM = T, 
    useMBPM_PM = F, useDavies_BM = F, useDavies_ECM = F)
tempdf = data.frame(useDavies = T, useTCGA = F, useMBPM = F, useMBPM_BM = F, useMBPM_PM = F, 
    useDavies_BM = T, useDavies_ECM = F)
parameter_df = rbind(parameter_df, tempdf)
tempdf = data.frame(useDavies = T, useTCGA = F, useMBPM = F, useMBPM_BM = F, useMBPM_PM = F, 
    useDavies_BM = F, useDavies_ECM = T)
parameter_df = rbind(parameter_df, tempdf)
tempdf = data.frame(useDavies = F, useTCGA = T, useMBPM = F, useMBPM_BM = F, useMBPM_PM = F, 
    useDavies_BM = F, useDavies_ECM = F)
parameter_df = rbind(parameter_df, tempdf)

Davies_fraction_altered_arr_BM = c()
Davies_fraction_altered_arr_ECM = c()

pdf("cnv_freq_across_samples.pdf", width = 21, height = 7)

for (z in 1:dim(parameter_df)[1]) {
    useDavies = parameter_df$useDavies[z]
    useTCGA = parameter_df$useTCGA[z]
    useMBPM = parameter_df$use_MBPM[z]
    useDavies_BM = parameter_df$useDavies_BM[z]
    useDavies_ECM = parameter_df$useDavies_ECM[z]
    useMBPM_BM = parameter_df$useMBPM_BM[z]
    useMBPM_PM = parameter_df$useMBPM_PM[z]

    # define folder location of .seg files, and tags to use in naming output pdf
    # files
    if (useDavies) {
        segfolder = "Davies_processed_WES"
        firsttitle = "davies"
        if (useDavies_BM) {
            secondtitle = "bm"
        } else if (useDavies_ECM) {
            secondtitle = "ecm"
        }
    } else if (useTCGA) {
        segfolder = "gdc_download_skcm_wes"
        firsttitle = "tcga"
        secondtitle = "ecm"
    } else if (useMBPM) {
        segfolder = "/data"
        firsttitle = "mbpm"
        if (useMBPM_BM) {
            secondtitle = "bm"
        } else if (useMBPM_PM) {
            secondtitle = "pm"
        }
    }

    seg_subfolders_or_files = list.files(segfolder)
    # in case of using MBPM data, define names of samples as stored on Amazon s3
    if (useMBPM) {
        if (useMBPM_BM) {
            seg_subfolders_or_files = c("patient5_sn", "bi006-skcm-snseq", "bi007-skcm-snseq", 
                "ns001-skcm-snseq", "ns002-skcm-snseq", "ns003-skcm-snseq", "ns006-skcm-snseq", 
                "ns008-skcm-snseq", "pa021-skcm-snseq", "pa045-skcm-snseq", "pa097-skcm-snseq", 
                "pa106-skcm-snseq", "pa123-skcm-snseq", "pa137-skcm-snseq", "pa144-skcm-snseq", 
                "pa145-skcm-snseq", "pa162-skcm-snseq")
        } else if (useMBPM_PM) {
            seg_subfolders_or_files = c("MPM01_sn", "MPM02_sn", "MPM04_sn", "MPM05_sn")
        }
    }

    all_samples_table = data.frame(sampleID = character(), chrom = integer(), start.pos = integer(), 
        end.pos = integer(), mean = double())
    for (i in 1:length(seg_subfolders_or_files)) {

        # extract sample name and phenotype (BM or ECM) from Davies file name
        if (useDavies) {
            afile = seg_subfolders_or_files[i]
            titleid_startidx = str_locate(afile, "ExomeCN.")[2] + 1
            titleid_endidx = str_locate(afile, ".seg")[1] - 1
            titleid = substr(afile, titleid_startidx, titleid_endidx)
            phenotype = davies_keyfile$phenotype[davies_keyfile$title == titleid]
        }

        # extract sample name and sample site (BM or ECM) from files under each tcga
        # subfolder
        if (useTCGA) {
            if (seg_subfolders_or_files[i] != "MANIFEST.txt" && seg_subfolders_or_files[i] != 
                "gdc_download_20210824_181450.338478.tar.gz") {
                afiles = list.files(paste0(segfolder, "/", seg_subfolders_or_files[i]))
                afile = ""
                for (j in 1:length(afiles)) {
                  if (length(grep("seg.txt", afiles[j])) != 0) {
                    afile = afiles[j]
                  }
                }
                aliquotid_startidx = str_locate(afile, "TCGA-SKCM.")[2] + 1
                aliquotid_endidx = str_locate(afile, ".ascat2")[1] - 1
                aliquotid = substr(afile, aliquotid_startidx, aliquotid_endidx)
                case_id = tcga_aliquot$case_id[tcga_aliquot$aliquot_id == aliquotid]
                sample_site = unique(tcga_clinical$site_of_resection_or_biopsy[tcga_clinical$case_id == 
                  case_id])
            }
        }

        # download infercnv result files for each file in MBPM dataset
        if (useMBPM) {
            if (useMBPM_BM) {
                system(paste0("aws s3 cp s3://snrna-seq/inferCNV/", seg_subfolders_or_files[i], 
                  "/inferCNV_cellbender_final_", seg_subfolders_or_files[i], "/run.final.infercnv_obj ", 
                  segfolder, "/infercnv_obj.rds"))
            } else if (useMBPM_PM) {
                system(paste0("aws s3 cp s3://snrna-seq/inferCNV/", seg_subfolders_or_files[i], 
                  "/inferCNV_final_", seg_subfolders_or_files[i], "/run.final.infercnv_obj ", 
                  segfolder, "/infercnv_obj.rds"))
            }
            infercnv_obj = readRDS(paste0(segfolder, "/infercnv_obj.rds"))
            afile = "infercnv_obj.rds"
        }

        # readin and append seg file data for each sample to all_samples_table for Davies
        # dataset, calculate fraction of genome with copy number < -.4 or > .4, and store
        # in array
        if (afile != "" && afile != "key.file.to.share.txt" && seg_subfolders_or_files[i] != 
            "MANIFEST.txt" && seg_subfolders_or_files[i] != "gdc_download_20210824_181450.338478.tar.gz") {
            if (useDavies) {
                if ((useDavies_ECM && phenotype == "ECM") || (useDavies_BM && phenotype == 
                  "BM")) {
                  atable = read.table(paste0(segfolder, "/", afile), header = T, 
                    sep = "\t", quote = NULL)
                  atable$sampleID = titleid
                  atable$chrom2 = atable$chrom
                  atable$chrom = NULL
                  atable$chrom = atable$chrom2
                  atable$chrom2 = NULL
                  atable$start.pos = atable$loc.start
                  atable$end.pos = atable$loc.end
                  atable$mean = atable$seg.mean
                  atable = subset(atable, chrom != "X")
                  atable = subset(atable, chrom != "Y")
                  atable$loc.start = NULL
                  atable$loc.end = NULL
                  atable$num.mark = NULL
                  atable$seg.mean = NULL
                  atable$ID = NULL
                  all_samples_table = rbind(all_samples_table, atable)

                  all_bases_altered = sum((atable$end.pos - atable$start.pos)[abs(atable$mean) > 
                    0.4])
                  fraction_altered = all_bases_altered/human_genome_length
                  if (useDavies_BM) {
                    Davies_fraction_altered_arr_BM = append(Davies_fraction_altered_arr_BM, 
                      fraction_altered)
                  } else {
                    Davies_fraction_altered_arr_ECM = append(Davies_fraction_altered_arr_ECM, 
                      fraction_altered)
                  }
                }
            } else if (useTCGA) {
                if (sample_site != "Brain, NOS" && length(grep("Skin", sample_site)) == 
                  0) {
                  atable = read.table(paste0(segfolder, "/", seg_subfolders_or_files[i], 
                    "/", afile), header = T, sep = "\t", quote = NULL)
                  atable$sampleID = atable$GDC_Aliquot
                  atable$chrom = substring(atable$Chromosome, 4)
                  atable$start.pos = atable$Start
                  atable$end.pos = atable$End
                  atable$mean = atable$Copy_Number
                  atable = subset(atable, chrom != "X")
                  atable = subset(atable, chrom != "Y")
                  atable$chrom = as.integer(atable$chrom)
                  atable$GDC_Aliquot = NULL
                  atable$Chromosome = NULL
                  atable$Start = NULL
                  atable$End = NULL
                  atable$Copy_Number = NULL
                  atable$Major_Copy_Number = NULL
                  atable$Minor_Copy_Number = NULL
                  all_samples_table = rbind(all_samples_table, atable)
                }
            } else if (useMBPM) {
                rowmediansarr = rowMedians(infercnv_obj@expr.data[, infercnv_obj@observation_grouped_cell_indices$tumor])
                atable = data.frame(sampleID = rep(seg_subfolders_or_files[i], length(rowmediansarr)), 
                  chrom = substring(infercnv_obj@gene_order$chr, 4), start.pos = infercnv_obj@gene_order$start, 
                  end.pos = infercnv_obj@gene_order$stop, mean = rowmediansarr)
                all_samples_table = rbind(all_samples_table, atable)
            }
        }

        if (useMBPM) {
            system(paste0("rm ", segfolder, "/infercnv_obj.rds"))
        }
    }

    # set cutoffs for considering whether a segment is copy number altered or not for
    # MBPM data, take cutoffs to exclude copy number that occurs less than 100 times
    # in data
    if (useDavies) {
        cutoffs = list(c(-0.1, 0.1), c(-0.2, 0.2), c(-0.3, 0.3), c(-0.4, 0.4), c(-0.5, 
            0.5))
    } else if (useTCGA) {
        cutoffs = list(c(1.9, 4.1))
    } else if (useMBPM) {
        all_samples_tablemean = table(all_samples_table$mean)
        all_cnvs = as.double(names(all_samples_tablemean))
        common_cnvs = as.double(names(all_samples_tablemean[all_samples_tablemean > 
            100]))
        posdists = all_cnvs - common_cnvs[2]
        posdists = posdists[posdists > 0]
        negdists = all_cnvs - common_cnvs[1]
        negdists = abs(negdists[negdists < 0])
        leastdist = min(c(posdists, negdists))
        cutoffs = list(c(min(common_cnvs) - leastdist/2, max(common_cnvs) + leastdist/2))
    }

    for (i in 1:length(cutoffs)) {
        gain_cutoff = cutoffs[[i]][2]
        loss_cutoff = cutoffs[[i]][1]
        print(plotFreq(segments = all_samples_table, thres.gain = gain_cutoff, thres.loss = loss_cutoff, 
            ylim = c(-100, 100), main = paste0(firsttitle, "_cnv_freq_", secondtitle, 
                " cutoffs: ", loss_cutoff, " ", gain_cutoff)))
    }
}

dev.off()

# create boxplot of altered genome fraction in Davies BM vs. ECM, at -0.4 to 0.4
# threshold
pdf("Davies_BM_vs_ECM_threshold_4_unpaired_same_shapes.pdf")
BM_arr = Davies_fraction_altered_arr_BM
ECM_arr = Davies_fraction_altered_arr_ECM

anunpairedpval = wilcox.test(BM_arr, ECM_arr)$p.val

unpaireddf1 = data.frame(fraction = BM_arr, sample_type = "BM")
unpaireddf2 = data.frame(fraction = ECM_arr, sample_type = "ECM")
unpaireddf = rbind(unpaireddf1, unpaireddf2)
print(unpaireddf %>% ggplot(aes(sample_type, fraction, fill = sample_type)) + geom_boxplot(aes(fill = sample_type)) + 
    geom_point(aes(fill = sample_type), shape = 21, size = 2) + ylim(0, 1) + ggtitle(paste0("Unpaired comparison\nDavies WES CNV Threshold of ", 
    altered_thresholds_arr[i], "\nWilcoxon unpaired p-val: ", anunpairedpval)) + 
    theme_set(theme_bw()))

dev.off()
