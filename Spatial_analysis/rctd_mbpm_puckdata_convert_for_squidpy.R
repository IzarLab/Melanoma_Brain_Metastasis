#!/usr/bin/env Rscript

### title: Load results of rctd pipe line for all spatial data, and convert their cell type assignments into txt format
### author: Yiping Wang date: 03/29/2022

pats = c("MBM05_rep1_slide","MBM06_slide","MBM07_slide","MBM08_slide","MBM11_rep1_slide","MBM18_slide","MBM13_slide","MPM08_pre_slide","MPM08_on_slide","MPM08_on_later_slide","MPM10_slide","MPM06_slide","MBM05_rep2_slide","MBM11_rep2_slide","puck5","puck6final","puck7_20_feature_threshold","puck8_20_feature_threshold")
for (pat in pats) {
  myRCTD = readRDS(paste0("/data/",pat,"_rctd_main.rds"))
  tempdf = data.frame(barcode=rownames(myRCTD@results$results_df))
  myRCTD@results$results_df = cbind(tempdf, myRCTD@results$results_df)
  write.table(myRCTD@results$results_df, paste0("/data/rctd_mbpm_puckdata_convert_for_squidpy/",pat,"_rctd_main.txt"), sep="\t", row.names=F, col.names=T, quote=FALSE)
}
system("aws s3 sync /data/rctd_mbpm_puckdata_convert_for_squidpy s3://snrna-seq/MBPM/newpuckdata/rctd_mbpm_puckdata_convert_for_squidpy")