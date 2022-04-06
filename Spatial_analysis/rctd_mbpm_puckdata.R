library(spacexr)
library(Matrix)
library(stringr)
library(Seurat)

### title: Use RCTD pipeline to assign cell type identities to SlideSeq samples
### author: Yiping Wang date: 03/29/2022

datadir = "/data"

#load in single-nuclei reference, filter for only single nuclei data, and filter out low-quality cells
system(paste0("aws s3 cp s3://snrna-seq/MBPM/data_MBPM_scn_v4.rds ",datadir,"/data_MBPM_scn_v4.rds"))
data_MBPM = readRDS(paste0(datadir,"/data_MBPM_scn_v4.rds"))
DefaultAssay(data_MBPM) = "RNA"
data_MBPM$placeholder = !(data_MBPM$cell_type_int %in% c("Low-quality cells","Doublets","Contamination","Undetermined"))
data_MBPM = subset(data_MBPM, placeholder)
data_MBPM = subset(data_MBPM, sequencing=="Single nuclei")

counts = data_MBPM@assays$RNA@counts
genenames = rownames(data_MBPM)
barcodes = colnames(data_MBPM)
rownames(counts) = genenames
colnames(counts) = barcodes

#extract cell_type_main information from single nuclei data object, create Reference object for rctd
cell_types = data_MBPM$cell_type_main
#cell_types = data_MBPM$cell_type_int
data_MBPM$cell_type_broad = "Non-immune"
data_MBPM$cell_type_broad[data_MBPM$cell_type_int %in% c("B cells","CD4+ T cells","CD8+ T cells","Dendritic cells","Mast cells","MDM","Monocytes","NK cells","Plasma cells","Tregs")] = "Immune"
data_MBPM$cell_type_broad[data_MBPM$cell_type_int %in% c("Tumor cells")] = "Tumor"
Idents(data_MBPM) = data_MBPM$cell_type_broad
#cell_types = data_MBPM$cell_type_broad
cell_types = str_replace_all(cell_types,"/","_")
names(cell_types) = barcodes
cell_types = factor(cell_types)

reference <- Reference(counts, cell_types)

pats = c("MBM05_rep1_slide","MBM06_slide","MBM07_slide","MBM08_slide","MBM11_rep1_slide","MBM18_slide","MBM13_slide","MPM08_pre_slide","MPM08_on_slide","MPM08_on_later_slide","MPM10_slide","MPM06_slide","MBM05_rep2_slide","MBM11_rep2_slide","puck5","puck6final","puck7_20_feature_threshold","puck8_20_feature_threshold")
#pats = c("puck6final","puck7_20_feature_threshold","puck8_20_feature_threshold")
#pats = c("MBM11_rep2_slide")
for (pat in pats) {
  # download rds object for each puck
  if (str_starts(pat,"puck")) {
    system(paste0("aws s3 cp s3://slide-seq/tangramrds/",pat,".rds ",datadir,"/",pat,".rds"))
  }
  else
  {
    system(paste0("aws s3 cp s3://snrna-seq/MBPM/newpuckdata/",pat,".rds ",datadir,"/",pat,".rds"))
  }
  data_pat = readRDS(paste0(datadir,"/",pat,".rds"))
  system(paste0("rm ",datadir,"/",pat,".rds"))

  #Extract coordinates and counts information for each puck
  DefaultAssay(data_pat) = "Spatial"
  barcodes = colnames(data_pat)
  genenames = rownames(data_pat)
  coords = data_pat$image@coordinates[colnames(data_pat),]
  if (pat=="puck5" || pat=="puck7_20_feature_threshold" || pat=="puck8_20_feature_threshold")
  {
    coords = coords[,c("x","y")]
  }
  if (pat=="puck6final")
  {
    #coords = data.frame(x=coords$xcoord,y=coords$ycoord)
    coords$x = coords$xcoord
    coords$y = coords$ycoord
    coords$xcoord = NULL
    coords$ycoord = NULL
  }
  counts = data_pat@assays$Spatial@counts
  rownames(counts) = genenames
  colnames(counts) = barcodes

  #run RCTD using single nuclei reference, and puck count and coordinate information
  #save in RDS file
  #output summary pdf figures
  puck <- SpatialRNA(coords, counts)
  myRCTD <- create.RCTD(puck, reference, max_cores = 8)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')

  saveRDS(myRCTD,paste0("/data/",pat,"_rctd_main.rds"))

  results <- myRCTD@results
  # normalize the cell type proportions to sum to 1.
  norm_weights = normalize_weights(results$weights) 
  cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
  spatialRNA <- myRCTD@spatialRNA
  resultsdir <- paste0('spatial/rctd_mbpm_puckdata/non_exclusive_main_sigs/',pat) ## you may change this to a more accessible directory on your computer.
  dir.create(resultsdir)

  # make the plots 
  # Plots the confident weights for each cell type as in full_mode (saved as 
  # 'results/cell_type_weights_unthreshold.pdf')
  plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights) 
  # Plots all weights for each cell type as in full_mode. (saved as 
  # 'results/cell_type_weights.pdf')
  plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights) 
  # Plots the weights for each cell type as in doublet_mode. (saved as 
  # 'results/cell_type_weights_doublets.pdf')
  plot_weights_doublet(cell_type_names, spatialRNA, resultsdir, results$weights_doublet, 
		       results$results_df) 
  # Plots the number of confident pixels of each cell type in 'full_mode'. (saved as 
  # 'results/cell_type_occur.pdf')
  plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)

  # makes a map of all cell types, (saved as 
  # 'results/all_cell_types.pdf')
  plot_all_cell_types(results$results_df, spatialRNA@coords, cell_type_names, resultsdir)

  # doublets
  #obtain a dataframe of only doublets
  doublets <- results$results_df[results$results_df$spot_class == "doublet_certain",] 
  # Plots all doublets in space (saved as 
  # 'results/all_doublets.pdf')
  plot_doublets(spatialRNA, doublets, resultsdir, cell_type_names)

  # Plots all doublets in space for each cell type (saved as 
  # 'results/all_doublets_type.pdf')
  plot_doublets_type(spatialRNA, doublets, resultsdir, cell_type_names) 
  # a table of frequency of doublet pairs 
  doub_occur <- table(doublets$second_type, doublets$first_type) 
  # Plots a stacked bar plot of doublet ocurrences (saved as 
  # 'results/doublet_stacked_bar.pdf')

  plot_doub_occur_stack(doub_occur, resultsdir, cell_type_names)
}