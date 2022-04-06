#!/usr/bin/env Rscript

#### Make inferCNV plots for selected samples

library(Seurat)
library(infercnv)
library(ggplot2)
library(pheatmap)
library(grid)
library(rlist)
library(stringr)

makeIndividualPlots = TRUE
makeCombinedPlot = FALSE
if (makeCombinedPlot) {
  #patslist = list(c("MBM05_sn","MBM06_sn","MBM07_sn","MBM08_sn","MBM09_sn","MBM10_sn","MBM11_sn","MBM12_sn","MBM13_sn","MBM14_sn","MBM15_sn","MBM16_sn","MBM17_sn","MBM18_sn","MBM19_sn","MBM20_sn"))
  #patslist = list(c("MBM01_sc","MBM02_sc","MBM03_sc","MBM04_sc","MBM05_sc"))
  patslist = list(c("MPM01_sn","MPM02_sn","MPM04_sn","MPM05_sn",'MPM06_sn','MPM07_sn','MPM08_sn','MPM09_sn','MPM10_sn','MPM11_sn'))
  #patslist = list(c("MBM05_sc","MBM05_sn"))
  cbpatslist = patslist
  s3folderlist = c("snrna-seq")
  cancertissues = c("malignant")
  outputnames = c("MPM_sn")
}
if (makeIndividualPlots) {
  #patslist = list(c("MBM05_sn"),c("MBM06_sn"),c("MBM07_sn"),c("MBM08_sn"),c("MBM09_sn"),c("MBM10_sn"),c("MBM11_sn"),c("MBM12_sn"),c("MBM13_sn"),c("MBM14_sn"),c("MBM15_sn"),c("MBM16_sn"),c("MBM17_sn"),c("MBM18_sn"),c("MBM19_sn"),c("MBM20_sn"))
  #patslist = list(c("MBM01_sc"),c("MBM02_sc"),c("MBM03_sc"),c("MBM04_sc"),c("MBM05_sc"))
  patslist = list(c("MPM01_sn"),c("MPM02_sn"),c("MPM04_sn"),c("MPM05_sn"),c("MPM06_sn"),c("MPM07_sn"),c("MPM08_sn"),c("MPM09_sn"),c("MPM10_sn"),c("MPM11_sn"))
  #patslist = list(c("MBM05_sc"),c("MBM05_sn"))
  cbpatslist = patslist
  s3folderlist = rep("snrna-seq",length(patslist))
  cancertissues = rep("malignant",length(patslist))
  outputnames = unlist(patslist)
}

selectstep = 1
use_one_cell_type = TRUE
one_cell_type = "malignant"
sort_by_cell_type = FALSE
use_subcluster_cell_types = FALSE
if (use_subcluster_cell_types) {
  sort_by_cell_type = TRUE
}
useRollingAverage = FALSE

# remove_ref_group = T

for (largeindex in 1:length(patslist))
{
  pats = patslist[[largeindex]]
  cbpats = cbpatslist[[largeindex]]
  s3folder = s3folderlist[largeindex]
  for (i in 1:length(pats))
  {
    pat = pats[i]
    system(paste0("aws s3 cp s3://",s3folder,"/inferCNV/",pat,"/inferCNV_subcluster_",pat,"/run.final.infercnv_obj ",pat,"_infercnv.rds --quiet"))
    infercnv_obj_temp = readRDS(paste0("",pat,"_infercnv.rds"))
    system(paste0("aws s3 cp s3://",s3folder,"/inferCNV/",pat,"/inferCNV_subcluster_",pat,"/data_",pat,"_cnv.rds ",pat,"_cb.rds --quiet"))
    seu = readRDS(paste0("",pat,"_cb.rds"))

    for (prev_name in names(infercnv_obj_temp@tumor_subclusters$subclusters))
    {
      infercnv_obj_temp@tumor_subclusters$subclusters[[prev_name]] = NULL
    }

    malignant_idxs = match(colnames(seu)[seu$malignant=="malignant"],colnames(infercnv_obj_temp@expr.data))
    malignant_idxs_names = colnames(seu)[seu$malignant=="malignant"]
    infercnv_obj_temp@tumor_subclusters$subclusters[[paste0("malignant_",pat)]] = list()
    infercnv_obj_temp@tumor_subclusters$subclusters[[paste0("malignant_",pat)]][["malignant.1"]] = malignant_idxs
    names(infercnv_obj_temp@tumor_subclusters$subclusters$malignant[["malignant.1"]]) = malignant_idxs_names

    infercnv_obj_temp@observation_grouped_cell_indices$malignant = malignant_idxs

    saveRDS(infercnv_obj_temp,paste0(pat,"_infercnv.rds"))
  }
}

for (largeindex in 1:length(patslist))
{

  pats = patslist[[largeindex]]
  cbpats = cbpatslist[[largeindex]]

  for (cancertissue in cancertissues)
  {
    one_cell_type = cancertissue
    all_gene_names = c()
    all_cell_names = c()
    for (patidx in 1:length(pats))
    {
      pat = pats[patidx]
      cbpat = cbpats[patidx]

      infercnv_obj_temp = readRDS(paste0("",pat,"_infercnv.rds"))
      colnames(infercnv_obj_temp@expr.data) = paste0(colnames(infercnv_obj_temp@expr.data),"_",pat)
      if (use_one_cell_type)
      {
	eval(parse(text=paste0("all_adipocyte_idxs = infercnv_obj_temp@observation_grouped_cell_indices$\"",one_cell_type,"\"")))
	infercnv_obj_temp@expr.data = infercnv_obj_temp@expr.data[,all_adipocyte_idxs]
	for (aname in names(infercnv_obj_temp@observation_grouped_cell_indices))
	{
	  if (aname!=one_cell_type)
	  {
	    infercnv_obj_temp@observation_grouped_cell_indices[[aname]] = NULL
	    infercnv_obj_temp@tumor_subclusters$subclusters[[aname]] = NULL
	    infercnv_obj_temp@tumor_subclusters$hc[[aname]] = NULL
	  }
	  else
	  {
	    infercnv_obj_temp@observation_grouped_cell_indices[[aname]] = NULL
	    infercnv_obj_temp@observation_grouped_cell_indices[[pat]] = 1:dim(infercnv_obj_temp@expr.data)[2]
	    for (subcluster in names(infercnv_obj_temp@tumor_subclusters$subclusters[[aname]]))
	    {
	      infercnv_obj_temp@tumor_subclusters$subclusters[[aname]][[subcluster]] = match(paste0(names(infercnv_obj_temp@tumor_subclusters$subclusters[[aname]][[subcluster]]),"_",pat),colnames(infercnv_obj_temp@expr.data))
	    }
	  }
	}
	if (length(infercnv_obj_temp@reference_grouped_cell_indices)==1)
	{
	  infercnv_obj_temp@reference_grouped_cell_indices[[1]] = c()
	}
      }

      all_gene_names = union(all_gene_names, rownames(infercnv_obj_temp@expr.data))
      all_cell_names = c(all_cell_names, colnames(infercnv_obj_temp@expr.data))

    }

    all_expr_data = matrix(0, length(all_gene_names), length(all_cell_names))
    all_gene_order = data.frame(chr = rep("",length(all_gene_names)))
    all_gene_order$start = 1
    all_gene_order$stop = 1
    all_reference_grouped_cell_indices = c()
    all_observation_grouped_cell_indices = list()
    all_cell_type_bped_main = data.frame(cell_type_bped_main = rep("",length(all_cell_names)))
    cluster_cnv_profiles = matrix(0, length(all_gene_names), 0)

    rownames(all_expr_data) = all_gene_names
    colnames(all_expr_data) = all_cell_names
    rownames(all_gene_order) = all_gene_names
    rownames(all_cell_type_bped_main) = all_cell_names
    rownames(cluster_cnv_profiles) = all_gene_names

    if (sort_by_cell_type)
    {
      for (patidx in 1:length(pats))
      {
        pat = pats[patidx]
	cbpat = cbpats[patidx]

        infercnv_obj_temp = readRDS(paste0("",pat,"_infercnv.rds"))
	if (use_one_cell_type)
	{
	  eval(parse(text=paste0("all_adipocyte_idxs = infercnv_obj_temp@observation_grouped_cell_indices$\"",one_cell_type,"\"")))
	  infercnv_obj_temp@expr.data = infercnv_obj_temp@expr.data[,all_adipocyte_idxs]
	  colnames(infercnv_obj_temp@expr.data) = paste0(colnames(infercnv_obj_temp@expr.data),"_",pat)
	  for (aname in names(infercnv_obj_temp@observation_grouped_cell_indices))
	  {
	    if (aname!=one_cell_type)
	    {
	      infercnv_obj_temp@observation_grouped_cell_indices[[aname]] = NULL
	      infercnv_obj_temp@tumor_subclusters$subclusters[[aname]] = NULL
	      infercnv_obj_temp@tumor_subclusters$hc[[aname]] = NULL
	    }
	    else
	    {
	      infercnv_obj_temp@observation_grouped_cell_indices[[aname]] = NULL
	      infercnv_obj_temp@observation_grouped_cell_indices[[pat]] = 1:dim(infercnv_obj_temp@expr.data)[2]
	      for (subcluster in names(infercnv_obj_temp@tumor_subclusters$subclusters[[aname]]))
	      {
		infercnv_obj_temp@tumor_subclusters$subclusters[[aname]][[subcluster]] = match(paste0(names(infercnv_obj_temp@tumor_subclusters$subclusters[[aname]][[subcluster]]),"_",pat),colnames(infercnv_obj_temp@expr.data))
	      }
	    }
	  }
	  if (length(infercnv_obj_temp@reference_grouped_cell_indices)==1)
	  {
	    infercnv_obj_temp@reference_grouped_cell_indices[[1]] = c()
	  }
	}
	orig_obj_temp = readRDS(paste0("",pat,"_cb.rds"))
	orig_obj_temp_colnames = paste0(colnames(orig_obj_temp),"_",pat)

	if (use_subcluster_cell_types)
	{
	  for (subcluster in names(infercnv_obj_temp@tumor_subclusters$subclusters))
	  {
	    for (subcluster1 in names(infercnv_obj_temp@tumor_subclusters$subclusters[[subcluster]]))
	    {
	      print(subcluster1)
	      if (subcluster!="tcell" && subcluster!="tcell_and_bcell")
	      {
		all_cell_type_bped_main$cell_type_bped_main[infercnv_obj_temp@tumor_subclusters$subclusters[[subcluster]][[subcluster1]]] = subcluster1
	      }
	    }
	  }
	}
	else
	{
	  matchidxs = match(orig_obj_temp_colnames,all_cell_names)
	  matchidxs = matchidxs[!is.na(matchidxs)]
	  all_cell_type_bped_main$cell_type_bped_main[matchidxs] = orig_obj_temp$celltype_bped_main
	}
	#browser()

      }

      all_cell_names = all_cell_names[order(all_cell_type_bped_main$cell_type_bped_main)]
      all_cell_type_bped_main = data.frame(cell_type_bped_main = all_cell_type_bped_main$cell_type_bped_main[order(all_cell_type_bped_main$cell_type_bped_main)])
      colnames(all_expr_data) = all_cell_names
      rownames(all_cell_type_bped_main) = all_cell_names
    }

    for (pat in pats)
    {
      infercnv_obj_temp = readRDS(paste0("",pat,"_infercnv.rds"))
      colnames(infercnv_obj_temp@expr.data) = paste0(colnames(infercnv_obj_temp@expr.data),"_",pat)

      if (use_one_cell_type)
      {
	eval(parse(text=paste0("all_adipocyte_idxs = infercnv_obj_temp@observation_grouped_cell_indices$\"",one_cell_type,"\"")))
	infercnv_obj_temp@expr.data = infercnv_obj_temp@expr.data[,all_adipocyte_idxs]
	for (aname in names(infercnv_obj_temp@observation_grouped_cell_indices))
	{
	  if (aname!=one_cell_type)
	  {
	    infercnv_obj_temp@observation_grouped_cell_indices[[aname]] = NULL
	    infercnv_obj_temp@tumor_subclusters$subclusters[[aname]] = NULL
	    infercnv_obj_temp@tumor_subclusters$hc[[aname]] = NULL
	  }
	  else
	  {
	    infercnv_obj_temp@observation_grouped_cell_indices[[aname]] = NULL
	    infercnv_obj_temp@observation_grouped_cell_indices[[pat]] = 1:dim(infercnv_obj_temp@expr.data)[2]
	    for (subcluster in names(infercnv_obj_temp@tumor_subclusters$subclusters[[aname]]))
	    {
	      infercnv_obj_temp@tumor_subclusters$subclusters[[aname]][[subcluster]] = match(paste0(names(infercnv_obj_temp@tumor_subclusters$subclusters[[aname]][[subcluster]]),"_",pat),colnames(infercnv_obj_temp@expr.data))
	    }
	  }
	}
	if (length(infercnv_obj_temp@reference_grouped_cell_indices)==1)
	{
	  infercnv_obj_temp@reference_grouped_cell_indices = list(c())
	}
      }

      if (!(use_one_cell_type) || use_subcluster_cell_types)
      {
	for (aname in names(infercnv_obj_temp@observation_grouped_cell_indices))
	{
	  infercnv_obj_temp@observation_grouped_cell_indices[[aname]] = NULL
	}
	for (a_cell_type in unique(all_cell_type_bped_main$cell_type_bped_main))
	{
	  infercnv_obj_temp@observation_grouped_cell_indices[[a_cell_type]] = match(rownames(all_cell_type_bped_main)[all_cell_type_bped_main$cell_type_bped_main==a_cell_type],colnames(infercnv_obj_temp@expr.data))
	}
      }

      all_expr_data[match(rownames(infercnv_obj_temp@expr.data),all_gene_names),match(colnames(infercnv_obj_temp@expr.data),all_cell_names)] = infercnv_obj_temp@expr.data

      all_reference_grouped_cell_indices = c(all_reference_grouped_cell_indices,match(colnames(infercnv_obj_temp@expr.data)[infercnv_obj_temp@reference_grouped_cell_indices[[1]]],all_cell_names))

      for (celltype in names(infercnv_obj_temp@observation_grouped_cell_indices))
      {
	if (!(celltype %in% names(all_observation_grouped_cell_indices)))
	{
	  eval(parse(text=paste0("all_observation_grouped_cell_indices = list.append(all_observation_grouped_cell_indices, \"",celltype,"\" = c())")))
	  print(length(all_observation_grouped_cell_indices))
	}
	eval(parse(text=paste0("an_observation_grouped_cell_indices = all_observation_grouped_cell_indices$`",celltype,"`")))
	eval(parse(text=paste0("an_observation_grouped_cell_indices_temp = infercnv_obj_temp@observation_grouped_cell_indices$`",celltype,"`")))
	an_observation_grouped_cell_indices = c(an_observation_grouped_cell_indices,match(colnames(infercnv_obj_temp@expr.data)[an_observation_grouped_cell_indices_temp],all_cell_names))
	eval(parse(text=paste0("all_observation_grouped_cell_indices[[\"",celltype,"\"]] = an_observation_grouped_cell_indices")))
      }

      all_gene_order$chr[match(rownames(infercnv_obj_temp@gene_order),all_gene_names)] = infercnv_obj_temp@gene_order$chr
      all_gene_order$start[match(rownames(infercnv_obj_temp@gene_order),all_gene_names)] = infercnv_obj_temp@gene_order$start
      all_gene_order$stop[match(rownames(infercnv_obj_temp@gene_order),all_gene_names)] = infercnv_obj_temp@gene_order$stop
    }

    all_expr_data = all_expr_data[order(all_gene_order$chr,all_gene_order$start),]
    all_gene_order = all_gene_order[order(as.numeric(all_gene_order$chr),as.numeric(all_gene_order$start)),]

    infercnv_obj = infercnv_obj_temp
    infercnv_obj@expr.data = all_expr_data
    infercnv_obj@gene_order = all_gene_order
    infercnv_obj@reference_grouped_cell_indices = list(all_reference_grouped_cell_indices)
    infercnv_obj@observation_grouped_cell_indices = all_observation_grouped_cell_indices
    infercnv_obj@tumor_subclusters = NULL

    common_gene_names = all_gene_names
    for (pat2 in pats)
    {
      infercnv_obj_temp2 = readRDS(paste0("",pat2,"_infercnv.rds"))
      common_gene_names = intersect(common_gene_names,rownames(infercnv_obj_temp2@gene_order))
    }

    infercnv_obj@expr.data = infercnv_obj@expr.data[match(common_gene_names,rownames(infercnv_obj@expr.data)),]
    infercnv_obj@gene_order = infercnv_obj@gene_order[match(common_gene_names,rownames(infercnv_obj@gene_order)),]

    selectidxs = seq(1,dim(infercnv_obj@expr.data)[2],selectstep)
    dfselect = data.frame(origidxs = selectidxs, newidxs = 1:length(selectidxs))
    infercnv_obj@expr.data = infercnv_obj@expr.data[,selectidxs]
    infercnv_obj@reference_grouped_cell_indices[[1]] = dfselect$newidxs[match(infercnv_obj@reference_grouped_cell_indices[[1]][infercnv_obj@reference_grouped_cell_indices[[1]] %in% selectidxs],dfselect$origidxs)]
    observation_grouped_cell_indices_temp = infercnv_obj@observation_grouped_cell_indices
    for (aname in names(infercnv_obj@observation_grouped_cell_indices))
    {
      testsub = dfselect$newidxs[match(infercnv_obj@observation_grouped_cell_indices[[aname]][infercnv_obj@observation_grouped_cell_indices[[aname]] %in% selectidxs],dfselect$origidxs)]
      #testsub = infercnv_obj@observation_grouped_cell_indices[[aname]][dfselect$newidxs[match(infercnv_obj@observation_grouped_cell_indices[[aname]] %in% selectidxs,dfselect$origidxs)]]
      print(length(testsub))
      print(aname)
    }
    for (aname in names(infercnv_obj@observation_grouped_cell_indices))
    {
      testsub = dfselect$newidxs[match(infercnv_obj@observation_grouped_cell_indices[[aname]][infercnv_obj@observation_grouped_cell_indices[[aname]] %in% selectidxs],dfselect$origidxs)]
      if (length(testsub)<2)
      {
	observation_grouped_cell_indices_temp[[aname]] = NULL
      }
      else
      {
	observation_grouped_cell_indices_temp[[aname]] = testsub
      }
    }
    infercnv_obj@observation_grouped_cell_indices = observation_grouped_cell_indices_temp

    obs_annotations_groups = rep(-1, length(colnames(infercnv_obj@expr.data)))
    names(obs_annotations_groups) = colnames(infercnv_obj@expr.data)
    obs_index_groupings = infercnv_obj@observation_grouped_cell_indices
    counter <- 1
    for (obs_index_group in obs_index_groupings) {
	obs_annotations_groups[ obs_index_group ] <- counter
	counter <- counter + 1
    }
    obs_annotations_groups[infercnv_obj@reference_grouped_cell_indices[[1]]] = counter+1
    infercnv_obj@observation_grouped_cell_indices$negidxs = which(obs_annotations_groups==-1)

    #fudges for single cell type plotting only
    if (use_one_cell_type)
    {
      #print(length(infercnv_obj@reference_grouped_cell_indices[[1]])==0)
      #print(sum(is.na(infercnv_obj@reference_grouped_cell_indices[[1]]))>0)
      #print(length(infercnv_obj@reference_grouped_cell_indices[[1]])==0 || sum(is.na(infercnv_obj@reference_grouped_cell_indices[[1]]))>0)
      if (length(infercnv_obj@reference_grouped_cell_indices[[1]])==0 || sum(is.na(infercnv_obj@reference_grouped_cell_indices[[1]]))>0)
      {
	infercnv_obj@reference_grouped_cell_indices[[1]] = c(1)
      }
    }
    if (length(infercnv_obj@observation_grouped_cell_indices$negidxs)==0)
    {
      infercnv_obj@observation_grouped_cell_indices$negidxs = NULL
    }

    chrs = unique(infercnv_obj@gene_order$chr)
    expr.data.temp = matrix(0,dim(infercnv_obj@expr.data)[1],dim(infercnv_obj@expr.data)[2])
    rownames(expr.data.temp) = rownames(infercnv_obj@expr.data)
    colnames(expr.data.temp) = colnames(infercnv_obj@expr.data)
    #gene_order_temp = data.frame(chr=character(), start = integer(), end = integer())
    #rownames(gene_order_temp) = rownames(infercnv_obj@gene_order)

    if (useRollingAverage)
    {
      for (z1 in 1:dim(infercnv_obj@expr.data)[2])
      {
	for (chr in chrs)
	{
	  expr.data.temp[infercnv_obj@gene_order$chr==chr,z1] = rollmean(infercnv_obj@expr.data[infercnv_obj@gene_order$chr==chr,z1],k=11,fill=NA)
	}
      }
      infercnv_obj@gene_order = infercnv_obj@gene_order[!is.na(expr.data.temp[,1]),]
      #expr.data.temp = na.omit(expr.data.temp)
      #infercnv_obj@expr.data = expr.data.temp
      infercnv_obj@expr.data = expr.data.temp[!is.na(expr.data.temp[,1]),]
      print("roll")
    }

    #lineardata = infercnv_obj@expr.data[1:(dim(infercnv_obj@expr.data)[1]*dim(infercnv_obj@expr.data)[2])]

    oldmedian = median(median(infercnv_obj@expr.data))
    infercnv_obj@expr.data[infercnv_obj@expr.data>.88 & infercnv_obj@expr.data<1.12] = oldmedian

    # if (remove_ref_group)
    # {
    #   if (length(infercnv_obj@reference_grouped_cell_indices[[1]])>=10)
    #   {
    #     infercnv_obj@reference_grouped_cell_indices[[1]] = infercnv_obj@reference_grouped_cell_indices[[1]][1:10]
    #   }
    # }

    #nonsense = nonsense+1
    source("misc/plot_cnv_fullcode.R")
    plot_cnv(infercnv_obj)
    if (use_one_cell_type)
    {
      system(paste0("mv infercnv.pdf ",outputnames[largeindex],"_",str_replace_all(one_cell_type," ","_"),"_infercnv.pdf"))
    }
    else
    {
      system(paste0("mv infercnv.pdf ",outputnames[largeindex],"_infercnv.pdf"))
    }
  }
}

