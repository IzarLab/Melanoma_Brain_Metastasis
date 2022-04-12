show.tree = function(tree, names) {
  named.tree = tree
  for (i in 1:dim(tree)[1]) {
    for (j in 1:dim(tree)[2]) {
      if (tree[i,j] < 0) {
        t = -1*tree[i,j]
        named.tree[i,j] <- names[t]
      }
    }      
  }
  return(named.tree)
}

SelectIntegrationFeaturesWeighted = function (object.list, w= NULL, nfeatures = 800, assay = NULL, verbose = TRUE, fvf.nfeatures = 2000, ...) 
{
  if(is.null(w)) {
    w <- rep(1,length(object.list))
    names(w) <- names(object.list)
  }
  print(w)
  
  
  if (!is.null(x = assay)) {
    if (length(x = assay) != length(x = object.list)) {
      stop("If specifying the assay, please specify one assay per object in the object.list")
    }
    for (ii in length(x = object.list)) {
      DefaultAssay(object = object.list[[ii]]) <- assay[ii]
    }
  }
  else {
    assay <- sapply(X = object.list, FUN = DefaultAssay)
  }
  for (ii in 1:length(x = object.list)) {
    if (length(x = VariableFeatures(object = object.list[[ii]])) == 
        0) {
      if (verbose) {
        message(paste0("No variable features found for object", 
                       ii, " in the object.list. Running FindVariableFeatures ..."))
      }
      object.list[[ii]] <- FindVariableFeatures(object = object.list[[ii]], 
                                                nfeatures = fvf.nfeatures, verbose = verbose, 
                                                ...)
    }
  }
  var.features <- unname(obj = unlist(x = lapply(X = 1:length(x = object.list), 
                                                 FUN = function(x) VariableFeatures(object = object.list[[x]], 
                                                                                    assay = assay[x]))))
  
  var.features <- sort(x = table(var.features), decreasing = TRUE) # sort by number of datasets in which gene is present
  #var.features <- unique(var.features) # sort by number of datasets in which gene is present
  
  var.features.mtx <- matrix(nrow=length(var.features),ncol = length(object.list))
  for (i in 1:length(x = object.list)) {
    var.features.mtx[,i] <-  as.numeric(names(var.features) %in% VariableFeatures(object = object.list[[i]], assay = assay[i]))
  }
  rownames(var.features.mtx) <- names(var.features)
  var.features.mtx.score <- as.numeric(var.features.mtx %*% w)
  names(var.features.mtx.score) <- names(var.features)
  
  var.features.mtx <- var.features.mtx[order(var.features.mtx.score,decreasing = T),]
  var.features <- sort(var.features.mtx.score,decreasing = T)
  
  #keep features present in all datasets
  for (i in 1:length(x = object.list)) {
    var.features <- var.features[names(x = var.features) %in% 
                                   rownames(x = object.list[[i]][[assay[i]]])]
  }
  
  
  tie.val <- var.features[min(nfeatures, length(x = var.features))]
  features <- names(x = var.features[which(x = var.features > 
                                             tie.val)])
  if (length(x = features) > 0) {
    feature.ranks <- sapply(X = features, FUN = function(x) {
      ranks <- sapply(X = object.list, FUN = function(y) {
        vf <- VariableFeatures(object = y)
        if (x %in% vf) {
          return(which(x = x == vf)) # returns the (variable feature/vst.variance.standardized) rank for gene x in dataset y
          
        }
        #endOfRank <- length(vf)+1
        #return(endOfRank)
        return(NULL)
      })
      median(x = unlist(x = ranks)) # calculates the median of the ranks for gene x across all datasets 
      #print(paste(ranks))
      #print(median(x = unlist(x = ranks)))
      #weighted.geomean(x = unlist(x = ranks),w=w[names(unlist(x = ranks))]) # calculates the weigthed geometrical mean of the ranks for gene x across all datasets 
    })
    
    features <- names(x = sort(x = feature.ranks))
  }
  features.tie <- var.features[which(x = var.features == tie.val)]
  tie.ranks <- sapply(X = names(x = features.tie), FUN = function(x) {
    ranks <- sapply(X = object.list, FUN = function(y) {
      vf <- VariableFeatures(object = y)
      if (x %in% vf) {
        return(which(x = x == vf))
      }
      #return(nfeatures+1)
      #return(length(vf)+1)
      return(NULL)
      #endOfRank <- length(vf)+1
      #return(endOfRank)
    })
    median(x = unlist(x = ranks))
    #print(paste(ranks))
    #print(median(x = unlist(x = ranks)))
    #print(paste(length(unlist(ranks)),length(w)))
    #weighted.geomean(x = unlist(x = ranks),w=w[names(unlist(x = ranks))]) # calculates the weigthed geometrical mean of the ranks for gene x across all datasets 
    
  })
  features <- c(features, names(x = head(x = sort(x = tie.ranks), 
                                         nfeatures - length(x = features))))
  return(list(features,var.features.mtx,var.features.mtx.score))
}

### stacked violins
## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          cols=NA,
                          ...) {
  #print(cols)
  p <- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  
  if(!is.na(cols)) { p <- p + scale_fill_manual(values=cols) }
  
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


make.pca <- function(object, pca.dim=50, dataSign=1) {
  
  varfeat <- object@assays$integrated@var.features
  
  refdata <- sign(dataSign)*data.frame(t(as.matrix(object@assays$integrated@data[varfeat,])))
  refdata <- refdata[, sort(colnames(refdata))]
  
  pca <- prcomp(refdata, rank. = pca.dim, scale. = TRUE, center = TRUE, retx=TRUE)
  
  #  outliers <- apply(pca$x, 2, function(x) which( abs(x - mean(x)) > (20 * sd(x)) ))
  #  plot(pca$x[,1],pca$x[,8])
  #  points(pca$x[outliers$PC8,1], pca$x[outliers$PC8,8], col="red", pch=19)
  
  return(pca)
}


make.umap <- function(object, pca.dim=10, n.neighbors=30, min.dist=0.3, metric="cosine", n.epochs=200, umap.dim=2, seed=42, pca.obj=NULL) {
  
  require(umap)
  
  varfeat <- object@assays$integrated@var.features
  
  refdata <- data.frame(t(as.matrix(object@assays$integrated@data[varfeat,])))
  refdata <- refdata[, sort(colnames(refdata))]
  
  if (is.null(pca.obj)) {
    pca <- prcomp(refdata, rank. = pca.dim, scale. = TRUE, center = TRUE, retx=TRUE)
    # Calculate PCA for the reference map
    #NB: this is the same as pca$x <- check consistency with the formula below
    ref.pca <- scale(refdata, pca$center, pca$scale) %*% pca$rotation
  } else {
    if (pca.dim > dim(pca.obj$x)[2]) {
      stop(paste0("Error. Requested more PCA dimensions (", pca.dim ,") than available (", dim(pca.obj$x)[2],")"))
    }
    ref.pca <- pca.obj$x[,1:pca.dim]
  }
  umap.config <- umap.defaults
  umap.config$n_neighbors = n.neighbors
  umap.config$min_dist = min.dist
  umap.config$metric = metric
  umap.config$n_epochs = n.epochs
  umap.config$n_components = umap.dim
  umap.config$random_state = seed
  umap.config$transform_state = seed
  
  ref.umap <- umap(ref.pca, config=umap.config)
  if (umap.dim==2) {
    colnames(ref.umap$layout) <- c("UMAP_1","UMAP_2")
  } else if (umap.dim==3) {
    colnames(ref.umap$layout) <- c("UMAP_1","UMAP_2","UMAP_3")
  }  
  return(ref.umap) 
}

apply.pca.obj.2 <- function(query, query.assay="RNA", pca.obj) {
  
  newdata <- data.frame(t(as.matrix(query@assays[[query.assay]]@data)))
  newdata <- newdata[ , order(names(newdata))]
  
  genes.use <-  sort(intersect(colnames(newdata), names(pca.obj$center)))
  
  newdata.var <- newdata[, genes.use]
  center.use <- pca.obj$center[genes.use]
  scale.use <- pca.obj$scale[genes.use]
  rotation.use <- pca.obj$rotation[genes.use,]
  
  npca <- scale(newdata.var, center.use, scale.use) %*% rotation.use
  
  return(npca)
}
