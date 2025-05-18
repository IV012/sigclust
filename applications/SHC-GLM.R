source("99_SigClust_GLM.R")

test_split <- function(data,ids1,ids2,var.genes,num_PCs,
                        batch,alpha_level,cores,posthoc=FALSE){
  cell1s <- data[ids1,]
  cell2s <- data[ids2,]
  gm.x <- rbind(cell1s,cell2s)
  plot(gm.x)
  batch <- batch[c(ids1,ids2)]
  label <- c(rep(1, length(ids1)), rep(2, length(ids2)))
  
  if (length(ids1) + length(ids2) <= 50){
    return(1)
  }else if (min(length(ids1)/length(ids2), length(ids2)/length(ids1)) <= 1/4){
    return(1)
  }
  else{
    gm.d <- dist(gm.x)
    gm.x <- cmdscale(gm.d, k=num_PCs)
    return(SigClust_GLM(gm.x, method="nonparametric", batch=batch, label=label)$pval)
  }
}

SHC.GLM <- function(data,batch=NULL,alpha=0.05,num_features=2500,
                  num_PCs=30,parallel=TRUE,cores=2) {
  if (!parallel) {
    cores <- 1
  }
  if (is.null(batch)) {
    batch <- rep("1",ncol(data))
  }
  if (is.factor(batch)|is.numeric(batch)) {
    warning("Converting batch vector to character")
    batch <- as.character(batch)
  }
  names(batch) <- colnames(data)
  if (is.null(colnames(data))) {
    warning("Assigning cell names automatically to the columns")
    colnames(data) <- paste0('cell',1:ncol(data))
  }
  
  # Get variable features
  num_features <- min(num_features, dim(data)[1])
  dev <- scry::devianceFeatureSelection(data)
  var.genes <- rownames(data)[order(dev,decreasing=TRUE)[1:num_features]]
  
  # Dimension reduction and clustering
  gm.x <- dev_reduce_dimension(t(data[var.genes,]),30)[[2]]
  gm.d <- dist(gm.x)
  hcl <- fastcluster::hclust(gm.d,method='ward.D')
  dend <- as.dendrogram(hcl)
  
  # Traverse the tree and store clusters and q-FWERs
  dends_to_test <- list(dend)
  clusters <- list()
  node0 <- NULL
  counter <- 0
  parents <- list('root')
  
  while (length(dends_to_test)>0) {
    # Identify the two-way split
    cuts <- dendextend::cutree(dends_to_test[[1]],k=2,order_clusters_as_data=FALSE)
    leaves <- dendextend::get_leaves_attr(dends_to_test[[1]],'label')
    ids1 <- leaves[cuts==1]
    ids2 <- leaves[cuts==2]
    alpha_level <- alpha*((length(leaves)-1)/(ncol(data)-1))
    
    # Remove any cells from batches with poor representation
    tab <- table(batch[c(ids1,ids2)],cuts)
    to.keep <- rownames(tab)[which(matrixStats::rowMins(tab)>20)]
    ids1 <- ids1[batch[ids1]%in%to.keep]
    ids2 <- ids2[batch[ids2]%in%to.keep]
    # Get p-value of the split
    if (length(to.keep)>0) {
      test <- test_split(gm.x,ids1,ids2,var.genes,num_PCs,
                         batch,alpha_level,cores,posthoc=FALSE)
      
    } else {
      test <- 1
    }
    
    # Compare to significance threshold
    if (test < alpha_level) {
      # If significant: add left and right branches to testing stack
      left <- dends_to_test[[1]][[1]]
      right <- dends_to_test[[1]][[2]]
      dends_to_test[[length(dends_to_test)+1]] <- left
      dends_to_test[[length(dends_to_test)+1]] <- right
      
      # Compute q-FWER
      if (is.null(node0)) {
        node0 <- data.tree::Node$new(paste0('Node 0: ',
                                 min(round(test*(ncol(data)-1)/(length(leaves)-1),2),1)))
      } else {
        do.call("<-",list(paste0('node',counter),
                          eval(parse(text=parents[[1]]))$AddChild(paste0('Node ',counter,': ',
                                                                         min(round(test*(ncol(data)-1)/(length(leaves)-1),2),1)))))
      }
      parents[[length(parents)+1]] <- paste0('node',counter)
      parents[[length(parents)+1]] <- paste0('node',counter)
      counter <- counter + 1
    } else {
      # If not significant: create cluster
      clusters[[length(clusters)+1]] <- leaves
      
      # Compute q-FWER
      if (is.null(node0)) {
        node0 <- data.tree::Node$new(paste0('Cluster 0: ',
                                 min(round(test*(ncol(data)-1)/(length(leaves)-1),2),1)))
      } else {
        do.call("<-",list(paste0('Cluster',length(clusters)),
                          eval(parse(text=parents[[1]]))$AddChild(paste0('Cluster ',length(clusters),': ',
                                                                         min(round(test*(ncol(data)-1)/(length(leaves)-1),2),1)))))
      }
    }
    dends_to_test[[1]] <- NULL
    parents[[1]] <- NULL
  }
  
  # Produce vector of cluster labels
  cluster_labels <- rep(0,ncol(data))
  names(cluster_labels) <- colnames(data)
  for (i in 1:length(clusters)) {
    cluster_labels[clusters[[i]]] <- i
  }
  
  return(list(cluster_labels,node0))
}