library(parallel)

# Compute Poisson deviances
poisson_dev_batch <- function(y,x) {
  if (is.null(x)) {
    n <- Matrix::colSums(y)
    pis <- Matrix::rowSums(y)/sum(y)
    mu <- crossprod(array(pis,dim=c(1,length(pis))),array(n,dim=c(1,length(n))))
    d <- 2 * (y * log(ifelse(y == 0, 1, y/mu)) - (y - mu))
    d[d<0] <- 0
    
    return(sqrt(d)*ifelse(y>mu,1,-1))
  } else {
    y1 <- lapply(unique(x),function(i) y[,x==i,drop=F])
    n <- lapply(y1,Matrix::colSums)
    pis <- lapply(y1,function(data) Matrix::rowSums(data)/sum(data))
    mu <- lapply(1:length(y1),function(ind)
      crossprod(array(pis[[ind]],dim=c(1,length(pis[[ind]]))),
                array(n[[ind]],dim=c(1,length(n[[ind]])))))
    d <- lapply(1:length(y1),function(ind)
      2 * (y1[[ind]] * log(ifelse(y1[[ind]] == 0, 1, y1[[ind]]/mu[[ind]])) -
             (y1[[ind]] - mu[[ind]])))
    
    res <- array(0,dim=dim(y))
    rownames(res) <- rownames(y)
    colnames(res) <- colnames(y)
    for (ind in 1:length(y1)) {
      d[[ind]][d[[ind]]<0] <- 0
      res[,x==unique(x)[ind]] <- as.matrix(sqrt(d[[ind]])*
                                             ifelse(y1[[ind]]>mu[[ind]],1,-1))
    }
    
    return(res)
  }
}

# Compute Poisson dispersion statistics
poisson_dispersion_stats <- function(y){
  n <- Matrix::colSums(y)
  pis <- Matrix::rowSums(y)/sum(y)
  mu <- crossprod(array(pis,dim=c(1,length(pis))),array(n,dim=c(1,length(n))))
  y2 <- (y - mu)^2 / mu
  
  disp <- Matrix::rowSums(y2)/ncol(y2)
  
  if (!'matrix'%in%class(y2)) {
    y2 <- as.matrix(y2)
  }
  
  return(sqrt(ncol(y))*(disp-1)/sqrt(matrixStats::rowVars(y2)))
}

# Perform dimension reduction
reduce_dimension <- function(y,x,num_PCs) {
  pdev <- poisson_dev_batch(y,x)
  pdev <- t(scale(Matrix::t(pdev),scale=F))
  PCs <- RSpectra::eigs_sym(as.matrix(tcrossprod(pdev)),k=num_PCs)
  projection <- t(crossprod(PCs$vectors,pdev))
  
  return(list(PCs, projection))
}

# Compute expected sum of squares from dimension reduction scores
compute_ess <- function(redduc) {
  sum((rowSums(sweep(redduc,2,colMeans(redduc),'-')^2)))
}

# Compute test statistic
ward_linkage <- function(redduc,labels) {
  ess1 <- compute_ess(redduc[labels==1,])
  ess2 <- compute_ess(redduc[labels==2,])
  ess <- compute_ess(redduc)
  return((ess-(ess1+ess2))/length(labels))
}


# Fit model for one batch
fit_model_batch <- function(y,on_genes,num_PCs) {
  # Compute sample moments of the on genes
  on_counts <- Matrix::t(y[on_genes,])
  cov <- cov(as.matrix(on_counts))
  means <- Matrix::colMeans(on_counts)
  
  # Use method of moments to estimate parameters
  sigmas <- log(((diag(cov)-means)/means^2)+1)
  mus <- log(means)-0.5*sigmas
  mus.sum <- tcrossprod(array(mus,dim=c(length(mus),1)),
                        array(1,dim=c(length(mus),1)))+
    tcrossprod(array(1,dim=c(length(mus),1)),array(mus,dim=c(length(mus),1)))
  sigmas.sum <- tcrossprod(array(sigmas,dim=c(length(sigmas),1)),
                           array(1,dim=c(length(sigmas),1)))+
    tcrossprod(array(1,dim=c(length(sigmas),1)),
               array(sigmas,dim=c(length(sigmas),1)))
  rhos <- suppressWarnings(log(cov/(exp(mus.sum+0.5*sigmas.sum))+1))
  rhos[is.na(rhos)|is.infinite(rhos)] <- (-10)
  diag(rhos) <- sigmas
  
  # Make the covariance matrix positive-definite
  on_cov_eigs <- RSpectra::eigs_sym(as.matrix(rhos),k=min(c(nrow(rhos)-1,num_PCs)))
  num_pos <- sum(on_cov_eigs$values>0)
  on_cov.sub <- on_cov_eigs$vectors[,1:num_pos]%*%
    sqrt(diag(on_cov_eigs$values[1:num_pos]))
  on_cov <- tcrossprod(on_cov.sub)
  diag(on_cov) <- diag(rhos)
  on_cov <- sfsmisc::posdefify(on_cov)
  on_cov.sqrt <- t(chol(on_cov))
  
  return(list(Matrix::rowMeans(y),mus,on_cov.sqrt))
}

# Fit model
fit_model <- function(y,on_genes,x,num_PCs) {
  on_means <- list()
  on_cov.sqrt <- list()
  lambdas <- list()
  
  for (b in unique(x)) {
    params <- fit_model_batch(y[,x==b],on_genes=on_genes,num_PCs=num_PCs)
    lambdas[[as.character(b)]] <- params[[1]]
    on_means[[as.character(b)]] <- params[[2]]
    on_cov.sqrt[[as.character(b)]] <- params[[3]]
  }
  
  return(list(lambdas,on_means,on_cov.sqrt))
}

# Generate one null sample
generate_null <- function(y,params,on_genes,x) {
  lambdas <- params[[1]]
  on_means <- params[[2]]
  on_cov.sqrt <- params[[3]]
  
  null <- array(0,dim=dim(y))
  rownames(null) <- rownames(y)
  
  for (b in unique(x)) {
    num_gen <- min(sum(x==b),1000)
    names(lambdas[[as.character(b)]]) <- rownames(y)
    null[-on_genes,which(x==b)[1:num_gen]] <-
      array(rpois(num_gen*(nrow(null)-length(on_genes)),
                  lambdas[[as.character(b)]][-on_genes]),
            dim=c(nrow(null)-length(on_genes),num_gen))
    Y <- exp(sweep(on_cov.sqrt[[as.character(b)]]%*%
                     array(rnorm(num_gen*length(on_genes)),
                           dim=c(length(on_genes),num_gen)),1,
                   on_means[[as.character(b)]],'+'))
    null[on_genes,which(x==b)[1:num_gen]] <-
      array(rpois(length(Y),Y),dim=dim(Y))
  }
  
  return(list(null[,colSums(null)>0],x[colSums(null)>0]))
}

# Generate one null sample and compute test statistic
generate_null_statistic <- function(y,params,on_genes,x,num_PCs,
                                    gm,labs,posthoc) {
  null_set <- generate_null(y,params,on_genes,x)
  null <- null_set[[1]]
  batch_null <- null_set[[2]]
  
  if (!posthoc) {
    null_gm <- reduce_dimension(null,batch_null,num_PCs)[[2]]
    null_gm.d <- dist(null_gm)
    hc2 <- cutree(hclust(null_gm.d,method='ward.D'),2)
  } else {
    null_gm <- reduce_dimension(null,batch_null,num_PCs)[[2]]
    pdev <- poisson_dev_batch(null,batch_null)
    pdev <- t(scale(Matrix::t(pdev),scale=F))
    gm2 <- t(crossprod(gm[[1]]$vectors,pdev))
    knns <- queryKNN(gm[[2]], gm2,
                     k=15, BNPARAM=AnnoyParam())$index
    neighbor.labels <- apply(knns,1,function(x) labs[x])
    hc2 <- unlist(apply(neighbor.labels,2,function(x)
      sample(names(table(x))[as.vector(table(x))==max(table(x))],1)))
    
    if (length(unique(hc2))==1) {
      null_gm.d <- dist(null_gm)
      hc2 <- cutree(hclust(null_gm.d,method='ward.D'),2)
    }
  }
  
  Qclust2 <- sapply(unique(batch_null),function(b)
    if (length(unique(hc2[batch_null==b]))==2 &
        min(table(hc2[batch_null==b]))>=2) {
      ward_linkage(null_gm[batch_null==b,],hc2[batch_null==b])
    } else {
      0
    })
  names(Qclust2) <- unique(batch_null)
  return(median(Qclust2))
}

# Test one split
test_split <- function(data,ids1,ids2,num_PCs,batch,
                       alpha_level,cores,posthoc) {
  # Re-order data
  cell1s <- data[,ids1]
  cell2s <- data[,ids2]
  true <- cbind(cell1s,cell2s)
  batch <- batch[c(ids1,ids2)]
  
  # Re-run dimension reduction and calculate test statistic
  gm <- reduce_dimension(true,batch,num_PCs)
  gm_sub.x <- gm[[2]]
  labs <- c(rep(1,length(ids1)),rep(2,length(ids2)))
  Qclust <- sapply(unique(batch),function(b)
    ward_linkage(gm_sub.x[batch==b,],labs[batch==b]))
  names(Qclust) <- unique(batch)
  stat <- median(Qclust)
  
  # Determine the set of "on" genes
  phi_stat <- poisson_dispersion_stats(true)
  check_means <- matrixStats::rowMins(sapply(unique(batch),function(b)
    Matrix::rowSums(true[,batch==b])))
  on_genes <- which(pnorm(phi_stat,lower.tail=F)<0.05&check_means!=0)
  
  # Fit model
  params <- fit_model(true,on_genes,batch,num_PCs)
  
  # Generate null distribution of test statistics
  Qclusts2_1 <- mclapply(1:10,function(i) {
    generate_null_statistic(true,params,on_genes,batch,
                            num_PCs,gm,labs,posthoc)
  },mc.cores=cores)
  
  # Quit early if p-value is much smaller or much larger than alpha level
  Qclusts2 <- unlist(Qclusts2_1)
  fit <- MASS::fitdistr(Qclusts2,'normal')
  pval <- 1-pnorm(stat,mean=fit$estimate[1],sd=fit$estimate[2])
  if (pval < 0.1*alpha_level | pval > 10*alpha_level) {
    return(pval)
  }
  
  # Otherwise, keep going
  Qclusts2_2 <- mclapply(11:50,function(i) {
    generate_null_statistic(true,params,on_genes,
                            batch,num_PCs,gm,labs,posthoc)
  },mc.cores=cores)
  
  # Compute smoothed p-value
  Qclusts2 <- c(Qclusts2,unlist(Qclusts2_2))
  fit <- fitdistr(Qclusts2,'normal')
  return(1-pnorm(stat,mean=fit$estimate[1],sd=fit$estimate[2]))
}

# Full clustering pipeline with built-in hypothesis testing
scSHC <- function(data,batch=NULL,alpha=0.05,num_features=2500,
                  num_PCs=30,parallel=T,cores=2) {
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
    colnames(data) <- paste0('cell',1:ncol(data))
  }
  
  # Dimension reduction and clustering
  gm.x <- reduce_dimension(data,batch,num_PCs)[[2]]
  gm.d <- dist(gm.x)
  hcl <- fastcluster::hclust(gm.d,method='ward.D')
  
  cuts <- dendextend::cutree(hcl,k=2,order_clusters_as_data=F)
  ids1 <- which(cuts==1)
  ids2 <- which(cuts==2)

  test <- test_split(data,ids1,ids2, num_PCs,
                     batch,alpha,cores,posthoc=F)

  return(test)
}


