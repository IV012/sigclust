library(MASS)
source("99_Unconditional_Simple_Test.R")


#### Internal Function to Estimate the Eigenvalues ####
.vareigen<-function(x,n,p,icovest){	
  #check the dimension of x to
  #match n and p
  if(dim(x)[1]==n & dim(x)[2]==p){
    mad1<-mad(x)
    simbackvar<-mad1^2
    # Jan. 23, 07; replace eigen by
    #svd to save memory
    xcov<-cov(x)
    xeig<-eigen(xcov, symmetric=TRUE, only.values =TRUE)
    veigval<-xeig$values
    vsimeigval<-xeig$values
    
    if(icovest==1){
      taub = 0
      tauu <- .sigclustcovest(veigval,simbackvar)$tau
      etau = (tauu-taub)/100
      ids = rep(0,100)
      for(i in 1:100){
        taus = taub + (i-1)*etau
        eigval.temp <- veigval - taus
        eigval.temp[eigval.temp < simbackvar] <- simbackvar
        ids[i] <- eigval.temp[1]/sum(eigval.temp)
      }
      tau <- taub + (which.max(ids)-1)*etau
      vsimeigval <- veigval - tau
      vsimeigval[vsimeigval<simbackvar] <- simbackvar
      #vsimeigval <- .sigclustcovest(veigval,simbackvar)$veigvest
    }
    
    if(icovest==2){
      vsimeigval[veigval<0] <- 0
    }
    
    if(icovest==3){  # Use original background noise thresholded estimate
      # (from Liu, et al, JASA paper)
      vsimeigval[veigval<simbackvar] <- simbackvar
    }
    list(veigval=veigval, simbackvar=simbackvar, vsimeigval=vsimeigval)
  }else{
    print("Wrong size of matrix x!")
    return(0)
  }
}

#### Internal Function to Perform 2-means Clustering ####
.cluster<-function(x, n, p,k=2){
  clust <- kmeans(x,2,iter.max = 1000)
  withinsum <- sum(clust$withinss)
  meanp <- colMeans(x)
  tx <- t(x)
  txdiff <- tx - meanp
  totalsum <- sum(txdiff^2)
  cindex <- withinsum/totalsum
  list(clust=clust, cindex=cindex)
}

#### Internal Function to Simulate the Null Distribution (Normal) ####
.simnull<-function(vsimeigval, n, p,k){
  simnorm<-matrix(0, n, p)
  for(i in 1:n){
    simnorm[i,]<-rnorm(p, sd=sqrt(vsimeigval))
  }
  simclust<-.cluster(simnorm, n, p, k )
  list(cindex=simclust$cindex)
}

#### Internal Function to Calculate the CI ####
.CI_Calculation <- function(x,Index,k){
  p <- ncol(x)
  if(p>1){
    meanp<-colMeans(x)
    tx<-t(x)
    txdiff<-tx-meanp
    totalsum<-sum(txdiff^2)  
    Subtotalsum <- 0
    for(i in 1:k){
      SubData <- x[Index==i,]
      # take care of the cluster of one single data point
      if(length(SubData) == p){
        meanp <- SubData
      }else{
        meanp<-colMeans(SubData) 
      }
      tx <- t(SubData)
      txdiff <- tx - meanp
      Subtotalsum <- Subtotalsum + sum(txdiff^2)  
    }
  }
  if(p==1){
    meanp<-mean(x)
    txdiff<-x-meanp
    totalsum<-sum(txdiff^2)
    Subtotalsum <- 0
    for(i in 1:k){
      SubData <- x[Index==i,]
      meanp <- mean(SubData)
      txdiff <- SubData - meanp
      Subtotalsum <- Subtotalsum + sum(txdiff^2)  
    }
  }
  CI <- Subtotalsum/totalsum
  return(CI)
}


#### SigClust for Count Data ####
SigClust_GLM <- function(x, nsim=100, nrep=1, label=NULL, k=2,
                         method = c("rift", "normal", "nonparametric")[1]){
  n <- dim(x)[1]
  p <- dim(x)[2]
  
  if (n <= 4){
    stop("Too few samples to calculate the cluster significance.")
  }
  
  x <- as.matrix(x)
  if (is.null(label)){
    xclust <- .cluster(x, n, p, k)
    for(i in 1:nrep){
      clust.temp <- .cluster(x, n, p, k)
      if(clust.temp$cindex < xclust$cindex)
        xclust <- clust.temp
    }
    label <- xclust$clust$cluster
  }
  if(is.null(batch)){
    bactch <- rep(1, n)
  }
  if (method == "normal"){
    xvareigen<-.vareigen(x, n, p, icovest)
    
    # calculate the simulated cluster index
    simcindex <- rep(0,nsim)
    for(i in 1:nsim){
      xsim <- .simnull(xvareigen$vsimeigval,n,p)
      simcindex[i] <- xsim$cindex
    }
    
    cindex <- .CI_Calculation(x, label, k)
    index <- (simcindex<=cindex)
    mindex <- mean(simcindex)
    sindex <- sd(simcindex)
    pval <- sum(index)/nsim
    pvalnorm <- pnorm(cindexlab, mindex, sindex)
  }else if (method == "nonparametric"){
    Ri <- numeric(0)
    for(batch_idx in unique(batch)){
      Ri <- c(Ri, Ri.crossfit(x[batch==batch_idx,], label[batch==batch_idx]))
    }
    log.truncation.ratio = 0
    asympt = KL.Test1.mixture(Ri, log.truncation.ratio, 0.01)
    simple = KL.Test2.mixture.exact(Ri, 0.05)
    simple_median = KL.Test2.mixture.median(Ri, log.truncation.ratio, 0.05)
    pvalnorm = asympt$pvalue
    pval = simple$pvalue
    pvalmed = simple_median$pvalue
  }else if (method == "rift"){
    idx <- sample(1:n,size = n,replace = F)
    Ri <- numeric(0)
    mus = rbind(apply(x[label == 1, ], 2, mean, na.rm=TRUE),
                apply(x[label == 2, ], 2, mean, na.rm=TRUE))
    i <- 0
    fold <- idx[(floor(n/10*i)+1):floor(n/10*(i+1))]
    train.x <- x[-fold,]
    test.x <- x[fold,]
    
    model1 = densityMclust(data = train.x, G = 1, modelNames = "EEI", plot=FALSE)
    # model2 = densityMclust(data = train.x, G = 2, modelNames = "EEI", plot=FALSE)
    model2 = fitting.mixture.model(data = train.x, G = 2, modelNames = "EEI",
                                    starting.means = mus,
                                    starting.variance = diag(rep(1, ncol(x))))
    
    Ri = c(Ri, Ri.compute(test.x, model1, model2))
    print(mean(Ri))
    log.truncation.ratio = 0
    simple = KL.Test2.mixture.exact(Ri, 0.05)
    simple_median = KL.Test2.mixture.median(Ri, log.truncation.ratio, 0.05)
    pval = simple$pvalue
    pvalnorm = simple_median$pvalue
  }else{
    stop("Method Not Support!")
  }
  
  return(list(pval=pval,
              pvalnorm=pvalnorm,
              method=method,
              target=target,
              family=family))
}

Ri.crossfit <- function(x, label, nfolds=5){
  idx_1 <- sample(which(label == 1))
  idx_2 <- sample(which(label == 2))
  
  folds_1 <- split(idx_1, cut(seq_along(idx_1), nfolds, labels = FALSE))
  folds_2 <- split(idx_2, cut(seq_along(idx_2), nfolds, labels = FALSE))
  
  Ri <- numeric(0)
  for(i in 1:nfolds){
    fold <- c(folds_1[[i]], folds_2[[i]])
    train.x <- x[-fold,]
    test.x <- x[fold,]
    tr.label <- label[-fold]
    
    n1 <- sum(tr.label == 1)
    n2 <- sum(tr.label == 2)
    n.tr <- n1 + n2
    
    # under the null distribution
    mu.null <- apply(train.x, 2, mean)
    cov.null <- cov(train.x)
    
    # under the alternative
    mu.c1 <- apply(train.x[tr.label == 1, ], 2, mean, na.rm=TRUE)
    mu.c2 <- apply(train.x[tr.label == 2, ], 2, mean, na.rm=TRUE)
    cov.c1 <- cov(train.x[tr.label == 1, ])
    cov.c2 <- cov(train.x[tr.label == 2, ])
    
    # Log likelihood under the null hypothesis
    log_likelihood_null <- dmvnorm(test.x, mean = mu.null, sigma = cov.null, log = TRUE)
    
    # Log likelihood under the alternative hypothesis
    log_likelihood_alt <- log(
      (n1 / n.tr) * dmvnorm(test.x, mean = mu.c1, sigma = cov.c1) +
        (n2 / n.tr) * dmvnorm(test.x, mean = mu.c2, sigma = cov.c2)
    )
    
    Ri <- c(Ri, log_likelihood_alt - log_likelihood_null)
  }
  return(Ri)
}