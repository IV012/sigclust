# coding for MDS-based SigClust
library(MASS)

# goal:
# 1. classic MDS-based SigClust
# 2. LDA projection based SigClust

# used functions

# eigenvalues estimation based on the original data
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


.sigclustcovest <- function(vsampeigv,sig2b){
  d <- length(vsampeigv)
  #Check have some eigenvalues < sig2b
  vtaucand <- vsampeigv - sig2b
  #find threshold to preserve power
  
  which <- which(vtaucand<=0)
  icut <- which[1] - 1
  powertail <- sum(vsampeigv[(icut+1):d])
  power2shift <- sig2b*(d-icut) - powertail
  vi <- c(1:icut)
  vcumtaucand <- sort(cumsum(sort(vtaucand[vi])),decreasing=TRUE)
  vpowershifted <- (vi-1)*vtaucand[vi] + vcumtaucand
  flag <- vpowershifted < power2shift
  if(sum(flag)==0){
    itau <- 0;
  }else{
    which <- which(flag>0)
    itau <- which[1]
  }
  if(itau==1){
    powerprop <- power2shift/vpowershifted[1]
    tau <- powerprop*vtaucand[1]
  }else if(itau==0){
    powerprop <- power2shift/vpowershifted[icut] 
    tau <- powerprop*vtaucand[icut] 
  }else{
    powerprop <- (power2shift-vpowershifted[itau])/
      (vpowershifted[itau-1]-vpowershifted[itau]) 
    tau <- vtaucand[itau] + powerprop*(vtaucand[itau-1] - vtaucand[itau]) 
  }
  veigvest <- vsampeigv - rep(tau,length(vsampeigv)) 
  flag <- veigvest > sig2b 
  veigvest <- flag*veigvest + (1-flag)*(sig2b*rep(1,d))
  list(veigvest=veigvest,tau=tau)
}


# only used for k=2
# Calculate the cluster index if labels are unknown
.cluster<-function(x, n, p,k=2){
  if(n>1){
    x<-as.matrix(x)
    #check the dimension of x to match
    #n and p
    if(dim(x)[1]==n & dim(x)[2]==p){	#run 2-means on x
      # cluster index for p>=2
  
      clust<-kmeans(x,2,iter.max = 1000)
      withinsum<-sum(clust$withinss)
      meanp<-colMeans(x)
      tx<-t(x)
      txdiff<-tx-meanp
      totalsum<-sum(txdiff^2)
      cindex<-withinsum/totalsum
      list(clust=clust, cindex=cindex)
    }else{
      print("Wrong size of matrix x!")
      return(0)
    }
  }else{
    print("Only one sample left, no need for clustering!")
    return(0)
  }
}


# Calculate the cluster index based on the given class label for general
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



.simnull<-function(vsimeigval, n, p,k){
  simnorm<-matrix(0, n, p)
  for(i in 1:n){
    simnorm[i,]<-rnorm(p, sd=sqrt(vsimeigval))
  }
  simclust<-.cluster(simnorm, n, p, k )
  list(cindex=simclust$cindex)
}



# main function
SigClust_MDS <- function(x,nsim,method=1,nrep=1,labflag=0,label=0,k=2,
                         covflag=0,xvareigen=0,icovest=1,dir=NULL){
  
  # x: data input
  # nsim: number of simulations
  # method: 1 original method 2: modified method
  # nrep: number of times for k-means operations
  # label: given or not 
  # k: number of clusters
  
  ## check the dimension of x to match n and p
  
  n <- dim(x)[1]
  p <- dim(x)[2] # p is the dimension of MDS embeddings
  
  if(n>1){

    x<-as.matrix(x)
    # estimate the labels if not give
    if(labflag==0){
      ##get the cluster index based on 2-means clustering results
      xclust<-.cluster(x,n,p,k)
      for(i in 1:nrep){
        clust.temp <- .cluster(x,n,p,k)
        if(clust.temp$cindex<xclust$cindex)
          xclust <- clust.temp
        xcindex <- xclust$cindex
        label <-xclust$clust$cluster
      }
    }
        
    if(method == 1){ # original method

      # estimation of eigenvalues
      # xcov <- cov(x)
      # xeig <- eigen(xcov, symmetric=TRUE, only.values =TRUE)
      # veigval <- xeig$values
      xvareigen<-.vareigen(x,n,p,icovest)
      
      # calculate the simulated cluster index
      simcindex <- rep(0,nsim)
      for(i in 1:nsim){
        xsim <- .simnull(xvareigen$vsimeigval,n,p)
        simcindex[i] <- xsim$cindex
      }
      
      # output cluster index
      if(labflag==0){
        index <- (simcindex<=xclust$cindex)
        mindex <- mean(simcindex)
        sindex <- sd(simcindex)
        pval <- sum(index)/nsim
        pvalnorm <- pnorm(xclust$cindex,mindex,sindex)
      }
      if(labflag==1){
        #calculate CI using given label
        cindex <- .CI_Calculation(x,label,k)
        index<-(simcindex<=cindex)
        mindex<-mean(simcindex)
        sindex<-sd(simcindex)
        pval<-sum(index)/nsim
        pvalnorm<-pnorm(cindexlab,mindex,sindex)
      }
     
    }
    
    if(method == 2){
      # modified method: only applicable to the case k=2
      if(k == 2){ 
        
      # column-wise cluster index
      xcindex <- .column_wise_CI(x,label)
      
      # estimation of eigenvalues
      xcov <- cov(x)
      xeig <- eigen(xcov, symmetric=TRUE, only.values =TRUE)
      veigval <- xeig$values
      
      simcindex_list <- rep(0,nsim)
      for(i in 1:nsim){
        simnorm <- matrix(rnorm(n*p, sd=sqrt(veigval)), n, p, byrow=TRUE)
        # column-wise cluster index
        label <- kmeans(simnorm,2)$cluster
        simcindex_list[i] <- .column_wise_CI(simnorm,label)
      }
      
      index <- (simcindex_list<=xcindex)
      mindex <- mean(simcindex_list)
      sindex <- sd(simcindex_list)
      pval <- sum(index)/nsim
      pvalnorm <- pnorm(xcindex,mindex,sindex)
    }else{
        print("The modified method is only for the case where k=2")
        return(0)
      }
    }
    
    return(list(pval=pval,
                pvalnorm=pvalnorm))  
    
  }else{
    print("Only one sample left, no need for clustering!")
    return(0)
  }
  
}


# main function for SigClust-True-MDS
SigClust_True_MDS <- function(x,nsim,Sigma,method=1,nrep=1,labflag=0,label=0,k=2,
                         covflag=0,xvareigen=0,icovest=1,dir=NULL){
  
  # x: data input
  # nsim: number of simulations
  # method: 1 original method 2: modified method
  # nrep: number of times for k-means operations
  # label: given or not 
  # k: number of clusters
  
  ## check the dimension of x to match n and p
  
  n <- dim(x)[1]
  r <- dim(x)[2] # p should be a small number
  p <- length(Sigma)
  
  if(n>1){
    
    x<-as.matrix(x)
    # estimate the labels if not give
    if(labflag==0){
      ##get the cluster index based on 2-means clustering results
      xclust<-.cluster(x,n,r,k=2)
      for(i in 1:nrep){
        clust.temp <- .cluster(x,n,r,k)
        if(clust.temp$cindex<xclust$cindex)
          xclust <- clust.temp
        xcindex <- xclust$cindex
        label <-xclust$clust$cluster
      }
    }
    
    if(method == 1){ # original method
      
      ##calculate the simulated cluster index
      # data generation under the null using given Sigma
      simcindex <- rep(0,nsim)
      for(i in 1:nsim){
        simnorm <- matrix(0, n, p)
        for(k in 1:n){
          simnorm[k,]<-rnorm(p, sd=sqrt(Sigma))
        }
        dissimilarity_matrix <- dist(simnorm, method = "euclidean", diag = T, upper = T, p = 2)
        mds_matrix <- cmdscale(dissimilarity_matrix,k = r)
        xclust <- .cluster(mds_matrix, n, r)
        simcindex[i] <- xclust$cindex
      }
      
      # output cluster index
      if(labflag==0){
        index <- (simcindex<=xclust$cindex)
        mindex <- mean(simcindex)
        sindex <- sd(simcindex)
        pval <- sum(index)/nsim
        pvalnorm <- pnorm(xclust$cindex,mindex,sindex)
      }
      if(labflag==1){
        #calculate CI using given label
        cindex <- .CI_Calculation(x,label,k)
        index<-(simcindex<=cindex)
        mindex<-mean(simcindex)
        sindex<-sd(simcindex)
        pval<-sum(index)/nsim
        pvalnorm<-pnorm(cindexlab,mindex,sindex)
      }
    }
    if(method == 2){
      # modified method
      if(k == 2){
        
        # column-wise cluster index
        xcindex <- .column_wise_CI(x,label)
        
        # data generation under null using given sigma
        simcindex<-rep(0,nsim)
        for(i in 1:nsim){
          simnorm <- matrix(rnorm(n*p, sd=sqrt(Sigma)), n, p, byrow=TRUE)
          dissimilarity_matrix <- dist(simnorm, method = "euclidean", diag = T, upper = T, p = 2)
          mds_matrix <- cmdscale(dissimilarity_matrix,k = r)
          # column wise test
          label <- kmeans(mds_matrix,2)$cluster
          simcindex[i] <- .column_wise_CI(mds_matrix,label)
        }
        
        index <- (simcindex<=xcindex)
        mindex <- mean(simcindex)
        sindex <- sd(simcindex)
        pval <- sum(index)/nsim
        pvalnorm <- pnorm(xcindex,mindex,sindex)
      }else{
        print("The modified method is only for the case where k=2")
        return(0)
      }
    }
    return(list(pval=pval,
                pvalnorm=pvalnorm))  
    
  }else{
    print("Only one sample left, no need for clustering!")
    return(0)
  }
  
}

.column_wise_CI <- function(x,label){
  
  n <- nrow(x)
  r <- ncol(x) # r is the dimension of MDS embeddings
  # column-wise cluster index
  xcindex <- 1
  for(col_index in 1:r){
    xclust <- .cluster(x[,col_index], n, 1, k=2)  
    xcindex <- min(xcindex, xclust$cindex)
  }
  
  # LDA projection of MDS matrix
  LDAModel <- lda(label~., data = data.frame(x))
  ProjX = x%*%(LDAModel$scaling)  # cluster index of three datasets
  
  xclust <- .cluster(ProjX, n, 1, k =2)  
  
  # take the minimum of three cluster indexes
  xcindex <- min(xcindex,xclust$cindex)
  return(xcindex)
}



# main function for SigClust-True
sigclust_true = function(x,n_sim,Sigma){
  n = nrow(x)
  p = ncol(x)
  xclust = .cluster(x, n, p)
  xcindex<-xclust$cindex
  
  simcindex<-rep(0,n_sim)
  for(i in 1:n_sim){
    xsim<-.simnull(Sigma, n, p)
    simcindex[i]<-xsim$cindex
  }
  
  index<-(simcindex<=xcindex)
  mindex<-mean(simcindex)
  sindex<-sd(simcindex)
  pval<-sum(index)/n_sim
  pvalnorm<-pnorm(xcindex,mindex,sindex)
  return(list(pval=pval,
              pvalnorm=pvalnorm))  
}


# supplementary code for generalized SigClust methods

Sigclust_True_k <- function(x,n_sim,Index,k,icovest){
  n <- nrow(x)
  p <- ncol(x)
  xcindex <- CI_Calculation_K(x,Index,k)
  
  # eigenvalues estimation
  if(p > 1){
    eigens <- .vareigen(x,n,p,icovest)$vsimeigval
  }
  if(p==1){
    eigens <- 1
  }
  
  simcindex<-rep(0,n_sim)
  for(i in 1:n_sim){
    xsim <- .simnull_k(eigens, n, p, k)
    simcindex[i] <- xsim$cindex
  } 
  
  index <- (simcindex<=xcindex)
  mindex <- mean(simcindex)
  sindex <- sd(simcindex)
  pval <- sum(index)/n_sim
  pvalnorm <- pnorm(xcindex,mindex,sindex)
  return(list(pval=pval,
              pvalnorm=pvalnorm))  
}



# Calculate the cluster index based on the given class label for general
CI_Calculation_K <- function(x,Index,k){
  p <- ncol(x)
  if(p>1){
    meanp<-colMeans(x)
    tx<-t(x)
    txdiff<-tx-meanp
    totalsum<-sum(txdiff^2)  
    Subtotalsum <- 0
    for(i in 1:k){
      SubData <- x[Index==i,]
      meanp<-colMeans(SubData)
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


.simnull_k<-function(vsimeigval, n, p, k){
  simnorm<-matrix(0, n, p)
  for(i in 1:n){
    simnorm[i,]<-rnorm(p, sd=sqrt(vsimeigval))
  }
  simclust<-.cluster_k(simnorm, n, p, k)
  list(cindex=simclust$cindex)
}


.cluster_k<-function(x, n, p, k){
  if(n>1){
    x<-as.matrix(x)
    #check the dimension of x to match
    #n and p
    if(dim(x)[1]==n & dim(x)[2]==p){	#run 2-means on x
      clust<-kmeans(x,k,iter.max = 1000)
      withinsum<-sum(clust$withinss)
      if(p>1){
        meanp<-colMeans(x)
        tx<-t(x)
        txdiff<-tx-meanp
        totalsum<-sum(txdiff^2)  
      }
      if(p==1){
        meanp<-mean(x)
        txdiff<-x-meanp
        totalsum<-sum(txdiff^2)
      }
      cindex<-withinsum/totalsum
      list(clust=clust, cindex=cindex)
    }else{
      print("Wrong size of matrix x!")
      return(0)
    }
  }else{
    print("Only one sample left, no need for clustering!")
    return(0)
  }
}


