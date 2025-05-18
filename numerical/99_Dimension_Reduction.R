##### Deviance Based Dimension Reduction ####

# Perform dimension reduction
dev_reduce_dimension <- function(y, num_PCs, family=c("poisson", "binary")[1]) {
  if (family == "binary"){
    pdev <- binary_dev_batch(y)
  }else{
    pdev <- poisson_dev_batch(y)
  }
  pdev <- scale(pdev, scale=F)
  PCs <- RSpectra::eigs_sym(as.matrix(crossprod(pdev)),k=num_PCs)
  projection <- pdev %*% PCs$vectors
  
  return(list(PCs, projection))
}

# Compute Poisson deviances
poisson_dev_batch <- function(y) {
  n <- Matrix::rowSums(y)
  pis <- Matrix::colSums(y)/sum(y)
  mu <- crossprod(array(pis,dim=c(1,length(pis))),array(n,dim=c(1,length(n))))
  mu <- t(mu)
  d <- 2 * (y * log(ifelse(y == 0, 1, y/mu)) - (y - mu))
  d[d<0] <- 0
  
  return(sqrt(d)*ifelse(y>mu,1,-1))
}

binary_dev_batch <- function(y) {
  n <- dim(y)[1]
  pis <- Matrix::colSums(y)/n
  mu <- crossprod(array(pis,dim=c(1,length(pis))),array(1,dim=c(1,n)))
  mu <- t(mu)
  d <- 2 * (y * log(ifelse(y == 0, 1, y/mu)) + 
              (1-y) * log(ifelse(y == 1, 1, (1-y)/(1-mu))))
  d[d<0] <- 0
  return(sqrt(d)*sign(y-mu))
}

#### Spectral Clustering Using Laplacian Matrix ####
kernpca <- function(X){
  S <- as.matrix(dist(X))
  D <- matrix(0, nrow=nrow(X), ncol = nrow(X))
  
  for (i in 1:nrow(X)) {
    index <- order(S[i,])[2:11]
    D[i,][index] <- 1 
  }
  D = D + t(D) 
  D[ D == 2 ] = 1
  degrees = colSums(D) 
  n = nrow(D)
  
  laplacian = ( diag(n) - diag(degrees^(-1/2)) %*% D %*% diag(degrees^(-1/2)) )
  eigenvectors = eigen(laplacian, symmetric = TRUE)
  n = nrow(laplacian)
  eigenvectors = eigenvectors$vectors[,(n - 2):(n - 1)]
  
  return(eigenvectors)
}

