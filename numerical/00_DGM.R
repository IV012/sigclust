#### Simulation Studies ####
#!/usr/bin/env Rscript
library("optparse")

option_list = list(
  make_option(c("-s", "--nsim"), type="numeric", default=100,
              help="number of simulations", metavar="numeric"),
  make_option(c("-d", "--data"), type="character", default="poisson",
              help="generative data distribution (binary, poisson, poislog, multinom)",
              metavar="character"),
  make_option(c("-n", "--n"), type="numeric", default=1000,
              help="number of samples", metavar="numeric"),
  make_option(c("-p", "--d"), type="numeric", default=1000,
              help="number of variables", metavar="numeric"),
  make_option(c("-a", "--a"), type="numeric", default=0,
              help="difference of clusters", metavar="numeric"),
  make_option(c("-r", "--snr"), type="numeric", default=0.1,
              help="signal noise ratio", metavar = "numeric")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#### Distributions Used in Simulation Studies ####

generate.binary <- function(n=1000, d=1000, a=FALSE, sigma=2, seed=1, snr=0.1){
  set.seed(seed)
  mu <- abs(runif(d, 0, 1))
  mu1 <- mu2 <- mu; 
  if(a){
    diff <- sample(d, d*snr)
    mu1[diff] <- abs(runif(diff, 0, 1))
  }
  
  set.seed(seed+1)
  X <- lapply(mu1, function(x) rbinom(n/2, 1, x))
  X <- matrix(unlist(X), n/2, d)
  set.seed(seed+2)
  Y <- lapply(mu2, function(x) rbinom(n/2, 1, x))
  Y <- matrix(unlist(Y), n/2, d)
  
  return(rbind(X, Y))
}

generate.poisson <- function(n=1000, d=1000, a=0, seed=1, snr=0.1, spike=FALSE){
  diff <- sample(d, d*snr)
  set.seed(seed)
  mu <- exp(rnorm(d, mean = 0, sd=1))
  mu1 <- mu2 <- mu; mu1[diff] <- mu1[diff]*exp(rnorm(length(diff), mean=0, sd=sqrt(a)))
  
  set.seed(seed+1)
  X <- lapply(mu1, function(x) rpois(n/2, x))
  X <- matrix(unlist(X), n/2, d)
  set.seed(seed+2)
  Y <- lapply(mu2, function(x) rpois(n/2, x))
  Y <- matrix(unlist(Y), n/2, d)
  
  return(rbind(X, Y))
}

generate.poislog <- function(n=1000, d=1000, a=0, seed=1, snr=0.1, spike=FALSE){
  # d must be divided by 10, n must be divided by 2
  diff <- sample(d, d*snr)
  set.seed(seed)
  mu <- matrix(exp(rnorm(d*n)),nrow=d)
  mu[diff ,1:(n/2)] <- mu[diff ,1:(n/2)] * exp(rnorm(d/10, sd=a))
  Y <- t(matrix(rpois(prod(dim(mu)),mu), nrow=nrow(mu)))
  return(Y)
}

generate.multinom <- function(n=5000, d=1000, a=1, seed=1, nclust=2, snr=0.1){
  set.seed(seed)
  ncells <- n / 2
  ngenes <- d
  ngenes_informative<-ngenes*snr
  
  # simulate two batches with different depths
  batch<-rep(1:2, each = nclust*ncells/2)
  ncounts <- rpois(ncells*nclust, lambda = 1000*batch)
  
  # generate profiles for clusters
  profiles_informative <- replicate(nclust, exp(rnorm(ngenes_informative, sd=a)))
  profiles_const<-matrix(ncol=nclust,rep(exp(rnorm(ngenes-ngenes_informative, sd=2)),nclust))
  profiles <- rbind(profiles_informative,profiles_const)
  
  # generate cluster labels
  clust <- sample(rep(1:nclust, each = ncells))
  
  # generate single-cell transcriptomes 
  counts <- sapply(seq_along(clust), function(i){
    rmultinom(1, ncounts[i], prob = profiles[,clust[i]])
  })
  rownames(counts) <- paste("gene", seq(nrow(counts)), sep = "_")
  colnames(counts) <- paste("cell", seq(ncol(counts)), sep = "_")
  
  # clean up rows
  Y <- counts[rowSums(counts) > 0, ]
  return(t(Y))
}

generator_list <- list(generate.binary, generate.poisson, generate.poislog, generate.multinom)
names(generator_list) <- c("binary", "poisson", "poislog", "multinom")