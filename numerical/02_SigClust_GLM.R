suppressPackageStartupMessages(library("sigclust2"))
library(filelock)
source("SigClust-MDS.R")
source("99_SigClust_GLM.R")
source("00_DGM.R")
library(glmpca)
library(parallel)

P_RIFT <- rep(NA, opt$nsim)
P_MRIFT <- rep(NA, opt$nsim)
file_name <- paste0("./result/", opt$data, "_output.csv")

sim.glm <- function(i){
  # Data Generation
  generator <- generator_list[[opt$data]]
  X <- generator(opt$n, opt$d, a=opt$a, seed=2024*i)
  filter <- apply(X, 2, sd) > 0
  X <- X[, filter]
  
  # Dimension Reduction
  if (opt$data == "binary"){
    Y <- glmpca(Y=t(X), L=min(opt$n/5, 30), fam="binom", sz=1, optimizer = "fisher")$factors
  }else if (opt$data == "multinom"){
    Y <- glmpca(Y=t(X), L=min(opt$n/5, 30), fam="binom", sz=NULL, optimizer = "fisher")$factors
  }else{
    Y <- glmpca(Y=t(X), L=min(opt$n/5, 30), fam="poi", optimizer = "fisher")$factors
  }

  dissimilarity_matrix <- dist(Y, method = "euclidean", diag = T, upper = T, p = 2)
  mds_matrix <- scale(cmdscale(dissimilarity_matrix,k = 2))
  # plot(mds_matrix)
  clust <- hclust(dissimilarity_matrix)
  labels <- cutree(clust, k = 2)
  if (sum(labels==2) < 0.1 * dim(X)[1]) {
    labels <- sample(1:2, dim(X)[1], replace=TRUE)
  }
  
  # Cluster Significance
  output <- SigClust_GLM(mds_matrix, nsim = 100, method="rift", label=labels)
  return(c(output$pval, output$pvalnorm))
}
output <- parallel::mclapply(seq(opt$nsim), sim.glm, mc.cores = 4)

i <- 1
for (pvals in output){
  P_RIFT[i] <- pvals[1]
  P_MRIFT[i] <- pvals[2]
  i <- i + 1
}

result <- data.frame(
  Distribution = opt$data,
  n = opt$n,
  d = opt$d,
  a = opt$a,
  Method = c(rep("GLM-RIFT", opt$nsim), rep("GLM-MRIFT", opt$nsim)),
  pval = c(P_RIFT, P_MRIFT)
)

lock_name <- paste0(file_name, ".lock")
lock <- lock(lock_name)
if (file.exists(file_name)) {
  current <- read.csv(file_name)
  result <- rbind(current, result)
}
write.csv(result, file=file_name, row.names = FALSE)