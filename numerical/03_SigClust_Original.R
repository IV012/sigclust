# BiocManager::install("pkimes/sigclust2") # used for all SigClust-based methods

suppressPackageStartupMessages(library("sigclust2"))
library(filelock)
source("SigClust-MDS.R")
source("99_Dimension_Reduction.R")
source("00_DGM.R")

file_name <- paste0("./result/", opt$data, "_output.csv")

generator <- generator_list[[opt$data]]
X <- generator(n=opt$n, d=opt$d, a=opt$a, seed=2024*opt$nsim, snr=opt$snr)
filter <- apply(X, 2, sd) > 0
X <- X[, filter]

# Dimension Reduction
# dissimilarity_matrix <- dist(X, method = "euclidean", diag = T, upper = T, p = 2)
# clust <- hclust(dissimilarity_matrix)
label <- kmeans(X, centers = 2, iter.max=1000, algorithm = "Forgy")$cluster

try(P_Soft <- sigclust(x=X, icovest=1, labels=label)$p_norm)
try(P_Hard <- sigclust(x=X, icovest=3, labels=label)$p_norm)

# mds_matrix <- cmdscale(dissimilarity_matrix,k = 2)
try(mds_matrix <- prcomp(x=X)$x[, 1:2])
try(output <- SigClust_MDS(mds_matrix, nsim = 100, label=label, method=2))
try(P_MDS <- output$pvalnorm)


result <- data.frame(
  Distribution = opt$data,
  n = opt$n,
  d = opt$d,
  a = opt$a,
  Method = c("SigClust-Hard", "SigClust-Soft", "SigClust-MDS"),
  pval = c(P_Hard, P_Soft, P_MDS),
  snr = opt$snr
)

lock_name <- paste0(file_name, ".lock")
lock <- lock(lock_name)
if (file.exists(file_name)) {
  current <- read.csv(file_name)
  result <- rbind(current, result)
}
write.csv(result, file=file_name, row.names = FALSE)

