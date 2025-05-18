# BiocManager::install("pkimes/sigclust2") # used for all SigClust-based methods

suppressPackageStartupMessages(library("sigclust2"))
library(filelock)
source("SigClust-MDS.R")
source("99_SigClust_GLM.R")
source("99_Dimension_Reduction.R")
source("00_DGM.R")

file_name <- paste0("./result/", opt$data, "_output.csv")

# Data Generation
generator <- generator_list[[opt$data]]
X <- generator(n=opt$n, d=opt$d, a=opt$a, seed=2024*opt$nsim, snr=opt$snr)
filter <- apply(X, 2, sd) > 0
X <- X[, filter]

# Dimension Reduction
Y <- dev_reduce_dimension(X, min(opt$n-1, 200), family=opt$data)[[2]]
dissimilarity_matrix <- dist(Y, method = "euclidean", diag = T, upper = T, p = 2)
mds_matrix <- cmdscale(dissimilarity_matrix,k = 2)
# plot(mds_matrix[, 1], mds_matrix[, 2])
# Cluster Significance
output <- SigClust_GLM(mds_matrix, nsim = 100, method="nonparametric")
P_RIFT <- output$pval
P_MRIFT <- output$pvalnorm

result <- data.frame(
  Distribution = opt$data,
  n = opt$n,
  d = opt$d,
  a = opt$a,
  Method = c("DEV-RIFT", "DEV-MRIFT"),
  pval = c(P_RIFT, P_MRIFT),
  snr = opt$snr
)

print(result)

lock_name <- paste0(file_name, ".lock")
lock <- lock(lock_name)
if (file.exists(file_name)) {
  current <- read.csv(file_name)
  result <- rbind(current, result)
}
write.csv(result, file=file_name, row.names = FALSE)

