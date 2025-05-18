# R Implementation for Statistical Significance of Clustering (SigClust)

## SigClust for the Count Data (SigClust-DEV)

### Dependency

- R packages: `sigclust2`, `scSHC`, `optparse`, `glmpca`, `mixtools`, `mclust`, `ggplot2`, `grid`, `gridExtra`.

- System: Windows, Mac, and Linux. Note that `scSHC` uses `mcapply` for parallel computing, which is not supported by Windows.

### General Implementation

The implementation of SigClust-DEV can be found in the folder `./numerical`.
The main files are `99_SigClust_GLM.R` and `99_Dimension_Reduction.R`, while we also include other files for reproducing the simulation results.

The function for deviance-based dimension reduction is `dev_reduce_dimension` within `99_Dimension_Reduction.R`, with the following arguments:

- `y`: the input data of an $\(n \times p\)$ matrix.

- `num_PCs`: the number of $k$ principle components to be used.

- `family`: The exponential family to be fitted. Either `poisson` or `binary`. Default `poisson`.

The output of `dev_reduce_dimension` includes `PCs` of an $p$ by $k$ loading matrix, and `projection` of a $n$ by $k$ projected data.

The function for cluster‑significance testing is `SigClust_GLM` within `SigClust_GLM.R`, with the following arguments:

- `x`: the input data of an \(n \times p\) matrix.  
- `nsim`: the number of Monte‑Carlo simulations used to approximate the null distribution. Default `100`.  
- `nrep`: the number of repeated k‑means initialisations when determining the observed cluster index. Default `1`.  
- `label`: an optional vector of length \(n\) containing existing two‑cluster assignments; if `NULL`, clusters are estimated internally. Default `NULL`.  
- `k`: the number of clusters; currently only `2` is supported. Default `2`.  
- `batch`: an optional vector of batch labels of length \(n\) for stratified tests; if `NULL`, all observations are treated as one batch. Default `NULL`.  
- `method`: the null‑distribution type, chosen from `"rift"` or `"normal"`. Default `"rift"`.

The output of  `SigClust_GLM` includes `pval`, `pvalnorm`, and `method`. For the best power, we encourage users to use `pvalnorm` by default.

Users may follow the sample code to test the clustering significance for a count matrix.

```r
# load the files for SigClust-DEV implementation
source("99_SigClust_GLM.R")
source("99_Dimension_Reduction.R")

# generate the simulation dataset
source("00_DGM.R")
X <- generator(n=opt$n, d=opt$d, a=opt$a, seed=2024*i)

# Dimension Reduction
Y <- dev_reduce_dimension(X, k=2, family=opt$data)[[2]]
# Cluster Significance
SigClust_GLM(mds_matrix, method="rift")$pvalnorm
```

### Numerical Studies

R scripts in the folder `./numerical` starting with 01 to 04 were used for simulation (i.e., Figure 1-3). Taking `01_SigClust_DEV.R` as an example, you may use `-h` and `--help` to know how to set the simulation settings.

```bash
Rscript ./numerical/01_SigClust_DEV.R --help
Rscript ./numerical/01_SigClust_DEV.R --n=1000 --data="binary" --a=0 # 100 Simulations for 1000*1000 Binary Data will be performed
```

### Real Data Applications

**Hydra scRNA data** can be found at the [data portal](https://portals.broadinstitute.org/single_cell/study/stem-cell-differentiation-trajectories-in-hydra-resolved-at-single-cell-resolution). We follow the [guideline](https://github.com/cejuliano/hydra_single_cell) here to process the data.

**The EHR datasets** are confidential due to privacy issues. Hierachical clustering with SigClust-DEV can be found in `./applications/SHC-GLM.R`.

## SigClust for General Non-Gaussian Data

To be implemented.
