# R Implementation for Statistical Significance of Clustering (SigClust)

## SigClust for the Count Data

### Dependency

- R packages: `sigclust2`, `scSHC`, `optparse`, `glmpca`, `mixtools`, `mclust`, `ggplot2`, `grid`, `gridExtra`.

- System: Windows, Mac, and Linux. Note that `scSHC` uses `mcapply` for parallel computing, which is not supported by Windows.

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
