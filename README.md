[![](https://cranlogs.r-pkg.org/badges/sensobol)](https://cran.r-project.org/package=sensobol)
[![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/grand-total/sensobol?color=blue)](https://r-pkg.org/pkg/sensobol)
[![](https://www.r-pkg.org/badges/version/sensobol?color=orange)](https://cran.r-project.org/package=sensobol)
 
 
# sensobol: an R package to compute variance-based sensitivity indices

The ``R`` package ``sensobol`` provides several functions to conduct variance-based uncertainty and sensitivity analysis, from the estimation of sensitivity indices to the visual representation of the results. It implements several state-of-the-art first and total-order estimators and allows the computation of up to third-order effects, as well as of the approximation error, in a swift and user-friendly way.

## Installation
To install the stable version on [CRAN](https://CRAN.R-project.org/package=sensobol), use

```r
install.packages("sensobol")
```
To install the development version, use devtools:

``` r
install.packages("devtools") # if you have not installed devtools package already
devtools::install_github("arnaldpuy/sensobol", build_vignettes = TRUE)
```

## Example

This brief example shows how to compute Sobol' indices. For a more detailed explanation of the package functions, check the vignette.

``` r
## Load the package:
library(sensobol)

## Define the base sample size and the parameters
N <- 2 ^ 8
params <- paste("X", 1:3, sep = "")

## Create sample matrix to compute first and total-order indices:
mat <- sobol_matrices(N = N, params = params)

## Compute the model output (using the Ishigami test function):
Y <- ishigami_Fun(mat)

## Compute and bootstrap the Sobol' indices:
ind <- sobol_indices(Y = Y, N = N, params = params)
```

## Citation

Please use the following citation if you use `sensobol` in your publications:

```r
A. Puy, S. Lo Piano, A. Saltelli, S. A. Levin (2021). sensobol: Computation of
  Variance-Based Sensitivity Indices. arxiv:2101.10103.
```

A BibTex entry for LaTex users is:

```r
@Manual{,
    title = {{sensobol}: {C}omputation of Variance-Based Sensitivity Indices},
    author = {Arnald Puy and Samuele Lo Piano and Andrea Satelli and Simon A. Levin},
    journal = {arxiv:2101.10103},
    year = {2021},
    url = {https://github.com/arnaldpuy/sensobol},
  }
```
