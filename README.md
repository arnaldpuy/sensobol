 [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/arnaldpuy/sensobol?branch=master&svg=true)](https://ci.appveyor.com/project/arnaldpuy/sensobol) [![Travis build status](https://travis-ci.org/arnaldpuy/sensobol.svg?branch=master)](https://travis-ci.org/arnaldpuy/sensobol) [![Codecov test coverage](https://codecov.io/gh/arnaldpuy/sensobol/branch/master/graph/badge.svg)](https://codecov.io/gh/arnaldpuy/sensobol?branch=master)
 
 
# sensobol

The goal of `sensobol` is to provide a set of functions to swiftly compute and visualize up to third-order Sobol' sensitivity indices. The functions allow to: 
- Create the sample matrices for the model evaluation.
- Compute and bootstrap up to third-order effects.
- Assess the approximation error of Sobol' indices.
- Plot the model uncertainty and the Sobol' indices.

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
N <- 1000
params <- paste("X", 1:3, sep = "")

## Create sample matrix to compute first and total-order indices:
mat <- sobol_matrices(N = N, params = params)

## Compute the model output (using the Ishigami test function):
Y <- ishigami_Mapply(mat)

## Compute and bootstrap the Sobol' indices:
sens <- sobol_indices(Y = Y, N = N, params = params, boot = TRUE, R = 100)
```

## Citation

Please use the following citation if you use `sensobol` in your publications:

```r
Arnald Puy (2019). sensobol: Computation of High-Order Sobol' Sensitivity Indices. R package
  version 0.2.2 http://github.com/arnaldpuy/sensobol
```

A BibTex entry for LaTex users is:

```r
@Manual{,
    title = {sensobol: Computation of High-Order Sobol' Sensitivity Indices},
    author = {Arnald Puy},
    year = {2019},
    note = {R package version 0.2.0},
    url = {http://github.com/arnaldpuy/sensobol},
  }
```
