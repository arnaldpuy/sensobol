[![](https://cranlogs.r-pkg.org/badges/sensobol)](https://cran.r-project.org/package=sensobol)
[![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/grand-total/sensobol?color=blue)](https://r-pkg.org/pkg/sensobol)
[![](https://www.r-pkg.org/badges/version/sensobol?color=orange)](https://cran.r-project.org/package=sensobol)
 
 
# sensobol: an R package to compute variance-based sensitivity indices

The `R` package `sensobol` provides functions to conduct variance-based
uncertainty and sensitivity analysis, from the estimation of the sensitivity
indices to the visual representation of the results. It implements several
state-of-the-art first and total-order estimators, computes effects up to the
fourth order, supports the analysis of grouped (e.g. correlated) inputs, and
offers randomised quasi-Monte Carlo designs — all in a swift and user-friendly
way.

## Installation

To install the stable version on [CRAN](https://CRAN.R-project.org/package=sensobol), use

```r
install.packages("sensobol")
```

To install the development version, use devtools:

``` r
install.packages("devtools") # if you have not installed devtools already
devtools::install_github("arnaldpuy/sensobol", build_vignettes = TRUE)
```

## Features

### Sensitivity indices

* **First and total-order indices** for every model input.
* **Higher-order effects:** second-, third- and fourth-order interaction
  indices, selected through `order = "first" / "second" / "third" / "fourth"`
  in `sobol_matrices()` and `sobol_indices()`.
* **Grouped inputs:** compute indices for *groups* of parameters rather than
  for individual parameters via the `groups` argument. This is the standard
  device for handling correlated inputs by moving them together.
* **Bootstrap confidence intervals** for all indices (`boot = TRUE`), with a
  choice of interval method (`norm`, `basic`, `percent`, `bca`).
* **Dummy parameter** index (`sobol_dummy()`) to gauge the numerical
  approximation error / noise floor of the estimates.
* **VARS-TO** total-order estimates through the star-based VARS method
  (`vars_matrices()`, `vars_to()`).
* **S-ersatz discrepancy** measure as an alternative, given-data sensitivity
  measure (`discrepancy_ersatz()`).
* Support for models defined as **systems of differential equations**
  (`sobol_ode()`).

### Estimators

First-order (`first =`) and total-order (`total =`) estimators available in
`sobol_indices()`:

| First-order (`first`) | Reference | Total-order (`total`) | Reference |
|---|---|---|---|
| `"sobol"` | Sobol' (1993) | `"jansen"` | Jansen (1999) |
| `"saltelli"` | Saltelli et al. (2010) | `"sobol"` | Sobol' (2001) |
| `"jansen"` | Jansen (1999) | `"homma"` | Homma & Saltelli (1996) |
| `"azzini"` | Azzini et al. (2020) | `"janon"` | Janon et al. (2014) |
| `"owen"` | Owen (2013) | `"glen"` | Glen & Isaacs (2012) |
| `"martinez"` | Martinez (2011) | `"azzini"` | Azzini et al. (2020) |
| `"mauntz"` | Saltelli et al. (2010) | `"saltelli"` | Saltelli et al. (2008) |
| | | `"owen"` | Owen (2013) |

Each estimator requires a specific sampling design (`matrices` argument of
`sobol_matrices()`). `sobol_indices()` checks the correspondence and raises an
informative error if an estimator is paired with an incompatible design. See
Table 3 of the vignette for the full estimator/design correspondence.

The `"owen"` first-order estimator is bias-corrected at `S_i = 0`, which makes
it well suited to discriminating genuinely non-influential inputs from
inputs with small but non-zero effects.

### Sampling

* **Sample matrices** with `sobol_matrices()` using either Sobol' quasi-random
  numbers (`type = "QRN"`, default), a Latin Hypercube Sampling design
  (`type = "LHS"`) or plain random numbers (`type = "R"`).
* **Randomised quasi-Monte Carlo (RQMC)** through the `scrambling` argument:
  - `"none"` (default): the deterministic Sobol' sequence, as in earlier
    releases.
  - `"shift"`: a Cranley-Patterson digital shift (no extra dependency).
  - `"owen"`: an in-house Sobol' generator (Joe-Kuo direction numbers) with
    hash-based Owen scrambling, independent of `randtoolbox` and supporting
    up to 250 dimensions.

  Use `seed` to make a scrambled design reproducible. Randomisation yields
  unbiased estimators and allows confidence intervals to be obtained from
  independent replications.

### Test functions and visualisation

* **Benchmark functions:** `ishigami_Fun()`, `sobol_Fun()`,
  `bratley1988_Fun()`, `bratley1992_Fun()`, `oakley_Fun()`, and a random
  `metafunction()` generator for large-scale benchmarking.
* **Plots:** `plot()` for the sensitivity indices of a `sensobol` object,
  `plot_uncertainty()` for the model-output distribution, `plot_scatter()` for
  model input–output scatterplots, and `plot_multiscatter()` for pairwise
  input combinations coloured by the output.

## Example

This brief example shows how to compute Sobol' indices. For a more detailed
explanation of the package functions, check the vignette.

``` r
## Load the package:
library(sensobol)

## Define the base sample size and the parameters:
N <- 2 ^ 10
params <- paste("X", 1:3, sep = "")

## Create the sample matrix to compute first and total-order indices:
mat <- sobol_matrices(N = N, params = params)

## Compute the model output (using the Ishigami test function):
Y <- ishigami_Fun(mat)

## Compute and bootstrap the Sobol' indices:
ind <- sobol_indices(Y = Y, N = N, params = params, boot = TRUE, R = 100)
ind

## Plot the indices:
plot(ind)
```

A few common variations:

``` r
## Up to second-order effects:
mat <- sobol_matrices(N = N, params = params, order = "second")
Y   <- ishigami_Fun(mat)
ind <- sobol_indices(Y = Y, N = N, params = params, order = "second")

## Grouped inputs (treat X2 and X3 as a single correlated group):
groups <- list(g1 = "X1", g2 = c("X2", "X3"))
mat <- sobol_matrices(N = N, params = params, groups = groups)
Y   <- ishigami_Fun(mat)
ind <- sobol_indices(Y = Y, N = N, params = params, groups = groups)

## A different estimator pairing on the four-matrix design:
mat <- sobol_matrices(N = N, params = params, matrices = c("A", "B", "AB", "BA"))
Y   <- ishigami_Fun(mat)
ind <- sobol_indices(Y = Y, N = N, params = params,
                     matrices = c("A", "B", "AB", "BA"),
                     first = "owen", total = "owen")

## Randomised QMC with reproducible Owen scrambling:
mat <- sobol_matrices(N = N, params = params, scrambling = "owen", seed = 1)
```

## Citation

Please use the following citation if you use `sensobol` in your publications:

```r
A. Puy, S. Lo Piano, A. Saltelli, S. A. Levin (2022). sensobol: Computation of 
Variance-Based Sensitivity Indices. Journal of Statistical Software 102(5), 
1-37. doi:10.18637/jss.v102.i05.
```

A BibTex entry for LaTex users is:

```r
@article{,
author = {Puy, Arnald and {Lo Piano}, Samuele and Saltelli, Andrea and Levin, Simon A.},
journal = {Journal of Statistical Software},
title = {{sensobol: an R package to compute variance-based sensitivity indices}},
doi = {10.18637/jss.v102.i05},
volume = {102}, 
number = {5},
pages = {1--37},
year = {2022}
}
```
