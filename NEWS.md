# sensobol 1.0.0

* This is a major package upgrade.

* The package now includes four first-order (Sobol', Salteli, Jansen and Azzini) and eight total-order (Jansen, Sobol', Homma, Saltelli, Janon, Glen, Azzini and VARS-TO) sensitivity estimators.

* The sample matrix can be constructed either with Sobol' quasi-random numbers, a latin hypercube design or random numbers.

* Several functions to plot the results of the uncertainty and sensitivity analysis have been incorporated: the function `plot_scatter` plots the model inputs against the output, whereas the function `plot_multiscatter` plots x_i against x_j and maps the resulting coordinate to its respective model output value.

* The package can now be used in models with either a scalar or a multivariate output.

# sensobol 0.2.2

* This release prepares the users for a major upcoming improvement of the package.

* Added warnings to the function `sobol_matrices`. The arguments `n` and `k` will
be substituted by `N` and `params` in the next release of the package. The 
arguments `second` and `third` will be substituted by the argument `order` in the next
release of the package.

* Added warnings to the functions `plot_uncertainty`, `plot_scatter`, `sobol_dummy` and `sobol_indices`. The argument `n` will be substituted by `N` in the next release of the package. 

* The function `sobol_ci` is deprecated. The computation of confidence intervals
will be done directly by `sobol_indices` in the next version of the package.

* The function `sobol_ci_dummy` is deprecated. The computation of confidence intervals
for dummy parameters will be done directly by `sobol_dummy` in the next version of the package.

* The function `sobol_replicas` is deprecated and will be removed in the next version
of the package.

* The test function `ishigami_Mapply` would be renamed as `ishigami_Fun` in the next
release of the package.

# sensobol 0.2.1

* Corrected a bug in the `plot_uncertainty` function. Now the function
demands to add the initial sample size `n` to visualize the model output uncertainty.

* Corrected some functions to adapt them to data.table 1.12.4.

* Added a new function (`sobol_replicas`) to easily extract the bootstrapped Sobol' indices.

# sensobol 0.2.0

* New option in the `sobol_matrices` function: the option `cluster` allows to create Sobol' matrices for clusters of parameters.

* New test functions added: 
  - Bratley & Fox 1988.
  - Bratley et al. 1992
  - Oakley & O'Hagan 2004
  
* Added references to all test functions.

# sensobol 0.1.1

* The vignette rendered wrongly in the previous version;
now the issue is corrected.

* Corrected the following note, found in the CRAN check: 
  - Namespace in Imports field not imported from: ‘parallel’.
 
* Some functions that were exported to the R package manual 
are now internal.

# sensobol 0.1.0

* Added a NEWS.md file to track changes to the package.

