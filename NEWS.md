# sensobol 1.1.4

* Added the function "discrepancy_ersatz" to use the S-ersatz discrepancy as a 
sensitivity measure. The S-ersatz is published in A. Puy, P. Roy and A. Saltelli. 2023.
Discrepancy measures for sensitivity analysis. arXiv: 2206.13470).

# sensobol 1.1.3

* Added the function "load_packages" to load (and install in case they are not already in the local R setup) all required libraries in just one call.

# sensobol 1.1.2

* Fixed a bug on sobol_matrices that caused the function to fail when a single parameter was passed into the "params" argument.

* Added a line in the help page of the sobol_indices and sobol_dummy functions informing that the Y argument does not accept NA or NaN values in the model output.

# sensobol 1.1.1

* The CITATION file has changed. Citations to sensobol should now refer to the paper published in the Journal of Statistical Software: 

Arnald Puy, Samuele Lo Piano, Andrea Saltelli, and Simon A. Levin. sensobol: an R package to compute variance-based sensitivity indices. Journal of Statistical Software 102.5 (2022), pp. 1-37. doi: 10.18637/jss.v102.i05

* The metafunction has been extended with a piecewise large, piecewise small and an oscillation function.

# sensobol 1.1.0

* The package now includes support to calculate and plot up to fourth-order effects.

# sensobol 1.0.4

* Fixed redundant links in the bibliography in .Rd files.

# sensobol 1.0.3

* Added a help page for the package.

* Corrected a bug in sobol_ode.

# sensobol 1.0.2

* The function `plot_sensobol`has been eliminated. Now the Sobol' indices can be plot by a single call to `plot`.

* Third-order indices are now referred to as `Sijl` rather than `Sijk`.

* The function `sobol_convergence`has been added to check the convergence of
Sobol' indices across sub-samples of the model output.

# sensobol 1.0.1

* The output of `sobol_indices` is now an object of class `sensobol`. Besides including the indices, it also informs on the sum of first-order effects, the estimators used in the computation and the total number of model runs.

* The output of `vars_to` is now an object of class `vars`. Besides including the indices, it also informs on the number of stars and the h value used.

* The function `plot_sobol` is deprecated and will be removed from future versions. Now the output of `sobol_indices`can be printed with a call to `plot`.

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

