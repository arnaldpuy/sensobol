#' It allows to rapidly compute, bootstrap and plot up to third-order Sobol'-based sensitivity
#' indices using several state-of-the-art first and total-order estimators. Sobol' indices
#' can be computed either for models that yield a scalar as a model output or for systems of
#' differential equations. The package also provides a suit of benchmark tests functions
#' and several options to obtain publication-ready figures of the model output uncertainty
#' and sensitivity-related analysis.
#'
#' A comprehensive empirical study of several total-order estimators included in sensobol can be found in \insertCite{Puyz;textual}{sensobol}.
#'
#' @author Arnald Puy (\email{arnald.puy@pm.me})
#'
#'
#' **Maintainer**: Arnald Puy (\email{arnald.puy@pm.me})
#'
#' @references
#' \insertAllCited{}
#'
#' @name sensobol-package
#' @aliases sensobol
#' @useDynLib sensobol, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @title sensobol: Computation of Variance-Based Sensitivity Indices
#' @keywords sensitivity uncertainty modeling
"_PACKAGE"



