
# FUNCTION TO CUT BY SIZE TO BE USED IN VARS-TO
##################################################################################

CutBySize <- function(m, block.size, nb = ceiling(m / block.size)) {
  int <- m / nb
  upper <- round(1:nb * int)
  lower <- c(1, upper[-nb] + 1)
  size <- c(upper[1], diff(upper))
  cbind(lower, upper, size)
}

# COMPUTATION OF VARS-TO
##################################################################################

#' Computation of VARS Total order index (VARS-TO)
#'
#' It computes VARS-TO following \insertCite{Razavi2016a;textual}{sensobol}.
#'
#' @param Y A numeric vector with the model output obtained from the matrix created with
#' \code{\link{vars_matrices}}.
#' @param star.centers Positive integer, number of star centers.
#' @param params Character vector with the name of the model inputs.
#' @param h Distance between pairs.
#' @param method Type of computation. If \code{method = "all.step"}, all pairs of points with values
#' \eqn{\Delta h, 2\Delta h, 3\Delta h,...} are used in each dimension. If \code{method = "one.step"},
#' only the pairs \eqn{\Delta h} away are used. The default is \code{method = "all.step"}.
#'
#' @details VARS is based on variogram analysis to characterize the spatial structure and variability
#' of a given model output across the input space \insertCite{Razavi2016a}{sensobol}. Variance-
#' based total-order effects can be computed as by-products of the VARS framework. The total-order index
#' is related to the variogram \eqn{\gamma(.)} and co-variogram \eqn{C(.)} functions by the
#' following equation:
#'
#' \deqn{T_i = \frac{\gamma (h_i) + E \left [C_{\mathbf{x}_{\sim i}} (h_i) \right]}{\hat{V}(y)} }
#'
#' where \eqn{x^*_{\sim i}} is a vector of all \eqn{k} factors except \eqn{x_i}.
#'
#' @return A \code{data.table} with the VARS-TO indices of each parameter.
#'
#' @importFrom magrittr "%>%"
#' @importFrom rlang ":="
#' @export
#'
#' @references
#'
#' \insertAllCited{}
#'
#' @examples
#' # Define settings
#' star.centers <- 10; params <- paste("X", 1:3, sep = ""); h <- 0.1
#'
#' # Create STAR-VARS
#' mat <- vars_matrices(star.centers = star.centers, params = params, h = h)
#'
#' # Run model
#' y <- sensobol::ishigami_Fun(mat)
#'
#' # Compute VARS-TO
#' ind <- vars_to(Y = y, star.centers = star.centers, params = params, h = h)
vars_to <- function(Y, star.centers, params, h, method = "all.step") {
  parameters <- NULL

  # REORGANIZE THE POINTS ------------------------------------------------

  n.cross.points <- length(params) * ((1 / h) - 1) + 1
  index.centers <- seq(1, length(Y), n.cross.points)
  mat <- matrix(Y[-index.centers], ncol = star.centers)
  indices <- CutBySize(nrow(mat), nb = length(params))
  out <- list()
  for(i in 1:nrow(indices)) {
    out[[i]] <- mat[indices[i, "lower"]:indices[i, "upper"], ]
  }

  # EXTRACT PAIRS OF POINTS SEPARATED h ----------------------------------

  if(method == "one.step") {
    d <- lapply(1:length(params), function(x)
      lapply(1:ncol(out[[x]]), function(j) {
        da <- c(out[[x]][, j][1],
                rep(out[[x]][, j][-c(1, length(out[[x]][, j]))], each = 2),
                out[[x]][, j][length(out[[x]][, j])])
      }))

  # EXTRACT PAIRS OF POINTS SEPARATED h, 2h, 3h, ... -------------------

  } else if(method == "all.step") {
    d <- lapply(1:length(params), function(x)
      lapply(1:ncol(out[[x]]), function(j) {
        da <- c(utils::combn(out[[x]][, j], 2))
      }))
  } else {
    stop("method should be either one.step or all.step")
  }
  out <- lapply(d, function(x)
    lapply(x, function(y) matrix(y, nrow = length(y) / 2, byrow = TRUE)))

  # COMPUTATION OF THE VARIOGRAM -----------------------------------------

  variogr <- lapply(out, function(x) lapply(x, function(y)
    mean(0.5 * (y[, 1] - y[, 2]) ^ 2)))
  variogr <- lapply(variogr, function(x) do.call(rbind, x))
  variogr <- unlist(lapply(variogr, mean))

  # COMPUTATION OF THE COVARIOGRAM ---------------------------------------

  covariogr <- lapply(out, function(x)
    lapply(x, function(y) stats::cov(y[, 1], y[, 2])))
  covariogr <- unlist(lapply(covariogr, function(x) Rfast::colmeans(do.call(rbind, x))))

  # VARS-TO --------------------------------------------------------------

  VY <- var(Y[index.centers])
  Ti <- (variogr + covariogr) / VY
  output <- data.table::data.table(Ti)
  output[, parameters:= params]
  return(output)
}
