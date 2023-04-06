
# DISCREPANCY ERSATZ ###########################################################

# S-ersatz function ------------------------------------------------------------
s_ersatz <- function(mat) {

  N <- nrow(mat)

  s <- ceiling(sqrt(N))

  # Create the zero matrix
  mat_zeroes <- matrix(0, s, s)

  # Compute index for x_i
  m <- ceiling(mat[, 1] * s)

  # Compute index for y
  x <- mat[, 2]
  n_norm <- (x-min(x))/(max(x)-min(x)) # Scale y to 0, 1
  n <- ceiling(n_norm * s)

  # Turn y==0 to y == 1
  n <- ifelse(n == 0, 1, n)

  # Merge and identify which cells are occupied by points
  ind <- cbind(m, n)
  mat_zeroes[ind] <- 1

  # Compute discrepancy
  out <- sum(mat_zeroes==1) / N

  return(out)
}

# Wrap-up ----------------------------------------------------------------------

#' Computation of the S-ersatz discrepancy.
#'
#' It allows to use the S-ersatz discrepancy measure by \insertCite{puy2023_discrepancy;textual}{sensobol}
#' as a sensitivity measure.
#'
#' @param mat A numeric matrix created with \code{\link{sobol_matrices}} and \code{matrices = "A"},
#' where each column represents an uncertain model input and each row a model simulation.
#' @param Y A numeric vector with the model output obtained from the matrix created with
#' \code{\link{sobol_matrices}}. The numeric vector should not contain any NA or NaN values.
#' @param params A character vector with the name of the model inputs.
#'
#' @importFrom Rdpack reprompt
#' @importFrom rlang :=
#'
#' @references
#' \insertAllCited{}
#'
#' @return A \code{data.table} object.
#' @export
#'
#' @details It is recommended to define \code{mat} using a power of 2 as a sample size.
#'
#' @examples
#' # Define settings
#' N <- 2^9; params <- paste("X", 1:8, sep = "")
#'
#' # Create sample matrix
#' mat <- sobol_matrices(N = N, params = params, matrices = "A")
#'
#' # Compute the Sobol' G function
#' Y <- sobol_Fun(mat)
#'
#' # Compute the S-ersatz discrepancy values
#' ind <- discrepancy_ersatz(mat = mat, Y = Y, params = params)
discrepancy_ersatz <- function(mat, Y, params) {

  parameters <- NULL

  value <- sapply(1:ncol(mat), function(j) {

    design <- cbind(mat[, j], Y)
    value <- s_ersatz(mat = design)

  })
  out <- data.table::data.table(value)[, parameters:= params]
  return(out)
}
