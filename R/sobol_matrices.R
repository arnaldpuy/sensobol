
# CREATION OF THE AB, BA, CB matrices
##################################################################################

scrambled_sobol <- function(matrices, A, B, C, order) {
  first <- 1:ncol(A)
  N <- nrow(A)

  # VECTORS WITH THE COLUMNS
  # -----------------------------------------------------------------

  if(order == "first") {
    loop <- first

  } else if (order == "second") {
    second <- c(first, utils::combn(1:ncol(A), 2, simplify = FALSE))
    loop <- second

  } else if (order == "third") {
    second <- c(first, utils::combn(1:ncol(A), 2, simplify = FALSE))
    third <- c(second, utils::combn(1:ncol(A), 3, simplify = FALSE))
    loop <- third

  } else {
    stop("order should be either first, second or third")
  }

  # CHECK WHICH MATRICES HAVE BEEN SELECTED
  # -----------------------------------------------------------------

  AB.mat <- "AB" %in% matrices
  BA.mat <- "BA" %in% matrices
  CB.mat <- "CB" %in% matrices

  # CONSTRUCT AB, BA MATRICES, ETC
  # -----------------------------------------------------------------

  if (AB.mat == TRUE) {
    X <- rbind(A, B)

    for(i in loop) {
      AB <- A
      AB[, i] <- B[, i]
      X <- rbind(X, AB)
    }
    AB <- X[(2 * N + 1):nrow(X), ]

  } else if (AB.mat == FALSE) {
    AB <- NULL
  }

  if (BA.mat == TRUE) {
    W <- rbind(A, B)
    for (i in loop) {
      BA <- B
      BA[, i] <- A[, i]
      W <- rbind(W, BA)
    }
    BA <- W[(2 * N + 1) : nrow(W), ]

  } else if (BA.mat == FALSE) {
    BA <- NULL
  }

  if (CB.mat == TRUE) {
    Z <- rbind(A, B)
    for(i in loop) {
      CB <- C
      CB[, i] <- B[, i]
      Z <- rbind(Z, CB)
    }
    CB <- Z[(2 * N + 1) : nrow(Z), ]

  } else if (CB.mat == FALSE) {
    CB <- NULL
  }

  # MERGE AND OUTPUT
  # -----------------------------------------------------------------

  final <- rbind(AB, BA, CB)
  return(final)
}

# FUNCTION TO CREATE THE SAMPLE MATRICES
##################################################################################

#' Creation of the sample matrices
#'
#' It creates the sample matrices to compute Sobol' first and total-order indices.
#' If needed, it also creates the sample matrices required to compute second and
#' third-order indices.
#'
#' @param matrices Character vector with the required matrices. The default
#' is \code{matrices = c("A", "B", "AB")}.
#' @param N Positive integer, initial sample size of the base sample matrix.
#' @param params Character vector with the name of the model inputs.
#' @param order One of "first", "second" or "third" to create a matrix to
#' compute first, second or up to third-order Sobol indices. The default is
#' \code{order = "first"}.
#' @param type Approach to construct the sample matrix. Options are:
#' * \code{type = "QRN"} (default): It uses \insertCite{Sobol1967;textual}{sensobol} Quasi-Random Numbers.
#' through a call to the function \code{\link{sobol}} of the \code{randtoolbox} package.
#' * \code{type = "LHS"}: It uses a Latin Hypercube Sampling Design
#' \insertCite{McKay1979}{sensobol} through a call
#' to the function \code{\link{randomLHS}} of the \code{lhs} package.
#' * \code{type = "R"}: It uses random numbers.
#' @param ... Further arguments in \code{\link{sobol}}.
#' @return A numeric matrix where each column is a model input distributed in (0,1) and each row
#' a sampling point.
#' @export

#' @details Before calling \code{sobol_matrices}, the user must decide which estimators
#' will be used to compute first and total-order indices, for this option conditions
#' the design of the sample matrix and therefore the argument \code{matrices}.
#' See Table 3 in the vignette for further details on the specific sampling designs required by
#' the estimators.
#'
#' The user can select one of the following sampling designs:
#' * \eqn{\mathbf{A}}, \eqn{\mathbf{B}}, \eqn{\mathbf{A}_B^{(i)}}.
#' * \eqn{\mathbf{A}}, \eqn{\mathbf{B}}, \eqn{\mathbf{B}_A^{(i)}}.
#' * \eqn{\mathbf{A}}, \eqn{\mathbf{B}}, \eqn{\mathbf{A}_B^{(i)}}, \eqn{\mathbf{B}_A^{(i)}}.
#'
#' If \code{order = "first"}, the function creates an \eqn{(N, 2k)} matrix according to the approach defined by
#' \code{type}, where the leftmost and the rightmost \eqn{k} columns are respectively allocated
#' to the \eqn{\mathbf{A}} and the \eqn{\mathbf{B}} matrix. Depending on the sampling design, it
#' also creates \eqn{k} \eqn{\mathbf{A}_B^{(i)}} (\eqn{\mathbf{B}_A^{(i)}}) matrices, where all
#' columns come from \eqn{\mathbf{A}} (\eqn{\mathbf{B}}) except the \eqn{i}-th, which comes from
#' \eqn{\mathbf{B}} (\eqn{\mathbf{A}}). All matrices are returned row-binded.
#'
#' If \code{order = "second"}, \eqn{\frac{k!}{2!(k-2)!}} extra \eqn{(N, k)} \eqn{\mathbf{A}_B^{(ij)}}
#' (\eqn{\mathbf{B}_A^{(ij)}}) matrices are created, where all columns come from \eqn{\mathbf{A}}
#' (\eqn{\mathbf{B}}) except the \eqn{i}-th and \eqn{j}-th, which come from \eqn{\mathbf{B}}
#' (\eqn{\mathbf{A}}). These matrices allow the computation of second-order effects, and are row-bound
#' to those created for first and total-order indices.
#'
#' If \code{order = "third"}, \eqn{\frac{k!}{3!(k-3)!}} extra \eqn{(N, k)} \eqn{\mathbf{A}_B^{(ijl)}}
#' (\eqn{\mathbf{B}_A^{(ijl)}}) matrices are bound below those created
#' for the computation of second-order effects. In these matrices, all columns come from \eqn{\mathbf{A}}
#' (\eqn{\mathbf{B}}) except the \eqn{i}-th, the \eqn{j}-th and the \eqn{l}-th, which come from \eqn{\mathbf{B}}
#' (\eqn{\mathbf{A}}). These matrices are needed to compute third-order effects, and are row-bound below
#' those created for second-order effects.
#'
#' All columns are distributed in (0,1). If the uncertainty in some parameter(s) is better described with
#' another distribution, the user should apply the required quantile inverse transformation to the column of
#' interest once the sample matrix is produced.
#'
#' @importFrom Rdpack reprompt
#'
#' @references
#'
#' \insertAllCited{}
#'
#' @examples
#' # Define settings
#' N <- 100; params <- paste("X", 1:10, sep = ""); order <-  "third"
#'
#' # Create sample matrix using Sobol' Quasi Random Numbers.
#' mat <- sobol_matrices(N = N, params = params, order = order)
#'
#' # Let's assume that the uncertainty in X3 is better described
#' # with a normal distribution with mean 0 and standard deviation 1:
#' mat[, 3] <- qnorm(mat[, 3], 0, 1)
sobol_matrices <- function(matrices = c("A", "B", "AB"),
                           N, params, order = "first",
                           type = "QRN", ...) {
  k <- length(params)
  n.matrices <- ifelse(any(stringr::str_detect(matrices, "C")) == FALSE, 2, 3)

  # SELECTION OF THE TYPE OF SAMPLE MATRIX
  # -----------------------------------------------------------------

  if (type == "QRN") {
    df <- randtoolbox::sobol(n = N, dim = k * n.matrices, ...)

  } else if (type == "R") {
    df <- replicate(k * n.matrices, stats::runif(N))

  } else if (type == "LHS") {
    df <- lhs::randomLHS(N, n.matrices * k)

  } else {
    stop ("method should be either QRN, R or LHS")
  }

  # CONSTRUCTION OF A, B AND C MATICES
  # -----------------------------------------------------------------

  A <- df[, 1:k]
  B <- df[, (k + 1) : (k * 2)]

  if (n.matrices == 3) {
    C <- df[, ((k * 2) + 1):(k * 3)]

  } else {
    C <- NULL
  }

  # CONSTRUCTION OF AB, BA MATRICES, ETC
  # -----------------------------------------------------------------

  out <- scrambled_sobol(matrices = matrices,
                         A = A, B = B, C = C,
                         order = order)
  A.mat <- "A" %in% matrices
  B.mat <- "B" %in% matrices
  C.mat <- "C" %in% matrices

  # SET NULL TO THOSE NOT USED
  # -----------------------------------------------------------------

  if (A.mat == FALSE) {
    A <- NULL
  }

  if (B.mat == FALSE) {
    B <- NULL
  }

  if (C.mat == FALSE) {
    C <- NULL
  }

  # BIND AND OUTPUT
  # -----------------------------------------------------------------

  final <- rbind(A, B, C, out)
  colnames(final) <- params
  return(final)
}
