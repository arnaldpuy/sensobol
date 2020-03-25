
scrambled_sobol <- function(matrices, A, B, C, order, cluster) {
  first <- 1:ncol(A)
  N <- nrow(A)
  if(order == "first") {
    loop <- first
  } else if(order == "second") {
    second <- c(first, utils::combn(1:ncol(A), 2, simplify = FALSE))
    loop <- second
  } else if(order == "third") {
    second <- c(first, utils::combn(1:ncol(A), 2, simplify = FALSE))
    third <- c(second, utils::combn(1:ncol(A), 3, simplify = FALSE))
    loop <- third
  } else {
    stop("order should be either first, second or third")
  }
  if(is.null(cluster) == FALSE) {
    loop <- cluster
  }
  AB.mat <- "AB" %in% matrices
  BA.mat <- "BA" %in% matrices
  CB.mat <- "CB" %in% matrices
  if(AB.mat == TRUE) {
    X <- rbind(A, B)
    for(i in loop) {
      AB <- A
      AB[, i] <- B[, i]
      X <- rbind(X, AB)
    }
    AB <- X[(2 * N + 1):nrow(X), ]
  } else if(AB.mat == FALSE) {
    AB <- NULL
  }
  if(BA.mat == TRUE) {
    W <- rbind(A, B)
    for(i in loop) {
      BA <- B
      BA[, i] <- A[, i]
      W <- rbind(W, BA)
    }
    BA <- W[(2 * N + 1) : nrow(W), ]
  } else if(BA.mat == FALSE) {
    BA <- NULL
  }
  if(CB.mat == TRUE) {
    Z <- rbind(A, B)
    for(i in loop) {
      CB <- C
      CB[, i] <- B[, i]
      Z <- rbind(Z, CB)
    }
    CB <- Z[(2 * N + 1) : nrow(Z), ]
  } else if(CB.mat == FALSE) {
    CB <- NULL
  }
  final <- rbind(AB, BA, CB)
  return(final)
}

#' Creation of the sample matrices
#'
#' It creates the sample matrices to compute Sobol' first and total-order indices.
#' If needed, it also creates the sample matrices required to compute second and
#' third-order indices. It uses \insertCite{Sobol1967;textual}{sensobol} quasi-random number sequences.
#'
#' @param matrices Vector with the required matrices. The default
#' is \code{matrices = c("A", "B", "AB")}.
#' @param N Integer, initial ample size of the base Sobol' matrix.
#' @param params Vector with the name of the model inputs.
#' @param order One of "first", "second" or "third" to create a matrix to
#' compute first, second or up to third-order Sobol' indices. The default is
#' \code{order = "first"}.
#' @param cluster List of vectors including the model inputs that form part of the
#' cluster/s. The default is \code{cluster = NULL}

#' @seealso Check the function \code{\link{sobol}} in the package \code{randtoolbox}
#' to see how the Sobol' quasi-random number sequences are constructed.
#' @return A matrix.
#' @export

#' @details The function generates a \eqn{(N, 2k)} matrix using Sobol' quasi-random
#'    number sequences, for \eqn{i=1,2,...,k} parameters. The first \emph{k}-matrix is
#'    the \strong{A} matrix and the remaining \emph{k}-matrix, the \strong{B}
#'    matrix. If \code{matrices} includes an "AB" ("BA") matrix and \code{order = "first"},
#'    the function also generates \emph{k} \strong{A}_B^{j} (\strong{B}_A^{j})
#'    matrices, where all collumns come from \strong{A} (\strong{B}) except the
#'    \emph{j}-th, which comes from \strong{B} (\strong{A}).
#'
#'  If \code{matrices} includes an "AB" ("BA") matrix and \code{order = "second"}
#'    or \code{order = "third"}, the function will also create \eqn{n} extra
#'    matrices, where \eqn{n = k! / 2!(k-2)!} or \eqn{n = k! /3!(k-3)!} respectively.
#'    These are needed to compute second or third-order Sobol' indices.
#'
#'  The argument \code{matrices} can include either a single element
#'   (i.e. \code{matrices = "A"}, " \code{matrices = "B"},  \code{matrices = "AB"}
#'   or \code{matrices = "BA"}), or a vector with a combination of matrices;
#'   for instance, \code{matrices = c("A", "AB")}, \code{matrices = c("A", "AB", "BA")},
#'   etc. The selection of the estimator to compute first and total-order indices
#'   will define which matrices are required:

#' * If the estimator of choice in \code{\link{sobol_indices}} is the Azzini estimator,
#'   \code{matrices = c("A", "B", "AB", "BA")}.
#'
#' * If the estimator in \code{\link{sobol_indices}} is the \insertCite{Janon2014;textual}{sensobol},
#'   \code{matrices = c("A", "AB", "BA")}.
#'
#' * If the estimator in \code{\link{sobol_indices}} for first-order indices is either the
#'   \insertCite{Saltelli2010a;textual}{sensobol} or the \insertCite{Jansen1999;textual}{sensobol},
#'   and the estimator for total order indices is either the
#'   \insertCite{Homma1996;textual}{sensobol}, the
#'   \insertCite{Jansen1999;textual}{sensobol} or the
#'   \insertCite{Sobol2001;textual}{sensobol}, \code{matrices = c("A", "B", "AB")}.
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
#' # Create sample matrix
#' mat <- sobol_matrices(N = N, params = params, order = order)
sobol_matrices <- function(matrices = c("A", "B", "AB"),
                           N, params, order = "first",
                           cluster = NULL) {
  k <- length(params)
  n.matrices <- ifelse(any(stringr::str_detect(matrices, "C")) == FALSE, 2, 3)
  df <- randtoolbox::sobol(n = N, dim = k * n.matrices)
  A <- df[, 1:k]
  B <- df[, (k + 1) : (k * 2)]
  if(n.matrices == 3) {
    C <- df[, ((k * 2) + 1):(k * 3)]
  } else {
    C <- NULL
  }
  out <- scrambled_sobol(matrices = matrices,
                         A = A, B = B, C = C,
                         order = order, cluster = cluster)
  A.mat <- "A" %in% matrices
  B.mat <- "B" %in% matrices
  C.mat <- "C" %in% matrices
  if(A.mat == FALSE) {
    A <- NULL
  }
  if(B.mat == FALSE) {
    B <- NULL
  }
  if(C.mat == FALSE) {
    C <- NULL
  }
  final <- rbind(A, B, C, out)
  colnames(final) <- params
  return(final)
}
