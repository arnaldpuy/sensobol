
scrambled_sobol <- function(matrices, A, B, C, order) {
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
#' third-order indices.
#'
#' @param matrices Vector with the required matrices. The default
#' is \code{matrices = c("A", "B", "AB")}.
#' @param N Integer, initial ample size of the base sample matrix.
#' @param params Character vector with the name of the model inputs.
#' @param order One of "first", "second" or "third" to create a matrix to
#' compute first, second or up to third-order Sobol indices. The default is
#' \code{order = "first"}.
#' @param type Approach to construct the sample matrix. Options are:
#' * \code{type = "QRN"}: It uses Sobol' Quasi-Random Numbers \insertCite{Sobol1967}{sensobol}
#' through a call to the function \code{\link{sobol}} of the \code{randtoolbox} package.
#' * \code{type = "LHS"}: It uses a Latin Hypercube Sampling Design
#' \insertCite{McKay1979}{sensobol} through a call
#' to the function \code{\link{randomLHS}} of the \code{lhs} package.
#' * \code{type = "R"}: It uses random numbers.
#' @param ... Further arguments in \code{\link{sobol}}.
#' @return A matrix where each column is a model input and each row a sampling point.
#' @export

#' @details Before calling \code{sobol_matrices}, the user must decide which estimators
#' will be used to compute first and total-order indices, for this option conditions
#' the design of the sample matrix and therefore the argument \code{matrices}.
#' See the vignette for further details on the specific sampling designs required by
#' the estimators.

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
                           type = "QRN", ...) {
  k <- length(params)
  n.matrices <- ifelse(any(stringr::str_detect(matrices, "C")) == FALSE, 2, 3)
  if(type == "QRN") {
    df <- randtoolbox::sobol(n = N, dim = k * n.matrices, ...)
  } else if(type == "R") {
    df <- replicate(k * n.matrices, stats::runif(N))
  } else if(type == "LHS") {
    df <- lhs::randomLHS(N, n.matrices * k)
  } else {
    stop("method should be either QRN, R or LHS")
  }
  A <- df[, 1:k]
  B <- df[, (k + 1) : (k * 2)]
  if(n.matrices == 3) {
    C <- df[, ((k * 2) + 1):(k * 3)]
  } else {
    C <- NULL
  }
  out <- scrambled_sobol(matrices = matrices,
                         A = A, B = B, C = C,
                         order = order)
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
