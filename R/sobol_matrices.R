
# RESOLUTION OF GROUPS INTO COLUMN INDICES
##################################################################################

#' @keywords internal
#' @noRd
resolve_groups <- function(params, groups) {

  # No groups: each parameter is its own (single-element) group
  # -----------------------------------------------------------------

  if (is.null(groups)) {
    g <- as.list(seq_along(params))
    names(g) <- params
    return(g)
  }

  # Vector form -> normalize to named list of character vectors
  # -----------------------------------------------------------------

  if (is.atomic(groups) && !is.list(groups)) {
    if (length(groups) != length(params))
      stop("'groups' as a vector must have the same length as 'params'.")
    grp <- as.character(groups)
    groups <- split(params, factor(grp, levels = unique(grp)))
  }

  # List form: validate strict partition
  # -----------------------------------------------------------------

  if (!is.list(groups))
    stop("'groups' must be either NULL, a character vector, or a named list.")

  if (is.null(names(groups)) || any(!nzchar(names(groups))) ||
      anyDuplicated(names(groups)))
    stop("'groups' must be a named list with unique, non-empty group labels.")

  if (!all(vapply(groups, is.character, logical(1))))
    stop("All entries of 'groups' must be character vectors of parameter names.")

  flat <- unlist(groups, use.names = FALSE)

  if (anyDuplicated(flat))
    stop(paste("'groups' must be a strict partition of 'params':",
               "each parameter must appear in exactly one group."))

  if (!setequal(flat, params))
    stop(paste("'groups' must be a strict partition of 'params':",
               "every parameter in 'params' must appear in exactly one group."))

  lapply(groups, function(g) match(g, params))
}

# CREATION OF THE AB, BA, CB matrices
##################################################################################

scrambled_sobol <- function(matrices, A, B, C, order, group_idx) {
  G <- length(group_idx)
  N <- nrow(A)

  # Vectors with the (combinations of) group indices
  # -----------------------------------------------------------------

  first <- as.list(seq_len(G))

  if (order == "first") {
    loop <- first

  } else if (order == "second") {
    loop <- c(first, utils::combn(seq_len(G), 2, simplify = FALSE))

  } else if (order == "third") {
    loop <- c(first,
              utils::combn(seq_len(G), 2, simplify = FALSE),
              utils::combn(seq_len(G), 3, simplify = FALSE))

  } else if (order == "fourth") {
    loop <- c(first,
              utils::combn(seq_len(G), 2, simplify = FALSE),
              utils::combn(seq_len(G), 3, simplify = FALSE),
              utils::combn(seq_len(G), 4, simplify = FALSE))

  } else {

    stop("order should be either first, second, third or fourth")
  }

  # Check which matrices have been selected
  # -----------------------------------------------------------------

  AB.mat <- "AB" %in% matrices
  BA.mat <- "BA" %in% matrices
  CB.mat <- "CB" %in% matrices

  # Construct Ab, Ba matrices, etc
  # -----------------------------------------------------------------

  if (AB.mat == TRUE) {
    X <- rbind(A, B)

    for (gset in loop) {
      cols <- unlist(group_idx[gset], use.names = FALSE)
      AB <- A
      AB[, cols] <- B[, cols]
      X <- rbind(X, AB)
    }
    AB <- X[(2 * N + 1):nrow(X), , drop = FALSE]

  } else if (AB.mat == FALSE) {
    AB <- NULL
  }

  if (BA.mat == TRUE) {
    W <- rbind(A, B)

    for (gset in loop) {
      cols <- unlist(group_idx[gset], use.names = FALSE)
      BA <- B
      BA[, cols] <- A[, cols]
      W <- rbind(W, BA)
    }
    BA <- W[(2 * N + 1) : nrow(W), , drop = FALSE]

  } else if (BA.mat == FALSE) {
    BA <- NULL
  }

  if (CB.mat == TRUE) {
    Z <- rbind(A, B)

    for (gset in loop) {
      cols <- unlist(group_idx[gset], use.names = FALSE)
      CB <- C
      CB[, cols] <- B[, cols]
      Z <- rbind(Z, CB)
    }
    CB <- Z[(2 * N + 1) : nrow(Z), , drop = FALSE]

  } else if (CB.mat == FALSE) {
    CB <- NULL
  }

  # Merge and output
  # -----------------------------------------------------------------

  final <- rbind(AB, BA, CB)
  return(final)
}

# FUNCTION TO CREATE THE SAMPLE MATRICES
##################################################################################

#' Creation of the sample matrices
#'
#' It creates the sample matrices to compute Sobol' first and total-order indices.
#' If needed, it also creates the sample matrices required to compute second,
#' third and fourth-order indices. Indices can be computed for individual parameters
#' or for groups of parameters (e.g. to handle correlated inputs by moving them
#' together).
#'
#' @param matrices Character vector with the required matrices. The default
#' is \code{matrices = c("A", "B", "AB")}.
#' @param N Positive integer, initial sample size of the base sample matrix.
#' @param params Character vector with the name of the model inputs.
#' @param order One of "first", "second", "third" or "fourth" to create a matrix to
#' compute first, second, third or up to fourth-order Sobol' indices. The default is
#' \code{order = "first"}. When \code{groups} is supplied, the order refers to
#' interactions between groups, not between individual parameters.
#' @param type Approach to construct the sample matrix. Options are:
#' * \code{type = "QRN"} (default): It uses \insertCite{Sobol1967;textual}{sensobol} Quasi-Random Numbers.
#' through a call to the function \code{\link[randtoolbox]{sobol}} of the \code{randtoolbox} package.
#' * \code{type = "LHS"}: It uses a Latin Hypercube Sampling Design
#' \insertCite{McKay1979}{sensobol} through a call
#' to the function \code{\link[lhs]{randomLHS}} of the \code{lhs} package.
#' * \code{type = "R"}: It uses random numbers.
#' @param groups Optional grouping of the parameters into a strict partition.
#' Either \code{NULL} (the default; each parameter is its own group, i.e.
#' the standard Sobol' design), a named list of character vectors of parameter
#' names (e.g. \code{groups = list(g1 = "X1", g2 = c("X2", "X3"))}), or a
#' character vector of length \code{length(params)} aligned with \code{params}
#' that labels each parameter with its group (e.g.
#' \code{groups = c("g1", "g2", "g2")}). Every parameter in \code{params} must
#' appear in exactly one group. When \code{groups} is supplied, conditional
#' matrices swap all columns belonging to a group simultaneously, which is the
#' standard device used to handle correlated inputs in Sobol'-based sensitivity
#' analysis by moving them together.
#' @param ... Further arguments in \code{\link[randtoolbox]{sobol}}.
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
#' The same process applies to create the matrices to compute fourth-order effects.
#'
#' When \code{groups} is supplied, the construction is identical except that the
#' index \eqn{i} (or \eqn{(i,j)}, \eqn{(i,j,l)}, \eqn{(i,j,l,m)}) ranges over
#' groups rather than individual parameters: a single conditional matrix
#' \eqn{\mathbf{A}_B^{(g)}} swaps all columns of group \eqn{g} from
#' \eqn{\mathbf{A}} to \eqn{\mathbf{B}} simultaneously. The number of conditional
#' matrices appended therefore depends on the number of groups \eqn{G} rather
#' than on \eqn{k}: \eqn{G} for first-order, \eqn{\binom{G}{2}} for second-order,
#' and so on. Downstream, \code{\link{sobol_indices}} interprets each conditional
#' block as a group-wise effect.
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
#'
#' # Grouped design: treat X2 and X3 as a single group (e.g. because they
#' # are correlated and should be moved together).
#' params <- paste("X", 1:4, sep = "")
#' groups <- list(g1 = "X1", g2 = c("X2", "X3"), g3 = "X4")
#' mat_g <- sobol_matrices(N = 100, params = params, groups = groups)
sobol_matrices <- function(matrices = c("A", "B", "AB"),
                           N, params, order = "first",
                           type = "QRN", groups = NULL, ...) {

  if (!is.numeric(N) || length(N) != 1 || N < 1 || N != floor(N))
    stop("'N' must be a single positive integer.")

  k <- length(params)
  n.matrices <- ifelse(any(stringr::str_detect(matrices, "C")) == FALSE, 2, 3)

  # Resolve and validate groups (NULL -> one group per parameter)
  # -----------------------------------------------------------------

  group_idx <- resolve_groups(params, groups)

  # Selection of the sample matrix type
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

  # Construction of A, B, ... matrices
  # -----------------------------------------------------------------

  A <- df[, 1:k, drop = FALSE]
  B <- df[, (k + 1) : (k * 2), drop = FALSE]

  if (n.matrices == 3) {
    C <- df[, ((k * 2) + 1):(k * 3), drop = FALSE]

  } else {
    C <- NULL
  }

  # Construction of Ab, Ba, ... matrices
  # -----------------------------------------------------------------

  out <- scrambled_sobol(matrices = matrices,
                         A = A, B = B, C = C,
                         order = order,
                         group_idx = group_idx)
  A.mat <- "A" %in% matrices
  B.mat <- "B" %in% matrices
  C.mat <- "C" %in% matrices

  # Set NULL to those not used
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

  # Bind and output
  # -----------------------------------------------------------------

  final <- rbind(A, B, C, out)
  colnames(final) <- params
  return(final)
}
