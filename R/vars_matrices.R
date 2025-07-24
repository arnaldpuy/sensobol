
# CREATE STAR-VARS MATRIX
##################################################################################

#' STAR-VARS sampling strategy
#'
#' It creates the STAR-VARS matrix needed to compute VARS-TO following \insertCite{Razavi2016a;textual}{sensobol}.
#'
#' @param star.centers Positive integer, number of star centers.
#' @param params Character vector with the name of the model inputs.
#' @param h Distance between pairs. The user should select between 0.001, 0.002, 0.005, 0.01,
#' 0.02, 0.05, 0.1, 0.2. Default is \code{h = 0.1}.
#' @param type Approach to construct the STAR-VARS. Options are:
#' * \code{type = "QRN"}: It uses \insertCite{Sobol1967;textual}{sensobol} Quasi-Random Numbers
#' through a call to the function \code{\link[randtoolbox]{sobol}} of the \code{randtoolbox} package.
#' * \code{type = "R"}: It uses random numbers.
#' @param ... Further arguments in \code{\link[randtoolbox]{sobol}}.
#'
#' @return A matrix where each column is a model input and each row a sampling point.
#' @export
#'
#' @details The user randomly selects \eqn{N_{star}} points across the factor space using
#' either Sobol' Quasi Random Numbers (\code{type = "QRN"}) or random numbers (\code{type = "R"}).
#' These are the \emph{star centres} and their location can be denoted as
#' \eqn{\mathbf{s}_v = s_{v_1},...,s_{v_i}, ..., s_{v_k}}, where \eqn{v=1,2,...,N_{star}}.
#' Then, for each star centre, the function generates a cross section of equally spaced points
#' \eqn{\Delta h} apart for each of the \eqn{k} model inputs, including and passing through the
#' star centre. The cross section is produced by fixing \eqn{\mathbf{s}_{v_{\sim i}}} and varying \eqn{s_i}.
#' Finally, for each factor all pairs of points with \eqn{h} values of \eqn{\Delta h, 2\Delta h, 3\Delta h}
#' and so on are extracted. The total computational cost of this design is
#' \eqn{N_t=N_{star} (k (\frac{1}{\Delta h} - 1) + 1)}.
#'
#' @references
#'
#' \insertAllCited{}
#'
#' @examples
#' # Define settings
#' star.centers <- 10; params <- paste("X", 1:5, sep = ""); h <- 0.1
#'
#' # Create STAR-VARS
#' mat <- vars_matrices(star.centers = star.centers, params = params, h = h)
vars_matrices <- function(star.centers, params, h = 0.1, type = "QRN",...) {
  out <- center <- sections <- A <- B <- AB <- X <- out <- list()
  h.values <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2)

  # Ensure that a proper h value has been selected
  # ----------------------------------------------------------------

  if(h %in% h.values == FALSE) {
    stop("Revise the selection of h")
  }

  # Define type of matrix
  # -----------------------------------------------------------------

  if(type == "QRN") {
    mat <- randtoolbox::sobol(n = star.centers, dim = length(params),...)

  } else if(type == "R") {
    mat <- replicate(length(params), stats::runif(star.centers))

  } else {
    stop("Type should be either QRN, R or LHS")
  }

  # Construct STAR_VARS
  # -----------------------------------------------------------------

  for(i in 1:nrow(mat)) {
    center[[i]] <- mat[i, ]
    sections[[i]] <- sapply(center[[i]], function(x) {
      all <- seq(x %% h, 1, h)
      non.zeros <- all[all!= 0]
    })
    B[[i]] <- sapply(1:ncol(mat), function(x)
      sections[[i]][, x][!sections[[i]][, x] %in% center[[i]][x]])
    A[[i]] <- matrix(center[[i]], nrow = nrow(B[[i]]),
                     ncol = length(center[[i]]), byrow = TRUE)
    X[[i]] <- rbind(A[[i]], B[[i]])

    for(j in 1:ncol(A[[i]])) {
      AB[[i]] <- A[[i]]
      AB[[i]][, j] <- B[[i]][, j]
      X[[i]] <- rbind(X[[i]], AB[[i]])
    }
    AB[[i]] <- X[[i]][(2 * nrow(B[[i]]) + 1):nrow(X[[i]]), ]
    out[[i]] <- rbind(unname(center[[i]]), AB[[i]])
  }
  output <- do.call(rbind, out)
  output[output == 1] <- 0.999
  return(output)
}
