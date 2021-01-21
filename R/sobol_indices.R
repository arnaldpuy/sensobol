
# Sobol' indices for a dummy parameter
sobol_dummy_boot <- function(d, i, N, params, boot) {
  k <- length(params)
  if (boot == TRUE) {
    m <- d[i, ]
  } else if (boot == FALSE) {
    m <- d
  }

  # COMPUTATION OF E(Y), V(Y), SI AND TI
  # -----------------------------------------------------------------

  f0 <- (1 / N) * sum(m[, 1] * m[, 2])
  VY <- 1 / (2 * N - 1) * sum(m[, 1]^2 + m[, 2]^2) - f0
  Si <- (1 / (N - 1) * sum(m[, 1] * m[, 2]) - f0) / VY
  STi <- 1 - (1 / (N - 1) * sum(m[, 2] * m[, 2]) - f0) / VY
  return(c(Si, STi))
}

#' Computation of Sobol' indices for a dummy parameter
#'
#' This function computes first and total-order Sobol' indices for a dummy
#' parameter following the formulae shown
#' in \insertCite{KhorashadiZadeh2017;textual}{sensobol}.
#'
#' @param Y A numeric vector with the model output obtained from the matrix created with
#' \code{\link{sobol_matrices}}.
#' @param N Positive integer, the initial sample size of the base sample matrix created with \code{\link{sobol_matrices}}.
#' @param params A character vector with the name of the model inputs.
#' @param boot Logical. If TRUE, the function bootstraps the Sobol' indices. If FALSE, it provides point
#' estimates. Default is \code{boot = FALSE}.
#' @param R Positive integer, number of bootstrap replicas.
#' @param parallel The type of parallel operation to be used (if any).
#' If missing, the default is taken from the option "boot.parallel"
#' (and if that is not set, "no"). For more information, check the
#' \code{parallel} option in the \code{boot} function of the \code{\link{boot}} package.
#' @param ncpus Positive integer: number of processes to be used in parallel operation:
#' typically one would chose this to the number of available CPUs.
#' Check the \code{ncpus} option in the \code{boot} function of the \code{\link{boot}} package.
#' @param conf Confidence intervals, number between 0 and 1. Default is \code{conf = 0.95}.
#' @param type Method to compute the confidence intervals. Default is \code{type = "norm"}.
#' Check the \code{type} option in the \code{boot} function of the \code{\link{boot}} package.
#'
#' @importFrom Rdpack reprompt
#' @importFrom rlang :=
#' @references
#' \insertAllCited{}
#'
#' @return A \code{data.table} object.
#' @export
#'
#' @examples
#' # Define settings
#' N <- 100; params <- paste("X", 1:3, sep = ""); R <- 10
#'
#' # Create sample matrix
#' mat <- sobol_matrices(N = N, params = params)
#'
#' # Compute Ishigami function
#' Y <- ishigami_Fun(mat)
#'
#' # Compute and bootstrap Sobol' indices for dummy parameter
#' ind.dummy <- sobol_dummy(Y = Y, N = N, params = params, boot = TRUE, R = R)

sobol_dummy <- function(Y, N, params, boot = FALSE, R = NULL, parallel = "no",
                        ncpus = 1, conf = 0.95, type = "norm") {
  k <- length(params)
  d <- matrix(Y, nrow = N)

  # COMPUTE WITH BOOTSTRAP
  # -----------------------------------------------------------------

  if (boot == TRUE) {
    out <- boot::boot(data = d, statistic = sobol_dummy_boot,
                      R = R, N = N, params = params,
                      parallel = parallel, ncpus = ncpus,
                      boot = boot)
    out <- bootstats(out, conf = conf, type = type)

  # COMPUTE WITHOUT BOOTSTRAP
  # -----------------------------------------------------------------

  } else if (boot == FALSE) {
    tmp <- sobol_dummy_boot(d = d, N = N, params = params, boot = FALSE)
    out <- data.table::data.table(tmp)
    data.table::setnames(out, "tmp", "original")
  }
  sensitivity <- c("Si", "Ti")
  parameters <- "dummy"
  out <- cbind(out, sensitivity, parameters)

  # TRANSFORM NEGATIVE VALUES TO ZEROES
  # -----------------------------------------------------------------

  colNames <- colnames(out)

  if (any(grepl("high.ci", colNames)) == TRUE) {
    cols_transform <- c("original", "low.ci", "high.ci")

  } else {
    cols_transform <- "original"
  }

  for(j in cols_transform){
    data.table::set(out, i = which(out[[j]] < 0), j = j, value = 0)
  }

  return(out)
}


sobol_boot <- function(d, i, N, params, matrices, R, first, total, order, boot) {

  # STOPPING RULE TO CHECK CONCORDANCE BETWEEN ESTIMATORS AND MATRIX
  # -------------------------------------------------------------------

  ms <- "Revise the correspondence between the matrices and the estimators"

  if (isTRUE(all.equal(matrices, c("A", "B", "AB")))) {
    if (!first == "saltelli" & !first == "jansen" |
       !total == "jansen" & !total == "sobol" & !total == "homma" &
       !total == "janon" & !total == "glen") {
      stop(ms)
    }

  } else if (isTRUE(all.equal(matrices, c("A", "B", "BA")))) {
    if (!first == "sobol"| !total == "saltelli") {
      stop(ms)
    }

  } else if (isTRUE(all.equal(matrices, c("A", "B", "AB", "BA")))) {

    if (!first == "azzini" | !total == "azzini" &
       !total == "jansen" & !total == "sobol" & !total == "homma" &
       !total == "janon" & !total == "glen" & !total == "saltelli") {

      if (!total == "azzini" | !first == "saltelli" & !first == "jansen" &
         !first == "azzini" & !first == "sobol") {

        stop(ms)
      }
    }
  }

  # -------------------------------------

  k <- length(params)
  if (boot == TRUE) {
    m <- d[i, ]
  } else if (boot == FALSE) {
    m <- d
  }
  if (order == "second") {
    k <- length(params) + length(utils::combn(params, 2, simplify = FALSE))

  } else if (order == "third") {
    k <- length(params) + length(utils::combn(params, 2, simplify = FALSE)) +
      length(utils::combn(params, 3, simplify = FALSE))
  }

  # DEFINE VECTORS BASED ON SAMPLE DESIGN
  # ------------------------------------------------------------------

  if (isTRUE(all.equal(matrices, c("A", "B", "AB")))) {
    Y_A <- m[, 1]
    Y_B <- m[, 2]
    Y_AB <- m[, -c(1, 2)]

  } else if (isTRUE(all.equal(matrices, c("A", "B", "BA")))) {
    Y_A <- m[, 1]
    Y_B <- m[, 2]
    Y_BA <- m[, -c(1, 2)]

  } else if (isTRUE(all.equal(matrices, c("A", "B", "AB", "BA")))) {
    Y_A <- m[, 1]
    Y_B <- m[, 2]
    Y_AB <- m[, 3:(k + 2)]
    Y_BA <- m[, (k + 3):ncol(m)]

  } # A warning might be needed here
  if (isTRUE(all.equal(matrices, c("A", "B", "AB"))) |
     isTRUE(all.equal(matrices, c("A", "B", "BA")))) {
    f0 <- 1 / (2 * N) * sum(Y_A + Y_B)
    VY <- 1 / (2 * N - 1) * sum((Y_A - f0)^2 + (Y_B - f0)^2)
  }

  # DEFINE FIRST-ORDER ESTIMATORS
  # --------------------------------------------------------------------

  # Define variance for estimators with A, B, AB; or A, B, BA matrices
  if (first == "saltelli" | first == "jansen" | first == "sobol") {
    f0 <- 1 / (2 * N) * sum(Y_A + Y_B)
    VY <- 1 / (2 * N - 1) * sum((Y_A - f0)^2 + (Y_B - f0)^2)
  }

  # ----------------------------------
  if (first == "sobol") {
    Vi <- 1 / N * Rfast::colsums(Y_A * Y_BA - f0^2)

  } else if (first == "saltelli") {
    Vi <- 1 / N * Rfast::colsums(Y_B * (Y_AB - Y_A))

  } else if (first == "jansen") {
    Vi <- VY - 1 / (2 * N) * Rfast::colsums((Y_B - Y_AB)^2)

  } else if (first == "azzini") {
    VY <- Rfast::colsums((Y_A - Y_B)^2 + (Y_BA - Y_AB)^2)
    Vi <- (2 * Rfast::colsums((Y_BA - Y_B) * (Y_A - Y_AB)))

  } else {
    stop("first should be sobol, saltelli, jansen or azzini")
  }

  if (first == "azzini") {
    Si <- Vi[1:length(params)] / VY[1:length(params)]

  } else {
    Si <- Vi[1:length(params)] / VY
  }

  # DEFINE TOTAL-ORDER ESTIMATORS
  # --------------------------------------------------------------------

  # Define variance for estimators with A, B, AB; or A, B, BA matrices
  if (total == "azzini" | total == "jansen" | total == "sobol" |
     total == "homma" | total == "janon" | total == "glen" | total == "saltelli") {
    f0 <- 1 / (2 * N) * sum(Y_A + Y_B)
    VY <- 1 / (2 * N - 1) * sum((Y_A - f0)^2 + (Y_B - f0)^2)
  }

  # ----------------------------------
  if (total == "jansen") {
    Ti <- (1 / (2 * N) * Rfast::colsums((Y_A - Y_AB)^2)) / VY

  } else if (total == "sobol") {
    Ti <- ((1 / N) * Rfast::colsums(Y_A * (Y_A - Y_AB))) / VY

  } else if (total == "homma") {
    Ti <- (VY - (1 / N) * Rfast::colsums(Y_A * Y_AB) + f0^2) / VY

  } else if (total == "saltelli") {
    Ti <- 1 - ((1 / N * Rfast::colsums(Y_B * Y_BA - f0^2)) / VY)

  } else if (total == "janon") {
    Ti <- 1 - (1 / N * Rfast::colsums(Y_A * Y_AB) -
                 (1/ N * Rfast::colsums((Y_A + Y_AB) / 2))^2) /
      (1 / N * Rfast::colsums((Y_A ^ 2 + Y_AB^2) / 2) -
         (1/ N * Rfast::colsums((Y_A + Y_AB) / 2))^2)

  } else if (total == "glen") {
    Ti <- 1 - (1 / (N - 1) *
                 Rfast::colsums(((Y_A - mean(Y_A)) * (Y_AB - Rfast::colmeans(Y_AB))) /
                                  sqrt(stats::var(Y_A) * Rfast::colVars(Y_AB))))

  } else if (total == "azzini") {
    Ti <- Rfast::colsums((Y_B - Y_BA)^2 + (Y_A - Y_AB)^2) /
      Rfast::colsums((Y_A - Y_B)^2 + (Y_BA - Y_AB)^2)

  } else {
    stop("total should be jansen, sobol, homma saltelli, janon, glen or azzini")
  }
  Ti <- Ti[1:length(params)]

  # DEFINE COMPUTATION OF SECOND-ORDER INDICES
  # ---------------------------------------------------------------------

  if (order == "second" | order == "third") {
    com2 <- utils::combn(1:length(params), 2, simplify = FALSE)
    mat2 <- t(mapply(c, Vi[(length(params) + 1):(length(params) + length(com2))],
                     lapply(com2, function(x) Vi[x])))
    Vij <- apply(mat2, 1, function(x) Reduce("-", x))

    if (first == "azzini") {
      VY <- Rfast::colsums((Y_A - Y_B)^2 + (Y_BA - Y_AB)^2)
      Sij <- Vij / VY[(length(params) + 1):(length(params) + ncol(utils::combn(1:length(params), 2)))]

    } else {
      Sij <- Vij / VY
    }

  } else {
    Sij <- NULL
  }

  # DEFINE COMPUTATION OF THIRD-ORDER INDICES
  # ---------------------------------------------------------------------

  if (order == "third") {
    tmp <- do.call(rbind, com2)
    Vij.vec <- as.numeric(paste(tmp[, 1], tmp[, 2], sep = ""))
    Vij.named <- Vij
    names(Vij.named) <- Vij.vec
    com3 <- utils::combn(1:length(params), 3, simplify = FALSE)
    Vi.only <- do.call(rbind, lapply(com3, function(x) Vi[x])) # Extract Vi, Vj, Vk
    Vijk.only <- utils::tail(Vi, length(com3)) # Extract Vijk
    tmp3 <- do.call(rbind, com3)
    first.pairwise <- lapply(paste(tmp3[, 1], tmp3[, 2], sep = ""), function(x) Vij.named[x])
    second.pairwise <- lapply(paste(tmp3[, 1], tmp3[, 3], sep = ""), function(x) Vij.named[x])
    third.pairwise <- lapply(paste(tmp3[, 2], tmp3[, 3], sep = ""), function(x) Vij.named[x])
    Vij.only <- t(mapply(cbind, first.pairwise, second.pairwise, third.pairwise))
    mat3 <- cbind(Vijk.only, Vij.only, Vi.only)
    Vijk <- apply(mat3, 1, function(x) Reduce("-", x))

    if (first == "azzini") {
      Sijk <- Vijk / utils::tail(VY, length(utils::combn(params, 3, simplify = FALSE)))

    } else {
      Sijk <- Vijk / VY
    }

  } else {
    Sijk <- NULL
  }
  return(c(Si, Ti, Sij, Sijk))
}


bootstats <- function(b, conf = conf, type = type) {
  p <- length(b$t0)
  lab <- c("original", "bias", "std.error", "low.ci", "high.ci")
  tmp <- as.data.frame(matrix(nrow = p,
                              ncol = length(lab),
                              dimnames = list(NULL, lab)))
  for (i in 1:p) {
    # original estimation, bias, standard deviation
    tmp[i, "original"] <- b$t0[i]
    tmp[i, "bias"] <- mean(b$t[, i]) - b$t0[i]
    tmp[i, "std.error"] <- stats::sd(b$t[, i])
    # confidence interval

    if (type == "norm") {
      ci <- boot::boot.ci(b, index = i, type = "norm", conf = conf)

      if (!is.null(ci)) {
        tmp[i, "low.ci"] <- ci$norm[2]
        tmp[i, "high.ci"] <- ci$norm[3]
      }

    } else if (type == "basic") {
      ci <- boot::boot.ci(b, index = i, type = "basic", conf = conf)
      if (!is.null(ci)) {
        tmp[i, "low.ci"] <- ci$basic[4]
        tmp[i, "high.ci"] <- ci$basic[5]
      }

    } else if (type == "percent") {
      ci <- boot::boot.ci(b, index = i, type = "perc", conf = conf)

      if (!is.null(ci)) {
        tmp[i, "low.ci"] <- ci$percent[4]
        tmp[i, "high.ci"] <- ci$percent[5]
      }

    } else if (type == "bca") {
      ci <- boot::boot.ci(b, index = i, conf = conf)

      if (!is.null(ci)) {
        tmp[i, "low.ci"] <- ci$bca[4]
        tmp[i, "high.ci"] <- ci$bca[5]
      }
    }
  }
  return(tmp)
}


#' Computation of Sobol' indices
#'
#' It allows to compute Sobol' indices up to the third order using state-of-the-art estimators.
#'
#'@param matrices Character vector with the required matrices. The default is \code{matrices = c("A", "B", "AB")}.
#' See \code{\link{sobol_matrices}}.
#' @param Y  numeric vector with the model output obtained from the matrix created with
#' \code{\link{sobol_matrices}}.
#' @param N Positive integer, the initial sample size of the base sample matrix created with \code{\link{sobol_matrices}}.
#' @param params Character vector with the name of the model inputs.
#' @param first Estimator to compute first-order indices. Options are:
#' * \code{first = "saltelli"} \insertCite{Saltelli2010a}{sensobol}.
#' * \code{first = "jansen"} \insertCite{Jansen1999}{sensobol}.
#' * \code{first = "sobol"}  \insertCite{Sobol1993}{sensobol}.
#' * \code{first = "azzini"} \insertCite{Azzini2020}{sensobol}.
#' @param total Estimator to compute total-order indices. Options are:
#' * \code{total = "jansen"} \insertCite{Jansen1999}{sensobol}.
#' * \code{total = "sobol"} \insertCite{Sobol2001}{sensobol}.
#' * \code{total = "homma"} \insertCite{Homma1996}{sensobol}.
#' * \code{total = "janon"} \insertCite{Janon2014}{sensobol}.
#' * \code{total = "glen"} \insertCite{Glen2012}{sensobol}.
#' * \code{total = "azzini"} \insertCite{Azzini2020}{sensobol}.
#' * \code{total = "saltelli"} \insertCite{Saltelli2008}{sensobol}.
#' @param order Whether to compute "first", "second", or "third" -order Sobol' indices. Default
#' is \code{order = "first"}.
#' @param boot Logical. If TRUE, the function bootstraps the Sobol' indices. If FALSE, it provides point
#' estimates. Default is \code{boot = FALSE}.
#' @param R Positive integer, number of bootstrap replicas. Default is NULL.
#' @param parallel The type of parallel operation to be used (if any).
#' If missing, the default is taken from the option "boot.parallel"
#' (and if that is not set, "no"). For more information, check the
#' \code{parallel} option in the \code{boot} function of the \code{\link{boot}} package.
#' @param ncpus Positive integer: number of processes to be used in parallel operation:
#' typically one would chose this to the number of available CPUs.
#' Check the \code{ncpus} option in the \code{boot} function of the \code{\link{boot}} package.
#' @param conf Confidence interval if \code{boot = TRUE}. Number between 0 and 1. Default is \code{conf = 0.95}.
#' @param type Method to compute the confidence interval if \code{boot = TRUE}. Default is "norm".
#' Check the \code{type} option in the \code{boot} function of the \code{\link{boot}} package.
#' @importFrom rlang ":="
#' @importFrom Rdpack reprompt
#' @importFrom stats var
#' @references
#' \insertAllCited{}
#'
#' @return A \code{data.table} object.
#' @seealso Check the function \code{\link{boot}} for further details on the bootstrapping
#' with regards to the methods available for the computation of confidence intervals in the \code{type} argument.
#' @export
#'
#' @details Any first and total-order estimator can be combined with the appropriate sampling design.
#' Check Table 3 of the vignette for a summary of all possible combinations, and Tables 1 and 2 for a
#' mathematical description of the estimators. If the analyst mismatches estimators and sampling designs,
#' the function will generate an error and urge to redefine the sample matrices or the estimators.
#'
#' For all estimators except \insertCite{Azzini2020;textual}{sensobol}'s and \insertCite{Janon2014;textual}{sensobol}'s,
#' \code{sobol_indices()} calculates the sample mean as \deqn{\hat{f}_0=\frac{1}{2N} \sum_{v=1}^{N}(f(\mathbf{A})_v + f(\mathbf{B})_v)\,,}
#' where \eqn{N} is the row dimension of the base sample matrix, and the unconditional sample variance as
#'
#' \deqn{\hat{V}(y) = \frac{1}{2N-1} \sum{v=1}^{N} ((f(\mathbf{A})_v - \hat{f})^2 + (f(\mathbf{B})_v - \hat{f})^2)\,,}
#' where \eqn{f(\mathbf{A})_v} (\eqn{f(\mathbf{B})_v}) indicates the model output \eqn{y} obtained after running the model \eqn{f}
#' in the \eqn{v}-th row of the \eqn{\mathbf{A}} (\eqn{\mathbf{B}}) matrix.
#'
#' For the Azzini estimator,
#' \deqn{\hat{V}(y) = \sum_{v=1}^{N} (f(\mathbf{A})_v - f(\mathbf{B})_v)^2 + (f(\mathbf{B}_A^{(i)})_v - f(\mathbf{A}_B^{(i)})_v) ^ 2}
#'
#' and for the Janon estimator,
#' \deqn{\hat{V}(y)=\frac{1}{N} \sum_{v=1}^{N} \frac{f(\mathbf{A})_v^2 + f(\mathbf{A}_B^{(i)})_v^2}{2}-f_0^2}
#'
#'where \eqn{f(\mathbf{A}_B^{(i)})_v} (\eqn{f(\mathbf{B}_A^{(i)})_v}) is the model output obtained after running the model \eqn{f} in
#'the \eqn{v}-th row of an \eqn{\mathbf{A}_B^{(i)})_v} (\eqn{\mathbf{B}_A^{(i)})_v}) matrix, where all columns come from \eqn{\mathbf{A}} (\eqn{\mathbf{B}})
#'except the \eqn{i}-th, which comes from \eqn{\mathbf{B}} (\eqn{\mathbf{A}}).
#'
#' @examples
#' # Define settings
#' N <- 1000; params <- paste("X", 1:3, sep = ""); R <- 10
#'
#' # Create sample matrix
#' mat <- sobol_matrices(N = N, params = params)
#'
#' # Compute Ishigami function
#' Y <- ishigami_Fun(mat)
#'
#' # Compute and bootstrap Sobol' indices
#' ind <- sobol_indices(Y = Y, N = N, params = params, boot = TRUE, R = R)
sobol_indices <- function(matrices = c("A", "B", "AB"), Y, N, params,
                          first = "saltelli", total = "jansen",
                          order = "first", boot = FALSE, R = NULL,
                          parallel = "no", ncpus = 1, conf = 0.95, type = "norm") {

  # CHECK CONCORDANCE BETWEEN BOOT AND R ARGUMENTS
  # ---------------------------------------------------------------------

  if (boot == FALSE & is.null(R) == FALSE | boot == TRUE & is.null(R) == TRUE) {
    stop("Bootstrapping requires boot = TRUE and an integer in R")
  }

  # DEFINE PARAMETERS
  # ----------------------------------------------------------------------

  sensitivity <- parameters <- NULL
  k <- length(params)
  d <- matrix(Y, nrow = N)

  # FUNCTION WHEN BOOT = FALSE
  # -----------------------------------------------------------------------

  if (boot == FALSE) {
    tmp <- sobol_boot(d = d, N = N, params = params, first = first, total = total,
                      order = order, boot = FALSE, matrices = matrices)
    out <- data.table::data.table(tmp)
    data.table::setnames(out, "tmp", "original")

    # FUNCTION WHEN BOOT = TRUE
    # -----------------------------------------------------------------------

  } else if (boot == TRUE) {
    tmp <- boot::boot(data = d, statistic = sobol_boot, R = R, N = N, params = params,
                      first = first, total = total, order = order, matrices = matrices,
                      parallel = parallel, ncpus = ncpus, boot = TRUE)
    out <- data.table::data.table(bootstats(tmp, conf = conf, type = type))

  } else {
    stop("boot has to be TRUE or FALSE")
  }

  # VECTORS OF PARAMETERS AND SENSITIVITY INDICES WHEN ORDER = FIRST
  # -----------------------------------------------------------------------

  if (order == "first") {
    parameters <- c(rep(params, times = 2))
    sensitivity <- c(rep(c("Si", "Ti"), each = k))

    # VECTORS OF PARAMETERS AND SENSITIVITY INDICES WHEN ORDER = SECOND
    # -----------------------------------------------------------------------

  } else if (order == "second") {
    vector.second <- unlist(lapply(utils::combn(params, 2, simplify = FALSE), function(x)
      paste0(x, collapse = ".")))
    parameters <- c(c(rep(params, times = 2)), vector.second)
    sensitivity <- c(rep(c("Si", "Ti"), each = length(params)),
                     rep("Sij", times = length(vector.second)))

    # VECTORS OF PARAMETERS AND SENSITIVITY INDICES WHEN ORDER = THIRD
    # -----------------------------------------------------------------------

  } else if (order == "third") {
    vector.second <- unlist(lapply(utils::combn(params, 2, simplify = FALSE), function(x)
      paste0(x, collapse = ".")))
    parameters <- c(c(rep(params, times = 2)), vector.second)
    vector.third <- unlist(lapply(utils::combn(params, 3, simplify = FALSE), function(x)
      paste0(x, collapse = ".")))
    parameters <- c(parameters, vector.third)
    sensitivity <- c(rep(c("Si", "Ti"), each = k),
                     rep("Sij", times = length(vector.second)),
                     rep("Sijk", times = length(vector.third)))

  } else {
    stop("order has to be first, second or third")
  }
  out <- cbind(out, sensitivity, parameters)
  return(out)
}



