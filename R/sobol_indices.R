
# Sobol' indices for a dummy parameter
sobol_dummy_boot <- function(d, i, N, params, boot) {
  k <- length(params)
  if(boot == TRUE) {
    m <- d[i, ]
  } else if(boot == FALSE) {
    m <- d
  }
  f0 <- (1 / N) * sum(m[, 1] * m[, 2])
  VY <- 1 / (2 * N - 1) * sum(m[, 1] ^ 2 + m[, 2] ^ 2) - f0
  Si <- (1 / (N - 1) * sum(m[, 1] * m[, 2]) - f0) / VY
  STi <- 1 - (1 / (N - 1) * sum(m[, 2] * m[, 2]) - f0) / VY
  return(c(Si, STi))
}

#' Computation of Sobol' indices for a dummy parameter
#'
#' This function computes first and total-order Sobol' indices for a dummy
#' parameter following the formulas shown
#' in \insertCite{KhorashadiZadeh2017;textual}{sensobol}.
#'
#' @param Y Numeric vector, model output.
#' @param N Integer, base sample size of the sample matrix created with \code{\link{sobol_matrices}}.
#' @param params Vector with the name of the model inputs.
#' @param boot Logical. Default is \code{boot = FALSE}.
#' @param R Integer, number of bootstrap replicas.
#' @param parallel The type of parallel operation to be used (if any).
#' If missing, the default is taken from the option "boot.parallel"
#' (and if that is not set, "no"). For more information, check the
#' \code{parallel} option in the \code{boot} function of the \code{\link{boot}} package.
#' @param ncpus Integer: number of processes to be used in parallel operation:
#' typically one would chose this to the number of available CPUs.
#' Check the \code{ncpus} option in the \code{boot} function of the \code{\link{boot}} package.
#' @param conf Number between 0 and 1. Default is \code{conf = 0.95}.
#' @param type Method to compute the confidence intervals. Default is \code{type = "norm"}.
#' Check the \code{type} option in the \code{boot} function of the \code{\link{boot}} package.
#'
#' @importFrom Rdpack reprompt
#' @importFrom rlang :=
#' @references
#' \insertAllCited{}
#'
#' @return A data.table object.
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
  if(boot == TRUE) {
    out <- boot::boot(data = d, statistic = sobol_dummy_boot,
                      R = R, N = N, params = params,
                      parallel = parallel, ncpus = ncpus,
                      boot = boot)
    out <- bootstats(out, conf = conf, type = type)
  } else if(boot == FALSE) {
    tmp <- sobol_dummy_boot(d = d, N = N, params = params, boot = FALSE)
    out <- data.table::data.table(tmp)
    data.table::setnames(out, "tmp", "original")
  }
  sensitivity <- c("Si", "Ti")
  parameters <- "dummy"
  out <- cbind(out, sensitivity, parameters)
  return(out)
}

sobol_boot <- function(d, i, N, params, R, first, total, order, boot) {
  if(first == "janon" & !total == "janon" | !first == "janon" & total == "janon") {
    stop("janon should be used simultaneously in first and total indices
         with an A, AB and BA matrices")
  }
  if(first == "azzini" & !total == "azzini" | !first == "azzini" & total == "azzini") {
    stop("azzini should be used simultaneously in first and total indices
         with an A, B, AB and BA matrices")
  }
  k <- length(params)
  if(boot == TRUE) {
    m <- d[i, ]
  } else if(boot == FALSE) {
    m <- d
  }
  if(order == "second") {
    k <- length(params) + length(utils::combn(params, 2, simplify = FALSE))
  } else if(order == "third") {
    k <- length(params) + length(utils::combn(params, 2, simplify = FALSE)) +
      length(utils::combn(params, 3, simplify = FALSE))
  }
  # DEFINE FIRST-ORDER ESTIMATORS -----------------------
  if(first == "jansen" | first == "saltelli" &
     total == "jansen" | total == "sobol" | total == "homma") {
    Y_A <- m[, 1]
    Y_B <- m[, 2]
    Y_AB <- m[, -c(1, 2)]
    f0 <- 1 / (2 * N) * sum(Y_A + Y_B)
    VY <- 1 / (2 * N - 1) * sum((Y_A - f0) ^ 2 + (Y_B - f0) ^ 2)
    if(first == "jansen") {
      Vi <- VY - 1 / (2 * N) * Rfast::colsums((Y_B - Y_AB) ^ 2)
      Si <- Vi[1:length(params)] / VY
    } else if(first == "saltelli") {
      Vi <- 1 / N * Rfast::colsums(Y_B * (Y_AB - Y_A))
      Si <- Vi[1:length(params)] / VY
    }
  }
  if(first == "janon" | total == "janon") {
    Y_A <- m[, 1]
    Y_AB <- m[, 2:(k + 1)]
    Y_BA <- m[,(k + 2):ncol(m)]
    Vi <- 1 / N * Rfast::colsums(Y_A * Y_BA - (1 / (2 * N) * Rfast::colsums(Y_A + Y_BA)) ^ 2)
    VY <- (1 / (2 * N) * Rfast::colsums(Y_A ^ 2 + Y_BA ^ 2) -
             (1 / (2 * N) * Rfast::colsums(Y_A + Y_BA)) ^ 2)
    Si <- Vi[1:length(params)] / VY[1:length(params)]
    STi <- 1 - (1 / N * Rfast::colsums(Y_A * Y_AB) -
                  (1/ N * Rfast::colsums((Y_A + Y_AB) / 2)) ^ 2) /
      (1 / N * Rfast::colsums((Y_A ^ 2 + Y_AB ^ 2) / 2) -
         (1/ N * Rfast::colsums((Y_A + Y_AB) / 2)) ^ 2)
    STi <- STi[1:length(params)]
  }
  if(first == "azzini" | total == "azzini") {
    Y_A <- m[, 1]
    Y_B <- m[, 2]
    Y_AB <- m[, 3:(3 + k - 1)]
    Y_BA <- m[, (ncol(m) - k + 1):ncol(m)]
    Vi <- 1 / N * Rfast::colsums(Y_A * Y_BA) - ((1 / N) * sum(Y_A * Y_B)) +
      (1 / N) * Rfast::colsums(Y_B * Y_AB) - ((1 / N) * Rfast::colsums(Y_AB * Y_BA))
    VY <- (1 / (2 * N ) * (sum(Y_B * (Y_B - Y_A)) + sum(Y_A * (Y_A - Y_B)) +
                             Rfast::colsums(Y_BA * (Y_BA - Y_AB)) +
                             Rfast::colsums(Y_AB * (Y_AB - Y_BA))))
    Si <- Vi[1:length(params)] / VY[1:length(params)]
    STi <- 1 - abs(Rfast::colsums((Y_A - Y_BA) * (Y_B - Y_AB)) /
                     (1 / 2 * Rfast::colsums((Y_A - Y_B) ^ 2 + (Y_AB - Y_BA) ^ 2)))
    STi <- STi[1:length(params)]
  }
  # DEFINE TOTAL-ORDER ESTIMATORS FOR THE REST ----------------
  if(total == "jansen") {
    STi <- (1 / (2 * N) * Rfast::colsums((Y_A - Y_AB) ^ 2)) / VY
  } else if(total == "homma") {
    STi <- (VY - (1 / N) * Rfast::colsums(Y_A * Y_AB) + f0 ^ 2) / VY
  } else if(total == "sobol") {
    STi <- ((1 / N) * Rfast::colsums(Y_A * (Y_A - Y_AB))) / VY
  }
  STi <- STi[1:length(params)]
  if(order == "second" | order == "third") {
    com2 <- utils::combn(1:length(params), 2, simplify = FALSE)
    mat2 <- t(mapply(c, Vi[(length(params) + 1):(length(params) + length(com2))], lapply(com2, function(x) Vi[x])))
    Vij <- apply(mat2, 1, function(x) Reduce("-", x))
    if(first == "janon" | first == "azzini") {
      Sij <- Vij / VY[(length(params) + 1):(length(VY) - length(utils::combn(params, 3, simplify = FALSE)))]
    } else {
      Sij <- Vij / VY
    }
  }
  if(order == "first") {
    Sij <- NULL
  }
  if(order == "third") {
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
    if(first =="janon" | first == "azzini") {
      Sijk <- Vijk / utils::tail(VY, length(utils::combn(params, 3, simplify = FALSE)))
    } else {
      Sijk <- Vijk / VY
    }
  } else {
    Sijk <- NULL
  }
  return(c(Si, STi, Sij, Sijk))
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
#' @param Y Numeric vector, model output.
#' @param N Integer, base sample size of the sample matrix created with \code{\link{sobol_matrices}}.
#' @param params Vector with the name of the model inputs.
#' @param first Estimator used to compute first order indices. Options are:
#' * \code{first = "saltelli"} \insertCite{Saltelli2010a}{sensobol} (Default).
#' * \code{first = "jansen"} \insertCite{Jansen1999}{sensobol}.
#' * \code{first = "janon"}  \insertCite{Janon2014}{sensobol}.
#' * \code{first = "azzini"}.
#' @param total Estimator used to compute total order indices. Options are:
#' * \code{total = "jansen"} \insertCite{Jansen1999}{sensobol} (Default).
#' * \code{total = "sobol"} \insertCite{Sobol2001}{sensobol}.
#' * \code{total = "homma"} \insertCite{Homma1996}{sensobol}.
#' * \code{total = "janon"} \insertCite{Janon2014}{sensobol}.
#' * \code{total = "azzini"}.
#' @param order Whether to compute "first", "second", or "third" order Sobol' indices. Default
#' is \code{order = "first"}.
#' @param boot Logical. If TRUE, bootstraps the Sobol' indices. If FALSE, it provides point
#' estimates. Default is \code{boot = FALSE}.
#' @param R Integer, number of bootstrap replicas. Default is NULL.
#' @param parallel The type of parallel operation to be used (if any).
#' If missing, the default is taken from the option "boot.parallel"
#' (and if that is not set, "no"). For more information, check the
#' \code{parallel} option in the \code{boot} function of the \code{\link{boot}} package.
#' @param ncpus Integer: number of processes to be used in parallel operation:
#' typically one would chose this to the number of available CPUs.
#' Check the \code{ncpus} option in the \code{boot} function of the \code{\link{boot}} package.
#' @param conf Confidence interval if \code{boot = TRUE}. Number between 0 and 1. Default is \code{conf = 0.95}.
#' @param type Method to compute the confidence interval if \code{boot = TRUE}. Default is "norm".
#' Check the \code{type} option in the \code{boot} function of the \code{\link{boot}} package.
#' @importFrom rlang ":="
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited{}
#'
#' @return A data.table object.
#' @seealso Check the function \code{\link{boot}} for further details on the bootstrapping
#' with regards to the methods available for the computation of confidence intervals in \code{type}.
#' @export
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
sobol_indices <- function(Y, N, params, first = "saltelli", total = "jansen",
                          order = "first", boot = FALSE, R = NULL, parallel = "no",
                          ncpus = 1, conf = 0.95, type = "norm") {
  if(boot == FALSE & is.null(R) == FALSE | boot == TRUE & is.null(R) == TRUE) {
    stop("Bootstrapping requires boot = TRUE and an integer in R")
  }
  sensitivity <- parameters <- NULL
  k <- length(params)
  d <- matrix(Y, nrow = N)
  if(boot == FALSE) {
    tmp <- sobol_boot(d = d, N = N, params = params, first = first, total = total,
                      order = order, boot = FALSE)
    out <- data.table::data.table(tmp)
    data.table::setnames(out, "tmp", "original")
  } else if(boot == TRUE) {
    tmp <- boot::boot(data = d, statistic = sobol_boot, R = R, N = N, params = params,
                      first = first, total = total, order = order,
                      parallel = parallel, ncpus = ncpus, boot = TRUE)
    out <- data.table::data.table(bootstats(tmp, conf = conf, type = type))
  } else {
    stop("boot has to be TRUE or FALSE")
  }
  if(order == "first") {
    parameters <- c(rep(params, times = 2))
    sensitivity <- c(rep(c("Si", "Ti"), each = k))
  } else if(order == "second") {
    vector.second <- unlist(lapply(utils::combn(params, 2, simplify = FALSE), function(x)
      paste0(x, collapse = ".")))
    parameters <- c(c(rep(params, times = 2)), vector.second)
    sensitivity <- c(rep(c("Si", "Ti"), each = length(params)),
                 rep("Sij", times = length(vector.second)))
  } else if(order == "third") {
    vector.second <- unlist(lapply(utils::combn(params, 2, simplify = FALSE), function(x)
      paste0(x, collapse = ".")))
    parameters <- c(c(rep(params, times = 2)), vector.second)
    vector.third <- unlist(lapply(utils::combn(params, 3, simplify = FALSE), function(x)
      paste0(x, collapse = ".")))
    parameters <- c(parameters, vector.third)
    sensitivity <- c(rep(c("Si", "Ti"), each = k),
                 rep("Sij", times = length(vector.second)),
                 rep("Sijk", times = length(vector.third)))
  }
  out <- cbind(out, sensitivity, parameters)
  return(out)
}


