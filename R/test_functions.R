
# SOBOL' G FUNCTION
##################################################################################

#' Sobol' G function
#'
#' It implements the \insertCite{Sobol1998;textual}{sensobol} G function.
#'
#' @param X A data frame or numeric matrix.
#'
#' @return A numeric vector with the model output.
#' @export
#'
#' @references
#' \insertAllCited{}
#'
#' @details The function requires eight model inputs and reads as
#' \deqn{y=\prod_{i=1}^{k} \frac{|4 x_i - 2| + a_i}{1 + a_i}\,,}
#' where \eqn{k=8}, \eqn{x_i\sim\mathcal{U}(0,1)} and \eqn{a=(0, 1, 4.5, 9, 99, 99, 99, 99)}.
#'
#' @examples
#' # Define settings
#' N <- 100; params <- paste("X", 1:8, sep = "")
#'
#' # Create sample matrix
#' mat <- sobol_matrices(N = N, params = params)
#'
#' # Compute Sobol' G
#' Y <- sobol_Fun(mat)
sobol_Fun <- function(X) {
  a <- c(0, 1, 4.5, 9, 99, 99, 99, 99)
  y <- 1

  for (j in 1:8) {
    y <- y * (abs(4 * X[, j] - 2) + a[j])/(1 + a[j])
  }
  return(y)
}

# ISHIGAMI FUNCTION
##################################################################################

ishigami <- function(X1, X2, X3) {
  A <- 2
  B <- 1
  sin(X1) + A * sin(X2) ^ 2 + B * X3 ^ 4 * sin(X1)
}

#' Ishigami function
#'
#' It implements the \insertCite{Ishigami1990;textual}{sensobol} function.
#'
#' @param X A data frame or numeric matrix where each column is a model input and each
#' row a sample point.
#'
#' @return A numeric vector with the model output.
#' @export
#'
#' @references
#' \insertAllCited{}
#'
#' @details The function requires 3 model inputs and reads as
#' \deqn{y=\sin(x_1) +a \sin(x_2) ^ 2 + b x_3 ^4 \sin(x_1)\,,}
#' where \eqn{a=2}, \eqn{b=1} and \eqn{(x_1,x_2,x_3)\sim\mathcal{U}(-\pi, +\pi)}. The
#' transformation of the distribution of the model inputs from \eqn{U(0, 1)} to
#' \eqn{U(-\pi, +\pi)}) is conducted internally.
#'
#' @examples
#' # Define settings
#' N <- 100; params <- paste("X", 1:3, sep = "")
#'
#' # Create sample matrix
#' mat <- sobol_matrices(N = N, params = params)
#'
#' # Compute Ishigami function
#' Y <- ishigami_Fun(mat)
ishigami_Fun <- function(X) {

  X <- apply(X, 2, function(x) x * (pi + pi) - pi)

  return(mapply(ishigami,
                X[, 1],
                X[, 2],
                X[, 3]))
}

# BRATLEY ET AL 1992 FUNCTION
##################################################################################

#' Bratley, Fox and Niederreiter (1992) function.
#'
#' It implements the \insertCite{Bratley1992;textual}{sensobol} function.
#'
#' @param X A data frame or numeric matrix where each column is a model input and each
#' row a sample point.
#'
#' @return A numeric vector with the model output.
#' @export
#' @references
#' \insertAllCited{}
#'
#' @details The function requires \eqn{k} model inputs and reads as:
#' \deqn{y=\sum_{i=1}^{k}(-1)^i\prod_{j=1}^{i}x_j\,,}
#' where \eqn{x_i\sim\mathcal{U}(0,1)}.
#'
#' @examples
#' # Define settings (test with k = 10)
#' N <- 100; params <- paste("X", 1:10, sep = "")
#'
#' # Create sample matrix
#' mat <- sobol_matrices(N = N, params = params)
#'
#' # Compute Bratley et al. (1992) function
#' Y <- bratley1992_Fun(mat)
bratley1992_Fun <- function(X) {
  # Preallocate
  mat <- tmp <- vector(mode = "list", length = nrow(X))
  Y <- vector(mode = "numeric", length = nrow(X))

  for (i in 1:nrow(X)) {
    mat[[i]] <- matrix(rep(X[i, ], times = ncol(X)),
                       nrow = ncol(X),
                       ncol = ncol(X),
                       byrow = TRUE)
    mat[[i]][upper.tri(mat[[i]])] <- 1
    tmp[[i]] <- matrixStats::rowProds(mat[[i]])
    Y[[i]] <- sum(tmp[[i]] * (-1) ^ (1:ncol(X)))
  }
  return(Y)
}

# BRATLEY AND FOX 1988 FUNCTION
##################################################################################

#' Bratley and Fox (1988) function
#'
#' It implements the \insertCite{Bratley1988;textual}{sensobol} function.
#'
#' @param X A data frame or numeric matrix where each column is a model input and each
#' row a sample point.
#'
#' @return A numeric vector with the model output.
#'
#' @details The function requires \eqn{k} model inputs and reads as follows:
#' \deqn{y=\prod_{i=1}^{k} |4x_i - 2 |\,,}
#' where \eqn{x_i\sim\mathcal{U}(0,1)}.
#'
#' @export
#'
#' @examples
#' # Define settings (test with k = 10)
#' N <- 100; params <- paste("X", 1:10, sep = "")
#'
#' # Create sample matrix
#' mat <- sobol_matrices(N = N, params = params)
#'
#' # Compute Bratley and Fox (1988) function
#' Y <- bratley1988_Fun(mat)
bratley1988_Fun <- function(X) {
  y <- 1

  for (j in 1:ncol(X)) {
    y <- y * (abs(4 * X[, j] - 2))}
  return(y)
}

# OAKLEY AND O'HAGAN FUNCTION
##################################################################################

#' Oakley & O'Hagan (2004) function
#'
#' It implements the \insertCite{Oakley2004;textual}{sensobol} function.
#'
#' @param X A data frame or numeric matrix where each column is a model input and each
#' row a sample point.
#'
#' @return A numeric vector with the model output.
#'
#' @export
#' @references
#' \insertAllCited{}
#'
#' @details The function requires 15 model inputs and reads as
#' \deqn{y=\mathbf{a}_1^T \bm{x} + \mathbf{a}_2 ^ T \sin(\mathbf{x}) + \mathbf{a}_3 ^ T \cos(\mathbf{x}) + \mathbf{x}^T \mathbf{M}\mathbf{x}\,,}
#' where \eqn{\mathbf{x}=x_1,x_2,...,x_k}, \eqn{k=15}, and values
#' for \eqn{\mathbf{a}^T_i,i=1,2,3} and \eqn{\mathbf{M}} are defined by \insertCite{Oakley2004;textual}{sensobol}. The
#' transformation of the distribution of the model inputs from \eqn{U(0, 1)} to
#' \eqn{N(0, 1)}) is conducted internally.
#'
#' @examples
#' # Define settings
#' N <- 100; params <- paste("X", 1:15, sep = "")
#'
#' # Create sample matrix
#' mat <- sobol_matrices(N = N, params = params)
#'
#' # Compute Oakley and O'Hagan (2004) function
#' Y <- oakley_Fun(mat)
oakley_Fun <- function(X) {

  a1 <- c(0.0118, 0.0456, 0.2297, 0.0393, 0.1177,
          0.3865, 0.3897, 0.6061, 0.6159, 0.4005,
          1.0741, 1.1474, 0.788, 1.1242, 1.1982)

  a2 <- c(0.4341, 0.0887, 0.0512, 0.3233, 0.1489,
          1.036, 0.9892, 0.9672, 0.8977, 0.8083,
          1.8426, 2.4712, 2.3946, 2.0045, 2.2621)

  a3 <- c(0.1044, 0.2057, 0.0774, 0.273, 0.1253,
          0.7526, 0.857, 1.0331, 0.8388, 0.797,
          2.2145, 2.0382, 2.4004, 2.0541, 1.9845)

  M <- c(-0.022482886, -0.18501666, 0.13418263, 0.36867264, 0.17172785, 0.13651143, -0.44034404,
         -0.081422854, 0.71321025, -0.44361072, 0.50383394, -0.024101458, -0.045939684, 0.21666181,
         0.055887417, 0.2565963, 0.053792287, 0.25800381, 0.23795905, -0.59125756, -0.081627077,
         -0.28749073, 0.41581639, 0.49752241, 0.083893165, -0.11056683, 0.033222351, -0.13979497,
         -0.031020556, -0.22318721, -0.055999811, 0.19542252, 0.095529005, -0.2862653, -0.14441303,
         0.22369356, 0.14527412, 0.28998481, 0.2310501, -0.31929879, -0.29039128, -0.20956898, 0.43139047,
         0.024429152, 0.044904409, 0.66448103, 0.43069872, 0.29924645, -0.16202441, -0.31479544,
         -0.39026802, 0.17679822, 0.057952663, 0.17230342, 0.13466011, -0.3527524, 0.25146896, -0.018810529,
         0.36482392, -0.32504618, -0.121278, 0.12463327, 0.10656519, 0.046562296, -0.21678617, 0.19492172,
         -0.065521126, 0.024404669, -0.09682886, 0.19366196, 0.33354757, 0.31295994, -0.083615456,
         -0.25342082, 0.37325717, -0.2837623, -0.32820154, -0.10496068, -0.22073452, -0.13708154,
         -0.14426375, -0.11503319, 0.22424151, -0.030395022, -0.51505615, 0.017254978, 0.038957118,
         0.36069184, 0.30902452, 0.050030193, -0.077875893, 0.003745656, 0.88685604, -0.26590028,
         -0.079325357, -0.042734919, -0.18653782, -0.35604718, -0.17497421, 0.088699956, 0.40025886,
         -0.055979693, 0.13724479, 0.21485613, -0.011265799, -0.09229473, 0.59209563, 0.031338285,
         -0.033080861, -0.24308858, -0.099798547, 0.034460195, 0.095119813, -0.3380162, 0.0063860024,
         -0.61207299, 0.081325416, 0.88683114, 0.14254905, 0.14776204, -0.13189434, 0.52878496, 0.12652391,
         0.045113625, 0.58373514, 0.37291503, 0.11395325, -0.29479222, -0.57014085, 0.46291592, -0.094050179,
         0.13959097, -0.38607402, -0.4489706, -0.14602419, 0.058107658, -0.32289338, 0.093139162,
         0.072427234, -0.56919401, 0.52554237, 0.23656926, -0.011782016, 0.071820601, 0.078277291,
         -0.13355752, 0.22722721, 0.14369455, -0.45198935, -0.55574794, 0.66145875, 0.34633299, 0.14098019,
         0.51882591, -0.28019898, -0.1603226, -0.068413337, -0.20428242, 0.069672173, 0.23112577,
         -0.044368579, -0.16455425, 0.21620977, 0.0042702105, -0.087399014, 0.31599556, -0.027551859,
         0.13434254, 0.13497371, 0.05400568, -0.17374789, 0.17525393, 0.060258929, -0.17914162, -0.31056619,
         -0.25358691, 0.025847535, -0.43006001, -0.62266361, -0.033996882, -0.29038151, 0.03410127,
         0.034903413, -0.12121764, 0.026030714, -0.33546274, -0.41424111, 0.05324838, -0.27099455,
         -0.026251302, 0.41024137, 0.26636349, 0.15582891, -0.18666254, 0.019895831, -0.24388652,
         -0.44098852, 0.012618825, 0.24945112, 0.071101888, 0.24623792, 0.17484502, 0.0085286769,
         0.2514707, -0.14659862, -0.08462515, 0.36931333, -0.29955293, 0.1104436, -0.75690139, 0.041494323,
         -0.25980564, 0.46402128, -0.36112127, -0.94980789, -0.16504063, 0.0030943325, 0.052792942,
         0.22523648, 0.38390366, 0.45562427, -0.18631744, 0.0082333995, 0.16670803, 0.16045688)

  M <- matrix(M, 15, 15, byrow = TRUE)

  Y <- vector()
  # transformation to normal distribution
  X <- apply(X, 2, function(x) stats::qnorm(x))

  for (i in 1:nrow(X)) {
    mat <- matrix(X[i, ])
    Y[i] <- a1 %*% mat + a2 %*% sin(mat) + a3 %*% cos(mat) + t(mat) %*% M %*% mat
  }
  return(Y)
}

# BECKER 2020 METAFUNCTION
##################################################################################

# Define the list of functions
function_list <- list(
  Linear = function(x) x,
  Quadratic = function(x) x ^ 2,
  Cubic = function(x) x ^ 3,
  Exponential = function(x) exp(1) ^ x / (exp(1) - 1),
  Periodic = function(x) sin(2 * pi * x) / 2,
  Discontinuous = function(x) ifelse(x > 0.5, 1, 0),
  Non.monotonic = function(x) 4 * (x - 0.5) ^ 2,
  Inverse = function(x) (10 - 1 / 1.1) ^ -1 * (x + 0.1) ^ - 1,
  No.effect = function(x) x * 0,
  Trigonometric = function(x) cos(x),
  Piecewise.large = function(x) ((-1) ^ as.integer(4 * x) *
                                   (0.125 - (x %% 0.25)) + 0.125),
  Piecewise.small = function(x) ((-1) ^ as.integer(32 * x) *
                                   (0.03125 - 2 * (x %% 0.03125)) + 0.03125) / 2,
  Oscillation = function(x) x ^ 2 - 0.2 * cos(7 * pi * x)
)

#' Random metafunction based on \insertCite{Becker2020;textual}{sensobol}'s metafunction.
#'
#' @param data A numeric matrix where each column is a model input and each row a sampling point.
#' @param k_2 Numeric value indicating the fraction of active pairwise interactions (between 0 and 1).
#' Default is \code{k_2 = 0.5}.
#' @param k_3 Numeric value indicating the fraction of active three-wise interactions
#' (between 0 and 1). Default is \code{k_2 = 0.2}.
#' @param epsilon Integer value. It fixes the seed for the random number generator.
#' The default is \code{epsilon = NULL}.
#' @return A numeric vector with the function output.
#' @export
#'
#' @importFrom Rdpack reprompt
#'
#' @details The metafunction randomly combines the following functions in a metafunction of dimension \eqn{k}:
#' * \eqn{f(x) = x ^ 3} (cubic).
#' * \eqn{f(x) = 1~\mbox{if}(x > 0.5), 0~\mbox{otherwise}} (discontinuous).
#' * \eqn{f(x) = \frac{e ^ x}{e - 1}} (exponential).
#' * \eqn{f(x) = \frac{10 - 1}{1.1} ^ {-1} (x + 0.1) ^ {-1}} (inverse).
#' * \eqn{f(x) = x} (linear)
#' * \eqn{f(x) = 0} (no effect).
#' * \eqn{f(x) = 4(x - 0.5) ^  2} (non-monotonic).
#' * \eqn{f(x) = \frac{\sin (2 \pi x)}{2}} (periodic).
#' * \eqn{f(x) = x ^ 2} (quadratic).
#' * \eqn{f(x) = \cos(x)} (trigonometric).
#'
#' It is constructed as follows:
#'
#' \deqn{y=\sum_{i=1}^{k}\alpha_i f^{u_i}(x_i) \\
#'  + \sum_{i=1}^{k_2}\beta_i f^{u_{V_{i,1}}}(x_{V_{i,1}}) f^{u_{V_{i,2}}} (x_{V_{i,2}}) \\
#' + \sum_{i=1}^{k_3}\gamma_i f^{u_{W_{i,1}}}(x_{W_{i,1}}) f^{u_{W_{i,2}}}(x_{W_{i,2}}) f^{u_{W_{i,3}}} (x_{W_{i,3}})}
#'
#' where \eqn{k} is the model dimensionality, \eqn{u} is a \eqn{k}-length vector formed by randomly
#' sampling with replacement the ten functions mentioned above, \eqn{V} and \eqn{W} are two matrices specifying the
#' number of pairwise and three-wise interactions given the model dimensionality,
#' and \eqn{\mathbf{\alpha}, \mathbf{\beta}, \mathbf{\gamma}} are three
#' vectors of length \eqn{k} generated by sampling from a mixture of two normal distributions
#' \eqn{\Psi=0.3\mathcal{N}(0, 5) + 0.7\mathcal{N}(0, 0.5)}.
#' See \insertCite{Puyj;textual}{sensobol} and \insertCite{Becker2020;textual}{sensobol} for a full
#' mathematical description of the metafunction approach.
#'
#' @references
#' \insertAllCited{}
#' @examples
#' # Define settings (number of model inputs = 86)
#' N <- 100; params <- paste("X", 1:86, sep = "")
#'
#' # Create sample matrix
#' mat <- sobol_matrices(N = N, params = params)
#'
#' # Compute metafunction
#' Y <- metafunction(mat)
metafunction <- function(data, k_2 = 0.5, k_3 = 0.2, epsilon = NULL) {
  k <- ncol(data)
  set.seed(epsilon)
  all_functions <- sample(names(function_list), k, replace = TRUE)
  set.seed(epsilon)
  components <- sample(1:2, prob = c(0.7, 0.3), size = 200,
                       replace = TRUE)
  mus <- c(0, 0)
  sds <- sqrt(c(0.5, 5))
  set.seed(epsilon)
  coefficients <- stats::rnorm(100) * sds[components] + mus[components]
  set.seed(epsilon)
  coefD1 <- sample(coefficients, k)
  d2 <- t(utils::combn(1:k, 2))
  set.seed(epsilon)
  d2M <- d2[sample(nrow(d2), size = ceiling(k * k_2), replace = FALSE),
  ]
  sample.size.d2M <- ifelse(is.vector(d2M) == TRUE, 1, nrow(d2M))
  set.seed(epsilon)
  coefD2 <- sample(coefficients, sample.size.d2M, replace = TRUE)
  d3 <- t(utils::combn(1:k, 3))
  set.seed(epsilon)
  size.d3 <- ifelse(nrow(d3) == 1, 1, ceiling(k * k_3))
  set.seed(epsilon)
  d3M <- d3[sample(nrow(d3), size = size.d3, replace = FALSE),
  ]
  sample.size <- ifelse(is.vector(d3M) == TRUE, 1, nrow(d3M))
  set.seed(epsilon)
  coefD3 <- sample(coefficients, sample.size, replace = TRUE)
  output <- sapply(seq_along(all_functions), function(x) function_list[[all_functions[x]]](data[,
                                                                                                x]))
  y1 <- Rfast::rowsums(mmult(output, coefD1))
  if (is.vector(d2M) == TRUE) {
    y2 <- sum(output[, d2M[1]] * output[, d2M[2]] * coefD2)
  }
  else {
    y2 <- Rfast::rowsums(mmult(output[, d2M[, 1]] * output[,
                                                           d2M[, 2]], coefD2))
  }
  if (is.vector(d3M) == TRUE) {
    y3 <- sum(output[, d3M[1]] * output[, d3M[2]] * output[,
                                                           d3M[3]] * coefD3)
  }
  else {
    y3 <- Rfast::rowsums(mmult(output[, d3M[, 1]] * output[,
                                                           d3M[, 2]] * output[, d3M[, 3]], coefD3))
  }
  Y <- y1 + y2 + y3
  return(Y)
}


