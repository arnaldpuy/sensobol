
# ODE SOLVER FOR SENSOBOL
##################################################################################

#' Wrapper around \code{deSolve} \code{\link{ode}}.
#'
#' It solves a system of ordinary differential equations and extracts the model output
#' at the selected times.
#'
#' @param d Character vector with the name of the model inputs.
#' @param times Numeric vector with the time sequences at which the model output is wanted.
#' @param state Initial values of the state variables.
#' @param func An R function as defined by \code{\link{ode}}.
#' @param ... Additional arguments passed to \code{\link{ode}}.
#'
#' @return A matrix with the output values.
#' @export
#'
#' @examples
#' # Define the model: the Lotka-Volterra system of equations
#' lotka_volterra_fun <- function(t, state, parameters) {
#' with(as.list(c(state, parameters)), {
#'   dX <- r * X * (1 - X / K) - alpha * X * Y
#'   dY <- -m * Y + theta * X * Y
#'   list(c(dX, dY))
#'  })
#'  }
#'
#' # Define the settings of the sensitivity analysis
#' N <- 2 ^ 5 # Sample size of sample matrix
#' params <- c("r", "alpha", "m", "theta", "K", "X", "Y") # Parameters
#'
#' # Define the timesteps
#' times <- seq(5, 15, 5)
#'
#' # Construct the sample matrix
#' mat <- sobol_matrices(N = N, params = params)
#'
#' # Transform to appropriate distributions
#' mat[, "r"] <- qunif(mat[, "r"], 0.8, 1.8)
#' mat[, "alpha"] <- qunif(mat[, "alpha"], 0.2, 1)
#' mat[, "m"] <- qunif(mat[, "m"], 0.6, 1)
#' mat[, "theta"] <- qunif(mat[, "theta"], 0.05, 0.15)
#' mat[, "K"] <- qunif(mat[, "K"], 47, 53)
#' mat[, "X"] <- floor(mat[, "X"] * (15 - 8 + 1) + 8)
#' mat[, "Y"] <- floor(mat[, "Y"] * (2 - 6 + 1) + 6)
#'
#' # Run the model
#' y <- list()
#' for (j in 1:length(times)) {
#'  for (i in 1:nrow(mat)) {
#'    y[[j]] <- sobol_ode(d = mat[i, ],
#'                        times = seq(0, j, 1),
#'                        state = c(X = mat[[i, "X"]], Y = mat[[i, "Y"]]),
#'                        func = lotka_volterra_fun)
#'  }
#'}

sobol_ode <- function(d, times, state, func,...) {
  out <- deSolve::ode(y = state, times = times, func = func, parms = d,...)
  return(out[times[length(times)], names(state)])
}

