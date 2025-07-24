
# ODE SOLVER FOR SENSOBOL
##################################################################################

#' Wrapper around \code{deSolve} \code{\link[deSolve]{ode}}.
#'
#' It solves a system of ordinary differential equations and extracts the model output
#' at the selected times.
#'
#' @param d Character vector with the name of the model inputs.
#' @param times Time sequence as defined by \code{\link[deSolve]{ode}}.
#' @param timeOutput Numeric vector determining the time steps at which
#' the output is wanted.
#' @param state Initial values of the state variables.
#' @param func An R function as defined by \code{\link[deSolve]{ode}}.
#' @param ... Additional arguments passed to \code{\link[deSolve]{ode}}.
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
#' # Define the times
#'  times <- seq(5, 20, 1)
#'
#'  # Define the times at which the output is wanted
#'  timeOutput <- c(10, 15)
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
#' for (i in 1:nrow(mat)) {
#'   y[[i]] <- sobol_ode(d = mat[i, ],
#'                      times = times,
#'                      timeOutput = timeOutput,
#'                      state = c(X = mat[[i, "X"]], Y = mat[[i, "Y"]]),
#'                      func = lotka_volterra_fun)
#'}

sobol_ode <- function(d, times, timeOutput, state, func,...) {

  out <- deSolve::ode(y = state, times = times, func = func, parms = d,...)

  check.timeOutput <- all(timeOutput %in% out[, "time"] )

  if(isFALSE(check.timeOutput) == TRUE) {
    stop("all values in timeOutput should be time steps in time")

  } else {
    output <- out[out[, "time"] %in% timeOutput, ]
  }

  return(output)
}

