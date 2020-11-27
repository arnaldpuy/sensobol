


#' Solver for ordinary differential equations
#'
#' @param params Vector with the name of the model inputs.
#' @param times Vector with the time sequences at which the model output is wanted.
#' @param state Initial values of the state variables.
#' @param func An R function as defined by \code{ode}.
#' @param ... Additional arguments passed to \code{ode}.
#'
#' @import foreach
#' @importFrom foreach %do%
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
#' R <- 100 # Number of bootstrap replicas
#'
#' # Construct the sample matrix
#' mat <- sobol_matrices(N = N, params = params)
#'
#' # Define the timesteps
#' times <- seq(5, 15, 5)
#'
#' # Run the model
sobol_ode <- function(params, times, state, func,...) {
  out <- deSolve::ode(y = state, times = times, func = func, parms = params,...)
  return(out[times[length(times)], names(state)])
}

