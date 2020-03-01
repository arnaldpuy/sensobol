

#' ODE solver in a Monte Carlo setting
#'
#'
#'
#' @param data A sample matrix created with the \code{\link{sobol_matrices}} function.
#' @param times Time points at which we desire the model output.
#' @param dt Time increase.
#' @param state Value of the state variables.
#' @param func Function as required by the \code{\link[deSolve]{ode}} function.
#'
#' @importFrom foreach "%dopar%"
#' @importFrom foreach "%:%"
#' @return A matrix
#'
#' @export
#'
#' @examples
#' # Define settings
#'  N <- 2 ^ 3 # Sample size of sample matrix
#' params <- c("r_b", "K", "beta", "alpha", "r_s", "K_s", "K_e", "r_e", "P", "T_e")
#'
#' # Create sample matrix
#' dt <- sobol_matrices(N = N, params = params)
#'
#' # Transform to appropriate distributions
#' dt[, "r_b"] <- qunif(dt[, "r_b"], 1.52, 1.6)
#' dt[, "K"] <- qunif(dt[, "K"], 100, 355)
#' dt[, "beta"] <- qunif(dt[, "beta"], 20000, 43200)
#' dt[, "alpha"] <- qunif(dt[, "alpha"], 1, 2)
#' dt[, "r_s"] <- qunif(dt[, "r_s"], 0.095, 0.15)
#' dt[, "K_s"] <- qunif(dt[, "K_s"], 24000, 25440)
#' dt[, "K_e"] <- qunif(dt[, "K_e"], 1, 1.2)
#' dt[, "r_e"] <- qunif(dt[, "r_e"], 0.92, 1)
#' dt[, "P"] <- qunif(dt[, "P"], 0.0015, 0.00195)
#' dt[, "T_e"] <- qunif(dt[, "T_e"], 0.7, 0.9)

# Define timesteps
#' times <- c(1, seq(25, 50, 25))
#'
#' # Run function
#' \donttest{sobol_ode(data = dt, dt = 1, state = state, func = budworm_diff)}

sobol_ode <- function(data, times, dt, state, func) {
  ode_fun <- function(d, times, state, func) {
    out <- deSolve::ode(y = state, times = times, func = func, parms = d)
    return(out[times[length(times)], names(state)])
  }
  i <- j <- NULL
  n.cores <- parallel::makeCluster(floor(parallel::detectCores() * 0.75))
  doParallel::registerDoParallel(n.cores)
  out <- foreach::foreach(j = times,
                          .combine = "rbind") %:%
    foreach::foreach(i=1:nrow(data),
                     .combine = "rbind",
                     .packages = "deSolve") %dopar%
    {
      ode_fun(d = data[i, ],
              times = seq(0, j, dt),
              state = state,
              func = func)
    }
  parallel::stopCluster(n.cores)
  out <- cbind(out, rep(times, each = nrow(data)))
  return(out)
}

