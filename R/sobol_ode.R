


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

