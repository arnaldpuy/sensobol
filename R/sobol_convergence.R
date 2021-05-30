
# NOW CHECK CONVERGENCE -------------------

# Create function to swiftly check convergence
sobol_sample <- function(matrices, Y, N, sub.sample, params, first, total,
                         order = order, seed = seed, ...) {

  dt <- data.table::data.table(matrix(Y, nrow = N))
  set.seed(seed)
  sample.Y <- unlist(dt[sample(.N, sub.sample)], use.names = FALSE)

  ind <- sobol_indices(matrices = matrices, Y = sample.Y, N = sub.sample, params = params,
                       first = first, total = total, order = order, ...)
  return(ind)
}

# Create function to plot convergence
plot_sobol_convergence <- function(dt) {

  Cost <- original <- sensitivity <- NULL

  gg <- ggplot2::ggplot(dt, ggplot2::aes(Cost, original, color = sensitivity)) +
    ggplot2::geom_line() +
    ggplot2::labs(x = "Total number of model runs", y = "Sobol' indices") +
    ggplot2::facet_wrap(~parameters) +
    theme_AP() +
    ggplot2::theme(legend.position = "top")
  return(gg)
}

add_ribbon <- function() {

  low.ci <- high.ci <- sensitivity <- NULL

  ggplot2::geom_ribbon(ggplot2::aes(ymin = low.ci, ymax = high.ci, fill = sensitivity),
                       alpha = 0.1, linetype = 0)
}



#' Check convergence of Sobol' indices.
#'
#' It checks the convergence of Sobol' indices on different sub-samples of the model output-.
#'
#' @param matrices Character vector with the required matrices. The default is \code{matrices = c("A", "B", "AB")}.
#' See \code{\link{sobol_matrices}}.
#' @param Y Numeric vector with the model output obtained from the matrix created with
#' \code{\link{sobol_matrices}}.
#' @param N Positive integer, the initial sample size of the base sample matrix created with \code{\link{sobol_matrices}}.
#' @param sub.sample Numeric vector with the sub-samples of the model output at which to check convergence.
#' @param params Character vector with the name of the model inputs.
#' @param first Estimator to compute first-order indices. Check options in \code{\link{sobol_indices}}.
#' @param total Estimator to compute total-order indices. Check options in \code{\link{sobol_indices}}.
#' @param order Whether to plot convergence for "second" or "third" order indices.
#' @param seed Whether to compute "first", "second", or "third" -order Sobol' indices. Default
#' is \code{order = "first"}.
#' @param plot.order Whether to plot convergence for "second" or "third"-order indices.
#' @param ... Further arguments in \code{\link{sobol_indices}}.
#'
#' @return A list with the results and the plots
#' @export
#'
#' @examples
#' # Define settings
#' matrices <- c("A", "B", "AB")
#' params <- paste("X", 1:3, sep = "")
#' N <- 2^10
#' first <- "saltelli"
#' total <- "jansen"
#' order <- "second"
#'
#' # Create sample matrix
#' mat <- sobol_matrices(N = N, params = params, order = order)
#'
#' # Compute Ishigami function
#' Y <- ishigami_Fun(mat)
#'
#' # Check convergence at specific sample sizes
#' sub.sample <- seq(100, N, 500) # Define sub-samples
#' sobol_convergence(matrices = matrices, Y = Y, N = N, sub.sample = sub.sample,
#' params = params, first = first, total = total, order = order, plot.order = order)

sobol_convergence <- function(matrices, Y, N, sub.sample, params, first, total,
                              order = order, seed = 666, plot.order, ...) {

  sensitivity <- NULL

  if (sub.sample[length(sub.sample)] > N) {
    stop("The last sub-sample should be equal or smaller than N")
  }

  Cost <- indices <- list()
  ind.list <- lapply(sub.sample, function(n)
    sobol_sample(matrices = matrices, Y = Y, N = N, sub.sample = n,
                 params = params, first = first, total = total, order = order, seed = seed,...))

  for (i in 1:length(ind.list)) {

    Cost[[i]] <- ind.list[[i]]$C
    indices[[i]] <- ind.list[[i]]$results
  }
  names(indices) <- Cost

  final.dt <- data.table::rbindlist(indices, idcol = "Cost")[, Cost:= as.numeric(Cost)]
  check.boot <- any(grepl("high.ci", colnames(final.dt)))
  first.order <- final.dt[sensitivity %in% c("Si", "Ti")]
  gg <- plot_sobol_convergence(first.order)

  if (check.boot == TRUE) {
    gg <- gg + add_ribbon()

  }

  output <- list(final.dt, gg)
  names(output) <- c("convergence", "plot.first")

  if (order == "second" & plot.order == "second" |
      order == "third" & plot.order == "second" |
      order == "third" & plot.order == "third") {
    second.order <- final.dt[sensitivity == "Sij"]
    gg2 <- plot_sobol_convergence(second.order) +
      ggplot2::scale_color_manual(values = "#7CAE00", name = "Sensitivity") +
      ggplot2::scale_fill_manual(values = "#7CAE00", name = "Sensitivity")

    if (check.boot == TRUE) {
      gg2 <- gg2 + add_ribbon()

    }

    output <- list(final.dt, gg, gg2)
    names(output) <- c("convergence", "plot.first", "plot.second")

  }

  if (order == "third" & plot.order == "third") {
    third.order <- final.dt[sensitivity == "Sijk"]
    gg3 <- plot_sobol_convergence(third.order) +
      ggplot2::scale_color_manual(values = "#C77CFF", name = "Sensitivity") +
      ggplot2::scale_fill_manual(values = "#C77CFF", name = "Sensitivity")

    if (check.boot == TRUE) {
      gg3 <- gg3 + add_ribbon()

    }

    output <- list(final.dt, gg, gg2, gg3)
    names(output) <- c("convergence", "plot.first", "plot.second", "plot.third")

  } else if (order == "first" & plot.order == "second" |
             order == "first" & plot.order == "third" |
             order == "second" & plot.order == "third") {

    stop("you need to compute high-order indices to plot convergence
         for high-order indices")
  }


  return(output)
}
