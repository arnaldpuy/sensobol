
theme_AP <- function() {
  ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          legend.background = ggplot2::element_rect(fill = "transparent",
                                           color = NA),
          legend.key = ggplot2::element_rect(fill = "transparent", color = NA),
          strip.background = ggplot2::element_rect(fill = "white"))
}

# PLOT SOBOL' FIRST AND TOTAL-ORDER INDICES -----------------------------------

#' Plot Sobol' indices.
#'
#' @param data The output of \code{sobol_indices}.
#' @param order If \code{order = "first"}, it plots first and total effects.
#' If \code{order = "second"}, it plots second-order effects. If \code{order = "third"}, it plots
#' third-order effects. Default is \code{order = "first"}
#'
#' @return A ggplot object.
#' @import ggplot2
#' @export
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
#'
#' # Plot Sobol' indices
#' plot_sobol(data = ind)
plot_sobol <- function(data, order = "first") {
  sensitivity <- low.ci <- high.ci <- parameters <- original <- NULL
  if(order == "first") {
    p <- data[sensitivity == "Si" | sensitivity == "Ti"]
    gg <- ggplot2::ggplot(p, ggplot2::aes(parameters, original, fill = sensitivity)) +
      ggplot2::geom_bar(stat = "identity",
               position = ggplot2::position_dodge(0.6),
               color = "black") +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = low.ci,
                        ymax = high.ci),
                    position = ggplot2::position_dodge(0.6)) +
      ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
      ggplot2::labs(x = "",
           y = "Sobol' index") +
      ggplot2::scale_fill_discrete(name = "Sobol' indices",
                          labels = c(expression(S[italic(i)]),
                                     expression(T[italic(i)]))) +
      theme_AP()
  } else if(!order == "first") {
    if(order == "second") {
      plot.type <- "Sij"
    } else if(order == "third") {
      plot.type <- "Sijk"
    } else {
      stop("Order should be either first, second or third")
    }
    p <- data[sensitivity == plot.type]
    gg <- ggplot2::ggplot(p, ggplot2::aes(stats::reorder(parameters, original),
                                 original)) +
      ggplot2::geom_point() +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = low.ci,
                        ymax = high.ci)) +
      ggplot2::labs(x = "",
           y = "Sobol' index") +
      ggplot2::geom_hline(yintercept = 0,
                 lty = 2,
                 color = "red") +
      theme_AP()
  }
  return(gg)
}

# PLOT MODEL OUTPUT UNCERTAINTY -----------------------------------------------

#' Plot model output uncertainty
#'
#' It creates an histogram with the model output distribution.
#'
#' @param Y A numeric vector with the model output.
#' @param N An integer with the initial sample size of the base matrix, defined in \code{\link{sobol_matrices}}.
#'
#' @return a ggplot2 object.
#' @import ggplot2
#' @export
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
#' # Plot uncertainty
#' plot_uncertainty(Y = Y, N = N)

plot_uncertainty <- function(Y, N = NULL) {
  if(is.vector(Y) == FALSE) {
    stop("Y should be a vector")
  }
  if(is.null(N) == TRUE) {
    stop("The size of the base sample matrix N should be specified")
  }
  Y <- Y[1:(2 * N)]
  df <- data.frame(Y)
  gg <- ggplot2::ggplot(df, ggplot2::aes(Y)) +
    ggplot2::geom_histogram(color = "black",
                   fill = "white") +
    ggplot2::labs(x = "y",
         y = "Count") +
    theme_AP()
  return(gg)
}


