
# PLOT SOBOL' FIRST AND TOTAL-ORDER INDICES -----------------------------------

#' Plot Sobol' indices.
#'
#' @param data The output of \code{sobol_indices}.
#' @param order If \code{order = "first"}, it plots first and total effects.
#' If \code{order = "second"}, it plots second-order effects. If \code{order = "third"}, it plots
#' third-order effects. Default is \code{order = "first"}.
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
    p <- data.table::data.table(data)[sensitivity == "Si" | sensitivity == "Ti"]
    gg <- ggplot2::ggplot(p, aes(parameters, original, fill = sensitivity)) +
      geom_bar(stat = "identity",
               position = position_dodge(0.6),
               color = "black") +
      geom_errorbar(aes(ymin = low.ci,
                        ymax = high.ci),
                    position = position_dodge(0.6)) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
      labs(x = "",
           y = "Sobol' index") +
      scale_fill_discrete(name = "Sobol' indices",
                          labels = c(expression(S[italic(i)]),
                                     expression(T[italic(i)]))) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.background = element_rect(fill = "transparent",
                                             color = NA),
            legend.key = element_rect(fill = "transparent",
                                      color = NA),
            legend.position = "top",
            strip.text.y = element_text(size = 6))
  } else if(!order == "first") {
    if(order == "second") {
      plot.type <- "Sij"
    } else if(order == "third") {
      plot.type <- "Sijk"
    } else {
      stop("Order should be either first, second or third")
    }
    p <- data.table::data.table(data)[sensitivity == plot.type]
    gg <- ggplot2::ggplot(p, aes(stats::reorder(parameters, original),
                                 original)) +
      geom_point() +
      geom_errorbar(aes(ymin = low.ci,
                        ymax = high.ci)) +
      theme_bw() +
      labs(x = "",
           y = "Sobol' index") +
      geom_hline(yintercept = 0,
                 lty = 2,
                 color = "red") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.background = element_rect(fill = "transparent",
                                             color = NA),
            legend.key = element_rect(fill = "transparent",
                                      color = NA),
            axis.text.x = element_text(angle = 45,
                                       hjust = 1))
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
  gg <- ggplot2::ggplot(df, aes(Y)) +
    geom_histogram(color = "black",
                   fill = "white") +
    labs(x = "Y",
         y = "Count") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.background = element_rect(fill = "transparent",
                                           color = NA),
          legend.key = element_rect(fill = "transparent",
                                    color = NA))
  return(gg)
}


