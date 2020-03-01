
# PLOT SOBOL' FIRST AND TOTAL-ORDER INDICES -----------------------------------

#' Plot Sobol' first and total-order indices
#'
#' @param data The output of \code{sobol_indices}.
#' @param dummy The output of the \code{sobol_dummy} function. If supplied and
#' \code{type = 1}, the plot includes two horizontal dashed lines showing the
#' highest confidence interval for the first and total-order effects of the
#' dummy parameter.
#' @param type An integer. If \code{type = 1}, it plots first and total effects.
#' If \code{type = 2}, it plots second-order effects. If \code{type = 3}, it plots
#' third-order effects. Default is \code{type = 1}.
#'
#' @return A ggplot object.
#' @import ggplot2
#' @export
#'
#' @examples
#' # Define settings
#' N <- 100; params <- paste("X", 1:3, sep = ""); R <- 10
#' # Create sample matrix
#' mat <- sobol_matrices(N = N, params = params)
#' # Compute Ishigami function
#' Y <- ishigami_Fun(mat)
#' # Compute and bootstrap Sobol' indices
#' ind <- sobol_indices(Y = Y, N = N, params = params, boot = TRUE, R = R)
#' # Plot Sobol' indices
#' plot_sobol(data = ind)

plot_sobol <- function(data, dummy = NULL, type = 1) {
  sensitivity <- low.ci <- high.ci <- parameters <- original <- NULL
  if(type == 1) {
    if(is.null(dummy) == FALSE) {
      plot.dummy <- geom_rect(data = dummy,
                              aes(ymin = 0,
                                  ymax = high.ci,
                                  xmin = -Inf,
                                  xmax = Inf,
                                  fill = sensitivity),
                              alpha = 0.2,
                              inherit.aes = FALSE)
    } else {
      plot.dummy <- NULL
    }
    p <- data[sensitivity == "Si" | sensitivity == "STi"]
    gg <- ggplot2::ggplot(p, aes(parameters, original,
                                 fill = sensitivity)) +
      geom_bar(stat = "identity",
               position = position_dodge(0.6),
               color = "black") +
      plot.dummy +
      geom_errorbar(aes(ymin = low.ci,
                        ymax = high.ci),
                    position = position_dodge(0.6)) +
      scale_fill_discrete(name = "Sobol' indices",
                          labels = c(expression(S[italic(i)]),
                                     expression(S[italic(T[i])]))) +
      labs(x = "",
           y = "Variance") +
      theme_bw() +
      theme(legend.position = "top",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.background = element_rect(fill = "transparent",
                                             color = NA),
            legend.key = element_rect(fill = "transparent",
                                      color = NA))
  } else if(!type == 1) {
    if(type == 2) {
      plot.type <- "Sij"
    } else if(type == 3) {
      plot.type <- "Sijk"
    } else {
      stop("Type should be either 1, 2 or 3")
    }
    p <- data[sensitivity == plot.type]
    gg <- ggplot2::ggplot(p, aes(stats::reorder(parameters, original),
                                 original)) +
      geom_point() +
      geom_errorbar(aes(ymin = low.ci,
                        ymax = high.ci)) +
      theme_bw() +
      labs(x = "",
           y = "Variance") +
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
#' @param n An integer with the initial sample size.
#'
#' @return a ggplot2 object.
#' @import ggplot2
#' @export
#'
#' @examples
#' # Define settings:

plot_uncertainty <- function(Y, n = NULL) {
  warning("The argument n will be substituted by N in the next release of the package")
  if(is.vector(Y) == FALSE) {
    stop("Y should be a vector")
  }
  if(is.null(n) == TRUE) {
    stop("n should be specified")
  }
  Y <- Y[1:(2*n)]
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


