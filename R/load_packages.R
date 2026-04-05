
# FUNCTION TO INSTALL AND LOAD REQUIRED PACKAGES
##################################################################################

#' Load R packages.
#'
#' The function loads R packages. If a package is not installed,
#' the function stops with an informative error message.
#'
#' @param x A character vector with the name of the packages to load.
#'
#' @return NULL
#' @export
#'
#' @examples
#' # Load packages:
#' \dontrun{load_packages(c("tidyverse", "data.table"))}
load_packages <- function(x) {

  for (i in x) {

    if (!requireNamespace(i, quietly = TRUE)) {
      stop(paste0("Package '", i, "' is required but not installed. ",
                  "Install it with: install.packages('", i, "')"))
    }
    library(i, character.only = TRUE)
  }
}
