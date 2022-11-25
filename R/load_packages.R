
# FUNCTION TO INSTALL AND LOAD REQUIRED PACKAGES
##################################################################################

#' Load (and install) R packages.
#'
#' The function loads R packages. If the packages are not already
#' in the local system, the function also downloads, installs and loads them.
#'
#' @param x A character vector with the name of the packages to load.
#'
#' @importFrom utils install.packages
#'
#' @return NULL
#' @export
#'
#' @examples
#' # Load packages:
#' \dontrun{load_packages(c("tidyverse", "data.table"))}
load_packages <- function(x) {

  for (i in x) {

    if (!require(i, character.only = TRUE)) {

      install.packages(i, dependencies = TRUE)

      library(i, character.only = TRUE)
    }
  }
}
