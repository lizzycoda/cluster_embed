install_if_missing <- function(pkg, version = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!is.null(version)) {
      if (!requireNamespace("remotes", quietly = TRUE)) {
        install.packages("remotes")
      }
      remotes::install_version(pkg, version = version, repos = "http://cran.us.r-project.org")
    } else {
      install.packages(pkg)
    }
  }
}

packages <- list(
  ggplot2   = "3.5.1",
  readr     = "2.5.1",
  dplyr     = "1.1.4",
  tidyr     = "1.3.1",
  Rtsne     = "0.17",
  MASS      = "7.3-61",
  vegan     = "2.6-10",
  FreeSortR = "1.3",
  igraph    = "2.1.3",
  dbscan    = "1.2.2",
  proxy     = "0.4-27",
  torch     = "0.14.2"
)


for (pkg in names(packages)) {
  install_if_missing(pkg, packages[[pkg]])
}

message("All packages installed (or already available).")
