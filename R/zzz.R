.onAttach <- function(libname, pkgname) {
  packageStartupMessage("MCPA!")
}

.onLoad <- function(libname, pkgname) {
  library(ggplot2)
  library(progress)
  library(PropCIs)
  library(dplyr)
  library(broom)
  op <- options()
  op.devtools <- list(
    devtools.path = "~/R-dev",
    devtools.install.args = "",
    devtools.name = "mcpa",
    devtools.desc.author = 'person("Julian", "Jara-Ettinger",
      "julian.jara-ettinger@yale.edu", role = c("auth", "cre"))',
    devtools.desc.license = "MIT",
    devtools.desc.suggests = NULL,
    devtools.desc = list()
  )
  toset <- !(names(op.devtools) %in% names(op))
  if(any(toset)) options(op.devtools[toset])
  invisible()
}
