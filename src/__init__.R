# This script is ran before every local script, doing housekeeping tasks in the
# docker.

options(
    repos = c("http://cran.us.r-project.org")
)

if (!interactive()) {
  # This is wrapped in an interactive call to not run it when sourcing this
  # while working in rstudio and the like.
  ROOT <- "/GATTACA" # The root in the docker
} else {
  ROOT <- getwd()
}

setwd(ROOT)

suppressMessages(library(tidyverse))
