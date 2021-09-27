# This script loads the virtual env (renv) so it has to be sourced at the
# beginning of every docker call. It also sets the working dir to here, so we
# can give relative paths in relative safety.

# It also sets all the variables that are needed in the other scripts.
# This assumes that the renv is already on and the wd is the root folder.

# This will assure that the `renv` env is active. For instance if this
# script was started in `--vanilla` mode (as it should).
if (!interactive()) {
  # This is wrapped in an interactive call to not run it when sourcing this
  # while working in rstudio and the like.
  ROOT <- "/GATTACA" # The root in the docker
  source(file.path(ROOT, "renv", "activate.R")) 
} else {
  ROOT <- getwd()
}

setwd(ROOT)

suppressMessages(library(tidyverse))
