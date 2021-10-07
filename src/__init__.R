# This script is ran before every local script, doing housekeeping tasks in the
# docker.

options(
    repos = c("http://cran.us.r-project.org")
)

# This NEEDS to be futile.logger THEN logger, as logger needs to overwrite
# stuff from futile.logger.
suppressMessages(library(futile.logger))
suppressMessages(library(logger))

setup_file_logging <- function (log_dir, log_name = NULL) {
  # Setup logging facilities to save to file.
  if (!is.null(log_name) && log_name == "NULL") {
    log_name <- NULL
  }
  
  if (is.null(log_name)) {
    start.timedate <- gsub(" ", "_", date())
    log_name <- paste0("GATTACA_", start.timedate, ".log")
  }

  log_path <- file.path(log_dir, log_name)
  file.create(log_path)
  log_appender(appender_tee(log_path))
}

if (!interactive()) {
  # This is wrapped in an interactive call to not run it when sourcing this
  # while working in rstudio and the like.
  ROOT <- "/GATTACA" # The root in the docker
} else {
  ROOT <- getwd()
  log_threshold(TRACE)
}

setwd(ROOT)

suppressMessages(library(tidyverse))
