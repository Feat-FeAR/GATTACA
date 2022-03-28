# ------------------------------------------------------------------------------
#
# This file runs all tests in the package. It is a substitute to the canonical
# functions that run tests (like `devtools::test()`)
#
# ------------------------------------------------------------------------------

log_threshold(DEBUG)

suppressMessages(library(testthat))

# This tells `__init__` to not reset the log level.
options(running_tests = TRUE)

#' Load test data from the resources folder
#'
#' The way to use this is to first make the output you want, then use
#' `save()` to save it to a file, and then use this function to reload it.
get_test_data <- function(name) {
  name <- paste0(name, ".testdata")
  load(file = file.path(ROOT, "src", "resources", "tests", name), envir = globalenv())
}


#' Remove all items in the global environment
#'
#' @param except A vector of strings with item names *not* to remove.
#'
#' @author MrHedmad
clear_global_env <- function(except = NULL) {
  all_things <- ls(all = TRUE, envir = globalenv())
  all_things <- all_things[!all_things %in% except]
  rm(list = all_things, envir = globalenv())
}

all_test_files <- list.files(file.path(ROOT, "tests"))
all_test_files <- all_test_files[startsWith(all_test_files, "test_")]
# Prevent infinite loops
all_test_files <- all_test_files[! all_test_files == "test_all.R"]

#' Prepare the stage for a test.
#'
#' This unloads everything, removes all items in the env and returns to the ROOT
#' so that the new test can be run cleanly.
#'
#' This is a bad way to do this, we could do stuff like `on.exit` etc but
#' I learned about their existance after this, so I can't be bothered to change.
#'
#' @author MrHedmad
prepare_env <- function() {
  # The `file` exception is for the loop.
  clear_global_env(except = c(
    "prepare_env", "ROOT", "all_test_files", "get_test_data", "file",
    "clear_global_env"
  ))
  setwd(ROOT)
}


log$info("Starting tests... Suppressing logging...")
log_threshold(ERROR)
prepare_env()
for (file in all_test_files) {
  source(file.path(ROOT, "src", "__init__.R"))
  testthat::test_file(file.path(ROOT, "tests", file))
  prepare_env()
}
# Reset the envir back
clear_global_env(except = "ROOT")
source(file.path(ROOT, "src", "__init__.R"))

log_threshold(DEBUG)

log$info("Finished running tests.")

warnings()
