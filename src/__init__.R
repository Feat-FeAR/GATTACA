# This script is ran before every local script, doing housekeeping tasks in the
# docker. The name is borrowed from Python.

options(
    repos = list(
      CRAN = c("https://cloud.r-project.org/"),
      BioCsoft = c("https://bioconductor.org/packages/3.14/bioc"),
      BioCann = c("https://bioconductor.org/packages/3.14/data/annotation"),
      BioCexp = c("https://bioconductor.org/packages/3.14/data/experiment")
    )
)

# This NEEDS to be futile.logger THEN logger, as logger needs to overwrite
# stuff from futile.logger.
suppressMessages(library(futile.logger))
suppressMessages(library(logger))

log_layout(layout_glue_generator(format = '{time} <{fn}> [{level}]: {msg}'))

if (!interactive()) {
  # This is wrapped in an interactive call to not run it when sourcing this
  # while working in rstudio and the like.
  ROOT <- "/GATTACA" # The root in the docker
} else {
  ROOT <- getwd()
  if (! getOption("running_tests", FALSE) == TRUE) {
    log_threshold(TRACE)
  }
}

setup_file_logging <- function (log_dir, log_name = NULL) {
  # Setup logging facilities to save to file.
  if (log_name == "NULL" | log_name == "") { log_name <- NULL }

  if (is.null(log_name)) {
    start.timedate <- gsub(" ", "_", date())
    log_name <- paste0("GATTACA_", start.timedate, ".log")
  }

  log_path <- file.path(log_dir, log_name)
  file.create(log_path)
  log_appender(appender_tee(log_path))

  # I save the set log path for later
  options(gattaca.log.path = log_path)

  # Setup the data log
  data_log_path <- gsub(pattern = "\\.log", replacement = ".data.log", log_path)
  if (data_log_path == log_path) {
    # The user did not end the log with `.log`. WE FORCE IT - MUHAHAHA
    data_log_path <- paste0(data_log_path, ".data")
  }
  file.create(data_log_path)
  options(gattaca.datalog.path = data_log_path)

  # Log the version of the docker to the file log. For posterity.
  log_info(paste("Version of current GATTACA container:", readLines(file.path(ROOT, "VERSION"))))
}

setwd(ROOT)

# Load all packages used by the tool. The order is important as R allows to
# overwrite functions (!!).

# These two are backbone packages.
suppressMessages(library(tidyverse))
suppressMessages(library(progress))

#' Gracefully load a series of packages, as to not spam the console.
#'
#' @param packages A vector of strings with the package names to load.
#'
#' @author MrHedmad
graceful_load <- function(packages) {
  log_info("Loading required packages...")
  log_debug("Loading packages:", paste(packages, collapse = ", "))

  pb <- progress_bar$new(
    format = "Loading... [:bar] :percent (:eta)",
    total = length(packages), clear = FALSE, width= 80)
  pb$tick(0)
  for (i in seq_along(packages)) {
    package <- packages[i]
    log_debug("Loading package: ", package)

    suppressMessages(library(package, character.only = TRUE))
    pb$tick()
  }
  log_debug("Finished loading packages")
}

graceful_load(c(
  "preprocessCore",   # Interarray Normalization by Quantile-Quantile Algorithm
  "PCAtools",         # Principal Component Analysis
  "limma",            # Empirical Bayes Method for Differential Expression
  "RankProd",         # Rank Product Method for Differential Expression
  "VennDiagram",      # Venn Diagrams
  "EnhancedVolcano",  # Volcano Plots
  "gplots",           # Heatmap with extensions - heatmap.2()
  "yaml",             # Yaml file parsing
  "UpSetR",           # UpSet Plots
  "oligo",            # Prepare agilent things
  "affycoretools",    # Filter out the affy control probes
  "reshape2",         # Reshaping functions
  "AnnotationDbi",    # Base for the annotations
  "org.Hs.eg.db"      # Base reader for all SQL packages
))

source(file.path(ROOT, "src", "tools", "tools.R"))

# Refactor conflicting functions to use them with the native pipe operator.
# The fact that I need to do this is idiotic. But, alas, this is R.
dp_select <- dplyr::select
