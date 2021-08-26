# Header Info ------------------------------------------------------------------
#
# Affymetrix CEL file -> Expression matrix
#
# a FeAR R-script
# based on: "Homer MiniTutorial - Processing Affymetrix Gene Expression Arrays"
# see pdf for more info
#
# General Script for CEL-file normalization through RMA algorithm
#
#   - CEL file loading
#   - RMA normalization
#       1 - Background correcting
#       2 - Normalizing
#       3 - Calculating expression (ProbeSet Summarization)
#   - Remove invalid probes
#   - Save expression matrix
# ------------------------------------------------------------------------------


thisFile <- function() {
  #' Returns the location of the file from where this is run from.
  #' 
  #' @returns A string with the current path. The path separators might
  #'   be arbitrary. 
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}

# Note for devs: this will only work on r-studio if used in interactive
# mode.
if (!interactive()) {
  root <- dirname(thisFile())
} else {
  if (!require("rstudioapi")) {
    install.packages("rstudioapi")
  }
  root <- dirname(rstudioapi::getSourceEditorContext()$path)
}

## NOTE:: This file is one folder deeper than the root, so i need to go one
# up. This is fragile, but I cannot be bothered to make a better implementation.
root <- dirname(root)

# This will assure that the `renv` env is active. For instance if this
# script was started in `--vanilla` mode (as it should).
tryCatch(
  {
    cat("Attempting to activate the existing `renv` environment...")
    source(file.path(root, "renv", "activate.R"))
    cat("...OK\n")
  },
  error = function(err) {
    stop("No `renv` project found. Did you run the installation script?")
  }
)

# Load the helper functions.
source(file.path(root, "STALKER_Functions.R"))   # Collection of custom functions


affy2expression <- function(
    input.folder, output.file, log_name = NULL,
    use.affy = FALSE, remove.controls = FALSE
    ) {
  # ---- Loading packages ----
  # As a rule: try oligo first, then if it doesn't work, try affy
  # For 3' IVT Affymetrix Arrays (older ones; e.g. Mouse 430 and Human U133 series)
  # For GeneChip/Exon Affymetrix Arrays (newer ones + popular old platforms)
  if (use.affy) { library(affy) } else { library(oligo) }
  library(limma)          # For plotMD() function
  library(openxlsx)       # Read, Write, and Edit .xlsx (Excel) Files
  library(logger)
  
  # Setup logging facilities
  output.dir <- dirname(output.file)
  start.timedate <- gsub(" ", "_", date())
  log.target <- if (is.null(log_name)) {
    file.path(output.dir, paste0("affy2expression_", start.timedate, ".log"))
  } else {
    log_name
  }
  file.create(log.target)
  log_appender(appender_tee(log.target))
  
  # ---- Load .CEL files ----
  log_info("Looking for .CEL files...")
  setwd(input.folder)
  
  celFiles = list.celfiles()
  paste0(
      "Found ", length(celFiles), "celfiles: ", paste0(celFiles, collapse = ", ")
    ) |>
    log_info()
  
  log_info("Parsing found files to raw data...")
  # This should also install-and-load the platform design package (e.g. pd.hg.u133a)
  affyRaw = read.celfiles(celFiles)
  
  log_info("Running RMA normalization...")
  eset = rma(affyRaw)
  
  if (remove.controls) {
    log_info("Removing control probes...")
    # Remove Affymetrix control probes
    probes.before = dim(eset)[1]
    # ^ anchor to match the start of string (see regular expressions)
    eset = eset[-grep("^AFFX", rownames(eset)),]
    probes.fter = dim(eset)[1]
    discarded = probes.before - probes.fter
    log_info(paste0(discarded, " Affymetrix control probes have been discarded"))
    
    # Check for missing values (NA) and NaN entries
    if (any(is.na(exprs((eset)))) || any(is.nan(exprs((eset))))) {
      log_warn("Detectem some missing values in the dataset. Has somenthing gone terribly wrong?")
    }
  }
  log_info("Saving Expression Matrix to file...")
  write.exprs(eset, file = output.file)
  log_info("affy2expression finished")
}

