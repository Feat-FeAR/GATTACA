# This script loads the virtual env (renv) so it has to be sourced at the
# beginning of every docker call. It also sets the working dir to here, so we
# can give relative paths in relative safety.

# Why not call it directly? Because I want to make it platform independent.

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

ROOT <- dirname(thisFile())

# This will assure that the `renv` env is active. For instance if this
# script was started in `--vanilla` mode (as it should).
tryCatch(
  {
    cat("Attempting to activate the existing `renv` environment...")
    source(file.path(ROOT, "renv", "activate.R"))
    cat("...OK\n")
  },
  error = function(err) {
    stop("No `renv` project found. Did you run the installation script?")
  }
)

setwd(ROOT)
