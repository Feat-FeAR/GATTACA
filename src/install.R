# Running in vanilla mode can cause the CRAN mirrors not to be set.
# This fixes it.

options(
  repos = list(CRAN="http://cran.rstudio.com/"),
  configure.args = list(
    preprocessCore = "--disable-threading"
  )
)

if (!requireNamespace("renv")) {
  print("Installing `renv` package...")
  if (!requireNamespace("remotes")) {
    print("Installing `remotes` package...")
    install.packages("remotes")
    print("...OK")
  }
  remotes::install_github("rstudio/renv")
  print("...OK")
}

thisFile <- function() {
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

here <- dirname(thisFile())
setwd(here)

tryCatch(
  {
    cat("Attempting to activate the existing `renv` environment...")
    source("./renv/activate.R")
    cat("...OK")
  },
  error = function(err) {
    cat("/nNo existing project found. Making a new one...")
    renv::init(bare = TRUE)
    cat("...OK")
  }
)
print("Restoring from renv.lock file... [[THIS COULD TAKE A VERY LONG TIME!!]]")
renv::restore()
print("Restore complete.")
