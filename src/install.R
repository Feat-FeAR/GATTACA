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

old.wd <- getwd()
new.wd <- dirname(thisFile())

setwd(new.wd)
tryCatch(
  {
    print("Attempting to activate the existing `renv` environment...")
    source("./renv/activate.R")
    print("...OK")
  },
  error = function(err) {
    print("No existing project found. Making a new one...")
    renv::init(bare = TRUE)
    print("...OK")
  }
)
print("Restoring from renv.lock file... [[THIS COULD TAKE A VERY LONG TIME!!]]")
renv::restore()
print("Restore complete.")

setwd(old.wd)