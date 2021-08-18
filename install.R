if (!requireNamespace("here")) {
  print("Installing `here` package...")
  install.packages("here")
  print("...OK")
}

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

