# Running in vanilla mode can cause the CRAN mirrors not to be set.
# This fixes it.

if (interactive()) {
  stop("I think you are sourcing this file manually. This isn't supposed to be done.")
}

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

print("Making new renv...")
renv::init(bare = TRUE)
print("Restoring from renv.lock file... [[THIS COULD TAKE A VERY LONG TIME!!]]")
renv::restore()
print("Restore complete.")
