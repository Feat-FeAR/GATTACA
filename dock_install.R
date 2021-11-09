
options(
  repos = c("http://cran.us.r-project.org"),
  configure.args = list(
    preprocessCore = "--disable-threading"
  )
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  print("Installing Bioconductor")
  install.packages("BiocManager")
}

ppackages <- readLines("/GATTACA/r_requirements.txt")

# Filter out comments and blank lines
..is_valid_package <- function(x) {
  if (x == "" | startsWith(x, "#")) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

packages <- trimws(packages, which = "both")
packages <- packages[sapply(packages, ..is_valid_package)]

tot_packages <- length(packages)
current_package <- 1

print(paste(
  "Installing", tot_packages, "packages:", paste(packages, collapse = ", ")
))
for (package in packages) {
  print(paste("Installing package", current_package, "of", tot_packages, ":", package))

  if (startsWith(package, "bioc::")) {
    package <- gsub("^bioc::", "", package)
    BiocManager::install(package, update = FALSE)
  } else {
    install.packages(package)
  }
  current_package <- current_package + 1
}
