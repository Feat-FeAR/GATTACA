# Some pointers for devs

Thank you for wanting to contribute to this repository. This file contains some guidelines for new submissions to the project, as well as well as a guide on how to work on the project locally.

## Environment polluters in the docker
R makes .Rproj and more session files that are automatically loaded upon startup. This is bad as they contain variables that migth have an effect on the analysis. This means two things:
    - Always run Rscript in the docker using the `--vanilla` option, so that these files are not loaded or saved.
    - Always `.gitignore` and `.dockerignore` these files (`*/.Rprofile` and `*/*.Rproj`) so that they do not propagate to the remote repository.
    - Remember to add the `-rm`

## Working and testing locally
It is useful to work interactively on the source code. To do this, you will
need the required libraries installed. The file
`/src/resources/r-requirements.txt` keeps a (hopefully up-to-date) list of the packages needed by
GATTACA to run. To install them, run something akin to this:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

packages <- read.table("/GATTACA/src/resources/r_requirements.txt")[, 1]

for (package in packages) {
  if (startsWith(package, "bioc::")) {
    package <- gsub("^bioc::", "", package)
    BiocManager::install(package)
  } else {
    install.packages(package)
  }
}
```

Afterwards, remember to `source` the `__init__.R` script before starting to work locally.

Even better, ignore all of this, and rebuild the docker each time you need to test the package. See the next section.

## Rebuilding the Docker
Docker saves intermediate containers in the cache, so that changes in the dockerfile will not (often) cause the container to rebuild completely.

What does this mean for us? It means that the very long installation process doesn't need to be restarted if the dockerfile section related to the installation of required packages is not changed.

So, feel free to change any other file without having to worry about lengthy recompilation. This also means that tests can (and should) be run in the container and not in interactive mode (like Rstudio).

Docker builds made during development should always be tagged with `bleeding`.

The usual workflow to work locally while rebuilding the docker instance is:
    - Make changes to the source files in `/src`;
    - Rebuild the docker with `docker build -t cmalabscience/gattaca:bleeding`;
    - Run `GATTACA -v bleeding ...` when testing.

