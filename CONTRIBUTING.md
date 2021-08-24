# Some pointers for devs

## Environment polluters
R makes .Rproj and more session files that are automatically loaded upon startup. This is bad as they contain `renv` activation lines that may be sourced before dependencies are loaded. This means two things:
    - Always run R and Rscript using the `--vanilla` option, so that these files are not loaded or saved.
    - Always `.gitignore` and `.dockerignore` these files (`*/.Rprofile` and `*/*.Rproj`) so that they do not propagate.
Together we can stop the spread of these files and the contamination of more environments!

