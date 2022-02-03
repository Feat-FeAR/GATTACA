library(logger)
library(progress)

data <- read.csv("./biocpackages.csv")
sub_data <- data

sub_data <- sub_data[endsWith(sub_data$package_name, ".db"),]
key <- grepl("affymetrix", sub_data$info, ignore.case = TRUE) | grepl("agilent", sub_data$info, ignore.case = TRUE)
sub_data <- sub_data[key,]

key2 <- grepl("h", sub_data$package_name, ignore.case = TRUE)
sub_data <- sub_data[key2,]

key3 <- grepl("mouse", sub_data$package_name, ignore.case = TRUE) |
  grepl("rat", sub_data$package_name, ignore.case = TRUE) |
  grepl("chicken", sub_data$package_name, ignore.case = TRUE) |
  grepl("drosophila", sub_data$package_name, ignore.case = TRUE) |
  grepl("zebrafish", sub_data$package_name, ignore.case = TRUE) |
  grepl("probeset", sub_data$package_name, ignore.case = TRUE) |
  grepl("cdna", sub_data$info, ignore.case = TRUE)
sub_data <- sub_data[!key3,]

print(sort(sub_data$package_name))

# Manual removal of dupes
sub_data <- sub_data[! sub_data$package_name %in% c(
  "ath1121501.db" #"hgu133a.db", "hgu133b.db"
), ]

#' Get the possible annotation names of a specific database.
#'
#' This also installs and loads the database package from Bioconductor.
#'
#' @param db_namespace The name of the library containing the db
#'
#' @returns A string vector with the names of the available annotations.
#'
#' @author MrHedmad
get_db_names <- function(db_namespace) {
  suppressWarnings({
    stopifnot("Invalid database name - cannot be empty"=db_namespace==character(0))
    if (!require(db_namespace, character.only = TRUE)) {
      BiocManager::install(db_namespace, update = FALSE)
      suppressPackageStartupMessages(library(db_namespace, character.only = TRUE))
    }
  })
  possibilities <- ls(paste0("package:", db_namespace))
  db_name = gsub("\\.db", "", db_namespace)
  possibilities <- sapply(
    possibilities, gsub,
    pattern = db_name, replacement = ""
  )

  # Clean out the things that start with _ or . as they are functions and junk
  possibilities <- sapply(
    possibilities, gsub,
    pattern = "^[_\\.].*", replacement = "")

  possibilities <- possibilities[possibilities != ""]
  names(possibilities) <- NULL
  return(possibilities)
}


#' Get specified annotations from a certain db.
#'
#' Possible db_names:
#'   hgu133a.db (Affymetrix Human Genome U133 Set (A))
#'   hgu133b.db (Affymetrix Human Genome U133 Set (B))
#'   hgu133plus2.db (Affymetrix Human Genome HG-U133 Plus 2.0 Array)
#'   HsAgilentDesign026652.db (Agilent-026652 Whole Human Genome Microarray 4x44K v2)
#'   hugene10sttranscriptcluster.db (Affymetrix Human Gene 1.0-ST Array)
#'
#' All bioconductor databases can use the following possible selections:
#'   ACCNUM, CHR, CHRLOC, CHRLOCEND, ENSEMBL, ENTREZID, ENZYME, GENENAME, GO,
#'   MAP, OMIM, PATH, PMID, REFSEQ, SYMBOL, UNIPROT
#'
#' @param db_name The id of the db. See above and the CHIP_TO_DB object.
#' @param selections A vector of strings with the annotations to get from the db.
#'   The available annotations can be seen above or by using the `get_db_names`
#'   function
#' @returns A data.frame with probe IDs as rownames and one column per annotation.
#'
#' @author MrHedmad
get_remote_annotations <- function(
  db_name, selections = c("ACCNUM", "SYMBOL", "GENENAME")
) {
  library(purrr)
  # get_db_names also loads the db in memory, so I don't do it here.
  possible_selections <- get_db_names(db_name)
  log_info("Checking selections...")
  if (!all(selections %in% possible_selections)) {
    stop(
      paste0(
        "Invalid selection(s): ",
        paste(selections[!selections %in% possible_selections], collapse = ", "),
        ". Possible selections: ",
        paste(possible_selections, collapse = ", "),
        "."
      )
    )
  }

  # Load the data
  log_info("Loading the annotation data...")
  db_clean_name <- gsub("\\.db", "", db_name)
  data <- as.list(rep(NA, times = length(selections)))
  names(data) <- selections

  for (selection in selections) {
    # This raises deprecation warnings for no reason.
    suppressWarnings(
      { data[selection] <- get(paste0(db_clean_name, selection)) }
    )
  }

  merge_genedata <- function(x, y) {
    return(
      merge(x, y, by = "probe_id", all = TRUE)
    )
  }

  collapse_to_str <- function(inputlist) {
    sapply(
      contents(inputlist),
      function(x) {
        suppressWarnings({if (is.na(x)){x} else {paste(x, collapse = " /// ")}})}
    )
  }

  dataframetize <- function(named_vector, name) {
    probe_id <- names(named_vector)
    data <- data.frame()[1:length(named_vector), ]
    data[[name]] <- named_vector
    data$probe_id <- probe_id
    return(data)
  }

  # We need dataframes to manipulate
  log_info("Collapsing annotations...")
  data <- map(data, collapse_to_str)
  log_info("Casting annotations...")
  container <- as.list(rep(NA, length(data)))
  names(container) <- names(data)
  for (i in seq_along(data)) {
    container[[i]] <- dataframetize(data[[i]], name = names(data)[i])
  }
  data <- container
  rm(container)
  log_info("Collapsing annotations again...")
  data <- purrr::reduce(data, merge_genedata)

  # The merging functions need the probe_ids as rownames
  log_info("Cleaning up rownames...")
  row.names(data) <- data$probe_id
  data$probe_id <- NULL

  return(data)
}


get_all_data <- function(bin, packages) {
  write("ID,SYMBOL,GENENAME,ENSEMBL,package_name,version", bin)
  pb <- progress_bar$new(format = "[:bar] :current/:total (:percent)", total = length(packages))
  pb$tick(0)
  for (package in packages) {
    packdata <- get_remote_annotations(package, selections = c("SYMBOL", "GENENAME", "ENSEMBL"))
    packdata$package_name <- rep(package, length(packdata$SYMBOL))
    packdata$version <- rep(packageVersion(package), length(packdata$SYMBOL))

    write.table(packdata, bin, append=TRUE, sep=",", col.names = FALSE)
    pb$tick()
  }
}

get_all_data(bin="./annotdataraw.csv", packages = sub_data$package_name)
