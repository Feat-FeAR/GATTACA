#!/usr/bin/env Rscript

# License ----------------------------------------------------------------------
# MIT License
#
# Copyright (c) 2021 Feat-FeAR
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
# ------------------------------------------------------------------------------


# This list maps shortcodes representing chips to their respective databases.
CHIP_TO_DB <- list(
  # Affymetrix Human Genome U133 Set (A)
  "hgu133a" = "hgu133a.db",
  # Affymetrix Human Genome U133 Set (B)
  "hgu133b" = "hgu133b.db",
  # Affymetrix Human Genome HG-U133 Plus 2.0 Array
  "hgu133plus2" = "hgu133plus2.db",
  # Agilent-026652 Whole Human Genome Microarray 4x44K v2
  "HsAgilentDesign026652" = "HsAgilentDesign026652.db",
  # Affymetrix Human Gene 1.0-ST Array
  "hugene10st" = "hugene10sttranscriptcluster.db"
)


#' Merge expression matrices with annotations and sort them.
#' 
#' Fuses by row names.
#' 
#' @param gene.stat The table of genes, usually a DEG summary-statistic top-table
#'   or an expression matrix.
#' @param annotation the matrix containing the annotation data
#' @param sort.by the name or index of the column used to sort the final data set
#' 
#' @author FeAR
merge_annotations <- function(gene.stat, annotation, sort.by = 1) {
  # 'merge' function to merge two matrix-like objects horizontally
  # and cast to data frame (right outer join).
  # NOTICE: both gene.stat and annotation are supposed to have the Probe_IDs
  # as row names
  joined = merge(
    annotation, gene.stat,
    by.x = "row.names", by.y = "row.names",
    all.y = TRUE
  )
  # The merge has to convert the row names to a column. This reverts it.
  rownames(joined) = joined[,1]
  gene.stat = joined[,-1]
  
  # Re-sort the data frame by the content of 'sort.by' column
  # ('sort.by' can be either a number or a column name)
  gene.stat = gene.stat[order(gene.stat[,sort.by]),]
  
  return(gene.stat)
}


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
    if (!require(db_namespace, character.only = TRUE)) {
      install.packages(paste0("bioc::", db_namespace))
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
  library(logger)
  library(purrr)
  # get_db_names also loads the db in memory, so I don't do it here.
  possible_selections <- get_db_names(db_name)
  log_info("Checking selections...")
  if (!all(selections %in% possible_selections)) {
    stop(
      paste0(
        "Invalid selection(s): ",
        paste(selection[!selections %in% possible_selections], collapse = ", "),
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
      {
        data[selection] <- get(paste0(db_clean_name, selection))
      }
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

#' Annotate a dataframe with data from a chip.
#' 
#' See the `get_remote_annotations` function for info on the possible chip_ids
#' and selections.
#' 
#' @param expression_set A data.frame with row.names the probe ids.
#' @param chip_id A str representing the chip used by the experiment.
#' @param selections A vector of str with valid annotation names to use to
#'   annotate the data.
#'   
#' @returns A dataframe with `selections` additional columns containing the
#'   selections.
annotate_data <- function(expression_set, chip_id, selections) {
  log_info("Checking chip selection...")
  database_name <- CHIP_TO_DB[[chip_id]]
  if (is.null(database_name)) {
    stop(paste0("Invalid chip id: ", chip_id))
  }
  
  log_info("Finding annotations...")
  annotations <- get_remote_annotations(database_name, selections = selections)
  
  # Check how many probes have no annotations
  missing_probes <- sum(
    ! row.names(expression_set) %in% row.names(annotations)
  )
  missing_perc <- missing_probes / length(expression_set) * 100
  if (missing_probes > 0) {
    log_warn(paste0(
      "There are ", missing_probes, " probes with no annotations. ",
      round(missing_perc, 4), "% of total."
    ))
  }
  
  # TODO : Maybe print out the number of NAs in the annotations for each
  # type of annotation?
  
  log_info("Merging annotations with the data...")
  merged_data <- merge_annotations(expression_set, annotations)
  
  return(merged_data)
}

#' Annotate a file containing expression data with annotations from a db.
#' 
#' The file needs to be in .csv format with a column named "probe_id" representing
#' probe ids. The output file is similarly formatted.
#' 
#' See the `get_remote_annotations` function for info on the possible chip_ids
#' and selections.
#' 
#' @param expression_data_path A str with the path to the csv file containing
#'   the data to be annotated.
#' @param output_path A full path to the output file.
#' @param chip_id A str representing the chip used by the experiment.
#' @param selections A vector of str with valid annotation names to use to
#'   annotate the data.
#' @param log_nale Name of the logfile written in the same directory as the 
#'   output file.
annotate_to_file <- function(
  expression_data_path, output_path,
  chip_id, selections,
  log_name = NULL
) {
  paste0(
    "Call: (exprpath/outputpath/chip/selections/logname):\n",
    paste(expression_data_path, output_path, chip_id, paste(selections, collapse = ", "), log_name,
    sep = " :: ")
  ) |>
    print()
  
  # Setup logging facilities
  output.dir <- dirname(output_path)
  start.timedate <- gsub(" ", "_", date())
  
  library(logger)
  
  log.target <- if (is.null(log_name)) {
    file.path(output.dir, paste0("GATTACA_annotator_", start.timedate, ".log"))
  } else {
    file.path(output.dir, log_name)
  }
  file.create(log.target)
  log_appender(appender_tee(log.target))
  
  log_info("Loading input data...")
  expression_set <- read.csv(file = expression_data_path, row.names = "probe_id")
 
  log_info("Annotating data...")
  annotated_set <- annotate_data(expression_set, chip_id, selections)
  
  log_info("Extracting probe ids...")
  annotated_set$probe_id <- row.names(annotated_set)
  
  write.csv(annotated_set, file = output_path, row.names = FALSE)
  log_info("Written annotations.")
}
