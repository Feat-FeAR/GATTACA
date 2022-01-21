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

log_debug("Sourcing the 'annotations.R' file.")


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


#' Annotate a dataframe with data from GEO.
#'
#' @param expression_set A data.frame with row.names the probe ids.
#' @param chip_id The GEO accession for the platform ID.
#'
#' @return A dataframe with additional annotation columns.
annotate_data <- function(expression_set, chip_id, selections) {
  log_info("Finding annotations...")
  annotations <- get_annotations_from_GEO(chip_id)

  # As annotations is a data.frame, both the is.na and the if make warnings
  # i suppress them for a tiny little bit
  options(warn=-1)
  if (is.na(annotations)){
    stop("No valid annotations found. Cannot proceed.")
  }
  options(warn=1)

  # Check the matching degree between array and annotation layouts
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

  # Print out the number of NAs in the annotations for each type of annotation
  notMap = matrix(0, nrow = 2, ncol = dim(annotations)[2],
                  dimnames = list(c("NA entries","%"), colnames(annotations)))
  for (i in colnames(annotations)) {
    notMap[1,i] = sum(is.na(annotations[, i]))
    notMap[2,i] = round(notMap[1,i]/dim(annotations)[1]*1e2, digits = 2)
  }
  log_info(paste("Missing annotations:", get.print.str(notMap), sep = "\n"))

  log_info("Merging annotations with the data...")
  merged_data <- merge_annotations(expression_set, annotations)

  return(merged_data)
}

#' Annotate a file containing expression data with annotations from GEO.
#'
#' The file needs to be in .csv format with a column named "probe_id" representing
#' probe ids. The output file is similarly formatted.
#'
#' @param expression_data_path A str with the path to the csv file containing
#'   the data to be annotated.
#' @param output_path A full path to the output file.
#' @param chip_id The GEO id of the platform to source annotations from.
#' @param selections A vector of str with valid annotation names to use to
#'   annotate the data.
annotate_to_file <- function(expression_data_path, output_path, chip_id) {
  log_info("Loading input data...")
  expression_set <- read_expression_data(expression_data_path)

  log_info("Annotating data...")
  annotated_set <- annotate_data(expression_set, chip_id, selections)

  write_expression_data(annotated_set, output_path)
  log_info("Written annotations.")
}


#' Get annotations directly from GEO with GEOquery
#'
#' This gets annotations directly from GEO with GEOquery. The data is then
#' unpacked, info logged, and returned as a data.frame ready to be merged.
#'
#' The probe_id column is named just "ID".
#'
#' @param geo_id The Geo ID to download from. Enforced to be a GEO Platform ID
#'          (starting with "GPL").
#' @return A data.frame with the annotations, and at least the "ID" column.
get_annotations_from_GEO <- function(geo_id) {
  log_info("Getting annotations for platform ID ", as.character(geo_id), "...")

  if (! startsWith(geo_id, "GPL")) {
    log_error("Invalid GEO Platform ID. Must begin with the 'GPL' identifier.")
    return(NA)
  }

  annotations <- tryCatch(
    getGEO(geo_id),
    error=\(e){
      log_error("Cannot find an annotation for GEO id ", as.character(geo_id), ".")
      return(NA)
    },
    warning=\(w){
    log_warn(w)
    }
  )

  if (! is(annotations, "GPL")) {
    log_error("Downloaded wrong type of data. Expected GPL object.")
    return(NA)
  }

  log_info("Done downloading from GEO.")

  # There is a specific slot but I don't trust it 100%
  column_names <- colnames(annotations@dataTable@table)
  log_info("GPL-INFO - Platform name: ", annotations@header$title)
  log_info("GPL-INFO - Manufacturer: ", annotations@header$manufacturer)
  log_info("GPL-INFO - Available data: ", paste(column_names, collapse = ", "))

  # Check if we have at least the ID column
  if (! "ID" %in% column_names) {
    log_error("This GPL does not come with the ID column. Cannot parse.")
    return(NA)
  }
  data <- annotations@dataTable@table
  rownames(data) <- data$ID
  data$ID <- NULL

  return(data)
}

