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

# Header Info ------------------------------------------------------------------
#
# From CEL file to Expression matrix
# For Affymetrix 3'IVT OR Gene/Exon ST Arrays
#
# a FeAR R-script - 14-Sep-2021
#
# NOTE:
# On both Gene 1.0 ST and Exon 1.0 ST Arrays a probe set is more or less an exon.
# Doing the analysis on probe set level means to analyze signals for each exon.
# On the contrary, a 'transcript cluster' contains all the probe sets of a gene,
# and therefore can be used to measure gene expressions. Accordingly, to perform
# analyses on gene level, always use target="core" as summarization target in
# the following functions from oligo package and "Transcript Cluster Annotations"
# packages that contain all annotations of the gene each probe set belongs to.
#
# General Script for CEL-file normalization through RMA algorithm
#
#   - CEL file loading
#   - RMA normalization
#       1 - Background correcting
#       2 - Normalizing
#       3 - Calculating expression (ProbeSet Summarization)
#   - Remove invalid probes
#   - Save expression matrix
#
# ------------------------------------------------------------------------------

# Load the helper functions.
source(file.path(ROOT, "src", "STALKER_Functions.R"))
source(file.path(ROOT, "src", "annotator.R"))

affy2expression <- function(
    input.folder, output.file, remove.controls = TRUE,
    plot.width = 16, plot.height = 9, use.pdf = TRUE,
    n_plots = Inf
    ) {
  paste0(
    "Call (input.folder, output.file, remove.controls): ",
    paste(input.folder, output.file, remove.controls, sep = " :: ")
  ) |>
    log_debug()
  # ---- Loading packages ----
  graceful_load(c(
    "oligo",
    "limma",
    "affycoretools",
    "reshape2"
  ))

  # Set options for printPlots
  options(
    "scriptName" = "prepaffy",
    "save.PNG.plot" = !use.pdf,
    "save.PDF.plot" = use.pdf,
    "plot.width" = plot.width,
    "plot.height" = plot.height
  )

  # ---- Load .CEL files ----
  log_info("Looking for .CEL files...")
  output.folder <- dirname(output.file)
  setwd(input.folder)

  celFiles = list.celfiles()
  paste0(
    "Found ", length(celFiles), " .CEL files: ", paste0(celFiles, collapse = ", ")
    ) |>
    log_info()

  log_info("Parsing found files to raw data...")
  # This should also install-and-load the platform design package (e.g. pd.hg.u133a)
  affyRaw = read.celfiles(celFiles)

  # We can now move to where the plots will be saved.
  setwd(output.folder)

  # Check Platform and set the exon.probes flag
  # This is done as newer platforms need different processing than old ones.
  #
  # The parameter 'target' (only for Gene ST and Exon ST) defines the degree
  # of summarization, "core" is the default option which use transcript
  # clusters containing "safely" annotated genes. For summarization on the
  # exon level (not recommended for Gene arrays), use "probeset".
  platform <- affyRaw@annotation
  log_info(paste0("Detected platform: ", platform))
  if (platform %in% c("pd.hg.u133a", "pd.hg.u133b", "pd.hg.u133.plus.2")) {
    exon.probes = FALSE
  } else if (platform %in% c("pd.hugene.1.0.st.v1")) {
    exon.probes = TRUE
  } else {
    stop("Invalid or unsupported platform")
  }

  # Making plots for quality control
  log_info("Making MA plots before normalization...")
  if (exon.probes) {
    unnormalized_data <- oligo::rma(
      affyRaw, target = "core", normalize = FALSE, background = FALSE
    )
  } else {
    unnormalized_data <- oligo::rma(
      affyRaw, normalize = FALSE, background = FALSE
    )
  }

  unnormalized_data |> exprs() |> as.data.frame() -> unnormalized_data
  ma.plots <- get_better_mas(
    unnormalized_data,
    title = "Unnormalized MA plot - {x} vs Median of other samples"
  )

  if (n_plots != Inf) {
    stopifnot(
      "Invalid amount of plots to display"={is.wholenumber(n_plots)},
      "Number of plots to display is too high."={length(ma.plots) > n_plots}
    )
    ma.plots <- ma.plots[1:n_plots]
  }

  for (i in seq_along(ma.plots)) {
    maplot <- ma.plots[[i]]
    printPlots(\() { suppressMessages(print(maplot)) }, paste(i, "-", maplot$labels$title))
  }

  log_info("Making overall boxplot")
  p <- function(){
    bplot <- ggplot(data = melt(unnormalized_data), aes(y = value, x = variable)) +
      geom_boxplot(outlier.alpha = 0.5, outlier.size = 1) +
      theme_bw() +
      scale_x_discrete(
        labels = 1:length(unnormalized_data)
      ) +
      ylab("Expression") + xlab("Sample") +
      ggtitle("Unnormalized Boxplots")
    print(bplot)
  }
  printPlots(suppressMessages(p), "Unnormalized Boxplots")

  rm(unnormalized_data)
  # Done making plots

  log_info("Running RMA normalization...")
  if (exon.probes) {
    expression_set = oligo::rma(affyRaw, target = "core")
  } else {
    expression_set = oligo::rma(affyRaw)
  }

  log_info(paste0("Dataset dimensions: ", dim(expression_set)[1], " probe sets x ",
                  dim(expression_set)[2], " samples"))

  if (remove.controls) {
    log_info("Removing control probes...")
    # Remove Affymetrix control probes
    probes.before = dim(expression_set)[1]
    if (exon.probes) {
      expression_set = getMainProbes(expression_set, level = "core")
    } else {
      # ^ anchor to match the start of string (see regular expressions)
      ctrl.index = grep("^AFFX", rownames(expression_set))
      # To prevent an integer(0) from deleting all probes
      if (length(ctrl.index) > 0) {
        expression_set = expression_set[-ctrl.index,]
      }
    }
    probes.after = dim(expression_set)[1]
    discarded = probes.before - probes.after
    discarded.percent = round(discarded / probes.before * 100, 4)
    paste0(
      discarded, " Affymetrix control probes have been discarded, ",
      discarded.percent, "% of the total."
    ) |>
      log_info()
    log_info(paste0("Dataset dimensions: ", probes.after, " probe sets x ",
                    dim(expression_set)[2], " samples"))
    # Check for missing values (NA) and NaN entries
    if (any(is.na(exprs((expression_set)))) || any(is.nan(exprs((expression_set))))) {
      log_warn("Detected some missing values in the dataset. Has something gone terribly wrong?")
    }
  }

  expression_set |> exprs() |> as.data.frame() -> transposed

  ma.plots <- get_better_mas(transposed, title = "Normalized MA plot - {x} vs Median of other samples")

  if (n_plots != Inf) {
    ma.plots <- ma.plots[1:n_plots]
  }

  for (i in seq_along(ma.plots)) {
    maplot <- ma.plots[[i]]
    printPlots(\() { suppressMessages(print(maplot)) }, paste(i, "-", maplot$labels$title))
  }

  log_info("Making overall boxplot")
  p <- function(){
    bplot <- ggplot(data = melt(transposed), aes(y = value, x = variable)) +
      geom_boxplot(outlier.alpha = 0.5, outlier.size = 1) +
      theme_bw() +
      scale_x_discrete(
        labels = 1:length(transposed)
      ) +
      ylab("Expression") + xlab("Sample") +
      ggtitle("Normalized Boxplots")
    print(bplot)
  }
  printPlots(suppressMessages(p), "Normalized Boxplots")

  log_info("Saving Expression Matrix to file...")
  write_expression_data(transposed, output.file)
  log_info("affy2expression finished")
}
