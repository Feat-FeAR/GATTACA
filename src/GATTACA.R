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
# GATTACA v3
# General Algorithm for The Transcriptional Analysis by one-Channel Arrays
#
# a FeAR R-script
#
# Pipeline for one-Color (HD) Microarrays
# Data are supposed to be already background-subtracted, log2-transformed, and
# interarray-normalized
#
# ------------------------------------------------------------------------------

# ---- Package Loading ----
# Load required packages and source external scripts

source(file.path(ROOT, "src", "STALKER_Functions.R"))   # Collection of custom functions

graceful_load(c(
  "preprocessCore",   # Interarray Normalization by Quantile-Quantile Algorithm
  "rafalib",          # Bland Altman Plots (aka MA Plots)
  "PCAtools",         # Principal Component Analysis
  "genefilter",       # Expression Gene Filtering
  "limma",            # Empirical Bayes Method for Differential Expression
  "RankProd",         # Rank Product Method for Differential Expression
  "VennDiagram",      # Venn Diagrams
  "EnhancedVolcano",  # Volcano Plots
  "gplots",           # Heatmap with extensions - heatmap.2() 
  "ggplot2",          # Box Plot and Bar Chart with Jitter (already loaded by PCAtools)
  "RColorBrewer",     # Color Palette for R - display.brewer.all()
  "yaml",             # Yaml file parsing
  "logger"           # Well, logging
))

# The annotation DB loads after dplyr so the `select` function gets overwritten.
# Why does R allow this?
dp_select <- dplyr::select

#' Run groupwise filtering on expression data.
#' 
#' The expression data *must* be named with colnames starting with the same
#' groups `groups`, for example by using `make.names(groups)`.
#' 
#' The data is filtered by keeping genes that are expressed more than
#' `expression_threshold` in at least `min_groupwise_presence` samples in
#' at least one group.
#' 
#' @param expression_set A data.frame with the data to be filtered.
#' @param groups A vector with the group names that correspond to the columns
#'   in the `expression_set`.
#' @param min_groupwise_presence The minimum groupwise presence used to filter.
#' @param expression_threshold The log2 expression treshold to use to consider
#'   a gene to be sufficiently expressed in a sample.
#'
#' @returns A data.frame with the filtered data.
filter_expression_data <- function(
  expression_set, groups, min_groupwise_presence, expression_threshold
) {
  log_info(
    "Filtering a ", ncol(expression_set), " cols by ", nrow(expression_set),
    " rows expression set."
  )

  original_dims <- dim(expression_set)

  report <- paste(
    "Initial dims: ", ncol(expression_set), " cols by ", original_dims[1], " rows",
    "\nMin presence per group: ", min_groupwise_presence,
    "\nExpression threshold (log2): ", expression_threshold
    )

  log_debug(
    "Raw colnames:", paste(colnames(expression_set), collapse = ", "),
    "\nGroups:", paste(groups, collapse = ", ")
  )

  group_sizes <- table(groups)
  min_presences <- ceiling(group_sizes * min_groupwise_presence)

  table_report <- data.frame(
    groups = names(group_sizes),
    sizes = as.vector(group_sizes),
    min_presence = as.vector(round(group_sizes * min_groupwise_presence, 2)),
    actual_min_presence = as.vector(min_presences)
  )

  report <- paste(report, get.print.str(table_report), sep = "\n")

  # To filter, we do three things. First, we look at the data and split it
  # groupwise. Secondly, we apply to each split, row-wise, a function that
  # gives us TRUE if at least `min_presence` columns have expression over
  # `expression_treshold` expression. We save these vectors into a data frame,
  # the `truth_dataframe`, with a col for each group. Finally, we collapse
  # this dataframe to a single truth key, with one entry per row, by applying
  # `any` row-wise, thus keeping an entry if the function described above is
  # true in at least one group.
  truth_dataframe <- data.frame(row.names = rownames(expression_set))

  for (group in unique(groups)) {
    # I cannot use "filter" as we first need to compare all cols to see if
    # the rowwise test passes, and only then we can filter.
    filter_fun <- function(frame) {
      kOverA_filter <- kOverA(min_presences[[group]], expression_threshold)
      frame |> apply(MARGIN = 1, FUN = kOverA_filter) ->
        filter_key
      return(filter_key)
    }
    
    expression_set |> dp_select(starts_with(group)) |> filter_fun() ->
      truth_dataframe[[group]]
  }

  truth_dataframe |> apply(2, sum) -> nr_kept_genes
  report <- paste(
    report, "--------------\nRetained genes:", get.print.str(nr_kept_genes),
    "Percentage:", get.print.str(nr_kept_genes / original_dims[1] * 100),
    sep = "\n"
  )

  truth_dataframe |> apply(1, any) -> truth_key
  expression_set <- expression_set[truth_key, ]

  report <- paste(
    report,
    paste("-----------\nFiltered dims:", ncol(expression_set), "cols by ", nrow(expression_set), "rows"),
    paste("Percentage retained genes:", nrow(expression_set) / original_dims[1] * 100),
    sep = "\n"
    )
  
  log_info(report)
  return(expression_set)
}


GATTACA <- function(options.path, input.file, output.dir) {
  # This function is so long, a description wouldn't fit here.
  # Refer to the project's README.

  # ---- Option parsing ----
  opts <- yaml.load_file(options.path)
  
  # ---- Making static functions ----
  pitstop <- pitstop.maker(opts$general$slowmode)
  printdata <- printif.maker(opts$general$show_data_snippets, topleft.head)
  
  log_info("Parsed options.")
  
  # ---- Move to output.dir ----
  setwd(output.dir)

  # ---- Variable setup ----
  # To reuse as-is the script, I unpack the variables from the yaml file here
  # TODO : Improvement - These should be checked here for basic validity, so
  #        we crash sooner rather than later. 
  log_info("Inputting variables...")

  user_colours <- opts$design$group_colors
  
  # Log2 expression threshold
  min_log2_expression = opts$design$filters$log2_expression
  
  # Fold Change Threshold
  thrFC = opts$design$filters$fold_change
  
  # Flags for script-wide IFs
  write_data_to_disk = !opts$switches$dryrun # The ! is important.
  if (opts$switches$dryrun & opts$general$save_png) {
    log_warn("save_png forced to be FALSE as this is a dryrun")
    opts$general$save_png <- FALSE
  }
  if (opts$switches$dryrun & opts$general$save_pdf) {
    log_warn("save_pdf forced to be FALSE as this is a dryrun")
    opts$general$save_pdf <- FALSE
  }
  
  # Global options suitable for STALKER_Functions
  options(
    scriptName = "GATTACA",
    save.PNG.plot = opts$general$save_png,
    save.PDF.plot = opts$general$save_pdf,
    use.annotations = !is.null(opts$general$annotation_chip_id),
    plot.width = opts$general$plot_width,
    plot.height = opts$general$plot_height,
    png_ppi = opts$general$png_resolution,
    enumerate.plots = opts$general$enumerate_plots
  )
  
  if (getOption("use.annotations")) {
    # If we need to use annotations, we need to load the annotation data
    source(file.path(ROOT, "src", "annotator.R"))
    annotation_data <- get_remote_annotations(
      CHIP_TO_DB[[opts$general$annotation_chip_id]], "SYMBOL"
    )
    log_info("Loaded annotations")
  } else {
    log_info("No annotations loaded")
  }

  # ---- Data Loading ----
  # Gene Expression Matrix - log2-Intensity-Values
  log_info("Loading data...")

  expression_set <- read.csv(input.file, header = TRUE, row.names = "probe_id")
  log_info(paste(
    "Loaded expression_set dimensions:", paste0(dim(expression_set), collapse = ", ")
  ))
  printdata(expression_set)

  pitstop("Finished loading data.")
  
  # ---- Load the experimental design ----
  log_info("Loading experimental design...")
  experimental_design <- design_parser(opts$design$experimental_design)
  
  # This is a list with as `groups` the groups and as `pairings` the 
  # ID for pairing
  experimental_design <- split_design(experimental_design)
  
  # Check the design
  if (length(experimental_design$groups) != ncol(expression_set)) {
    stop(paste0(
      "The number of patients (", ncol(expression_set),
      ") does not match the number of design groups (", length(experimental_design$groups), ")"
    ))
  }
  
  # Test if we are in paired or unpaired mode
  ..pairing_NAs <- is.na(experimental_design$pairings)
  # We need either ALL or NONE patient NAs
  if (all(..pairing_NAs)) {
    log_info("No sample pairing detected. Running in unpaired mode.")
    paired_mode <- FALSE
  } else if (!all(..pairing_NAs)) {
    log_info("Found sample pairings. Running in paired mode.")
    paired_mode <- TRUE
  } else {
    stop("Some samples have pairing data and some do not. Cannot proceed with pairing ambiguity.")
  }
  
  log_info("Experimental desing loaded.")

  # Tidy Sample Names According to the Experimental Design
  # Create a new vector containing tidy group names
  log_info("Making tidy group names using the experimental design...")
  unique_simple_groups <- unique(experimental_design$groups)
  unique_groups <- make.unique(experimental_design$groups, sep = "_")
  ..old_colnames <- colnames(expression_set)
  colnames(expression_set) <- unique_groups
  
  if (length(user_colours) < length(unique_simple_groups)) {
    stop("Too few colors in \'user_colours\' vector!")
  }
  # Bind colours to groups
  names(user_colours[1:length(unique_simple_groups)]) <- unique_simple_groups
  
  # Save Correspondences Table so we can check it later
  ..corrTable = cbind(..old_colnames, colnames(expression_set)) # Cast to matrix
  colnames(..corrTable) = c("Original_ID", "Group_Name")
  printdata(..corrTable)

  if (write_data_to_disk) {
    write.table(
      ..corrTable,
      "Corresp Table.tsv", sep = "\t",
      col.names = TRUE, row.names = TRUE, quote = FALSE
    )
    log_info(paste("'correspondence_table.tsv' has been saved in", output.dir))
  }
  
  # Check the contrasts
  raw_contrasts <- opts$design$contrasts
  
  unpack_contrasts <- function(contrasts_list) {
    lapply(contrasts_list, strsplit, split = "-") |> 
      lapply(unlist) ->
      unpacked
    return(unpacked)
  }
  
  raw_contrasts |> unpack_contrasts() -> group_contrasts
  
  # Check if the contrasts groups are actually in the data
  ..check <- unlist(group_contrasts) %in% unique_simple_groups
  if (! all(..check)) {
    stop(paste0(
      "Some groups defined in the contrasts are not present in the data. ",
      "Conflicting groups: ", paste0(group_contrasts[!..check])
    ))
  }
  
  log_info("Design loaded and approved.")
  pitstop("")


  # ---- Preliminary Boxplots ----
  # Distribution Inspection to Spot Non-Logarithmic Data
  printPlots(\() {boxplot(expression_set, las = 2)}, "Raw boxpot")
  printPlots(\() {plotDensities(expression_set, legend = FALSE)}, "Raw density")
  

  # ---- Normalization ----
  # After-RMA 2nd Quantile Normalization
  if (opts$switches$renormalize) {
    expression_set <- qq_normalize(expression_set)
  }


  # ---- MA-Plot & Box-Plot ----
  # Normalization Final Check with Figure Production
  printPlots(function() {
    boxplot(
      expression_set,
      las = 2, col = user_colours,
      main = "Expression values per sample", ylab = "log2 (intesity)"
      )
    },
    "Final Boxplot"
  )
  printPlots(function() {
    plotDensities(expression_set, legend = FALSE, main = "Expression values per sample")
  }, "Final Density")
  
  pitstop("Maybe check the plots and come back?")
  
  # MA-Plot for bias detection
  # From limma package: array 'column' vs the average of all the arrays other than that
  printPlots(function(){plotMD(expression_set, column = 1)}, "Mean-Difference Plot")
  
  for (combo in combn(unique_simple_groups, 2, simplify = FALSE)) {
    expression_set |> group_by(across(starts_with(combo[[1]]))) |>
      rowMeans() -> group1
    expression_set |> group_by(across(starts_with(combo[[2]]))) |>
      rowMeans() -> group2

    p <- function() {
      maplot(
        group1, group2, 
        xlab = "A (Average log-expression)", ylab = "M (Expression log-ratio)",
        n = 5e4,
        curve.add = TRUE, curve.col = user_colours[2], curve.lwd = 1.5, curve.n = 1e4,
        pch = 20, cex = 0.1
      )
      title(main = paste(combo[[1]], "vs", combo[[2]]))
      abline(h = 0, col = user_colours[1], lty = 2) # lty = line type
      abline(h = c(1,-1), col = user_colours[1])
      abline(v = min_log2_expression, col = user_colours[1]) # Platform-specific log2-expression threshold
    }
    printPlots(p, paste0("MA-Plot", combo[[1]], "_vs_", combo[[2]]))
  }


  # ---- Clustering ----
  log_info("Starting sample-wise hierarchical clustering for batch-effect detection...")
  # Matrix Transpose t() is used because dist() computes the distances between
  # the ROWS of a matrix
  # Distance matrix (NOTE: t(expression_set) is coerced to matrix)
  expression_set |> t() |> dist() |> hclust(method = "ward.D2") ->
    hierarchical_clusters

  printPlots(\(){plot(hierarchical_clusters)}, "Dendrogram")
  # Desired number of clusters
  # TODO : Should this be an option?
  nr_shown_clusters = 6

  p <- function() {
    plot(hierarchical_clusters)
    # Red borders around the kNum clusters 
    rect.hclust(
      hierarchical_clusters, k = nr_shown_clusters, border = user_colours[2]
    )
  }
  printPlots(p, "Dendrogram and Clusters")
  
  log_info("Finished with hierarchical clustering.")
  pitstop("Maybe check the plots and come back?")
  

  # ---- PCA ----
  # Performed on Samples for Batch-Effect Detection
  log_info("Performing PCA to detect batch sampling...")
  # Bundle some metadata in the PCA object for later.
  # Strictly enforced that rownames(metadata) == colnames(expression_set)
  metadata = data.frame(
    groups = experimental_design$groups,
    row.names = unique_groups
  )

  # Do the PCA (centering the data before performing PCA, by default)
  PCA_object = pca(expression_set, metadata = metadata)
  log_info("Finished running PCA. Plotting results...")
  # Plot results
  printPlots(\(){print(screeplot(PCA_object))}, "Scree Plot")

  p <- function() {
    suppressMessages(print(biplot(PCA_object, colby = "groups", colkey = user_colours)))
  }
  printPlots(p, "PCA")
  p <- function() {
    suppressMessages(print(pairsplot(PCA_object, colby = "groups", colkey = user_colours)))
  }
  printPlots(p, "PCA Pairs")
  pitstop(paste0(
    "Take a look at the PCA results.",
    " If there are some batched samples, remove them and re-run GATTACA."
    ))
  
  # ---- SD vs Mean Plot ----
  # Poisson Hypothesis Check
  log_info("Using Poisson to produce a SD vs Mean plot...")
  
  # Un-log intensity values
  log_info("Returning to linear intensities...")
  unlogged_expression_set = 2 ^ expression_set
  printdata(unlogged_expression_set)
  pitstop("Did the 'unlogging' mess anything up?")
  
  # Store values using matrices
  log_info("Making matrices...")
  ..nr_unique_simple_groups <- length(unique_simple_groups)
  
  matrix(nrow = nrow(expression_set), ncol = ..nr_unique_simple_groups + 1) |>
    as.data.frame() -> ..means_frame
  colnames(..means_frame) <- c(unique_simple_groups, "Global")
  
  # Make a copy in memory
  ..sd_matrix <- data.frame(..means_frame) 

  ..correlation_vector <- as.vector(rep(NA, length(unique_simple_groups) + 1))
  names(..correlation_vector) <- c(unique_simple_groups, "Global")
  
  # Statistics for each group...
  log_info("Calculating groupwise statistics...")
  for (group in unique_simple_groups) {
    unlogged_expression_set |> dp_select(starts_with(group)) |>
      rowMeans(na.rm = TRUE) ->
      ..means_frame[[group]]
    
    unlogged_expression_set |> dp_select(starts_with(group)) |>
      apply(1, sd, na.rm = TRUE) ->
      ..sd_matrix[[group]]
    
    cor(..means_frame[[group]], ..sd_matrix[[group]]) ->
      ..correlation_vector[[group]]
  }
  
  # ...and for the whole experiment
  log_info("Calculating global statistics...")

  ..means_frame[["Global"]] = rowMeans(unlogged_expression_set, na.rm = TRUE)
  ..sd_matrix[["Global"]] = apply(unlogged_expression_set, 1, sd, na.rm = TRUE)
  ..correlation_vector[["Global"]] = cor(..means_frame[,"Global"], ..sd_matrix[,"Global"])
  
  # Scatter plot
  log_info("Making plots...")
  p <- function() {
    par(mfrow = c(1, length(unique_simple_groups)+1))
    X.max = max(..means_frame, na.rm = TRUE)
    Y.max = max(..sd_matrix, na.rm = TRUE)
    for (group in unique_simple_groups) {
      
      plot(..means_frame[[group]], ..sd_matrix[[group]],
           xlab = "Mean", ylab = "SD",
           xlim = c(0, X.max), ylim = c(0, Y.max),
           pch = 20, cex = 0.1)
      title(main = group)
      mtext(side = 3, paste("Corr =", toString(round(..sd_matrix[[group]], digits = 5))))
    }
    par(mfrow = c(1, 1))
  }
  printPlots(p, "SD_vs_Mean Plot") # Save just the 'Global' one
  
  # ---- Filtering ----

  expression_set <- filter_expression_data(
    expression_set, groups = experimental_design$groups,
    min_groupwise_presence = opts$design$filters$min_groupwise_presence,
    expression_threshold = min_log2_expression
  )

  log_info("Done filtering.")


  # ---- DE by Limma ----
  if (opts$switches$limma) {
    log_info("Running DEA with `limma`...")
    pitstop("")

    # Differential Expression Assessment and Figure Production
    log_info("Running differential expression analysis in paired mode.")

    log_info("Making limma design matrix...")
    
    if (paired_mode) {
      ..factor_pairings <- as.factor(experimental_design$pairings)
      ..pairings_levels <- levels(..factor_pairings)
      ..pairings_levels <- ..pairings_levels[1:length(..pairings_levels)-1]
      ..pairings_levels <- make.names(..pairings_levels)
      ..groups_levels <- sort(unique_simple_groups)
      ..limma_design <- model.matrix(
        ~ 0 + experimental_design$groups + ..factor_pairings
      )
      colnames(..limma_design) <- c(..groups_levels, ..pairings_levels)
    } else {
      ..limma_design <- model.matrix(~ 0 + experimental_design$groups)
      colnames(..limma_design) <- unique_simple_groups
    }

    log_info("Design matrix:\n", get.print.str(..limma_design))

    log_info("Making contrasts matrix...")
    makeContrasts(
      contrasts = raw_contrasts,
      levels = ..limma_design
    ) -> ..contrast_matrix

    log_info("Contrasts matrix:\n", get.print.str(..contrast_matrix))

    log_info("Computing contrasts...")
    lmFit(expression_set, ..limma_design) -> ..limma_fit
    ..limma_fit |> contrasts.fit(..contrast_matrix) |> eBayes() -> ..limma_Bayes

    # Print Results (Top-Ten genes) for all the contrasts of interest
    # This is only useful in slowmode.
    if (opts$general$slowmode) {
      for (i in seq_along(raw_contrasts)) {
        # 'print' because automatic printing is turned off in loops (and functions)
        # `cat` is also OK here as we are in slowmode.
        cat("\nDEG Top-List for contrast: ", raw_contrasts[i], "\n", sep = "")
        print(topTable(..limma_Bayes, coef = i, adjust.method = "BH", sort.by = "B"))
        cat("\n") # The one time I like a bit more space.
      }
    }
    pitstop("Next: Save all DEGs.")

    log_info("Getting Differential expressions...")
    progress <- txtProgressBar(min = 0, max = length(raw_contrasts))
    DEGs.limma = list() # Create an empty list
    for (i in seq_along(raw_contrasts)) {
      DEGs.limma[[i]] = topTable(
        ..limma_Bayes, coef = i, number = Inf,
        adjust.method = "BH", sort.by = "B"
      ) # this is a list of Data Frames
      setTxtProgressBar(progress, i)
    }
    close(progress)

    # Save full DEG Tables
    if (write_data_to_disk) {
      log_info("Saving Differential expression tables...")
      progress <- txtProgressBar(min = 0, max = length(raw_contrasts))
      ..limma_output_expression <- DEGs.limma[[i]]
      ..limma_output_expression$probe_id <- rownames(..limma_output_expression)
      rownames(..limma_output_expression) <- NULL
      for (i in seq_along(raw_contrasts)) {
        degTabName = paste("Limma - DEG Table ", raw_contrasts[i], ".csv", sep = "")
        write.csv(..limma_output_expression, degTabName, row.names = FALSE, quote = FALSE)
        setTxtProgressBar(progress, i)
      }
      close(progress)
      log_info(paste0("Saved tables in ", output.dir))
    }
    
    # Summary of DEGs (you can change the Log2-Fold-Change Threshold lfc...)
    results.limma = decideTests(
      ..limma_Bayes, adjust.method = "BH", p.value = 0.05,
      lfc = thrFC
    )
    p <- function() {
      vennDiagram(results.limma)
    }
    printPlots(p, "Limma Venn")
    
    # Show Hyperparameters
    d0 = ..limma_Bayes$df.prior           # prior degrees of freedom
    dg = mean(..limma_fit$df.residual)    # original degrees of freedom
    hyp = cbind(
      c(..limma_Bayes$s2.prior,     # prior variance
      mean(..limma_fit$sigma^2),    # mean sample residual variance
      mean(..limma_Bayes$s2.post),  # mean posterior residual variance
      d0, dg, d0/(d0+dg))   # Shrinkage degree
    )
    rownames(hyp) = c(
      "Prior Var",
      "Sample Residual Mean Var",
      "Posterior Residual Mean Var",
      "Prior df",
      "Original df",
      "Shrinkage degree"
    )
    colnames(hyp) = "Hyperparameters"
    log_info(paste0("Hyperparameters:: \n", get.print.str(hyp)))

    pitstop("Finished running DEA.")

    # ---- Limma Plot ----
    log_info("Making DEG plots...")
    # MA-Plots with significant DEGs and Volcano plots
    
    # Find Axis Limits
    log_info("Finding axis limits...")
    max.M.value = 0
    for (i in seq_along(raw_contrasts)) {
      temp = max(abs(DEGs.limma[[i]]$logFC))
      if (temp > max.M.value) {
        max.M.value = temp
      }
    }
    max.A.value = 0
    for (i in seq_along(raw_contrasts)) {
      temp = max(DEGs.limma[[i]]$AveExpr)
      if (temp > max.A.value) {
        max.A.value = temp
      }
    }
    min.A.value = Inf
    for (i in seq_along(raw_contrasts)) {
      temp = min(DEGs.limma[[i]]$AveExpr)
      if (temp < min.A.value) {
        min.A.value = temp
      }
    }
    min.P.value = 1
    for (i in seq_along(raw_contrasts)) {
      temp = min(DEGs.limma[[i]]$P.Value)
      if (temp < min.P.value) {
        min.P.value = temp
      }
    }
    
    # MA-Plot with DEGs
    log_info("Making Limma MA plots...")
    # I cannot place a progress bar here as `printPlots` logs to stdout.
    for (i in seq_along(raw_contrasts)) {
      # Mark in red/blue all the up-/down- regulated genes (+1/-1 in 'results.limma' matrix)
      p <- function() {
        plotMD(
          ..limma_Bayes, status = results.limma[,i],
          values = c(1,-1), hl.col = user_colours[c(2,1)],
          xlab = "A (Average log-expression)",
          ylab = "M (log2-Fold-Change)"
        )
        abline(h = 0, col = user_colours[1], lty = 2) # lty = line type
        abline(h = c(thrFC,-thrFC), col = user_colours[1])
        abline(v = min_log2_expression, col = user_colours[1]) # Platform-specific log2-expression threshold
      }
      printPlots(p, paste("MA-Plot with Limma DEGs ", raw_contrasts[i], sep = ""))
    }
    
    # Volcano Plots
    log_info("Making Limma Volcano plots...")
    for (i in seq_along(raw_contrasts)) {
      tot.DEG = sum(DEGs.limma[[i]]$adj.P.Val < 0.05) # Total number of significant DEGs (without any FC cutoff)
      high.DEG = min(c(5,tot.DEG)) # To highlight no more than 5 genes per plot
      
      # Significance Threshold
      # TODO : Discuss if this is right
      # Find that p-value corresponding to BH-adj.p-value ~ 0.05 (or Bonferroni point when tot.DEG = 0)
      thrP = (0.05/nrow(expression_set))*(tot.DEG + 1)
      # TODO : Implement this?
      # Alternative approach - Suitable also for correction methods other than BH
      # WARNING: DEG list has to be sorted by p-value or B statistics!
      #if (tot.DEG > 0) {
      #  thrP = DEGs.limma[[i]][tot.DEG + 1, "P.Value"] # This is the p-value of the first non-DEG
      #} else {
      #  thrP = (0.05/filtSize) # Bonferroni point
      #}
      # Check the threshold p-value
      log_info(
        paste0(
          "\n", raw_contrasts[i], " - Threshold Report:\n",
          "--------------------------------------\n",
          "  p-value threshold  =  ", thrP, "\n",
          "  -log10(p-value)    =  ", -log10(thrP), "\n",
          "  Gene Ranking       =  ", tot.DEG, ":", tot.DEG + 1, "\n\n", sep = ""
        )
      )
      
      # Enhanced Volcano Plot
      if (getOption("use.annotations")) {
        ..annotation_data <- merge_annotations(DEGs.limma[[i]], annotation_data)
        ..volcano_labels <- ..annotation_data$SYMBOL
      } else {
        ..volcano_labels = rownames(DEGs.limma[[i]])
      }
      p <- function() {
        # NOTICE: When in a for loop, you have to explicitly print your
        # resulting EnhancedVolcano object
        suppressWarnings(
          # This prints warnings as - I think - the internal implementation
          # uses xlim and ylim but ggplot2 ignores them.
          print(
            EnhancedVolcano(
              DEGs.limma[[i]],
              x = "logFC", y = "P.Value",
              pCutoff = thrP, FCcutoff = thrFC,
              pointSize = 1,
              col = c("black", "black", "black", user_colours[2]),
              lab = ..volcano_labels,
              #selectLab = myLabels[1:high.DEG],
              labSize = 4,
              title = raw_contrasts[i],
              subtitle = "Limma",
              legendPosition = "none"
            )
          )
        )
      }
      printPlots(p, paste("Volcano with Limma DEGs ", raw_contrasts[i], sep = ""))

    }
  }
  
  # ---- DE by RankProduct ----
  if (opts$switches$rankproduct) {
    log_info("Running differential expression analysis with rankproduct...")
    log_info("Running in unpaired mode.")
    # Differential Expression Assessment and Figure Production
    results.RP = matrix(
      data = 0,
      nrow = nrow(expression_set), ncol = length(raw_contrasts)
    )
    rownames(results.RP) = rownames(expression_set)
    colnames(results.RP) = raw_contrasts

    for (i in seq_along(raw_contrasts)) {
      log_info(paste0("Running analysis ", i, " of ", length(raw_contrasts)))
      
      log_info("Finding control and case groups...")
      ..unpacked_groups = strsplit(raw_contrasts[i], split = "-", fixed = TRUE)[[1]]
      
      ..control_group <- ..unpacked_groups[1]
      ..treated_group <- ..unpacked_groups[2]
      
      expression_set |> dp_select(
        starts_with(..control_group) | starts_with(..treated_group)
        ) -> ..sub_expression_set
      
      # The above filter rearranges the columns. This crashes later as the
      # column order is needed to match the pairings for subtraction.
      # So we return to the original one here.
      ..original_col_order <- colnames(expression_set)[
        colnames(expression_set) %in% colnames(..sub_expression_set)
      ] 
      ..sub_expression_set <- ..sub_expression_set[, ..original_col_order]
      
      # Make the cl array, representing class labels of the sample
      colnames(..sub_expression_set) |> startsWith(..treated_group) |> as.numeric() ->
        ..rp_class_labels
      
      if (paired_mode) {
        # To make a paired run, we need to subtract the matching paired sets.
        # First, we need to get the pairings that survived the filtering
        colnames(expression_set) %in% colnames(..sub_expression_set) ->
          ..surviving_cols
        ..sub_parings <- experimental_design$pairings[..surviving_cols]

        # We make an empty dataframe to hold the subtraction result
        ..subtracted_expression_set <- data.frame(row.names = rownames(..sub_expression_set))
        for (pairing in unique(..sub_parings)) {

          ..paired_expression_set <- ..sub_expression_set[..sub_parings == pairing]
          
          ..paired_expression_set |> dp_select(starts_with(..control_group)) ->
            ..control_set
          ..paired_expression_set |> dp_select(starts_with(..treated_group)) ->
            ..treated_set
          
          if (ncol(..control_set) != 1) {
            log_warn(
              paste0(
                "There are more than 1 '", ..control_group, 
                "' paired samples (",
                paste(colnames(..control_set), collapse = ", "), ").",
                " The average of them will be used instead."
              )
            )
            ..control_set <- rowMeans(..control_set)
            colnames(..control_set) <- paste0("averaged_", ..control_group)
          } else if (ncol(..control_set) == 0) {
            stop(
              paste0(
                "There are no paired columns for group '", ..control_group,
                "' for pairing ", pairing
              )
            )
          }
          
          if (ncol(..treated_set) != 1) {
            log_warn(
              paste0(
                "There are more than 1 '", ..treated_group, 
                "' paired samples (",
                paste(colnames(..treated_set), collapse = ", "), ").",
                " The average of them will be used instead."
              )
            )
            ..treated_set <- rowMeans(..treated_set)
            colnames(..treated_set) <- paste0("averaged_", ..treated_group)
          } else if (ncol(..treated_set) == 0) {
            stop(
              paste0(
                "There are no paired columns for group '", ..treated_group,
                "' for pairing ", pairing
              )
            )
          }
          
          # Now, control and treated sets can be subtracted as they are one
          # and only one column.
          ..subtracted_partial_set <- ..control_set - ..treated_set
          
          ..subtracted_expression_set <- merge(
            ..subtracted_expression_set, ..subtracted_partial_set,
            all = TRUE, by = 0
          )
          # the merge makes the rownames into their own col
          rownames(..subtracted_expression_set) <- ..subtracted_expression_set$Row.names
          ..subtracted_expression_set$Row.names <- NULL
        }
        log_info("Finished pairing dataset.")
        ..sub_expression_set <- ..subtracted_expression_set
        ..rp_class_labels <- rep(1, ncol(..subtracted_expression_set))
      }

      # invisible(capture.output()) is to suppress automatic output to console
      # WARNING: therein <- (instead of =) is mandatory for assignment!
      log_info("Running RankProduct...")
      RP.out <- RP.advance(
        ..sub_expression_set,
        ..rp_class_labels, origin = rep(1, ncol(..sub_expression_set)),
        gene.names = rownames(expression_set),
        logged = TRUE, na.rm = FALSE, plot = FALSE, rand = 123
      )
      cat("\n")
      p <- function() {
        invisible(capture.output(plotRP(RP.out, cutoff = 0.05)))
      }
      printPlots(p, paste("Rankprod ", raw_contrasts[i]))
      
      # Compute full DEG Tables (returns a list of 2 matrices, not data frames)
      log_info("Computing DEG table...")
      invisible(
        capture.output(
          # The Inf cutoff is to print all genes
          DEGs.RP <- topGene(RP.out, logged = TRUE, logbase = 2, cutoff = Inf)
        )
      )
      for (j in 1:2) {
         # Invert FC to get Case vs Ctrl and take the log2 values
        DEGs.RP[[j]][,3] = log2(1/DEGs.RP[[j]][,3])
        colnames(DEGs.RP[[j]])[3] = "Log2FC" # Correct column name
      }
      
      # Print Results (Top-Ten genes) for all the contrasts of interest
      log_info("Getting top contrasts...")
      log_info(paste0("DEG Top-List for contrast: ", raw_contrasts[i]))
      tops = rbind(DEGs.RP$Table1[1:10,], DEGs.RP$Table2[1:10,])
      printdata(tops) # just on-Screen
      pitstop("")
      
      # Save full DEG Tables
      if (write_data_to_disk) {
        log_info("Saving full DEG tables...")
        upDegTabName = paste0("RP_Up - DEG Table ", raw_contrasts[i], ".csv")
        dwnDegTabName = paste0("RP_Down - DEG Table ", raw_contrasts[i], ".csv")
        ..rpd_DEGs_up <- DEGs.RP[[1]]
        ..rpd_DEGs_up$probe_id <- rownames(..rpd_DEGs_up)
        rownames(..rpd_DEGs_up) <- NULL
        write.csv(..rpd_DEGs_up, upDegTabName, row.names = FALSE, quote = FALSE)
        ..rpd_DEGs_down <- DEGs.RP[[2]]
        ..rpd_DEGs_down$probe_id <- rownames(..rpd_DEGs_down)
        rownames(..rpd_DEGs_down) <- NULL
        write.csv(..rpd_DEGs_down, dwnDegTabName, row.names = FALSE, quote = FALSE)
        log_info(paste0(
          upDegTabName, "' and '", dwnDegTabName,
          "' have been saved in ", output.dir
        ))
      }
      
      # Fetch indexes of significant DEGs
      log_info("Finding DEGs of interest...")
      ups.Index = which(DEGs.RP$Table1[,4] < 0.05 & DEGs.RP$Table1[,3] > thrFC)
      dwn.Index = which(DEGs.RP$Table2[,4] < 0.05 & DEGs.RP$Table2[,3] < -thrFC)
      log_info(
        paste(
          "Found", length(ups.Index), "upregulated and",
          length(dwn.Index), "downregulated genes."
        )
      )
      
      # Here 'results.RP' is filled in analogy to 'results.limma':
      # (-1,0,+1) for (Down, NotSig, Up)-regulated genes
      # Accessing matrix rows by row name: matrix[rowname, col_i]
      results.RP[rownames(DEGs.RP$Table1)[ups.Index],i] = 1
      results.RP[rownames(DEGs.RP$Table2)[dwn.Index],i] = -1
      
    }
    
    log_info("Calculating summary statistics...")
    summary.RP = rbind(
      colSums(results.RP == -1), # Cast to matrix
      colSums(results.RP == 0),
      colSums(results.RP == 1)
    )
    rownames(summary.RP) = c("Down", "NotSig", "Up")
    print(summary.RP)
    
    log_info("Finished DEA.")
    pitstop("")
  }
  
  if (opts$switches$limma & opts$switches$rankproduct) {
    log_info("Making comparison plots between limma and rankproduct...")
    # Venn diagrams of DEGs from Limma and RP methods
    
    # To suppress 'venn.diagram()' logging
    futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
    
    # Plot Venn diagrams
    log_info("Plotting Venn diagrams...")
    for (i in seq_along(raw_contrasts)) {
      for (j in c(1, -1)) {

        DEG.id.limma = rownames(results.limma)[which(results.limma[,i] == j)]
        DEG.id.RP = rownames(results.RP)[which(results.RP[,i] == j)]
        
        if (j == 1) {
          venn.sub = "UP-regulated DEGs"
        } else {
          venn.sub = "DOWN-regulated DEGs"
        }
        
        if (length(DEG.id.limma) == 0 & length(DEG.id.RP) == 0) {
          log_warn("Both the sets are empty in the contrast: ", raw_contrasts[i], " (", venn.sub, ")")
          next # Skip the current iteration of the for loop without terminating it
        }
        
        venn.plot = venn.diagram(
          x = list(DEG.id.limma, DEG.id.RP),
          filename = NULL, # to print just on screen
          force.unique = TRUE,
          main = raw_contrasts[i], main.cex = 2, main.fontface = "bold", main.fontfamily = "sans", # Title
          sub = venn.sub, sub.fontfamily = "sans", # Subtitle
          lwd = 2, lty = "blank", fill = user_colours[1:2], # circles
          cex = 2, fontface = "bold", fontfamily = "sans", # numbers
          category.names = c("Limma", "Rank Product"), # names
          cat.cex = 2,
          cat.fontface = "bold",
          cat.default.pos = "outer",
          cat.pos = c(-150, 150),
          cat.dist = c(0.055, 0.055),
          cat.fontfamily = "sans"
        )
        
        # Create a new canvas and draw the Venn
        p <- function() {
          grid.newpage()
          grid.draw(venn.plot)
        }
        printPlots(
          p,
          paste0(
            "Comparison Venn ", raw_contrasts[i], "_",
            strsplit(venn.sub, split = "-")[[1]][1]
          )
        )
      }
    }
  }
  
  log_info("GATTACA finished")
}
