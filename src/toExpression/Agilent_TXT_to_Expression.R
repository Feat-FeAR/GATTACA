# Header Info ------------------------------------------------------------------
#
# Agilent TXT raw data --to--> expression matrix
#
# a FeAR R-script - 27-Dec-2021
#
# NOTE
# read.maimages() function from limma requires as its first argument a data frame
# containing (at least) a column called 'FileName' with the names of the raw-data
# files to be used for the analysis. For Agilent arrays such a list is usually
# given in the form of a tab-delimited file named "Targets.txt", possibly
# containing other kind of information about the experimental design. However,
# this script ignores any Targets.txt file found within the input directory and
# builds its own target list run-time.
#
# Script outline:
#   - Raw data loading
#   - 'normexp' background correction
#   - Quantile-Quantile interarray normalization
#   - Negative control probe evaluation
#   - Invalid probes removal
#   - Add annotation
#   - Expression matrix saving
#

log_debug("Sourcing the 'Agilent_TXT_to_Expression.R' file.")

agil2expression <- function (
  input_dir, output_file,
  grep_pattern="*.(txt|TXT)",
  remove_controls = TRUE,
  plot.width = 16,
  plot.height = 9,
  use.pdf = TRUE,
  n_plots = Inf
) {

  set.seed(1) # This module uses random sampling

  # Options for printPlots
  # Set options for printPlots
  options(
    "scriptName" = "prepagil",
    "save.PNG.plot" = !use.pdf,
    "save.PDF.plot" = use.pdf,
    "plot.width" = plot.width,
    "plot.height" = plot.height
  )

  # Inputting data
  output_dir <- dirname(output_file)
  setwd(output_dir)

  log_info("Finding input files matching pattern...")
  raw_files <- list.files(path = input_dir, pattern = grep_pattern)
  
  # Remove possible "Targets.txt" file from row_files list
  target.index <- grep("targets.txt", raw_files, ignore.case = TRUE)
  if (length(target.index) > 0) {
    raw_files = raw_files[-target.index]
  }
  
  log_info(
    paste("Found", length(raw_files), "input files:", paste(raw_files, collapse = ", "))
  )
  log_info("Reading in input files...")
  expression_data <- read.maimages(
    files = file.path(input_dir, raw_files),
    source = "agilent.median",
    green.only = TRUE
  )

  # Print MA plots to diagnose the data.
  print_data <- as.data.frame(expression_data$E)
  colnames(print_data) <- raw_files
  print_data <- log(print_data, 2)
  if (nrow(print_data) > 100000) {
    log_info("Reducing raw dataset...")
    print_data <- reservoir_sample(print_data, 100000)
    print_data <- as.data.frame(print_data)
  }
  ma.plots <- get_better_mas(print_data, title = "Raw probes MA plot - {x} vs Median of other samples")

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
    bplot <- ggplot(data = melt(print_data), aes(y = value, x = variable)) +
      geom_boxplot(outlier.alpha = 0.5, outlier.size = 1) +
      theme_bw() +
      scale_x_discrete(
        labels = 1:length(print_data)
      ) +
      ylab("Expression") + xlab("Sample") +
      ggtitle("Unnormalized Boxplots")
    print(bplot)
  }
  printPlots(suppressMessages(p), "Unnormalized Boxplots")

  rm(print_data)

  log_info("Finished inputting data.")

  # Normalization
  log_info("Running Normalization")
  log_info("Background correcting...")
  expression_data <- limma::backgroundCorrect(expression_data, method = "normexp", offset = 50)

  # This step also log2s the data
  log_info("Running interarray normalization...")
  expression_data = limma::normalizeBetweenArrays(expression_data, method = "quantile")

  expression_set <- expression_data$E
  log_info(paste("Final dataset dimensions - Cols:", ncol(expression_set), "Rows:", nrow(expression_set)))

  # Check the negative control probes to estimate the log2_intensity of
  # unhybridized spots, to be used as the threshold value for filtering in GATTACA
  neg.ctrl = expression_data$genes$ControlType == -1
  neg.id = unique(expression_data$genes$ProbeName[neg.ctrl])
  
  log_info(paste0(sum(neg.ctrl), " Negative-Control probes have been found, corresponding to ",
                  length(neg.id), " unique probe(s) [", neg.id, "]. Mean value of ",
                  length(expression_set[neg.ctrl,]), " unhybridized spots = ",
                  mean(expression_set[neg.ctrl,])))
  
  # Print MA plots to diagnose the data.
  colnames(expression_set) <- make.names(raw_files)
  ma.plots <- get_better_mas(as.data.frame(expression_set), title = "Normalized MA plot - {x} vs Median of other samples")

  if (n_plots != Inf) {
    ma.plots <- ma.plots[1:n_plots]
  }

  for (i in seq_along(ma.plots)) {
    maplot <- ma.plots[[i]]
    printPlots(\() { suppressMessages(print(maplot)) }, paste(i, "-", maplot$labels$title))
  }

  log_info("Making overall boxplot")
  p <- function(){
    bplot <- ggplot(data = melt(as.data.frame(expression_set)), aes(y = value, x = variable)) +
      geom_boxplot(outlier.alpha = 0.5, outlier.size = 1) +
      theme_bw() +
      scale_x_discrete(
        labels = 1:length(as.data.frame(expression_set))
      ) +
      ylab("Expression") + xlab("Sample") +
      ggtitle("Normalized Boxplots")
    print(bplot)
  }
  printPlots(suppressMessages(p), "Normalized Boxplots")

  if (remove_controls) {
    log_info("Filtering out control probes...")
    control_probes <- abs(expression_data$genes$ControlType) == 1
    control_probes_percent <- round(sum(control_probes) / nrow(expression_set) * 100, 4)
    log_info(paste(
      "Found", sum(control_probes), "control probes,", control_probes_percent,
      "% of total probes. Removing them..."
    ))

    expression_data <- expression_data[!control_probes, ]
    expression_set <- expression_data$E
    log_info(paste("Filtered dataset dimensions - Cols:", ncol(expression_set), "Rows:", nrow(expression_set)))
  }

  log_info("Finding replicate probes and collapsing them...")
  # Replace value of replicate probes with their mean - Probe_ID is used to identify the replicates
  expression_data = avereps(expression_data,  ID = expression_data$genes$ProbeName)

  expression_set <- expression_data$E
  log_info(paste("Final dataset dimensions - Cols:", ncol(expression_set), "Rows:", nrow(expression_set)))

  expression_set <- as.data.frame(expression_set)
  colnames(expression_set) <- make.names(raw_files)

  log_info("Saving output file...")

  write_expression_data(expression_set, output_file)

}
