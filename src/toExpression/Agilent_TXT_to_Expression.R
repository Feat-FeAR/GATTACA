# Header Info ------------------------------------------------------------------
#
# Agilent TXT file + Target File -> Expression matrix
#
# a FeAR R-script - 24-Mar-2021
# based on: "Mannheimia MiniTutorial - Processing Agilent Arrays"
# see pdf for more info
#
# NOTE:
# The "target file" is just a tab-delimited text file created by the user,
# containing the experimental design and featuring a column called 'FileName'.
#
# General Script for Agilent-file normalization through normexp+quantile algorithm
#
#   - TXT-target file loading
#   - Raw data reading
#   - 'normexp' background correcting
#   - Quantile-Quantile interarray normalization
#   - Remove invalid probes
#   - Add annotation
#   - Save expression matrix
#

log_debug("Sourcing the 'Agilent_TXT_to_Expression.R' file.")

agil2expression <- function (
  input_dir, output_file, analysis_program,
  grep_pattern="*.txt",
  offset = 0,
  remove_controls = TRUE,
  plot.width = 16,
  plot.height = 9,
  use.pdf = TRUE,
  n_plots = Inf
) {

  set.seed(1) # This module uses random sampling

  paste0(
    "Call (input_dir, output_file, analysis_program, grep_pattern, offset, remove_controls, plot.width, plot.heigth, use.pdf): ",
    paste(input_dir, output_file, analysis_program, grep_pattern, offset, remove_controls, plot.width, plot.height, use.pdf, sep = " :: ")
  ) |>
    log_debug()

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

  log_info(
    paste("Found", length(raw_files), "input files:", paste(raw_files, collapse = ", "))
  )
  log_info("Reading in input files...")
  expression_data <- read.maimages(
    files = file.path(input_dir, raw_files),
    source = analysis_program,
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
  expression_data <- backgroundCorrect(expression_data, method = "normexp", offset = offset)

  # This step also log2s the data
  log_info("Running interarray normalization...")
  expression_data = normalizeBetweenArrays(expression_data, method = "quantile")

  expression_set <- expression_data$E
  log_info(paste("Final dataset dimensions - Cols:", ncol(expression_set), "Rows:", nrow(expression_set)))

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
