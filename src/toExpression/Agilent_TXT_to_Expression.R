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


# ---- Load Required libraries ----

source(file.path(ROOT, "src", "STALKER_Functions.R"))

graceful_load(c("limma"))

agil2expression <- function (
  input_dir, output_file, analysis_program,
  grep_pattern="*.txt",
  offset = 0,
  remove_controls = TRUE
) {
  # Inputting data
  output_dir <- dirname(output_file)
  
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
