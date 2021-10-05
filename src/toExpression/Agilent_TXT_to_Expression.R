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
library(logger)
library(limma)

source(file.path(ROOT, "src", "STALKER_Functions.R"))


agil2expression <- function (
  input_dir, output_file, analysis_program,
  grep_pattern="*.txt", green_only = FALSE,
  offset = 0,
  remove_controls = TRUE,
  log_name = NULL
) {
  # Setup logging facilities
  output_dir <- dirname(output_file)
  start.timedate <- gsub(" ", "_", date())
  
  # I don't know if log_name was passed as a string or NULL, so in the call
  # made by entry I have to wrap the input in "" to make it a valid string.
  # This causes NULL to become "NULL", and therefore I have to do this badness
  if (!is.null(log_name)) {
    log_name <- if (log_name == "NULL") {NULL} else {log_name}
  }
  log.target <- if (is.null(log_name)) {
    file.path(output_dir, paste0("affy2expression_", start.timedate, ".log"))
  } else {
    file.path(output_dir, log_name)
  }
  file.create(log.target)
  log_appender(appender_tee(log.target))
  
  
  # Inputting data
  output_dir <- dirname(output_file)
  
  log_info("Finding input files matching pattern...")
  raw_files <- list.files(path = input_dir, pattern = grep_pattern)
  
  log_info(
    paste("Found", length(raw_files), "input files:", paste(raw_files, collapse = ", "))
  )
  log_info("Reading in input files...")
  raw_data <- read.maimages(
    files = file.path(input_dir, raw_files),
    source = analysis_program,
    green.only = green_only
  )
  
  log_info("Finished inputting data.")
  
  # Normalization
  log_info("Running Normalization")
  log_info("Background correcting...")
  raw_BGcorrected <- backgroundCorrect(raw_data, method = "normexp", offset = offset)
  
  log_info("Running interarray normalization...")
  raw_BGandNormalized = normalizeBetweenArrays(raw_BGcorrected, method = "quantile")
  
  expression_set <- raw_BGandNormalized$E
  log_info(paste("Final dataset dimensions - Cols:", ncol(expression_set), "Rows:", nrow(expression_set)))
  
  if (remove_controls) {
    log_info("Filtering out control probes...")
    control_probes <- abs(raw_BGandNormalized$genes$ControlType) == 1
    control_probes_percent <- round(sum(control_probes) / nrow(expression_set) * 100, 4)
    log_info(paste(
      "Found", sum(control_probes), "control probes,", control_probes_percent,
      "% of total probes. Removing them..."
    ))
    
    filtered_data <- raw_BGandNormalized[!control_probes, ]
    expression_set <- filtered_data$E
    log_info(paste("Filtered dataset dimensions - Cols:", ncol(expression_set), "Rows:", nrow(expression_set)))
  }
  
  log_info("Finding replicate probes and collapsing them...")
  # Replace value of replicate probes with their mean - Probe_ID is used to identify the replicates
  filtered_data = avereps(filtered_data,  ID = filtered_data$genes$ProbeName)
  
  expression_set <- filtered_data$E
  log_info(paste("Final dataset dimensions - Cols:", ncol(expression_set), "Rows:", nrow(expression_set)))
  
  expression_set <- as.data.frame(expression_set)
  colnames(expression_set) <- make.names(raw_files)
  
  log_info("Saving output file...")
  
  expression_set$probe_id <- rownames(expression_set)
  rownames(expression_set) <- NULL
  
  write.csv(expression_set, file = output_file, quote = FALSE, row.names = FALSE)
  
}
