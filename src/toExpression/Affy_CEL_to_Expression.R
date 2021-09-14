# Header Info ------------------------------------------------------------------
#
# Affymetrix CEL files -> Expression matrix
#
# a FeAR R-script
# based on: "Homer MiniTutorial - Processing Affymetrix Gene Expression Arrays"
# see pdf for more info
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
# ------------------------------------------------------------------------------

# Load the helper functions.
source(file.path(ROOT, "STALKER_Functions.R"))   # Collection of custom functions
source(file.path(ROOT, "annotator.R"))

affy2expression <- function(
    input.folder, output.file, log_name = NULL,
    use.affy = FALSE, remove.controls = TRUE
    ) {
  # ---- Loading packages ----
  library(oligo)
  library(limma)          # For plotMD() function
  library(openxlsx)       # Read, Write, and Edit .xlsx (Excel) Files
  library(logger)
  library(affycoretools)

  print(paste0(
    "Call: (in/out/log/affy/rm) ",
    input.folder, output.file, log_name,
    use.affy, remove.controls, sep = " :: "
  ))
  
  # Setup logging facilities
  output.dir <- dirname(output.file)
  start.timedate <- gsub(" ", "_", date())
  
  # I don't know if log_name was passed as a string or NULL, so in the call
  # made by entry I have to wrap the input in "" to make it a valid string.
  # This causes NULL to become "NULL", and therefore I have to do this badness
  log_name <- if (log_name == "NULL") {NULL} else {log_name}
  log.target <- if (is.null(log_name)) {
    file.path(output.dir, paste0("affy2expression_", start.timedate, ".log"))
  } else {
    file.path(output.dir, log_name)
  }
  file.create(log.target)
  log_appender(appender_tee(log.target))
  
  # ---- Load .CEL files ----
  log_info("Looking for .CEL files...")
  setwd(input.folder)
  
  celFiles = list.celfiles()
  paste0(
      "Found ", length(celFiles), " celfiles: ", paste0(celFiles, collapse = ", ")
    ) |>
    log_info()
  
  log_info("Parsing found files to raw data...")
  # This should also install-and-load the platform design package (e.g. pd.hg.u133a)
  affyRaw = read.celfiles(celFiles)
  
  # Check Platform and set the exon.probes flag
  # This is done as newer platforms need different processing than old ones.
  platform <- affyRaw@annotation
  log_info(paste0("Detected platform: ", platform))
  if (platform %in% c("pd.hg.u133a", "pd.hg.u133b", "pd.hg.u133.plus.2")) {
    exon.probes = FALSE
  } else if (platform %in% c("pd.hugene.1.0.st.v1")) {
    exon.probes = TRUE
  } else {
    stop("Invalid or unsupported platform")
  }
  
  if (exon.probes) {
    log_info("Running RMA normalization...")
    # The parameter 'target' (only for Gene ST and Exon ST) defines the degree
    # of summarization, "core" is the default option which use transcript
    # clusters containing "safely" annotated genes. For summarization on the
    # exon level (not recommended for Gene arrays), use "probeset".
    expression_set = rma(affyRaw, target = "core")
  } else {
    expression_set = rma(affyRaw)
  }
  
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
      discarded.percent, "% of total."
    ) |>
      log_info()
    # Check for missing values (NA) and NaN entries
    if (any(is.na(exprs((expression_set)))) || any(is.nan(exprs((expression_set))))) {
      log_warn("Detected some missing values in the dataset. Has something gone terribly wrong?")
    }
  }
  
  # Transposing this object is a pain. But `write.exprs` actually does what we
  # need. So I do this badness.
  log_info("Transposing and extracting probe ids...")
  temp <- tempfile()
  on.exit(unlink(temp)) # This assures that the tempfile is deleted no matter what
  write.exprs(expression_set, file = temp)

  transposed <- read.table(file = temp, row.names = 1)
  
  transposed$probe_id <- row.names(transposed)
  row.names(transposed) <- NULL
  
  log_info("Saving Expression Matrix to file...")
  write.csv(transposed, file = output.file, row.names = FALSE)
  log_info("affy2expression finished")
}
