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


thisFile <- function() {
  #' Returns the location of the file from where this is run from.
  #' 
  #' @returns A string with the current path. The path separators might
  #'   be arbitrary. 
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}

# Note for devs: this will only work on r-studio if used in interactive
# mode.
if (!interactive()) {
  root <- dirname(thisFile())
} else {
  if (!require("rstudioapi")) {
    install.packages("rstudioapi")
  }
  root <- dirname(rstudioapi::getSourceEditorContext()$path)
}

# This will assure that the `renv` env is active. For instance if this
# script was started in `--vanilla` mode (as it should).
tryCatch(
  {
    cat("Attempting to activate the existing `renv` environment...")
    source(file.path(root, "renv", "activate.R"))
    cat("...OK\n")
  },
  error = function(err) {
    stop("No `renv` project found. Did you run the installation script?")
  }
)

# ---- Package Loading ----
# Load required packages and source external scripts

library(preprocessCore)   # Interarray Normalization by Quantile-Quantile Algorithm
library(rafalib)          # Bland Altman Plots (aka MA Plots)
library(PCAtools)         # Principal Component Analysis
library(genefilter)       # Expression Gene Filtering
library(limma)            # Empirical Bayes Method for Differential Expression
library(RankProd)         # Rank Product Method for Differential Expression
library(VennDiagram)      # Venn Diagrams
library(openxlsx)         # Reading, Writing, and Editing of .xlsx (Excel) Files
library(EnhancedVolcano)  # Volcano Plots
library(gplots)           # Heatmap with extensions - heatmap.2() 
library(ggplot2)          # Box Plot and Bar Chart with Jitter (already loaded by PCAtools)
library(RColorBrewer)     # Color Palette for R - display.brewer.all()
library(yaml)             # Yaml file parsing
library(logger)           # Well, logging

source(file.path(root, "STALKER_Functions.R"))   # Collection of custom functions


GATTACA <- function(options.path, input.file, output.dir) {
  # This function is so long, a description wouldn't fit here.
  # Refer to the project's README.

  # ---- Option parsing ----
  opts <- yaml.load_file(options.path)
  

  # ---- Making static functions ----
  pitstop <- pitstop.maker(opts$general$slowmode)
  printdata <- printif.maker(opts$general$show_data_snippets, topleft.head)
  
  # Setup logging facilities
  start.timedate <- gsub(" ", "_", date())
  log.target <- if (is.null(opts$general$log_name)) {
    file.path(output.dir, paste0("GATTACA_", start.timedate, ".log"))
  } else {
    file.path(output.dir, opts$general$log_name)
  }
  file.create(log.target)
  log_appender(appender_tee(log.target))
  
  log_info("Parsed options.")
  
  # ---- Prepare Annotations ----
  log_info("Grabbing annotations from database...")
  if (opts$chip$use_remote) {
    # NOTE : This is why languages should have a `match` statement...
    # Using a remote annotation database
    if (opts$chip$type == "Affymetrix Human Genome U133 Set (A)") {
      log_info("Using Affymetrix Human Genome U133 Set (A)")
      library(hgu133a.db)
      db.data <- list(
        acc = hgu133aACCNUM,
        symb = hgu133aSYMBOL,
        gname = hgu133aGENENAME
      )
    } else if (opts$chip$type == "Affymetrix Human Genome U133 Set (B)") {
      log_info("Using Affymetrix Human Genome U133 Set (B)")
      library(hgu133b.db)
      db.data <- list(
        acc = hgu133bACCNUM,
        symb = hgu133bSYMBOL,
        gname = hgu13bGENENAME
      )
    } else if (opts$chip$type == "Affymetrix Human Genome HG-U133 Plus 2.0 Array") {
      log_info("Using Affymetrix Human Genome HG-U133 Plus 2.0 Array")
      library(hgu133plus2.db)
      db.data <- list(
        acc = hgu133plus2ACCNUM,
        symb = hgu133plus2SYMBOL,
        gname = hgu133plus2GENENAME
      )
    } else if (opts$chip$type == "Agilent-026652 Whole Human Genome Microarray 4x44K v2") {
      log_info("Using Agilent-026652 Whole Human Genome Microarray 4x44K v2")
      library(HsAgilentDesign026652.db)
      db.data <- list(
        acc = HsAgilentDesign026652ACCNUM,
        symb = HsAgilentDesign026652SYMBOL,
        gname = HsAgilentDesign026652GENENAME
      )
    } else if (opts$chip$type == "Affymetrix Human Gene 1.0-ST Array") {
      log_info("Using Affymetrix Human Gene 1.0-ST Array")
      library(hugene10sttranscriptcluster.db)
      db.data <- list(
        acc = hugene10sttranscriptclusterACCNUM,
        symb = hugene10sttranscriptclusterSYMBOL,
        gname = hugene10sttranscriptclusterGENENAME
      )
    }

    annotation <- data.frame(
      Accession = sapply(contents(db.data$acc), paste, collapse = ", "),
      GeneSymbol = sapply(contents(db.data$symb), paste, collapse = ", "),
      Description = sapply(contents(db.data$gname), paste, collapse = ", ")
    )

    missing.annotations <- sum(annotation[,2] == "NA")

  } else {
    # Using local annotation database
    log_info("Using a local database: `", opts$chip$database_path, "`")
    tryCatch(
      {
        log_info("Attempting to load the database...")
        annotation <- read.xlsx(
          opts$chip$database_path,
          colNames = TRUE, rowNames = TRUE,
          sep.names = "_")
        log_info("Done loading. Attempting to use database...")
        annotation <-
          annotation[,c("Representative_Public_ID", "Gene_Symbol", "Gene_Title")]
        missing.annotations <- sum(annotation[,2] == "---")
      }
    )
  }
  log_info("Done loading annotations.")
  missing.percentage <- round(missing.annotations/dim(annotation)[1]*1e2, digits = 2)
  
  ## Reached milestone: loaded annotations
  log_info(
    paste0(
      missing.annotations, " unannotated genes (", missing.percentage, " %).\n"
    )
  )
  pitstop("Finished loading annotations.")
  

  # ---- Variable setup ----
  # To reuse as-is the script, I unpack the variables from the yaml file here
  # TODO : Improvement - These should be checked here for basic validity, so
  #   we crash sooner rather than later. 
  log_info("Inputting variables...")
  myFolder <- output.dir
  myFile <- input.file
  
  # TODO : Remove these.
  # Row offset, including the header (How many rows above the first number?)
  rowOffset = 1
  # Column containing (unique) gene identifiers
  colWithID = 1
  # Column offset, including row names (How many columns before the first number?)
  colOffset = 1
  
  # Experimental Design - Group Names (start with control condition)
  groups <- opts$design$experimental_groups
  # Experimental Design vector
  design <- opts$design$experimental_design
  # The design vector will be later on tested for elongation (Compact design) 
  # and correctness (both design modes)
  
  myColors <- opts$design$group_colors
  if (length(myColors) < length(groups)) {
    stop("Too few colors in \'myColors\' vector!")
  }
  
  # Log2 expression threshold
  thr0 = opts$design$filters$log2_expression
  
  # Fold Change Threshold
  thrFC = opts$design$filters$fold_change
  
  # Flags for script-wide IFs
  saveOut = !opts$switches$dryrun # The ! is important.
  logConversion = opts$switches$log2_transform
  secondNorm = opts$switches$renormalize
  
  # Global options suitable for STALKER_Functions
  options(
    scriptName = "GATTACA",
    save.PNG.plot = opts$general$save_png,
    save.PDF.plot = opts$general$save_pdf,
    append.annot = opts$general$append_annotation,
    plot.width = opts$general$plot_width,
    plot.height = opts$general$plot_height
  )
  
  if (getOption("append.annot")) {
    # If 'annotation' is not defined an error message will be displayed
    dim.annot = dim(annotation)
    log_info(
      paste0(
        "A ", dim.annot[1], " x ", dim.annot[2],
        " annotation dataframe has been loaded."
      )
    )
  } else {
    annotation = NULL
    log_info("No annotation loaded.")
  }
  
  ### TODO: PLEASE REFACTOR ME!
  annot <- annotation
  rm(annotation)


  # ---- Data Loading ----
  # Gene Expression Matrix - log2-Intensity-Values
  log_info("Loading data...")
  setwd(myFolder)
  
  # Read raw data
  log_info("Loading raw data...")
  dataset = read.table(myFile, header = FALSE, sep = "\t", dec = ".")
  d = dim(dataset)
  log_info(paste("Raw dataset dimensions:", paste0(d, collapse = ", ")))
  printdata(dataset)
  
  # Extract column headings (sample identifiers)
  log_info("Extracting sample identifiers...")
  header  = read.table(myFile, nrows = 1, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  dataset = read.table(myFile, skip = rowOffset, header = FALSE, sep = "\t", dec = ".")
  colnames(dataset) = unlist(header)
  d = dim(dataset)
  log_info(paste("Headed dataset dimensions:", paste0(d, collapse = ", ")))
  printdata(dataset)
  
  # Extract row names (gene identifiers)
  log_info("Extracting gene identifiers...")
  rownames(dataset) = dataset[,colWithID]
  dataset = dataset[,(colOffset+1):d[2]]
  header  = header[,(colOffset+1):d[2]]
  d = dim(dataset)
  log_info(paste("Final dataset dimensions:", paste0(d, collapse = ", ")))
  printdata(dataset)

  # NOTICE:
  # The rows of the final expression matrix could be 1 more than annot matrix rows in case of Agilent arrays
  # because of the presence of the NegativeControl probe (-)3xSLv1
  pitstop("Finished loading data.")
  

  # ---- Rename Samples ----
  # Tidy Sample Names According to the Experimental Design
  log_info("Checking experimental design...")
  sampleName = vector() # To declare an empty vector
  m = length(groups)
  badMsg = "Bad Experimental \'design\' vector!\n\n"
  
  # Design check
  if (length(design) != m && length(design) != d[2]) {
    stop(badMsg)
  }
  if (length(design) == m) {
    design = rep(c(1:m), design) # Convert from Compact to Full Design mode
    log_info("Compact design mode approved")
  } else if (length(design) == d[2]) {
    if (max(design) > m || min(design) < 1 || sum((design %% 1) != 0)) {
      stop(badMsg)
    } else {
      log_info("Full design mode approved")
    }
  }
  
  # Create a new vector containing tidy group names
  log_info("Making tidy group names using the experimental design...")
  for (i in 1:m) {
    index = which(design == i)
    grpName = paste0(groups[i], "_", c(1:length(index)))
    sampleName[index] = grpName
  }
  
  # Correspondences Table
  corrTable = cbind(unlist(header), sampleName) # Cast to matrix
  colnames(corrTable) = c("Original_ID", "R_Name")
  rownames(corrTable) = c(1:d[2])
  printdata(corrTable)

  if (saveOut) {
    write.table(
      corrTable,
      "Corresp Table.txt", sep = "\t",
      col.names = TRUE, row.names = TRUE, quote = FALSE
    )
    log_info(paste("'Corresp Table.txt' has been saved in", myFolder))
  }
  
  colnames(dataset) = sampleName
  log_info("Desing matrix loaded.")
  printdata(dataset)
  pitstop("Done loading the design matrix. Do the snippets of data above look ok?")


  # ---- Log Transformation ----
  # Distribution Inspection to Spot Non-Logarithmic Data
  p <- function() {
    boxplot(dataset, las = 2) # las = style of axis labels
  }
  printPlots(p, "Raw boxpot")
  p <- function() {
    plotDensities(dataset, legend = FALSE) # From limma package
  }
  printPlots(p, "Raw density")
  
  # log2-transform if intensities are in linear scale
  if (logConversion) {
    log_info("Log-transforming the intensities...")

    # Are there any negative values?
    negVal = length(which(dataset < 0))
    if (negVal > 0) {
      stop(negVal, " negative expression values detected. Cannot log-normalize.")
    }

    # Are there any zero values?
    zeroVal = length(which(dataset == 0))
    # Value imputation
    if (zeroVal > 0) {
      log_warn(zeroVal, " null values found. They will get removed.")
      # Non-zero minimum value
      imputation = min(dataset[which(dataset > 0, arr.ind = TRUE)])
      #imputation = NaN # Alternatively
      dataset[which(dataset == 0, arr.ind = TRUE)] = imputation
      log_warn(zeroVal, " null values have been imputated with", imputation)
    }
    
    log_info("Log transforming...")
    dataset = log2(dataset)
    log_info("Transformation complete.")
    printdata(dataset)
    pitstop("Done log transformation. Continue if the outputs above look ok.")
  }
  

  # ---- Normalization ----
  # After-RMA 2nd Quantile Normalization
  if (secondNorm) {
    log_info("Renormalizing data...")
    normData = normalize.quantiles(as.matrix(dataset)) # From preprocessCore package
    normData = as.data.frame(normData)
    rownames(normData) = rownames(dataset)
    colnames(normData) = colnames(dataset)
    
    # Inspect...
    printdata(normData)
    p <- function () {
      boxplot(normData, las = 2) # las = style of axis labels
    }
    printPlots(p, "Renormalized Boxplot")
    p <- function () {
      plotDensities(normData, legend = FALSE) # From limma package
    }
    printPlots(p, "Renormalized Densities")
    pitstop("Do the normalized boxplots and data look ok?")
    
    # ...if ok, then reassign
    dataset = normData
    log_info("Renormalization complete.")
  }
  

  # ---- MA-Plot & Box-Plot ----
  # Normalization Final Check with Figure Production
  p <- function() {
    boxplot(
      dataset,
      las = 2, col = myColors[design],
      main = "Expression values per sample", ylab = "log2 (intesity)"
    ) # las = style of axis labels
  }
  printPlots(p, "Final Boxplots")
  p <- function() {
    plotDensities(dataset, legend = FALSE, main = "Expression values per sample") # From limma package
  }
  printPlots(p, "Final Density")
  
  pitstop("Maybe check the plots and come back?")
  
  # MA-Plot for bias detection
  # From limma package: array 'column' vs the average of all the arrays other than that
  p <- function() {
    plotMD(dataset, column = 1)
  }
  printPlots(p, "Mean-Difference Plot")
  for (i in 1:(m-1)) { # All the possible combinations of two groups
    for (j in (i+1):m) {
      
      arr1 = rowMeans(dataset[,which(design == i)])
      arr2 = rowMeans(dataset[,which(design == j)])
      
      # From rafalib package
      p <- function() {
        maplot(
          arr1, arr2, 
          xlab = "A (Average log-expression)", ylab = "M (Expression log-ratio)",
          n = 5e4,
          curve.add = TRUE, curve.col = myColors[2], curve.lwd = 1.5, curve.n = 1e4,
          pch = 20, cex = 0.1
        )
        # As an alternative, from limma package (without trend curve)
        #plotMD(cbind(arr2, arr1),
        #       xlab = "A (Average log-expression)",
        #       ylab = "M (Expression log-ratio)")
        title(main = paste(groups[j], "vs", groups[i]))
        abline(h = 0, col = myColors[1], lty = 2) # lty = line type
        abline(h = c(1,-1), col = myColors[1])
        abline(v = thr0, col = myColors[1]) # Platform-specific log2-expression threshold
      }
      printPlots(p, paste("MA-Plot ", groups[j], "_vs_", groups[i], sep = ""))
    }
  }


  # ---- Clustering ----
  # Sample-wise Hierarchical Clustering for Batch-Effect Detection
  log_info("Starting sample-wise hierarchical clustering for batch-effect detection...")
  # Matrix Transpose t() is used because dist() computes the distances between
  # the ROWS of a matrix
  # Distance matrix (NOTE: t(dataset) is coerced to matrix)
  sampleDist = dist(t(dataset), method = "euclidean")
  hc = hclust(sampleDist, method = "ward.D")
  p <- function() {
    plot(hc) # Display Dendrogram
  }
  printPlots(p, "Dendrogram")
  # Desired number of clusters
  kNum = 6
  
  clust = cutree(hc, k = kNum) # Cut tree into kNum clusters
  clustList = list() # Create an empty list
  for (i in 1:kNum) {
    clustList[[i]] = cbind(clust[which(clust == i)]) # list of matrices
  }

  p <- function() {
    plot(hc)
    print(rect.hclust(hc, k = kNum, border = myColors[2])) # Red borders around the kNum clusters 
  }
  printPlots(p, "Dendrogram and Clusters")
  
  log_info("Finished with hierarchical clustering.")
  pitstop("Maybe check the plots and come back?")
  

  # ---- PCA ----
  # Performed on Samples for Batch-Effect Detection
  log_info("Performing PCA to detect batch sampling...")
  # Strictly enforced that rownames(metadata) == colnames(dataset)
  myMetadata = data.frame(row.names = colnames(dataset))
  # Add column 'Group' to metadata dataframe
  myMetadata$Group = rep(NA, d[2])
  for (i in 1:m) {
    myMetadata$Group[which(design == i)] = groups[i]
  }
  # Do the PCA (centering the data before performing PCA, by default)
  pcaOut = pca(dataset, metadata = myMetadata)
  log_info("Finished running PCA. Plotting results...")
  # Plot results
  p <- function() {
    print(screeplot(pcaOut))
  }
  printPlots(p, "Scree Plot")
  
  colMAtch = myColors[1:m] # Vector for color matching
  
  names(colMAtch) = groups
  p <- function() {
    suppressMessages(print(biplot(pcaOut, colby = "Group", colkey = colMAtch)))
  }
  printPlots(p, "PCA")
  p <- function() {
    suppressMessages(print(pairsplot(pcaOut, colby = "Group", colkey = colMAtch)))
  }
  printPlots(p, "PCA Pairs")

  # Possibly remove some 'batched' samples, by sample name, e.g.:
  # toBeRemoved = c("TG_1","WT_2","TGFK_1","Ab_5")
  # TODO: implement this.
  toBeRemoved = c()
  if (length(toBeRemoved) > 0) {
    for (i in 1:length(toBeRemoved)) {
      rem.Index = which(colnames(dataset) == toBeRemoved[i])
      dataset = dataset[,-rem.Index]
      design = design[-rem.Index]
      sampleName = sampleName[-rem.Index]
    }
    log_info(paste0(length(toBeRemoved), " samples have been removed."))
    d = dim(dataset)
    paste0("Sub-dataset dimensions: ", d)
    printdata(dataset)
  }
  
  pitstop("Take a look at the PCA results.")
  
  # NOTE - NOTE - NOTE - NOTE - NOTE - NOTE - NOTE - NOTE - NOTE - NOTE - NOTE - NOTE
  # You can rerun this section to get the PC representation without the batched samples
  # TODO : Wrap this step in a fun, and place a general option to redo this to 
  # make new plots. Remember to have the fun rename the plots so that they don't
  # overwrite each other.
  # TODO : In interactive mode, allow the interactive selection of which samples 
  # to discard. NOTE :: This lowers reproducibility. I suggest **just having the
  # user either set in the conf which samples to remove**, or not allow the
  # option in the first place.

  
  # ---- SD vs Mean Plot ----
  # Poisson Hypothesis Check
  log_info("Using Poisson to produce a SD vs Mean plot...")
  
  # Un-log intensity values
  log_info("Returning to linear intensities...")
  unlogged.dataset = 2^dataset
  printdata(unlogged.dataset)
  pitstop("Did the 'unlogging' mess anything up?")
  
  # Store values using matrices
  log_info("Making matrices...")
  meansArr = matrix(nrow = d[1], ncol = m+1)
  SDsArr = matrix(nrow = d[1], ncol = m+1)
  corrArr = vector()
  
  # Statistics for each group...
  log_info("Calculating groupwise statistics...")
  for (i in 1:m) {
    meansArr[,i] = rowMeans(unlogged.dataset[,which(design == i)], na.rm = TRUE) # Within-group mean
    SDsArr[,i] = apply(unlogged.dataset[,which(design == i)], 1, sd, na.rm = TRUE) # Within-group SD
    corrArr[i] = cor(meansArr[,i], SDsArr[,i]) # Mean-SD Correlation
  }
  # ...and for the whole experiment
  log_info("Calculating overall statistics...")
  meansArr[,m+1] = rowMeans(unlogged.dataset, na.rm = TRUE)
  SDsArr[,m+1] = apply(unlogged.dataset, 1, sd, na.rm = TRUE)
  corrArr[m+1] = cor(meansArr[,m+1], SDsArr[,m+1])
  meansArr[is.nan(meansArr)] = NA # To convert NaN in NA when some of the original groups are actually missing
  
  # Scatter plot
  log_info("Making plots...")
  p <- function() {
    par(mfrow = c(1, m+1)) # Optional, for sub-plotting
    X.max = max(meansArr, na.rm = TRUE)
    Y.max = max(SDsArr, na.rm = TRUE)
    for (i in 1:(m+1)) {
      
      plot(meansArr[,i], SDsArr[,i],
           xlab = "Mean", ylab = "SD",
           xlim = c(0, X.max), ylim = c(0, Y.max),
           pch = 20, cex = 0.1)
      
      if (i <= m) {
        title(main = groups[i])
      }
      else {
        title(main = "Global")
      }
      mtext(side = 3, paste("Corr =", toString(round(corrArr[i], digits = 5))))
    }
  }
  printPlots(p, "SD_vs_Mean Plot") # Save just the 'Global' one
  par(mfrow = c(1, 1)) # To reset sub-plotting
  
  
  # ---- Filtering ----
  # Eliminate Low-Intensity Probes
  log_info("Startig filtering steps...")
  # Minimum gene presence per group - Default=0.80
  kFac = opts$design$filters$min_groupwise_presence 
  grSize = vector()
  for (i in 1:m) {
    grSize[i] = sum(design == i)
  }
  kSize = ceiling(grSize * kFac) # 'ceiling' to be more stringent than 'round'
  ## NOTE :: I have NO idea why it doesen't work without the ending [1]...
  log_info(paste0("Filtering with a k Factor of ", kFac, " and k Size of ", kSize)[1])
  
  # Filtering Table
  log_info("Making filtering report ...")
  # Cast to matrix
  filTable = cbind(grSize, grSize*kFac, kSize, round(100*(kSize/grSize), 1))
  colnames(filTable) = c("Group_Size", "Min_Presence", "Rounded", "Actual %")
  rownames(filTable) = groups
  mgppg = paste0("Minimum gene presence per group = ", kFac*100, "% of the samples.")
  log_info(mgppg)
  printdata(filTable)

  if (saveOut) {
    log_info("Saving filtering report...")
    write(
      paste0("Filtering Summary File\n======================\n"),
      "Filtering Report.txt"
    )
    write(
      paste0("Presence threshold thr0 = ", thr0),
      "Filtering Report.txt",
      append = TRUE
    )
    write(mgppg, "Filtering Report.txt", append = TRUE)
    suppressWarnings(
      write.table(
        filTable, "Filtering Report.txt", sep = "\t", append = TRUE,
        col.names = TRUE, row.names = TRUE, quote = FALSE
      )
    )
    log_info(paste0("'Filtering Report.txt' has been saved in", myFolder))
  }
  pitstop("Does the filtering report look ok? If so, continue.")

  log_info("Running Filtering...")
  # Above thr0 value in at least kSize samples out of grSize (kFac*100 % - Default=80%)...
  presenceIndx = matrix(nrow = d[1], ncol = m) # Store values using matrices
  for (i in 1:m) {
    f1 = kOverA(kSize[i], thr0)
    ffun1 = filterfun(f1)
    presenceIndx[,i] = genefilter(dataset[,which(design == i)], ffun1)
  }
  retained = as.matrix(colSums(presenceIndx)) # Cast to matrix
  colnames(retained) = c("Retained_Genes")
  rownames(retained) = groups
  if (saveOut) {
    write("\n", "Filtering Report.txt", append = TRUE)
    suppressWarnings(
      write.table(
        retained, "Filtering Report.txt", append = TRUE,
        sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE
      )
    )
  }
  
  #...in at least 1 group (row-wise logical OR)
  unionSet = (rowSums(presenceIndx) > 0)
  
  filtSize = sum(unionSet)
  report = paste0(
    "\nFiltering Report:\n=================\n",
    "  ", d[1] - filtSize, " genes have been filtered out (",
    round(((d[1] - filtSize)/d[1])*1e2, digits = 2), "%)\n",
    "  ", filtSize, " retained genes out of ", d[1], " (",
    round((filtSize/d[1])*1e2, digits = 2), "%)\n\n"
  )
  cat(report)
  if (saveOut) {
    write(paste("\n", report, sep = ""), "Filtering Report.txt", append = TRUE)
  }
  
  dataset = dataset[which(unionSet),]
  log_info("Done filtering...")
  

  # ---- Contrast Definition ----
  log_info("Selecting contrasts...")
  # Reassignment to get 'filtSize' defined even in the case of no filtering
  filtSize = dim(dataset)[1]
  
  counter = 1
  myContr = vector()
  for (i in 1:(m-1)) { # All the possible combinations of two groups
    for (j in (i+1):m) {
      # Character vector specifying contrasts
      myContr[counter] = paste(groups[j], "-", groups[i], sep = "")
      counter = counter + 1
    }
  }
  
  myContr = as.matrix(myContr)
  rownames(myContr) = c(1:length(myContr))
  
  # Check with the user if the contrasts are OK
  paste0(1:length(myContr), ": ", myContr[1,]) |>
    paste0(collapse = " / ") ->
    available.contrasts
  log_info(paste0("Available contrasts: ", available.contrasts))
  selected.contrasts <- opts$design$contrasts
  log_info(paste0("Selected contrasts: ", selected.contrasts))
  pitstop("Check the contrasts before continuing.")
  # Check if the contrasts are ok
  log_info("Cheking contrasts...")
  nr.contrasts <- length(myContr[1,])
  if (length(selected.contrasts) > nr.contrasts || any(selected.contrasts > nr.contrasts)) {
    stop(paste0("Invalid contrast selection: ", selected.contrasts))
  }
  log_info("Selecting contrasts...")
  myContr = as.matrix(myContr[selected.contrasts,])
  rownames(myContr) = c(1:length(myContr))
  log_info("Contrast set.")
  pitstop("")
  
  # ---- DE by Limma ----
  if (opts$switches$limma$run_DEA) {
    log_info("Running DEA with `limma`...")
    pitstop("")
    if (opts$switches$limma$run_paired_DEA == FALSE) {
      # ---- Unpaired analysis ----
      # Differential Expression Assessment and Figure Production
      log_info("Running differential expression analysis in unpaired mode.")
      log_info("Making limma desing matrix...")
      limmaDesign = matrix(data = 0, nrow = d[2], ncol = m)
      for (i in 1:m) {
        limmaDesign[which(design == i),i] = 1
      }
      colnames(limmaDesign) = groups
      
      fit = lmFit(dataset, limmaDesign)
      
      log_info("Making contrasts matrix...")
      contrast.matrix = makeContrasts(
        contrasts = myContr,
        levels = limmaDesign)
      
      log_info("Computing contrasts...")
      fit2 = contrasts.fit(fit, contrast.matrix)
      efit2 = eBayes(fit2)
      
      # Print Results (Top-Ten genes) for all the contrasts of interest
      # This is only useful in slowmode.
      if (opts$general$slowmode) {
        for (i in 1:length(myContr)) {
          # 'print' because automatic printing is turned off in loops (and functions)
          cat("\nDEG Top-List for contrast: ", myContr[i], "\n", sep = "")
          print(topTable(efit2, coef = i, adjust.method = "BH", sort.by = "B"))
          cat("\n") # The one time I like a bit more space.
        }
      }
      pitstop("Compute all DEGs? ")

      # Compute full DEG Tables
      log_info("Computing Differential expressions...")
      progress <- txtProgressBar(min = 0, max = length(myContr))
      DEGs.limma = list() # Create an empty list
      for (i in 1:length(myContr)) {
        DEGs.limma[[i]] = topTable(
          efit2, coef = i, number = filtSize,
          adjust.method = "BH", sort.by = "B"
        ) # this is a list of Data Frames
        DEGs.limma[[i]] = appendAnnotation(DEGs.limma[[i]], annot, sort.by = "adj.P.Val")
        setTxtProgressBar(progress, i)
      }
      close(progress)
      
      # Save full DEG Tables
      if (saveOut) {
        log_info("Saving Differential expression tables...")
        progress <- txtProgressBar(min = 0, max = length(myContr))
        for (i in 1:length(myContr)) {
          degTabName = paste("Limma - DEG Table ", myContr[i], ".txt", sep = "")
          write.table(DEGs.limma[[i]], degTabName, sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
          setTxtProgressBar(progress, i)
        }
        close(progress)
        log_info(paste0("Saved tables in ", myFolder))
      }
      
      # Summary of DEGs (you can change the Log2-Fold-Change Threshold lfc...)
      results.limma = decideTests(efit2, adjust.method = "BH", p.value = 0.05, lfc = thrFC)
      p <- function() {
        vennDiagram(results.limma)
      }
      printPlots(p, "Limma Venn")
      
      # Show Hyperparameters
      d0 = efit2$df.prior           # prior degrees of freedom
      dg = mean(fit$df.residual)    # original degrees of freedom
      hyp = cbind(
        c(efit2$s2.prior,     # prior variance
        mean(fit$sigma^2),    # mean sample residual variance
        mean(efit2$s2.post),  # mean posterior residual variance
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
      
      # Save significant DEG list in Excel format with annotations
      log_info("Saving DEGs into xlsx files. Generating tables...")
      progress <- txtProgressBar(min = 0, max = length(myContr))
      all.limma.sigs = list() # Create an empty list
      for (i in 1:length(myContr)) {
        all.limma.sigs[[i]] = topTable(
          efit2, coef = i, number = filtSize,
          adjust.method = "BH", sort.by = "B",
          p.value = 0.05, lfc = thrFC
        ) # list of Data Frames
        setTxtProgressBar(progress, i)
      }
      close(progress)

      log_info("Saving tables...")
      progress <- txtProgressBar(min = 0, max = length(myContr))
      for (i in 1:length(myContr)) {
        if (saveOut & dim(all.limma.sigs[[i]])[1] > 0) {
          all.limma.sigs[[i]] = appendAnnotation(all.limma.sigs[[i]], annot, sort.by = "adj.P.Val")
          write.xlsx(
            all.limma.sigs[[i]],
            paste("Significant Genes by limma - ", myContr[i], ".xlsx", sep = ""),
            colNames = TRUE, rowNames = TRUE, sheetName = myContr[i],
            keepNA = TRUE, firstRow = TRUE, overwrite = TRUE
          ) # Freezes the first row!
        }
        setTxtProgressBar(progress, i)
      }
      close(progress)
      pitstop("Finishen running DEA.")
    } else {
      # ---- Paired Analysis ----
      # Differential Expression Assessment in Paired-Sample Design
      # (for two-group comparison only)
      log_info("Running differential expression in paired mode.")
      # Reduce dataset to the sole two groups to be compared by paired test
      log_info("Reducing dataset to only data of interest...")
      # TODO : These tests should be placed at the start of the script!
      if (length(myContr) > 1) {
        stop("Too many contrasts loaded - Paired test is for two-group comparison only!")
      }
      pg = strsplit(myContr, "-", fixed = TRUE) # List of the two groups to pair
      index = which(grepl(pg[[1]][1], sampleName) | grepl(pg[[1]][2], sampleName))
      dataset = dataset[,index]
      design = design[index]
      sampleName = sampleName[index]
      d = dim(dataset)
      log_info(paste("Sub-dataset dimensions: ", paste(d, collapse = ", ")))
      printdata(dataset)
      
      log_info("Getting patient pairings...")
      # Define pairing based on sampleName
      # Always to be checked !!
      patient.ID = factor(as.numeric(substring(sampleName, regexpr("_", sampleName) + 1)))
      #patient.ID = factor(c(1,1,2,2,3,3,4,5,6,4,5,6)) # ...or go manually
      log_info(paste0("Patient IDs are: ", paste0(patient.ID, collapse = ", ")))
      pitstop("Are the patient IDs ok?")
      cond = factor(groups[design])
      
      log_info("Making limma design formula...")
      limmaDesign = model.matrix(~patient.ID + cond)
      
      log_info("Fitting model...")
      fit = lmFit(dataset, limmaDesign)
      efit2 = eBayes(fit)
      
      # Print Results (Top-Ten genes)
      log_info("Getting results...")
      coeff = tail(colnames(limmaDesign), n = 1) # Return the last element of the vector
      log_info(paste0("DEG Top-List for contrast: ", myContr, " (design column: ", coeff, ")"))
      topTable(efit2, coef = coeff, adjust.method = "BH", sort.by = "B") # just on-Screen
      
      pitstop("")
      
      # Compute full DEG Tables
      log_info("Getting Differential expressions...")
      DEGs.limma = list() # Create an empty list (for compatibility with independent-sample test)
      DEGs.limma[[1]] = topTable(
        efit2, coef = coeff, number = filtSize,
        adjust.method = "BH", sort.by = "B"
      ) # list of Data Frames
      DEGs.limma[[1]] = appendAnnotation(DEGs.limma[[1]], annot, sort.by = "adj.P.Val")
      
      # Save full DEG Tables
      if (saveOut) {
        log_info("Saving Differential expression table...")
        degTabName = paste("Limma - DEG Table ", myContr, " - Paired.txt", sep = "")
        write.table(DEGs.limma[[1]], degTabName, sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
        cat("\n'", degTabName, "' has been saved in ", myFolder, "\n\n", sep = "")
      }
      
      # BUG:: This call probaby distorts the data to make the wrong MA plot
      # described in issue #7.
      results.limma = decideTests(efit2, adjust.method = "BH", p.value = 0.05, lfc = thrFC)
      
      # Show Hyperparameters
      d0 = efit2$df.prior           # prior degrees of freedom
      dg = mean(fit$df.residual)    # original degrees of freedom
      hyp = cbind(
        c(efit2$s2.prior,       # prior variance
        mean(fit$sigma^2),      # mean sample residual variance
        mean(efit2$s2.post),    # mean posterior residual variance
        d0, dg, d0/(d0+dg))     # Shrinkage degree
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
      # TODO : The formatting of the hyperparameters is terrible. But i don't 
      # want to lose anymore time.
      log_info(paste0("Hyperparameters:: ", get.print.str(hyp)))
      # Get rid of this print when you fix above.
      print(hyp)
      
      # Save significant DEG list in Excel format with annotations
      log_info("Saving DEGs into a xlsx file. Generating table...")
      all.limma.sigs = list() # Create an empty list
      all.limma.sigs[[1]] = topTable(
        efit2, coef = coeff, number = filtSize,
        adjust.method = "BH", sort.by = "B",
        p.value = 0.05, lfc = thrFC
      ) # list of Data Frames
      
      if (saveOut & dim(all.limma.sigs[[1]])[1] > 0) {
        log_info("Saving table...")
        all.limma.sigs[[1]] = appendAnnotation(all.limma.sigs[[1]], annot, sort.by = "adj.P.Val")
        write.xlsx(
          all.limma.sigs[[1]],
          paste0("Significant Genes by limma - ", myContr, " - Paired.xlsx"),
          colNames = TRUE, rowNames = TRUE, sheetName = myContr,
          keepNA = TRUE, firstRow = TRUE, overwrite = TRUE
        ) # Freezes the first row!
      }
      
      pitstop("Finished running DEA with limma.")
    }
    
    # ---- Limma Plot ----
    log_info("Making DEG plots...")
    # MA-Plots with significant DEGs and Volcano plots
    
    # Find Axis Limits
    log_info("Finding axis limits...")
    max.M.value = 0
    for (i in 1:length(myContr)) {
      temp = max(abs(DEGs.limma[[i]]$logFC))
      if (temp > max.M.value) {
        max.M.value = temp
      }
    }
    max.A.value = 0
    for (i in 1:length(myContr)) {
      temp = max(DEGs.limma[[i]]$AveExpr)
      if (temp > max.A.value) {
        max.A.value = temp
      }
    }
    min.A.value = Inf
    for (i in 1:length(myContr)) {
      temp = min(DEGs.limma[[i]]$AveExpr)
      if (temp < min.A.value) {
        min.A.value = temp
      }
    }
    min.P.value = 1
    for (i in 1:length(myContr)) {
      temp = min(DEGs.limma[[i]]$P.Value)
      if (temp < min.P.value) {
        min.P.value = temp
      }
    }
    
    # MA-Plot with DEGs
    log_info("Making Limma MA plots...")
    # I cannot place a progress bar here as `printPlots` logs to stdout.
    for (i in 1:length(myContr)) {
      # Mark in red/blue all the up-/down- regulated genes (+1/-1 in 'results.limma' matrix)
      p <- function() {
        plotMD(
          efit2, coef = i, status = results.limma[,i],
          values = c(1,-1), hl.col = myColors[c(2,1)],
          xlim = c(min.A.value, max.A.value), ylim = c(-max.M.value, max.M.value),
          xlab = "A (Average log-expression)",
          ylab = "M (log2-Fold-Change)"
        )
        abline(h = 0, col = myColors[1], lty = 2) # lty = line type
        abline(h = c(thrFC,-thrFC), col = myColors[1])
        abline(v = thr0, col = myColors[1]) # Platform-specific log2-expression threshold
      }
      printPlots(p, paste("MA-Plot with Limma DEGs ", myContr[i], sep = ""))
    }
    
    # Volcano Plots
    log_info("Making Limma Volcano plots...")
    for (i in 1:length(myContr)) {
      tot.DEG = sum(DEGs.limma[[i]]$adj.P.Val < 0.05) # Total number of significant DEGs (without any FC cutoff)
      high.DEG = min(c(5,tot.DEG)) # To highlight no more than 5 genes per plot
      
      # Significance Threshold
      # Find that p-value corresponding to BH-adj.p-value ~ 0.05 (or Bonferroni point when tot.DEG = 0)
      thrP = (0.05/filtSize)*(tot.DEG + 1)
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
          "\n", myContr[i], " - Threshold Report:\n",
          "--------------------------------------\n",
          "  p-value threshold  =  ", thrP, "\n",
          "  -log10(p-value)    =  ", -log10(thrP), "\n",
          "  Gene Ranking       =  ", tot.DEG, ":", tot.DEG + 1, "\n\n", sep = ""
        )
      )
      
      # Check within the DEG list
      # WARNING: DEG list has to be sorted by p-value or B statistics!
      # It should be the alpha-crossing-point for BH-adj.p-values column...
      # ...and a thrP-containing interval for un-adjusted p-values column.
      if (tot.DEG > 0) {
        print(DEGs.limma[[i]][tot.DEG:(tot.DEG + 1), c("logFC","AveExpr","t","P.Value","adj.P.Val","B")])
      } else {
        print(DEGs.limma[[i]][1, c("logFC","AveExpr","t","P.Value","adj.P.Val","B")])
      }
      # TODO: A pitstop here maybe?
      
      # Enhanced Volcano Plot
      if (getOption("append.annot")) {
        # Case-insensitive search of Gene Symbol column
        myLabels = DEGs.limma[[i]][,grep("symbol",tolower(colnames(DEGs.limma[[i]])))]
      } else {
        myLabels = rownames(DEGs.limma[[i]])
      }
      p <- function() {
        # NOTICE: When in a for loop, you have to explicitly print your
        # resulting EnhancedVolcano object
        print(
          EnhancedVolcano(
            DEGs.limma[[i]],
            x = "logFC", y = "P.Value", ylim = c(0, -log10(min.P.value)),
            pCutoff = thrP, FCcutoff = thrFC,
            pointSize = 1,
            col = c("black", "black", "black", myColors[2]),
            lab = myLabels,
            #selectLab = myLabels[1:high.DEG],
            labSize = 4,
            title = myContr[i],
            subtitle = "Limma",
            legendPosition = "none"
          )
        )
      }
      printPlots(p, paste("Volcano with Limma DEGs ", myContr[i], sep = ""))
      
      # TODO :: Implement this?
      if (FALSE) {
        # Alternative Volcano Plot From limma package
        volcanoplot(efit2, coef = i, style = "p-value",
                    highlight = high.DEG, names = rownames(dataset), hl.col = myColors[1],
                    xlim = c(-max.M.value, max.M.value), ylim = c(0, -log10(min.P.value)),
                    xlab = "log2-Fold-Change")
        title(main = myContr[i])
        abline(v = c(thrFC,-thrFC), col = myColors[2])
        abline(h = -log10(thrP), col = myColors[2])
        abline(h = -log10(0.05), col = myColors[2], lty = 2)
        printPlots(paste("Volcano with Limma DEGs ", myContr[i], sep = ""))
      }
    }
  }
  
  # ---- DE by RankProduct ----
  if (opts$switches$rankproduct$run_DEA) {
    log_info("Running Differential expression analysis with rankproduct...")
    if (opts$switches$rankproduct$run_paired_DEA == FALSE) {
      log_info("Running in paired mode.")
      # Differential Expression Assessment and Figure Production
      results.RP = matrix(data = 0, nrow = filtSize, ncol = length(myContr))
      rownames(results.RP) = rownames(dataset)
      colnames(results.RP) = myContr
      for (i in 1:length(myContr)) {
        log_info(paste0("Running analysis ", i, "of", length(myContr)))
        
        log_info("Finding control and case groups...")
        contr.groups = strsplit(myContr[i], split = "-", fixed = TRUE)
        
        case.id = which(groups == contr.groups[[1]][1])
        case.index = which(design == case.id)
        ctrl.id = which(groups == contr.groups[[1]][2])
        ctrl.index = which(design == ctrl.id)
        
        sub.dataset = dataset[,c(ctrl.index,case.index)]
        cl = c(rep(0,length(ctrl.index)), rep(1,length(case.index)))
        
        # Single Origin
        log_info("Setting origins...")
        origin = rep(1,length(cl))
        # Multiple origins
        # TODO : Implement this?
        # For-looping would require a length(myContr)-long list... to be done
        #origin = c(1,1,1,1,2,3,1,1,1,1,2,3,2,2)
        
        # invisible(capture.output()) is to suppress automatic output to console
        # WARNING: therein <- (instead of =) is mandatory for assignment!
        log_info("Running RankProduct...")
        invisible(
          capture.output(
            RP.out <- RP.advance(
              sub.dataset, cl, origin, gene.names = rownames(dataset),
              logged = TRUE, na.rm = FALSE, plot = FALSE, rand = 123
            )
          )
        )
        invisible(capture.output(plotRP(RP.out, cutoff = 0.05)))
        
        # Compute full DEG Tables (returns a list of 2 matrices, not data frames)
        log_info("Computing DEG table...")
        invisible(
          capture.output(
            DEGs.RP <- topGene(RP.out, logged = TRUE, logbase = 2, num.gene = filtSize)
          )
        )
        for (j in 1:2) {
           # Invert FC to get Case vs Ctrl and take the log2 values
          DEGs.RP[[j]][,3] = log2(1/DEGs.RP[[j]][,3])
          colnames(DEGs.RP[[j]])[3] = "Log2FC" # Correct column name
        }
        
        # Print Results (Top-Ten genes) for all the contrasts of interest
        log_info("Getting top contrasts...")
        log_info(paste0("DEG Top-List for contrast: ", myContr[i]))
        tops = rbind(DEGs.RP$Table1[1:10,], DEGs.RP$Table2[1:10,])
        # Sort by PFP (~FDR) and take just the 'absolute' Top-Ten
        tops = tops[order(tops[,4])[1:10],]
        print(tops) # just on-Screen
        pitstop("")
        
        # Save full DEG Tables
        if (saveOut) {
          log_info("Saving full DEG tables...")
          upDegTabName = paste0("RP_Up - DEG Table ", myContr[i], ".txt")
          dwnDegTabName = paste0("RP_Down - DEG Table ", myContr[i], ".txt")
          write.table(
            DEGs.RP[[1]], upDegTabName,
            sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE
          )
          write.table(
            DEGs.RP[[2]], dwnDegTabName,
            sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE
          )
          log_info(
            paste0(
              upDegTabName, "' and '", dwnDegTabName,
              "' have been saved in ", myFolder
            )
          )
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
        
        # Save significant DEG list in Excel format with annotations
        if (saveOut & (length(ups.Index) > 1 | length(dwn.Index) > 1)) {
          log_info("Saving annotated DEGs in excel format...")
          # To join two objects vertically
          all.RP.sigs = rbind(DEGs.RP$Table1[ups.Index,], DEGs.RP$Table2[dwn.Index,])
          all.RP.sigs = appendAnnotation(all.RP.sigs, annot, sort.by = "pfp")
          write.xlsx(
            all.RP.sigs, paste0("Significant Genes by RP - ", myContr[i], ".xlsx"),
            colNames = TRUE, rowNames = TRUE, sheetName = myContr[i],
            keepNA = TRUE, firstRow = TRUE, overwrite = TRUE
          ) # Freezes the first row!
        }
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
    } else {
      # ---- DE by RP - Paired ----
      # One-class Analysis of the Paired (log)-Differences
      log_info("Running in paired mode.")

      log_info("Reducing dataset to values of interest...")
      if (length(myContr) > 1) {
        stop("Too many contrasts loaded - Paired test is for two-group comparison only!\n\n")
      }
      pg = strsplit(myContr, "-", fixed = TRUE) # List of the two groups to pair
      index = which(grepl(pg[[1]][1], sampleName) | grepl(pg[[1]][2], sampleName))
      dataset = dataset[,index]
      design = design[index]
      sampleName = sampleName[index]
      d = dim(dataset)
      cat("\nSub-dataset dimensions:", d, "\n\n", sep = " ")
      printdata(dataset)
    
      log_info("Calculating results matrix...")
      results.RP = matrix(data = 0, nrow = filtSize, ncol = 1)
      rownames(results.RP) = rownames(dataset)
      colnames(results.RP) = myContr
      
      log_info("Finding pairs based on sample names...")
      case.index = which(grepl(pg[[1]][1], sampleName))
      ctrl.index = which(grepl(pg[[1]][2], sampleName))
      # ...or go manually
      # TODO : Implement this?
      #case.index = c(1,2,3,5,7)
      #ctrl.index = c(8,9,10,4,6)
      
      # Check the re-ordering before pairing
      cat("Pairings:/n")
      print(cbind(sampleName[case.index], sampleName[ctrl.index]))
      pitstop("Do the pairings look ok?")
      
      # Pair Samples
      log_info("Pairing samples...")
      paired.dataset = dataset[,case.index] - dataset[,ctrl.index]
      cl = rep(1,dim(paired.dataset)[2])
      
      # Single Origin
      log_info("Setting origins...")
      origin = rep(1, length(paired.dataset))
      # Multiple origins
      # TODO : Implement this?
      #origin = c(1,1,1,1,2,2)
      
      # invisible(capture.output()) is to suppress automatic output to console
      # WARNING: therein <- (instead of =) is mandatory for assignment!
      log_info("Running rankproduct...")
      invisible(
        capture.output(
          RP.out <- RP.advance(
            paired.dataset, cl, origin, gene.names = rownames(dataset),
            logged = TRUE, na.rm = FALSE, plot = FALSE, rand = 123
          )
        )
      )
      p <- function() {
        invisible(capture.output(plotRP(RP.out, cutoff = 0.05)))
      }
      printPlots(p, "Rankproduct plot")
      
      # Compute full DEG Tables (it returns a list of 2 matrices, not data frames)
      log_info("Computing DEG tables...")
      invisible(
        capture.output(
          DEGs.RP <- topGene(RP.out, logged = TRUE, logbase = 2, num.gene = filtSize)
        )
      )
      for (j in 1:2) {
        DEGs.RP[[j]][,3] = log2(DEGs.RP[[j]][,3]) # Take the log2 values
        colnames(DEGs.RP[[j]])[3] = "Log2FC" # Correct column name
      }
      
      # Print Results (Top-Ten genes) for the contrast of interest
      cat("DEG Top-List for paired contrast: ", myContr, "\n", sep = "")
      tops = rbind(DEGs.RP$Table1[1:10,], DEGs.RP$Table2[1:10,])
      # Sort by PFP (~FDR) and take just the 'absolute' Top-Ten
      tops = tops[order(tops[,4])[1:10],]
      print(tops) # just on-Screen
      pitstop("")
      
      # Save full DEG Tables
      if (saveOut) {
        log_info("Saving full DEG tables...")
        upDegTabName = paste("RP_Up - DEG Table ", myContr, ".txt", sep = "")
        dwnDegTabName = paste("RP_Down - DEG Table ", myContr, ".txt", sep = "")
        write.table(DEGs.RP[[1]], upDegTabName, sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
        write.table(DEGs.RP[[2]], dwnDegTabName, sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
        log_info(
          paste0(
            upDegTabName, "' and '", dwnDegTabName, "' have been saved in ", myFolder
          )
        )
      }
      
      # Fetch indexes of significant DEGs
      log_info("Finding DEGs of interest...")
      dwn.Index = which(DEGs.RP$Table1[,4] < 0.05 & DEGs.RP$Table1[,3] < -thrFC)
      ups.Index = which(DEGs.RP$Table2[,4] < 0.05 & DEGs.RP$Table2[,3] > thrFC)
      log_info(
        paste(
          "Found", length(ups.Index), "upregulated and",
          length(dwn.Index), "downregulated genes."
        )
      )
      
      # Here 'results.RP' is filled in analogy to 'results.limma':
      # (-1,0,+1) for (Down, NotSig, Up)-regulated genes
      # Accessing matrix rows by row name: matrix[rowname, col_i]
      results.RP[rownames(DEGs.RP$Table1)[dwn.Index],1] = -1
      results.RP[rownames(DEGs.RP$Table2)[ups.Index],1] = 1
      
      # Save significant DEG list in Excel format with annotations
      if (saveOut & (length(ups.Index) > 1 | length(dwn.Index) > 1)) {
        log_info("Saving annotated DEGs in excel format...")
        # To join two objects vertically use `rbind`
        all.RP.sigs = rbind(DEGs.RP$Table1[dwn.Index,], DEGs.RP$Table2[ups.Index,])
        all.RP.sigs = appendAnnotation(all.RP.sigs, annot, sort.by = "pfp")
        write.xlsx(
          all.RP.sigs,
          paste("Significant Genes by RP - ", myContr, ".xlsx", sep = ""),
          colNames = TRUE, rowNames = TRUE, sheetName = myContr,
          keepNA = TRUE, firstRow = TRUE
        ) # Freezes the first row!
      }
      
      log_info("Calculating summary statistics...")
      summary.RP = rbind(
        colSums(results.RP == -1), # Cast to matrix
        colSums(results.RP == 0),
        colSums(results.RP == 1)
      )
      rownames(summary.RP) = c("Down", "NotSig", "Up")
      print(summary.RP)
      
      log_info("Finished running DEA.")
      pitstop("")
    }
  }
  
  if (opts$switches$limma$run_DEA & opts$switches$rankproduct$run_DEA) {
    log_info("Making comparison plots between limma and rankproduct...")
    # DEG Comparison -------------------------------------------------------------------------------------------------------
    # Venn diagrams of DEGs from Limma and RP methods
    
    # To suppress 'venn.diagram()' logging
    futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
    
    # Plot Venn diagrams
    log_info("Plotting Venn diagrams...")
    for (i in 1:length(myContr)) {
      for (j in c(1, -1)) {

        DEG.id.limma = rownames(results.limma)[which(results.limma[,i] == j)]
        DEG.id.RP = rownames(results.RP)[which(results.RP[,i] == j)]
        
        if (j == 1) {
          venn.sub = "UP-regulated DEGs"
        } else {
          venn.sub = "DOWN-regulated DEGs"
        }
        
        if (length(DEG.id.limma) == 0 & length(DEG.id.RP) == 0) {
          cat("\nBoth the sets are empty in the contrast: ", myContr[i], " (", venn.sub, ")\n\n", sep = "")
          next # Skip the current iteration of the for loop without terminating it
        }
        
        venn.plot = venn.diagram(
          x = list(DEG.id.limma, DEG.id.RP),
          filename = NULL, # to print just on screen
          force.unique = TRUE,
          main = myContr[i], main.cex = 2, main.fontface = "bold", main.fontfamily = "sans", # Title
          sub = venn.sub, sub.fontfamily = "sans", # Subtitle
          lwd = 2, lty = "blank", fill = myColors[1:2], # circles
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
            "Comparison Venn ", myContr[i], "_",
            strsplit(venn.sub, split = "-")[[1]][1]
          )
        )
      }
    }
  }
  
  log_info("GATTACA finished")
}
