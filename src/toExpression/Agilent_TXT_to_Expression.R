# Header Info ----------------------------------------------------------------------------------------------------------
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





# Load Packages --------------------------------------------------------------------------------------------------------

library(limma)      # Empirical Bayes Method for Differential Expression (here used for low-level analysis only)
#library(openxlsx)   # Read, Write, and Edit .xlsx (Excel) Files





# Variable Definition --------------------------------------------------------------------------------------------------

system.root = "D:\\UniUPO Drive\\" # SilverLife @ Home
system.root = "D:\\Drive UniUPO\\" # SkyLake2   @ DBIOS

k = 4 # Experiment to analyze
endName = c("Exp 2017-10-27", # 1
            "Exp 2017-11-04", # 2
            "Exp 2018-01-12", # 3
            "Exp 2018-01-15", # 4
            "Exp 2018-05-18", # 5
            "Exp 2018-06-07", # 6
            "Exp 2018-10-12") # 7

save.PNG.plot = TRUE
save.PDF.plot = FALSE





# Function Definition --------------------------------------------------------------------------------------------------
# User-Defined Functions

# Save a graphical output to 'Output Figures' sub-directory
#   figureName  = output file name (without extension)
#   PNG.bool    = T to print the currently displayed figure in PNG format
#   PDF.bool    = T to print the currently displayed figure in PDF format
printPlots = function(figureName, PNG.bool = save.PNG.plot, PDF.bool = save.PDF.plot)
{
  figSubFolder = "Output Figures"
  fullName = file.path(figSubFolder, figureName, fsep = .Platform$file.sep) # OS-independent path separator
  
  if (!file.exists(figSubFolder) & (PNG.bool | PDF.bool)) {
    dir.create(figSubFolder)
    cat("\nNew folder '", figSubFolder, "' has been created in the PWD...\n\n", sep = "")
  }
  if (PNG.bool) { # invisible(capture.output()) to suppress automatic output to console
    invisible(capture.output(
      dev.print(device = png, filename = paste(fullName, ".png", sep = ""), width = 820, height = 600)))
  }
  if (PDF.bool) {
    invisible(capture.output(
      dev.print(device = pdf, paste(fullName, ".pdf", sep = ""))))
  }
}





# Load TXt files and read raw data -------------------------------------------------------------------------------------

myFolder = paste(system.root, "WORKS\\202x - Article - Colon\\Data\\1 - Raw Data\\", endName[k], sep = "")
setwd(myFolder)

# Load target files - 'readTargets()' will by default look for the 'FileName' column in the specified file
targets = readTargets("Targets.txt", row.names = "SampleNumber")
targets
m = dim(targets)[1]

# Convert the data to an EListRaw object (data object for single-channel arrays)
# Specify green.only=TRUE for one-channel data and retain information about background via gIsWellAboveBG
raw = read.maimages(targets, source = 'agilent.median', green.only = TRUE, other.columns = 'gIsWellAboveBG')

# Generate QC plots
boxplot(log2(as.matrix(raw)), names = c(1:m), las = 1, # las = style of axis labels
        main = "Raw Expression", xlab = "Sample Index", ylab = "log2(intensity)") # Raw data to be logged
printPlots("1 - BoxPlot - Raw")
plotDensities(log2(as.matrix(raw)), legend = FALSE)
printPlots("1 - Densities - Raw")
# MA-plots from limma package (without trend curve)
samp = c(11,12) # One sample vs one sample (by sample indexes)
plotMD(log2(as.matrix(raw)[,samp]), main = paste("Raw Expression: ", samp[1], " vs ", samp[2], sep = ""),
       xlab = "A (Average log-expression)", ylab = "M (Expression log-ratio)")
printPlots(paste("1 - MAPlot - Raw", samp[1], "vs", samp[2], sep = ""))





# Background subtraction and interarray normalization ------------------------------------------------------------------

# Subtract the background signal
ofst = 50 # Offset to shrunk towards zero log-ratios at the lower intensities (see 'anti-fanning' in MA-Plot)
raw_BGcorrected = backgroundCorrect(raw, method = "normexp", offset = ofst)

# Generate QC plots
boxplot(log2(as.matrix(raw_BGcorrected)), names = c(1:m), las = 1, # las = style of axis labels
        main = "Unnormalized Expression", xlab = "Sample Index", ylab = "log2(intensity)") # Data still to be logged
printPlots("2 - BoxPlot - BGcorrected")
plotDensities(log2(as.matrix(raw_BGcorrected)), legend = FALSE)
printPlots("2 - Densities - BGcorrected")
samp = c(11,12) # One sample vs one sample (by sample indexes)
plotMD(log2(as.matrix(raw_BGcorrected)[,samp]),
       main = paste("Unnormalized Expression: ", samp[1], " vs ", samp[2], sep = ""),
       xlab = "A (Average log-expression)", ylab = "M (Expression log-ratio)")
printPlots(paste("2 - MAPlot - BGcorrected", samp[1], "vs", samp[2], sep = ""))

# Quantile-normalize and log-transform the data
raw_BGandNormalized = normalizeBetweenArrays(raw_BGcorrected, method = "quantile")

# Generate QC plots
boxplot(as.matrix(raw_BGandNormalized), names = c(1:m), las = 1, # las = style of axis labels
        main = "Normalized Expression", xlab = "Sample Index", ylab = "log2(intensity)") # Already logged data
printPlots("3 - BoxPlot - BGandNormalized")
plotDensities(as.matrix(raw_BGandNormalized), legend = FALSE)
printPlots("3 - Densities - BGandNormalized")
samp = c(11,12) # One sample vs one sample (by sample indexes)
plotMD(as.matrix(raw_BGandNormalized)[,samp],
       main = paste("Normalized Expression: ", samp[1], " vs ", samp[2], sep = ""),
       xlab = "A (Average log-expression)", ylab = "M (Expression log-ratio)")
printPlots(paste("3 - MAPlot - BGandNormalized", samp[1], "vs", samp[2], sep = ""))

# Extract the expression matrix from EListRaw object
expressionMatrix = raw_BGandNormalized$E
d = dim(expressionMatrix)
cat("\nDataset dimensions:", d, "\n\n", sep = " ")
# Create a new vector containing tidy group names
if(sum(colnames(targets) == "BarCode") == 1) {
  grpName = paste("Sample_", targets$BarCode, ".", c(1:m), sep = "")
} else {
  grpName = paste("Sample_", c(1:m), sep = "")
}
colnames(expressionMatrix) = grpName
expressionMatrix[1:10,]





# Optional - Remove invalid probes -------------------------------------------------------------------------------------

# Filter out control probes (and those that fail the 'gIsWellAboveBG' check)
ctrlProbes = abs(raw_BGandNormalized$genes$ControlType) == 1 # abs() to include also the NegativeControl probe (-)3xSLv1
cat("\n", sum(ctrlProbes), " control probes detected\n\n", sep = "")
#as.data.frame(raw_BGandNormalized)[ctrlProbes, 7] # Get names of Control Probes detected
isExpr = rep(TRUE, dim(raw_BGandNormalized)[1]) # No second filter - Default
if (FALSE) {
  # Well above BG in at least 50% of the samples
  isExpr = rowSums(raw_BGandNormalized$other$gIsWellAboveBG > 0) >= dim(targets)[1]/2
  sum(isExpr)
}
raw_BGandNormalized_filt = raw_BGandNormalized[!ctrlProbes & isExpr, ]

# Replace value of replicate probes with their mean - Probe_ID is used to identify the replicates
raw_BGandNormalized_filt_mean = avereps(raw_BGandNormalized_filt,  ID = raw_BGandNormalized_filt$genes$ProbeName)

# Generate QC plots
pre = dim(raw_BGandNormalized_filt)
post = dim(raw_BGandNormalized_filt_mean)
cat("\n", (pre - post)[1], " replicate probes have been averaged\n\n", sep = "")
cat("\nFinal dataset dimensions:", post, "\n\n", sep = " ")
samp = c(11,12) # One sample vs one sample (by sample indexes)
plotMD(as.matrix(raw_BGandNormalized_filt_mean)[,samp],
       main = paste("Normalized Expression - Filtered: ", samp[1], " vs ", samp[2], sep = ""),
       xlab = "A (Average log-expression)", ylab = "M (Expression log-ratio)")
abline(h = 0, col = "cornflowerblue", lty = 2) # lty = line type
abline(h = c(1,-1), col = "red3")
abline(v = 7, col = "cornflowerblue") # Platform-specific log2-expression threshold
printPlots(paste("4 - MAPlot - Final", samp[1], "vs", samp[2], sep = ""))

# Extract the expression matrix from EListRaw object
expressionMatrix = raw_BGandNormalized_filt_mean$E
colnames(expressionMatrix) = grpName
expressionMatrix[1:10,]





# Adding Gene Annotation -----------------------------------------------------------------------------------------------
# TO BE DONE





# Save expression matrix -----------------------------------------------------------------------------------------------

# Fix headings for subsequent import
expressionMatrix = as.data.frame(expressionMatrix) # Cast to Data Frame
expressionMatrix = cbind(rownames(expressionMatrix), expressionMatrix) # Add Probe_IDs as the first column
colnames(expressionMatrix)[1] = "ProbeID"

# Save the matrix without annotation
write.table(expressionMatrix, file = paste("Quantile-norm_logExpr - ", endName[k] ,".txt", sep = ""),
            sep = "\t", col.names = TRUE, row.names = FALSE)
cat("\nExpression matrix (without annotations) has been saved in\n", myFolder, "\n\n", sep = "")

# Save the matrix with annotations
#write.table(joined, file = "RMA-normalized_logExpression_U133Plus2_Annotated.txt",
#            sep = "\t", col.names = TRUE, row.names = FALSE)
#cat("\nExpression matrix (with annotations) has been saved in\n", myFolder, "\n\n", sep = "")

