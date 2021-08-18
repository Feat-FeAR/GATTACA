# Header Info ----------------------------------------------------------------------------------------------------------
#
# GATTACA v2 - General Algorithm for The Transcriptional Analysis by one-Channel Arrays
#
# a FeAR R-script - 04-Aug-2021
#
# Pipeline for one-Color (HD) Microarrays
# Data are supposed to be already background-subtracted, log2-transformed, and interarray-normalized
#

thisFile <- function() {
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

## Note for devs: this will only work on r-studio
if (!interactive()) {
  here <- head(strsplit(thisFile(), "/")[[1]], -1)
  here <- paste0(here, collapse = "/")
} else {
  if (!require("rstudioapi")) {
    install.packages("rstudioapi")
  }
  here <- dirname(rstudioapi::getSourceEditorContext()$path)
}


tryCatch(
  {
    cat("Attempting to activate the existing `renv` environment...")
    source(paste0(here, "/renv/activate.R"))
    cat("...OK\n")
  },
  error = function(err) {
    stop("No `renv` project found. Did you run the installation script?")
  }
)

# * Package Loading ----------------------------------------------------------------------------------------------------
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

source(paste0(here, "/STALKER_Functions.R"))   # Collection of custom functions


# Prepare Annotations --------------------------------------------------------------------------------------------------
# Choose your platform!

system.root = "E:\\UniTo Drive\\"            # SilverLife @ Home
system.root = "D:\\UniTo Drive\\"            # SkyLake @ DBIOS
system.root = "C:\\Users\\aleph\\OneDrive\\" # SurfaceFeAR

# Affymetrix Human Genome U133 Set (A)
      # Remote annotation database (chip hgu133a) - [2475 missing GeneSymbols]
      library(hgu133a.db)
      #hgu133a() # List of available annotations
      annot = data.frame(Accession   = sapply(contents(hgu133aACCNUM), paste, collapse = ", "),
                         GeneSymbol  = sapply(contents(hgu133aSYMBOL), paste, collapse = ", "),
                         Description = sapply(contents(hgu133aGENENAME), paste, collapse = ", "))
      miss = sum(annot[,2] == "NA")
      
      # Local annotation database - [1228 missing GeneSymbols]
      setwd(paste(system.root, "Coding\\R scripts\\Annotations\\HG-U133A-na36-annot-csv", sep = ""))
      annot = read.xlsx("HG-U133A.na36.annot.xlsx", colNames = TRUE, rowNames = TRUE, sep.names = "_") # As Data Frame
      #as.matrix(colnames(annot)) # List of available annotations
      annot = annot[,c("Representative_Public_ID", "Gene_Symbol", "Gene_Title")]
      miss = sum(annot[,2] == "---")

# Affymetrix Human Genome U133 Set (B)
      # Remote annotation database (chip hgu133b) - [7176 missing GeneSymbols]
      library(hgu133b.db)
      #hgu133b() # List of available annotations
      annot = data.frame(Accession   = sapply(contents(hgu133bACCNUM), paste, collapse = ", "),
                         GeneSymbol  = sapply(contents(hgu133bSYMBOL), paste, collapse = ", "),
                         Description = sapply(contents(hgu133bGENENAME), paste, collapse = ", "))
      miss = sum(annot[,2] == "NA")
      
      # Local annotation database - [5875 missing GeneSymbols]
      setwd(paste(system.root, "Coding\\R scripts\\Annotations\\HG-U133B-na36-annot-csv", sep = ""))
      annot = read.xlsx("HG-U133B.na36.annot.xlsx", colNames = TRUE, rowNames = TRUE, sep.names = "_") # As Data Frame
      #as.matrix(colnames(annot)) # List of available annotations
      annot = annot[,c("Representative_Public_ID", "Gene_Symbol", "Gene_Title")]
      miss = sum(annot[,2] == "---")

# Affymetrix Human Genome HG-U133 Plus 2.0 Array
      # Remote annotation database (chip hgu133plus2) - [12,770 missing GeneSymbols]
      library(hgu133plus2.db)
      #hgu133plus2() # List of available annotations
      annot = data.frame(Accession   = sapply(contents(hgu133plus2ACCNUM), paste, collapse = ", "),
                         GeneSymbol  = sapply(contents(hgu133plus2SYMBOL), paste, collapse = ", "),
                         Description = sapply(contents(hgu133plus2GENENAME), paste, collapse = ", "))
      miss = sum(annot[,2] == "NA")
      
      # Local annotation database - [9619 missing GeneSymbols]
      setwd(paste(system.root, "Coding\\R scripts\\Annotations\\HG-U133_Plus_2-na36-annot-csv", sep = ""))
      annot = read.xlsx("HG-U133_Plus_2.na36.annot.xlsx", colNames = TRUE, rowNames = TRUE, sep.names = "_") # As Data Frame
      #as.matrix(colnames(annot)) # List of available annotations
      annot = annot[,c("Representative_Public_ID", "Gene_Symbol", "Gene_Title")]
      miss = sum(annot[,2] == "---")

# Agilent-028005 SurePrint G3 Mouse GE 8x60K
      # Local annotation database - [21865 missing GeneSymbols]
      setwd(paste(system.root, "Coding\\R scripts\\Annotations\\Agilent-028005_SurePrintG3MouseGE8x60K", sep = ""))
      annot = read.xlsx("028005_D_AA_20181026.xlsx", colNames = TRUE, rowNames = TRUE, sep.names = "_") # As Data Frame
      #as.matrix(colnames(annot)) # List of available annotations
      annot = annot[,c("GeneSymbol", "GeneName", "Description")]
      miss = sum(is.na(annot[,1]))

# Agilent-026652 Whole Human Genome Microarray 4x44K v2
      # Remote annotation database (chip HsAgilentDesign026652) - [6076 missing GeneSymbols]
      library(HsAgilentDesign026652.db)
      #HsAgilentDesign026652() # List of available annotations
      annot = data.frame(Accession   = sapply(contents(HsAgilentDesign026652ACCNUM), paste, collapse = ", "),
                         GeneSymbol  = sapply(contents(HsAgilentDesign026652SYMBOL), paste, collapse = ", "),
                         Description = sapply(contents(HsAgilentDesign026652GENENAME), paste, collapse = ", "))
      miss = sum(annot[,2] == "NA")
      
      # Local annotation database - [7035 missing GeneSymbols]
      setwd(paste(system.root, "Coding\\R scripts\\Annotations\\Agilent-026652_WholeHumanGenome4x44Kv2", sep = ""))
      annot = read.xlsx("GPL13497_noParents.an.xlsx", colNames = TRUE, rowNames = TRUE, sep.names = "_") # As Data Frame
      as.matrix(colnames(annot)) # List of available annotations
      annot = annot[,c("GeneSymbols", "GeneNames")]
      miss = sum(is.na(annot[,1]))

cat("\n", miss, " unannotated genes (", round(miss/dim(annot)[1]*1e2, digits = 2), " %) \n\n", sep = "")





# * Variable Definition ------------------------------------------------------------------------------------------------
# User-Defined Experiment-Specific Variables

myFolder = paste(system.root, "WORKS\\0010 - Ongoing\\target_path", sep = "")
myFile = "Expression_Matrix_FileName.txt"

rowOffset = 1 # Row offset, including the header (How many rows above the first number?)
colWithID = 1 # Column containing (unique) gene identifiers
colOffset = 2 # Column offset, including row names (How many columns before the first number?)

groups = c("WT","Ab","AbFK","TG","TGFK") # Experimental Design - Group Names (start with control condition)
design = c(5,5,5,5,4) # Experimental Design - Ordered Group Size (compact Design Mode)

# ...or build a custom Experimental Design vector (Full Design Mode)
# Use ALL natural numbers from 1 to m=length(groups) - According to 'groups' ordering
design = c(3,3,1,1,1,1,1,4,4,4,5,5,2,2,2,2,2,3,3,3,4,4,5,5)

myColors = c("cornflowerblue","firebrick3","olivedrab3","darkgoldenrod1","purple","magenta3")
if (length(myColors) < length(groups)) {
  cat("\n")
  stop("Too few colors in \'myColors\' vector!\n\n")
}

# Filter: log2-expression conventional threshold (to be checked from case to case)
# Agilent
    # Check the un-hybridized (-)3xSLv1 NegativeControl probe as a plausible value:
    # e.g. mean(as.numeric(dataset["(-)3xSLv1",])))
    thr0 = 6
# Affymetrix
    thr0 = 4

# Fold Change Threshold. Usually, |log2FC|>0.5 OR |log2FC|>1
thrFC = 0.5

# Flags for script-wide IF
saveOut = TRUE
logConversion = FALSE
secondNorm = FALSE

# Global options suitable for STALKER_Functions
# Erase them by typing: options(scriptName = NULL, ...)
options(scriptName = "GATTACA",
        save.PNG.plot = TRUE,
        save.PDF.plot = FALSE,
        append.annot = TRUE)

if (getOption("append.annot")) {
  dim.annot = dim(annot) # If 'annot' is not defined an error message will be displayed
  cat("\nA ", dim.annot[1], " x ", dim.annot[2], " annotation dataframe has been loaded\n\n", sep = "")
} else {
  annot = NULL
  cat("\nNo annotation loaded\n\n", sep = "")
}





# * Data Loading -------------------------------------------------------------------------------------------------------
# Gene Expression Matrix - log2-Intensity-Values

setwd(myFolder)

# Read raw data
dataset = read.table(myFile, header = FALSE, sep = "\t", dec = ".")
d = dim(dataset)
cat("\nRaw dataset dimensions:", d, "\n\n", sep = " ")
dataset[1:10,1:5]

# Extract column headings (sample identifiers)
header  = read.table(myFile, nrows = 1, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
dataset = read.table(myFile, skip = rowOffset, header = FALSE, sep = "\t", dec = ".")
colnames(dataset) = unlist(header)
d = dim(dataset)
cat("\nHeaded dataset dimensions:", d, "\n\n", sep = " ")
dataset[1:10,1:5]

# Extract row names (gene identifiers)
rownames(dataset) = dataset[,colWithID]
dataset = dataset[,(colOffset+1):d[2]]
header  = header[,(colOffset+1):d[2]]
d = dim(dataset)
cat("\nFinal dataset dimensions:", d, "\n\n", sep = " ")
dataset[1:10,1:5]

# NOTICE:
# The rows of the final expression matrix could be 1 more than annot matrix rows in case of Agilent arrays
# because of the presence of the NegativeControl probe (-)3xSLv1





# * Rename Samples  ----------------------------------------------------------------------------------------------------
# Tidy Sample Names According to the Experimental Design

sampleName = vector() # To declare an empty vector
m = length(groups)
badMsg = "Bad Experimental \'design\' vector!\n\n"

# Design check
if (length(design) != m && length(design) != d[2]) {
  cat("\n")
  stop(badMsg)
}
if (length(design) == m) {
  design = rep(c(1:m), design) # Convert from Compact to Full Design mode
  cat("\nCompact \'design\' mode approved\n\n")
} else if (length(design) == d[2]) {
  if (max(design) > m || min(design) < 1 || sum((design %% 1) != 0)) { # %% is R modulus operator (for integer checking)
    cat("\n")
    stop(badMsg)
  } else {
    cat("\nFull \'design\' mode approved\n\n")
  }
}

# Create a new vector containing tidy group names
for (i in 1:m) {
  index = which(design == i)
  grpName = paste(groups[i], "_", c(1:length(index)), sep = "") # Concatenate Strings with no separator
  sampleName[index] = grpName
}

# Correspondences Table
corrTable = cbind(unlist(header), sampleName) # Cast to matrix
colnames(corrTable) = c("Original_ID", "R_Name")
rownames(corrTable) = c(1:d[2])
corrTable
if (saveOut) {
  write.table(corrTable, "Corresp Table.txt", sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
  cat("\n'Corresp Table.txt' has been saved in", myFolder, "\n\n")
}

colnames(dataset) = sampleName
dataset[1:10,1:5]





# Log Transformation ---------------------------------------------------------------------------------------------------
# Distribution Inspection to Spot Non-Logarithmic Data

boxplot(dataset, las = 2) # las = style of axis labels
plotDensities(dataset, legend = FALSE) # From limma package

# log2-transform iff intensities are in linear scale
if (logConversion) {
  
  # Are there any negative values?
  negVal = length(which(dataset < 0))
  if (negVal > 0) {
    cat("\n")
    stop(negVal, " negative expression values detected!\n\n")
  }
  
  # Are there any zero values?
  zeroVal = length(which(dataset == 0))
  cat("\n", zeroVal, "null values found\n\n")
  
  # Value imputation
  if (zeroVal > 0) {
    # Non-zero minimum value
    imputation = min(dataset[which(dataset > 0, arr.ind = TRUE)])
    #imputation = NaN # Alternatively
    
    dataset[which(dataset == 0, arr.ind = TRUE)] = imputation
    cat("\n", zeroVal, "null values have been imputated with", imputation, "\n\n")
  }
  
  length(which(dataset == 0))
  min(dataset)
  
  dataset = log2(dataset)
}





# Normalization --------------------------------------------------------------------------------------------------------
# After-RMA 2nd Quantile Normalization

if (secondNorm) {
  normData = normalize.quantiles(as.matrix(dataset)) # From preprocessCore package
  normData = as.data.frame(normData)
  rownames(normData) = rownames(dataset)
  colnames(normData) = colnames(dataset)
  
  # Inspect...
  normData[1:10,1:5]
  boxplot(normData, las = 2) # las = style of axis labels
  plotDensities(normData, legend = FALSE) # From limma package
  
  # ...if ok, then reassign
  dataset = normData
}





# MA-Plot & Box-Plot ---------------------------------------------------------------------------------------------------
# Normalization Final Check with Figure Production

boxplot(dataset, las = 2, col = myColors[design],
        main = "Expression values per sample", ylab = "log2 (intesity)") # las = style of axis labels
printPlots("1 - Boxplots")
plotDensities(dataset, legend = FALSE, main = "Expression values per sample") # From limma package
printPlots("2 - Density")

# MA-Plot for bias detection
plotMD(dataset, column = 1) # From limma package: array 'column' vs the average of all the arrays other than that
for (i in 1:(m-1)) { # All the possible combinations of two groups
  for (j in (i+1):m) {
    
    arr1 = rowMeans(dataset[,which(design == i)])
    arr2 = rowMeans(dataset[,which(design == j)])
    
    # From rafalib package
    maplot(arr1, arr2, 
           xlab = "A (Average log-expression)",
           ylab = "M (Expression log-ratio)",
           n = 5e4,
           curve.add = TRUE, curve.col = myColors[2], curve.lwd = 1.5, curve.n = 1e4,
           pch = 20, cex = 0.1)
    
    # As an alternative, from limma package (without trend curve)
    #plotMD(cbind(arr2, arr1),
    #       xlab = "A (Average log-expression)",
    #       ylab = "M (Expression log-ratio)")
    
    title(main = paste(groups[j], "vs", groups[i]))
    abline(h = 0, col = myColors[1], lty = 2) # lty = line type
    abline(h = c(1,-1), col = myColors[1])
    abline(v = thr0, col = myColors[1]) # Platform-specific log2-expression threshold
    printPlots(paste("3 - MA-Plot ", groups[j], "_vs_", groups[i], sep = ""))
  }
}





# Clustering -----------------------------------------------------------------------------------------------------------
# Sample-wise Hierarchical Clustering for Batch-Effect Detection

# Matrix Transpose t() is used because dist() computes the distances between the ROWS of a matrix
sampleDist = dist(t(dataset), method = "euclidean") # Distance matrix (NOTE: t(dataset) is coerced to matrix)
hc = hclust(sampleDist, method = "ward.D")
plot(hc) # Display Dendrogram
printPlots("4 - Dendrogram")
# Desired number of clusters
kNum = 6

clust = cutree(hc, k = kNum) # Cut tree into kNum clusters
clustList = list() # Create an empty list
for (i in 1:kNum) {
  clustList[[i]] = cbind(clust[which(clust == i)]) # list of matrices
}
clustList

rect.hclust(hc, k = kNum, border = myColors[2]) # Red borders around the kNum clusters 
printPlots("4 - Dendrogram and Clusters")





# PCA ------------------------------------------------------------------------------------------------------------------
# Performed on Samples for Batch-Effect Detection

myMetadata = data.frame(row.names = colnames(dataset)) # Strictly enforced that rownames(metadata) == colnames(dataset)

# Add column 'Group' to metadata dataframe ($ allows accessing/creating dataframe columns by name)
myMetadata$Group = rep(NA, d[2])
for (i in 1:m) {
  myMetadata$Group[which(design == i)] = groups[i]
}

# Do the PCA (centering the data before performing PCA, by default)
pcaOut = pca(dataset, metadata = myMetadata)

# Plot results
screeplot(pcaOut)
printPlots("5 - Scree Plot")

colMAtch = myColors[1:m] # Vector for color matching
names(colMAtch) = groups
biplot(pcaOut, colby = "Group", colkey = colMAtch)
printPlots("6 - PCA")
pairsplot(pcaOut, colby = "Group", colkey = colMAtch)
printPlots("6 - PCA Pairs")

# Possibly remove some 'batched' samples, by sample name, e.g.: toBeRemoved = c("TG_1","WT_2","TGFK_1","Ab_5")
toBeRemoved = c()
if (length(toBeRemoved) > 0) {
  for (i in 1:length(toBeRemoved)) {
    rem.Index = which(colnames(dataset) == toBeRemoved[i])
    dataset = dataset[,-rem.Index]
    design = design[-rem.Index]
    sampleName = sampleName[-rem.Index]
  }
  cat("\n", length(toBeRemoved), " samples have been removed\n", sep = "")
  d = dim(dataset)
  cat("\nSub-dataset dimensions:", d, "\n\n", sep = " ")
  dataset[1:2,]
}

# NOTE - NOTE - NOTE - NOTE - NOTE - NOTE - NOTE - NOTE - NOTE - NOTE - NOTE - NOTE
# You can rerun this section to get the PC representation without the batched samples





# SD vs Mean Plot ------------------------------------------------------------------------------------------------------
# Poisson Hypothesis Check

# Un-log intensity values
unlogged.dataset = 2^dataset
unlogged.dataset[1:10,1:5]

# Store values using matrices
meansArr = matrix(nrow = d[1], ncol = m+1)
SDsArr = matrix(nrow = d[1], ncol = m+1)
corrArr = vector()

# Statistics for each group...
for (i in 1:m) {
  meansArr[,i] = rowMeans(unlogged.dataset[,which(design == i)], na.rm = TRUE) # Within-group mean
  SDsArr[,i] = apply(unlogged.dataset[,which(design == i)], 1, sd, na.rm = TRUE) # Within-group SD
  corrArr[i] = cor(meansArr[,i], SDsArr[,i]) # Mean-SD Correlation
}
# ...and for the whole experiment
meansArr[,m+1] = rowMeans(unlogged.dataset, na.rm = TRUE)
SDsArr[,m+1] = apply(unlogged.dataset, 1, sd, na.rm = TRUE)
corrArr[m+1] = cor(meansArr[,m+1], SDsArr[,m+1])
meansArr[is.nan(meansArr)] = NA # To convert NaN in NA when some of the original groups are actually missing

# Scatter plot
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
printPlots("7 - SD_vs_Mean Plot") # Save just the 'Global' one
par(mfrow = c(1, 1)) # To reset sub-plotting





# Filtering ------------------------------------------------------------------------------------------------------------
# Eliminate Low-Intensity Probes

kFac = 0.75 # Minimum gene presence per group - Default=0.80
grSize = vector()
for (i in 1:m) {
  grSize[i] = sum(design == i)
}
kSize = ceiling(grSize * kFac) # 'ceiling' to be more stringent than 'round'

# Filtering Table
filTable = cbind(grSize, grSize*kFac, kSize, round(100*(kSize/grSize), 1)) # Cast to matrix
colnames(filTable) = c("Group_Size", "Min_Presence", "Rounded", "Actual %")
rownames(filTable) = groups
mgppg = paste("\nMinimum gene presence per group = ", kFac*100, "% of the samples\n\n", sep = "")
cat(mgppg)
filTable
if (saveOut) {
  write(paste("Filtering Summary File\n",
              "======================", sep = ""), "Filtering Report.txt")
  write(paste("\nPresence threshold thr0 = ", thr0, sep = ""), "Filtering Report.txt", append = TRUE)
  write(mgppg, "Filtering Report.txt", append = TRUE)
  suppressWarnings(write.table(filTable, "Filtering Report.txt", sep = "\t", append = TRUE,
                               col.names = TRUE, row.names = TRUE, quote = FALSE))
  cat("\n'Filtering Report.txt' has been saved in", myFolder, "\n\n")
}

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
retained
if (saveOut) {
  write("\n", "Filtering Report.txt", append = TRUE)
  suppressWarnings(write.table(retained, "Filtering Report.txt", append = TRUE,
                               sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE))
}

#...in at least 1 group (row-wise logical OR)
unionSet = (rowSums(presenceIndx) > 0)

filtSize = sum(unionSet)
report = paste("\nFiltering Report:\n",
               "-----------------\n",
               "  ", d[1] - filtSize, " genes have been filtered out (",
               round(((d[1] - filtSize)/d[1])*1e2, digits = 2), "%)\n",
               "  ", filtSize, " retained genes out of ", d[1], " (",
               round((filtSize/d[1])*1e2, digits = 2), "%)\n\n", sep = "")
cat(report)
if (saveOut) {
  write(paste("\n", report, sep = ""), "Filtering Report.txt", append = TRUE)
}

dataset = dataset[which(unionSet),]





# * Contrast Definition ------------------------------------------------------------------------------------------------

filtSize = dim(dataset)[1] # Reassignment to get 'filtSize' defined even in the case of no filtering

counter = 1
myContr = vector()
for (i in 1:(m-1)) { # All the possible combinations of two groups
  for (j in (i+1):m) {
    myContr[counter] = paste(groups[j], "-", groups[i], sep = "") # Character vector specifying contrasts
    counter = counter + 1
  }
}

myContr = as.matrix(myContr) # Cast to matrix
rownames(myContr) = c(1:length(myContr))
myContr
# Modify here to retain (and reorder) just the contrasts of interest
myContr = as.matrix(myContr[c(1,6,3),])
rownames(myContr) = c(1:length(myContr))
myContr





# * DE by Limma --------------------------------------------------------------------------------------------------------
# Differential Expression Assessment and Figure Production

limmaDesign = matrix(data = 0, nrow = d[2], ncol = m)
for (i in 1:m) {
  limmaDesign[which(design == i),i] = 1
}
colnames(limmaDesign) = groups
limmaDesign

fit = lmFit(dataset, limmaDesign)

contrast.matrix = makeContrasts(
  contrasts = myContr,
  levels = limmaDesign)
contrast.matrix

fit2 = contrasts.fit(fit, contrast.matrix)
efit2 = eBayes(fit2)

# Print Results (Top-Ten genes) for all the contrasts of interest
for (i in 1:length(myContr)) {
  # 'print' because automatic printing is turned off in loops (and functions)
  cat("\nDEG Top-List for contrast: ", myContr[i], "\n", sep = "")
  print(topTable(efit2, coef = i, adjust.method = "BH", sort.by = "B")) # just on-Screen
}

# Compute full DEG Tables
DEGs.limma = list() # Create an empty list
for (i in 1:length(myContr)) {
  DEGs.limma[[i]] = topTable(efit2, coef = i, number = filtSize,
                             adjust.method = "BH", sort.by = "B") # list of Data Frames
  DEGs.limma[[i]] = appendAnnotation(DEGs.limma[[i]], annot, sort.by = "adj.P.Val")
}

# Save full DEG Tables
if (saveOut) {
  for (i in 1:length(myContr)) {
    degTabName = paste("Limma - DEG Table ", myContr[i], ".txt", sep = "")
    write.table(DEGs.limma[[i]], degTabName, sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
    cat("\n'", degTabName, "' has been saved in ", myFolder, "\n\n", sep = "")
  }
}

# Summary of DEGs (you can change the Log2-Fold-Change Threshold lfc...)
results.limma = decideTests(efit2, adjust.method = "BH", p.value = 0.05, lfc = thrFC)
summary(results.limma)
vennDiagram(results.limma)
printPlots("8 - Limma Venn")

# Show Hyperparameters
d0 = efit2$df.prior           # prior degrees of freedom
dg = mean(fit$df.residual)    # original degrees of freedom
hyp = cbind(c(efit2$s2.prior, # prior variance
      mean(fit$sigma^2),      # mean sample residual variance
      mean(efit2$s2.post),    # mean posterior residual variance
      d0, dg, d0/(d0+dg)))    # Shrinkage degree
rownames(hyp) = c("Prior Var",
                  "Sample Residual Mean Var",
                  "Posterior Residual Mean Var",
                  "Prior df",
                  "Original df",
                  "Shrinkage degree")
colnames(hyp) = "Hyperparameters"
hyp

# Save significant DEG list in Excel format with annotations
all.limma.sigs = list() # Create an empty list
for (i in 1:length(myContr)) {
  all.limma.sigs[[i]] = topTable(efit2, coef = i, number = filtSize,
                                 adjust.method = "BH", sort.by = "B",
                                 p.value = 0.05, lfc = thrFC) # list of Data Frames
}
for (i in 1:length(myContr)) {
  if (saveOut & dim(all.limma.sigs[[i]])[1] > 0) {
    all.limma.sigs[[i]] = appendAnnotation(all.limma.sigs[[i]], annot, sort.by = "adj.P.Val")
    write.xlsx(all.limma.sigs[[i]], paste("Significant Genes by limma - ", myContr[i], ".xlsx", sep = ""),
               colNames = TRUE, rowNames = TRUE, sheetName = myContr[i],
               keepNA = TRUE, firstRow = TRUE) # Freezes the first row!
  }
}





# ~ DE by Limma - Paired -------------------------------------------------------------------------------------------------
# Differential Expression Assessment in Paired-Sample Design (for two-group comparison only)

# Reduce dataset to the sole two groups to be compared by paired test
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
dataset[1:2,]

# Define pairing based on sampleName
patient.ID = factor(as.numeric(substring(sampleName, regexpr("_", sampleName) + 1))) # Always to be checked !!
#patient.ID = factor(c(1,1,2,2,3,3,4,5,6,4,5,6)) # ...or go manually
cond = factor(groups[design])

limmaDesign = model.matrix(~patient.ID + cond)
limmaDesign

fit = lmFit(dataset, limmaDesign)
efit2 = eBayes(fit)

# Print Results (Top-Ten genes)
coeff = tail(colnames(limmaDesign), n = 1) # Return the last element of the vector
cat("\nDEG Top-List for contrast: ", myContr, " (design column: ", coeff, ")\n", sep = "")
topTable(efit2, coef = coeff, adjust.method = "BH", sort.by = "B") # just on-Screen

# Compute full DEG Tables
DEGs.limma = list() # Create an empty list (for compatibility with independent-sample test)
DEGs.limma[[1]] = topTable(efit2, coef = coeff, number = filtSize,
                           adjust.method = "BH", sort.by = "B") # list of Data Frames
DEGs.limma[[1]] = appendAnnotation(DEGs.limma[[1]], annot, sort.by = "adj.P.Val")

# Save full DEG Tables
if (saveOut) {
  degTabName = paste("Limma - DEG Table ", myContr, " - Paired.txt", sep = "")
  write.table(DEGs.limma[[1]], degTabName, sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
  cat("\n'", degTabName, "' has been saved in ", myFolder, "\n\n", sep = "")
}

# Summary of DEGs for paired test ...to be implemented)

# Show Hyperparameters
d0 = efit2$df.prior           # prior degrees of freedom
dg = mean(fit$df.residual)    # original degrees of freedom
hyp = cbind(c(efit2$s2.prior, # prior variance
      mean(fit$sigma^2),      # mean sample residual variance
      mean(efit2$s2.post),    # mean posterior residual variance
      d0, dg, d0/(d0+dg)))    # Shrinkage degree
rownames(hyp) = c("Prior Var",
                  "Sample Residual Mean Var",
                  "Posterior Residual Mean Var",
                  "Prior df",
                  "Original df",
                  "Shrinkage degree")
colnames(hyp) = "Hyperparameters"
hyp

# Save significant DEG list in Excel format with annotations
all.limma.sigs = list() # Create an empty list
all.limma.sigs[[1]] = topTable(efit2, coef = coeff, number = filtSize,
                               adjust.method = "BH", sort.by = "B",
                               p.value = 0.05, lfc = thrFC) # list of Data Frames

if (saveOut & dim(all.limma.sigs[[1]])[1] > 0) {
  all.limma.sigs[[1]] = appendAnnotation(all.limma.sigs[[1]], annot, sort.by = "adj.P.Val")
  write.xlsx(all.limma.sigs[[1]], paste("Significant Genes by limma - ", myContr, " - Paired.xlsx", sep = ""),
             colNames = TRUE, rowNames = TRUE, sheetName = myContr,
             keepNA = TRUE, firstRow = TRUE) # Freezes the first row!
}





# Limma Plot -----------------------------------------------------------------------------------------------------------
# MA-Plots with significant DEGs and Volcano plots

# Find Axis Limits
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
for (i in 1:length(myContr)) {
  # Mark in red/blue all the up-/down- regulated genes (+1/-1 in 'results.limma' matrix)
  plotMD(efit2, coef = i, status = results.limma[,i], values = c(1,-1), hl.col = myColors[c(2,1)],
         xlim = c(min.A.value, max.A.value), ylim = c(-max.M.value, max.M.value),
         xlab = "A (Average log-expression)",
         ylab = "M (log2-Fold-Change)")
  abline(h = 0, col = myColors[1], lty = 2) # lty = line type
  abline(h = c(thrFC,-thrFC), col = myColors[1])
  abline(v = thr0, col = myColors[1]) # Platform-specific log2-expression threshold
  printPlots(paste("9 - MA-Plot with Limma DEGs ", myContr[i], sep = ""))
}

# Volcano Plots
for (i in 1:length(myContr)) {
  
  tot.DEG = sum(DEGs.limma[[i]]$adj.P.Val < 0.05) # Total number of significant DEGs (without any FC cutoff)
  high.DEG = min(c(5,tot.DEG)) # To highlight no more than 5 genes per plot
  
  # Significance Threshold
  # Find that p-value corresponding to BH-adj.p-value ~ 0.05 (or Bonferroni point when tot.DEG = 0)
  thrP = (0.05/filtSize)*(tot.DEG + 1)
  
  # Alternative approach - Suitable also for correction methods other than BH
  # WARNING: DEG list has to be sorted by p-value or B statistics!
  #if (tot.DEG > 0) {
  #  thrP = DEGs.limma[[i]][tot.DEG + 1, "P.Value"] # This is the p-value of the first non-DEG
  #} else {
  #  thrP = (0.05/filtSize) # Bonferroni point
  #}
  
  # Check the threshold p-value
  cat("\n", myContr[i], " - Threshold Report:\n",
      "--------------------------------------\n",
      "  p-value threshold  =  ", thrP, "\n",
      "  -log10(p-value)    =  ", -log10(thrP), "\n",
      "  Gene Ranking       =  ", tot.DEG, ":", tot.DEG + 1, "\n\n", sep = "")
  
  # Check within the DEG list - WARNING: DEG list has to be sorted by p-value or B statistics!
  # It should be the alpha-crossing-point for BH-adj.p-values column...
  # ...and a thrP-containing interval for un-adjusted p-values column.
  if (tot.DEG > 0) {
    print(DEGs.limma[[i]][tot.DEG:(tot.DEG + 1), c("logFC","AveExpr","t","P.Value","adj.P.Val","B")])
    cat("\n")
  } else {
    print(DEGs.limma[[i]][1, c("logFC","AveExpr","t","P.Value","adj.P.Val","B")])
    cat("\n")
  }
  
  # Enhanced Volcano Plot
  if (getOption("append.annot")) {
    # Case-insensitive search of Gene Symbol column
    myLabels = DEGs.limma[[i]][,grep("symbol",tolower(colnames(DEGs.limma[[i]])))]
  } else {
    myLabels = rownames(DEGs.limma[[i]])
  }
  print( # NOTICE: When in a for loop, you have to explicitly print your resulting EnhancedVolcano object
    EnhancedVolcano(DEGs.limma[[i]],
                    x = "logFC",
                    y = "P.Value",
                    ylim = c(0, -log10(min.P.value)),
                    pCutoff = thrP,
                    FCcutoff = thrFC,
                    pointSize = 1,
                    col = c("black", "black", "black", myColors[2]),
                    lab = myLabels,
                    #selectLab = myLabels[1:high.DEG],
                    labSize = 4,
                    title = myContr[i],
                    subtitle = "Limma",
                    legendPosition = "none"))
  printPlots(paste("10 - Volcano with Limma DEGs ", myContr[i], sep = ""))
  
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
    printPlots(paste("10 - Volcano with Limma DEGs ", myContr[i], sep = ""))
  }
}





# * DE by RankProduct --------------------------------------------------------------------------------------------------
# Differential Expression Assessment and Figure Production

results.RP = matrix(data = 0, nrow = filtSize, ncol = length(myContr))
rownames(results.RP) = rownames(dataset)
colnames(results.RP) = myContr
for (i in 1:length(myContr)) {
  
  contr.groups = strsplit(myContr[i], split = "-", fixed = TRUE)
  
  case.id = which(groups == contr.groups[[1]][1])
  case.index = which(design == case.id)
  ctrl.id = which(groups == contr.groups[[1]][2])
  ctrl.index = which(design == ctrl.id)
  
  sub.dataset = dataset[,c(ctrl.index,case.index)]
  cl = c(rep(0,length(ctrl.index)), rep(1,length(case.index)))
  
  # Single Origin
  origin = rep(1,length(cl))
  # Multiple origins
  #origin = c(1,1,1,1,2,3,1,1,1,1,2,3,2,2) # For-looping would require a length(myContr)-long list... to be done
  
  # invisible(capture.output()) is to suppress automatic output to console
  # WARNING: therein <- (instead of =) is mandatory for assignment!
  invisible(capture.output(RP.out <- RP.advance(sub.dataset, cl, origin, gene.names = rownames(dataset),
                                                logged = TRUE, na.rm = FALSE, plot = FALSE, rand = 123)))
  invisible(capture.output(plotRP(RP.out, cutoff = 0.05)))
  
  # Compute full DEG Tables (returns a list of 2 matrices, not data frames)
  invisible(capture.output(DEGs.RP <- topGene(RP.out, logged = TRUE, logbase = 2, num.gene = filtSize)))
  for (j in 1:2) {
    DEGs.RP[[j]][,3] = log2(1/DEGs.RP[[j]][,3]) # Invert FC to get Case vs Ctrl and take the log2 values
    colnames(DEGs.RP[[j]])[3] = "Log2FC" # Correct column name
  }
  
  # Print Results (Top-Ten genes) for all the contrasts of interest
  cat("\nDEG Top-List for contrast: ", myContr[i], "\n", sep = "")
  tops = rbind(DEGs.RP$Table1[1:10,], DEGs.RP$Table2[1:10,])
  tops = tops[order(tops[,4])[1:10],] # Sort by PFP (~FDR) and take just the 'absolute' Top-Ten
  print(tops) # just on-Screen
  
  # Save full DEG Tables
  if (saveOut) {
    upDegTabName = paste("RP_Up - DEG Table ", myContr[i], ".txt", sep = "")
    dwnDegTabName = paste("RP_Down - DEG Table ", myContr[i], ".txt", sep = "")
    write.table(DEGs.RP[[1]], upDegTabName, sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
    write.table(DEGs.RP[[2]], dwnDegTabName, sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
    cat("\n'", upDegTabName, "' and '", dwnDegTabName, "' \nhave been saved in ", myFolder, "\n\n", sep = "")
  }
  
  # Fetch indexes of significant DEGs
  ups.Index = which(DEGs.RP$Table1[,4] < 0.05 & DEGs.RP$Table1[,3] > thrFC)
  dwn.Index = which(DEGs.RP$Table2[,4] < 0.05 & DEGs.RP$Table2[,3] < -thrFC)
  
  # Here 'results.RP' is filled in analogy to 'results.limma': (-1,0,+1) for (Down, NotSig, Up)-regulated genes
  # Accessing matrix rows by row name: matrix[rowname, col_i]
  results.RP[rownames(DEGs.RP$Table1)[ups.Index],i] = 1
  results.RP[rownames(DEGs.RP$Table2)[dwn.Index],i] = -1
  
  # Save significant DEG list in Excel format with annotations
  if (saveOut & (length(ups.Index) > 1 | length(dwn.Index) > 1)) {
    all.RP.sigs = rbind(DEGs.RP$Table1[ups.Index,], DEGs.RP$Table2[dwn.Index,]) # To join two objects vertically
    all.RP.sigs = appendAnnotation(all.RP.sigs, annot, sort.by = "pfp")
    write.xlsx(all.RP.sigs, paste("Significant Genes by RP - ", myContr[i], ".xlsx", sep = ""),
               colNames = TRUE, rowNames = TRUE, sheetName = myContr[i],
               keepNA = TRUE, firstRow = TRUE) # Freezes the first row!
  }
}

summary.RP = rbind(colSums(results.RP == -1), # Cast to matrix
                   colSums(results.RP == 0),
                   colSums(results.RP == 1))
rownames(summary.RP) = c("Down", "NotSig", "Up")
summary.RP





# ~ DE by RP - Paired --------------------------------------------------------------------------------------------------
# One-class Analysis of the Paired (log)-Differences

# Reduce dataset to the sole two groups to be compared by paired test
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
dataset[1:2,]

# Prepare Results matrix
results.RP = matrix(data = 0, nrow = filtSize, ncol = 1)
rownames(results.RP) = rownames(dataset)
colnames(results.RP) = myContr

# Define pairing based on sampleName
case.index = which(grepl(pg[[1]][1], sampleName))
ctrl.index = which(grepl(pg[[1]][2], sampleName))
# ...or go manually
#case.index = c(1,2,3,5,7)
#ctrl.index = c(8,9,10,4,6)

# Check the re-ordering before pairing
cbind(sampleName[case.index], sampleName[ctrl.index])

# Pair Samples
paired.dataset = dataset[,case.index] - dataset[,ctrl.index]
cl = rep(1,dim(paired.dataset)[2])

# Single Origin
origin = rep(1, length(paired.dataset))
# Multiple origins
#origin = c(1,1,1,1,2,2)

# invisible(capture.output()) is to suppress automatic output to console
# WARNING: therein <- (instead of =) is mandatory for assignment!
invisible(capture.output(RP.out <- RP.advance(paired.dataset, cl, origin, gene.names = rownames(dataset),
                                              logged = TRUE, na.rm = FALSE, plot = FALSE, rand = 123)))
invisible(capture.output(plotRP(RP.out, cutoff = 0.05)))

# Compute full DEG Tables (it returns a list of 2 matrices, not data frames)
invisible(capture.output(DEGs.RP <- topGene(RP.out, logged = TRUE, logbase = 2, num.gene = filtSize)))
for (j in 1:2) {
  DEGs.RP[[j]][,3] = log2(DEGs.RP[[j]][,3]) # Take the log2 values
  colnames(DEGs.RP[[j]])[3] = "Log2FC" # Correct column name
}

# Print Results (Top-Ten genes) for the contrast of interest
cat("\nDEG Top-List for paired contrast: ", myContr, "\n", sep = "")
tops = rbind(DEGs.RP$Table1[1:10,], DEGs.RP$Table2[1:10,])
tops = tops[order(tops[,4])[1:10],] # Sort by PFP (~FDR) and take just the 'absolute' Top-Ten
print(tops) # just on-Screen

# Save full DEG Tables
if (saveOut) {
  upDegTabName = paste("RP_Up - DEG Table ", myContr, ".txt", sep = "")
  dwnDegTabName = paste("RP_Down - DEG Table ", myContr, ".txt", sep = "")
  write.table(DEGs.RP[[1]], upDegTabName, sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
  write.table(DEGs.RP[[2]], dwnDegTabName, sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
  cat("\n'", upDegTabName, "' and '", dwnDegTabName, "' \nhave been saved in ", myFolder, "\n\n", sep = "")
}

# Fetch indexes of significant DEGs
dwn.Index = which(DEGs.RP$Table1[,4] < 0.05 & DEGs.RP$Table1[,3] < -thrFC)
ups.Index = which(DEGs.RP$Table2[,4] < 0.05 & DEGs.RP$Table2[,3] > thrFC)

# Here 'results.RP' is filled in analogy to 'results.limma': (-1,0,+1) for (Down, NotSig, Up)-regulated genes
# Accessing matrix rows by row name: matrix[rowname, col_i]
results.RP[rownames(DEGs.RP$Table1)[dwn.Index],1] = -1
results.RP[rownames(DEGs.RP$Table2)[ups.Index],1] = 1

# Save significant DEG list in Excel format with annotations
if (saveOut & (length(ups.Index) > 1 | length(dwn.Index) > 1)) {
  all.RP.sigs = rbind(DEGs.RP$Table1[dwn.Index,], DEGs.RP$Table2[ups.Index,]) # To join two objects vertically
  all.RP.sigs = appendAnnotation(all.RP.sigs, annot, sort.by = "pfp")
  write.xlsx(all.RP.sigs, paste("Significant Genes by RP - ", myContr, ".xlsx", sep = ""),
             colNames = TRUE, rowNames = TRUE, sheetName = myContr,
             keepNA = TRUE, firstRow = TRUE) # Freezes the first row!
}

summary.RP = rbind(colSums(results.RP == -1), # Cast to matrix
                   colSums(results.RP == 0),
                   colSums(results.RP == 1))
rownames(summary.RP) = c("Down", "NotSig", "Up")
summary.RP





# RP Plot --------------------------------------------------------------------------------------------------------------
# Volcano plot of the last contrast

# NOTICE_1: While matrices are constrained to a single data type, data frames can feature columns of different types.
# NOTICE_2: While duplicated row (and column) names are allowed in a matrix, they are not allowed in a data frame.
# 
# With these two points in mind, the following steps allow merging up- and down-RP DEGs in just one table,
# retaining for each duplicate probe only the one with the lower p-value.

dub.DEGs.RP = rbind(DEGs.RP[[1]], DEGs.RP[[2]]) # Join matrices vertically (duplicated row names allowed)
probe.IDs = rownames(dub.DEGs.RP) # Store row names in a saving variable
rownames(dub.DEGs.RP) = NULL # Remove row names from matrix (so it can be coerced to a data frame)
dub.DEGs.RP = as.data.frame(dub.DEGs.RP) # Cast to data frame
dub.DEGs.RP = cbind(probe.IDs, dub.DEGs.RP) # Reinsert Probe IDs as a standard entry (to be used as grouping variable)

# Remove duplicates based on P.value column condition
uni.DEGs.RP = dub.DEGs.RP[with(dub.DEGs.RP, ave(P.value, probe.IDs, FUN = min) == P.value),]
# ...that is a more readable short expression for the extended:
# uni.DEGs.RP = dub.DEGs.RP[ave(dub.DEGs.RP$P.value, dub.DEGs.RP$probe.IDs, FUN = min) == dub.DEGs.RP$P.value,]

# Second round, in case of replicate probe.IDs with the same P.value
uni.DEGs.RP = uni.DEGs.RP[with(uni.DEGs.RP, ave(pfp, probe.IDs, FUN = min) == pfp),]

# Check for duplicates
sum(duplicated(dub.DEGs.RP$probe.IDs)) # It must be dim(DEGs.RP[[i]])[1]
sum(duplicated(uni.DEGs.RP$probe.IDs)) # It must be 0

# Restore Probe_IDs as row names
rownames(uni.DEGs.RP) = uni.DEGs.RP$probe.IDs
uni.DEGs.RP = uni.DEGs.RP[,-1]

# Append annotation and sort data frame by pfp (adj.p-value)
uni.DEGs.RP = appendAnnotation(uni.DEGs.RP, annot, sort.by = "pfp")

min.P.value = min(uni.DEGs.RP$P.value) # Find y-axis limit
tot.DEG = sum(uni.DEGs.RP$pfp < 0.05) # Total number of significant DEGs (without any FC cutoff)
high.DEG = min(c(5,tot.DEG)) # To highlight no more than 5 genes per plot

# Significance Threshold
# Find that p-value corresponding to BH-adj.p-value ~ 0.05 (or Bonferroni point when tot.DEG = 0)
# WARNING: DEG list has to be sorted by p-value (or pfp)!
if (tot.DEG > 0) {
  thrP = uni.DEGs.RP[tot.DEG + 1, "P.value"] # This is the p-value of the first non-DEG
} else {
  thrP = (0.05/filtSize) # Bonferroni point
}

# Check the threshold p-value
cat("\n", myContr[length(myContr)], " - Threshold Report:\n",
    "--------------------------------------\n",
    "  p-value threshold  =  ", thrP, "\n",
    "  -log10(p-value)    =  ", -log10(thrP), "\n",
    "  Gene Ranking       =  ", tot.DEG, ":", tot.DEG + 1, "\n\n", sep = "")

# Check within the DEG list - WARNING: DEG list has to be sorted by p-value (or pfp)!
# It should be the alpha-crossing-point for adj.p-values column
if (tot.DEG > 0) {
  print(uni.DEGs.RP[tot.DEG:(tot.DEG + 1), c("gene.index", "RP/Rsum", "Log2FC", "pfp", "P.value")])
  cat("\n")
} else {
  print(uni.DEGs.RP[1, c("gene.index", "RP/Rsum", "Log2FC", "pfp", "P.value")])
  cat("\n")
}

# Enhanced Volcano Plot
if (getOption("append.annot")) {
  # Case-insensitive search of Gene Symbol column
  myLabels = uni.DEGs.RP[,grep("symbol",tolower(colnames(uni.DEGs.RP)))]
} else {
  myLabels = rownames(uni.DEGs.RP)
}
print( # NOTICE: When in a for loop, you have to explicitly print your resulting EnhancedVolcano object
  EnhancedVolcano(uni.DEGs.RP,
                  x = "Log2FC",
                  y = "P.value",
                  ylim = c(0, -log10(min.P.value)),
                  pCutoff = thrP,
                  FCcutoff = thrFC,
                  pointSize = 1,
                  col = c("black", "black", "black", myColors[2]),
                  lab = myLabels,
                  #selectLab = myLabels[1:high.DEG],
                  labSize = 4,
                  title = myContr[length(myContr)],
                  subtitle = "Rank Product",
                  legendPosition = "none"))
printPlots(paste("11 - Volcano with RP DEGs ", myContr[length(myContr)], sep = ""))





# DEG Comparison -------------------------------------------------------------------------------------------------------
# Venn diagrams of DEGs from Limma and RP methods

# To suppress 'venn.diagram()' logging
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

# Plot Venn diagrams
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
    
    venn.plot = venn.diagram(x = list(DEG.id.limma, DEG.id.RP),
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
                             cat.fontfamily = "sans")
    
    # Create a new canvas and draw the Venn
    grid.newpage()
    grid.draw(venn.plot)
    printPlots(paste("12 - Comparison Venn ", myContr[i], "_", strsplit(venn.sub, split = "-")[[1]][1], sep = ""))
  }
}





# Single Gene View -----------------------------------------------------------------------------------------------------
# Plot single-gene comparison charts (boxes/bars and jitter dots)

# Character Vector containing the Probe_IDs of the genes of interest
goi = c("A_51_P309854",
        "A_55_P2130129",
        "A_51_P506733",
        "A_55_P2007243")

# Plot single-gene data points
ct = "MS"
if (getOption("append.annot")) {
  singleGeneView(dataset, groups, design, goi, chart.type = ct, ann = annot)
} else {
  singleGeneView(dataset, groups, design, goi, chart.type = ct)
}

# Provide a quantitative readout in console
for (i in 1:length(goi)) {
  
  # Prepare
  geneStats.tit = paste("\nSingle-Gene Summary Statistics",
                        "\n==============================\n", sep = "")
  if (getOption("append.annot")) {
    geneStats.head = paste("\nProbe_ID:\t", goi[i], "\n",
                          "Gene Symbol:\t", annot[goi[i], grepl("Symbol", colnames(annot))], "\n",
                          "Gene Name:\t", annot[goi[i], grepl("Name", colnames(annot))], "\n", sep = "")
  } else {
    geneStats.head = paste("\nProbe_ID:\t", goi[i], "\n", sep = "")
  }
  geneStats.num = descStat1G(dataset[goi[i],], groups, design, 4)
  
  # Print on screen
  if (i == 1) cat(geneStats.tit, "\n", sep = "") # One-line IF
  cat(geneStats.head, "\n", sep = "")
  print(geneStats.num)
  cat("\n\n")
  
  # Print on file
  if (saveOut) {
    if (i == 1) write(geneStats.tit, "Single Gene Stats.txt") # One-line IF
    write(geneStats.head, "Single Gene Stats.txt", append = TRUE)
    
    suppressWarnings(write.table(geneStats.num, "Single Gene Stats.txt", append = TRUE, quote = FALSE,
                                 row.names = TRUE, col.names = TRUE, sep = "\t"))
    write("\n", "Single Gene Stats.txt", append = TRUE)
  }
}








####
#
# P2RX Focus --- TO BE CORRECTED and GENERALIZED
#  - usare grepl per "Symbol" come sopra
#  - inserire un IF nel caso di assenza di annotazioni
#
####

# In case of limma
for (i in 0:8) { # 0 and 8 are just "negative controls"
  cat(paste("P2RX", i, " probed: ", sep = ""))
  if (length(grep(paste("P2RX", i, sep = ""), annot[,"Gene_Symbol"])) != 0) {
    cat("YES\n")
  } else {
    cat("NO\n")
  }
}

i = 1
P2X.table = DEGs.limma[[i]][grep("P2RX", DEGs.limma[[i]][,"Gene_Symbol"]),]
P2X.table

write.xlsx(P2X.table, paste("P2RX panel - ", myContr[i], ".xlsx", sep = ""),
           colNames = TRUE, rowNames = TRUE, sheetName = myContr[i],
           keepNA = TRUE, firstRow = TRUE) # Freezes the first row!



# In case of RP
for (i in 0:8) { # 0 and 8 are just "negative controls"
  cat(paste("P2RX", i, " probed: ", sep = ""))
  if (length(grep(paste("P2RX", i, sep = ""), annot[,"GeneSymbol"])) != 0) {
    cat("YES\n")
  } else {
    cat("NO\n")
  }
}

for (i in 1:2) {
  DEGs.RP.annot = appendAnnotation(DEGs.RP[[i]], annot, sort.by = "pfp")
  P2X.table = DEGs.RP.annot[grep("P2RX", DEGs.RP.annot[,"GeneSymbol"]),]
  print(P2X.table)
  write.xlsx(P2X.table, paste("P2RX panel_", i, ".xlsx", sep = ""),
             colNames = TRUE, rowNames = TRUE, sheetName = myContr,
             keepNA = TRUE, firstRow = TRUE) # Freezes the first row!
}



