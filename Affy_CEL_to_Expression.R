# Header Info ----------------------------------------------------------------------------------------------------------
#
# Affymetrix CEL file -> Expression matrix
#
# a FeAR R-script - 13-Mar-2021
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
#   - Add annotation
#   - Save expression matrix
#





# Load Packages --------------------------------------------------------------------------------------------------------

# As a rule: try oligo first, then if it doesn't work, try affy
library(affy)  # For 3' IVT Affymetrix Arrays (older ones; e.g. Mouse 430 and Human U133 series)
library(oligo) # For GeneChip/Exon Affymetrix Arrays (newer ones + popular old platforms)

# Select Annotation Packages
library(hgu133a.db)     # Affymetrix Human Genome U133 A Set annotation data (chip hgu133a)
library(hgu133b.db)     # Affymetrix Human Genome U133 B Set annotation data (chip hgu133b)
library(hgu133plus2.db) # Affymetrix Human Genome U133 Plus 2.0 Array annotation data (chip hgu133plus2)

library(limma)          # For plotMD() function
library(openxlsx)       # Read, Write, and Edit .xlsx (Excel) Files





# Variable Definition --------------------------------------------------------------------------------------------------

system.root = "D:\\UniUPO Drive\\" # SilverLife @ Home
system.root = "D:\\Drive UniUPO\\" # SkyLake2   @ DBIOS





# Load CEL files and read raw data -------------------------------------------------------------------------------------

myFolder = "D:\\Dropbox\\WORKS in PROGRESS\\PDAC Bioinfo\\3 - Badea 2008\\Data\\CEL_files"
myFolder = "C:\\Users\\FeAR\\Desktop\\CEL_files_U133A"
setwd(myFolder)

celFiles = list.celfiles()
celFiles

# Parser to CEL files
affyRaw = read.celfiles(celFiles) # This should also install-and-load the platform design package (e.g. pd.hg.u133a)

# Plots before normalization
boxplot(affyRaw, names = c(1:length(celFiles)), xlab = "Sample Index",
        main = "Unnormalized Expression", ylab = "log2(intensity)", transfo = log2) # Raw data to be logged
MAplot(affyRaw, which = c(7,11)) # which	= integer: sample index describing which samples are to be plotted





# RMA normalization algorithm ------------------------------------------------------------------------------------------

eset = rma(affyRaw)

exprs(eset)[1:10,1:10]

# Plots after normalization
boxplot(eset, names = c(1:length(celFiles)), xlab = "Sample Index",
        main = "RMA-Normalized Expression", ylab = "log2(intensity)", transfo = identity) # Already logged (by RMA)
MAplot(eset, which = c(7,11)) # Single samples vs pseudo-median reference chip
samp = c(7,11) # One sample vs one sample
plotMD(exprs(eset)[,samp], main = paste("RMA-Normalized Expression: ", samp[1], " vs ", samp[2], sep = ""),
       xlab = "A (Average log-expression)", ylab = "M (Expression log-ratio)")





# Optional - Remove invalid probes -------------------------------------------------------------------------------------

# Remove Affymetrix control probes
probes.before = dim(eset)[1]
eset = eset[-grep("^AFFX", rownames(eset)),] # ^ anchor to match the start of string (see regular expressions)
probes.fter = dim(eset)[1]
discarded = probes.before - probes.fter
cat("\n", discarded, "Affymetrix control probes have been discarded\n\n", sep = " ")

# Check for missing values (NA) and NaN entries
any(is.na(exprs((eset))))
any(is.nan(exprs((eset))))
# Optional - Remove/handle them
# to be implemented...





# Adding Gene Annotation -----------------------------------------------------------------------------------------------

myFrame = data.frame(exprs(eset))

# List of the available annotations
#hgu133a()
#hgu133b()
hgu133plus2()

annot = data.frame(Accession   = sapply(contents(hgu133plus2ACCNUM), paste, collapse = ", "),
                   GeneSymbol  = sapply(contents(hgu133plus2SYMBOL), paste, collapse = ", "),
                   Description = sapply(contents(hgu133plus2GENENAME), paste, collapse = ", "))
naSymb = "NA"

# [OR]

# Local Annotation File
setwd(paste(system.root, "Coding\\R scripts\\Annotations\\HG-U133_Plus_2-na36-annot-csv", sep = ""))
annot = read.xlsx("HG-U133_Plus_2.na36.annot.xlsx", colNames = TRUE, rowNames = TRUE, sep.names = "_") # As Data Frame
as.matrix(colnames(annot)) # List of the available annotations
annot = annot[,c("Representative_Public_ID", "Gene_Symbol", "Gene_Title")]
naSymb = "---"
setwd(myFolder)

# Check missing annotations
nc = dim(annot)[2]
notMap = matrix(0, nrow = 1, ncol = nc, dimnames = list("Not Mapped", colnames(annot)))
for (i in 1:nc) {
  notMap[,i] = length(which(annot[,i] == naSymb))
}
notMap

# To merge two data frames horizontally by one or more common key variables:
#  - inner join (default): Return only the rows that have matching keys in both the tables (~ intersection)
#  - outer join (all = T): Return all rows from both the tables, joining the records that have matching (~ union)
#  - left outer (all.x = T): Return all rows from the left table, and any rows with matching keys from the right table
#  - right outer (all.y = T): Return all rows from the right table, and any rows with matching keys from the left table
#  - cross join (by = NULL): Return the Cartesian product
joined = merge(annot, myFrame, by.x = "row.names", by.y = "row.names", all.y = TRUE) # Right outer join

colnames(joined)[1] = "ProbeID"
joined[1:10,1:6]





# Save expression matrix -----------------------------------------------------------------------------------------------

# Without annotation
write.exprs(eset, file = "RMA-normalized_logExpression_U133Plus2.txt")
cat("\nExpression matrix (without annotations) has been saved in\n", myFolder, "\n\n", sep = "")

# With annotations
write.table(joined, file = "RMA-normalized_logExpression_U133Plus2_Annotated.txt",
            sep = "\t", col.names = TRUE, row.names = FALSE)
cat("\nExpression matrix (with annotations) has been saved in\n", myFolder, "\n\n", sep = "")




