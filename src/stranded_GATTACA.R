### STRANDED GATTACA PIECES - Will not run on their own.

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
