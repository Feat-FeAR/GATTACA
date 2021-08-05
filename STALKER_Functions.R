# Header Info ----------------------------------------------------------------------------------------------------------
#
# STALKER Functions
#
# A Collection of R functions to be used with:
#     GATTACA
#
# a FeAR R-script - 04-Aug-2021
#





# Save a graphical output to '<folderPrefix> Figures' sub-directory
#
#   figureName    = name of the output file (without extension)
#   folderPrefix  = prefix for naming the saving subfolder (default to the name of the parent script sourcing this one)
#   PNG           = boolean: T to print the currently displayed figure in PNG format
#   PDF           = boolean: T to print the currently displayed figure in PDF format
#
printPlots = function(figureName, folderPrefix = getOption("scriptName"),
                      PNG = getOption("save.PNG.plot"), PDF = getOption("save.PDF.plot"))
{
  flag = FALSE # A dummy flag to insert a couple of 'new lines' in case of WARNINGs
  
  # Check argument values
  # NOTE: considering that getOption("...") returns NULL for undefined arguments, IFs are evaluated only when:
  # the corresponding global option is not defined
  #  AND
  # no argument is passed runtime
  if (is.null(PNG)) { 
    PNG = TRUE
    flag = TRUE
    cat("\nWARNING: \'save.PNG.plot\' option defaulted to TRUE")
  }
  if (is.null(PDF)) {
    PDF = TRUE
    flag = TRUE
    cat("\nWARNING: \'save.PDF.plot\' option defaulted to TRUE")
  }
  if (is.null(folderPrefix)) {
    figSubFolder = "Figures"
  } else {
    figSubFolder = paste(folderPrefix, "Figures", sep = " ")
  }
  
  fullName = file.path(figSubFolder, figureName, fsep = .Platform$file.sep) # OS-independent path separator
  
  if (!file.exists(figSubFolder) && (PNG || PDF)) {
    dir.create(figSubFolder)
    flag = TRUE
    cat("\nNew folder '", figSubFolder, "' has been created in the current WD", sep = "")
  }
  if (PNG) { # invisible(capture.output()) to suppress automatic output to console
    invisible(capture.output(
      dev.print(device = png, filename = paste(fullName, ".png", sep = ""), width = 820, height = 600)))
  }
  if (PDF) {
    invisible(capture.output(
      dev.print(device = pdf, paste(fullName, ".pdf", sep = ""))))
  }
  if (flag) {
    cat("\n\n")
  }
}





# Append annotation to DEG-statistic top-table and sort (do nothing if do.the.job == FALSE)
#
#   gene.stat   = the table of genes, usually a DEG summary-statistic top-table (or an expression matrix)
#   ann         = the matrix containing the annotation data
#   do.the.job  = F to skip the appending task by global settings, without the need for an external IF
#   sort.by     = the name or index of the column used to sort the final data set
#
appendAnnotation = function(gene.stat, ann, do.the.job = getOption("append.annot"), sort.by = 1)
{
  # Check argument values
  if (is.null(do.the.job)) {
    do.the.job = TRUE
    cat("\nWARNING: \'append.annot\' option defaulted to TRUE\n\n")
  }
  
  if (do.the.job) {
    # 'merge' function to merge two matrix-like objects horizontally and cast to data frame (right outer join)
    # NOTICE: both gene.stat and ann are supposed to have the Probe_IDs as rownames
    joined = merge(ann, gene.stat, by.x = "row.names", by.y = "row.names", all.y = TRUE)
    rownames(joined) = joined[,1]
    gene.stat = joined[,-1]
    
    # Re-sort the data frame by the content of 'sort.by' column ('sort.by' can be either a number or a column name)
    gene.stat = gene.stat[order(gene.stat[,sort.by]),]
  }
  
  return(gene.stat)
}





# Return basics descriptive statistics of a single gene, by group label
#
#   gene  = Numeric vector or single-row data frame from gene expression matrix
#   gr    = Group names
#   des   = Experimental design (full design mode vector)
#
descStat1G = function(gene, gr, des)
{
  # Define a new empty data frame
  stat.frame = data.frame(GROUP = character(),
                          n = integer(),
                          MEAN = double(),
                          VAR = double(),
                          SD = double(),
                          SEM = double(),
                          stringsAsFactors = FALSE)
  
  for (i in 1:length(gr)) {
    
    n.gene = as.numeric(gene[des == i]) # Downcast to numeric vector
    
    stat.frame[i,1] = gr[i]
    stat.frame[i,2] = sum(des == i)
    stat.frame[i,3] = mean(n.gene)
    stat.frame[i,4] = var(n.gene)
    stat.frame[i,5] = sd(n.gene)
    stat.frame[i,6] = sd(n.gene)/sqrt(sum(des == i)) # SEM
  }
  
  return(stat.frame)
}





# Plot single gene comparison chart
#
#   exp.mat     = Expression matrix (as data frame)
#   gr          = Group names
#   des         = Experimental design (full design mode vector)
#   gois        = Genes of interest by probe (char vector)
#   chart.type  = "BP" (Box Plot), "BC" (Bar Chart), or "MS" (Mean & SEM)
#
singleGeneView = function(exp.mat, gr, des, gois, chart.type = "BP")
{
  geo = switch(chart.type,
               "BP" = "point",
               "BC" = "bar",
               "MS" = "crossbar")
  
  for (i in 1:length(gois)) {
    
    var.expr = as.numeric(exp.mat[gois[i],]) # Downcast to vector
    var.groups = gr[des]
    sgex = data.frame(var.expr, var.groups) # Single Gene Expression Data Frame
    sgs = descStat1G(exp.mat[gois[i],], gr, des) # Single Gene Summary Data Frame
    
    if (chart.type == "BP") {
      
      print( # NOTICE: When in a for loop, you have to explicitly print your resulting ggplot object
        ggplot(data = sgex, aes(var.groups, var.expr)) +
          theme_bw(base_size = 15, base_rect_size = 1.5) +
          xlab("Group") + # In the following functions, when data=NULL (default), the data is inherited from ggplot()
          ylab("log2 Expression") +
          ggtitle(label = "Box Plot with Jitter", subtitle = paste("Probe ID: ", gois[i], sep = "")) +
          geom_boxplot(width = 0.5, size = 0.5, notch = TRUE, outlier.shape = NA) +
          stat_summary(fun = "mean", geom = geo, color = "red3", size = 2) +
          geom_jitter(position = position_jitter(width = 0.1, height = 0, seed = 123), size = 1.5))
      
    } else if (chart.type == "BC" | chart.type == "MS") {
      
      print(
        ggplot(data = sgex, aes(var.groups, var.expr)) +
          theme_bw(base_size = 15, base_rect_size = 1.5) +
          xlab("Group") +
          ylab("log2 Expression") +
          ggtitle(label = "Box Plot with Jitter", subtitle = paste("Probe ID: ", gois[i], sep = "")) +
          stat_summary(fun = "mean", geom = geo, color = "black", size = 0.5, width = 0.2) +
          # Recommended alternative for bar charts in ggplot2:
          #geom_bar(data = sgs, aes(GROUP, MEAN), stat = "identity", color = "black", size = 0.5, width = 0.2) +
          geom_errorbar(data = sgs, aes(GROUP, MEAN, ymin = MEAN - SEM, ymax = MEAN + SEM), size = 1, width = 0.1) + 
          geom_jitter(position = position_jitter(width = 0.1, height = 0, seed = 123), size = 1.5))
      
    } else {
      
      cat("\n")
      stop("Invalid chart.type!\n\n")
      
    }
  }
}