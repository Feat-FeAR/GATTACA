# Header Info ------------------------------------------------------------------
#
# STALKER Functions
#
# A Collection of R functions to be used with:
#     GATTACA
#
# a FeAR R-script
#
# ------------------------------------------------------------------------------


#' Save a graphical output to '<folderPrefix> Figures' sub-directory.
#' 
#' Automatically makes the output folder if not there.
#' 
#' @param plotfun A function that resolves to printing a plot to an open
#'   device. This takes a function and not a plot object as some plotting
#'   facilities (notably the default one) cannot print plot objects conveniently.
#' @param figureName name of the output file (without extension)
#' @param folderPrefixprefix for naming the saving subfolder
#'   (defaults to `scriptName` option)
#' @param PNG If true, print the plot to PNG format. Defaults to `save.PNG.plot`
#'   option or TRUE (if NULL).
#' @param PDF If true, print the plot to PDF format. Defaults to `save.PDF.plot`
#'   option or TRUE (if NULL).
#' @param plot.width Width of the plot. Defaults to `plot.width` option or 820.
#' @param plot.height Height of the plot. Defaults to `plot.height` option or 600.
#' 
#' @author FeAR
printPlots = function(
  plotfun,
  figureName,
  folderPrefix = getOption("scriptName"),
  PNG = getOption("save.PNG.plot"), PDF = getOption("save.PDF.plot"),
  plot.width = getOption("plot.width"), plot.height = getOption("plot.height")
  ) {
  library(logger)
  flag = FALSE # A dummy flag to insert a couple of 'new lines' in case of WARNINGs

  # Check argument values
  # NOTE: considering that getOption("...") returns NULL for undefined arguments,
  #       IFs are evaluated only when:
  #           the corresponding global option is not defined
  #             AND
  #           no argument is passed runtime
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
  
  if (is.null(plot.width)) { 
    plot.width = 820
  }
  if (is.null(plot.height)) { 
    plot.height = 600
  }
  
  fullName = file.path(figSubFolder, figureName, fsep = .Platform$file.sep) # OS-independent path separator
  
  if (!file.exists(figSubFolder) && (PNG || PDF)) {
    dir.create(figSubFolder)
    flag = TRUE
    log_info("New folder '", figSubFolder, "' has been created in the current WD", sep = "")
  }
  if (PNG) { # invisible(capture.output()) to suppress automatic output to console
    log_info(paste0("Saving ", figureName, ".png ..."))
    png(
      file = paste(fullName, ".png", sep = ""),
      width = plot.width, height = plot.height
    )
    plotfun()
    dev.off()
  }
  if (PDF) {
    log_info(paste0("Saving ", figureName, ".pdf ..."))
    pdf(
      file = paste(fullName, ".pdf", sep = ""),
      width = plot.width, height = plot.height
    )
    plotfun()
    dev.off()
  }
  if (flag) {
    cat("\n\n")
  }
}


#' Append annotation to DEG-statistic top-table and sort
#' 
#' @param gene.stat The table of genes, usually a DEG summary-statistic top-table
#'   or an expression matrix.
#' @param ann the matrix containing the annotation data
#' @param do.the.job If false, don't do anything. Defaults to `append.annot`
#'   option or TRUE if the option in NULL.
#' @param sort.by the name or index of the column used to sort the final data set
#' 
#' @author FeAR
appendAnnotation = function(gene.stat, ann, do.the.job = getOption("append.annot"), sort.by = 1)
{
  #' 
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


#' Return basics descriptive statistics of a single gene, by group label
#' 
#' @param gene Numeric vector or single-row data frame from gene expression matrix
#' @param gr Group names
#' @param des Experimental design (full design mode vector)
#' @param prec Decimal precision. Defaults to 4.
descStat1G = function(gene, gr, des, prec = 4)
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
    stat.frame[i,3] = round(mean(n.gene), digits = prec)
    stat.frame[i,4] = round(var(n.gene), digits = prec)
    stat.frame[i,5] = round(sd(n.gene), digits = prec)
    stat.frame[i,6] = round(sd(n.gene)/sqrt(sum(des == i)), digits = prec) # SEM
  }
  
  return(stat.frame)
}


#' Plot single gene comparison chart
#' 
#' @param exp.mat Expression matrix (as data frame)
#' @param gr Group names
#' @param des Experimental design (full design mode vector)
#' @param gois Genes of interest by probe (char vector)
#' @param chart.type "BP" (Box Plot), "BC" (Bar Chart), or "MS" (Mean & SEM)
#'   Defaults to "BP".
#' @param ann Optional annotation data frame.
#' 
#' @author FeAR
singleGeneView = function(exp.mat, gr, des, gois, chart.type = "BP", ann = NULL)
{
  geo = switch(
    chart.type,
    "BP" = "point",
    "BC" = "bar",
    "MS" = "crossbar"
  )
  
  for (i in 1:length(gois)) {
    
    var.expr = as.numeric(exp.mat[gois[i],]) # Downcast to vector
    var.groups = gr[des]
    sgex = data.frame(var.expr, var.groups) # Single Gene Expression Data Frame
    sgs = descStat1G(exp.mat[gois[i],], gr, des, 6) # Single Gene Summary Data Frame
    
    if (is.null(ann)) {
      gene.symb = ""
    } else {
      gene.symb = paste(ann[gois[i], grepl("Symbol", colnames(ann))], " - ", sep = "")
    }
    
    if (chart.type == "BP") {
      
      print( # NOTICE: When in a for loop, you have to explicitly print your resulting ggplot object
        ggplot(data = sgex, aes(var.groups, var.expr)) +
          theme_bw(base_size = 15, base_rect_size = 1.5) +
          xlab("Group") + # In the following functions, when data=NULL (default), the data is inherited from ggplot()
          ylab("log2 Expression") +
          ggtitle(label = "Box Plot with Jitter", subtitle = paste(gene.symb, "Probe ID: ", gois[i], sep = "")) +
          geom_boxplot(width = 0.5, size = 0.5, notch = TRUE, outlier.shape = NA) +
          stat_summary(fun = "mean", geom = geo, color = "red3", size = 2) +
          geom_jitter(position = position_jitter(width = 0.1, height = 0, seed = 123), size = 1.5))
      
    } else if (chart.type == "BC" | chart.type == "MS") {
      
      print(
        ggplot(data = sgex, aes(var.groups, var.expr)) +
          theme_bw(base_size = 15, base_rect_size = 1.5) +
          xlab("Group") +
          ylab("log2 Expression") +
          ggtitle(label = "Mean & SEM Plot with Jitter", subtitle = paste(gene.symb, "Probe ID: ", gois[i], sep = "")) +
          stat_summary(fun = "mean", geom = geo, color = "black", size = 0.5, width = 0.2) +
          # Recommended alternative for bar charts in ggplot2:
          #geom_bar(data = sgs, aes(GROUP, MEAN), stat = "identity", color = "black", size = 0.5, width = 0.2) +
          geom_errorbar(data = sgs, aes(GROUP, MEAN, ymin = MEAN - SEM, ymax = MEAN + SEM), size = 1, width = 0.1) + 
          geom_jitter(position = position_jitter(width = 0.1, height = 0, seed = 123), size = 1.5))
      
    } else {
      
      cat("\n")
      stop("Invalid chart.type!\n\n")
    }
    
    printPlots(paste("SingleGene Plot - ", chart.type, " - ", gois[i], sep = ""))
  }
}


# https://stackoverflow.com/questions/14469522/stop-an-r-program-without-error
#' Stop a program anywhere, but with no errors.
stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}


#' Make a pistop function.
#' 
#' A pitstop function prints out a prompt that the user can reply to
#' to either continue or silently stop the program.
#' 
#' @param check A boolean that governs if all generated pitstops run (`TRUE`)
#'   or not (`FALSE`). This is useful to globally skip all pitstops.
#' 
#' @returns A `pitstop` function that takes a message as prompt.
#' 
#' @author Hedmad
pitstop.maker <- function(check) {
  pitstop <- function(message) {
    if (check) {
      readline(prompt=paste0(message, " Continue? [yes/NO]: ")) |>
        tolower() ->
        response
      if (response %in% c("yes", "ye", "y")) {
        return()
      } else {
        stop_quietly()
      }
    }
  }
  return(pitstop)
}


#' Print out the topleft part of the input data. Useful to inspect the data
#' 
#' @param data The data to be printed. Has to be subsetted with [] and processed
#'   by `head()`. This works for dataframes and matrices.
#'   
#' @author Hedmad
topleft.head <- function(data) {
  y <- dim(data)[2]
  floored.y <- floor(y * 0.25)
  if (floored.y == 1 & y >= 2) {floored.y <- 2}
  
  return(head(data[, c(1:min(5, floored.y))]))
}


#' Make a printif. A printif prints only if a global check passes.
#' 
#' @param check A boolean that governs if all generated printifs run (`TRUE`)
#'   or not (`FALSE`). This is useful to globally skip all prints.
#' @param applied.fun A function applied to the input of all printifs.
#'   Useful if all inputs need to be preprocessed in the same way.
#'   
#' @author Hedmad
printif.maker <- function(check, applied.fun = function(x) {x}) {
  printif <- function(incoming) {
    if (check) {
      # TODO :: this cannot be suppressed by `suppressMessages`
      print(applied.fun(incoming))
    }
  }
}


#' Capture the output of `print(data)` to a string.
#' 
#' Only supports data that can be printed.
#' 
#' @param data Object to be printed.
#' 
#' @returns A single string.
#' 
#' @author Hedmad.
get.print.str <- function(data) {
  return(paste0(capture.output(print(data)), collapse = "\n"))
}
