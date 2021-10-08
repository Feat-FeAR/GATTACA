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
  folderPrefix = getOption("scriptName", ""),
  PNG = getOption("save.PNG.plot", TRUE), PDF = getOption("save.PDF.plot", TRUE),
  plot.width = getOption("plot.width", 16), plot.height = getOption("plot.height", 9),
  png_ppi = getOption("png_ppi")
  ) {
  library(logger)

  # Check argument values
  # NOTE: considering that getOption("...") returns NULL for undefined arguments,
  #       IFs are evaluated only when:
  #           the corresponding global option is not defined
  #             AND
  #           no argument is passed runtime
  if (is.null(PNG)) { 
    PNG = TRUE
    log_warn("save.PNG.plot option defaulted to TRUE")
  }
  if (is.null(PDF)) {
    PDF = TRUE
    log_warn("save.PDF.plot option defaulted to TRUE")
  }

  figSubFolder = paste(folderPrefix, "Figures", sep = " ")
  
  fullName = file.path(figSubFolder, figureName)
  
  if (!file.exists(figSubFolder) && (PNG || PDF)) {
    dir.create(figSubFolder)
    log_info("New folder '", figSubFolder, "' has been created in the current WD", sep = "")
  }
  if (PNG) { # invisible(capture.output()) to suppress automatic output to console
    log_info(paste0("Saving ", figureName, ".png ..."))
    png(
      file = paste(fullName, ".png", sep = ""),
      width = plot.width*png_ppi, height = plot.height*png_ppi
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
appendAnnotation = function(
  gene.stat, ann,
  do.the.job = getOption("append.annot"),
  sort.by = 1
  )
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
  stat.frame = data.frame(
    GROUP = character(),
    n = integer(),
    MEAN = double(),
    VAR = double(),
    SD = double(),
    SEM = double(),
    stringsAsFactors = FALSE
  )
  
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
#' @author MrHedmad
pitstop.maker <- function(check) {
  pitstop <- function(message) {
    if (check) {
      cat(paste0(message, " Continue? [yes/NO]: "))
      readLines("stdin", n=1) |> tolower() -> response
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
#' @author MrHedmad
topleft.head <- function(data) {
  y <- ncol(data)

  return(head(data[, c(1:min(5, max(2, y)))]))
}


#' Make a printif. A printif prints only if a global check passes.
#' 
#' @param check A boolean that governs if all generated printifs run (`TRUE`)
#'   or not (`FALSE`). This is useful to globally skip all prints.
#' @param applied.fun A function applied to the input of all printifs.
#'   Useful if all inputs need to be preprocessed in the same way.
#'   
#' @author MrHedmad
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
#' @author MrHedmad
get.print.str <- function(data) {
  return(paste0(capture.output(print(data)), collapse = "\n"))
}


#' Finds all the chars in `possibilities` that are also in `input`,
#' returning a vector of unique chars.
#' 
#' @param str1 The first string
#' @param str2 The second string
#' 
#' @author MrHedmad
str_intersection <- function(str1, str2) {
  intersection <- c()
  str1 <- strsplit(str1, "")[[1]]
  str2 <- strsplit(str2, "")[[1]]

  for (char in str1) {
    if (char %in% str2 && ! char %in% intersection) {
      intersection <- c(intersection, char)
    }
  }

  return(intersection)
}

#' Find all numbers in a string and return them in a vector.
#' 
#' @author MrHedmad
find_nums <- function(chars) {
  library(stringr)
  pattern <- "([0-9])+"
  if (length(grep(pattern, chars)) == 0) {
    # There are no numbers here...
    return(NULL) 
  }
  
  captures <- str_extract_all(chars, pattern, simplify = TRUE)
  numbers <- sapply(captures, as.numeric)
  
  return(numbers)
}


#' Removes all characters in `to.remove` from `original` returning a string.
#'
#' @author MrHedmad
subtract.str <- function(to.remove, original) {
  result <- gsub(paste0("[", to.remove, "]"), "", original)
  return(result)
}

#' Expand labels with the intercalated strategy.
#' 
#' @param core A string of comma separated labels to expand
#' @param times The number of times to expand
#' @param initial.num The initial number to use to expand wildcards.
#' 
#' @returns The expanded labels as a single string
#' 
#' @author MrHedmad
expand_intercalated <- function(core, times, initial.num) {
  result <- ""
  rep.number <- initial.num
  for (i in rep(0, times)) {
    uniquestr <- gsub("\\*", rep.number, core)
    result <- paste(result, uniquestr, sep = ",")
    rep.number <- rep.number + 1
  }
  # The result starts with an extra , as there is a paste("", ..., sep = ",")
  # This gets rid of it.
  result <- paste0(strsplit(result, "")[[1]][-1], collapse = "")
  
  return(result)
}


#' Expand labels with the ordered strategy.
#' 
#' @param core A string of comma separated labels to expand
#' @param times The number of times to expand
#' @param initial.num The initial number to use to expand wildcards.
#' 
#' @returns The expanded labels as a single string
#' 
#' @author MrHedmad
expand_ordered <- function(core, times, initial.num) {
  result <- ""
  rep.number <- initial.num
  for (entry in strsplit(core, ",")[[1]]) {

    for (i in rep(0, times)){
      uniquestr <- gsub("\\*", rep.number, entry)
      result <- paste(result, uniquestr, sep = ",")
      rep.number <- rep.number + 1
    }

    rep.number <- initial.num
  }
  
  # The result starts with an extra , as there is a paste("", ..., sep = ",")
  # This gets rid of it.
  result <- paste0(strsplit(result, "")[[1]][-1], collapse = "")
  
  return(result)
}


#' Get the capture groups from a regexec and regmatches call.
#' 
#' @param patters A RegEx pattern to use.
#' @param string The string to test against
#' 
#' @returns A vector with the first match found in first position and the
#' capture groups in the rest of the vector.
#' 
#' @author MrHedmad
get_captures <- function(pattern, string) {
  matches <- regexec(pattern, string, perl = TRUE)
  captures <- regmatches(string, matches)
  
  return(captures[[1]])
}


#' Parser for experimental designs
#' Each sample is marked with a letter (or series of letters) representing
#' conditions. They can also be marked with arbitrary numbers representing
#' sample pairs for paired designs. The letter and optional numbers are known as
#' labels, and are separated by commas (all spaces are ignored).
#'
#' The parser respects already expanded labels
#'   > design_parser("a, b, c") # Unpaired
#'     [1] "a, b, c"
#'   > design_parser("a1, b2, c3") # Paired
#'     [1] "a1, b2, c3"
#' To avoid high repetition, the parser supports two types of pattern expansion:
#' Round brackets represent intercalated expansions, where labels inside the
#' brackets are repeated in the order inside the brackets for the number
#' of times specified:
#'   > design_parser("(a, b):3") # Unpaired 'intercalated'
#'     [1] "a, b, a, b, a, b"
#' Square brackets represent ordered expansions, where each label in the brackets
#' is repeated on its own the number of times specified:
#'   > design_parser("[a, b]:3") # Unpaired 'ordered'
#'     [1] "a, a, a, b, b, b"
#'   > design_parser("[a, b, c]:3") # Unpaired 'ordered'
#'     [1] "a, a, a, b, b, b, c, c, c"
#' The * wildcard inside the brackets will be replaced with unique numbers, so
#' that expansion is actually useful for paired designs:
#'   > design_parser("(a*, b*):3") # Paired 'intercalated'
#'     [1] "a1, b1, a2, b2, a3, b3"
#'   > design_parser("[a*, b*]:3") # Paired 'ordered'
#'     [1] "a1, a2, a3, b1, b2, b3"
#'   > design_parser("[a*, b*, c*]:2") # Paired 'ordered'
#'     [1] "a1, a2, b1, b2, c1, c2"
#' The wildcard is assured to not collide with any already used numbers in the
#' pattern, starting one number after the largest number found:
#'   > design_parser("a3, b7, (a*, b*):3") # Paired 'intercalated'
#'     [1] "a3, b7, a8, b8, a9, b9, a10, b10"
#'   > design_parser("a3, b7, [a*, b*]:3") # Paired 'ordered'
#'     [1] "a3, b7, a8, a9, a10, b8, b9, b10"
#'  The various patterns can be mixed and matched:
#'   > design_parser("(a*, b*):3, a3, b6, [a*]:2")
#'     [1] "a7, b7, a8, b8, a9, b9, a3, b6, a10, a11"
#'  Note: Both types of expansion work identically if only one label is specified
#'  in the brackets.
#'
#'  @param rawstr The raw design string to be parsed.
#'  
#'  @author MrHedmad
design_parser <- function(rawstr) {
  rawstr <- gsub(" ", "", rawstr)
  
  # The `:x` modifiers pollute finding the max patient number. This removes
  # them.
  str_no_times = gsub(":[0-9]+", "", rawstr)
  used.numbers <- find_nums(str_no_times)
  initial.num <- if (is.null(used.numbers)) {1} else {max(used.numbers) + 1}
  
  # Split the initial string by commas *not* in parenthesis.
  splitted <- strsplit(rawstr, ",(?![^(^\\[]*[\\)\\]])", perl = TRUE)[[1]]
  

  # I put a lot of comments here as there are many nested ifs and its hard
  # to know what is doing what.
  intercalated.pattern <- "\\((.*?)\\):([0-9]+)"
  ordered.pattern <- "\\[(.*?)\\]:([0-9]+)"
  
  result <- ""
  
  for (unexpanded in splitted) {
    # There is some code duplication here that I leave for clarity as
    # this bit is complex enough already...
    if (length(grep(intercalated.pattern, unexpanded)) == 1) {
      captures <- get_captures(intercalated.pattern, unexpanded)
      expanded <- expand_intercalated(
        core = captures[2],
        times = as.integer(captures[3]),
        initial.num = initial.num
        )
      # As we've used as.integer(captures[3]) many nums, we need to increase
      # this number that much for the other expansions.
      initial.num <- initial.num + as.integer(captures[3])
      
      result <- paste0(result, expanded, sep = ",")
      next
    }
    if (length(grep(ordered.pattern, unexpanded)) == 1) {
      captures <- get_captures(ordered.pattern, unexpanded)
      expanded <- expand_ordered(
        core = captures[2],
        times = as.integer(captures[3]),
        initial.num = initial.num
      )
      # As we've used as.integer(captures[3]) many nums, we need to increase
      # this number that much for the other expansions.
      initial.num <- initial.num + as.integer(captures[3])
      
      result <- paste0(result, expanded, sep = ",")
      next
    }
    
    # If we get here then the unexpanded pattern needs no expanding
    result <- paste0(result, unexpanded, sep = ",")
  }

  # The result starts with an extra , as there is a paste("", ..., sep = ",")
  # This gets rid of it.
  result <- gsub('.{1}$', '', result)
  result <- gsub(",", ", ", result) # Some extra space to be easy on the eyes
  
  # The new design just needs to be splitted again...
  result <- strsplit(result, ", ")[[1]]
  return(result)
}


#' Split a parsed design into the groups and the pairings vectors.
#' 
#' Essentially, the character (a-z) portion of each entry will be included in the 
#' groups vector, while the numerical portion the same for the pairings vector.
#' 
#' @param experimental_design An experimental design vector (such as one parsed
#'   by `design_parser`).
#' @returns A list with the groups vector of str in slot $groups and the pairings
#'   vector of ints in the $pairings vector.
#'   
#' @author MrHedmad
split_design <- function(experimental_design) {
  library(purrr)

  get_group <- function(x) {
    # Get rid of any numbers
    x <- gsub("[0-9]+", "", x, perl = TRUE)
    return(x)
  }

  get_pairings <- function(x) {
    # Get rid of any letters
    x <- gsub("[a-zA-Z]+", "", x, perl = TRUE)
    if (x == "") {
      return(NA)
    }
    return(as.integer(x))
  }
  # This complains about some comparison, but it seems to be a false positive.
  result <- suppressWarnings(
      list(
      groups = unlist(map(experimental_design, get_group)),
      pairings = unlist(map(experimental_design, get_pairings))
    )
  )
  return(result)
}


#' Makes the elements of a character vector unique by appending sequence numbers
#' to duplicates. Just like make.unique {base}, but starting from 1.
#' 
#' @param exp_design The input character vector
#' @param sep The name-number separator to be used
#' 
#' @returns A list with the groups vector of str in slot $groups and the pairings
#'   vector of ints in the $pairings vector.
#'   
#' @author FeAR
make.unique_from1 <- function(exp_design, sep = ".") {
  result <- vector(mode = "character", length = length(exp_design))
  unique_labels <- unique(exp_design)
  
  for (grp in unique_labels) {
    index = which(exp_design == grp)
    new.names = paste0(grp, sep, c(1:length(index)))
    result[index] = new.names
  }
  return(result)
}
  
  
  
