# Header Info ------------------------------------------------------------------
#
# This file hosts some tool functions that are used throughout the package.
#
# It also sources other files like it, for easy access.
#
# ------------------------------------------------------------------------------

log_debug("Sourcing the 'tools.R' file.")

# Source the other tool files.
# Note: The other tool files are dependent on each other, so you should source
# all of them at once (such as by sourcing this file).
source(file.path(ROOT, "src", "tools", "plots.R"))
source(file.path(ROOT, "src", "tools", "design_parser.R"))
source(file.path(ROOT, "src", "tools", "annotations.R"))


# https://stackoverflow.com/questions/14469522/stop-an-r-program-without-error
#' Stop a program anywhere, but with no errors.
stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}


#' Make a pitstop function.
#'
#' A pitstop function prints out a prompt that the user can reply to
#' to either continue or silently stop the program.
#'
#' @param check A bool that governs if all generated pitstops run (`TRUE`)
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


#' Print out the topleft part of the input data. Useful to inspect data that
#' is both tall and wide
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


#' Create a function to push to a data log saved in getOption("gattaca.datalog.path").
#'
#' The function creates the path to the data log when making a new
#' function, and also overwrites the contents of identically-named files.
#'
#' A data log is a weird log where snippets of data are saved for future
#' inspection. The log therefore has to log multi-line strings without any
#' log levels. This function makes functions that can push multi-line strings
#' to such logs.
#'
#' If the data log path is not found, this function logs a warning and returns
#' a function that does nothing.
#'
#' @param prefix The prefix to give to every log entry, either as a string or
#'   a function that returns a string. "<m>" Is replaced with the contents of
#'   the optional "message" parameter passed to the resulting log pusher.
#' @param suffix The suffix to give to every log entry, either as a string or
#'   a function that returns a string.
#' @param padding If set, prefixes every line pushed to the log (detected by
#'   newline characters) with the padding. Useful if some indentation is
#'   wanted for the pushed data.
#'
#' @returns A function that takes a mandatory `string_data` parameter with the
#'   data to log, and an optional `message` to describe the data.
make_data_log_pusher <- function(
    prefix = \(){paste("[", date(), "] <m> >>>>", sep = "")},
    suffix = "<<<<",
    padding = "    "
  ) {

  configured_push <- function(string_data, message = ""){
    path <- getOption("gattaca.datalog.path")

    if (is.null(path)) {
      log_warn("No datalog path was found. Skipping istantiating the data log.")
      return(\(...){})
    }
    dir <- dirname(path)
    dir.create(dir, showWarnings = FALSE)

    sprefix <- if (is.function(prefix)) {prefix()} else {prefix}
    sprefix <- gsub("<m>", message, sprefix)
    ssuffix <- if (is.function(suffix)) {suffix()} else {suffix}

    string_data <- paste0(padding, gsub("\\n", paste0("\n", padding), string_data))

    pushed_string <- paste(sprefix, string_data, ssuffix, sep = "\n")

    # This function behaves like `tee`
    cat(pushed_string)
    write(pushed_string, file = path, append = TRUE)
  }
}

push_to_data_log <- make_data_log_pusher()

log_data <- function(data, message = "", shorten = TRUE, is.string = FALSE) {
      if (is.string) {
        sdata <- data
      } else {
        sdata <- if (shorten) {get.print.str(topleft.head(data))} else {get.print.str((data))}
      }
      push_to_data_log(sdata, message = message)
    }

#' Finds all the chars in `str1` that are also in `str2`,
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

  names(numbers) <- NULL

  return(numbers)
}


#' Removes all characters in `to.remove` from `original` returning a string.
#'
#' @author MrHedmad
subtract.str <- function(to.remove, original) {
  if (to.remove == "") {
    return(original)
  }
  result <- gsub(paste0("[", to.remove, "]"), "", original)
  return(result)
}


#' This function saves expression data in the correct format, handling
#' moving columns around and things.
#'
#' Expression data is any dataframe with rownames saved to the csv as `probe_id`,
#' so can technically hold anything. Is mainly used in this package to save
#' either DEG tables or expression sets.
#'
#' @param expression_data The expression data to save. A data.frame with
#'   probe ids as rownames.
#' @param target Path to where the file will be saved.
#' @param verbose Should the function write to the log file?
#'
#' @author MrHedmad
write_expression_data <- function (
  expression_data, target = "expression_data.csv",
  verbose = TRUE
) {
  expression_data$probe_id <- rownames(expression_data)
  rownames(expression_data) <- NULL

  expression_data %>% dplyr::select("probe_id", everything()) -> expression_data

  if (verbose) {
    log_info(paste(
      "Saving a", ncol(expression_data) - 1, "cols by",
      nrow(expression_data), "rows expression dataset"
    ))
  }

  write.csv(expression_data, target, row.names = FALSE, quote = TRUE)
}


#' Read and extract row names from a .csv containing expression data, such as
#' one saved by the `write_expression_data` function.
#'
#' @param target Path to the file to read.
#'
#' @returns A data.frame with  the loaded data, with probe ids as row names.
#'
#' @author MrHedmad
read_expression_data <- function (target, verbose = TRUE) {
  expression_data <- read.csv(target)
  rownames(expression_data) <- expression_data$probe_id
  expression_data$probe_id <- NULL

  if (verbose) {
    log_info(paste(
      "Loaded a", ncol(expression_data), "cols by",
      nrow(expression_data), "rows expression dataset from '",
      target, "'"
    ))
  }

  return(expression_data)
}


#' Run quantile-quantile normalization on an expression data set.
#'
#' The qq normalization is really powerful, but may distort results.
#'
#' @param expression_set The set to normalize.
#'
#' @returns The normalized expression dataset, with the same number and names
#'   of columns and rows.
#'
#' @author MrHedmad
qq_normalize <- function(expression_set) {
  library(preprocessCore)

  log_info("Running quantile-quantile normalization...")
  expression_set |> as.matrix() |> normalize.quantiles() |> as.data.frame() ->
    normalized_data

  rownames(normalized_data) <- rownames(expression_set)
  colnames(normalized_data) <- colnames(expression_set)

  return(normalized_data)
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


#' Get the row wise median of a data frame.
#'
#' The rownames are preserved. The colname is replaced to be "Median".
#'
#' @param data_set The dataframe to process.
#'
#' @author Hedmad
get_median <- function(data_set) {
  new_set <- apply(data_set, MARGIN = 1, FUN = median)
  new_set |> as.data.frame() -> new_set
  rownames(new_set) <- rownames(data_set)
  colnames(new_set) <- "Median"
  return(new_set)
}


## From: https://github.com/r-lib/fs/blob/master/R/sanitize.R
#' Sanitize a filename so that it is valid to write to.
#'
#' @param filename The filename to sanitize.
#' @param replacement The replacement string to use to replace illegal
#'   characters. If it is an illegal character itself, the output will not
#'   be sanitized.
#'
#' @returns The sanitized filename
path_sanitize <- function(filename, replacement = "") {
  illegal <- "[/\\?<>\\:*|\":]"
  control <- "[[:cntrl:]]"
  reserved <- "^[.]+$"
  windows_reserved <- "^(con|prn|aux|nul|com[0-9]|lpt[0-9])([.].*)?$"
  windows_trailing <- "[. ]+$"

  filename <- gsub(illegal, replacement, filename)
  filename <- gsub(control, replacement, filename)
  filename <- gsub(reserved, replacement, filename)
  filename <- gsub(windows_reserved, replacement, filename, ignore.case = TRUE)
  filename <- gsub(windows_trailing, replacement, filename)

  # TODO: this substr should really be unicode aware, so it doesn't chop a
  # multibyte code point in half.
  filename <- substr(filename, 1, 255)
  if (replacement == "") {
    return(filename)
  }
  path_sanitize(filename, "")
}



#' Sample data with the reservoir sampling algorithm.
#'
#' See https://en.wikipedia.org/wiki/Reservoir_sampling
#' This is an < O(n) way to sample a dataset while keeping the same probability
#' of inclusion to every entry to k/n, where k is the sample size and n is
#' the len of the data.
#'
#' @param data Either a data.frame or a vector of data to subset.
#' @param k The sample size
#'
#' @returns A sample of the same type of the input data
#'
#' @author Hedmad
reservoir_sample <- function(data, k) {
  if (is.data.frame(data)) {
    data_size <- nrow(data)
    stopifnot("Sampling size must be lower than data length"=data_size > k)
  } else {
    data_size <- length(data)
    stopifnot("Sampling size must be lower than data length"=data_size > k)
  }

  reservoir <- 1:k
  i = k
  pick_prob <- exp(log(runif(1))/1)

  while (i <= data_size) {
    i <- i + floor(log(runif(1))/log(1 - pick_prob)) + 1
    if (i <= data_size) {
      reservoir[floor(runif(1, 1, k+1))] <- i
      pick_prob <- pick_prob * exp(log(runif(1)) / k)
    }
  }
  if (is.data.frame(data)) {
    return(data[reservoir,])
  }
  return(data[reservoir])
}


#' Compute the binned mean of a vector.
#'
#' Bins the data in the vector into equally-sized bins following the order of
#' bin_by. Then calculates the mean value in each bin, and returns the resulting
#' vector. If `bin_by` is NULL, uses the data instead to calculate bins.
#'
#' May result in `NA`s if there are some bins with no numbers inside of them.
#' Passing na.rm removes them when computing the means.
#'
#' @param x A numeric vector with the data to bin
#' @param bin_by Either `NULL` or a numeric vector of the same size of `x` that
#'   will be used to bin by.
#' @param n_bins The number of bins to divide the data in.
#' @param na.rm If `TRUE`, removes
#'
#' @returns A numeric vector of length `n_bins` with the binned means.
#'
#' @author Hedmad
bin_mean <- function(x, bin_by = x, n_bins = 10, na.rm = FALSE) {
  stopifnot(
    "n_bins must be a positive integer" = {n_bins > 0 & is.wholenumber(n_bins)},
    "x cannot be of length 0" = {length(x) != 0},
    "x must be numeric" = {is.numeric(x)},
    "bin_by must be numeric" = {is.numeric(bin_by)},
    "the lengths of bin_by and x must be identical" = {length(bin_by) == length(x)},
    "n_bins must be shorter or equal to the length of bin_by" = {n_bins <= length(bin_by)}
  )

  bin_thresholds <- seq(
    from = range(bin_by)[1], to = range(bin_by)[2], length.out = n_bins + 1
  )
  last_thr <- bin_thresholds[1]
  data <- rep(NULL, n_bins)
  for (i in seq_along(bin_thresholds[-1])) {
    current_thr <- bin_thresholds[-1][i]

    target_vals <- which(bin_by >= last_thr & bin_by <= current_thr)
    data[i] <- mean(x[target_vals], na.rm = na.rm)
    last_thr <- current_thr
  }
  data[is.nan(data)] <- NA
  return(data)
}


#' Test if the input is a whole number (see `?is.integer`)
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  return(abs(x - round(x)) < tol)
}


#' This function returns a corrected alpha threshold for the Benjamini-Hochberg
#' procedure, to be used instead of the inputted one to filter by P-value
#' instead of by the adjusted p-value.
#'
#' @param adj_pvals A vector with the BH-adjusted p-values.
#' @param alpha_threshold The desired alpha threshold to correct.
#'
#' @returns The new p-value threshold
#'
#' @author Hedmad, FeAR
find_BH_critical_p <- function(adj_pvals, alpha_threshold = 0.05) {
  stopifnot(
    "length of adj_pvals must be greater than 0" = {length(adj_pvals) > 0},
    "alpha_threshold must not be negative" = {alpha_threshold >= 0}
  )
  if (any(adj_pvals < 0)) {
    warning("Found some adjusted p-values which are less than 0.")
  }
  if (any(adj_pvals == 0)) {
    warning("Found some adjusted p-values which are 0.")
  }
  if (any(adj_pvals > 1)) {
    warning("Found some adjusted p-values which are greater than 1.")
  }
  total_hits = sum(adj_pvals < alpha_threshold)
  corrected_alpha = (alpha_threshold / length(adj_pvals)) * (total_hits + 1)
  return(corrected_alpha)
}


#' Return TRUE if all items in the vector are identical, false otherwise.
#'
#' @param vector The vector to test.
#'
#' @author Hedmad
all_identical <- function(vector) {
  if (length(vector) == 1) {
    return(TRUE)
  }
  return(all(vector == vector[1]))
}

#' The kOverA function from GeneFilter. We need just this one, so loading
#' the package is overkill.
kOverA <- function (k, A = 100, na.rm = TRUE) {
  function(x) {
    if (na.rm)
      x <- x[!is.na(x)]
    sum(x > A) >= k
  }
}