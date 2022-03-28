# The logging libs around are overkill. I just need two types of files:
# A data log, for writing data, and a normal log, to write output info.
# The log needs to be always saved in '/GATTACA/target/', and always
# sent to the console.

LOGLEVELS <- list(
    "debug" = 10,
    "info" = 20,
    "warning" = 30,
    "error" = 40
)

.set_logname <- function(name) {
    # Set the log name to `name`. The path and extensions of _data.log and .log
    # are automatically set.
    options("GATTACA_log" = paste0("/GATTACA/logs/", name, ".log"))
    options("GATTACA_data_log" = paste0("/GATTACA/logs/", name, "_data.log"))
}

.set_loglevel <- function(level, target) {
    # Levels can be "debug" < "info" < "warning" < "error".
    # Only messages that match that level or above are logged.
    if (! level %in% names(LOGLEVELS)) {stop(paste("Invalid loglevel", level))}
    if (! target %in% c("stdout", "file")) {stop(paste("Invalid target", target))}

    # R FUCKING SUCKS BALLS - The paste0 function does not get evaluated, so it
    # tries to set it as the name, and fails. So i have to do this shit
    path <- paste0("GATTACA_log_level_", target)
    newopt <- list(LOGLEVELS[level])
    names(newopt) <- path
    options(newopt)
}

.set_loglevel_file <- function(level) {.set_loglevel(level, "file")}
.set_loglevel_stdout <- function(level) {.set_loglevel(level, "stdout")}

.prep_data <- function(arg_list, padding = "    ") {
    string <- ""
    for (i in seq_along(arg_list)){
        string_data <- get.print.str(arg_list[[i]])
        padded <- paste0(padding, gsub("\\n", paste0("\n", padding), string_data))
        string <- paste0(string, "[", names(arg_list)[i], "] >>>\n", padded, "\n<<<\n")
    }
    return(string)
}

.push_to_log <- function(..., level, kind = "standard") {
    args <- list(...)

    if (! kind %in% c("standard", "data")) {stop(paste0("Invalid kind ", kind))}

    slevel <- names(LOGLEVELS)[unlist(LOGLEVELS) >= level][1]

    if (kind == "standard") {
        string <- paste0(unlist(args), collapse = "")
        msg <- paste0(date(), " [", slevel, "]: ", string, "\n")
    } else {
        msg <- .prep_data(args)
    }

    if (level >= getOption("GATTACA_log_level_stdout", default = 20)){
        cat(msg)
    }

    if (level >= getOption("GATTACA_log_level_file", default = 20)) {
        if (kind == "standard") {
            file <- getOption("GATTACA_log", default = stop("Unset log name. Cannot log messages to file."))
        } else {
            file <- getOption("GATTACA_data_log", default = stop("Unset log name. Cannot log messages to file."))
        }
        cat(msg, file = file, append = TRUE)
    }
}

log <- list(
    set_name = .set_logname,
    set_console_level = .set_loglevel_stdout,
    set_file_level = .set_loglevel_file,

    debug = \(...){.push_to_log(..., level = LOGLEVELS[["debug"]])},
    info = \(...){.push_to_log(..., level = LOGLEVELS[["info"]])},
    warn = \(...){.push_to_log(..., level = LOGLEVELS[["warning"]])},
    error = \(...){.push_to_log(..., level = LOGLEVELS[["error"]])},
    data = \(...){.push_to_log(..., level = LOGLEVELS[["info"]], kind = "data")}
)